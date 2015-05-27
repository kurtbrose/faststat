'''
faststat is a *streaming*, *light-weight* statistics library designed for embedding
in other Python applications.  *Steaming* means
that faststat operates on data points as they arrive, without needing to store
previous data points.  *Light-weight* means that statistics do not take up a great
deal of CPU or RAM.  Adding a data point to a stat object is a 0.5 - 3 microsecond operation.  
Each stat object takes about 4kiB of memory.
'''
import array
import random
import math
import collections
import time
import functools
import json
import weakref

import _faststat
import cache
import format


class Sample(object):
    '''
    This class implements Reservoir Sampling to keep a random sample of an infinite stream.
    See http://gregable.com/2007/10/reservoir-sampling.html for one description.

    This class is kept separate from the other stats, because its memory usage is far greater.
    '''
    def __init__(self, sample_size=2**14, type='f'):
        self.sample = array.array(type)
        self.sample_size = sample_size
        self.num_vals = 0

    def add(self, val):
        if self.num_vals < self.sample_size:
            self.sample.append(val)
        else:
            pos = random.randint(0, self.num_vals)
            if pos < self.sample_size:
                self.sample[pos] = val
        self.num_vals += 1

DEFAULT_PERCENTILES = (0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
EXPO_AVGS = (1.0/2, 1.0/4, 1.0/8, 1.0/16, 1.0/32, 1.0/64)


# keep buckets for intervals in size from 100ns to ~14 hours
TIME_BUCKETS = sum( 
    [(1*10**x, 2*10**x, 5*10**x) for x in range(2, 13)], ())
# useful buckets for unsigned integers up to 64 bits
UINT_BUCKETS = (1, 2, 3, 4, 5, 6, 7, 8, 9) + sum(
    [(1*10**x, 2*10**x, 5*10**x) for x in range(1, 20)], ())
# useful buckets for signed integers up to 64 bits
INT_BUCKETS = tuple(reversed([-e for e in UINT_BUCKETS[:-3]])) + (0,) + UINT_BUCKETS[:-3]

DEFAULT_BUCKETS = (0, 1e-5, 1e-4, 1e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
    10, 20, 50, 100, 200, 500, 1e3, 2e3, 5e3, 1e4, 1e5, 1e6)
DEFAULT_BUCKETS = tuple(reversed([-e for e in DEFAULT_BUCKETS])) + (0,) + DEFAULT_BUCKETS

ONE_MIN_NS = int(60e9)
ONE_HOUR_NS = 60 * ONE_MIN_NS

WINDOW_COUNTS = [(64, ONE_MIN_NS), (32, ONE_HOUR_NS)]

class _BaseStats(object):
    'base class to avoid repeating code'
    def __init__(self, buckets, lastN, percentiles, interval, expo_avgs,
            window_counts, num_top):
        buckets = buckets + (float("inf"),)
        lastN = int(2**math.ceil(math.log(lastN)/math.log(2)))
        num_top = int(2**math.ceil(math.log(lastN)/math.log(2)))
        if interval:
            self.interval = Interval(window_counts=())
            interval = self.interval._stats
        else:
            interval = None
        self._stats = _faststat.Stats(buckets, lastN, percentiles, interval,
            expo_avgs, window_counts, num_top)

    @property
    def variance(self):
        if self.n < 2:
            return float("nan")
        return self.m2 / (self.n - 1)

    @property
    def trimean(self):
        p = self.get_percentiles()
        return (p[0.25] + 2 * p[0.5] + p[0.75]) / 4

    @property
    def skewness(self):
        if not self.m2:
            return float("nan")
        return self.n ** 0.5 * self.m3 / self.m2 ** 1.5

    @property
    def kurtosis(self):
        if not self.m2:
            return float("nan")
        return self.n * self.m4 / self.m2 ** 2 - 3

    @property
    def percentiles(self):
        return self.get_percentiles()

    @property
    def buckets(self):
        return self.get_buckets()

    @property
    def geometric_mean(self):
        'nth root of product of data points'
        return math.exp(self.sum_of_logs / self.n)

    @property
    def harmonic_mean(self):
        'inverse of mean of inverses of data points'
        return self.n / self.sum_of_inv

    @property
    def window_median(self):
        if self.num_prev:
            prev = sorted([self.get_prev(i)[1] for i in range(self.num_prev)])
            return prev[len(prev)/2]

    @property
    def expo_avgs(self):
        return self.get_expo_avgs()

    @property
    def lag_avgs(self):
        '''
        same data as expo_avgs, but with keys as the average age
        of the data -- assuming evenly spaced data points -- rather
        than decay rates
        '''
        if not self.interval:
            return
        interval = self.interval.mean
        return dict([(interval/alpha, val) 
            for alpha, val in self.get_expo_avgs().items()])

    def __repr__(self):
        p = self.percentiles
        if self.n < len(p):
            quartiles = "(n too small)"
        else:
            quartiles = (_sigfigs(p.get(0.25, -1)), 
                _sigfigs(p.get(0.5, -1)), _sigfigs(p.get(0.75, -1)))
        return '<faststat.{0} n={1} mean={2} quartiles={3}>'.format(
            type(self).__name__, self.n, _sigfigs(self.mean), quartiles)


class Stats(_BaseStats):
    '''
    Call add(value) to add a data point.
    '''
    def __init__(self, buckets=DEFAULT_BUCKETS, lastN=64, percentiles=DEFAULT_PERCENTILES,
            interval=True):
        super(Stats, self).__init__(buckets, lastN, percentiles, interval, EXPO_AVGS,
            WINDOW_COUNTS, num_top=64)
        self.add = self._stats.add

    def __getattr__(self, name):
        if name not in ("tick", "end"):
            return getattr(self._stats, name)


class _TimeStats(_BaseStats):
    def get_percentiles(self):
        data = self._stats.get_percentiles()
        for k in data:
            data[k] = ntd_float(data[k])
        return data

    @property
    def mean(self):
        return ntd_float(self._stats.mean)

    @property
    def max(self):
        return ntd_float(self._stats.max)

    @property
    def min(self):
        return ntd_float(self._stats.min)


class Interval(_TimeStats):
    '''
    Call tick() to register occurrences.
    Note that calling tick() N times results in N-1 data points.
    '''
    def __init__(self, buckets=TIME_BUCKETS, lastN=64, percentiles=DEFAULT_PERCENTILES,
            window_counts=WINDOW_COUNTS):
        super(Interval, self).__init__(buckets, lastN, percentiles, False, (), (), num_top=64)
        self.tick = self._stats.tick

    def __getattr__(self, name):
        if name not in ("add", "end"):
            return getattr(self._stats, name)


class Duration(_TimeStats):
    '''
    Represents statistics for a duration.
    Call end(start_time_nanos) to add a data point.
    '''
    def __init__(self, buckets=TIME_BUCKETS, lastN=64, percentiles=DEFAULT_PERCENTILES,
            interval=True):
        super(Duration, self).__init__(buckets, lastN, percentiles, interval, EXPO_AVGS,
            WINDOW_COUNTS, num_top=64)
        self.end = self._stats.end

    def __getattr__(self, name):
        if name not in ("add", "tick"):
            return getattr(self._stats, name)


class Markov(object):
    '''
    Represents the states of a Markov process.  The transition between states are
    modeled as Intervals, and the time spent in a given state is modeled as
    Durations.
    '''
    def __init__(self):
        self.state_durations = collections.defaultdict(Duration)
        self.transition_intervals = collections.defaultdict(Interval)
        self.transitor_states = collections.defaultdict(int)
        self._weakref_holder = {}
        self.state_counts = collections.defaultdict(functools.partial(Stats, interval=False))

    def _transition(self, nxt, cur=None, since=None):
        '''
        Register that a transition has taken place.
        nxt is an identifier for the state being entered.
        cur is an identifier for the state being left.
        since is the time at which the previous state was entered.
        '''
        self.transition_intervals[(cur, nxt)].tick()
        if since:
            self.state_durations[cur].end(since)

    def make_transitor(self, state):
        '''
        Creates and returns a new Markov.Transitor, in the passed state.
        '''
        return Markov.Transitor(self, state)

    def _cleanup(self, ref):
        'cleanup after a transitor weakref fires'
        self.transitor_states[self._weakref_holder[ref]] -= 1
        del self._weakref_holder[ref]

    class Transitor(object):
        '''
        An extremely light-weight object that simply tracks a current
        state and the time of the last transition.
        '''
        def __init__(self, markov, state):
            self.markov = markov
            self.state = state
            self.markov._transition(state)
            self.markov.transitor_states[state] += 1
            state_count = self.markov.transitor_states[state]
            self.markov.state_counts[state].add(state_count)
            self.weakref = weakref.ref(self, markov._cleanup)
            self.markov._weakref_holder[self.weakref] = state
            self.last_transition = nanotime()

        def transition(self, state):
            '''
            Notify the parent Markov stats object of a transition
            from the current state to the passed state.
            '''
            self.markov.transitor_states[self.state] -= 1
            self.markov.transitor_states[state] += 1
            old_state_count = self.markov.transitor_states[self.state]
            new_state_count = self.markov.transitor_states[state]
            self.markov.state_counts[self.state].add(old_state_count)
            self.markov.state_counts[state].add(new_state_count)
            self.markov._weakref_holder[self.weakref] = state
            self.markov._transition(state, self.state, self.last_transition)
            self.last_transition, self.state = nanotime(), state

        def __repr__(self):
            return '<faststat.Markov.Transitor({0})>'.format(self.state)


class PathStats(object):
    '''
    Represents a set of paths taken.  Unlike Markov, these states remember
    their "history".

    Because there are likely to be many more paths than states, the Path
    is more aggressively memory optimized, and also employs a SegmentedCache
    internally to limit the number of unique path durations that will be
    stored.
    '''
    def __init__(self, maxsize=2048):
        self.path_stats = cache.SegmentedCache(maxsize)
        self._weakref_path_map = {}

    def make_walker(self, start):
        'returns a walker object that tracks a path'
        return PathStats.Walker(self, start)

    def _commit(self, ref):
        'commit a walkers data after it is collected'
        path_times = self._weakref_path_map[ref]
        path_times.append(nanotime())
        del self._weakref_path_map[ref]
        path = tuple(path_times[1::2])
        times = path_times[::2]
        if path not in self.path_stats:
            # tuple to save a tiny bit of memory
            self.path_stats[path] = tuple([
                Duration(interval=False) for i in range(len(path))])
        path_stats = self.path_stats[path]
        for i in range(1, len(times)):
            path_stats[i - 1]._stats.add(times[i] - times[i - 1])

    def pformat(self, prefix=()):
        '''
        Makes a pretty ASCII format of the data, suitable for
        displaying in a console or saving to a text file.
        Returns a list of lines.
        '''
        nan = float("nan")

        def sformat(segment, stat):
            FMT = "n={0}, mean={1}, p50/95={2}/{3}, max={4}"
            line_segs = [segment]
            for s in [stat]:
                p = s.get_percentiles()
                p50, p95 = p.get(0.50, nan), p.get(0.95, nan)
                line_segs.append(FMT.format(s.n, s.mean, p50, p95, s.max))
            return '{0}: {1}'.format(*line_segs)

        lines = []
        for path in sorted(self.path_stats.keys()):
            lines.append('=====================')
            for seg, stat in zip(path, self.path_stats[path]):
                lines.append(sformat(seg, stat))
        return lines

    def _finished_segment(self, path, since, start):
        if path not in self.state_stats:
            self.state_stats[path] = (
                Duration(interval=False), Duration(interval=False))
        dur, offset = self.state_stats[path]
        dur.end(since)
        offset.end(start)

    class Walker(object):
        '''
        A light-weight object that tracks a current path and the time of
        the last transition.  Similar to Tranistor for Markov.
        '''
        def __init__(self, pathstats, segment="NEW"):
            self.pathstats = pathstats
            self._commiter = weakref.ref(self, self.pathstats._commit)
            self.path = self.pathstats._weakref_path_map[self._commiter] = []
            self.push(segment)

        def push(self, segment):
            '''
            pushes a new segment onto the path, closing out the previous segment
            '''
            self.path.append(nanotime())
            self.path.append(segment)
            self.curseg = segment

        def pop(self):
            self.push(PathStats.POP)

        def branch(self):
            child = PathStats.Walker(
                self.pathstats, PathStats.BRANCH_C)
            self.push(PathStats.BRANCH_P)
            return child

        def join(self, walker):
            self.push((PathStats.JOIN, tuple(walker.path)))
            walker.push(PathStats.JOINED)

    BRANCH_P, BRANCH_C, JOIN, JOINED, POP = "BRANCH_P", "BRANCH_C", "JOIN", "JOINED", "POP"

    def __repr__(self):
        return "<PathStats npaths={0}>".format(len(self.state_stats))


class ntd_float(float):
    'a float which represents the difference of two timestamps in nanoseconds'
    def __repr__(self):
        return format.si_format(self / 1e9, "s")

    def __format__(self, format_spec):
        return repr(self)


TimeSeries = functools.partial(Stats, interval=False)
nanotime = _faststat.nanotime


def _sigfigs(n, sigfigs=3):
    'helper function to round a number to significant figures'
    n = float(n)
    if n == 0 or math.isnan(n):  # avoid math domain errors
        return n
    return round(n, -int(math.floor(math.log10(abs(n))) - sigfigs + 1))


def merge_moments(m_a, m_a2, m_a3, m_a4, n_a, m_b, m_b2, m_b3, m_b4, n_b):
    '''
    Merge moments of two samples A and B.
    parameters are 
    m_a, ..., m_a4 = first through fourth moment of sample A
    n_a = size of sample A
    m_b, ..., m_b4 = first through fourth moment of sample B
    n_b = size of sample B
    '''
    delta = m_b - m_a
    delta_2 = delta * delta
    delta_3 = delta * delta_2
    delta_4 = delta * delta_3
    n_x = n_a + n_b
    m_x = m_a + delta * n_b / n_x
    m_x2 = m_a2 + m_b2 + delta_2 * n_a * n_b / n_x
    m_x3 = m_a3 + m_b3 + delta_3 * n_a * n_b * (n_a - n_b) + 3 * delta * (n_a * m_2b - n_b * m_2a) / n_x
    m_x4 = (m_a4 + m_b4 + delta_4 * (n_a * n_b * (n_a * n_a - n_a * n_b + n_b * n_b)) / (n_x ** 3) +
            6 * delta_2 * (n_a * n_a * m_b2 + n_b * n_b * m_a2) / (n_x ** 2) +
            4 * delta * (n_a * m_b3 - n_b * m_a3) / n_x )
    return m_x, m_x2, m_x3, m_x4, n_x
