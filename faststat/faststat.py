'''
Intended as a light-weight general purpose stats collector.

Has several functions:
1-
Computation of higher order statistics based on
http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics

2-
Keeping a list of the last several entries seen.

3-
Keeping a list of bucket counts of amount seen less than or equal to each amount.

These functions are all rolled into a single class, because at this low level
the function call overhead is significant.

Also a C-implementation is provided which improves the performance by a factor
of ~15.

Intended for use in long-running applications, where there will be many
more values added than stats queries made.
'''
import array
import random
import math
import collections
import time
import functools

import _faststat


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


class Interval(_BaseStats):
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


class Duration(_BaseStats):
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
        self.state_durations  = collections.defaultdict(faststat.Duration)
        self.transition_intervals = collections.defaultdict(faststat.Interval)

    def transition(self, nxt, cur=None, since=None):
        '''
        Register that a transition has taken place.
        nxt is an identifier for the state being entered.
        cur is an identifier for the state being left.
        since is the time at which the previous state was entered.
        '''
        self.transition_intervals[(cur, nxt)].tick()
        if since:
            self.state_durations[cur].end(since)


TimeSeries = functools.partial(Stats, interval=False)


nanotime = _faststat.nanotime


def _sigfigs(n, sigfigs=3):
    'helper function to round a number to significant figures'
    n = float(n)
    if n == 0 or math.isnan(n):  # avoid math domain errors
        return n
    return round(n, -int(math.floor(math.log10(abs(n))) - sigfigs + 1))
