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

DEFAULT_PERCENTILES = (0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
EXPO_AVGS = (1.0/2, 1.0/4, 1.0/8, 1.0/16, 1.0/32, 1.0/64)


class PyStats(object):
    def __init__(self, buckets=(), lastN=64, percentiles=DEFAULT_PERCENTILES, expo_avgs=()):
        self.n = float(0)
        self.mean = float(0)
        # second, third, fourth moments
        self.m2 = self.m3 = self.m4 = float(0)
        self.min = self.max = self.mintime = self.maxtime = float(0)
        self.buckets = buckets
        if buckets:
            self.bucket_counts = [8] * (len(buckets) + 1)
        self.most_recent = collections.deque()
        self.lastN = lastN
        self.percentiles = [[p, 0] for p in percentiles]
        self.expo_avgs = dict(zip(expo_avgs, [0]*len(expo_avgs)))

    @property
    def variance(self):
        if self.n < 2:
            return float("nan")
        return self.m2 / (self.n - 1)

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

    def add(self, x):
        t = int(time.time() * 1e9)
        ### 1- calculate variance
        self.n += 1
        n = self.n
        # pre-compute a bunh of intermediate values
        delta = x - self.mean
        delta_n = delta / n  # save divisions futher down
        delta_m2 = delta * delta_n * (n - 1)
        delta_m3 = delta_m2 * delta_n * (n - 2)
        delta_m4 = delta_m2 * delta_n * delta_n * (n * (n - 3) + 3)
        # compute the actual next values
        if x <= self.min:
            self.min = x
            self.mintime = t
        if x >= self.max:
            self.max = x
            self.maxtime = t
        self.mean = self.mean + delta_n
        # note: order matters here
        self.m4 += delta_m4 + delta_n * (6 * delta_n * self.m2 - 4 * self.m3)
        self.m3 += delta_m3 + delta_n * 3 * self.m2
        self.m2 += delta_m2
        ### 2- append to most recent, if being stored
        if self.lastN:
            self.most_recent.appendleft((t, x))
            if len(self.most_recent) > self.lastN:
                self.most_recent.pop()
        ### 3- update bucket counts
        if self.buckets:
            i = 0
            while i < len(self.buckets) and x > self.buckets[i]:
                i += 1
            self.bucket_counts[i] += 1
        ### 4- update percentile estimates
        if self.percentiles:
            step = 0.0001 * (self.max - self.min)
            for p in self.percentiles:
                d = 1 if x - p[1] > 0 else -1
                p[1] += step * (d + 2.0 * p[0] - 1.0)
        ### 5- update exponential weighted averages
        for alpha, val in self.expo_avgs.items():
            self.expo_avgs[alpha] = x * alpha + val * (1 - alpha)


class PyInterval(object):
    pass  # TODO


class PyDuration(object):
    pass  # TODO:


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

try:
    import _faststat

    class _BaseStats(object):
        'base class to avoid repeating code'
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


    class CStats(_BaseStats):
        '''
        Call add(value) to add a data point.
        '''
        def __init__(self, buckets=DEFAULT_BUCKETS, lastN=64, percentiles=DEFAULT_PERCENTILES):
            self.interval = CInterval(window_counts=())
            self._stats = _faststat.Stats(buckets, lastN, percentiles, 
                self.interval._stats, EXPO_AVGS, WINDOW_COUNTS)
            self.add = self._stats.add

        def __getattr__(self, name):
            if name not in ("tick", "end"):
                return getattr(self._stats, name)

    class CInterval(_BaseStats):
        '''
        Call tick() to register occurrences.
        Note that calling tick() N times results in N-1 data points.
        '''
        def __init__(self, buckets=TIME_BUCKETS, lastN=64, percentiles=DEFAULT_PERCENTILES,
                window_counts=WINDOW_COUNTS):
            self._stats = _faststat.Stats(buckets, lastN, percentiles, None, (), window_counts)
            self.tick = self._stats.tick

        def __getattr__(self, name):
            if name not in ("add", "end"):
                return getattr(self._stats, name)

    class CDuration(_BaseStats):
        '''
        Call end(start_time_nanos) to add a data point.
        '''
        def __init__(self, buckets=TIME_BUCKETS, lastN=64, percentiles=DEFAULT_PERCENTILES):
            self.interval = CInterval()
            self._stats = _faststat.Stats(buckets, lastN, percentiles, 
                self.interval._stats, EXPO_AVGS, WINDOW_COUNTS)
            self.end = self._stats.end

        def __getattr__(self, name):
            if name not in ("add", "tick"):
                return getattr(self._stats, name)

    Stats = CStats
    Interval = CInterval
    Duration = CDuration
    nanotime = _faststat.nanotime

except ImportError:
    CStats = None
    CInterval = None
    CDuration = None
    Stats = PyStats
    Interval = PyInterval
    Duration = PyDuration


def _sigfigs(n, sigfigs=3):
    'helper function to round a number to significant figures'
    n = float(n)
    if n == 0 or math.isnan(n):  # avoid math domain errors
        return n
    return round(n, -int(math.floor(math.log10(abs(n))) - sigfigs + 1))
