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


class PyStats(object):
    def __init__(self, buckets=(), lastN=64, percentiles=DEFAULT_PERCENTILES):
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


class PyInterval(object):
    pass  # TODO


class PyDuration(object):
    pass  # TODO:


try:
    import _faststat

    # keep buckets for intervals in size from 100ns to ~14 hours
    TIME_BUCKETS = sum( 
        [(1*10**-x, 2*10**-x, 5*10**-x) for x in (7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4)], ())

    class _BaseStats(object):
        'base class to avoid repeating code'
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

        @property
        def percentiles(self):
            return self.get_percentiles()

        @property
        def buckets(self):
            return self.get_buckets()


    class CStats(_BaseStats):
        '''
        Call add(value) to add a data point.
        '''
        def __init__(self, buckets=(), lastN=64, percentiles=DEFAULT_PERCENTILES):
            self.interval = CInterval()
            self._stats = _faststat.Stats(buckets, lastN, percentiles, self.interval._stats)
            self.add = self._stats.add

        def __getattr__(self, name):
            if name not in ("tick", "end"):
                return getattr(self._stats, name)

    class CInterval(_BaseStats):
        '''
        Call tick() to register occurrences.
        Note that calling tick() N times results in N-1 data points.
        '''
        def __init__(self, buckets=TIME_BUCKETS, lastN=64, percentiles=DEFAULT_PERCENTILES):
            self._stats = _faststat.Stats(buckets, lastN, percentiles, None)
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
            self._stats = _faststat.Stats(buckets, lastN, percentiles, self.interval._stats)
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
