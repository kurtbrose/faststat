'''
Computation of higher order statistics based on
http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics

Intended for use in long-running applications.
'''

try:
    import _faststat

    class Stats(object):
        def __init__(self):
            self._stats = _faststat.Stats()

        @property
        def variance(self):
            return self.m2 / (self.n - 1)

        @property
        def skewness(self):
            return self.n ** 0.5 * self.m3 / self.m2 ** 1.5

        @property
        def kurtosis(self):
            return self.n * self.m4 / self.m2 ** 2 - 3

        def __getattr__(self, name):
            return getattr(self._stats, name)

except ImportError:
    class Stats(object):
        def __init__(self):
            self.n = float(0)
            self.mean = float(0)
            # second, third, fourth moments
            self.m2 = self.m3 = self.m4 = float(0)
            self.min = self.max = float(0)

            @property
            def variance(self):
                return self.m2 / (self.n - 1)

            @property
            def skewness(self):
                return self.n ** 0.5 * self.m3 / self.m2 ** 1.5

            @property
            def kurtosis(self):
                return self.n * self.m4 / self.m2 ** 2 - 3

            def add(self, x):
                ### calculate variance
                self.n += 1
                n = self.n
                # pre-compute a bunh of intermediate values
                delta = x - self.mean
                delta_n = delta / n  # save divisions futher down
                delta_m2 = delta * delta_n * (n - 1)
                delta_m3 = delta_m2 * delta_n * (n - 2)
                delta_m4 = delta_m2 * delta_n * delta_n * (n * n - 3 * n + 3)
                # compute the actual next values
                self.min = x if x < self.min else self.min
                self.max = x if x > self.max else self.max
                self.mean = self.mean + delta_n
                # note: order matters here
                self.m4 += delta_m4 + delta_n * (6 * delta_n * self.m2 - 4 * self.m3)
                self.m3 += delta_m3 + delta_n * 3 * self.m2
                self.m2 += delta_m2


def test():
    import random
    import time
    random.seed(10)  # make test repeatable
    data = [random.normalvariate(1.0, 1.0) for i in range(int(1e6))]
    stats = Stats()
    start = time.time()
    for d in data:
        stats.add(d)
    print time.time() - start, "microseconds per point"
    print "mean (should be 1)", stats.mean
    print "kurtosis (should be 0)", stats.kurtosis
    print "variance (should be 1)", stats.variance
    print "skewness (should be 0)", stats.skewness

if __name__ == "__main__":
    test()
