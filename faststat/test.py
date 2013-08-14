'''
test implementation from wikipedia
http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics
'''
import random
import time
from faststat import Stats, PyStats


def online_variance(data):
    n = 0
    mean = 0
    M2 = 0

    for x in data:
        n = n + 1
        delta = x - mean
        mean = mean + delta/n
        M2 = M2 + delta*(x - mean)

    variance = M2/(n - 1)
    return variance


def online_kurtosis(data):
    n = 0
    mean = 0
    M2 = 0
    M3 = 0
    M4 = 0

    for x in data:
        n1 = n
        n = n + 1
        delta = x - mean
        delta_n = delta / n
        delta_n2 = delta_n * delta_n
        term1 = delta * delta_n * n1
        mean = mean + delta_n
        M4 = M4 + term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3
        M3 = M3 + term1 * delta_n * (n - 2) - 3 * delta_n * M2
        M2 = M2 + term1

    kurtosis = (n*M4) / (M2*M2) - 3
    return kurtosis


def test(py=False):
    random.seed(10)  # make test repeatable
    data = [random.normalvariate(1.0, 1.0) for i in range(int(1e6))]
    if py:
        stats = PyStats()
    else:
        stats = Stats()
    start = time.time()
    for d in data:
        stats.add(d)
    print time.time() - start, "microseconds per point"
    print "mean (should be 1)", stats.mean
    print "kurtosis / reference kurtosis", stats.kurtosis, online_kurtosis(data)
    print "variance / reference variance", stats.variance, online_variance(data)
    print "skewness (should be 0)", stats.skewness
    print "max, min", stats.max, stats.min
    print "m2, m3, m4", stats.m2, stats.m3, stats.m4

if __name__ == "__main__":
    test()
    test(True)
