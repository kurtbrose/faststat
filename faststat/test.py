'''
test implementation from wikipedia
http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics
'''
import random
import time
from faststat import Stats


def p2_parabolic(l_v, l_n, c_v, c_n, r_v, r_n, d):
    return (c_v + (d / (r_n - l_n)) * (
        (c_n - l_n + d) * (r_v - c_v) / (r_n - c_n) +
        (r_n - c_n - d) * (c_v - l_v) / (c_n - l_n)))


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


def test():
    random.seed(103)  # make test repeatable
    data = [random.normalvariate(1.0, 1.0) for i in range(int(1e6))]
    stats = Stats()
    start = time.time()
    for d in data:
        stats.add(d)
    print time.time() - start, "microseconds per point"
    '''
    print "mean (should be 1)", stats.mean
    print "kurtosis / reference kurtosis", stats.kurtosis, online_kurtosis(data)
    print "variance / reference variance", stats.variance, online_variance(data)
    print "skewness (should be 0)", stats.skewness
    print "max, min, mintime, maxtime", stats.max, stats.min, stats.mintime, stats.maxtime
    print "m2, m3, m4", stats.m2, stats.m3, stats.m4
    print "geometric_mean, harmonic mean", stats.geometric_mean, stats.harmonic_mean
    print "interval.min, interval.geometric_mean, interval.harmonic_mean",
    print stats.interval.min, stats.interval.geometric_mean, stats.interval.harmonic_mean
    print "expo_avgs (should be 1)", stats.expo_avgs
    print "window_counts", stats.get_window_counts()
    print "top 10", sorted(stats.get_topN())[-10:]
    return stats
    '''
    import pprint
    pprint.pprint(stats.get_qdigest_state())


if __name__ == "__main__":
    test()