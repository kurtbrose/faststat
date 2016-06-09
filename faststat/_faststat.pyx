from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_New
from libc.math cimport log


cdef extern from "faststat_core.c":
    ctypedef struct faststat_P2Percentile:
        unsigned short percentile
        double val
        unsigned int n

    ctypedef struct faststat_Bucket:
        float max
        unsigned int count

    ctypedef struct faststat_DataPoint:
        double value
        unsigned long long nanotime

    void faststat_P2Percentile_update(
        faststat_P2Percentile *array, int num_percentiles, unsigned long long n,
        double max, double min, double x)

    void faststat_Bucket_update(faststat_Bucket *array, double x)

    void faststat_DataPoint_update_topN(
        faststat_DataPoint *topN, unsigned int num_top, double x, unsigned long long t)

    void faststat_WindowCount_update(
        faststat_WindowCount *window_counts, unsigned int num_window_counts, unsigned long long t)

    float NAN


ctypedef struct ExpoAvg:
    double val
    double alpha


cdef class Stats:
    cdef:
        unsigned long long n
        double mean, min, max, m2, m3, m4
        double sum_of_logs, sum_of_inv
        unsigned long long mintime, maxtime, lasttime
        unsigned int num_percentiles
        faststat_P2Percentile *percentiles
        unsigned int num_buckets
        faststat_Bucket *buckets
        unsigned int num_expo_avgs
        ExpoAvg *expo_avgs
        double window_sum
        unsigned int num_prev  # must be a power of 2
        faststat_DataPoint *lastN
        unsigned int num_top  # must be a power of 2
        faststat_DataPoint *topN
        unsigned int num_window_counts
        faststat_WindowCount *window_counts
        void *alloc = NULL

    def __cinit__(self, list buckets, list percentiles, Stats interval, 
            list expo_avgs, list window_counts, int num_top):
        cdef:
            size_t to_alloc = 0
            void *alloc = NULL
            int percentile_offset, buckets_offset

        # figure out size of memory block needed and allocate all at once
        percentile_offset = to_alloc
        if percentiles:
            to_alloc += round_size(sizeof(faststat_P2Percentile) * len(percentiles))
        buckets_offset = to_alloc
        if buckets:
            if buckets[-1] != float('+inf'):
                raise ValueError('last bucket must be positive infinity; found: {0}'.format(
                    buckets[-1]))
            to_alloc += round_size(sizeof(faststat_Bucket) * len(buckets))




        alloc = PyMem_Malloc(to_alloc * 4)
        if alloc == NULL:
            raise MemoryError  # could not allocate faststat memory
        self.alloc = alloc

    def add(self, double x not None, unsigned long long t not None):
        faststat_WindowCount_update(self.window_counts, self.num_window_counts, self.lasttime, t)
        self.lasttime = t
        self.n += 1
        if self.n == 1:
            self.min = self.max = x
            self.mintime = self.maxtime = self.lasttime
        if x <= self.min:
            self.mintime = self.lasttime
            self.min = x
        if x >= self.max:
            self.maxtime = self.lasttime
            self.max = x
        self.sum_of_logs += log(x) if x > 0 else NAN
        self.sum_of_inv += 1 / x if x > 0 else NAN
        self._update_moments(x)
        faststat_P2Percentile_update(self.percentiles, self.num_percentiles, self.n, self.max, self.min, x)
        faststat_Bucket_update(self.buckets)
        self._update_expo_avgs(x)
        faststat_DataPoint_update_topN(self.topN, self.num_top, x, t)
        self._update_lastN(x)

    cdef void _update_moments(self, double x):
        '''
        update mean, and 2nd, 3rd, and 4th moments
        '''
        cdef double n, delta, delta_n, delta_m2, delta_m3, delta_m4
        n = self.n  # assign to double to ensure compiler keeps all intermediate results as floats
        delta = x - self.mean
        delta_n = delta / n
        delta_m2 = delta * delta_n * (n - 1)
        delta_m3 = delta_m2 * delta_n * (n - 2)
        delta_m4 = delta_m2 * delta_n * delta_n * (n * (n - 3) + 3)
        # compute updated values
        self.mean = self.mean + delta_n
        # note: order matters here
        self.m4 += delta_m4 + delta_n * (6 * delta_n * self.m2 - 4 * self.m3)
        self.m3 += delta_m3 + delta_n * 3 * self.m2
        self.m2 += delta_m2

    cdef void _update_expo_avgs(self, double x):
        cdef:
            unsigned int i
            double val, alpha
        for i in range(self.num_expo_avgs):
            val = self.expo_avgs[i].val
            alpha = self.expo_avgs[i].alpha
            # this equation ensures no "gain"
            self.expo_avgs[i].val = x * alpha + val * (1 - alpha)

    cdef void _update_lastN(self, double x):
        cdef unsigned int offset
        if self.num_prev == 0:
            return
        # bit-wise or here is much faster than modulo, and correct
        # as long as num_prev is a power of 2
        offset = (self.n - 1) & (self.num_prev - 1)
        self.window_sum -= self.lastN[offset].value
        self.window_sum += x
        self.lastN[offset].value = x
        self.lastN[offset].nanotime = self.lasttime


    property n:
        'number of points seen'
        def __get__(self): return self.n

    property mean:
        def __get__(self): return self.mean

    property sum_of_logs:
        def __get__(self): return self.sum_of_logs

    def __dealloc__(self):
        if self.alloc:
            PyMem_Free(self.alloc)


cdef size_t round_size(size_t size):
    'round up allocation size to an even multiple of 8 bytes'
    if not size % 8:  # already an even multiple of 8
        return size
    return size + 8 - (size % 8)
