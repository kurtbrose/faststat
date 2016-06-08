from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_New


cdef extern from "faststat_core.c":
    ctypedef faststat_Stats:
        pass
    int faststat_P2Percentile_size


cdef class Stats:
    cdef:
        faststat_Stats stats_struct
        void *extra_heap = NULL

    def __cinit__(self, list buckets, list percentiles, Stats interval, 
            list expo_avgs, list window_counts, int num_top):
        cdef:
            size_t to_alloc = 0
            void *alloc = NULL
            int percentile_offset

        if percentiles:
            percentile_offset = to_alloc
            to_alloc += round_size(faststat_P2Percentile_size * len(percentiles))


        alloc = PyMem_Malloc(to_alloc * 4)
        if alloc == NULL:
            raise MemoryError  # could not allocate faststat memory
        self.extra_heap = alloc

        faststat_Stats_init(
            &(self.stats_struct),
            alloc + percentile_offset, len(percentiles),
            alloc + )


    def __dealloc__(self):
        if self.extra_heap:
            PyMem_Free(self.extra_heap)


cdef size_t round_size(size_t size):
    'round up allocation size to an even multiple of 8 bytes'
    if not size % 8:  # already an even multiple of 8
        return size
    return size + 8 - (size % 8)
