#include <stdio.h>
#include <math.h>

//define gettime() which returns integral nanoseconds since epoch
//as a 64 bit integer for a variety of platforms
#ifdef _WIN32
//windows has its own brand of time keeping functions
#include <windows.h>
#define DELTA_EPOCH_IN_SECS  11644473600ULL
//difference between Jan 1, 1601 and Jan 1, 1970 (unix epoch)

static unsigned long long _nanotime(void) {
    FILETIME ft;
    ULARGE_INTEGER result;
    GetSystemTimeAsFileTime(&ft); //returns time in 100ns intervals since Jan 1, 1601
    result.HighPart = ft.dwHighDateTime;
    result.LowPart = ft.dwLowDateTime;
    result.QuadPart -= DELTA_EPOCH_IN_SECS * 10000000ULL; // 1000 (ms) * 1000 (us) * 10 (100ns)
    return result.QuadPart * 100;
}

// for old versions of MSVC which do not include NAN macro
#ifndef NAN
    static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
    #define NAN (*(const float *) __nan)
#endif

#elif defined linux
//linux has clock_gettime(CLOCK_REALTIME) which is ns since epoch -- perfect
#include <time.h>

static unsigned long long _nanotime(void) {
    struct timespec ts;
    if(clock_gettime(CLOCK_REALTIME, &ts) == -1) {
        return 0;
    }
    return 1000000000ULL * ts.tv_sec + ts.tv_nsec;
}

#else
//for those oddballs like OSX and BSD, fall back on gettimeofday() which is at least microseconds
#include <time.h>

static unsigned long long _nanotime(void) {
    struct timeval tv;
    if(gettimeofday(&tv, NULL) == -1) {
        return 0;
    }
    return tv.tv_sec * 1000000000ULL + tv.tv_usec * 1000ULL;
}

#endif

static unsigned long long nanotime_override = 0;

static unsigned long long nanotime(void) {
    if(nanotime_override) {
        return nanotime_override;
    }
    return _nanotime();
}

//percentile point for usage in P2 algorithm
typedef struct {
    unsigned short percentile;  //divide by 0xFFFF to get a float between 0 and 1
    double val;  //estimate of current percentile value
    unsigned int n;  //estimate of how many values were less than this
} faststat_P2Percentile;


typedef struct {
    float max;
    unsigned int count;
} faststat_Bucket;


// keeping this size a power of 2 makes pointer arithmetic in heap more efficient
typedef struct {
    double value;
    unsigned long long nanotime;
} faststat_DataPoint;

// represents an exponential moving average
// aka a low pass filter, infinite impulse response filter
typedef struct {
    double val;
    double alpha;
} faststat_ExpoAvg;


// represents a count for a given interval,
// aligned on unix epoch
typedef struct {
    unsigned short num_windows;  // number of counts -- MUST BE A POWER OF 2
    unsigned long long window_size_nanosecs;  // size of each window in seconds
    unsigned int *counts;  // counts for the previous num_windows intervals
} faststat_WindowCount;


// for representing a normally distributed variable
typedef struct faststat_Stats_struct {
    unsigned long long n;
    double mean, min, max, m2, m3, m4;
    double sum_of_logs, sum_of_inv;  // for geometric and harmonic mean
    unsigned long long mintime, maxtime, lasttime;
    unsigned int num_percentiles;
    faststat_P2Percentile *percentiles;
    unsigned int num_buckets;
    faststat_Bucket *buckets; // last bucket MUST BE +inf
    unsigned int num_expo_avgs;
    faststat_ExpoAvg *expo_avgs;
    double window_avg;
    unsigned int num_prev; // MUST BE A POWER OF 2!
    faststat_DataPoint *lastN;
    unsigned int num_top; // MUST BE A POWER OF 2!
    faststat_DataPoint *topN;
    unsigned int num_window_counts;
    //window counts must be sorted by window_size, to
    //make handling code cleaner/smaller
    faststat_WindowCount *window_counts;
    struct faststat_Stats_struct *interval;
} faststat_Stats;


int faststat_Stats_init(faststat_Stats *self) {
    *self = {0};
}


// round up to the nearest power of 2
// http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
// this is faster than BSR instruction on most x86
unsigned int round_up_pwr_2(unsigned int v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;

}


#define faststat_P2Percentile_size sizeof(faststat_P2Percentile)
