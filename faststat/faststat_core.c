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


// be careful; if-condition must properly terminate when max == +inf, even for nan and +inf
#define OFFSET(n)  if(!(x >= array[i+n].max)) { array[i+n].count++; break; } 

void faststat_Bucket_update(faststat_Bucket *array, double x) {
    unsigned int i;
    for(i=0; ; i += 16) {
        OFFSET( 0) OFFSET( 1) OFFSET( 2) OFFSET( 3) 
        OFFSET( 4) OFFSET( 5) OFFSET( 6) OFFSET( 7) 
        OFFSET( 8) OFFSET( 9) OFFSET(10) OFFSET(11) 
        OFFSET(12) OFFSET(13) OFFSET(14) OFFSET(15)
    }
}


//helper for _update_percentiles
static void _p2_update_point(double l_v, double l_n, faststat_P2Percentile *cur,
                            double r_v, double r_n, unsigned long long n) {
    int d;
    double percentile, new_val, c_v, c_n, diff;
    percentile = ((double)cur->percentile) / 0x10000;
    c_n = cur->n;
    diff = (n - 1) * percentile + 1 - c_n;
    // clamp d at +/- 1
    if(diff >= 1) {
        d = 1;
    } else if(diff <= -1) {
        d = -1;
    } else {
        return;
    }
    c_v = cur->val;
    if(l_n < c_n + d && c_n + d < r_n) {  // try updating estimate with parabolic
        new_val = c_v + (d / (r_n - l_n)) * ( 
            (c_n - l_n + d) * (r_v - c_v) / (r_n - c_n) +
            (r_n - c_n - d) * (c_v - l_v) / (c_n - l_n));
        if(l_v >= new_val || r_v <= new_val) {  // fall back on linear
            if(d == 1) {
                new_val = c_v + (r_v - c_v) / (r_n - c_n);
            } else {  // d == -1
                new_val = c_v - (l_v - c_v) / (l_n - c_n);
            }
        }
        cur->val = new_val;
        cur->n += d;
    }
}


static void _insert_percentile_sorted(
        faststat_P2Percentile *array, int num_percentiles,
        unsigned long long n, double x) {
    unsigned long long num, i;  // prevent loss of precision compiler warning
    double tmp;
    num = n < num_percentiles ? n : num_percentiles;
    for(i = 0; i < num-1; i++) { // insert in sorted order
        if(x < percentiles[i].val) {
            tmp = x;
            x = percentiles[i].val;
            percentiles[i].val = tmp;
        }
    }
    percentiles[num-1].val = x;
}


//update percentiles using piece-wise parametric algorithm
void faststat_P2Percentile_update(
        faststat_P2Percentile *array, int num_percentiles, unsigned long long n,
        double max, double min, double x) {
    unsigned int i;
    faststat_P2Percentile *right, *left, *cur, *prev, *nxt;
    right = &(array[num_percentiles - 1]);
    left = array;
    if(!(right->n < n) ) { // just insert until n > num_percentiles
        _insert_percentile_sorted(array, num_percentiles, n, x);
        return;
    }
    //right-most is stopping case; handle first
    if(x < right->val && right->n + 1 < n) {
        right->n++;
    }
    //handle the rest of the points
    prev = right;
    for(i = num_percentiles-2; ; i--) {
        cur = &(percentiles[i]);
        if(x < cur->val && cur->n + 1 < prev->n) {
            cur->n++;
        }
        prev = cur;
        if(i == 0) { //making i unsigned fixes some warnings
            break;
        }
    }
    //left-most point is a special case
    nxt = &(percentiles[1]);
    _p2_update_point(min, 0, left, nxt->val, nxt->n, n);
    cur = left;
    for(i=1; i < num_percentiles - 1; i++) {
        prev = cur;
        cur = nxt;
        nxt = &(percentiles[i+1]);
        _p2_update_point(prev->val, prev->n, cur, nxt->val, nxt->n, n);
    }
    _p2_update_point(cur->val, cur->n, right, (double)max, (double)n, n);
} 

//this algorithm deterministically favors storing newer values over older values
// expects topN to be a 1-based index to make 
void faststat_DataPoint_update_topN(
        faststat_DataPoint *topN, unsigned int num_top, double x, unsigned long long t) {
    faststat_DataPoint *cur, *left, *right, *end, *min, *topN;
    unsigned int cur_i, left_i, right_i;
    if(x < topN->value) {
        return;
    }
    topN -= 1; // switch to 1-based indexing to save additions navigating down heap
    // (left child = 2 * index instead of 2 * index + 1; profiling found this to make a difference)
    // replace the smallest element of the topN with the new point
    topN[1].value = x;
    topN[1].nanotime = t;
    // restore the heap condition
    cur = topN + 1;  // use pointers instead of array indices
    cur_i = 1;
    left = topN + 2;
    left_i = 2;
    right = topN + 3;
    right_i = 3;
    end = topN + 1 + self->num_top;
    while(right < end) {
        if(left->value == right->value) {
            min = left->nanotime > right->nanotime ? left : right;
        } else {
            min = left->value < right->value ? left : right;
        }
        if(cur->value < min->value) {
            break;
        }
        // swap cur with min of left, right
        x = min->value;
        t = min->nanotime;
        min->value = cur->value;
        min->nanotime = cur->nanotime;
        cur->value = x;
        cur->nanotime = t;
        // set up for the next layer of the heap
        cur = min;
        cur_i = min == left ? left_i : right_i;
        left_i = cur_i * 2;
        right_i = left_i + 1;
        left = topN + left_i;
        right = topN + right_i;
    }
}


// re-zero all of the windows which have been "missed" between self->lasttime and t
static void _rezero_window_counts(
        faststat_WindowCount *window_counts, unsigned int num_window_counts,
        unsigned long long lasttime, unsigned long long t)) {
    faststat_WindowCount *cur;
    unsigned int i;
    unsigned long long j, last_window, cur_window;
    for(i = 0; i < num_window_counts; i++) {
        cur = &(window_counts[i]);
        last_window = lasttime / cur->window_size_nanosecs;
        cur_window = t / cur->window_size_nanosecs;
        if(last_window == cur_window) {
            break;  // because window_counts are sorted by window_size_nanosecs,
        }           // the next window cannot miss unless the current one does
        if(cur_window - last_window >= cur->num_windows) {
            memset(cur->counts, 0, sizeof(*cur->counts) * cur->num_windows);
            continue;  // if the entire array is getting zero'd, just use memset
        }
        // TODO: convert this to a memset instead of a loop (perhaps)
        for(j = last_window + 1; j <= cur_window; j++) {
            cur->counts[j & (cur->num_windows - 1)] = 0;
        }  // zero out all the "missed" windows
    }
}


void faststat_WindowCount_update(
        faststat_WindowCount *window_counts, unsigned int num_window_counts,
        unsigned long long lasttime, unsigned long long t) {
    faststat_WindowCount *cur;
    unsigned int i, cur_count;
    _rezero_window_counts(window_countx, num_window_counts, lasttime, t);
    // step 2 -- increment current counts
    for(i = 0; i < num_window_counts; i++) {
        cur = &(window_counts[i]);
        // use the current time as the index into the circular array to save some memory
        cur_count = (t / cur->window_size_nanosecs) & (cur->num_windows - 1);
        ++(cur->counts[cur_count]);
    }
}
