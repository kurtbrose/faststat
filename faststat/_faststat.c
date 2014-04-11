#include <Python.h>
#include <structmember.h>
#include <pymem.h>
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
    PyObject_HEAD
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

/*
typedef struct {
    unsigned int n;
    unsigned int num_prev;

} faststat_StatsGroup;
*/

char* NEW_ARGS[] = {"buckets", "lastN", "percentiles", "interval", "expo_avgs", 
    "window_counts", "num_top", NULL};


static PyObject* faststat_Stats_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    faststat_Stats *self;
    PyObject *buckets, *percentiles, *interval, *expo_avgs, *window_counts, *cur;
    int num_prev, num_buckets, num_percentiles, num_expo_avgs, num_window_counts, num_top;
    int i, total, offset;
    double temp;
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "OiOOOOi", NEW_ARGS, 
            &buckets, &num_prev, &percentiles, &interval, &expo_avgs, &window_counts, &num_top)) {
        return NULL;
    }

    buckets = PySequence_Fast(buckets, "expected a sequence");
    percentiles = PySequence_Fast(percentiles, "expected a sequence");
    expo_avgs = PySequence_Fast(expo_avgs, "expected a sequence");
    window_counts = PySequence_Fast(window_counts, "expected a sequence");
    if(!buckets || !percentiles || !expo_avgs || !window_counts) { 
        // TODO: decref on buckets and percentiles
        return NULL;
    }
    num_buckets = (int)PySequence_Fast_GET_SIZE(buckets);
    num_percentiles = (int)PySequence_Fast_GET_SIZE(percentiles);
    num_expo_avgs = (int)PySequence_Fast_GET_SIZE(expo_avgs);
    num_window_counts = (int)PySequence_Fast_GET_SIZE(window_counts);

    self = (faststat_Stats*)type->tp_alloc(type, 0);
    if(self != NULL) {
        self->interval = NULL;
        self->n = 0;
        self->mean = self->m2 = self->m3 = self->m4 = self->min = self->max = 0;
        self->sum_of_logs = self->sum_of_inv = 0;
        self->mintime = self->maxtime = self->lasttime = 0;
        self->num_percentiles = num_percentiles;
        if(interval != Py_None ) {
            self->interval = (faststat_Stats*)interval; // WARNING: incompatible pointer type..
        } else {                 // TODO: figure out a better test of type here
            self->interval = NULL;
        }
        if(num_percentiles) {
            self->percentiles = PyMem_New(faststat_P2Percentile, num_percentiles);
            for(i=0; i<num_percentiles; i++) {
                temp = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(percentiles, i));
                self->percentiles[i].percentile = (unsigned short)(temp * 0x10000);
                self->percentiles[i].val = 0;
                self->percentiles[i].n = i + 1;
            }
        } else {
            self->percentiles = NULL;
        }
        self->num_buckets = num_buckets;
        if(num_buckets) {
            self->buckets = PyMem_New(faststat_Bucket, num_buckets);
            for(i=0; i<num_buckets; i++) {
                self->buckets[i].count = 0;
                self->buckets[i].max = (float)PyFloat_AsDouble(PySequence_Fast_GET_ITEM(buckets, i));
                // don't bother checking for error; let it raise later
            }
        } else {
            self->buckets = NULL;
        }
        self->num_expo_avgs = num_expo_avgs;
        if(num_expo_avgs) {
            self->expo_avgs = PyMem_New(faststat_ExpoAvg, num_expo_avgs);
            for(i=0; i<num_expo_avgs; i++) {
                self->expo_avgs[i].val = 0;
                self->expo_avgs[i].alpha = (double)PyFloat_AsDouble(PySequence_Fast_GET_ITEM(expo_avgs, i));
            }
        } else {
            self->expo_avgs = NULL;
        }
        self->num_prev = num_prev;
        if(num_prev) {
            self->lastN = PyMem_New(faststat_DataPoint, num_prev);
            for(i=0; i<num_prev; i++) {
                self->lastN[i].value = 0;
                self->lastN[i].nanotime = 0;
            }
        } else {
            self->lastN = NULL;
        }
        if(num_top == 0) {
            num_top = 1;
        }
        self->num_top = num_top;
        self->topN = PyMem_New(faststat_DataPoint, num_prev);
        self->topN -= 1; //use 1 based indexing

        self->num_window_counts = num_window_counts;
        if(num_window_counts) {
            self->window_counts = PyMem_New(faststat_WindowCount, num_window_counts);
            PyList_Sort(window_counts);
            total = 0;
            for(i=0; i<num_window_counts; i++) {
                cur = PySequence_Fast_GET_ITEM(window_counts, i);
                if(!PyTuple_Check(cur)) {
                    continue;
                }
                PyArg_ParseTuple(cur, "HK",
                    &(self->window_counts[i].num_windows),
                    &(self->window_counts[i].window_size_nanosecs));
                total += self->window_counts[i].num_windows;
            }
            // allocate all of the window counts as one contiguous block
            self->window_counts[0].counts = PyMem_New(unsigned int, total);
            memset(self->window_counts[0].counts, 0, sizeof(unsigned int) * total);
            offset = self->window_counts[0].num_windows;
            for(i=1; i<num_window_counts; i++) {
                self->window_counts[i].counts = self->window_counts[0].counts + offset;
                offset += self->window_counts[i].num_windows;
            }
        } else {
            self->window_counts = NULL;
        }
    }

    if(PyErr_Occurred()) {
       Py_DECREF(self);
       return NULL;
    }

    return (PyObject*) self;
}


static void faststat_Stats_dealloc(faststat_Stats* self) {
    if(self->percentiles) {
        PyMem_Del(self->percentiles);
    }
    if(self->buckets) {
        PyMem_Del(self->buckets);
    }
    if(self->expo_avgs) {
        PyMem_Del(self->expo_avgs);
    }
    if(self->lastN) {
        PyMem_Del(self->lastN);
    }
    if(self->window_counts) {
        // see constructor; all window_counts are allocated as one chunk
        PyMem_Del(self->window_counts[0].counts);
        PyMem_Del(self->window_counts);
    }
}


static PyMemberDef faststat_Stats_members[] = {
    {"n", T_UINT, offsetof(faststat_Stats, n), READONLY, "numder of points"},
    {"mean", T_DOUBLE, offsetof(faststat_Stats, mean), READONLY, "mean"},
    {"min", T_DOUBLE, offsetof(faststat_Stats, min), READONLY, "min"},
    {"max", T_DOUBLE, offsetof(faststat_Stats, max), READONLY, "max"},
    {"sum_of_logs", T_DOUBLE, offsetof(faststat_Stats, sum_of_logs), READONLY, 
        "sum of logs of values, provided all values were positive definite, else NaN (userful for geometric mean)"},
    {"sum_of_inv", T_DOUBLE, offsetof(faststat_Stats, sum_of_inv), READONLY,
        "sum of inverses of values, provided all were positive definite, else NaN (useful for harmonic mean)"},    
    {"lasttime", T_ULONGLONG, offsetof(faststat_Stats, lasttime), READONLY,
                      "time (in nanoseconds since epoch) of last call to add"},
    {"mintime", T_ULONGLONG, offsetof(faststat_Stats, mintime), READONLY, 
                      "time (in nanoseconds since epoch) that minimum value was seen"},
    {"maxtime", T_ULONGLONG, offsetof(faststat_Stats, maxtime), READONLY,
                      "time (in nanoseconds since epoch) that maximum value was seen"},
    {"m2", T_DOUBLE, offsetof(faststat_Stats, m2), READONLY, "m2"},
    {"m3", T_DOUBLE, offsetof(faststat_Stats, m3), READONLY, "m3"},
    {"m4", T_DOUBLE, offsetof(faststat_Stats, m4), READONLY, "m4"},
    {"interval", T_OBJECT, offsetof(faststat_Stats, interval), READONLY, "interval"},
    {"window_avg", T_DOUBLE, offsetof(faststat_Stats, window_avg), READONLY, 
        "average of stored most recent data points"},
    {"num_prev", T_ULONG, offsetof(faststat_Stats, num_prev), READONLY,
        "number of most recent data points stored (accessible via get_prev() )"},
    {NULL}
};


//update mean, and second third and fourth moments
static void _update_moments(faststat_Stats *self, double x) {
    double n, delta, delta_n, delta_m2, delta_m3, delta_m4;
    n = (double)self->n; // note: math with 32 bit ints can cause problems
    //pre-compute a bunch of intermediate values
    delta = x - self->mean;
    delta_n = delta / n;
    delta_m2 = delta * delta_n * (n - 1);
    delta_m3 = delta_m2 * delta_n * (n - 2);
    delta_m4 = delta_m2 * delta_n * delta_n * (n * (n - 3) + 3);
    //compute updated values
    self->mean = self->mean + delta_n;
    //note: order matters here
    self->m4 += delta_m4 + delta_n * (6 * delta_n * self->m2 - 4 * self->m3);
    self->m3 += delta_m3 + delta_n * 3 * self->m2;
    self->m2 += delta_m2;
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


static void _insert_percentile_sorted(faststat_Stats *self, double x) {
    unsigned long long num, i;  // prevent loss of precision compiler warning
    double tmp;
    num = self->n < self->num_percentiles ? self->n : self->num_percentiles;
    for(i = 0; i < num-1; i++) { //insert in sorted order
        if(x < self->percentiles[i].val) {
            tmp = x;
            x = self->percentiles[i].val;
            self->percentiles[i].val = tmp;
        }
    }
    self->percentiles[num-1].val = x;
}


//update percentiles using piece-wise parametric algorithm
static void _update_percentiles(faststat_Stats *self, double x) {
    unsigned int i;
    //double percentile; //TODO: remove me
    faststat_P2Percentile *right, *left, *cur, *prev, *nxt;
    right = &(self->percentiles[self->num_percentiles-1]);
    left = &(self->percentiles[0]);
    if(!(right->n < self->n) ) { // just insert until self->n > self->num_percentiles
        _insert_percentile_sorted(self, x);
        return;
    }
    //right-most is stopping case; handle first
    if(x < right->val && right->n + 1 < self->n) {
        right->n++;
    }
    //handle the rest of the points
    prev = right;
    for(i = self->num_percentiles-2; ; i--) {
        cur = &(self->percentiles[i]);
        if(x < cur->val && cur->n + 1 < prev->n) {
            cur->n++;
        }
        prev = cur;
        if(i == 0) { //making i unsigned fixes some warnings
            break;
        }
    }
    //left-most point is a special case
    nxt = &(self->percentiles[1]);
    _p2_update_point(self->min, 0, left, nxt->val, nxt->n, self->n);
    cur = left;
    for(i=1; i < self->num_percentiles - 1; i++) {
        prev = cur;
        cur = nxt;
        nxt = &(self->percentiles[i+1]);
        _p2_update_point(prev->val, prev->n, cur, nxt->val, nxt->n, self->n);
    }
    _p2_update_point(cur->val, cur->n, right, (double)self->max, (double)self->n, self->n);
} 

// be careful; if-condition must properly terminate when max == +inf, even for nan and +inf
#define OFFSET(n)  if(!(x >= self->buckets[i+n].max)) { self->buckets[i+n].count++; break; } 

static void _update_buckets(faststat_Stats *self, double x) {
    unsigned int i;
    for(i=0; ; i += 16) {
        OFFSET( 0) OFFSET( 1) OFFSET( 2) OFFSET( 3) 
        OFFSET( 4) OFFSET( 5) OFFSET( 6) OFFSET( 7) 
        OFFSET( 8) OFFSET( 9) OFFSET(10) OFFSET(11) 
        OFFSET(12) OFFSET(13) OFFSET(14) OFFSET(15)
    }
}

#undef OFFSET

static void _update_lastN(faststat_Stats *self, double x) {
    unsigned int offset;
    if(self->num_prev == 0) { return; }
    offset = (self->n - 1) & (self->num_prev - 1);
    self->window_avg -= self->lastN[offset].value / (1.0 * self->num_prev);
    self->window_avg += x / (1.0 * self->num_prev);
    self->lastN[offset].value = x;
    self->lastN[offset].nanotime = self->lasttime;
}


//this algorithm deterministically favors storing newer values over older values
static void _update_topN(faststat_Stats *self, double x, unsigned long long t) {
    faststat_DataPoint *cur, *left, *right, *end, *min, *topN;
    unsigned int cur_i, left_i, right_i;
    // uses one based indexing to save additions when navigating down heap
    topN = self->topN;
    if(x < topN[1].value) {
        return;
    }
    // replace the smallest element of the topN with the new point
    topN[1].value = x;
    topN[1].nanotime = t;
    // restore the heap condition
    cur = topN + 1; //use pointers instead of array indices
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


static void _update_expo_avgs(faststat_Stats *self, double x) {
    unsigned int i;
    double val, alpha;
    for(i = 0; i < self->num_expo_avgs; i++) {
        val = self->expo_avgs[i].val;
        alpha = self->expo_avgs[i].alpha;
        // this equation ensures no "gain"
        val = x * alpha + val * (1 - alpha);
        self->expo_avgs[i].val = val;
    }
}

// re-zero all of the windows which have been "missed" between self->lasttime and t
static void _rezero_window_counts(faststat_Stats *self, unsigned long long t) {
    faststat_WindowCount *cur;
    unsigned int i;
    unsigned long long j, last_window, cur_window;
    for(i = 0; i < self->num_window_counts; i++) {
        cur = &(self->window_counts[i]);
        last_window = self->lasttime / cur->window_size_nanosecs;
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


static void _update_window_counts(faststat_Stats *self, unsigned long long t) {
    faststat_WindowCount *cur;
    unsigned int i, cur_count;
    _rezero_window_counts(self, t);
    // step 2 -- increment current counts
    for(i = 0; i < self->num_window_counts; i++) {
        cur = &(self->window_counts[i]);
        // use the current time as the index into the circular array to save some memory
        cur_count = (t / cur->window_size_nanosecs) & (cur->num_windows - 1);
        ++(cur->counts[cur_count]);
    }
}


static void _add(faststat_Stats *self, double x, unsigned long long t) {
    //update extremely basic values: number, min, and max
    self->lasttime = t;
    self->n++;
    if(self->n == 1) {
        self->min = self->max = x;
        self->mintime = self->maxtime = self->lasttime;
    }
    if(x <= self->min) {
        self->mintime = self->lasttime;
        self->min = x;
    }
    if(x >= self->max) {
        self->maxtime = self->lasttime;
        self->max = x;
    }
    // TODO: any platforms not support the NAN macro?
    self->sum_of_logs += x > 0 ? log(x) : NAN;
    self->sum_of_inv += x > 0 ? 1 / x : NAN;
    _update_moments(self, x);
    _update_percentiles(self, x);
    _update_buckets(self, x);
    _update_expo_avgs(self, x);
    _update_lastN(self, x);
    _update_window_counts(self, t);
    _update_topN(self, x, t);
}


static PyObject* faststat_Stats_add(faststat_Stats *self, PyObject *args) {
    //visual studios hates in-line variable declaration
    double x;
    unsigned long long t;
    x = 0;
    t = 0;
    if(PyArg_ParseTuple(args, "d", &x)) {
        t = nanotime();
        if(self->interval && self->lasttime) {
            unsigned long long interval;
            interval = t - self->lasttime;
            if(interval == 0) {
                interval = 1;
            }
            // ensure interval is at least 1 nanosecond to not mess up
            // harmonic and geometric mean (1 ns is noise on this scale)
            _add(self->interval, (double)(interval), t);
        }
        _add(self, x, t);
    }
    if(PyErr_Occurred()) { return NULL; }
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject* faststat_Stats_end(faststat_Stats *self, PyObject *args) {
    unsigned long long end;
    unsigned long long start;
    end = start = 0;
    if(PyArg_ParseTuple(args, "K", &start)) {
        end = nanotime();
        if(self->interval && self->lasttime) {
            _add(self->interval, (double)(end - self->lasttime), end);
        }
        _add(self, (double)(end - start), end);
    }
    if(PyErr_Occurred()) { return NULL; }
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject* faststat_Stats_tick(faststat_Stats *self, PyObject *args) {
    //tricky part is how to handle the first tick
    // weird part will be that calling tick N times results in N-1 data points
    unsigned long long t;
    t = nanotime();
    if(self->lasttime) {
        _add(self, (double)(t - self->lasttime), t);
    } else {
        self->lasttime = t;
    }
    if(PyErr_Occurred()) { return NULL; }
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject* faststat_Stats_get_percentiles(faststat_Stats* self, PyObject *args) {
    PyObject *p_dict;
    faststat_P2Percentile *cur;
    double cur_val;
    unsigned int i;
    p_dict = PyDict_New();
    for(i = 0; i < self->num_percentiles; i++) {
        cur = &(self->percentiles[i]);
        cur_val = ((double)cur->percentile) / 0x10000;
        cur_val = floor(10000 * cur_val + 0.5) / 10000;  
        //re-round to handle slop from being 16 bit number
        // (note: windows math.h does not include round; use floor)
        PyDict_SetItem(
            p_dict, 
            PyFloat_FromDouble(cur_val), 
            PyFloat_FromDouble(cur->val));
    }
    if(PyErr_Occurred()) { 
        Py_DECREF(p_dict);
        return NULL; 
    }
    return p_dict;
}


static PyObject* faststat_Stats_get_buckets(faststat_Stats* self, PyObject *args) {
    PyObject *b_dict;
    faststat_Bucket *cur;
    unsigned int i;
    unsigned long long leftover;
    leftover = self->n;
    b_dict = PyDict_New();
    for(i = 0; i < self->num_buckets; i++) {
        cur = &(self->buckets[i]);
        leftover -= cur->count;
        PyDict_SetItem(
            b_dict,
            PyFloat_FromDouble(cur->max),
            PyLong_FromUnsignedLongLong(cur->count));
    }
    PyDict_SetItem(b_dict, Py_None, PyLong_FromUnsignedLongLong(leftover));
    if(PyErr_Occurred()) { 
        Py_DECREF(b_dict);
        return NULL; 
    }
    return b_dict;
}


static PyObject* faststat_Stats_get_expoavgs(faststat_Stats *self, PyObject *args) {
    PyObject *b_dict;
    faststat_ExpoAvg *cur;
    unsigned int i;
    b_dict = PyDict_New();
    for(i = 0; i < self->num_expo_avgs; i++) {
        cur = &(self->expo_avgs[i]);
        PyDict_SetItem(b_dict,
            PyFloat_FromDouble(cur->alpha),
            PyFloat_FromDouble(cur->val));
    }
    return b_dict;
}


static PyObject* faststat_Stats_get_prev(faststat_Stats *self, PyObject *args) {
    int offset;
    double val;
    PyObject *tuple, *pyval, *pytime;
    unsigned long long nanotime;
    if(self->num_prev == 0) {
        Py_INCREF(Py_None);
        return Py_None;
    }

    offset = 0;
    if(PyArg_ParseTuple(args, "i", &offset)) {
        offset = ((self->n - 1)  + (self->num_prev - offset)) & (self->num_prev - 1);
        val = self->lastN[offset].value;
        nanotime = self->lastN[offset].nanotime;
        pyval = PyFloat_FromDouble(val);
        pytime = PyLong_FromUnsignedLongLong(nanotime);
        if(pyval != NULL && pytime != NULL) {
            tuple = PyTuple_Pack(2, pytime, pyval);
            if(tuple != NULL) {
                return tuple;
            }
        }
    }
    if(PyErr_Occurred()) { return NULL; }
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject* faststat_Stats_get_topN(faststat_Stats *self, PyObject *args) {
    PyObject *ret;
    unsigned int i;
    ret = PyList_New(self->num_top);
    for(i=0; i<self->num_top; i++) {
        PyList_SetItem(ret, i, Py_BuildValue(
            "(dK)", self->topN[i].value, self->topN[i].nanotime));
    }
    if(PyErr_Occurred()) { 
        Py_DECREF(ret);
        return NULL; 
    }
    return ret;
}


static PyObject* faststat_Stats_get_window_counts(faststat_Stats *self, PyObject *args) {
    unsigned long long t;
    PyObject *window_count_dict, *cur_items;
    faststat_WindowCount *cur;
    unsigned long long i, j, cur_window;
    t = nanotime();
    _rezero_window_counts(self, t);
    window_count_dict = PyDict_New();
    for(i = 0; i < self->num_window_counts; i++) {
        cur = self->window_counts + i;
        cur_items = PyTuple_New(cur->num_windows);
        cur_window = t / cur->window_size_nanosecs;
        for(j = 0; j < cur->num_windows; j++) {
            PyTuple_SetItem(
                cur_items, (Py_ssize_t)j, 
                PyLong_FromUnsignedLong(
                    cur->counts[(j + cur_window) & (cur->num_windows - 1)]));
        }
        PyDict_SetItem(
            window_count_dict,
            PyLong_FromUnsignedLongLong(cur->window_size_nanosecs),
            cur_items);
    }
    if(PyErr_Occurred()) { 
        Py_DECREF(window_count_dict);
        return NULL; 
    }
    return window_count_dict;
}


static PyMethodDef faststat_Stats_methods[] = {
    {"add", (PyCFunction)faststat_Stats_add, METH_VARARGS, "add a data point"},
    {"end", (PyCFunction)faststat_Stats_end, METH_VARARGS, 
        "add a duration data point, whose start time is passed"},
    {"tick", (PyCFunction)faststat_Stats_tick, METH_NOARGS, 
        "add an interval data point between now and the last tick"},
    {"get_percentiles", (PyCFunction)faststat_Stats_get_percentiles, METH_NOARGS, 
        "construct percentiles dictionary"},
    {"get_buckets", (PyCFunction)faststat_Stats_get_buckets, METH_NOARGS,
        "construct buckets dictionary"},
    {"get_expo_avgs", (PyCFunction)faststat_Stats_get_expoavgs, METH_NOARGS,
        "get a dictionary of decay rates to previous averages"},
    {"get_prev", (PyCFunction)faststat_Stats_get_prev, METH_VARARGS,
        "get the nth previous sample"},
    {"get_topN", (PyCFunction)faststat_Stats_get_topN, METH_NOARGS,
        "get the highest values"},
    {"get_window_counts", (PyCFunction)faststat_Stats_get_window_counts, METH_NOARGS,
        "get a dictionary of window intervals to window counts"},
    {NULL}
};


static PyTypeObject faststat_StatsType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "_faststat.Stats",          /*tp_name*/
    sizeof(faststat_Stats),    /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)faststat_Stats_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    "online stats collector",  /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    faststat_Stats_methods,    /* tp_methods */
    faststat_Stats_members,    /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    faststat_Stats_new,        /* tp_new */
};


static PyObject* pynanotime(PyObject *_) {
    PyObject *result;
    result = PyLong_FromUnsignedLongLong(nanotime());
    if(PyErr_Occurred()) { return NULL; }
    return result;
}

static PyObject* pynanotime_override(PyObject *_, PyObject *args) {
    unsigned long long t;
    if(PyArg_ParseTuple(args, "K", &t)) {
        nanotime_override = t;
    }
    if(PyErr_Occurred()) { return NULL; }
    Py_INCREF(Py_None);
    return Py_None;
}


static PyMethodDef module_methods[] = { 
    {"nanotime", (PyCFunction)pynanotime, METH_NOARGS, 
        "get integral nanoseconds since unix epoch"},
    {"_nanotime_override", (PyCFunction)pynanotime_override, METH_VARARGS,
        "override time seen by all faststat operations, useful for testing time based algoritmhs"},
    {NULL} };


#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC init_faststat(void) {
    PyObject *module;

    if(PyType_Ready(&faststat_StatsType) < 0)
        return;

    module = Py_InitModule3("_faststat", module_methods, "fast statistics");

    if(module == NULL)
        return;

    Py_INCREF(&faststat_StatsType);
    PyModule_AddObject(module, "Stats", (PyObject*)&faststat_StatsType);
    nanotime_override = 0;
}
