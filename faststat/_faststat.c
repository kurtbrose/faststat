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

static unsigned long long nanotime() {
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

static unsigned long long nanotime() {
    struct timespec ts;
    if(clock_gettime(CLOCK_REALTIME, &ts) == -1) {
        return 0;
    }
    return 1000000000ULL * ts.tv_sec + ts.tv_nsec;
}

#else
//for those oddballs like OSX and BSD, fall back on gettimeofday() which is at least microseconds
#include <time.h>

static unsigned long long nanotime() {
    struct timeval tv;
    if(gettimeofday(&tv, NULL) == -1) {
        return 0;
    }
    return tv.tv_sec * 1000000000ULL + tv.tv_usec * 1000ULL;
}

#endif

//percentile point for usage in P2 algorithm
typedef struct {
    unsigned short percentile;  //divide by 0xFFFF to get a float between 0 and 1
    double val;  //estimate of current percentile value
    double n;  //estimate of how many values were less than this
} faststat_P2Percentile;  // 2 + 8 + 8 = 18 (probably 20 with padding)


typedef struct {
    float max;
    double count;
} faststat_Bucket;  // 8


typedef struct {
    double value;
    double weight;
    unsigned long long nanotime;
} faststat_PrevSample;  // 24

// for representing a normally distributed variable
typedef struct faststat_Stats_struct {
    PyObject_HEAD
    unsigned long long n;  // 8
    double total_weight;  // 8
    double mean, min, max, m2, m3, m4;  // 6 * 8 = 48
    unsigned long long mintime, maxtime, lasttime;  // 3 * 8 = 24
    unsigned int num_percentiles;  // 4
    faststat_P2Percentile *percentiles;  // n * 20
    unsigned int num_buckets;  // 4
    faststat_Bucket *buckets;  // m * 8
    unsigned int num_prev;  // 4
    faststat_PrevSample *lastN;  // p * 16
    struct faststat_Stats_struct *interval;  // x2
} faststat_Stats;  // total = 2 * (8 + 8 + 48 + 24 + 4 + 20n + 4 + 8m + 4 + 16p)
// plausible size estimate: 3k !!  TODO: prioritize memory usage, make size reporting function

/*
typedef struct {
    unsigned int n;
    unsigned int num_prev;

} faststat_StatsGroup;
*/

char* NEW_ARGS[5] = {"buckets", "lastN", "percentiles", "interval", NULL};


static PyObject* faststat_Stats_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    faststat_Stats *self;
    PyObject *buckets, *percentiles, *interval;
    int num_prev, num_buckets, num_percentiles;
    int i;
    double temp;
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "OiOO", NEW_ARGS, 
                                 &buckets, &num_prev, &percentiles, &interval)) {
        return NULL;
    }

    buckets = PySequence_Fast(buckets, "expected a sequence");
    percentiles = PySequence_Fast(percentiles, "expected a sequence");
    if(!buckets || !percentiles) { // TODO: decref on buckets and percentiles
        return NULL;
    }
    num_buckets = (int)PySequence_Fast_GET_SIZE(buckets);
    num_percentiles = (int)PySequence_Fast_GET_SIZE(percentiles);

    self = (faststat_Stats*)type->tp_alloc(type, 0);
    if(self != NULL) {
        self->interval = NULL;
        self->n = 0;
        self->total_weight = 0;
        self->mean = self->m2 = self->m3 = self->m4 = self->min = self->max = 0;
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
                self->percentiles[i].n = 0;
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
        self->num_prev = num_prev;
        if(num_prev) {
            self->lastN = PyMem_New(faststat_PrevSample, num_prev);
            for(i=0; i<num_prev; i++) {
                self->lastN[i].value = 0;
                self->lastN[i].nanotime = 0;
            }
        } else {
            self->lastN = NULL;
        }
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
    if(self->lastN) {
        PyMem_Del(self->lastN);
    }
}


static PyMemberDef faststat_Stats_members[] = {
    {"n", T_UINT, offsetof(faststat_Stats, n), READONLY, "numder of points"},
    {"mean", T_DOUBLE, offsetof(faststat_Stats, mean), READONLY, "mean"},
    {"min", T_DOUBLE, offsetof(faststat_Stats, min), READONLY, "min"},
    {"max", T_DOUBLE, offsetof(faststat_Stats, max), READONLY, "max"},
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


#define MIN_DIFF (1.001)
#define NOT_TOO_CLOSE(larger, smaller) (fabs((larger) / (smaller)) > MIN_DIFF)


//helper for _update_percentiles
static void _p2_update_point(double l_v, double l_n, faststat_P2Percentile *cur,
                            double r_v, double r_n, double total_weight) {
    double percentile, new_val, c_v, c_n, diff, new_n;
    percentile = ((double)cur->percentile) / 0x10000;
    c_n = cur->n;
    new_n = total_weight * percentile;
    if(! ( NOT_TOO_CLOSE(new_n, c_n) || NOT_TOO_CLOSE(c_n, new_n) ) ) {
        return; // avoid massive cancellation when computing difference
    }
    diff = new_n - c_n;

    c_v = cur->val;
    if( NOT_TOO_CLOSE(c_n + diff, l_n) && NOT_TOO_CLOSE(r_n, c_n + diff)) {  // try updating estimate with parabolic
        new_val = c_v + (diff / (r_n - l_n)) * ( 
            (c_n - l_n + diff) * (r_v - c_v) / (r_n - c_n) +
            (r_n - c_n - diff) * (c_v - l_v) / (c_n - l_n));
        if( NOT_TOO_CLOSE(new_val, l_v) || NOT_TOO_CLOSE(r_v, new_val)) {  // fall back on linear
            new_val += c_v + diff * (r_v - c_v) / (r_n - c_n);
        }
        printf("val: %f -> %f  ; n: %f -> %f\n", cur->val, new_val, cur->n, new_n);
        if(new_val != new_val) {

        } else {
            cur->val = new_val;
            cur->n = new_n;
        }
    }
}


static void _insert_percentile_sorted(faststat_Stats *self, double x, double weight) {
    unsigned long long num, i;  // prevent loss of precision compiler warning
    double tmp;
    num = self->n < self->num_percentiles ? self->n : self->num_percentiles;
    for(i = 0; i < num-1; i++) { //insert in sorted order
        if(x < self->percentiles[i].val) {
            tmp = x;
            x = self->percentiles[i].val;
            self->percentiles[i].val = tmp;
            tmp = weight;
            weight = self->percentiles[i].n;
            self->percentiles[i].n = tmp;
        }
    }
    self->percentiles[num-1].val = x;
    self->percentiles[num-1].n = weight;
}


static void _init_weights(faststat_Stats *self) {
    int i;
    double sofar;
    sofar = 0;
    for(i = 0; i < self->num_percentiles; i++) {
        sofar += self->percentiles[i].n;
        self->percentiles[i].n += sofar;
    }
}


//update percentiles using piece-wise parametric algorithm
static void _update_percentiles(faststat_Stats *self, double x, double weight) {
    unsigned int i;
    //double percentile; //TODO: remove me
    faststat_P2Percentile *right, *left, *cur, *prev, *nxt;
    right = &(self->percentiles[self->num_percentiles-1]);
    left = &(self->percentiles[0]);
    if( self->n <= self->num_percentiles ) { // just insert until self->n > self->num_percentiles
        _insert_percentile_sorted(self, x, weight);
        if( self->n == self->num_percentiles) {
            _init_weights(self);
        }
        return;
    }
    //right-most is stopping case; handle first
    if(x < right->val && NOT_TOO_CLOSE(self->total_weight, right->n + weight)) {
        printf("++++");
        right->n += weight;
    }
    //handle the rest of the points
    prev = right;
    for(i = self->num_percentiles-2; ; i--) {
        cur = &(self->percentiles[i]);
        if(x < cur->val && NOT_TOO_CLOSE(prev->n, cur->n + weight)) {
            printf("++++");
            cur->n += weight;
        }
        prev = cur;
        if(i == 0) { //making i unsigned fixes some warnings
            break;
        }
    }
    //left-most point is a special case
    nxt = &(self->percentiles[1]);
    _p2_update_point(self->min, 0, left, nxt->val, nxt->n, self->total_weight);
    cur = left;
    for(i=1; i < self->num_percentiles - 1; i++) {
        prev = cur;
        cur = nxt;
        nxt = &(self->percentiles[i+1]);
        _p2_update_point(prev->val, prev->n, cur, nxt->val, nxt->n, self->total_weight);
    }
    _p2_update_point(cur->val, cur->n, right, self->max, self->n, self->total_weight);
} 


static void _update_buckets(faststat_Stats *self, double x) {
    unsigned int i;
    for(i=0; i < self->num_buckets; i++) {
        if(x < self->buckets[i].max) {
            self->buckets[i].count++;
            break;
        }
    }
}


static void _update_lastN(faststat_Stats *self, double x) {
    unsigned int offset;
    if(self->num_prev == 0) { return; }
    offset = (self->n - 1) % self->num_prev;
    self->lastN[offset].value = x;
    self->lastN[offset].nanotime = self->lasttime;
}


static void _add(faststat_Stats *self, double x, unsigned long long t) {
    //update extremely basic values: number, min, and max
    self->lasttime = t;
    self->n++;
    self->total_weight += 1.0;
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
    _update_moments(self, x);
    _update_percentiles(self, x, 1.0);
    _update_buckets(self, x);
    _update_lastN(self, x);
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
            _add(self->interval, (double)(t - self->lasttime), t);
        }
        _add(self, x, t);
    }
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
    Py_INCREF(p_dict);
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
    Py_INCREF(b_dict);
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
        offset = ((self->n - 1)  + (self->num_prev - offset)) % self->num_prev;
        val = self->lastN[offset].value;
        nanotime = self->lastN[offset].nanotime;
        pyval = PyFloat_FromDouble(val);
        pytime = PyLong_FromUnsignedLongLong(nanotime);
        if(pyval != NULL && pytime != NULL) {
            tuple = PyTuple_Pack(2, pytime, pyval);
            if(tuple != NULL) {
                Py_INCREF(pyval);
                Py_INCREF(pytime);
                Py_INCREF(tuple);
                return tuple;
            }
        }
    }
    Py_INCREF(Py_None);
    return Py_None;
}


static PyMethodDef faststat_Stats_methods[] = {
    {"add", (PyCFunction)faststat_Stats_add, METH_VARARGS, "add a data point"},
    {"end", (PyCFunction)faststat_Stats_end, METH_VARARGS, "add a duration data point, whose start time is passed"},
    {"tick", (PyCFunction)faststat_Stats_tick, METH_NOARGS, "add an interval data point between now and the last tick"},
    {"get_percentiles", (PyCFunction)faststat_Stats_get_percentiles, METH_NOARGS, 
                "construct percentiles dictionary"},
    {"get_buckets", (PyCFunction)faststat_Stats_get_buckets, METH_NOARGS,
                "construct buckets dictionary"},
    {"get_prev", (PyCFunction)faststat_Stats_get_prev, METH_VARARGS,
                "get the nth previous sample"},
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


PyObject* pynanotime(PyObject* args) {
    PyObject *result;
    result = PyLong_FromUnsignedLongLong(nanotime());
    Py_INCREF(result);
    return result;
}


static PyMethodDef module_methods[] = { 
    {"nanotime", (PyCFunction)pynanotime, METH_NOARGS, 
                "get integral nanoseconds since unix epoch"},
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
}
