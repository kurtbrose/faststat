#include <Python.h>
#include <structmember.h>
#include <pymem.h>
#include <stdio.h>


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


typedef struct {
    double value;
    double timestamp;
} faststat_PrevSample;


typedef struct {
    PyObject_HEAD
    unsigned int n;
    double mean, m2, m3, m4, min, max;
    int num_percentiles;
    faststat_P2Percentile* percentiles;
    int num_buckets;
    faststat_Bucket* buckets;
    int num_prev;
    faststat_PrevSample* lastN;
} faststat_Stats;


char* NEW_ARGS[3] = {"buckets", "lastN", "percentiles"};

static PyObject* faststat_Stats_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    faststat_Stats *self;
    PyObject *buckets, *percentiles;
    int num_prev, num_buckets, num_percentiles;
    int i;
    double temp;
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "OdO", NEW_ARGS, 
                                 &buckets, &num_prev, &percentiles)) {
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
        self->n = 0;
        self->mean = self->m2 = self->m3 = self->m4 = self->min = self->max = 0;
        self->num_percentiles = num_percentiles;
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
                self->buckets[i].max = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(buckets, i));
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
                self->lastN[i].timestamp = 0;
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
    {"n", T_UINT, offsetof(faststat_Stats, n), 0, "numder of points"},
    {"mean", T_DOUBLE, offsetof(faststat_Stats, mean), 0, "mean"},
    {"min", T_DOUBLE, offsetof(faststat_Stats, min), 0, "min"},
    {"max", T_DOUBLE, offsetof(faststat_Stats, max), 0, "max"},
    {"m2", T_DOUBLE, offsetof(faststat_Stats, m2), 0, "m2"},
    {"m3", T_DOUBLE, offsetof(faststat_Stats, m3), 0, "m3"},
    {"m4", T_DOUBLE, offsetof(faststat_Stats, m4), 0, "m4"},
    {NULL}
};


//update mean, and second third and fourth moments
static void _update_moments(faststat_Stats *self, double x) {
    double n, delta, delta_n, delta_m2, delta_m3, delta_m4;
    n = self->n; // note: math with 32 bit ints can cause problems
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
                            double r_v, double r_n, unsigned int n) {
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
    int num, i;
    double tmp;
    num = self->n < self->num_percentiles ? self->n : self->num_percentiles;
    for(i = 0; i < num-1; i++) { //insert in sorted order
        if(x < self->percentiles[i].val) {
            tmp = x;
            x = self->percentiles[i].val;
            self->percentiles[i].val = x;
        }
    }
    self->percentiles[num-1].val = x;
}


//update percentiles using piece-wise parametric algorithm
static void _update_percentiles(faststat_Stats *self, double x) {
    int i;
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
    for(i = self->num_percentiles-2; i >= 0; i--) {
        cur = &(self->percentiles[i]);
        if(x < cur->val && cur->n + 1 < prev->n) {
            cur->n++;
        }
        prev = cur;
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
    _p2_update_point(cur->val, cur->n, right, self->max, self->n, self->n);
} 


static void _update_buckets(faststat_Stats *self, double x) {

}


static PyObject* faststat_Stats_add(faststat_Stats *self, PyObject *args) {
    //visual studios hates in-line variable declaration
    double x;
    x = 0;
    if(PyArg_ParseTuple(args, "d", &x)) {
        //update extremely basic values: number, min, and max
        self->n++;
        if(self->n == 1) {
            self->min = self->max = x;
        }
        self->min = x < self->min ? x : self->min;
        self->max = x > self->max ? x : self->max;
        _update_moments(self, x);
        _update_percentiles(self, x);
        _update_buckets(self, x);
    }
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject* faststat_Stats_get_percentiles(faststat_Stats* self, PyObject *args) {
    PyObject *p_dict;
    faststat_P2Percentile *cur;
    double cur_val;
    int i;
    p_dict = PyDict_New();
    for(i = 0; i < self->num_percentiles; i++) {
        cur = &(self->percentiles[i]);
        cur_val = ((double)cur->percentile) / 0x10000;
        cur_val = round(10000 * cur_val) / 10000;  //re-round to handle slop from being 16 bit number
        PyDict_SetItem(
            p_dict, 
            PyFloat_FromDouble(cur_val), 
            PyFloat_FromDouble(cur->val));
    }
    Py_INCREF(p_dict);
    return p_dict;
}


static PyMethodDef faststat_Stats_methods[] = {
    {"add", (PyCFunction)faststat_Stats_add, METH_VARARGS, "add a data point"},
    {"get_percentiles", (PyCFunction)faststat_Stats_get_percentiles, METH_NOARGS, 
                "construct percentiles dictionary"},
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


static PyMethodDef module_methods[] = { {NULL} };


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
