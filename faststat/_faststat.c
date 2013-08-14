#include <Python.h>
#include <structmember.h>

typedef struct {
    PyObject_HEAD
    int n;
    double mean, m2, m3, m4, min, max;
} faststat_Stats;

static PyObject* faststat_Stats_new(PyTypeObject *type, PyObject* args, PyObject *kwds) {
    faststat_Stats* self;

    self = (faststat_Stats*)type->tp_alloc(type, 0);
    if(self != NULL) {
        self->n = 0;
        self->mean = self->m2 = self->m3 = self->m4 = self->min = self->max = 0;
    }

    return (PyObject*) self;
}

static PyMemberDef faststat_Stats_members[] = {
    {"n", T_INT, offsetof(faststat_Stats, n), 0, "numder of points"},
    {"mean", T_DOUBLE, offsetof(faststat_Stats, mean), 0, "mean"},
    {"min", T_DOUBLE, offsetof(faststat_Stats, min), 0, "min"},
    {"max", T_DOUBLE, offsetof(faststat_Stats, max), 0, "max"},
    {"m2", T_DOUBLE, offsetof(faststat_Stats, m2), 0, "m2"},
    {"m3", T_DOUBLE, offsetof(faststat_Stats, m3), 0, "m3"},
    {"m4", T_DOUBLE, offsetof(faststat_Stats, m4), 0, "m4"},
    {NULL}
};

static PyObject* faststat_Stats_add(faststat_Stats *self, PyObject *args) {
    double x = 0;
    if(PyArg_ParseTuple(args, "d", &x)) {
        self->n++;
        //pre-compute a bunch of intermediate values
        double n = self->n; // note: math with 32 bit ints can cause problems
        double delta = x - self->mean;
        double delta_n = delta / n;
        double delta_m2 = delta * delta_n * (n - 1);
        double delta_m3 = delta_m2 * delta_n * (n - 2);
        double delta_m4 = delta_m2 * delta_n * delta_n * (n * (n - 3) + 3);
        //compute updated values
        self->min = x < self->min ? x : self->min;
        self->max = x > self->max ? x : self->max;
        self->mean = self->mean + delta_n;
        //note: order matters here
        self->m4 += delta_m4 + delta_n * (6 * delta_n * self->m2 - 4 * self->m3);
        self->m3 += delta_m3 + delta_n * 3 * self->m2;
        self->m2 += delta_m2;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef faststat_Stats_methods[] = {
    {"add", (PyCFunction)faststat_Stats_add, METH_VARARGS, "add a data point"},
    {NULL}
};

static PyTypeObject faststat_StatsType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "_faststat.Stats",          /*tp_name*/
    sizeof(faststat_Stats),    /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    0,                         /*tp_dealloc*/
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
