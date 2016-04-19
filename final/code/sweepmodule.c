// addapted from http://blog.debao.me/2013/04/my-first-c-extension-to-numpy/

#include "Python.h"
#include "numpy/arrayobject.h"

static PyObject*
run_sweep (PyObject *dummy, PyObject *args)
{
    PyObject *arg1=NULL;
    PyObject *arr1=NULL;
    int nd;

    if (!PyArg_ParseTuple(args, "O", &arg1))
        return NULL;

    arr1 = PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_IN_ARRAY);
    if (arr1 == NULL)
        return NULL;

    nd = PyArray_NDIM(arr1);   //number of dimensions

    Py_DECREF(arr1);

    return PyInt_FromLong(nd);
}

static struct PyMethodDef methods[] = {
    {"run_sweep", run_sweep, METH_VARARGS, "descript of run_sweep"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initsweep (void)
{
    (void)Py_InitModule("sweep", methods);
    import_array();
}
