// adapted from http://blog.debao.me/2013/04/my-first-c-extension-to-numpy/
// and http://docs.scipy.org/doc/numpy-1.10.1/user/c-info.how-to-extend.html

#include "Python.h"
#include "numpy/arrayobject.h"

/*
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
*/

/*
static PyObject *
run_sweep(PyObject *dummy, PyObject *args)
{
    PyObject *arg1=NULL, *arg2=NULL, *arg3=NULL, *out=NULL;
    PyObject *oarr=NULL;

    if (!PyArg_ParseTuple(args, "OOOO!", &arg1, &arg2, &arg3,
        &PyArray_Type, &out)) return NULL;

    oarr = PyArray_FROM_OTF(out, NPY_DOUBLE, NPY_INOUT_ARRAY);
    if (oarr == NULL) goto fail;

     code that makes use of arguments
     You will probably need at least
       nd = PyArray_NDIM(<..>)    -- number of dimensions
       dims = PyArray_DIMS(<..>)  -- npy_intp array of length nd
                                     showing length in each dim.
       dptr = (double *)PyArray_DATA(<..>) -- pointer to data.

       If an error occurs goto fail.



    PyArray_GETPTR3 (E, i, j)



    Py_DECREF(oarr);
    Py_INCREF(Py_None);
    return Py_None;

 fail:
    PyArray_XDECREF_ERR(oarr);
    return NULL;
}
*/

static PyObject*
run_sweep (PyObject *dummy, PyObject *args)
{
    PyObject *arg1=NULL, *arg2=NULL, *out=NULL;
    PyArrayObject *arr1=NULL, *arr2=NULL, *oarr=NULL;

    if (!PyArg_ParseTuple(args, "OOO!", &arg1, &arg2,
        &PyArray_Type, &out)) return NULL;

    arr1 = (PyArrayObject*)PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_IN_ARRAY);
    if (arr1 == NULL) return NULL;
    arr2 = (PyArrayObject*)PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_IN_ARRAY);
    if (arr2 == NULL) goto fail;
    oarr = (PyArrayObject*)PyArray_FROM_OTF(out, NPY_DOUBLE, NPY_INOUT_ARRAY);
    if (oarr == NULL) goto fail;

    /*vv* code that makes use of arguments *vv*/

    int nd = PyArray_NDIM(arr1);   //number of dimensions
    npy_intp *shape = PyArray_DIMS(arr1);  // npy_intp array of length nd showing length in each dim.

    int i;
    int j;
    for (i=0; i<oarr->dimensions[0]; ++i) {
        for (j=0; j<oarr->dimensions[1]; ++j) {
            double *v = (double*)PyArray_GETPTR2(oarr, i, j);
            *v = *v * 2;
        }
    }
    /*^^* code that makes use of arguments *^^*/

    Py_DECREF(arr1);
    Py_DECREF(arr2);
    Py_DECREF(oarr);
    Py_INCREF(Py_None);
    return Py_None;

 fail:
    Py_XDECREF(arr1);
    Py_XDECREF(arr2);
    PyArray_XDECREF_ERR(oarr);
    return NULL;
}


static struct PyMethodDef methods[] = {
    {"run_sweep", run_sweep, METH_VARARGS, "description of run_sweep"},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

PyMODINIT_FUNC
initsweep (void)
{
    (void)Py_InitModule("sweep", methods);
    import_array();
}
