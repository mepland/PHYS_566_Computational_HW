// adapted from http://blog.debao.me/2013/04/my-first-c-extension-to-numpy/
// and http://docs.scipy.org/doc/numpy-1.10.1/user/c-info.how-to-extend.html

// compile and install module with:
// python module_setup_for_sweep.py install
// or TODO
// python module_setup_for_sweep.py build_ext --inplace


/*
test code

import sweep
import numpy as np
a = np.array([[1.1, 2.2, 3.3], [1.2, 1.3, 1.4]])
a = np.array([[1.1, 1.2, 1.3], [2.1, 2.2, 2.3], [3.1, 3.2, 3.3]])

sweep.run_sweep(10.0, 0, a)
print a

sweep.run_sweep(10.0, 1, a)
print a

exit()

*/

#include "Python.h"
#include "numpy/arrayobject.h"

static PyObject*
run_sweep (PyObject *dummy, PyObject *args)
{
    double T;
    int NN_type;

    PyObject *out=NULL;
    PyArrayObject *oarr=NULL;

    if (!PyArg_ParseTuple(args, "diO!", &T, &NN_type,
        &PyArray_Type, &out)) return NULL;

    oarr = (PyArrayObject*)PyArray_FROM_OTF(out, NPY_DOUBLE, NPY_INOUT_ARRAY);
    if (oarr == NULL) goto fail;

    /*vv* code that makes use of arguments *vv*/

    int i;
    int j;
    for (i=0; i < oarr->dimensions[0]; ++i) {
        for (j=0; j < oarr->dimensions[1]; ++j) {
            double *v = (double*)PyArray_GETPTR2(oarr, i, j);
            if(NN_type == 0){
              *v = *v + T;
            }
            else if(NN_type == 1){
              *v = *v + 2*T;
            }

        }
    }
    /*^^* code that makes use of arguments *^^*/


    Py_DECREF(oarr);
    Py_INCREF(Py_None);
    return Py_None;

 fail:
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
