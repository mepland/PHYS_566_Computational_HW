// adapted from
// http://blog.debao.me/2013/04/my-first-c-extension-to-numpy/
// http://docs.scipy.org/doc/numpy-1.10.1/user/c-info.how-to-extend.html
// http://stackoverflow.com/questions/24661141/argument-passing-in-python-module-written-in-c

// compile and install module with:
// python module_setup_for_sweepMod.py install
// or
// python module_setup_for_sweepMod.py build_ext --inplace

#include "Python.h"
#include "numpy/arrayobject.h"
#include <math.h>
#include <stdlib.h>

// make a mod function that works with negative numbers
// from http://stackoverflow.com/questions/11720656/modulo-operation-with-negative-numbers

int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}

// sweepMod.run_sweep(T, kB, J, NN_type, seed, world_grid)

static PyObject*
run_sweep (PyObject *dummy, PyObject *args)
{
    double T, kB, J;
    int NN_type, seed;

    PyObject *out=NULL;
    PyArrayObject *oarr=NULL;

    if (!PyArg_ParseTuple(args, "dddiiO!", &T, &kB, &J, &NN_type, &seed,
        &PyArray_Type, &out)) return NULL;

    oarr = (PyArrayObject*)PyArray_FROM_OTF(out, NPY_DOUBLE, NPY_INOUT_ARRAY);
    if (oarr == NULL) goto fail;



    /*vv* code that makes use of arguments *vv*/

    int i, j, k; // declare for loop variables, have to do it here, limited by old compiler

    // make arrays for NN points, first 4 are Von Neumann, last 4 add the corners for Moore
    int NN_x[] = { 1,  0, -1,  0,  1, -1, -1,  1};
    int NN_y[] = { 0,  1,  0, -1,  1,  1, -1, -1};

    int NN_size = 0;

    if(NN_type == 0) NN_size = 4; // Von Neumann
    else if(NN_type == 1) NN_size = 8; // Von Neumann


    // setup the rng seed, NOTE MUST BE CAREFUL TO GIVE GOOD SEEDS
    srand((unsigned)seed);


    // loop over the grid

    int n = oarr->dimensions[0]; // find n from array

    for(i=0; i < n; ++i){
      for(j=0; j < n; ++j){
        double *this_point = (double*)PyArray_GETPTR2(oarr, i, j); // get i, j elements value

        // compute DeltaE

        // find E original
        double E_NN_original = 0.0;
        for(k=0; k < NN_size; ++k){
          double *NN_point1 = (double*)PyArray_GETPTR2(oarr, mod((i+NN_x[k]), n), mod((j+NN_y[k]), n) ); // get kth NN elements value
          E_NN_original += -J*(*this_point)*(*NN_point1);
        } // end k loop

        // test flip the spin
        *this_point = -(*this_point);

        // find E flipped
        double E_NN_flipped = 0.0;
        for(k=0; k < NN_size; ++k){
          double *NN_point2 = (double*)PyArray_GETPTR2(oarr, mod((i+NN_x[k]), n), mod((j+NN_y[k]), n) ); // get kth NN elements value
          E_NN_flipped += -J*(*this_point)*(*NN_point2);
        } // end k loop

        // actually compute DeltaE
        double DeltaE = E_NN_flipped - E_NN_original;


        // if DeltaE <= 0 always keep the spin, ie do nothing as it has already been flipped
        if( DeltaE > 0.0){ // if here, keep the spin with probability p = exp(-DeltaE/kB*T)

          double p = exp( -DeltaE/(kB*T) );
          double r = (double)rand() / (double)((unsigned)RAND_MAX + 1); // [0, 1) rng

           // if r < p keep the spin, ie do nothing as it has already been flipped
           if( r >= p ){

             // flip the spin back to its original position
             *this_point = -(*this_point);

           } // end r >= p if statement

        } // end DeltaE > 0.0 if statement

      } // end j loop
    } // end i loop

    /*^^* code that makes use of arguments *^^*/


    Py_DECREF(oarr);
    Py_INCREF(Py_None);
    return Py_None;

 fail:
    PyArray_XDECREF_ERR(oarr);
    return NULL;
} // end run_sweep()

// sweepMod.find_E(J, NN_type, world_grid)

static PyObject*
find_E (PyObject *dummy, PyObject *args)
{
    double J;
    int NN_type;

    PyObject *in=NULL;
    PyArrayObject *iarr=NULL;

    if (!PyArg_ParseTuple(args, "diO!", &J, &NN_type,
        &PyArray_Type, &in)) return NULL;

    iarr = (PyArrayObject*)PyArray_FROM_OTF(in, NPY_DOUBLE, NPY_ENSURECOPY); // use NPY_ENSURECOPY to work off a copy of the in array, safer and prevents errors
    if (iarr == NULL) return NULL;


    /*vv* code that makes use of arguments *vv*/

    int i, j, k; // declare for loop variables, have to do it here, limited by old compiler

    // make arrays for NN points, first 4 are Von Neumann, last 4 add the corners for Moore
    int NN_x[] = { 1,  0, -1,  0,  1, -1, -1,  1};
    int NN_y[] = { 0,  1,  0, -1,  1,  1, -1, -1};

    int NN_size = 0;

    if(NN_type == 0) NN_size = 4; // Von Neumann
    else if(NN_type == 1) NN_size = 8; // Von Neumann

    // setup E
    double E = 0.0;

    // loop over the grid

    int n = iarr->dimensions[0]; // find n from array

    for(i=0; i < n; ++i){
      for(j=0; j < n; ++j){
        double *this_point = (double*)PyArray_GETPTR2(iarr, i, j); // get i, j elements value

        // compute the contribution to E from this spin 
        for(k=0; k < NN_size; ++k){
          double *NN_point = (double*)PyArray_GETPTR2(iarr, mod((i+NN_x[k]), n), mod((j+NN_y[k]), n) ); // get kth NN elements value
          E += -J*(*this_point)*(*NN_point);
        } // end k loop

      } // end j loop
    } // end i loop

    /*^^* code that makes use of arguments *^^*/

    PyObject *ret = Py_BuildValue("d", E);
    return ret;

} // end find_E()


static struct PyMethodDef methods[] = {
    {"run_sweep", run_sweep, METH_VARARGS, "description of run_sweep"},
    {"find_E", find_E, METH_VARARGS, "description of find_E"},
    {NULL, NULL, 0, NULL} /* Sentinel */
};




PyMODINIT_FUNC
initsweepMod (void)
{
    (void)Py_InitModule("sweepMod", methods);
    import_array();
}
