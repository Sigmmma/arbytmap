#include <stdio.h>
#include "Python.h"
#include "abstract.h"
#include "longobject.h"
#include "modsupport.h"
#include "object.h"


static PyObject *py_swizzle_array(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*:swizzle_array",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5]))
        return Py_None;

    swizzle_char_array(
        bufs[0].len, bufs[0].buf, bufs[1].len, bufs[1].buf,
        bufs[2].len, bufs[2].buf, bufs[3].len, bufs[3].buf,
        bufs[4].len, bufs[4].buf, bufs[5].len, bufs[5].buf,
        bufs[5].itemsize, 1);

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]);
    PyBuffer_Release(&bufs[3]);
    PyBuffer_Release(&bufs[4]);
    PyBuffer_Release(&bufs[5]);

    /* Convert from a C integer value to a Python integer instance */
    return Py_None;
}

static PyObject *py_unswizzle_array(PyObject *self, PyObject *args) {
    Py_buffer bufs[2];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*:swizzle_array",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5]))
        return Py_None;

    swizzle_char_array(
        bufs[0].len, bufs[0].buf, bufs[1].len, bufs[1].buf,
        bufs[2].len, bufs[2].buf, bufs[3].len, bufs[3].buf,
        bufs[4].len, bufs[4].buf, bufs[5].len, bufs[5].buf,
        bufs[5].itemsize, 0);

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);

    /* Convert from a C integer value to a Python integer instance */
    return Py_None;
}

/* A list of all the methods defined by this module.
"METH_VARGS" tells Python how to call the handler.
The {NULL, NULL} entry indicates the end of the method definitions */
static PyMethodDef bitmap_io_ext_methods[] = {
    //{"pad_24bit_array", py_pad_24bit_array, METH_VARARGS, ""},
    //{"pad_48bit_array", py_pad_48bit_array, METH_VARARGS, ""},
    //{"unpad_24bit_array", py_unpad_24bit_array, METH_VARARGS, ""},
    //{"unpad_48bit_array", py_unpad_48bit_array, METH_VARARGS, ""},
    //{"uncompress_rle", py_uncompress_rle, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}      /* sentinel */
};

/* When Python imports a C module named 'X' it loads the
module then looks for a method named "init"+X and calls it.*/
static struct PyModuleDef bitmap_io_ext_module = {
    PyModuleDef_HEAD_INIT,
    "bitmap_io_ext",
    "A set of C functions to replace certain speed intensive bitmap io functions",
    -1,
    bitmap_io_ext_methods,
};

PyMODINIT_FUNC PyInit_bitmap_io_ext(void) {
    return PyModule_Create(&bitmap_io_ext_module);
}
