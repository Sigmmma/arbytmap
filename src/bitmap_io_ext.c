#include <stdio.h>
#include "Python.h"
#include "abstract.h"
#include "longobject.h"
#include "modsupport.h"
#include "object.h"


static void pad_24bit_array(
    Py_buffer *padded_pix_buf, Py_buffer *unpadded_pix_buf)
{
    unsigned long long i=0, j=0, max_i=0;

    unsigned long (*padded_pix)[];
    unsigned char (*unpadded_pix)[];

    padded_pix   = (unsigned long(*)[]) padded_pix_buf->buf;
    unpadded_pix = (unsigned char(*)[]) unpadded_pix_buf->buf;

    max_i = (padded_pix_buf->len)/4;

    for (i=0; i < max_i; i++) {
        j = i*3;
        (*padded_pix)[i] = ((*unpadded_pix)[j] +
                           ((*unpadded_pix)[j+1]<<8) +
                           ((*unpadded_pix)[j+2]<<16));
    }
}

static void unpad_24bit_array(
    Py_buffer *unpadded_pix_buf, Py_buffer *padded_pix_buf)
{
    unsigned long long i=0, j=0, k=0, max_i=0;
    unsigned long pixel=0;

    unsigned char (*unpadded_pix)[];
    unsigned char (*padded_pix_8)[];
    unsigned long (*padded_pix_32)[];

    unpadded_pix  = (unsigned char(*)[]) unpadded_pix_buf->buf;
    padded_pix_8  = (unsigned char(*)[]) padded_pix_buf->buf;
    padded_pix_32 = (unsigned long(*)[]) padded_pix_buf->buf;

    max_i = (unpadded_pix_buf->len)/3;

    if ( padded_pix_buf->itemsize == 4) {
        for (i=0; i < max_i; i++) {
            j=i*3;
            pixel = (*padded_pix_32)[i];
            (*unpadded_pix)[j] = pixel&0xFF; j++;
            (*unpadded_pix)[j] = (pixel>>8)&0xFF; j++;
            (*unpadded_pix)[j] = (pixel>>16)&0xFF;
        }
    } else {
        for (i=0; i < max_i; i++) {
            /*
            If the itemsize isnt 4, the pixels havent been
            packed, and since the channel order in Arbytmap
            is ARGB, we need to skip each pixels first byte
            */
            j=i*3; k=i*4 + 1;
            (*unpadded_pix)[j] = (*padded_pix_8)[k]; j++; k++;
            (*unpadded_pix)[j] = (*padded_pix_8)[k]; j++; k++;
            (*unpadded_pix)[j] = (*padded_pix_8)[k];
        }
    }
}

static void pad_48bit_array(
    Py_buffer *padded_pix_buf, Py_buffer *unpadded_pix_buf)
{
    unsigned long long i=0, j=0, max_i=0;

    unsigned long long (*padded_pix)[];
    unsigned short (*unpadded_pix)[];

    padded_pix   = (unsigned long long(*)[]) padded_pix_buf->buf;
    unpadded_pix = (unsigned short(*)[]) unpadded_pix_buf->buf;

    max_i = (padded_pix_buf->len)/8;

    for (i=0; i < max_i; i++) {
        j = i*3;
        (*padded_pix)[i] = ((unsigned long long)(*unpadded_pix)[j] +
                           ((unsigned long long)(*unpadded_pix)[j+1]<<16) +
                           ((unsigned long long)(*unpadded_pix)[j+2]<<32));
    }
}

static void unpad_48bit_array(
    Py_buffer *unpadded_pix_buf, Py_buffer *padded_pix_buf)
{
    unsigned long long i=0, j=0, k=0, max_i=0, pixel=0;

    unsigned short (*unpadded_pix)[];
    unsigned short (*padded_pix_16)[];
    unsigned long long(*padded_pix_64)[];

    unpadded_pix  = (unsigned short(*)[]) unpadded_pix_buf->buf;
    padded_pix_16 = (unsigned short(*)[]) padded_pix_buf->buf;
    padded_pix_64 = (unsigned long long(*)[]) padded_pix_buf->buf;

    max_i = (unpadded_pix_buf->len)/6;

    if ( padded_pix_buf->itemsize == 8) {
        for (i=0; i < max_i; i++) {
            j=i*3;
            pixel = (*padded_pix_64)[i];
            (*unpadded_pix)[j] = pixel&0xFFff; j++;
            (*unpadded_pix)[j] = (pixel>>16)&0xFFff; j++;
            (*unpadded_pix)[j] = (pixel>>32)&0xFFff;
        }
    } else {
        for (i=0; i < max_i; i++) {
            /*
            If the itemsize isnt 8, the pixels havent been
            packed, and since the channel order in Arbytmap
            is ARGB, we need to skip each pixels first 2 bytes
            */
            j=i*3; k=i*4 + 1;
            (*unpadded_pix)[j] = (*padded_pix_16)[k]; j++; k++;
            (*unpadded_pix)[j] = (*padded_pix_16)[k]; j++; k++;
            (*unpadded_pix)[j] = (*padded_pix_16)[k];
        }
    }
}

static PyObject *py_pad_24bit_array(PyObject *self, PyObject *args) {
    Py_buffer bufs[2];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*:pad_24bit_array", &bufs[0], &bufs[1]))
        return Py_None;

    pad_24bit_array(&bufs[0], &bufs[1]);

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);

    return Py_None;
}

static PyObject *py_pad_48bit_array(PyObject *self, PyObject *args) {
    Py_buffer bufs[2];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*:pad_48bit_array", &bufs[0], &bufs[1]))
        return Py_None;

    pad_48bit_array(&bufs[0], &bufs[1]);

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);

    return Py_None;
}

static PyObject *py_unpad_24bit_array(PyObject *self, PyObject *args) {
    Py_buffer bufs[2];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*:unpad_24bit_array", &bufs[0], &bufs[1]))
        return Py_None;

    unpad_24bit_array(&bufs[0], &bufs[1]);

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);

    return Py_None;
}

static PyObject *py_unpad_48bit_array(PyObject *self, PyObject *args) {
    Py_buffer bufs[2];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*:unpad_48bit_array", &bufs[0], &bufs[1]))
        return Py_None;

    unpad_48bit_array(&bufs[0], &bufs[1]);

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);

    return Py_None;
}

/* A list of all the methods defined by this module.
"METH_VARGS" tells Python how to call the handler.
The {NULL, NULL} entry indicates the end of the method definitions */
static PyMethodDef bitmap_io_ext_methods[] = {
    {"pad_24bit_array", py_pad_24bit_array, METH_VARARGS, ""},
    {"pad_48bit_array", py_pad_48bit_array, METH_VARARGS, ""},
    {"unpad_24bit_array", py_unpad_24bit_array, METH_VARARGS, ""},
    {"unpad_48bit_array", py_unpad_48bit_array, METH_VARARGS, ""},
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
