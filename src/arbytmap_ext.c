#include "shared.h"


static void depalettize_bitmap_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *unpacked_idx_buf,
    Py_buffer *unpacked_pal_buf, sint8 ucc)
{
    uint64 i=0, j=0, max_i=0;
    uint8 *unpacked_pix;
    uint8 *unpacked_pal;
    uint8 *unpacked_idx;

    unpacked_pix = (uint8*)unpacked_pix_buf->buf;
    unpacked_pal = (uint8*)unpacked_pal_buf->buf;
    unpacked_idx = (uint8*)unpacked_idx_buf->buf;
    max_i = unpacked_pix_buf->len;

    if (ucc == 4) {
        for (i=0; i < max_i; i += 4) {
            j = unpacked_idx[i>>2]<<2;
            unpacked_pix[i]   = unpacked_pal[j];
            unpacked_pix[i+1] = unpacked_pal[j+1];
            unpacked_pix[i+2] = unpacked_pal[j+2];
            unpacked_pix[i+3] = unpacked_pal[j+3];
        }
    } else if (ucc == 2) {
        for (i=0; i < max_i; i += 2) {
            j = unpacked_idx[i>>1]<<1;
            unpacked_pix[i]   = unpacked_pal[j];
            unpacked_pix[i+1] = unpacked_pal[j+1];
        }
    } else {
        for (i=0; i < max_i; i++) {
            unpacked_pix[i] = unpacked_pal[unpacked_idx[i]];
        }
    }
}

static void downsample_bitmap_8(
    Py_buffer *downsamp_pix_buf, Py_buffer *swizzled_pix_buf,
    uint32 merge_count, sint8 ucc)
{
    uint64 i=0, j=0, max_i=0, swizz_i=0;
    uint64 merge_val=0;
    uint8 *downsamp_pix;
    uint8 *swizzled_pix;

    downsamp_pix = (uint8*)downsamp_pix_buf->buf;
    swizzled_pix = (uint8*)swizzled_pix_buf->buf;
    max_i = downsamp_pix_buf->len;

    for (i=0; i < max_i; i++) {
        merge_val = 0;
        for (j=0; j < merge_count; j++) {
            merge_val += swizzled_pix[swizz_i];
            swizz_i++;
        }
        downsamp_pix[i] = (merge_val/merge_count)&0xFF;
    }
}

/*
    Deep color versions of the above functions
*/

static void depalettize_bitmap_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *unpacked_idx_buf,
    Py_buffer *unpacked_pal_buf, sint8 ucc)
{
    uint64 i=0, j=0, max_i=0;
    uint16 *unpacked_pix;
    uint16 *unpacked_pal;
    uint8 *unpacked_idx;

    unpacked_pix = (uint16*)unpacked_pix_buf->buf;
    unpacked_pal = (uint16*)unpacked_pal_buf->buf;
    unpacked_idx = (uint8*)unpacked_idx_buf->buf;
    max_i = unpacked_pix_buf->len;

    if (ucc == 4) {
        for (i=0; i < max_i; i += 4) {
            j = unpacked_idx[i>>2]<<2;
            unpacked_pix[i]   = unpacked_pal[j];
            unpacked_pix[i+1] = unpacked_pal[j+1];
            unpacked_pix[i+2] = unpacked_pal[j+2];
            unpacked_pix[i+3] = unpacked_pal[j+3];
        }
    } else if (ucc == 2) {
        for (i=0; i < max_i; i += 2) {
            j = unpacked_idx[i>>1]<<1;
            unpacked_pix[i]   = unpacked_pal[j];
            unpacked_pix[i+1] = unpacked_pal[j+1];
        }
    } else {
        for (i=0; i < max_i; i++) {
            unpacked_pix[i] = unpacked_pal[unpacked_idx[i]];
        }
    }
}

static void downsample_bitmap_16(
    Py_buffer *downsamp_pix_buf, Py_buffer *swizzled_pix_buf,
    uint32 merge_count, sint8 ucc)
{
    uint64 i=0, j=0, max_i=0, swizz_i=0;
    uint64 merge_val=0;
    uint16 *downsamp_pix;
    uint16 *swizzled_pix;

    downsamp_pix = (uint16*)downsamp_pix_buf->buf;
    swizzled_pix = (uint16*)swizzled_pix_buf->buf;
    max_i = downsamp_pix_buf->len;

    for (i=0; i < max_i; i++) {
        merge_val = 0;
        for (j=0; j < merge_count; j++) {
            merge_val += swizzled_pix[swizz_i];
            swizz_i++;
        }
        downsamp_pix[i] = (merge_val/merge_count)&0xFF;
    }
}


static PyObject *py_depalettize_bitmap(PyObject *self, PyObject *args) {
    Py_buffer bufs[3];
    sint8 ucc, i;

    // Get the pointers to each of the array objects and channel count
    if (!PyArg_ParseTuple(args, "w*w*w*b:depalettize_bitmap",
        &bufs[0], &bufs[1], &bufs[2], &ucc))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[0].itemsize == 2) {
        depalettize_bitmap_16(&bufs[0], &bufs[1], &bufs[2], ucc);
    } else {
        depalettize_bitmap_8(&bufs[0], &bufs[1], &bufs[2], ucc);
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_downsample_bitmap(PyObject *self, PyObject *args) {
    Py_buffer bufs[2];
    uint32 merge_count;
    sint8 ucc, i;

    // Get the pointers to each of the array objects and channel count
    if (!PyArg_ParseTuple(args, "w*w*kb:downsample_bitmap",
        &bufs[0], &bufs[1], &merge_count, &ucc)) {
        return Py_BuildValue("");  // return Py_None while incrementing it
    }

    if (bufs[0].itemsize == 2) {
        downsample_bitmap_16(&bufs[0], &bufs[1], merge_count, ucc);
    } else {
        downsample_bitmap_8(&bufs[0], &bufs[1], merge_count, ucc);
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_populate_scaler_array(PyObject *self, PyObject *args) {
    Py_buffer buf;
    double scale;
    uint64 i=0, max_i=0;
    uint8  *char_arr;
    uint16 *short_arr;

    if (!PyArg_ParseTuple(args, "w*d:populate_scaler_array", &buf, &scale)) {
        return Py_BuildValue("");  // return Py_None while incrementing it
    }

    max_i = buf.len;
    if (buf.itemsize == 2) {
        short_arr = (uint16 *)buf.buf;
        max_i /= 2;
        for (i=0; i < max_i; i++)
            short_arr[i] = (uint16)(i*scale + 0.5);
    } else {
        char_arr = (uint8  *)buf.buf;
        for (i=0; i < max_i; i++)
            char_arr[i] = (uint8)(i*scale + 0.5);
    }

    // Release the buffer objects
    PyBuffer_Release(&buf);

    return Py_BuildValue("");  // return Py_None while incrementing it
}

/* A list of all the methods defined by this module.
"METH_VARGS" tells Python how to call the handler.
The {NULL, NULL} entry indicates the end of the method definitions */
static PyMethodDef arbytmap_ext_methods[] = {
    {"depalettize_bitmap", py_depalettize_bitmap, METH_VARARGS, ""},
    {"downsample_bitmap", py_downsample_bitmap, METH_VARARGS, ""},
    {"populate_scaler_array", py_populate_scaler_array, METH_VARARGS, ""},
    //{"downsample_bitmap_gamma", py_downsample_bitmap_gamma, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}      /* sentinel */
};

/* When Python imports a C module named 'X' it loads the
module then looks for a method named "init"+X and calls it.*/
static struct PyModuleDef arbytmap_ext_module = {
    PyModuleDef_HEAD_INIT,
    "arbytmap_ext",
    "A set of C functions to replace certain speed intensive Arbytmap functions",
    -1,
    arbytmap_ext_methods,
};

PyMODINIT_FUNC PyInit_arbytmap_ext(void) {
    return PyModule_Create(&arbytmap_ext_module);
}
