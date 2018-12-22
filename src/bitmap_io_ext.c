#include "shared.h"


static void pad_24bit_array(
    Py_buffer *padded_pix_buf, Py_buffer *unpadded_pix_buf)
{
    uint64 i=0, j=0, max_i=0;

    uint32 *padded_pix;
    uint8 *unpadded_pix;

    padded_pix   = (uint32*) padded_pix_buf->buf;
    unpadded_pix = (uint8*) unpadded_pix_buf->buf;

    max_i = (padded_pix_buf->len)/4;

    for (i=0; i < max_i; i++) {
        j = i*3;
        padded_pix[i] = (unpadded_pix[j] + (unpadded_pix[j+1]<<8) + (unpadded_pix[j+2]<<16));
    }
}

static void unpad_24bit_array(
    Py_buffer *unpadded_pix_buf, Py_buffer *padded_pix_buf)
{
    uint64 i=0, j=0, k=0, max_i=0;
    uint32 pixel=0;

    uint8 *unpadded_pix;
    uint8 *padded_pix_8;
    uint32 *padded_pix_32;

    unpadded_pix  = (uint8*) unpadded_pix_buf->buf;
    padded_pix_8  = (uint8*) padded_pix_buf->buf;
    padded_pix_32 = (uint32*) padded_pix_buf->buf;

    max_i = (unpadded_pix_buf->len)/3;

    if (padded_pix_buf->itemsize == 4) {
        for (i=0; i < max_i; i++) {
            j=i*3;
            pixel = padded_pix_32[i];
            unpadded_pix[j] = pixel&0xFF; j++;
            unpadded_pix[j] = (pixel>>8)&0xFF; j++;
            unpadded_pix[j] = (pixel>>16)&0xFF;
        }
    } else {
        for (i=0; i < max_i; i++) {
            /*
            If the itemsize isnt 4, the pixels havent been
            packed, and since the channel order in Arbytmap
            is ARGB, we need to skip each pixels first byte
            */
            j=i*3; k=i*4 + 1;
            unpadded_pix[j] = padded_pix_8[k]; j++; k++;
            unpadded_pix[j] = padded_pix_8[k]; j++; k++;
            unpadded_pix[j] = padded_pix_8[k];
        }
    }
}

static void pad_48bit_array(
    Py_buffer *padded_pix_buf, Py_buffer *unpadded_pix_buf)
{
    uint64 i=0, j=0, max_i=0;

    uint64 *padded_pix;
    uint16 *unpadded_pix;

    padded_pix   = (uint64*) padded_pix_buf->buf;
    unpadded_pix = (uint16*) unpadded_pix_buf->buf;

    max_i = (padded_pix_buf->len)/8;

    for (i=0; i < max_i; i++) {
        j = i*3;
        padded_pix[i] = ((uint64)unpadded_pix[j] +
                        ((uint64)unpadded_pix[j+1]<<16) +
                        ((uint64)unpadded_pix[j+2]<<32));
    }
}

static void unpad_48bit_array(
    Py_buffer *unpadded_pix_buf, Py_buffer *padded_pix_buf)
{
    uint64 i=0, j=0, k=0, max_i=0, pixel=0;

    uint16 *unpadded_pix;
    uint16 *padded_pix_16;
    uint64 *padded_pix_64;

    unpadded_pix  = (uint16*) unpadded_pix_buf->buf;
    padded_pix_16 = (uint16*) padded_pix_buf->buf;
    padded_pix_64 = (uint64*) padded_pix_buf->buf;

    max_i = (unpadded_pix_buf->len)/6;

    if (padded_pix_buf->itemsize == 8) {
        for (i=0; i < max_i; i++) {
            j=i*3;
            pixel = padded_pix_64[i];
            unpadded_pix[j] = pixel&0xFFff; j++;
            unpadded_pix[j] = (pixel>>16)&0xFFff; j++;
            unpadded_pix[j] = (pixel>>32)&0xFFff;
        }
    } else {
        for (i=0; i < max_i; i++) {
            /*
            If the itemsize isnt 8, the pixels havent been
            packed, and since the channel order in Arbytmap
            is ARGB, we need to skip each pixels first 2 bytes
            */
            j=i*3; k=i*4 + 1;
            unpadded_pix[j] = padded_pix_16[k]; j++; k++;
            unpadded_pix[j] = padded_pix_16[k]; j++; k++;
            unpadded_pix[j] = padded_pix_16[k];
        }
    }
}

static PyObject *py_pad_24bit_array(PyObject *self, PyObject *args) {
    Py_buffer bufs[2];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*:pad_24bit_array", &bufs[0], &bufs[1])) {
        return Py_BuildValue("");  // return Py_None while incrementing it
    }

    pad_24bit_array(&bufs[0], &bufs[1]);

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_pad_48bit_array(PyObject *self, PyObject *args) {
    Py_buffer bufs[2];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*:pad_48bit_array", &bufs[0], &bufs[1])) {
        return Py_BuildValue("");  // return Py_None while incrementing it
    }

    pad_48bit_array(&bufs[0], &bufs[1]);

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_unpad_24bit_array(PyObject *self, PyObject *args) {
    Py_buffer bufs[2];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*:unpad_24bit_array", &bufs[0], &bufs[1])) {
        return Py_BuildValue("");  // return Py_None while incrementing it
    }

    unpad_24bit_array(&bufs[0], &bufs[1]);

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_unpad_48bit_array(PyObject *self, PyObject *args) {
    Py_buffer bufs[2];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*:unpad_48bit_array", &bufs[0], &bufs[1])) {
        return Py_BuildValue("");  // return Py_None while incrementing it
    }

    unpad_48bit_array(&bufs[0], &bufs[1]);

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_crop_pixel_data(PyObject *self, PyObject *args) {
    Py_buffer bufs[2];
    uint8  *pixels, *new_pixels;
    sint64 z_stride, y_stride, x_stride, src_i, dst_i,
        src_y_skip0, dst_y_skip0, src_x_skip0, dst_x_skip0,
        src_y_skip1, dst_y_skip1, src_x_skip1, dst_x_skip1,
        max_src_i, max_dst_i, y, z;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*KKKKKKKKKKKKK:crop_pixel_data",
        &bufs[0], &bufs[1], &z_stride, &y_stride, &x_stride, &src_i, &dst_i,
        &src_y_skip0, &dst_y_skip0, &src_x_skip0, &dst_x_skip0,
        &src_y_skip1, &dst_y_skip1, &src_x_skip1, &dst_x_skip1)) {
        return Py_BuildValue("");  // return Py_None while incrementing it
    }

    src_x_skip1 += x_stride;
    dst_x_skip1 += x_stride;

    max_src_i = src_i + (src_y_skip0 + (src_x_skip0 + src_x_skip1) * y_stride + src_y_skip1) * z_stride;
    max_dst_i = dst_i + (dst_y_skip0 + (dst_x_skip0 + dst_x_skip1) * y_stride + dst_y_skip1) * z_stride;

    /*PySys_FormatStdout("%d  %d    %d  %d\n", max_src_i, max_dst_i, 
                       bufs[0].len, bufs[1].len);*/
    if (max_src_i < bufs[0].len || max_dst_i < bufs[1].len) {
        PyBuffer_Release(&bufs[0]);
        PyBuffer_Release(&bufs[1]);
        PySys_FormatStdout("Invalid offsets supplied to bitmap_io_ext.crop_pixel_data\n");
        return Py_BuildValue("");  // return Py_None while incrementing it
    }

    pixels = (uint8*)bufs[0].buf;
    new_pixels = (uint8*)bufs[1].buf;
    for (z = 0; z < z_stride; z++) {
        src_i += src_y_skip0;
        dst_i += dst_y_skip0;
        for (y = 0; y < y_stride; y++) {
            src_i += src_x_skip0;
            dst_i += dst_x_skip0;
            /*if (src_i + x_stride > bufs[0].len || dst_i + x_stride > bufs[1].len) {
            PySys_FormatStdout("FUCK\n");
            break;
            }*/
            memcpy(&new_pixels[dst_i], &pixels[src_i], x_stride);
            src_i += src_x_skip1;
            dst_i += dst_x_skip1;
        }
        src_i += src_y_skip1;
        dst_i += dst_y_skip1;
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_swap_array_items(PyObject *self, PyObject *args) {
    Py_buffer bufs[2];
    uint8  *pixels, *temp_pix;
    sint16 *chan_map, step;
    uint64 max_i, i, j;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*:swap_array_items", &bufs[0], &bufs[1])) {
        return Py_BuildValue("");  // return Py_None while incrementing it
    }

    pixels   = (uint8  *)bufs[0].buf;
    chan_map = (sint16 *)bufs[1].buf;
    max_i = bufs[0].len;
    step  = (sint16)(bufs[1].len / bufs[1].itemsize);
    temp_pix = malloc(step);
    if (temp_pix == NULL) {
        PyBuffer_Release(&bufs[0]);
        PyBuffer_Release(&bufs[1]);
        PySys_FormatStdout("Could not allocate pixel swap buffer in bitmap_io_ext.swap_array_items\n");
        return Py_BuildValue("");  // return Py_None while incrementing it
    }

    for (i = 0; i < max_i; i += step) {
        memcpy(temp_pix, &(pixels[i]), step);
        for (j = 0; j < step; j++) {
            if (chan_map[j] < 0) {
                pixels[i + j] = 0;
            } else {
                pixels[i + j] = temp_pix[chan_map[j]];
            }
        }
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);
    free(temp_pix);

    return Py_BuildValue("");  // return Py_None while incrementing it
}

// This was GOING to be a function to quickly make an array object and
// return it to the caller without having to copy it. Guess I cant do that.
/*
static PyObject *py_make_array(PyObject *self, PyObject *args) {
    PyObject *obj;
    Py_buffer array_buffer;
    Py_ssize_t len, itemsize;
    sint32 format_len;
    sint8 *format;
    sint8 *buf;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "u#k:make_array", format, &format_len, &len))
        return Py_BuildValue("");  // return Py_None while incrementing it

    buf = calloc(len, 1);
    if (format_len != 1 || buf == NULL) {
        return Py_None;
    }

    switch (*format) {
        case 'b': case 'B':
            itemsize = sizeof(sint8);
            break;
        case 'h': case 'H':
            itemsize = sizeof(sint16);
            break;
        case 'u':
            itemsize = sizeof(Py_UNICODE);
            break;
        case 'i': case 'I':
            itemsize = sizeof(sint32);
            break;
        case 'l': case 'L':
            itemsize = sizeof(sint32);
            break;
        case 'q': case 'Q':
            itemsize = sizeof(sint64);
            break;
        case 'f':
            itemsize = sizeof(float);
            break;
        case 'd':
            itemsize = sizeof(double);
            break;
    }

    array_buffer.buf = buf;
    array_buffer.obj = obj;
    array_buffer.len = len;
    array_buffer.readonly = (sint32)0;
    array_buffer.itemsize = itemsize;
    array_buffer.format = format;

    // need to somehow turn obj into a valid PyObject

    return obj;
}
*/
/* A list of all the methods defined by this module.
"METH_VARGS" tells Python how to call the handler.
The {NULL, NULL} entry indicates the end of the method definitions */
static PyMethodDef bitmap_io_ext_methods[] = {
    {"pad_24bit_array", py_pad_24bit_array, METH_VARARGS, ""},
    {"pad_48bit_array", py_pad_48bit_array, METH_VARARGS, ""},
    {"unpad_24bit_array", py_unpad_24bit_array, METH_VARARGS, ""},
    {"unpad_48bit_array", py_unpad_48bit_array, METH_VARARGS, ""},
    {"swap_array_items", py_swap_array_items, METH_VARARGS, ""},
    {"crop_pixel_data", py_crop_pixel_data, METH_VARARGS, ""},
    //{"make_array", py_make_array, METH_VARARGS, ""},
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
