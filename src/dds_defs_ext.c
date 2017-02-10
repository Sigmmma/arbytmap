#include <stdio.h>
#include <math.h>
#include "Python.h"
#include "abstract.h"
#include "longobject.h"
#include "modsupport.h"
#include "object.h"


static void unpack_dxt1_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *r_scale, Py_buffer *g_scale, Py_buffer *b_scale,
    char pix_per_tex, char chan0, char chan1, char chan2, char chan3,
    unsigned short unpack_max)
{
    Py_ssize_t packed_tex_size;
    unsigned char (*unpacked_pix)[];
    unsigned char (*r_scale)[], (*g_scale)[], (*b_scale)[];

    char ucc = 4;  // unpacked channel count
    char chans_per_tex = ucc*pix_per_tex;
    unsigned long long i=0, j=0, max_i=0, pxl_i=0;

    unpacked_pix = (unsigned char(*)[])unpacked_pix_buf->buf;
    a_scale = (unsigned char(*)[])a_scale_buf->buf;
    r_scale = (unsigned char(*)[])r_scale_buf->buf;
    g_scale = (unsigned char(*)[])g_scale_buf->buf;
    b_scale = (unsigned char(*)[])b_scale_buf->buf;

    packed_tex_size = packed_tex_buf->itemsize;
    packed_tex = (unsigned long(*)[]) packed_tex_buf->buf;

    //loop through each texel
    for i in range(len(packed_tex)/2):
        pxl_i = i*chans_per_tex
        j = i*2

        /*if the format DXT1 then the two entries in the array
        are the colors and the color indexing in that order.*/
        color0 = packed[j] & 0xFFff
        color1 = (packed[j] >> 16) & 0xFFff
        color_idx = packed[j+1]

        //unpack the colors
        c_0[1] = r_scale[(color0>>11) & 31]
        c_0[2] = g_scale[(color0>>5) & 63]
        c_0[3] = b_scale[(color0) & 31]

        c_1[1] = r_scale[(color1>>11) & 31]
        c_1[2] = g_scale[(color1>>5) & 63]
        c_1[3] = b_scale[(color1) & 31]

        //if the first color is a larger integer
        //then color key transparency is NOT used
        if color0 > color1:
            c_2[1] = (c_0[1]*2 + c_1[1])/3
            c_2[2] = (c_0[2]*2 + c_1[2])/3
            c_2[3] = (c_0[3]*2 + c_1[3])/3
            colors[3] = [
                unpack_max,
                (c_0[1] + 2*c_1[1])/3,
                (c_0[2] + 2*c_1[2])/3,
                (c_0[3] + 2*c_1[3])/3]
        else:
            c_2[1] = (c_0[1]+c_1[1])/2
            c_2[2] = (c_0[2]+c_1[2])/2
            c_2[3] = (c_0[3]+c_1[3])/2
            colors[3] = transparent

        for j in pixel_indices:
            color = colors[(color_idx >> (j*2))&3]
            off = j*ucc + pxl_i
            unpacked[off + chan0] = color[0]
            unpacked[off + chan1] = color[1]
            unpacked[off + chan2] = color[2]
            unpacked[off + chan3] = color[3]

}

static void unpack_dxt1_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *r_scale, Py_buffer *g_scale, Py_buffer *b_scale,
    char pixels_per_texel, char chan0, char chan1, char chan2, char chan3,
    unsigned short unpack_max)
{

}


static PyObject *py_unpack_dxt1(PyObject *self, PyObject *args) {
    Py_buffer bufs[5];
    char pix_per_tex;
    char chans[4];
    unsigned short unpack_max;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*bbbbbH:unpack_dxt1",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &pix_per_tex,
        &chans[0], &chans[1], &chans[2], &chans[3], &unpack_max))
        return Py_None;

    if (bufs[0].itemsize == 2) {
        unpack_dxt1_16(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], pix_per_tex,
            chans[0], chans[1], chans[2], chans[3], unpack_max);
    } else {
        unpack_dxt1_8(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], pix_per_tex,
            chans[0], chans[1], chans[2], chans[3], unpack_max);
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);

    return Py_None;
}

/* A list of all the methods defined by this module.
"METH_VARGS" tells Python how to call the handler.
The {NULL, NULL} entry indicates the end of the method definitions.*/
static PyMethodDef dds_defs_ext_methods[] = {
    {"unpack_dxt1", py_unpack_dxt1, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}      /* sentinel */
};

/* When Python imports a C module named 'X' it loads the
module then looks for a method named "init"+X and calls it.*/
static struct PyModuleDef dds_defs_ext_module = {
    PyModuleDef_HEAD_INIT,
    "dds_defs_ext",
    "A set of C functions to replace certain speed intensive dds functions",
    -1,
    dds_defs_ext_methods,
};

PyMODINIT_FUNC PyInit_dds_defs_ext(void) {
    return PyModule_Create(&dds_defs_ext_module);
}
