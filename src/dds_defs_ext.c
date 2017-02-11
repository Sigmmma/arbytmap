#include <stdio.h>
#include <math.h>
#include "Python.h"
#include "abstract.h"
#include "longobject.h"
#include "modsupport.h"
#include "object.h"

#define READ_DXT_COLORS(x, unpack_max)\
    color0 = (*packed_tex)[x]&0xFFff;\
    color1 = (*packed_tex)[x]>>16;\
    color_idx = (*packed_tex)[x+1];\
    \
    c_0[1] = (*r_scale)[(color0>>11) & 31];\
    c_0[2] = (*g_scale)[(color0>>5) & 63];\
    c_0[3] = (*b_scale)[(color0) & 31];\
    \
    c_1[1] = (*r_scale)[(color1>>11) & 31];\
    c_1[2] = (*g_scale)[(color1>>5) & 63];\
    c_1[3] = (*b_scale)[(color1) & 31];\
    \
    if (color0 > color1) {\
        c_2[1] = (c_0[1]*2 + c_1[1])/3;\
        c_2[2] = (c_0[2]*2 + c_1[2])/3;\
        c_2[3] = (c_0[3]*2 + c_1[3])/3;\
        \
        colors[3] = &c_3;\
        c_3[0] = unpack_max;\
        c_3[1] = (c_0[1] + 2*c_1[1])/3;\
        c_3[2] = (c_0[2] + 2*c_1[2])/3;\
        c_3[3] = (c_0[3] + 2*c_1[3])/3;\
    } else {\
        c_2[1] = (c_0[1]+c_1[1])/2;\
        c_2[2] = (c_0[2]+c_1[2])/2;\
        c_2[3] = (c_0[3]+c_1[3])/2;\
        colors[3] = &transparent;\
    }

#define UNPACK_DXT_COLORS()\
    color = colors[(color_idx>>(j<<1))&3];\
    off = j*ucc + pxl_i;\
    (*unpacked_pix)[off + chan1] = (*color)[1];\
    (*unpacked_pix)[off + chan2] = (*color)[2];\
    (*unpacked_pix)[off + chan3] = (*color)[3];

#define READ_DXT5_ALPHA()\
    alpha0 = (*packed_tex)[j]&0xFF;\
    alpha1 = ((*packed_tex)[j]>>8)&0xFF;\
    a_lookup[0] = alpha0;\
    a_lookup[1] = alpha1;\
    alpha_idx = (((unsigned long long)(*packed_tex)[j+1]<<16) +\
                 ((*packed_tex)[j]>>16));\
    if (alpha0 > alpha1) {\
        a_lookup[2] = (*a_scale)[(alpha0*6 + alpha1)/7];\
        a_lookup[3] = (*a_scale)[(alpha0*5 + alpha1*2)/7];\
        a_lookup[4] = (*a_scale)[(alpha0*4 + alpha1*3)/7];\
        a_lookup[5] = (*a_scale)[(alpha0*3 + alpha1*4)/7];\
        a_lookup[6] = (*a_scale)[(alpha0*2 + alpha1*5)/7];\
        a_lookup[7] = (*a_scale)[(alpha0   + alpha1*6)/7];\
    } else {\
        a_lookup[2] = (*a_scale)[(alpha0*4 + alpha1)/5];\
        a_lookup[3] = (*a_scale)[(alpha0*3 + alpha1*2)/5];\
        a_lookup[4] = (*a_scale)[(alpha0*2 + alpha1*3)/5];\
        a_lookup[5] = (*a_scale)[(alpha0   + alpha1*4)/5];\
        a_lookup[6] = (*a_scale)[0];\
        a_lookup[7] = (*a_scale)[255];\
    }


static void unpack_dxt1_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *r_scale_buf, Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    char pix_per_tex, char chan0, char chan1, char chan2, char chan3,
    unsigned char unpack_max)
{
    Py_ssize_t packed_tex_size = 4;
    unsigned char (*unpacked_pix)[];
    unsigned long (*packed_tex)[];
    unsigned char (*r_scale)[], (*g_scale)[], (*b_scale)[];
    unsigned short color0, color1;
    unsigned long color_idx;

    char ucc=4;  // unpacked channel count
    char chans_per_tex = ucc*pix_per_tex;
    unsigned long long i=0, j=0, max_i=0, pxl_i=0, off=0;
    unsigned char c_0[4]={unpack_max,0,0,0}, c_1[4]={unpack_max,0,0,0};
    unsigned char c_2[4]={unpack_max,0,0,0}, c_3[4]={unpack_max,0,0,0};
    unsigned char transparent[4]={0,0,0,0};
    unsigned char (*color)[4];
    unsigned char (*colors[4])[4];

    unpacked_pix = (unsigned char(*)[])unpacked_pix_buf->buf;
    r_scale = (unsigned char(*)[])r_scale_buf->buf;
    g_scale = (unsigned char(*)[])g_scale_buf->buf;
    b_scale = (unsigned char(*)[])b_scale_buf->buf;

    colors[0] = &c_0; colors[1] = &c_1; colors[2] = &c_2; colors[3] = &c_3;
    packed_tex = (unsigned long(*)[]) packed_tex_buf->buf;
    max_i = (unpacked_pix_buf->len)/chans_per_tex;

    //loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        j = i<<1;

        /*If the format DXT1 then the two entries in the array
        are the colors and the color indexing in that order.
        Also, if the first color is a larger integer
        then color key transparency is NOT used.*/
        READ_DXT_COLORS(j, unpack_max);

        for (j=0; j<pix_per_tex; j++) {
            UNPACK_DXT_COLORS();
            (*unpacked_pix)[off + chan0] = (*color)[0];
        }
    }
}


static void unpack_dxt2_3_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    char pix_per_tex, char chan0, char chan1, char chan2, char chan3)
{
    Py_ssize_t packed_tex_size = 4;
    unsigned char (*unpacked_pix)[];
    unsigned long (*packed_tex)[];
    unsigned char (*a_scale)[], (*r_scale)[], (*g_scale)[], (*b_scale)[];
    unsigned short color0, color1;
    unsigned long color_idx;

    char ucc=4;  // unpacked channel count
    char chans_per_tex = ucc*pix_per_tex;
    unsigned long long i=0, j=0, max_i=0, pxl_i=0, off=0, alpha=0;
    unsigned char c_0[4]={0,0,0,0}, c_1[4]={0,0,0,0};
    unsigned char c_2[4]={0,0,0,0}, c_3[4]={0,0,0,0};
    unsigned char transparent[4]={0,0,0,0};
    unsigned char (*color)[4];
    unsigned char (*colors[4])[4];

    unpacked_pix = (unsigned char(*)[])unpacked_pix_buf->buf;
    a_scale = (unsigned char(*)[])a_scale_buf->buf;
    r_scale = (unsigned char(*)[])r_scale_buf->buf;
    g_scale = (unsigned char(*)[])g_scale_buf->buf;
    b_scale = (unsigned char(*)[])b_scale_buf->buf;

    colors[0] = &c_0; colors[1] = &c_1; colors[2] = &c_2; colors[3] = &c_3;
    packed_tex = (unsigned long(*)[]) packed_tex_buf->buf;
    max_i = (unpacked_pix_buf->len)/chans_per_tex;

    //loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        j = i<<2;

        alpha = ((unsigned long long)(*packed_tex)[j+1]<<32) + (*packed_tex)[j];
        READ_DXT_COLORS(j+2, 0);

        for (j=0; j<pix_per_tex; j++) {
            UNPACK_DXT_COLORS();
            (*unpacked_pix)[off + chan0] = (*a_scale)[(alpha>>(j<<2))&15];
        }
    }
}


static void unpack_dxt4_5_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    char pix_per_tex, char chan0, char chan1, char chan2, char chan3)
{
    Py_ssize_t packed_tex_size = 4;
    unsigned char (*unpacked_pix)[];
    unsigned long (*packed_tex)[];
    unsigned char (*a_scale)[], (*r_scale)[], (*g_scale)[], (*b_scale)[];
    unsigned short color0, color1;
    unsigned long color_idx;

    char ucc=4;  // unpacked channel count
    char chans_per_tex = ucc*pix_per_tex;
    unsigned long long i=0, j=0, max_i=0, pxl_i=0, off=0, alpha_idx=0;
    unsigned char c_0[4]={0,0,0,0}, c_1[4]={0,0,0,0};
    unsigned char c_2[4]={0,0,0,0}, c_3[4]={0,0,0,0};
    unsigned char transparent[4]={0,0,0,0};
    unsigned char (*color)[4];
    unsigned char (*colors[4])[4];
    unsigned char alpha0=0, alpha1=0, a_lookup[8];

    unpacked_pix = (unsigned char(*)[])unpacked_pix_buf->buf;
    a_scale = (unsigned char(*)[])a_scale_buf->buf;
    r_scale = (unsigned char(*)[])r_scale_buf->buf;
    g_scale = (unsigned char(*)[])g_scale_buf->buf;
    b_scale = (unsigned char(*)[])b_scale_buf->buf;

    colors[0] = &c_0; colors[1] = &c_1; colors[2] = &c_2; colors[3] = &c_3;
    packed_tex = (unsigned long(*)[]) packed_tex_buf->buf;
    max_i = (unpacked_pix_buf->len)/chans_per_tex;

    //loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        j = i<<2;

        READ_DXT5_ALPHA();
        READ_DXT_COLORS(j+2, 0);

        for (j=0; j<pix_per_tex; j++) {
            UNPACK_DXT_COLORS();
            (*unpacked_pix)[off + chan0] = a_lookup[(alpha_idx>>(j<<3))&7];
        }
    }
}


/*
    Deep color versions of the above functions
*/


static void unpack_dxt1_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *r_scale_buf, Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    char pix_per_tex, char chan0, char chan1, char chan2, char chan3,
    unsigned short unpack_max)
{
    Py_ssize_t packed_tex_size = 4;
    unsigned short (*unpacked_pix)[];
    unsigned long (*packed_tex)[];
    unsigned short (*r_scale)[], (*g_scale)[], (*b_scale)[];
    unsigned short color0, color1;
    unsigned long color_idx;

    char ucc=4;  // unpacked channel count
    char chans_per_tex = ucc*pix_per_tex;
    unsigned long long i=0, j=0, max_i=0, pxl_i=0, off=0;
    unsigned short c_0[4]={unpack_max,0,0,0}, c_1[4]={unpack_max,0,0,0};
    unsigned short c_2[4]={unpack_max,0,0,0}, c_3[4]={unpack_max,0,0,0};
    unsigned short transparent[4]={0,0,0,0};
    unsigned short (*color)[4];
    unsigned short (*colors[4])[4];

    unpacked_pix = (unsigned short(*)[])unpacked_pix_buf->buf;
    r_scale = (unsigned short(*)[])r_scale_buf->buf;
    g_scale = (unsigned short(*)[])g_scale_buf->buf;
    b_scale = (unsigned short(*)[])b_scale_buf->buf;

    colors[0] = &c_0; colors[1] = &c_1; colors[2] = &c_2; colors[3] = &c_3;
    packed_tex = (unsigned long(*)[]) packed_tex_buf->buf;
    max_i = (unpacked_pix_buf->len)/chans_per_tex;

    //loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        j = i<<1;

        /*If the format DXT1 then the two entries in the array
        are the colors and the color indexing in that order.
        Also, if the first color is a larger integer
        then color key transparency is NOT used.*/
        READ_DXT_COLORS(j, unpack_max);

        for (j=0; j<pix_per_tex; j++) {
            UNPACK_DXT_COLORS();
            (*unpacked_pix)[off + chan0] = (*color)[0];
        }
    }
}

static void unpack_dxt2_3_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    char pix_per_tex, char chan0, char chan1, char chan2, char chan3)
{
    Py_ssize_t packed_tex_size = 4;
    unsigned short (*unpacked_pix)[];
    unsigned long (*packed_tex)[];
    unsigned short (*a_scale)[], (*r_scale)[], (*g_scale)[], (*b_scale)[];
    unsigned short color0, color1;
    unsigned long color_idx;

    char ucc=4;  // unpacked channel count
    char chans_per_tex = ucc*pix_per_tex;
    unsigned long long i=0, j=0, max_i=0, pxl_i=0, off=0, alpha=0;
    unsigned short c_0[4]={0,0,0,0}, c_1[4]={0,0,0,0};
    unsigned short c_2[4]={0,0,0,0}, c_3[4]={0,0,0,0};
    unsigned short transparent[4]={0,0,0,0};
    unsigned short (*color)[4];
    unsigned short (*colors[4])[4];

    unpacked_pix = (unsigned short(*)[])unpacked_pix_buf->buf;
    a_scale = (unsigned short(*)[])a_scale_buf->buf;
    r_scale = (unsigned short(*)[])r_scale_buf->buf;
    g_scale = (unsigned short(*)[])g_scale_buf->buf;
    b_scale = (unsigned short(*)[])b_scale_buf->buf;

    colors[0] = &c_0; colors[1] = &c_1; colors[2] = &c_2; colors[3] = &c_3;
    packed_tex = (unsigned long(*)[]) packed_tex_buf->buf;
    max_i = (unpacked_pix_buf->len)/chans_per_tex;

    //loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        j = i<<2;

        alpha = ((unsigned long long)(*packed_tex)[j+1]<<32) + (*packed_tex)[j];
        READ_DXT_COLORS(j+2, 0);

        for (j=0; j<pix_per_tex; j++) {
            UNPACK_DXT_COLORS();
            (*unpacked_pix)[off + chan0] = (*a_scale)[(alpha>>(j<<2))&15];
        }
    }
}


static void unpack_dxt4_5_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    char pix_per_tex, char chan0, char chan1, char chan2, char chan3)
{
    Py_ssize_t packed_tex_size = 4;
    unsigned short (*unpacked_pix)[];
    unsigned long (*packed_tex)[];
    unsigned short (*a_scale)[], (*r_scale)[], (*g_scale)[], (*b_scale)[];
    unsigned short color0, color1;
    unsigned long color_idx;

    char ucc=4;  // unpacked channel count
    char chans_per_tex = ucc*pix_per_tex;
    unsigned long long i=0, j=0, max_i=0, pxl_i=0, off=0, alpha_idx=0;
    unsigned short c_0[4]={0,0,0,0}, c_1[4]={0,0,0,0};
    unsigned short c_2[4]={0,0,0,0}, c_3[4]={0,0,0,0};
    unsigned short transparent[4]={0,0,0,0};
    unsigned short (*color)[4];
    unsigned short (*colors[4])[4];
    unsigned short alpha0=0, alpha1=0, a_lookup[8];

    unpacked_pix = (unsigned short(*)[])unpacked_pix_buf->buf;
    a_scale = (unsigned short(*)[])a_scale_buf->buf;
    r_scale = (unsigned short(*)[])r_scale_buf->buf;
    g_scale = (unsigned short(*)[])g_scale_buf->buf;
    b_scale = (unsigned short(*)[])b_scale_buf->buf;

    colors[0] = &c_0; colors[1] = &c_1; colors[2] = &c_2; colors[3] = &c_3;
    packed_tex = (unsigned long(*)[]) packed_tex_buf->buf;
    max_i = (unpacked_pix_buf->len)/chans_per_tex;

    //loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        j = i<<2;

        READ_DXT5_ALPHA();
        READ_DXT_COLORS(j+2, 0);

        for (j=0; j<pix_per_tex; j++) {
            UNPACK_DXT_COLORS();
            (*unpacked_pix)[off + chan0] = a_lookup[(alpha_idx>>(j<<3))&7];
        }
    }
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
            chans[0], chans[1], chans[2], chans[3], unpack_max&0xFF);
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]);
    PyBuffer_Release(&bufs[3]);
    PyBuffer_Release(&bufs[4]);

    return Py_None;
}


static PyObject *py_unpack_dxt2_3(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];
    char pix_per_tex;
    char chans[4];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*bbbbb:unpack_dxt2_3",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        &pix_per_tex, &chans[0], &chans[1], &chans[2], &chans[3]))
        return Py_None;

    if (bufs[0].itemsize == 2) {
        unpack_dxt2_3_16(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
            pix_per_tex, chans[0], chans[1], chans[2], chans[3]);
    } else {
        unpack_dxt2_3_8(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
            pix_per_tex, chans[0], chans[1], chans[2], chans[3]);
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]);
    PyBuffer_Release(&bufs[3]);
    PyBuffer_Release(&bufs[4]);
    PyBuffer_Release(&bufs[5]);

    return Py_None;
}


static PyObject *py_unpack_dxt4_5(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];
    char pix_per_tex;
    char chans[4];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*bbbbb:unpack_dxt4_5",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        &pix_per_tex, &chans[0], &chans[1], &chans[2], &chans[3]))
        return Py_None;

    if (bufs[0].itemsize == 2) {
        unpack_dxt4_5_16(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
            pix_per_tex, chans[0], chans[1], chans[2], chans[3]);
    } else {
        unpack_dxt4_5_8(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
            pix_per_tex, chans[0], chans[1], chans[2], chans[3]);
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]);
    PyBuffer_Release(&bufs[3]);
    PyBuffer_Release(&bufs[4]);
    PyBuffer_Release(&bufs[5]);

    return Py_None;
}

/* A list of all the methods defined by this module.
"METH_VARGS" tells Python how to call the handler.
The {NULL, NULL} entry indicates the end of the method definitions.*/
static PyMethodDef dds_defs_ext_methods[] = {
    {"unpack_dxt1", py_unpack_dxt1, METH_VARARGS, ""},
    {"unpack_dxt2_3", py_unpack_dxt2_3, METH_VARARGS, ""},
    {"unpack_dxt4_5", py_unpack_dxt4_5, METH_VARARGS, ""},
    //{"unpack_dxt5a", py_unpack_dxt5a, METH_VARARGS, ""},
    //{"unpack_dxn", py_unpack_dxn, METH_VARARGS, ""},
    //{"unpack_ctx1", py_unpack_ctx1, METH_VARARGS, ""},
    //{"unpack_u8v8", py_unpack_u8v8, METH_VARARGS, ""},
    //{"pack_dxt1", py_pack_dxt1, METH_VARARGS, ""},
    //{"pack_dxt2_3", py_pack_dxt2_3, METH_VARARGS, ""},
    //{"pack_dxt4_5", py_pack_dxt4_5, METH_VARARGS, ""},
    //{"pack_dxt5a", py_pack_dxt5a, METH_VARARGS, ""},
    //{"pack_dxn", py_pack_dxn, METH_VARARGS, ""},
    //{"pack_ctx1", py_pack_ctx1, METH_VARARGS, ""},
    //{"pack_u8v8", py_pack_u8v8, METH_VARARGS, ""},
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
