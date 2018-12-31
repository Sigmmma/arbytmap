#include "shared.h"

static uint64 get_dxgi_tiled_address(uint32 offset, uint32 tile_ct_x, uint32 tile_len) 
{
    uint32 aligned_width = (tile_ct_x + 31) & 0xFFffFFe0;

    uint32 log_bpp = (tile_len >> 2) + ((tile_len >> 1) >> (tile_len >> 2));
    uint32 offset_b = offset << log_bpp;
    uint32 offset_t = (((offset_b & 0xFFffF000) >> 3) +
                       ((offset_b & 0x700) >> 2)) + (offset_b & 0x3F);
    uint32 offset_m = offset_t >> (7 + log_bpp);

    // x position values
    uint32 tile_x = (((offset_t >> (5 + log_bpp)) & 2) + (offset_b >> 6)) & 3;
    uint32 macro_x = (((offset_m % (aligned_width >> 5)) << 2) + tile_x) << 3;
    uint32 micro_x = ((((offset_t >> 1) & 0xFFFFFFF0) +
                      (offset_t & 0xF)) & ((tile_len << 3) - 1)) >> log_bpp;

    // y position values
    uint32 tile_y = ((offset_t >> (6 + log_bpp)) & 1) + ((offset_b & 0x800) >> 10);
    uint32 macro_y = (((offset_m / (aligned_width >> 5)) << 2) + tile_y) << 3;
    uint32 micro_y = (((offset_t & (((tile_len << 6) - 1) & 0xFFffFFe0)) +
                      ((offset_t & 15) << 1)) >> (3 + log_bpp)) & 0xFFffFFfe;

    uint64 x_address = (macro_x + micro_x) & 0xFFffFFff;
    uint64 y_address = macro_y + micro_y + ((offset_t & 0x10) >> 4);

    return x_address + (y_address << 32);
}


static PyObject *py_dxgi_tile_array(PyObject *self, PyObject *args) {
    Py_buffer bufs[2], *tiled_buf, *untiled_buf;
    int tile_mode;
    uint8 *tiled, *untiled;
    uint32 i, j, x, y;
    uint32 offset, t_idx, u_idx;
    uint32 x_chunks, y_chunks;
    uint32 width, height, depth;
    uint32 b_width, b_height, b_size;
    uint64 xy_address;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "pw*w*IIIIII:dxgi_tile_array",
        &tile_mode, &bufs[0], &bufs[1],
        &width, &height, &depth, &b_width, &b_height, &b_size)) {
        return Py_BuildValue("");  // return Py_None while incrementing it
    }

    if (tile_mode) {
        untiled_buf = &bufs[0];
        tiled_buf   = &bufs[1];
    } else {
        untiled_buf = &bufs[1];
        tiled_buf   = &bufs[0];
    }

    untiled = (uint8 *)untiled_buf->buf;
    tiled = (uint8 *)tiled_buf->buf;

    x_chunks = width / b_width;
    y_chunks = height / b_height;
    /*PySys_FormatStdout("UNTILED: %d  %d\n", untiled_buf->len, untiled_buf->itemsize);
    PySys_FormatStdout("  %d  %d  %d  %d\n", width, height, depth, b_size);
    PySys_FormatStdout("TILED: %d  %d\n", tiled_buf->len, tiled_buf->itemsize);
    PySys_FormatStdout("  %d  %d  %d  %d\n", b_width, b_height, depth, b_size);*/

    if (((x_chunks * y_chunks * depth * b_size) > bufs[0].len) ||
        ((x_chunks * y_chunks * depth * b_size) > bufs[1].len)) {
        RELEASE_PY_BUFFER_ARRAY(bufs, i)
        PySys_FormatStdout("Invalid offsets supplied to tiler_ext.dxgi_tile_array\n");
        return Py_BuildValue("");  // return Py_None while incrementing it
    }

    if (tile_mode) {
        for (i = 0; i < y_chunks; i++) {
            for (j = 0; j < x_chunks; j++) {
                offset = i * x_chunks + j;
                xy_address = get_dxgi_tiled_address(offset, x_chunks, b_size);
                x = (uint32)xy_address;
                y = (uint32)(xy_address >> 32);
                t_idx = offset * b_size;
                u_idx = (y * x_chunks + x) * b_size;
                if ((t_idx + b_size > tiled_buf->len) ||
                    (u_idx + b_size > untiled_buf->len)) {
                    PySys_FormatStdout("Calculated address outside buffer: x=%d\ty=%d\n", x, y);
                    i = y_chunks;
                    break;
                }
                memcpy(&tiled[t_idx], &untiled[u_idx], (size_t)b_size);
            }
        }
    } else {
        for (i = 0; i < y_chunks; i++) {
            for (j = 0; j < x_chunks; j++) {
                offset = i * x_chunks + j;
                xy_address = get_dxgi_tiled_address(offset, x_chunks, b_size);
                x = (uint32)xy_address;
                y = (uint32)(xy_address >> 32);
                t_idx = offset * b_size;
                u_idx = (y * x_chunks + x) * b_size;
                if ((t_idx + b_size > tiled_buf->len) ||
                    (u_idx + b_size > untiled_buf->len)) {
                    PySys_FormatStdout("Calculated address outside buffer: x=%d\ty=%d\n", x, y);
                    i = y_chunks;
                    break;
                }
                memcpy(&untiled[u_idx], &tiled[t_idx], (size_t)b_size);
            }
        }
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

/* A list of all the methods defined by this module.
"METH_VARGS" tells Python how to call the handler.
The {NULL, NULL} entry indicates the end of the method definitions */
static PyMethodDef tiler_ext_methods[] = {
    {"dxgi_tile_array", py_dxgi_tile_array, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}      /* sentinel */
};

/* When Python imports a C module named 'X' it loads the
module then looks for a method named "init"+X and calls it.*/
static struct PyModuleDef tiler_ext_module = {
    PyModuleDef_HEAD_INIT,
    "tiler_ext",
    "A set of C functions to replace certain speed intensive tiler functions",
    -1,
    tiler_ext_methods,
};

PyMODINIT_FUNC PyInit_tiler_ext(void) {
    return PyModule_Create(&tiler_ext_module);
}
