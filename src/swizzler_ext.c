#include "shared.h"

static void swizzle_char_array(
    Py_ssize_t c_size, sint8 *c_offs_sint8s,
    Py_ssize_t x_size, sint8 *x_offs_sint8s,
    Py_ssize_t y_size, sint8 *y_offs_sint8s,
    Py_ssize_t z_size, sint8 *z_offs_sint8s,
    Py_ssize_t swizz_size,   sint8 *swizz_arr_sint8s,
    Py_ssize_t unswizz_size, sint8 *unswizz_arr_sint8s,
    Py_ssize_t stride, sint8 swizz, char *method_name)
{
    sint32 i=0, zi=0, yi=0, xi=0, ci=0;
    uint32 z=0, y=0, x=0, c=0, yz=0, xyz=0, cxyz=0;
    uint32 max_i=0;
    uint32 *c_offs = (uint32 *)c_offs_sint8s;
    uint32 *x_offs = (uint32 *)x_offs_sint8s;
    uint32 *y_offs = (uint32 *)y_offs_sint8s;
    uint32 *z_offs = (uint32 *)z_offs_sint8s;
    uint8 *swizz_arr = (uint8*)swizz_arr_sint8s;
    uint8 *unswizz_arr = (uint8*)unswizz_arr_sint8s;
    c_size /= 4; x_size /= 4; y_size /= 4; z_size /= 4;
    swizz_size /= stride; unswizz_size /= stride;
    max_i = (uint32)(c_size * x_size * y_size * z_size);

    if ((max_i != swizz_size) || (max_i != unswizz_size)) {
        PySys_FormatStdout("Pixel buffers supplied to swizzler.%s are not the expected size.\n", method_name);
        return;
    }

    if (swizz) {
        if (stride == 8) {
            for (zi=0; zi < z_size; zi++) {
                z = z_offs[zi];
                for (yi=0; yi < y_size; yi++) {
                    yz = z + y_offs[yi];
                    for (xi=0; xi < x_size; xi++) {
                        xyz = yz + x_offs[xi];
                        for (ci=0; ci < c_size; ci++) {
                            cxyz = (xyz + c_offs[ci])<<3;
                            swizz_arr[cxyz] = unswizz_arr[i]; i++; cxyz++;
                            swizz_arr[cxyz] = unswizz_arr[i]; i++; cxyz++;
                            swizz_arr[cxyz] = unswizz_arr[i]; i++; cxyz++;
                            swizz_arr[cxyz] = unswizz_arr[i]; i++; cxyz++;
                            swizz_arr[cxyz] = unswizz_arr[i]; i++; cxyz++;
                            swizz_arr[cxyz] = unswizz_arr[i]; i++; cxyz++;
                            swizz_arr[cxyz] = unswizz_arr[i]; i++; cxyz++;
                            swizz_arr[cxyz] = unswizz_arr[i]; i++;
                        }
                    }
                }
            }
        } else if (stride == 4) {
            for (zi=0; zi < z_size; zi++) {
                z = z_offs[zi];
                for (yi=0; yi < y_size; yi++) {
                    yz = z + y_offs[yi];
                    for (xi=0; xi < x_size; xi++) {
                        xyz = yz + x_offs[xi];
                        for (ci=0; ci < c_size; ci++) {
                            cxyz = (xyz + c_offs[ci])<<2;
                            swizz_arr[cxyz] = unswizz_arr[i]; i++; cxyz++;
                            swizz_arr[cxyz] = unswizz_arr[i]; i++; cxyz++;
                            swizz_arr[cxyz] = unswizz_arr[i]; i++; cxyz++;
                            swizz_arr[cxyz] = unswizz_arr[i]; i++;
                        }
                    }
                }
            }
        } else if (stride == 2) {
            for (zi=0; zi < z_size; zi++) {
                z = z_offs[zi];
                for (yi=0; yi < y_size; yi++) {
                    yz = z + y_offs[yi];
                    for (xi=0; xi < x_size; xi++) {
                        xyz = yz + x_offs[xi];
                        for (ci=0; ci < c_size; ci++) {
                            cxyz = (xyz + c_offs[ci])<<1;
                            swizz_arr[cxyz] = unswizz_arr[i]; i++; cxyz++;
                            swizz_arr[cxyz] = unswizz_arr[i]; i++;
                        }
                    }
                }
            }
        } else {
            for (zi=0; zi < z_size; zi++) {
                z = z_offs[zi];
                for (yi=0; yi < y_size; yi++) {
                    yz = z + y_offs[yi];
                    for (xi=0; xi < x_size; xi++) {
                        xyz = yz + x_offs[xi];
                        for (ci=0; ci < c_size; ci++) {
                            swizz_arr[xyz + c_offs[ci]] = unswizz_arr[i];
                            i++;
                        }
                    }
                }
            }
        }
    } else {
        if (stride == 8) {
            for (zi=0; zi < z_size; zi++) {
                z = z_offs[zi];
                for (yi=0; yi < y_size; yi++) {
                    yz = z + y_offs[yi];
                    for (xi=0; xi < x_size; xi++) {
                        xyz = yz + x_offs[xi];
                        for (ci=0; ci < c_size; ci++) {
                            cxyz = (xyz + c_offs[ci])<<3;
                            swizz_arr[i] = unswizz_arr[cxyz]; i++; cxyz++;
                            swizz_arr[i] = unswizz_arr[cxyz]; i++; cxyz++;
                            swizz_arr[i] = unswizz_arr[cxyz]; i++; cxyz++;
                            swizz_arr[i] = unswizz_arr[cxyz]; i++; cxyz++;
                            swizz_arr[i] = unswizz_arr[cxyz]; i++; cxyz++;
                            swizz_arr[i] = unswizz_arr[cxyz]; i++; cxyz++;
                            swizz_arr[i] = unswizz_arr[cxyz]; i++; cxyz++;
                            swizz_arr[i] = unswizz_arr[cxyz]; i++;
                        }
                    }
                }
            }
        } else if (stride == 4) {
            for (zi=0; zi < z_size; zi++) {
                z = z_offs[zi];
                for (yi=0; yi < y_size; yi++) {
                    yz = z + y_offs[yi];
                    for (xi=0; xi < x_size; xi++) {
                        xyz = yz + x_offs[xi];
                        for (ci=0; ci < c_size; ci++) {
                            cxyz = (xyz + c_offs[ci])<<2;
                            swizz_arr[i] = unswizz_arr[cxyz]; i++; cxyz++;
                            swizz_arr[i] = unswizz_arr[cxyz]; i++; cxyz++;
                            swizz_arr[i] = unswizz_arr[cxyz]; i++; cxyz++;
                            swizz_arr[i] = unswizz_arr[cxyz]; i++;
                        }
                    }
                }
            }
        } else if (stride == 2) {
            for (zi=0; zi < z_size; zi++) {
                z = z_offs[zi];
                for (yi=0; yi < y_size; yi++) {
                    yz = z + y_offs[yi];
                    for (xi=0; xi < x_size; xi++) {
                        xyz = yz + x_offs[xi];
                        for (ci=0; ci < c_size; ci++) {
                            cxyz = (xyz + c_offs[ci])<<1;
                            swizz_arr[i] = unswizz_arr[cxyz]; i++; cxyz++;
                            swizz_arr[i] = unswizz_arr[cxyz]; i++;
                        }
                    }
                }
            }
        } else {
            for (zi=0; zi < z_size; zi++) {
                z = z_offs[zi];
                for (yi=0; yi < y_size; yi++) {
                    yz = z + y_offs[yi];
                    for (xi=0; xi < x_size; xi++) {
                        xyz = yz + x_offs[xi];
                        for (ci=0; ci < c_size; ci++) {
                            swizz_arr[i] = unswizz_arr[xyz + c_offs[ci]];
                            i++;
                        }
                    }
                }
            }
        }
    }
}


static PyObject *py_swizzle_array(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];
    sint8 i;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*:swizzle_array",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    swizzle_char_array(
        bufs[0].len, bufs[0].buf, bufs[1].len, bufs[1].buf,
        bufs[2].len, bufs[2].buf, bufs[3].len, bufs[3].buf,
        bufs[4].len, bufs[4].buf, bufs[5].len, bufs[5].buf,
        bufs[5].itemsize, 1, "swizzle_array");

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_unswizzle_array(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];
    sint8 i;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*:unswizzle_array",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    swizzle_char_array(
        bufs[0].len, bufs[0].buf, bufs[1].len, bufs[1].buf,
        bufs[2].len, bufs[2].buf, bufs[3].len, bufs[3].buf,
        bufs[4].len, bufs[4].buf, bufs[5].len, bufs[5].buf,
        bufs[5].itemsize, 0, "unswizzle_array");

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

/* A list of all the methods defined by this module.
"METH_VARGS" tells Python how to call the handler.
The {NULL, NULL} entry indicates the end of the method definitions */
static PyMethodDef swizzler_ext_methods[] = {
    {"swizzle_array", py_swizzle_array, METH_VARARGS, ""},
    {"unswizzle_array", py_unswizzle_array, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}      /* sentinel */
};

/* When Python imports a C module named 'X' it loads the
module then looks for a method named "init"+X and calls it.*/
static struct PyModuleDef swizzler_ext_module = {
    PyModuleDef_HEAD_INIT,
    "swizzler_ext",
    "A set of C functions to replace certain speed intensive swizzler functions",
    -1,
    swizzler_ext_methods,
};

PyMODINIT_FUNC PyInit_swizzler_ext(void) {
    return PyModule_Create(&swizzler_ext_module);
}
