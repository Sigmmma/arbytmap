#include "shared.h"

//short-hand macros for packer routines

//unpacking individual components
#define PACK_ARGB_A_MERGE(typ) ((typ)a_scale[(unpacked_pix[i<<2]+a_rnd)/a_div]<<a_shift)
#define PACK_ARGB_R_MERGE(typ) ((typ)r_scale[(unpacked_pix[(i<<2)+1]+r_rnd)/r_div]<<r_shift)
#define PACK_ARGB_G_MERGE(typ) ((typ)g_scale[(unpacked_pix[(i<<2)+2]+g_rnd)/g_div]<<g_shift)
#define PACK_ARGB_B_MERGE(typ) ((typ)b_scale[(unpacked_pix[(i<<2)+3]+b_rnd)/b_div]<<b_shift)

#define PACK_AI_A_MERGE(typ)   ((typ)a_scale[(unpacked_pix[i<<1]+a_rnd)/a_div]<<a_shift)
#define PACK_AI_I_MERGE(typ)   ((typ)i_scale[(unpacked_pix[(i<<1)+1]+i_rnd)/i_div]<<i_shift)

#define PACK_ARGB_A(typ) ((typ)a_scale[unpacked_pix[(i<<2)]]<<a_shift)
#define PACK_ARGB_R(typ) ((typ)r_scale[unpacked_pix[(i<<2)+1]]<<r_shift)
#define PACK_ARGB_G(typ) ((typ)g_scale[unpacked_pix[(i<<2)+2]]<<g_shift)
#define PACK_ARGB_B(typ) ((typ)b_scale[unpacked_pix[(i<<2)+3]]<<b_shift)

#define PACK_AI_A(typ)   ((typ)a_scale[unpacked_pix[(i<<1)]]<<a_shift)
#define PACK_AI_I(typ)   ((typ)i_scale[unpacked_pix[(i<<1)+1]]<<i_shift)


#define PACK_ARGB_A_MERGE_16(typ) ((typ)a_scale_16[(unpacked_pix[i<<2]+a_rnd)/a_div]<<a_shift)
#define PACK_ARGB_R_MERGE_16(typ) ((typ)r_scale_16[(unpacked_pix[(i<<2)+1]+r_rnd)/r_div]<<r_shift)
#define PACK_ARGB_G_MERGE_16(typ) ((typ)g_scale_16[(unpacked_pix[(i<<2)+2]+g_rnd)/g_div]<<g_shift)
#define PACK_ARGB_B_MERGE_16(typ) ((typ)b_scale_16[(unpacked_pix[(i<<2)+3]+b_rnd)/b_div]<<b_shift)

#define PACK_AI_A_MERGE_16(typ)   ((typ)a_scale_16[(unpacked_pix[i<<1]+a_rnd)/a_div]<<a_shift)
#define PACK_AI_I_MERGE_16(typ)   ((typ)i_scale_16[(unpacked_pix[(i<<1)+1]+i_rnd)/i_div]<<i_shift)

#define PACK_ARGB_A_16(typ) ((typ)a_scale_16[unpacked_pix[(i<<2)]]<<a_shift)
#define PACK_ARGB_R_16(typ) ((typ)r_scale_16[unpacked_pix[(i<<2)+1]]<<r_shift)
#define PACK_ARGB_G_16(typ) ((typ)g_scale_16[unpacked_pix[(i<<2)+2]]<<g_shift)
#define PACK_ARGB_B_16(typ) ((typ)b_scale_16[unpacked_pix[(i<<2)+3]]<<b_shift)

#define PACK_AI_A_16(typ)   ((typ)a_scale_16[unpacked_pix[(i<<1)]]<<a_shift)
#define PACK_AI_I_16(typ)   ((typ)i_scale_16[unpacked_pix[(i<<1)+1]]<<i_shift)

//unpacking ARGB and AI in different ways
#define PACK_ARGB_MERGE(typ) (\
PACK_ARGB_A_MERGE(typ) + PACK_ARGB_R_MERGE(typ) + PACK_ARGB_G_MERGE(typ) + PACK_ARGB_B_MERGE(typ))
#define PACK_AI_MERGE(typ) (PACK_AI_A_MERGE(typ) + PACK_AI_I_MERGE(typ))

#define PACK_ARGB(typ) (PACK_ARGB_A(typ) | PACK_ARGB_R(typ) | PACK_ARGB_G(typ) | PACK_ARGB_B(typ))
#define PACK_AI(typ) (PACK_AI_A(typ) | PACK_AI_I(typ))

//deep color versions of the above
#define PACK_ARGB_MERGE_16(typ) (\
PACK_ARGB_A_MERGE_16(typ) + PACK_ARGB_R_MERGE_16(typ) +\
PACK_ARGB_G_MERGE_16(typ) + PACK_ARGB_B_MERGE_16(typ))
#define PACK_AI_MERGE_16(typ) (PACK_AI_A_MERGE_16(typ) + PACK_AI_I_MERGE_16(typ))

#define PACK_ARGB_16(typ) (\
PACK_ARGB_A_16(typ) | PACK_ARGB_R_16(typ) | PACK_ARGB_G_16(typ) | PACK_ARGB_B_16(typ))
#define PACK_AI_16(typ) (PACK_AI_A_16(typ) | PACK_AI_I_16(typ))


static void pack_raw_4_channel_merge_8bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint64 a_div, sint64 r_div, sint64 g_div, sint64 b_div,
    sint8 a_shift,  sint8 r_shift,  sint8 g_shift,  sint8 b_shift)
{
    Py_ssize_t packed_pix_size;
    uint8 *unpacked_pix;
    uint8 *a_scale, *r_scale, *g_scale, *b_scale;
    uint64 i=0, max_i=0;
    sint64 a_rnd, r_rnd, g_rnd, b_rnd;

    uint8  *packed_pix_8;
    uint16 *packed_pix_16;
    uint32 *packed_pix_32;
    uint64 *packed_pix_64;

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8 = (uint8*)packed_pix_buf->buf;
    packed_pix_16 = (uint16*)packed_pix_buf->buf;
    packed_pix_32 = (uint32*)packed_pix_buf->buf;
    packed_pix_64 = (uint64*)packed_pix_buf->buf;

    unpacked_pix = (uint8*)unpacked_pix_buf->buf;
    a_scale = (uint8*)a_scale_buf->buf;
    r_scale = (uint8*)r_scale_buf->buf;
    g_scale = (uint8*)g_scale_buf->buf;
    b_scale = (uint8*)b_scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;
    a_rnd=a_div/2; r_rnd=r_div/2; g_rnd=g_div/2; b_rnd=b_div/2;

    if (packed_pix_size == 8) {
        for (i=0; i < max_i; i++) {
            packed_pix_64[i] = PACK_ARGB_MERGE(uint64)&0xFFffFFffFFffFFffULL;
        }
    } else if (packed_pix_size == 4) {
        for (i=0; i < max_i; i++) {
            packed_pix_32[i] = PACK_ARGB_MERGE(uint32)&0xFFffFFffUL;
        }
    } else if (packed_pix_size == 2) {
        for (i=0; i < max_i; i++) {
            packed_pix_16[i] = PACK_ARGB_MERGE(uint16)&0xFFff;
        }
    } else {
        for (i=0; i < max_i; i++) {
            packed_pix_8[i] = PACK_ARGB_MERGE(uint8)&0xFF;
        }
    }
}

static void pack_raw_2_channel_merge_8bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *i_scale_buf,
    sint64 a_div, sint64 i_div, sint8 a_shift,  sint8 i_shift)
{
    Py_ssize_t packed_pix_size;
    uint8 *unpacked_pix;
    uint8 *a_scale, *i_scale;
    uint64 i=0, max_i=0;
    sint64 a_rnd=a_div/2, i_rnd=i_div/2;

    uint8  *packed_pix_8;
    uint16 *packed_pix_16;
    uint32 *packed_pix_32;
    uint64 *packed_pix_64;

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8 = (uint8*)packed_pix_buf->buf;
    packed_pix_16 = (uint16*)packed_pix_buf->buf;
    packed_pix_32 = (uint32*)packed_pix_buf->buf;
    packed_pix_64 = (uint64*)packed_pix_buf->buf;

    unpacked_pix = (uint8*)unpacked_pix_buf->buf;
    a_scale = (uint8*)a_scale_buf->buf;
    i_scale = (uint8*)i_scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 8) {
        for (i=0; i < max_i; i++) {
            packed_pix_64[i] = PACK_AI_MERGE(uint64)&0xFFffFFffFFffFFffULL;
        }
    } else if (packed_pix_size == 4) {
        for (i=0; i < max_i; i++) {
            packed_pix_32[i] = PACK_AI_MERGE(uint32)&0xFFffFFffUL;
        }
    } else if (packed_pix_size == 2) {
        for (i=0; i < max_i; i++) {
            packed_pix_16[i] = PACK_AI_MERGE(uint16)&0xFFff;
        }
    } else {
        for (i=0; i < max_i; i++) {
            packed_pix_8[i] = PACK_AI_MERGE(uint8)&0xFF;
        }
    }
}

static void pack_raw_4_channel_8bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 a_shift,  sint8 r_shift,  sint8 g_shift,  sint8 b_shift)
{
    Py_ssize_t packed_pix_size;
    uint8 *unpacked_pix;
    uint8 *a_scale, *r_scale, *g_scale, *b_scale;
    uint64 i=0, max_i=0;

    uint8  *packed_pix_8;
    uint16 *packed_pix_16;
    uint32 *packed_pix_32;
    uint64 *packed_pix_64;

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8 = (uint8*)packed_pix_buf->buf;
    packed_pix_16 = (uint16*)packed_pix_buf->buf;
    packed_pix_32 = (uint32*)packed_pix_buf->buf;
    packed_pix_64 = (uint64*)packed_pix_buf->buf;

    unpacked_pix = (uint8*)unpacked_pix_buf->buf;
    a_scale = (uint8*)a_scale_buf->buf;
    r_scale = (uint8*)r_scale_buf->buf;
    g_scale = (uint8*)g_scale_buf->buf;
    b_scale = (uint8*)b_scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 8) {
        for (i=0; i < max_i; i++) {
            packed_pix_64[i] = PACK_ARGB(uint64)&0xFFffFFffFFffFFffULL;
        }
    } else if (packed_pix_size == 4) {
        for (i=0; i < max_i; i++) {
            packed_pix_32[i] = PACK_ARGB(uint32)&0xFFffFFffUL;
        }
    } else if (packed_pix_size == 2) {
        for (i=0; i < max_i; i++) {
            packed_pix_16[i] = PACK_ARGB(uint16)&0xFFff;
        }
    } else {
        for (i=0; i < max_i; i++) {
            packed_pix_8[i] = PACK_ARGB(uint8)&0xFF;
        }
    }
}

static void pack_raw_2_channel_8bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *i_scale_buf,
    sint8 a_shift, sint8 i_shift)
{
    Py_ssize_t packed_pix_size;
    uint8 *unpacked_pix;
    uint8 *a_scale, *i_scale;
    uint64 i=0, max_i=0;

    uint8  *packed_pix_8;
    uint16 *packed_pix_16;
    uint32 *packed_pix_32;
    uint64 *packed_pix_64;

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8 = (uint8*)packed_pix_buf->buf;
    packed_pix_16 = (uint16*)packed_pix_buf->buf;
    packed_pix_32 = (uint32*)packed_pix_buf->buf;
    packed_pix_64 = (uint64*)packed_pix_buf->buf;

    unpacked_pix = (uint8*)unpacked_pix_buf->buf;
    a_scale = (uint8*)a_scale_buf->buf;
    i_scale = (uint8*)i_scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 8) {
        for (i=0; i < max_i; i++) {
            packed_pix_64[i] = PACK_AI(uint64)&0xFFffFFffFFffFFffULL;
        }
    } else if (packed_pix_size == 4) {
        for (i=0; i < max_i; i++) {
            packed_pix_32[i] = PACK_AI(uint32)&0xFFffFFffUL;
        }
    } else if (packed_pix_size == 2) {
        for (i=0; i < max_i; i++) {
            packed_pix_16[i] = PACK_AI(uint16)&0xFFff;
        }
    } else {
        for (i=0; i < max_i; i++) {
            packed_pix_8[i] = PACK_AI(uint8)&0xFF;
        }
    }
}

static void pack_raw_1_channel_8bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *scale_buf, sint8 shift)
{
    Py_ssize_t packed_pix_size;
    uint8 *unpacked_pix;
    uint8 *scale;
    uint64 i=0, max_i=0;

    uint8  *packed_pix_8;
    uint16 *packed_pix_16;
    uint32 *packed_pix_32;
    uint64 *packed_pix_64;

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8 = (uint8*)packed_pix_buf->buf;
    packed_pix_16 = (uint16*)packed_pix_buf->buf;
    packed_pix_32 = (uint32*)packed_pix_buf->buf;
    packed_pix_64 = (uint64*)packed_pix_buf->buf;

    unpacked_pix = (uint8*)unpacked_pix_buf->buf;
    scale = (uint8*)scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 8) {
        for (i=0; i < max_i; i++) {
            packed_pix_64[i] = (uint64)(scale[unpacked_pix[i]]<<shift);
        }
    } else if (packed_pix_size == 4) {
        for (i=0; i < max_i; i++) {
            packed_pix_32[i] = scale[unpacked_pix[i]]<<shift;
        }
    } else if (packed_pix_size == 2) {
        for (i=0; i < max_i; i++) {
            packed_pix_16[i] = scale[unpacked_pix[i]]<<shift;
        }
    } else {
        for (i=0; i < max_i; i++) {
            packed_pix_8[i] = scale[unpacked_pix[i]]<<shift;
        }
    }
}


/*
    Deep color versions of the above functions
*/





static void pack_raw_4_channel_merge_16bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint64 a_div, sint64 r_div, sint64 g_div, sint64 b_div,
    sint8 a_shift,  sint8 r_shift,  sint8 g_shift,  sint8 b_shift)
{
    Py_ssize_t packed_pix_size;
    uint16 *unpacked_pix;
    uint8 *a_scale, *r_scale, *g_scale, *b_scale;
    uint16 *a_scale_16, *r_scale_16, *g_scale_16, *b_scale_16;
    uint64 i=0, max_i=0;
    sint64 a_rnd=a_div/2, r_rnd=r_div/2, g_rnd=g_div/2, b_rnd=b_div/2;

    uint8  *packed_pix_8;
    uint16 *packed_pix_16;
    uint32 *packed_pix_32;
    uint64 *packed_pix_64;

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8 = (uint8*)packed_pix_buf->buf;
    packed_pix_16 = (uint16*)packed_pix_buf->buf;
    packed_pix_32 = (uint32*)packed_pix_buf->buf;
    packed_pix_64 = (uint64*)packed_pix_buf->buf;

    unpacked_pix = (uint16*)unpacked_pix_buf->buf;
    a_scale = (uint8*)a_scale_buf->buf;
    r_scale = (uint8*)r_scale_buf->buf;
    g_scale = (uint8*)g_scale_buf->buf;
    b_scale = (uint8*)b_scale_buf->buf;
    a_scale_16 = (uint16*)a_scale_buf->buf;
    r_scale_16 = (uint16*)r_scale_buf->buf;
    g_scale_16 = (uint16*)g_scale_buf->buf;
    b_scale_16 = (uint16*)b_scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;
    if (packed_pix_size == 8) {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
                packed_pix_64[i] = PACK_ARGB_MERGE_16(uint64)&0xFFffFFffFFffFFffULL;
            }
        } else {
            for (i=0; i < max_i; i++) {
                packed_pix_64[i] = PACK_ARGB_MERGE(uint64)&0xFFffFFffFFffFFffULL;
            }
        }
    } else if (packed_pix_size == 4) {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
                packed_pix_32[i] = PACK_ARGB_MERGE_16(uint32)&0xFFffFFffUL;
            }
        } else {
            for (i=0; i < max_i; i++) {
                packed_pix_32[i] = PACK_ARGB_MERGE(uint32)&0xFFffFFffUL;
            }
        }
    } else if (packed_pix_size == 2) {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
                packed_pix_16[i] = PACK_ARGB_MERGE_16(uint16)&0xFFff;
            }
        } else {
            for (i=0; i < max_i; i++) {
                packed_pix_16[i] = PACK_ARGB_MERGE(uint16)&0xFFff;
            }
        }
    } else {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
                packed_pix_8[i] = PACK_ARGB_MERGE_16(uint8)&0xFF;
            }
        } else {
            for (i=0; i < max_i; i++) {
                packed_pix_8[i] = PACK_ARGB_MERGE(uint8)&0xFF;
            }
        }
    }
}

static void pack_raw_2_channel_merge_16bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *i_scale_buf,
    sint64 a_div, sint64 i_div, sint8 a_shift,  sint8 i_shift)
{
    Py_ssize_t packed_pix_size;
    uint16 *unpacked_pix;
    uint8  *a_scale,  *i_scale;
    uint16 *a_scale_16, *i_scale_16;
    uint64 i=0, max_i=0;
    sint64 a_rnd=a_div/2, i_rnd=i_div/2;

    uint8  *packed_pix_8;
    uint16 *packed_pix_16;
    uint32 *packed_pix_32;
    uint64 *packed_pix_64;

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8 = (uint8*)packed_pix_buf->buf;
    packed_pix_16 = (uint16*)packed_pix_buf->buf;
    packed_pix_32 = (uint32*)packed_pix_buf->buf;
    packed_pix_64 = (uint64*)packed_pix_buf->buf;

    unpacked_pix = (uint16*)unpacked_pix_buf->buf;
    a_scale = (uint8*)a_scale_buf->buf;
    i_scale = (uint8*)i_scale_buf->buf;
    a_scale_16 = (uint16*)a_scale_buf->buf;
    i_scale_16 = (uint16*)i_scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 8) {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
                packed_pix_64[i] = PACK_AI_MERGE_16(uint64)&0xFFffFFffFFffFFffULL;
            }
        } else {
            for (i=0; i < max_i; i++) {
                packed_pix_64[i] = PACK_AI_MERGE(uint64)&0xFFffFFffFFffFFffULL;
            }
        }
    } else if (packed_pix_size == 4) {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
                packed_pix_32[i] = PACK_AI_MERGE_16(uint32)&0xFFffFFffUL;
            }
        } else {
            for (i=0; i < max_i; i++) {
                packed_pix_32[i] = PACK_AI_MERGE(uint32)&0xFFffFFffUL;
            }
        }
    } else if (packed_pix_size == 2) {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
                packed_pix_16[i] = PACK_AI_MERGE_16(uint16)&0xFFff;
            }
        } else {
            for (i=0; i < max_i; i++) {
                packed_pix_16[i] = PACK_AI_MERGE(uint16)&0xFFff;
            }
        }
    } else {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
                packed_pix_8[i] = PACK_AI_MERGE_16(uint8)&0xFF;
            }
        } else {
            for (i=0; i < max_i; i++) {
                packed_pix_8[i] = PACK_AI_MERGE(uint8)&0xFF;
            }
        }
    }
}


static void pack_raw_4_channel_16bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 a_shift,  sint8 r_shift,  sint8 g_shift,  sint8 b_shift)
{
    Py_ssize_t packed_pix_size;
    uint16 *unpacked_pix;
    uint8 *a_scale, *r_scale, *g_scale, *b_scale;
    uint16 *a_scale_16, *r_scale_16, *g_scale_16, *b_scale_16;
    uint64 i=0, max_i=0;

    uint8  *packed_pix_8;
    uint16 *packed_pix_16;
    uint32 *packed_pix_32;
    uint64 *packed_pix_64;

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8 = (uint8*)packed_pix_buf->buf;
    packed_pix_16 = (uint16*)packed_pix_buf->buf;
    packed_pix_32 = (uint32*)packed_pix_buf->buf;
    packed_pix_64 = (uint64*)packed_pix_buf->buf;

    a_scale = (uint8*)a_scale_buf->buf;
    r_scale = (uint8*)r_scale_buf->buf;
    g_scale = (uint8*)g_scale_buf->buf;
    b_scale = (uint8*)b_scale_buf->buf;
    a_scale_16 = (uint16*)a_scale_buf->buf;
    r_scale_16 = (uint16*)r_scale_buf->buf;
    g_scale_16 = (uint16*)g_scale_buf->buf;
    b_scale_16 = (uint16*)b_scale_buf->buf;
    unpacked_pix = (uint16*)unpacked_pix_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 8) {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
                packed_pix_64[i] = PACK_ARGB_16(uint64)&0xFFffFFffFFffFFffULL;
            }
        } else {
            for (i=0; i < max_i; i++) {
                packed_pix_64[i] = PACK_ARGB(uint64)&0xFFffFFffFFffFFffULL;
            }
        }
    } else if (packed_pix_size == 4) {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
                packed_pix_32[i] = PACK_ARGB_16(uint32)&0xFFffFFffUL;
            }
        } else {
            for (i=0; i < max_i; i++) {
                packed_pix_32[i] = PACK_ARGB(uint32)&0xFFffFFffUL;
            }
        }
    } else if (packed_pix_size == 2) {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
                packed_pix_16[i] = PACK_ARGB_16(uint16)&0xFFff;
            }
        } else {
            for (i=0; i < max_i; i++) {
                packed_pix_16[i] = PACK_ARGB(uint16)&0xFFff;
            }
        }
    } else {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
                packed_pix_8[i] = PACK_ARGB_16(uint8)&0xFFff;
            }
        } else {
            for (i=0; i < max_i; i++) {
                packed_pix_8[i] = PACK_ARGB(uint8)&0xFF;
            }
        }
    }
}

static void pack_raw_2_channel_16bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *i_scale_buf,
    sint8 a_shift, sint8 i_shift)
{
    Py_ssize_t packed_pix_size;
    uint16 *unpacked_pix;
    uint8  *a_scale,  *i_scale;
    uint16 *a_scale_16, *i_scale_16;
    uint64 i=0, max_i=0;

    uint8  *packed_pix_8;
    uint16 *packed_pix_16;
    uint32 *packed_pix_32;
    uint64 *packed_pix_64;

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8 = (uint8*)packed_pix_buf->buf;
    packed_pix_16 = (uint16*)packed_pix_buf->buf;
    packed_pix_32 = (uint32*)packed_pix_buf->buf;
    packed_pix_64 = (uint64*)packed_pix_buf->buf;

    unpacked_pix = (uint16*)unpacked_pix_buf->buf;
    a_scale = (uint8*)a_scale_buf->buf;
    i_scale = (uint8*)i_scale_buf->buf;
    a_scale_16 = (uint16*)a_scale_buf->buf;
    i_scale_16 = (uint16*)i_scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 8) {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
               packed_pix_64[i] = PACK_AI_16(uint64)&0xFFffFFffFFffFFffULL;
            }
        } else {
            for (i=0; i < max_i; i++) {
               packed_pix_64[i] = PACK_AI(uint64)&0xFFffFFffFFffFFffULL;
            }
        }
    } else if (packed_pix_size == 4) {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
               packed_pix_32[i] = PACK_AI_16(uint32)&0xFFffFFffUL;
            }
        } else {
            for (i=0; i < max_i; i++) {
               packed_pix_32[i] = PACK_AI(uint32)&0xFFffFFffUL;
            }
        }
    } else if (packed_pix_size == 2) {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
               packed_pix_16[i] = PACK_AI_16(uint16)&0xFFff;
            }
        } else {
            for (i=0; i < max_i; i++) {
               packed_pix_16[i] = PACK_AI(uint16)&0xFFff;
            }
        }
    } else {
        if (a_scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
               packed_pix_8[i] = PACK_AI_16(uint8)&0xFF;
            }
        } else {
            for (i=0; i < max_i; i++) {
               packed_pix_8[i] = PACK_AI(uint8)&0xFF;
            }
        }
    }
}

static void pack_raw_1_channel_16bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *scale_buf, sint8 shift)
{
    Py_ssize_t packed_pix_size;
    uint16 *unpacked_pix;
    uint8  *scale;
    uint16 *scale_16;
    uint64 i=0, max_i=0;

    uint8  *packed_pix_8;
    uint16 *packed_pix_16;
    uint32 *packed_pix_32;
    uint64 *packed_pix_64;

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8 = (uint8*)packed_pix_buf->buf;
    packed_pix_16 = (uint16*)packed_pix_buf->buf;
    packed_pix_32 = (uint32*)packed_pix_buf->buf;
    packed_pix_64 = (uint64*)packed_pix_buf->buf;

    unpacked_pix = (uint16*)unpacked_pix_buf->buf;
    scale = (uint8*)scale_buf->buf;
    scale_16 = (uint16*)scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 8) {
        if (scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
                packed_pix_64[i] = ((uint64)scale_16[unpacked_pix[i]]<<shift);
            }
        } else {
            for (i=0; i < max_i; i++) {
                packed_pix_64[i] = ((uint64)scale[unpacked_pix[i]]<<shift);
            }
        }
    } else if (packed_pix_size == 4) {
        if (scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
                packed_pix_32[i] = scale_16[unpacked_pix[i]]<<shift;
            }
        } else {
            for (i=0; i < max_i; i++) {
                packed_pix_32[i] = scale[unpacked_pix[i]]<<shift;
            }
        }
    } else if (packed_pix_size == 2) {
        if (scale_buf->itemsize == 2) {
            for (i=0; i < max_i; i++) {
                packed_pix_16[i] = scale_16[unpacked_pix[i]]<<shift;
            }
        } else {
            for (i=0; i < max_i; i++) {
                packed_pix_16[i] = scale[unpacked_pix[i]]<<shift;
            }
        }
    } else {
        for (i=0; i < max_i; i++) {
            packed_pix_8[i] = scale[unpacked_pix[i]]<<shift;
        }
    }
}

static void pack_indexing(
    Py_buffer *packed_indexing_buf, Py_buffer *unpacked_indexing_buf,
    sint8 indexing_size)
{
    //THIS FUNCTION IS CURRENTLY UNTESTED.
    uint8 *packed_indexing;
    uint8 *unpacked_indexing;
    uint64 i=0, j=0, max_i=0;

    packed_indexing = (uint8*) packed_indexing_buf->buf;
    unpacked_indexing = (uint8*) unpacked_indexing_buf->buf;

    max_i = packed_indexing_buf->len;

    if (indexing_size == 4) {
        for (i=0; i < max_i; i++) {
            j = i<<1;
            packed_indexing[i] = ( unpacked_indexing[j] +
                                  (unpacked_indexing[j+1]<<4))&255;
        }
    } else if (indexing_size == 2) {
        for (i=0; i < max_i; i++) {
            j = i<<2;
            packed_indexing[i] = ( unpacked_indexing[j] +
                                  (unpacked_indexing[j+1]<<2)+
                                  (unpacked_indexing[j+2]<<4)+
                                  (unpacked_indexing[j+3]<<6))&255;
        }
    } else {
        for (i=0; i < max_i; i++) {
            j = i<<3;
            packed_indexing[i] = ( unpacked_indexing[j] +
                                  (unpacked_indexing[j+1]<<1)+
                                  (unpacked_indexing[j+2]<<2)+
                                  (unpacked_indexing[j+3]<<3)+
                                  (unpacked_indexing[j+4]<<4)+
                                  (unpacked_indexing[j+5]<<5)+
                                  (unpacked_indexing[j+6]<<6)+
                                  (unpacked_indexing[j+7]<<7))&255;
        }
    }
}


static PyObject *py_pack_raw_4_channel_merge(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];
    sint64 divs[4];
    sint8 shifts[4], i;

    // Get the pointers to each of the arrays, divisors, and shifts
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*LLLLbbbb:pack_raw_4_channel_merge",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        &divs[0], &divs[1], &divs[2], &divs[3],
        &shifts[0], &shifts[1], &shifts[2], &shifts[3]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[1].itemsize == 2) {
        pack_raw_4_channel_merge_16bpp(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
            divs[0], divs[1], divs[2], divs[3],
            shifts[0], shifts[1], shifts[2], shifts[3]);
    } else {
        pack_raw_4_channel_merge_8bpp(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
            divs[0], divs[1], divs[2], divs[3],
            shifts[0], shifts[1], shifts[2], shifts[3]);
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_pack_raw_2_channel_merge(PyObject *self, PyObject *args) {
    Py_buffer bufs[4];
    sint64 divs[2];
    sint8 shifts[2], i;

    // Get the pointers to each of the arrays, divisors, and shifts
    if (!PyArg_ParseTuple(args, "w*w*w*w*LLbb:pack_raw_2_channel_merge",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3],
        &divs[0], &divs[1], &shifts[0], &shifts[1]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[1].itemsize == 2) {
        pack_raw_2_channel_merge_16bpp(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3],
            divs[0], divs[1], shifts[0], shifts[1]);
    } else {
        pack_raw_2_channel_merge_8bpp(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3],
            divs[0], divs[1], shifts[0], shifts[1]);
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_pack_raw_4_channel(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];
    sint8 shifts[4], i;

    // Get the pointers to each of the arrays, masks, and shifts
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*bbbb:pack_raw_4_channel",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        &shifts[0], &shifts[1], &shifts[2], &shifts[3]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[1].itemsize == 2) {
        pack_raw_4_channel_16bpp(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
            shifts[0], shifts[1], shifts[2], shifts[3]);
    } else {
        pack_raw_4_channel_8bpp(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
            shifts[0], shifts[1], shifts[2], shifts[3]);
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_pack_raw_2_channel(PyObject *self, PyObject *args) {
    Py_buffer bufs[4];
    sint8 shifts[2], i;

    // Get the pointers to each of the arrays, masks, and shifts
    if (!PyArg_ParseTuple(args, "w*w*w*w*bb:pack_raw_2_channel",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &shifts[0], &shifts[1]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[1].itemsize == 2) {
        pack_raw_2_channel_16bpp(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], shifts[0], shifts[1]);
    } else {
        pack_raw_2_channel_8bpp(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], shifts[0], shifts[1]);
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_pack_raw_1_channel(PyObject *self, PyObject *args) {
    Py_buffer bufs[3];
    sint8 shift, i;

    // Get the pointers to each of the arrays, mask, and shift
    if (!PyArg_ParseTuple(args, "w*w*w*b:pack_raw_1_channel",
        &bufs[0], &bufs[1], &bufs[2], &shift))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[1].itemsize == 2) {
        pack_raw_1_channel_16bpp(&bufs[0], &bufs[1], &bufs[2], shift);
    } else {
        pack_raw_1_channel_8bpp(&bufs[0], &bufs[1], &bufs[2], shift);
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_pack_indexing(PyObject *self, PyObject *args) {
    Py_buffer bufs[2];
    sint8 indexing_size, i;

    // Get the pointers to each of the arrays and indexing size
    if (!PyArg_ParseTuple(args, "w*w*b:pack_indexing",
        &bufs[0], &bufs[1], &indexing_size))
        return Py_BuildValue("");  // return Py_None while incrementing it

    pack_indexing(&bufs[0], &bufs[1], indexing_size);

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

/* A list of all the methods defined by this module.
"METH_VARGS" tells Python how to call the handler.
The {NULL, NULL} entry indicates the end of the method definitions */
static PyMethodDef raw_packer_ext_methods[] = {
    {"pack_raw_4_channel", py_pack_raw_4_channel, METH_VARARGS, ""},
    {"pack_raw_2_channel", py_pack_raw_2_channel, METH_VARARGS, ""},
    {"pack_raw_1_channel", py_pack_raw_1_channel, METH_VARARGS, ""},
    {"pack_raw_4_channel_merge", py_pack_raw_4_channel_merge, METH_VARARGS, ""},
    {"pack_raw_2_channel_merge", py_pack_raw_2_channel_merge, METH_VARARGS, ""},
    {"pack_indexing", py_pack_indexing, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}      /* sentinel */
};

/* When Python imports a C module named 'X' it loads the
module then looks for a method named "init"+X and calls it.*/
static struct PyModuleDef raw_packer_ext_module = {
    PyModuleDef_HEAD_INIT,
    "raw_packer_ext",
    "A set of C functions to replace certain speed intensive pixel packer functions",
    -1,
    raw_packer_ext_methods,
};

PyMODINIT_FUNC PyInit_raw_packer_ext(void) {
    return PyModule_Create(&raw_packer_ext_module);
}
