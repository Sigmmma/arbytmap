#include <stdio.h>
#include "Python.h"
#include "abstract.h"
#include "longobject.h"
#include "modsupport.h"
#include "object.h"

static void pack_raw_4_channel_merge_8bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    long long a_div, long long r_div, long long g_div, long long b_div,
    char a_shift,  char r_shift,  char g_shift,  char b_shift)
{
    Py_ssize_t packed_pix_size;
    unsigned char (*unpacked_pix)[];
    unsigned char (*a_scale)[], (*r_scale)[], (*g_scale)[], (*b_scale)[];
    unsigned long long i=0, max_i=0;
    long long a_rnd, r_rnd, g_rnd, b_rnd;

    unsigned char  (*packed_pix_8)[];
    unsigned short (*packed_pix_16)[];
    unsigned long  (*packed_pix_32)[];
    unsigned long long (*packed_pix_64)[];

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8  = (unsigned char(*)[]) packed_pix_buf->buf;
    packed_pix_16 = (unsigned short(*)[])packed_pix_buf->buf;
    packed_pix_32 = (unsigned long(*)[]) packed_pix_buf->buf;
    packed_pix_64 = (unsigned long long(*)[])packed_pix_buf->buf;

    unpacked_pix = (unsigned char(*)[])unpacked_pix_buf->buf;
    a_scale = (unsigned char(*)[])a_scale_buf->buf;
    r_scale = (unsigned char(*)[])r_scale_buf->buf;
    g_scale = (unsigned char(*)[])g_scale_buf->buf;
    b_scale = (unsigned char(*)[])b_scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;
    a_rnd=a_div/2; r_rnd=r_div/2; g_rnd=g_div/2; b_rnd=b_div/2;

    if (packed_pix_size == 8) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_64)[i] = (
                ((*a_scale)[((*unpacked_pix)[(i<<2)]+a_rnd)/a_div]<<a_shift) +
                ((*r_scale)[((*unpacked_pix)[(i<<2)+1]+r_rnd)/r_div]<<r_shift) +
                ((*g_scale)[((*unpacked_pix)[(i<<2)+2]+g_rnd)/g_div]<<g_shift) +
                ((*b_scale)[((*unpacked_pix)[(i<<2)+3]+b_rnd)/b_div]<<b_shift))&0xFFffFFffFFffFFffLL;
        }
    } else if (packed_pix_size == 4) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_32)[i] = (
                ((*a_scale)[((*unpacked_pix)[(i<<2)]+a_rnd)/a_div]<<a_shift) +
                ((*r_scale)[((*unpacked_pix)[(i<<2)+1]+r_rnd)/r_div]<<r_shift) +
                ((*g_scale)[((*unpacked_pix)[(i<<2)+2]+g_rnd)/g_div]<<g_shift) +
                ((*b_scale)[((*unpacked_pix)[(i<<2)+3]+b_rnd)/b_div]<<b_shift))&0xFFffFFffL;
        }
    } else if (packed_pix_size == 2) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_16)[i] = (
                ((*a_scale)[((*unpacked_pix)[(i<<2)]+a_rnd)/a_div]<<a_shift) +
                ((*r_scale)[((*unpacked_pix)[(i<<2)+1]+r_rnd)/r_div]<<r_shift) +
                ((*g_scale)[((*unpacked_pix)[(i<<2)+2]+g_rnd)/g_div]<<g_shift) +
                ((*b_scale)[((*unpacked_pix)[(i<<2)+3]+b_rnd)/b_div]<<b_shift))&0xFFff;
        }
    } else {
        for (i=0; i < max_i; i++) {
            (*packed_pix_8)[i] = (
                ((*a_scale)[((*unpacked_pix)[(i<<2)]+a_rnd)/a_div]<<a_shift) +
                ((*r_scale)[((*unpacked_pix)[(i<<2)+1]+r_rnd)/r_div]<<r_shift) +
                ((*g_scale)[((*unpacked_pix)[(i<<2)+2]+g_rnd)/g_div]<<g_shift) +
                ((*b_scale)[((*unpacked_pix)[(i<<2)+3]+b_rnd)/b_div]<<b_shift))&0xFF;
        }
    }
}

static void pack_raw_2_channel_merge_8bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *i_scale_buf,
    long long a_div, long long i_div, char a_shift,  char i_shift)
{
    Py_ssize_t packed_pix_size;
    unsigned char (*unpacked_pix)[];
    unsigned char (*a_scale)[], (*i_scale)[];
    unsigned long long i=0, max_i=0;
    long long a_rnd=a_div/2, i_rnd=i_div/2;

    unsigned char  (*packed_pix_8)[];
    unsigned short (*packed_pix_16)[];
    unsigned long  (*packed_pix_32)[];

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8  = (unsigned char(*)[]) packed_pix_buf->buf;
    packed_pix_16 = (unsigned short(*)[])packed_pix_buf->buf;
    packed_pix_32 = (unsigned long(*)[]) packed_pix_buf->buf;

    unpacked_pix = (unsigned char(*)[])unpacked_pix_buf->buf;
    a_scale = (unsigned char(*)[])a_scale_buf->buf;
    i_scale = (unsigned char(*)[])i_scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 4) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_32)[i] = (
                ((*a_scale)[((*unpacked_pix)[(i<<1)]+a_rnd)/a_div]<<a_shift) +
                ((*i_scale)[((*unpacked_pix)[(i<<1)+1]+i_rnd)/i_div]<<i_shift))&0xFFffFFffL;
        }
    } else if (packed_pix_size == 2) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_16)[i] = (
                ((*a_scale)[((*unpacked_pix)[(i<<1)]+a_rnd)/a_div]<<a_shift) +
                ((*i_scale)[((*unpacked_pix)[(i<<1)+1]+i_rnd)/i_div]<<i_shift))&0xFFff;
        }
    } else {
        for (i=0; i < max_i; i++) {
            (*packed_pix_8)[i] = (
                ((*a_scale)[((*unpacked_pix)[(i<<1)]+a_rnd)/a_div]<<a_shift) +
                ((*i_scale)[((*unpacked_pix)[(i<<1)+1]+i_rnd)/i_div]<<i_shift))&0xFF;
        }
    }
}

static void pack_raw_4_channel_8bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    char a_shift,  char r_shift,  char g_shift,  char b_shift)
{
    Py_ssize_t packed_pix_size;
    unsigned char (*unpacked_pix)[];
    unsigned char (*a_scale)[], (*r_scale)[], (*g_scale)[], (*b_scale)[];
    unsigned long long i=0, max_i=0;

    unsigned char  (*packed_pix_8)[];
    unsigned short (*packed_pix_16)[];
    unsigned long  (*packed_pix_32)[];
    unsigned long long (*packed_pix_64)[];

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8  = (unsigned char(*)[]) packed_pix_buf->buf;
    packed_pix_16 = (unsigned short(*)[])packed_pix_buf->buf;
    packed_pix_32 = (unsigned long(*)[]) packed_pix_buf->buf;
    packed_pix_64 = (unsigned long long(*)[])packed_pix_buf->buf;

    unpacked_pix = (unsigned char(*)[])unpacked_pix_buf->buf;
    a_scale = (unsigned char(*)[])a_scale_buf->buf;
    r_scale = (unsigned char(*)[])r_scale_buf->buf;
    g_scale = (unsigned char(*)[])g_scale_buf->buf;
    b_scale = (unsigned char(*)[])b_scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 8) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_64)[i] = (
                ((*a_scale)[(*unpacked_pix)[(i<<2)]]<<a_shift) +
                ((*r_scale)[(*unpacked_pix)[(i<<2)+1]]<<r_shift) +
                ((*g_scale)[(*unpacked_pix)[(i<<2)+2]]<<g_shift) +
                ((*b_scale)[(*unpacked_pix)[(i<<2)+3]]<<b_shift))&0xFFffFFffFFffFFffLL;
        }
    } else if (packed_pix_size == 4) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_32)[i] = (
                ((*a_scale)[(*unpacked_pix)[(i<<2)]]<<a_shift) +
                ((*r_scale)[(*unpacked_pix)[(i<<2)+1]]<<r_shift) +
                ((*g_scale)[(*unpacked_pix)[(i<<2)+2]]<<g_shift) +
                ((*b_scale)[(*unpacked_pix)[(i<<2)+3]]<<b_shift))&0xFFffFFffL;
        }
    } else if (packed_pix_size == 2) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_16)[i] = (
                ((*a_scale)[(*unpacked_pix)[(i<<2)]]<<a_shift) +
                ((*r_scale)[(*unpacked_pix)[(i<<2)+1]]<<r_shift) +
                ((*g_scale)[(*unpacked_pix)[(i<<2)+2]]<<g_shift) +
                ((*b_scale)[(*unpacked_pix)[(i<<2)+3]]<<b_shift))&0xFFff;
        }
    } else {
        for (i=0; i < max_i; i++) {
            (*packed_pix_8)[i] = (
                ((*a_scale)[(*unpacked_pix)[(i<<2)]]<<a_shift) +
                ((*r_scale)[(*unpacked_pix)[(i<<2)+1]]<<r_shift) +
                ((*g_scale)[(*unpacked_pix)[(i<<2)+2]]<<g_shift) +
                ((*b_scale)[(*unpacked_pix)[(i<<2)+3]]<<b_shift))&0xFF;
        }
    }
}

static void pack_raw_2_channel_8bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *i_scale_buf,
    char a_shift, char i_shift)
{
    Py_ssize_t packed_pix_size;
    unsigned char (*unpacked_pix)[];
    unsigned char (*a_scale)[], (*i_scale)[];
    unsigned long long i=0, max_i=0;

    unsigned char  (*packed_pix_8)[];
    unsigned short (*packed_pix_16)[];
    unsigned long  (*packed_pix_32)[];

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8  = (unsigned char(*)[]) packed_pix_buf->buf;
    packed_pix_16 = (unsigned short(*)[])packed_pix_buf->buf;
    packed_pix_32 = (unsigned long(*)[]) packed_pix_buf->buf;

    unpacked_pix = (unsigned char(*)[])unpacked_pix_buf->buf;
    a_scale = (unsigned char(*)[])a_scale_buf->buf;
    i_scale = (unsigned char(*)[])i_scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 4) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_32)[i] = (
                ((*a_scale)[(*unpacked_pix)[(i<<1)]]<<a_shift) +
                ((*i_scale)[(*unpacked_pix)[(i<<1)+1]]<<i_shift))&0xFFffFFffL;
        }
    } else if (packed_pix_size == 2) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_16)[i] = (
                ((*a_scale)[(*unpacked_pix)[(i<<1)]]<<a_shift) +
                ((*i_scale)[(*unpacked_pix)[(i<<1)+1]]<<i_shift))&0xFFff;
        }
    } else {
        for (i=0; i < max_i; i++) {
            (*packed_pix_8)[i] = (
                ((*a_scale)[(*unpacked_pix)[(i<<1)]]<<a_shift) +
                ((*i_scale)[(*unpacked_pix)[(i<<1)+1]]<<i_shift))&0xFF;
        }
    }
}

static void pack_raw_1_channel_8bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf, Py_buffer *scale_buf)
{
    Py_ssize_t packed_pix_size;
    unsigned char (*unpacked_pix)[];
    unsigned char (*scale)[];
    unsigned long long i=0, max_i=0;

    unsigned char  (*packed_pix_8)[];
    unsigned short (*packed_pix_16)[];

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8  = (unsigned char(*)[]) packed_pix_buf->buf;
    packed_pix_16 = (unsigned short(*)[])packed_pix_buf->buf;

    unpacked_pix = (unsigned char(*)[])unpacked_pix_buf->buf;
    scale = (unsigned char(*)[])scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 2) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_16)[i] = (*scale)[(*unpacked_pix)[i]];
        }
    } else {
        for (i=0; i < max_i; i++) {
            (*packed_pix_8)[i] = (*scale)[(*unpacked_pix)[i]];
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
    long long a_div, long long r_div, long long g_div, long long b_div,
    char a_shift,  char r_shift,  char g_shift,  char b_shift)
{
    Py_ssize_t packed_pix_size;
    unsigned short (*unpacked_pix)[];
    unsigned short (*a_scale)[], (*r_scale)[], (*g_scale)[], (*b_scale)[];
    unsigned long long i=0, max_i=0;
    long long a_rnd=a_div/2, r_rnd=r_div/2, g_rnd=g_div/2, b_rnd=b_div/2;

    unsigned char  (*packed_pix_8)[];
    unsigned short (*packed_pix_16)[];
    unsigned long  (*packed_pix_32)[];
    unsigned long long (*packed_pix_64)[];

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8  = (unsigned char(*)[]) packed_pix_buf->buf;
    packed_pix_16 = (unsigned short(*)[])packed_pix_buf->buf;
    packed_pix_32 = (unsigned long(*)[]) packed_pix_buf->buf;
    packed_pix_64 = (unsigned long long(*)[])packed_pix_buf->buf;

    unpacked_pix = (unsigned short(*)[])unpacked_pix_buf->buf;
    a_scale = (unsigned short(*)[])a_scale_buf->buf;
    r_scale = (unsigned short(*)[])r_scale_buf->buf;
    g_scale = (unsigned short(*)[])g_scale_buf->buf;
    b_scale = (unsigned short(*)[])b_scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 8) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_64)[i] = (
                ((*a_scale)[((*unpacked_pix)[(i<<2)]+a_rnd)/a_div]<<a_shift) +
                ((*r_scale)[((*unpacked_pix)[(i<<2)+1]+r_rnd)/r_div]<<r_shift) +
                ((*g_scale)[((*unpacked_pix)[(i<<2)+2]+g_rnd)/g_div]<<g_shift) +
                ((*b_scale)[((*unpacked_pix)[(i<<2)+3]+b_rnd)/b_div]<<b_shift))&0xFFffFFffFFffFFffLL;
        }
    } else if (packed_pix_size == 4) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_32)[i] = (
                ((*a_scale)[((*unpacked_pix)[(i<<2)]+a_rnd)/a_div]<<a_shift) +
                ((*r_scale)[((*unpacked_pix)[(i<<2)+1]+r_rnd)/r_div]<<r_shift) +
                ((*g_scale)[((*unpacked_pix)[(i<<2)+2]+g_rnd)/g_div]<<g_shift) +
                ((*b_scale)[((*unpacked_pix)[(i<<2)+3]+b_rnd)/b_div]<<b_shift))&0xFFffFFffL;
        }
    } else if (packed_pix_size == 2) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_16)[i] = (
                ((*a_scale)[((*unpacked_pix)[(i<<2)]+a_rnd)/a_div]<<a_shift) +
                ((*r_scale)[((*unpacked_pix)[(i<<2)+1]+r_rnd)/r_div]<<r_shift) +
                ((*g_scale)[((*unpacked_pix)[(i<<2)+2]+g_rnd)/g_div]<<g_shift) +
                ((*b_scale)[((*unpacked_pix)[(i<<2)+3]+b_rnd)/b_div]<<b_shift))&0xFFff;
        }
    } else {
        for (i=0; i < max_i; i++) {
            (*packed_pix_8)[i] = (
                ((*a_scale)[((*unpacked_pix)[(i<<2)]+a_rnd)/a_div]<<a_shift) +
                ((*r_scale)[((*unpacked_pix)[(i<<2)+1]+r_rnd)/r_div]<<r_shift) +
                ((*g_scale)[((*unpacked_pix)[(i<<2)+2]+g_rnd)/g_div]<<g_shift) +
                ((*b_scale)[((*unpacked_pix)[(i<<2)+3]+b_rnd)/b_div]<<b_shift))&0xFF;
        }
    }
}

static void pack_raw_2_channel_merge_16bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *i_scale_buf,
    long long a_div, long long i_div, char a_shift,  char i_shift)
{
    Py_ssize_t packed_pix_size;
    unsigned short (*unpacked_pix)[];
    unsigned short (*a_scale)[], (*i_scale)[];
    unsigned long long i=0, max_i=0;
    long long a_rnd=a_div/2, i_rnd=i_div/2;

    unsigned char  (*packed_pix_8)[];
    unsigned short (*packed_pix_16)[];
    unsigned long  (*packed_pix_32)[];

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8  = (unsigned char(*)[]) packed_pix_buf->buf;
    packed_pix_16 = (unsigned short(*)[])packed_pix_buf->buf;
    packed_pix_32 = (unsigned long(*)[]) packed_pix_buf->buf;

    unpacked_pix = (unsigned short(*)[])unpacked_pix_buf->buf;
    a_scale = (unsigned short(*)[])a_scale_buf->buf;
    i_scale = (unsigned short(*)[])i_scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 4) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_32)[i] = (
                ((*a_scale)[((*unpacked_pix)[(i<<1)]+a_rnd)/a_div]<<a_shift) +
                ((*i_scale)[((*unpacked_pix)[(i<<1)+1]+i_rnd)/i_div]<<i_shift))&0xFFffFFffL;
        }
    } else if (packed_pix_size == 2) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_16)[i] = (
                ((*a_scale)[((*unpacked_pix)[(i<<1)]+a_rnd)/a_div]<<a_shift) +
                ((*i_scale)[((*unpacked_pix)[(i<<1)+1]+i_rnd)/i_div]<<i_shift))&0xFFff;
        }
    } else {
        for (i=0; i < max_i; i++) {
            (*packed_pix_8)[i] = (
                ((*a_scale)[((*unpacked_pix)[(i<<1)]+a_rnd)/a_div]<<a_shift) +
                ((*i_scale)[((*unpacked_pix)[(i<<1)+1]+i_rnd)/i_div]<<i_shift))&0xFF;
        }
    }
}


static void pack_raw_4_channel_16bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    char a_shift,  char r_shift,  char g_shift,  char b_shift)
{
    Py_ssize_t packed_pix_size;
    unsigned short (*unpacked_pix)[];
    unsigned short (*a_scale)[], (*r_scale)[], (*g_scale)[], (*b_scale)[];
    unsigned long long i=0, max_i=0;

    unsigned char  (*packed_pix_8)[];
    unsigned short (*packed_pix_16)[];
    unsigned long  (*packed_pix_32)[];
    unsigned long long (*packed_pix_64)[];

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8  = (unsigned char(*)[]) packed_pix_buf->buf;
    packed_pix_16 = (unsigned short(*)[])packed_pix_buf->buf;
    packed_pix_32 = (unsigned long(*)[]) packed_pix_buf->buf;
    packed_pix_64 = (unsigned long long(*)[])packed_pix_buf->buf;

    unpacked_pix = (unsigned short(*)[])unpacked_pix_buf->buf;
    a_scale = (unsigned short(*)[])a_scale_buf->buf;
    r_scale = (unsigned short(*)[])r_scale_buf->buf;
    g_scale = (unsigned short(*)[])g_scale_buf->buf;
    b_scale = (unsigned short(*)[])b_scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 8) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_64)[i] = (
                ((*a_scale)[(*unpacked_pix)[(i<<2)]]<<a_shift) +
                ((*r_scale)[(*unpacked_pix)[(i<<2)+1]]<<r_shift) +
                ((*g_scale)[(*unpacked_pix)[(i<<2)+2]]<<g_shift) +
                ((*b_scale)[(*unpacked_pix)[(i<<2)+3]]<<b_shift))&0xFFffFFffFFffFFffLL;
        }
    } else if (packed_pix_size == 4) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_32)[i] = (
                ((*a_scale)[(*unpacked_pix)[(i<<2)]]<<a_shift) +
                ((*r_scale)[(*unpacked_pix)[(i<<2)+1]]<<r_shift) +
                ((*g_scale)[(*unpacked_pix)[(i<<2)+2]]<<g_shift) +
                ((*b_scale)[(*unpacked_pix)[(i<<2)+3]]<<b_shift))&0xFFffFFffL;
        }
    } else if (packed_pix_size == 2) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_16)[i] = (
                ((*a_scale)[(*unpacked_pix)[(i<<2)]]<<a_shift) +
                ((*r_scale)[(*unpacked_pix)[(i<<2)+1]]<<r_shift) +
                ((*g_scale)[(*unpacked_pix)[(i<<2)+2]]<<g_shift) +
                ((*b_scale)[(*unpacked_pix)[(i<<2)+3]]<<b_shift))&0xFFff;
        }
    } else {
        for (i=0; i < max_i; i++) {
            (*packed_pix_8)[i] = (
                ((*a_scale)[(*unpacked_pix)[(i<<2)]]<<a_shift) +
                ((*r_scale)[(*unpacked_pix)[(i<<2)+1]]<<r_shift) +
                ((*g_scale)[(*unpacked_pix)[(i<<2)+2]]<<g_shift) +
                ((*b_scale)[(*unpacked_pix)[(i<<2)+3]]<<b_shift))&0xFF;
        }
    }
}

static void pack_raw_2_channel_16bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *i_scale_buf,
    char a_shift, char i_shift)
{
    Py_ssize_t packed_pix_size;
    unsigned short (*unpacked_pix)[];
    unsigned short (*a_scale)[], (*i_scale)[];
    unsigned long long i=0, max_i=0;

    unsigned char  (*packed_pix_8)[];
    unsigned short (*packed_pix_16)[];
    unsigned long  (*packed_pix_32)[];

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_8  = (unsigned char(*)[]) packed_pix_buf->buf;
    packed_pix_16 = (unsigned short(*)[])packed_pix_buf->buf;
    packed_pix_32 = (unsigned long(*)[]) packed_pix_buf->buf;

    unpacked_pix = (unsigned short(*)[])unpacked_pix_buf->buf;
    a_scale = (unsigned short(*)[])a_scale_buf->buf;
    i_scale = (unsigned short(*)[])i_scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    if (packed_pix_size == 4) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_32)[i] = (
                ((*a_scale)[(*unpacked_pix)[(i<<1)]]<<a_shift) +
                ((*i_scale)[(*unpacked_pix)[(i<<1)+1]]<<i_shift))&0xFFffFFffL;
        }
    } else if (packed_pix_size == 2) {
        for (i=0; i < max_i; i++) {
            (*packed_pix_16)[i] = (
                ((*a_scale)[(*unpacked_pix)[(i<<1)]]<<a_shift) +
                ((*i_scale)[(*unpacked_pix)[(i<<1)+1]]<<i_shift))&0xFFff;
        }
    } else {
        for (i=0; i < max_i; i++) {
            (*packed_pix_8)[i] = (
                ((*a_scale)[(*unpacked_pix)[(i<<1)]]<<a_shift) +
                ((*i_scale)[(*unpacked_pix)[(i<<1)+1]]<<i_shift))&0xFF;
        }
    }
}

static void pack_raw_1_channel_16bpp(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf, Py_buffer *scale_buf)
{
    Py_ssize_t packed_pix_size;
    unsigned short (*unpacked_pix)[];
    unsigned short (*scale)[];
    unsigned long long i=0, max_i=0;
    unsigned short (*packed_pix_16)[];

    packed_pix_size = packed_pix_buf->itemsize;

    packed_pix_16 = (unsigned short(*)[])packed_pix_buf->buf;
    unpacked_pix = (unsigned short(*)[])unpacked_pix_buf->buf;
    scale = (unsigned short(*)[])scale_buf->buf;

    max_i = (packed_pix_buf->len)/packed_pix_size;

    for (i=0; i < max_i; i++) {
        (*packed_pix_16)[i] = (*scale)[(*unpacked_pix)[i]];
    }
}

static void pack_indexing(
    Py_buffer *packed_indexing_buf, Py_buffer *unpacked_indexing_buf,
    char indexing_size)
{
    //THIS FUNCTION IS CURRENTLY UNTESTED.
    unsigned char (*packed_indexing)[];
    unsigned char (*unpacked_indexing)[];
    unsigned long long i=0, max_i=0;

    packed_indexing = (unsigned char(*)[]) packed_indexing_buf->buf;
    unpacked_indexing = (unsigned char(*)[]) unpacked_indexing_buf->buf;

    max_i = packed_indexing_buf->len;

    if (indexing_size == 4) {
        for (i=0; i < max_i; i++) {
            (*packed_indexing)[i] = ( (*unpacked_indexing)[(i<<1)] +
                                     ((*unpacked_indexing)[(i<<1)+1]<<4))&255;
        }
    } else if (indexing_size == 2) {
        for (i=0; i < max_i; i++) {
            (*packed_indexing)[i] = ( (*unpacked_indexing)[(i<<2)] +
                                     ((*unpacked_indexing)[(i<<2)+1]<<2)+
                                     ((*unpacked_indexing)[(i<<2)+2]<<4)+
                                     ((*unpacked_indexing)[(i<<2)+3]<<6))&255;
        }
    } else {
        for (i=0; i < max_i; i++) {
            (*packed_indexing)[i] = ( (*unpacked_indexing)[(i<<3)] +
                                     ((*unpacked_indexing)[(i<<3)+1]<<1)+
                                     ((*unpacked_indexing)[(i<<3)+2]<<2)+
                                     ((*unpacked_indexing)[(i<<3)+3]<<3)+
                                     ((*unpacked_indexing)[(i<<3)+4]<<4)+
                                     ((*unpacked_indexing)[(i<<3)+5]<<5)+
                                     ((*unpacked_indexing)[(i<<3)+6]<<6)+
                                     ((*unpacked_indexing)[(i<<3)+7]<<7))&255;
        }
    }
}


static PyObject *py_pack_raw_4_channel_merge(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];
    long long divs[4];
    char shifts[4];

    // Get the pointers to each of the arrays, divisors, and shifts
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*LLLLbbbb:pack_raw_4_channel_merge",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        &divs[0], &divs[1], &divs[2], &divs[3],
        &shifts[0], &shifts[1], &shifts[2], &shifts[3]))
        return Py_None;

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
    PyBuffer_Release(&bufs[0]); PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]); PyBuffer_Release(&bufs[3]);
    PyBuffer_Release(&bufs[4]); PyBuffer_Release(&bufs[5]);

    return Py_None;
}

static PyObject *py_pack_raw_2_channel_merge(PyObject *self, PyObject *args) {
    Py_buffer bufs[4];
    long long divs[2];
    char shifts[2];

    // Get the pointers to each of the arrays, divisors, and shifts
    if (!PyArg_ParseTuple(args, "w*w*w*w*LLbb:pack_raw_2_channel_merge",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3],
        &divs[0], &divs[1], &shifts[0], &shifts[1]))
        return Py_None;

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
    PyBuffer_Release(&bufs[0]); PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]); PyBuffer_Release(&bufs[3]);

    return Py_None;
}

static PyObject *py_pack_raw_4_channel(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];
    char shifts[4];

    // Get the pointers to each of the arrays, masks, and shifts
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*bbbb:pack_raw_4_channel",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        &shifts[0], &shifts[1], &shifts[2], &shifts[3]))
        return Py_None;

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
    PyBuffer_Release(&bufs[0]); PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]); PyBuffer_Release(&bufs[3]);
    PyBuffer_Release(&bufs[4]); PyBuffer_Release(&bufs[5]);

    return Py_None;
}

static PyObject *py_pack_raw_2_channel(PyObject *self, PyObject *args) {
    Py_buffer bufs[4];
    char shifts[2];

    // Get the pointers to each of the arrays, masks, and shifts
    if (!PyArg_ParseTuple(args, "w*w*w*w*bb:pack_raw_2_channel",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &shifts[0], &shifts[1]))
        return Py_None;

    if (bufs[1].itemsize == 2) {
        pack_raw_2_channel_16bpp(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], shifts[0], shifts[1]);
    } else {
        pack_raw_2_channel_8bpp(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], shifts[0], shifts[1]);
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]); PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]); PyBuffer_Release(&bufs[3]);

    return Py_None;
}

static PyObject *py_pack_raw_1_channel(PyObject *self, PyObject *args) {
    Py_buffer bufs[3];
    char shift;

    // Get the pointers to each of the arrays, mask, and shift
    if (!PyArg_ParseTuple(args, "w*w*w*:pack_raw_1_channel",
        &bufs[0], &bufs[1], &bufs[2], &shift))
        return Py_None;

    if (bufs[1].itemsize == 2) {
        pack_raw_1_channel_16bpp(&bufs[0], &bufs[1], &bufs[2]);
    } else {
        pack_raw_1_channel_8bpp(&bufs[0], &bufs[1], &bufs[2]);
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]);

    return Py_None;
}

static PyObject *py_pack_indexing(PyObject *self, PyObject *args) {
    Py_buffer bufs[2];
    char indexing_size;

    // Get the pointers to each of the arrays and indexing size
    if (!PyArg_ParseTuple(args, "w*w*b:pack_indexing",
        &bufs[0], &bufs[1], &indexing_size))
        return Py_None;

    pack_indexing(&bufs[0], &bufs[1], indexing_size);

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);

    return Py_None;
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
