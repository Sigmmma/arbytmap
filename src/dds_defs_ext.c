#include <math.h>
#include "shared.h"

// to prevent loss of precision, cast these to ULL's before squaring
#define SQ(x) ((uint64)x)*((uint64)x)

#define READ_DXT_COLORS(x)\
    color0 = packed_tex[x]&0xFFff;\
    color1 = packed_tex[x]>>16;\
    color_idx = packed_tex[x+1];\
    \
    c_0[1] = r_scale[(color0>>11) & 31];\
    c_0[2] = g_scale[(color0>>5)  & 63];\
    c_0[3] = b_scale[ color0      & 31];\
    \
    c_1[1] = r_scale[(color1>>11) & 31];\
    c_1[2] = g_scale[(color1>>5)  & 63];\
    c_1[3] = b_scale[ color1      & 31];\
    \
    if (color0 > color1) {\
        c_2[1] = (c_0[1]*2 + c_1[1])/3;\
        c_2[2] = (c_0[2]*2 + c_1[2])/3;\
        c_2[3] = (c_0[3]*2 + c_1[3])/3;\
        \
        colors[3] = c_3;\
        c_3[1] = (c_0[1] + 2*c_1[1])/3;\
        c_3[2] = (c_0[2] + 2*c_1[2])/3;\
        c_3[3] = (c_0[3] + 2*c_1[3])/3;\
    } else {\
        c_2[1] = (c_0[1]+c_1[1])/2;\
        c_2[2] = (c_0[2]+c_1[2])/2;\
        c_2[3] = (c_0[3]+c_1[3])/2;\
        colors[3] = transparent;\
    }

#define UNPACK_DXT_COLORS()\
    color = colors[(color_idx>>(j<<1))&3];\
    off = j*ucc + pxl_i;\
    if (chan1 >= 0) unpacked_pix[off + 1] = color[chan1];\
    if (chan2 >= 0) unpacked_pix[off + 2] = color[chan2];\
    if (chan3 >= 0) unpacked_pix[off + 3] = color[chan3];

#define UNPACK_DXT2345_COLORS()\
    color = colors[(color_idx>>(j<<1))&3];\
    off = j*ucc + pxl_i;\
    if      (chan0 == 0) unpacked_pix[off] = a_scale[a];\
    else if (chan0 >  0) unpacked_pix[off] = color[chan0];\
    else                 unpacked_pix[off] = unpack_max;\
    \
    if      (chan1 == 0) unpacked_pix[off + 1] = a_scale[a];\
    else if (chan1 >  0) unpacked_pix[off + 1] = color[chan1];\
    else                 unpacked_pix[off + 1] = unpack_max;\
    \
    if      (chan2 == 0) unpacked_pix[off + 2] = a_scale[a];\
    else if (chan2 >  0) unpacked_pix[off + 2] = color[chan2];\
    else                 unpacked_pix[off + 2] = unpack_max;\
    \
    if      (chan3 == 0) unpacked_pix[off + 3] = a_scale[a];\
    else if (chan3 >  0) unpacked_pix[off + 3] = color[chan3];\
    else                 unpacked_pix[off + 3] = unpack_max;

#define READ_DXT5_A(lookup, idx, j, val0, val1)\
    lookup[0] = val0 = packed_tex[j]&0xFF;\
    lookup[1] = val1 = (packed_tex[j]>>8)&0xFF;\
    idx = (((uint64)packed_tex[(j)+1])<<16) + (packed_tex[(j)]>>16);\
    if (val0 > val1) {\
        lookup[2] = (val0*6 + val1)/7;\
        lookup[3] = (val0*5 + val1*2)/7;\
        lookup[4] = (val0*4 + val1*3)/7;\
        lookup[5] = (val0*3 + val1*4)/7;\
        lookup[6] = (val0*2 + val1*5)/7;\
        lookup[7] = (val0   + val1*6)/7;\
    } else {\
        lookup[2] = (val0*4 + val1)/5;\
        lookup[3] = (val0*3 + val1*2)/5;\
        lookup[4] = (val0*2 + val1*3)/5;\
        lookup[5] = (val0   + val1*4)/5;\
        lookup[6] = 0;\
        lookup[7] = 255;\
    }

#define CALC_B_NORMALIZE(r, g, b, x, y, mask, sign_mask, typ)\
    if (r&sign_mask) {x = r&mask;} else {x = mask - r;}\
    if (g&sign_mask) {y = g&mask;} else {y = mask - g;}\
    d = (double)(mask*mask - x*x - y*y);\
    if (d > 0) {\
        b = (typ)(sqrt(d)) + sign_mask;\
    } else {\
        n_len = sqrt(mask*mask - d)/mask;\
        x = (typ)(x/n_len);\
        y = (typ)(y/n_len);\
        \
        if (r&sign_mask) {r = x + sign_mask;} else {r = mask - x;}\
        if (g&sign_mask) {g = y + sign_mask;} else {g = mask - y;}\
        b = sign_mask;\
    }

#define CALC_W_NORMALIZE(u, v, w, max_len, mask, sign_mask, typ)\
    /*This function takes a ones signed u and v that were cast as unsigned\
      integers so it can preform its own ones signed operations.*/\
    if (u&sign_mask) {u -= mask;}\
    if (v&sign_mask) {v -= mask;}\
    d = (double)(max_len*max_len - u*u - v*v);\
    if (d > 0) {\
        w = (typ)(sqrt(d));\
    } else {\
        n_len = sqrt(max_len*max_len - d)/max_len;\
        u = (typ)(u/n_len);\
        v = (typ)(v/n_len);\
        w = 0;\
    }

#define PICK_DXT_PALETTE_DIST()\
    for (j = 0; j < chans_per_tex; j+=ucc) {\
        r = unpacked[r_pxl_i + j];\
        g = unpacked[g_pxl_i + j];\
        b = unpacked[b_pxl_i + j];\
        for (k = j + ucc; k < chans_per_tex; k+=ucc) {\
            /* if this will be solid black we wont care about its color*/\
            if (make_alpha && (unpacked[pxl_i + k] < a_cutoff)) continue;\
            dist1 = (SQ(r - unpacked[r_pxl_i + k])+\
                     SQ(g - unpacked[g_pxl_i + k])+\
                     SQ(b - unpacked[b_pxl_i + k]));\
            if (dist1 >= dist0) {\
                dist0 = dist1;\
                c_0i = j;\
                c_1i = k;\
            }\
        }\
    }

#define PACK_DXT_COLOR(c, c_i, packed_c)\
    c[1] = unpacked[r_pxl_i + c_i];\
    c[2] = unpacked[g_pxl_i + c_i];\
    c[3] = unpacked[b_pxl_i + c_i];\
    packed_c = ((r_scale[c[1]]<<11) + (g_scale[c[2]]<<5) + b_scale[c[3]])&0xFFff;

#define SWAP_COLORS()\
    tmp[1] = c_0[1]; tmp[2] = c_0[2]; tmp[3] = c_0[3];\
    c_0[1] = c_1[1]; c_0[2] = c_1[2]; c_0[3] = c_1[3];\
    c_1[1] = tmp[1]; c_1[2] = tmp[2]; c_1[3] = tmp[3];\
    tmp_color = color0;\
    color0    = color1;\
    color1    = tmp_color;

#define CALC_DXT_IDX_NO_ALPHA(idx)\
    c_2[1] = (c_0[1]*2 + c_1[1])/3; c_3[1] = (c_0[1] + c_1[1]*2)/3;\
    c_2[2] = (c_0[2]*2 + c_1[2])/3; c_3[2] = (c_0[2] + c_1[2]*2)/3;\
    c_2[3] = (c_0[3]*2 + c_1[3])/3; c_3[3] = (c_0[3] + c_1[3]*2)/3;\
    for (j=0; j<chans_per_tex; j+=ucc) {\
        r = unpacked[r_pxl_i+j];\
        g = unpacked[g_pxl_i+j];\
        b = unpacked[b_pxl_i+j];\
        dist0 = SQ(r-c_0[1]) + SQ(g-c_0[2]) + SQ(b-c_0[3]);\
        dist2 = SQ(r-c_2[1]) + SQ(g-c_2[2]) + SQ(b-c_2[3]);\
        dist3 = SQ(r-c_3[1]) + SQ(g-c_3[2]) + SQ(b-c_3[3]);\
        dist1 = SQ(r-c_1[1]) + SQ(g-c_1[2]) + SQ(b-c_1[3]);\
        if (dist1 < dist0) {\
            if (dist1 < dist3) {\
                idx += 1<<(j>>1);\
            } else {\
                idx += 3<<(j>>1);\
            }\
        } else if (dist2 < dist0) {\
            idx += 2<<(j>>1);\
        }\
    }

#define CALC_DXT_IDX_ALPHA(idx)\
    c_2[1] = (c_0[1] + c_1[1])/2;\
    c_2[2] = (c_0[2] + c_1[2])/2;\
    c_2[3] = (c_0[3] + c_1[3])/2;\
    for (j=0; j<chans_per_tex; j+=ucc) {\
        if (unpacked[pxl_i + j] < a_cutoff) {\
            idx += 3<<(j>>1);\
            continue;\
        }\
        r = unpacked[r_pxl_i+j];\
        g = unpacked[g_pxl_i+j];\
        b = unpacked[b_pxl_i+j];\
        dist0 = SQ(r-c_0[1]) + SQ(g-c_0[2]) + SQ(b-c_0[3]);\
        dist2 = SQ(r-c_2[1]) + SQ(g-c_2[2]) + SQ(b-c_2[3]);\
        dist1 = SQ(r-c_1[1]) + SQ(g-c_1[2]) + SQ(b-c_1[3]);\
        if ((dist1 < dist0) && (dist1 < dist2)) {\
            idx += 1<<(j>>1);\
        } else if (dist2 < dist0) {\
            idx += 2<<(j>>1);\
        }\
    }

#define PICK_DXT5_ALPHA_DIST(a0, a1, a_pxl_i, a_tmp, a_scale)\
    a0 = a1 = a_scale[unpacked[a_pxl_i]];\
    for (j = ucc; j < chans_per_tex; j += ucc) {\
        a_tmp = a_scale[unpacked[a_pxl_i + j]];\
        if (a_tmp > a0) a0 = a_tmp;\
        if (a_tmp < a1) a1 = a_tmp;\
    }

/*there are 4 interpolated colors in PICK_DXT5_ALPHA_IDX_0_255 mode
0 = a0               1 = a1
2 = (6*a0 +   a1)/7  3 = (5*a0 + 2*a1)/7
4 = (4*a0 + 3*a1)/7  5 = (3*a0 + 4*a1)/7
6 = (2*a0 + 5*a1)/7  7 = (  a0 + 6*a1)/7*/
#define PICK_DXT5_ALPHA_IDX(a0, a1, a_idx, a_pxl_i, a_tmp, a_dif, a_scale)\
    a_dif = a0 - a1;\
    for (j = 0; j<chans_per_tex; j += ucc) {\
        a_tmp = ((a_scale[unpacked[a_pxl_i + j]] - a1) * 7 + (a_dif >> 1)) / a_dif;\
        if (a_tmp == 0) {\
            a_idx += 1ULL << ((j * 3) / ucc);\
        } else if (a_tmp < 7) {\
            a_idx += ((uint64)(8 - a_tmp)) << ((j * 3) / ucc);\
        }\
    }


/*there are 4 interpolated colors in PICK_DXT5_ALPHA_IDX_0_255 mode
0 =  a0              1 = a1
2 = (4*a0 +   a1)/5  3 = (3*a0 + 2*a1)/5
4 = (2*a0 + 3*a1)/5  5 = (  a0 + 4*a1)/5
6 =  0               7 = 255*/
#define PICK_DXT5_ALPHA_IDX_0_255(a0, a1, a_idx, a_pxl_i, a_tmp, a_dif, a_scale)\
    a0 = 255; a1 = 0;\
    for (j = 0; j<chans_per_tex; j += ucc) {\
        a_tmp = a_scale[unpacked[a_pxl_i + j]];\
        if ((a_tmp < a0) && (a_tmp != 0))   a0 = a_tmp;\
        if ((a_tmp > a1) && (a_tmp != 255)) a1 = a_tmp;\
    }\
    a_dif = a1 - a0;\
    if (!a1) {\
        a0 = a_dif = 0;\
        a1 = 255;\
    }\
    for (j = 0; j<chans_per_tex; j += ucc) {\
        a_tmp = a_scale[unpacked[a_pxl_i + j]];\
        if (a_tmp == 0) {\
            a_idx += 6ULL << ((j * 3) / ucc);\
        } else if (a_tmp == 255) {\
            a_idx += 7ULL << ((j * 3) / ucc);\
        } else if (a_dif) {\
            a_tmp = ((a_tmp - a0) * 5 + (a_dif >> 1)) / a_dif;\
            if (a_tmp == 5) {\
                a_idx += 1ULL << ((j * 3) / ucc);\
            } else if (a_tmp > 0) {\
                a_idx += ((uint64)(a_tmp + 1)) << ((j * 3) / ucc);\
            }\
        }\
    }

static void unpack_dxt1_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *r_scale_buf, Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint16 chan0, sint16 chan1, sint16 chan2, sint16 chan3)
{
    uint8 *unpacked_pix;
    uint32 *packed_tex;
    uint8 *r_scale, *g_scale, *b_scale;
    uint16 color0, color1;
    uint32 color_idx;

    sint8 ucc=4;  // unpacked channel count
    sint8 chans_per_tex = ucc*pix_per_tex;
    uint64 i, j, max_i, pxl_i, off;
    uint8 c_0[4]={255,0,0,0}, c_1[4]={255,0,0,0},
                  c_2[4]={255,0,0,0}, c_3[4]={255,0,0,0};
    uint8 transparent[4]={0,0,0,0};
    uint8 *color, *colors[4];

    unpacked_pix = (uint8*)unpacked_pix_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;
    b_scale = (uint8 *)b_scale_buf->buf;

    colors[0] = c_0; colors[1] = c_1; colors[2] = c_2; colors[3] = c_3;
    packed_tex = (uint32 *) packed_tex_buf->buf;
    max_i = (unpacked_pix_buf->len)/chans_per_tex;

    //loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        j = i<<1;

        /*If the format DXT1 then the two entries in the array
        are the colors and the color indexing in that order.
        Also, if the first color is a larger integer
        then color key transparency is NOT used.*/
        READ_DXT_COLORS(j);

        for (j=0; j<pix_per_tex; j++) {
            UNPACK_DXT_COLORS();
            if (chan0 >= 0) unpacked_pix[off] = color[chan0];
            else            unpacked_pix[off] = 0xFF;
        }
    }
}


static void unpack_dxt2_3_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint16 chan0, sint16 chan1, sint16 chan2, sint16 chan3)
{
    uint8 *unpacked_pix;
    uint32 *packed_tex;
    uint8 *a_scale, *r_scale, *g_scale, *b_scale;
    uint16 color0, color1;
    uint32 color_idx;

    sint8 ucc=4;  // unpacked channel count
    sint8 chans_per_tex = ucc*pix_per_tex;
    uint64 i, j, max_i, pxl_i, off, alpha;
    uint8 c_0[4]={0,0,0,0}, c_1[4]={0,0,0,0},
                  c_2[4]={0,0,0,0}, c_3[4]={0,0,0,0};
    uint8 transparent[4]={0,0,0,0};
    uint8 *color, *colors[4], a, unpack_max = 0xFF;

    unpacked_pix = (uint8 *)unpacked_pix_buf->buf;
    a_scale = (uint8 *)a_scale_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;
    b_scale = (uint8 *)b_scale_buf->buf;

    colors[0] = c_0; colors[1] = c_1; colors[2] = c_2; colors[3] = c_3;
    packed_tex = (uint32 *) packed_tex_buf->buf;
    max_i = (unpacked_pix_buf->len)/chans_per_tex;

    //loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        j = i<<2;

        alpha = ((uint64)packed_tex[j+1]<<32) + packed_tex[j];
        READ_DXT_COLORS(j+2);

        for (j=0; j<pix_per_tex; j++) {
            a = (alpha>>(j<<2))&15;
            UNPACK_DXT2345_COLORS();
        }
    }
}


static void unpack_dxt4_5_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint16 chan0, sint16 chan1, sint16 chan2, sint16 chan3)
{
    uint8 *unpacked_pix;
    uint32 *packed_tex;
    uint8 *a_scale, *r_scale, *g_scale, *b_scale;
    uint16 color0, color1;
    uint32 color_idx;

    sint8 ucc=4;  // unpacked channel count
    sint8 chans_per_tex = ucc*pix_per_tex;
    uint64 i, j, max_i, pxl_i, off, alpha_idx;
    uint8 c_0[4]={0,0,0,0}, c_1[4]={0,0,0,0},
                  c_2[4]={0,0,0,0}, c_3[4]={0,0,0,0};
    uint8 transparent[4]={0,0,0,0}, alpha0, alpha1, a_lookup[8];
    uint8 *color, *colors[4], a, unpack_max=0xFF;

    unpacked_pix = (uint8 *)unpacked_pix_buf->buf;
    a_scale = (uint8 *)a_scale_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;
    b_scale = (uint8 *)b_scale_buf->buf;

    colors[0] = c_0; colors[1] = c_1; colors[2] = c_2; colors[3] = c_3;
    packed_tex = (uint32 *) packed_tex_buf->buf;
    max_i = (unpacked_pix_buf->len)/chans_per_tex;

    //loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        j = i<<2;

        READ_DXT5_A(a_lookup, alpha_idx, j, alpha0, alpha1);
        READ_DXT_COLORS(j+2);
        for (j=0; j<pix_per_tex; j++) {
            a = a_lookup[(alpha_idx>>(3*j))&7];
            UNPACK_DXT2345_COLORS();
        }
    }
}


static void unpack_dxt5a_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *scale0_buf, Py_buffer *scale1_buf,
    Py_buffer *scale2_buf, Py_buffer *scale3_buf,
    sint8 ucc, sint8 scc, sint8 pix_per_tex, sint16 *chans)
{
    uint8 *unpacked_pix;
    uint32 *packed_tex;
    uint8 *scales[4], *scale;

    sint16 chans_per_tex=ucc*pix_per_tex, chan;
    uint64 i, j, max_i, pxl_i, idx;
    uint8 val0, val1, lookup[8];

    unpacked_pix = (uint8 *)unpacked_pix_buf->buf;
    scales[0] = (uint8 *)scale0_buf->buf;
    scales[1] = (uint8 *)scale1_buf->buf;
    scales[2] = (uint8 *)scale2_buf->buf;
    scales[3] = (uint8 *)scale3_buf->buf;

    packed_tex = (uint32 *) packed_tex_buf->buf;
    max_i = unpacked_pix_buf->len/pix_per_tex;

    //loop through each texel
    for (i=0; i < max_i; i++) {
        //chan = chans[i%ucc];
        chan = i%ucc;
        pxl_i = (i/scc)*chans_per_tex + chan;
        scale = scales[chan];

        READ_DXT5_A(lookup, idx, i << 1, val0, val1);
        for (j=0; j<pix_per_tex; j++) {
            unpacked_pix[pxl_i + j*ucc] = scale[lookup[(idx>>(j*3))&7]];
        }
    }
}


static void unpack_dxn_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *r_scale_buf, Py_buffer *g_scale_buf,
    sint8 ucc, sint16 pix_per_tex, sint16 chan1, sint16 chan2, sint16 chan3)
{
    uint8 *unpacked_pix;
    uint32 *packed_tex;
    uint8 *r_scale, *g_scale;

    sint8 chans_per_tex=ucc*pix_per_tex;
    uint64 i, j, max_i, pxl_i;
    uint64 r_idx, g_idx, r_pxl_i, g_pxl_i, b_pxl_i;
    uint8 r0, r1, r_lookup[8], g0, g1, g_lookup[8];
    uint8 x, y, r, g, b, colors[4]={0,0,0,0};
    double d, n_len;

    unpacked_pix = (uint8 *)unpacked_pix_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;

    packed_tex = (uint32 *) packed_tex_buf->buf;
    max_i = (unpacked_pix_buf->len)/chans_per_tex;

    //loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        r_pxl_i = pxl_i + 1;
        g_pxl_i = pxl_i + 2;
        b_pxl_i = pxl_i + 3;
        j = i<<2;

        READ_DXT5_A(g_lookup, g_idx, j, g0, g1);
        READ_DXT5_A(r_lookup, r_idx, j + 2, r0, r1);

        for (j=0; j<pix_per_tex; j++) {
            g = y = g_scale[g_lookup[(g_idx>>(j*3))&7]];
            r = x = r_scale[r_lookup[(r_idx>>(j*3))&7]];
            CALC_B_NORMALIZE(r, g, b, x, y, 0x7F, 0x80, uint8);
            colors[1] = r; colors[2] = g; colors[3] = b;
            if (chan1 >= 0) unpacked_pix[r_pxl_i + j*ucc] = colors[chan1];
            if (chan2 >= 0) unpacked_pix[g_pxl_i + j*ucc] = colors[chan2];
            if (chan3 >= 0) unpacked_pix[b_pxl_i + j*ucc] = colors[chan3];
        }
    }
}


static void unpack_ctx1_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *r_scale_buf, Py_buffer *g_scale_buf,
    sint8 ucc, sint16 pix_per_tex, sint16 chan1, sint16 chan2, sint16 chan3)
{
    uint8 *unpacked_pix;
    uint32 *packed_tex;
    uint8 *r_scale, *g_scale;
    uint32 color_idx;

    sint8 chans_per_tex = ucc*pix_per_tex;
    uint64 i, j, max_i, pxl_i, off;
    uint8 c_0[4] = {255,0,0,0}, c_1[4] = {255,0,0,0},
                  c_2[4] = {255,0,0,0}, c_3[4] = {255,0,0,0};
    uint8 *color, *colors[4], x, y;
    double d, n_len;
    colors[0] = c_0; colors[1] = c_1; colors[2] = c_2; colors[3] = c_3;

    packed_tex = (uint32 *)packed_tex_buf->buf;
    unpacked_pix = (uint8*)unpacked_pix_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;

    max_i = (unpacked_pix_buf->len) / chans_per_tex;

    //loop through each texel
    for (i = 0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        j = i << 1;

        c_0[1] = r_scale[packed_tex[j] & 255];
        c_0[2] = g_scale[(packed_tex[j] >> 8) & 255];
        CALC_B_NORMALIZE(c_0[1], c_0[2], c_0[3], x, y, 0x7F, 0x80, uint8);

        c_1[1] = r_scale[(packed_tex[j] >> 16) & 255];
        c_1[2] = g_scale[(packed_tex[j] >> 24) & 255];
        CALC_B_NORMALIZE(c_1[1], c_1[2], c_1[3], x, y, 0x7F, 0x80, uint8);

        color_idx = packed_tex[j + 1];

        c_2[1] = (c_0[1] * 2 + c_1[1]) / 3;
        c_2[2] = (c_0[2] * 2 + c_1[2]) / 3;
        c_2[3] = (c_0[3] * 2 + c_1[3]) / 3;

        c_3[1] = (c_0[1] + 2 * c_1[1]) / 3;
        c_3[2] = (c_0[2] + 2 * c_1[2]) / 3;
        c_3[3] = (c_0[3] + 2 * c_1[3]) / 3;
        for (j = 0; j<pix_per_tex; j++) {
            UNPACK_DXT_COLORS();
        }
    }
}


static void unpack_v16u16_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *r_scale_buf, Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 ucc, sint16 r_chan, sint16 g_chan, sint16 b_chan)
{
    uint32 *packed_pix;
    uint8 *unpacked_pix;
    uint8 *r_scale, *g_scale, *b_scale;

    uint64 i, max_i, pxl_i;
    sint32 u, v, w;
    double d, n_len;

    packed_pix   = (uint32 *)packed_pix_buf->buf;
    unpacked_pix = (uint8 *)unpacked_pix_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;
    b_scale = (uint8 *)g_scale_buf->buf;

    max_i = (unpacked_pix_buf->len)/ucc;
    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i*ucc;

        u = packed_pix[i] & 0xFFff;
        v = packed_pix[i] >> 16;
        CALC_W_NORMALIZE(u, v, w, 0x7Fff, 0xFFff, 0x8000, sint32);
        unpacked_pix[pxl_i + r_chan] = r_scale[u + 0x8000];
        unpacked_pix[pxl_i + g_chan] = g_scale[v + 0x8000];
        unpacked_pix[pxl_i + b_chan] = b_scale[w + 0x8000];
    }
}


static void unpack_v8u8_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *r_scale_buf, Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 ucc, sint16 chan1, sint16 chan2, sint16 chan3)
{
    uint16 *packed_pix;
    uint8  *unpacked_pix;
    uint8 *r_scale, *g_scale, *b_scale, colors[4]={0,0,0,0};

    uint64 i, max_i, pxl_i;
    sint16 u, v, w;
    double d, n_len;

    packed_pix   = (uint16 *)packed_pix_buf->buf;
    unpacked_pix = (uint8 *)unpacked_pix_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;
    b_scale = (uint8 *)g_scale_buf->buf;

    max_i = (unpacked_pix_buf->len)/ucc;
    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i*ucc;

        u = packed_pix[i] & 0xFF;
        v = packed_pix[i] >> 8;
        CALC_W_NORMALIZE(u, v, w, 0x7F, 0xFF, 0x80, sint16);
        colors[1] = r_scale[u + 0x80];
        colors[2] = g_scale[v + 0x80];
        colors[3] = b_scale[w + 0x80];
        if (chan1 >= 0) unpacked_pix[pxl_i + 1] = colors[chan1];
        if (chan2 >= 0) unpacked_pix[pxl_i + 2] = colors[chan2];
        if (chan3 >= 0) unpacked_pix[pxl_i + 3] = colors[chan3];
    }
}


static void pack_dxt1_8(
    Py_buffer *packed_tex_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *r_scale_buf, Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 can_have_alpha, uint8 a_cutoff)
{
    uint32 *packed;
    uint8 *unpacked;
    uint8 *r_scale, *g_scale, *b_scale;

    sint8 ucc=4;  // unpacked channel count
    sint8 chans_per_tex = ucc*pix_per_tex, make_alpha=0;
    uint64 c_0i, c_1i, i, j, k, max_i;
    uint64 pxl_i, r_pxl_i, g_pxl_i, b_pxl_i;
    uint8 c_0[4], c_1[4], c_2[4], c_3[4], tmp[4], r, g, b;
    uint16 color0, color1, tmp_color;
    uint32 idx, dist0, dist1, dist2, dist3;

    packed   = (uint32 *)packed_tex_buf->buf;
    unpacked = (uint8 *)unpacked_pix_buf->buf;
    r_scale  = (uint8 *)r_scale_buf->buf;
    g_scale  = (uint8 *)g_scale_buf->buf;
    b_scale  = (uint8 *)b_scale_buf->buf;

    max_i = (packed_tex_buf->len)/8; // 8 bytes per texel

    //loop through each texel
    for (i=0; i < max_i; i++) {
        idx = dist0 = c_0i = c_1i = 0;
        pxl_i = i*chans_per_tex;
        r_pxl_i = pxl_i + 1;
        g_pxl_i = pxl_i + 2;
        b_pxl_i = pxl_i + 3;

        // figure out if we are using color key transparency for this pixel
        // by seeing if any of the alpha values are below the cutoff bias
        if (can_have_alpha) {
            make_alpha = 0;
            for (j = 0; j<chans_per_tex; j += ucc) {
                if (unpacked[pxl_i + j] < a_cutoff) {
                    make_alpha = 1;
                    break;
                }
            }
        }
        // compare distance between all pixels and find the two furthest apart
        // (we are actually comparing the area of the distance as it's faster)
        PICK_DXT_PALETTE_DIST();

        // store furthest apart colors for use
        PACK_DXT_COLOR(c_0, c_0i, color0);
        PACK_DXT_COLOR(c_1, c_1i, color1);

        if (!make_alpha && (color0 == color1)) {
            ;// do nothing except save the colors to the array
        } else {
            //  make_alpha == no_alpha_for_texel
            if (make_alpha == (color0 > color1)) {SWAP_COLORS();}

            if (color0 > color1) {CALC_DXT_IDX_NO_ALPHA(idx);}
            else                 {CALC_DXT_IDX_ALPHA(idx);}
        }
        packed[i << 1] = (color1 << 16) + color0;
        packed[(i << 1) + 1] = idx;
    }
}


static void pack_dxt2_3_8(
    Py_buffer *packed_tex_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf, sint8 pix_per_tex)
{
    uint32 *packed;
    uint8 *unpacked;
    uint8 *a_scale, *r_scale, *g_scale, *b_scale;

    sint8 ucc=4;  // unpacked channel count
    sint8 chans_per_tex = ucc*pix_per_tex, make_alpha = 0, a_cutoff = 0;
    uint64 c_0i, c_1i, i, j, k, max_i, alpha;
    uint64 pxl_i, r_pxl_i, g_pxl_i, b_pxl_i;
    uint8 c_0[4], c_1[4], c_2[4], c_3[4], tmp[4], r, g, b;
    uint16 color0, color1, tmp_color;
    uint32 idx, dist0, dist1, dist2, dist3;

    packed = (uint32 *) packed_tex_buf->buf;
    unpacked = (uint8 *)unpacked_pix_buf->buf;
    a_scale = (uint8 *)a_scale_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;
    b_scale = (uint8 *)b_scale_buf->buf;

    max_i = (packed_tex_buf->len)/16; // 16 bytes per texel

    //loop through each texel
    for (i=0; i < max_i; i++) {
        idx = dist0 = alpha = c_0i = c_1i = 0;
        pxl_i = i*chans_per_tex;
        r_pxl_i = pxl_i+1;
        g_pxl_i = pxl_i+2;
        b_pxl_i = pxl_i+3;

        // calculate the alpha
        for (j=0; j<chans_per_tex; j+=ucc) {
            alpha += ((uint64)a_scale[unpacked[pxl_i+j]])<<j;
        }
        // compare distance between all pixels and find the two furthest apart
        // (we are actually comparing the area of the distance as it's faster)
        PICK_DXT_PALETTE_DIST();

        // store furthest apart colors for use
        PACK_DXT_COLOR(c_0, c_0i, color0);
        PACK_DXT_COLOR(c_1, c_1i, color1);

        if (color0 != color1) {
            if (color0 < color1) {SWAP_COLORS();}
            CALC_DXT_IDX_NO_ALPHA(idx);
        }
        packed[i << 2]       = alpha & 0xFFffFFff;
        packed[(i << 2) + 1] = alpha >> 32;
        packed[(i << 2) + 2] = (color1 << 16) + color0;
        packed[(i << 2) + 3] = idx;
    }
}


static void pack_dxt4_5_8(
    Py_buffer *packed_tex_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf, sint8 pix_per_tex)
{
    uint32 *packed;
    uint8 *unpacked;
    uint8 *a_scale, *r_scale, *g_scale, *b_scale;

    sint8 ucc=4;  // unpacked channel count
    sint8 chans_per_tex = ucc*pix_per_tex, make_alpha = 0, a_cutoff = 0;
    uint64 c_0i, c_1i, i, j, k, max_i;
    uint64 a_idx, pxl_i, r_pxl_i, g_pxl_i, b_pxl_i;
    uint8 c_0[4], c_1[4], c_2[4], c_3[4], tmp[4];
    uint8 a0, a1, a_tmp, a_dif, r, g, b;
    uint16 color0, color1, tmp_color;
    uint32 idx, dist0, dist1, dist2, dist3;

    packed = (uint32 *) packed_tex_buf->buf;
    unpacked = (uint8 *)unpacked_pix_buf->buf;
    a_scale = (uint8 *)a_scale_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;
    b_scale = (uint8 *)b_scale_buf->buf;

    max_i = (packed_tex_buf->len)/16; // 16 bytes per texel

    //loop through each texel
    for (i=0; i < max_i; i++) {
        idx = dist0 = c_0i = c_1i = a_idx = 0;
        pxl_i = i*chans_per_tex;
        r_pxl_i = pxl_i+1;
        g_pxl_i = pxl_i+2;
        b_pxl_i = pxl_i+3;

        // calculate the min and max alpha values
        PICK_DXT5_ALPHA_DIST(a0, a1, pxl_i, a_tmp, a_scale);

        // calculate the alpha indexing
        if (a0 == a1) {
            ; // values are the same; do nothing
        } else if (a1 && (a0 != 255)) {
            PICK_DXT5_ALPHA_IDX(a0, a1, a_idx, pxl_i, a_tmp, a_dif, a_scale);
        } else {
            // max is 255 or min is 0
            PICK_DXT5_ALPHA_IDX_0_255(a0, a1, a_idx, pxl_i, a_tmp, a_dif, a_scale);
        }

        // compare distance between all pixels and find the two furthest apart
        // (we are actually comparing the area of the distance as it's faster)
        PICK_DXT_PALETTE_DIST();

        // store furthest apart colors for use
        PACK_DXT_COLOR(c_0, c_0i, color0);
        PACK_DXT_COLOR(c_1, c_1i, color1);

        if (color0 != color1) {
            if (color0 < color1) { SWAP_COLORS(); }
            CALC_DXT_IDX_NO_ALPHA(idx);
        }
        packed[i << 2] = ((a_idx & 0xFFff) << 16) + (a1 << 8) + a0;
        packed[(i << 2) + 1] = (a_idx >> 16) & 0xFFffFFff;
        packed[(i << 2) + 2] = (color1 << 16) + color0;
        packed[(i << 2) + 3] = idx;
    }
}


static void pack_v8u8_8(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *u_scale_buf, Py_buffer *v_scale_buf,
    sint8 ucc, sint16 u_chan, sint16 v_chan)
{
    uint16 *packed_pix;
    uint8  *unpacked_pix;
    uint8  *u_scale, *v_scale;
    uint64 i, max_i, pxl_i;

    packed_pix   = (uint16 *)packed_pix_buf->buf;
    unpacked_pix = (uint8 *)unpacked_pix_buf->buf;
    u_scale      = (uint8 *)u_scale_buf->buf;
    v_scale      = (uint8 *)v_scale_buf->buf;

    max_i = (packed_pix_buf->len)/2;
    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i*ucc;
        packed_pix[i] = ((((uint32)v_scale[unpacked_pix[pxl_i + v_chan]]) << 8) +
                                          u_scale[unpacked_pix[pxl_i + u_chan]]) ^ 0x8080;
    }
}


static void pack_v16u16_8(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *u_scale_buf, Py_buffer *v_scale_buf,
    sint8 ucc, sint16 u_chan, sint16 v_chan)
{
    uint32 *packed_pix;
    uint8 *unpacked_pix;
    uint16 *u_scale, *v_scale;
    uint64 i, max_i, pxl_i;

    packed_pix   = (uint32  *)packed_pix_buf->buf;
    unpacked_pix = (uint8  *)unpacked_pix_buf->buf;
    u_scale      = (uint16 *)u_scale_buf->buf;
    v_scale      = (uint16 *)v_scale_buf->buf;

    max_i = (packed_pix_buf->len)/4;
    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i*ucc;
        packed_pix[i] = ((((uint32)v_scale[unpacked_pix[pxl_i + v_chan]]) << 16) +
                                          u_scale[unpacked_pix[pxl_i + u_chan]]) ^ 0x80008000UL;
    }
}


/*
    Deep color versions of the above functions
*/


static void unpack_dxt1_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *r_scale_buf, Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 chan0, sint8 chan1, sint8 chan2, sint8 chan3)
{
    uint16 *unpacked_pix;
    uint32 *packed_tex;
    uint16 *r_scale, *g_scale, *b_scale;
    uint16 color0, color1;
    uint32 color_idx;

    sint8 ucc=4;  // unpacked channel count
    sint8 chans_per_tex = ucc*pix_per_tex;
    uint64 i, j, max_i, pxl_i, off;
    uint16 c_0[4]={65535,0,0,0}, c_1[4]={65535,0,0,0},
                   c_2[4]={65535,0,0,0}, c_3[4]={65535,0,0,0};
    uint16 transparent[4]={0,0,0,0};
    uint16 *color, *colors[4];

    unpacked_pix = (uint16 *)unpacked_pix_buf->buf;
    r_scale = (uint16 *)r_scale_buf->buf;
    g_scale = (uint16 *)g_scale_buf->buf;
    b_scale = (uint16 *)b_scale_buf->buf;

    colors[0] = c_0; colors[1] = c_1; colors[2] = c_2; colors[3] = c_3;
    packed_tex = (uint32 *)packed_tex_buf->buf;
    max_i = (unpacked_pix_buf->len)/(2*chans_per_tex);

    //loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        j = i<<1;

        /*If the format DXT1 then the two entries in the array
        are the colors and the color indexing in that order.
        Also, if the first color is a larger integer
        then color key transparency is NOT used.*/
        READ_DXT_COLORS(j);

        for (j=0; j<pix_per_tex; j++) {
            UNPACK_DXT_COLORS();
            if (chan0 >= 0) unpacked_pix[off + chan0] = color[0];
            else            unpacked_pix[off] = 0xFFff;
        }
    }
}

static void unpack_dxt2_3_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 chan0, sint8 chan1, sint8 chan2, sint8 chan3)
{
    uint16 *unpacked_pix;
    uint32 *packed_tex;
    uint16 *a_scale, *r_scale, *g_scale, *b_scale;
    uint16 color0, color1;
    uint32 color_idx;

    sint8 ucc=4;  // unpacked channel count
    sint8 chans_per_tex = ucc*pix_per_tex;
    uint64 i, j, max_i, pxl_i, off, alpha;
    uint16 c_0[4]={0,0,0,0}, c_1[4]={0,0,0,0},
                   c_2[4]={0,0,0,0}, c_3[4]={0,0,0,0};
    uint16 transparent[4]={0,0,0,0};
    uint16 *color, *colors[4], a, unpack_max=0xFFff;

    unpacked_pix = (uint16 *)unpacked_pix_buf->buf;
    a_scale = (uint16 *)a_scale_buf->buf;
    r_scale = (uint16 *)r_scale_buf->buf;
    g_scale = (uint16 *)g_scale_buf->buf;
    b_scale = (uint16 *)b_scale_buf->buf;

    colors[0] = c_0; colors[1] = c_1; colors[2] = c_2; colors[3] = c_3;
    packed_tex = (uint32 *) packed_tex_buf->buf;
    max_i = (unpacked_pix_buf->len)/(2*chans_per_tex);

    //loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        j = i<<2;

        alpha = ((uint64)packed_tex[j+1]<<32) + packed_tex[j];
        READ_DXT_COLORS(j+2);

        for (j=0; j<pix_per_tex; j++) {
            a = (alpha>>(j<<2))&15;
            UNPACK_DXT2345_COLORS();
        }
    }
}


static void unpack_dxt4_5_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 chan0, sint8 chan1, sint8 chan2, sint8 chan3)
{
    uint16 *unpacked_pix;
    uint32 *packed_tex;
    uint16 *a_scale, *r_scale, *g_scale, *b_scale;
    uint16 color0, color1;
    uint32 color_idx;

    sint8 ucc=4;  // unpacked channel count
    sint8 chans_per_tex = ucc*pix_per_tex;
    uint64 i, j, max_i, pxl_i, off, alpha_idx;
    uint16 c_0[4]={0,0,0,0}, c_1[4]={0,0,0,0},
                   c_2[4]={0,0,0,0}, c_3[4]={0,0,0,0};
    uint16 transparent[4]={0,0,0,0};
    uint16 *color, *colors[4], a, unpack_max=0xFFff;
    uint8 alpha0, alpha1, a_lookup[8];

    unpacked_pix = (uint16 *)unpacked_pix_buf->buf;
    a_scale = (uint16 *)a_scale_buf->buf;
    r_scale = (uint16 *)r_scale_buf->buf;
    g_scale = (uint16 *)g_scale_buf->buf;
    b_scale = (uint16 *)b_scale_buf->buf;

    colors[0] = c_0; colors[1] = c_1; colors[2] = c_2; colors[3] = c_3;
    packed_tex = (uint32 *) packed_tex_buf->buf;
    max_i = (unpacked_pix_buf->len)/(2*chans_per_tex);

    //loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        j = i<<2;

        READ_DXT5_A(a_lookup, alpha_idx, j, alpha0, alpha1);
        READ_DXT_COLORS(j+2);
        for (j=0; j<pix_per_tex; j++) {
            a = a_lookup[(alpha_idx>>(3*j))&7];
            UNPACK_DXT2345_COLORS();
        }
    }
}


static void unpack_dxt5a_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *scale0_buf, Py_buffer *scale1_buf,
    Py_buffer *scale2_buf, Py_buffer *scale3_buf,
    sint8 ucc, sint8 scc, sint8 pix_per_tex, sint16 *chans)
{
    uint16 *unpacked_pix;
    uint32  *packed_tex;
    uint16 *scales[4], *scale;

    sint16 chans_per_tex=ucc*pix_per_tex, chan;
    uint64 i, j, max_i, pxl_i, idx;
    uint8 val0, val1, lookup[8];

    unpacked_pix = (uint16 *)unpacked_pix_buf->buf;
    scales[0] = (uint16 *)scale0_buf->buf;
    scales[1] = (uint16 *)scale1_buf->buf;
    scales[2] = (uint16 *)scale2_buf->buf;
    scales[3] = (uint16 *)scale3_buf->buf;

    packed_tex = (uint32 *)packed_tex_buf->buf;
    max_i = unpacked_pix_buf->len/(2*pix_per_tex);

    //loop through each texel
    for (i=0; i < max_i; i++) {
        //chan = chans[i%ucc];
        chan = i%ucc;
        pxl_i = (i/scc)*chans_per_tex + chan;
        scale = scales[chan];

        READ_DXT5_A(lookup, idx, i << 1, val0, val1);
        for (j=0; j<pix_per_tex; j++) {
            unpacked_pix[pxl_i + j*ucc] = scale[lookup[(idx>>(j*3))&7]];
        }
    }
}


static void unpack_dxn_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *r_scale_buf, Py_buffer *g_scale_buf,
    sint8 ucc, sint8 pix_per_tex, sint8 chan1, sint8 chan2, sint8 chan3)
{
    uint16 *unpacked_pix;
    uint32 *packed_tex;
    uint16 *r_scale, *g_scale;

    sint8 chans_per_tex=ucc*pix_per_tex;
    uint64 i, j, max_i, pxl_i;
    uint64 r_idx, g_idx, r_pxl_i, g_pxl_i, b_pxl_i;
    uint16 r0, r1, r_lookup[8], g0, g1, g_lookup[8];
    uint16  x, y, r, g, b, colors[4]={0,0,0,0};
    double d, n_len;

    unpacked_pix = (uint16 *)unpacked_pix_buf->buf;
    r_scale      = (uint16 *)r_scale_buf->buf;
    g_scale      = (uint16 *)g_scale_buf->buf;

    packed_tex = (uint32 *)packed_tex_buf->buf;
    max_i = (unpacked_pix_buf->len)/(2*chans_per_tex);

    //loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        r_pxl_i = pxl_i + 1;
        g_pxl_i = pxl_i + 2;
        b_pxl_i = pxl_i + 3;
        j = i<<2;

        READ_DXT5_A(g_lookup, g_idx, j, g0, g1);
        READ_DXT5_A(r_lookup, r_idx, j + 2, r0, r1);

        for (j=0; j<pix_per_tex; j++) {
            g = y = g_scale[g_lookup[(g_idx>>(j*3))&7]];
            r = x = r_scale[r_lookup[(r_idx>>(j*3))&7]];
            CALC_B_NORMALIZE(r, g, b, x, y, 0x7Fff, 0x8000, uint16);
            colors[1] = r; colors[2] = g; colors[3] = b;
            if (chan1 >= 0) unpacked_pix[r_pxl_i + j*ucc] = colors[chan1];
            if (chan2 >= 0) unpacked_pix[g_pxl_i + j*ucc] = colors[chan2];
            if (chan3 >= 0) unpacked_pix[b_pxl_i + j*ucc] = colors[chan3];
        }
    }
}


static void unpack_ctx1_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_tex_buf,
    Py_buffer *r_scale_buf, Py_buffer *g_scale_buf,
    sint8 ucc, sint8 pix_per_tex, sint8 chan1, sint8 chan2, sint8 chan3)
{
    uint16 *unpacked_pix;
    uint32  *packed_tex;
    uint16 *r_scale, *g_scale;
    uint32 color_idx;

    sint8 chans_per_tex = ucc*pix_per_tex;
    uint64 i, j, max_i, pxl_i, off;
    uint16 c_0[4] = {65535,0,0,0}, c_1[4] = {65535,0,0,0},
                   c_2[4] = {65535,0,0,0}, c_3[4] = {65535,0,0,0};
    uint16 *color, *colors[4], x, y;
    double d, n_len;

    unpacked_pix = (uint16*)unpacked_pix_buf->buf;
    r_scale = (uint16 *)r_scale_buf->buf;
    g_scale = (uint16 *)g_scale_buf->buf;

    colors[0] = c_0; colors[1] = c_1; colors[2] = c_2; colors[3] = c_3;
    packed_tex = (uint32 *)packed_tex_buf->buf;
    max_i = (unpacked_pix_buf->len) / (2*chans_per_tex);

    //loop through each texel
    for (i = 0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        j = i << 1;

        c_0[1] = r_scale[packed_tex[j] & 255];
        c_0[2] = g_scale[(packed_tex[j] >> 8) & 255];
        CALC_B_NORMALIZE(c_0[1], c_0[2], c_0[3], x, y, 0x7Fff, 0x8000, uint16);

        c_1[1] = r_scale[(packed_tex[j] >> 16) & 255];
        c_1[2] = g_scale[(packed_tex[j] >> 24) & 255];
        CALC_B_NORMALIZE(c_1[1], c_1[2], c_1[3], x, y, 0x7Fff, 0x8000, uint16);

        color_idx = packed_tex[j + 1];

        c_2[1] = (c_0[1] * 2 + c_1[1]) / 3;
        c_2[2] = (c_0[2] * 2 + c_1[2]) / 3;
        c_2[3] = (c_0[3] * 2 + c_1[3]) / 3;

        c_3[1] = (c_0[1] + 2 * c_1[1]) / 3;
        c_3[2] = (c_0[2] + 2 * c_1[2]) / 3;
        c_3[3] = (c_0[3] + 2 * c_1[3]) / 3;
        for (j = 0; j<pix_per_tex; j++) {
            UNPACK_DXT_COLORS();
        }
    }
}


static void unpack_v8u8_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *r_scale_buf, Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 ucc, sint8 chan1, sint8 chan2, sint8 chan3)
{
    uint16 *packed_pix;
    uint16 *unpacked_pix;
    uint16 *r_scale, *g_scale, *b_scale, colors[4]={0,0,0,0};

    uint64 i, max_i, pxl_i;
    sint16 u, v, w;
    double d, n_len;

    packed_pix   = (uint16 *)packed_pix_buf->buf;
    unpacked_pix = (uint16 *)unpacked_pix_buf->buf;
    r_scale = (uint16 *)r_scale_buf->buf;
    g_scale = (uint16 *)g_scale_buf->buf;
    b_scale = (uint16 *)g_scale_buf->buf;

    max_i = (unpacked_pix_buf->len)/(ucc*2);

    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i*ucc;

        u = packed_pix[i] & 0xFF;
        v = packed_pix[i] >> 8;
        CALC_W_NORMALIZE(u, v, w, 0x7F, 0xFF, 0x80, sint16);
        colors[1] = r_scale[u + 0x80]; 
        colors[2] = g_scale[v + 0x80]; 
        colors[3] = b_scale[w + 0x80];
        if (chan1 >= 0) unpacked_pix[pxl_i + 1] = colors[chan1];
        if (chan2 >= 0) unpacked_pix[pxl_i + 2] = colors[chan2];
        if (chan3 >= 0) unpacked_pix[pxl_i + 3] = colors[chan3];
    }
}


static void unpack_v16u16_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *r_scale_buf, Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 ucc, sint8 chan1, sint8 chan2, sint8 chan3)
{
    uint32  *packed_pix;
    uint16 *unpacked_pix;
    uint16 *r_scale, *g_scale, *b_scale, colors[4]={0,0,0,0};

    uint64 i, max_i, pxl_i;
    sint32 u, v, w;
    double d, n_len;

    packed_pix   = (uint32 *)packed_pix_buf->buf;
    unpacked_pix = (uint16 *)unpacked_pix_buf->buf;
    r_scale = (uint16 *)r_scale_buf->buf;
    g_scale = (uint16 *)g_scale_buf->buf;
    b_scale = (uint16 *)g_scale_buf->buf;

    max_i = (unpacked_pix_buf->len)/(ucc*2);
    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i*ucc;

        u = packed_pix[i] & 0xFFff;
        v = packed_pix[i] >> 16;
        CALC_W_NORMALIZE(u, v, w, 0x7Fff, 0xFFff, 0x8000, sint32);
        colors[1] = r_scale[u + 0x8000]; 
        colors[2] = g_scale[v + 0x8000]; 
        colors[3] = b_scale[w + 0x8000];
        if (chan1 >= 0) unpacked_pix[pxl_i + 1] = colors[chan1];
        if (chan2 >= 0) unpacked_pix[pxl_i + 2] = colors[chan2];
        if (chan3 >= 0) unpacked_pix[pxl_i + 3] = colors[chan3];
    }
}


static void pack_dxt1_16(
    Py_buffer *packed_tex_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *r_scale_buf, Py_buffer *g_scale_buf,  Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 can_have_alpha, uint16 a_cutoff)
{
    uint32 *packed;
    uint16 *unpacked;
    uint8 *r_scale, *g_scale, *b_scale;

    sint8 ucc=4;  // unpacked channel count
    sint8 chans_per_tex=ucc*pix_per_tex, make_alpha=0;
    uint64 c_0i, c_1i, i, j, k, max_i;
    uint64 pxl_i, r_pxl_i, g_pxl_i, b_pxl_i;
    uint16 c_0[4], c_1[4], c_2[4], c_3[4], tmp[4], r, g, b;
    uint16 color0, color1, tmp_color;
    uint32 idx;
    uint64 dist0, dist1, dist2, dist3;

    packed = (uint32 *) packed_tex_buf->buf;
    unpacked = (uint16 *)unpacked_pix_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;
    b_scale = (uint8 *)b_scale_buf->buf;

    max_i = (packed_tex_buf->len)/8; // 8 bytes per texel

    //loop through each texel
    for (i=0; i < max_i; i++) {
        idx = dist0 = c_0i = c_1i = 0;
        pxl_i = i*chans_per_tex;
        r_pxl_i = pxl_i + 1;
        g_pxl_i = pxl_i + 2;
        b_pxl_i = pxl_i + 3;

        // figure out if we are using color key transparency for this pixel
        // by seeing if any of the alpha values are below the cutoff bias
        if (can_have_alpha) {
            make_alpha = 0;
            for (j = 0; j<chans_per_tex; j += ucc) {
                if (unpacked[pxl_i + j] < a_cutoff) {
                    make_alpha = 1;
                    break;
                }
            }
        }
        // compare distance between all pixels and find the two furthest apart
        // (we are actually comparing the area of the distance as it's faster)
        PICK_DXT_PALETTE_DIST();

        // store furthest apart colors for use
        PACK_DXT_COLOR(c_0, c_0i, color0);
        PACK_DXT_COLOR(c_1, c_1i, color1);

        if (!make_alpha && (color0 == color1)) {
            ;// do nothing except save the colors to the array
        } else {
            //  make_alpha == no_alpha_for_texel
            if (make_alpha == (color0 > color1)) {SWAP_COLORS();}

            if (color0 > color1) {CALC_DXT_IDX_NO_ALPHA(idx);}
            else                 {CALC_DXT_IDX_ALPHA(idx);}
        }
        packed[i << 1] = (color1 << 16) + color0;
        packed[(i << 1) + 1] = idx;
    }
}


static void pack_dxt2_3_16(
    Py_buffer *packed_tex_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf, sint8 pix_per_tex)
{
    uint32 *packed;
    uint16 *unpacked;
    uint8 *a_scale, *r_scale, *g_scale, *b_scale;

    sint8 ucc=4;  // unpacked channel count
    sint8 chans_per_tex = ucc*pix_per_tex, make_alpha = 0, a_cutoff = 0;
    uint64 c_0i, c_1i, i, j, k, max_i, alpha;
    uint64 pxl_i, r_pxl_i, g_pxl_i, b_pxl_i;
    uint16 c_0[4], c_1[4], c_2[4], c_3[4], tmp[4], r, g, b;
    uint16 color0, color1, tmp_color;
    uint32 idx;
    uint64 dist0, dist1, dist2, dist3;

    packed = (uint32 *) packed_tex_buf->buf;
    unpacked = (uint16 *)unpacked_pix_buf->buf;
    a_scale = (uint8 *)a_scale_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;
    b_scale = (uint8 *)b_scale_buf->buf;

    max_i = (packed_tex_buf->len)/16; // 8 bytes per texel

    //loop through each texel
    for (i=0; i < max_i; i++) {
        idx = dist0 = alpha = c_0i = c_1i = 0;
        pxl_i = i*chans_per_tex;
        r_pxl_i = pxl_i+1;
        g_pxl_i = pxl_i+2;
        b_pxl_i = pxl_i+3;

        // calculate the alpha
        for (j=0; j<chans_per_tex; j+=ucc) {
            alpha += ((uint64)a_scale[unpacked[pxl_i+j]])<<j;
        }
        // compare distance between all pixels and find the two furthest apart
        // (we are actually comparing the area of the distance as it's faster)
        PICK_DXT_PALETTE_DIST();

        // store furthest apart colors for use
        PACK_DXT_COLOR(c_0, c_0i, color0);
        PACK_DXT_COLOR(c_1, c_1i, color1);

        if (color0 != color1) {
            if (color0 < color1) { SWAP_COLORS(); }
            CALC_DXT_IDX_NO_ALPHA(idx);
        }
        packed[i << 2]       = alpha & 0xFFffFFff;
        packed[(i << 2) + 1] = alpha >> 32;
        packed[(i << 2) + 2] = (color1 << 16) + color0;
        packed[(i << 2) + 3] = idx;
    }
}

static void pack_dxt4_5_16(
    Py_buffer *packed_tex_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf, sint8 pix_per_tex)
{
    uint32 *packed;
    uint16 *unpacked;
    uint8 *a_scale, *r_scale, *g_scale, *b_scale;

    sint8 ucc=4;  // unpacked channel count
    sint8 chans_per_tex = ucc*pix_per_tex, make_alpha = 0, a_cutoff = 0;
    uint64 c_0i, c_1i, i, j, k, max_i;
    uint64 a_idx, pxl_i, r_pxl_i, g_pxl_i, b_pxl_i;
    uint16 c_0[4], c_1[4], c_2[4], c_3[4], tmp[4], r, g, b;
    uint16 a0, a1, a_tmp, a_dif;
    uint16 color0, color1, tmp_color;
    uint32 idx;
    uint64 dist0, dist1, dist2, dist3;

    packed = (uint32 *) packed_tex_buf->buf;
    unpacked = (uint16 *)unpacked_pix_buf->buf;
    a_scale = (uint8 *)a_scale_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;
    b_scale = (uint8 *)b_scale_buf->buf;

    max_i = (packed_tex_buf->len)/16; // 16 bytes per texel

    //loop through each texel
    for (i=0; i < max_i; i++) {
        idx = dist0 = c_0i = c_1i = a_idx = 0;
        pxl_i = i*chans_per_tex;
        r_pxl_i = pxl_i+1;
        g_pxl_i = pxl_i+2;
        b_pxl_i = pxl_i+3;

        // calculate the min and max alpha values
        PICK_DXT5_ALPHA_DIST(a0, a1, pxl_i, a_tmp, a_scale);

        // calculate the alpha indexing
        if (a0 == a1) {
            ; // values are the same; do nothing
        } else if (a1 && (a0 != 255)) {
            PICK_DXT5_ALPHA_IDX(a0, a1, a_idx, pxl_i, a_tmp, a_dif, a_scale);
        } else {
            // max is 255 or min is 0
            PICK_DXT5_ALPHA_IDX_0_255(a0, a1, a_idx, pxl_i, a_tmp, a_dif, a_scale);
        }

        // compare distance between all pixels and find the two furthest apart
        // (we are actually comparing the area of the distance as it's faster)
        PICK_DXT_PALETTE_DIST();

        // store furthest apart colors for use
        PACK_DXT_COLOR(c_0, c_0i, color0);
        PACK_DXT_COLOR(c_1, c_1i, color1);

        if (color0 != color1) {
            if (color0 < color1) {SWAP_COLORS();}
            CALC_DXT_IDX_NO_ALPHA(idx);
        }
        packed[i << 2] = ((a_idx & 0xFFff) << 16) + (a1 << 8) + a0;
        packed[(i << 2) + 1] = (a_idx >> 16) & 0xFFffFFff;
        packed[(i << 2) + 2] = (color1 << 16) + color0;
        packed[(i << 2) + 3] = idx;
    }
}


static void pack_v8u8_16(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *u_scale_buf, Py_buffer *v_scale_buf,
    sint8 ucc, sint8 u_chan, sint8 v_chan)
{
    uint16 *packed_pix;
    uint16 *unpacked_pix;
    uint8  *u_scale, *v_scale;
    uint64 i, max_i, pxl_i;

    packed_pix   = (uint16 *)packed_pix_buf->buf;
    unpacked_pix = (uint16 *)unpacked_pix_buf->buf;
    u_scale      = (uint8  *)u_scale_buf->buf;
    v_scale      = (uint8  *)v_scale_buf->buf;

    max_i = (packed_pix_buf->len)/2;
    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i*ucc;
        packed_pix[i] = ((((uint16)v_scale[unpacked_pix[pxl_i + v_chan]]) << 8) +
                                           u_scale[unpacked_pix[pxl_i + u_chan]]) ^ 0x8080;
    }
}


static void pack_v16u16_16(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
    Py_buffer *u_scale_buf, Py_buffer *v_scale_buf,
    sint8 ucc, sint8 u_chan, sint8 v_chan)
{
    uint32  *packed_pix;
    uint16 *unpacked_pix;
    uint16 *u_scale, *v_scale;
    uint64 i, max_i, pxl_i;

    packed_pix   = (uint32  *)packed_pix_buf->buf;
    unpacked_pix = (uint16 *)unpacked_pix_buf->buf;
    u_scale      = (uint16 *)u_scale_buf->buf;
    v_scale      = (uint16 *)v_scale_buf->buf;

    max_i = (packed_pix_buf->len)/4;
    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i*ucc;
        packed_pix[i] = ((((uint32)v_scale[unpacked_pix[pxl_i + v_chan]]) << 16) +
                                          u_scale[unpacked_pix[pxl_i + u_chan]]) ^ 0x80008000UL;
    }
}


static PyObject *py_unpack_dxt1(PyObject *self, PyObject *args) {
    Py_buffer bufs[5];
    sint8 pix_per_tex;
    sint16 chans[4];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*bhhhh:unpack_dxt1",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &pix_per_tex,
        &chans[0], &chans[1], &chans[2], &chans[3]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[0].itemsize == 2) {
        unpack_dxt1_16(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], pix_per_tex,
            chans[0], chans[1], chans[2], chans[3]);
    } else {
        unpack_dxt1_8(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], pix_per_tex,
            chans[0], chans[1], chans[2], chans[3]);
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]);
    PyBuffer_Release(&bufs[3]);
    PyBuffer_Release(&bufs[4]);

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_unpack_dxt2_3(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];
    sint8 pix_per_tex;
    sint16 chans[4];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*bhhhh:unpack_dxt2_3",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        &pix_per_tex, &chans[0], &chans[1], &chans[2], &chans[3]))
        return Py_BuildValue("");  // return Py_None while incrementing it

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

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_unpack_dxt4_5(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];
    sint8 pix_per_tex;
    sint16 chans[4];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*bhhhh:unpack_dxt4_5",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        &pix_per_tex, &chans[0], &chans[1], &chans[2], &chans[3]))
        return Py_BuildValue("");  // return Py_None while incrementing it

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

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_unpack_dxt5a(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];
    sint8 pix_per_tex, ucc, scc;
    sint16 chans[4];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*bbbhhhh:unpack_dxt5a",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        &ucc, &scc, &pix_per_tex, &chans[0], &chans[1], &chans[2], &chans[3]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[0].itemsize == 2) {
        unpack_dxt5a_16(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
            ucc, scc, pix_per_tex, chans);
    } else {
        unpack_dxt5a_8(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
            ucc, scc, pix_per_tex, chans);
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]);
    PyBuffer_Release(&bufs[3]);
    PyBuffer_Release(&bufs[4]);
    PyBuffer_Release(&bufs[5]);

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_unpack_dxn(PyObject *self, PyObject *args) {
    Py_buffer bufs[4];
    sint8 pix_per_tex, ucc;
    sint16 chans[3];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*bbhhh:unpack_dxn",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3],
        &ucc, &pix_per_tex, &chans[0], &chans[1], &chans[2]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[0].itemsize == 2) {
        unpack_dxn_16(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3],
            ucc, pix_per_tex, chans[0], chans[1], chans[2]);
    } else {
        unpack_dxn_8(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3],
            ucc, pix_per_tex, chans[0], chans[1], chans[2]);
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]);
    PyBuffer_Release(&bufs[3]);

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_unpack_ctx1(PyObject *self, PyObject *args) {
    Py_buffer bufs[4];
    sint8 pix_per_tex, ucc;
    sint16 chans[3];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*bbhhh:unpack_ctx1",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3],
        &ucc, &pix_per_tex, &chans[0], &chans[1], &chans[2]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[0].itemsize == 2) {
        unpack_ctx1_16(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3],
            ucc, pix_per_tex, chans[0], chans[1], chans[2]);
    } else {
        unpack_ctx1_8(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3],
            ucc, pix_per_tex, chans[0], chans[1], chans[2]);
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]);
    PyBuffer_Release(&bufs[3]);

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_unpack_vu(PyObject *self, PyObject *args) {
    Py_buffer bufs[5];
    sint8 ucc;
    sint16 chans[3];

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*bhhh:unpack_vu",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4],
        &ucc, &chans[0], &chans[1], &chans[2]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[1].itemsize == 2) {  // unpacking v8u8
        if (bufs[0].itemsize == 2) {  // to a16r16g16b16
            unpack_v8u8_16(
                &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4],
                ucc, chans[0], chans[1], chans[2]);
        } else {  // to a8r8g8b8
            unpack_v8u8_8(
                &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4],
                ucc, chans[0], chans[1], chans[2]);
        }
    } else if (bufs[1].itemsize == 4) {  // unpacking v16u16
        if (bufs[0].itemsize == 2) {  // to a16r16g16b16
            unpack_v16u16_16(
                &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4],
                ucc, chans[0], chans[1], chans[2]);
        } else {  // to a8r8g8b8
            unpack_v16u16_8(
                &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4],
                ucc, chans[0], chans[1], chans[2]);
        }
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]);
    PyBuffer_Release(&bufs[3]);
    PyBuffer_Release(&bufs[4]);

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_pack_dxt1(PyObject *self, PyObject *args) {
    Py_buffer bufs[5];
    sint8 pix_per_tex, can_have_alpha;
    uint16 a_cutoff;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*bbH:pack_dxt1",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4],
        &pix_per_tex, &can_have_alpha, &a_cutoff))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[1].itemsize == 2) {
        pack_dxt1_16(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4],
            pix_per_tex, can_have_alpha, a_cutoff);
    } else {
        pack_dxt1_8(
            &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4],
            pix_per_tex, can_have_alpha, a_cutoff&0xFF);
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]);
    PyBuffer_Release(&bufs[3]);
    PyBuffer_Release(&bufs[4]);

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_pack_dxt2_3(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];
    sint8 pix_per_tex;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*b:pack_dxt2_3",
        &bufs[0], &bufs[1],
        &bufs[2], &bufs[3], &bufs[4], &bufs[5], &pix_per_tex))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[1].itemsize == 2) {
        pack_dxt2_3_16(
            &bufs[0], &bufs[1],
            &bufs[2], &bufs[3], &bufs[4], &bufs[5], pix_per_tex);
    } else {
        pack_dxt2_3_8(
            &bufs[0], &bufs[1],
            &bufs[2], &bufs[3], &bufs[4], &bufs[5], pix_per_tex);
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]);
    PyBuffer_Release(&bufs[3]);
    PyBuffer_Release(&bufs[4]);
    PyBuffer_Release(&bufs[5]);

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_pack_dxt4_5(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];
    sint8 pix_per_tex;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*b:pack_dxt4_5",
        &bufs[0], &bufs[1],
        &bufs[2], &bufs[3], &bufs[4], &bufs[5], &pix_per_tex))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[1].itemsize == 2) {
        pack_dxt4_5_16(
            &bufs[0], &bufs[1],
            &bufs[2], &bufs[3], &bufs[4], &bufs[5], pix_per_tex);
    } else {
        pack_dxt4_5_8(
            &bufs[0], &bufs[1],
            &bufs[2], &bufs[3], &bufs[4], &bufs[5], pix_per_tex);
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]);
    PyBuffer_Release(&bufs[3]);
    PyBuffer_Release(&bufs[4]);
    PyBuffer_Release(&bufs[5]);

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_pack_vu(PyObject *self, PyObject *args) {
    Py_buffer bufs[4];
    sint16 u_chan, v_chan, ucc;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*bhh:pack_vu",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &ucc, &u_chan, &v_chan))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[1].itemsize == 2) {  // packing a16r16g16b16
        if (bufs[0].itemsize == 4) {  // to v16u16
            pack_v16u16_16(&bufs[0], &bufs[1], &bufs[2], &bufs[3], ucc, u_chan, v_chan);
        } else {  // to v8u8
            pack_v8u8_16(&bufs[0], &bufs[1], &bufs[2], &bufs[3], ucc, u_chan, v_chan);
        }
    } else if (bufs[1].itemsize == 1) {  // packing a8r8g8b8
        if (bufs[0].itemsize == 4) {  // to v16u16
            pack_v16u16_8(&bufs[0], &bufs[1], &bufs[2], &bufs[3], ucc, u_chan, v_chan);
        } else {  // to v8u8
            pack_v8u8_8(&bufs[0], &bufs[1], &bufs[2], &bufs[3], ucc, u_chan, v_chan);
        }
    }

    // Release the buffer objects
    PyBuffer_Release(&bufs[0]);
    PyBuffer_Release(&bufs[1]);
    PyBuffer_Release(&bufs[2]);
    PyBuffer_Release(&bufs[3]);

    return Py_BuildValue("");  // return Py_None while incrementing it
}


/* A list of all the methods defined by this module.
"METH_VARGS" tells Python how to call the handler.
The {NULL, NULL} entry indicates the end of the method definitions.*/
static PyMethodDef dds_defs_ext_methods[] = {
    {"unpack_dxt1", py_unpack_dxt1, METH_VARARGS, ""},
    {"unpack_dxt2_3", py_unpack_dxt2_3, METH_VARARGS, ""},
    {"unpack_dxt4_5", py_unpack_dxt4_5, METH_VARARGS, ""},
    {"unpack_dxt5a", py_unpack_dxt5a, METH_VARARGS, ""},
    {"unpack_dxn", py_unpack_dxn, METH_VARARGS, ""},
    {"unpack_ctx1", py_unpack_ctx1, METH_VARARGS, ""},
    {"unpack_vu", py_unpack_vu, METH_VARARGS, ""},

    {"pack_dxt1", py_pack_dxt1, METH_VARARGS, ""},
    {"pack_dxt2_3", py_pack_dxt2_3, METH_VARARGS, ""},
    {"pack_dxt4_5", py_pack_dxt4_5, METH_VARARGS, ""},
    //{"pack_dxt5a", py_pack_dxt5a, METH_VARARGS, ""},
    //{"pack_dxn", py_pack_dxn, METH_VARARGS, ""},
    //{"pack_ctx1", py_pack_ctx1, METH_VARARGS, ""},
    {"pack_vu", py_pack_vu, METH_VARARGS, ""},
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
