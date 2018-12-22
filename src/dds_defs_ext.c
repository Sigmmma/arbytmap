#pragma warning(disable: 4101)
#include <math.h>
#include "shared.h"

// to prevent loss of precision, cast these to ULL's before squaring
#define SQ(x) ((uint64)x)*((uint64)x)

#define DEFINE_UNPACK_VARIABLES(packed_type, src_unpacked_type, dst_unpacked_type)\
    packed_type *packed_pix     = (packed_type *) packed_pix_buf->buf;\
    dst_unpacked_type *unpacked_pix = (dst_unpacked_type *) unpacked_pix_buf->buf;\
    dst_unpacked_type *scales[4]={(dst_unpacked_type *)a_scale_buf->buf,\
                                  (dst_unpacked_type *)r_scale_buf->buf,\
                                  (dst_unpacked_type *)g_scale_buf->buf,\
                                  (dst_unpacked_type *)b_scale_buf->buf};\
    dst_unpacked_type *scale = (dst_unpacked_type *)unpacked_pix_buf->buf;\
    src_unpacked_type src_unpacked_max = ((1 << (8 * sizeof(src_unpacked_type))) - 1);\
    dst_unpacked_type dst_unpacked_max = ((1 << (8 * sizeof(dst_unpacked_type))) - 1);\
    int src_chan, dst_chan;\
    uint64 i, j, pxl_i, max_i;

#define DEFINE_DXT_UNPACK_VARIABLES()\
    int chans_per_tex=ucc*pix_per_tex, txl_pxl_i;\
    uint64 idx;

#define DEFINE_DXT_COLOR_UNPACK_VARIABLES()\
    uint16 color0, color1;\
    uint32 color_idx;\
    uint8 c_0[4]={255,0,0,0}, c_1[4]={255,0,0,0},\
          c_2[4]={255,0,0,0}, c_3[4]={255,0,0,0},\
          transparent[4]={0,0,0,0};\
    uint8 *color, *colors[4] = {c_0, c_1, c_2, c_3};

#define READ_DXT_COLORS(i)\
    color0 = packed_pix[(i)]&0xFFff;\
    color1 = packed_pix[(i)]>>16;\
    color_idx = packed_pix[(i)+1];\
    \
    c_0[1] = (((color0>>11) & 31)*255 + 15)/31;\
    c_0[2] = (((color0>>5)  & 63)*255 + 31)/63;\
    c_0[3] = (( color0      & 31)*255 + 15)/31;\
    \
    c_1[1] = (((color1>>11) & 31)*255 + 15)/31;\
    c_1[2] = (((color1>>5)  & 63)*255 + 31)/63;\
    c_1[3] = (( color1      & 31)*255 + 15)/31;\
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

#define UNPACK_DXT_COLORS(dst_chan)\
    for (dst_chan = 0; dst_chan < ucc; dst_chan++) {\
        src_chan = chan_map[dst_chan];\
        if (src_chan < 0 && dst_chan == 0)\
            unpacked_pix[pxl_i] = dst_unpacked_max;\
        else if (src_chan >= 0)\
            unpacked_pix[pxl_i] = scales[dst_chan][color[src_chan]];\
        pxl_i++;\
    }

#define UNPACK_DXT2345_COLORS(dst_chan)\
    for (dst_chan = 0; dst_chan < ucc; dst_chan++) {\
        if (chan_map[dst_chan] < 0 && dst_chan == 0)\
            unpacked_pix[pxl_i] = dst_unpacked_max;\
        else if (chan_map[dst_chan] > 0)\
            unpacked_pix[pxl_i] = scales[dst_chan][color[chan_map[dst_chan]]];\
        else if (chan_map[dst_chan] == 0)\
            unpacked_pix[pxl_i] = scales[dst_chan][a];\
        pxl_i++;\
    }

#define READ_DXT5_A(lookup, idx, i, val0, val1)\
    lookup[0] = val0 = packed_pix[i]&0xFF;\
    lookup[1] = val1 = (packed_pix[i]>>8)&0xFF;\
    idx = (((uint64)packed_pix[(i)+1])<<16) | (packed_pix[(i)]>>16);\
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

#define _CALC_B_NORMALIZE(r, g, b, x, y, mask, sign_mask, typ)\
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

#define CALC_B_NORMALIZE(r, g, b, x, y, unpacked_max, typ)\
_CALC_B_NORMALIZE(r, g, b, x, y, (unpacked_max >> 1), ((unpacked_max + 1) >> 1), typ)

#define _CALC_W_NORMALIZE(u, v, w, max_len, mask, sign_mask, typ)\
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

#define CALC_W_NORMALIZE(u, v, w, unpacked_max, typ)\
_CALC_W_NORMALIZE(u, v, w, (unpacked_max>> 1), unpacked_max, ((unpacked_max + 1) >> 1), typ)

#define PICK_DXT_PALETTE_DIST(i, j)\
    for (i = 0; i < chans_per_tex; i+=ucc) {\
        r = unpacked[r_pxl_i + i];\
        g = unpacked[g_pxl_i + i];\
        b = unpacked[b_pxl_i + i];\
        for (j = i + ucc; j < chans_per_tex; j+=ucc) {\
            /* if this will be solid black we wont care about its color*/\
            if (make_alpha && (unpacked[pxl_i + j] < a_cutoff)) continue;\
            dist1 = (SQ(r - unpacked[r_pxl_i + j])+\
                     SQ(g - unpacked[g_pxl_i + j])+\
                     SQ(b - unpacked[b_pxl_i + j]));\
            if (dist1 >= dist0) {\
                dist0 = dist1;\
                c_0i = i;\
                c_1i = j;\
            }\
        }\
    }

#define PACK_DXT_COLOR(c, i, packed_c)\
    c[1] = unpacked[r_pxl_i + i];\
    c[2] = unpacked[g_pxl_i + i];\
    c[3] = unpacked[b_pxl_i + i];\
    packed_c = ((((r_scale[c[1]] * 31 + 15UL) / 255)<<11) |\
                (((g_scale[c[2]] * 63 + 31UL) / 255)<<5) |\
                 ((b_scale[c[3]] * 31 + 15UL) / 255)) & 0xFFff;

#define SWAP_COLORS()\
    tmp[1] = c_0[1]; tmp[2] = c_0[2]; tmp[3] = c_0[3];\
    c_0[1] = c_1[1]; c_0[2] = c_1[2]; c_0[3] = c_1[3];\
    c_1[1] = tmp[1]; c_1[2] = tmp[2]; c_1[3] = tmp[3];\
    tmp_color = color0;\
    color0    = color1;\
    color1    = tmp_color;

#define CALC_DXT_IDX_NO_ALPHA(idx, i)\
    c_2[1] = (c_0[1]*2 + c_1[1])/3; c_3[1] = (c_0[1] + c_1[1]*2)/3;\
    c_2[2] = (c_0[2]*2 + c_1[2])/3; c_3[2] = (c_0[2] + c_1[2]*2)/3;\
    c_2[3] = (c_0[3]*2 + c_1[3])/3; c_3[3] = (c_0[3] + c_1[3]*2)/3;\
    for (i=0; i<chans_per_tex; i+=ucc) {\
        r = unpacked[r_pxl_i+i];\
        g = unpacked[g_pxl_i+i];\
        b = unpacked[b_pxl_i+i];\
        dist0 = SQ(r-c_0[1]) + SQ(g-c_0[2]) + SQ(b-c_0[3]);\
        dist2 = SQ(r-c_2[1]) + SQ(g-c_2[2]) + SQ(b-c_2[3]);\
        dist3 = SQ(r-c_3[1]) + SQ(g-c_3[2]) + SQ(b-c_3[3]);\
        dist1 = SQ(r-c_1[1]) + SQ(g-c_1[2]) + SQ(b-c_1[3]);\
        if (dist1 < dist0) {\
            if (dist1 < dist3) {\
                idx |= 1<<(i>>1);\
            } else {\
                idx |= 3<<(i>>1);\
            }\
        } else if (dist2 < dist0) {\
            idx |= 2<<(i>>1);\
        }\
    }

#define CALC_DXT_IDX_ALPHA(idx, i)\
    c_2[1] = (c_0[1] + c_1[1])/2;\
    c_2[2] = (c_0[2] + c_1[2])/2;\
    c_2[3] = (c_0[3] + c_1[3])/2;\
    for (i=0; i<chans_per_tex; i+=ucc) {\
        if (unpacked[pxl_i + i] < a_cutoff) {\
            idx |= 3<<(i>>1);\
            continue;\
        }\
        r = unpacked[r_pxl_i+i];\
        g = unpacked[g_pxl_i+i];\
        b = unpacked[b_pxl_i+i];\
        dist0 = SQ(r-c_0[1]) + SQ(g-c_0[2]) + SQ(b-c_0[3]);\
        dist2 = SQ(r-c_2[1]) + SQ(g-c_2[2]) + SQ(b-c_2[3]);\
        dist1 = SQ(r-c_1[1]) + SQ(g-c_1[2]) + SQ(b-c_1[3]);\
        if ((dist1 < dist0) && (dist1 < dist2)) {\
            idx |= 1<<(i>>1);\
        } else if (dist2 < dist0) {\
            idx |= 2<<(i>>1);\
        }\
    }

#define PICK_DXT5_ALPHA_DIST(a0, a1, a_pxl_i, a_tmp, a_scale, i)\
    a0 = a1 = a_scale[unpacked[a_pxl_i]];\
    for (i = ucc; i < chans_per_tex; i += ucc) {\
        a_tmp = a_scale[unpacked[a_pxl_i + i]];\
        if (a_tmp > a0) a0 = a_tmp;\
        if (a_tmp < a1) a1 = a_tmp;\
    }

/*there are 4 interpolated colors in PICK_DXT5_ALPHA_IDX_0_255 mode
0 = a0               1 = a1
2 = (6*a0 +   a1)/7  3 = (5*a0 + 2*a1)/7
4 = (4*a0 + 3*a1)/7  5 = (3*a0 + 4*a1)/7
6 = (2*a0 + 5*a1)/7  7 = (  a0 + 6*a1)/7*/
#define PICK_DXT5_ALPHA_IDX(a0, a1, a_idx, a_pxl_i, a_tmp, a_dif, a_scale, i)\
    a_dif = a0 - a1;\
    for (i = 0; i<chans_per_tex; i += ucc) {\
        a_tmp = ((a_scale[unpacked[a_pxl_i + i]] - a1) * 7 + (a_dif >> 1)) / a_dif;\
        if (a_tmp == 0) {\
            a_idx |= 1ULL << ((i * 3) / ucc);\
        } else if (a_tmp < 7) {\
            a_idx |= ((uint64)(8 - a_tmp)) << ((i * 3) / ucc);\
        }\
    }


/*there are 4 interpolated colors in PICK_DXT5_ALPHA_IDX_0_255 mode
0 =  a0              1 = a1
2 = (4*a0 +   a1)/5  3 = (3*a0 + 2*a1)/5
4 = (2*a0 + 3*a1)/5  5 = (  a0 + 4*a1)/5
6 =  0               7 = 255*/
#define PICK_DXT5_ALPHA_IDX_0_255(a0, a1, a_idx, a_pxl_i, a_tmp, a_dif, a_scale, i)\
    a0 = 255; a1 = 0;\
    for (i = 0; i<chans_per_tex; i += ucc) {\
        a_tmp = a_scale[unpacked[a_pxl_i + i]];\
        if ((a_tmp < a0) && (a_tmp != 0))   a0 = a_tmp;\
        if ((a_tmp > a1) && (a_tmp != 255)) a1 = a_tmp;\
    }\
    a_dif = a1 - a0;\
    if (!a1) {\
        a0 = a_dif = 0;\
        a1 = 255;\
    }\
    for (i = 0; i<chans_per_tex; i += ucc) {\
        a_tmp = a_scale[unpacked[a_pxl_i + i]];\
        if (a_tmp == 0) {\
            a_idx |= 6ULL << ((i * 3) / ucc);\
        } else if (a_tmp == 255) {\
            a_idx |= 7ULL << ((i * 3) / ucc);\
        } else if (a_dif) {\
            a_tmp = ((a_tmp - a0) * 5 + (a_dif >> 1)) / a_dif;\
            if (a_tmp == 5) {\
                a_idx |= 1ULL << ((i * 3) / ucc);\
            } else if (a_tmp > 0) {\
                a_idx |= ((uint64)(a_tmp + 1)) << ((i * 3) / ucc);\
            }\
        }\
    }

static void unpack_dxt1_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 ucc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint8, uint8)
    DEFINE_DXT_UNPACK_VARIABLES()
    DEFINE_DXT_COLOR_UNPACK_VARIABLES()

    max_i = unpacked_pix_buf->len / chans_per_tex;

    // loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i * chans_per_tex;
        j = i << 1;

        /*If the format DXT1 then the two entries in the array
        are the colors and the color indexing in that order.
        Also, if the first color is a larger integer
        then color key transparency is NOT used.*/
        READ_DXT_COLORS(j);

        for (txl_pxl_i = 0; txl_pxl_i < pix_per_tex; txl_pxl_i++) {
            color = colors[(color_idx >> (txl_pxl_i << 1)) & 3];
            UNPACK_DXT_COLORS(dst_chan);
        }
    }
}


static void unpack_dxt2_3_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 ucc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint8, uint8)
    DEFINE_DXT_UNPACK_VARIABLES()
    DEFINE_DXT_COLOR_UNPACK_VARIABLES()

    uint64 alpha;
    uint8 a;

    max_i = unpacked_pix_buf->len / chans_per_tex;

    // loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i * chans_per_tex;
        j = i << 2;

        alpha = ((uint64)packed_pix[j + 1] << 32) + packed_pix[j];
        READ_DXT_COLORS(j + 2);

        for (txl_pxl_i = 0; txl_pxl_i < pix_per_tex; txl_pxl_i++) {
            color = colors[(color_idx >> (txl_pxl_i << 1)) & 3];
            a = (((alpha >> (txl_pxl_i << 2)) & 15) * 255) / 15;
            UNPACK_DXT2345_COLORS(dst_chan);
        }
    }
}


static void unpack_dxt4_5_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 ucc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint8, uint8)
    DEFINE_DXT_UNPACK_VARIABLES()
    DEFINE_DXT_COLOR_UNPACK_VARIABLES()

    uint64 alpha_idx;
    uint8 a, alpha0, alpha1, a_lookup[8];

    max_i = unpacked_pix_buf->len / chans_per_tex;

    // loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i * chans_per_tex;
        j = i << 2;

        READ_DXT5_A(a_lookup, alpha_idx, j, alpha0, alpha1);
        READ_DXT_COLORS(j + 2);
        for (txl_pxl_i = 0; txl_pxl_i < pix_per_tex; txl_pxl_i++) {
            color = colors[(color_idx >> (txl_pxl_i << 1)) & 3];
            a = a_lookup[(alpha_idx >> (3 * txl_pxl_i)) & 7];
            UNPACK_DXT2345_COLORS(dst_chan);
        }
    }
}

static void unpack_dxt3a_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 ucc, sint8 scc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint8, uint8)
    DEFINE_DXT_UNPACK_VARIABLES()

    uint64 alpha;

    max_i = unpacked_pix_buf->len / chans_per_tex;

    //loop through each destination channel
    for (dst_chan = 0; dst_chan < ucc; dst_chan++) {
        scale = scales[dst_chan];
        src_chan = chan_map[dst_chan];
        pxl_i = dst_chan;

        if (src_chan < 0) {
            // not reading anything for this destination channel.
            // either leave it full black, or set it to full white.
            if (dst_chan == 0)
                // set alpha to full white for the entire image
                for (i = 0; i < (uint64)(unpacked_pix_buf->len); i += ucc)
                    unpacked_pix[i] = dst_unpacked_max;
            continue;
        }

        // loop through each texel
        for (i=0; i < max_i; i++) {
            j = (i * scc + src_chan) << 1;
            alpha = ((uint64)packed_pix[j + 1] << 32) + packed_pix[j];

            for (txl_pxl_i = 0; txl_pxl_i<pix_per_tex; txl_pxl_i++) {
                unpacked_pix[pxl_i] = scale[((alpha & 0xf) * 255) / 15];
                alpha = alpha >> 4;
                pxl_i += ucc;
            }
        }
    }
}

static void unpack_dxt5a_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 ucc, sint8 scc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint8, uint8)
    DEFINE_DXT_UNPACK_VARIABLES()
    uint8 val0, val1, lookup[8];

    max_i = unpacked_pix_buf->len / chans_per_tex;

    //loop through each destination channel
    for (dst_chan = 0; dst_chan < ucc; dst_chan++) {
        scale = scales[dst_chan];
        src_chan = chan_map[dst_chan];
        pxl_i = dst_chan;

        if (src_chan < 0) {
            // not reading anything for this destination channel.
            // either leave it full black, or set it to full white.
            if (dst_chan == 0)
                // set alpha to full white for the entire image
                for (i = 0; i < (uint64)(unpacked_pix_buf->len); i += ucc)
                    unpacked_pix[i] = dst_unpacked_max;
            continue;
        }

        // loop through each texel
        for (i=0; i < max_i; i++) {
            j = (i * scc + src_chan) << 1;
            READ_DXT5_A(lookup, idx, j, val0, val1);
            for (txl_pxl_i = 0; txl_pxl_i<pix_per_tex; txl_pxl_i++) {
                unpacked_pix[pxl_i] = scale[lookup[(idx >> (3 * txl_pxl_i)) & 7]];
                pxl_i += ucc;
            }
        }
    }
}


static void unpack_dxn_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 ucc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint8, uint8)
    DEFINE_DXT_UNPACK_VARIABLES()

    uint8 r0, r1, r_lookup[8], g0, g1, g_lookup[8];
    uint8 x, y, r, g, b, color[4]={0,0,0,0};
    uint64 r_idx, g_idx;
    double d, n_len;

    max_i = unpacked_pix_buf->len / chans_per_tex;
    pxl_i = 0;

    // loop through each texel
    for (i=0; i < max_i; i++) {
        j = i << 2;

        READ_DXT5_A(g_lookup, g_idx, j,     g0, g1);
        READ_DXT5_A(r_lookup, r_idx, j + 2, r0, r1);

        for (txl_pxl_i = 0; txl_pxl_i<pix_per_tex; txl_pxl_i++) {
            g = y = g_lookup[(g_idx >> (txl_pxl_i * 3)) & 7];
            r = x = r_lookup[(r_idx >> (txl_pxl_i * 3)) & 7];
            CALC_B_NORMALIZE(r, g, b, x, y, src_unpacked_max, uint8);
            color[1] = r; color[2] = g; color[3] = b;
            //loop through each destination channel
            for (dst_chan = 0; dst_chan < ucc; dst_chan++) {
                if (chan_map[dst_chan] >= 0)
                    unpacked_pix[pxl_i] = scales[dst_chan][color[chan_map[dst_chan]]];
                else if (dst_chan == 0)
                    unpacked_pix[pxl_i] = dst_unpacked_max;
                pxl_i++;
            }
        }
    }
}


static void unpack_ctx1_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 ucc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint8, uint8)
    DEFINE_DXT_UNPACK_VARIABLES()
    DEFINE_DXT_COLOR_UNPACK_VARIABLES()
    uint8 x, y, a = 0;  // set a to 0 as the scale should only have 1 value in it
    double d, n_len;

    max_i = unpacked_pix_buf->len / chans_per_tex;

    // loop through each texel
    for (i = 0; i < max_i; i++) {
        pxl_i = i * chans_per_tex;
        j = i << 1;

        c_0[1] = packed_pix[j] & 255;
        c_0[2] = (packed_pix[j] >> 8) & 255;
        CALC_B_NORMALIZE(c_0[1], c_0[2], c_0[3], x, y, src_unpacked_max, uint8);

        c_1[1] = (packed_pix[j] >> 16) & 255;
        c_1[2] = (packed_pix[j] >> 24) & 255;
        CALC_B_NORMALIZE(c_1[1], c_1[2], c_1[3], x, y, src_unpacked_max, uint8);

        color_idx = packed_pix[j + 1];

        c_2[1] = (c_0[1] * 2 + c_1[1]) / 3;
        c_2[2] = (c_0[2] * 2 + c_1[2]) / 3;
        c_2[3] = (c_0[3] * 2 + c_1[3]) / 3;

        c_3[1] = (c_0[1] + 2 * c_1[1]) / 3;
        c_3[2] = (c_0[2] + 2 * c_1[2]) / 3;
        c_3[3] = (c_0[3] + 2 * c_1[3]) / 3;
        for (txl_pxl_i = 0; txl_pxl_i<pix_per_tex; txl_pxl_i++) {
            color = colors[(color_idx >> (txl_pxl_i << 1)) & 3];
            UNPACK_DXT2345_COLORS(dst_chan);
        }
    }
}


static void unpack_v16u16_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 ucc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint16, uint8)
    uint16 color[4] = { 0,0,0,0 };
    sint16 u, v, w;
    double d, n_len;

    max_i = unpacked_pix_buf->len / ucc;
    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i*ucc;

        u = packed_pix[i] & 0xFFff;
        v = packed_pix[i] >> 16;
        CALC_W_NORMALIZE(u, v, w, src_unpacked_max, sint32);
        color[1] = u + ((src_unpacked_max + 1) >> 1);
        color[2] = v + ((src_unpacked_max + 1) >> 1);
        color[3] = w + ((src_unpacked_max + 1) >> 1);

        //loop through each destination channel
        for (dst_chan = 0; dst_chan < ucc; dst_chan++) {
            if (chan_map[dst_chan] >= 0)
                unpacked_pix[pxl_i] = scales[dst_chan][color[chan_map[dst_chan]]];
            else if (dst_chan == 0)
                unpacked_pix[pxl_i] = dst_unpacked_max;
            pxl_i++;
        }
    }
}


static void unpack_v8u8_8(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 ucc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint16, uint8, uint8)
    uint8 color[4] = { 0,0,0,0 };
    sint16 u, v, w;
    double d, n_len;

    max_i = unpacked_pix_buf->len / ucc;

    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i * ucc;

        u = packed_pix[i] & 0xFF;
        v = packed_pix[i] >> 8;
        CALC_W_NORMALIZE(u, v, w, src_unpacked_max, sint16);
        color[1] = u + ((src_unpacked_max + 1) >> 1);
        color[2] = v + ((src_unpacked_max + 1) >> 1);
        color[3] = w + ((src_unpacked_max + 1) >> 1);

        //loop through each destination channel
        for (dst_chan = 0; dst_chan < ucc; dst_chan++) {
            if (chan_map[dst_chan] >= 0)
                unpacked_pix[pxl_i] = scales[dst_chan][color[chan_map[dst_chan]]];
            else if (dst_chan == 0)
                unpacked_pix[pxl_i] = dst_unpacked_max;
            pxl_i++;
        }
    }
}


static void pack_dxt1_8(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
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

    packed   = (uint32 *)packed_pix_buf->buf;
    unpacked = (uint8 *)unpacked_pix_buf->buf;
    r_scale  = (uint8 *)r_scale_buf->buf;
    g_scale  = (uint8 *)g_scale_buf->buf;
    b_scale  = (uint8 *)b_scale_buf->buf;

    max_i = packed_pix_buf->len / 8; // 8 bytes per texel

    // loop through each texel
    for (i=0; i < max_i; i++) {
        c_0i = c_1i = idx = dist0 = 0;
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
        PICK_DXT_PALETTE_DIST(j, k);

        // store furthest apart colors for use
        PACK_DXT_COLOR(c_0, c_0i, color0);
        PACK_DXT_COLOR(c_1, c_1i, color1);

        if (make_alpha || (color0 != color1)) {
            //  make_alpha == no_alpha_for_texel
            if (make_alpha == (color0 > color1)) {
                SWAP_COLORS()
            }

            if (color0 > color1) {
                CALC_DXT_IDX_NO_ALPHA(idx, j)
            } else {
                CALC_DXT_IDX_ALPHA(idx, j)
            }
        }
        packed[i << 1] = (color1 << 16) | color0;
        packed[(i << 1) + 1] = idx;
    }
}


static void pack_dxt2_3_8(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
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

    packed = (uint32 *) packed_pix_buf->buf;
    unpacked = (uint8 *)unpacked_pix_buf->buf;
    a_scale = (uint8 *)a_scale_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;
    b_scale = (uint8 *)b_scale_buf->buf;

    max_i = packed_pix_buf->len / 16; // 16 bytes per texel

    // loop through each texel
    for (i=0; i < max_i; i++) {
        alpha = c_0i = c_1i = idx = dist0 = 0;
        pxl_i = i*chans_per_tex;
        r_pxl_i = pxl_i+1;
        g_pxl_i = pxl_i+2;
        b_pxl_i = pxl_i+3;

        // calculate the alpha
        for (j=0; j<chans_per_tex; j+=ucc) {
            alpha |= ((uint64)((a_scale[unpacked[pxl_i + j]] * 15 + 7) / 255)) << j;
        }
        // compare distance between all pixels and find the two furthest apart
        // (we are actually comparing the area of the distance as it's faster)
        PICK_DXT_PALETTE_DIST(j, k);

        // store furthest apart colors for use
        PACK_DXT_COLOR(c_0, c_0i, color0);
        PACK_DXT_COLOR(c_1, c_1i, color1);

        if (color0 != color1) {
            if (color0 < color1) {
                SWAP_COLORS()
            }
            CALC_DXT_IDX_NO_ALPHA(idx, j)
        }
        packed[i << 2]       = alpha & 0xFFffFFff;
        packed[(i << 2) + 1] = alpha >> 32;
        packed[(i << 2) + 2] = (color1 << 16) | color0;
        packed[(i << 2) + 3] = idx;
    }
}


static void pack_dxt4_5_8(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
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

    packed = (uint32 *) packed_pix_buf->buf;
    unpacked = (uint8 *)unpacked_pix_buf->buf;
    a_scale = (uint8 *)a_scale_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;
    b_scale = (uint8 *)b_scale_buf->buf;

    max_i = packed_pix_buf->len / 16; // 16 bytes per texel

    // loop through each texel
    for (i=0; i < max_i; i++) {
        c_0i = c_1i = a_idx = idx = dist0 = 0;
        pxl_i = i*chans_per_tex;
        r_pxl_i = pxl_i+1;
        g_pxl_i = pxl_i+2;
        b_pxl_i = pxl_i+3;

        // calculate the min and max alpha values
        PICK_DXT5_ALPHA_DIST(a0, a1, pxl_i, a_tmp, a_scale, j)

        // calculate the alpha indexing
        if (a0 == a1) {
            ; // values are the same; do nothing
        } else if (a1 && (a0 != 255)) {
            PICK_DXT5_ALPHA_IDX(a0, a1, a_idx, pxl_i, a_tmp, a_dif, a_scale, j)
        } else {
            // max is 255 or min is 0
            PICK_DXT5_ALPHA_IDX_0_255(a0, a1, a_idx, pxl_i, a_tmp, a_dif, a_scale, j)
        }

        // compare distance between all pixels and find the two furthest apart
        // (we are actually comparing the area of the distance as it's faster)
        PICK_DXT_PALETTE_DIST(j, k);

        // store furthest apart colors for use
        PACK_DXT_COLOR(c_0, c_0i, color0);
        PACK_DXT_COLOR(c_1, c_1i, color1);

        if (color0 != color1) {
            if (color0 < color1) {
                SWAP_COLORS()
            }
            CALC_DXT_IDX_NO_ALPHA(idx, j)
        }
        packed[i << 2] = ((a_idx & 0xFFff) << 16) | (a1 << 8) | a0;
        packed[(i << 2) + 1] = (a_idx >> 16) & 0xFFffFFff;
        packed[(i << 2) + 2] = (color1 << 16) | color0;
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

    max_i = packed_pix_buf->len / 2;  // 2 bytes per pixel
    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i*ucc;
        packed_pix[i] = ((((uint32)v_scale[unpacked_pix[pxl_i + v_chan]]) << 8) |
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

    max_i = packed_pix_buf->len / 4;  // 4 bytes per pixel
    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i*ucc;
        packed_pix[i] = ((((uint32)v_scale[unpacked_pix[pxl_i + v_chan]]) << 16) |
                                   u_scale[unpacked_pix[pxl_i + u_chan]]) ^ 0x80008000UL;
    }
}


/*
    Deep color versions of the above functions
*/


static void unpack_dxt1_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 ucc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint8, uint16)
    DEFINE_DXT_UNPACK_VARIABLES()
    DEFINE_DXT_COLOR_UNPACK_VARIABLES()

    max_i = unpacked_pix_buf->len / (chans_per_tex * 2);

    // loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i * chans_per_tex;
        j = i << 1;

        /*If the format DXT1 then the two entries in the array
        are the colors and the color indexing in that order.
        Also, if the first color is a larger integer
        then color key transparency is NOT used.*/
        READ_DXT_COLORS(j);

        for (txl_pxl_i = 0; txl_pxl_i < pix_per_tex; txl_pxl_i++) {
            color = colors[(color_idx >> (txl_pxl_i << 1)) & 3];
            UNPACK_DXT_COLORS(dst_chan);
        }
    }
}

static void unpack_dxt2_3_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 ucc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint8, uint16)
    DEFINE_DXT_UNPACK_VARIABLES()
    DEFINE_DXT_COLOR_UNPACK_VARIABLES()

    uint64 alpha;
    uint8 a;

    max_i = unpacked_pix_buf->len / (chans_per_tex * 2);

    // loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i * chans_per_tex;
        j = i << 2;

        alpha = ((uint64)packed_pix[j + 1] << 32) + packed_pix[j];
        READ_DXT_COLORS(j + 2);

        for (txl_pxl_i = 0; txl_pxl_i < pix_per_tex; txl_pxl_i++) {
            a = (((alpha >> (txl_pxl_i << 2)) & 15) * 255) / 15;
            color = colors[(color_idx >> (txl_pxl_i << 1)) & 3];
            UNPACK_DXT2345_COLORS(dst_chan);
        }
    }
}


static void unpack_dxt4_5_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 ucc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint8, uint16)
    DEFINE_DXT_UNPACK_VARIABLES()
    DEFINE_DXT_COLOR_UNPACK_VARIABLES()

    uint64 alpha_idx;
    uint8 a, alpha0, alpha1, a_lookup[8];

    max_i = unpacked_pix_buf->len / (chans_per_tex * 2);

    // loop through each texel
    for (i=0; i < max_i; i++) {
        pxl_i = i * chans_per_tex;
        j = i << 2;

        READ_DXT5_A(a_lookup, alpha_idx, j, alpha0, alpha1);
        READ_DXT_COLORS(j + 2);
        for (txl_pxl_i = 0; txl_pxl_i < pix_per_tex; txl_pxl_i++) {
            a = a_lookup[(alpha_idx >> (3 * txl_pxl_i)) & 7];
            color = colors[(color_idx >> (txl_pxl_i << 1)) & 3];
            UNPACK_DXT2345_COLORS(dst_chan);
        }
    }
}

static void unpack_dxt3a_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 ucc, sint8 scc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint8, uint16)
    DEFINE_DXT_UNPACK_VARIABLES()

    uint64 alpha;

    max_i = unpacked_pix_buf->len / (chans_per_tex * 2);

    //loop through each destination channel
    for (dst_chan = 0; dst_chan < ucc; dst_chan++) {
        scale = scales[dst_chan];
        src_chan = chan_map[dst_chan];
        pxl_i = dst_chan;

        if (src_chan < 0) {
            // not reading anything for this destination channel.
            // either leave it full black, or set it to full white.
            if (dst_chan == 0)
                // set alpha to full white for the entire image
                for (i = 0; i < (uint64)(unpacked_pix_buf->len / 2); i += ucc)
                    unpacked_pix[i] = dst_unpacked_max;
            continue;
        }

        // loop through each texel
        for (i=0; i < max_i; i++) {
            j = (i * scc + src_chan) << 1;
            alpha = ((uint64)packed_pix[j + 1] << 32) + packed_pix[j];

            for (txl_pxl_i = 0; txl_pxl_i<pix_per_tex; txl_pxl_i++) {
                unpacked_pix[pxl_i] = scale[((alpha & 0xf) * 255) / 15];
                alpha = alpha >> 4;
                pxl_i += ucc;
            }
        }
    }
}


static void unpack_dxt5a_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 ucc, sint8 scc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint8, uint16)
    DEFINE_DXT_UNPACK_VARIABLES()
    uint8 val0, val1, lookup[8];

    max_i = unpacked_pix_buf->len / (chans_per_tex * 2);

    //loop through each destination channel
    for (dst_chan = 0; dst_chan < ucc; dst_chan++) {
        scale = scales[dst_chan];
        src_chan = chan_map[dst_chan];
        pxl_i = dst_chan;

        if (src_chan < 0) {
            // not reading anything for this destination channel.
            // either leave it full black, or set it to full white.
            if (dst_chan == 0)
                // set alpha to full white for the entire image
                for (i = 0; i < (uint64)(unpacked_pix_buf->len / 2); i += ucc)
                    unpacked_pix[i] = dst_unpacked_max;
            continue;
        }

        // loop through each texel
        for (i=0; i < max_i; i++) {
            j = (i * scc + src_chan) << 1;
            READ_DXT5_A(lookup, idx, j, val0, val1);
            for (txl_pxl_i = 0; txl_pxl_i<pix_per_tex; txl_pxl_i++) {
                unpacked_pix[pxl_i] = scale[lookup[(idx >> (3 * txl_pxl_i)) & 7]];
                pxl_i += ucc;
            }
        }
    }
}


static void unpack_dxn_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 ucc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint8, uint16)
    DEFINE_DXT_UNPACK_VARIABLES()

    uint8 r0, r1, r_lookup[8], g0, g1, g_lookup[8];
    uint8 x, y, r, g, b, color[4]={0,0,0,0};
    uint64 r_idx, g_idx;
    double d, n_len;

    max_i = unpacked_pix_buf->len / (chans_per_tex * 2);
    pxl_i = 0;

    // loop through each texel
    for (i=0; i < max_i; i++) {
        j = i << 2;

        READ_DXT5_A(g_lookup, g_idx, j,     g0, g1);
        READ_DXT5_A(r_lookup, r_idx, j + 2, r0, r1);

        for (txl_pxl_i = 0; txl_pxl_i<pix_per_tex; txl_pxl_i++) {
            g = y = g_lookup[(g_idx >> (txl_pxl_i * 3)) & 7];
            r = x = r_lookup[(r_idx >> (txl_pxl_i * 3)) & 7];
            CALC_B_NORMALIZE(r, g, b, x, y, src_unpacked_max, uint8);
            color[1] = r; color[2] = g; color[3] = b;
            //loop through each destination channel
            for (dst_chan = 0; dst_chan < ucc; dst_chan++) {
                if (chan_map[dst_chan] >= 0)
                    unpacked_pix[pxl_i] = scales[dst_chan][color[chan_map[dst_chan]]];
                else if (dst_chan == 0)
                    unpacked_pix[pxl_i] = dst_unpacked_max;
                pxl_i++;
            }
        }
    }
}


static void unpack_ctx1_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 pix_per_tex, sint8 ucc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint8, uint16)
    DEFINE_DXT_UNPACK_VARIABLES()
    DEFINE_DXT_COLOR_UNPACK_VARIABLES()
    uint8 x, y, a = 0;  // set a to 0 as the scale should only have 1 value in it
    double d, n_len;

    max_i = unpacked_pix_buf->len / (chans_per_tex * 2);

    // loop through each texel
    for (i = 0; i < max_i; i++) {
        pxl_i = i*chans_per_tex;
        j = i << 1;

        c_0[1] = packed_pix[j] & 255;
        c_0[2] = (packed_pix[j] >> 8) & 255;
        CALC_B_NORMALIZE(c_0[1], c_0[2], c_0[3], x, y, src_unpacked_max, uint8);

        c_1[1] = (packed_pix[j] >> 16) & 255;
        c_1[2] = (packed_pix[j] >> 24) & 255;
        CALC_B_NORMALIZE(c_1[1], c_1[2], c_1[3], x, y, src_unpacked_max, uint8);

        color_idx = packed_pix[j + 1];

        c_2[1] = (c_0[1] * 2 + c_1[1]) / 3;
        c_2[2] = (c_0[2] * 2 + c_1[2]) / 3;
        c_2[3] = (c_0[3] * 2 + c_1[3]) / 3;

        c_3[1] = (c_0[1] + 2 * c_1[1]) / 3;
        c_3[2] = (c_0[2] + 2 * c_1[2]) / 3;
        c_3[3] = (c_0[3] + 2 * c_1[3]) / 3;
        for (txl_pxl_i = 0; txl_pxl_i<pix_per_tex; txl_pxl_i++) {
            color = colors[(color_idx >> (txl_pxl_i << 1)) & 3];
            UNPACK_DXT2345_COLORS(dst_chan);
        }
    }
}


static void unpack_v8u8_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 ucc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint16, uint8, uint16)
    uint8 color[4] = { 0,0,0,0 };
    sint16 u, v, w;
    double d, n_len;

    max_i = unpacked_pix_buf->len / (ucc * 2);
    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i*ucc;

        u = packed_pix[i] & src_unpacked_max;
        v = packed_pix[i] >> 8;
        CALC_W_NORMALIZE(u, v, w, src_unpacked_max, sint8);
        color[1] = u + ((src_unpacked_max + 1) >> 1);
        color[2] = v + ((src_unpacked_max + 1) >> 1);
        color[3] = w + ((src_unpacked_max + 1) >> 1);

        //loop through each destination channel
        for (dst_chan = 0; dst_chan < ucc; dst_chan++) {
            if (chan_map[dst_chan] >= 0)
                unpacked_pix[pxl_i] = scales[dst_chan][color[chan_map[dst_chan]]];
            else if (dst_chan == 0)
                unpacked_pix[pxl_i] = dst_unpacked_max;
            pxl_i++;
        }
    }
}


static void unpack_v16u16_16(
    Py_buffer *unpacked_pix_buf, Py_buffer *packed_pix_buf,
    Py_buffer *a_scale_buf, Py_buffer *r_scale_buf,
    Py_buffer *g_scale_buf, Py_buffer *b_scale_buf,
    sint8 ucc, sint8 *chan_map)
{
    DEFINE_UNPACK_VARIABLES(uint32, uint16, uint16)
    uint16 color[4] = { 0,0,0,0 };
    sint32 u, v, w;
    double d, n_len;

    max_i = unpacked_pix_buf->len / (ucc * 2);
    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i*ucc;

        u = packed_pix[i] & src_unpacked_max;
        v = packed_pix[i] >> 16;
        CALC_W_NORMALIZE(u, v, w, src_unpacked_max, sint16);
        color[1] = u + ((src_unpacked_max + 1) >> 1);
        color[2] = v + ((src_unpacked_max + 1) >> 1);
        color[3] = w + ((src_unpacked_max + 1) >> 1);

        //loop through each destination channel
        for (dst_chan = 0; dst_chan < ucc; dst_chan++) {
            if (chan_map[dst_chan] >= 0)
                unpacked_pix[pxl_i] = scales[dst_chan][color[chan_map[dst_chan]]];
            else if (dst_chan == 0)
                unpacked_pix[pxl_i] = dst_unpacked_max;
            pxl_i++;
        }
    }
}


static void pack_dxt1_16(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
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

    packed = (uint32 *) packed_pix_buf->buf;
    unpacked = (uint16 *)unpacked_pix_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;
    b_scale = (uint8 *)b_scale_buf->buf;

    max_i = packed_pix_buf->len / 8; // 8 bytes per texel

    // loop through each texel
    for (i=0; i < max_i; i++) {
        dist0 = c_0i = c_1i = idx = 0;
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
        PICK_DXT_PALETTE_DIST(j, k);

        // store furthest apart colors for use
        PACK_DXT_COLOR(c_0, c_0i, color0);
        PACK_DXT_COLOR(c_1, c_1i, color1);

        if (make_alpha || (color0 != color1)) {
            //  make_alpha == no_alpha_for_texel
            if (make_alpha == (color0 > color1)) {
                SWAP_COLORS()
            }

            if (color0 > color1) {
                CALC_DXT_IDX_NO_ALPHA(idx, j)
            } else {
                CALC_DXT_IDX_ALPHA(idx, j)
            }
        }
        packed[i << 1] = (color1 << 16) | color0;
        packed[(i << 1) + 1] = idx;
    }
}


static void pack_dxt2_3_16(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
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

    packed = (uint32 *) packed_pix_buf->buf;
    unpacked = (uint16 *)unpacked_pix_buf->buf;
    a_scale = (uint8 *)a_scale_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;
    b_scale = (uint8 *)b_scale_buf->buf;

    max_i = packed_pix_buf->len / 16; // 8 bytes per texel

    // loop through each texel
    for (i=0; i < max_i; i++) {
        dist0 = alpha = c_0i = c_1i = idx = 0;
        pxl_i = i*chans_per_tex;
        r_pxl_i = pxl_i+1;
        g_pxl_i = pxl_i+2;
        b_pxl_i = pxl_i+3;

        // calculate the alpha
        for (j=0; j<chans_per_tex; j+=ucc) {
            alpha |= ((uint64)((a_scale[unpacked[pxl_i+j]] * 15 + 7) / 255))<<j;
        }
        // compare distance between all pixels and find the two furthest apart
        // (we are actually comparing the area of the distance as it's faster)
        PICK_DXT_PALETTE_DIST(j, k);

        // store furthest apart colors for use
        PACK_DXT_COLOR(c_0, c_0i, color0);
        PACK_DXT_COLOR(c_1, c_1i, color1);

        if (color0 != color1) {
            if (color0 < color1) {
                SWAP_COLORS()
            }
            CALC_DXT_IDX_NO_ALPHA(idx, j)
        }
        packed[i << 2]       = alpha & 0xFFffFFff;
        packed[(i << 2) + 1] = alpha >> 32;
        packed[(i << 2) + 2] = (color1 << 16) | color0;
        packed[(i << 2) + 3] = idx;
    }
}

static void pack_dxt4_5_16(
    Py_buffer *packed_pix_buf, Py_buffer *unpacked_pix_buf,
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

    packed = (uint32 *) packed_pix_buf->buf;
    unpacked = (uint16 *)unpacked_pix_buf->buf;
    a_scale = (uint8 *)a_scale_buf->buf;
    r_scale = (uint8 *)r_scale_buf->buf;
    g_scale = (uint8 *)g_scale_buf->buf;
    b_scale = (uint8 *)b_scale_buf->buf;

    max_i = packed_pix_buf->len / 16; // 16 bytes per texel

    // loop through each texel
    for (i=0; i < max_i; i++) {
        dist0 = c_0i = c_1i = a_idx = idx = 0;
        pxl_i = i*chans_per_tex;
        r_pxl_i = pxl_i+1;
        g_pxl_i = pxl_i+2;
        b_pxl_i = pxl_i+3;

        // calculate the min and max alpha values
        PICK_DXT5_ALPHA_DIST(a0, a1, pxl_i, a_tmp, a_scale, j);

        // calculate the alpha indexing
        if (a0 == a1) {
            ; // values are the same; do nothing
        } else if (a1 && (a0 != 255)) {
            PICK_DXT5_ALPHA_IDX(a0, a1, a_idx, pxl_i, a_tmp, a_dif, a_scale, j);
        } else {
            // max is 255 or min is 0
            PICK_DXT5_ALPHA_IDX_0_255(a0, a1, a_idx, pxl_i, a_tmp, a_dif, a_scale, j);
        }

        // compare distance between all pixels and find the two furthest apart
        // (we are actually comparing the area of the distance as it's faster)
        PICK_DXT_PALETTE_DIST(j, k);

        // store furthest apart colors for use
        PACK_DXT_COLOR(c_0, c_0i, color0);
        PACK_DXT_COLOR(c_1, c_1i, color1);

        if (color0 != color1) {
            if (color0 < color1) {
                SWAP_COLORS()
            }
            CALC_DXT_IDX_NO_ALPHA(idx, j);
        }
        packed[i << 2] = ((a_idx & 0xFFff) << 16) | (a1 << 8) | a0;
        packed[(i << 2) + 1] = (a_idx >> 16) & 0xFFffFFff;
        packed[(i << 2) + 2] = (color1 << 16) | color0;
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

    max_i = packed_pix_buf->len / 2;  // 2 bytes per pixel
    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i*ucc;
        packed_pix[i] = ((((uint16)v_scale[unpacked_pix[pxl_i + v_chan]]) << 8) |
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

    packed_pix   = (uint32 *)packed_pix_buf->buf;
    unpacked_pix = (uint16 *)unpacked_pix_buf->buf;
    u_scale      = (uint16 *)u_scale_buf->buf;
    v_scale      = (uint16 *)v_scale_buf->buf;

    max_i = packed_pix_buf->len / 4;  // 4 bytes per pixel
    //loop through each pixel
    for (i=0; i < max_i; i++) {
        pxl_i = i*ucc;
        packed_pix[i] = ((((uint32)v_scale[unpacked_pix[pxl_i + v_chan]]) << 16) |
                                   u_scale[unpacked_pix[pxl_i + u_chan]]) ^ 0x80008000UL;
    }
}

static PyObject *py_unpack_dxt1(PyObject *self, PyObject *args) {
    Py_buffer bufs[7];
    sint8 pix_per_tex, ucc, i;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*bbw*:unpack_dxt1",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5], &pix_per_tex, &ucc, &bufs[6]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[0].itemsize == 2) {
        unpack_dxt1_16(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        pix_per_tex, ucc, bufs[6].buf);
    } else {
        unpack_dxt1_8(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        pix_per_tex, ucc, bufs[6].buf);
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_unpack_dxt2_3(PyObject *self, PyObject *args) {
    Py_buffer bufs[7];
    sint8 pix_per_tex, ucc, i;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*bbw*:unpack_dxt2_3",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5], &pix_per_tex, &ucc, &bufs[6]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[0].itemsize == 2) {
        unpack_dxt2_3_16(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        pix_per_tex, ucc, bufs[6].buf);
    } else {
        unpack_dxt2_3_8(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        pix_per_tex, ucc, bufs[6].buf);
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_unpack_dxt4_5(PyObject *self, PyObject *args) {
    Py_buffer bufs[7];
    sint8 pix_per_tex, ucc, i;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*bbw*:unpack_dxt4_5",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5], &pix_per_tex, &ucc, &bufs[6]))
        return Py_BuildValue("");  // return Py_None while incrementing it


    if (bufs[0].itemsize == 2) {
        unpack_dxt4_5_16(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        pix_per_tex, ucc, bufs[6].buf);
    } else {
        unpack_dxt4_5_8(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        pix_per_tex, ucc, bufs[6].buf);
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_unpack_dxt3a(PyObject *self, PyObject *args) {
    Py_buffer bufs[7];
    sint8 pix_per_tex, ucc, scc, i;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*bbbw*:unpack_dxt3a",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5], &pix_per_tex, &ucc, &scc, &bufs[6]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[0].itemsize == 2) {
        unpack_dxt3a_16(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
                        pix_per_tex, ucc, scc, bufs[6].buf);
    } else {
        unpack_dxt3a_8(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
                       pix_per_tex, ucc, scc, bufs[6].buf);
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}


static PyObject *py_unpack_dxt5a(PyObject *self, PyObject *args) {
    Py_buffer bufs[7];
    sint8 pix_per_tex, ucc, scc, i;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*bbbw*:unpack_dxt5a",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5], &pix_per_tex, &ucc, &scc, &bufs[6]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[0].itemsize == 2) {
        unpack_dxt5a_16(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
                        pix_per_tex, ucc, scc, bufs[6].buf);
    } else {
        unpack_dxt5a_8(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
                       pix_per_tex, ucc, scc, bufs[6].buf);
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_unpack_dxn(PyObject *self, PyObject *args) {
    Py_buffer bufs[7];
    sint8 pix_per_tex, ucc, i;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*bbw*:unpack_dxn",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5], &pix_per_tex, &ucc, &bufs[6]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[0].itemsize == 2) {
        unpack_dxn_16(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        pix_per_tex, ucc, bufs[6].buf);
    } else {
        unpack_dxn_8(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        pix_per_tex, ucc, bufs[6].buf);
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_unpack_ctx1(PyObject *self, PyObject *args) {
    Py_buffer bufs[7];
    sint8 pix_per_tex, ucc, i;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*bbw*:unpack_ctx1",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5], &pix_per_tex, &ucc, &bufs[6]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[0].itemsize == 2) {
        unpack_ctx1_16(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        pix_per_tex, ucc, bufs[6].buf);
    } else {
        unpack_ctx1_8(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
        pix_per_tex, ucc, bufs[6].buf);
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_unpack_vu(PyObject *self, PyObject *args) {
    Py_buffer bufs[7];
    sint8 ucc, i;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*w*w*bw*:unpack_vu",
        &bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5], &ucc, &bufs[6]))
        return Py_BuildValue("");  // return Py_None while incrementing it

    if (bufs[1].itemsize == 2) {  // unpacking v8u8
        if (bufs[0].itemsize == 2) {  // to a16r16g16b16
            unpack_v8u8_16(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
            ucc, bufs[6].buf);
        } else {  // to a8r8g8b8
            unpack_v8u8_8(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
            ucc, bufs[6].buf);
        }
    } else if (bufs[1].itemsize == 4) {  // unpacking v16u16
        if (bufs[0].itemsize == 2) {  // to a16r16g16b16
            unpack_v16u16_16(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
            ucc, bufs[6].buf);
        } else {  // to a8r8g8b8
            unpack_v16u16_8(&bufs[0], &bufs[1], &bufs[2], &bufs[3], &bufs[4], &bufs[5],
            ucc, bufs[6].buf);
        }
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_pack_dxt1(PyObject *self, PyObject *args) {
    Py_buffer bufs[5];
    sint8 pix_per_tex, can_have_alpha, i;
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
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_pack_dxt2_3(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];
    sint8 pix_per_tex, i;

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
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_pack_dxt4_5(PyObject *self, PyObject *args) {
    Py_buffer bufs[6];
    sint8 pix_per_tex, i;

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
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_pack_vu(PyObject *self, PyObject *args) {
    Py_buffer bufs[4];
    sint8 u_chan, v_chan, ucc, i;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*w*w*bbb:pack_vu",
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
    RELEASE_PY_BUFFER_ARRAY(bufs, i)

    return Py_BuildValue("");  // return Py_None while incrementing it
}

static PyObject *py_dxt_swizzle(PyObject *self, PyObject *args) {
    Py_buffer bufs[2];
    int swizzle;
    uint8 ucc;
    sint64 txl_stride, single_txl_stride;
    sint64 txl_ct_y, txl_ct_x, txl_w, txl_h;

    uint8 *src_pixels, *dst_pixels;
    sint64 tx_y_idx, tx_idx, y_idx;
    sint64 i = 0, j = 0, k, a, i_tx, i_tx_y, i_tx_yx;

    // Get the pointers to each of the array objects
    if (!PyArg_ParseTuple(args, "w*w*pbKKKK:dxt_swizzle",
        &bufs[0], &bufs[1], &swizzle, &ucc, &txl_ct_y, &txl_ct_x, &txl_w, &txl_h))
        return Py_BuildValue("");  // return Py_None while incrementing it

    src_pixels = (uint8*)bufs[0].buf;
    dst_pixels = (uint8*)bufs[1].buf;

    single_txl_stride = txl_w * ucc * bufs[0].itemsize;
    txl_stride = txl_h * txl_ct_x * single_txl_stride;

    if (bufs[0].len != bufs[1].len || bufs[0].itemsize != bufs[1].itemsize) {
        RELEASE_PY_BUFFER_ARRAY(bufs, i)
        PySys_FormatStdout("Pixel buffers supplied to dds_defs_ext.dxt_swizzle are not the same size.\n");
        return Py_BuildValue("");  // return Py_None while incrementing it
    } else if (txl_stride * txl_ct_y > bufs[0].len) {
        RELEASE_PY_BUFFER_ARRAY(bufs, i)
        PySys_FormatStdout("Pixel buffers supplied to dds_defs_ext.dxt_swizzle are smaller than required.\n");
        return Py_BuildValue("");  // return Py_None while incrementing it
    }

    for (tx_y_idx = 0; tx_y_idx < txl_ct_y; tx_y_idx++) {
        i_tx = i;
        if (swizzle) {
            for (tx_idx = 0; tx_idx < txl_ct_x; tx_idx++) {
                i_tx_y = i_tx;
                for (y_idx = 0; y_idx < txl_h; y_idx++) {
                    i_tx_yx = i_tx_y;
                    for (k = 0; k < single_txl_stride; k++) {
                        dst_pixels[j] = src_pixels[i_tx_yx];
                        j++;
                        i_tx_yx++;
                    }
                    i_tx_y += txl_ct_x * single_txl_stride;
                }
                i_tx += single_txl_stride;
            }
        } else {
            for (tx_idx = 0; tx_idx < txl_ct_x; tx_idx++) {
                i_tx_y = i_tx;
                for (y_idx = 0; y_idx < txl_h; y_idx++) {
                    i_tx_yx = i_tx_y;
                    for (k = 0; k < single_txl_stride; k++) {
                        dst_pixels[i_tx_yx] = src_pixels[j];
                        j++;
                        i_tx_yx++;
                    }
                    i_tx_y += txl_ct_x * single_txl_stride;
                }
                i_tx += single_txl_stride;
            }
        }
        i += txl_stride;
    }

    // Release the buffer objects
    RELEASE_PY_BUFFER_ARRAY(bufs, a)

    return Py_BuildValue("");  // return Py_None while incrementing it
}


/* A list of all the methods defined by this module.
"METH_VARGS" tells Python how to call the handler.
The {NULL, NULL} entry indicates the end of the method definitions.*/
static PyMethodDef dds_defs_ext_methods[] = {
    {"unpack_dxt1", py_unpack_dxt1, METH_VARARGS, ""},
    {"unpack_dxt2_3", py_unpack_dxt2_3, METH_VARARGS, ""},
    {"unpack_dxt4_5", py_unpack_dxt4_5, METH_VARARGS, ""},
    {"unpack_dxt3a", py_unpack_dxt3a, METH_VARARGS, ""},
    {"unpack_dxt5a", py_unpack_dxt5a, METH_VARARGS, ""},
    //{"unpack_dxt3a1111", py_unpack_dxt3a1111, METH_VARARGS, ""},
    {"unpack_dxn", py_unpack_dxn, METH_VARARGS, ""},
    {"unpack_ctx1", py_unpack_ctx1, METH_VARARGS, ""},
    {"unpack_vu", py_unpack_vu, METH_VARARGS, ""},

    {"pack_dxt1", py_pack_dxt1, METH_VARARGS, ""},
    {"pack_dxt2_3", py_pack_dxt2_3, METH_VARARGS, ""},
    {"pack_dxt4_5", py_pack_dxt4_5, METH_VARARGS, ""},
    //{"pack_dxt3a", py_pack_dxt3a, METH_VARARGS, ""},
    //{"pack_dxt5a", py_pack_dxt5a, METH_VARARGS, ""},
    //{"pack_dxt3a1111", py_pack_dxt3a1111, METH_VARARGS, ""},
    //{"pack_dxn", py_pack_dxn, METH_VARARGS, ""},
    //{"pack_ctx1", py_pack_ctx1, METH_VARARGS, ""},
    {"pack_vu", py_pack_vu, METH_VARARGS, ""},

    {"dxt_swizzle", py_dxt_swizzle, METH_VARARGS, ""},
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
