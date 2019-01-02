from array import array
from math import sqrt

#this will be the reference to the bitmap convertor module.
#once the module loads this will become the reference to it.
ab = None

try:
    from arbytmap.ext import dds_defs_ext
    fast_dds_defs = True
except Exception:
    fast_dds_defs = False


def get_texel_pixel_count(width, height):
    return min(width, 4) * min(height, 4)


def initialize():
    """FOR DXT FORMATS, ALPHA CHANNELS ARE TREATED SPECIALLY,
    BUT ARE EXPLICITELY PLACED HERE TO MAKE SURE THEY DONT
    CAUSE CHANNEL MAP SWAPPING PROBLEMS"""

    ab.FORMAT_DXT1 = "DXT1"
    ab.FORMAT_DXT2 = "DXT2"
    ab.FORMAT_DXT3 = "DXT3"
    ab.FORMAT_DXT4 = "DXT4"
    ab.FORMAT_DXT5 = "DXT5"

    # uses only the alpha channel of dxt3
    ab.FORMAT_DXT3A  = "DXT3A"
    ab.FORMAT_DXT3Y  = "DXT3Y"
    ab.FORMAT_DXT3AY = "DXT3AY"
    # uses only the alpha channel of dxt3, and each bit is
    # used as an stencil mask for each of the ARGB channels.
    # this format is basically A1R1G1B1 with a dxt texel swizzle
    ab.FORMAT_DXT3A1111 = "DXT3A1111"   #NOT YET IMPLEMENTED

    ab.FORMAT_DXT5A  = "DXT5A"
    ab.FORMAT_DXT5Y  = "DXT5Y"
    ab.FORMAT_DXT5AY = "DXT5AY"

    # normal map formats
    ab.FORMAT_DXT5NM = "DXT5NM"         #NOT YET IMPLEMENTED
    ab.FORMAT_DXN    = "DXN"
    ab.FORMAT_CTX1   = "CTX1"
    ab.FORMAT_V8U8   = "V8U8"
    ab.FORMAT_V16U16 = "V16U16"
    ab.FORMAT_R8G8   = "R8G8"
    ab.FORMAT_R16G16 = "R16G16"
    ab.FORMAT_G8B8   = "G8B8"
    ab.FORMAT_G16B16 = "G16B16"

    combine = lambda base, **main: {
        k: (base[k] if k not in main else main[k]) for k in
        set(main.keys()).union(base.keys())}

    dxt_specs = dict(
        compressed=True, dds_format=True, raw_format=False,
        packed_size_calc=dxt_packed_size_calc,
        packed_width_calc=packed_dxt_dimension_calc,
        packed_height_calc=packed_dxt_dimension_calc,
        packed_typecode='I', packed_field_sizes=(2, ),
        block_width=4, block_height=4, 
        )

    ab.register_format(ab.FORMAT_DXT1, 1, **combine(
        dxt_specs, bpp=4, depths=(8, 8, 8, 8),
        unpacker=unpack_dxt1, packer=pack_dxt1))

    for fmt in (ab.FORMAT_DXT2, ab.FORMAT_DXT3):
        ab.register_format(fmt, 1, **combine(
            dxt_specs, bpp=8, depths=(8, 8, 8, 8),
            premultiplied=(fmt == ab.FORMAT_DXT2),
            unpacker=unpack_dxt2_3, packer=pack_dxt2_3))

    for fmt in (ab.FORMAT_DXT4, ab.FORMAT_DXT5):
        ab.register_format(fmt, 1, **combine(
            dxt_specs, bpp=8, depths=(8, 8, 8, 8),
            premultiplied=(fmt == ab.FORMAT_DXT4),
            unpacker=unpack_dxt4_5, packer=pack_dxt4_5))

    for fmt in (ab.FORMAT_DXT3A, ab.FORMAT_DXT3Y):
        ab.register_format(fmt, 1, **combine(
            dxt_specs, bpp=4, depths=(8,)),
            unpacker=unpack_dxt3a, packer=pack_dxt3a)

    for fmt in (ab.FORMAT_DXT5A, ab.FORMAT_DXT5Y):
        ab.register_format(fmt, 1, **combine(
            dxt_specs, bpp=4, depths=(8,)),
            unpacker=unpack_dxt5a, packer=pack_dxt5a)

    ab.register_format(ab.FORMAT_DXT3AY, 1, **combine(
        dxt_specs, bpp=8, depths=(8, 8),
        unpacker=unpack_dxt3a, packer=pack_dxt3a))

    ab.register_format(ab.FORMAT_DXT5AY, 1, **combine(
        dxt_specs, bpp=8, depths=(8, 8),
        unpacker=unpack_dxt5a, packer=pack_dxt5a))

    ab.register_format(ab.FORMAT_DXN, 1, **combine(
        dxt_specs, bpp=8, depths=(8, 8, 8),
        unpacker=unpack_dxn, packer=pack_dxn))

    ab.register_format(ab.FORMAT_CTX1, 1, **combine(
        dxt_specs, bpp=4, depths=(8, 8, 8),
        unpacker=unpack_ctx1, packer=pack_ctx1))

    ab.register_format(ab.FORMAT_V8U8, 1, bpp=16, dds_format=True,
                       unpacker=unpack_v8u8, packer=pack_v8u8,
                       depths=(8,8,8), offsets=(0,8,0),
                       masks=(0, 0xFF, 0xFF), packed_field_sizes=(2, ))

    ab.register_format(ab.FORMAT_V16U16, 1, bpp=32, dds_format=True,
                       unpacker=unpack_v16u16, packer=pack_v16u16,
                       depths=(16,16,16), offsets=(0,16,0),
                       masks=(0, 0xFFff, 0xFFff), packed_field_sizes=(4, ))

    ab.register_format(ab.FORMAT_R8G8, 1, bpp=16, dds_format=True,
                       unpacker=unpack_r8g8, packer=pack_r8g8,
                       depths=(8,8,8), offsets=(0,8,0),
                       masks=(0, 0xFF, 0xFF), packed_field_sizes=(2, ))

    ab.register_format(ab.FORMAT_R16G16, 1, bpp=32, dds_format=True,
                       unpacker=unpack_r16g16, packer=pack_r16g16,
                       depths=(16,16,16), offsets=(0,16,0),
                       masks=(0, 0xFFff, 0xFFff), packed_field_sizes=(4, ))

    ab.register_format(ab.FORMAT_G8B8, 1, bpp=16, dds_format=True,
                       unpacker=unpack_g8b8, packer=pack_g8b8,
                       depths=(8,8,8), offsets=(8,0,0),
                       masks=(0xFF, 0xFF, 0), packed_field_sizes=(2, ))

    ab.register_format(ab.FORMAT_G16B16, 1, bpp=32, dds_format=True,
                       unpacker=unpack_g16b16, packer=pack_g16b16,
                       depths=(16,16,16), offsets=(16,0,0),
                       masks=(0xFFff, 0xFFff, 0), packed_field_sizes=(4, ))
    

def _dxt_swizzle(src_pixels, orig_width, orig_height, channel_ct, swizz=False):
    width, height = clip_dxt_dimensions(orig_width, orig_height)
    txl_ct_x = 1 if width  < 4 else width  // 4
    txl_ct_y = 1 if height < 4 else height // 4

    txl_w = 4 if txl_ct_x > 1 else orig_width
    txl_h = 4 if txl_ct_y > 1 else orig_height

    assert len(src_pixels) % channel_ct == 0
    dst_pixels = ab.bitmap_io.make_array(src_pixels.typecode, len(src_pixels))

    # 4 channels per pixel, 16 pixels per texel
    if fast_dds_defs:
        dds_defs_ext.dxt_swizzle(
            src_pixels, dst_pixels, swizz, channel_ct,
            txl_ct_y, txl_ct_x, txl_w, txl_h)
    else:
        txl_stride = txl_h * width * channel_ct
        tx_block_offs = tuple(range(0, width * channel_ct, txl_w * channel_ct))
        y_block_offs  = tuple(y * width * channel_ct for y in range(txl_h))
        x_block_offs  = tuple(range(0, txl_w * channel_ct, channel_ct))
        c_block_offs  = tuple(range(channel_ct))
        i = j = 0
        for tx_y in range(txl_ct_y):
            if swizz:
                for tx in tx_block_offs:
                    i_tx = i + tx
                    for y in y_block_offs:
                        i_tx_y = i_tx + y
                        for x in x_block_offs:
                            i_tx_yx = i_tx_y + x
                            for c in c_block_offs:
                                dst_pixels[j] = src_pixels[i_tx_yx + c]
                                j += 1
            else:
                for tx in tx_block_offs:
                    i_tx = i + tx
                    for y in y_block_offs:
                        i_tx_y = i_tx + y
                        for x in x_block_offs:
                            i_tx_yx = i_tx_y + x
                            for c in c_block_offs:
                                dst_pixels[i_tx_yx + c] = src_pixels[j]
                                j += 1

            i += txl_stride

    return dst_pixels


def unswizzle_dxt(pixels, orig_width, orig_height, channel_ct):
    return _dxt_swizzle(pixels, orig_width, orig_height, channel_ct, False)


def swizzle_dxt(pixels, orig_width, orig_height, channel_ct):
    return _dxt_swizzle(pixels, orig_width, orig_height, channel_ct, True)


def dxt_packed_size_calc(fmt, width, height, depth=1):
    width, height = clip_dxt_dimensions(width, height)
    return (ab.BITS_PER_PIXEL[fmt] * height * width * depth)//8


def packed_dxt_dimension_calc(dim, mip_level, tiled=False):
    dim = dim >> mip_level
    if dim <= 4: return 4
    return dim + (4 - (dim % 4)) % 4


def clip_dxt_dimensions(width, height):
    return (packed_dxt_dimension_calc(width, 0),
            packed_dxt_dimension_calc(height, 0))


def unpack_dxt1(arby, bitmap_index, width, height, depth=1):
    packed = arby.texture_block[bitmap_index]
    assert packed.typecode == 'I'

    unpack_code = arby._UNPACK_ARRAY_CODE
    unpack_size = ab.PIXEL_ENCODING_SIZES[unpack_code]
    unpack_max = (1<<(unpack_size*8)) - 1

    ucc = arby.unpacked_channel_count
    width, height, depth = ab.clip_dimensions(width, height, depth)
    texel_width, texel_height, _ = ab.clip_dimensions(width//4, height//4)

    pixels_per_texel = (width//texel_width)*(height//texel_height)
    dxt_width, dxt_height = clip_dxt_dimensions(width, height)
    unpacked = ab.bitmap_io.make_array(unpack_code, dxt_width*dxt_height*ucc)

    upscales = list(arby.channel_upscalers)
    chan_map = list(arby.channel_mapping)

    while len(upscales) < 4: upscales.append(array(upscales[0].typecode, [0]))
    while len(chan_map) < 4: chan_map.append(-1)

    if fast_dds_defs:
        a_scale, r_scale, g_scale, b_scale = upscales[: 4]
        dds_defs_ext.unpack_dxt1(
            unpacked, packed, a_scale, r_scale, g_scale, b_scale,
            pixels_per_texel, ucc, array("b", chan_map[: 4]))
    else:
        channels_per_texel = ucc*pixels_per_texel
        pixel_indices = range(pixels_per_texel)
        upscales = tuple(tuple(scale) for scale in upscales)

        #create the arrays to hold the color channel data
        c_0 = [255,0,0,0]
        c_1 = [255,0,0,0]
        c_2 = [255,0,0,0]
        c_3 = [255,0,0,0]
        transparent = [0,0,0,0]

        #stores the colors in a way we can easily access them
        colors = [c_0, c_1, c_2, c_3]

        #loop through each texel
        for i in range(len(packed)//2):
            pxl_i = i*channels_per_texel
            j = i*2

            #if the format DXT1 then the two entries in the array
            #are the colors and the color indexing in that order.
            color0 = packed[j] & 65535
            color1 = (packed[j] >> 16) & 65535
            color_idx = packed[j+1]

            #unpack the colors
            c_0[1] = (((color0>>11) & 31)*255 + 15)//31
            c_1[1] = (((color1>>11) & 31)*255 + 15)//31
            c_0[2] = (((color0>>5) & 63)*255 + 31)//63
            c_1[2] = (((color1>>5) & 63)*255 + 31)//63
            c_1[3] = ((color1 & 31)*255 + 15)//31
            c_0[3] = ((color0 & 31)*255 + 15)//31

            #if the first color is a larger integer
            #then color key transparency is NOT used
            if color0 > color1:
                c_2[1] = (c_0[1]*2 + c_1[1])//3
                c_2[2] = (c_0[2]*2 + c_1[2])//3
                c_2[3] = (c_0[3]*2 + c_1[3])//3
                c_3[:] = [255, (c_0[1] + 2*c_1[1])//3,
                          (c_0[2] + 2*c_1[2])//3, (c_0[3] + 2*c_1[3])//3]
            else:
                c_2[1] = (c_0[1]+c_1[1])//2
                c_2[2] = (c_0[2]+c_1[2])//2
                c_2[3] = (c_0[3]+c_1[3])//2
                c_3[:] = transparent

            for j in pixel_indices:
                color = colors[(color_idx >> (j*2))&3]
                off = j*ucc + pxl_i
                dst_chan = 0
                for src_chan in chan_map:
                    if src_chan < 0 and dst_chan == 0:
                        # alpha and not reading alpha. set to full white
                        unpacked[off] = unpack_max
                    elif src_chan >= 0:
                        unpacked[off + dst_chan] = upscales[dst_chan][color[src_chan]]
                    dst_chan += 1

    return unswizzle_dxt(unpacked, width, height * depth, ucc)


def unpack_dxt2_3(arby, bitmap_index, width, height, depth=1):
    packed = arby.texture_block[bitmap_index]
    assert packed.typecode == 'I'

    unpack_code = arby._UNPACK_ARRAY_CODE
    unpack_size = ab.PIXEL_ENCODING_SIZES[unpack_code]
    unpack_max = (1<<(unpack_size*8)) - 1

    ucc = arby.unpacked_channel_count
    width, height, depth = ab.clip_dimensions(width, height, depth)
    texel_width, texel_height, _ = ab.clip_dimensions(width//4, height//4)

    pixels_per_texel = (width//texel_width)*(height//texel_height)

    #create a new array to hold the pixels after we unpack them
    dxt_width, dxt_height = clip_dxt_dimensions(width, height)
    unpacked = ab.bitmap_io.make_array(unpack_code, dxt_width*dxt_height*ucc)

    upscales = list(arby.channel_upscalers)
    chan_map = list(arby.channel_mapping)

    while len(upscales) < 4: upscales.append(array(upscales[0].typecode, [0]))
    while len(chan_map) < 4: chan_map.append(-1)

    if fast_dds_defs:
        a_scale, r_scale, g_scale, b_scale = upscales[: 4]
        dds_defs_ext.unpack_dxt2_3(
            unpacked, packed, a_scale, r_scale, g_scale, b_scale,
            pixels_per_texel, ucc, array("b", chan_map[: 4]))
    else:
        channels_per_texel = ucc*pixels_per_texel
        pixel_indices = range(pixels_per_texel)
        upscales = tuple(tuple(scale) for scale in upscales)

        #create the arrays to hold the color channel data
        c_0 = [255,0,0,0]
        c_1 = [255,0,0,0]
        c_2 = [255,0,0,0]
        c_3 = [255,0,0,0]

        #stores the colors in a way we can easily access them
        colors = [c_0, c_1, c_2, c_3]

        #loop through each texel
        for i in range(len(packed)//4):
            pxl_i = i*channels_per_texel
            j = i*4

            #DXT2/3 is much simpler than DXT4/5
            alpha = (packed[j+1]<<32) | packed[j]
            color0 = packed[j+2] & 65535
            color1 = (packed[j+2] >> 16) & 65535
            color_idx = packed[j+3]

            if color0 < color1:
                color0, color1 = color1, color0

            #unpack the colors
            c_0[1] = (((color0>>11) & 31)*255 + 15)//31
            c_1[1] = (((color1>>11) & 31)*255 + 15)//31
            c_0[2] = (((color0>>5) & 63)*255 + 31)//63
            c_1[2] = (((color1>>5) & 63)*255 + 31)//63
            c_1[3] = ((color1 & 31)*255 + 15)//31
            c_0[3] = ((color0 & 31)*255 + 15)//31

            c_2[1] = (c_0[1]*2 + c_1[1])//3
            c_2[2] = (c_0[2]*2 + c_1[2])//3
            c_2[3] = (c_0[3]*2 + c_1[3])//3

            c_3[1] = (c_0[1] + c_1[1]*2)//3
            c_3[2] = (c_0[2] + c_1[2]*2)//3
            c_3[3] = (c_0[3] + c_1[3]*2)//3

            for j in pixel_indices:
                color = colors[(color_idx >> (j*2))&3]
                off = j*ucc + pxl_i
                a = (((alpha >> (j*4)) & 15)*255)//15

                dst_chan = 0
                for src_chan in chan_map:
                    if src_chan < 0 and dst_chan == 0:
                        # alpha and not reading alpha. set to full white
                        unpacked[off] = unpack_max
                    elif src_chan > 0:
                        unpacked[off + dst_chan] = upscales[dst_chan][color[src_chan]]
                    elif src_chan == 0:
                        unpacked[off + dst_chan] = upscales[dst_chan][a]
                    dst_chan += 1

    return unswizzle_dxt(unpacked, width, height * depth, ucc)


def unpack_dxt4_5(arby, bitmap_index, width, height, depth=1):
    packed = arby.texture_block[bitmap_index]
    assert packed.typecode == 'I'

    unpack_code = arby._UNPACK_ARRAY_CODE
    unpack_size = ab.PIXEL_ENCODING_SIZES[unpack_code]
    unpack_max = (1<<(unpack_size*8)) - 1

    ucc = arby.unpacked_channel_count
    width, height, depth = ab.clip_dimensions(width, height, depth)
    texel_width, texel_height, _ = ab.clip_dimensions(width//4, height//4)

    pixels_per_texel = (width//texel_width)*(height//texel_height)

    #create a new array to hold the pixels after we unpack them
    dxt_width, dxt_height = clip_dxt_dimensions(width, height)
    unpacked = ab.bitmap_io.make_array(unpack_code, dxt_width*dxt_height*ucc)

    upscales = list(arby.channel_upscalers)
    chan_map = list(arby.channel_mapping)

    while len(upscales) < 4: upscales.append(array(upscales[0].typecode, [0]))
    while len(chan_map) < 4: chan_map.append(-1)

    if fast_dds_defs:
        a_scale, r_scale, g_scale, b_scale = upscales[: 4]
        dds_defs_ext.unpack_dxt4_5(
            unpacked, packed, a_scale, r_scale, g_scale, b_scale,
            pixels_per_texel, ucc, array("b", chan_map[: 4]))
    else:
        a_lookup = [0,0,0,0,0,0,0,0]

        channels_per_texel = ucc*pixels_per_texel
        pixel_indices = range(pixels_per_texel)
        upscales = tuple(tuple(scale) for scale in upscales)

        #create the arrays to hold the color channel data
        c_0 = [255,0,0,0]
        c_1 = [255,0,0,0]
        c_2 = [255,0,0,0]
        c_3 = [255,0,0,0]

        #stores the colors in a way we can easily access them
        colors = [c_0, c_1, c_2, c_3]

        #loop through each texel
        for i in range(len(packed)//4):
            pxl_i = i*channels_per_texel
            j = i*4

            a_lookup[0] = alpha0 = packed[j] & 255
            a_lookup[1] = alpha1 = (packed[j] >> 8) & 255
            alpha_idx = ((packed[j]>>16) & 65535) | (packed[j+1] << 16)

            #depending on which alpha is larger
            #the indexing is calculated differently
            if alpha0 > alpha1:
                a_lookup[2] = (alpha0*6 + alpha1)//7
                a_lookup[3] = (alpha0*5 + alpha1*2)//7
                a_lookup[4] = (alpha0*4 + alpha1*3)//7
                a_lookup[5] = (alpha0*3 + alpha1*4)//7
                a_lookup[6] = (alpha0*2 + alpha1*5)//7
                a_lookup[7] = (alpha0   + alpha1*6)//7
            else:
                a_lookup[2] = (alpha0*4 + alpha1)//5
                a_lookup[3] = (alpha0*3 + alpha1*2)//5
                a_lookup[4] = (alpha0*2 + alpha1*3)//5
                a_lookup[5] = (alpha0   + alpha1*4)//5
                a_lookup[6] = 0
                a_lookup[7] = 255

            #half of the first array entry in DXT4/5 format is both
            #alpha values and the first third of the indexing
            color0 = packed[j+2] & 65535
            color1 = (packed[j+2]>>16) & 65535
            color_idx = packed[j+3]

            if color0 < color1:
                color0, color1 = color1, color0

            #unpack the colors
            c_0[1] = (((color0>>11) & 31)*255 + 15)//31
            c_1[1] = (((color1>>11) & 31)*255 + 15)//31
            c_0[2] = (((color0>>5) & 63)*255 + 31)//63
            c_1[2] = (((color1>>5) & 63)*255 + 31)//63
            c_1[3] = ((color1 & 31)*255 + 15)//31
            c_0[3] = ((color0 & 31)*255 + 15)//31

            c_2[1] = (c_0[1]*2 + c_1[1])//3
            c_2[2] = (c_0[2]*2 + c_1[2])//3
            c_2[3] = (c_0[3]*2 + c_1[3])//3

            c_3[1] = (c_0[1] + c_1[1]*2)//3
            c_3[2] = (c_0[2] + c_1[2]*2)//3
            c_3[3] = (c_0[3] + c_1[3]*2)//3

            for j in pixel_indices:
                color = colors[(color_idx >> (j*2))&3]
                off = j*ucc + pxl_i
                a = a_lookup[(alpha_idx >> (j*3))&7]

                dst_chan = 0
                for src_chan in chan_map:
                    if src_chan < 0 and dst_chan == 0:
                        # alpha and not reading alpha. set to full white
                        unpacked[off] = unpack_max
                    elif src_chan > 0:
                        unpacked[off + dst_chan] = upscales[dst_chan][color[src_chan]]
                    elif src_chan == 0:
                        unpacked[off + dst_chan] = upscales[dst_chan][a]
                    dst_chan += 1

    return unswizzle_dxt(unpacked, width, height * depth, ucc)


def unpack_dxt3a(arby, bitmap_index, width, height, depth=1):
    packed = arby.texture_block[bitmap_index]
    assert packed.typecode == 'I'

    unpack_code = arby._UNPACK_ARRAY_CODE
    unpack_size = ab.PIXEL_ENCODING_SIZES[unpack_code]
    unpack_max = (1<<(unpack_size*8)) - 1

    ucc = arby.unpacked_channel_count
    scc = arby.source_channel_count
    width, height, depth = ab.clip_dimensions(width, height, depth)
    texel_width, texel_height, _ = ab.clip_dimensions(width//4, height//4)

    pixels_per_texel = (width//texel_width)*(height//texel_height)

    #create a new array to hold the pixels after we unpack them
    dxt_width, dxt_height = clip_dxt_dimensions(width, height)
    unpacked = ab.bitmap_io.make_array(unpack_code, dxt_width*dxt_height*ucc)

    upscales = list(arby.channel_upscalers)
    chan_map = list(arby.channel_mapping)

    while len(upscales) < 4: upscales.append(array(upscales[0].typecode, [0]))
    while len(chan_map) < 4: chan_map.append(-1)

    if fast_dds_defs:
        a_scale, r_scale, g_scale, b_scale = upscales[: 4]
        dds_defs_ext.unpack_dxt3a(
            unpacked, packed, a_scale, r_scale, g_scale, b_scale,
            pixels_per_texel, ucc, scc, array("b", chan_map[: 4]))
    else:
        channels_per_texel = ucc*pixels_per_texel
        pixel_indices = range(pixels_per_texel)
        upscales = tuple(tuple(scale) for scale in upscales)
        #loop through each texel
        for dst_chan in range(ucc):
            scale = upscales[dst_chan]
            src_chan = chan_map[dst_chan]

            if src_chan < 0:
                # not reading anything for this destination channel.
                # either leave it full black, or set it to full white.
                if dst_chan == 0:
                    # set alpha to full white
                    for off in range(0, len(unpacked), ucc):
                        unpacked[off] = unpack_max
                continue

            pxl_i = dst_chan
            for i in range(2 * src_chan, len(packed), 2 * scc):
                alpha = (packed[i+1]<<32) | packed[i]
                for j in pixel_indices:
                    unpacked[pxl_i] = scale[(((alpha>>(j*4)) & 15)*255)//15]
                    pxl_i += ucc

    return unswizzle_dxt(unpacked, width, height * depth, ucc)


def unpack_dxt5a(arby, bitmap_index, width, height, depth=1):
    packed = arby.texture_block[bitmap_index]
    assert packed.typecode == 'I'

    unpack_code = arby._UNPACK_ARRAY_CODE
    unpack_size = ab.PIXEL_ENCODING_SIZES[unpack_code]
    unpack_max = (1<<(unpack_size*8)) - 1

    ucc = arby.unpacked_channel_count
    scc = arby.source_channel_count
    width, height, depth = ab.clip_dimensions(width, height, depth)
    texel_width, texel_height, _ = ab.clip_dimensions(width//4, height//4)

    pixels_per_texel = (width//texel_width)*(height//texel_height)

    #create a new array to hold the pixels after we unpack them
    dxt_width, dxt_height = clip_dxt_dimensions(width, height)
    unpacked = ab.bitmap_io.make_array(unpack_code, dxt_width*dxt_height*ucc)

    upscales = list(arby.channel_upscalers)
    chan_map = list(arby.channel_mapping)

    while len(upscales) < 4: upscales.append(array(upscales[0].typecode, [0]))
    while len(chan_map) < 4: chan_map.append(-1)

    if fast_dds_defs:
        a_scale, r_scale, g_scale, b_scale = upscales[: 4]
        dds_defs_ext.unpack_dxt5a(
            unpacked, packed, a_scale, r_scale, g_scale, b_scale,
            pixels_per_texel, ucc, scc, array("b", chan_map[: 4]))
    else:
        lookup = [0,0,0,0,0,0,0,0]

        channels_per_texel = ucc*pixels_per_texel
        pixel_indices = range(pixels_per_texel)
        upscales = tuple(tuple(scale) for scale in upscales)
        #loop through each texel
        for dst_chan in range(ucc):
            scale = upscales[dst_chan]
            src_chan = chan_map[dst_chan]

            if src_chan < 0:
                # not reading anything for this destination channel.
                # either leave it full black, or set it to full white.
                if dst_chan == 0:
                    # set alpha to full white
                    for off in range(0, len(unpacked), ucc):
                        unpacked[off] = unpack_max
                continue

            pxl_i = dst_chan
            for i in range(2 * src_chan, len(packed), 2 * scc):
                lookup[0] = val0 = packed[i] & 255
                lookup[1] = val1 = (packed[i] >> 8) & 255
                idx = ((packed[i]>>16) & 65535) | (packed[i+1] << 16)

                # depending on which value is larger
                # the indexing is calculated differently
                if val0 > val1:
                    lookup[2] = (val0*6 + val1)//7
                    lookup[3] = (val0*5 + val1*2)//7
                    lookup[4] = (val0*4 + val1*3)//7
                    lookup[5] = (val0*3 + val1*4)//7
                    lookup[6] = (val0*2 + val1*5)//7
                    lookup[7] = (val0   + val1*6)//7
                else:
                    lookup[2] = (val0*4 + val1)//5
                    lookup[3] = (val0*3 + val1*2)//5
                    lookup[4] = (val0*2 + val1*3)//5
                    lookup[5] = (val0   + val1*4)//5
                    lookup[6] = 0
                    lookup[7] = 255

                for j in pixel_indices:
                    unpacked[pxl_i] = scale[lookup[(idx >> (j*3))&7]]
                    pxl_i += ucc

    return unswizzle_dxt(unpacked, width, height * depth, ucc)


def unpack_dxn(arby, bitmap_index, width, height, depth=1):
    packed = arby.texture_block[bitmap_index]
    assert packed.typecode == 'I'

    unpack_code = arby._UNPACK_ARRAY_CODE
    unpack_size = ab.PIXEL_ENCODING_SIZES[unpack_code]
    unpack_max = (1<<(unpack_size*8)) - 1

    zero_point = sign_mask = 0x80
    mask = sign_mask - 1
    mask_sq = mask**2

    ucc = arby.unpacked_channel_count
    width, height, depth = ab.clip_dimensions(width, height, depth)
    texel_width, texel_height, _ = ab.clip_dimensions(width//4, height//4)

    pixels_per_texel = (width//texel_width)*(height//texel_height)
    channels_per_texel = ucc*pixels_per_texel

    #create a new array to hold the pixels after we unpack them
    dxt_width, dxt_height = clip_dxt_dimensions(width, height)
    unpacked = ab.bitmap_io.make_array(unpack_code, dxt_width*dxt_height*ucc)

    upscales = list(arby.channel_upscalers)
    chan_map = list(arby.channel_mapping)

    while len(upscales) < 4: upscales.append(array(upscales[0].typecode, [0]))
    while len(chan_map) < 4: chan_map.append(-1)

    if fast_dds_defs:
        a_scale, r_scale, g_scale, b_scale = upscales[: 4]
        dds_defs_ext.unpack_dxn(
            unpacked, packed, a_scale, r_scale, g_scale, b_scale,
            pixels_per_texel, ucc, array("b", chan_map[: 4]))
    else:
        # convert to tuples for faster access
        upscales = tuple(tuple(scale) for scale in upscales)
        pixel_indices = range(pixels_per_texel)
        r_lookup = [0,0,0,0,0,0,0,0]
        g_lookup = [0,0,0,0,0,0,0,0]

        #loop through each texel
        for i in range(len(packed)//4):
            pxl_i = i*channels_per_texel
            j = i*4

            g_lookup[0] = g0 = packed[j]&255
            g_lookup[1] = g1 = (packed[j]>>8)&255
            g_idx = ((packed[j]>>16)&65535) + (packed[j+1]<<16)

            r_lookup[0] = r0 = packed[j+2]&255
            r_lookup[1] = r1 = (packed[j+2]>>8)&255
            r_idx = ((packed[j+2]>>16)&65535) + (packed[j+3]<<16)

            #depending on which alpha value is larger
            #the indexing is calculated differently
            if g0 > g1:
                g_lookup[2] = (g0*6 + g1  )//7
                g_lookup[3] = (g0*5 + g1*2)//7
                g_lookup[4] = (g0*4 + g1*3)//7
                g_lookup[5] = (g0*3 + g1*4)//7
                g_lookup[6] = (g0*2 + g1*5)//7
                g_lookup[7] = (g0   + g1*6)//7
            else:
                g_lookup[2] = (g0*4 + g1  )//5
                g_lookup[3] = (g0*3 + g1*2)//5
                g_lookup[4] = (g0*2 + g1*3)//5
                g_lookup[5] = (g0   + g1*4)//5
                g_lookup[6] = 0
                g_lookup[7] = 255

            if r0 > r1:
                r_lookup[2] = (r0*6 + r1  )//7
                r_lookup[3] = (r0*5 + r1*2)//7
                r_lookup[4] = (r0*4 + r1*3)//7
                r_lookup[5] = (r0*3 + r1*4)//7
                r_lookup[6] = (r0*2 + r1*5)//7
                r_lookup[7] = (r0   + r1*6)//7
            else:
                r_lookup[2] = (r0*4 + r1  )//5
                r_lookup[3] = (r0*3 + r1*2)//5
                r_lookup[4] = (r0*2 + r1*3)//5
                r_lookup[5] = (r0   + r1*4)//5
                r_lookup[6] = 0
                r_lookup[7] = 255

            for k in pixel_indices:
                g = y = g_lookup[(g_idx>>(k*3))&7]
                r = x = r_lookup[(r_idx>>(k*3))&7]

                off = k*ucc + pxl_i

                # we're normalizing the coordinates
                # here, not just unpacking them
                x = r&mask if r&sign_mask else zero_point - r
                y = g&mask if g&sign_mask else zero_point - g

                d = mask_sq - x**2 - y**2
                if d > 0:
                    b = int(sqrt(d)) + zero_point
                else:
                    b = zero_point
                    n_len = sqrt(mask_sq - d)/mask
                    x = int(x/n_len)
                    y = int(y/n_len)

                    r = x + zero_point if r&sign_mask else zero_point - x
                    g = y + zero_point if g&sign_mask else zero_point - y

                color = [0, r, g, b]
                dst_chan = 0
                for src_chan in chan_map:
                    if src_chan <= 0 or dst_chan == 0:
                        # set alpha to full white
                        unpacked[off] = unpack_max
                    elif src_chan >= 0:
                        unpacked[off + dst_chan] = upscales[dst_chan][color[src_chan]]
                    dst_chan += 1

    return unswizzle_dxt(unpacked, width, height * depth, ucc)


def unpack_ctx1(arby, bitmap_index, width, height, depth=1):
    packed = arby.texture_block[bitmap_index]
    assert packed.typecode == 'I'

    unpack_code = arby._UNPACK_ARRAY_CODE
    unpack_size = ab.PIXEL_ENCODING_SIZES[unpack_code]
    unpack_max = (1<<(unpack_size*8)) - 1

    zero_point = sign_mask = 0x80
    mask = sign_mask - 1
    mask_sq = mask**2

    ucc = arby.unpacked_channel_count
    width, height, depth = ab.clip_dimensions(width, height, depth)
    texel_width, texel_height, _ = ab.clip_dimensions(width//4, height//4)

    pixels_per_texel = (width//texel_width)*(height//texel_height)
    channels_per_texel = ucc*pixels_per_texel

    pixel_indices = range(pixels_per_texel)

    #create a new array to hold the pixels after we unpack them
    dxt_width, dxt_height = clip_dxt_dimensions(width, height)
    unpacked = ab.bitmap_io.make_array(unpack_code, dxt_width*dxt_height*ucc)

    upscales = list(arby.channel_upscalers)
    chan_map = list(arby.channel_mapping)

    while len(upscales) < 4: upscales.append(array(upscales[0].typecode, [0]))
    while len(chan_map) < 4: chan_map.append(-1)

    if fast_dds_defs:
        a_scale, r_scale, g_scale, b_scale = upscales[: 4]
        dds_defs_ext.unpack_ctx1(
            unpacked, packed, a_scale, r_scale, g_scale, b_scale,
            pixels_per_texel, ucc, array("b", chan_map[: 4]))
    else:
        #create the arrays to hold the color channel data
        c_0 = [0,0,0,0]
        c_1 = [0,0,0,0]
        c_2 = [0,0,0,0]
        c_3 = [0,0,0,0]

        #stores the colors in a way we can easily access them
        colors = [c_0, c_1, c_2, c_3]

        # convert to tuples for faster access
        upscales = tuple(tuple(scale) for scale in upscales)

        #loop through each texel
        for i in range(len(packed)//2):
            j = i*2
            pxl_i = i*channels_per_texel

            values = packed[j]
            idx = packed[j+1]

            # unpack the colors
            c_0[1] = x0 = r0 = (values) & 255
            c_0[2] = y0 = g0 = (values>>8) & 255
            c_1[1] = x1 = r1 = (values>>16) & 255
            c_1[2] = y1 = g1 = (values>>24) & 255

            #calculate the z-components
            # we're normalizing the coordinates here, not just unpacking them
            x0 = x0&mask if x0&sign_mask else zero_point - x0
            y0 = y0&mask if y0&sign_mask else zero_point - y0
            x1 = x1&mask if x1&sign_mask else zero_point - x1
            y1 = y1&mask if y1&sign_mask else zero_point - y1

            d = mask_sq - x0**2 - y0**2
            if d > 0:
                b0 = int(sqrt(d)) + zero_point
            else:
                b0 = zero_point
                n_len = sqrt(mask_sq - d)/mask
                x0 = int(x0/n_len)
                y0 = int(y0/n_len)

                r0 = x0 + zero_point if r0&sign_mask else zero_point - x0
                g0 = y0 + zero_point if g0&sign_mask else zero_point - y0

            d = mask_sq - x1**2 - y1**2
            if d > 0:
                b1 = int(sqrt(d)) + zero_point
            else:
                b1 = zero_point
                n_len = sqrt(mask_sq - d)/mask
                x1 = int(x1/n_len)
                y1 = int(y1/n_len)

                r1 = x1 + zero_point if r1&sign_mask else zero_point - x1
                g1 = y1 + zero_point if g1&sign_mask else zero_point - y1

            # store the normalized colors
            c_0[1] = r0; c_1[1] = r1
            c_0[2] = g0; c_1[2] = g1
            c_0[3] = b0; c_1[3] = b1

            # calculate the in-between colors
            c_2[1] = (c_0[1]*2 + c_1[1])//3
            c_2[2] = (c_0[2]*2 + c_1[2])//3
            c_2[3] = (c_0[3]*2 + c_1[3])//3

            c_3[1] = (c_0[1] + c_1[1]*2)//3
            c_3[2] = (c_0[2] + c_1[2]*2)//3
            c_3[3] = (c_0[3] + c_1[3]*2)//3

            for k in pixel_indices:
                color = colors[(idx >> (k*2))&3]
                off = k*ucc + pxl_i
                dst_chan = 0
                for src_chan in chan_map:
                    if src_chan <= 0 or dst_chan == 0:
                        # set alpha to full white
                        unpacked[off + dst_chan] = unpack_max
                    elif src_chan >= 0:
                        unpacked[off + dst_chan] = upscales[dst_chan][color[src_chan]]
                    dst_chan += 1

    return unswizzle_dxt(unpacked, width, height * depth, ucc)


def unpack_v8u8(arby, bitmap_index, width, height, depth=1):
    return unpack_vu(arby, bitmap_index, width, height, depth, 8)


def unpack_v16u16(arby, bitmap_index, width, height, depth=1):
    return unpack_vu(arby, bitmap_index, width, height, depth, 16)


def unpack_vu(arby, bitmap_index, width, height, depth=1, bpc=8):
    packed = arby.texture_block[bitmap_index]

    #create a new array to hold the pixels after we unpack them
    unpack_code = arby._UNPACK_ARRAY_CODE
    unpack_size = ab.PIXEL_ENCODING_SIZES[unpack_code]
    unpack_max = (1<<(unpack_size*8)) - 1

    ucc = arby.unpacked_channel_count
    bytes_per_pixel = unpack_size*ucc
    unpacked = ab.bitmap_io.make_array(
        unpack_code, width*height, bytes_per_pixel)

    upscales = list(arby.channel_upscalers)
    chan_map = list(arby.channel_mapping)

    while len(upscales) < 4: upscales.append(array(upscales[0].typecode, [0]))
    while len(chan_map) < 4: chan_map.append(-1)

    if fast_dds_defs:
        a_scale, r_scale, g_scale, b_scale = upscales[: 4]
        dds_defs_ext.unpack_vu(
            unpacked, packed, a_scale, r_scale, g_scale, b_scale,
            ucc, array("b", chan_map[: 4]))
        return unpacked

    sign_mask = 1 << (bpc - 1)   # == 128 for 8bpc
    chan_mask = (1 << bpc) - 1   # == 255 for 8bpc
    dist_max  = (sign_mask - 1)  # == 127 for 8bpc
    dist_max_sq = dist_max**2    # == 16129 for 8bpc

    # convert to tuples for faster access
    upscales = tuple(tuple(scale) for scale in upscales)
    for i in range(0, len(packed)):
        # RGB normal maps use unsigned chars, which maps to:
        #     [0, 255] -> [-1, 1]
        # V8U8 uses signed chars, which maps(as unsigned chars) to:
        #     [0, 127] -> [+0, 1]    and    [128, 255] -> [-1, -0]
        # Ones compliment is used here to simplify math and to allow
        # all components to have a zero point and to make both sides
        # of the zero point have an equal numbers of points.

        off = ucc*i
        u = packed[i]&chan_mask
        v = (packed[i]>>bpc)&chan_mask
        if u&sign_mask: u -= chan_mask
        if v&sign_mask: v -= chan_mask

        # we're normalizing the coordinates here, not just unpacking them
        d = dist_max_sq - u**2 - v**2
        if d > 0:
            w = int(sqrt(d))
        else:
            n_len = sqrt(dist_max_sq - d)/dist_max
            u = int(u/n_len)
            v = int(v/n_len)
            w = 0

        colors = [0, u + sign_mask, v + sign_mask, w + sign_mask]
        dst_chan = 0
        for src_chan in chan_map:
            if src_chan < 0 and dst_chan == 0:
                # alpha and not reading alpha. set to full white
                unpacked[off] = unpack_max
            elif src_chan >= 0:
                unpacked[off + dst_chan] = upscales[dst_chan][colors[src_chan]]
            dst_chan += 1

    return unpacked


def unpack_r8g8(arby, bitmap_index, width, height, depth=1):
    return unpack_rg(arby, bitmap_index, width, height, depth, 8)


def unpack_r16g16(arby, bitmap_index, width, height, depth=1):
    return unpack_rg(arby, bitmap_index, width, height, depth, 16)


def unpack_rg(arby, bitmap_index, width, height, depth=1, bpc=8):
    packed = arby.texture_block[bitmap_index]

    #create a new array to hold the pixels after we unpack them
    unpack_code = arby._UNPACK_ARRAY_CODE
    unpack_size = ab.PIXEL_ENCODING_SIZES[unpack_code]
    unpack_max = (1<<(unpack_size*8)) - 1

    ucc = arby.unpacked_channel_count
    bytes_per_pixel = unpack_size*ucc
    unpacked = ab.bitmap_io.make_array(
        unpack_code, width*height, bytes_per_pixel)

    upscales = list(arby.channel_upscalers)
    chan_map = list(arby.channel_mapping)

    while len(upscales) < 4: upscales.append(array(upscales[0].typecode, [0]))
    while len(chan_map) < 4: chan_map.append(-1)

    if False and fast_dds_defs:
        # NOT IMPLEMENTED YET
        a_scale, r_scale, g_scale, b_scale = upscales[: 4]
        dds_defs_ext.unpack_gr(
            unpacked, packed, a_scale, r_scale, g_scale, b_scale,
            ucc, array("b", chan_map[: 4]))
        return unpacked

    sign_mask = 1 << (bpc - 1)   # == 128 for 8bpc
    chan_mask = (1 << bpc) - 1   # == 255 for 8bpc
    dist_max  = (sign_mask - 1)  # == 127 for 8bpc
    dist_max_sq = dist_max**2    # == 16129 for 8bpc

    # convert to tuples for faster access
    upscales = tuple(tuple(scale) for scale in upscales)
    for i in range(0, len(packed)):
        off = ucc*i
        u = ((packed[i]>>bpc)&chan_mask) - dist_max
        v = (packed[i]&chan_mask) - dist_max
        if u < 0: u += 1
        if v < 0: v += 1

        # we're normalizing the coordinates here, not just unpacking them
        d = dist_max_sq - u**2 - v**2
        if d > 0:
            w = int(sqrt(d))
        else:
            n_len = sqrt(dist_max_sq - d)/dist_max
            u = int(u/n_len)
            v = int(v/n_len)
            w = 0

        colors = [0, u + sign_mask, v + sign_mask, w + sign_mask]
        dst_chan = 0
        for src_chan in chan_map:
            if src_chan < 0 and dst_chan == 0:
                # alpha and not reading alpha. set to full white
                unpacked[off] = unpack_max
            elif src_chan >= 0:
                unpacked[off + dst_chan] = upscales[dst_chan][colors[src_chan]]
            dst_chan += 1

    return unpacked


def unpack_g8b8(arby, bitmap_index, width, height, depth=1):
    return unpack_gb(arby, bitmap_index, width, height, depth, 8)


def unpack_g16b16(arby, bitmap_index, width, height, depth=1):
    return unpack_gb(arby, bitmap_index, width, height, depth, 16)


def unpack_gb(arby, bitmap_index, width, height, depth=1, bpc=8):
    packed = arby.texture_block[bitmap_index]

    #create a new array to hold the pixels after we unpack them
    unpack_code = arby._UNPACK_ARRAY_CODE
    unpack_size = ab.PIXEL_ENCODING_SIZES[unpack_code]
    unpack_max = (1<<(unpack_size*8)) - 1

    ucc = arby.unpacked_channel_count
    bytes_per_pixel = unpack_size*ucc
    unpacked = ab.bitmap_io.make_array(
        unpack_code, width*height, bytes_per_pixel)

    upscales = list(arby.channel_upscalers)
    chan_map = list(arby.channel_mapping)

    while len(upscales) < 4: upscales.append(array(upscales[0].typecode, [0]))
    while len(chan_map) < 4: chan_map.append(-1)

    if False and fast_dds_defs:
        # NOT IMPLEMENTED YET
        a_scale, r_scale, g_scale, b_scale = upscales[: 4]
        dds_defs_ext.unpack_gb(
            unpacked, packed, a_scale, r_scale, g_scale, b_scale,
            ucc, array("b", chan_map[: 4]))
        return unpacked

    sign_mask = 1 << (bpc - 1)   # == 128 for 8bpc
    chan_mask = (1 << bpc) - 1   # == 255 for 8bpc
    dist_max  = (sign_mask - 1)  # == 127 for 8bpc
    dist_max_sq = dist_max**2    # == 16129 for 8bpc

    # convert to tuples for faster access
    upscales = tuple(tuple(scale) for scale in upscales)
    for i in range(0, len(packed)):
        off = ucc*i
        v = ((packed[i]>>bpc)&chan_mask) - dist_max
        w = (packed[i]&chan_mask) - dist_max
        if v < 0: v += 1
        if w < 0: w += 1

        # we're normalizing the coordinates here, not just unpacking them
        d = dist_max_sq - v**2 - w**2
        if d > 0:
            u = int(sqrt(d))
        else:
            n_len = sqrt(dist_max_sq - d)/dist_max
            v = int(v/n_len)
            w = int(w/n_len)
            u = 0

        colors = [0, u + sign_mask, v + sign_mask, w + sign_mask]
        dst_chan = 0
        for src_chan in chan_map:
            if src_chan < 0 and dst_chan == 0:
                # alpha and not reading alpha. set to full white
                unpacked[off] = unpack_max
            elif src_chan >= 0:
                unpacked[off + dst_chan] = upscales[dst_chan][colors[src_chan]]
            dst_chan += 1

    return unpacked


########################################
'''######## PACKING ROUTINES ########'''
########################################


def pack_dxt1(arby, unpacked, width, height, depth=1):
    ucc, bpt = arby.unpacked_channel_count, 8
    width, height, depth = ab.clip_dimensions(width, height, depth)
    dxt_width, dxt_height = clip_dxt_dimensions(width, height)
    texel_width, texel_height, _ = ab.clip_dimensions(dxt_width//4, dxt_height//4)
    pixels_per_texel = get_texel_pixel_count(width, height)
    channels_per_texel = ucc*pixels_per_texel
    can_have_alpha = arby.color_key_transparency
    a_cutoff = arby.one_bit_bias

    _, r_scale, g_scale, b_scale = arby.channel_downscalers
    repacked = ab.bitmap_io.make_array("I", texel_width*texel_height, bpt)
    unpacked = swizzle_dxt(unpacked, width, height * depth, ucc)

    if fast_dds_defs:
        dds_defs_ext.pack_dxt1(
            repacked, unpacked, r_scale, g_scale, b_scale,
            pixels_per_texel, can_have_alpha, a_cutoff)
        return repacked

    #this is the indexing for each pixel in each texel
    #values are multiplied by 4 to account for the channels
    pixel_indices = range(0, channels_per_texel, ucc)
    make_alpha = False
    c_2 = [0,0,0,0]
    c_3 = [0,0,0,0]

    #shorthand names
    rpa = repacked
    upa = unpacked

    # convert to tuples for faster access
    r_scale, g_scale, b_scale = tuple(r_scale), tuple(g_scale), tuple(b_scale)

    #loop for each texel
    for txl_i in range(0, len(repacked), 2):
        dist0 = dist1 = c_0i = c_1i = idx = 0

        pxl_i = (txl_i//2)*channels_per_texel
        r_pxl_i = pxl_i + 1
        g_pxl_i = pxl_i + 2
        b_pxl_i = pxl_i + 3

        # compare distance between all pixels and find the two furthest apart
        # (we are actually comparing the area of the distance as it's faster)
        for i in pixel_indices:
            r = upa[r_pxl_i + i]
            g = upa[g_pxl_i + i]
            b = upa[b_pxl_i + i]
            for j in pixel_indices:
                if j <= i: continue
                dist1 = ((r - upa[r_pxl_i + j])**2+
                         (g - upa[g_pxl_i + j])**2+
                         (b - upa[b_pxl_i + j])**2)
                if dist1 > dist0:
                    dist0 = dist1
                    c_0i = i
                    c_1i = j

        # store furthest apart colors for use
        c_0 = upa[pxl_i + c_0i: pxl_i + c_0i + 4]
        c_1 = upa[pxl_i + c_1i: pxl_i + c_1i + 4]

        # quantize the colors down to 16 bit color and repack
        color0 = ((((r_scale[c_0[1]]*31+15)//255)<<11) |
                  (((g_scale[c_0[2]]*63+31)//255)<<5) |
                  (b_scale[c_0[3]]*31+15)//255)
        color1 = ((((r_scale[c_1[1]]*31+15)//255)<<11) |
                  (((g_scale[c_1[2]]*63+31)//255)<<5) |
                  (b_scale[c_1[3]]*31+15)//255)

        # figure out if we are using color key transparency for this pixel
        #by seeing if any of the alpha values are below the cutoff bias
        if can_have_alpha:
            make_alpha = False
            for i in pixel_indices:
                if upa[pxl_i+i] < a_cutoff:
                    make_alpha = True
                    break

        if color0 == color1 and not make_alpha:
            rpa[txl_i] = (color1<<16) | color0
            continue

        # if the current color selection doesn't match what we want then
        # we reverse which color is which (if we are using transparency then
        # the first color as an integer must be smaller or equal to the second)
        if make_alpha == (color0 > color1):
            c_0, c_1 = c_1, c_0
            rpa[txl_i] = (color0<<16) | color1
        else:
            rpa[txl_i] = (color1<<16) | color0

        # calculate the intermediate colors
        #If the target format is DXT2/3/4/5 then no CK transparency is used.
        #CK mode will only be selected if both colors are the same.
        #If both colors are the same then it is fine to run non-CK
        #calculation on it since it will default to index zero.
        #That is why the DXT3/5 calculation is in this part only
        if rpa[txl_i]&65535 > rpa[txl_i]>>16:
            c_2[1] = (c_0[1]*2 + c_1[1])//3
            c_2[2] = (c_0[2]*2 + c_1[2])//3
            c_2[3] = (c_0[3]*2 + c_1[3])//3

            c_3[1] = (c_0[1] + c_1[1]*2)//3
            c_3[2] = (c_0[2] + c_1[2]*2)//3
            c_3[3] = (c_0[3] + c_1[3]*2)//3

            # calculate each pixel's closest match
            # and assign it the proper index
            for i in pixel_indices:
                r = upa[r_pxl_i+i]
                g = upa[g_pxl_i+i]
                b = upa[b_pxl_i+i]
                dists = ((r-c_0[1])**2 + (g-c_0[2])**2 + (b-c_0[3])**2,
                         (r-c_1[1])**2 + (g-c_1[2])**2 + (b-c_1[3])**2,
                         (r-c_2[1])**2 + (g-c_2[2])**2 + (b-c_2[3])**2,
                         (r-c_3[1])**2 + (g-c_3[2])**2 + (b-c_3[3])**2)

                idx += dists.index(min(dists))<<(i>>1)

            rpa[txl_i+1] = idx
            continue

        c_2[1] = (c_0[1]+c_1[1])//2
        c_2[2] = (c_0[2]+c_1[2])//2
        c_2[3] = (c_0[3]+c_1[3])//2
        #here, c_3 represents zero color and fully transparent

        #calculate each pixel's closest match and assign it the proper index
        for i in pixel_indices:
            if upa[pxl_i+i] < a_cutoff:
                idx += 3<<(i>>1)
                continue
            r = upa[r_pxl_i+i]
            g = upa[g_pxl_i+i]
            b = upa[b_pxl_i+i]
            dists = ((r-c_0[1])**2 + (g-c_0[2])**2 + (b-c_0[3])**2,
                     (r-c_1[1])**2 + (g-c_1[2])**2 + (b-c_1[3])**2,
                     (r-c_2[1])**2 + (g-c_2[2])**2 + (b-c_2[3])**2)

            idx += dists.index(min(dists))<<(i>>1)

        rpa[txl_i+1] = idx

    return repacked


def pack_dxt2_3(arby, unpacked, width, height, depth=1):
    ucc, bpt = arby.unpacked_channel_count, 16
    ucc = arby.unpacked_channel_count
    width, height, depth = ab.clip_dimensions(width, height, depth)
    dxt_width, dxt_height = clip_dxt_dimensions(width, height)
    texel_width, texel_height, _ = ab.clip_dimensions(dxt_width//4, dxt_height//4)
    pixels_per_texel = get_texel_pixel_count(width, height)
    channels_per_texel = ucc*pixels_per_texel

    a_scale, r_scale, g_scale, b_scale = arby.channel_downscalers
    repacked = ab.bitmap_io.make_array("I", texel_width*texel_height, bpt)
    unpacked = swizzle_dxt(unpacked, width, height * depth, ucc)

    if fast_dds_defs:
        dds_defs_ext.pack_dxt2_3(
            repacked, unpacked,
            a_scale, r_scale, g_scale, b_scale, pixels_per_texel)
        return repacked

    # convert to tuples for faster access
    a_scale, r_scale, g_scale, b_scale = tuple(a_scale), tuple(r_scale),\
                                         tuple(g_scale), tuple(b_scale)

    #this is the indexing for each pixel in each texel
    #values are multiplied by 4 to account for the channels
    pixel_indices = range(0, channels_per_texel, ucc)
    c_2 = [0,0,0,0]
    c_3 = [0,0,0,0]

    #shorthand names
    rpa = repacked
    upa = unpacked

    #loop for each texel
    for txl_i in range(0, len(repacked), 4):
        dist0 = dist1 = c_0i = c_1i = 0

        pxl_i = (txl_i//4)*channels_per_texel
        r_pxl_i = pxl_i + 1
        g_pxl_i = pxl_i + 2
        b_pxl_i = pxl_i + 3

        '''CALCULATE THE ALPHA'''
        # calculate alpha channel for DXT 2/3
        # coincidentally, the number of channels(4) matches the number of
        # bits in the alpha(4), so the shift is the same as the channel index
        alpha = sum(((a_scale[upa[pxl_i+i]]*15 + 7)//255) << i
                    for i in pixel_indices)

        rpa[txl_i]   = alpha&0xFFffFFff
        rpa[txl_i+1] = alpha>>32

        # CALCULATE THE COLORS
        # compare distance between all pixels and find the two furthest apart
        # (we are actually comparing the area of the distance as it's faster)
        for i in pixel_indices:
            r = upa[i + r_pxl_i]
            g = upa[i + g_pxl_i]
            b = upa[i + b_pxl_i]
            for j in pixel_indices:
                if j <= i: continue
                dist1 = ((r - upa[r_pxl_i + j])**2+
                         (g - upa[g_pxl_i + j])**2+
                         (b - upa[b_pxl_i + j])**2)
                if dist1 > dist0:
                    dist0 = dist1
                    c_0i = i
                    c_1i = j

        # store furthest apart colors for use
        c_0 = upa[pxl_i + c_0i: pxl_i + c_0i + 4]
        c_1 = upa[pxl_i + c_1i: pxl_i + c_1i + 4]

        # quantize the colors down to 16 bit color and repack
        color0 = ((((r_scale[c_0[1]]*31+15)//255)<<11) |
                  (((g_scale[c_0[2]]*63+31)//255)<<5) |
                  (b_scale[c_0[3]]*31+15)//255)
        color1 = ((((r_scale[c_1[1]]*31+15)//255)<<11) |
                  (((g_scale[c_1[2]]*63+31)//255)<<5) |
                  (b_scale[c_1[3]]*31+15)//255)

        if color0 != color1:
            # if the current color selection doesn't match what
            # we want then we reverse which color is which
            if color0 < color1:
                c_0, c_1 = c_1, c_0
                color0, color1 = color1, color0

            idx = 0
            c_2[1] = (c_0[1]*2 + c_1[1])//3
            c_2[2] = (c_0[2]*2 + c_1[2])//3
            c_2[3] = (c_0[3]*2 + c_1[3])//3

            c_3[1] = (c_0[1] + c_1[1]*2)//3
            c_3[2] = (c_0[2] + c_1[2]*2)//3
            c_3[3] = (c_0[3] + c_1[3]*2)//3

            # calculate each pixel's closest match
            # and assign it the proper index
            for i in pixel_indices:
                r = upa[r_pxl_i+i]
                g = upa[g_pxl_i+i]
                b = upa[b_pxl_i+i]
                dists = ((r-c_0[1])**2 + (g-c_0[2])**2 + (b-c_0[3])**2,
                         (r-c_1[1])**2 + (g-c_1[2])**2 + (b-c_1[3])**2,
                         (r-c_2[1])**2 + (g-c_2[2])**2 + (b-c_2[3])**2,
                         (r-c_3[1])**2 + (g-c_3[2])**2 + (b-c_3[3])**2)

                idx += dists.index(min(dists))<<(i>>1)

            rpa[txl_i+3] = idx
        rpa[txl_i+2] = (color1<<16) | color0

    return repacked


def pack_dxt4_5(arby, unpacked, width, height, depth=1):
    ucc, bpt = arby.unpacked_channel_count, 16
    ucc = arby.unpacked_channel_count
    width, height, depth = ab.clip_dimensions(width, height, depth)
    dxt_width, dxt_height = clip_dxt_dimensions(width, height)
    texel_width, texel_height, _ = ab.clip_dimensions(dxt_width//4, dxt_height//4)
    pixels_per_texel = get_texel_pixel_count(width, height)
    channels_per_texel = ucc*pixels_per_texel

    a_scale, r_scale, g_scale, b_scale = arby.channel_downscalers
    repacked = ab.bitmap_io.make_array("I", texel_width*texel_height, bpt)
    unpacked = swizzle_dxt(unpacked, width, height * depth, ucc)

    if fast_dds_defs:
        dds_defs_ext.pack_dxt4_5(
            repacked, unpacked, a_scale, r_scale, g_scale, b_scale,
            pixels_per_texel)
        return repacked

    # convert to tuples for faster access
    a_scale, r_scale, g_scale, b_scale = tuple(a_scale), tuple(r_scale),\
                                         tuple(g_scale), tuple(b_scale)

    #this is the indexing for each pixel in each texel
    #values are multiplied by 4 to account for the channels
    pixel_indices = range(0, channels_per_texel, ucc)
    c_0 = [0,0,0,0]
    c_1 = [0,0,0,0]
    c_2 = [0,0,0,0]
    c_3 = [0,0,0,0]

    #shorthand names
    rpa = repacked
    upa = unpacked

    #loop for each texel
    for txl_i in range(0, len(repacked), 4):
        dist0 = dist1 = c_0i = c_1i = alpha_idx = 0

        #cache so it doesn't have to keep being calculated
        pxl_i = (txl_i//4)*channels_per_texel
        r_pxl_i = pxl_i + 1
        g_pxl_i = pxl_i + 2
        b_pxl_i = pxl_i + 3

        # CALCULATE THE ALPHA
        #find the most extreme values
        alpha_vals = tuple(map(lambda i: a_scale[upa[pxl_i+i]], pixel_indices))
        alpha0 = max(alpha_vals)
        alpha1 = min(alpha_vals)

        if alpha0 == alpha1:
            # if they are the same number then
            # the indexing can stay at all zero
            pass
        elif alpha1 and alpha0 != 255:
            # if the most extreme values are NOT 0 or
            # 255, use them as the interpolation values
            # In this mode, value_0 must be greater than value_1
            alpha_dif = alpha0 - alpha1
            half_dif = alpha_dif//2
            # calculate and store which interpolated
            # index each alpha value is closest to
            for i in range(len(alpha_vals)):
                # 0 = c_0                 1 = c_1
                # 2 = (6*c_0 + c_1)//7    3 = (5*c_0 + 2*c_1)//7
                # 4 = (4*c_0 + 3*c_1)//7  5 = (3*c_0 + 4*c_1)//7
                # 6 = (2*c_0 + 5*c_1)//7  7 = (c_0 + 6*c_1)//7

                # calculate how far between both colors
                # that the value is as a 0 to 7 int
                tmp = ((alpha_vals[i] - alpha1)*7 + half_dif)//alpha_dif
                if tmp == 0:
                    alpha_idx |= 1<<(i*3)
                elif tmp < 7:
                    # Because the colors are stored in opposite
                    # order, we need to invert the index
                    alpha_idx |= (8-tmp)<<(i*3)
        else:
            # In this mode, value_0 must be less than or equal to value_1
            # if the most extreme values ARE 0 and 255 though, then
            # we need to calculate the second most extreme values
            alpha0 = 255
            alpha1 = 0
            for val in alpha_vals:
                # store if lowest int so far
                if val < alpha0 and val:        alpha0 = val
                # store if greatest int so far
                if val > alpha1 and val != 255: alpha1 = val

            if alpha1:
                alpha_dif = alpha1 - alpha0
            else:
                alpha0 = alpha_dif = 0
                alpha1 = 255

            half_dif = alpha_dif//2

            # calculate and store which interpolated
            # index each alpha value is closest to
            for i in range(len(alpha_vals)):
                # there are 4 interpolated colors in this mode
                # 0 =  c_0                1 = c_1
                # 2 = (4*c_0 + c_1)//5    3 = (3*c_0 + 2*c_1)//5
                # 4 = (2*c_0 + 3*c_1)//5  5 = (c_0 + 4*c_1)//5
                # 6 =  0                  7 = 255
                comp = alpha_vals[i]
                if comp == 0:
                    # if the value is 0 we set it to index 6
                    alpha_idx |= 6<<(i*3)
                elif comp == 255:
                    # if the value is 255 we set it to index 7
                    alpha_idx |= 7<<(i*3)
                elif alpha_dif:
                    # calculate how far between both colors
                    # that the value is as a 0 to 5 int
                    tmp = ((comp - alpha0)*5 + half_dif)//alpha_dif
                    if tmp == 5:
                        alpha_idx |= 1<<(i*3)
                    elif tmp > 0:
                        alpha_idx |= (tmp+1)<<(i*3)

        rpa[txl_i]   = ((alpha_idx<<16) + (alpha1<<8) + alpha0)&0xFFffFFff
        rpa[txl_i+1] = alpha_idx>>16

        # CALCULATE THE COLORS
        # compare distance between all pixels and find the two furthest apart
        # (we are actually comparing the area of the distance as it's faster)
        for i in pixel_indices:
            r = upa[r_pxl_i + i]
            g = upa[g_pxl_i + i]
            b = upa[b_pxl_i + i]
            for j in pixel_indices:
                if j <= i: continue
                dist1 = ((r - upa[r_pxl_i + j])**2+
                         (g - upa[g_pxl_i + j])**2+
                         (b - upa[b_pxl_i + j])**2)
                if dist1 > dist0:
                    dist0 = dist1
                    c_0i = i
                    c_1i = j

        # store furthest apart colors for use
        c_0 = upa[pxl_i + c_0i: pxl_i + c_0i + 4]
        c_1 = upa[pxl_i + c_1i: pxl_i + c_1i + 4]

        # quantize the colors down to 16 bit color and repack
        color0 = ((((r_scale[c_0[1]]*31+15)//255)<<11) |
                  (((g_scale[c_0[2]]*63+31)//255)<<5) |
                  (b_scale[c_0[3]]*31+15)//255)
        color1 = ((((r_scale[c_1[1]]*31+15)//255)<<11) |
                  (((g_scale[c_1[2]]*63+31)//255)<<5) |
                  (b_scale[c_1[3]]*31+15)//255)

        if color0 != color1:
            # if the current color selection doesn't match what
            # we want then we reverse which color is which
            if color0 < color1:
                c_0, c_1 = c_1, c_0
                color0, color1 = color1, color0

            idx = 0
            c_2[1] = (c_0[1]*2 + c_1[1])//3
            c_2[2] = (c_0[2]*2 + c_1[2])//3
            c_2[3] = (c_0[3]*2 + c_1[3])//3

            c_3[1] = (c_0[1] + c_1[1]*2)//3
            c_3[2] = (c_0[2] + c_1[2]*2)//3
            c_3[3] = (c_0[3] + c_1[3]*2)//3

            # calculate each pixel's closest match
            # and assign it the proper index
            for i in pixel_indices:
                r = upa[i + r_pxl_i]
                g = upa[i + g_pxl_i]
                b = upa[i + b_pxl_i]
                dists = ((r-c_0[1])**2 + (g-c_0[2])**2 + (b-c_0[3])**2,
                         (r-c_1[1])**2 + (g-c_1[2])**2 + (b-c_1[3])**2,
                         (r-c_2[1])**2 + (g-c_2[2])**2 + (b-c_2[3])**2,
                         (r-c_3[1])**2 + (g-c_3[2])**2 + (b-c_3[3])**2)

                idx += dists.index(min(dists))<<(i>>1)

            rpa[txl_i+3] = idx
        rpa[txl_i+2] = (color1<<16) | color0

    return repacked


def pack_dxt3a(arby, unpacked, width, height, depth=1):
    width, height, depth = ab.clip_dimensions(width, height, depth)
    #this is how many texels wide/tall the texture is
    dxt_width, dxt_height = clip_dxt_dimensions(width, height)
    texel_width, texel_height, _ = ab.clip_dimensions(dxt_width//4, dxt_height//4)

    #create a new array to hold the texels after we repack them
    ucc = arby.unpacked_channel_count
    assert arby.target_channel_count == ucc
    bpt = ucc*8

    scales = list(arby.channel_downscalers)
    repacked = ab.bitmap_io.make_array("I", texel_width*texel_height, bpt)
    unpacked = swizzle_dxt(unpacked, width, height * depth, ucc)

    pixels_per_texel   = get_texel_pixel_count(width, height)
    channels_per_texel = ucc*pixels_per_texel
    pixel_indices = range(0, channels_per_texel, ucc)

    if False and fast_dds_defs:
        # NOT IMPLEMENTED
        dds_defs_ext.pack_dxt3a(repacked, unpacked, pixels_per_texel, *scales)
        return repacked

    #shorthand names
    rpa = repacked
    upa = unpacked

    # convert to tuples for faster access
    for i in range(len(scales)):
        scales[i] = tuple(scales[i])

    #loop for each texel
    for txl_i in range(0, len(repacked), 2):
        #cache so it doesn't have to keep being calculated
        pxl_i = (txl_i//(2*ucc))*channels_per_texel
        chan = (txl_i//2)%ucc
        scale = scales[chan]

        # CALCULATE THE ALPHA
        alpha = a_shift = 0
        for i in pixel_indices:
            alpha |= ((scale[upa[pxl_i + i]]*15 + 7)//255) << a_shift
            a_shift += 4

        rpa[txl_i]   = alpha&0xFFffFFff
        rpa[txl_i+1] = alpha>>32

    return repacked


def pack_dxt5a(arby, unpacked, width, height, depth=1):
    width, height, depth = ab.clip_dimensions(width, height, depth)
    #this is how many texels wide/tall the texture is
    dxt_width, dxt_height = clip_dxt_dimensions(width, height)
    texel_width, texel_height, _ = ab.clip_dimensions(dxt_width//4, dxt_height//4)

    #create a new array to hold the texels after we repack them
    ucc = arby.unpacked_channel_count
    assert arby.target_channel_count == ucc
    bpt = ucc*8

    scales = list(arby.channel_downscalers)
    repacked = ab.bitmap_io.make_array("I", texel_width*texel_height, bpt)
    unpacked = swizzle_dxt(unpacked, width, height * depth, ucc)

    pixels_per_texel   = get_texel_pixel_count(width, height)
    channels_per_texel = ucc*pixels_per_texel
    pixel_indices = range(0, channels_per_texel, ucc)

    if False and fast_dds_defs:
        # NOT IMPLEMENTED
        dds_defs_ext.pack_dxt5a(repacked, unpacked, pixels_per_texel, *scales)
        return repacked

    #shorthand names
    rpa = repacked
    upa = unpacked

    # convert to tuples for faster access
    for i in range(len(scales)):
        scales[i] = tuple(scales[i])

    #loop for each texel
    for txl_i in range(0, len(repacked), 2):
        #cache so it doesn't have to keep being calculated
        pxl_i = (txl_i//(2*ucc))*channels_per_texel
        chan = (txl_i//2)%ucc
        idx = 0
        scale = scales[chan]

        vals = tuple(map(lambda i: scale[upa[pxl_i+i+chan]], pixel_indices))
        val0 = max(vals)
        val1 = min(vals)

        if val0 == val1:
            # if they are the same number then
            # the indexing can stay at all zero
            pass
        elif val1 and val0 != 255:
            # if the most extreme values are NOT 0 or
            # 255, use them as the interpolation values
            # In this mode, value_0 must be greater than value_1
            dif = val0 - val1
            half_dif = dif//2
            # calculate and store which interpolated
            # index each value is closest to
            for i in range(len(vals)):
                # 0 = c_0                 1 = c_1
                # 2 = (6*c_0 + c_1)//7    3 = (5*c_0 + 2*c_1)//7
                # 4 = (4*c_0 + 3*c_1)//7  5 = (3*c_0 + 4*c_1)//7
                # 6 = (2*c_0 + 5*c_1)//7  7 = (c_0 + 6*c_1)//7

                # calculate how far between both colors
                # that the value is as a 0 to 7 int
                tmp = ((vals[i] - val1)*7 + half_dif)//dif
                if tmp == 0:
                    idx |= 1<<(i*3)
                elif tmp < 7:
                    # Because the colors are stored in opposite
                    # order, we need to invert the index
                    idx |= (8-tmp)<<(i*3)
        else:
            # In this mode, value_0 must be less than or equal to value_1
            # if the most extreme values ARE 0 and 255 though, then
            # we need to calculate the second most extreme values
            val0 = 255
            val1 = 0
            for val in vals:
                # store if lowest int so far
                if val < val0 and val:        val0 = val
                # store if greatest int so far
                if val > val1 and val != 255: val1 = val

            if val1:
                dif = val1 - val0
            else:
                val0 = dif = 0
                val1 = 255

            half_dif = dif//2

            # calculate and store which interpolated
            # index each value is closest to
            for i in range(len(vals)):
                # there are 4 interpolated colors in this mode
                # 0 =  c_0                1 = c_1
                # 2 = (4*c_0 + c_1)//5    3 = (3*c_0 + 2*c_1)//5
                # 4 = (2*c_0 + 3*c_1)//5  5 = (c_0 + 4*c_1)//5
                # 6 =  0                  7 = 255
                comp = vals[i]
                if comp == 0:
                    # if the value is 0 we set it to index 6
                    idx |= 6<<(i*3)
                elif comp == 255:
                    # if the value is 255 we set it to index 7
                    idx |= 7<<(i*3)
                elif dif:
                    # calculate how far between both colors
                    # that the value is as a 0 to 5 int
                    tmp = ((comp - val0)*5 + half_dif)//dif
                    if tmp == 5:
                        idx |= 1<<(i*3)
                    elif tmp > 0:
                        idx |= (tmp+1)<<(i*3)

        rpa[txl_i] = ((idx<<16) | (val1<<8) | val0)&0xFFffFFff
        rpa[txl_i+1] = idx>>16

    return repacked


def pack_dxn(arby, unpacked, width, height, depth=1):
    width, height, depth = ab.clip_dimensions(width, height, depth)
    dxt_width, dxt_height = clip_dxt_dimensions(width, height)
    texel_width, texel_height, _ = ab.clip_dimensions(dxt_width//4, dxt_height//4)

    #create a new array to hold the texels after we repack them
    bpt = 16
    ucc = arby.unpacked_channel_count

    scales = list(arby.channel_downscalers)
    repacked = ab.bitmap_io.make_array("I", texel_width*texel_height, bpt)
    unpacked = swizzle_dxt(unpacked, width, height * depth, ucc)

    pixels_per_texel   = get_texel_pixel_count(width, height)
    channels_per_texel = ucc*pixels_per_texel
    pixel_indices = range(0, channels_per_texel, ucc)

    if False and fast_dds_defs:
        # NOT IMPLEMENTED
        dds_defs_ext.pack_dxn(repacked, unpacked, pixels_per_texel, *scales)
        return repacked

    #shorthand names
    rpa = repacked
    upa = unpacked

    # convert to tuples for faster access
    for i in range(len(scales)):
        scales[i] = tuple(scales[i])

    #loop for each texel
    for txl_i in range(0, len(repacked), 2):
        #cache so it doesn't have to keep being calculated
        pxl_i = (txl_i>>2)*channels_per_texel
        idx = 0

        # figure out if we're packing red or green(1=red, 2=green)
        chan = (((txl_i>>1)+1)%2)+1
        scale = scales[chan]

        vals = tuple(map(lambda i: scale[upa[pxl_i+i+chan]], pixel_indices))
        val0 = max(vals)
        val1 = min(vals)

        if val0 == val1:
            # if they are the same number then
            # the indexing can stay at all zero
            pass
        elif val1 and val0 != 255:
            # if the most extreme values are NOT 0 or
            # 255, use them as the interpolation values
            # In this mode, value_0 must be greater than value_1
            dif = val0 - val1
            half_dif = dif//2
            # calculate and store which interpolated
            # index each value is closest to
            for i in range(len(vals)):
                # 0 = c_0                 1 = c_1
                # 2 = (6*c_0 + c_1)//7    3 = (5*c_0 + 2*c_1)//7
                # 4 = (4*c_0 + 3*c_1)//7  5 = (3*c_0 + 4*c_1)//7
                # 6 = (2*c_0 + 5*c_1)//7  7 = (c_0 + 6*c_1)//7

                # calculate how far between both colors
                # that the value is as a 0 to 7 int
                tmp = ((vals[i] - val1)*7 + half_dif)//dif
                if tmp == 0:
                    idx |= 1<<(i*3)
                elif tmp < 7:
                    # Because the colors are stored in opposite
                    # order, we need to invert the index
                    idx |= (8-tmp)<<(i*3)
        else:
            # In this mode, value_0 must be less than or equal to value_1
            # if the most extreme values ARE 0 and 255 though, then
            # we need to calculate the second most extreme values
            val0 = 255
            val1 = 0
            for val in vals:
                # store if lowest int so far
                if val < val0 and val:        val0 = val
                # store if greatest int so far
                if val > val1 and val != 255: val1 = val

            if val1:
                dif = val1 - val0
            else:
                val0 = dif = 0
                val1 = 255

            half_dif = dif//2

            # calculate and store which interpolated
            # index each value is closest to
            for i in range(len(vals)):
                # there are 4 interpolated colors in this mode
                # 0 =  c_0                1 = c_1
                # 2 = (4*c_0 + c_1)//5    3 = (3*c_0 + 2*c_1)//5
                # 4 = (2*c_0 + 3*c_1)//5  5 = (c_0 + 4*c_1)//5
                # 6 =  0                  7 = 255
                comp = vals[i]
                if comp == 0:
                    # if the value is 0 we set it to index 6
                    idx |= 6<<(i*3)
                elif comp == 255:
                    # if the value is 255 we set it to index 7
                    idx |= 7<<(i*3)
                elif dif:
                    # calculate how far between both colors
                    # that the value is as a 0 to 5 int
                    tmp = ((comp - val0)*5 + half_dif)//dif
                    if tmp == 5:
                        idx |= 1<<(i*3)
                    elif tmp > 0:
                        idx |= (tmp+1)<<(i*3)

        rpa[txl_i]   = ((idx<<16) | (val1<<8) | val0)&0xFFffFFff
        rpa[txl_i+1] = idx>>16

    return repacked


def pack_ctx1(arby, unpacked, width, height, depth=1):
    width, height, depth = ab.clip_dimensions(width, height, depth)
    dxt_width, dxt_height = clip_dxt_dimensions(width, height)
    texel_width, texel_height, _ = ab.clip_dimensions(dxt_width//4, dxt_height//4)

    #create a new array to hold the texels after we repack them
    bpt = 8
    ucc = arby.unpacked_channel_count
    repacked = ab.bitmap_io.make_array("I", texel_width*texel_height, bpt)
    unpacked = swizzle_dxt(unpacked, width, height * depth, ucc)

    _, r_scale, g_scale, __ = arby.channel_downscalers

    pixels_per_texel   = get_texel_pixel_count(width, height)
    channels_per_texel = ucc*pixels_per_texel
    pixel_indices = range(0, channels_per_texel, ucc)

    if False and fast_dds_defs:
        # NOT IMPLEMENTED
        dds_defs_ext.pack_ctx1(repacked, unpacked, r_scale, g_scale,
                               pixels_per_texel)
        return repacked

    #shorthand names
    rpa = repacked
    upa = unpacked

    # convert to tuples for faster access
    r_scale, g_scale = tuple(r_scale), tuple(g_scale)

    #loop for each texel
    for txl_i in range(0, len(repacked), 2):
        dist0 = dist1 = c_0i = c_1i = idx = 0
        xy_0 = [0,0,0,0]
        xy_1 = [0,0,0,0]
        xy_2 = [0,0,0,0]
        xy_3 = [0,0,0,0]

        #cache so it doesn't have to keep being calculated
        pxl_i = (txl_i//2)*channels_per_texel
        r_pxl_i = pxl_i + 1
        g_pxl_i = pxl_i + 2

        # compare distance between all pixels and find the two furthest apart
        #(we are actually comparing the area of the distance as it's faster)
        for i in pixel_indices:
            for j in pixel_indices:
                if j <= i: continue
                dist1 = ((upa[r_pxl_i + i] - upa[r_pxl_i + j])**2 +
                         (upa[g_pxl_i + i] - upa[g_pxl_i + j])**2)
                if dist1 > dist0:
                    dist0 = dist1
                    c_0i = i
                    c_1i = j

        # store furthest apart colors for use
        xy_0[0] = r_scale[upa[r_pxl_i + c_0i]]
        xy_0[1] = g_scale[upa[g_pxl_i + c_0i]]

        xy_1[0] = r_scale[upa[r_pxl_i + c_1i]]
        xy_1[1] = g_scale[upa[g_pxl_i + c_1i]]

        color0 = xy_0[0] | (xy_0[1]<<8)
        color1 = xy_1[0] | (xy_1[1]<<8)

        rpa[txl_i] = color0 | (color1<<16)
        if color0 != color1:
            # calculate the intermediate colors
            xy_2[0] = (xy_0[0]*2 + xy_1[0])//3
            xy_2[1] = (xy_0[1]*2 + xy_1[1])//3

            xy_3[0] = (xy_0[0] + xy_1[0]*2)//3
            xy_3[1] = (xy_0[1] + xy_1[1]*2)//3

            # calculate each pixel's closest match
            # and assign it the proper index
            for i in pixel_indices:
                x = r_scale[upa[r_pxl_i + i]]
                y = g_scale[upa[g_pxl_i + i]]
                dist0 = (x-xy_0[0])**2 + (y-xy_0[1])**2
                dist1 = (x-xy_1[0])**2 + (y-xy_1[1])**2

                # add appropriate indexing value to array
                if dist0 <= dist1: #closer to color 0
                    if dist0 > (x-xy_2[0])**2 + (y-xy_2[1])**2:
                        #closest to color 2
                        idx |= 2<<(i//2)
                elif dist1 < (x-xy_3[0])**2 + (y-xy_3[1])**2:
                    #closest to color 1
                    idx |= 1<<(i//2)
                else: #closest to color 3
                    idx |= 3<<(i//2)

        rpa[txl_i+1] = idx

    return repacked


def pack_v8u8(arby, unpacked, width, height, depth=1):
    return pack_vu(arby, unpacked, width, height, depth, 8)


def pack_v16u16(arby, unpacked, width, height, depth=1):
    return pack_vu(arby, unpacked, width, height, depth, 16)


def pack_vu(arby, unpacked, width, height, depth=1, bpc=8):
    ucc = arby.unpacked_channel_count
    if ucc < 2:
        raise TypeError("Cannot convert image with less than 2 channels "
                        "to V%sU%s." % (bpc, bpc))

    bytes_per_pixel = (bpc * 2)//8
    typecode = ab.INVERSE_PIXEL_ENCODING_SIZES[bytes_per_pixel]
    packed = ab.bitmap_io.make_array(typecode, len(unpacked)//ucc)
    _, u_scale, v_scale, __ = arby.channel_downscalers
    if ucc == 2:
        chan0, chan1 = 0, 1
    else:
        chan0, chan1 = 1, 2

    if fast_dds_defs:
        dds_defs_ext.pack_vu(packed, unpacked, u_scale, v_scale,
                             ucc, chan0, chan1)
        return packed

    # convert to tuples for faster access
    u_scale, v_scale = tuple(u_scale), tuple(v_scale)

    sign_mask = 1 << (bpc - 1)
    sign_mask = sign_mask + (sign_mask << bpc)
    for i in range(0, len(unpacked), ucc):
        # RGB normal maps use unsigned chars, which maps to:
        #     [0, 255] -> [-1, 1]
        # V8U8 uses signed chars, which maps(as unsigned chars) to:
        #     [0, 127] -> [+0, 1]    and    [128, 255] -> [-1, -0]
        # Ones compliment is used here to simplify math and to allow
        # all components to have a zero point and to make both sides
        # of the zero point have an equal numbers of points.
        packed[i//ucc] = (((v_scale[unpacked[i + chan1]]<<bpc) |
                            u_scale[unpacked[i + chan0]])^sign_mask)

    return packed


def pack_r8g8(arby, unpacked, width, height, depth=1):
    return pack_rg(arby, unpacked, width, height, depth, 8)


def pack_r16g16(arby, unpacked, width, height, depth=1):
    return pack_rg(arby, unpacked, width, height, depth, 16)


def pack_rg(arby, unpacked, width, height, depth=1, bpc=8):
    ucc = arby.unpacked_channel_count
    if ucc < 2:
        raise TypeError("Cannot convert image with less than 2 channels "
                        "to R%sG%s." % (bpc, bpc))

    bytes_per_pixel = (bpc * 2)//8
    typecode = ab.INVERSE_PIXEL_ENCODING_SIZES[bytes_per_pixel]
    packed = ab.bitmap_io.make_array(typecode, len(unpacked)//ucc)
    _, r_scale, g_scale, __ = arby.channel_downscalers
    if ucc == 2:
        chan0, chan1 = 0, 1
    else:
        chan0, chan1 = 1, 2

    if False and fast_dds_defs:
        # NOT IMPLEMENTED YET
        dds_defs_ext.pack_rg(packed, unpacked, r_scale, g_scale,
                             ucc, chan0, chan1)
        return packed

    # convert to tuples for faster access
    r_scale, g_scale = tuple(r_scale), tuple(g_scale)
    for i in range(0, len(unpacked), ucc):
        packed[i//ucc] = ((g_scale[unpacked[i + chan1]]<<bpc) |
                           r_scale[unpacked[i + chan0]])

    return packed


def pack_g8b8(arby, unpacked, width, height, depth=1):
    return pack_gb(arby, unpacked, width, height, depth, 8)


def pack_g16b16(arby, unpacked, width, height, depth=1):
    return pack_gb(arby, unpacked, width, height, depth, 16)


def pack_gb(arby, unpacked, width, height, depth=1, bpc=8):
    ucc = arby.unpacked_channel_count
    if ucc < 2:
        raise TypeError("Cannot convert image with less than 2 channels "
                        "to G%sB%s." % (bpc, bpc))

    bytes_per_pixel = (bpc * 2)//8
    typecode = ab.INVERSE_PIXEL_ENCODING_SIZES[bytes_per_pixel]
    packed = ab.bitmap_io.make_array(typecode, len(unpacked)//ucc)
    _, __, g_scale, b_scale = arby.channel_downscalers
    if ucc == 2:
        chan0, chan1 = 0, 1
    else:
        chan0, chan1 = 1, 2

    if False and fast_dds_defs:
        # NOT IMPLEMENTED YET
        dds_defs_ext.pack_gb(packed, unpacked, g_scale, b_scale,
                             ucc, chan0, chan1)
        return packed

    # convert to tuples for faster access
    g_scale, b_scale = tuple(g_scale), tuple(b_scale)
    for i in range(0, len(unpacked), ucc):
        packed[i//ucc] = ((b_scale[unpacked[i + chan1]]<<bpc) |
                           g_scale[unpacked[i + chan0]])

    return packed
