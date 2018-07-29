import os
import time
import mmap

from array import array
from copy import deepcopy
from math import ceil, log
from struct import pack_into, unpack_from, unpack
from traceback import format_exc
tga_def = dds_def = png_def = None

from arbytmap.constants import *
try:
    from supyr_struct.defs.bitmaps.dds import dds_def
    from supyr_struct.defs.bitmaps.tga import tga_def
    from supyr_struct.defs.bitmaps.objs.png import pad_idat_data
    from supyr_struct.defs.bitmaps.png import png_def
    from supyr_struct.defs.util import fcc
except Exception:
    print("SupyrStruct was not loaded. Cannot load or save non-raw images.")

try:
    from arbytmap.ext import bitmap_io_ext
    fast_bitmap_io = True
except:
    fast_bitmap_io = False


#this will be the reference to the bitmap convertor module.
#once the module loads this will become the reference to it.
ab = None


def load_from_dds_file(convertor, input_path, ext, **kwargs):
    """Loads a DDS file into the convertor."""
    dds_file = dds_def.build(filepath="%s.%s" % (input_path, ext))

    try:
        head = dds_file.data.header
        fmt_head  = head.dds_pixelformat
        fmt_flags = fmt_head.flags
        err = ""

        if fmt_head.four_cc.enum_name == "DX10":
            err += "CANNOT LOAD DX10 DDS FILES.\n"

        if head.caps2.volume and head.caps2.cubemap:
            err += ("ERROR: DDS HEADER INVALID. TEXTURE " +
                    "SPECIFIED AS BOTH CUBEMAP AND VOLUMETRIC.\n")

        mipmap_count = max(head.mipmap_count - 1, 0)
        typ = ab.TYPE_2D
        sub_bitmap_count = 1
        if head.caps2.volume:
            typ = ab.TYPE_3D
        elif head.caps2.cubemap:
            typ = ab.TYPE_CUBEMAP
            sub_bitmap_count = sum(bool(head.caps2[n]) for name in
                                   ("pos_x", "pos_y", "pos_z",
                                    "neg_x", "neg_y", "neg_z"))

        fmt = None
        bitdepths = set()
        for mask in (fmt_head.r_bitmask, fmt_head.g_bitmask,
                     fmt_head.b_bitmask, fmt_head.a_bitmask):
            bitdepths.add(sum((mask>>i)&1 for i in range(32)))

        if fmt_flags.four_cc:
            # the texture has a compression method
            fmt = fmt_head.four_cc.enum_name
            if fmt == "CxV8U8":
                fmt = ab.FORMAT_V8U8
            elif fmt.startswith("LIN_"):
                fmt = fmt.lstrip("LIN_")

            if   fmt == ab.FORMAT_DXT3A:
                if   fmt_flags.alpha_only: pass
                elif fmt_flags.has_alpha:  fmt = ab.FORMAT_DXT3AY
                elif fmt_flags.luminance:  fmt = ab.FORMAT_DXT3Y
            elif fmt == ab.FORMAT_DXT5A:
                if   fmt_flags.alpha_only: pass
                elif fmt_flags.has_alpha:  fmt = ab.FORMAT_DXT5AY
                elif fmt_flags.luminance:  fmt = ab.FORMAT_DXT5Y

            if fmt not in ab.VALID_FORMATS:
                fmt = None
        elif fmt_flags.rgb_space:
            if fmt_head.rgb_bitcount == 8:
                if   bitdepths == set((0, 2, 3)): fmt = ab.FORMAT_R3G3B2
            elif fmt_head.rgb_bitcount in (15, 16):
                if   bitdepths == set((2, 3, 8)): fmt = ab.FORMAT_A8R3G3B2
                elif bitdepths == set((0, 5, 6)): fmt = ab.FORMAT_R5G6B5
                elif bitdepths == set((1, 5)):    fmt = ab.FORMAT_A1R5G5B5
                elif bitdepths == set((0, 5)):    fmt = ab.FORMAT_A1R5G5B5
                elif bitdepths == set((4,  )):    fmt = ab.FORMAT_A4R4G4B4
            elif fmt_head.rgb_bitcount == 24:  fmt = ab.FORMAT_R8G8B8
            elif fmt_head.rgb_bitcount == 32:
                if   bitdepths == set((0, 8)): fmt = ab.FORMAT_X8R8G8B8
                elif bitdepths == set((8,  )): fmt = ab.FORMAT_A8R8G8B8
        elif fmt_flags.alpha_only:
            if   fmt_head.rgb_bitcount == 8:  fmt = ab.FORMAT_A8
            elif fmt_head.rgb_bitcount == 16: fmt = ab.FORMAT_A16
        elif fmt_flags.has_alpha:
            if   fmt_head.rgb_bitcount == 16: fmt = ab.FORMAT_A8L8
            elif fmt_head.rgb_bitcount == 32: fmt = ab.FORMAT_A16L16
        elif fmt_flags.luminance:
            if   fmt_head.rgb_bitcount == 8:  fmt = ab.FORMAT_L8
            elif fmt_head.rgb_bitcount == 16: fmt = ab.FORMAT_L16
        elif fmt_flags.vu_space:
            if   fmt_head.rgb_bitcount == 16: fmt = ab.FORMAT_V8U8
            elif fmt_head.rgb_bitcount == 32: fmt = ab.FORMAT_V16U16
        elif fmt_flags.yuv_space:  fmt = ab.FORMAT_Y8U8V8

        if fmt is None:
            err += "UNABLE TO DETERMINE DDS FORMAT. FAILED TO LOAD TEXTURE.\n"

        if err:
            print(err)
            return

        chan_ct = ab.CHANNEL_COUNTS[fmt]
        channel_map = None
        if not fmt_flags.four_cc and chan_ct > 2:
            channel_map = (fmt_head.a_bitmask, fmt_head.r_bitmask,
                           fmt_head.g_bitmask, fmt_head.b_bitmask)[: chan_ct]

        tex_info = {"width": head.width, "height": head.height,
                    "depth": max(head.depth, 1), "texture_type": typ,
                    "filepath": dds_file.filepath, "format": fmt,
                    "mipmap_count": mipmap_count,
                    "sub_bitmap_count": sub_bitmap_count,
                    "channel_mapping": channel_map}

        temp = []
        #loop over each mipmap and cube face
        #and turn them into pixel arrays
        dds_data, off = dds_file.data.pixel_data, 0
        for sb in range(sub_bitmap_count):
            for m in range(mipmap_count + 1):
                dims = ab.get_mipmap_dimensions(
                    head.width, head.height, head.depth, m)
                off = bitmap_bytes_to_array(dds_data, off, temp, fmt, *dims)

        # rearrange the images so they are sorted by [mip][bitmap]
        tex_block = [None]*len(temp)
        for m in range(mipmap_count + 1):
            for sb in range(sub_bitmap_count):
                tex_block[m*sub_bitmap_count + sb] = temp[
                    sb*(mipmap_count + 1) + m]
        convertor.load_new_texture(texture_block=tex_block,
                                   texture_info=tex_info)
    except:
        print(format_exc())


def load_from_tga_file(convertor, input_path, ext, **kwargs):
    """Loads a TGA file into the convertor."""
    tga_file = tga_def.build(filepath="%s.%s" % (input_path, ext))

    head = tga_file.data.header
    image_desc = head.image_descriptor
    pixels  = tga_file.data.pixels_wrapper.pixels
    palette = tga_file.data.color_table
    cm_bpp = head.color_map_depth
    bpp    = head.bpp
    alpha_depth = image_desc.alpha_bit_count
    if cm_bpp == 15: cm_bpp, alpha_depth = 16, 1
    if    bpp == 15: bpp,    alpha_depth = 16, 1

    #do another check to make sure image is color mapped
    color_mapped = head.image_type.format.enum_name == "color_mapped_rgb"

    tex_info = {"width":head.width, "height":head.height, "depth":1,
                "texture_type":"2D", "mipmap_count":0,
                "sub_bitmap_count":1, "filepath":input_path}

    err = ""
    #figure out what color format we've got
    bpp_test = cm_bpp if color_mapped else bpp

    fmt_name = head.image_type.format.enum_name
    fmt = None
    if bpp_test == 1 or fmt_name == "bw_1_bit":
        fmt = ab.FORMAT_A1
        err += "Unable to load black and white 1-bit color Targa images."
    elif bpp_test == 8:
        if fmt_name == "unmapped_rgb": pass
        elif alpha_depth == 0: fmt = ab.FORMAT_L8
        elif alpha_depth == 4: fmt = ab.FORMAT_A4L4
        elif alpha_depth == 8: fmt = ab.FORMAT_A8
    elif bpp_test == 16:
        if   alpha_depth == 8: fmt = ab.FORMAT_A8L8
        elif alpha_depth == 1: fmt = ab.FORMAT_A1R5G5B5
    elif bpp_test == 24:       fmt = ab.FORMAT_R8G8B8
    elif bpp_test == 32:
        if   alpha_depth == 0: fmt = ab.FORMAT_X8R8G8B8
        elif alpha_depth == 8: fmt = ab.FORMAT_A8R8G8B8

    if fmt not in ab.VALID_FORMATS:
        err += ("Unable to load %sbit color Targa images." % bpp)

    tex_info["format"] = fmt
    if image_desc.interleaving.data:
        err += "Unable to load Targa images with interleaved pixels."

    if err:
        print(err)
        return

    tex_block = []
    if image_desc.screen_origin.enum_name == "lower_left":
        tga_file.flip_image_origin()
        pixels = tga_file.data.pixels_wrapper.pixels

    typecode = ab.PACKED_TYPECODES[tex_info["format"]]
    if color_mapped:
        if cm_bpp == 24:
            palette = pad_24bit_array(palette)
        elif cm_bpp == 48:
            palette = pad_48bit_array(palette)
        else:
            palette = array(typecode, palette)

        #if the color map doesn't start at zero
        #then we need to shift around the palette
        if head.color_map_origin:
            palette = (palette[head.color_map_origin: ] +
                       palette[: head.color_map_origin])

        tex_info.update(palette=[palette], palettize=1, indexing_size=bpp)
        pixel_array = array("B", pixels)
    elif bpp == 24:
        pixel_array = pad_24bit_array(pixels)
    elif bpp == 48:
        pixel_array = pad_48bit_array(pixels)
    else:
        pixel_array = array(typecode, pixels)

    tex_block.append(pixel_array)
    convertor.load_new_texture(texture_block=tex_block,
                               texture_info=tex_info)


def save_to_rawdata_file(convertor, output_path, ext, **kwargs):
    """Saves the currently loaded texture to a raw file.
    The file has no header and in most cases wont be able
    to be directly opened be applications."""

    final_output_path = output_path
    if not os.path.exists(os.path.dirname(final_output_path)):
        os.makedirs(os.path.dirname(output_path))

    filenames = []
    if convertor.is_palettized():
        print("Cannot save palettized images to RAW files.")
        return filenames

    fmt = convertor.format
    tex_block = convertor.texture_block
    sub_bitmap_ct = convertor.sub_bitmap_count
    overwrite      = kwargs.get("overwrite", True)
    mip_levels     = kwargs.get("mip_levels", (0, ))
    bitmap_indexes = kwargs.get("bitmap_indexes", "all")

    if bitmap_indexes == "all":
        bitmap_indexes = range(sub_bitmap_ct)
    elif isinstance(bitmap_indexes, int):
        bitmap_indexes = (bitmap_indexes, )

    if mip_levels == "all":
        mip_levels = range(convertor.mipmap_count + 1)
    elif isinstance(mip_levels, int):
        mip_levels = (mip_levels, )

    for m in mip_levels:
        width  = max(convertor.width  // (1<<m), 1)
        height = max(convertor.height // (1<<m), 1)
        depth  = max(convertor.depth  // (1<<m), 1)

        mip_output_path = output_path
        if len(mip_levels) > 1:
            mip_output_path = "%s_mip%s" % (mip_output_path, m)

        for sb in bitmap_indexes:
            index = sb + m*sub_bitmap_ct
            if index >= len(tex_block):
                continue

            final_output_path = mip_output_path
            if sub_bitmap_ct > 1:
                final_output_path = "%s_tex%s" % (final_output_path, sb)

            final_output_path = "%s.%s" % (final_output_path, ext)

            if not overwrite and os.path.exists(final_output_path):
                continue

            with open(final_output_path, 'wb+') as raw_file:
                pixel_array = tex_block[index]
                if not convertor.packed:
                    pixel_array = convertor.pack(
                        pixel_array, width, height, depth)
                    if pixel_array is None:
                        print("ERROR: UNABLE TO PACK IMAGE DATA.\n")
                        continue

                if ab.BITS_PER_PIXEL[fmt] == 24:
                    raw_file.write(unpad_24bit_array(pixel_array))
                elif ab.BITS_PER_PIXEL[fmt] == 48:
                    raw_file.write(unpad_48bit_array(pixel_array))
                else:
                    raw_file.write(pixel_array)
            filenames.append(final_output_path)

    return filenames


def save_to_dds_file(convertor, output_path, ext, **kwargs):
    """Saves the currently loaded texture to a DDS file"""
    typ = convertor.texture_type
    fmt = convertor.format

    swizzle_mode = kwargs.pop("swizzle_mode", convertor.swizzled)
    channel_map = kwargs.pop("channel_mapping", None)
    if channel_map is not None or convertor.swizzled != swizzle_mode:
        conv_cpy = deepcopy(convertor)
        conv_cpy.load_new_conversion_settings(
            swizzle_mode=swizzle_mode, channel_mapping=channel_map)
        conv_cpy.convert_texture()
        return conv_cpy.save_to_file(output_path="%s.%s" % (output_path, ext),
                                     **kwargs)

    if fmt in (ab.FORMAT_A16R16G16B16, ab.FORMAT_R16G16B16):
        print("ERROR: CANNOT SAVE %s TO DDS.\nCANCELLING DDS SAVE." %
              fmt)
        return []

    dds_file = dds_def.build()
    dds_file.data.pixel_data = b''
    dds_file.filepath = "%s.%s" % (output_path, ext)
    if not kwargs.get("overwrite", True) and os.path.exists(dds_file.filepath):
        return []

    mip_levels     = kwargs.get("mip_levels", "all")
    bitmap_indexes = kwargs.get("bitmap_indexes", "all")

    if bitmap_indexes == "all":
        bitmap_indexes = range(convertor.sub_bitmap_count)
    elif isinstance(bitmap_indexes, int):
        bitmap_indexes = (bitmap_indexes, )

    if mip_levels == "all":
        mip_levels = range(convertor.mipmap_count + 1)
    elif isinstance(mip_levels, int):
        mip_levels = (mip_levels, )

    head = dds_file.data.header
    flags = head.flags
    fmt_head  = head.dds_pixelformat
    fmt_flags = fmt_head.flags

    w, h, d = convertor.width, convertor.height, convertor.depth
    bpp = ab.BITS_PER_PIXEL[fmt]
    masks = ab.CHANNEL_MASKS[fmt]
    offsets = ab.CHANNEL_OFFSETS[fmt]
    channel_count = ab.CHANNEL_COUNTS[fmt]
    if fmt in ab.THREE_CHANNEL_FORMATS:
        channel_count = 3

    palette_unpacker  = convertor.palette_unpacker
    indexing_unpacker = convertor.indexing_unpacker
    depalettizer      = convertor.depalettize_bitmap

    flags.linearsize = True
    fmt_flags.four_cc = True
    line_size = (w if w > 0 else 1) * bpp

    # compute this for the block compressed formats
    head.pitch_or_linearsize = (line_size * 4) // (8 * 16)

    if fmt in (ab.FORMAT_DXT3A, ab.FORMAT_DXT3Y, ab.FORMAT_DXT3AY,
               ab.FORMAT_DXT5A, ab.FORMAT_DXT5Y, ab.FORMAT_DXT5AY,
               ab.FORMAT_CTX1):
        fmt_name = "LIN_%s" % fmt
        if fmt in (ab.FORMAT_DXT3AY, ab.FORMAT_DXT3Y, ab.FORMAT_DXT3A,
                   ab.FORMAT_DXT5AY, ab.FORMAT_DXT5Y, ab.FORMAT_DXT5A):
            fmt_name = fmt_name[:len("LIN_DXT_")] + "A"
            fmt_flags.luminance  = "Y" in fmt
            fmt_flags.alpha_only = not fmt_flags.luminance
            fmt_flags.has_alpha  = "A" in fmt and not fmt_flags.alpha_only

        fmt_head.four_cc.set_to(fmt_name)
    elif "DXT" in fmt.upper() or fmt in (ab.FORMAT_DXN, ):
        fmt_head.four_cc.set_to(fmt)
    elif fmt in (ab.FORMAT_A1, ab.FORMAT_A4, ab.FORMAT_L4, ab.FORMAT_A4L4,
                 ab.FORMAT_L5V5U5, ab.FORMAT_X8L8V8U8, ab.FORMAT_Q8L8V8U8,
                 ab.FORMAT_Q8W8V8U8, ab.FORMAT_Q16W16V16U16,
                 ab.FORMAT_R8G8, ab.FORMAT_R16G16, ab.FORMAT_A2W10V10U10,
                 ab.FORMAT_A2B10G10R10, ab.FORMAT_A2R10G10B10,
                 ab.FORMAT_A16B16G16R16):
        flags.linearsize = False
        flags.pitch = True
        fmt_head.four_cc.set_to(fmt)
    else:
        # non-fourcc format
        flags.linearsize = False
        flags.pitch = True
        head.pitch_or_linearsize = (line_size + 7) // 8

        fmt_flags.four_cc = False
        fmt_head.rgb_bitcount = bpp
        if fmt in (ab.FORMAT_V8U8, ab.FORMAT_V16U16):
            fmt_flags.vu_space = True
            fmt_head.r_bitmask = masks[1] << offsets[1]
            fmt_head.g_bitmask = masks[2] << offsets[2]
        elif channel_count >= 3:
            if fmt == ab.FORMAT_Y8U8V8:
                fmt_flags.yuv_space = True
            else:
                fmt_flags.rgb_space = True
                fmt_flags.has_alpha = channel_count > 3
                if fmt_flags.has_alpha:
                    fmt_head.a_bitmask = masks[0] << offsets[0]

            fmt_head.r_bitmask = masks[1] << offsets[1]
            fmt_head.g_bitmask = masks[2] << offsets[2]
            fmt_head.b_bitmask = masks[3] << offsets[3]
        elif channel_count == 2:
            fmt_head.a_bitmask = masks[0] << offsets[0]
            fmt_head.r_bitmask = masks[1] << offsets[1]
            fmt_flags.has_alpha = True
            fmt_flags.luminance = True
        elif fmt == ab.FORMAT_A8:
            fmt_head.a_bitmask = masks[0] << offsets[0]
            fmt_flags.alpha_only = True
        else:
            fmt_head.r_bitmask = masks[0] << offsets[0]
            fmt_flags.luminance = True

    head.width, head.height, head.depth = w, h, d
    if typ == ab.TYPE_3D:
        head.caps.complex = True
        head.caps2.volume = flags.depth = True
    elif typ == ab.TYPE_CUBEMAP:
        head.caps.complex = True
        head.caps2.cubemap = True
        for name in ("pos_x", "pos_y", "pos_z",
                     "neg_x", "neg_y", "neg_z")[: len(bitmap_indexes)]:
            head.caps2[name] = True

    head.mipmap_count = len(mip_levels) - 1
    if convertor.mipmap_count:
        head.caps.complex = True
        head.caps.mipmaps = flags.mipmaps = True
    if convertor.photoshop_compatability:
        head.mipmap_count += 1

    #write each of the pixel arrays into the bitmap
    for sb in bitmap_indexes:
        #write each of the pixel arrays into the bitmap
        for m in mip_levels:
            # get the index of the bitmap we'll be working with
            i = m*convertor.sub_bitmap_count + sb
            pixels = convertor.texture_block[i]
            w, h, d = ab.get_mipmap_dimensions(
                convertor.width, convertor.height, convertor.depth, m)

            if convertor.is_palettized():
                pal = convertor.palette[i]
                if convertor.palette_packed:
                    pal = palette_unpacker(pal)
                if convertor.packed:
                    pixels = indexing_unpacker(pixels)

                pixels = convertor.pack_raw(depalettizer(pal, pixels))
            elif not convertor.packed:
                pixels = convertor.pack(pixels, w, h, d)

            if pixels is None:
                print("ERROR: UNABLE TO PACK IMAGE DATA.\nCANCELLING WRITE.")
                return []

            if bpp == 24:
                pixels = unpad_24bit_array(pixels)
            dds_file.data.pixel_data += pixels

    dds_file.serialize(temp=False, backup=False, calc_pointers=False)
    return [dds_file.filepath]


def save_to_tga_file(convertor, output_path, ext, **kwargs):
    """Saves the currently loaded texture to a TGA file"""
    fmt = convertor.format
    filenames = []

    make_copy = fmt not in (
        ab.FORMAT_A1, ab.FORMAT_L8, ab.FORMAT_A8,
        ab.FORMAT_A1R5G5B5, ab.FORMAT_R8G8B8,
        ab.FORMAT_X8R8G8B8, ab.FORMAT_A8R8G8B8)

    swizzle_mode = kwargs.pop("swizzle_mode", convertor.swizzled)
    if ("channel_mapping" in kwargs or (fmt != ab.FORMAT_A8R8G8B8 and make_copy)
            or convertor.swizzled != swizzle_mode):
        conv_cpy = deepcopy(convertor)
        # TODO: optimize this so the only textures loaded in and converted
        # are the mip_levels and sub_bitmaps that were requested to be saved
        conv_cpy.load_new_conversion_settings(
            target_format=ab.FORMAT_A8R8G8B8, swizzle_mode=swizzle_mode,
            channel_mapping=kwargs.pop("channel_mapping", None))
        conv_cpy.convert_texture()
        return conv_cpy.save_to_file(output_path="%s.%s" % (output_path, ext),
                                     **kwargs)

    channel_count = ab.CHANNEL_COUNTS[fmt]

    tga_file = tga_def.build()
    head = tga_file.data.header
    image_desc = head.image_descriptor
    image_desc.screen_origin.set_to("upper_left")
    if convertor.is_palettized():
        head.has_color_map.set_to("yes")
        head.image_type.format.set_to("color_mapped_rgb")
        head.color_map_length = 2**convertor.indexing_size
        head.color_map_depth = ab.BITS_PER_PIXEL[fmt]

        head.bpp = 8
        if convertor.target_indexing_size > 8:
            head.bpp = convertor.indexing_size
    else:
        head.bpp = ab.BITS_PER_PIXEL[fmt]
        head.image_type.format.set_to("bw_8_bit")
        if channel_count > 1:
            image_desc.alpha_bit_count = ab.CHANNEL_DEPTHS[fmt][0]
        if channel_count > 2:
            head.image_type.format.set_to("unmapped_rgb")

    final_output_path = output_path
    pals = convertor.palette
    tex_block = convertor.texture_block
    sub_bitmap_ct = convertor.sub_bitmap_count
    overwrite      = kwargs.get("overwrite", True)
    mip_levels     = kwargs.get("mip_levels", (0, ))
    bitmap_indexes = kwargs.get("bitmap_indexes", "all")

    if bitmap_indexes == "all":
        bitmap_indexes = range(sub_bitmap_ct)
    elif isinstance(bitmap_indexes, int):
        bitmap_indexes = (bitmap_indexes, )

    if mip_levels == "all":
        mip_levels = range(convertor.mipmap_count + 1)
    elif isinstance(mip_levels, int):
        mip_levels = (mip_levels, )

    for m in mip_levels:
        width  = max(convertor.width  // (1<<m), 1)
        height = max(convertor.height // (1<<m), 1)
        depth  = max(convertor.depth  // (1<<m), 1)
        head.width  = width
        head.height = height*depth
        mip_output_path = output_path
        if len(mip_levels) > 1:
            mip_output_path = "%s_mip%s" % (mip_output_path, m)

        for sb in bitmap_indexes:
            index = sb + m*sub_bitmap_ct
            if index >= len(tex_block):
                continue

            final_output_path = mip_output_path
            if sub_bitmap_ct > 1:
                final_output_path = "%s_tex%s" % (final_output_path, sb)

            tga_file.filepath = "%s.%s" % (final_output_path, ext)
            if not overwrite and os.path.exists(tga_file.filepath):
                continue

            if convertor.is_palettized():
                pal = pals[index]
                idx = tex_block[index]
                if not convertor.palette_packed:
                    pal = convertor.palette_packer(pal)

                '''need to pack the indexing and make sure it's 8-bit
                   since TGA doesn't support less than 8 bit indexing'''

                if not convertor.packed:
                    idx = convertor.indexing_packer(idx)
                elif convertor.indexing_size < 8:
                    temp = convertor.target_indexing_size
                    convertor.target_indexing_size = 8
                    try:
                        idx = convertor.indexing_packer(
                            convertor.indexing_unpacker(idx))
                    except Exception as e:
                        convertor.target_indexing_size = temp
                        raise e
                    finally:
                        convertor.target_indexing_size = temp

                tga_file.data.color_table = bytes(pal)
                tga_file.data.pixels_wrapper.pixels = idx
            else:
                pixel_array = tex_block[index]
                if not convertor.packed:
                    pixel_array = convertor.pack(pixel_array, width, height, 0)
                    if pixel_array is None:
                        print("ERROR: UNABLE TO PACK IMAGE DATA.\n"+
                              "CANCELLING TGA SAVE.")
                        return filenames

                if ab.BITS_PER_PIXEL[fmt] == 24:
                    pixel_array = unpad_24bit_array(pixel_array)

                tga_file.data.pixels_wrapper.pixels = pixel_array

            tga_file.serialize(temp=False, backup=False, calc_pointers=False)
            filenames.append(tga_file.filepath)

    return filenames


def save_to_png_file(convertor, output_path, ext, **kwargs):
    """Saves the currently loaded texture to a PNG file"""
    fmt = convertor.format
    palettized = convertor.is_palettized()
    if fmt in ab.THREE_CHANNEL_FORMATS:
        channel_count = 3
    else:
        channel_count = ab.CHANNEL_COUNTS[fmt]

    if channel_count <= 2:
        valid_depths = (1, 2, 4, 8, 16)
    else:
        valid_depths = (8, 16)

    bit_depth = fmt_depth = max(ab.CHANNEL_DEPTHS[fmt] + (1, ))
    target_depth = 1<<int(ceil(log(bit_depth, 2)))
    keep_alpha = kwargs.get("keep_alpha", channel_count <= 2)

    filenames = []

    if channel_count >= 2:
        # png doesnt allow 2 channel greyscale, so convert them to 4 channel
        if keep_alpha and channel_count == 4:
            if target_depth > 8: fmt_to_save_as = ab.FORMAT_A16R16G16B16
            else:                fmt_to_save_as = ab.FORMAT_A8R8G8B8
        else:
            if target_depth > 8: fmt_to_save_as = ab.FORMAT_R16G16B16
            else:                fmt_to_save_as = ab.FORMAT_R8G8B8
    elif fmt == ab.FORMAT_A16: fmt_to_save_as = fmt = ab.FORMAT_L16
    elif fmt == ab.FORMAT_A8:  fmt_to_save_as = fmt = ab.FORMAT_L8
    elif fmt == ab.FORMAT_A4:  fmt_to_save_as = fmt = ab.FORMAT_L4
    elif fmt == ab.FORMAT_A2:  fmt_to_save_as = fmt = ab.FORMAT_L2
    elif fmt == ab.FORMAT_A1:  fmt_to_save_as = fmt = ab.FORMAT_L1
    elif bit_depth > 8:      fmt_to_save_as = ab.FORMAT_L16
    elif target_depth == 8:  fmt_to_save_as = ab.FORMAT_L8
    elif target_depth == 4:  fmt_to_save_as = ab.FORMAT_L4
    elif target_depth == 2:  fmt_to_save_as = ab.FORMAT_L2
    elif target_depth == 1:  fmt_to_save_as = ab.FORMAT_L1

    swizzle_mode = kwargs.pop("swizzle_mode", convertor.swizzled)
    channel_map = kwargs.pop("channel_mapping", None)
    if (bit_depth not in valid_depths or channel_map is not None or
            fmt != fmt_to_save_as or convertor.swizzled != swizzle_mode):
        conv_cpy = deepcopy(convertor)
        # TODO: optimize this so the only textures loaded in and converted
        # are the mip_levels and sub_bitmaps that were requested to be saved
        if target_depth > 8:
            conv_cpy.set_deep_color_mode(True)

        conv_cpy.load_new_conversion_settings(
            target_format=fmt_to_save_as, swizzle_mode=convertor.swizzle_mode,
            channel_mapping=channel_map)
        if not conv_cpy.convert_texture():
            return []
        return conv_cpy.save_to_file(output_path="%s.%s" % (output_path, ext),
                                     **kwargs)

    png_file = png_def.build()
    png_file.data.chunks.append(case="IHDR")
    head = png_file.data.chunks[-1]
    plte_chunk = None
    if palettized:
        bit_depth = convertor.indexing_size
        color_type = "indexed_color"
        png_file.data.chunks.append(case="sRGB")
        png_file.data.chunks.append(case="PLTE")
        plte_chunk = png_file.data.chunks[-1]
        if channel_count == 4:
            png_file.data.chunks.append(case="tRNS")
            trns_chunk = png_file.data.chunks[-1]
    elif channel_count > 2:
        color_type = "truecolor"
        if channel_count == 4:
            color_type = "truecolor_with_alpha"

        png_file.data.chunks.append(case="sRGB")
    elif channel_count == 2:
        color_type = "greyscale_with_alpha"
    else:
        color_type = "greyscale"

    head.bit_depth = bit_depth
    head.color_type.set_to(color_type)

    png_file.data.chunks.append(case="IDAT")
    idat_chunk = png_file.data.chunks[-1]
    png_file.data.chunks.append(case="IEND")

    tex_block = convertor.texture_block
    pals = convertor.palette
    sub_bitmap_ct = convertor.sub_bitmap_count
    overwrite      = kwargs.get("overwrite", True)
    mip_levels     = kwargs.get("mip_levels", (0, ))
    bitmap_indexes = kwargs.get("bitmap_indexes", "all")
    png_compress_level = kwargs.get("png_compress_level", None)

    if bitmap_indexes == "all":
        bitmap_indexes = range(sub_bitmap_ct)
    elif isinstance(bitmap_indexes, int):
        bitmap_indexes = (bitmap_indexes, )

    if mip_levels == "all":
        mip_levels = range(convertor.mipmap_count + 1)
    elif isinstance(mip_levels, int):
        mip_levels = (mip_levels, )

    for m in mip_levels:
        width  = max(convertor.width  // (1<<m), 1)
        height = max(convertor.height // (1<<m), 1)
        depth  = max(convertor.depth  // (1<<m), 1)
        head.width  = width
        head.height = height*depth
        mip_output_path = output_path
        if len(mip_levels) > 1:
            mip_output_path = "%s_mip%s" % (mip_output_path, m)

        for sb in bitmap_indexes:
            index = sb + m*sub_bitmap_ct
            if index >= len(tex_block):
                continue

            final_output_path = mip_output_path
            if sub_bitmap_ct > 1:
                final_output_path = "%s_tex%s" % (final_output_path, sb)

            png_file.filepath = "%s.%s" % (final_output_path, ext)
            if not overwrite and os.path.exists(png_file.filepath):
                continue

            stride = width * head.bit_depth
            pix = tex_block[index]
            if palettized:
                pal = pals[index]
                if not convertor.packed:
                    pix = convertor.indexing_packer(pix)

                if channel_count == 4:
                    if convertor.palette_packed:
                        pal = convertor.palette_unpacker(pal)
                    alpha_pal = array(
                        "B", (pal[i] for i in range(0, len(pal), 4)))
                    trns_chunk.palette = alpha_pal
                    old_pal = pal
                    pal = bytearray(len(pal)*3//4)
                    for i in range(len(pal)//3):
                        j = i*4
                        i = i*3
                        pal[i: i+3] = old_pal[j+1: j+4]
                    plte_chunk.data = pal
                else:
                    if not convertor.palette_packed:
                        pal = convertor.palette_packer(pal)
                    plte_chunk.data = bytes(pal)

            else:
                stride *= channel_count
                if channel_count == 1 and fmt_depth == 16:
                    stride = stride//2

                if not convertor.packed:
                    pix = convertor.pack_raw(pix)

                if channel_count <= 2:
                    pix = bytearray(pix)
                    if fmt_depth == 16:
                        stride *= 2  # no idea why this is needed, but it is...
                        if channel_count != 1:
                            swap_channels(pix, (1, 0))
                elif channel_count == 4:
                    pix = bytearray(pix)
                    if   fmt_depth == 8:
                        swap_channels(pix, (2, 1, 0, 3))
                    elif fmt_depth == 16:
                        swap_channels(pix, (5, 4, 3, 2, 1, 0, 7, 6))
                elif fmt_depth == 8:
                    pix = bytearray(unpad_24bit_array(pix))
                    swap_channels(pix, (2, 1, 0))
                elif fmt_depth == 16:
                    pix = bytearray(unpad_48bit_array(pix))
                    swap_channels(pix, (5, 4, 3, 2, 1, 0))
            png_file.set_chunk_data(
                idat_chunk, pad_idat_data(pix, stride//8),
                png_compress_level=png_compress_level)
            png_file.serialize(temp=False, backup=False, calc_pointers=False)
            filenames.append(png_file.filepath)

    return filenames


def get_pixel_bytes_size(fmt, width, height, depth=1):
    if ab.PACKED_SIZE_CALCS.get(fmt):
        return ab.PACKED_SIZE_CALCS[fmt](fmt, width, height, depth)

    width, height, depth = ab.clip_dimensions(width, height, depth)
    return (ab.BITS_PER_PIXEL[fmt] * height * width * depth)//8


def make_array(typecode, item_ct, item_size=None, fill=0):
    # it would be nice to be able to make an array of w/e size
    # without having to create a bytearray first and throw it away
    if item_size is None:
        item_size = PIXEL_ENCODING_SIZES.get(typecode, 1)
    return array(typecode, bytes([fill])*item_ct*item_size)


def crop_pixel_data(pix, chan_ct, width, height, depth,
                    x0=0, new_width=-1, y0=0, new_height=-1, z0=0, new_depth=-1):
    if new_width  < 0: new_width  = width
    if new_height < 0: new_height = height
    if new_depth  < 0: new_depth  = depth

    new_pix = make_array(
        pix.typecode,
        new_width*new_height*new_depth*chan_ct,
        pix.itemsize)

    if len(pix) == 0:
        return new_pix

    src_x_skip0, src_y_skip0, src_z_skip0 = max(x0, 0), max(y0, 0), max(z0, 0)

    src_x_skip1 = max(width  - src_x_skip0 - new_width,  0)
    src_y_skip1 = max(height - src_y_skip0 - new_height, 0)
    src_z_skip1 = max(depth  - src_z_skip0 - new_depth,  0)

    x_stride = width  - src_x_skip0 - src_x_skip1
    y_stride = height - src_y_skip0 - src_y_skip1
    z_stride = depth  - src_z_skip0 - src_z_skip1

    if 0 in (x_stride, y_stride, z_stride):
        return new_pix

    dst_x_skip0, dst_y_skip0, dst_z_skip0 = max(-x0, 0), max(-y0, 0), max(-z0, 0)

    dst_x_skip1 = max(new_width  - dst_x_skip0 - x_stride, 0)
    dst_y_skip1 = max(new_height - dst_y_skip0 - y_stride, 0)

    src_z_skip0 *= width * height * chan_ct
    src_y_skip0 *= width * chan_ct
    src_y_skip1 *= width * chan_ct
    dst_z_skip0 *= new_width * new_height * chan_ct
    dst_y_skip0 *= new_width * chan_ct
    dst_y_skip1 *= new_width * chan_ct

    src_x_skip0 *= chan_ct
    dst_x_skip0 *= chan_ct
    src_x_skip1 *= chan_ct
    dst_x_skip1 *= chan_ct
    x_stride *= chan_ct
    
    if fast_bitmap_io:
        bitmap_io_ext.crop_pixel_data(
            pix, new_pix, z_stride, y_stride, x_stride,
            src_z_skip0, dst_z_skip0,
            src_y_skip0, dst_y_skip0, src_x_skip0, dst_x_skip0,
            src_y_skip1, dst_y_skip1, src_x_skip1, dst_x_skip1)
    else:
        src_i = src_z_skip0
        dst_i = dst_z_skip0
        src_x_skip1 += x_stride
        dst_x_skip1 += x_stride
        for z in range(z_stride):
            src_i += src_y_skip0
            dst_i += dst_y_skip0
            for y in range(y_stride):
                src_i += src_x_skip0
                dst_i += dst_x_skip0
                new_pix[dst_i: dst_i + x_stride] = pix[src_i: src_i + x_stride]
                src_i += src_x_skip1
                dst_i += dst_x_skip1

            src_i += src_y_skip1
            dst_i += dst_y_skip1

    return new_pix


def swap_channels(pix, channel_map):
    step = len(channel_map)
    channel_map = tuple(channel_map)
    src_map     = tuple(range(step))
    if channel_map == src_map:
        return

    assert set(channel_map) == set(src_map)

    if fast_bitmap_io:
        return bitmap_io_ext.swap_channels(pix, array("H", channel_map))

    for i in range(0, len(pix), step):
        orig = pix[i: i + step]
        for j in src_map:
            pix[i + j] = orig[channel_map[j]]


def bitmap_bytes_to_array(rawdata, offset, texture_block, fmt,
                          width, height, depth=1, bitmap_size=None, **kwargs):
    """This function will create an array of pixels of width*height*depth from
    an iterable, sliceable, object, and append it to the supplied texture_block.
    This function will return the offset of the end of the pixel data so that
    textures following the current one can be found."""
    #get the texture encoding
    encoding = ab.PACKED_TYPECODES[fmt]

    pixel_size = ab.PIXEL_ENCODING_SIZES[encoding]

    #get how many bytes the texture is going to be if it wasnt provided
    if bitmap_size is None:
        bitmap_size = bitmap_data_end = get_pixel_bytes_size(
            fmt, width, height, depth)
    bitmap_data_end = bitmap_size

    '''24 bit images are handled a bit differently since lots of
    things work on whole powers of 2. "2" can not be raised to an
    integer power to yield "24", whereas it can be for 8, 16, and 32.
    To fix this, the bitmap will be padded with an alpha channel on
    loading and ignored on saving. This will bring the 24 bit image
    up to 32 bit and make everything work just fine.'''
    if ab.BITS_PER_PIXEL[fmt] == 24:
        pixel_array = pad_24bit_array(rawdata[offset: offset + bitmap_size])
    elif ab.BITS_PER_PIXEL[fmt] == 48:
        pixel_array = pad_48bit_array(rawdata[offset: offset + bitmap_size])
    else:
        pixel_array = array(encoding, rawdata[offset: offset + bitmap_size])

    #if not enough pixel data was supplied, extra will be added
    if len(pixel_array)*pixel_size < bitmap_size:
        #print("WARNING: PIXEL DATA SUPPLIED DID NOT MEET "+
        #      "THE SIZE EXPECTED. PADDING WITH ZEROS.")
        pixel_array.extend(
            make_array(pixel_array.typecode,
                       bitmap_size - len(pixel_array)*pixel_size, 1))

    #add the pixel array to the current texture block
    texture_block.append(pixel_array)
    return offset + bitmap_data_end


def bitmap_palette_to_array(rawdata, offset, palette_block, fmt, palette_count):
    return bitmap_bytes_to_array(rawdata, offset, palette_block,
                                 fmt, palette_count, 1)


def bitmap_indexing_to_array(rawdata, offset, indexing_block,
                             width, height, depth=1):
    indexing_block.append(array("B", rawdata[offset:offset+width*height*depth]))
    return offset + width*height*depth


def pad_24bit_array(unpadded):
    if not hasattr(unpadded, 'typecode'):
        unpadded = array("B", unpadded)
    elif unpadded.typecode != 'B':
        raise TypeError(
            "Bad typecode for unpadded 24bit array. Expected B, got %s" %
            unpadded.typecode)

    if fast_bitmap_io:
        padded = make_array("I", len(unpadded)//3)
        bitmap_io_ext.pad_24bit_array(padded, unpadded)
        return padded

    return array(
        "I", map(lambda x:(
            unpadded[x] + (unpadded[x+1]<<8)+ (unpadded[x+2]<<16)),
                 range(0, len(unpadded), 3)))


def pad_48bit_array(unpadded):
    if not hasattr(unpadded, 'typecode'):
        unpadded = array("B", unpadded)
    elif unpadded.typecode != 'B':
        raise TypeError(
            "Bad typecode for unpadded 24bit array. Expected B, got %s" %
            unpadded.typecode)

    if fast_bitmap_io:
        padded = make_array("Q", len(unpadded)//3)
        bitmap_io_ext.pad_48bit_array(padded, unpadded)
        return padded

    return array(
        "Q", map(lambda x:(
            unpadded[x] + (unpadded[x+1]<<16)+ (unpadded[x+2]<<32)),
                 range(0, len(unpadded), 3)))


def unpad_24bit_array(padded):
    """given a 24BPP pixel data array that has been padded to
    32BPP, this will return an unpadded, unpacked, array copy.
    The endianness of the data will be little."""

    if padded.typecode == "I":
        # pixels have been packed
        unpadded = make_array("B", len(padded), 3)
        if fast_bitmap_io:
            bitmap_io_ext.unpad_24bit_array(unpadded, padded)
        else:
            for i in range(len(padded)):
                unpadded[i*3]   = padded[i]&255
                unpadded[i*3+1] = (padded[i]>>8)&255
                unpadded[i*3+2] = (padded[i]>>16)&255
    elif padded.typecode == "B":
        # pixels have NOT been packed

        # Because they havent been packed, it should be assumed
        # the channel order is the default one, namely ARGB.
        # Since we are removing the alpha channel, remove
        # the first byte from each pixel
        unpadded = make_array("B", len(padded)//4, 3)
        if fast_bitmap_io:
            bitmap_io_ext.unpad_24bit_array(unpadded, padded)
        else:
            for i in range(len(padded)//4):
                unpadded[i*3]   = padded[i*4+1]
                unpadded[i*3+1] = padded[i*4+2]
                unpadded[i*3+2] = padded[i*4+3]
    else:
        raise TypeError(
            "Bad typecode for padded 24bit array. Expected B or I, got %s" %
            padded.typecode)

    return unpadded


def unpad_48bit_array(padded):
    """given a 48BPP pixel data array that has been padded to
    64BPP, this will return an unpadded, unpacked, array copy.
    The endianness of the data will be little."""

    if padded.typecode == "Q":
        # pixels have been packed
        unpadded = make_array("H", len(padded), 6)
        if fast_bitmap_io:
            bitmap_io_ext.unpad_48bit_array(unpadded, padded)
        else:
            for i in range(len(padded)):
                j = i*3
                unpadded[j]   = padded[i]&65535
                unpadded[j+1] = (padded[i]>>16)&65535
                unpadded[j+2] = (padded[i]>>32)&65535
    elif padded.typecode == "H":
        # pixels have NOT been packed

        # Because they havent been packed, it should be assumed
        # the channel order is the default one, namely ARGB.
        # Since we are removing the alpha channel, remove
        # the first two bytes from each pixel
        unpadded = make_array("H", len(padded)//4, 6)
        if fast_bitmap_io:
            bitmap_io_ext.unpad_48bit_array(unpadded, padded)
        else:
            for i in range(len(padded)//4):
                j = i*3
                i *= 4
                unpadded[j]   = padded[i+1]
                unpadded[j+1] = padded[i+2]
                unpadded[j+2] = padded[i+3]
    else:
        raise TypeError(
            "Bad typecode for padded 48bit array. Expected H or Q, got %s" %
            padded.typecode)

    return unpadded


file_writers = {"raw":save_to_rawdata_file}
file_readers = {}

if tga_def is not None:
    file_writers["tga"] = save_to_tga_file
    file_readers["tga"] = load_from_tga_file

if dds_def is not None:
    file_writers["dds"] = save_to_dds_file
    file_readers["dds"] = load_from_dds_file

if png_def is not None:
    file_writers["png"] = save_to_png_file
