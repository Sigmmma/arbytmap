import os
import time
import mmap

from struct import pack_into, unpack_from, unpack
from array import array
from traceback import format_exc
from constants import *
try:
    from supyr_struct.defs.bitmaps.tga import tga_def
    from supyr_struct.defs.bitmaps.dds import dds_def
    from supyr_struct.defs.util import fcc
except Exception:
    print("Bitmap IO was not loaded. Cannot load or save images.")
    tga_def = dds_def = None

try:
    try:
        from .ext import bitmap_io_ext
    except Exception:
        from ext import bitmap_io_ext
    fast_bitmap_io = True
except Exception:
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
                elif bitdepths == set((0, 4)):    fmt = ab.FORMAT_A4R4G4B4
            elif fmt_head.rgb_bitcount == 24:  fmt = ab.FORMAT_R8G8B8
            elif fmt_head.rgb_bitcount == 32:
                if   bitdepths == set((0, 8)): fmt = ab.FORMAT_X8R8G8B8
                elif bitdepths == set((8,  )): fmt = ab.FORMAT_A8R8G8B8
        elif fmt_flags.alpha_only: fmt = ab.FORMAT_A8
        elif fmt_flags.has_alpha:  fmt = ab.FORMAT_A8L8
        elif fmt_flags.luminance:  fmt = ab.FORMAT_L8
        elif fmt_flags.yuv_space:  fmt = ab.FORMAT_Y8U8V8
        elif fmt_flags.vu_space:
            if   fmt_head.rgb_bitcount == 16: fmt = ab.FORMAT_V8U8
            elif fmt_head.rgb_bitcount == 32: fmt = ab.FORMAT_V16U16

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
                dims = ab.clip_dimensions(
                    head.width>>m, head.height>>m, max(head.depth, 1)>>m, fmt)
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

    if bpp_test == 1 or head.image_type.format.enum_name == "bw_1_bit":
        #tex_info["format"] = ab.FORMAT_A1
        err += "Unable to load black and white 1-bit color Targa images."
    elif bpp_test == 8:
        tex_info["format"] = {
            True:  ab.FORMAT_A8,
            False: ab.FORMAT_L8}.get(alpha_depth == 8)
    elif bpp_test == 16:
        tex_info["format"] = {
            0: ab.FORMAT_R5G6B5,   1: ab.FORMAT_A1R5G5B5,
            4: ab.FORMAT_A4R4G4B4, 8: ab.FORMAT_A8L8,
            }.get(alpha_depth)
    elif bpp_test == 24:
        tex_info["format"] = ab.FORMAT_R8G8B8
    elif bpp_test == 32:
        tex_info["format"] = {
            True:  ab.FORMAT_A8R8G8B8,
            False: ab.FORMAT_X8R8G8B8}.get(alpha_depth == 8)

    if tex_info.get("format") is None:
        err += ("Unable to load %sbit color Targa images." % bpp)

    if image_desc.interleaving.data:
        err += "Unable to load Targa images with interleaved pixels."
        
    if err:
        print(err)
        return

    tex_block = []
    if image_desc.screen_origin.enum_name == "lower_left":
        tga_file.flip_image_origin()

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
        tex_block.append(array("B", pixels))
    elif bpp == 24:
        tex_block.append(pad_24bit_array(pixels))
    elif bpp == 48:
        tex_block.append(pad_48bit_array(pixels))
    else:
        pixel_array = array(typecode, pixels)

    convertor.load_new_texture(texture_block=tex_block,
                               texture_info=tex_info)


def save_to_rawdata_file(convertor, output_path, ext, **kwargs):
    """Saves the currently loaded texture to a raw file.
    The file has no header and in most cases wont be able
    to be directly opened be applications."""

    final_output_path = output_path

    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
        
    if convertor.is_palettized():
        print("Cannot save palettized images to RAW files.")
        return

    for sb in range(convertor.sub_bitmap_count):
        sb_output_path = output_path
        w, h, d = conv.width, conv.height, conv.depth
        if convertor.sub_bitmap_count > 1:
            sb_output_path = "%s_tex%s" % (sb_output_path, sb)
                
        #write each of the pixel arrays into the bitmap
        for m in range(convertor.mipmap_count + 1):
            i = m*convertor.sub_bitmap_count + sb
            
            final_output_path = sb_output_path
            if convertor.mipmap_count:
                final_output_path = "%s_mip%s" % (final_output_path, m)
        
            with open(final_output_path + ext, 'w+b') as raw_file:
                pixel_array = convertor.texture_block[i]
                if convertor.packed:
                    pixel_array = convertor.pack(pixel_array, w, h, d)
                    w, h, d = ab.clip_dimensions(w//2, h//2, d//2)
                    if pixel_array is None:
                        print("ERROR: UNABLE TO PACK IMAGE DATA.\n"+
                              "CANCELLING DDS SAVE.")
                        return()
                
                if ab.BITS_PER_PIXEL[convertor.format] == 24:
                    raw_file.write(unpad_24bit_array(pixel_array))
                elif ab.BITS_PER_PIXEL[convertor.format] == 48:
                    raw_file.write(unpad_48bit_array(pixel_array))
                else:
                    raw_file.write(pixel_array)


def save_to_dds_file(convertor, output_path, ext, **kwargs):
    """Saves the currently loaded texture to a DDS file"""
    typ = convertor.texture_type
    fmt = convertor.format

    if fmt in (ab.FORMAT_A16R16G16B16, ab.FORMAT_R16G16B16):
        print("ERROR: CANNOT SAVE '%s' FORMAT TO DDS.\nCANCELLING DDS SAVE." %
              fmt)
        return

    dds_file = dds_def.build()
    dds_file.data.pixel_data = b''
    dds_file.filepath = "%s.%s" % (output_path, ext)
    if not kwargs.get("overwrite", True) and os.path.exists(dds_file.filepath):
        return

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
    min_w, _, __ = ab.clip_dimensions(w, h, d, fmt)
    line_size = min_w * bpp

    # compute this for the block compressed formats
    head.pitch_or_linearsize = (line_size * 4) // (8 * 16)

    if fmt in (ab.FORMAT_DXT3A, ab.FORMAT_DXT3Y, ab.FORMAT_DXT3AY,
               ab.FORMAT_DXT5A, ab.FORMAT_DXT5Y, ab.FORMAT_DXT5AY,
               ab.FORMAT_CTX1):
        fmt_name = "LIN_%s" % fmt
        if fmt in (ab.FORMAT_DXT3AY, ab.FORMAT_DXT3Y, ab.FORMAT_DXT3A,
                   ab.FORMAT_DXT5AY, ab.FORMAT_DXT5Y, ab.FORMAT_DXT5A):
            fmt_name = fmt_name[:len("LIN_DXT_")] + "A"
            fmt_flags.luminance = "Y" in fmt
            fmt_flags.alpha_only = not fmt_flags.luminance
            fmt_flags.has_alpha  = "A" in fmt and not fmt_flags.alpha_only

        fmt_head.four_cc.set_to(fmt_name)
    elif "DXT" in fmt.upper() or fmt in (ab.FORMAT_DXN, ):
        fmt_head.four_cc.set_to(fmt)
    elif fmt in (ab.FORMAT_A1, ab.FORMAT_A4, ab.FORMAT_L4, ab.FORMAT_A4L4,
                 ab.FORMAT_L5V5U5, ab.FORMAT_X8L8V8U8, ab.FORMAT_Q8L8V8U8,
                 ab.FORMAT_Q8W8V8U8, ab.FORMAT_Q16W16V16U16,
                 ab.FORMAT_G8R8, ab.FORMAT_G16R16, ab.FORMAT_A2W10V10U10,
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
                     "neg_x", "neg_y", "neg_z")[: convertor.sub_bitmap_count]:
            head.caps2[name] = True

    head.mipmap_count = convertor.mipmap_count
    if convertor.mipmap_count:
        head.caps.complex = True
        head.caps.mipmaps = flags.mipmaps = True
    if convertor.photoshop_compatability:
        head.mipmap_count += 1
    dds_file.pprint(printout=True)
    #write each of the pixel arrays into the bitmap
    for sb in range(convertor.sub_bitmap_count):
        #write each of the pixel arrays into the bitmap
        for m in range(convertor.mipmap_count + 1):
            # get the index of the bitmap we'll be working with
            i = m*convertor.sub_bitmap_count + sb
            pixels = convertor.texture_block[i]

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
                return

            if bpp == 24:
                pixels = unpad_24bit_array(pixels)
            dds_file.data.pixel_data += pixels

        w, h, d = ab.clip_dimensions(w//2, h//2, d//2)

    dds_file.serialize(temp=False, backup=False, calc_pointers=False)


def save_to_tga_file(convertor, output_path, ext, **kwargs):
    """Saves the currently loaded texture to a TGA file"""
    conv = convertor
    fmt = conv.format
    
    if (fmt in (ab.FORMAT_R5G6B5, ab.FORMAT_A4R4G4B4,
                ab.FORMAT_A8L8, ab.FORMAT_V8U8) or fmt in ab.COMPRESSED_FORMATS):
        print("CANNOT EXTRACT THIS FORMAT TO TGA. EXTRACTING TO DDS INSTEAD.")
        save_to_dds_file(conv, output_path, "dds", **kwargs)
        return
    elif ab.BITS_PER_PIXEL[fmt] > 32:
        print("ERROR: CANNOT SAVE BITMAP OF HIGHER THAN 32 BIT "+
              "COLOR DEPTH TO DDS.\nCANCELLING TGA SAVE.")
        return

    channel_count = ab.CHANNEL_COUNTS[fmt]

    tga_file = tga_def.build()
    head = tga_file.data.header
    image_desc = head.image_descriptor

    head.width  = conv.width
    head.height = conv.height*conv.depth
    image_desc.screen_origin.set_to("upper_left")
    if conv.is_palettized():
        head.has_color_map.set_to("yes")
        head.image_type.format.set_to("color_mapped_rgb")
        head.color_map_length = 2**conv.indexing_size
        head.color_map_depth = ab.BITS_PER_PIXEL[fmt]

        head.bpp = 8
        if conv.target_indexing_size > 8:
            head.bpp = conv.indexing_size
    else:
        head.bpp = ab.BITS_PER_PIXEL[fmt]
        head.image_type.format.set_to("bw_8_bit")
        if channel_count > 1:
            image_desc.alpha_bit_count = CHANNEL_DEPTHS[fmt][0]
        if channel_count > 2:
            head.image_type.format.set_to("unmapped_rgb")

    final_output_path = output_path
    width = conv.width
    height = conv.height
    pals = conv.palette
    tex_block = conv.texture_block
    mip_count = convertor.mipmap_count+1
    overwrite = kwargs.get("overwrite", True)

    for sb in range(conv.sub_bitmap_count):
        if conv.sub_bitmap_count > 1:
            final_output_path = "%s_tex%s" % (output_path, sb)

        tga_file.filepath = "%s.%s" % (final_output_path, ext)
        if not overwrite and os.path.exists(tga_file.filepath):
            continue

        if conv.is_palettized():
            pal = pals[sb]
            idx = tex_block[sb]
            if conv.palette_packed:
                pal = conv.palette_packer(pal)
                
            '''need to pack the indexing and make sure it's 8-bit
               since TGA doesn't support less than 8 bit indexing'''

            if not conv.packed:
                idx = conv.indexing_packer(idx)
            elif conv.indexing_size < 8:
                temp = conv.target_indexing_size
                conv.target_indexing_size = 8
                idx = conv.indexing_packer(conv.indexing_unpacker(idx))
                conv.target_indexing_size = temp
                    
            tga_file.color_table = pal
            tga_file.pixels_wrapper.pixels = idx
        else:
            pixel_array = tex_block[sb]
            if not conv.packed:
                pixel_array = conv.pack(pixel_array, width, height, 0)
                width, height, _ = ab.clip_dimensions(width//2, height//2)
                if pixel_array is None:
                    print("ERROR: UNABLE TO PACK IMAGE DATA.\n"+
                          "CANCELLING TGA SAVE.")
                    return

            if ab.BITS_PER_PIXEL[fmt] == 24:
                pixel_array = unpad_24bit_array(pixel_array)
            
            tga_file.pixels_wrapper.pixels = pixel_array

        tga_file.serialize(temp=False, backup=False, calc_pointers=False)


def get_pixel_bytes_size(fmt, width, height, depth=1):
    pixel_size = ab.PIXEL_ENCODING_SIZES[ab.PACKED_TYPECODES[fmt]]

    #make sure the dimensions for the format are correct
    width, height, depth = ab.clip_dimensions(width, height, depth, fmt)
    
    bitmap_size = ab.pixel_count_to_array_length(
        height*width*depth, pixel_size, fmt)*pixel_size

    return bitmap_size


def make_array(typecode, item_ct, item_size=None, fill=0):
    # it would be nice to be able to make an array of w/e size
    # without having to create a bytearray first and throw it away
    if item_size is None:
        item_size = PIXEL_ENCODING_SIZES.get(typecode, 1)
    return array(typecode, bytes([fill])*item_ct*item_size)


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
        print("WARNING: PIXEL DATA SUPPLIED DID NOT MEET "+
              "THE SIZE EXPECTED. PADDING WITH ZEROS.")
        pixel_array.extend(
            make_array(pixel_array.typecode, bitmap_size - len(pixel_array), 1))
    
    #add the pixel array to the current texture block
    texture_block.append(pixel_array)
    return offset + bitmap_data_end


def bitmap_palette_to_array(rawdata, offset, palette_block, fmt, palette_count):
    return bitmap_bytes_to_array(rawdata, offset, palette_block,
                                 fmt, palette_count, 1)


def bitmap_indexing_to_array(rawdata, offset, indexing_block,
                             width, height, depth=1):
    """This function will create an array of pixels of width*height*depth from
       an iterable, sliceable, object. Since indexing never seems to be more
       than 8 bit, we won't worry about higher bit counts. Appends indexing
       array to supplied indexing_block and returns the end offset
    """
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
        padded = make_array("L", len(unpadded)//3)
        bitmap_io_ext.pad_24bit_array(padded, unpadded)
        return padded

    return array(
        "L", map(lambda x:(
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
    
    if padded.typecode == "L":
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
            "Bad typecode for padded 24bit array. Expected H or Q, got %s" %
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
