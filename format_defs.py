from array import array
from math import ceil, log
from traceback import format_exc

from arbytmap.constants import *

"""TEXTURE TYPES"""
TYPE_2D = "2D"
TYPE_3D = "3D"
TYPE_CUBEMAP = "CUBE"

"""TEXTURE FORMATS"""
FORMAT_A1 = "A1"  #NOT YET IMPLEMENTED
FORMAT_L1 = "L1"  #NOT YET IMPLEMENTED
FORMAT_A2 = "A2"  #NOT YET IMPLEMENTED
FORMAT_L2 = "L2"  #NOT YET IMPLEMENTED
FORMAT_A4 = "A4"  #NOT YET IMPLEMENTED
FORMAT_L4 = "L4"  #NOT YET IMPLEMENTED
FORMAT_A8 = "A8"
FORMAT_L8 = "L8"
FORMAT_AL8 = "AL8"
FORMAT_A16 = "A16"
FORMAT_L16 = "L16"
FORMAT_A2L2 = "A2L2"  #NOT YET IMPLEMENTED
FORMAT_A4L4 = "A4L4"
FORMAT_A8L8 = "A8L8"
FORMAT_A16L16 = "A16L16"
FORMAT_R3G3B2 = "R3G3B2"
FORMAT_R5G6B5 = "R5G6B5"
FORMAT_R8G8B8 = "R8G8B8"
FORMAT_A8R3G3B2 = "A8R3G3B2"
FORMAT_A4R4G4B4 = "A4R4G4B4"
FORMAT_A1R5G5B5 = "A1R5G5B5"
FORMAT_X8R8G8B8 = "X8R8G8B8"
FORMAT_A8R8G8B8 = "A8R8G8B8"

FORMAT_L5V5U5   = "L5V5U5"  # NOT PROPERLY IMPLEMENTED
FORMAT_Y8U8V8   = "Y8U8V8"  # NOT PROPERLY IMPLEMENTED
FORMAT_X8L8V8U8 = "X8L8V8U8"  # NOT PROPERLY IMPLEMENTED
FORMAT_Q8L8V8U8 = "Q8L8V8U8"  # NOT PROPERLY IMPLEMENTED

FORMAT_Q8W8V8U8 = "Q8W8V8U8"  # NOT PROPERLY IMPLEMENTED

#DEEP COLOR SUPPORT
FORMAT_A2B10G10R10  = "A2B10G10R10"
FORMAT_A2R10G10B10  = "A2R10G10B10"
FORMAT_A2W10V10U10  = "A2W10V10U10"  # NOT PROPERLY IMPLEMENTED

FORMAT_R16G16B16    = "R16G16B16"
FORMAT_A16B16G16R16 = "A16B16G16R16"
FORMAT_A16R16G16B16 = "A16R16G16B16"
FORMAT_Q16W16V16U16 = "Q16W16V16U16"  # NOT PROPERLY IMPLEMENTED


#FLOATING POINT SUPPORT
#NONE OF THESE ARE IMPLEMENTED YET
FORMAT_R16F = "R16F"                  #Will require Numpy
FORMAT_R32F = "R32F"
FORMAT_R16G16B16F    = "R16G16B16F"   #Will require Numpy
FORMAT_A16R16G16B16F = "A16R16G16B16F"
FORMAT_R32G32B32F    = "R32G32B32F"
FORMAT_A32R32G32B32F = "A32R32G32B32F"


DEFAULT_UNPACK_FORMAT = FORMAT_A8R8G8B8

VALID_FORMATS = set()
RAW_FORMATS = set()
COMPRESSED_FORMATS = set()
DDS_FORMATS = set()
THREE_CHANNEL_FORMATS = set()
PREMULTIPLIED_FORMATS = set()

UNPACKERS = {}
PACKERS = {}

BITS_PER_BLOCK = {}
BLOCK_WIDTHS = {}
BLOCK_HEIGHTS = {}

PACKED_WIDTH_CALCS = {}
PACKED_HEIGHT_CALCS = {}
PACKED_DEPTH_CALCS = {}
PACKED_SIZE_CALCS = {}
PACKED_FIELD_SIZES = {}

SUB_BITMAP_COUNTS = {TYPE_2D:1, TYPE_3D:1, TYPE_CUBEMAP:6}

# this is how many BITS(NOT BYTES) each format's pixels take up
BITS_PER_PIXEL = {}

# this is the typecode that each format's pixel data array will use
PACKED_TYPECODES = {}

# the number of channels possible in each format,
# regardless of whether or not they are raw formats
CHANNEL_COUNTS = {}

# this is how many bits the depth of each channel is for each raw format
CHANNEL_DEPTHS = {}

# this is the mask of each channel in each format
CHANNEL_MASKS = {}

# this is how far right the channel is shifted when
# unpacked and left when repacked
CHANNEL_OFFSETS = {}

ALL_FORMAT_COLLECTIONS = {
    "VALID_FORMAT":VALID_FORMATS, "BITS_PER_PIXEL":BITS_PER_PIXEL,
    "RAW_FORMAT":RAW_FORMATS, "THREE_CHANNEL_FORMAT":THREE_CHANNEL_FORMATS,
    "PREMULTIPLIED_FORMATS": PREMULTIPLIED_FORMATS,
    "COMPRESSED_FORMAT":COMPRESSED_FORMATS, "DDS_FORMAT":DDS_FORMATS,
    "PACKED_TYPECODES":PACKED_TYPECODES, "UNPACKER":UNPACKERS,
    "CHANNEL_COUNT":CHANNEL_COUNTS, "PACKER":PACKERS,
    "CHANNEL_OFFSETS":CHANNEL_OFFSETS, "PACKED_SIZE_CALCS": PACKED_SIZE_CALCS,
    "CHANNEL_MASKS":CHANNEL_MASKS, "CHANNEL_DEPTHS":CHANNEL_DEPTHS,
    "PACKED_WIDTH_CALCS": PACKED_WIDTH_CALCS,
    "PACKED_HEIGHT_CALCS": PACKED_HEIGHT_CALCS,
    "PACKED_DEPTH_CALCS": PACKED_DEPTH_CALCS,
    "PACKED_FIELD_SIZES": PACKED_FIELD_SIZES,
    "BITS_PER_BLOCK": BITS_PER_BLOCK,
    "BLOCK_WIDTHS": BLOCK_WIDTHS, "BLOCK_HEIGHTS": BLOCK_HEIGHTS,
    }

def register_format(format_id, silence_redefine=False, **kwargs):
    """
    Registers a texture format based on the provided information.
    """
    try:
        depths  = tuple(kwargs.get("depths", ()))
        offsets = tuple(kwargs.get("offsets", ()))
        masks   = tuple(kwargs.get("masks", ()))
        bpp     = kwargs.get("bpp", sum(depths))
        channel_count   = kwargs.get("channel_count")
        packed_typecode = kwargs.get("packed_typecode")

        if format_id is None:
            raise TypeError(
                "No identifier supplied for format.\n" +
                "This must be a hashable type, such as an int or str.")
        elif format_id in VALID_FORMATS:
            if silence_redefine:
                return

            raise TypeError((
                "Cannot add '%s' format definition to Arbytmap as " +
                "that format identifier is already in use.") % format_id)
        elif not bpp or not depths:
            raise TypeError((
                "Cannot define '%s' format without a given bits per " +
                "pixel or specifying each channels bit depths.") % format_id)
        elif max(depths) > 16:
            raise TypeError(
                "Cannot define '%s' as the max bits per channel is 16." %
                format_id)
        elif bpp > 64:
            raise TypeError(
                "Cannot define '%s' as the max bits per pixel is 64." %
                format_id)

        if not packed_typecode:
            packed_typecode = INVERSE_PIXEL_ENCODING_SIZES[(bpp + 7)//8]

        if not channel_count:
            channel_count = len(depths)

        if channel_count > 4:
            raise TypeError("Cannot create format '%s' with more than " +
                            "four channels." % format_id)

        if channel_count == 3:
            THREE_CHANNEL_FORMATS.add(format_id)
            channel_count = 4

        if not offsets:
            offsets = tuple(sum(depths[:i]) * (depths[i] != 0)
                            for i in range(len(depths)))

        if not masks:
            masks = tuple((2**depths[i] - 1) for i in range(len(depths)))

        if format_id in THREE_CHANNEL_FORMATS:
            if len(depths)  == 3: depths  += (0, )
            if len(offsets) == 3: offsets += (0, )
            if len(masks)   == 3: masks   += (0, )

        if kwargs.get("raw_format",    1): RAW_FORMATS.add(format_id)
        if kwargs.get("dds_format",    0): DDS_FORMATS.add(format_id)
        if kwargs.get("compressed",    0): COMPRESSED_FORMATS.add(format_id)
        if kwargs.get("premultiplied", 0): PREMULTIPLIED_FORMATS.add(format_id)

        if kwargs.get("unpacker"):
            UNPACKERS[format_id] = kwargs["unpacker"]
        if kwargs.get("packer"):
            PACKERS[format_id]   = kwargs["packer"]
        if kwargs.get("packed_size_calc"):
            PACKED_SIZE_CALCS[format_id] = kwargs["packed_size_calc"]

        BLOCK_WIDTHS[format_id]  = kwargs.get("block_width", 1)
        BLOCK_HEIGHTS[format_id] = kwargs.get("block_height", 1)
        bpb = kwargs.get("bpb", bpp * (BLOCK_WIDTHS[format_id] *
                                       BLOCK_HEIGHTS[format_id]))
            
        PACKED_WIDTH_CALCS[format_id] = kwargs.get("packed_width_calc",
                                                   packed_dimension_calc)
        PACKED_HEIGHT_CALCS[format_id] = kwargs.get("packed_height_calc",
                                                    packed_dimension_calc)
        PACKED_DEPTH_CALCS[format_id] = kwargs.get("packed_depth_calc",
                                                   packed_dimension_calc)
        PACKED_FIELD_SIZES[format_id] = kwargs.get("packed_field_sizes",
                                                   (bpp // 8, ))


        BITS_PER_PIXEL[format_id]   = bpp
        BITS_PER_BLOCK[format_id]   = bpb
        PACKED_TYPECODES[format_id] = packed_typecode
        CHANNEL_COUNTS[format_id]   = channel_count
        # we unpack to ARGB, not BGRA. Reverse the tuples
        CHANNEL_MASKS[format_id]    = tuple(masks[::-1])
        CHANNEL_DEPTHS[format_id]   = tuple(depths[::-1])
        CHANNEL_OFFSETS[format_id]  = tuple(offsets[::-1])

        VALID_FORMATS.add(format_id)
    except:
        print("Error occurred while trying to define new texture format.")
        print(format_exc())


def unregister_format(format_id):
    """
    Removed the specified format from all collections that define it.
    """
    for val in ALL_FORMAT_COLLECTIONS.values():
        if format_id not in val:
            continue
        elif isinstance(val, dict):
            del val[format_id]
        else:
            val.pop(val.index(format_id))


def print_format(format_id, printout=True):
    out_str = "%s Format Definition:\n" % format_id
    for key in sorted(ALL_FORMAT_COLLECTIONS.keys()):
        val = ALL_FORMAT_COLLECTIONS[key]
        if not isinstance(val, dict):
            out_str += '    %s: %s\n' % (key, format_id in val)
        elif format_id in val:
            out_str += '    %s: %s\n' % (key, val[format_id])

    out_str += '\n'
    if printout:
        print(out_str)
    return out_str


def get_mipmap_dimensions(width, height, depth, mip):
    return clip_dimensions(width>>mip, height>>mip, depth>>mip)


def packed_dimension_calc(dim, mip_level, tiled):
    dim = max(1, dim >> mip_level)
    if tiled:
        dim = 2**int(ceil(log(dim, 2.0)))
    return dim


def clip_dimensions(width, height, depth=1):
    return max(1, width), max(1, height), max(1, depth)


# Need to implement unpacking more than 1 pixel
# per byte for the formats i've commented out.
#register_format(FORMAT_A1, 1, depths=(1,))
#register_format(FORMAT_L1, 1, depths=(1,))
#register_format(FORMAT_A2, 1, depths=(2,))
#register_format(FORMAT_L2, 1, depths=(2,))
#register_format(FORMAT_A4, 1, depths=(4,))
#register_format(FORMAT_L4, 1, depths=(4,))
register_format(FORMAT_A8, 1, depths=(8,))
register_format(FORMAT_L8, 1, depths=(8,))
register_format(FORMAT_AL8, 1, depths=(8,8),
                offsets=(0,0), masks=(255, 255), bpp=8)
register_format(FORMAT_A16, 1, depths=(16,))
register_format(FORMAT_L16, 1, depths=(16,))
#register_format(FORMAT_A2L2, 1, depths=(2,2))
register_format(FORMAT_A4L4, 1, depths=(4,4))
register_format(FORMAT_A8L8, 1, depths=(8,8))
register_format(FORMAT_A16L16, 1, depths=(16,16))

register_format(FORMAT_R3G3B2,   1, depths=(2,3,3))
register_format(FORMAT_R5G6B5,   1, depths=(5,6,5))
register_format(FORMAT_R8G8B8,   1, depths=(8,8,8))
register_format(FORMAT_A1R5G5B5, 1, depths=(5,5,5,1))
register_format(FORMAT_A4R4G4B4, 1, depths=(4,4,4,4))
register_format(FORMAT_A8R3G3B2, 1, depths=(2,3,3,8))
register_format(FORMAT_X8R8G8B8, 1, depths=(8,8,8), bpp=32)
register_format(FORMAT_A8R8G8B8, 1, depths=(8,8,8,8))
register_format(FORMAT_R16G16B16, 1, depths=(16,16,16))
register_format(FORMAT_A16R16G16B16, 1, depths=(16,16,16,16))
register_format(FORMAT_A16B16G16R16, 1, depths=(16,16,16,16),
                offsets=(32,16,0,48), masks=(65535, 65535, 65535, 65535))
register_format(FORMAT_A2R10G10B10,  1, depths=(10,10,10,2))
register_format(FORMAT_A2B10G10R10,  1, depths=(10,10,10,2),
                offsets=(20,10,0,30), masks=(1023, 1023, 1023, 3))

# will need a converter to go between RGB and YUV color spaces
register_format(FORMAT_Y8U8V8, 1, depths=(8,8,8))

# will need a converter to go between RGB and LUV color spaces
register_format(FORMAT_L5V5U5,   1, depths=(5,5,5))
register_format(FORMAT_X8L8V8U8, 1, depths=(8,8,8), bpp=32)
register_format(FORMAT_Q8L8V8U8, 1, depths=(8,8,8,8))

register_format(FORMAT_Q8W8V8U8, 1, depths=(8,8,8,8))

register_format(FORMAT_A2W10V10U10,  1, depths=(10,10,10,2))
register_format(FORMAT_Q16W16V16U16, 1, depths=(16,16,16,16))
