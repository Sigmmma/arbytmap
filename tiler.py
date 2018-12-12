from array import array
from math import ceil, log
from traceback import format_exc
import time

import arbytmap as ab

try:
    from arbytmap.ext import tiler_ext
    fast_tiler = True
except Exception:
    fast_tiler = False


def get_tiling_info(fmt):
    return (ab.format_defs.BLOCK_WIDTHS[fmt],
            ab.format_defs.BLOCK_HEIGHTS[fmt],
            ab.format_defs.BITS_PER_BLOCK[fmt] // 8)


class Tiler():
    converter = None
    tile_method = ""
    _methods = {}

    def __init__(self, **kwargs):
        self.converter = kwargs.get("converter")
        self.tile_method = kwargs.get("tile_method", "DEFAULT")

        if self.tile_method not in self._methods:
            raise TypeError("Unknown tiler method '%'" %
                            kwargs.get("tile_method"))

    def add_method(*args):
        if not args:
            return
        elif isinstance(args[0], Tiler):
            args = args[1:]

        tile_method, tile_func = args
        if tile_method in Tiler._methods:
            raise ValueError("Method '%s' already exists." % tile_method)
        Tiler._methods[tile_method] = tile_func

    def tile_texture(self, force=False, delete_old=True):
        conv = self.converter
        if conv is None:
            return False

        tex_block = conv.texture_block

        if tex_block is None:
            print("ERROR: NO TEXTURE LOADED. CANNOT PREFORM " +
                  "TILING OPERATION WITHOUT A LOADED TEXTURE")
            return False

        mode = conv.tile_mode
        if force:
            mode = not conv.tiled

        if mode == conv.tiled:
            return True

        # used to keep track of which pixel array we are reading
        width, height, depth = conv.width, conv.height, conv.depth
        if not conv.packed:
            raise TypeError("Cannot tile/untile unpacked texture.")

        i = 0
        for m in range(conv.mipmap_count + 1):
            m_width  = conv.get_packed_width(width, m)
            m_height = conv.get_packed_height(height, m)
            m_depth  = conv.get_packed_depth(depth, m)

            for s in range(conv.sub_bitmap_count):
                # get the pixel array to be tiled/untiled
                pixels = tex_block[i]

                # make the new array to place the tiled data into
                if isinstance(pixels, array):
                    tiled = array(pixels.typecode, pixels)
                elif isinstance(pixels, bytearray):
                    tiled = bytearray(len(pixels))
                else:
                    raise TypeError(
                        'Pixel array is not the proper type. Expected ' +
                        'array.array or bytearray, got %s' % type(pixels))

                self._tile_block(mode, pixels, tiled, conv.format,
                                 m_width, m_height, m_depth)

                # replace the old pixels with the new tiled one
                tex_block[i] = tiled

                # delete the old pixel array
                if delete_old:
                    # delete the old pixel array
                    del pixels[:]

                i += 1

        # now that we're done (un)tiling
        # the bitmap we invert the boolean
        conv.tiled = not conv.tiled

        # no errors occurred so we return a success
        return True

    def _tile_block(self, mode, pixels, modified_pixels,
                    fmt, width, height, depth):
        b_width, b_height, b_size = get_tiling_info(fmt)

        if mode:
            tiled, untiled = modified_pixels, pixels
        else:
            tiled, untiled = pixels, modified_pixels

        tiler = self._methods.get(self.tile_method)
        if fast_tiler:
            if mode:
                tiler_ext.tile_array()
            else:
                untiler_ext.tile_array()
            return
        elif tiler is None:
            raise TypeError("Unknown tile method '%s'" % self.tile_method)

        dst_b_size = b_size
        if isinstance(untiled, array):
            assert tiled.itemsize == untiled.itemsize
            dst_b_size = dst_b_size // untiled.itemsize

        x_chunks = width // b_width
        y_chunks = height // b_height
        #print("%sx%s  %sx%s  %sx%s %sx%s" % (width, height,
        #                                     x_chunks, y_chunks,
        #                                     b_width, b_height,
        #                                     b_size, untiled.itemsize))
        if mode:
            for i in range(y_chunks):
                offset = i * x_chunks
                for j in range(x_chunks):
                    x, y = tiler(offset + j, x_chunks, b_size)
                    t_idx = (i * x_chunks + j) * dst_b_size
                    u_idx = (y * x_chunks + x) * dst_b_size
                    tiled[t_idx: t_idx + dst_b_size] = untiled[
                        u_idx: u_idx + dst_b_size]
        else:
            for i in range(y_chunks):
                offset = i * x_chunks
                for j in range(x_chunks):
                    x, y = tiler(offset + j, x_chunks, b_size)
                    t_idx = (i * x_chunks + j) * dst_b_size
                    u_idx = (y * x_chunks + x) * dst_b_size
                    untiled[u_idx: u_idx + dst_b_size] = tiled[
                        t_idx: t_idx + dst_b_size]


def get_dxgi_tiled_address(offset, tile_ct_x, tile_len):
    aligned_width = (tile_ct_x + 31) & 0xFFffFFe0

    log_bpp = (tile_len >> 2) + ((tile_len >> 1) >> (tile_len >> 2))
    offset_b = offset << log_bpp
    offset_t = (((offset_b & 0xFFffF000) >> 3) +
               ((offset_b & 0x700) >> 2)) + (offset_b & 0x3F)
    offset_m = offset_t >> (7 + log_bpp)

    # x position values
    tile_x  = (((offset_t >> (5 + log_bpp)) & 2) + (offset_b >> 6)) & 3
    macro_x = (((offset_m % (aligned_width >> 5)) << 2) + tile_x) << 3
    micro_x = ((((offset_t >> 1) & 0xFFFFFFF0) +
               (offset_t & 0xF)) & ((tile_len << 3) - 1)) >> log_bpp

    # y position values
    tile_y  = ((offset_t >> (6 + log_bpp)) & 1) + ((offset_b & 0x800) >> 10)
    macro_y = (((offset_m // (aligned_width >> 5)) << 2) + tile_y) << 3
    micro_y = (((offset_t & (((tile_len << 6) - 1) & 0xFFffFFe0)) +
                ((offset_t & 15) << 1)) >> (3 + log_bpp)) & 0xFFffFFfe

    return macro_x + micro_x, macro_y + micro_y + ((offset_t & 0x10) >> 4)


Tiler.add_method("DEFAULT", get_dxgi_tiled_address)
Tiler.add_method("DXGI", get_dxgi_tiled_address)
