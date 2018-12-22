import os, time
import cProfile
from traceback import format_exc
import arbytmap as ab
curr_dir = os.path.abspath(os.curdir)
start = time.time()

bitmap_test = ab.Arbytmap()
def convert_chain(last_fmt, formats, **kwargs):
    ck_trans   = kwargs.pop("ck_trans", False)
    keep_alpha = kwargs.pop("keep_alpha", False)
    for fmt in formats:
        kw = {}
        src_chan_ct = ab.CHANNEL_COUNTS[last_fmt]
        tar_chan_ct = ab.CHANNEL_COUNTS[fmt]
        if src_chan_ct <= tar_chan_ct:
            if src_chan_ct == tar_chan_ct:
                pass
            elif tar_chan_ct == 2:
                if "A" in last_fmt:
                    kw["channel_mapping"] = ab.A_TO_AL
                else:
                    kw["channel_mapping"] = ab.L_TO_AL
            elif src_chan_ct == 2:
                kw["channel_mapping"] = ab.AL_TO_ARGB
            elif "A" in last_fmt:
                kw["channel_mapping"] = ab.A_TO_ARGB
            else:
                kw["channel_mapping"] = ab.L_TO_ARGB
        elif src_chan_ct == 4:
            if tar_chan_ct == 2:
                kw["channel_merge_mapping"] = ab.M_ARGB_TO_AL
            elif "A" in fmt:
                kw["channel_merge_mapping"] = ab.M_ARGB_TO_A
            else:
                kw["channel_merge_mapping"] = ab.M_ARGB_TO_L
        elif src_chan_ct == 2:
            if "A" in fmt:
                kw["channel_mapping"] = ab.ANYTHING_TO_A
            else:
                kw["channel_mapping"] = ab.AL_TO_L

        try:
            bitmap_test.load_from_file(
                input_path=curr_dir + "\\test_files\\disc_%s.dds" % last_fmt)
            bitmap_test.load_new_conversion_settings(
                target_format=fmt, color_key_transparency=ck_trans, **kw)
            start = time.time()
            bitmap_test.convert_texture()
            bitmap_test.save_to_file(
                output_path=curr_dir + "\\test_files\\disc_%s.dds" % fmt)
            bitmap_test.save_to_file(
                output_path=curr_dir + "\\test_files\\disc_%s.png" % fmt,
                keep_alpha=keep_alpha)
            last_fmt = fmt
            print("Completed %s in %.4f seconds" % (fmt, time.time()-start))
        except Exception:
            print(format_exc())
            print("Failed %s" % fmt)


def run_test(print_formats=False, deep_color=False,
             keep_alpha=False, ck_trans=False):
    print(("fast_arbytmap      = %s\n"
           "fast_raw_packers   = %s\n"
           "fast_raw_unpackers = %s\n"
           "fast_bitmap_io     = %s\n"
           "fast_swizzler      = %s\n"
           "fast_dds_defs      = %s\n"
           "fast_tiler         = %s\n")
          % (ab.fast_arbytmap, ab.fast_raw_packer, ab.fast_raw_unpacker,
             ab.bitmap_io.fast_bitmap_io, ab.swizzler.fast_swizzler,
             ab.dds_defs.fast_dds_defs, ab.tiler.fast_tiler))

    if print_formats:
        for fmt in sorted(ab.VALID_FORMATS):
            ab.print_format(fmt, True)

    bitmap_test.set_deep_color_mode(deep_color)
    '''
'''
    input('Press "Enter" to begin conversion 1')
    convert_chain(ab.FORMAT_A8R8G8B8, (
        ab.FORMAT_R8G8B8, ab.FORMAT_R5G6B5, ab.FORMAT_R3G3B2),
                  ck_trans=ck_trans, keep_alpha=keep_alpha)

    input('Press "Enter" to begin conversion 2')
    convert_chain(ab.FORMAT_A8R8G8B8, (
        ab.FORMAT_A1R5G5B5, ab.FORMAT_A4R4G4B4, ab.FORMAT_A8R3G3B2),
                  ck_trans=ck_trans, keep_alpha=keep_alpha)

    input('Press "Enter" to begin conversion 3')
    convert_chain(ab.FORMAT_A8R8G8B8, (
        ab.FORMAT_DXT5, ab.FORMAT_DXT4,
        ab.FORMAT_DXT3, ab.FORMAT_DXT2, ab.FORMAT_DXT1),
                  ck_trans=ck_trans, keep_alpha=keep_alpha)

    input('Press "Enter" to begin conversion 4')
    convert_chain(ab.FORMAT_A8R8G8B8, (
        ab.FORMAT_V16U16,
        ab.FORMAT_A16B16G16R16, ab.FORMAT_A2R10G10B10, ab.FORMAT_A2B10G10R10,
        ab.FORMAT_V8U8, ab.FORMAT_DXN, ab.FORMAT_CTX1),
                  ck_trans=ck_trans, keep_alpha=keep_alpha)

    input('Press "Enter" to begin conversion 5')
    convert_chain(ab.FORMAT_A8R8G8B8, (
        ab.FORMAT_A8L8,
        ab.FORMAT_A16, ab.FORMAT_A8,
        ab.FORMAT_DXT5A, ab.FORMAT_DXT3A,
        #ab.FORMAT_A4, ab.FORMAT_A2, ab.FORMAT_A1,
        ), ck_trans=ck_trans, keep_alpha=keep_alpha)

    input('Press "Enter" to begin conversion 6')
    convert_chain(ab.FORMAT_A8R8G8B8, (
        ab.FORMAT_L8,
        ab.FORMAT_DXT5Y, ab.FORMAT_DXT3Y,
        ), ck_trans=ck_trans, keep_alpha=keep_alpha)

    input('Press "Enter" to begin conversion 6')
    convert_chain(ab.FORMAT_A8R8G8B8, (ab.FORMAT_A16L16, ab.FORMAT_A8L8,
                                       ab.FORMAT_DXT5AY, ab.FORMAT_DXT3AY,
                                       ab.FORMAT_A4L4),
                  ck_trans=ck_trans, keep_alpha=keep_alpha)

    input('Press "Enter" to begin conversion 7')
    convert_chain(ab.FORMAT_A8R8G8B8, (ab.FORMAT_L8, ab.FORMAT_DXT5Y,
                                       ab.FORMAT_L8, ab.FORMAT_L16),
                  ck_trans=ck_trans, keep_alpha=keep_alpha)
    
if __name__ == "__main__":
    try:
        run_test(print_formats=False, deep_color=True,
                 keep_alpha=False, ck_trans=False)
    except:
        print(format_exc())

    input("Finished")
