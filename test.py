import os, time
import cProfile
curr_dir = os.path.abspath(os.curdir)
start = time.time()

def convert_chain(last_fmt, formats):
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

        bitmap_test.load_from_file(
            input_path=curr_dir + "\\test_files\\disc_%s.dds" % last_fmt)
        bitmap_test.load_new_conversion_settings(target_format=fmt,
                                                 color_key_transparency=True,
                                                 **kw)
        start = time.time()
        bitmap_test.convert_texture()
        bitmap_test.save_to_file(
            output_path=curr_dir + "\\test_files\\disc_%s.dds" % fmt)
        last_fmt = fmt
        print("Completed %s in %s seconds" % (fmt, time.time()-start))

from traceback import format_exc
try:
    import arbytmap as ab
    print(("fast_arbytmap      = %s\n"
           "fast_raw_packers   = %s\n"
           "fast_raw_unpackers = %s\n"
           "fast_bitmap_io     = %s\n"
           "fast_swizzler      = %s\n"
           "fast_dds_defs      = %s\n")
          % (ab.fast_arbytmap, ab.fast_raw_packer, ab.fast_raw_unpacker,
             ab.fast_bitmap_io, ab.swizzler.fast_swizzler,
             ab.dds_defs.fast_dds_defs))

    #for fmt in sorted(ab.VALID_FORMATS):
    #    ab.print_format(fmt, True)
    #input()

    bitmap_test = ab.Arbytmap()
    bitmap_test.set_deep_color_mode(True)
    '''

    input('Press "Enter" to continue conversion 1')
    convert_chain(ab.FORMAT_A8R8G8B8, (
        ab.FORMAT_X8R8G8B8, ab.FORMAT_R8G8B8, ab.FORMAT_R5G6B5,
        ab.FORMAT_A1R5G5B5, ab.FORMAT_A4R4G4B4, ab.FORMAT_R3G3B2,
        ab.FORMAT_R8G8B8))
 
    input('Press "Enter" to continue conversion 2')
    convert_chain(ab.FORMAT_A8R8G8B8, (
        ab.FORMAT_A8L8, ab.FORMAT_A8R3G3B2, ab.FORMAT_A4R4G4B4,
        ab.FORMAT_A4L4, ab.FORMAT_A1R5G5B5, ab.FORMAT_L8))

    input('Press "Enter" to continue conversion 3')
    convert_chain(ab.FORMAT_A8R8G8B8, (
        ab.FORMAT_DXT5, ab.FORMAT_DXT4, ab.FORMAT_DXT3,
        ab.FORMAT_DXT2, ab.FORMAT_DXT1))
'''
    convert_chain(ab.FORMAT_A8R8G8B8, (ab.FORMAT_A16B16G16R16,))
    input('Press "Enter" to continue conversion 4')
    convert_chain(ab.FORMAT_A8R8G8B8, (
        ab.FORMAT_V16U16, ab.FORMAT_A16B16G16R16,
        ab.FORMAT_A2R10G10B10, ab.FORMAT_A2B10G10R10, ab.FORMAT_V8U8,
        ab.FORMAT_DXN, ab.FORMAT_CTX1, ab.FORMAT_R8G8B8))

    input('Press "Enter" to continue conversion 5')
    convert_chain(ab.FORMAT_A8L8, (ab.FORMAT_DXT5AY, ab.FORMAT_A4L4,
                                   ab.FORMAT_A8R8G8B8))

    input('Press "Enter" to continue conversion 6')
    convert_chain(ab.FORMAT_A8R8G8B8, (ab.FORMAT_L8, ab.FORMAT_DXT5Y,
                                       ab.FORMAT_L8))
except:
    print(format_exc())

input("Finished")
