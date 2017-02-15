import os, time
curr_dir = os.path.abspath(os.curdir)
start = time.time()

from traceback import format_exc
try:
    import arbytmap as ab

    bitmap_test = ab.Arbytmap()
    #bitmap_test.set_deep_color_mode(True)
    #for f in sorted(bc.VALID_FORMATS):
    #    bc.print_format(f)
    
    bitmap_test.load_from_file(input_path="dxt5.dds")
    bitmap_test.load_new_conversion_settings(target_format=ab.FORMAT_A8R8G8B8)
    bitmap_test.print_info(1,1,1,0,0)
    print('Press "Enter" to continue conversion'); input()
    start = time.time()
    bitmap_test.convert_texture()
    bitmap_test.save_to_file(output_path = curr_dir + "\\out.dds")
except:
    print(format_exc())

print("Completed in %s seconds" % (time.time()-start))
input()
