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
    
    bitmap_test.load_from_file(input_path="DXT5_Alpha.dds")
    bitmap_test.load_new_conversion_settings(target_format=ab.FORMAT_A8Y8)
    print('Press "Enter" to continue conversion'); input()
    start = time.time()
    bitmap_test.convert_texture()
    bitmap_test.save_to_file(output_path = curr_dir + "\\Test.tga")
except:
    print(format_exc())

print("Completed in", str(time.time()-start).split('.')[0], "seconds.")
input()
