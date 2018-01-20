"""COLOR CHANNEL ORDERS"""
#these are different orders that the color channels can be in.
#rather than make multiple types for each format with 3 or more
#channels, we define the different ways the channels can be ordered
"""channel values are in little endian: least significant byte first.
   So for example, in C_ORDER_BGRA, B is the first byte, and if the
   pixel were read as a 32 bit integer, B's bits would be values 0-255
"""
C_ORDER_ARGB = "ARGB"  # <---THE ORDER Arbytmap UNPACKS TO
C_ORDER_ABGR = "ABGR"
C_ORDER_RGBA = "RGBA"
C_ORDER_BGRA = "BGRA"  # <---DEFAULT MOST IMAGE FORMATS PACK TO

#Default pixel storing order for most image formats is little endian BGRA
'''The format Arbytmap will unpack pixels to will ALWAYS be ARGB.
   This does NOT apply to 1 and 2 channel formats though, so they are
   instead always unpacked to AY or A.
   In the future I intend to make it so you can specify which format
   to load from/save to, but the internal format will always be ARGB.
'''
C_ORDER_DEFAULT = C_ORDER_ARGB

PIXEL_ENCODING_SIZES = {"B":1, "H":2, "I":4, "Q":8, "b":1, "h":2, "i":4, "q":8}
INVERSE_PIXEL_ENCODING_SIZES = {
    0:"B", 1:"B",
    2:"H",
    3:"I", 4:"I",
    5:"Q", 6:"Q", 7:"Q", 8:"Q"}

# if a channel has this in it's divisor it
# will be erased when the bitmap is repacked.
# use 62 bits so adding the rounding amount doesnt
# cause the int to go over the max value of a 64bit int.
CHANNEL_ERASE_DIVISOR = 2**62

# this is the default amount of bits per pixel for palettized bitmaps
DEFAULT_INDEXING_SIZE = 8

#HALF FLOAT "f" WILL REQUIRE Numpy
PIXEL_ENCODING_SIZES_F = {"f":2, "F":4}
INVERSE_PIXEL_ENCODING_SIZES_F = {2:"f", 4:"F"}



"""##################"""
### CHANNEL MAPPINGS ###
"""##################"""

# the way channel mappings work is that each index is one channel.
# in their standard form the value at each index should be the
# number of that index. Ex:(0, 1, 2, 3)

# to remove channels you should create a mapping with exactly
# how many channels you want to have and have the value of each
# index be the channel that you want to place there.
# Ex: (A, R, G, B) to (G, B) would use the mapping (2, 3)

# to switch channels around you would create a mapping with one
# index for each channel in the target format. the value at each
# index will be the index of the channel you want from the source format.
# Ex: (A, R, G, B) to (B, G, R, A) would use the mapping (3, 2, 1, 0)

# if you want a blank channel to be made then set the value at that index to -1
# Ex: converting A8 to A8R8G8B8 = (0, -1, -1, -1)


"""these channel mappings are used to swap ALPHA AND
INTENSITY, but ONLY if the source bitmap is A8Y8"""
#          ( A, L)
AL_TO_LA = ( 1, 0)
A_TO_AL  = ( 0,-1)
L_TO_AL  = (-1, 0)

"""these channel mappings are used to convert different
formats to Y8 and A8. these are also used for converting to AY8.
just use the one that preserves the channel you want to keep"""
#               (A)
ANYTHING_TO_A = (0,)
#         (L)
AL_TO_L = (1,)

"""these channel mappings are to convert
A8, Y8, AY, and YA to A8R8G8B8 and X8R8G8B8"""
#            ( A,  R,  G,  B)
A_TO_ARGB  = ( 0, -1, -1, -1)
L_TO_ARGB  = (-1,  0,  0,  0)

AL_TO_ARGB = ( 0,  1,  1,  1)
LA_TO_ARGB = ( 1,  0,  0,  0)


"""########################"""
### CHANNEL MERGE MAPPINGS ###
"""########################"""

# why merge mappings are used is that if the target format has
# less channels than the source then we either need to remove
# channels or merge them together. we can remove them with the
# above channel mappings, but for things like RGB to monochrome
# we need to merge the pixels together to get the average intensity.

#merge mapping are the length of the source's channel count. each index is
#which channel in the target format to merge the channel from the source into.
#Ex: merging ARGB's 4 channel into A8Y8 would be (0, 1, 1, 1)

#              ( A,  R,  G,  B )
M_ARGB_TO_AL = ( 0,  1,  1,  1 )
M_ARGB_TO_LA = ( 1,  0,  0,  0 )
M_ARGB_TO_L  = ( -1, 0,  0,  0 )
M_ARGB_TO_A  = ( 0, -1, -1, -1 )
