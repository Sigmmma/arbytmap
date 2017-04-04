# Arbytmap

## What is this repository for?

* Arbytmap is a power-of-2 bitmap conversion and manipulation module. Arbytmap is intended for use in converting bitmaps from one format to another, half-rezzing and generating mipmaps, swizzling bitmaps, and writing converted bitmaps to files. Arbytmap is currently undergoing a large scale cleanup, as much of it was written years ago when I was a much worse programmer.

* Many of the conversion functions have accelerator modules written in C, so this module(if properly compiled and installed) can reach speeds seen in lower level languages.

* Raw formats currently supported by this module are as follows:
```A8, Y8, AY8, A8Y8, R3G3B2, 5G6B5, R8G8B8, Y8U8V8, A1R5G5B5, A4R4G4B4, X8R8G8B8, A8R8G8B8, R16G16B16, A16R16G16B16```

* DXT formats currently supported by this module are as follows:
```DXT1/2/3/4/5, DXN, DXT5A, DXT5Y, DXT5AY, CTX1, U8V8```

## Todo

* Completely redo the dds and tga reading/writing system(use [supyr_struct](https://bitbucket.org/moses_of_egypt/supyr_struct) for handling creating and reading the files).

* Make/finish C functions for working with the remaining formats.

* Clean up and standardize the interface for loading bitmaps and conversion settings into the Arbytmap class.

* Cleanup/redo pretty much everything that looks/functions bad.

* Anything else I can think of(I have the flu right now, so I'm having a hard time writing this).

## Who do I talk to?

* Devin Bobadilla (Author of arbytmap) mosesbobadilla@gmail.com