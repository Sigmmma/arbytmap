# Arbytmap

## What is this repository for?

* Current Version: 0.6.0

* Arbytmap is a power-of-2 bitmap conversion and manipulation module. Arbytmap is intended for use in converting bitmaps from one format to another, half-rezzing and generating mipmaps, swizzling bitmaps, and writing converted bitmaps to files. Arbytmap is currently undergoing a large scale cleanup, as much of it was written years ago when I was a much worse programmer.

* Many of the conversion functions have accelerator modules written in C, so this module(if properly compiled and installed) can reach speeds seen in lower level langauges.


## Todo

* Completely redo the dds and tga reading/writing system(use supyr_struct for handling creating and reading the files).

* Make/finish C functions for working with the remaining formats.

* Clean up and standardize the interface for loading bitmaps and conversion settings into the Arbytmap class.

* Clean/redo up pretty much everything that looks/functions bad.

* Anything else I can think of(I have the flu right now, so I'm having a hard time writing this).

## Who do I talk to?

* [Devin Bobadilla](mosesbobadilla@gmail.com)(Author of arbytmap)