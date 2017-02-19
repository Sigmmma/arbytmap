#!/usr/bin/env python
from os.path import dirname, join
try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

curr_dir = dirname(__file__)

#               YYYY.MM.DD
release_date = "2017.02.15"
version = (0, 6, 0)

try:
    try:
        long_desc = open(join(curr_dir, "readme.rst")).read()
    except Exception:
        long_desc = open(join(curr_dir, "readme.md")).read()
except Exception:
    long_desc = 'Could not read long description from readme.'

setup(
    name="arbytmap",
    description='A power-of-2 texture manipulation module for python 3.',
    long_description=long_desc,
    version="0.6.0",
    url='http://bitbucket.org/Moses_of_Egypt/arbytmap',
    author='Devin Bobadilla',
    author_email='MosesBobadilla@gmail.com',
    license='MIT',
    packages=[
        'arbytmap',
        'arbytmap.ext',
        ],
    ext_modules = [
        Extension("arbytmap.ext.arbytmap_ext", ["arbytmap\\src\\arbytmap_ext.c"]),
        Extension("arbytmap.ext.bitmap_io_ext", ["arbytmap\\src\\bitmap_io_ext.c"]),
        Extension("arbytmap.ext.dds_defs_ext", ["arbytmap\\src\\dds_defs_ext.c"]),
        Extension("arbytmap.ext.raw_packer_ext", ["arbytmap\\src\\raw_packer_ext.c"]),
        Extension("arbytmap.ext.raw_unpacker_ext", ["arbytmap\\src\\raw_unpacker_ext.c"]),
        Extension("arbytmap.ext.swizzler_ext", ["arbytmap\\src\\swizzler_ext.c"])
        ],
    platforms=["POSIX", "Windows"],
    keywords="arbytmap, texture, bitmap, converter, image, editing",
    install_requires=[],
    requires=[],
    provides=['arbytmap'],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        ],
    zip_safe=False,
    )
