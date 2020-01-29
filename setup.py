#!/usr/bin/env python
import sys
from traceback import format_exc

try:
    from setuptools import setup, Extension, Command
except ImportError:
    from distutils.core import setup, Extension, Command
from distutils.command.build_ext import build_ext
from distutils.errors import CCompilerError, DistutilsExecError, \
     DistutilsPlatformError

import arbytmap


is_pypy = hasattr(sys, 'pypy_translation_info')
ext_errors = None
if sys.platform == 'win32':
   ext_errors = (CCompilerError, DistutilsExecError, DistutilsPlatformError,
                 IOError, ValueError)

class BuildFailed(Exception):
    pass

class ve_build_ext(build_ext):
    # This class allows C extension building to fail.

    def run(self):
        try:
            build_ext.run(self)
        except DistutilsPlatformError:
            raise BuildFailed()

    def build_extension(self, ext):
        if ext_errors:
            try:
                build_ext.build_extension(self, ext)
            except ext_errors:
                raise BuildFailed()
        else:
            build_ext.build_extension(self, ext)


long_desc = open(curr_dir, "README.md").read()

setup_kwargs = dict(
    name="arbytmap",
    description='A texture manipulation module for python 3.',
    long_description=long_desc,
    long_description_content_type='text/markdown',
    version='%s.%s.%s' % arbytmap.__version__,
    url='https://github.com/MosesofEgypt/arbytmap',
    author='Devin Bobadilla',
    author_email='MosesBobadilla@gmail.com',
    license='MIT',
    packages=[
        'arbytmap',
        'arbytmap.ext',
        ],
    ext_modules = [
        Extension("arbytmap.ext.arbytmap_ext",     ["arbytmap/src/arbytmap_ext.c"]),
        Extension("arbytmap.ext.bitmap_io_ext",    ["arbytmap/src/bitmap_io_ext.c"]),
        Extension("arbytmap.ext.dds_defs_ext",     ["arbytmap/src/dds_defs_ext.c"]),
        Extension("arbytmap.ext.raw_packer_ext",   ["arbytmap/src/raw_packer_ext.c"]),
        Extension("arbytmap.ext.raw_unpacker_ext", ["arbytmap/src/raw_unpacker_ext.c"]),
        Extension("arbytmap.ext.swizzler_ext",     ["arbytmap/src/swizzler_ext.c"])
        ],
    package_data={
        'arbytmap': ["src/*", '*.txt', '*.md', '*.rst'],
        },
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
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Multimedia :: Graphics :: Graphics Conversion",
        "Programming Language :: C",
        ],
    zip_safe=False,
    cmdclass=dict(build_ext=ve_build_ext)
    )

success = False
kwargs = dict(setup_kwargs)
if not is_pypy:
    try:
        setup(**kwargs)
        success = True
    except BuildFailed:
        print(format_exc())
        print('*' * 80)
        print("WARNING: The C accelerator modules could not be compiled.\n"
              "Attempting to install without accelerators now.\n"
              "Any errors that occurred are printed above.")
        print('*' * 80)

if not success:
    kwargs.pop('ext_modules')
    setup(**kwargs)
    print("Installation successful.")
