try:
    from setuptools import setup, Extension, Command
except ImportError:
    from distutils.core import setup, Extension, Command

import arbytmap

long_desc = open("README.md").read()

setup(
    name="arbytmap",
    description='A texture manipulation module for python 3.',
    long_description=long_desc,
    long_description_content_type='text/markdown',
    version='%s.%s.%s' % arbytmap.__version__,
    url='https://github.com/Sigmmma/arbytmap',
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
        'arbytmap': ["src/*", '*.[tT][xX][tT]', '*.[mM][dD]'],
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
    )
