from distutils.core import setup, Extension

_faststat = Extension('_faststat', sources=['faststat/_faststat.c'])

setup(
    name='faststat',
    version='0.1',
    description='fast online statistics collection',
    ext_modules=[_faststat])
