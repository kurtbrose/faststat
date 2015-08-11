import platform
from distutils.core import setup, Extension

extra_compile_args = []
if platform.system() == 'Windows':
    extra_compile_args = ['/MT',]
libraries = []
if platform.system() not in ('Windows', 'Darwin'):
    libraries = ['rt',]

_faststat = Extension('_faststat', sources=['faststat/_faststat.c'], 
                      libraries=libraries, extra_compile_args=extra_compile_args)

setup(
    name='faststat',
    version='0.3.1',
    author="Kurt Rose",
    author_email="kurt@kurtrose.com",
    description='fast online statistics collection',
    license="MIT",
    url="http://github.com/doublereedkurt/faststat",
    long_description='...',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
    ],
    packages=['faststat'],
    ext_modules=[_faststat])


