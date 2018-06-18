# OpenLibm

[![Travis](https://travis-ci.org/JuliaMath/openlibm.svg?branch=master)](https://travis-ci.org/JuliaMath/openlibm)
[![AppVeyor](https://ci.appveyor.com/api/projects/status/sia04r4089rr19uc/branch/master?svg=true)](https://ci.appveyor.com/project/ararslan/openlibm-19152/branch/master)

[OpenLibm](http://www.openlibm.org) is an effort to have a high quality, portable, standalone
C mathematical library ([`libm`](http://en.wikipedia.org/wiki/libm)).
It can be used standalone in applications and programming language
implementations.

The project was born out of a need to have a good `libm` for the
[Julia programming langage](http://www.julialang.org) that worked
consistently across compilers and operating systems, and in 32-bit and
64-bit environments.

## Platform support

OpenLibm builds on Linux, Mac OS X, Windows, FreeBSD, OpenBSD, and DragonFly BSD.
It builds with both GCC and clang. Although largely tested and widely
used on x86 architectures, OpenLibm also supports ARM and
PowerPC.

## Build instructions

1. Use GNU Make to build OpenLibm. This is `make` on most systems, but `gmake` on BSDs.
2. Use `make USEGCC=1` to build with GCC. This is the default on
   Linux and Windows.
3. Use `make USECLANG=1` to build with clang. This is the default on OS X, FreeBSD,
   and OpenBSD.
4. Architectures are auto-detected. Use `make ARCH=i386` to force a
   build for i386. Other supported architectures are i486, i586, and
   i686. GCC 4.8 is the minimum requirement for correct codegen on
   older 32-bit architectures.

## Acknowledgements

PowerPC support for openlibm was graciously sponsored by IBM.
