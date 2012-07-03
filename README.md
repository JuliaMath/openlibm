OpenLIBM is an effort to have a high quality independent LIBM
library. It is meant to be used standalone in applications and
programming language implementations, and perhaps even as a reference
for LIBM implementations in OSes.

OpenLIBM builds on Linux and Mac OS X, and with little effort, 
should build on FreeBSD as well. It builds with both, GCC and clang.

The OpenLIBM code derives from the FreeBSD msun implementation, which
in turn derives from FDLIBM 5.3.

Build instructions:

`make` or `make USEGCC=1` to build with GCC.
`make USECLANG=1` to build with clang.
