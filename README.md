## OpenLIBM v0.2

OpenLIBM is an effort to have a high quality standalone LIBM
library. It is meant to be used standalone in applications and
programming language implementations.

OpenLIBM builds on Linux, Mac OS X, and Windows, and with little effort, 
should build on FreeBSD as well. It builds with both, GCC and clang.

The OpenLIBM code derives from the FreeBSD msun implementation, which
in turn derives from FDLIBM 5.3. As a result, it has a number of fixes and
updates that have accumulated over the years in msun, and also optimized
assembly versions of many functions.

### Build instructions:

1. `make` or `make USEGCC=1` to build with GCC.
2. `make USECLANG=1` to build with clang.
