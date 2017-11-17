#!/bin/bash

set -e
set -x

FFLAGS="-Wall -Wextra -Wimplicit-interface -fPIC -g -fcheck=all -fbacktrace"
FFLAGS="-Wall -Wextra -Wimplicit-interface -fPIC -O3 -march=native -ffast-math -funroll-loops"

CFLAGS="-fno-gnu89-inline -fno-builtin -fPIC -m64 -std=c99 -Wall -O3"
CFLAGS="$CFLAGS -march=native -funroll-loops"
#LTO="-flto"
gcc $LTO $CFLAGS -Wno-implicit-function-declaration -c redux.c -o redux.o

#CFLAGS="$CFLAGS -ffast-math"

gcc $LTO $CFLAGS -I../include -DASSEMBLER -D__BSD_VISIBLE -Wno-implicit-function-declaration -c ../src/s_exp2.c -o s_exp2.c.o
gcc $LTO $CFLAGS -I../include -DASSEMBLER -D__BSD_VISIBLE -Wno-implicit-function-declaration -c ../src/e_pow.c -o e_pow.c.o
gfortran $LTO $FFLAGS -c test_exp2.f90 -o test_exp2.o
gfortran $LTO $FFLAGS test_exp2.o s_exp2.c.o e_pow.c.o redux.o -o test_exp2

./test_exp2
