#!/bin/bash

set -e
set -x

cd ..

gcc -fno-gnu89-inline -fno-builtin -O3 -fPIC -m64 -std=c99 -Wall -Iinclude -DASSEMBLER -D__BSD_VISIBLE -Wno-implicit-function-declaration -c src/s_exp2.c -o src/s_exp2.c.o
