for i in *.c; do clang -D__BSD_VISIBLE -Wno-implicit-function-declaration -c -O0 -I. -I../include -I../ld128 $i; done 2>err
