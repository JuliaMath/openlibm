for i in *.c; do gcc -D__BSD_VISIBLE -c -O0 -I. -I../include -I../ld128 $i; done 2>err
