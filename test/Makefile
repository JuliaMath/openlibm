OPENLIBM_HOME=$(abspath ..)
include ../Make.inc

all: test-double test-float test-double-system test-float-system

test-double: test-double.c libm-test.c
	$(CC) -g $< -D__BSD_VISIBLE -I ../include -I../src ../libopenlibm.a -o $@

test-float: test-float.c libm-test.c
	$(CC) -g $< -D__BSD_VISIBLE -I ../include -I../src ../libopenlibm.a -o $@

test-double-system: test-double.c libm-test.c
	$(CC) -g $< -DSYS_MATH_H -lm -o $@

test-float-system: test-float.c libm-test.c
	$(CC) -g $< -DSYS_MATH_H -lm -o $@

clean:
	rm -f *~ *# test-double test-float test-double-system test-float-system

