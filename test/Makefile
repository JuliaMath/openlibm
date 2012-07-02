OPENLIBM_HOME=$(abspath ..)
include ../Make.inc

all: test-double test-float

test-double: test-double.c
	gcc $< -D__BSD_VISIBLE -I ../include -I../src ../libopenlibm.a -o $@

test-float: test-float.c
	gcc $< -D__BSD_VISIBLE -I ../include -I../src ../libopenlibm.a -o $@

clean:
	rm -f *~ *# test-double test-float
