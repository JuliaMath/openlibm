OPENLIBM_HOME=$(abspath ..)
include ../Make.inc

all: test-double test-float # test-double-system test-float-system

bench: bench-syslibm bench-openlibm

test-double: test-double.c libm-test.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_add_TARGET_$(ARCH)) $(LDFLAGS) -g $@.c -D__BSD_VISIBLE -I ../include -I../src ../libopenlibm.a -o $@

test-float: test-float.c libm-test.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_add_TARGET_$(ARCH)) $(LDFLAGS) -g $@.c -D__BSD_VISIBLE -I ../include -I../src ../libopenlibm.a -o $@

test-double-system: test-double.c libm-test.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_add_TARGET_$(ARCH)) $(LDFLAGS) -g $< -DSYS_MATH_H -lm -o $@

test-float-system: test-float.c libm-test.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_add_TARGET_$(ARCH)) $(LDFLAGS) -g $< -DSYS_MATH_H -lm -o $@

bench-openlibm: libm-bench.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_add_TARGET_$(ARCH)) $(LDFLAGS) -O2 $< ../libopenlibm.a -o $@

bench-syslibm: libm-bench.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_add_TARGET_$(ARCH)) $(LDFLAGS) -O2 $< -lm -o $@

clean:
	rm -fr test-double test-float test-double-system test-float-system bench-openlibm bench-syslibm *.dSYM
