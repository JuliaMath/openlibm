OPENLIBM_HOME=$(abspath ..)
include ../Make.inc

# Set rpath of tests to builddir for loading shared library
OPENLIBM_LIB = -L.. -lopenlibm
ifneq ($(OS),WINNT)
ifneq ($(OS),Darwin)
OPENLIBM_LIB += -Wl,-rpath=$(OPENLIBM_HOME)
endif
else # WINNT
CFLAGS_add += -DIMPORT_EXPORTS
endif

all: test-double test-float # test-double-system test-float-system

# Per-issue regression tests: every test/regression/*.c is a self-contained
# program named after the issue it covers. New files are picked up automatically
# (no Makefile edit needed); the top-level Makefile runs them via `make test`.
REGRESSION_SRCS := $(sort $(wildcard regression/*.c))
REGRESSION_BINS := $(REGRESSION_SRCS:.c=)

regression-build: $(REGRESSION_BINS)

regression/%: regression/%.c regression/regress-util.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_add) $(LDFLAGS) $(LDFLAGS_arch) $< -D__BSD_VISIBLE -I ../include -I../src $(OPENLIBM_LIB) -o $@

bench: bench-syslibm bench-openlibm

test-double: test-double.c libm-test.c libm-test-ulps.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_add) $(LDFLAGS) $(LDFLAGS_arch) $@.c -D__BSD_VISIBLE -I ../include -I../src $(OPENLIBM_LIB) -o $@

test-float: test-float.c libm-test.c libm-test-ulps.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_add) $(LDFLAGS) $(LDFLAGS_arch) $@.c -D__BSD_VISIBLE -I ../include -I../src $(OPENLIBM_LIB) -o $@

test-double-system: test-double.c libm-test.c libm-test-ulps.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_add) $(LDFLAGS) $(LDFLAGS_arch) $< -DSYS_MATH_H -lm -o $@

test-float-system: test-float.c libm-test.c libm-test-ulps.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_add) $(LDFLAGS) $(LDFLAGS_arch) $< -DSYS_MATH_H -lm -o $@

bench-openlibm: libm-bench.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_add) $(LDFLAGS) $(LDFLAGS_arch) $< $(OPENLIBM_LIB) -o $@

bench-syslibm: libm-bench.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_add) $(LDFLAGS) $(LDFLAGS_arch) $< -lm -o $@

clean:
	rm -fr test-double test-float test-double-system test-float-system bench-openlibm bench-syslibm *.dSYM
	rm -fr $(REGRESSION_BINS) $(addsuffix .exe,$(REGRESSION_BINS)) regression/*.dSYM

.PHONY: all bench regression-build clean
