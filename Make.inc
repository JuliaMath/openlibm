# -*- mode: makefile-gmake -*-

OS = $(shell uname)
ARCH = $(shell uname -m)

FC = gfortran
FFLAGS += -O3

USEGCC = 1
USECLANG = 0

CFLAGS += -std=c99 -Wall -O3 -I$(OPENLIBM_HOME) -I$(OPENLIBM_HOME)/include -I$(OPENLIBM_HOME)/ld80 -I$(OPENLIBM_HOME)/$(ARCH) -I$(OPENLIBM_HOME)/src -DASSEMBLER -D__BSD_VISIBLE -Wno-implicit-function-declaration

ifeq ($(USECLANG),1)
USEGCC = 0
CC = clang
CFLAGS += -fno-builtin
endif

ifeq ($(USEGCC),1)
CC = gcc
CFLAGS += -fno-gnu89-inline
endif
AR = ar

default: all

%.c.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

%.f.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.S.o: %.S
	$(CC) $(filter -m% -B% -I% -D%,$(CFLAGS)) -c $< -o $@

clean:
	rm -fr *.o *.c.o *.S.o *~ test-double test-float test-double-system test-float-system *.dSYM

# OS-specific stuff
ifeq ($(ARCH),i386)
override ARCH = i387
endif

ifeq ($(ARCH),x86_64)
override ARCH = amd64
endif 

ifeq ($(ARCH),i686)
override ARCH = i387
endif

ifeq ($(OS), MINGW32_NT-6.1)
OS=WINNT
CFLAGS+=-nodefaultlibs
FFLAGS+=-nodefaultlibs
endif

ifeq ($(OS), Linux)
SHLIB_EXT = so
CFLAGS+=-fPIC
FFLAGS+=-fPIC
endif

ifeq ($(OS), FreeBSD)
SHLIB_EXT = so
CFLAGS+=-fPIC
FFLAGS+=-fPIC
endif

ifeq ($(OS), Darwin)
SHLIB_EXT = dylib
CFLAGS+=-fPIC
FFLAGS+=-fPIC
endif

ifeq ($(OS), WINNT)
SHLIB_EXT = dll
endif
