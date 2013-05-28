# -*- mode: makefile-gmake -*-

OS := $(shell uname)

FC = gfortran
FFLAGS += -O3

USEGCC = 1
USECLANG = 0


ifeq ($(USECLANG),1)
USEGCC = 0
CC = clang
CFLAGS_add += -fno-builtin
endif

ifeq ($(USEGCC),1)
CC = gcc
CFLAGS_add += -fno-gnu89-inline
endif
AR = ar

ARCH := $(shell $(CC) -dumpmachine | sed "s/\([^-]*\).*$$/\1/")

CFLAGS_add += -std=c99 -Wall -O3 -I$(OPENLIBM_HOME) -I$(OPENLIBM_HOME)/include -I$(OPENLIBM_HOME)/ld80 -I$(OPENLIBM_HOME)/$(ARCH) -I$(OPENLIBM_HOME)/src -DASSEMBLER -D__BSD_VISIBLE -Wno-implicit-function-declaration

default: all

%.c.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS_add) $(CFLAGS) -c $< -o $@

%.f.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.S.o: %.S
	$(CC) $(SFLAGS) $(filter -m% -B% -I% -D%,$(CFLAGS_add)) -c $< -o $@

clean:
	rm -fr *.o *.c.o *.S.o *~ test-double test-float test-double-system test-float-system *.dSYM

# OS-specific stuff
ifeq ($(ARCH),i386)
override ARCH := i387
endif

ifeq ($(ARCH),x86_64)
override ARCH := amd64
endif 

ifeq ($(ARCH),i686)
override ARCH := i387
endif

ifneq (,$(findstring MINGW,$(OS)))
override OS=WINNT
endif
#keep these if statements these separate
ifeq ($(OS), WINNT)
CFLAGS_add+=-nodefaultlibs
FFLAGS+=-nodefaultlibs
endif

ifeq ($(OS), Linux)
SHLIB_EXT = so
CFLAGS_add+=-fPIC
FFLAGS+=-fPIC
endif

ifeq ($(OS), FreeBSD)
SHLIB_EXT = so
CFLAGS_add+=-fPIC
FFLAGS+=-fPIC
endif

ifeq ($(OS), Darwin)
SHLIB_EXT = dylib
CFLAGS_add+=-fPIC
FFLAGS+=-fPIC
endif

ifeq ($(OS), WINNT)
SHLIB_EXT = dll
endif
