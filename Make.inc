# -*- mode: makefile-gmake -*-

OS := $(shell uname)
# Do not forget to bump SOMINOR when changing VERSION,
# and SOMAJOR when breaking ABI in a backward-incompatible way
VERSION = 0.5
SOMAJOR = 1
SOMINOR = 1
DESTDIR =
prefix = /usr/local
bindir = $(prefix)/bin
libdir = $(prefix)/lib
includedir = $(prefix)/include

USEGCC = 1
USECLANG = 0

ifeq ($(OS), Darwin)
USEGCC = 0
USECLANG = 1
endif

AR = ar

ifeq ($(USECLANG),1)
USEGCC = 0
CC = clang
CFLAGS_add += -fno-builtin
endif

ifeq ($(USEGCC),1)
CC = gcc
CFLAGS_add += -fno-gnu89-inline
endif

ARCH := $(shell $(CC) -dumpmachine | sed "s/\([^-]*\).*$$/\1/")
ifeq ($(ARCH),mingw32)
$(error "the mingw32 compiler you are using fails the openblas testsuite. please see the Julia README.windows.md document for a replacement")
endif

CFLAGS_add += -std=c99 -Wall -O3 -I$(OPENLIBM_HOME) -I$(OPENLIBM_HOME)/include -I$(OPENLIBM_HOME)/ld80 -I$(OPENLIBM_HOME)/$(ARCH) -I$(OPENLIBM_HOME)/src -DASSEMBLER -D__BSD_VISIBLE -Wno-implicit-function-declaration

default: all

%.c.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS_add) $(CFLAGS) -c $< -o $@

%.S.o: %.S
	$(CC) $(SFLAGS) $(filter -m% -B% -I% -D%,$(CFLAGS_add)) -c $< -o $@

# OS-specific stuff
ifeq ($(ARCH),i386)
override ARCH := i387
endif
ifeq ($(ARCH),i486)
override ARCH := i387
endif
ifeq ($(ARCH),i586)
override ARCH := i387
endif
ifeq ($(ARCH),i686)
override ARCH := i387
endif

ifeq ($(ARCH),x86_64)
override ARCH := amd64
endif

ifneq (,$(findstring MINGW,$(OS)))
override OS=WINNT
endif

#keep these if statements separate
ifeq ($(OS), WINNT)
SHLIB_EXT = dll
SONAME_FLAG = -soname
CFLAGS_add += -nodefaultlibs
shlibdir = $(bindir)
else
shlibdir = $(libdir)
endif

ifeq ($(OS), Linux)
SHLIB_EXT = so
SONAME_FLAG = -soname
CFLAGS_add+=-fPIC
endif

ifeq ($(OS), FreeBSD)
SHLIB_EXT = so
SONAME_FLAG = -soname
CFLAGS_add+=-fPIC
endif

ifeq ($(OS), Darwin)
SHLIB_EXT = dylib
SONAME_FLAG = -install_name
CFLAGS_add+=-fPIC
endif
