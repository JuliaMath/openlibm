# -*- mode: makefile-gmake -*-

OS := $(shell uname)
# Do not forget to bump SOMINOR when changing VERSION,
# and SOMAJOR when breaking ABI in a backward-incompatible way
VERSION = 0.5-dev
SOMAJOR = 2
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

ifeq ($(OS), FreeBSD)
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

ARCH ?= $(shell $(CC) -dumpmachine | sed "s/\([^-]*\).*$$/\1/")
ARCH_ORIGIN := $(origin ARCH)

ifeq ($(ARCH),mingw32)
$(error "the mingw32 compiler you are using fails the openblas testsuite. please see the Julia README.windows.md document for a replacement")
endif

CFLAGS_add += -std=c99 -Wall -O1 -Wno-implicit-function-declaration -DASSEMBLER -D__BSD_VISIBLE -I$(OPENLIBM_HOME) -I$(OPENLIBM_HOME)/include -I$(OPENLIBM_HOME)/ld80 -I$(OPENLIBM_HOME)/$(ARCH) -I$(OPENLIBM_HOME)/src

default: all

%.c.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_add) -c $< -o $@

%.S.o: %.S
	$(CC) $(CPPFLAGS) $(SFLAGS) $(SFLAGS_add) $(filter -m% -B% -I% -D%,$(CFLAGS_add)) -c $< -o $@

# OS-specific stuff
REAL_ARCH := $(ARCH)
ifeq ($(findstring arm,$(ARCH)),arm)
override ARCH := arm
endif
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
  override CFLAGS_add += -nodefaultlibs
  shlibdir = $(bindir)
 else
  ifeq ($(OS), Darwin)
    SHLIB_EXT = dylib
    SONAME_FLAG = -install_name
  else
    SHLIB_EXT = so
    SONAME_FLAG = -soname
  endif
  override CFLAGS_add += -fPIC
  shlibdir = $(libdir)
endif

# The target specific FLAGS_add
ifeq ($(ARCH_ORIGIN),file)
CFLAGS_add_TARGET_$(ARCH) :=
SFLAGS_add_TARGET_$(ARCH) :=
LDFLAGS_add_TARGET_$(ARCH) :=
else
ifeq ($(ARCH),i387)
CFLAGS_add_TARGET_$(ARCH)  := -m32 -march=$(REAL_ARCH)
SFLAGS_add_TARGET_$(ARCH)  := -m32 -march=$(REAL_ARCH)
LDFLAGS_add_TARGET_$(ARCH) := -m32 -march=$(REAL_ARCH)
endif
CFLAGS_add_TARGET_x86_64  := -m64
SFLAGS_add_TARGET_x86_64  := -m64
LDFLAGS_add_TARGET_x86_64 := -m64
# Arm
ifeq ($(ARCH),arm)
ifneq ($(REAL_ARCH),arm)
CFLAGS_add_TARGET_$(ARCH)  := -march=$(REAL_ARCH)
SFLAGS_add_TARGET_$(ARCH)  := -march=$(REAL_ARCH)
LDFLAGS_add_TARGET_$(ARCH) := -march=$(REAL_ARCH)
else
$(error No known generic arm cflags. Please specify arch type)
endif
endif
endif

# Actually finish setting the FLAGS_add
CFLAGS_add += $(CFLAGS_add_TARGET_$(ARCH))
LDFLAGS_add += $(LDFLAGS_add_TARGET_$(ARCH))
SFLAGS_add += $(SFLAGS_add_TARGET_$(ARCH))
