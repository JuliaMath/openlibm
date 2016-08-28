OPENLIBM_HOME=$(abspath .)
include ./Make.inc

SUBDIRS = src $(ARCH) bsdsrc
ifneq ($(ARCH), arm)
ifneq ($(ARCH), powerpc)
SUBDIRS += ld80
endif
endif

define INC_template
TEST=test
override CUR_SRCS = $(1)_SRCS
include $(1)/Make.files
SRCS += $$(addprefix $(1)/,$$($(1)_SRCS))
endef

DIR=test

$(foreach dir,$(SUBDIRS),$(eval $(call INC_template,$(dir))))

DUPLICATE_NAMES = $(filter $(patsubst %.S,%,$($(ARCH)_SRCS)),$(patsubst %.c,%,$(src_SRCS)))
DUPLICATE_SRCS = $(addsuffix .c,$(DUPLICATE_NAMES))

OBJS =  $(patsubst %.f,%.f.o,\
	$(patsubst %.S,%.S.o,\
	$(patsubst %.c,%.c.o,$(filter-out $(addprefix src/,$(DUPLICATE_SRCS)),$(SRCS)))))

# If we're on windows, don't do versioned shared libraries. Also, generate an import library
# for the DLL. If we're on OSX, put the version number before the .dylib.  Otherwise,
# put it after.
ifeq ($(OS), WINNT)
OLM_MAJOR_MINOR_SHLIB_EXT := $(SHLIB_EXT)
LDFLAGS_add += -Wl,--out-implib,libopenlibm.$(OLM_MAJOR_MINOR_SHLIB_EXT).a
else
ifeq ($(OS), Darwin)
OLM_MAJOR_MINOR_SHLIB_EXT := $(SOMAJOR).$(SOMINOR).$(SHLIB_EXT)
OLM_MAJOR_SHLIB_EXT := $(SOMAJOR).$(SHLIB_EXT)
else
OLM_MAJOR_MINOR_SHLIB_EXT := $(SHLIB_EXT).$(SOMAJOR).$(SOMINOR)
OLM_MAJOR_SHLIB_EXT := $(SHLIB_EXT).$(SOMAJOR)
endif
endif

.PHONY: all check test clean distclean \
	install install-static install-shared install-pkgconfig install-headers

all: libopenlibm.a libopenlibm.$(OLM_MAJOR_MINOR_SHLIB_EXT)

check test: test/test-double test/test-float
	test/test-double
	test/test-float

libopenlibm.a: $(OBJS)
	$(AR) -rcs libopenlibm.a $(OBJS)

libopenlibm.$(OLM_MAJOR_MINOR_SHLIB_EXT): $(OBJS)
	$(CC) -shared $(OBJS) $(LDFLAGS) $(LDFLAGS_add) -Wl,$(SONAME_FLAG),libopenlibm.$(OLM_MAJOR_SHLIB_EXT) -o $@
ifneq ($(OS),WINNT)
	ln -sf $@ libopenlibm.$(OLM_MAJOR_SHLIB_EXT)
	ln -sf $@ libopenlibm.$(SHLIB_EXT)
endif

test/test-double: libopenlibm.$(OLM_MAJOR_MINOR_SHLIB_EXT)
	$(MAKE) -C test test-double

test/test-float: libopenlibm.$(OLM_MAJOR_MINOR_SHLIB_EXT)
	$(MAKE) -C test test-float

clean:
	rm -f amd64/*.o arm/*.o bsdsrc/*.o i387/*.o ld128/*.o ld80/*.o src/*.o
	rm -f libopenlibm.a libopenlibm.*$(SHLIB_EXT)*
	$(MAKE) -C test clean

openlibm.pc: openlibm.pc.in Make.inc Makefile
	echo "prefix=${prefix}" > openlibm.pc
	echo "version=${VERSION}" >> openlibm.pc
	cat openlibm.pc.in >> openlibm.pc

install-static: libopenlibm.a
	mkdir -p $(DESTDIR)$(libdir)
	cp -f -a libopenlibm.a $(DESTDIR)$(libdir)/

install-shared: libopenlibm.$(OLM_MAJOR_MINOR_SHLIB_EXT)
	mkdir -p $(DESTDIR)$(shlibdir)
	cp -f -a libopenlibm.*$(SHLIB_EXT)* $(DESTDIR)$(shlibdir)/

install-pkgconfig: openlibm.pc
	mkdir -p $(DESTDIR)$(pkgconfigdir)
	cp -f -a openlibm.pc $(DESTDIR)$(pkgconfigdir)/

install-headers:
	mkdir -p $(DESTDIR)$(includedir)/openlibm
	cp -f -a include/*.h $(DESTDIR)$(includedir)/openlibm
	cp -f -a src/*.h $(DESTDIR)$(includedir)/openlibm

install: install-static install-shared install-pkgconfig install-headers
