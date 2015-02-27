OPENLIBM_HOME=$(abspath .)
include ./Make.inc

SUBDIRS = src $(ARCH) bsdsrc
ifneq ($(ARCH), arm)
SUBDIRS += ld80
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

all: libopenlibm.a libopenlibm.$(SHLIB_EXT) 
	$(MAKE) -C test
libopenlibm.a: $(OBJS)  
	$(AR) -rcs libopenlibm.a $(OBJS)
libopenlibm.$(SHLIB_EXT): $(OBJS)
ifeq ($(OS),WINNT)
	$(CC) -shared $(OBJS) $(LDFLAGS) $(LDFLAGS_add) -Wl,$(SONAME_FLAG),libopenlibm.$(SHLIB_EXT) -o libopenlibm.$(SHLIB_EXT)
else
	$(CC) -shared $(OBJS) $(LDFLAGS) $(LDFLAGS_add) -Wl,$(SONAME_FLAG),libopenlibm.$(SHLIB_EXT).$(SOMAJOR) -o libopenlibm.$(SHLIB_EXT).$(SOMAJOR).$(SOMINOR)
	@-ln -sf libopenlibm.$(SHLIB_EXT).$(SOMAJOR).$(SOMINOR) libopenlibm.$(SHLIB_EXT).$(SOMAJOR)
	@-ln -sf libopenlibm.$(SHLIB_EXT).$(SOMAJOR).$(SOMINOR) libopenlibm.$(SHLIB_EXT)
endif

clean:
	@for dir in $(SUBDIRS) .; do \
		rm -fr $$dir/*.o $$dir/*.a $$dir/*.$(SHLIB_EXT)*; \
	done
	@rm -f test/test-double test/test-float

distclean:
	-rm -f $(OBJS) *.a *.$(SHLIB_EXT) libopenlibm.*
	-$(MAKE) -C test clean

openlibm.pc: openlibm.pc.in Make.inc Makefile
	echo "prefix=${prefix}" > openlibm.pc
	echo "version=${VERSION}" >> openlibm.pc
	cat openlibm.pc.in >> openlibm.pc

install: all openlibm.pc
	mkdir -p $(DESTDIR)$(shlibdir)
	mkdir -p $(DESTDIR)$(libdir)/pkgconfig
	mkdir -p $(DESTDIR)$(includedir)/openlibm
	cp -f -a libopenlibm.$(SHLIB_EXT)* $(DESTDIR)$(shlibdir)/
	cp -f -a libopenlibm.a $(DESTDIR)$(libdir)/
	cp -f -a include/openlibm*.h $(DESTDIR)$(includedir)/
	cp -f -a openlibm.pc $(DESTDIR)$(libdir)/pkgconfig/
