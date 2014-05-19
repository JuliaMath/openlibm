OPENLIBM_HOME=$(abspath .)
include ./Make.inc

SUBDIRS = src ld80 $(ARCH) bsdsrc

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
	$(CC) -shared $(OBJS) $(LDFLAGS) -Wl,$(SONAME_FLAG),libopenlibm.$(SHLIB_EXT) -o libopenlibm.$(SHLIB_EXT)
else
	$(CC) -shared $(OBJS) $(LDFLAGS) -Wl,$(SONAME_FLAG),libopenlibm.$(SHLIB_EXT).$(VERSION) -o libopenlibm.$(SHLIB_EXT).$(VERSION)
	@-ln -sf libopenlibm.$(SHLIB_EXT).$(VERSION) libopenlibm.$(SHLIB_EXT).$(word 1,$(VERSION_SPLIT)).$(word 2,$(VERSION_SPLIT))
	@-ln -sf libopenlibm.$(SHLIB_EXT).$(VERSION) libopenlibm.$(SHLIB_EXT).$(word 1,$(VERSION_SPLIT))
	@-ln -sf libopenlibm.$(SHLIB_EXT).$(VERSION) libopenlibm.$(SHLIB_EXT)
endif

clean:
	@for dir in $(SUBDIRS) .; do \
		rm -fr $$dir/*.o $$dir/*.a $$dir/*.$(SHLIB_EXT)*; \
	done

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
	cp -a libopenlibm.$(SHLIB_EXT)* $(DESTDIR)$(shlibdir)/
	cp -a libopenlibm.a $(DESTDIR)$(libdir)/
	cp -a src/openlibm.h $(DESTDIR)$(includedir)/
	cp -a openlibm.pc $(DESTDIR)$(libdir)/pkgconfig/
ifneq ($(wildcard $(ARCH)/bsd_asm.h),)
	cp -a $(ARCH)/bsd_asm.h $(DESTDIR)$(includedir)/openlibm/
endif
ifneq ($(wildcard $(ARCH)/bsd_cdefs.h),)
	cp -a $(ARCH)/bsd_cdefs.h $(DESTDIR)$(includedir)/openlibm/
endif
