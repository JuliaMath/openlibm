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
	$(FC) -shared $(OBJS) $(LDFLAGS) -o libopenlibm.$(SHLIB_EXT).$(VERSION)
	ln -s libopenlibm.$(SHLIB_EXT).$(VERSION) libopenlibm.$(SHLIB_EXT).$(word 1,$(VERSION_SPLIT)).$(word 2,$(VERSION_SPLIT))
	ln -s libopenlibm.$(SHLIB_EXT).$(VERSION) libopenlibm.$(SHLIB_EXT).$(word 1,$(VERSION_SPLIT))
	ln -s libopenlibm.$(SHLIB_EXT).$(VERSION) libopenlibm.$(SHLIB_EXT)


clean:
	@for dir in $(SUBDIRS) .; do \
		rm -fr $$dir/*.o $$dir/*.a $$dir/*.$(SHLIB_EXT)*; \
	done

distclean:
	rm -f $(OBJS) *.a *.$(SHLIB_EXT)
	$(MAKE) -C test clean

install: all
	mkdir -p $(DESTDIR)$(libdir)
	mkdir -p $(DESTDIR)$(PREFIX)/include
	cp -a libopenlibm.$(SHLIB_EXT)* libopenlibm.a $(DESTDIR)$(libdir)/
	cp -a src/openlibm.h $(DESTDIR)$(PREFIX)/include/
