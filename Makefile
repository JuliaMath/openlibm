include ./Make.inc

all:
	$(MAKE) -C src all
	$(MAKE) -C ld128 all
	$(MAKE) -C bsdsrc all
	$(QUIET_LINK)ar -rcs libopenlibm.a src/*.o ld128/*.o bsdsrc/*.o
	$(QUIET_LINK)$(CC) -shared -fPIC src/*.o ld128/*.o bsdsrc/*.o -o libopenlibm.$(SHLIB_EXT)

cleanall:
	$(MAKE) -C src clean
	$(MAKE) -C ld128 clean
	$(MAKE) -C bsdsrc clean
	rm -f *.a *.$(SHLIB_EXT)
