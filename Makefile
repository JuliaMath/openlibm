include ./Make.inc

all:
	$(MAKE) -C src all
	$(MAKE) -C ld128 all
	$(QUIET_LINK)ar -rcs libopenlibm.a src/*.o ld128/*.o
	$(QUIET_LINK)$(CC) -shared -fPIC src/*.o ld128/*.o -o libopenlibm.$(SHLIB_EXT)

cleanall:
	$(MAKE) -C src clean
	$(MAKE) -C ld128 clean
	rm -f *.a *.$(SHLIB_EXT)
