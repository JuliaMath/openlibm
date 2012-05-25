include ./Make.inc

all:
	$(MAKE) -C src all
	$(MAKE) -C ld128 all
	$(MAKE) -C bsdsrc all
	ar -rcs libopenlibm.a src/*.c.o ld128/*.c.o bsdsrc/*.c.o
	$(CC) -shared -fPIC src/*.c.o ld128/*.c.o bsdsrc/*.c.o -o libopenlibm.$(SHLIB_EXT)

cleanall:
	$(MAKE) -C src clean
	$(MAKE) -C ld128 clean
	$(MAKE) -C bsdsrc clean
	rm -f *.a *.$(SHLIB_EXT)
