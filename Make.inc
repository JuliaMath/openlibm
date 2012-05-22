CC=gcc
CFLAGS=-Wall -O2 -I. -I../include -I../ld128 -I../src -D__BSD_VISIBLE -Wno-implicit-function-declaration

default: all

%.o: %.c
	$(QUIET_CC)$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o *~

# OS-specific stuff
OS = $(shell uname)
ARCH = $(shell uname -m)


ifeq ($(OS), MINGW32_NT-6.1)
OS=WINNT
endif

ifeq ($(OS), Linux)
SHLIB_EXT = so
CFLAGS+=-fPIC
endif

ifeq ($(OS), FreeBSD)
SHLIB_EXT = so
CFLAGS+=-fPIC
endif

ifeq ($(OS), Darwin)
SHLIB_EXT = dylib
CFLAGS+=-std=c99 -fPIC
endif

ifeq ($(OS), WINNT)
SHLIB_EXT = dll
endif

# Colors for make

CCCOLOR="\033[34m"
LINKCOLOR="\033[34;1m"
SRCCOLOR="\033[33m"
BINCOLOR="\033[37;1m"
MAKECOLOR="\033[32;1m"
ENDCOLOR="\033[0m"

QUIET_CC = @printf '    %b %b\n' $(CCCOLOR)CC$(ENDCOLOR) $(SRCCOLOR)$@$(ENDCOLOR);
QUIET_LINK = @printf '    %b %b\n' $(LINKCOLOR)LINK$(ENDCOLOR) $(BINCOLOR)$@$(ENDCOLOR);
