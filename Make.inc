OS = $(shell uname)
ARCH = $(shell uname -m)

USEGCC = 1
USECLANG = 0

ifeq ($(USECLANG),1)
USEGCC = 0
CC = clang
CFLAGS= -std=c99 -Wall -O4 -I$(OPENLIBM_HOME) -I$(OPENLIBM_HOME)/include -I$(OPENLIBM_HOME)/ld80 -I$(OPENLIBM_HOME)/$(ARCH) -I$(OPENLIBM_HOME)/src -DASSEMBLER -D__BSD_VISIBLE -Wno-implicit-function-declaration
endif

ifeq ($(USEGCC),1)
CC = gcc
CFLAGS= -std=c99 -Wall -O3 -I$(OPENLIBM_HOME) -I$(OPENLIBM_HOME)/include -I$(OPENLIBM_HOME)/ld80 -I$(OPENLIBM_HOME)/$(ARCH) -I$(OPENLIBM_HOME)/src -DASSEMBLER -D__BSD_VISIBLE -Wno-implicit-function-declaration
endif

default: all

%.c.o: %.c
	$(QUIET_CC)$(CC) $(CFLAGS) -c $< -o $@

%.S.o: %.S
	$(QUIET_ASM)$(CC) $(filter -f% -m% -B% -I% -D%,$(CFLAGS)) -c $< -o $@

clean:
	rm -f *.o *.c.o *.S.o *~

# OS-specific stuff
ifeq ($(ARCH),i386)
override ARCH = i387
endif

ifeq ($(ARCH),x86_64)
override ARCH = amd64
endif 

ifeq ($(ARCH),i686)
override ARCH = i387
endif

ifeq ($(OS), MINGW32_NT-6.1)
OS=WINNT
CFLAGS+=-nodefaultlibs
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
CFLAGS+=-fPIC
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
QUIET_ASM = @printf '    %b %b\n' $(CCCOLOR)ASM$(ENDCOLOR) $(SRCCOLOR)$@$(ENDCOLOR);
QUIET_LINK = @printf '    %b %b\n' $(LINKCOLOR)LINK$(ENDCOLOR) $(BINCOLOR)$@$(ENDCOLOR);
