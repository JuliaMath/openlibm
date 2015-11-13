#!/bin/sh
set -eux

case "$TARGET" in
host)
	uname -a
	export LOADER=
	make ;;
arm32)
	sudo bash -c 'echo >> /etc/apt/sources.list "deb http://archive.ubuntu.com/ubuntu/ trusty main restricted universe"'
	sudo apt-get update
	sudo apt-get -y install gcc-4.7-arm-linux-gnueabihf qemu binfmt-support
	make CC="arm-linux-gnueabihf-gcc-4.7"
	export LD_LIBRARY_PATH=/usr/arm-linux-gnueabihf/lib
	#export LOADER=/usr/arm-linux-gnueabihf/lib/ld-linux-armhf.so.3
	export LOADER="echo TESTS DISABLED ON ARM"
	;;
*)
	echo 'Unknown TARGET!'
	exit 1
	;;
esac

$LOADER make check

make clean && git status --ignored --porcelain && test -z "$(git status --ignored --porcelain)"
