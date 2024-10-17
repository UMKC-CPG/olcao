all:

install:
	cd ../../build/release
	export FC=hdf5
	cmake3 ../..
	make install
	unpackOLCAODB

clean:

distclean: clean

.PHONY: all install clean distclean
