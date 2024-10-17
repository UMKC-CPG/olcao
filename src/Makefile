all:

install:
	cd ../build/release
	export FC=hdf5
	use -p cmake-3.6.3
	cmake ../..
	make install
	unpackOLCAODB

clean:

distclean: clean

.PHONY: all install clean distclean
