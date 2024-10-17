all:

install:
	. ../.olcao/olcaorc && \
	cd ../build/release && \
	export FC=h5fc && \
	export PATH=/apps/share64/debian7/cmake/cmake-3.6.3/bin:${PATH} && \
	cmake ../.. && \
	make install && \
	unpackOLCAODB

clean:

distclean: clean

.PHONY: all install clean distclean
