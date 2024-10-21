all:

install:
	. ../.olcao/olcaorc && \
	cd ../build/release && \
	export FC=h5fc && \
	. /etc/environ.sh && use -e -r cmake-3.6.3 && \
	cmake ../.. && \
	make install && \
	unpackOLCAODB

clean:
	rm -rf ../build/release/*

distclean: clean
	rm -rf ../bin
	rm -rf ../share

.PHONY: all install clean distclean
