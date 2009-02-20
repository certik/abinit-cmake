#
# Makefile for the NetCDF library embedded in ABINIT
#

include ../../config.mk

all_targets all xmlf90: install
	@echo "$(xmlf90_pkg_name) is now ready for use."

uncompress: uncompress-stamp
	@echo "$(xmlf90_pkg_name) has been uncompressed."

configure: configure-stamp
	@echo "$(xmlf90_pkg_name) has been configured."

build: build-stamp
	@echo "$(xmlf90_pkg_name) has been built."

install: install-stamp
	@echo "$(xmlf90_pkg_name) has been installed in tmp."

uncompress-stamp:
	gzip -cd $(abinit_tardir)/$(xmlf90_pkg_name).tar.gz | tar xf -
	touch uncompress-stamp

configure-stamp: uncompress
	-mkdir tmp
	cp build.mk $(xmlf90_pkg_name)/macros/fortran.mk
	touch configure-stamp

build-stamp: configure
	cd $(xmlf90_pkg_name) && FLIB_ROOT="$(abinit_builddir)/plugins/xmlf90/$(xmlf90_pkg_name)/macros" /bin/sh build.sh
	touch build-stamp

install-stamp: build
	cp $(xmlf90_pkg_name)/macros/lib/*.a .
	-cp $(xmlf90_pkg_name)/macros/modules/*.mod .
	touch install-stamp

clean:
	rm -f dummy
