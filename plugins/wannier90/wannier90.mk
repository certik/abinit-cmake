#
# Makefile for the NetCDF library embedded in ABINIT
#

include ../../config.mk

all_targets all wannier90: install
	@echo "$(wannier90_pkg_name) is now ready for use."

uncompress: uncompress-stamp
	@echo "$(wannier90_pkg_name) has been uncompressed."

configure: configure-stamp
	@echo "$(wannier90_pkg_name) has been configured."

build: build-stamp
	@echo "$(wannier90_pkg_name) has been built."

install: install-stamp
	@echo "$(wannier90_pkg_name) has been installed in tmp."

uncompress-stamp:
	gzip -cd $(abinit_tardir)/$(wannier90_pkg_name).tar.gz | tar xf -
	cd $(wannier90_pkg_name) && patch -p1 < $(abinit_srcdir)/plugins/wannier90/$(wannier90_pkg_name)-0001.patch
	touch uncompress-stamp

configure-stamp: uncompress
	-mkdir tmp
	echo "F90 = $(FC)" > $(wannier90_pkg_name)/make.sys
	echo "FCOPTS = $(CPPFLAGS_WANNIER90) $(FCFLAGS_FREEFORM) $(FCFLAGS_WANNIER90)" >> $(wannier90_pkg_name)/make.sys
	echo "LDOPTS = $(CPPFLAGS_WANNIER90) $(FCFLAGS_FREEFORM) $(FCFLAGS_WANNIER90) $(FC_LDFLAGS)" >> $(wannier90_pkg_name)/make.sys
	echo "LIBS = $(lib_linalg_libs) $(FCLIBS_WANNIER90)" >> $(wannier90_pkg_name)/make.sys
	touch configure-stamp

build-stamp: configure
	cd $(wannier90_pkg_name) && $(MAKE) wannier lib
	touch build-stamp

install-stamp: build
	cp $(wannier90_pkg_name)/wannier90.x $(wannier90_pkg_name)/libwannier.a .
	touch install-stamp

clean:
	rm -f dummy
