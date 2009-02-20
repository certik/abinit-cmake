#
# Makefile for the linear algebra library embedded in ABINIT
#

include ../../config.mk

string_f_pkg_name   = string_f-2006
string_f_pkg_string = LibString_F 2006 Octopus SVN checkout
libxc_pkg_name      = libxc-0.9
libxc_pkg_string    = LibXC 0.9 Upstream Release

all_targets all etsf_xc: install
	@echo "$(etsf_xc_pkg_name) is now ready for use."

uncompress: uncompress-stamp
	@echo "$(etsf_xc_pkg_name) has been uncompressed."

configure: configure-stamp
	@echo "$(etsf_xc_pkg_name) has been configured."

build: build-stamp
	@echo "$(etsf_xc_pkg_name) has been built."

install: install-stamp
	@echo "$(etsf_xc_pkg_name) has been installed."

uncompress-stamp:
	gzip -cd $(abinit_tardir)/$(etsf_xc_pkg_name).tar.gz | tar xf -
	gzip -cd $(string_f_pkg_name).tar.gz | tar xf -
	gzip -cd $(libxc_pkg_name).tar.gz | tar xf -
	touch uncompress-stamp

string_f-config: uncompress
	-mkdir tmp
	cd $(string_f_pkg_name) && \
	 CPP="$(CPP)" \
	 CPPFLAGS="$(CPPFLAGS_ETSF_XC)" \
	 CC="$(CC)" \
	 CFLAGS="$(CFLAGS_ETSF_XC)" \
	 ./configure --prefix=$(PWD)/tmp

string_f-build: string_f-config
	cd $(string_f_pkg_name) && $(MAKE)

string_f-install: string_f-build
	cd $(string_f_pkg_name) && $(MAKE) install

build-stamp: configure

configure-stamp: uncompress string_f-install
	cd $(libxc_pkg_name) && \
	 CPP="$(CPP)" \
	 CPPFLAGS="$(CPPFLAGS_ETSF_XC)" \
	 CC="$(CC)" \
	 CFLAGS="$(CFLAGS_ETSF_XC)" \
	 FC="$(FC)" \
	 FCFLAGS="$(FCFLAGS_FREEFORM) $(FCFLAGS_ETSF_XC)" \
	 ./configure \
	  --prefix="$(abinit_builddir)/plugins/etsf_xc/tmp" \
	  --with-string_f="$(abinit_builddir)/plugins/etsf_xc/$(string_f_pkg_name)" \
	  --enable-fortran
	touch configure-stamp

build-stamp: configure
	cd $(libxc_pkg_name) && $(MAKE)
	touch build-stamp

install-stamp: build
	cd $(libxc_pkg_name) && $(MAKE) install
	-cp $(libxc_pkg_name)/src/xc_types.mod \
	 $(libxc_pkg_name)/src/libxc.mod tmp/include
	touch install-stamp

clean:
	rm -rf $(string_f_pkg_name) $(libxc_pkg_name)
	rm -rf $(string_f_pkg_name).tar.gz $(libxc_pkg_name).tar.gz
