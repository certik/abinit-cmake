#
# Makefile for the XML Fortran 90 library embedded in ABINIT
#

include ../../config.mk

all_targets all fox: install
	@echo "$(fox_pkg_name) is now ready for use."

uncompress: uncompress-stamp
	@echo "$(fox_pkg_name) has been uncompressed."

configure: configure-stamp
	@echo "$(fox_pkg_name) has been configured."

build: build-stamp
	@echo "$(fox_pkg_name) has been built."

install: install-stamp
	@echo "$(fox_pkg_name) has been installed in tmp."

uncompress-stamp:
	gzip -cd $(abinit_tardir)/$(fox_pkg_name).tar.gz | tar xf -
	touch uncompress-stamp

configure-stamp: uncompress
	-mkdir tmp
	cd $(fox_pkg_name) && \
	 CPP="$(CPP)" \
	 CPPFLAGS="$(CPPFLAGS_FOX)" \
	 CC="$(CC)" \
	 CFLAGS="$(CFLAGS_FOX)" \
	 CXX="$(CXX)" \
	 CXXFLAGS="$(CXXFLAGS_FOX)" \
	 F77="$(FC)" \
	 FFLAGS="$(FCFLAGS_FIXEDFORM) $(FCFLAGS_FOX)" \
	 F90="$(FC)" \
	 F90FLAGS="$(FCFLAGS_FREEFORM) $(FCFLAGS_FOX)" \
	 FC="$(FC)" \
	 FCFLAGS="$(FCFLAGS_FREEFORM) $(FCFLAGS_FOX)" \
	 ./configure \
	  --prefix=$(PWD)/tmp \
	  --enable-wcml
	touch configure-stamp

build-stamp: configure
	cd $(fox_pkg_name) && $(MAKE)
	touch build-stamp

install-stamp: build
	cd $(fox_pkg_name) && $(MAKE) install
	touch install-stamp

clean:
	rm -f dummy 
