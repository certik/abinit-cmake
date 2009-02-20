#
# Makefile for the NetCDF library embedded in ABINIT
#

include ../../config.mk

all_targets all netcdf: install
	@echo "$(netcdf_pkg_name) is now ready for use."

uncompress: uncompress-stamp
	@echo "$(netcdf_pkg_name) has been uncompressed."

configure: configure-stamp
	@echo "$(netcdf_pkg_name) has been configured."

build: build-stamp
	@echo "$(netcdf_pkg_name) has been built."

install: install-stamp
	@echo "$(netcdf_pkg_name) has been installed in tmp."

uncompress-stamp:
	gzip -cd $(abinit_tardir)/$(netcdf_pkg_name).tar.gz | tar xf -
	touch uncompress-stamp

configure-stamp: uncompress
	-mkdir tmp
	cd $(netcdf_pkg_name) && \
	 CPP="$(CPP)" \
	 CPPFLAGS="$(CPPFLAGS_NETCDF)" \
	 CC="$(CC)" \
	 CFLAGS="$(CFLAGS_NETCDF)" \
	 CXX="$(CXX)" \
	 CXXFLAGS="$(CXXFLAGS_NETCDF)" \
	 F77="$(FC)" \
	 FFLAGS="$(FCFLAGS_FIXEDFORM) $(FCFLAGS_NETCDF)" \
	 F90="$(FC)" \
	 F90FLAGS="$(FCFLAGS_FREEFORM) $(FCFLAGS_NETCDF)" \
	 FC="$(FC)" \
	 FCFLAGS="$(FCFLAGS_FREEFORM) $(FCFLAGS_NETCDF)" \
	 ./configure \
	  --prefix=$(PWD)/tmp \
	  $(CONFIGOPT_NETCDF)
	touch configure-stamp

build-stamp: configure
	cd $(netcdf_pkg_name) && $(MAKE)
	touch build-stamp

install-stamp: build
	cd $(netcdf_pkg_name) && $(MAKE) -i install
	touch install-stamp

clean:
	rm -f dummy
