#
# Makefile for the ETSF_IO library embedded in ABINIT
#

include ../../config.mk

all_targets all etsf_io: install
	@echo "$(etsf_io_pkg_name) is now ready for use."

uncompress: uncompress-stamp
	@echo "$(etsf_io_pkg_name) has been uncompressed."

configure: configure-stamp
	@echo "$(etsf_io_pkg_name) has been configured."

build: build-stamp
	@echo "$(etsf_io_pkg_name) has been built."

install: install-stamp
	@echo "$(etsf_io_pkg_name) has been installed in tmp."

uncompress-stamp:
	gzip -cd $(abinit_tardir)/$(etsf_io_pkg_name).tar.gz | tar xf -
	touch uncompress-stamp

configure-stamp: uncompress
	-mkdir tmp
	cd $(etsf_io_pkg_name) && \
	 F77="$(FC)" \
	 FFLAGS="$(FCFLAGS_FIXEDFORM) $(FCFLAGS_ETSF_IO)" \
	 F90="$(FC)" \
	 F90FLAGS="$(FCFLAGS_FREEFORM) $(FCFLAGS_ETSF_IO)" \
	 FC="$(FC)" \
	 FCFLAGS="$(FCFLAGS_FREEFORM) $(FCFLAGS_ETSF_IO)" \
	 ./configure --prefix=$(PWD)/tmp --with-netcdf-module-path=`echo "$(lib_netcdf_includes)" | sed "s/-I//"` --with-netcdf-ldflags=`echo "$(lib_netcdf_libs)" | sed "s/-lnetcdf//"` --disable-build-tests --with-moduledir=$(PWD)/tmp/include
	touch configure-stamp

build-stamp: configure
	cd $(etsf_io_pkg_name) && $(MAKE)
	touch build-stamp

install-stamp: build
	cd $(etsf_io_pkg_name) && $(MAKE) -i install
	touch install-stamp

clean:
	rm -f dummy
