#
# Makefile for the BigDFT library embedded in ABINIT
#

include ../../config.mk

all_targets all bigdft: install
	@echo "$(bigdft_pkg_name) is now ready for use."

uncompress: uncompress-stamp
	@echo "$(bigdft_pkg_name) has been uncompressed."

configure: configure-stamp
	@echo "$(bigdft_pkg_name) has been configured."

build: build-stamp
	@echo "$(bigdft_pkg_name) has been built."

install: install-stamp
	@echo "$(bigdft_pkg_name) has been installed in tmp."

uncompress-stamp:
	gzip -cd $(abinit_tardir)/$(bigdft_pkg_name).tar.gz | tar xf -
	touch uncompress-stamp

# In the configure libxc are enable temporary for the old poisson solver
# lib to compile. In will be removed in future version of BigDFT.
configure-stamp: uncompress
	-mkdir tmp
	cd $(bigdft_pkg_name) && \
	 F90="$(FC)" \
	 F90_LDFLAGS="$(FC_LDFLAGS_EXTRA)" \
	 F90FLAGS="$(FCFLAGS_FREEFORM) $(FCFLAGS_BIGDFT)" \
	 FC="$(FC)" \
	 FC_LDFLAGS="$(FC_LDFLAGS_EXTRA)" \
	 FCFLAGS="$(FCFLAGS_FREEFORM) $(FCFLAGS_BIGDFT)"\
	 ./configure \
	  --prefix=$(PWD)/tmp \
	  --enable-mpi="$(enable_mpi)" \
	  --disable-binaries \
	  --disable-libxc \
	  --enable-libpoissonsolver \
	  --with-moduledir=$(PWD)/tmp/include
	touch configure-stamp

build-stamp: configure
	cd $(bigdft_pkg_name) && $(MAKE)
	touch build-stamp

install-stamp: build
	-cd $(bigdft_pkg_name) && $(MAKE) -i install
	touch install-stamp

clean:
	rm -f dummy
