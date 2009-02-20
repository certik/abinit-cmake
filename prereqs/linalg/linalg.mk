#
# Makefile for the linear algebra library embedded in ABINIT
#

include ../../config.mk

all_targets all linalg: install
	@echo "$(linalg_pkg_name) is now ready for use."

uncompress: uncompress-stamp
	@echo "$(linalg_pkg_name) has been uncompressed."

configure: configure-stamp
	@echo "$(linalg_pkg_name) has been configured."

build: build-stamp
	@echo "$(linalg_pkg_name) has been built."

install: install-stamp
	@echo "$(linalg_pkg_name) has been installed."

uncompress-stamp:
	gzip -cd $(abinit_srcdir)/prereqs/linalg/$(linalg_pkg_name).tar.gz | tar xf -
	touch uncompress-stamp

configure-stamp: uncompress
	touch configure-stamp

build-stamp: configure
	cd blas && $(MAKE) FC="$(FC)" FCFLAGS="$(FCFLAGS_64BITS) $(FCFLAGS_FIXEDFORM) $(FCFLAGS_LINALG)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" RANLIB="$(RANLIB)"
	cd lapack && $(MAKE) FC="$(FC)" FCFLAGS="$(FCFLAGS_64BITS) $(FCFLAGS_FIXEDFORM) $(FCFLAGS_LINALG)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" RANLIB="$(RANLIB)"
	touch build-stamp

install-stamp: build
	$(MKDIR_P) tmp/lib
	cp blas/libblas.a lapack/liblapack.a tmp/lib
	touch install-stamp

clean:
	rm -rf blas lapack
	rm -f uncompress-stamp configure-stamp build-stamp install-stamp
