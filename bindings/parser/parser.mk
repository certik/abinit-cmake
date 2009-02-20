#
# Makefile for the parser bindings of ABINIT
#

include ../../config.mk

PYTHON_CPPFLAGS = `python-config --includes`

all: libabinis.a python

absources:
	cd $(top_builddir)/plugins/bigdft && $(MAKE)
	cd $(top_builddir)/src/defs && $(MAKE)
	cd $(top_builddir)/src/lib00numeric && $(MAKE)
	cd $(top_builddir)/src/01manage_mpi && $(MAKE)
	cd $(top_builddir)/src/11util && $(MAKE)
	cd $(top_builddir)/src/12parser && $(MAKE)
	cd $(top_builddir)/src/12geometry && $(MAKE)
	cd $(top_builddir)/src/13xml && $(MAKE)
	cd $(top_builddir)/src/13recipspace && $(MAKE)
	cd $(top_builddir)/src/13iovars && $(MAKE)

libabinit_tmpdir = tmp-libabinis-objects
libabinis.a: absources libbindings.a 
	test -e "$(libabinit_tmpdir)" || $(MKDIR_P) $(libabinit_tmpdir)
	cd $(libabinit_tmpdir) && $(AR) x ../$(top_builddir)/src/lib00numeric/liblib00numeric.a
	cd $(libabinit_tmpdir) && $(AR) x ../$(top_builddir)/src/defs/libdefs.a
	cd $(libabinit_tmpdir) && $(AR) x ../$(top_builddir)/src/01manage_mpi/lib01manage_mpis.a
	cd $(libabinit_tmpdir) && $(AR) x ../$(top_builddir)/src/11util/lib11util.a
	cd $(libabinit_tmpdir) && $(AR) x ../$(top_builddir)/src/12parser/lib12parser.a
	cd $(libabinit_tmpdir) && $(AR) x ../$(top_builddir)/src/12geometry/lib12geometry.a
	cd $(libabinit_tmpdir) && $(AR) x ../$(top_builddir)/src/13xml/lib13xml.a
	cd $(libabinit_tmpdir) && $(AR) x ../$(top_builddir)/src/13recipspace/lib13recipspace.a
	cd $(libabinit_tmpdir) && $(AR) x ../$(top_builddir)/src/13iovars/lib13iovars.a
	cd $(libabinit_tmpdir) && $(AR) x ../$(top_builddir)/bindings/parser/libbindings.a
	$(AR) $(ARFLAGS) libabinis.a $(libabinit_tmpdir)/*.$(OBJEXT)
	rm -rf $(libabinit_tmpdir)


dtset_c.h dtset_c.static.h dtset_py.h dtset_f90.inc ab_dtset_f90_get.f90: dtset.pickle parse_dtset.py
	$(srcdir)/parse_dtset.py

ab_dtset_c.o: ab_dtset_c.c ab_dtset_c.h dtset_c.h dtset_c.static.h
	$(CC) $(CPPFLAGS) -I. $(CFLAGS) -c $(srcdir)/ab_dtset_c.c

ab_dtset_f90.o: ab_dtset_f90.f90 dtset_f90.inc ab_dtset_f90_get.f90
	$(FC) $(CPPFLAGS) -I. -I$(top_builddir)/src/defs $(FCFLAGS) -c $(srcdir)/ab_dtset_f90.f90

dtsetinit.o: dtsetinit.F90
	$(FC) $(CPPFLAGS) -I. -I$(top_builddir)/src/defs $(FCFLAGS) -c $(srcdir)/dtsetinit.F90

libbindings.a: ab_dtset_c.o ab_dtset_f90.o dtsetinit.o
	$(AR) $(ARFLAGS) libbindings.a ab_dtset_c.o ab_dtset_f90.o dtsetinit.o
	$(RANLIB) libbindings.a

clean:
	rm -f *.o *.mod libbindings.a libabinis.a
	rm -f abinit.so
	rm -f dtset_c.h dtset_c.static.h dtset_py.h dtset_f90.inc ab_dtset_f90_get.f90


ab_dtset_py.o: dtset_py.h ab_dtset_c.h
	$(CC) $(CPPFLAGS) $(PYTHON_CPPFLAGS) -I. $(CFLAGS) -c $(srcdir)/ab_dtset_py.c


python: ab_dtset_py.o
	$(FC) -shared -o abinit.so ab_dtset_py.o -L. -labinis



examples: example-c example-f90

example-c: libabinis.a example-c.c fallbacks.o
	$(FC) -o example-c -I. $(srcdir)/example-c.c libabinis.a fallbacks.o

example-f90: libabinis.a example-f90.f90 fallbacks.o
	$(FC) -I. -I$(top_builddir)/src/defs $(FCFLAGS) -o example-f90 $(srcdir)/example-f90.f90 libabinis.a fallbacks.o

fallbacks.o: fallbacks.f90
	$(FC) -I. -I$(top_builddir)/src/defs $(FCFLAGS) -c $(srcdir)/fallbacks.f90
