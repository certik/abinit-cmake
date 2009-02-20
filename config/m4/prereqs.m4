# -*- Autoconf -*-
#
# Copyright (c) 2005-2008 ABINIT Group (Yann Pouillon)
# All rights reserved.
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for the libraries required ABINIT
#



# ABI_PREREQ_FFT()
# ----------------
#
# Sets all variables needed to handle an external FFT library.
#
AC_DEFUN([ABI_PREREQ_FFT],
[dnl Initial setup
 lib_fft_includes=""
 lib_fft_libs=""
 
 dnl Define preprocessing options
 if test "${enable_fftw}" = "yes"; then
  AC_MSG_WARN([FFTW support is still under development])
  AC_DEFINE(HAVE_FFTW,1,[Define to 1 if you want to use the FFTW library])

  AC_MSG_CHECKING([whether to use the threaded version of FFTW])
  AC_MSG_RESULT(${enable_fftw_threads})
  if test "${enable_fftw_threads}" = "yes"; then
   AC_DEFINE(HAVE_FFTW_THREADS,1,[Define to 1 if you want to use the threaded FFTW library])
  fi

  lib_fft_includes="${with_fftw_includes}"
  lib_fft_libs="${with_fftw_libs}"
 fi

 dnl Output result
 AC_MSG_CHECKING([whether to use the FFTW library])
 AC_MSG_RESULT(${enable_fftw})

 dnl Substitute variables needed for the use of the library
 AC_SUBST(lib_fft_includes)
 AC_SUBST(lib_fft_libs)
]) # ABI_PREREQ_FFT



# ABI_PREREQ_LINALG()
# -------------------
#
# Sets all variables needed to handle the LINALG external library.
#
AC_DEFUN([ABI_PREREQ_LINALG],
[dnl Initial setup
 lib_linalg_includes=""
 lib_linalg_libs=""

 dnl Define variables needed to build the library
 if test -z "${CPPFLAGS_LINALG}"; then
  CPPFLAGS_LINALG="${CPPFLAGS}"
 fi
 AC_SUBST(CPPFLAGS_LINALG)
 if test -z "${CFLAGS_LINALG}"; then
  CFLAGS_LINALG="${CFLAGS}"
 fi
 AC_SUBST(CFLAGS_LINALG)
 if test -z "${CXXFLAGS_LINALG}"; then
  CXXFLAGS_LINALG="${CXXFLAGS}"
 fi
 AC_SUBST(CXXFLAGS_LINALG)
 if test -z "${FCFLAGS_LINALG}"; then
  FCFLAGS_LINALG="${FCFLAGS}"
 fi
 AC_SUBST(FCFLAGS_LINALG)

 dnl Add optimizations for Fortran flags
 FCFLAGS_LINALG="${FCFLAGS_LINALG} ${fcflags_opt_linalg}"

 dnl Set type from command-line option
 linalg_type="${with_linalg_type}"

 dnl Check whether library option has been specified
 if test "${with_linalg_libs}" = ""; then
  lib_linalg_includes=""
  lib_linalg_libs="-L\$(abinit_builddir)/prereqs/linalg -llapack -lblas"
  build_linalg="yes"
  test "${linalg_type}" = "" && linalg_type="abinit"
 else
  lib_linalg_includes="${with_linalg_includes}"
  lib_linalg_libs="${with_linalg_libs}"
  build_linalg="no"
  test "${linalg_type}" = "" && linalg_type="external"
 fi

 dnl Apply tricks if wanted
 if test "${enable_tricks}" = "yes"; then
  ABI_TRICKS_LINALG(${linalg_type})
  if test "${linalg_tricks_bypass}" = "yes"; then
   build_linalg="no"
  fi
 fi

 dnl Package information to export
 dnl linalg_pkg_name="@PKG_NAME@"
 dnl linalg_pkg_string="@PKG_DESC@"

 dnl Output result
 AC_MSG_CHECKING([whether to build the LINALG library])
 AC_MSG_RESULT([${build_linalg}])

 dnl Substitute variables needed to build the library
 dnl AC_SUBST(linalg_pkg_name)
 dnl AC_SUBST(linalg_pkg_string)

 dnl Substitute variables needed for the use of the library
 AC_SUBST(lib_linalg_includes)
 AC_SUBST(lib_linalg_libs)
 AC_SUBST(build_linalg)
 AC_SUBST(linalg_type)

 dnl Inform Automake
 AM_CONDITIONAL(DO_BUILD_LINALG,test "${build_linalg}" = "yes")
]) # ABI_PREREQ_LINALG
