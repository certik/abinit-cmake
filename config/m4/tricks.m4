# -*- Autoconf -*-
#
# Copyright (c) 2006-2008 ABINIT Group (Yann Pouillon)
# All rights reserved.
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Tricks for compilers and external libraries
#

 ##############################################################################

#
# Compilers
#



# ABI_TRICKS_AR(SYSTEM)
# ---------------------
#
# Applies archiver tricks and workarounds depending on the operating system.
#
AC_DEFUN([ABI_TRICKS_AR],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl For some mysterious reason, the archiver is not called properly
 dnl anymore. The following line is a workaround.
 test "${ARFLAGS_CMD}" = "" && ARFLAGS_CMD="rc"

 case "$1" in

  aix*)
   ARFLAGS_64BITS="-X 64"
   ;;

 esac
]) # ABI_TRICKS_AR



# ABI_TRICKS_CPP(SYSTEM)
# ----------------------
#
# Applies C preprocessor tricks and workarounds depending on the
# operating system.
#
AC_DEFUN([ABI_TRICKS_CPP],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl Find the true C preprocessor, needed by the wrapper
 AC_MSG_CHECKING([for the true C preprocessor])
 TRUE_CPP=""
 TRUE_CPPFLAGS=""
 if test -x "/lib/cpp"; then
  TRUE_CPP="/lib/cpp"
 fi
 if test "${TRUE_CPP}" = ""; then
  AC_PATH_PROG(TRUE_CPP,cpp,/usr/bin/cpp)
 fi
 AC_MSG_RESULT([${TRUE_CPP}])

 dnl Add command-line options to preprocessor calls
 AC_MSG_CHECKING([for C preprocessor options])
 case "$1" in

  aix*)
   TRUE_CPPFLAGS="-P"
   ;;

  *)
   TRUE_CPPFLAGS="-P -std=c89"
   ;;

 esac
 AC_MSG_RESULT([${TRUE_CPPFLAGS}])
]) # ABI_TRICKS_CPP



# ABI_TRICKS_CC(COMPILER, VERSION)
# ---------------------------------
#
# Applies tricks and workarounds depending on C compiler type and
# version.
#
AC_DEFUN([ABI_TRICKS_CC],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
 m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl

 if test "$1" != "UNKNOWN"; then
  AC_MSG_NOTICE([applying C compiler tricks (type: $1, version: $2)])
 fi

 case "$1" in

  compaq)
   if test "${ac_cv_c_bigendian}" = "yes"; then
    CFLAGS_EXTRA="${CFLAGS_EXTRA} -convert big_endian"
   fi
   ;;

  gnu)
   case "${target}" in
    mips*)
     CFLAGS_64BITS="-mabi=64"
     ;;
    *)
     CFLAGS_64BITS="-m64"
     ;;
   esac
   ;;

  ibm)
   CFLAGS_64BITS="-q64"
   ;;

  intel)
   CC_LDFLAGS_EXTRA="${CC_LDFLAGS_EXTRA} -static-libcxa -i-static"
   case "$2" in 
    9.0|9.1|10.0|10.1)
     CFLAGS_EXTRA="${CFLAGS_EXTRA} -vec-report0"
     ;;
   esac
   ;;

  pathscale)
    CFLAGS_64BITS="-m64"
   ;;

 esac
]) # ABI_TRICKS_CC



# ABI_TRICKS_CXX(COMPILER, VERSION)
# ----------------------------------
#
# Applies tricks and workarounds depending on C++ compiler type and
# version.
#
AC_DEFUN([ABI_TRICKS_CXX],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
 m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl

 if test "$1" != "UNKNOWN"; then
  AC_MSG_NOTICE([applying C++ compiler tricks (type: $1, version: $2)])
 fi

 case "$1" in

  compaq)
   if test "${ac_cv_c_bigendian}" = "yes"; then
    CXXFLAGS_EXTRA="${CXXFLAGS_EXTRA} -convert big_endian"
   fi
   ;;

  gnu)
   case "${target}" in
    mips*)
     CXXFLAGS_64BITS="-mabi=64"
     ;;
    *)
     CXXFLAGS_64BITS="-m64"
     ;;
   esac
   ;;

  ibm)
   CXXFLAGS_64BITS="-q64"
   ;;

  intel)
   CXX_LDFLAGS_EXTRA="${CXX_LDFLAGS_EXTRA} -static-libcxa -i-static"
   case "$2" in 
    9.0|9.1|10.0|10.1)
     CXXFLAGS_EXTRA="${CXXFLAGS_EXTRA} -vec-report0"
     ;;
   esac
   ;;

  pathscale)
   CXXFLAGS_64BITS="-m64"
   ;;

 esac
]) # ABI_TRICKS_CXX



# ABI_TRICKS_FC(COMPILER, VERSION)
# ---------------------------------
#
# Applies tricks and workarounds depending on Fortran compiler type and
# version.
#
AC_DEFUN([ABI_TRICKS_FC],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
 m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl

 if test "$1" != "UNKNOWN"; then
  AC_MSG_NOTICE([applying Fortran compiler tricks (type: $1, version: $2)])
 fi

 case "$1" in

  absoft)
   fc_wrap="yes"
   ;;

  compaq)
   if test "${ac_cv_c_bigendian}" = "yes"; then
    FCFLAGS_EXTRA="${FCFLAGS_EXTRA} -convert big_endian"
   fi
   ;;

  fujitsu)
   FCFLAGS_FREEFORM="-Free"
   FCFLAGS_FIXEDFORM="-Fixed"
   FCFLAGS_EXTRA="${FCFLAGS_EXTRA} -Am -Ee -Ep"
   fc_wrap="yes"
   ;;

  g95|gnu)
   FCFLAGS_64BITS="-m64"
   FCFLAGS_EXTRA="${FCFLAGS_EXTRA} -fno-second-underscore"
   ;;

  hitachi)
   FCFLAGS_EXTRA="-hf95 -nosave -nohugeary"
   ;;

  ibm)
   FCFLAGS_64BITS="-q64"
   FCFLAGS_FREEFORM="-qsuffix=cpp=F90:f=f90 -qfree=f90"
   FCFLAGS_FIXEDFORM="-qsuffix=cpp=F:f=f -qfixed"
   fc_wrap="yes"
   ;;

  intel)
   case "$2" in
    7.*)
     FCFLAGS_EXTRA="${FCFLAGS_EXTRA} -132"
     FC_LDFLAGS_EXTRA="${FC_LDFLAGS_EXTRA} -Vaxlib"
     ;;

    8.0)
     AC_MSG_WARN([IFORT 8.0 is severely buggy and not able to produce reliable ABINIT binaries])
     AC_MSG_ERROR([Please use another version of IFORT, or another compiler],80)
     ;;

    8.1)
     FCFLAGS_EXTRA="${FCFLAGS_EXTRA} -extend_source"
     FC_LDFLAGS_EXTRA="${FC_LDFLAGS_EXTRA} -static-libcxa"
     ;;

    9.0|9.1|10.0|10.1)
     FCFLAGS_EXTRA="${FCFLAGS_EXTRA} -extend_source"
     FC_LDFLAGS_EXTRA="${FC_LDFLAGS_EXTRA} -i-static -static-libcxa"
     case "${build_cpu}" in
      ia64)
       ;;
      *)
       FCFLAGS_EXTRA="${FCFLAGS_EXTRA} -vec-report0"
       ;;
     esac
     ;;

   esac
   ;;

  mipspro)
   FCFLAGS_64BITS="-64"
   if test "${abi_cpu_type}" != ""; then
    FCFLAGS_EXTRA="${FCFLAGS_EXTRA} -${abi_cpu_type}"
   fi
   fc_wrap="yes"
   ;;

  pathscale)
   FCFLAGS_64BITS="-m64"
   FCFLAGS_EXTRA="${FCFLAGS_EXTRA} -extend-source -fno-second-underscore"
   ;;

  pgi)
   FCFLAGS_EXTRA="${FCFLAGS_EXTRA} -Mextend"
   case "$2" in
    6.0)
     FCFLAGS_EXTRA="${FCFLAGS_EXTRA} -Msave"
     ;;
   esac
   ;;

  sun)
   mkdir -p "${abinit_builddir}/tmp-modules"
   FCFLAGS_EXTRA='-moddir=$(abinit_builddir)/tmp-modules -M$(abinit_builddir)/tmp-modules'
   ;;

 esac
]) # ABI_TRICKS_FC



 ##############################################################################

#
# External libraries
#



# ABI_TRICKS_BIGDFT()
# -------------------
#
# Applies tricks and workarounds to have the BigDFT library correctly
# linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_BIGDFT],
[
 AC_MSG_NOTICE([applying BigDFT tricks (not needed yet)])
]) # ABI_TRICKS_BIGDFT



# ABI_TRICKS_ETSF_IO()
# --------------------
#
# Applies tricks and workarounds to have the ETSF I/O library correctly
# linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_ETSF_IO],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 if test "$1" != "UNKNOWN"; then
  AC_MSG_NOTICE([applying ETSF_IO tricks])
 fi

 case "$1" in

  sun)
   lib_netcdf_includes='-I$(abinit_builddir)/tmp-modules'
   ;;

 esac
]) # ABI_TRICKS_ETSF_IO



# ABI_TRICKS_ETSF_XC()
# --------------------
#
# Applies tricks and workarounds to have the ETSF_XC library correctly
# linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_ETSF_XC],
[
 AC_MSG_NOTICE([applying ETSF_XC tricks (not needed yet)])
]) # ABI_TRICKS_ETSF_XC



# ABI_TRICKS_FFTW()
# -----------------
#
# Applies tricks and workarounds to have the FFTW library correctly
# linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_FFTW],
[
 AC_MSG_NOTICE([applying FFTW tricks (not needed yet)])
]) # ABI_TRICKS_FFTW



# ABI_TRICKS_FOX()
# ----------------
#
# Applies tricks and workarounds to have the FoX library correctly
# linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_FOX],
[
 AC_MSG_NOTICE([applying FoX tricks (not needed yet)])
]) # ABI_TRICKS_FOX



# ABI_TRICKS_LINALG(TYPE)
# -----------------------
#
# Applies tricks and workarounds to have the optimized linear algebra
# libraries correctly linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_LINALG],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 if test "$1" != ""; then
  AC_MSG_NOTICE([applying linear algebra tricks (type: $1)])
 fi

 case "$1" in

  abinit)
   linalg_tricks_bypass="no"
   ;;

  atlas)
   FC_LDFLAGS_EXTRA="${FC_LDFLAGS_EXTRA} -llapack -lblas"
   linalg_tricks_bypass="yes"
   ;;

  essl)
   FCFLAGS_EXTRA="${FCFLAGS_EXTRA} -qessl"
   FC_LDFLAGS_EXTRA="${FC_LDFLAGS_EXTRA} -lessl"
   linalg_tricks_bypass="yes"
   ;;

  acml|asl|cxml|mkl|mlib|sgimath|sunperf)
   AC_MSG_WARN([tricks not yet implemented for $1])
   ;;

 esac

 AC_SUBST(linalg_tricks_bypass)
]) # ABI_TRICKS_LINALG



# ABI_TRICKS_NETCDF(COMPILER, VERSION)
# ------------------------------------
#
# Applies tricks and workarounds to have the optimized linear algebra
# libraries correctly linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_NETCDF],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
 m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl

 if test "$1" != "UNKNOWN"; then
  AC_MSG_NOTICE([applying compiler-specific NetCDF tricks (type: $1, version: $2)])
 fi

 CONFIGOPT_NETCDF="${CONFIGOPT_NETCDF} --disable-cxx"

 case "$1" in

  compaq)
   FCFLAGS_NETCDF="${FCFLAGS_NETCDF} -assume no2underscores"
   ;;

  g95)
   CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -Df2cFortran"
   FCFLAGS_NETCDF="${FCFLAGS_NETCDF} -fsecond-underscore"
   ;;

  gcc)
   FCFLAGS_NETCDF="${FCFLAGS_NETCDF} -fsecond-underscore"
   ;;

  ibm)
   CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -DIBMR2Fortran"
   ;;

  intel)
   CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -DpgiFortran"
   case "$2" in
    7.1|8.0|8.1|9.0|9.1)
     FCFLAGS_NETCDF="${FCFLAGS_NETCDF} -mp"
     ;;
    *)
     FCFLAGS_NETCDF="${FCFLAGS_NETCDF} -ip -no-prec-div"
     ;;
   esac
   ;;

  pathscale)
   CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -Df2cFortran"
   FCFLAGS_NETCDF="${FCFLAGS_NETCDF} -fsecond-underscore"
   ;;

  pgi)
   CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -DpgiFortran"
   ;;

  sun)
   CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -DsunFortran"
   ;;

 esac
]) # ABI_TRICKS_NETCDF



# ABI_TRICKS_WANNIER90(COMPILER, VERSION)
# ---------------------------------------
#
# Applies tricks and workarounds to have the Wannier90 bindings correctly
# linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_WANNIER90],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
 m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl

 if test "$1" != "UNKNOWN"; then
  AC_MSG_NOTICE([applying compiler-specific Wannier90 tricks (type: $1, version: $2)])
 fi

 case "$1" in

  intel)
   FCLIBS_WANNIER90="${FCLIBS_WANNIER90} -lsvml"
   ;;

 esac
]) # ABI_TRICKS_WANNIER90



# ABI_TRICKS_XMLF90()
# -------------------
#
# Applies tricks and workarounds to have the XMLF90 library correctly
# linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_XMLF90],
[
 AC_MSG_NOTICE([applying XMLF90 tricks (not needed yet)])
]) # ABI_TRICKS_XMLF90
