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
# Library detection support
#



# ABI_CHECK_LIBS_BLAS_LAPACK()
# ----------------------------
#
# Checks whether the BLAS and LAPACK libraries are present on specific
# platforms, in order to avoid building them along with ABINIT.
#
AC_DEFUN([ABI_CHECK_LIBS_BLAS_LAPACK],
[case "${target}" in

  *fujitsu*)
    dnl Check for ssl2vp
    AC_CHECK_LIB([ssl2vp],[main],[ssl2vp_lib="ssl2vp"],[ssl2vp_lib=""])
    if test "${ssl2vp_lib}" != ""; then
     AC_MSG_NOTICE([using BLAS/LAPACK routines from SSL2VP library])
     AC_DEFINE(HAVE_SSL2VP,1,[Define to 1 if you have the SSL2VP library])
     HAVE_SSL2VP=1
     blas_lib=""
     blas_libs=""
     blas_includes=""
     lapack_lib=""
     lapack_libs="-lssl2vp"
     lapack_includes=""
     build_blas="no"
     build_lapack="no"
    fi
    ;;

  *irix*)
    dnl Check for COMPLIB.SGIMATH
    AC_CHECK_LIB([complib.sgimath],[main],[complib_sgimath_lib="complib.sgimath"],[complib_sgimath_lib=""])
    if test "${complib_sgimath_lib}" != ""; then
     AC_MSG_NOTICE([using BLAS/LAPACK routines from COMPLIB.SGIMATH library])
     AC_DEFINE(HAVE_COMPLIB_SGIMATH,1,[Define to 1 if you have the COMPLIB.SGIMATH library])
     HAVE_COMPLIB_SGIMATH=1
     blas_lib=""
     blas_libs=""
     blas_includes=""
     lapack_lib=""
     lapack_libs="-lcomplib.sgimath"
     lapack_includes=""
     build_blas="no"
     build_lapack="no"
    fi
    ;;

  *hp*)
    dnl Check for VECLIB
    AC_CHECK_LIB([veclib],[main],[veclib_lib="veclib"],[veclib_lib=""])
    if test "${veclib_lib}" != ""; then
     AC_MSG_NOTICE([using BLAS/LAPACK routines from COMPLIB.SGIMATH library])
     AC_DEFINE(HAVE_VECLIB,1,[Define to 1 if you have the VECLIB library])
     HAVE_VECLIB=1
     blas_lib=""
     blas_libs=""
     blas_includes=""
     lapack_lib=""
     lapack_libs="-lveclib"
     lapack_includes=""
     build_blas="no"
     build_lapack="no"
    fi
    ;;

  *linux*)
    dnl Check for ATLAS
    AC_CHECK_LIB([atlas],[main],[atlas_lib="atlas"],[atlas_lib=""])
    if test "${atlas_lib}" != ""; then
     AC_MSG_NOTICE([using BLAS/LAPACK routines from ATLAS library])
     AC_DEFINE(HAVE_ATLAS,1,[Define to 1 if you have the ATLAS library])
     blas_lib=""
     blas_libs="-lblas"
     blas_includes=""
     lapack_lib=""
     lapack_libs="-llapack"
     lapack_includes=""
     build_blas="no"
     build_lapack="no"
    fi
    dnl Check for MKL
    AC_CHECK_LIB([mkl],[main],[mkl_lib="mkl"],[mkl_lib=""])
    if test "${mkl_lib}" != ""; then
     AC_MSG_NOTICE([using BLAS/LAPACK routines from MKL library])
     AC_DEFINE(HAVE_MKL,1,[Define to 1 if you have the MKL library])
     blas_lib=""
     blas_libs=""
     blas_includes=""
     lapack_lib=""
     lapack_libs="-lmkl"
     lapack_includes=""
     build_blas="no"
     build_lapack="no"
    fi
    ;;

  *nec*)
    dnl Check for ASL
    AC_CHECK_LIB([asl],[main],[asl_lib="asl"],[asl_lib=""])
    if test "${asl_lib}" != ""; then
     AC_MSG_NOTICE([using BLAS/LAPACK routines from ASL library])
     AC_DEFINE(HAVE_ASL,1,[Define to 1 if you have the ASL library])
     HAVE_ASL=1
     blas_lib=""
     blas_libs=""
     blas_includes=""
     lapack_lib=""
     lapack_libs="-lasl"
     lapack_includes=""
     build_blas="no"
     build_lapack="no"
    fi
    ;;

  *sun*)
    dnl Check for MATHLIB
    AC_CHECK_LIB([mathlib],[main],[mathlib_lib="mathlib"],[mathlib_lib=""])
    if test "${mathlib_lib}" != ""; then
     AC_MSG_NOTICE([using BLAS/LAPACK routines from MATHLIB library])
     AC_DEFINE(HAVE_MATHLIB,1,[Define to 1 if you have the MATHLIB library])
     HAVE_MATHLIB=1
     blas_lib=""
     blas_libs=""
     blas_includes=""
     lapack_lib=""
     lapack_libs="-lmathlib"
     lapack_includes=""
     build_blas="no"
     build_lapack="no"
    fi
    ;;

esac
]) # ABI_CHECK_LIBS_BLAS_LAPACK



# ABI_CHECK_LIBS_FFT()
# --------------------
#
# Checks for FFT libraries on specific platforms.
#
# NOTE: This feature is not activated in ABINIT 5.0, and its implementation
#       is still very incomplete.
#
AC_DEFUN([ABI_CHECK_LIBS_FFT],
[dnl Look for some architecture-dependent implementations
 case "${target}" in

  *irix*)
    dnl Check for DFFTW_THREADS
    AC_CHECK_LIB([dfftw_threads],[main],[dfftw_threads_lib="dfftw_threads"],[dfftw_threads_lib=""])
    dnl Check for DFFTW
    AC_CHECK_LIB([dfftw],[main],[dfftw_lib="dfftw"],[dfftw_lib=""])
    if test "${dfftw_threads_lib}" != ""; then
     AC_DEFINE(HAVE_DFFTW_THREADS,1,[Define to 1 if you have the DFFTW_THREADS library])
     HAVE_DFFTW_THREADS=1
     fft_lib="${dfftw_threads_lib}"
     fft_libs="-ldfftw_threads"
     fft_includes=""
     build_fft="no"
    elif test "${dfftw_lib}" != ""; then
     AC_DEFINE(HAVE_DFFTW,1,[Define to 1 if you have the DFFTW library])
     HAVE_DFFTW=1
     fft_lib="${dfftw_lib}"
     fft_libs="-ldfftw"
     fft_includes=""
     build_fft="no"
    fi
    ;;

  *linux*)
   dnl Check for FFTW
   AC_CHECK_LIB([fftw],[main],[fftw_lib="fftw"],[fftw_lib=""])
   dnl Check for RT
   AC_CHECK_LIB([rt],[clock_gettime],[rt_lib="rt"],[rt_lib=""])
   dnl Check for PTHREAD
   AC_CHECK_LIB([pthread],[pthread_initialize],[pthread_lib="-lpthread"],[pthread_lib=""])
   dnl Configure environment
   if test "${fftw_lib}" != ""; then
    AC_DEFINE(HAVE_FFTW,1,[Define to 1 if you have the FFTW library])
    HAVE_FFTW=1
    fft_lib="fftw"
    fft_libs="-lfftw"
    fft_includes=""
    build_fft="no"
    if test "${rt_lib}" != ""; then
     AC_DEFINE(HAVE_RT,1,[Define to 1 if you have the RT library])
     HAVE_RT=1
     fft_libs="${fft_libs} -lrt"
    fi
    if test "${pthread_lib}" != ""; then
     AC_DEFINE(HAVE_PTHREAD,1,[Define to 1 if you have the PTHREAD library])
     HAVE_PTHREAD=1
     fft_libs="${fft_libs} -lpthread"
    fi
   fi
   ;;

esac
]) # ABI_CHECK_LIBS_FFT
