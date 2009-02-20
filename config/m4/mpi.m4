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
# MPI support for ABINIT
#



# _ABI_MPI_CHECK_GENERIC(PREFIX)
# ------------------------------
#
# Looks for a generic implementation of MPI, using the provided prefix.
#
AC_DEFUN([_ABI_MPI_CHECK_GENERIC],
[dnl Set default values
 mpi_generic_prefix="$1"
 mpi_generic_usable="no"
 mpi_generic_cc=""
 mpi_generic_cxx=""
 mpi_generic_fc=""
 mpi_generic_runner=""
 mpi_generic_cppflags=""
 mpi_generic_cflags=""
 mpi_generic_cxxflags=""
 mpi_generic_fcflags=""
 mpi_generic_ldflags=""

 AC_MSG_CHECKING([for a usable MPI implementation])

 dnl Check whether generic files are available
 if test "${mpi_generic_prefix}" != ""; then

  dnl Look for a library (might not be sufficient)
  if test -s "${mpi_generic_prefix}/lib/libmpi.a"; then
   mpi_generic_ldflags="-L${mpi_generic_prefix}/lib"
  fi

  dnl Look for a Fortran include file
  if test -s "${mpi_generic_prefix}/include/mpif.h"; then
   mpi_generic_cppflags="-I${mpi_generic_prefix}/include"
  fi

  dnl Look for a C compiler
  if test -x "${mpi_generic_prefix}/bin/mpicc"; then
   mpi_generic_cc="${mpi_generic_prefix}/bin/mpicc"
  fi

  dnl Look for C++ compiler
  if test -x "${mpi_generic_prefix}/bin/mpic++"; then
   mpi_generic_cxx="${mpi_generic_prefix}/bin/mpic++"
  fi

  dnl Look for a Fortran 90 compiler
  if test -x "${mpi_generic_prefix}/bin/mpif90"; then
   mpi_generic_fc="${mpi_generic_prefix}/bin/mpif90"
  fi

  dnl Look for a runner
  if test -x "${mpi_generic_prefix}/bin/mpirun"; then
   mpi_generic_runner="${mpi_generic_prefix}/bin/mpirun"
  fi
  if test "${mpi_generic_runner}" = ""; then
   if test -x "${mpi_generic_prefix}/bin/mpiexec"; then
    mpi_generic_runner="${mpi_generic_prefix}/bin/mpiexec"
   fi
  fi
  if test "${mpi_generic_runner}" = ""; then
   if test -x "${mpi_generic_prefix}/bin/dmpirun"; then
    mpi_generic_runner="${mpi_generic_prefix}/bin/dmpirun"
   fi
  fi
  if test "${mpi_generic_runner}" = ""; then
   if test -x "${mpi_generic_prefix}/bin/srun"; then
    mpi_generic_runner="${mpi_generic_prefix}/bin/srun"
   fi
  fi

 fi

 dnl Decide whether generic MPI implementation can be used or not
 if test "${mpi_generic_fc}" != "" -a "${mpi_generic_runner}" != ""; then
  mpi_generic_usable="yes"
 fi
 if test "${mpi_generic_cppflags}" != "" -a "${mpi_generic_ldflags}" != ""; then
  mpi_generic_usable="yes"
 fi

 dnl BEGIN DEBUG
 dnl AC_SUBST(mpi_generic_cc)
 dnl AC_SUBST(mpi_generic_cxx)
 dnl AC_SUBST(mpi_generic_fc)
 dnl AC_SUBST(mpi_generic_runner)
 dnl AC_SUBST(mpi_generic_cppflags)
 dnl AC_SUBST(mpi_generic_cflags)
 dnl AC_SUBST(mpi_generic_cxxflags)
 dnl AC_SUBST(mpi_generic_fcflags)
 dnl AC_SUBST(mpi_generic_ldflags)
 dnl END DEBUG

 AC_MSG_RESULT([${mpi_generic_usable}])
]) # _ABI_MPI_CHECK_GENERIC



# _ABI_MPI_CHECK_NATIVE()
# -----------------------
#
# Looks for a native implementation of MPI, i.e. checks if the Fortran
# compiler is able to produce MPI binaries.
#
AC_DEFUN([_ABI_MPI_CHECK_NATIVE],
[dnl Set default values
 mpi_native_usable="no"
 mpi_native_level=""
 mpi_native_cc=""
 mpi_native_cxx=""
 mpi_native_fc=""
 mpi_native_runner=""

 dnl Back-up build environment
 ABI_ENV_BACKUP

 dnl Try to compile a C MPI program
 AC_MSG_CHECKING([for a native C MPI support])

 CPPFLAGS="${CPPFLAGS} ${MPI_CPPFLAGS}"
 CFLAGS="${CFLAGS} ${MPI_CFLAGS}"
 LDFLAGS="${CC_LDFLAGS} ${MPI_CC_LDFLAGS}"
 LIBS="${CC_LIBS} ${MPI_CC_LIBS}"

 AC_LANG_PUSH([C])
 AC_LINK_IFELSE([AC_LANG_PROGRAM(
  [[#include <stdlib.h>
#include "mpi.h"]],
  [[
      MPI_Init(NULL,NULL);
  ]])], [mpi_native_cc="yes"], [mpi_native_cc="no"])
 AC_LANG_POP

 AC_MSG_RESULT([${mpi_native_cc}])

 dnl Try to compile a C++ MPI program
 AC_MSG_CHECKING([for a native C++ MPI support])

 CXXFLAGS="${CXXFLAGS} ${MPI_CXXFLAGS}"
 LDFLAGS="${CXX_LDFLAGS} ${MPI_CXX_LDFLAGS}"
 LIBS="${CXX_LIBS} ${MPI_CXX_LIBS}"

 AC_LANG_PUSH([C++])
 AC_LINK_IFELSE([AC_LANG_PROGRAM(
  [[@%:@include "mpi.h"]],
  [[
      MPI::Init()
  ]])], [mpi_native_cxx="yes"], [mpi_native_cxx="no"])
 AC_LANG_POP

 AC_MSG_RESULT([${mpi_native_cxx}])

 dnl Try to compile a Fortran MPI program
 AC_MSG_CHECKING([for a native Fortran MPI support])

 FCFLAGS="${FCFLAGS} ${MPI_FCFLAGS}"
 LDFLAGS="${FC_LDFLAGS} ${MPI_FC_LDFLAGS}"
 LIBS="${FC_LIBS} ${MPI_FC_LIBS}"

 AC_LANG_PUSH([Fortran])

 AC_LINK_IFELSE([AC_LANG_PROGRAM([],
  [[  
      include "mpif.h"
      integer :: ierr
      call mpi_init(ierr)
      call mpi_finalize(ierr)
  ]])], [mpi_native_fc="yes"], [mpi_native_fc="no"])

 AC_MSG_RESULT([${mpi_native_fc}])

 dnl Try to compile a Fortran MPI2 program
 if test "${mpi_native_fc}" = "yes"; then
  AC_MSG_CHECKING([for MPI standard level supported])

  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
   [[  
       use mpi
       integer :: ierr
       call mpi_init(ierr)
       call mpi_finalize(ierr)
   ]])], [mpi_native_level="2"], [mpi_native_level="1"])

  AC_MSG_RESULT([MPI-${mpi_native_level}])
 fi

 AC_LANG_POP

 dnl Restore build environment
 ABI_ENV_RESTORE

 dnl Set usability according to results
 if test "${mpi_native_fc}" = "yes"; then
  mpi_native_usable="yes"
 fi

 dnl Look for mpirun
 if test "${mpi_native_usable}" = "yes"; then
  AC_PATH_PROG(mpi_native_runner,mpirun)
  if test "${mpi_native_runner}" = ""; then
   AC_PATH_PROG(mpi_native_runner,mpiexec)
  fi
  if test "${mpi_native_runner}" = ""; then
   AC_PATH_PROG(mpi_native_runner,dmpirun)
  fi
  if test "${mpi_native_runner}" = ""; then
   AC_PATH_PROG(mpi_native_runner,srun)
  fi
 fi

 dnl BEGIN DEBUG
 dnl AC_SUBST(mpi_native_usable)
 dnl AC_SUBST(mpi_native_cc)
 dnl AC_SUBST(mpi_native_cxx)
 dnl AC_SUBST(mpi_native_fc)
 dnl AC_SUBST(mpi_native_runner)
 dnl END DEBUG
]) # _ABI_MPI_CHECK_NATIVE



# ABI_MPI_CHECK()
# ---------------
#
# Tries first to determine which kind of MPI implementation is available,
# then applies tricks and workarounds. The order of search is important,
# as some implementations are currently better supported than others.
#
# Note: ABI_MPI_INIT() should be called first.
#
AC_DEFUN([ABI_MPI_CHECK],
[dnl Check whether MPI support has been forced on or off
 if test "${enable_mpi}" = "yes" -o "${enable_mpi}" = "native"; then

  dnl Look for a native implementation
  if test "${mpi_usable}" != "yes"; then
   _ABI_MPI_CHECK_NATIVE
   mpi_usable="${mpi_native_usable}"
   if test "${mpi_native_usable}" = "yes"; then
    mpi_type="native"
    mpi_level="${mpi_native_level}"
    test "${MPI_RUNNER}" = "" -a "${mpi_native_runner}" != "" && \
     MPI_RUNNER="${mpi_native_runner}"
   fi
  fi

  dnl Enable MPI support if something usable has been found
  if test "${mpi_usable}" = "yes"; then
   enable_mpi="yes"
   if test "${MPI_RUNNER}" = ""; then
    AC_MSG_WARN([MPI runner not found])
   fi
  else
   enable_mpi="no"
  fi
 else
  if test "${enable_mpi}" = "manual"; then
   AC_MSG_NOTICE([using manual settings for MPI configuration])
   enable_mpi="yes"
   mpi_usable="unknown"
   mpi_type="unknown"
  else
   enable_mpi="no"
   MPI_CPPFLAGS=""
   MPI_CFLAGS=""
   MPI_CC_LDFLAGS=""
   MPI_CXXFLAGS=""
   MPI_CXX_LDFLAGS=""
   MPI_FCFLAGS=""
   MPI_FC_LDFLAGS=""
   MPI_RUNNER=""
  fi
 fi
]) # ABI_MPI_CHECK



# ABI_MPI_INIT()
# --------------
#
# Sets MPI default values.
#
AC_DEFUN([ABI_MPI_INIT],
[dnl Set default values
 mpi_usable="no"
 mpi_type="unknown"
 mpi_level="1"

 dnl Init prefix
 if test "${with_mpi_prefix}" != "" -a ! -d "${with_mpi_prefix}"; then
  with_mpi_prefix=""
 fi
 if test "${with_mpi_prefix}" = "yes" -o "${with_mpi_prefix}" = "no"; then
  mpi_prefix=""
 else
  mpi_prefix="${with_mpi_prefix}"
 fi

 dnl Init CPP flags
 if test "${with_mpi_cppflags}" = "yes" -o "${with_mpi_cppflags}" = "no"; then
  MPI_CPPFLAGS=""
 else
  MPI_CPPFLAGS="${with_mpi_cppflags}"
 fi

 dnl Init C flags
 if test "${with_mpi_cflags}" = "yes" -o "${with_mpi_cflags}" = "no"; then
  MPI_CFLAGS=""
 else
  MPI_CFLAGS="${with_mpi_cflags}"
 fi
 if test "${with_mpi_cc_ldflags}" = "yes" -o "${with_mpi_cc_ldflags}" = "no"; then
  MPI_CC_LDFLAGS=""
 else
  MPI_CC_LDFLAGS="${with_mpi_cc_ldflags}"
 fi
 if test "${with_mpi_cc_libs}" = "yes" -o "${with_mpi_cc_libs}" = "no"; then
  MPI_CC_LIBS=""
 else
  MPI_CC_LIBS="${with_mpi_cc_libs}"
 fi

 dnl Init C++ flags
 if test "${with_mpi_cxxflags}" = "yes" -o "${with_mpi_cxxflags}" = "no"; then
  MPI_CXXFLAGS=""
 else
  MPI_CXXFLAGS="${with_mpi_cxxflags}"
 fi
 if test "${with_mpi_cxx_ldflags}" = "yes" -o \
         "${with_mpi_cxx_ldflags}" = "no"; then
  MPI_CXX_LDFLAGS=""
 else
  MPI_CXX_LDFLAGS="${with_mpi_cxx_ldflags}"
 fi
 if test "${with_mpi_cxx_libs}" = "yes" -o "${with_mpi_cxx_libs}" = "no"; then
  MPI_CXX_LIBS=""
 else
  MPI_CXX_LIBS="${with_mpi_cxx_libs}"
 fi

 dnl Init Fortran flags
 if test "${with_mpi_fcflags}" = "yes" -o "${with_mpi_fcflags}" = "no"; then
  MPI_FCFLAGS=""
 else
  MPI_FCFLAGS="${with_mpi_fcflags}"
 fi
 if test "${with_mpi_fc_ldflags}" = "yes" -o \
         "${with_mpi_fc_ldflags}" = "no"; then
  MPI_FC_LDFLAGS=""
 else
  MPI_FC_LDFLAGS="${with_mpi_fc_ldflags}"
 fi
 if test "${with_mpi_fc_libs}" = "yes" -o "${with_mpi_fc_libs}" = "no"; then
  MPI_FC_LIBS=""
 else
  MPI_FC_LIBS="${with_mpi_fc_libs}"
 fi

 dnl Init runner
 if test "${with_mpi_runner}" = "yes" -o "${with_mpi_runner}" = "no"; then
  MPI_RUNNER=""
 else
  MPI_RUNNER="${with_mpi_runner}"
 fi

 dnl Activate MPI if a valid prefix is specified
 if test "${with_mpi_prefix}" != "" -a "${enable_mpi}" = ""; then
  AC_MSG_NOTICE([enabling MPI support])
  enable_mpi="yes"
 fi

 dnl Set MPI support to "no" if still unset
 if test "${enable_mpi}" = ""; then
  AC_MSG_NOTICE([disabling MPI support])
  enable_mpi="no"
 fi

 dnl Report MPI support status
 AC_MSG_CHECKING([for MPI support requested])
 AC_MSG_RESULT([${enable_mpi}])

 dnl Look for a generic implementation
 if test "${enable_mpi}" = "yes" -a "${mpi_usable}" != "yes"; then

  dnl Search within prefix
  if test "${mpi_prefix}" != ""; then
   AC_MSG_NOTICE([using MPI prefix ${mpi_prefix}])
   _ABI_MPI_CHECK_GENERIC([${mpi_prefix}])
  fi

  dnl Fall back into /usr
  if test "${mpi_generic_usable}" != "yes"; then
   AC_MSG_NOTICE([looking for MPI in /usr])
   _ABI_MPI_CHECK_GENERIC([/usr])
   test "${mpi_generic_usable}" = "yes" && mpi_prefix="/usr"
  fi

  dnl Try /usr/local as a last resort
  if test "${mpi_generic_usable}" != "yes"; then
   AC_MSG_NOTICE([looking for MPI in /usr/local])
   _ABI_MPI_CHECK_GENERIC([/usr/local])
   test "${mpi_generic_usable}" = "yes" && mpi_prefix="/usr/local"
  fi

  dnl Propagate results
  if test "${mpi_generic_usable}" = "yes"; then
   mpi_usable="${mpi_generic_usable}"
   mpi_type="generic"
   test "${mpi_generic_cc}" != "" && CC="${mpi_generic_cc}"
   test "${mpi_generic_cxx}" != "" && CXX="${mpi_generic_cxx}"
   test "${mpi_generic_fc}" != "" && FC="${mpi_generic_fc}"
   test "${mpi_generic_runner}" != "" && MPI_RUNNER="${mpi_generic_runner}"

   AC_MSG_CHECKING([the type of MPI implementation we have])

   dnl Check whether LAM files are available
   if test "${mpi_type}" = "generic"; then
    if test -s "${mpi_prefix}/lib/liblam.a"; then
     mpi_type="lam"
     mpi_level="1"
    fi
   fi

   dnl Check whether MPICH files are available
   if test "${mpi_type}" = "generic"; then
    if test -s "${mpi_prefix}/lib/libmpich.a" -a \
            -s "${mpi_prefix}/lib/libfmpich.a"; then
     mpi_type="mpich"
     mpi_level="1"
    fi
   fi

   dnl Check whether Open-MPI files are available
   if test "${mpi_type}" = "generic"; then
    if test -x "${mpi_prefix}/bin/opal_wrapper" -a \
            -x "${mpi_prefix}/bin/orted" -a \
            -x "${mpi_prefix}/bin/orterun"; then
     mpi_type="openmpi"
     mpi_level="2"
    fi
   fi

   AC_MSG_RESULT([${mpi_type}])
  else
   if test "${enable_mpi}" != "no"; then
    AC_MSG_NOTICE([check for a native MPI compiler support deferred])
   fi
  fi
 fi

 dnl Enable substitution
 AC_SUBST(mpi_usable)
 AC_SUBST(mpi_type)
 AC_SUBST(mpi_level)
 AC_SUBST(MPI_CPPFLAGS)
 AC_SUBST(MPI_CFLAGS)
 AC_SUBST(MPI_CC_LDFLAGS)
 AC_SUBST(MPI_CXXFLAGS)
 AC_SUBST(MPI_CXX_LDFLAGS)
 AC_SUBST(MPI_FCFLAGS)
 AC_SUBST(MPI_FC_LDFLAGS)
 AC_SUBST(MPI_RUNNER)
]) # ABI_MPI_INIT
