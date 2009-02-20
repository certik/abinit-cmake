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
# Initialization
#



# ABI_INIT_CPU_INFO(TARGET)
# -------------------------
#
# Sets architecture-related variables from the information given by the
# specified target. This is a helper for many other ABINIT macros, that
# should be called quite early in the configure script.
#
# At present, the variables set are:
#
#  * abi_cpu_model  : CPU model, if guessed;
#  * abi_cpu_64bits : whether the CPU is 64 bits or not.
#
AC_DEFUN([ABI_INIT_CPU_INFO],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 abi_cpu_model=""
 abi_cpu_bits=""
 abi_cpu_64bits=""

 case "$1" in

  alpha*|powerpc*)
   abi_cpu_model="${target_cpu}"
   abi_cpu_64bits=`echo "${abi_cpu_model}" | grep '64$'`
   if test "${abi_cpu_64bits}" = ""; then
    abi_cpu_64bits="no"
    abi_cpu_bits="32"
   else
    abi_cpu_64bits="yes"
    abi_cpu_bits="64"
   fi
   ;;

  i686-*linux*)
   dnl Athlon ?
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model=`cat /proc/cpuinfo | grep 'Athlon'`
    if test "${abi_cpu_model}" != ""; then
     abi_cpu_model="athlon"
     abi_cpu_64bits="no"
     abi_cpu_bits="32"
    fi
   fi
   dnl Pentium 3 ?
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model=`cat /proc/cpuinfo | grep 'Pentium III'`
    if test "${abi_cpu_model}" != ""; then
     abi_cpu_model="pentium3"
     abi_cpu_64bits="no"
     abi_cpu_bits="32"
    fi
   fi
   dnl Pentium 4 ?
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) Pentium(R) 4'`
    if test "${abi_cpu_model}" != ""; then
     abi_cpu_model="pentium4"
     abi_cpu_64bits="no"
     abi_cpu_bits="32"
    fi
   fi
   dnl Pentium 4M ?
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) Pentium(R) M'`
    if test "${abi_cpu_model}" != ""; then
     abi_cpu_model="pentium4"
     abi_cpu_64bits="no"
     abi_cpu_bits="32"
    fi
   fi
   dnl Centrino ?
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) CPU           T2400'`
    if test "${abi_cpu_model}" != ""; then
     abi_cpu_model="intel_centrino"
     abi_cpu_64bits="no"
     abi_cpu_bits="32"
    fi
   fi
   dnl Pentium CoreDuo ?
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) CPU           T2050'`
    if test "${abi_cpu_model}" != ""; then
     abi_cpu_model="intel_coreduo"
     abi_cpu_64bits="no"
     abi_cpu_bits="32"
    fi
   fi
   dnl Pentium Core2 ?
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) Core(TM)2 CPU'`
    if test "${abi_cpu_model}" != ""; then
     abi_cpu_model="intel_core2"
     abi_cpu_64bits="no"
     abi_cpu_bits="32"
    fi
   fi
   dnl Unknown
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model="unknown"
    abi_cpu_64bits="unknown"
    abi_cpu_bits="32"
   fi
   ;;

  ia64-*linux*)
   dnl Itanium 1 ?
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model=`cat /proc/cpuinfo | grep 'Itanium' | grep -v 'Itanium 2'`
    if test "${abi_cpu_model}" != ""; then
     abi_cpu_model="itanium1"
    fi
   fi
   dnl Itanium 2 ?
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model=`cat /proc/cpuinfo | grep 'Itanium 2'`
    if test "${abi_cpu_model}" != ""; then
     abi_cpu_model="itanium2"
    fi
   fi
   dnl Unknown
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model="unknown"
   fi
   dnl The processor is anyway 64-bit
   abi_cpu_64bits="yes"
   abi_cpu_bits="64"
   ;;

  mips*irix*)
   # Get processor type
   abi_cpu_model=`hinv 2> /dev/null | grep '^CPU: MIPS '`
   if test "${abi_cpu_model}" != ""; then
    abi_cpu_model=`echo "${abi_cpu_model}" | awk '{print tolower($3)}'`
   fi
   abi_cpu_64bits="yes"
   abi_cpu_bits="64"
   ;;

  x86_64-*linux*)
   dnl Athlon64 ?
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model=`cat /proc/cpuinfo | grep 'Athlon'`
    if test "${abi_cpu_model}" != ""; then
     abi_cpu_model="athlon64"
    fi
   fi
   dnl Opteron ?
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model=`cat /proc/cpuinfo | grep 'Opteron'`
    if test "${abi_cpu_model}" != ""; then
     abi_cpu_model="opteron"
    fi
   fi
   dnl Sempron ?
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model=`cat /proc/cpuinfo | grep 'Sempron'`
    if test "${abi_cpu_model}" != ""; then
     abi_cpu_model="athlon64"
    fi
   fi
   dnl Xeon ?
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) XEON(TM)'`
    if test "${abi_cpu_model}" != ""; then
     abi_cpu_model="xeon"
    fi
   fi
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) Xeon(R)'`
    if test "${abi_cpu_model}" != ""; then
     abi_cpu_model="xeon"
    fi
   fi
   dnl Unknown
   if test "${abi_cpu_model}" = ""; then
    abi_cpu_model="unknown"
   fi
   dnl The processor is anyway 64-bit
   abi_cpu_64bits="yes"
   abi_cpu_bits="64"
   ;;

 esac

 AC_SUBST(abi_cpu_model)
 AC_SUBST(abi_cpu_64bits)
]) # ABI_INIT_CPU_INFO



# ABI_INIT_HEADER()
# -----------------
#
# Initializes the contents of the header file produced by Autoheader.
#
AC_DEFUN([ABI_INIT_HEADER],
[dnl Set top of file ...
 AH_TOP([/*
 * Copyright (c) 2005-2008 ABINIT Group (Yann Pouillon)
 * All rights reserved.
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

/* ABINIT configuration */

#ifndef _ABINIT_CONFIG_H
#define _ABINIT_CONFIG_H

#ifdef __INTEL_COMPILER
#define FC_INTEL 1
#endif

])

 dnl ... as well as bottom
 AH_BOTTOM([#endif /* _ABINIT_CONFIG_H */])
]) # ABI_INIT_HEADER



# ABI_INIT_INSTALL_DIRS()
# -----------------------
#
# Sets installation directories.
#
AC_DEFUN([ABI_INIT_INSTALL_DIRS],
[dnl Set-up prefix
 if test "${prefix}" = "NONE"; then
  abinit_prefix="${ac_default_prefix}"
 else
  abinit_prefix="${prefix}"
 fi

 dnl Set-up all directory names
 abinit_bindir="${abinit_prefix}/abinit/${ABINIT_VERSION_BASE}/bin"
 abinit_chkdir="${abinit_prefix}/abinit/${ABINIT_VERSION_BASE}/tests"
 abinit_datdir="${abinit_prefix}/abinit"
 abinit_docdir="${abinit_prefix}/abinit/${ABINIT_VERSION_BASE}/doc"
 abinit_incdir="${abinit_prefix}/abinit/${ABINIT_VERSION_BASE}/include"
 abinit_libdir="${abinit_prefix}/abinit/${ABINIT_VERSION_BASE}/lib"
 abinit_mandir="${abinit_prefix}/abinit/man"
 abinit_rundir="${abinit_prefix}/abinit/bin"
 abinit_wwwdir="${abinit_prefix}/abinit/${ABINIT_VERSION_BASE}/www"

 dnl Substitute all variables
 AC_SUBST(abinit_prefix)
 AC_SUBST(abinit_bindir)
 AC_SUBST(abinit_chkdir)
 AC_SUBST(abinit_datdir)
 AC_SUBST(abinit_docdir)
 AC_SUBST(abinit_incdir)
 AC_SUBST(abinit_libdir)
 AC_SUBST(abinit_mandir)
 AC_SUBST(abinit_rundir)
 AC_SUBST(abinit_wwwdir)
]) # ABI_INIT_INSTALL_DIRS



# ABI_INIT_PREREQS()
# ------------------
#
# Sets parameters of prerequisite libraries.
#
AC_DEFUN([ABI_INIT_PREREQS],
[dnl LINALG setup
 linalg_supported_types="abinit acml asl atlas cxml essl external mkl mlib sgimath sunperf"
 linalg_pkg_name="lapack-abinit_5.6"
 linalg_pkg_string="An old, robust, version of the Lapack library (hacked by Yann Pouillon)"

 AC_SUBST(linalg_supported_types)
 AC_SUBST(linalg_pkg_name)
 AC_SUBST(linalg_pkg_string)
]) # ABI_INIT_PREREQS



# ABI_INIT_TARGET()
# -----------------
#
# Initializes the target name for the platform ABINIT is about to be built on.
#
# Note: to be called after the detection of the Fortran compiler type.
#
AC_DEFUN([ABI_INIT_TARGET],
[dnl Clean-up operating system name
 [abi_target_os=`echo ${target_os} | sed -e 's/-.*//'`]
 
 ABINIT_TARGET="${target_cpu}_${abi_target_os}_${fc_type}${fc_version}"
 AC_DEFINE_UNQUOTED(ABINIT_TARGET,"${ABINIT_TARGET}",
  [ABINIT target description])
 AC_SUBST(ABINIT_TARGET)
]) # ABI_INIT_TARGET



# ABI_INIT_VERSION()
# ------------------
#
# Sets all variables related to the current version of ABINIT.
#
AC_DEFUN([ABI_INIT_VERSION],
[dnl Get version from Autoconf
 ABINIT_VERSION="${PACKAGE_VERSION}"
 ABINIT_VERSION_MAJOR=`echo "${ABINIT_VERSION}" | cut -d. -s -f1`
 ABINIT_VERSION_MINOR=`echo "${ABINIT_VERSION}" | cut -d. -s -f2`
 ABINIT_VERSION_MICRO=`echo "${ABINIT_VERSION}" | cut -d. -s -f3`
 ABINIT_VERSION_MINOR=`echo "${ABINIT_VERSION_MINOR}" | sed -e 's/[a-z]//g'`
 if test "${ABINIT_VERSION_MICRO}" = ""; then
  ABINIT_VERSION_MICRO=`echo "${ABINIT_VERSION}" | cut -b4-`
 fi
 if test "${ABINIT_VERSION_MICRO}" = ""; then
  ABINIT_VERSION_MICRO="dev"
 fi
 ABINIT_VERSION_BUILD=`date '+%Y%m%d'`

 ABINIT_VERSION_BASE="${ABINIT_VERSION_MAJOR}.${ABINIT_VERSION_MINOR}"

 dnl Make numbers available to source files
 AC_DEFINE_UNQUOTED(ABINIT_VERSION,"${ABINIT_VERSION}",
  [ABINIT whole version number])
 AC_DEFINE_UNQUOTED(ABINIT_VERSION_MAJOR,"${ABINIT_VERSION_MAJOR}",
  [ABINIT major version number])
 AC_DEFINE_UNQUOTED(ABINIT_VERSION_MINOR,"${ABINIT_VERSION_MINOR}",
  [ABINIT minor version number])
 AC_DEFINE_UNQUOTED(ABINIT_VERSION_MICRO,"${ABINIT_VERSION_MICRO}",
  [ABINIT micro version number (patch level)])
 AC_DEFINE_UNQUOTED(ABINIT_VERSION_BUILD,"${ABINIT_VERSION_BUILD}",
  [ABINIT build date])
 AC_DEFINE_UNQUOTED(ABINIT_VERSION_BASE,"${ABINIT_VERSION_BASE}",
  [ABINIT base version number])
 AC_SUBST(ABINIT_VERSION)
 AC_SUBST(ABINIT_VERSION_MAJOR)
 AC_SUBST(ABINIT_VERSION_MINOR)
 AC_SUBST(ABINIT_VERSION_MICRO)
 AC_SUBST(ABINIT_VERSION_BUILD)
 AC_SUBST(ABINIT_VERSION_BASE)
]) # ABI_INIT_VERSION
