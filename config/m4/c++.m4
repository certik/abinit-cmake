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
# C++ compilers support
#



# _ABI_CHECK_CXX_COMPAQ(COMPILER)
# -------------------------------
#
# Checks whether the specified C++ compiler is the COMPAQ C++ compiler.
# If yes, tries to determine its version number and sets the cxx_type
# and cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CXX_COMPAQ],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Compaq C++ compiler])
 cxx_info_string=`$1 -V 2>&1`
 abi_result=`echo "${cxx_info_string}" | grep '^Compaq C++'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  cxx_info_string=""
  cxx_type="UNKNOWN"
  cxx_version="UNKNOWN"
 else
  AC_DEFINE([COMPAQ_CXX],1,[Define to 1 if you are using the COMPAQ C++ compiler])
  cxx_type="compaq"
  cxx_version=`echo "${cxx_info_string}" | grep '^Compiler Driver' | sed -e 's/Compiler Driver V//; s/-.*//'`
  if test "${cxx_version}" = "${cxx_info_string}"; then
   cxx_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CXX_COMPAQ



# _ABI_CHECK_CXX_GCC(COMPILER)
# ----------------------------
#
# Checks whether the specified C++ compiler is the GCC C++ compiler.
# If yes, tries to determine its version number and sets the cxx_type
# and cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CXX_GCC],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the GCC C++ compiler])
 cxx_info_string=`$1 --version 2>&1 | head -1`
 if test "${ac_cv_cxx_compiler_gnu}" != "yes"; then
  cxx_info_string=""
  cxx_type="UNKNOWN"
  cxx_version="UNKNOWN"
  abi_result="no"
 else
  AC_DEFINE([GCC_CXX],1,[Define to 1 if you are using the GCC C++ compiler])
  cxx_type="gcc"
  cxx_version=`echo ${cxx_info_string} | sed -e 's/.*(GCC) //; s/ .*//'`
  if test "${cxx_version}" = "${cxx_info_string}"; then
   abi_result=`echo "${cxx_info_string}" | grep ' '`
   if test "${abi_result}" != ""; then
    cxx_version="UNKNOWN"
   fi
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CXX_GCC



# _ABI_CHECK_CXX_IBM(COMPILER)
# ----------------------------
#
# Checks whether the specified C++ compiler is the IBM XL C++ compiler.
# If yes, tries to determine its version number and sets the cxx_type
# and cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CXX_IBM],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the IBM XL C++ compiler])
 cxx_info_string=`$1 -qversion 2>&1 | head -1`
 cxx_garbage=`$1 -qversion 2>&1 | wc -l | sed -e 's/ //g'`
 abi_result=`echo "${cxx_info_string}" | grep 'IBM(R) XL C/C++'`
 if test "${abi_result}" = ""; then
  abi_result=`echo "${cxx_info_string}" | grep 'C for AIX'`
 fi
 if test "${abi_result}" = ""; then
  abi_result="no"
  cxx_info_string=""
  cxx_type="UNKNOWN"
  cxx_version="UNKNOWN"
  if test "${cxx_garbage}" -gt 50; then
   AC_DEFINE([IBM_CXX],1,[Define to 1 if you are using the IBM XL C++ compiler])
   cxx_type="ibm"
   cxx_version="UNKNOWN"
   abi_result="yes"
  fi
 else
  AC_DEFINE([IBM_CXX],1,[Define to 1 if you are using the IBM XL C++ compiler])
  cxx_type="ibm"
  cxx_version=`echo "${cxx_info_string}" | sed -e 's/.* V//; s/ .*//'`
  if test "${cxx_version}" = "${cxx_info_string}"; then
   cxx_version=`echo "${cxx_info_string}" | sed -e 's/C for AIX version //'`
  fi
  if test "${cxx_version}" = "${cxx_info_string}"; then
   cxx_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CXX_IBM



# _ABI_CHECK_CXX_INTEL(COMPILER)
# ------------------------------
#
# Checks whether the specified C++ compiler is the Intel C++ compiler.
# If yes, tries to determine its version number and sets the cxx_type
# and cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CXX_INTEL],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Intel C++ compiler])
 cxx_info_string=`$1 -v -V 2>&1`
 abi_result=`echo "${cxx_info_string}" | head -n 1 | grep 'Intel(R) C++ Compiler'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  cxx_info_string=""
  cxx_type="UNKNOWN"
  cxx_version="UNKNOWN"
 else
  AC_DEFINE([INTEL_CXX],1,[Define to 1 if you are using the Intel C++ compiler])
  cxx_type="intel"
  cxx_version=`echo "${abi_result}" | sed -e 's/.*Version //; s/ .*//'`
  if test "${cxx_version}" = "${abi_result}"; then
   cxx_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CXX_INTEL



# _ABI_CHECK_CXX_PATHSCALE(COMPILER)
# ----------------------------------
#
# Checks whether the specified C++ compiler is the PathScale C++ compiler.
# If yes, tries to determine its version number and sets the cxx_type
# and cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CXX_PATHSCALE],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the PathScale C++ compiler])
 cxx_info_string=`$1 --version 2>&1`
 abi_result=`echo "${cxx_info_string}" | grep '^PathScale'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  cxx_info_string=""
  cxx_type="UNKNOWN"
  cxx_version="UNKNOWN"
 else
  AC_DEFINE([PATHSCALE_CXX],1,[Define to 1 if you are using the PathScale C++ compiler])
  cxx_type="pathscale"
  cxx_version=`echo "${abi_result}" | sed -e 's/.* Version //; s/ .*//'`
  if test "${cxx_version}" = "${abi_result}"; then
   cxx_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CXX_PATHSCALE



# _ABI_CHECK_CXX_PGI(COMPILER)
# ----------------------------
#
# Checks whether the specified C++ compiler is the Portland Group C++
# compiler. If yes, tries to determine its version number and sets the
# cxx_type and cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CXX_PGI],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Portland Group C++ compiler])
 cxx_info_string=`$1 -v -V 2>&1 | sed -e '/^$/d'`
 abi_result=`echo "${cxx_info_string}" | grep '^pgCC'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  cxx_info_string=""
  cxx_type="UNKNOWN"
  cxx_version="UNKNOWN"
 else
  AC_DEFINE([PGI_CXX],1,[Define to 1 if you are using the Portland Group C++ compiler])
  cxx_type="pgi"
  cxx_version=`echo "${abi_result}" | sed -e 's/.* //; s/-.*//'`
  if test "${cxx_version}" = "${abi_result}"; then
   cxx_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CXX_PGI



# _ABI_CHECK_CXX_SUN(COMPILER)
# ----------------------------
#
# Checks whether the specified C++ compiler is the Sun C++ compiler.
# If yes, tries to determine its version number and sets the
# cxx_type and cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CXX_SUN],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Sun C++ compiler])
 cxx_info_string=`$1 -V 2>&1 | head -n 1`
 abi_result=`echo "${cxx_info_string}" | grep 'Sun' | grep ' C++ '`
 if test "${abi_result}" = ""; then
  abi_result="no"
  cxx_info_string=""
  cxx_type="UNKNOWN"
  cxx_version="UNKNOWN"
 else
  AC_DEFINE([SUN_CXX],1,[Define to 1 if you are using the Sun C++ compiler])
  cxx_type="sun"
  cxx_version=`echo "${abi_result}" | sed -e 's/.* C++ //; s/ .*//'`
  if test "${cxx_version}" = "${abi_result}"; then
   cxx_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CXX_SUN



# ABI_PROG_CXX()
# --------------
#
# Tries to determine which type of C++ compiler is installed.
#
AC_DEFUN([ABI_PROG_CXX],
[dnl Init
 if test "${cxx_type}" = ""; then
  cxx_type="UNKNOWN"
 fi

 dnl Determine C++ compiler type (the order is important)
 AC_MSG_CHECKING([which type of C++ compiler we have])

 dnl Check Intel before GCC, in order not to confuse the detector
 if test "${cxx_type}" = "UNKNOWN"; then
  _ABI_CHECK_CXX_INTEL(${CXX})
 fi
 if test "${cxx_type}" = "UNKNOWN"; then
  _ABI_CHECK_CXX_GCC(${CXX})
 fi
 if test "${cxx_type}" = "UNKNOWN"; then
  _ABI_CHECK_CXX_COMPAQ(${CXX})
 fi
 if test "${cxx_type}" = "UNKNOWN"; then
  _ABI_CHECK_CXX_PATHSCALE(${CXX})
 fi
 if test "${cxx_type}" = "UNKNOWN"; then
  _ABI_CHECK_CXX_PGI(${CXX})
 fi
 if test "${cxx_type}" = "UNKNOWN"; then
  _ABI_CHECK_CXX_SUN(${CXX})
 fi
 if test "${cxx_type}" = "UNKNOWN"; then
  _ABI_CHECK_CXX_IBM(${CXX})
 fi

 dnl Fall back to generic when detection fails
 if test "${cxx_type}" = "UNKNOWN"; then
  cxx_type="generic"
  cxx_version="0.0"
 fi

 dnl Normalise C++ compiler version
 cxx_version=`echo ${cxx_version} | cut -d. -f1-2`

 dnl Display final result
 AC_MSG_RESULT([${cxx_type} ${cxx_version}])

 dnl Schedule compiler info for substitution
 AC_SUBST(cxx_type)
 AC_SUBST(cxx_version)
]) # ABI_PROG_CXX
