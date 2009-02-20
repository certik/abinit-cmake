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
# C compilers support
#



# _ABI_CHECK_CC_COMPAQ(COMPILER)
# ------------------------------
#
# Checks whether the specified C compiler is the COMPAQ C compiler.
# If yes, tries to determine its version number and sets the cc_type
# and cc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CC_COMPAQ],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Compaq C compiler])
 cc_info_string=`$1 -V 2>&1`
 abi_result=`echo "${cc_info_string}" | grep '^Compaq C '`
 if test "${abi_result}" = ""; then
  abi_result="no"
  cc_info_string=""
  cc_type="UNKNOWN"
  cc_version="UNKNOWN"
 else
  AC_DEFINE([COMPAQ_CC],1,[Define to 1 if you are using the COMPAQ C compiler])
  cc_type="compaq"
  cc_info_string=`$1 -V 2>&1 | grep '^Compiler Driver'`
  cc_version=`echo "${cc_info_string}" | sed -e 's/Compiler Driver V//; s/ .*//'`
  if test "${cc_version}" = "${cc_info_string}"; then
   cc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CC_COMPAQ



# _ABI_CHECK_CC_GCC(COMPILER)
# ---------------------------
#
# Checks whether the specified C compiler is the GCC C compiler.
# If yes, tries to determine its version number and sets the cc_type
# and cc_version variables accordingly.
#
# Note: This macro should be called after AC_PROG_CC.
#
AC_DEFUN([_ABI_CHECK_CC_GCC],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the GCC C compiler])
 cc_info_string=`$1 --version 2>&1 | grep '^gcc' | head -n 1`
 if test "${ac_cv_c_compiler_gnu}" != "yes"; then
  cc_type="UNKNOWN"
  cc_version="UNKNOWN"
  abi_result="no"
 else
  AC_DEFINE([GCC_CC],1,[Define to 1 if you are using the GCC C compiler])
  cc_type="gcc"
  cc_version=`echo ${cc_info_string} | sed -e 's/.*(GCC) //; s/ .*//'`
  if test "${cc_version}" = "${cc_info_string}"; then
   abi_result=`echo "${cc_info_string}" | grep ' '`
   if test "${abi_result}" != ""; then
    cc_version="UNKNOWN"
   fi
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CC_GCC



# _ABI_CHECK_CC_IBM(COMPILER)
# ---------------------------
#
# Checks whether the specified C compiler is the IBM XL C compiler.
# If yes, tries to determine its version number and sets the cc_type
# and cc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CC_IBM],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the IBM XL C compiler])
 cc_info_string=`$1 -qversion 2>&1 | head -n 1`
 cc_garbage=`$1 -qversion 2>&1 | wc -l | sed -e 's/ //g'`
 abi_result=`echo "${cc_info_string}" | grep 'IBM(R) XL C/C++'`
 if test "${abi_result}" = ""; then
  abi_result=`echo "${cc_info_string}" | grep 'C for AIX'`
 fi
 if test "${abi_result}" = ""; then
  abi_result="no"
  cc_info_string=""
  cc_type="UNKNOWN"
  cc_version="UNKNOWN"
  if test "${cc_garbage}" -gt 50; then
   AC_DEFINE([IBM_CC],1,[Define to 1 if you are using the IBM XL C compiler])
   cc_type="ibm"
   cc_version="UNKNOWN"
   abi_result="yes"
  fi
 else
  AC_DEFINE([IBM_CC],1,[Define to 1 if you are using the IBM XL C compiler])
  cc_type="ibm"
  cc_version=`echo "${abi_result}" | sed -e 's/.* V//; s/ .*//'`
  if test "${cc_version}" = "${abi_result}"; then
   cc_version=`echo "${abi_result}" | sed -e 's/C for AIX version //'`
  fi
  if test "${cc_version}" = "${abi_result}"; then
   cc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CC_IBM



# _ABI_CHECK_CC_INTEL(COMPILER)
# -----------------------------
#
# Checks whether the specified C compiler is the Intel C compiler.
# If yes, tries to determine its version number and sets the cc_type
# and cc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CC_INTEL],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Intel C compiler])
 cc_info_string=`$1 -V 2>&1`
 abi_result=`echo "${cc_info_string}" | head -n 1 | grep 'Intel(R) C Compiler'`
 if test "${abi_result}" = ""; then
  abi_result=`echo "${cc_info_string}" | head -n 1 | grep 'Intel(R) C++ Compiler'`
 fi
 if test "${abi_result}" = ""; then
  abi_result="no"
  cc_info_string=""
  cc_type="UNKNOWN"
  cc_version="UNKNOWN"
 else
  AC_DEFINE([INTEL_CC],1,[Define to 1 if you are using the Intel C compiler])
  cc_type="intel"
  cc_version=`echo "${abi_result}" | sed -e 's/.*Version //; s/ .*//'`
  if test "${cc_version}" = "${abi_result}"; then
   cc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CC_INTEL



# _ABI_CHECK_CC_PATHSCALE(COMPILER)
# ---------------------------------
#
# Checks whether the specified C compiler is the PathScale C compiler.
# If yes, tries to determine its version number and sets the cc_type
# and cc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CC_PATHSCALE],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the PathScale C compiler])
 cc_info_string=`$1 --version 2>&1`
 abi_result=`echo "${cc_info_string}" | grep '^PathScale'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  cc_info_string=""
  cc_type="UNKNOWN"
  cc_version="UNKNOWN"
 else
  AC_DEFINE([PATHSCALE_CC],1,[Define to 1 if you are using the PathScale C compiler])
  cc_type="pathscale"
  cc_version=`echo "${abi_result}" | sed -e 's/.* Version //; s/ .*//'`
  if test "${cc_version}" = "${abi_result}"; then
   cc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CC_PATHSCALE



# _ABI_CHECK_CC_PGI(COMPILER)
# ---------------------------
#
# Checks whether the specified C compiler is the Portland Group C compiler.
# If yes, tries to determine its version number and sets the cc_type
# and cc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CC_PGI],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the PGI C compiler])
 cc_info_string=`$1 -V 2>&1 | sed -e '/^$/d'`
 abi_result=`echo "${cc_info_string}" | grep '^pgcc'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  cc_info_string=""
  cc_type="UNKNOWN"
  cc_version="UNKNOWN"
 else
  AC_DEFINE([PGI_CC],1,[Define to 1 if you are using the Portland Group C compiler])
  cc_type="pgi"
  cc_version=`echo "${abi_result}" | sed -e 's/.* //; s/-.*//'`
  if test "${cc_version}" = "${abi_result}"; then
   cc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CC_PGI



# _ABI_CHECK_CC_SUN(COMPILER)
# ---------------------------
#
# Checks whether the specified C compiler is the Sun C compiler.
# If yes, tries to determine its version number and sets the cc_type
# and cc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CC_SUN],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Sun C compiler])
 cc_info_string=`$1 -V 2>&1 | head -n 1`
 abi_result=`echo "${cc_info_string}" | grep 'Sun' | grep ' C '`
 if test "${abi_result}" = ""; then
  abi_result="no"
  cc_info_string=""
  cc_type="UNKNOWN"
  cc_version="UNKNOWN"
 else
  AC_DEFINE([SUN_CC],1,[Define to 1 if you are using the Sun C compiler])
  cc_type="sun"
  cc_version=`echo "${abi_result}" | sed -e 's/.* C //; s/ .*//'`
  if test "${cc_version}" = "${abi_result}"; then
   cc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CC_SUN



# ABI_PROG_CC()
# -------------
#
# Tries to determine which type of C compiler is installed.
#
AC_DEFUN([ABI_PROG_CC],
[dnl Init
 if test "${cc_type}" = ""; then
  cc_type="UNKNOWN"
 fi

 dnl Determine C compiler type (the order is important)
 AC_MSG_CHECKING([which type of compiler we have])

 dnl Check Intel before GCC, in order not to confuse the detector
 if test "${cc_type}" = "UNKNOWN"; then
  _ABI_CHECK_CC_INTEL(${CC})
 fi
 if test "${cc_type}" = "UNKNOWN"; then
  _ABI_CHECK_CC_GCC(${CC})
 fi
 if test "${cc_type}" = "UNKNOWN"; then
  _ABI_CHECK_CC_COMPAQ(${CC})
 fi
 if test "${cc_type}" = "UNKNOWN"; then
  _ABI_CHECK_CC_PATHSCALE(${CC})
 fi
 if test "${cc_type}" = "UNKNOWN"; then
  _ABI_CHECK_CC_PGI(${CC})
 fi
 if test "${cc_type}" = "UNKNOWN"; then
  _ABI_CHECK_CC_SUN(${CC})
 fi
 if test "${cc_type}" = "UNKNOWN"; then
  _ABI_CHECK_CC_IBM(${CC})
 fi

 dnl Fall back to generic when detection fails
 if test "${cc_type}" = "UNKNOWN"; then
  cc_type="generic"
  cc_version="0.0"
 fi

 dnl Normalise C compiler version
 cc_version=`echo ${cc_version} | cut -d. -f1-2`

 dnl Display final result
 AC_MSG_RESULT([${cc_type} ${cc_version}])

 dnl Schedule compiler info for substitution
 AC_SUBST(cc_type)
 AC_SUBST(cc_version)
]) # ABI_PROG_CC
