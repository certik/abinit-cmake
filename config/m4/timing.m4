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
# Fortran timing subroutines
#



# _ABI_CHECK_FORTRAN_ETIME()
# --------------------------
#
# Checks whether the Fortran compiler supports the etime() subroutine.
#
AC_DEFUN([_ABI_CHECK_FORTRAN_ETIME],
[dnl Init
 fc_has_etime="no"

 AC_MSG_CHECKING([whether the Fortran compiler accepts etime()])

 dnl Try to compile a program calling etime()
 tmp_fc_srcext="${FC_SRCEXT}"
 FC_SRCEXT="F90"
 AC_LANG_PUSH([Fortran])
 AC_LINK_IFELSE([AC_LANG_PROGRAM([],
  [[
      call etime(1)
  ]])], [fc_has_etime="yes"])
 AC_LANG_POP()
 FC_SRCEXT="${tmp_fc_srcext}"
 unset tmp_fc_srcext

 if test "${fc_has_etime}" = "yes"; then
  AC_DEFINE([HAVE_FORTRAN_ETIME],1,
   [Define to 1 if your Fortran compiler supports etime()])
 fi

 AC_MSG_RESULT(${fc_has_etime})
]) # _ABI_CHECK_FORTRAN_ETIME



# ABI_FC_TIMING()
# ---------------
#
# Tries to determine which Fortran timing routines are available.
#
AC_DEFUN([ABI_FC_TIMING],
[dnl Init
 fc_timing="standard"

 dnl Look for etime() support
 if test "${fc_timing}" = "standard"; then
  _ABI_CHECK_FORTRAN_ETIME
 fi

 dnl Schedule info for substitution
 AC_SUBST(fc_timing)
]) # ABI_FC_TIMING
