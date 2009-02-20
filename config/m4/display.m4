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
# Miscellaneous macros
#



# ABI_MSG_END()
# -------------
#
# Prints a message at the end of the configure process.
#
AC_DEFUN([ABI_MSG_END],
[[cat <<EOF

Configuration complete.
You may now type "make" to build ABINIT.
(or "make -j4", "make multi" or "make multi_alt" on a SMP machine)

EOF
]]) # ABI_MSG_END



# ABI_MSG_FC_BUGGY(FC_TYPE)
# -------------------------
#
# Prints a message explaining why a compiler has been wrapped, or giving
# advice for its use.
#
AC_DEFUN([ABI_MSG_FC_BUGGY],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 case "$1" in

  absoft)
   ABI_MSG_NOTICE([fc-absoft],[About the ABSoft Fortran compiler])
   ;;

  ibm)
   ABI_MSG_NOTICE([fc-ibm],[About the IBM XL Fortran compiler])
   ;;

  intel)
   ABI_MSG_NOTICE([fc-intel],[About the Intel Fortran compiler])
   ;;

 esac
]) # ABI_MSG_FC_BUGGY



dnl ABI_MSG_NOTICE(FILE, TITLE)
dnl ---------------------------
dnl
dnl Print a framed message to attract users' attention to something.
dnl
AC_DEFUN([ABI_MSG_NOTICE],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl Init
 abi_msg_file="${abinit_srcdir}/config/messages/$1.msg"
 abi_msg_title="$2"
 test "${abi_msg_title}" = "" && abi_msg_title="IMPORTANT NOTE"

 dnl Format title
 abi_msg_spacer="                                                            "
 abi_msg_tmp1=`echo "${abi_msg_title}" | sed -e 's/./ /g'`
 abi_msg_tmp2=`echo "${abi_msg_tmp1}" | grep "${abi_msg_spacer}"`
 abi_msg_spacer=`echo "${abi_msg_spacer}" | sed -e "s/${abi_msg_tmp1}//"`
 test "${abi_msg_tmp2}" = "" || abi_msg_spacer=""
 abi_msg_title="${abi_msg_title}${abi_msg_spacer}"

 if test -s "${abi_msg_file}"; then

 dnl Print header
 echo ""
 echo "        +--------------------------------------------------------------+"
 echo "        | ${abi_msg_title} |"
 echo "        +--------------------------------------------------------------+"

 dnl Format and write message
 while read abi_msg_line; do
  abi_msg_line=`eval echo ${abi_msg_line}`
  abi_msg_spacer="                                                            "
  abi_msg_tmp1=`echo "${abi_msg_line}" | sed -e 's/./ /g'`
  abi_msg_tmp2=`echo "${abi_msg_tmp1}" | grep "${abi_msg_spacer}"`
  test "${abi_msg_tmp1}" = "" || \
   abi_msg_spacer=`echo "${abi_msg_spacer}" | sed -e "s/${abi_msg_tmp1}//"`
  test "${abi_msg_tmp2}" = "" || abi_msg_spacer=""
  echo "        | ${abi_msg_line}${abi_msg_spacer} |"
 done <"${abi_msg_file}"

 dnl Print footer
 echo "        +--------------------------------------------------------------+"
 echo ""

 else
  AC_MSG_WARN([message file ${abi_msg_file} not found])
 fi
]) dnl ABI_MSG_NOTICE


# ABI_MSG_SECTION(TITLE)
# ----------------------
#
# Prints a nice title for each section.
#
AC_DEFUN([ABI_MSG_SECTION],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 abi_sec_title="$1"

 dnl Calculate amount of space chars needed for pretty-printing
 abi_sec_spaces="                                                                      "
 abi_sec_tmp="${abi_sec_title}"
 while test "${abi_sec_tmp}" != ""; do
  abi_sec_spaces=`echo "${abi_sec_spaces}" | sed -e 's/^.//'`
  abi_sec_tmp=`echo "${abi_sec_tmp}" | sed -e 's/^.//'`
 done

 echo ""
 echo " =============================================================================="
 echo " === ${abi_sec_title}${abi_sec_spaces} ==="
 echo " =============================================================================="
 echo ""
]) # ABI_MSG_SECTION
