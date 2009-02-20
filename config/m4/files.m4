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
# File I/O
#



# ABI_LOAD_OPTIONS()
# ------------------
#
# Looks for options to the configure script in predefined locations:
#
#  * system-wide : /etc/abinit/build/hostname.ac
#  * per-user    : ~/.abinit/build/hostname.ac
#  * local       : <current_directory>/hostname.ac
#
# and eventually in a command-line-specified file. "hostname" is the
# name of the machine without the domain name. The last valid file is
# the one to be read.
#
# NOTE
#
#    The name choosen for the per-user the config file allows the peaceful
#    coexistence of several option sets on machines sharing the same home
#    directory.
#
AC_DEFUN([ABI_LOAD_OPTIONS],
[dnl Setup file names
 abi_hostname=`hostname | sed -e 's/\..*//'`
 abi_sys_options="/etc/abinit/build/${abi_hostname}.ac"
 abi_per_options="${HOME}/.abinit/build/${abi_hostname}.ac"
 abi_src_options="${abinit_srcdir}/${abi_hostname}.ac"
 abi_loc_options="./${abi_hostname}.ac"
 abi_cmd_options="${with_config_file}"
 abi_cfg_options=""

 dnl Select and read config file
 if test "${enable_config_file}" = "yes"; then
  for abi_options in "${abi_sys_options}" "${abi_per_options}" \
                     "${abi_src_options}" "${abi_loc_options}" \
                     "${abi_cmd_options}"; do
   if test -s "${abi_options}"; then
    abi_cfg_options="${abi_options}"
   fi
  done
  if test "${abi_cfg_options}" != ""; then
   AC_MSG_NOTICE([reading options from ${abi_cfg_options}])
   . ${abi_cfg_options}
  fi
 fi
]) # ABI_LOAD_OPTIONS



# ABI_LOAD_OPTFLAGS(SUFFIX, COMPILER, VERSION, ARCHITECTURE)
# ----------------------------------------------------------
#
# Looks for default optimization flags for a specified compiler and
# a specified version, running on a specified architecture. Load them
# if found.
#
AC_DEFUN([ABI_LOAD_OPTFLAGS],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
 m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl
 m4_if([$3], , [AC_FATAL([$0: missing argument 3])])dnl
 m4_if([$4], , [AC_FATAL([$0: missing argument 3])])dnl

 dnl Init
 abi_result=""
 abi_optflags_file=""

 dnl Explore all the possibilities
 for tmp_optflags_file in \
  "${ac_top_srcdir}/config/optflags/generic_$1/all/all.standard" \
  "${ac_top_srcdir}/config/optflags/$2_$1/all/all.standard" \
  "${ac_top_srcdir}/config/optflags/$2_$1/all/$4.standard" \
  "${ac_top_srcdir}/config/optflags/$2_$1/$3/all.standard" \
  "${ac_top_srcdir}/config/optflags/$2_$1/$3/$4.standard" \
  "${ac_top_srcdir}/config/optflags/generic_$1/all/all.${abi_optflags_mode}" \
  "${ac_top_srcdir}/config/optflags/$2_$1/all/all.${abi_optflags_mode}" \
  "${ac_top_srcdir}/config/optflags/$2_$1/all/$4.${abi_optflags_mode}" \
  "${ac_top_srcdir}/config/optflags/$2_$1/$3/all.${abi_optflags_mode}" \
  "${ac_top_srcdir}/config/optflags/$2_$1/$3/$4.${abi_optflags_mode}"; do

  if test -s "${tmp_optflags_file}"; then
   abi_optflags_file="${tmp_optflags_file}"
   abi_result=`echo "${abi_optflags_file}" | \
    sed -e 's,.*optflags/,,; s,^\([[^_]]*\)_\([[^/]]*\),\2: \1,'`
  fi
 done

 dnl Source the file
 #AC_MSG_NOTICE([checking ${abi_optflags_file}])
 if test "${abi_optflags_file}" != ""; then
  AC_MSG_NOTICE([applying optimisations for ${abi_result}])
  . "${abi_optflags_file}"
 else
  AC_MSG_WARN([could not find suitable optimisations])
 fi
]) # ABI_LOAD_OPTFLAGS
