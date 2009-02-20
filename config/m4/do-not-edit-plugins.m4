# -*- Autoconf -*-
#
# Copyright (c) 2005-2008 ABINIT Group (Yann Pouillon)
# All rights reserved.
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

# Generated by make-macros-plugins on 2009/02/14 09:58:06 +0000

#
# ABINIT plug-in support for the "configure" script
#

#
# IMPORTANT NOTE
#
# This file has been automatically generated by the make-macros-plugins
# script. If you try to edit it, your changes will systematically be
# overwritten.
#



# ABI_PLUGIN_BIGDFT()
# -------------------
#
# Sets all variables needed to handle the BIGDFT plug-in.
#
AC_DEFUN([ABI_PLUGIN_BIGDFT],
[dnl Initial setup
 lib_bigdft_includes=""
 lib_bigdft_libs=""
 build_bigdft="no"
 abi_plug_ready="no"
 abi_plug_tarball="no"
 abi_plug_bins=""
 abi_plug_incs=""
 abi_plug_libs=""

 dnl Define variables needed to build the library
 bigdft_pkg_name="bigdft-1.0-1"
 AC_SUBST(bigdft_pkg_name)
 bigdft_pkg_string="BigDFT library 1.0-1 (hacked by D. Caliste to add up-to-date module support)"
 AC_SUBST(bigdft_pkg_string)

 if test -z "${CPPFLAGS_BIGDFT}"; then
  CPPFLAGS_BIGDFT="${CPPFLAGS}"
 fi
 AC_SUBST(CPPFLAGS_BIGDFT)
 if test -z "${CFLAGS_BIGDFT}"; then
  CFLAGS_BIGDFT="${CFLAGS}"
 fi
 AC_SUBST(CFLAGS_BIGDFT)
 if test -z "${CXXFLAGS_BIGDFT}"; then
  CXXFLAGS_BIGDFT="${CXXFLAGS}"
 fi
 AC_SUBST(CXXFLAGS_BIGDFT)
 if test -z "${FCFLAGS_BIGDFT}"; then
  FCFLAGS_BIGDFT="${FCFLAGS}"
 fi
 AC_SUBST(FCFLAGS_BIGDFT)

 dnl Add optimizations for Fortran flags
 FCFLAGS_BIGDFT="${FCFLAGS_BIGDFT} ${fcflags_opt_bigdft}"

 dnl Check whether to activate plug-in
 if test "${enable_bigdft}" = "yes"; then

  dnl Check for already installed plug-ins
  if test "${with_plugins_prefix}" != ""; then
   AC_MSG_CHECKING([whether BIGDFT is ready for use])
   abi_plug_ready="yes"
   abi_plug_bins="${with_plugins_prefix}/bin"
   abi_plug_incs="-I${with_plugins_prefix}/include"
   abi_plug_libs="-L${with_plugins_prefix}/lib"
   
   if test -s "${with_plugins_prefix}/lib/libbigdft.a"; then
    abi_plug_libs="${abi_plug_libs} -lbigdft"
   else
    abi_plug_ready="no"
   fi

   if test -s "${with_plugins_prefix}/lib/libpoissonsolver.a"; then
    abi_plug_libs="${abi_plug_libs} -lpoissonsolver"
   else
    abi_plug_ready="no"
   fi

   AC_MSG_RESULT([${abi_plug_ready}])
  fi

  if test "${abi_plug_ready}" = "yes"; then
   with_bigdft_bins="${abi_plug_bins}"
   with_bigdft_includes="${abi_plug_incs}"
   with_bigdft_libs="${abi_plug_libs}"
  fi

  dnl Check whether command-line plug-in options have been specified
  if test "${with_bigdft_includes}" = "" -o "${with_bigdft_libs}" = ""; then

   dnl Check for a tarball repository
   AC_MSG_CHECKING([for a source tarball of BIGDFT])
   if test -s "${abinit_tardir}/bigdft-1.0-1.tar.gz"; then
    abi_plug_tarball="yes"
   fi
   AC_MSG_RESULT([${abi_plug_tarball}])

   dnl Get the package
   if test "${abi_plug_ready}" = "no"; then
    if test "${abi_plug_tarball}" = "no"; then
     AC_MSG_NOTICE([downloading BIGDFT - this may take a while])

     if test ! -s "${abinit_tardir}/bigdft-1.0-1.tar.gz"; then
      AC_MSG_CHECKING([availability of BIGDFT from URL 1])
      ${WGET} --timeout=15 --tries=1 -q -O \
       "${abinit_tardir}/bigdft-1.0-1.tar.gz" \
       'http://www.abinit.org/plugins/bigdft-1.0-1.tar.gz'
      test -s "${abinit_tardir}/bigdft-1.0-1.tar.gz" && abi_plug_tarball="yes"
      AC_MSG_RESULT([${abi_plug_tarball}])
     fi

    fi

    dnl Enable plug-in support only if the download was successful
    if test "${abi_plug_tarball}" = "yes"; then
     lib_bigdft_includes="-I\$(abinit_builddir)/plugins/bigdft"
     lib_bigdft_libs="-L\$(abinit_builddir)/plugins/bigdft -lbigdft -lpoissonsolver"
     build_bigdft="yes"
    else
     lib_bigdft_includes=""
     lib_bigdft_libs=""
     enable_bigdft="no"
     build_bigdft="no"
     AC_MSG_WARN([could not download BIGDFT plug-in tarball])
     AC_MSG_WARN([support for BIGDFT plug-in has been disabled])
    fi
   fi
  else
   lib_bigdft_includes="${with_bigdft_includes}"
   lib_bigdft_libs="${with_bigdft_libs}"
   build_bigdft="no"
  fi

 else
  enable_bigdft="no"
  build_bigdft="no"
 fi

 dnl Set variables required to build binaries
  dnl No binary for BIGDFT

 dnl Output results
 AC_MSG_CHECKING([whether to enable the BIGDFT plug-in])
 AC_MSG_RESULT([${enable_bigdft}])
 AC_MSG_CHECKING([whether to build the BIGDFT plug-in])
 AC_MSG_RESULT([${build_bigdft}])

 dnl Apply tricks if wanted
 if test "${build_bigdft}" = "yes" -a "${enable_tricks}" = "yes"; then
  ABI_TRICKS_BIGDFT(${fc_type},${fc_version})
  if test "${bigdft_tricks_bypass}" = "yes"; then
   build_bigdft="no"
  fi
 fi

 dnl Define preprocessing macro
 if test "${enable_bigdft}" = "yes"; then
  AC_DEFINE([HAVE_BIGDFT],1,[Define to 1 if you want support for BIGDFT])
 fi

 dnl Substitute variables needed for the use of the plug-in
 AC_SUBST(lib_bigdft_includes)
 AC_SUBST(lib_bigdft_libs)
 AC_SUBST(build_bigdft)

 dnl Inform Automake
 AM_CONDITIONAL(DO_ENABLE_BIGDFT,test "${enable_bigdft}" = "yes")
 AM_CONDITIONAL(DO_BUILD_BIGDFT,test "${build_bigdft}" = "yes")
]) # ABI_PLUGIN_BIGDFT



# ABI_PLUGIN_NETCDF()
# -------------------
#
# Sets all variables needed to handle the NETCDF plug-in.
#
AC_DEFUN([ABI_PLUGIN_NETCDF],
[dnl Initial setup
 lib_netcdf_includes=""
 lib_netcdf_libs=""
 build_netcdf="no"
 abi_plug_ready="no"
 abi_plug_tarball="no"
 abi_plug_bins=""
 abi_plug_incs=""
 abi_plug_libs=""

 dnl Define variables needed to build the library
 netcdf_pkg_name="netcdf-3.6.3"
 AC_SUBST(netcdf_pkg_name)
 netcdf_pkg_string="NetCDF library 3.6.3 (upstream release)"
 AC_SUBST(netcdf_pkg_string)

 if test -z "${CPPFLAGS_NETCDF}"; then
  CPPFLAGS_NETCDF="${CPPFLAGS}"
 fi
 AC_SUBST(CPPFLAGS_NETCDF)
 if test -z "${CFLAGS_NETCDF}"; then
  CFLAGS_NETCDF="${CFLAGS}"
 fi
 AC_SUBST(CFLAGS_NETCDF)
 if test -z "${CXXFLAGS_NETCDF}"; then
  CXXFLAGS_NETCDF="${CXXFLAGS}"
 fi
 AC_SUBST(CXXFLAGS_NETCDF)
 if test -z "${FCFLAGS_NETCDF}"; then
  FCFLAGS_NETCDF="${FCFLAGS}"
 fi
 AC_SUBST(FCFLAGS_NETCDF)

 dnl Add optimizations for Fortran flags
 FCFLAGS_NETCDF="${FCFLAGS_NETCDF} ${fcflags_opt_netcdf}"

 dnl Check whether to activate plug-in
 if test "${enable_netcdf}" = "yes"; then

  dnl Check for already installed plug-ins
  if test "${with_plugins_prefix}" != ""; then
   AC_MSG_CHECKING([whether NETCDF is ready for use])
   abi_plug_ready="yes"
   abi_plug_bins="${with_plugins_prefix}/bin"
   abi_plug_incs="-I${with_plugins_prefix}/include"
   abi_plug_libs="-L${with_plugins_prefix}/lib"
   
   if test ! -x "${with_plugins_prefix}/bin/ncdump"; then
    abi_plug_ready="no"
   fi

   if test ! -x "${with_plugins_prefix}/bin/ncgen"; then
    abi_plug_ready="no"
   fi

   if test -s "${with_plugins_prefix}/lib/libnetcdf.a"; then
    abi_plug_libs="${abi_plug_libs} -lnetcdf"
   else
    abi_plug_ready="no"
   fi

   AC_MSG_RESULT([${abi_plug_ready}])
  fi

  if test "${abi_plug_ready}" = "yes"; then
   with_netcdf_bins="${abi_plug_bins}"
   with_netcdf_includes="${abi_plug_incs}"
   with_netcdf_libs="${abi_plug_libs}"
  fi

  dnl Check whether command-line plug-in options have been specified
  if test "${with_netcdf_includes}" = "" -o "${with_netcdf_libs}" = ""; then

   dnl Check for a tarball repository
   AC_MSG_CHECKING([for a source tarball of NETCDF])
   if test -s "${abinit_tardir}/netcdf-3.6.3.tar.gz"; then
    abi_plug_tarball="yes"
   fi
   AC_MSG_RESULT([${abi_plug_tarball}])

   dnl Get the package
   if test "${abi_plug_ready}" = "no"; then
    if test "${abi_plug_tarball}" = "no"; then
     AC_MSG_NOTICE([downloading NETCDF - this may take a while])

     if test ! -s "${abinit_tardir}/netcdf-3.6.3.tar.gz"; then
      AC_MSG_CHECKING([availability of NETCDF from URL 1])
      ${WGET} --timeout=15 --tries=1 -q -O \
       "${abinit_tardir}/netcdf-3.6.3.tar.gz" \
       'http://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-3.6.3.tar.gz'
      test -s "${abinit_tardir}/netcdf-3.6.3.tar.gz" && abi_plug_tarball="yes"
      AC_MSG_RESULT([${abi_plug_tarball}])
     fi
     if test ! -s "${abinit_tardir}/netcdf-3.6.3.tar.gz"; then
      AC_MSG_CHECKING([availability of NETCDF from URL 2])
      ${WGET} --timeout=15 --tries=1 -q -O \
       "${abinit_tardir}/netcdf-3.6.3.tar.gz" \
       'http://www.abinit.org/plugins/netcdf-3.6.3.tar.gz'
      test -s "${abinit_tardir}/netcdf-3.6.3.tar.gz" && abi_plug_tarball="yes"
      AC_MSG_RESULT([${abi_plug_tarball}])
     fi

    fi

    dnl Enable plug-in support only if the download was successful
    if test "${abi_plug_tarball}" = "yes"; then
     lib_netcdf_includes="-I\$(abinit_builddir)/plugins/netcdf"
     lib_netcdf_libs="-L\$(abinit_builddir)/plugins/netcdf -lnetcdf"
     build_netcdf="yes"
    else
     lib_netcdf_includes=""
     lib_netcdf_libs=""
     enable_netcdf="no"
     build_netcdf="no"
     AC_MSG_WARN([could not download NETCDF plug-in tarball])
     AC_MSG_WARN([support for NETCDF plug-in has been disabled])
    fi
   fi
  else
   lib_netcdf_includes="${with_netcdf_includes}"
   lib_netcdf_libs="${with_netcdf_libs}"
   build_netcdf="no"
  fi

 else
  enable_netcdf="no"
  build_netcdf="no"
 fi

 dnl Set variables required to build binaries
  if test "${enable_netcdf}" = "yes"; then
  FCLIBS_NETCDF='$(FCLIBS)'
 else
  FCLIBS_NETCDF=''
 fi
 AC_SUBST(FCLIBS_NETCDF)

 dnl Output results
 AC_MSG_CHECKING([whether to enable the NETCDF plug-in])
 AC_MSG_RESULT([${enable_netcdf}])
 AC_MSG_CHECKING([whether to build the NETCDF plug-in])
 AC_MSG_RESULT([${build_netcdf}])

 dnl Apply tricks if wanted
 if test "${build_netcdf}" = "yes" -a "${enable_tricks}" = "yes"; then
  ABI_TRICKS_NETCDF(${fc_type},${fc_version})
  if test "${netcdf_tricks_bypass}" = "yes"; then
   build_netcdf="no"
  fi
 fi

 dnl Define preprocessing macro
 if test "${enable_netcdf}" = "yes"; then
  AC_DEFINE([HAVE_NETCDF],1,[Define to 1 if you want support for NETCDF])
 fi

 dnl Substitute variables needed for the use of the plug-in
 AC_SUBST(lib_netcdf_includes)
 AC_SUBST(lib_netcdf_libs)
 AC_SUBST(build_netcdf)

 dnl Inform Automake
 AM_CONDITIONAL(DO_ENABLE_NETCDF,test "${enable_netcdf}" = "yes")
 AM_CONDITIONAL(DO_BUILD_NETCDF,test "${build_netcdf}" = "yes")
]) # ABI_PLUGIN_NETCDF



# ABI_PLUGIN_ETSF_IO()
# -------------------
#
# Sets all variables needed to handle the ETSF_IO plug-in.
#
AC_DEFUN([ABI_PLUGIN_ETSF_IO],
[dnl Initial setup
 lib_etsf_io_includes=""
 lib_etsf_io_libs=""
 build_etsf_io="no"
 abi_plug_ready="no"
 abi_plug_tarball="no"
 abi_plug_bins=""
 abi_plug_incs=""
 abi_plug_libs=""

 dnl Define variables needed to build the library
 etsf_io_pkg_name="etsf_io-1.0.2"
 AC_SUBST(etsf_io_pkg_name)
 etsf_io_pkg_string="ETSF I/O library 1.0.2 (upstream release)"
 AC_SUBST(etsf_io_pkg_string)

 if test -z "${CPPFLAGS_ETSF_IO}"; then
  CPPFLAGS_ETSF_IO="${CPPFLAGS}"
 fi
 AC_SUBST(CPPFLAGS_ETSF_IO)
 if test -z "${CFLAGS_ETSF_IO}"; then
  CFLAGS_ETSF_IO="${CFLAGS}"
 fi
 AC_SUBST(CFLAGS_ETSF_IO)
 if test -z "${CXXFLAGS_ETSF_IO}"; then
  CXXFLAGS_ETSF_IO="${CXXFLAGS}"
 fi
 AC_SUBST(CXXFLAGS_ETSF_IO)
 if test -z "${FCFLAGS_ETSF_IO}"; then
  FCFLAGS_ETSF_IO="${FCFLAGS}"
 fi
 AC_SUBST(FCFLAGS_ETSF_IO)

 dnl Add optimizations for Fortran flags
 FCFLAGS_ETSF_IO="${FCFLAGS_ETSF_IO} ${fcflags_opt_etsf_io}"

 dnl Check whether to activate plug-in
 if test "${enable_etsf_io}" = "yes"; then

  dnl Check for already installed plug-ins
  if test "${with_plugins_prefix}" != ""; then
   AC_MSG_CHECKING([whether ETSF_IO is ready for use])
   abi_plug_ready="yes"
   abi_plug_bins="${with_plugins_prefix}/bin"
   abi_plug_incs="-I${with_plugins_prefix}/include"
   abi_plug_libs="-L${with_plugins_prefix}/lib"
   
   if test ! -x "${with_plugins_prefix}/bin/etsf_io"; then
    abi_plug_ready="no"
   fi

   if test -s "${with_plugins_prefix}/lib/libetsf_io.a"; then
    abi_plug_libs="${abi_plug_libs} -letsf_io"
   else
    abi_plug_ready="no"
   fi

   if test -s "${with_plugins_prefix}/lib/libetsf_io_utils.a"; then
    abi_plug_libs="${abi_plug_libs} -letsf_io_utils"
   else
    abi_plug_ready="no"
   fi

   if test -s "${with_plugins_prefix}/lib/libetsf_io_low_level.a"; then
    abi_plug_libs="${abi_plug_libs} -letsf_io_low_level"
   else
    abi_plug_ready="no"
   fi

   AC_MSG_RESULT([${abi_plug_ready}])
  fi

  if test "${abi_plug_ready}" = "yes"; then
   with_etsf_io_bins="${abi_plug_bins}"
   with_etsf_io_includes="${abi_plug_incs}"
   with_etsf_io_libs="${abi_plug_libs}"
  fi

  dnl Check whether command-line plug-in options have been specified
  if test "${with_etsf_io_includes}" = "" -o "${with_etsf_io_libs}" = ""; then

   dnl Check for a tarball repository
   AC_MSG_CHECKING([for a source tarball of ETSF_IO])
   if test -s "${abinit_tardir}/etsf_io-1.0.2.tar.gz"; then
    abi_plug_tarball="yes"
   fi
   AC_MSG_RESULT([${abi_plug_tarball}])

   dnl Get the package
   if test "${abi_plug_ready}" = "no"; then
    if test "${abi_plug_tarball}" = "no"; then
     AC_MSG_NOTICE([downloading ETSF_IO - this may take a while])

     if test ! -s "${abinit_tardir}/etsf_io-1.0.2.tar.gz"; then
      AC_MSG_CHECKING([availability of ETSF_IO from URL 1])
      ${WGET} --timeout=15 --tries=1 -q -O \
       "${abinit_tardir}/etsf_io-1.0.2.tar.gz" \
       'http://www.abinit.org/plugins/etsf_io-1.0.2.tar.gz'
      test -s "${abinit_tardir}/etsf_io-1.0.2.tar.gz" && abi_plug_tarball="yes"
      AC_MSG_RESULT([${abi_plug_tarball}])
     fi

    fi

    dnl Enable plug-in support only if the download was successful
    if test "${abi_plug_tarball}" = "yes"; then
     lib_etsf_io_includes="-I\$(abinit_builddir)/plugins/etsf_io"
     lib_etsf_io_libs="-L\$(abinit_builddir)/plugins/etsf_io -letsf_io -letsf_io_utils -letsf_io_low_level"
     build_etsf_io="yes"
    else
     lib_etsf_io_includes=""
     lib_etsf_io_libs=""
     enable_etsf_io="no"
     build_etsf_io="no"
     AC_MSG_WARN([could not download ETSF_IO plug-in tarball])
     AC_MSG_WARN([support for ETSF_IO plug-in has been disabled])
    fi
   fi
  else
   lib_etsf_io_includes="${with_etsf_io_includes}"
   lib_etsf_io_libs="${with_etsf_io_libs}"
   build_etsf_io="no"
  fi

 else
  enable_etsf_io="no"
  build_etsf_io="no"
 fi

 dnl Set variables required to build binaries
  if test "${enable_etsf_io}" = "yes"; then
  FCLIBS_ETSF_IO='$(FCLIBS)'
 else
  FCLIBS_ETSF_IO=''
 fi
 AC_SUBST(FCLIBS_ETSF_IO)

 dnl Output results
 AC_MSG_CHECKING([whether to enable the ETSF_IO plug-in])
 AC_MSG_RESULT([${enable_etsf_io}])
 AC_MSG_CHECKING([whether to build the ETSF_IO plug-in])
 AC_MSG_RESULT([${build_etsf_io}])

 dnl Apply tricks if wanted
 if test "${build_etsf_io}" = "yes" -a "${enable_tricks}" = "yes"; then
  ABI_TRICKS_ETSF_IO(${fc_type},${fc_version})
  if test "${etsf_io_tricks_bypass}" = "yes"; then
   build_etsf_io="no"
  fi
 fi

 dnl Define preprocessing macro
 if test "${enable_etsf_io}" = "yes"; then
  AC_DEFINE([HAVE_ETSF_IO],1,[Define to 1 if you want support for ETSF_IO])
 fi

 dnl Substitute variables needed for the use of the plug-in
 AC_SUBST(lib_etsf_io_includes)
 AC_SUBST(lib_etsf_io_libs)
 AC_SUBST(build_etsf_io)

 dnl Inform Automake
 AM_CONDITIONAL(DO_ENABLE_ETSF_IO,test "${enable_etsf_io}" = "yes")
 AM_CONDITIONAL(DO_BUILD_ETSF_IO,test "${build_etsf_io}" = "yes")
]) # ABI_PLUGIN_ETSF_IO



# ABI_PLUGIN_ETSF_XC()
# -------------------
#
# Sets all variables needed to handle the ETSF_XC plug-in.
#
AC_DEFUN([ABI_PLUGIN_ETSF_XC],
[dnl Initial setup
 lib_etsf_xc_includes=""
 lib_etsf_xc_libs=""
 build_etsf_xc="no"
 abi_plug_ready="no"
 abi_plug_tarball="no"
 abi_plug_bins=""
 abi_plug_incs=""
 abi_plug_libs=""

 dnl Define variables needed to build the library
 etsf_xc_pkg_name="etsf_xc-0.9"
 AC_SUBST(etsf_xc_pkg_name)
 etsf_xc_pkg_string="ETSF XC library 0.9 (hacked from the LibXC of Octopus)"
 AC_SUBST(etsf_xc_pkg_string)

 if test -z "${CPPFLAGS_ETSF_XC}"; then
  CPPFLAGS_ETSF_XC="${CPPFLAGS}"
 fi
 AC_SUBST(CPPFLAGS_ETSF_XC)
 if test -z "${CFLAGS_ETSF_XC}"; then
  CFLAGS_ETSF_XC="${CFLAGS}"
 fi
 AC_SUBST(CFLAGS_ETSF_XC)
 if test -z "${CXXFLAGS_ETSF_XC}"; then
  CXXFLAGS_ETSF_XC="${CXXFLAGS}"
 fi
 AC_SUBST(CXXFLAGS_ETSF_XC)
 if test -z "${FCFLAGS_ETSF_XC}"; then
  FCFLAGS_ETSF_XC="${FCFLAGS}"
 fi
 AC_SUBST(FCFLAGS_ETSF_XC)

 dnl Add optimizations for Fortran flags
 FCFLAGS_ETSF_XC="${FCFLAGS_ETSF_XC} ${fcflags_opt_etsf_xc}"

 dnl Check whether to activate plug-in
 if test "${enable_etsf_xc}" = "yes"; then

  dnl Check for already installed plug-ins
  if test "${with_plugins_prefix}" != ""; then
   AC_MSG_CHECKING([whether ETSF_XC is ready for use])
   abi_plug_ready="yes"
   abi_plug_bins="${with_plugins_prefix}/bin"
   abi_plug_incs="-I${with_plugins_prefix}/include"
   abi_plug_libs="-L${with_plugins_prefix}/lib"
   
   if test -s "${with_plugins_prefix}/lib/libxc.a"; then
    abi_plug_libs="${abi_plug_libs} -lxc"
   else
    abi_plug_ready="no"
   fi

   if test -s "${with_plugins_prefix}/lib/libstring_f.a"; then
    abi_plug_libs="${abi_plug_libs} -lstring_f"
   else
    abi_plug_ready="no"
   fi

   AC_MSG_RESULT([${abi_plug_ready}])
  fi

  if test "${abi_plug_ready}" = "yes"; then
   with_etsf_xc_bins="${abi_plug_bins}"
   with_etsf_xc_includes="${abi_plug_incs}"
   with_etsf_xc_libs="${abi_plug_libs}"
  fi

  dnl Check whether command-line plug-in options have been specified
  if test "${with_etsf_xc_includes}" = "" -o "${with_etsf_xc_libs}" = ""; then

   dnl Check for a tarball repository
   AC_MSG_CHECKING([for a source tarball of ETSF_XC])
   if test -s "${abinit_tardir}/etsf_xc-0.9.tar.gz"; then
    abi_plug_tarball="yes"
   fi
   AC_MSG_RESULT([${abi_plug_tarball}])

   dnl Get the package
   if test "${abi_plug_ready}" = "no"; then
    if test "${abi_plug_tarball}" = "no"; then
     AC_MSG_NOTICE([downloading ETSF_XC - this may take a while])

     if test ! -s "${abinit_tardir}/etsf_xc-0.9.tar.gz"; then
      AC_MSG_CHECKING([availability of ETSF_XC from URL 1])
      ${WGET} --timeout=15 --tries=1 -q -O \
       "${abinit_tardir}/etsf_xc-0.9.tar.gz" \
       'http://www.abinit.org/plugins/etsf_xc-0.9.tar.gz'
      test -s "${abinit_tardir}/etsf_xc-0.9.tar.gz" && abi_plug_tarball="yes"
      AC_MSG_RESULT([${abi_plug_tarball}])
     fi

    fi

    dnl Enable plug-in support only if the download was successful
    if test "${abi_plug_tarball}" = "yes"; then
     lib_etsf_xc_includes="-I\$(abinit_builddir)/plugins/etsf_xc"
     lib_etsf_xc_libs="-L\$(abinit_builddir)/plugins/etsf_xc -lxc -lstring_f"
     build_etsf_xc="yes"
    else
     lib_etsf_xc_includes=""
     lib_etsf_xc_libs=""
     enable_etsf_xc="no"
     build_etsf_xc="no"
     AC_MSG_WARN([could not download ETSF_XC plug-in tarball])
     AC_MSG_WARN([support for ETSF_XC plug-in has been disabled])
    fi
   fi
  else
   lib_etsf_xc_includes="${with_etsf_xc_includes}"
   lib_etsf_xc_libs="${with_etsf_xc_libs}"
   build_etsf_xc="no"
  fi

 else
  enable_etsf_xc="no"
  build_etsf_xc="no"
 fi

 dnl Set variables required to build binaries
  dnl No binary for ETSF_XC

 dnl Output results
 AC_MSG_CHECKING([whether to enable the ETSF_XC plug-in])
 AC_MSG_RESULT([${enable_etsf_xc}])
 AC_MSG_CHECKING([whether to build the ETSF_XC plug-in])
 AC_MSG_RESULT([${build_etsf_xc}])

 dnl Apply tricks if wanted
 if test "${build_etsf_xc}" = "yes" -a "${enable_tricks}" = "yes"; then
  ABI_TRICKS_ETSF_XC(${fc_type},${fc_version})
  if test "${etsf_xc_tricks_bypass}" = "yes"; then
   build_etsf_xc="no"
  fi
 fi

 dnl Define preprocessing macro
 if test "${enable_etsf_xc}" = "yes"; then
  AC_DEFINE([HAVE_ETSF_XC],1,[Define to 1 if you want support for ETSF_XC])
 fi

 dnl Substitute variables needed for the use of the plug-in
 AC_SUBST(lib_etsf_xc_includes)
 AC_SUBST(lib_etsf_xc_libs)
 AC_SUBST(build_etsf_xc)

 dnl Inform Automake
 AM_CONDITIONAL(DO_ENABLE_ETSF_XC,test "${enable_etsf_xc}" = "yes")
 AM_CONDITIONAL(DO_BUILD_ETSF_XC,test "${build_etsf_xc}" = "yes")
]) # ABI_PLUGIN_ETSF_XC



# ABI_PLUGIN_FOX()
# -------------------
#
# Sets all variables needed to handle the FOX plug-in.
#
AC_DEFUN([ABI_PLUGIN_FOX],
[dnl Initial setup
 lib_fox_includes=""
 lib_fox_libs=""
 build_fox="no"
 abi_plug_ready="no"
 abi_plug_tarball="no"
 abi_plug_bins=""
 abi_plug_incs=""
 abi_plug_libs=""

 dnl Define variables needed to build the library
 fox_pkg_name="FoX-4.0.1"
 AC_SUBST(fox_pkg_name)
 fox_pkg_string="FoX Fortran XML library 4.0.1 (upstream release)"
 AC_SUBST(fox_pkg_string)

 if test -z "${CPPFLAGS_FOX}"; then
  CPPFLAGS_FOX="${CPPFLAGS}"
 fi
 AC_SUBST(CPPFLAGS_FOX)
 if test -z "${CFLAGS_FOX}"; then
  CFLAGS_FOX="${CFLAGS}"
 fi
 AC_SUBST(CFLAGS_FOX)
 if test -z "${CXXFLAGS_FOX}"; then
  CXXFLAGS_FOX="${CXXFLAGS}"
 fi
 AC_SUBST(CXXFLAGS_FOX)
 if test -z "${FCFLAGS_FOX}"; then
  FCFLAGS_FOX="${FCFLAGS}"
 fi
 AC_SUBST(FCFLAGS_FOX)

 dnl Add optimizations for Fortran flags
 FCFLAGS_FOX="${FCFLAGS_FOX} ${fcflags_opt_fox}"

 dnl Check whether to activate plug-in
 if test "${enable_fox}" = "yes"; then

  dnl Check for already installed plug-ins
  if test "${with_plugins_prefix}" != ""; then
   AC_MSG_CHECKING([whether FOX is ready for use])
   abi_plug_ready="yes"
   abi_plug_bins="${with_plugins_prefix}/bin"
   abi_plug_incs="-I${with_plugins_prefix}/include"
   abi_plug_libs="-L${with_plugins_prefix}/lib"
   
   if test ! -x "${with_plugins_prefix}/bin/FoX-config"; then
    abi_plug_ready="no"
   fi

   if test -s "${with_plugins_prefix}/lib/libFoX_wxml.a"; then
    abi_plug_libs="${abi_plug_libs} -lFoX_wxml"
   else
    abi_plug_ready="no"
   fi

   if test -s "${with_plugins_prefix}/lib/libFoX_wcml.a"; then
    abi_plug_libs="${abi_plug_libs} -lFoX_wcml"
   else
    abi_plug_ready="no"
   fi

   if test -s "${with_plugins_prefix}/lib/libFoX_utils.a"; then
    abi_plug_libs="${abi_plug_libs} -lFoX_utils"
   else
    abi_plug_ready="no"
   fi

   if test -s "${with_plugins_prefix}/lib/libFoX_sax.a"; then
    abi_plug_libs="${abi_plug_libs} -lFoX_sax"
   else
    abi_plug_ready="no"
   fi

   if test -s "${with_plugins_prefix}/lib/libFoX_common.a"; then
    abi_plug_libs="${abi_plug_libs} -lFoX_common"
   else
    abi_plug_ready="no"
   fi

   if test -s "${with_plugins_prefix}/lib/libfsys.a"; then
    abi_plug_libs="${abi_plug_libs} -lfsys"
   else
    abi_plug_ready="no"
   fi

   AC_MSG_RESULT([${abi_plug_ready}])
  fi

  if test "${abi_plug_ready}" = "yes"; then
   with_fox_bins="${abi_plug_bins}"
   with_fox_includes="${abi_plug_incs}"
   with_fox_libs="${abi_plug_libs}"
  fi

  dnl Check whether command-line plug-in options have been specified
  if test "${with_fox_includes}" = "" -o "${with_fox_libs}" = ""; then

   dnl Check for a tarball repository
   AC_MSG_CHECKING([for a source tarball of FOX])
   if test -s "${abinit_tardir}/FoX-4.0.1.tar.gz"; then
    abi_plug_tarball="yes"
   fi
   AC_MSG_RESULT([${abi_plug_tarball}])

   dnl Get the package
   if test "${abi_plug_ready}" = "no"; then
    if test "${abi_plug_tarball}" = "no"; then
     AC_MSG_NOTICE([downloading FOX - this may take a while])

     if test ! -s "${abinit_tardir}/FoX-4.0.1.tar.gz"; then
      AC_MSG_CHECKING([availability of FOX from URL 1])
      ${WGET} --timeout=15 --tries=1 -q -O \
       "${abinit_tardir}/FoX-4.0.1.tar.gz" \
       'http://www.uszla.me.uk/software/source/FoX/FoX-4.0.1.tgz'
      test -s "${abinit_tardir}/FoX-4.0.1.tar.gz" && abi_plug_tarball="yes"
      AC_MSG_RESULT([${abi_plug_tarball}])
     fi
     if test ! -s "${abinit_tardir}/FoX-4.0.1.tar.gz"; then
      AC_MSG_CHECKING([availability of FOX from URL 2])
      ${WGET} --timeout=15 --tries=1 -q -O \
       "${abinit_tardir}/FoX-4.0.1.tar.gz" \
       'http://www.abinit.org/plugins/FoX-4.0.1.tar.gz'
      test -s "${abinit_tardir}/FoX-4.0.1.tar.gz" && abi_plug_tarball="yes"
      AC_MSG_RESULT([${abi_plug_tarball}])
     fi

    fi

    dnl Enable plug-in support only if the download was successful
    if test "${abi_plug_tarball}" = "yes"; then
     lib_fox_includes="-I\$(abinit_builddir)/plugins/fox"
     lib_fox_libs="-L\$(abinit_builddir)/plugins/fox -lFoX_wxml -lFoX_wcml -lFoX_utils -lFoX_sax -lFoX_common -lfsys"
     build_fox="yes"
    else
     lib_fox_includes=""
     lib_fox_libs=""
     enable_fox="no"
     build_fox="no"
     AC_MSG_WARN([could not download FOX plug-in tarball])
     AC_MSG_WARN([support for FOX plug-in has been disabled])
    fi
   fi
  else
   lib_fox_includes="${with_fox_includes}"
   lib_fox_libs="${with_fox_libs}"
   build_fox="no"
  fi

 else
  enable_fox="no"
  build_fox="no"
 fi

 dnl Set variables required to build binaries
  if test "${enable_fox}" = "yes"; then
  FCLIBS_FOX='$(FCLIBS)'
 else
  FCLIBS_FOX=''
 fi
 AC_SUBST(FCLIBS_FOX)

 dnl Output results
 AC_MSG_CHECKING([whether to enable the FOX plug-in])
 AC_MSG_RESULT([${enable_fox}])
 AC_MSG_CHECKING([whether to build the FOX plug-in])
 AC_MSG_RESULT([${build_fox}])

 dnl Apply tricks if wanted
 if test "${build_fox}" = "yes" -a "${enable_tricks}" = "yes"; then
  ABI_TRICKS_FOX(${fc_type},${fc_version})
  if test "${fox_tricks_bypass}" = "yes"; then
   build_fox="no"
  fi
 fi

 dnl Define preprocessing macro
 if test "${enable_fox}" = "yes"; then
  AC_DEFINE([HAVE_FOX],1,[Define to 1 if you want support for FOX])
 fi

 dnl Substitute variables needed for the use of the plug-in
 AC_SUBST(lib_fox_includes)
 AC_SUBST(lib_fox_libs)
 AC_SUBST(build_fox)

 dnl Inform Automake
 AM_CONDITIONAL(DO_ENABLE_FOX,test "${enable_fox}" = "yes")
 AM_CONDITIONAL(DO_BUILD_FOX,test "${build_fox}" = "yes")
]) # ABI_PLUGIN_FOX



# ABI_PLUGIN_WANNIER90()
# -------------------
#
# Sets all variables needed to handle the WANNIER90 plug-in.
#
AC_DEFUN([ABI_PLUGIN_WANNIER90],
[dnl Initial setup
 lib_wannier90_includes=""
 lib_wannier90_libs=""
 build_wannier90="no"
 abi_plug_ready="no"
 abi_plug_tarball="no"
 abi_plug_bins=""
 abi_plug_incs=""
 abi_plug_libs=""

 dnl Define variables needed to build the library
 wannier90_pkg_name="wannier90-1.1"
 AC_SUBST(wannier90_pkg_name)
 wannier90_pkg_string="Wannier90 program 1.1 (upstream release)"
 AC_SUBST(wannier90_pkg_string)

 if test -z "${CPPFLAGS_WANNIER90}"; then
  CPPFLAGS_WANNIER90="${CPPFLAGS}"
 fi
 AC_SUBST(CPPFLAGS_WANNIER90)
 if test -z "${CFLAGS_WANNIER90}"; then
  CFLAGS_WANNIER90="${CFLAGS}"
 fi
 AC_SUBST(CFLAGS_WANNIER90)
 if test -z "${CXXFLAGS_WANNIER90}"; then
  CXXFLAGS_WANNIER90="${CXXFLAGS}"
 fi
 AC_SUBST(CXXFLAGS_WANNIER90)
 if test -z "${FCFLAGS_WANNIER90}"; then
  FCFLAGS_WANNIER90="${FCFLAGS}"
 fi
 AC_SUBST(FCFLAGS_WANNIER90)

 dnl Add optimizations for Fortran flags
 FCFLAGS_WANNIER90="${FCFLAGS_WANNIER90} ${fcflags_opt_wannier90}"

 dnl Check whether to activate plug-in
 if test "${enable_wannier90}" = "yes"; then

  dnl Check for already installed plug-ins
  if test "${with_plugins_prefix}" != ""; then
   AC_MSG_CHECKING([whether WANNIER90 is ready for use])
   abi_plug_ready="yes"
   abi_plug_bins="${with_plugins_prefix}/bin"
   abi_plug_incs="-I${with_plugins_prefix}/include"
   abi_plug_libs="-L${with_plugins_prefix}/lib"
   
   if test ! -x "${with_plugins_prefix}/bin/wannier90.x"; then
    abi_plug_ready="no"
   fi

   if test -s "${with_plugins_prefix}/lib/libwannier.a"; then
    abi_plug_libs="${abi_plug_libs} -lwannier"
   else
    abi_plug_ready="no"
   fi

   AC_MSG_RESULT([${abi_plug_ready}])
  fi

  if test "${abi_plug_ready}" = "yes"; then
   with_wannier90_bins="${abi_plug_bins}"
   with_wannier90_includes="${abi_plug_incs}"
   with_wannier90_libs="${abi_plug_libs}"
  fi

  dnl Check whether command-line plug-in options have been specified
  if test "${with_wannier90_libs}" = ""; then

   dnl Check for a tarball repository
   AC_MSG_CHECKING([for a source tarball of WANNIER90])
   if test -s "${abinit_tardir}/wannier90-1.1.tar.gz"; then
    abi_plug_tarball="yes"
   fi
   AC_MSG_RESULT([${abi_plug_tarball}])

   dnl Get the package
   if test "${abi_plug_ready}" = "no"; then
    if test "${abi_plug_tarball}" = "no"; then
     AC_MSG_NOTICE([downloading WANNIER90 - this may take a while])

     if test ! -s "${abinit_tardir}/wannier90-1.1.tar.gz"; then
      AC_MSG_CHECKING([availability of WANNIER90 from URL 1])
      ${WGET} --timeout=15 --tries=1 -q -O \
       "${abinit_tardir}/wannier90-1.1.tar.gz" \
       'http://quasiamore.mit.edu/wannier/code/wannier90-1.1.tar.gz'
      test -s "${abinit_tardir}/wannier90-1.1.tar.gz" && abi_plug_tarball="yes"
      AC_MSG_RESULT([${abi_plug_tarball}])
     fi
     if test ! -s "${abinit_tardir}/wannier90-1.1.tar.gz"; then
      AC_MSG_CHECKING([availability of WANNIER90 from URL 2])
      ${WGET} --timeout=15 --tries=1 -q -O \
       "${abinit_tardir}/wannier90-1.1.tar.gz" \
       'http://www.abinit.org/plugins/wannier90-1.1.tar.gz'
      test -s "${abinit_tardir}/wannier90-1.1.tar.gz" && abi_plug_tarball="yes"
      AC_MSG_RESULT([${abi_plug_tarball}])
     fi

    fi

    dnl Enable plug-in support only if the download was successful
    if test "${abi_plug_tarball}" = "yes"; then
     lib_wannier90_includes=""
     lib_wannier90_libs="-L\$(abinit_builddir)/plugins/wannier90 -lwannier"
     build_wannier90="yes"
    else
     lib_wannier90_includes=""
     lib_wannier90_libs=""
     enable_wannier90="no"
     build_wannier90="no"
     AC_MSG_WARN([could not download WANNIER90 plug-in tarball])
     AC_MSG_WARN([support for WANNIER90 plug-in has been disabled])
    fi
   fi
  else
   lib_wannier90_includes="${with_wannier90_includes}"
   lib_wannier90_libs="${with_wannier90_libs}"
   build_wannier90="no"
  fi

 else
  enable_wannier90="no"
  build_wannier90="no"
 fi

 dnl Set variables required to build binaries
  if test "${enable_wannier90}" = "yes"; then
  FCLIBS_WANNIER90='$(FCLIBS)'
 else
  FCLIBS_WANNIER90=''
 fi
 AC_SUBST(FCLIBS_WANNIER90)

 dnl Output results
 AC_MSG_CHECKING([whether to enable the WANNIER90 plug-in])
 AC_MSG_RESULT([${enable_wannier90}])
 AC_MSG_CHECKING([whether to build the WANNIER90 plug-in])
 AC_MSG_RESULT([${build_wannier90}])

 dnl Apply tricks if wanted
 if test "${build_wannier90}" = "yes" -a "${enable_tricks}" = "yes"; then
  ABI_TRICKS_WANNIER90(${fc_type},${fc_version})
  if test "${wannier90_tricks_bypass}" = "yes"; then
   build_wannier90="no"
  fi
 fi

 dnl Define preprocessing macro
 if test "${enable_wannier90}" = "yes"; then
  AC_DEFINE([HAVE_WANNIER90],1,[Define to 1 if you want support for WANNIER90])
 fi

 dnl Substitute variables needed for the use of the plug-in
 AC_SUBST(lib_wannier90_includes)
 AC_SUBST(lib_wannier90_libs)
 AC_SUBST(build_wannier90)

 dnl Inform Automake
 AM_CONDITIONAL(DO_ENABLE_WANNIER90,test "${enable_wannier90}" = "yes")
 AM_CONDITIONAL(DO_BUILD_WANNIER90,test "${build_wannier90}" = "yes")
]) # ABI_PLUGIN_WANNIER90



# ABI_PLUGIN_XMLF90()
# -------------------
#
# Sets all variables needed to handle the XMLF90 plug-in.
#
AC_DEFUN([ABI_PLUGIN_XMLF90],
[dnl Initial setup
 lib_xmlf90_includes=""
 lib_xmlf90_libs=""
 build_xmlf90="no"
 abi_plug_ready="no"
 abi_plug_tarball="no"
 abi_plug_bins=""
 abi_plug_incs=""
 abi_plug_libs=""

 dnl Define variables needed to build the library
 xmlf90_pkg_name="xmlf90-1.2g"
 AC_SUBST(xmlf90_pkg_name)
 xmlf90_pkg_string="XML Fortran 90 Library 1.2g (upstream release)"
 AC_SUBST(xmlf90_pkg_string)

 if test -z "${CPPFLAGS_XMLF90}"; then
  CPPFLAGS_XMLF90="${CPPFLAGS}"
 fi
 AC_SUBST(CPPFLAGS_XMLF90)
 if test -z "${CFLAGS_XMLF90}"; then
  CFLAGS_XMLF90="${CFLAGS}"
 fi
 AC_SUBST(CFLAGS_XMLF90)
 if test -z "${CXXFLAGS_XMLF90}"; then
  CXXFLAGS_XMLF90="${CXXFLAGS}"
 fi
 AC_SUBST(CXXFLAGS_XMLF90)
 if test -z "${FCFLAGS_XMLF90}"; then
  FCFLAGS_XMLF90="${FCFLAGS}"
 fi
 AC_SUBST(FCFLAGS_XMLF90)

 dnl Add optimizations for Fortran flags
 FCFLAGS_XMLF90="${FCFLAGS_XMLF90} ${fcflags_opt_xmlf90}"

 dnl Check whether to activate plug-in
 if test "${enable_xmlf90}" = "yes"; then

  dnl Check for already installed plug-ins
  if test "${with_plugins_prefix}" != ""; then
   AC_MSG_CHECKING([whether XMLF90 is ready for use])
   abi_plug_ready="yes"
   abi_plug_bins="${with_plugins_prefix}/bin"
   abi_plug_incs="-I${with_plugins_prefix}/include"
   abi_plug_libs="-L${with_plugins_prefix}/lib"
   
   if test -s "${with_plugins_prefix}/lib/libflib.a"; then
    abi_plug_libs="${abi_plug_libs} -lflib"
   else
    abi_plug_ready="no"
   fi

   AC_MSG_RESULT([${abi_plug_ready}])
  fi

  if test "${abi_plug_ready}" = "yes"; then
   with_xmlf90_bins="${abi_plug_bins}"
   with_xmlf90_includes="${abi_plug_incs}"
   with_xmlf90_libs="${abi_plug_libs}"
  fi

  dnl Check whether command-line plug-in options have been specified
  if test "${with_xmlf90_includes}" = "" -o "${with_xmlf90_libs}" = ""; then

   dnl Check for a tarball repository
   AC_MSG_CHECKING([for a source tarball of XMLF90])
   if test -s "${abinit_tardir}/xmlf90-1.2g.tar.gz"; then
    abi_plug_tarball="yes"
   fi
   AC_MSG_RESULT([${abi_plug_tarball}])

   dnl Get the package
   if test "${abi_plug_ready}" = "no"; then
    if test "${abi_plug_tarball}" = "no"; then
     AC_MSG_NOTICE([downloading XMLF90 - this may take a while])

     if test ! -s "${abinit_tardir}/xmlf90-1.2g.tar.gz"; then
      AC_MSG_CHECKING([availability of XMLF90 from URL 1])
      ${WGET} --timeout=15 --tries=1 -q -O \
       "${abinit_tardir}/xmlf90-1.2g.tar.gz" \
       'http://www.abinit.org/plugins/xmlf90-1.2g.tar.gz'
      test -s "${abinit_tardir}/xmlf90-1.2g.tar.gz" && abi_plug_tarball="yes"
      AC_MSG_RESULT([${abi_plug_tarball}])
     fi

    fi

    dnl Enable plug-in support only if the download was successful
    if test "${abi_plug_tarball}" = "yes"; then
     lib_xmlf90_includes="-I\$(abinit_builddir)/plugins/xmlf90"
     lib_xmlf90_libs="-L\$(abinit_builddir)/plugins/xmlf90 -lflib"
     build_xmlf90="yes"
    else
     lib_xmlf90_includes=""
     lib_xmlf90_libs=""
     enable_xmlf90="no"
     build_xmlf90="no"
     AC_MSG_WARN([could not download XMLF90 plug-in tarball])
     AC_MSG_WARN([support for XMLF90 plug-in has been disabled])
    fi
   fi
  else
   lib_xmlf90_includes="${with_xmlf90_includes}"
   lib_xmlf90_libs="${with_xmlf90_libs}"
   build_xmlf90="no"
  fi

 else
  enable_xmlf90="no"
  build_xmlf90="no"
 fi

 dnl Set variables required to build binaries
  dnl No binary for XMLF90

 dnl Output results
 AC_MSG_CHECKING([whether to enable the XMLF90 plug-in])
 AC_MSG_RESULT([${enable_xmlf90}])
 AC_MSG_CHECKING([whether to build the XMLF90 plug-in])
 AC_MSG_RESULT([${build_xmlf90}])

 dnl Apply tricks if wanted
 if test "${build_xmlf90}" = "yes" -a "${enable_tricks}" = "yes"; then
  ABI_TRICKS_XMLF90(${fc_type},${fc_version})
  if test "${xmlf90_tricks_bypass}" = "yes"; then
   build_xmlf90="no"
  fi
 fi

 dnl Define preprocessing macro
 if test "${enable_xmlf90}" = "yes"; then
  AC_DEFINE([HAVE_XMLF90],1,[Define to 1 if you want support for XMLF90])
 fi

 dnl Substitute variables needed for the use of the plug-in
 AC_SUBST(lib_xmlf90_includes)
 AC_SUBST(lib_xmlf90_libs)
 AC_SUBST(build_xmlf90)

 dnl Inform Automake
 AM_CONDITIONAL(DO_ENABLE_XMLF90,test "${enable_xmlf90}" = "yes")
 AM_CONDITIONAL(DO_BUILD_XMLF90,test "${build_xmlf90}" = "yes")
]) # ABI_PLUGIN_XMLF90
