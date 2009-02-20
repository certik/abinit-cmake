#!/bin/sh
#
#    ABINIT 5 automatic build script
#    Copyright (C) 2006-2007 Yann Pouillon
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

#
# NOTE: YOU WILL NEED TO EDIT THIS SCRIPT BEFORE IT WORKS.
#

# Initialize environment
abinit_package="abinit-5.2dev"
abinit_hostname=`hostname | sed -e 's/\..*//'`
abinit_srcdir="/scratch/pouillon"
abinit_config="${abinit_srcdir}/${abinit_hostname}-custom.ac"

now=`date '+%Y/%m/%d %H:%M'`
cat <<EOF
Build of ${abinit_package} on `hostname`
Started ${now}

Step         Ret  Stdout  Stderr Duration
------------ --- ------- ------- --------
EOF

# Prepare temporary directory and config file
mkdir -p "${abinit_srcdir}" || exit 1
cd "${abinit_srcdir}"
if test ! -e "${abinit_config}"; then
 ln -s "${HOME}/.abinit/build/${abinit_hostname}.ac" "${abinit_config}"
fi

# Download ABINIT
start=`date '+%s'`
(cd ${abinit_srcdir} && rm -rf ${abinit_package} && \
  wget -q -O - http://www.abinit.org/pouillon/${abinit_package}.tar.gz | \
   gunzip -c - | tar xf - && \
  cd ${abinit_package} && mkdir tmp && cd tmp)
cd ${abinit_srcdir}/${abinit_package}/tmp
ret=${?}
stop=`date '+%s'`

# Touch all files to prevent clock skew issues
find "${abinit_srcdir}/${abinit_package}" -exec touch {} \;

# Print Download stats
osize=`ls -l download.stdout 2> /dev/null | awk '{print $5}'`
esize=`ls -l download.stderr 2> /dev/null | awk '{print $5}'`
osize=${osize:--1}
esize=${esize:--1}

printf '%-12s %3d %7d %7d %8d\n' DOWNLOAD ${ret} ${osize} ${esize} \
 $((${stop}-${start}))

# ---------------------------------------------------------------------------- #

# Configure ABINIT
start=`date '+%s'`
../configure --with-config-file="${abinit_config}" \
 > configure.stdout 2> configure.stderr
ret=${?}
stop=`date '+%s'`

# Print config stats
osize=`ls -l configure.stdout 2> /dev/null | awk '{print $5}'`
esize=`ls -l configure.stderr 2> /dev/null | awk '{print $5}'`
osize=${osize:--1}
esize=${esize:--1}

printf '%-12s %3d %7d %7d %8d\n' CONFIG ${ret} ${osize} ${esize} \
 $((${stop}-${start}))

# ---------------------------------------------------------------------------- #

# Build ABINIT
start=`date '+%s'`
make > make.stdout 2> make.stderr
ret=${?}
stop=`date '+%s'`

# print make stats
osize=`ls -l make.stdout 2> /dev/null | awk '{print $5}'`
esize=`ls -l make.stderr 2> /dev/null | awk '{print $5}'`
osize=${osize:--1}
esize=${esize:--1}

printf '%-12s %3d %7d %7d %8d\n' BUILD ${ret} ${osize} ${esize} \
 $((${stop}-${start}))

# ---------------------------------------------------------------------------- #

# Install ABINIT
start=`date '+%s'`
rm -rf ${abinit_srcdir}/abinit-install
make install DESTDIR=${abinit_srcdir}/abinit-install \
 > install.stdout 2> install.stderr
ret=${?}
stop=`date '+%s'`

# print install stats
osize=`ls -l install.stdout 2> /dev/null | awk '{print $5}'`
esize=`ls -l install.stderr 2> /dev/null | awk '{print $5}'`
osize=${osize:--1}
esize=${esize:--1}

printf '%-12s %3d %7d %7d %8d\n' INSTALL ${ret} ${osize} ${esize} \
 $((${stop}-${start}))

# ---------------------------------------------------------------------------- #

cd tests || exit 1

for series in 'in' fast cpu v1 v2 v3 v4 netcdf; do
 start=`date '+%s'`

 echo "" >> ../tests.stdout
 echo "${series} tests" >> ../tests.stdout
 echo "" >> ../tests.stderr
 echo "${series} tests" >> ../tests.stderr
# if test "${series}" = "v1" -o "${series}" = "v2" -o \
#         "${series}" = "v3" -o "${series}" = "v4"; then
#  make tests_${series} start=01 stop=15 >> ../tests.stdout 2>> ../tests.stderr
#  ret=${?}
# else
 make tests_${series} >> ../tests.stdout 2>> ../tests.stderr
 ret=${?}
# fi

 stop=`date '+%s'`

 # print tests stats
 osize=`ls -l ../tests.stdout 2> /dev/null | awk '{print $5}'`
 esize=`ls -l ../tests.stderr 2> /dev/null | awk '{print $5}'`
 osize=${osize:--1}
 esize=${esize:--1}

 printf '%-12s %3d %7d %7d %8d\n' ${series} ${ret} ${osize} ${esize} \
  $((${stop}-${start}))
done

now=`date '+%Y/%m/%d %H:%M'`
cat <<EOF

Finished ${now}
EOF
