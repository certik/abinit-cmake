#!/bin/sh
#
# Copyright (C) 2005-2008 ABINIT Group (Yann Pouillon)
# All rights reserved.
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

set -e

# Init
my_name="run-netcdf-tests"
my_cnffile="tests.env"

# Check config file
if test -s "${my_cnffile}"; then
 . "${my_cnffile}"
else
 echo "${my_name}: config file ${my_cnffile} not found - aborting now."
 exit 1
fi

# Finish init
my_outdir="${abinit_outdir}/netcdf/,,`hostname`_`date '+%Y%m%d'`"

mkdir -p ${my_outdir} && cd ${my_outdir} && ${abinit_bindir}/abinetcdf
exit ${?}
