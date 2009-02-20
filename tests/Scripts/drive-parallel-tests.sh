# sh script for running various tests of abinis and abinip
# Actually run different instances of the Run script.
# See the latter script for more information.

# Copyright (C) 2002-2008 ABINIT group (XG,LSi,YPouillon)
# This file is distributed under the terms of the
# GNU General Public License, see ~abinit/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
# For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
#
# Usage under sh-shell:
# ( drive-parallel-tests machine_name sets [cases] ) >& log_file
#
# For example :
# ( drive-parallel-tests sleepy A ) >& log_file
# ( drive-parallel-tests sleepy A-D ) >& log_file
# ( drive-parallel-tests sleepy A 0 ) >& log_file
# ( drive-parallel-tests sleepy A 1-2 ) >& log_file
#
# The list of allowed machine names is mentioned in the run-parallel-tests
# script. More information on the arguments is also provided in the same script.
#

set -e

# Init
my_name="drive-parallel-tests"
my_cnffile="tests.env"

# Check arguments
if test "${#}" -lt "2"; then
  echo "Usage: ${my_name} machine_name sets [cases]"
  echo ""
  echo "Two arguments must be provided giving machine name and mode"
  exit 0
fi

# Set-up environment
if test -s "${my_cnffile}"; then
 . ${my_cnffile}
else
 echo "${my_name}: ${my_cnffile} not found - aborting now"
 exit 1
fi

# Finish init
machine_name="${1}"
run_parallel_tests="${PERL} ${abinit_rundir}/run-parallel-tests.pl"

cd "${abinit_outdir}/paral"
(${run_parallel_tests} ${machine_name} ${2} ${3}) 2>&1 1> ,,log

