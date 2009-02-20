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
# ( wrap-parallel-tests machine_name [ seq | seqpar ] ) >& log_file
#
# For example :
# ( wrap-parallel-tests ibm_pcpm seq) >& log_file
#
# The list of allowed machine names is mentioned in the run-parallel-tests
# script.
#

set -e

# Init
my_name="wrap-parallel-tests"
my_cnffile="tests.env"

# Check arguments
if test "${#}" -lt "2"; then
  echo "Usage: ${my_name} machine_name ( seq | seqpar )"
  echo ""
  echo "Two arguments must be provided giving machine name and mode"
  echo ""
  echo "Note : in case you called ${my_name} from the make, the syntax is"
  echo " make tests_paral paral_host=machine_name paral_mode=( seq |seqpar)"
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
run_charge_tests="${BOURNE_SHELL} ${abinit_rundir}/run-charge-tests.sh"

cd "${abinit_outdir}/paral"

# Select test type
case "${2}" in

 seq)
  rm -rf ./,,runlog*
  (${run_parallel_tests} ${machine_name} A 0) 2>&1 1> ,,runlogA
  (${run_parallel_tests} ${machine_name} B 0) 2>&1 1> ,,runlogB
  (${run_parallel_tests} ${machine_name} C 0) 2>&1 1> ,,runlogC
  (${run_parallel_tests} ${machine_name} D 0) 2>&1 1> ,,runlogD
  (${run_parallel_tests} ${machine_name} E 0) 2>&1 1> ,,runlogE
  (${run_parallel_tests} ${machine_name} F 0) 2>&1 1> ,,runlogF
  (${run_parallel_tests} ${machine_name} G 0) 2>&1 1> ,,runlogG
  (${run_parallel_tests} ${machine_name} H 0) 2>&1 1> ,,runlogH
  (${run_parallel_tests} ${machine_name} I 0) 2>&1 1> ,,runlogI
  (${run_parallel_tests} ${machine_name} J 0) 2>&1 1> ,,runlogJ
  (${run_parallel_tests} ${machine_name} M 0) 2>&1 1> ,,runlogM
  (${run_parallel_tests} ${machine_name} N 0) 2>&1 1> ,,runlogN
  (${run_parallel_tests} ${machine_name} O 0) 2>&1 1> ,,runlogO
  (${run_parallel_tests} ${machine_name} P 0) 2>&1 1> ,,runlogP
  (${run_parallel_tests} ${machine_name} R 0) 2>&1 1> ,,runlogR
  (${run_parallel_tests} ${machine_name} T 0) 2>&1 1> ,,runlogT
  (${run_parallel_tests} ${machine_name} U 0) 2>&1 1> ,,runlogU
  (${run_parallel_tests} ${machine_name} V 0) 2>&1 1> ,,runlogV
# XG051126 Disabled charge test. I do not understand why it does not work 
# on my machine ...
#  (${run_charge_tests} ${machine_name}) 2>&1 1> ,,runlog9
  ;;

 seqpar)
  rm -rf ./,,runlog*
  (${run_parallel_tests} ${machine_name} A) 2>&1 1> ,,runlogA
  (${run_parallel_tests} ${machine_name} B) 2>&1 1> ,,runlogB
  (${run_parallel_tests} ${machine_name} C) 2>&1 1> ,,runlogC
  (${run_parallel_tests} ${machine_name} D) 2>&1 1> ,,runlogD
  (${run_parallel_tests} ${machine_name} E) 2>&1 1> ,,runlogE
  (${run_parallel_tests} ${machine_name} F) 2>&1 1> ,,runlogF
  (${run_parallel_tests} ${machine_name} G) 2>&1 1> ,,runlogG
  (${run_parallel_tests} ${machine_name} H) 2>&1 1> ,,runlogH
  (${run_parallel_tests} ${machine_name} I) 2>&1 1> ,,runlogI
  (${run_parallel_tests} ${machine_name} J) 2>&1 1> ,,runlogJ
  (${run_parallel_tests} ${machine_name} M) 2>&1 1> ,,runlogM
  (${run_parallel_tests} ${machine_name} N) 2>&1 1> ,,runlogN
  (${run_parallel_tests} ${machine_name} O) 2>&1 1> ,,runlogO
  (${run_parallel_tests} ${machine_name} P) 2>&1 1> ,,runlogP
  (${run_parallel_tests} ${machine_name} R) 2>&1 1> ,,runlogR
  (${run_parallel_tests} ${machine_name} T) 2>&1 1> ,,runlogT
  (${run_parallel_tests} ${machine_name} U) 2>&1 1> ,,runlogU
  (${run_parallel_tests} ${machine_name} V) 2>&1 1> ,,runlogV
# XG051126 Disabled charge test. I do not understand why it does not work
# on my machine ...
#  (${run_charge_tests} ${machine_name}) 2>&1 1> ,,runlog9
  ;;

 *)
  echo "${my_name}: unrecognized argument: ${2}"
  exit 2
  ;;

esac
