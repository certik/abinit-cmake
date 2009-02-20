#
# Wrapper for the standard tests of abinit
#

my_name="wrap-standard-tests"
my_cnffile="tests.env"

# Check arguments
if test "${#}" -lt "1"; then
  echo "Usage: ${my_name} machine_name [ start_test [ stop_test ] ]"
  echo ""
  exit 0
fi

# Set-up environment
if test -s "${my_cnffile}"; then
 . ${my_cnffile}
else
 echo "${my_name}: ${my_cnffile} not found - aborting now"
 exit 1
fi

# Save log if machine_name is chkinabi
if test "${1}" = "chkinabi"; then
 my_logfile='>& ,,chkinabi.log'
else
 my_logfile=''
fi

${PERL} Scripts/run-standard-tests.pl $@ ${my_logfile}
exit ${?}
