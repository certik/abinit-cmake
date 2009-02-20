#
# Copyright (C) 2005-2008 ABINIT Group (Yann Pouillon)
# All rights reserved.
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

# Stop at first error
set -e

# Init
. config.sh

# Check arguments
if test "${#}" = "0"; then
 echo "Usage: make-split-dist number_of_tarballs"
 echo ""
 exit 0
fi

start_time=`date '+%s'`

cat <<EOF
Make-split-dist report
=========================

EOF

ntar="${1}"

# Create dist directory
echo -n "Creating dist directory..."
make distdir > /dev/null
echo "done."



end_time=`date '+%s'`

cat <<EOF
-- 
Time elapsed : `awk "END{print ${end_time}-${start_time}}" < /dev/null` s

EOF
