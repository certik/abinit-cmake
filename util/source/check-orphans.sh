# This script looks for source files not present in object lists

[ -s configure.ac ] || exit 1

echo "Looking for orphan files"
echo "------------------------"
echo ""

for d in prereqs/[a-z]* plugins/[a-z]* src/defs src/[0-9]*
do
	orphans=""

	echo -n "Directory $d:"

	cd ${d} > /dev/null
	for f in `ls *.f *.F90 2> /dev/null`
	do
		file_lists="object_list"
		[ -s paral_list ] && file_lists="${file_lists} paral_list"
		grep "${f%\.[fF]*}" ${file_lists} > /dev/null || orphans="${orphans} ${f}"
	done
	cd - > /dev/null 2>&1

	echo "${orphans:- none}"
done
