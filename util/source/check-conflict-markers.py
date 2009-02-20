#!/usr/bin/env python

import re
import os
import sys

re_markers  = re.compile("^(<<<<<<< TREE|=======|>>>>>>> MERGE-SOURCE)$")

retval = 0

for root,dirs,files in os.walk("."):
  # Prune Bazaar subdirs
  if ".bzr" in dirs:
    dirs.remove(".bzr")

  # Display conflict markers found
  for item in files:
    chk_data = file("%s/%s" % (root,item),"r").readlines()
    chk_stat = False
    for line in chk_data:
      if ( re_markers.match(line) ):
        chk_stat = True
        retval = 1
    if ( chk_stat ):
      sys.stderr.write("Conflict markers in %s/%s\n" % (root,item))

sys.exit(retval)
