!
! Copyright (C) 2005-2008 ABINIT Group (Yann Pouillon)
! All rights reserved.
!
! This file is part of the ABINIT software package. For license information,
! please see the COPYING file in the top-level directory of the ABINIT source
! distribution.
!
! OpenVMS version of the run-basic-tests.sh
!
$ dir_old = f$env("default")
$ my_name="run-basic-tests"
$ my_cnffile="tests.env"
!
! Check arguments
$ if ( p2 .eqs. "" )
$ then
$   write sys$output "Usage: @''my_name' test_dir test_number"
$   exit
$ endif
! Finish init
$ test_dir= p1
$ test_number= p2
$ set default 'p1'
$ my_output="test''test_number'"
!
! Clean-up
$ delete 'my_output'.in;*/nolog
$ delete 'my_output'.out;*/nolog
$ delete 'my_output'.log;*/nolog
$ delete 'my_output'i_*.*;*/nolog
$ delete 'my_output'o_*.*;*/nolog
$ delete 'my_output'_*.*;*/nolog
$ copy [.Input]test'test_number'.in 'my_output'.in
!
! Write abinit.files for test
$ create 'my_output'.files
$ open/append filf 'my_output'.files
$ write filf "''my_output'.in"
$ write filf "''my_output'.out"
$ write filf "''my_output'i"
$ write filf "''my_output'o"
$ write filf "''my_output'"
$ if test_number .eq. 1
$ then
$   write filf "[-.Psps_for_tests]01h.pspgth"
$ else
$   if test_number .eq. 2
$   then
$     write filf "[-.Psps_for_tests]14si.pspgth"
$   else
$     if test_number .eq. 3
$     then
$       write filf "[-.Psps_for_tests]01h.pspgth"
$       write filf "[-.Psps_for_tests]04be.pspgth"
$     else
$       if test_number .eq. 4
$       then
$         write filf "[-.Psps_for_tests]70yb.pspnc"
$       else
$         if test_number .eq. 5
$         then
$           write filf "[-.Psps_for_tests]13al.pspgth"
$         else
$           if test_number .eq. 6
$           then
$             write filf "[-.Psps_for_tests]14si.xml"
$           else
$	      write sys$output "illegal test number"
$	      close filf
$	      set default 'dir_old'
$	      exit
$           endif
$         endif
$       endif
$     endif
$   endif
$ endif
$ close filf
$ pipe run [--]abinis.exe < 'my_output'.files > 'my_output'.log
$ type 'my_output'_STATUS.dat
$ set default 'dir_old'
$ exit
