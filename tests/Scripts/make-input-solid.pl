#
# Script for generating input files for solids, corresponding
# to all elements in the periodic table. 

# Copyright (C) 2001-2008 ABINIT group (RC,XG)
# This file is distributed under the terms of the
# GNU General Public License, see ~abinit/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
# For the initials of contributors, see ~abinit/doc/developers/contributors.txt .

# NOTE : under Unix, a Run.solid script will be automatically generated by
# the command  make perl  in the ~abinit directory.
#
# Usage :
# unix C-shell: Solid  (to be filled by Razvan)
#
# For example, under unix C-shell :
#  Solid  (to be filled by Razvan)
#

use Getopt::Long;

#preparatives:

#initialize the main variables (default)
$calcul = 0;
$step = 0;
$iscf = 5;
$dilatmx = 1.2;
$ecut = 5;
$ecutsm = 2.0;
$ionmov = 3;
$ixc = 1;
$kptopt = 1;
$ndtset = 0;
$ngkpt = 2;
$nstep = 50;
$ntime = 50;
$ntype = 1;
$occopt = 0;
$optcell = 2;
$toldfe = 0.00000001;
$toldff = 0.00001;
$tolmxf = 0.0001;
$tolvrs = 0;
$tsmear = 0.01;

# interpret the input data (ARGV)
if ($#ARGV < 1 ) {
  print " Please, provide at least the option pseudo and\n";
  print " the name of the pseudopotential file\n";
  print " then, the number of an element \n";
  print " for example : Solid -pseudo ../Psps_for_tests/01h.bare -test 1 \n";
  exit 16;
}

#check the input to find the script's options 
for ($i = 0; $i < $#ARGV; $i++) {
   if ($ARGV[$i] eq "pseudo") {
   $PSEUDO=$ARGV[$i+1];
    if (-r $PSEUDO) {
     eval("goto flagok");
    }
    else {
     print "Error 16. Pseudopotential file, $PSEUDO, does not exist! \n\n\n";
     exit 16;
    };
   };
   $calcul = $ARGV[$i+1] if $ARGV[$i] eq "calcul";
   $step = $ARGV[$i+1] if $ARGV[$i] eq "step";
flagok:
}
print $calcul,"\n";
print $ixc,"\n";

# open every file
open(PSP,"<$PSEUDO");
open(INP,">ab.in");

# read the header of the pseudopotential file

#analyze all the entry options
for ($i = 0; $i < $#ARGV; $i++) {
if ($ARGV[$i] eq "ecut") {
 $ecut = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "ecutsm") {
 $ecutsm = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "dilatmx") {
 $dilatmx = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "ionmov") {
 $ionmov = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "iscf") {
 $iscf = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "ndtset") {
 $ndtset = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "ngkpt") {
 $ngkpt = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "nstep") {
 $nstep = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "ntime") {
 $ntime = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "occopt") {
 $occopt = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "optcell") {
 $optcell = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "toldfe") {
 $toldfe = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "toldff") {
 $toldff = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "tolvrs") {
 $tolvrs = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "tolmxf") {
 $tolmxf = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "tsmear") {
 $tsmear = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "ixc") {
 $ixc = $ARGV[$i+1];
 eval("goto nexti");
 };
if ($ARGV[$i] eq "znucl") {
 $znucl = $ARGV[$i+1];
 eval("goto nexti");
 };
eval("goto nexti") if $ARGV[$i] eq "kptopt";
eval("goto nexti") if $ARGV[$i] eq "type";
eval("goto nexti") if $ARGV[$i] eq "ntype";
eval("goto nexti") if substr($ARGV[$i],0,3) eq "spg";
eval("goto nexti") if $ARGV[$i] eq "bravlatt";
eval("goto nexti") if $ARGV[$i] eq "angdeg";
eval("goto nexti") if $ARGV[$i] eq "rprim";
eval("goto nexti") if $ARGV[$i] eq "natom";
eval("goto nexti") if substr($ARGV[$i],0,1) eq "x";
eval("goto nexti") if substr($ARGV[$i],0,2) eq "rf";
eval("goto nexti") if $ARGV[$i] eq "acell";
eval("goto nexti") if $ARGV[$i] eq "angdeg";
eval("goto nexti") if $ARGV[$i] eq "rprim";
eval("goto nexti") if $ARGV[$i] eq "shiftk";
eval("goto nexti") if $ARGV[$i] eq "pseudo";
eval("goto nexti") if $ARGV[$i] eq "step";
eval("goto nexti") if $ARGV[$i] eq "calcul";
print INP "$ARGV[$i] $ARGV[$i+1] \n";

nexti:
$i++;
}


# write a header for the input file
print INP "# ABINIT input file generated using the solid.pl script \n";
print INP "# name of the pseudopotential file : $PSEUDO \n";

#we start to build the input file

#GS calculation
if ($calcul eq 0) {
 print INP "ecut $ecut\n";
 print INP "ngkpt $ngkpt $ngkpt $ngkpt 1\n";
 if ($tolvrs eq 0) {
  print INP sprintf "toldfe %15.12f \n", $toldfe;
 }
 else {
  print INP sprintf "tolvrs  %15.12f \n", $tolvrs;
 };
}
#RX calculation
elsif ($calcul eq 1) {
 print INP "ecut $ecut\n";
 print INP "ecutsm $ecutsm\n";
 print INP "dilatmx $dilatmx\n";
 print INP "ngkpt $ngkpt $ngkpt $ngkpt 1\n";
 print INP "ionmov $ionmov\n";
 print INP "optcell $optcell\n";
 print INP sprintf "toldff  %15.12f \n", $toldff;
 print INP sprintf "tolmxf  %15.12f \n", $tolmxf;
 print INP "ntime $ntime\n";
}
#CE calculation
elsif ($calcul eq 2) {
 $ndtset=5 if $ndtset eq 0;
 print INP "ndtset $ndtset\n";
 $step = 5 if ($step eq 0);
 for ($ii=1; $ii<$ndtset+1; $ii++) {
  $ecuti=$ecut+$ii*$step;  
  print INP "ecut$ii $ecuti\n";
 };
 print INP "ngkpt $ngkpt $ngkpt $ngkpt 1\n";
 if ($tolvrs eq 0) {
  print INP sprintf "toldfe %15.12f \n", $toldfe;
 }
 else {
  print INP sprintf "tolvrs  %15.12f \n", $tolvrs;
 };
}
#CK calculation
elsif ($calcul eq 3) {
 $ndtset=2 if $ndtset eq 0;
 print INP "ndtset $ndtset\n";
 $step = 2 if ($step eq 0);
 for ($ii=1; $ii<$ndtset+1; $ii++) {
  $ngkpti = $ngkpt + ($ii-1)*$step;
  print INP "ngkpt$ii $ngkpti $ngkpti $ngkpti 1\n";
 };
 print INP "ecut $ecut\n";
 if ($tolvrs eq 0) {
  print INP sprintf "toldfe %15.12f \n", $toldfe; 
 }
 else {
  print INP sprintf "tolvrs  %15.12f \n", $tolvrs;
 };
}
#ERROR - undefined calculation
else
{
print "Undefined type of calculation, $calcul\n\n\n";
exit 5;
};

#pseudopotential-dependent variables  
print INP "znucl $znucl \n";
print INP "ixc $ixc\n";

# common variables:
print INP "iscf $iscf\n";
print INP "kptopt $kptopt\n";
print INP "nstep $nstep\n";
print INP "ntype $ntype\n";
print INP sprintf "tsmear %5.3f \n", $tsmear;

#chemical element choice 
eval("goto Case_$znucl");
# previous command should jump to a specific element; 
# if it doesn't, the element is undefined
print "Error: test $znucl is not defined.";
exit 8;

#WARNING! 
#the values for acell are given in Angstroms for input in specvar and 
#in bohr when they are written directly in the INP file

SWITCH_ELEMENTS: {

Case_1:          #H
$occopt = 1 if $occopt eq 0;
&specvar(1,'hcp',3.650,5.830);
eval("goto Case_LAST");

Case_2:          #He
$occopt = 1 if $occopt eq 0;
&specvar(1,'hcp',3.57,5.83);
eval("goto Case_LAST");

Case_3:          #Li
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',4.404);
eval("goto Case_LAST");

Case_4:          #Be
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',2.2866,3.5833);
eval("goto Case_LAST");

Case_5:          #B
print INP "acell 16.5351 16.5351 9.5620\n";
print INP "angdeg 90 90 90\n";
print INP "natom 50\n";
print INP "natrd 5\n";
$occopt=1 if $occopt eq 0;
print INP "spgroup 134\n";
print INP "type 50*1\n";
print INP "xred \n";
print INP "0.3253      0.0883      0.3985\n";
print INP "0.2272      0.0805      0.0865\n";
print INP "0.1195      0.1195      0.378\n";
print INP "0.2425      0.2425      0.5815\n";
print INP "0.          0.          0.5\n";
eval("goto Case_LAST");

Case_6:          #C
$occopt = 1 if $occopt eq 0;
&specvar(1,'diam',3.56679);
eval("goto Case_LAST");

Case_7:          #N
print INP "acell 10.66561 10.66561 10.66561\n";
print INP "angdeg 90 90 90\n";
print INP "natom 8\n";
print INP "natrd 1\n";
$occopt=1 if $occopt eq 0;
print INP "spgroup 205\n";
print INP "type 1\n";
print INP "xred \n";
print INP "0.054 0.054 0.054\n";
eval("goto Case_LAST");

Case_8:          #O
print INP "acell 12.90683 12.90683 12.90683\n";
print INP "angdeg 90 90 90\n";
print INP "natom 8\n";
print INP "natrd 2\n";
$occopt=1 if $occopt eq 0;
print INP "spgroup 223\n";
print INP "type 50*1\n";
print INP "xred \n";
print INP "0.253 0.5 0.\n";
print INP "0. 0. 0.\n";
eval("goto Case_LAST");

Case_9:          #F
$occopt = 1 if $occopt eq 0;
print "not yet available!";
eval("goto Case_LAST");

Case_10:          #Ne
$occopt = 1 if $occopt eq 0;
&specvar(1,'fcc',4.429);
eval("goto Case_LAST");

Case_11:          #Na
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',4.225);
eval("goto Case_LAST");

Case_12:          #Mg
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',3.20927,5.21033);
eval("goto Case_LAST");

Case_13:          #Al
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',4.04958);
eval("goto Case_LAST");

Case_14:          #Si
$occopt = 1 if $occopt eq 0;
&specvar(1,'diam',5.43070);
eval("goto Case_LAST");

Case_15:          #P
print INP "acell 6.255 8.277 19.842\n";
print INP "angdeg 90 90 90\n";
print INP "chkprim 0";
print INP "natom 8\n";
print INP "natrd 1\n";
$occopt=1 if $occopt eq 0;
print INP "spgroup 64\n";
print INP "spgaxor 4\n";
print INP "type 8*1\n";
print INP "xred \n";
print INP "0.0 0.090 0.098\n";
#print INP "0.0 0.910 0.902\n";
#print INP "0.0 0.590 0.402\n";
#print INP "0.0 0.410 0.502\n";
#print INP "0.5 0.590 0.902\n";
#print INP "0.5 0.410 0.098\n";
#print INP "0.5 0.090 0.598\n";
#print INP "0.5 0.910 0.402\n";
eval("goto Case_LAST");

Case_16:          #S, rhombohedral
print INP "acell 12.20763 12.20763 12.20763\n";
print INP "angdeg 115.18 115.18 115.18\n";
print INP "natom 18\n";
print INP "natrd 1\n";
print INP "spgroup 148\n";
print INP "spgorig 2\n";
print INP "type 6*1\n";
$occopt=1 if $occopt eq 0;
print INP "xred 0.1454 0.1882 0.1055\n";
#print INP "     0.8118 0.9572 0.1055\n";
#print INP "     0.0428 0.8546 0.1055\n";
#print INP "     0.8546 0.8118 0.8945\n";
#print INP "     0.1882 0.0428 0.8945\n";
#print INP "     0.9572 0.8546 0.8945\n";
eval("goto Case_LAST");

Case_17:          #Cl
print INP "acell 11.79189 15.60913 8.46597\n";
print INP "angdeg 90 90 90\n";
print INP "chkprim 0\n";
print INP "natom 8\n";
print INP "natrd 1\n";
$occopt=1 if $occopt eq 0;
print INP "spgroup 64\n";
print INP "spgaxor 4\n";
print INP "type 8*1\n";
print INP "xred \n";
print INP "0.0 0.100 0.130\n";
eval("goto Case_LAST");

Case_18:          #Ar
$occopt = 1 if $occopt eq 0;
&specvar(1,'fcc',5.256);
eval("goto Case_LAST");

Case_19:          #K
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',5.225);
eval("goto Case_LAST");

Case_20:          #Ca
$occopt = 4 if $occopt eq 0;
&specvar(4,"fcc ",5.576);
eval("goto Case_LAST");

Case_21:          #Sc
$occopt = 4 if $occopt eq 0;
&specvar(4,"fcc ",4.541);
eval("goto Case_LAST");

Case_22:          #Ti
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',2.95,4.686);
eval("goto Case_LAST");

Case_23:          #V
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',3.24);
eval("goto Case_LAST");

Case_24:          #Cr
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',2.8839);
eval("goto Case_LAST");

Case_25:          #Mn
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',3.081);
eval("goto Case_LAST");
#structure stable above 1134 C up to the melting point 1245 C

Case_26:          #Fe
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',3.591);
eval("goto Case_LAST");

Case_27:          #Co
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',3.548);
eval("goto Case_LAST");

Case_28:          #Ni
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',2.65,4.33);
eval("goto Case_LAST");

Case_29:          #Cu
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',3.61496);
eval("goto Case_LAST");

Case_30:          #Zn
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',2.6648,4.9467);
eval("goto Case_LAST");

Case_31:          #Ga       
print INP "acell 8.523987 8.535325 14.446577\n";
print INP "angdeg 90 90 90\n";
print INP "chkprim 0\n";
print INP "natom 8\n";
print INP "natrd 1\n";
$occopt=1 if $occopt eq 0;
print INP "spgroup 64\n";
print INP "spgaxor 4\n";
print INP "type 8*1\n";
print INP "xred \n";
print INP "0.0 0.078500000000 0.152500000000\n";
#print INP "0.0 0.921500000000 0.847500000000\n";
#print INP "0.0 0.578500000000 0.347500000000\n";
#print INP "0.0 0.421500000000 0.652500000000\n";
#print INP "0.5 0.578500000000 0.847500000000\n";
#print INP "0.5 0.421500000000 0.152500000000\n";
#print INP "0.5 0.078500000000 0.652500000000\n";
#print INP "0.5 0.921500000000 0.347500000000\n";
eval("goto Case_LAST");

Case_32:          #Ge
$occopt = 1 if $occopt eq 0;
&specvar(1,'diam',5.65695);
eval("goto Case_LAST");

Case_33:          #As        
print INP "acell 7.80645 7.80645 7.80645\n";
print INP "angdeg 54.1 54.1 54.1\n";
print INP "natom 2\n";
print INP "type 2*1\n";
$occopt=1 if $occopt eq 0;
print INP "xred 0.226 0.226 0.226\n";
print INP "     0.774 0.774 0.774\n";
eval("goto Case_LAST");


Case_34:          #Se        
print INP "acell 8.230078 8.230078 9.353104\n";
print INP "angdeg 90 90 120\n";
print INP "natom 3\n";
print INP "natrd 1\n";
print INP "type 3*1\n";
$occopt=1 if $occopt eq 0;
print INP "spgroup 154\n";
print INP "xred 0.217 0.0   0.0\n";
#print INP "     0.783 0.783 0.333333333333\n";
#print INP "     0.0   0.217 0.666666666667\n";
eval("goto Case_LAST");

Case_35:          #Br
print INP "acell 12.6045 16.4784 8.4659\n";
print INP "angdeg 90 90 90\n";
print INP "chkprim 0\n";
print INP "natom 8\n";
print INP "natrd 1\n";
$occopt=1 if $occopt eq 0;
print INP "spgroup 64\n";
print INP "spgaxor 4\n";
print INP "type 8*1\n";
print INP "xred \n";
print INP "0.0 0.1100 0.1350\n";
eval("goto Case_LAST");

Case_36:          #Kr
$occopt = 1 if $occopt eq 0;
&specvar(1,'fcc',5.721);    #5.706
eval("goto Case_LAST");

Case_37:          #Rb
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',5.585);
eval("goto Case_LAST");

Case_38:          #Sr
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',6.0847);
eval("goto Case_LAST");

Case_39:          #Y
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',3.6474,5.7306);
eval("goto Case_LAST");

Case_40:          #Zr
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',3.232,5.147);
eval("goto Case_LAST");

Case_41:          #Nb
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',3.3004);
eval("goto Case_LAST");

Case_42:          #Mo
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',4.16);
eval("goto Case_LAST");

Case_43:          #Tc
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',2.735,4.388);
eval("goto Case_LAST");

Case_44:          #Ru
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',2.70389,4.28168);
eval("goto Case_LAST");

Case_45:          #Rh
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',3.8031);
eval("goto Case_LAST");

Case_46:          #Pd
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',3,8898);
eval("goto Case_LAST");

Case_47:          #Ag
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',4.0862);
eval("goto Case_LAST");

Case_48:          #Cd
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',2.97887,5.61765);
eval("goto Case_LAST");

Case_49:          #In        
print INP "acell 6.13027 6.13027 9.33147\n";
print INP "angdeg 90 90 90\n";
print INP "natom 2\n";
print INP "type 2*1\n";
$occopt=1 if $occopt eq 4;
print INP "xred 0.0 0.0 0.0\n";
print INP "     0.5 0.5 0.5\n";
print INP "chkprim 0\n";
eval("goto Case_LAST");

Case_50:          #Sn
$occopt = 1 if $occopt eq 0;
&specvar(1,'diam',6.4912);
eval("goto Case_LAST");

Case_51:          #Sb        
print INP "acell 8.51625803 8.51625803 8.51625803\n";
print INP "angdeg 57.6 57.6 57.6\n";
print INP "natom 2\n";
print INP "type 2*1\n";
$occopt=1 if $occopt eq 4;
print INP "xred 0.233 0.233 0.233\n";
print INP "     0.767 0.767 0.767\n";
eval("goto Case_LAST");

Case_52:          #Te        
print INP "acell 8.4034792 8.4034792 11.1775780\n";
print INP "angdeg 90 90 120\n";
print INP "natom 3\n";
print INP "natrd 1\n";
print INP "type 3*1\n";
$occopt=1 if $occopt eq 0;
print INP "spgroup 154\n";
print INP "xred 0.269 0.0 0.0\n";
#print INP "     0.731 0.731 0.333333333333\n";
#print INP "     0.0 0.731 0.666666666667\n";
eval("goto Case_LAST");

Case_53:          #I
print INP "acell 13.73844 18.50692 9.05186\n";
print INP "angdeg 90 90 90\n";
print INP "chkprim 0\n";
print INP "natom 8\n";
print INP "natrd 1\n";
$occopt=1 if $occopt eq 0;
print INP "spgroup 64\n";
print INP "spgaxor 4\n";
print INP "type 8*1\n";
print INP "xred \n";
print INP "0.0 0.1156 0.1493\n";
eval("goto Case_LAST");

Case_54:          #Xe
$occopt = 1 if $occopt eq 0;
&specvar(1,'fcc',6.197); #6.2023
eval("goto Case_LAST");

Case_55:          #Cs
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',6.045);
eval("goto Case_LAST");

Case_56:          #Ba
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',5.0);
eval("goto Case_LAST");

Case_57:          #La
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',5.296);
eval("goto Case_LAST");

Case_58:          #Ce
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',5.1612);
eval("goto Case_LAST");

Case_59:          #Pr
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',5.161);
eval("goto Case_LAST");

Case_60:          #Nd
print INP "acell 6.91243 6.91243 22.29725\n";
print INP "angdeg 90 90 120\n";
print INP "natom 4\n";
print INP "type 4*1\n";
$occopt=1 if $occopt eq 4;
print INP "xred\n";
print INP "0 0 0 \n";
print INP "0 0 0.5 \n";
print INP "0.333333333333 0.666666666667 0.75 \n";
print INP "0.666666666667 0.333333333333 0.25 \n";
eval("goto Case_LAST");

Case_61:          #Pm
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',3.651,5.8601);
eval("goto Case_LAST");

Case_62:          #Sm   
print INP "acell 16.999975 16.999975 16.999975\n";
print INP "angdeg 23.13 23.13 23.13\n";
print INP "natom 3\n";
print INP "natrd 2\n";
print INP "type 3*1\n";
$occopt=1 if $occopt eq 4;
print INP "spgroup 166\n";
print INP "spgorig 2\n";
print INP "xred 0.0 0.0 0.0\n";
print INP "    0.222 0.222 0.222\n";
#print INP "     0.778 0.778 0.778\n";
eval("goto Case_LAST");

Case_63:          #Eu
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',4.551);
eval("goto Case_LAST");

Case_64:          #Gd
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',3.629,5.796);
eval("goto Case_LAST");

Case_65:          #Tb
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',3,601,5.6936);
eval("goto Case_LAST");

Case_66:          #Dy
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',3.584,5.668);
eval("goto Case_LAST");

Case_67:          #Ho
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',3.5773,5.6158);
eval("goto Case_LAST");

Case_68:          #Er
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',3.55,5.59);
eval("goto Case_LAST");

Case_69:          #Tm
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',3.5375,5.5546);
eval("goto Case_LAST");

Case_70:          #Yb
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',5.4862);
eval("goto Case_LAST");

Case_71:          #Lu
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',3.5031,5.5509);
eval("goto Case_LAST");

Case_72:          #Hf
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',3.1967,5.0578);
eval("goto Case_LAST");

Case_73:          #Ta
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',3.3058);
eval("goto Case_LAST");

Case_74:          #W
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',3.16469);
eval("goto Case_LAST");

Case_75:          #Re
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',2.7608,4.4582);
eval("goto Case_LAST");

Case_76:          #Os
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',2.7352,4.319);
eval("goto Case_LAST");

Case_77:          #Ir
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',3.8394);
eval("goto Case_LAST");

Case_78:          #Pt
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',3.9231);
eval("goto Case_LAST");

Case_79:          #Au
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',4.07825);
eval("goto Case_LAST");

Case_80:          #Hg     
print INP "acell 5.643288 5.643288 5.643288\n";
print INP "angdeg 70.46 70.46 70.46\n";
print INP "natom 1\n";
print INP "type 1\n";
$occopt=1 if $occopt eq 4;
print INP "xred 0.0 0.0 0.0\n";
eval("goto Case_LAST");

Case_81:          #Tl
$occopt = 4 if $occopt eq 0;
&specvar(4,'hcp',3.438,5.478);
eval("goto Case_LAST");

Case_82:          #Pb
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',4.9505);
eval("goto Case_LAST");

Case_83:          #Bi   
print INP "acell 8.96845 8.96845 8.96845\n";
print INP "angdeg 57.14 57.14 57.14\n";
print INP "natom 2\n";
print INP "type 1 1\n";
$occopt=1 if $occopt eq 0;
print INP "xred 0.23806 0.23806 0.23806\n";
print INP "     0.76194 0.76194 0.76194\n";
eval("goto Case_LAST");

Case_84:          #Po
print INP "acell 6.34759 6.34759 6.34759\n";
print INP "angdeg 90 90 90\n";
print INP "natom 1\n";
print INP "type 1\n";
$occopt=1 if $occopt eq 0;
print INP "xred 0.0 0.0 0.0\n";
eval("goto Case_LAST");

Case_85:          #At
print "not yet available !";
eval("goto Case_LAST");

Case_86:          #Rn
print "not yet available !";
eval("goto Case_LAST");

Case_87:          #Fr
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',6.50);
eval("goto Case_LAST");

Case_88:          #Ra
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',5.148);
eval("goto Case_LAST");

Case_89:          #Ac
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',5.311);
eval("goto Case_LAST");

Case_90:          #Th
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',5.0843);
eval("goto Case_LAST");

Case_91:          #Pa  
print INP "acell 7.41717 7.41717 6.11893\n";
print INP "angdeg 90 90 90\n";
print INP "natom 2\n";
print INP "type 2*1\n";
$occopt=4 if $occopt eq 0;
print INP "xred 0.0 0.0 0.0\n";
print INP "     0.5 0.5 0.5\n";
print INP "chkprim 0\n";
eval("goto Case_LAST");

Case_92:          #U    
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',3.474);
eval("goto Case_LAST");

Case_93:          #Np
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',3.52);
eval("goto Case_LAST");

Case_94:          #Pu
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',3.638);
eval("goto Case_LAST");

Case_95:          #Am
$occopt = 4 if $occopt eq 0;
&specvar(4,'fcc',4.894);
eval("goto Case_LAST");

Case_96:          #Cm
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',3.747);
eval("goto Case_LAST");

Case_97:          #Bk
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',3.802);
eval("goto Case_LAST");

Case_98:          #Cf
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',3.857);
eval("goto Case_LAST");

Case_99:          #Es
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',3.913);
eval("goto Case_LAST");

Case_100:          #Fm
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',3.968);
eval("goto Case_LAST");

Case_101:          #Md
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',4.023);
eval("goto Case_LAST");

Case_102:          #No
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',4.078);
eval("goto Case_LAST");

Case_103:          #Lr
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',4.134);
eval("goto Case_LAST");

Case_104:          #Rf
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',4.189);
eval("goto Case_LAST");

Case_105:          #Ha
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',4.244);
eval("goto Case_LAST");

Case_106:          #Sg
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',4.300);
eval("goto Case_LAST");

Case_107:          #Ns
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',4.355);
eval("goto Case_LAST");

Case_108:          #Hs
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',4.410);
eval("goto Case_LAST");

Case_109:          #Mt
$occopt = 4 if $occopt eq 0;
&specvar(4,'bcc',4.465);
eval("goto Case_LAST");

Case_LAST:
} #end SWITCH_ELEMENTS




#build the "in" file
$term3=".GS" if $calcul eq 0;
$term3=".RX" if $calcul eq 1;
$term3=".CE" if $calcul eq 2;
$term3=".CK" if $calcul eq 3;
$term5=".LDA" if $ixc < 10;
$term5=".GGA" if $ixc > 10;
$inpfile="in.Z".$znucl.$term3.".EC".$ecut.$term5;

$status = system("mv ab.in $inpfile");

#build the *.files file
$ifile=substr($inpfile,3);
$filesfile=$inpfile."files";
open(FILES,">$filesfile");
print FILES "$inpfile\n";
print FILES "$ifile.output\n";
print FILES "$ifile.inp\n";
print FILES "$ifile.out\n";
print FILES "$ifile.tmp\n";
print FILES "$PSEUDO\n";
close(FILES);

$status = system("abinis < $filesfile > $ifile.log");

#****************************************
sub specvar {
      local($occopt,$strtyp,$cell1,$cell2) = @_;


print INP "occopt $occopt \n";

#select the structure type
# and then fill in the input file with the rprim, acell and xred values

#FCC
if ($strtyp eq 'fcc') {
 $a01=$cell1/0.529177249;
 print INP "acell $a01 $a01 $a01  \n";
 print INP "natom 1\n";
 print INP "type 1\n";
 print INP "rprim 0.0 0.5 0.5  0.5 0.0 0.5  0.5 0.5 0.0 \n";
 print INP "xred 0.0 0.0 0.0 \n";
}

#BCC
if ($strtyp eq 'bcc') {
 $a01=$cell1/0.529177249;
 print INP "acell $a01 $a01 $a01  \n";
 print INP "natom 1\n"; 
 print INP "type 1\n";
 print INP "rprim -0.5 0.5 0.5  0.5 -0.5 0.5  0.5 0.5 -0.5 \n";
 print INP "xred 0.0 0.0 0.0 \n";
}

#DIAMOND
if ($strtyp eq 'diam') {
 $a01=$cell1/0.529177249;
 $x1=$cell1/4;
 print INP "acell $a01 $a01 $a01  \n";
 print INP "natom 2\n";
 print INP "type 1 1\n";
 print INP "rprim 0.0 0.5 0.5  0.5 0.0 0.5  0.5 0.5 0.0 \n";
 print INP "xcart 0.0 0.0 0.0 \n";
 print INP "      $x1 $x1 $x1 \n";
}

#HCP
if ($strtyp eq 'hcp') {
 $a01=$cell1/0.529177249;
 $a02=$cell2/0.529177249;
 print INP "acell  $a01 $a01 $a02\n";
 print INP "natom 2\n";
 print INP "type 1 1\n";
 print INP "angdeg 90 90 120 \n"; 
 print INP "xred 0.333333333333 0.666666666667 0.25 \n";
 print INP "     0.666666666667 0.333333333333 0.75 \n";
}

}

