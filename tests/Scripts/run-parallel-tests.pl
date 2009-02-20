# Script for running various tests of abinip, the parallel version of abinit
# Before running this script, read README in this directory.

# Copyright (C) 1998-2008 ABINIT group (XG,LSi)
# This file is distributed under the terms of the
# GNU General Public License, see ~abinit/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
# For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
#
# Because of some execution problems, the script executes only one
# 'test' at a time (See the keyword 'HERE'), meaning that only one
# input file is considered, for different cases, from pure sequential
# to parallel with 10 processors.

# Also, only selected machine names are allowed, see the keyword 'MACH_LIST'
#
# Usage :
#  unix C-shell: Run machine_name sets [cases] >& log_file
#	 unix bash: Run machine_name sets [cases] > log_file 2>&1
#  DOS/Windows: perl Run.pl machine_name sets [cases] > log_file
#		sets or cases can be specified as a) a single set or case;
#			b) a range Start-Last or Start-end of sets or cases.
#
# For example :
# (Run ibm_pcpm A) >& log_file        run, for input file A,
#                                        all cases on the ibm_pcpm cluster
# (Run ibm_pcpm A 0) >& log_file      run, for input file A,
#                                        only case 0 on the ibm_pcpm cluster
# (Run t3e B 1-end) >& log_file     run, for input file B,
#                                        from case_1 to the end on a Cray-T3E
# (Run t3e C-D  1-2a) >& log_file     run, for input files C to D,
#                                        cases between 1 and 2a included,
#                                        on a Cray-T3E.
#*****************************************************************************

# Load ABINIT environment
$abinit_bindir = $ENV{'abinit_bindir'};
$abinit_rundir = $ENV{'abinit_rundir'};
$abinit_pspdir = $ENV{'abinit_pspdir'};
$abinit_inpdir = $ENV{'abinit_inpdir'};
$abinit_outdir = $ENV{'abinit_outdir'};
$PERL = $ENV{'PERL'};

$, = ' ';               # set output field separator
$\ = "\n";              # set output record separator
#
@MachList=('t3e','sgi','tux','vpp','hp','compaq_es40','ibm_44p','sr8k','mpd','sun','dummy','test','spinoza','fock','max','chpit','lemaitre','sleepy(-compilo)','pwrg4','mpich2','generic','chum(-compilo)', 'cea');
@CompiloList = ('g95','pgi','gcc41','gcc42','gcc43','pathscale','intel','intel10','open64','sun');
@Cases = ('0','1','2','4','10');
$debug = 2;
# Check whether a machine name is provided
if ($#ARGV < 1) {
  print 'At least, two arguments must be provided giving machine name and set number';
  print "Choose among :@MachList" if ($#ARGV < 0);
  exit 16;
	}
#
$UNIXSLASH = '/';	# unix subdirectory delimitor in file paths (octal 057)
$DOSSLASH = '\\';	# subdirectory delimitor in DOS file paths (\ escaped)
# try unix 'uname' command in the Bourne shell manner
$rc = system('uname > uname.tmp 2>&1');
unlink ('uname.tmp');
$OStype = $ENV{'OSTYPE'};	# OSTYPE environment on unix systems (HP/UX excepted)
if ($rc == 0) {
	$SLASH = $UNIXSLASH;		# subdirectory delimitor in file paths
	$COPY_COMMAND = 'cp -p';	# unix copy file command
	$DIFF_COMMAND = 'diff -b';	# unix difference command
	$ERASE_COMMAND = 'rm -f';	# unix delete file command
	$XSUFX = '';		# no special suffix for executable module
	if ($OStype eq 'cygwin32') {
# since perl for Windows NT is not a PGI Workstation command (cygwin32) but a
# DOS module, calling the fldiff perl script by means of the "system" function
# requires DOS conventions
		$PERLSLASH = $DOSSLASH;   # subdirectory delimitor as for DOS file paths
		$PLSUFX = '.pl';	# DOS suffix for perl script
		}
	else {
# for perl under normal unix systems:
		$PERLSLASH = $UNIXSLASH;	# subdirectory delimitor for perl file paths
		$PLSUFX = '.pl';		# no special suffix for perl script
		}
	}
# When 'uname' command is unknown, the return code is non zero. Under Windows it may
# be 1 or 256 depending of the Perl version and the command interpretor; check OS
# environmental variable
elsif (($rc == 1 || $rc == 256) && $ENV{'OS'} eq 'Windows_NT') {			
	$OStype = $ENV{'OS'};
	$SLASH = $DOSSLASH;		# subdirectory delimitor in DOS file paths
	$COPY_COMMAND = 'copy';		# DOS copy file command
	$DIFF_COMMAND = 'fc /w /l /n';	# DOS difference command
# since unlink <file*> fails on some versions of perl for DOS, let's use del
	$ERASE_COMMAND = 'del /q';	# DOS delete file command
	$XSUFX = '.exe';		# DOS suffix for executable module
	$PERLSLASH = $DOSSLASH;   # subdirectory delimitor in DOS file path
	$PLSUFX = '.pl';		# DOS suffix for perl script
	$PERL = 'perl';			# perl.exe MUST be accessible through the DOS PATH
	}
else {
	print "unrecognized Operating System $OStype";
	exit 99;
	}

$CODE_SEQ = &perlpath("$abinit_bindir/abinis$XSUFX");	# codename for sequential version
$CODE_PAR = &perlpath("$abinit_bindir/abinip$XSUFX");	# codename for parallel version
$PSPS = &transpath($abinit_pspdir);		# pseudopotential directory  
# Set a default for the second pseudopotential (the latter is needed
# only for the tests number E and H)
$PSP2 = &transpath("$PSPS/33as.SGS_mod");
$FLDIFF = &perlpath("$PERL $abinit_rundir/fldiff.pl");	# fldiff script in utilities directory
#*****************************************************************************
 
# Check whether the machine name is allowed, and initialize its characteristics
# It is the place to add new machines and characteristics
# (but do not forget to mention the new name in the MACH_LIST variable above)

# The following characteristics are adapted for machines where neither lam
# nor mpich have to be used, and where there is always one thread per PE.
$MACH = $ARGV[0];
$LAM = 0;
$MPICH = 0;
# List the set of cases for each machine (mask the second characteristics
# of tests, that may be 0, 1, 2, 4, 10 ; 0, 1 and 2 are always allowed)
$CASE0 = 1;
$CASE1 = 1;
$CASE2 = 1;
$CASE4 = 1;
$CASE10 = 1;
# However, the list of cases depend on the test
#-----------------------------------
# For the t3e in berkeley and Garching
if ($MACH eq 't3e') {
  $MPIRUN_NP = 'mpprun -n';
	}
#-----------------------------------
# For a SGI R1000
elsif ($MACH eq 'sgi') {
  $MPIRUN_NP = 'mpprun -np';
	}
#-----------------------------------
# For the cluster of DEC/Linux "tux" machine at the UCL.
# Note that on this cluster, the local 
# node (tux0) is used for the first thread. 
# This explains why the number of lines
# of the cluster files is smaller than usual.
# 
elsif ($MACH eq 'tux') {
  $MPICH = 1;
  $MPICH_DIR = '/usr/local/mpi/bin';
  $MPIRUN_NP = "$MPICH_DIR/mpirun -np";
  $MACH0='tux0';
  $MACH1='tux2';
  $MACH2='tux3';
  $MACH3='tux5';
  $MACH4='tux6';
  $MACH5='tux7';
  $MACH6='tux8';
  $MACH7='tux9';
  $MACH8='tux10';
  $MACH9='tux11';
  unlink 'cluster1','cluster2','cluster4','cluster10';
	open (FILE,">cluster1") || die "Unable to open FILE cluster1";
	print FILE $MACH0;
	close (FILE);
	open (FILE,">cluster2") || die "Unable to open FILE cluster2";
	print FILE $MACH1;
	close (FILE);
	open (FILE,">cluster4") || die "Unable to open FILE cluster4";
        print FILE $MACH1;
        print FILE $MACH2;
        print FILE $MACH3;
	close (FILE);
	open (FILE,">cluster10") || die "Unable to open FILE cluster10";
        print FILE $MACH1;
        print FILE $MACH2;
        print FILE $MACH3;
        print FILE $MACH4;
        print FILE $MACH5;
        print FILE $MACH6;
        print FILE $MACH7;
        print FILE $MACH8;
        print FILE $MACH9;
	close (FILE);
	}
#-----------------------------------
# For the cluster of DEC/Linux "boop" machine at the UCL.
# Note that on this cluster, the local
# node (tux0) is used for the first thread.
# This explains why the number of lines
# of the cluster files is smaller than usual. 
# 
elsif ($MACH eq 'boop') {
  $MPICH = 1;
  $MPICH_DIR = '/usr/local/mpi/bin';
  $MPIRUN_NP = "$MPICH_DIR/mpirun -np";
  $MACH1='tux0'; 
  $MACH2='boop1';
  $MACH3='boop2';
  $MACH4='boop3';
  unlink 'cluster1','cluster2','cluster4','cluster10';
        open (FILE,">cluster1") || die "Unable to open FILE cluster1";
        print FILE $MACH1;
        close (FILE);
        open (FILE,">cluster2") || die "Unable to open FILE cluster2";
        print FILE $MACH2;
        close (FILE);
        open (FILE,">cluster4") || die "Unable to open FILE cluster4";
        print FILE $MACH2;
        print FILE $MACH2;
        print FILE $MACH2;
        close (FILE);
        open (FILE,">cluster10") || die "Unable to open FILE cluster10";
        print FILE $MACH2;
        print FILE $MACH2;
        print FILE $MACH2;
        print FILE $MACH2;
        print FILE $MACH2;
        print FILE $MACH2;
        print FILE $MACH2;
        print FILE $MACH2;
        print FILE $MACH2;
        close (FILE);
        }
#-----------------------------------
# For the cluster of Intel/Linux machine at the UCL.
# The syntax for the execution of the parallel run is different
# from other machines (see the modifs related to "cox" throughout
# Warning : cox10 had a problem
elsif ($MACH eq 'cox') {
  $MPICH = 1;
  $MPICH_DIR = '/usr/local/mpi/bin';
  $MPIRUN_NP = "$MPICH_DIR/mpirun -np";
  $MACH1 = 'cox3.pcpm.ucl.ac.be';
  $MACH2 = 'cox11.pcpm.ucl.ac.be';
  $MACH3 = 'cox12.pcpm.ucl.ac.be';
  $MACH4 = 'cox14.pcpm.ucl.ac.be';
  $MACH5 = 'cox15.pcpm.ucl.ac.be';
  unlink 'cluster1','cluster2','cluster4','cluster10';
	open (FILE,">cluster1") || die "Unable to open FILE cluster1";
	print FILE '1';
	print FILE "$MACH1 2";
	close (FILE);
	open (FILE,">cluster2") || die "Unable to open FILE cluster2";
	print FILE '2';
	print FILE "$MACH1 2";
	print FILE "$MACH2 2";
	close (FILE);
	open (FILE,">cluster4") || die "Unable to open FILE cluster4";
	print FILE '4';
	print FILE "$MACH1 2";
	print FILE "$MACH2 2";
	print FILE "$MACH3 2";
	print FILE "$MACH4 2";
	close (FILE);
	open (FILE,">cluster10") || die "Unable to open FILE cluster10";
	print FILE '10';
	print FILE "$MACH1 2";
	print FILE "$MACH1 4";
	print FILE "$MACH2 2";
	print FILE "$MACH2 4";
	print FILE "$MACH3 2";
	print FILE "$MACH3 4";
	print FILE "$MACH4 2";
	print FILE "$MACH4 4";
	print FILE "$MACH5 2";
	print FILE "$MACH5 4";
	close (FILE);
	}
#-----------------------------------
# For the IBM 44P machine at the PCPM lab, using MPICH
elsif ($MACH eq 'ibm_44p') {
  $MPICH = 1;
  $MPICH_DIR = '/usr/local/mpi/bin';
  $MPIRUN_NP = "$MPICH_DIR/mpirun -np";
  $MACH1 = 'dirac';
  $CASE10 = 0;
  unlink 'cluster1','cluster2','cluster4';
	open (FILE,">cluster1") || die "Unable to open FILE cluster1";
	print FILE $MACH1;
	close (FILE);
	open (FILE,">cluster2") || die "Unable to open FILE cluster2";
	print FILE $MACH1;
	print FILE $MACH1;
	close (FILE);
	open (FILE,">cluster4") || die "Unable to open FILE cluster4";
		for ($i = 0; $i < 4; $i++) {
		print FILE $MACH1;
		}
	close (FILE);
	}
#-----------------------------------
# For the IBM Power 5  machine at the PCPM lab, using MPICH
elsif ($MACH eq 'fock') {
  $MPICH = 1;
  $MPICH_DIR = '/usr/local/mpich/bin';
  $MPIRUN_NP = "$MPICH_DIR/mpirun -np";
  $MACH1 = 'fock';
  $CASE10 = 0;
  unlink 'cluster1','cluster2','cluster4';
        open (FILE,">cluster1") || die "Unable to open FILE cluster1";
        print FILE $MACH1;
        close (FILE);
        open (FILE,">cluster2") || die "Unable to open FILE cluster2";
        print FILE $MACH1;
        print FILE $MACH1;
        close (FILE);
        open (FILE,">cluster4") || die "Unable to open FILE cluster4";
                for ($i = 0; $i < 4; $i++) {
                print FILE $MACH1;
                }
        close (FILE);
        }
#-----------------------------------
# For the Fujitsu VPP700 at Riken
elsif ($MACH eq 'vpp') {
  $MPIRUN_NP = ' -np';
  $CASE10 = 0;
	}
#-----------------------------------
# For the Compaq/DEC 4 processors in LLN
elsif ($MACH eq 'compaq_es40') {
  $MPIRUN_NP = ' dmpirun -np';
  $CASE10 = 0;
	}
#-----------------------------------
# For the HP integrity Itanium 4 processor in LLN
elsif ($MACH eq 'chpit') {
  $MPIRUN_NP = ' /usr/local/mpich/bin/mpirun -np';
  $CASE10 = 0;
        }
#-----------------------------------
# For the Sun V20 Opteron-based machine in LLN
elsif ($MACH eq 'lemaitre') {
  $MPIRUN_NP = '/usr/local/mpich-1.2.6-eth-intel9.1/bin/mpirun -np';
        }
#-----------------------------------
# For the HP exemplar, S-Class and N4000
elsif ($MACH eq 'hp') {
  $MPIRUN_NP = 'mpirun -j -np';
  $CASE10 = 0;
	}
#-----------------------------------
# For the Hitachi SR8000
elsif ($MACH eq 'sr8k' || $MACH eq 'mpd') {
  $MPIRUN_NP = 'mpiexec -n';
	}
# 
elsif ($MACH eq 'cea') {
  $MPIRUN_NP = 'cea_mprun -p parallele_128 -n ';
	}
#-----------------------------------
# For the Sun sunfire 2 proc
elsif ($MACH eq 'sun') {
  $MPIRUN_NP = '/opt/SUNWhpc/HPC5.0/bin/mprun -np';
  $CASE4 = 0;
  $CASE10 = 0;
	}
#-----------------------------------
# For a bi-processor SGI R14000
elsif ($MACH eq 'spinoza') {
  $MPIRUN_NP = 'mpirun -np';
  $CASE4 = 0;
  $CASE10 = 0;
        }
#-----------------------------------
# For the Intel PC XEON "dummy" 2 procs (with hyperthreading)
elsif ($MACH eq 'dummy') {
  $MPICH = 1;
  $MPICH_DIR = '/usr/local/mpi-pgi/bin';
  $MPIRUN_NP = "$MPICH_DIR/mpirun -np";
  $MACH0='dummy.pcpm.ucl.ac.be';
  $MACH1='dummy.pcpm.ucl.ac.be:4';
  $MACH2='dummy.pcpm.ucl.ac.be:2';
  $CASE10 = 0;
  unlink 'cluster1','cluster2','cluster4';
        open (FILE,">cluster1") || die "Unable to open FILE cluster1";
        print FILE $MACH1;
        close (FILE);
        open (FILE,">cluster2") || die "Unable to open FILE cluster2";
        print FILE $MACH2;
        close (FILE);
        open (FILE,">cluster4") || die "Unable to open FILE cluster4";
        print FILE $MACH1;
        close (FILE);
        }
#-----------------------------------
# For the Intel PC XEON "sleepy" 2 procs (with hyperthreading)
elsif ($MACH =~ /^sleepy-?(.*)/) {
  $MPICH = 1;
  $MPICH_DIR = "/usr/local/mpi-$1/bin";
  if ($1 eq '') { $MPICH_DIR = "/usr/local/mpi-ifc9.1-64/bin"; }
# verify if $MPICH_DIR exists
  if (! -e $MPICH_DIR ) {
     print "The $MPICH_DIR folder is not present...";
     print "Choose among : @CompiloList";
     exit 15;
  }
# Note that one uses the mpirun from PGI, but the one from IFC is OK also ...
  $MPIRUN_NP = "$MPICH_DIR/mpirun -np";
  $MACH0='sleepy.pcpm.ucl.ac.be';
  $MACH1='sleepy.pcpm.ucl.ac.be:4';
  $MACH2='sleepy.pcpm.ucl.ac.be:2';
  $CASE10 = 0;
  unlink 'cluster1','cluster2','cluster4';
        open (FILE,">cluster1") || die "Unable to open FILE cluster1";
        print FILE $MACH1;
        close (FILE);
        open (FILE,">cluster2") || die "Unable to open FILE cluster2";
        print FILE $MACH2;
        close (FILE);
        open (FILE,">cluster4") || die "Unable to open FILE cluster4";
        print FILE $MACH1;
        close (FILE);
}
#-----------------------------------
# For the AMD SUN GALAXY X4200  "chum" 2 x AMD 2220 Dual Cores with OpenMPI or MPICH2
elsif ($MACH =~ /^chum-?(.*)/) {
  $MPICH = 2;
  $OPENMPI = 1;
#  $MPICH_DIR = "/usr/local/mpich2-$1/bin";
  $MPICH_DIR = "/usr/local/openmpi-$1/bin";
  if ($1 eq '') { $MPICH_DIR = "/usr/local/openmpi-intel/bin"; }
  if ($1 eq 'sun') { $MPICH_DIR = "/usr/local/mpich2-sun/bin"; $OPENMPI = 0;}
#  if ($1 eq 'open64') { $MPICH_DIR = "/usr/local/mpich2-open64/bin"; $OPENMPI = 0;}
# verify if $MPICH_DIR exists
  if (! -e $MPICH_DIR ) {
     print "The $MPICH_DIR folder is not present...";
     print "Choose among : @CompiloList";
     exit 15;
  }
print "The $MPICH_DIR will be used...($OPENMPI)";
#
  $MPIRUN_NP = "$MPICH_DIR/mpirun -np";
  $MACH0='chum.pcpm.ucl.ac.be';
  $MACH1='chum.pcpm.ucl.ac.be:4';
  $MACH2='chum.pcpm.ucl.ac.be:2';
  $CASE10 = 0;
  unlink 'cluster1','cluster2','cluster4';
        open (FILE,">cluster1") || die "Unable to open FILE cluster1";
        print FILE $MACH1;
        close (FILE);
        open (FILE,">cluster2") || die "Unable to open FILE cluster2";
        print FILE $MACH2;
        close (FILE);
        open (FILE,">cluster4") || die "Unable to open FILE cluster4";
        print FILE $MACH1;
        close (FILE);
}
#-----------------------------------
# For the MacBookPRo with Intel 2 cores "pwrg4"
elsif ($MACH eq 'pwrg4') {
  $MPICH = 1;
  $MPICH_DIR = "/usr/local/mpich-1.2.7p1-gcc-4.3.0/bin"; 
  $MPIRUN_NP = "$MPICH_DIR/mpirun -np";
  $MACH0='pwrg4.pcpm.ucl.ac.be';
  $MACH1='pwrg4.pcpm.ucl.ac.be:4';
  $MACH2='pwrg4.pcpm.ucl.ac.be:2';
  $CASE10 = 0;
  unlink 'cluster1','cluster2','cluster4';
        open (FILE,">cluster1") || die "Unable to open FILE cluster1";
        print FILE $MACH1;
        close (FILE);
        open (FILE,">cluster2") || die "Unable to open FILE cluster2";
        print FILE $MACH2;
        close (FILE);
        open (FILE,">cluster4") || die "Unable to open FILE cluster4";
        print FILE $MACH1;
        close (FILE);
        }
#-----------------------------------
# For the cluster of PowerPC Mac OS X "max" with biprocs
# For the cluster of PowerPC Mac OS X "max" with biprocs
elsif ($MACH eq 'max') {
  $MPICH = 1;
  $MPICH_DIR = '/usr/local/mpich-1.2.7/bin';
  $MPIRUN_NP = "$MPICH_DIR/mpirun -np";
  $MACH0='max';
  $MACH1='max';
  $MACH2='maxi8';
  $MACH3='maxi8';
  $CASE10 = 0;
  unlink 'cluster1','cluster2','cluster4';
        open (FILE,">cluster1") || die "Unable to open FILE cluster1";
        print FILE $MACH0;
        close (FILE);
        open (FILE,">cluster2") || die "Unable to open FILE cluster2";
        print FILE $MACH0;
        print FILE $MACH1;
        close (FILE);
        open (FILE,">cluster4") || die "Unable to open FILE cluster4";
        print FILE $MACH0;
        print FILE $MACH1;
        print FILE $MACH2;
        print FILE $MACH3;
        close (FILE);
        }
#-----------------------------------
# For mpich2-1.0.5p4 (gmatteo@tux)
elsif ($MACH eq 'mpich2') {
 $MPIRUN_NP = '/opt/mpich2-1.0.5p4/bin/mpirun -np ';
        }
#-----------------------------------
# generic
# use the mpirun that is in the PATH
elsif ($MACH eq 'generic') {
 $MPIRUN_NP = 'mpirun -np ';
        }
#-----------------------------------
# This is a dummy machine name, for sequential tests only
# (One should check that this is only used in the sequential case)
elsif ($MACH eq 'test') {
 $MPIRUN_NP = ' ';
        }
#-----------------------------------
else {
  print "The machine name $MACH is not allowed";
  print "Choose among : @MachList";
  exit 12;
	}

print "Testing the Abinit code on the $MACH $OStype platform";
# Here stop the examination of machine names
# and the definition of their characteristics
#*******************************************************************************
$_ = $ARGV[1];
$hitnum = m/[A-W]/;	# check for upper case
$hitrng = m/([A-W])-(.*)/;	# check for range
if ($hitrng) {
	$FirstSet = $1;
	$LastSet = $2;
	$LastSet = '' if ($LastSet eq 'end');
	}
elsif ($hitnum) {	# if last set is omitted, run only first one
	$FirstSet = $ARGV[1];
	$LastSet = $ARGV[1];
	}
else {
	print "Invalid sets: $ARGV[1]";
	exit 24;
	}

# set label for default first and last test number
$FirstCase = '0';
$LastCase = '10';
$_ = $ARGV[2];
$hit = m/(.*)-(.*)/;	# check for range
	if ($hit) {
		$FirstCase = $1;
		$LastCase = $2;
		$LastCase = '10' if ($LastCase eq 'end');
		}
	elsif ($ARGV[2] ne '') {	# if last test is omitted, run only first one
		$FirstCase = $ARGV[2];
		$LastCase = $ARGV[2];
		}

print "Following tests will be run: sets $FirstSet to $LastSet, cases $FirstCase to $LastCase";

#################################################
# set a date flag to make a new directory today #
#################################################
($sec,$min,$hour,$mday,$ymon,$yyear,$wday,$yday,$isdst)=localtime(time);
$ymon++;	# ymon was 0-11
$yyear +=1900;	# yyear was relative to 1900
$YYYYMMDD = sprintf("%4.4d",$yyear).sprintf("%2.2d",$ymon).sprintf("%2.2d",$mday);

$WORK_DIR = ',,'.$MACH.'_'.$YYYYMMDD;
if (! -e $WORK_DIR || ! -d $WORK_DIR) {
	mkdir ($WORK_DIR,0755);		# Mode 0755 ignored under DOS-Windows
	}
else {
	print "Can not create directory, $WORK_DIR already exists";
	}

print "cd $WORK_DIR";
chdir ("$WORK_DIR");

# **************************************** 
# Loop on all the sets
for ($T1 = $FirstSet; $T1 <= $LastSet; $T1 ++) {
# Make preliminary operations that might differ for different tests

# HERE, the set of possible tests and corresponding input file
	if ($T1 eq 'A' || $T1 eq 'B') {
	  $PSP1 = "$PSPS/14si.psp";
	  }
#
	elsif ($T1 eq 'C') {
	  $PSP1 = "$PSPS/13al.pspgth";
	  $PSP2 = "$PSPS/33as.SGS_mod";
	  }
#
	elsif ($T1 eq 'D') {
	  $PSP1 = "$PSPS/42mo.pspnc";
          $CASE10=0;
	  }
#
	elsif ($T1 eq 'E' || $T1 eq 'H' || $T1 eq 'J') {
	  $PSP1 = "$PSPS/31ga.SGS_mod";
          $PSP2 = "$PSPS/33as.SGS_mod";
	  $OPTFLDIFF = '-medium';
	  }
#
	elsif ($T1 eq 'F') {
	  $PSP1 = "$PSPS/13al.pspgth";
          $PSP2 = "$PSPS/33as.SGS_mod";
		unlink 'tF.i_DS1_DEN','tF.i_DS2_DEN';
	  system("$COPY_COMMAND tC0.o_DEN tF.i_DS1_DEN"); 
	  system("$COPY_COMMAND tC0.o_DEN tF.i_DS2_DEN"); 
	  }
#
	elsif ($T1 eq 'G') {
	  $PSP1 = "$PSPS/13al.pspgth";
          $PSP2 = "$PSPS/33as.SGS_mod";
	  $OPTFLDIFF = '-medium';
	  unlink 'tG.i_WFK','tG.i_WFQ';
	  system("$COPY_COMMAND tC0.o_WFK tG.i_WFK");
	  system("$COPY_COMMAND tC0.o_WFK tG.i_WFQ");
	  }
#
	elsif ($T1 eq 'I') {
	  $PSP1 = "$PSPS/26fe.pspnc";
	  $OPTFLDIFF = '-easy';
          $CASE10=0;
	  }
#
        elsif ($T1 eq 'K') {
          $PSP1 = "$PSPS/31ga.pspnc";
          $PSP2 = "$PSPS/33as.pspnc";
          $OPTFLDIFF = '-medium';
          }
#
        elsif ($T1 eq 'L') {
          $PSP1 = "$PSPS/7n.pspnc";
          $PSP2 = "$PSPS/14si.pspnc";
          }
#
        elsif ($T1 eq 'M') {
          $PSP1 = "$PSPS/14si.pspnc";
          }
#
        elsif ($T1 eq 'N') {
          $PSP1 = "$PSPS/14si.pspnc";
          }
#
        elsif ($T1 eq 'O') {
          $PSP1 = "$PSPS/14si.pspnc";
          $CASE10=0;
          }
#This is to test MPI and paral_kgb, only with 4 procs in parallel
        elsif ($T1 eq 'P') {
          $PSP1 = "$PSPS/14si.phoney_mod";
          $OPTFLDIFF = '-easy';
          $CASE1= 0;
          $CASE2= 0;
          $CASE10=0;
          }
#This is to test MPI_IO norm conserving (with the MPI flag), only with 4 procs in parallel
        elsif ($T1 eq 'Q') {
          $PSP1 = "$PSPS/14si.phoney_mod";
	  $CASE0= 0;
	  $CASE1= 0;
	  $CASE2= 0;
	  $CASE10=0;
          }
#This is to test MPI and paral_kgb, only with 4 procs in parallel, paw case
        elsif ($T1 eq 'R') {
          $PSP1 = "$PSPS/6c_lda.paw";
          $OPTFLDIFF = '-easy';
          $CASE1= 0;
          $CASE2= 0;
          $CASE10=0;
          }
#This is to test MPI_IO PAW (with the MPI flag), only with 4 procs in parallel
        elsif ($T1 eq 'S') {
          $PSP1 = "$PSPS/6c_lda.paw";
	  $CASE0= 0;
	  $CASE1= 0;
	  $CASE2= 0;
	  $CASE10=0;
          }
#This is to test the recursion code, in sequential and with 4 procs in parallel
        elsif ($T1 eq 'T') {
          $PSP1 = "$PSPS/2he.psphgh";
	  $CASE1= 0;
	  $CASE2= 0;
	  $CASE10=0;
          }
#This is to test the KSS generation parallelized over k-points
        elsif ($T1 eq 'U') {
          $PSP1 = "$PSPS/14si.pspnc";
          $CASE10=0;
          }
#This is to test the SC-GW with advanced memory management in parallel
        elsif ($T1 eq 'V') {
          $PSP1 = "$PSPS/11na.pspnc";
          }
	$FLDreport = "fldiff.set$T1.report";
	unlink $FLDreport;
	$RUNfile = 'abfiles.run';
	if ($LAM == 1) {
		$cmd = "$LAMDIR/wipe ../cluster1";
		print $cmd if ($debug > 1);
		system($cmd);
		$cmd = "$LAMDIR/wipe ../cluster2";
		print $cmd if ($debug > 1);
		system($cmd);
		}
#This is to test MPI and paral_kgb, only with 4 procs in parallel, paw case
        elsif ($T1 eq 'W') {
          $PSP1 = "$PSPS/6c_lda.paw";
          $OPTFLDIFF = '-easy';
          $CASE1= 0;
          $CASE2= 0;
          $CASE10=0;
          }
	
# ****************************************
	foreach $case (@Cases) {
# Jump to first test to be run
		next if ($T2 eq '' && $case ne $FirstCase);
#
#Here, set characteristics of the test : the name of the test_case is $T1$T2 
		$T2 = $case;
		if ( ($T2 eq '0' && $CASE0 == 0) || ($T2 eq '1' && $CASE1 == 0) || 
                     ($T2 eq '2' && $CASE2 == 0) || ($T2 eq '4' && $CASE4 == 0) || 
		     ($T2 eq '10' && $CASE10 == 0) ) {
			last if ($T2 eq $LastCase); 
			next;
			}
#Set some variables according to the value of $T2
		$MODE = $T2 eq '0' ? 'SEQ' : 'PAR';
		$CLUST = "../cluster$T2";
		$NPAR = $T2;
		unlink $RUNfile;
		system ("$ERASE_COMMAND t$T1$T2*");	# unlink with metachar
#Initialize the files file
		open (FILE,">$RUNfile") || die "Unable to open FILE $RUNfile";
		print FILE "$abinit_inpdir/paral/Input/t$T1.in";
		print FILE "t$T1$T2.out";
		print FILE "t$T1.i";
		print FILE "t$T1$T2.o";
		print FILE "t$T1$T2";
		print FILE $PSP1;
		print FILE $PSP2;
		close (FILE);
#Run the job
		if ($MACH eq 'cox') {
			$MachFile = "--gm-f $CLUST";
			}
		elsif ($MPICH == 1 && $MODE eq 'PAR') { 
			$MachFile = "-machinefile $CLUST";
			}
		else {
			$MachFile = '';
			}
# if OPENMPI, start mpd
		if ($OPENMPI == 0) {
			$MPDSTAT=`$MPICH_DIR/mpdtrace`;
			if ( $MPDSTAT =~ !/^chum/ ) {
        			system("$MPICH_DIR/mpd &");
				}
			}

		if ($LAM == 1 && $MODE eq 'PAR') { 
			$cmd = "$LAMDIR/lamboot $CLUST ";
			print $cmd if ($debug > 1);
			system($cmd);
			}
#
		if ($MODE eq 'SEQ') {
			$CALL_CODE = ($T2 eq '0' && $MACH eq 'sr8k') ? "prun $CODE_SEQ" : $CODE_SEQ
			}
		else {
			$CALL_CODE = $MACH eq 'vpp' ? "$CODE_PAR $MPIRUN_NP $NPAR" : "$MPIRUN_NP $NPAR $MachFile $CODE_PAR"
			}
		$cmd = "$CALL_CODE < $RUNfile > t$T1$T2.log";
		print $cmd if ($debug > 1);
		system($cmd);
#
		if ($LAM == 1 && $MODE eq 'PAR') {
		$cmd = "$LAMDIR/wipe $CLUST";
		print $cmd if ($debug > 1);
		system($cmd);
		}
#Analyze the output
		$cmd = "$DIFF_COMMAND t$T1$T2.out $abinit_inpdir/paral/Refs/t$T1$NPAR.out > diff.t$T1$T2";
		print $cmd if ($debug > 1);
		system($cmd);
# append label with test number to report file
		open (FLDREP,">>$FLDreport") || die "Unable to open FLDREP for test $T1$T2";
		print FLDREP "\n","Case_$T1$T2";	# Case_NN
		close (FLDREP);
		$cmd = "$FLDIFF -ignore $OPTFLDIFF t$T1$T2.out $abinit_inpdir/paral/Refs/t$T1$NPAR.out Case_$T1$T2 >> $FLDreport";
		print $cmd if ($debug > 1);
		system($cmd);
#Eventually exit
		last if ($T2 eq $LastCase);
		}
# End of loop on $T1
	}
# ***************************************************************************
sub transpath {
	local($path) = @_;
#
# purpose: translate unix-like path of file to DOS-like according to $SLASH
# argument:
#	$path = path to be translated
# output: this subroutine acts as a function and returns the path
# according to host conventions
	
	$path =~ tr/\057/\\/ if ($SLASH eq '\\');
	return $path;
	}

# ****************************************
sub perlpath {
	local($path) = @_;
#
# purpose: translate unix-like path of file to DOS-like according to $PERLSLASH.
#   This is necessary when calling a DOS command like perl under PGI Workstation.
# argument:
#	$path = path to be translated
# output: this subroutine acts as a function and returns the path
# according to host conventions
	
	$path =~ tr/\057/\\/ if ($PERLSLASH eq '\\');
	return $path;
	}
