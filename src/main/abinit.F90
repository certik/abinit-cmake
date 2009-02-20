!{\src2tex{textfont=tt}}
!!****p* ABINIT/abinit
!! NAME
!! abinit
!!
!! FUNCTION
!! Main routine for conducting LDA calculations by self-consistent cycles.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, MKV, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! The new user is strongly adviced to read the
!! latest version of the file ~abinit/doc/users/new_user_guide.html
!! before trying to modify or even use the code.
!! Even experienced users of the code should also be careful in coding,
!! please read the latest version of the file ~abinit/doc/developers/rules_coding
!!
!! The present main routine drives the following operations :
!!
!! 1) Eventually initialize MPI
!! 2) Initialize overall timing of run
!! 3) Print greeting for interactive user and
!!    Read names of files (input, output, rootinput, rootoutput, roottemporaries),
!!    create the name of the status file, initialize the status subroutine.
!! 4) Open output file and print herald at top of output and log files
!! 5) Read the input file, and store the information in a long string of characters
!! 6) Take ndtset from the input string, then allocate
!!    the arrays whose dimensions depends only on ndtset and msym.
!! 7) Continue to analyze the input string, and allocate the remaining arrays.
!!    Also modulate the timing according to timopt.
!! 8) Finish to read the "file" file completely,
!!    and also initialize pspheads (the pseudopotential header information)
!! 9) Provide defaults for the variables that have not yet been initialized.
!! 10) Call the main input routine.
!! 11) Echo input data to output file and log file
!! 12) Perform additional checks on input data
!!  At this stage, all the information from the "files" file and "input" file
!!  have been read and checked.
!! ___________________________________________
!! 13) Perform main calculation  (call driver)
!! -------------------------------------------
!!
!! 14) Give final echo of coordinates, etc.
!! 15) Timing analysis
!! 16) Delete the status file, and, for build-in tests,
!!       analyse the correctness of results
!! 17) Write the final timing, close the output file, and write a final line
!!       to the log file
!! 18) Eventual cleaning of MPI run
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! PARENTS
!!
!! CHILDREN
!!      chkinp,chkvars,date_and_time,driver,dtsetfree,herald,indefo
!!      initpapichoice,intagm,invars0,invars1m,invars2m,iofn1,iofn2,leave_new
!!      leave_test,mpi_allreduce,mpi_comm_rank,mpi_comm_size,mpi_finalize
!!      mpi_init,outvars,papi_init,parsefile,status,testfi,timab,timana,timein
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

program abinit

 use defs_basis
 use defs_datatypes
 use defs_infos
#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12parser
 use interfaces_13iovars
 use interfaces_18seqpar
 use interfaces_21drive
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments -----------------------------------

!Local variables-------------------------------
!fix for Fortran interfaces
!scalars
 character(len=30) :: token
!no_abirules
!
!===============================================================================
!  abinit_version designate overall code version
!  mpw=maximum number of planewaves in basis sphere
!  unkg,unkgq,unkg1,unkg,unwff1,unwff2,unwffgs,unwfkq,
!  unwft1,unwft2,unylm,unylm1,unpaw,unpaw1 and unpawq...
!!  define input and output unit numbers.
!   These unit numbers are transferred
!   down to the adequate routines.
!   Other unit numbers (ab_in,ab_out,std_out,tmp_unit)
!   have been defined in defs_basis.f .
!  The array filnam is used for the name of input and output files,
!  and roots for generic input, output or temporary files.
!  Pseudopotential file names are set in iofn2, and are contained in pspheads.
!  The name filstat will be needed beyond gstate to check
!  the appearance of the "exit" flag, to make a hasty exit, as well as
!  in order to output the status of the computation.
!==============================================================================
! Declarations
! Define input and output unit numbers (do not forget, unit 5 and 6
! are standard input and output)
! Also, unit number 21, 22 and 23 are used in nstdy3, for the 3 dot
! wavefunctions. Others unit numbers will be used in the case
! of the variational and 2n+1 expressions.
! In defs_basis, one defines :
!  std_in=5, ab_in=5, std_out=6, ab_out=7, tmp_unit=9, tmp_unit2=10
 integer,parameter :: unchi0=42,unddb=16,unddk=15,unkg1=19,unkg=17,unkgq=18
 integer,parameter :: unpaw=26,unpaw1=27,unpawq=28,unpos=30
 integer,parameter :: unwff1=1,unwff2=2,unwffgs=3,unwfkq=4,unwft1=11
 integer,parameter :: unwft2=12,unwftgs=13,unwftkq=14,unylm=24,unylm1=25
 integer,parameter :: unkss=40,unscr=41,unqps=43
! Define "level of the routine", for debugging purposes
 integer,parameter :: level=1
 integer :: choice,dmatpuflag,idtset,iexit,ii,ios,iounit,istatr,istatshft,jdtset,lenstr,marr
 integer :: msym,mu,mxlpawu,mxmband,mxmband_upper,mxnatom,mxnatsph,mxnatpawu,mxncenter
 integer :: mxnatvshift,mxnconeq,mxnkptgw
 integer :: mxnkpt,mxnorb,mxnnos,mxnqptdm
 integer :: mxnspinor,mxnsppol,mxnsym,mxntypat,natom,ndtset,ndtset_alloc,nfft,nkpt,npsp
 integer :: nsppol,option,prtvol,timopt,papiopt,tread,usepaw
 integer :: intarr(1)
 integer,allocatable :: bravais_(:,:),mband_upper_(:),nband(:),npwtot(:)
 real(dp) :: cpui,etotal,walli,zion_max
 real(dp) :: dprarr(1),strten(6),tsec(2)
 real(dp),allocatable :: fred(:,:),xred(:,:)
 character(len=24) :: codename
 character(len=500) :: message
 character(len=strlen) :: string,string_raw
 character(len=fnlen) :: filstat
 character(len=fnlen) :: filnam(5)
!character(len=30) :: token
 type(datafiles_type) :: dtfil
 type(dataset_type),allocatable :: dtsets(:)
 type(MPI_type) :: mpi_enreg
 type(pspheader_type),allocatable :: pspheads(:)
 type(results_out_type),allocatable :: results_out(:)
            integer :: ierr
#if defined HAVE_NETCDF
 logical :: netcdf_input
#endif
 logical :: xml_output
 integer :: values(8)
 character(len=5) :: strzone
 character(len=8) :: strdat
 character(len=10) :: strtime
!no_abirules
#if defined MPI
          !Variables introduced for MPI version
           real(dp) :: tsec_s(2)
#endif

!******************************************************************
!BEGIN EXECUTABLE SECTION
!
!1) Eventually initialize MPI : one should write a separate routine -init_mpi_enreg- for doing that !!

!Default for sequential use
 mpi_enreg%world_comm=0
 mpi_enreg%world_group=0
 mpi_enreg%me=0
 mpi_enreg%nproc=1
 mpi_enreg%num_group_fft = 0 ! in some cases not initialized but referenced in xdef_comm.F90
 mpi_enreg%paral_compil=0
 mpi_enreg%paral_compil_mpio=0

!Initialize MPI
#if defined MPI
           call MPI_INIT(ierr)
           mpi_enreg%world_comm=MPI_COMM_WORLD
           mpi_enreg%world_group=MPI_GROUP_NULL
           call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_enreg%me,ierr)
           call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_enreg%nproc,ierr)
           write(6,*)' abinit : nproc,me=',mpi_enreg%nproc,mpi_enreg%me
           mpi_enreg%paral_compil=1
#endif

!Signal MPI I/O compilation has been activated
#if defined MPI_IO
           mpi_enreg%paral_compil_mpio=1
           if(mpi_enreg%paral_compil==0)then
            write(message,'(6a)') ch10,&
&            ' abinit : ERROR -',ch10,&
&            '  In order to use MPI_IO, you must compile with the MPI flag ',ch10,&
&            '  Action : recompile your code with different CPP flags.'
            call wrtout(06,message,'COLL')
            call leave_new('COLL')
           end if
#endif

!Initialize spaceComm, used in leave_test
 mpi_enreg%spaceComm=mpi_enreg%world_comm
!Initialize paral_compil_kpt, actually always equal to paral_compil
!(paral_compil_kpt should be suppressed after big cleaning)
 mpi_enreg%paral_compil_kpt=0
 if(mpi_enreg%paral_compil==1) mpi_enreg%paral_compil_kpt=1

!Other values of mpi_enreg are dataset dependent, and should NOT be initialized
!inside abinit.F90 .
!XG 071118 : At present several other values are
!initialized temporarily inside invars1.F90, FROM THE DTSET
!VALUES. In order to releave the present constraint of having mpi_enreg
!equal for all datasets, they should be reinitialized from the dtset values
!inside invars2m.F90 (where there is a loop over datasets, and finally,
!reinitialized from the dataset values inside each big routine called by driver,
!according to the kind of parallelisation that is needed there.
!One should have one init_mpi_dtset routine (or another name) per big routine (well, there is also
!the problem of TDDFT ...). Also, one should have a clean_mpi_dtset called at the end
!of each big routine, as well as invars1.F90 or invars2m.F90 .

!2) Initialize overall timing of run:
#ifdef HAVE_PAPI
 call papi_init()
#endif

 call timein(cpui,walli)

 call timab(1,0,tsec)

!Start to accumulate time for the entire run.
!The end of accumulation is in timana.f
 call timab(1,1,tsec)

!3) Print greeting for interactive user,
!read names of files (input, output, rootinput, rootoutput, roottemporaries),
!create the name of the status file, initialize the status subroutine.

 call timab(41,1,tsec)

 call iofn1(filnam,filstat,mpi_enreg)

 iexit=99
 call status(0,filstat,iexit,level,'first status  ')

!4) Open output file and print herald at top of output and log files

 if(mpi_enreg%me==0)then

  open (unit=ab_out,file=filnam(2),form='formatted',status='new')
  rewind (unit=ab_out)
  codename='ABINIT'//repeat(' ',18)
  call herald(codename,abinit_version,ab_out)
  call herald(codename,abinit_version,std_out)
! Write names of files
  write(message, '(a,a,a,a,a,a,a,a,a,a,a,a)' )&
&  '- input  file    -> ',trim(filnam(1)),ch10,&
&  '- output file    -> ',trim(filnam(2)),ch10,&
&  '- root for input  files -> ',trim(filnam(3)),ch10,&
&  '- root for output files -> ',trim(filnam(4)),ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')
 end if

!5) Read the file, stringify it and return the number of datasets.
 call status(0,filstat,iexit,level,'call instrng  ')
 call parsefile(filnam(1), lenstr, ndtset, string)

 ndtset_alloc=ndtset ; if(ndtset==0)ndtset_alloc=1
 allocate(dtsets(0:ndtset_alloc))
!Here, the default msym.
 msym=384

!7) Continue to analyze the input string, get upper dimensions,
!and allocate the remaining arrays.
!Also modulate the timing according to timopt.

 call status(0,filstat,iexit,level,'call invars0  ')

 call invars0(dtsets,ab_out,istatr,istatshft,lenstr,&
& mxnatom,mxntypat,ndtset,ndtset_alloc,npsp,string)

!DEBUG
!write(6,*)' abinit : npsp,mxntypat=',npsp,mxntypat
!ENDDEBUG


!Read timopt and pass it to timab
 timopt=1 ; marr=1
 if(mpi_enreg%paral_compil==1)timopt=0
 token = 'timopt'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) timopt=intarr(1)

!Read papiopt
 call  initpapichoice(string(1:lenstr),lenstr,papiopt)

 call timab(41,2,tsec)


 call timab(timopt,5,tsec)

!8) Finish to read the "file" file completely, as npsp is known,
!and also initialize pspheads, that contains the important information
!from the pseudopotential headers, as well as the psp filename

 call status(0,filstat,iexit,level,'call iofn2    ')

 call timab(42,1,tsec)

 allocate(pspheads(npsp))
 call iofn2(npsp,pspheads,mpi_enreg)

!If (all) pspcod are 7 then this is a PAW calculation
 usepaw=0;if(minval(abs(pspheads(1:npsp)%pspcod-7))==0) usepaw=1
 do idtset=0,ndtset_alloc
  dtsets(idtset)%usepaw=usepaw
 end do

!Take care of other dimensions, and part of the content of dtsets
!that is or might be needed early.
!zion_max=maxval(pspheads(1:npsp)%zionpsp) ! This might not work properly with HP compiler
 zion_max=pspheads(1)%zionpsp
 do ii=1,npsp
  zion_max=max(pspheads(ii)%zionpsp,zion_max)
 end do

 allocate(bravais_(11,0:ndtset_alloc))
 allocate(mband_upper_ (  0:ndtset_alloc))

 call invars1m(bravais_,dmatpuflag,dtsets,ab_out,lenstr,mband_upper_,mpi_enreg,&
& msym,mxlpawu,mxmband_upper,mxnatom,mxnatpawu,mxnatsph,mxnatvshift,mxncenter,mxnconeq,mxnkptgw,&
& mxnkpt,mxnorb,mxnnos,mxnqptdm,mxnspinor,mxnsppol,mxnsym,ndtset,ndtset_alloc,string,&
& zion_max)

!9) Provide defaults for the variables that have not yet been initialized.

!Be careful : at this fifth call of status, istatr taken
!from inprep will be saved definitively.
 call status(0,filstat,istatr,level,'call indefo   ')

 call indefo(dtsets,mpi_enreg,ndtset_alloc)

!If all the pseudopotentials have the same pspxc, override the default
!value for dtsets 1 to ndtset
 if(minval(abs((pspheads(1:npsp)%pspxc-pspheads(1)%pspxc)))==0)then
  dtsets(1:ndtset_alloc)%ixc=pspheads(1)%pspxc
 end if

 call timab(42,2,tsec)

!10) Call the main input routine.

 call status(0,filstat,istatshft,level,'call inmain   ')

 call timab(43,1,tsec)

 call invars2m(bravais_,dtsets,ab_out,lenstr,&
& mband_upper_,mpi_enreg,msym,ndtset,ndtset_alloc,npsp,pspheads,string)
 deallocate(bravais_, mband_upper_)

!mxmband=maxval(dtsets(1:ndtset_alloc)%mband) ! THis might not wotk with the HP compiler
 mxmband=dtsets(1)%mband
 do ii=1,ndtset_alloc
  mxmband=max(dtsets(ii)%mband,mxmband)
 end do

 call timab(43,2,tsec)

!11) Echo input data to output file and log file

 call status(0,filstat,iexit,level,'call outvars(1')

 call timab(44,1,tsec)

!For evolving variables, and results
 allocate(results_out(0:ndtset_alloc))

 do idtset=0,ndtset_alloc
  allocate(results_out(idtset)%fcart(3,mxnatom))
  allocate(results_out(idtset)%fred(3,mxnatom))
  allocate(results_out(idtset)%npwtot(mxnkpt))
  allocate(results_out(idtset)%occ(mxmband_upper*mxnkpt*mxnsppol))
  allocate(results_out(idtset)%vel(3,mxnatom))
  allocate(results_out(idtset)%xred(3,mxnatom))
  results_out(idtset)%natom      =dtsets(idtset)%natom
  results_out(idtset)%acell(:)   =dtsets(idtset)%acell_orig(:)
  results_out(idtset)%occ(:)     =dtsets(idtset)%occ_orig(:)
  results_out(idtset)%rprim(:,:) =dtsets(idtset)%rprim_orig(:,:)
  results_out(idtset)%rprimd(:,:)=dtsets(idtset)%rprimd_orig(:,:)
  results_out(idtset)%vel(:,:)   =dtsets(idtset)%vel_orig(:,:)
  results_out(idtset)%xred(:,:)  =dtsets(idtset)%xred_orig(:,:)
  results_out(idtset)%etotal     =zero
  results_out(idtset)%fred(:,:)  =zero
  results_out(idtset)%fcart(:,:) =zero
  results_out(idtset)%npwtot(:)  =0
  results_out(idtset)%strten(:)  =zero
 end do

 if(mpi_enreg%me==0) then

! Echo input to output file on unit ab_out, and to log file on unit 06 :
  choice=1
  do ii=1,2
   if(ii==1)iounit=ab_out
   if(ii==2)iounit=6

!  DEBUG
!  write(6,*)' abinit : before outvars, ndtset_alloc=',ndtset_alloc
!  ENDDEBUG

   call outvars (choice,dmatpuflag,dtsets,iounit,istatr,istatshft,mpi_enreg,&
&   mxlpawu,mxmband,mxnatom,mxnatpawu,mxnatsph,mxnatvshift,mxncenter,mxnconeq,mxnkptgw,mxnkpt,&
&   mxnorb,mxnnos,mxnqptdm,mxnspinor,mxnsppol,mxnsym,mxntypat,&
&   ndtset,ndtset_alloc,npsp,pspheads,results_out,timopt)
  end do

  if (dtsets(1)%outputXML == 1) then
!  ABINIT has been compiled with XML output support, then we open the
!  channel of the XML output file.
   open(unit = ab_xml_out, file = trim(filnam(4))//".xml", form = "formatted", &
&   action = "write")
   write(ab_xml_out, "(A)") '<?xml version="1.0" encoding="utf-8"?>'
!  DTD declaration : root element is "abinitRun", and DTD file is not public
!  but given in the distribution by the file abinitRun.dtd.
   write(ab_xml_out, "(A)") '<!DOCTYPE abinitRun SYSTEM "extras/post_processing/abinitRun.dtd">'
!  Creating the root markup.
   write(ab_xml_out, "(A)") '<abinitRun>'
!  Store the launch date for the timeInfo markup (at the end).
   call date_and_time(strdat,strtime,strzone,values)
   xml_output = .true.
  else
   xml_output = .false.
  end if

! End of me==0 section
 end if

!This synchronization is not strictly needed, but without it,
!there are problems with Tv1#93 in parallel, PGI compiler, on Intel/PC
 call leave_test(mpi_enreg)

!DEBUG
!write(6,*)' abinit : after outvars'
!ENDDEBUG

 call timab(44,2,tsec)

!12) Perform additional checks on input data

 call status(0,filstat,iexit,level,'call chkinp   ')

 call timab(45,1,tsec)

 call chkinp(dtsets,ab_out,mpi_enreg,ndtset,ndtset_alloc,npsp,pspheads)

!DEBUG
!write(6,*)' abinit : after chkinp'
!stop
!ENDDEBUG

!Check whether the string only contains valid keywords
 call chkvars(string)

!DEBUG
!write(6,*)' abinit : after chkvars'
!stop
!ENDDEBUG

!At this stage, all the information from the "files" file and "input" file
!have been read and checked.

 call timab(45,2,tsec)

!13) Perform main calculation
!The timing is done in gstate

 call status(0,filstat,iexit,level,'call driver   ')

 prtvol=dtsets(1)%prtvol
 if(prtvol==-level)then
  write(message,'(a1,a,a1,a,i1,a)') ch10,' abinit : before driver ',&
&  ch10,'  prtvol=-',level,', debugging mode => stop '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!Unit numbers are transferred to dtfil
 dtfil%unchi0 =unchi0
 dtfil%unddb  =unddb
 dtfil%unddk  =unddk
 dtfil%unkg   =unkg
 dtfil%unkgq  =unkgq
 dtfil%unkg1  =unkg1
 dtfil%unkss  =unkss
 dtfil%unqps  =unqps
 dtfil%unscr  =unscr
 dtfil%unwff1 =unwff1
 dtfil%unwff2 =unwff2
 dtfil%unwffgs=unwffgs
 dtfil%unwffkq=unwfkq
 dtfil%unwft1 =unwft1
 dtfil%unwft2 =unwft2
 dtfil%unwftgs=unwftgs
 dtfil%unwftkq=unwftkq
 dtfil%unylm  =unylm
 dtfil%unylm1 =unylm1
 dtfil%unpaw  =unpaw
 dtfil%unpaw1 =unpaw1
 dtfil%unpawq =unpawq
 dtfil%unpos  =unpos

 call driver(abinit_version,cpui,dtfil,dtsets,filnam,filstat,&
& mpi_enreg,ndtset,ndtset_alloc,npsp,pspheads,results_out,walli)
 call status(0,filstat,iexit,level,'after driver  ')

!14) Give final echo of coordinates, etc.

 call timab(46,1,tsec)

 write(message,'(a,a,a,62a,80a)') ch10,&
& '== END DATASET(S) ',('=',mu=1,62),ch10,('=',mu=1,80)
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,message,'COLL')

 if(mpi_enreg%me==0) then

! Echo input to output file on unit ab_out, and to log file on unit 06
  choice=2
  do ii=1,2
   if(ii==1)iounit=ab_out
   if(ii==2)iounit=6
   write(iounit,*)' '
   call outvars (choice,dmatpuflag,dtsets,iounit,istatr,istatshft,mpi_enreg,&
&   mxlpawu,mxmband,mxnatom,mxnatpawu,mxnatsph,mxnatvshift,mxncenter,mxnconeq,mxnkptgw,mxnkpt,&
&   mxnorb,mxnnos,mxnqptdm,mxnspinor,mxnsppol,mxnsym,mxntypat,&
&   ndtset,ndtset_alloc,npsp,pspheads,results_out,timopt)
   if(ii==2)write(06,*)' '
  end do

! DEBUG
! write(6,*)' abinit : after outvars '
! stop
! ENDDEBUG

 end if ! mpi_enreg%me==0

!In prevision of the next two calls, some variables need to be transfered.
!They concern the case ndtset<2, so take first value.
 natom=dtsets(1)%natom ; nkpt=dtsets(1)%nkpt ; nsppol=dtsets(1)%nsppol
 nfft=dtsets(1)%nfft   ; etotal=results_out(1)%etotal
 allocate(nband(nkpt*nsppol),npwtot(nkpt),fred(3,natom),xred(3,natom))
 fred(:,:)  =results_out(1)%fred(:,1:natom)
 nband(:)   =dtsets(1)%nband(1:nkpt*nsppol)
 npwtot(:)  =results_out(1)%npwtot(1:nkpt)
 strten(:)  =results_out(1)%strten(:)
 xred(:,:)  =results_out(1)%xred(:,1:natom)

 call timab(46,2,tsec)

!15) Timing analysis

 if(timopt/=0)then
  call status(0,filstat,iexit,level,'call timana   ')
  call timana (mpi_enreg, natom, nband, ndtset, nfft, nkpt, npwtot, nsppol, timopt, papiopt)
 else
#if defined MPI
            if(mpi_enreg%me==0)then
!            This is for the automatic tests
             write(ab_out,'(5a)')ch10,ch10,&
&             '- Timing analysis has been suppressed with timopt=0',ch10,ch10
            end if
#endif
 end if

!16) Delete the status file, and, for build-in tests,
!analyse the correctness of results.

!DEBUG
!write(6,*)' abinit : before testfi '
!stop
!ENDDEBUG

 if(ndtset==0)then
  call testfi(etotal,filnam,filstat,fred,natom,strten,xred)
 else
  open (tmp_unit,file=trim(filstat),form='formatted',status='unknown')
  close (tmp_unit,status='delete')
 end if

!One should have here the explicit deallocation of all arrays
 do idtset=0,ndtset_alloc
  deallocate(results_out(idtset)%fcart)
  deallocate(results_out(idtset)%fred)
  deallocate(results_out(idtset)%npwtot)
  deallocate(results_out(idtset)%occ)
  deallocate(results_out(idtset)%vel)
  deallocate(results_out(idtset)%xred)
 end do

 do idtset=0,ndtset_alloc
  call dtsetfree(dtsets(idtset))
 end do

 deallocate(dtsets,fred,nband,npwtot)
 deallocate(pspheads,results_out,xred)

!17) Write the final timing, close the output file, and write a final line
!to the log file

 call timein(tsec(1),tsec(2))
 tsec(1)=tsec(1)-cpui
 tsec(2)=tsec(2)-walli

#if defined MPI
           write(6, '(a,i4,a,f13.1,a,f13.1)' )&
&           ' Proc.',mpi_enreg%me,' individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
           if(mpi_enreg%me==0)then
            write(ab_out, '(a,a,a,i4,a,f13.1,a,f13.1)' ) &
&            '-',ch10,'- Proc.',mpi_enreg%me,' individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
           end if
           call MPI_ALLREDUCE(tsec,tsec_s,2,MPI_DOUBLE_PRECISION,MPI_SUM,&
&           MPI_COMM_WORLD,ierr)
           tsec=tsec_s

#endif

 write(message, '(a,80a,a,a,a)' ) ch10,('=',mu=1,80),ch10,ch10,&
& ' Calculation completed.'
 call wrtout(ab_out,message,'COLL')

 if(mpi_enreg%me==0)then
  write(ab_out, '(a,f13.1,a,f13.1)' ) &
&  '+Overall time at end (sec) : cpu=',tsec(1),'  wall=',tsec(2)
  close (unit=ab_out)

  if (xml_output) then
   write(ab_xml_out, "(A)") '  <!-- timeInfo markup : cpu and wall attributes are given in seconds. -->'
   write(ab_xml_out, "(A,I0,A,I0,A,I0,A)", advance = "NO") &
&   '  <timeInfo day="', values(3) , &
&   '" month="', values(2) ,'" year="', values(1) ,'" '
   write(message, *) tsec(1)
   write(ab_xml_out, "(A,A,A)", advance = "NO") 'cpu="', trim(message) ,'"'
   write(message, *) tsec(2)
   write(ab_xml_out, "(A,A,A)") ' wall="', trim(message) ,'" />'
   write(ab_xml_out, "(A)") "</abinitRun>"
!  ABINIT has been compiled with XML output support, then we close the
!  channel of the XML output file.
   close(unit = ab_xml_out)
  end if

  write(message, '(a,a)' ) ch10,' Calculation completed.'
  call wrtout(06,message,'COLL')
 end if

!18) Eventual cleaning of MPI run

#if defined MPI
           call MPI_FINALIZE(ierr)
#endif

 end program abinit
!!***
