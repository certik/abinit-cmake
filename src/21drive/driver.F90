!{\src2tex{textfont=tt}}
!!****f* ABINIT/driver
!! NAME
!! driver
!!
!! FUNCTION
!! Driver for ground state, response function, susceptibility, screening
!! and sigma calculations. The present routine drives the following operations.
!! An outer loop allows computation related to different data sets.
!! For each data set, either a GS calculation, a RF calculation,
!! a SUS calculation, a SCR calculation or a SIGMA calculation is made.
!! In both cases, the input variables are transferred in the proper variables,
!! selected big arrays are allocated, then the gstate, respfn or suscep subroutines are called.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG,MKV,MM,MT,FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! codvsn= code version
!! cpui=initial CPU time
!! dtsets(0:ndtset_alloc)=<type datasets_type>contains all input variables
!! filnam(5)=character strings giving file names
!! filstat=character strings giving name of status file
!! mpi_enreg=informations about MPI parallelization
!! ndtset=number of datasets
!! ndtset_alloc=number of datasets, corrected for allocation of at
!!               least one data set.
!! npsp=number of pseudopotentials
!! pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!! walli=initial wall clock time
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  Input/Output
!! dtfil=<type datafiles_type>infos about file names, file unit numbers
!!  (part of which were initialized previously)
!! results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including evolving variables
!!   Default values are set up in the calling routine
!!
!! NOTES
!! The array filnam is used for the name of input and output files,
!! and roots for generic input, output or temporary files.
!! Pseudopotential file names are set in pspini and pspatm,
!! using another name. The name filstat will be needed beyond gstate to check
!! the appearance of the "exit" flag, to make a hasty exit, as well as
!! in order to output the status of the computation.
!!
!! TODO
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      chkdilatmx,chkexi,getdim_nloc,gstate,leave_new,mkfilename,mkrdim
!!      nonlinear,pstate,respfn,screening,sigma,status,suscep,timab,wannier
!!      wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine driver(codvsn,cpui,dtfil,dtsets,filnam,filstat,&
& mpi_enreg,ndtset,ndtset_alloc,npsp,pspheads,results_out,walli)

 use defs_basis
 use defs_parameters
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_13iovars
 use interfaces_13paw
 use interfaces_13psp
 use interfaces_17suscep
 use interfaces_18seqpar
 use interfaces_21drive, except_this_one => driver
 use interfaces_21paral_md
 use interfaces_21rdm
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndtset,ndtset_alloc,npsp
 real(dp),intent(in) :: cpui,walli
 character(len=6),intent(in) :: codvsn
 character(len=fnlen),intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(inout) :: dtfil
!arrays
 character(len=fnlen),intent(in) :: filnam(5)
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(pspheader_type),intent(in) :: pspheads(npsp)
 type(results_out_type),intent(inout) :: results_out(0:ndtset_alloc)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=2
 integer,save :: dimekb_old=-1,lmnmax_old=-1,lnmax_old=-1,mqgridff_old=0
 integer,save :: mqgridvl_old=0,ntypat_old=-1,paw_size_old=-1,usepaw_old=-1
 integer :: ceksph,getcell,getocc,getvel,getxcart,getxred,idtset,iexit,iget,ii
 integer :: ilang,ipsp,ipspalch,ireadden,ireadwf,itypalch,itypat,jdtset
 integer :: jdtset_status,l_max,l_size_max,lmnmax,lmnmaxso,lnmax,lnmaxso,mband
 integer :: mdtset,mgfft,mk1mem,mkmem,mkqmem,mpsang,mpssoang,mpw,mtypalch,mu
 integer :: n1xccc,natom,nfft,nkpt,nspden,nsppol,nsym,openexit,paw_size,prtvol
 integer :: usepaw,will_read
 real(dp) :: etotal
 character(len=2) :: appen
 character(len=4) :: stringfile
 character(len=500) :: message
 character(len=9) :: stringvar
 character(len=fnlen) :: filchi0,fildens1in,fildensin,filkss,filqps,filscr
 character(len=fnlen) :: filvhain,fnamewff1,fnamewffddk,fnamewffk,fnamewffq
 type(dataset_type) :: dtset
 type(pawang_type) :: pawang
 type(pseudopotential_type) :: psps
 type(results_gs_type) :: results_gs
 type(vardims_type) :: abidims
!arrays
 integer :: mkmems(3)
 integer,allocatable :: jdtset_(:),npwtot(:)
 real(dp) :: acell(3),rprim(3,3),rprimd(3,3),rprimdget(3,3),strten(6),tsec(2)
 real(dp),allocatable :: occ(:),vel(:,:),xcart(:,:),xred(:,:),xredget(:,:)
 character(len=fnlen) :: filnam_ds(5)
 type(pawrad_type),allocatable :: pawrad(:)
 type(pawtab_type),allocatable :: pawtab(:)

!******************************************************************

!DEBUG ! Do not comment this line, for the time being XG020913
!write(6,*)' driver : enter '
!stop
!ENDDEBUG

 call timab(100,1,tsec)
 call status(0,filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 prtvol=dtsets(1)%prtvol
 if(prtvol==-level)then
  write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,&
&  ' driver : enter , debug mode '
  call wrtout(06,message,'COLL')
 end if

 mdtset=99

 if(ndtset>mdtset)then
  write(message, '(a,a,a,a,i2,a,i5,a)' )ch10,&
&  ' driver : BUG ',ch10,&
&  '  The maximal allowed ndtset is ',mdtset,&
&  ' while the input value is ',ndtset,'.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!DEBUG ! Do not comment this line for the time being XG020913
!write(6,*)' driver : before mtypalch '
!ENDDEBUG

!mtypalch=maxval(dtsets(1:ndtset_alloc)%ntypalch) ! Likely troubles with HP compiler
 mtypalch=dtsets(1)%ntypalch
 do ii=1,ndtset_alloc
  mtypalch=max(dtsets(ii)%ntypalch,mtypalch)
 end do

!DEBUG ! Do not comment this line for the time being XG020913
!write(6,*)' driver : before allocate, npsp= ',npsp
!ENDDEBUG

!Allocation of some arrays independent of the dataset
 allocate(psps%filpsp(npsp))
 allocate(psps%pspcod(npsp))
 allocate(psps%pspdat(npsp))
 allocate(psps%pspso(npsp))
 allocate(psps%pspxc(npsp))
 allocate(psps%title(npsp))
 allocate(psps%zionpsp(npsp))
 allocate(psps%znuclpsp(npsp))
 call psp2params_init(psps%gth_params, npsp)

 psps%filpsp(1:npsp)=pspheads(1:npsp)%filpsp
 psps%pspcod(1:npsp)=pspheads(1:npsp)%pspcod
 psps%pspdat(1:npsp)=pspheads(1:npsp)%pspdat
 psps%pspso(1:npsp)=pspheads(1:npsp)%pspso
 psps%pspxc(1:npsp)=pspheads(1:npsp)%pspxc
 psps%title(1:npsp)=pspheads(1:npsp)%title
 psps%zionpsp(1:npsp)=pspheads(1:npsp)%zionpsp
 psps%znuclpsp(1:npsp)=pspheads(1:npsp)%znuclpsp

 allocate(jdtset_(0:ndtset))
 if(ndtset/=0)then
  jdtset_(:)=dtsets(0:ndtset)%jdtset
 else
  jdtset_(0)=0
 end if

!*********************************************************************
!Big loop on datasets

!Do loop on idtset (allocate statements are present)
 do idtset=1,ndtset_alloc

  jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=1

  if(ndtset>=2)then
   jdtset_status=jdtset
  else
   jdtset_status=0
  end if

  call status(jdtset_status,filstat,iexit,level,'loop jdtset   ')

  write(message,'(a,80a,a,a,i2,a,66a,a)') ch10,&
&  ('=',mu=1,80),ch10,&
&  '== DATASET ',jdtset,' ',('=',mu=1,66),ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,message,'PERS')     ! PERS is choosen to make debugging easier

! Determine here wether the calculation is PAW
  usepaw  =0
  if (pspheads(1)%pspcod==7) usepaw=1  ! If paw, all pspcod necessarily are 7 (see iofn2)

! Copy input variables into a local dtset.
  call dtsetCopy(dtset, dtsets(idtset))

! DEBUG
! write(6,*)' driver : dtset%acell_orig(:)=',dtset%acell_orig(:)
! write(6,*)' driver : dtset%rprim_orig(:,:)=',dtset%rprim_orig(:,:)
! ENDDEBUG

! Set other values
  dtset%jdtset = jdtset
  dtset%ndtset = ndtset

! Copy input values
  acell(:)         = dtset%acell_orig(:)
  mkmems(1)        = dtset%mkmem
  mkmems(2)        = dtset%mkqmem
  mkmems(3)        = dtset%mk1mem
  psps%optnlxccc   = dtset%optnlxccc
  psps%mqgrid_ff   = dtset%mqgrid
  if (usepaw == 1) then
   psps%mqgrid_vl = dtset%mqgriddg
  else
   psps%mqgrid_vl = dtset%mqgrid
  end if
  rprim(:,:)       = dtset%rprim_orig(:,:)

! BEGIN TF_CHANGES
  mpi_enreg%paral_compil_respfn=dtset%paral_rf
  mpi_enreg%ngroup_respfn=dtset%ngroup_rf
! END TF_CHANGES

! Allocate arrays
  allocate(occ(dtset%mband*dtset%nkpt*dtset%nsppol))
  allocate(vel(3,dtset%natom) )
  allocate(xred(3,dtset%natom))

  occ   (:)   = dtset%occ_orig(1:dtset%mband*dtset%nkpt*dtset%nsppol)
  vel   (:,:) = dtset%vel_orig(:,1:dtset%natom)
  xred  (:,:) = dtset%xred_orig(:,1:dtset%natom)

! ****************************************************************************
! Treat the file names (get variables)


  filnam_ds(1:5)=filnam(1:5)

! If multi dataset mode, special treatment of filenames 3 and 4 (density and
! wavefunctions input and output, as well as other output files)
  if(ndtset>0)then
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   filnam_ds(3)=trim(filnam(3))//'_DS'//appen
   filnam_ds(4)=trim(filnam(4))//'_DS'//appen
!  DEBUG
!  write(6,*)' filnam_ds(3)',filnam_ds(3)
!  write(6,*)' filnam_ds(4)',filnam_ds(4)
!  ENDDEBUG
  end if

! According to getwfk and irdwfk, build _WFK file name, referred as fnamewffk
  stringfile='_WFK' ; stringvar='wfk'
  call mkfilename(filnam,fnamewffk,dtset%getwfk,idtset,dtset%irdwfk,jdtset_,&
&  ndtset,stringfile,stringvar,will_read)

  if(dtset%optdriver/=1)ireadwf=will_read
  if(ndtset/=0 .and. dtset%optdriver==1 .and. will_read==0)then
   write(message, '(a,a,a,a,a,a,a,a,i3,a,a,a,i3,a,i3,a,a,a)' )ch10,&
&   ' driver : ERROR -',ch10,&
&   '  At least one of the input variables irdwfk and getwfk ',ch10,&
&   '  must refer to a valid _WFK file, in the response function',ch10,&
&   '  case, while for idtset=',idtset,',',ch10,&
&   '  they are irdwfk=',dtset%irdwfk,', and getwfk=',dtset%getwfk,'.',ch10,&
&   '  Action : correct irdwfk or getwfk in your input file.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! Treatment of the other get wavefunction variable, if response function case
! or nonlinear case
  if ((dtset%optdriver==1).or.(dtset%optdriver==5)) then

!  According to getwfq and irdwfq, build _WFQ file name, referred as fnamewffq
   stringfile='_WFQ' ; stringvar='wfq'
   call mkfilename(filnam,fnamewffq,dtset%getwfq,idtset,dtset%irdwfq,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)
!  If fnamewffq is not initialized thanks to getwfq or irdwfq, use fnamewffk
   if(will_read==0)fnamewffq=fnamewffk

!  According to get1wf and ird1wf, build _1WF file name, referred as fnamewff1
   stringfile='_1WF' ; stringvar='1wf'
   call mkfilename(filnam,fnamewff1,dtset%get1wf,idtset,dtset%ird1wf,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)
   ireadwf=will_read

!  According to getddk and irdddk, build _1WF file name, referred as fnamewffddk
   stringfile='_1WF' ; stringvar='ddk'
   call mkfilename(filnam,fnamewffddk,dtset%getddk,idtset,dtset%irdddk,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)

  end if ! optdriver is 1 or 5

! According to getden, build _DEN file name, referred as fildensin
! A default is available if getden is 0
  stringfile='_DEN' ; stringvar='den'
  call mkfilename(filnam,fildensin,dtset%getden,idtset,0,jdtset_,&
&  ndtset,stringfile,stringvar,will_read)
  if(will_read==0)fildensin=trim(filnam_ds(3))//'_DEN'
  ireadden=will_read
! if (optdriver==0.and.ireadwf/=0) ireadden=0

! According to getden, build _VHA file name, referred as filvhain
! By default if positron=1
  stringfile='_VHA' ; stringvar='vha'
  call mkfilename(filnam,filvhain,0,idtset,0,jdtset_,&
&  ndtset,stringfile,stringvar,will_read)
  if(will_read==0)filvhain=trim(filnam_ds(3))//'_VHA'


! According to get1den, build _DEN file name, referred as fildens1in
! A default is available if get1den is 0
  stringfile='_DEN' ; stringvar='1den'
  call mkfilename(filnam,fildens1in,dtset%get1den,idtset,0,jdtset_,&
&  ndtset,stringfile,stringvar,will_read)
  if(will_read==0)fildens1in=trim(filnam_ds(3))//'_DEN'

  if (dtset%optdriver==4) then
!  According to getscr and irdscr, build _SCR file name, referred as filscr
!  A default is available if getscr is 0
   stringfile='_SCR' ; stringvar='scr'
   call mkfilename(filnam,filscr,dtset%getscr,idtset,dtset%irdscr,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)
   if(will_read==0)filscr=trim(filnam_ds(3))//'_SCR'

!  According to getsuscep and irdsuscep, build _SUS file name, referred as filchi0  
!  A default is available if getsuscep is 0
   stringfile='_SUS' ; stringvar='sus'
   call mkfilename(filnam,filchi0,dtset%getsuscep,idtset,dtset%irdsuscep,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)
   if(will_read==0)filchi0=TRIM(filnam_ds(3))//'_SUS'
  end if

  if ((dtset%optdriver==3).or.(dtset%optdriver==4)) then
!  According to getkss and irdkss, build _KSS file name, referred as filkss
!  A default is available if getkss is 0
   stringfile='_KSS' ; stringvar='kss'
   call mkfilename(filnam,filkss,dtset%getkss,idtset,dtset%irdkss,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)
   if(will_read==0)filkss=trim(filnam_ds(3))//'_KSS'

!  According to getqps and irdqps, build _QPS file name, referred as filqps
!  A default is available if getqps is 0
   stringfile='_QPS' ; stringvar='qps'
   call mkfilename(filnam,filqps,dtset%getqps,idtset,dtset%irdqps,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)
   if(will_read==0)filqps=trim(filnam_ds(3))//'_QPS'
  end if

  dtfil%ireadden      =ireadden
  dtfil%ireadwf       =ireadwf
  dtfil%filnam_ds(1:5)=filnam_ds(1:5)
  dtfil%filchi0       =filchi0    
  dtfil%filqps        =filqps
  dtfil%filscr        =filscr
  dtfil%fildensin     =fildensin
  dtfil%fildens1in    =fildens1in
  dtfil%filkss        =filkss
  dtfil%filstat       =filstat
  dtfil%filvhain      =filvhain
  dtfil%fnamewffk     =fnamewffk
  dtfil%fnamewffq     =fnamewffq
  dtfil%fnamewffddk   =fnamewffddk
  dtfil%fnamewff1     =fnamewff1

! ****************************************************************************
! Treat other get variables

! If multi dataset mode, and already the second dataset,
! treatment of other get variables.
  if( ndtset>1 .and. idtset>1 )then

   getocc=dtset%getocc
   getvel=dtset%getvel
   getxcart=dtset%getxcart
   getxred=dtset%getxred
   getcell=dtset%getcell
!  Should be in chkinp ; should also check that the number of atoms coincide.
   if(getxcart/=0 .and. getxred/=0) then
    write(message, '(a,a,a,a,a,a,i3,a,a,a,i3,a,i3,a,a,a)' )ch10,&
&    ' driver : ERROR -',ch10,&
&    '  The input variables getxcart and getxred cannot be',ch10,&
&    '  simultaneously non-zero, while for idtset=',idtset,',',ch10,&
&    '  they are getxcart=',getxcart,', and getxred=',getxred,'.',ch10,&
&    '  Action : correct getxcart or getxred in your input file.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

!  DEBUG
!  write(6,*)' ndtset,idtset,getxred=',ndtset,idtset,getxred
!  ENDDEBUG

   if(getocc>0 .or. (getocc<0 .and. idtset+getocc>0) )then
!   In case getocc is a negative number (so must add to idtset)
    if(getocc<0 .and. idtset+getocc>0) iget=idtset+getocc
    if(getocc>0)then
     do iget=1,idtset
      if( dtsets(iget)%jdtset==getocc )exit
     end do
     if(iget==idtset)then
!     The index of the dataset, from which the data ought to be taken,
!     does not correspond to a previous dataset.
      write(message, '(a,a,a,a,i3,a,a,i3,a,a,a,a)' )ch10,&
&      ' driver : ERROR -',ch10,&
&      '  The component number',idtset,' of the input variable getocc,',&
&      ' equal to',getocc,',',ch10,&
&      '  does not correspond to an existing index.',&
&      '  Action : correct getocc or jdtset in your input file.'
      call wrtout(6,message,'COLL')
      call leave_new('COLL')
     end if
    end if
    occ(:)=results_out(iget)%occ(1:dtset%mband*dtset%nkpt*dtset%nsppol)
    write(message, '(a,i3,a,a)' )&
&    ' driver : getocc/=0, take occ from output of dataset with index',&
&    dtsets(iget)%jdtset,'.',ch10
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
   end if

!  Getcell has to be done BEFORE getxcart
!  since acell and rprim will be used
   if(getcell>0 .or. (getcell<0 .and. idtset+getcell>0) )then
!   In case getocc is a negative number (so must add to idtset)
    if(getcell<0 .and. idtset+getcell>0) iget=idtset+getcell
    if(getcell>0)then
     do iget=1,idtset
      if( dtsets(iget)%jdtset==getcell )exit
     end do
     if(iget==idtset)then
!     The index of the dataset, from which the data ought to be taken,
!     does not correspond to a previous dataset.
      write(message, '(a,a,a,a,i3,a,a,i3,a,a,a,a)' )ch10,&
&      ' driver : ERROR -',ch10,&
&      '  The component number',idtset,' of the input variable getcell,',&
&      ' equal to',getcell,',',ch10,&
&      '  does not correspond to an existing index.',&
&      '  Action : correct getcell or jdtset in your input file.'
      call wrtout(6,message,'COLL')
      call leave_new('COLL')
     end if
    end if
    acell(:)=results_out(iget)%acell(:)
    rprim(:,:)=results_out(iget)%rprim(:,:)
    write(message, '(a,i3,a,a)' )&
&    ' driver : getcell/=0, take acell and rprim from output of dataset with index',&
&    dtsets(iget)%jdtset,'.',ch10
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')

!   Check that the new acell and rprim are consistent with the input dilatmx
    call mkrdim(acell,rprim,rprimd)
    call chkdilatmx(dtset%dilatmx,rprimd,dtset%rprimd_orig)

   end if

   if(getxred>0 .or. (getxred<0 .and. idtset+getxred>0) )then
!   In case getxred is a negative number (so must add to idtset)
    if(getxred<0 .and. idtset+getxred>0) iget=idtset+getxred
    if(getxred>0)then
     do iget=1,idtset
      if( dtsets(iget)%jdtset==getxred )exit
     end do
     if(iget==idtset)then
!     The index of the dataset, from which the data ought to be taken,
!     does not correspond to a previous dataset.
      write(message, '(a,a,a,a,i3,a,a,i3,a,a,a,a)' )ch10,&
&      ' driver : ERROR -',ch10,&
&      '  The component number',idtset,' of the input variable getxred,',&
&      ' equal to',getxred,',',ch10,&
&      '  does not correspond to an existing index.',&
&      '  Action : correct getxred or jdtset in your input file.'
      call wrtout(6,message,'COLL')
      call leave_new('COLL')
     end if
    end if
    xred(:,:)=results_out(iget)%xred(:,1:dtset%natom)
    write(message, '(a,i3,a,a)' )&
&    ' driver : getxred/=0, take xred from output of dataset with index',&
&    dtsets(iget)%jdtset,'.',ch10
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
   end if

   if(getxcart>0 .or. (getxcart<0 .and. idtset+getxcart>0) )then
!   In case getxcart is a negative number (so must add to idtset)
    if(getxcart<0 .and. idtset+getxcart>0) iget=idtset+getxcart
    if(getxcart>0)then
     do iget=1,idtset
      if( dtsets(iget)%jdtset==getxcart )exit
     end do
     if(iget==idtset)then
!     The index of the dataset, from which the data ought to be taken,
!     does not correspond to a previous dataset.
      write(message, '(a,a,a,a,i3,a,a,i3,a,a,a,a)' )ch10,&
&      ' driver : ERROR -',ch10,&
&      '  The component number',idtset,' of the input variable getxcart,',&
&      ' equal to',getxcart,',',ch10,&
&      '  does not correspond to an existing index.',&
&      '  Action : correct getxcart or jdtset in your input file.'
      call wrtout(6,message,'COLL')
      call leave_new('COLL')
     end if
    end if
!   Compute xcart of the previous dataset
    allocate(xcart(3,dtset%natom),xredget(3,dtset%natom))
    rprimdget(:,:)=results_out(iget)%rprimd(:,:)
    xredget (:,:)=results_out(iget)%xred(:,1:dtset%natom)
    call xredxcart(dtset%natom,1,rprimdget,xcart,xredget)
!   xcart from previous dataset is computed. Now, produce xred for the new dataset.
!   Now new acell and rprim ...
    call mkrdim(acell,rprim,rprimd)
    call xredxcart(dtset%natom,-1,rprimd,xcart,xred)
    deallocate(xcart,xredget)
    write(message, '(a,i3,a,a)' )&
&    ' driver : getxcart/=0, take xcart from output of dataset with index',&
&    dtsets(iget)%jdtset,'.',ch10
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
   end if

   if(getvel>0 .or. (getvel<0 .and. idtset+getvel>0) )then
!   In case getvel is a negative number (so must add to idtset)
    if(getvel<0 .and. idtset+getvel>0) iget=idtset+getvel
    if(getvel>0)then
     do iget=1,idtset
      if( dtsets(iget)%jdtset==getvel )exit
     end do
     if(iget==idtset)then
!     The index of the dataset, from which the data ought to be taken,
!     does not correspond to a previous dataset.
      write(message, '(a,a,a,a,i3,a,a,i3,a,a,a,a)' )ch10,&
&      ' driver : ERROR -',ch10,&
&      '  The component number',idtset,' of the input variable getvel,',&
&      ' equal to',getvel,',',ch10,&
&      '  does not correspond to an existing index.',&
&      '  Action : correct getvel or jdtset in your input file.'
      call wrtout(6,message,'COLL')
      call leave_new('COLL')
     end if
    end if
    vel(:,:)=results_out(iget)%vel(:,1:dtset%natom)
    write(message, '(a,i3,a,a)' )&
&    ' driver : getvel/=0, take vel from output of dataset with index',&
&    dtsets(iget)%jdtset,'.',ch10
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
   end if

  end if

! ****************************************************************************
! Treat the pseudopotentials : initialize the psps variable

  call status(jdtset_status,filstat,iexit,level,'init psps     ')

! Note that mpsang is the max of 1+lmax, with minimal value 1 (even for local psps, at present)
! mpsang=max(maxval(pspheads(1:npsp)%lmax)+1,1) ! might not work with HP compiler
! n1xccc=maxval(pspheads(1:npsp)%xccc)
  mpsang=1
  n1xccc=pspheads(1)%xccc
  do ii=1,npsp
   mpsang=max(pspheads(ii)%lmax+1,mpsang)
   n1xccc=max(pspheads(ii)%xccc,n1xccc)
  end do

! Determine the maximum number of projectors, for the set of pseudo atom
  call getdim_nloc(lmnmax,lmnmaxso,lnmax,lnmaxso,dtset%mixalch,npsp,dtset%npspalch,&
&  dtset%ntypat,dtset%ntypalch,pspheads)

  psps%mpsang   = mpsang
  psps%mtypalch = mtypalch
  psps%npsp     = npsp
  psps%npspalch = dtset%npspalch
  psps%ntypat   = dtset%ntypat
  psps%ntypalch = dtset%ntypalch
  psps%ntyppure = dtset%ntyppure
  psps%n1xccc   = n1xccc

! Set the flag for reciprocal space or real space calculations
  psps%vlspl_recipSpace = (dtset%icoulomb /= 1)
! changed by RShaltaf 
  psps%positron = dtset%positron
  psps%usepaw  =usepaw
  psps%useylm  =dtset%useylm

  allocate(psps%algalch(psps%ntypalch))
  allocate(psps%mixalch(psps%npspalch,psps%ntypalch))
  psps%algalch(1:psps%ntypalch)=dtset%algalch(1:psps%ntypalch)
  psps%mixalch(1:psps%npspalch,1:psps%ntypalch)=dtset%mixalch(1:psps%npspalch,1:psps%ntypalch)

! Set mpspso and psps%pspso
! Warning : mpspso might be different for each dataset.
  psps%mpspso=1
  do ipsp=1,dtset%npsp
   if(dtset%nspinor==1)then
    psps%pspso(ipsp)=0
   else
    if(dtset%so_psp(ipsp)/=1)then
     psps%pspso(ipsp)=dtset%so_psp(ipsp)
    else
     psps%pspso(ipsp)=pspheads(ipsp)%pspso
    end if
    if(psps%pspso(ipsp)/=0)psps%mpspso=2
   end if
!  Ideally the following line should not exist, but at present, the space has to be booked
   if(pspheads(ipsp)%pspso/=0)psps%mpspso=2
  end do

! Set mpssoang, lmnmax, lnmax
  if(psps%mpspso==1)then
   psps%mpssoang=psps%mpsang
   psps%lmnmax  =lmnmax
   psps%lnmax   =lnmax
  else
   psps%mpssoang=2*psps%mpsang-1
   psps%lmnmax=lmnmaxso
   psps%lnmax=lnmaxso
  end if
  if (psps%useylm==0) then
   psps%lmnmax=psps%lnmax
  end if

! Set dimekb
  if (psps%usepaw==0) then
   psps%dimekb=psps%lnmax
  else
   psps%dimekb=psps%lmnmax*(psps%lmnmax+1)/2
  end if

! The following arrays are often not deallocated before the end of the dtset loop
! and might keep their content from one dataset to the other,
! if the conditions are fulfilled
  if(dimekb_old/=psps%dimekb .or. ntypat_old/=dtset%ntypat .or. &
&  usepaw_old/=psps%usepaw) then
   if(idtset/=1)deallocate(psps%ekb)
   allocate(psps%ekb(psps%dimekb,dtset%ntypat*(1-psps%usepaw)))
   dimekb_old=psps%dimekb
  end if
  if(lmnmax_old/=psps%lmnmax .or. ntypat_old/=dtset%ntypat)then
   if(idtset/=1)deallocate(psps%indlmn)
   allocate(psps%indlmn(6,psps%lmnmax,dtset%ntypat))
   lmnmax_old=psps%lmnmax
  end if
  if(mqgridff_old/=psps%mqgrid_ff .or. lnmax_old/=psps%lnmax .or. ntypat_old/=dtset%ntypat)then
   if(idtset/=1)deallocate(psps%ffspl,psps%qgrid_ff)
   allocate(psps%ffspl(psps%mqgrid_ff,2,psps%lnmax,dtset%ntypat))
   allocate(psps%qgrid_ff(psps%mqgrid_ff))
   mqgridff_old=psps%mqgrid_ff
   lnmax_old=psps%lnmax
  end if
  if(mqgridvl_old/=psps%mqgrid_vl .or. ntypat_old/=dtset%ntypat)then
   if(idtset/=1)deallocate(psps%vlspl,psps%qgrid_vl)
   if (idtset/=1 .and. .not.psps%vlspl_recipSpace) then
    deallocate(psps%dvlspl)
   end if
   allocate(psps%vlspl(psps%mqgrid_vl,2,dtset%ntypat),psps%qgrid_vl(psps%mqgrid_vl))
   if (.not.psps%vlspl_recipSpace) then
    allocate(psps%dvlspl(psps%mqgrid_vl,2,dtset%ntypat))
   end if
   mqgridvl_old=psps%mqgrid_vl
  end if
  if(ntypat_old/=dtset%ntypat.or. usepaw_old/=psps%usepaw)then
   if(idtset/=1)deallocate(psps%xccc1d)
   allocate(psps%xccc1d(n1xccc*(1-psps%usepaw),6,dtset%ntypat))
   usepaw_old=psps%usepaw
  end if
  if(ntypat_old/=dtset%ntypat)then
   if(idtset/=1)deallocate(psps%xcccrc,psps%ziontypat,psps%znucltypat)
   allocate(psps%xcccrc(dtset%ntypat))
   allocate(psps%znucltypat(dtset%ntypat))
   allocate(psps%ziontypat(dtset%ntypat))
   ntypat_old=dtset%ntypat
!  The correct dimension of pawrad/tab is ntypat.
!  In case of paw, no alchemical psp is allowed, so npsp=ntypat
!  However, in case of alchemical psps, pawrad/tab(ipsp) is invoked in
!  pspini. So, in order to avoid any problem, declare pawrad/tab
!  at paw_size=max(ntypat,npsp).
   paw_size=0;if (psps%usepaw==1) paw_size=max(dtset%ntypat,npsp)
  end if
  psps%ziontypat(:)=dtset%ziontypat(:)

! ****************************************************************************
! PAW allocations.

  call status(jdtset_status,filstat,iexit,level,'PAW allocs    ')

  if (paw_size/=paw_size_old) then
   if(idtset/=1) then
    call pawalloc(dtset,idtset,mpsang,psps%mqgrid_vl,npsp,2,paw_size,paw_size_old,&
&    pawang,pawrad,pawtab,pspheads)
    deallocate(pawrad,pawtab)
   end if
   allocate(pawrad(paw_size),pawtab(paw_size))
  end if
  call pawalloc(dtset,idtset,mpsang,psps%mqgrid_vl,npsp,1,paw_size,paw_size_old,&
&  pawang,pawrad,pawtab,pspheads)
  paw_size_old=paw_size

! ****************************************************************************
! ETSF I/O

  call status(jdtset_status,filstat,iexit,level,'ETSF IO       ')

! Fill in abidims structure
  abidims%mband    = dtset%mband
  abidims%mproj    = 1                 ! FIXME
  abidims%mpsang   = mpsang
  abidims%mpw      = dtset%mpw
  abidims%natom    = dtset%natom
  abidims%natsph   = dtset%natsph
  abidims%nberry   = dtset%nberry
  abidims%nconeq   = dtset%nconeq
  abidims%nfft     = dtset%nfft
  abidims%nfreqsus = dtset%nfreqsus
  abidims%ngrid1   = dtset%ngfft(1)
  abidims%ngrid2   = dtset%ngfft(2)
  abidims%ngrid3   = dtset%ngfft(3)
  abidims%nkpt     = dtset%nkpt
  abidims%nkptgw   = dtset%nkptgw
  abidims%npsp     = npsp
  abidims%npspalch = dtset%npspalch
  abidims%nqptdm   = dtset%nqptdm
  abidims%nshiftk  = dtset%nshiftk
  abidims%nspden   = dtset%nspden
  abidims%nspinor  = dtset%nspinor
  abidims%nsppol   = dtset%nsppol
  abidims%nsym     = dtset%nsym
  abidims%ntypat   = dtset%ntypat
  abidims%ntypalch = dtset%ntypalch
  abidims%wfs_dim1 = -1
  abidims%wfs_dim2 = -1
  abidims%npw_tiny = -1


! ****************************************************************************
! At this stage, all the data needed for the treatment of one dataset
! have been transferred from multi-dataset arrays.

  iexit=0

! Smaller integer arrays :
  allocate(npwtot(dtset%nkpt))

  allocate(results_gs%fcart(3,dtset%natom),results_gs%fred(3,dtset%natom))
! Also allocate some other components of results_gs,
! even if they are not used in the present routine.
  allocate(results_gs%gresid(3,dtset%natom),results_gs%grewtn(3,dtset%natom))
  allocate(results_gs%grxc(3,dtset%natom),results_gs%synlgr(3,dtset%natom))

! Initialize these to zero (needed when hasty exit)
  etotal=zero
  strten(:)=zero
  call energies_init(results_gs%energies)
  results_gs%etotal      = zero
  results_gs%fcart(:,:)  = zero
  results_gs%fred(:,:)   = zero
  results_gs%gresid(:,:) = zero
  results_gs%grewtn(:,:) = zero
  results_gs%grxc(:,:)   = zero
  results_gs%strten(:)   = zero
  results_gs%synlgr(:,:) = zero

  select case(dtset%optdriver)

   case(0)

    call status(jdtset_status,filstat,iexit,level,'call gstate   ')

    if(dtset%parareel==0)then

     mpi_enreg%parareel=0
     call gstate(acell,codvsn,cpui,dtfil,dtset,iexit,&
&     mpi_enreg,&
&     npwtot,dtset%nspinor,&
&     occ,pawang,pawrad,pawtab,psps,results_gs,rprim,vel,walli,xred)

    else

     mpi_enreg%parareel=1
     mpi_enreg%paral_level=1
     call pstate(acell,codvsn,cpui,dtfil,dtset,iexit,&
&     dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,dtset%nfft,dtset%nkpt,&
&     npwtot,dtset%nspden,dtset%nspinor,dtset%nsppol,dtset%nsym,&
&     occ,pawrad,pawtab,psps,results_gs,rprim,vel,walli,xred)

    end if

    etotal=results_gs%etotal
    strten(1:6)=results_gs%strten(1:6)

    call status(jdtset_status,filstat,iexit,level,'after gstate  ')
   case(1)

    call status(jdtset_status,filstat,iexit,level,'call respfn   ')

    mpi_enreg%paral_level=2
    call respfn(codvsn,cpui,dtfil,dtset,etotal,iexit,mkmems,mpi_enreg,&
&    npwtot,dtset%nspinor,occ,pawang,pawrad,pawtab,psps,walli,xred)

    call status(jdtset_status,filstat,iexit,level,'after respfn  ')

   case(2)

    call status(jdtset_status,filstat,iexit,level,'call suscep   ')

    mpi_enreg%paral_level=2
    call suscep(dtfil,dtset,iexit,&
&    dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,dtset%nfft,dtset%nkpt,&
&    dtset%nspden,dtset%nspinor,dtset%nsppol,dtset%nsym,occ,xred,abidims,dtset%ngfft)

    call status(jdtset_status,filstat,iexit,level,'after suscep  ')

   case(3)

    call status(jdtset_status,filstat,iexit,level,'call screening')

    mpi_enreg%paral_level=2
    call screening(acell,codvsn,dtfil,dtset,iexit,mpi_enreg,pawang,pawrad,pawtab,psps,rprim,xred)

    call status(jdtset_status,filstat,iexit,level,'after screenin')

   case(4)

    call status(jdtset_status,filstat,iexit,level,'call sigma    ')

    mpi_enreg%paral_level=2
    call sigma(acell,codvsn,dtfil,dtset,iexit,mpi_enreg,pawang,pawrad,pawtab,psps,rprim,xred)

    call status(jdtset_status,filstat,iexit,level,'after sigma   ')

   case(5)

    call status(jdtset_status,filstat,iexit,level,'call nonlinear   ')

    mpi_enreg%paral_level=2
    call nonlinear(codvsn,dtfil,dtset,etotal,iexit,&
&    dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,dtset%nfft,dtset%nkpt,npwtot,dtset%nspden,&
&    dtset%nspinor,dtset%nsppol,dtset%nsym,occ,pawrad,pawtab,psps,xred)

    call status(jdtset_status,filstat,iexit,level,'after nonlinear  ')

   case(6)

    call status(jdtset_status,filstat,iexit,level,'call wannier   ')

    mpi_enreg%paral_level=2
    call wannier(dtfil,dtset,iexit,dtset%mband,mpi_enreg,dtset%nkpt,dtset%nsppol)

    call status(jdtset_status,filstat,iexit,level,'after nonlinear  ')

   case(7)

!   According to getkss, build _KSS file name, referred as filkss
!   A default is available if getkss is 0
    stringfile='_KSS' ; stringvar='kss'
    call mkfilename(filnam,filkss,dtset%getkss,idtset,dtset%irdkss,jdtset_,&
&    ndtset,stringfile,stringvar,will_read)
    if(will_read==0)filkss=trim(filnam_ds(3))//'_KSS'

    dtfil%filkss        =filkss

    call rdm(acell,dtfil,dtset,pawtab,mpi_enreg,rprim)

    case default

!   Error the choice is either 0 -> gstate, 1 -> respfn, 2 -> suscep,
!   3 -> screening, 4 -> sigma,  5 -> nonlinear, 6 -> wannier
    write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&    ' driver : ERROR -',ch10,&
&    '  The variable optdriver must be between 0 and 6, but was ',&
&    dtset%optdriver,ch10,&
&    '  This is not allowed.  ',ch10,&
&    '  Action : modify optdriver in the input file.'

    call wrtout(06,  message,'COLL')
    call leave_new('COLL')

  end select

! ****************************************************************************

! DEBUG
! write(6,*)' driver : after respfn'
! ENDDEBUG

! Transfer of multi dataset outputs from temporaries :
! acell, xred, occ rprim, and vel might be modified from their
! input values
! etotal, fcart, fred, and strten have been computed
! npwtot was already computed before, but is stored only now

  results_out(idtset)%acell(:)          =acell(:)
  results_out(idtset)%etotal            =etotal
  results_out(idtset)%rprim(:,:)        =rprim(:,:)
  call mkrdim(acell,rprim,rprimd)
  results_out(idtset)%rprimd(:,:)       =rprimd(:,:)
  results_out(idtset)%strten(:)         =strten(:)
  results_out(idtset)%fcart(1:3,1:dtset%natom)=results_gs%fcart(:,:)
  results_out(idtset)%fred(1:3,1:dtset%natom) =results_gs%fred(:,:)
  results_out(idtset)%npwtot(1:dtset%nkpt)       =npwtot(1:dtset%nkpt)
  results_out(idtset)%occ(1:dtset%mband*dtset%nkpt*dtset%nsppol)=occ(:)
  results_out(idtset)%vel(:,1:dtset%natom)    =vel(:,:)
  results_out(idtset)%xred(:,1:dtset%natom)   =xred(:,:)

! DEBUG
! write(6,*)' driver : 0'
! ENDDEBUG

  call dtsetFree(dtset)

! DEBUG
! write(6,*)' driver : 1'
! ENDDEBUG

  deallocate(psps%algalch,psps%mixalch)

  deallocate(occ)
  deallocate(vel,xred)
  deallocate(npwtot)
  deallocate(results_gs%fcart,results_gs%fred)
  deallocate(results_gs%gresid,results_gs%grewtn)
  deallocate(results_gs%grxc,results_gs%synlgr)

! DEBUG
! write(6,*)' driver : 2, iexit=',iexit
! ENDDEBUG

  if(iexit/=0)exit

! Check whether exiting was required by the user.
! If found then beat a hasty exit from time steps
  openexit=1 ; if(dtset%chkexit==0) openexit=0
! update for the end of parareel case : chkexi will crash otherwise
  if (mpi_enreg%parareel == 1) then
   mpi_enreg%parareel=0
  end if

  call chkexi(zero,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)

  if (iexit/=0)exit

! End do loop on idtset (allocate statements are present -
! an exit statement is present)
 end do

!*********************************************************************

!DEBUG
!write(6,*)' driver : 3'
!ENDDEBUG

 deallocate(psps%qgrid_ff)
 deallocate(psps%qgrid_vl)
 deallocate(psps%xcccrc)
 deallocate(psps%xccc1d)
 deallocate(psps%vlspl)
 if (.not.psps%vlspl_recipSpace) then
  deallocate(psps%dvlspl)
 end if
 deallocate(psps%ekb)
 deallocate(psps%indlmn)
 deallocate(psps%filpsp)
 deallocate(psps%pspcod)
 deallocate(psps%pspdat)
 deallocate(psps%pspso)
 deallocate(psps%pspxc)
 deallocate(psps%title)
 deallocate(psps%znuclpsp)
 deallocate(psps%znucltypat)
 deallocate(psps%zionpsp)
 deallocate(psps%ziontypat)
 deallocate(psps%ffspl)
 call psp2params_free(psps%gth_params)

!PAW deallocation
 call pawalloc(dtset,idtset,mpsang,psps%mqgrid_vl,npsp,3,paw_size,paw_size_old,&
& pawang,pawrad,pawtab,pspheads)
 deallocate(pawrad,pawtab)

 deallocate(jdtset_)

!DEBUG
!write(6,*)' driver : before exit '
!write(6, '(a,9f5.2)' ) ' driver : before exit  , rprim ',rprim (:,:)
!call flush(6)
!stop
!ENDDEBUG

!DEBUG
!write(6,*)' driver : return for test memory leak '
!return
!ENDDEBUG


 call status(0,filstat,iexit,level,'exit          ')
 call timab(100,2,tsec)

end subroutine driver
!!***
