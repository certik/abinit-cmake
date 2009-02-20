!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkinp
!! NAME
!! chkinp
!!
!! FUNCTION
!! Check consistency of input data against itself.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, MKV, DRH, MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  iout=unit number for output file
!!  mpi_enreg=informations about MPI parallelization
!!  ndtset=number of datasets
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!      one data set.
!!  npsp=number of pseudopotentials
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      chkdpr,chkgrp,chkint,chkorthsy,leave_new,metric,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine chkinp(dtsets,iout,mpi_enreg,ndtset,ndtset_alloc,npsp,pspheads)

 use defs_basis
 use defs_datatypes
#if defined HAVE_ETSF_IO
 use etsf_io
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_13iovars, except_this_one => chkinp
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,ndtset,ndtset_alloc,npsp
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables-------------------------------
!scalars
 logical :: twvl
 integer :: accesswff,awtr,bantot,berryopt,ceksph,dmatpuopt,dmatudiag,fftalg,frzfermi,fftgw
 integer :: gwpara,gwmem,gwcomp,ia,iatom,ib,iband,idtset,ierr,ii,ikpt,ilang,intxc,ionmov
 integer :: iprcch,iprcel,iprctfvw,ipsp,iscf,isppol,isym,itypat,ixc,jdtset
 integer :: jellslab,jj,kptopt,kssform,localrdwf,maxiatsph,mband,mffmem
 integer :: miniatsph,mk1mem,mkmem,mkqmem,mproj,mu,natom,natsph,nbandkss
 integer :: nbdblock,nberry,nfft,nfftdg,nkpt,nkpt_me,nloalg,npband,npfft,npkpt,npspalch,npwkss
 integer :: nqpt,nshiftk,nspden,nspinor,nsppol,nstep,nsym,ntypalch,ntypat,occopt
 integer :: optcell,optdriver,optforces,optstress,pawlcutd,pawlmix,pawmixdg
 integer :: pawnhatxc,pawnzlm,pawoptmix,pawprtdos,pawprtvol,pawspnorb,pawstgylm,pawusecp,pawxcdev,prepanl,prepgkk
 integer :: prtden,prtdensph,prtdos,prtdosm,prteig,prtfsurf,prtgeo,prtnabla,prtstm,prtvha
 integer :: prtvhxc,prtvol,prtvxc,prtwant,prtwf,response,rfstrs,spmeth,symchi
 integer :: symmorphi,symsigma,usedmatpu,useexexch,usepaw,usepawu,usewvl,useylm,wfoptalg
 real(dp) :: boxcutmin,charge,diecut,diemac,dosdeltae,ecut,ecuteps,ecutsigx,ecutsm,ecutwfn,exchmix,fixmom,kptnrm
 real(dp) :: norm,pawecutdg,residual,slabwsrad,slabzbeg,slabzend
 real(dp) :: stmbias,sumalch,sumocc,toldfe,toldff,tolrff,tolwfr,tsmear,ucvol
 real(dp) :: wtksum,wvl_hgrid,zatom
 character(len=500) :: message
!arrays
 integer :: bdberry(4),cond_values(3),ngfft(3),ngfftdg(3),nprojmax(0:3)
 integer :: rfatpol(2)
 integer :: kptrlatt(3,3)
 integer,allocatable :: iatsph(:),istwfk(:),lexexch(:),lpawu(:),nband(:)
 integer,allocatable :: so_psp(:),symafm(:),symrel(:,:,:),typat(:)
 real(dp) :: gmet(3,3),gprimd(3,3),prods(3,3),rmet(3,3),rmet_sym(3,3)
 real(dp) :: rprimd(3,3),rprimd_sym(3,3)
 real(dp),allocatable :: frac(:,:),jpawu(:),kpt(:,:),mixalch(:,:),occ(:)
 real(dp),allocatable :: ratsph(:),tnons(:,:),upawu(:),wtk(:),xred(:,:),znucl(:)
 real(dp),allocatable :: shiftk(:,:)
 character(len=9) :: cond_string(3)

! *************************************************************************

!DEBUG
!write(6,*)' chkinp : enter '
!stop
!ENDDEBUG

!Print machine precision (other machine parameters are computed
!in the dlamch function, see Lapack library)
 write(message,'(a,a,1p,e24.16)' ) ch10,&
& ' chkinp: machine precision is ',epsilon(0.0_dp)
 call wrtout(06,  message,'COLL')

!Some initialisations
 ierr=0
 cond_string(1:3)=' '
 cond_values(1:3)=(/0,0,0/)

!Do loop on idtset (allocate statements are present)
 do idtset=1,ndtset_alloc

! Copy input dataset values
  dmatpuopt=dtsets(idtset)%dmatpuopt
  dmatudiag=dtsets(idtset)%dmatudiag
  useexexch=dtsets(idtset)%useexexch
  fftalg   =dtsets(idtset)%ngfft(7)
  fftgw    =dtsets(idtset)%fftgw
  iprctfvw =dtsets(idtset)%iprctfvw
  kptnrm   =dtsets(idtset)%kptnrm
  mband    =dtsets(idtset)%mband
  mffmem   =dtsets(idtset)%mffmem
  mkmem    =dtsets(idtset)%mkmem
  mkqmem   =dtsets(idtset)%mkqmem
  mk1mem   =dtsets(idtset)%mk1mem
  natsph   =dtsets(idtset)%natsph
  natom    =dtsets(idtset)%natom
  nkpt     =dtsets(idtset)%nkpt
  nloalg   =dtsets(idtset)%nloalg(5)
  npband   =dtsets(idtset)%npband
  npfft    =dtsets(idtset)%npfft
  npkpt    =dtsets(idtset)%npkpt
  npspalch =dtsets(idtset)%npspalch
  nshiftk  =dtsets(idtset)%nshiftk
  nspden   =dtsets(idtset)%nspden
  nspinor  =dtsets(idtset)%nspinor
  nsppol   =dtsets(idtset)%nsppol
  nstep    =dtsets(idtset)%nstep
  nsym     =dtsets(idtset)%nsym
  ntypat   =dtsets(idtset)%ntypat
  ntypalch =dtsets(idtset)%ntypalch
  occopt   =dtsets(idtset)%occopt
  optdriver=dtsets(idtset)%optdriver
  optforces=dtsets(idtset)%optforces
  optstress=dtsets(idtset)%optstress
  prepanl  =dtsets(idtset)%prepanl
  prepgkk  =dtsets(idtset)%prepgkk
  prtden   =dtsets(idtset)%prtden
  prtdensph=dtsets(idtset)%prtdensph
  prteig   =dtsets(idtset)%prteig
  prtfsurf =dtsets(idtset)%prtfsurf
  prtnabla =dtsets(idtset)%prtnabla
  prtstm   =dtsets(idtset)%prtstm
  prtvha   =dtsets(idtset)%prtvha
  prtvhxc  =dtsets(idtset)%prtvhxc
  prtvol   =dtsets(idtset)%prtvol
  prtvxc   =dtsets(idtset)%prtvxc
  prtwant  =dtsets(idtset)%prtwant
  prtwf    =dtsets(idtset)%prtwf
  symmorphi=dtsets(idtset)%symmorphi
  symchi   =dtsets(idtset)%symchi
  symsigma =dtsets(idtset)%symsigma
  usedmatpu=dtsets(idtset)%usedmatpu
  usepaw   =dtsets(idtset)%usepaw
  usepawu  =dtsets(idtset)%usepawu
  usewvl   =dtsets(idtset)%usewvl
  rfatpol(1)=dtsets(idtset)%rfatpol(1)
  rfatpol(2)=dtsets(idtset)%rfatpol(2)
  rprimd(:,:)=dtsets(idtset)%rprimd_orig(:,:)
  ngfft(:) =dtsets(idtset)%ngfft(1:3)
  ngfftdg(:)=dtsets(idtset)%ngfftdg(1:3)

! Allocate arrays
  allocate(kpt(3,nkpt),iatsph(natsph),istwfk(nkpt))
  allocate(mixalch(npspalch,ntypalch),nband(nkpt*nsppol))
  allocate(occ(mband*nkpt*nsppol),shiftk(3,nshiftk),so_psp(npsp))
  allocate(ratsph(ntypat))
  allocate(symafm(nsym),symrel(3,3,nsym),tnons(3,nsym),typat(natom))
  allocate(wtk(nkpt),xred(3,natom),znucl(npsp) )
  allocate(upawu(ntypat),jpawu(ntypat),lpawu(ntypat),lexexch(ntypat))

  iatsph(:)  =dtsets(idtset)%iatsph(1:natsph)
  istwfk(:)  =dtsets(idtset)%istwfk(1:nkpt)
  jpawu  (:) =dtsets(idtset)%jpawu   (1:ntypat)
  kpt   (:,:)=dtsets(idtset)%kpt(1:3,1:nkpt)
  kptrlatt(:,:)=dtsets(idtset)%kptrlatt(1:3,1:3)
  lpawu (:)  =dtsets(idtset)%lpawu   (1:ntypat)
  lexexch (:)  =dtsets(idtset)%lexexch   (1:ntypat)
  nband (:)  =dtsets(idtset)%nband(1:nkpt*nsppol)
  occ   (:)  =dtsets(idtset)%occ_orig(1:mband*nkpt*nsppol)
  ratsph(:)  =dtsets(idtset)%ratsph(1:ntypat)
  so_psp(:)  =dtsets(idtset)%so_psp(1:npsp)
  symafm(:)  =dtsets(idtset)%symafm(1:nsym)
  symrel(:,:,:)=dtsets(idtset)%symrel(1:3,1:3,1:nsym)
  tnons (:,:)=dtsets(idtset)%tnons(1:3,1:nsym)
  typat (:)  =dtsets(idtset)%typat(1:natom)
  upawu (:)  =dtsets(idtset)%upawu  (1:ntypat)
  wtk   (:)  =dtsets(idtset)%wtk   (1:nkpt)
  xred  (:,:)=dtsets(idtset)%xred_orig(:,1:natom)
  znucl (:)  =dtsets(idtset)%znucl(1:npsp)
  if(npspalch>0.and.ntypalch>0)then
   mixalch(:,:)=dtsets(idtset)%mixalch(1:npspalch,1:ntypalch)
  end if

  jdtset=dtsets(idtset)%jdtset
  if(ndtset==0)jdtset=0

  if(jdtset/=0)then
   write(message, '(a,a,a,i2,a)' ) ch10,&
&   ' chkinp: Checking input parameters for consistency,',&
&   ' jdtset=',jdtset,'.'
  else
   write(message, '(a,a)' ) ch10,&
&   ' chkinp: Checking input parameters for consistency.'
  end if
  call wrtout(iout,message,'COLL')
  call wrtout(06,  message,'COLL')

  accesswff =dtsets(idtset)%accesswff
  awtr      =dtsets(idtset)%awtr
  boxcutmin =dtsets(idtset)%boxcutmin
  ceksph    =dtsets(idtset)%ceksph
  dosdeltae =dtsets(idtset)%dosdeltae
  frzfermi  =dtsets(idtset)%frzfermi
  gwcomp    =dtsets(idtset)%gwcomp
  gwpara    =dtsets(idtset)%gwpara
  gwmem     =dtsets(idtset)%gwmem
  ionmov    =dtsets(idtset)%ionmov
  intxc     =dtsets(idtset)%intxc
  iprcel    =dtsets(idtset)%iprcel
  iprcch    =dtsets(idtset)%iprcch
  iprctfvw  =dtsets(idtset)%iprctfvw
  iscf      =dtsets(idtset)%iscf
  ixc       =dtsets(idtset)%ixc
  nqpt      =dtsets(idtset)%nqpt
  optcell   =dtsets(idtset)%optcell
  kptopt    =dtsets(idtset)%kptopt
  localrdwf =dtsets(idtset)%localrdwf
  nberry    =dtsets(idtset)%nberry
  bdberry(1:4)=dtsets(idtset)%bdberry(1:4)
  nbandkss  =dtsets(idtset)%nbandkss
  npwkss    =dtsets(idtset)%npwkss
  berryopt  =dtsets(idtset)%berryopt
  wfoptalg  =dtsets(idtset)%wfoptalg
  nbdblock  =dtsets(idtset)%nbdblock
  kssform   =dtsets(idtset)%kssform
  useylm    =dtsets(idtset)%useylm

  charge   =dtsets(idtset)%charge
  ecut     =dtsets(idtset)%ecut
  tsmear   =dtsets(idtset)%tsmear
  ecuteps  =dtsets(idtset)%ecuteps
  ecutsigx =dtsets(idtset)%ecutsigx
  ecutsm   =dtsets(idtset)%ecutsm
  ecutwfn  =dtsets(idtset)%ecutwfn
  exchmix  =dtsets(idtset)%exchmix
  fixmom   =dtsets(idtset)%fixmom

  diecut   =dtsets(idtset)%diecut
  diemac   =dtsets(idtset)%diemac

  jellslab =dtsets(idtset)%jellslab
  slabwsrad=dtsets(idtset)%slabwsrad
  slabzbeg =dtsets(idtset)%slabzbeg
  slabzend =dtsets(idtset)%slabzend

  pawecutdg=dtsets(idtset)%pawecutdg
  pawlcutd =dtsets(idtset)%pawlcutd
  pawlmix  =dtsets(idtset)%pawlmix
  pawmixdg =dtsets(idtset)%pawmixdg
  pawnhatxc=dtsets(idtset)%pawnhatxc
  pawnzlm  =dtsets(idtset)%pawnzlm
  pawoptmix=dtsets(idtset)%pawoptmix
  pawprtdos=dtsets(idtset)%pawprtdos
  pawprtvol=dtsets(idtset)%pawprtvol
  pawspnorb=dtsets(idtset)%pawspnorb
  pawstgylm=dtsets(idtset)%pawstgylm
  pawusecp =dtsets(idtset)%pawusecp
  pawxcdev =dtsets(idtset)%pawxcdev

  prtdos   =dtsets(idtset)%prtdos
  prtdosm  =dtsets(idtset)%prtdosm
  prtgeo   =dtsets(idtset)%prtgeo
  prtvxc   =dtsets(idtset)%prtvxc

  rfstrs   =dtsets(idtset)%rfstrs

  spmeth   =dtsets(idtset)%spmeth
  stmbias  =dtsets(idtset)%stmbias
  toldfe   =dtsets(idtset)%toldfe
  toldff   =dtsets(idtset)%toldff
  tolrff   =dtsets(idtset)%tolrff
  tolwfr   =dtsets(idtset)%tolwfr

  response=0
  if(dtsets(idtset)%rfelfd/=0 .or. dtsets(idtset)%rfmgfd/=0 .or. &
&  dtsets(idtset)%rfphon/=0 .or. dtsets(idtset)%rfstrs/=0  )response=1

  call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

! check if nsym=1 in phonon calculation in finite electric field
  if( (response==1) .and. (dtsets(idtset)%berryopt==4) ) then
   if (dtsets(idtset)%nsym/=1) then
    write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '  For phonon calculation in finite electric field',ch10,&
&    '  nsym > 1 is not allowed currently.',ch10,&
&    ' Action : modify value of nsym to be 1 in input file.'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    ierr=ierr+1
   end if
  end if

! ** Here begins the checking section **************************************
! Check the values of variables, using alphabetical order

! accesswff
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'accesswff',accesswff,4,(/0,1,2,3/),0,0,iout)
  if (accesswff == 3) then
   do ikpt=1,nkpt
    if(istwfk(ikpt)/=1)then
!    ETSF_IO, current limitation to istwfk 1
     write(message,'(8a,I0,a)' ) ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '  When accesswff==3, all the components of istwfk must be 1.',ch10,&
&     '  Not yet programmed for time-reversal symmetry.',ch10,&
&     ' Action : modify value of istwfk to be 1 for all k points in input file.'
     call wrtout(iout,message,'COLL')
     call wrtout(6,message,'COLL')

!    The following line was wrong : istwfk is used previously, to dimension the
!    number of plane waves, so it cannot be modified here ...
!    (And moreover, no input variable should be modified in chkinp, that is
!    supposed only to do the checking).
!    dtsets(idtset)%istwfk(:) = 1

     ierr=ierr+1
     exit
    end if
   end do
  end if

  if(mpi_enreg%paral_compil_mpio==0)then
   cond_string(1)='MPIO flag in makefile_macros' ; cond_values(1)=0
!  Make sure that accesswff is 0 or 2 or 1
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'accesswff',accesswff,4,(/0,1,2,3/),0,0,iout)
  end if

! amu
! Check that atomic masses are > 0 if ionmov = 1
  if (ionmov==1) then
   do itypat=1,ntypat
    if (dtsets(idtset)%amu(itypat)<=0.0_dp) then
     write(message, '(a,a,a,a,i3,a,1p,e12.4,a,a,a,a,a)' ) ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '  For ionmov =1, amu(',itypat,' ) was input as',&
&     dtsets(idtset)%amu(itypat),' .',ch10,&
&     ' Input value must be > 0 for molecular dynamics.',&
&     ch10,' Action : modify value of amu in input file.'
     call wrtout(iout,message,'COLL')
     call wrtout(06,  message,'COLL')
     ierr=ierr+1
    end if
   end do
  end if

! bdberry
  if(berryopt>0 .and. berryopt/=4 .and. nberry>0)then
   do ii=1,2*nsppol
    if (bdberry(ii)<1) then
     write(message, '(a,a,a,a,i3,a,i5,a,a,a,a,a)' ) ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '  For berryopt>0, bdberry(',ii,' ) was input as',&
&     bdberry(ii),' .',ch10,&
&     ' Input value must be > 0 when berryopt>0.',&
&     ch10,' Action : modify value of bdberry in input file.'
     call wrtout(iout,message,'COLL')
     call wrtout(06,  message,'COLL')
     ierr=ierr+1
    end if
   end do
   if(bdberry(2)<bdberry(1))then
    write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '  For berryopt>0, bdberry(1) is larger than bdberry(2).',ch10,&
&    '  This is not allowed.',ch10,&
&    ' Action : modify value of bdberry in input file.'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    ierr=ierr+1
   end if
   if(nsppol==2 .and. bdberry(4)<bdberry(3))then
    write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '  For berryopt>0 and nsppol=2, bdberry(3) is larger than bdberry(4).',ch10,&
&    '  This is not allowed.',ch10,&
&    ' Action : modify value of bdberry in input file.'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    ierr=ierr+1
   end if
!  Make sure all nband(nkpt) are >= bdberry
   do isppol=1,nsppol
    do ikpt=1,nkpt
     if (nband(ikpt+(isppol-1)*nkpt)<=bdberry(2*isppol)) then
      cond_string(1)='ikpt'
      cond_values(1)=ikpt
      cond_string(2)='isppol'
      cond_values(2)=isppol
      cond_string(3)='nband'
      cond_values(3)=nband(ikpt+(isppol-1)*nkpt)
      call chkint(0,3,cond_string,cond_values,ierr,&
&      'bdberry',bdberry(2*isppol),1,(/nband(ikpt+(isppol-1)*nkpt)/),&
&      -1,nband(ikpt+(isppol-1)*nkpt),iout)
      if(ierr==1)exit
     end if
    end do
   end do
  end if

! HERE, should use chkint ...
! berryopt
  if ((berryopt < -3).or.(berryopt > 4)) then
   write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&   ' chkinp : ERROR -',ch10,&
&   '  berryopt is found to be ',berryopt,ch10,&
&   '  but the only allowed values are between -3 and 4.',ch10,&
&   '  Action: change berryopt in your input file'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! toldfe, occopt, kptopt, mkmem and nspinor in case berryopt > 0

  if (berryopt /= 0) then

   if ((toldfe < tiny(one)).and.(toldff < tiny(one)).and.&
&   (tolrff < tiny(one)).and.(berryopt == 4)) then
    write(message,'(a,a,a,a,a,a,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The values of toldfe, toldff, and tolrff are found to be zero.',ch10,&
&    '  This is not allowed in a Berry phase calculation of the',ch10,&
&    '  the electric field response (berryopt = 4).',ch10,&
&    '  Action : change toldfe, toldff, or tolrff in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

   if (mkmem == 0 .and. berryopt<0) then
    write(message,'(a,a,a,a,i3,a,a,a,a,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of mkmem is found to be ',mkmem,ch10,&
&    '  This is not allowed in a Berry phase calculation of the',ch10,&
&    '  polarization, the ddk or the electric field',ch10,&
&    '  response.',ch10,&
&    '  Action : change mkmem in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

   if (occopt /= 1 .and. berryopt<0) then
    write(message,'(a,a,a,a,i3,a,a,a,a,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of occopt is found to be ',occopt,ch10,&
&    '  This is not allowed in a Berry phase calculation of the',ch10,&
&    '  polarization, the ddk or the electric field',ch10,&
&    '  response.',ch10,&
&    '  Action : put occopt = 1 in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

   if (nspinor /= 1) then
    write(message,'(a,a,a,a,i3,a,a,a,a,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of nspinor is found to be ',nspinor,ch10,&
&    '  This is not allowed in a Berry phase calculation of the',ch10,&
&    '  polarization, the ddk or the electric field',ch10,&
&    '  response.',ch10,&
&    '  Action : put nspinor = 1 in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

  end if

! Non-linear response calculations

  if (optdriver == 5) then

   if (nspinor /= 1) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of nspinor is found to be ',nspinor,ch10,&
&    '  This is not allowed in a non-linear response calculation.',ch10,&
&    '  Action : put nspinor = 1 in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

   if (nsppol /= 1) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of nsppol is found to be ',nsppol,ch10,&
&    '  This is not allowed in a non-linear response calculation.',ch10,&
&    '  Action : put nspinor = 1 in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

   if (occopt /= 1) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of occopt is found to be ',occopt,ch10,&
&    '  This is not allowed in a non-linear response calculation.',ch10,&
&    '  Action : put occopt = 1 in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

   if (mkmem == 0) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of mkmem is found to be ',mkmem,ch10,&
&    '  This is not allowed in a non-linear response calculation.',ch10,&
&    '  Action : change mkmem in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

   if ((kptopt /= 2).and.(kptopt /= 3)) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of kptopt is found to be ',kptopt,ch10,&
&    '  This is not allowed in a non-linear response calculation.',ch10,&
&    '  Action : put kptopt = 2 or 3 in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

   if ((ixc /= 3).and.(ixc /= 7)) then
    write(message,'(a,a,a,a,i3,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of ixc is found to be ',ixc,ch10,&
&    '  A non-linear response calculation can only be performed for',&
&    ' ixc = 3 or 7.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

  end if     ! optdriver == 5

  if (prepanl == 1) then

   if ((ixc /= 3).and.(ixc /= 7)) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of ixc is found to be ',ixc,ch10,&
&    '  When you prepare a non-linear response calculation (prepanl=1),',ch10,&
&    '  you should use ixc = 3 or ixc = 7.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

   if (prtden /= 1) then
    write(message,'(a,a,a,a,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  When you prepare a non-linear response calculation (prepanl=1),',ch10,&
&    '  you should put prdten=1 in your input file in order to write',ch10,&
&    '  the first-order density changes to a disk file.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

  end if     ! prepanl

! boxcutmin
  if(response==1)then
   cond_string(1)='response' ; cond_values(1)=1
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'boxcutmin',boxcutmin,0,two,iout)
  end if

! ceksph
! Make sure that ceksph is 0. The value 1 is only allowed in newsp.
! This restriction could be removed in a future version, but the information
! should then be contained in the produced wavefunction file !
  if (ceksph/=0) then
   write(message, '(6a,i4,4a)' ) ch10,&
&   'chkinp : ERROR -',ch10,&
&   '  ceksph must be 0 when used in the code abinit.',ch10,&
&   '  Its value, in the input file, is',ceksph,ch10,&
&   '  1 is allowed in newsp. Other values are not allowed. ',ch10,&
&   '  Action : change ceksph to 0 in input file.'
   call wrtout(iout,message,'COLL')
   call wrtout(06,  message,'COLL')
   ierr=ierr+1
  end if

! diecut
  if(iscf==-1)then
   cond_string(1)='iscf' ; cond_values(1)=-1
!  Checks that presently diecut is 4*ecut
   if( abs(diecut-4._dp*ecut) > 1.0d-8 ) then
    write(message, '(a,a,a,a,a,a,es14.6,a,es14.6,a,a,a)' ) ch10,&
&    'chkinp : ERROR -',ch10,&
&    ' When iscf=-1, diecut MUST be 4*ecut, while it is found that',ch10,&
&    ' ecut=',ecut,' and diecut=',diecut,'.',ch10,&
&    ' Action : change one of these values in the input file.'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    ierr=ierr+1
   end if
  end if

! diemac
  call chkdpr(0,0,cond_string,cond_values,ierr,'diemac',diemac,1,0.01_dp,iout)

! dmatpuopt
  if (usepawu==1.or.usepawu==2) then
   cond_string(1)='usepawu' ; cond_values(1)=1
   call chkint(0,1,cond_string,cond_values,ierr,&
&   'dmatpuopt',dmatpuopt,3,(/1,2,3/),0,1,iout)
  end if

! dmatudiag
  if (usepawu==1.or.usepawu==2) then
   cond_string(1)='usepawu' ; cond_values(1)=1
   call chkint(0,1,cond_string,cond_values,ierr,&
&   'dmatudiag',dmatudiag,3,(/0,1,2/),0,1,iout)
  end if

! dosdeltae
  call chkdpr(0,0,cond_string,cond_values,ierr,'dosdeltae',dosdeltae,1,0.0_dp,iout)

! ecut
! With planewaves, one must use positive ecut
  if(usewvl==0)then
   if (abs(ecut+1._dp)<tol8) then
    write(message, '(6a)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '   The input keyword "ecut" is compulsory !',ch10,&
&    '  Action : add a value for "ecut" in the input file.'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    call leave_new('COLL')
   else
    cond_string(1)='usewvl' ; cond_values(1)=usewvl
    call chkdpr(1,1,cond_string,cond_values,ierr,&
&    'ecut',ecut,1,tol8,iout)
   end if
  end if

! pawecutdg (placed here to stop before ngfftdg)
  if (usepaw==1) then
   if (abs(pawecutdg+1._dp)<tol8) then
    write(message, '(6a)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '  The input keyword "pawecutdg" is compulsory when PAW is activated !',ch10,&
&    '  Action : add a value for "pawecutdg" in the input file.'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    call leave_new('COLL')
   else
    cond_string(1)='pawecutdg' ; cond_values(1)=ecut
    call chkdpr(1,0,cond_string,cond_values,ierr,&
&    'pawecutdg',pawecutdg,1,ecut,iout)
   end if
  end if

! ecuteps
  if(optdriver==3)then
   call chkdpr(0,0,cond_string,cond_values,ierr,'ecuteps',ecuteps,1,0.0_dp,iout)
   if(fftgw<20 .and. fftgw/=0)then
    if(ecutwfn<ecuteps)then
     write(message,'(4a,es16.6,a,es16.6,a,6a)')ch10,&
&     ' chkinp : ERROR -',ch10,&
&     '  The values of ecutwfn and ecuteps are ', ecutwfn,' and ',ecuteps,ch10,&
&     '  With fftgw lower than 20, one expect ecuteps to be smaller than ecutwfn.',ch10,&
&     '  Indeed, one is wasting memory without gaining CPU time or accuracy.',ch10,&
&     '  Action : use another value of fftgw (e.g. 21), or adjust ecutwfn with ecuteps.'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')
    end if
   end if
  end if

! ecutsigx
  if(optdriver==4)then
   call chkdpr(0,0,cond_string,cond_values,ierr,'ecutsigx',ecutsigx,1,0.0_dp,iout)
   if(fftgw<20)then
    if(ecutwfn<ecutsigx)then
     write(message,'(4a,es16.6,a,es16.6,a,6a)')ch10,&
&     ' chkinp : ERROR -',ch10,&
&     '  The values of ecutwfn and ecutsigx are ', ecutwfn,' and ',ecutsigx,ch10,&
&     '  With fftgw lower than 20, one expect ecutsigx to be smaller than ecutwfn.',ch10,&
&     '  Indeed, one is wasting memory without gaining CPU time or accuracy.',ch10,&
&     '  Action : use another value of fftgw (e.g. 21), or adjust ecutwfn with ecutsigx.'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')
    end if
   end if
  end if

! ecutsm
  call chkdpr(0,0,cond_string,cond_values,ierr,'ecutsm',ecutsm,1,0.0_dp,iout)
! With non-zero optcell, one must use non-zero ecutsm
  if(optcell/=0 )then
   cond_string(1)='optcell' ; cond_values(1)=optcell
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'ecutsm',ecutsm,1,tol8,iout)
  end if

! exchmix
  call chkdpr(0,0,cond_string,cond_values,ierr,'exchmix',exchmix,1,0.0_dp,iout)

! fixmom
  if(nsppol==1)then
   cond_string(1)='nsppol' ; cond_values(1)=1
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'fixmom',fixmom,0,-99.99_dp,iout)
  end if
  if(optdriver==1)then
   cond_string(1)='optdriver' ; cond_values(1)=1
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'fixmom',fixmom,0,-99.99_dp,iout)
  end if
  if(optdriver==2)then
   cond_string(1)='optdriver' ; cond_values(1)=2
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'fixmom',fixmom,0,-99.99_dp,iout)
  end if
  if(prtdos==1)then
   cond_string(1)='prtdos' ; cond_values(1)=1
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'fixmom',fixmom,0,-99.99_dp,iout)
  end if

! fftgw
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'fftgw',fftgw,8,(/00,01,10,11,20,21,30,31/),0,0,iout)

! frzfermi
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'frzfermi',frzfermi,2,(/0,1/),0,0,iout)

! gwpara
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'gwpara',gwpara,3,(/0,1,2/),0,0,iout)

! if(gwpara==2)then
! write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
! &   ' chkinp: ERROR -',ch10,&
! &   '  The value gwpara=2 is temporarily forbidden,',ch10,&
! &   '  except for expert users.',ch10,&
! &   ' Action : modify gwpara, or suppress the present test - if you are an expert user.'
! call wrtout(6,message,'COLL')
! call leave_new('COLL')
! endif

! gwmem
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'gwmem',gwmem,4,(/0,1,10,11/),0,0,iout)

! gwcomp
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'gwcomp',gwcomp,2,(/0,1/),0,0,iout)
! When usepaw/=0, gwcomp must be 0
  if(usepaw/=0)then
   cond_string(1)='usepaw' ; cond_values(1)=usepaw
!  Make sure that gwcomp==0
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'gwcompv',gwcomp,1,(/0/),0,0,iout)
  end if
  if(gwcomp/=0 .and. optdriver==3 .and. ( awtr /=1 .or. spmeth /=0 )) then
   write(message,'(8a)' ) ch10,&
&   ' chkinp: ERROR -',ch10,&
&   '  When gwcomp/=0, the Adler-Wiser formula with time-reversal should be used',ch10,&
&   '  Action : set awtr to 1 or/and spmeth to 0'
   call wrtout(iout,message,'COLL')
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! iatsph between 1 and natom
  maxiatsph=maxval(iatsph(:))
  miniatsph=minval(iatsph(:))
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'iatsph',miniatsph,1,(/1/),1,1,iout)
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'iatsph',maxiatsph,1,(/natom/),-1,natom,iout)

! icoulomb
  call chkint(1,1,cond_string,cond_values,ierr,&
&  'icoulomb',dtsets(idtset)%icoulomb,2,(/0,1/),0,0,iout)
  if(dtsets(idtset)%icoulomb == 1)then
   if (dtsets(idtset)%nspden > 2) then
    write(message,'(a,a,a,a,i3,a,a,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of nspden is found to be ', dtsets(idtset)%nspden, ch10,&
&    '  The computation of Hartree potential with free boundary counditions', ch10, &
&    '  is not allowed with spin (collinear or not).',ch10,&
&    '  Action : put nspden = 1 or 2 in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
  end if


! intxc
  if(iscf==-1)then
   cond_string(1)='iscf' ; cond_values(1)=-1
!  Make sure that intxc is 0
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'intxc',intxc,1,(/0/),0,0,iout)
  end if
! TEMPORARY
  if(optdriver==1)then
   cond_string(1)='optdriver' ; cond_values(1)=1
!  Make sure that intxc is 0
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'intxc',intxc,1,(/0/),0,0,iout)
  end if

! ionmov
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'ionmov',ionmov,13,(/0,1,2,3,4,5,6,7,8,9,12,13,20/),0,0,iout)
! When optcell/=0, ionmov must be 2 or 3
  if(optcell/=0)then
   cond_string(1)='optcell' ; cond_values(1)=optcell
!  Make sure that ionmov==2, 3 or 13
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'ionmov',ionmov,3,(/2,3,13/),0,0,iout)
  end if

! iprcch
!   if (abs(iprcch)>=1.and.iscf>=10.and.nspden==4) then
!    write(message,'(8a)')ch10,&
! &   ' chkinp : ERROR -',ch10,&
! &   '  When non-collinear magnetism is activated (nspden=4),',ch10,&
! &   '  iprcch/=0 is not compatible with SCF mixing on density (iscf>=10) !',ch10,&
! &   '  Action : choose SCF mixing on potential (iscf<10) or change iprcch value.'
!    call wrtout(6,message,'COLL')
!    ierr=ierr+1
!   end if
  if (iprcch<0.and.mod(iprcel,100)>=61.and.(iprcel<71.or.iprcel>79)) then
   cond_string(1)='iprcel';cond_values(1)=iprcel
   call chkint(0,2,cond_string,cond_values,ierr,&
&   'iprcch',iprcch,1,(/0/),1,0,iout)
  end if

! iprcel
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'iprcel',iprcel,1,(/0/),1,21,iout)
  if(nsppol==2 .and. (occopt>=3 .and. occopt<=7) )then
   cond_string(1)='nsppol' ; cond_values(1)=nsppol
   cond_string(2)='occopt' ; cond_values(2)=occopt
   call chkint(2,2,cond_string,cond_values,ierr,&
&   'iprcel',iprcel,1,(/0/),0,0,iout)
  end if

! iprctfvw
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'iprctfvw',iprctfvw,4,(/0,1,2,3/),0,0,iout)
  if((nsppol/=1).or.(iscf > 9))then
   cond_string(1)='nsppol' ; cond_values(1)=nsppol
   cond_string(2)='iscf' ; cond_values(2)=iscf
   call chkint(2,2,cond_string,cond_values,ierr,&
&   'iprctfvw',iprctfvw,1,(/0/),0,0,iout)
  end if


! iscf
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'iscf',iscf,18,(/-3,-2,-1,1,2,3,4,5,6,7,11,12,13,14,15,16,17,22/),0,0,iout)
! If ionmov==4, iscf must be 2, 12, 5 or 6.
  if(ionmov==4)then
   cond_string(1)='ionmov' ; cond_values(1)=4
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'iscf',iscf,4,(/2,12,5,6/),0,0,iout)
  end if
! If PAW, iscf cannot be -1, 11
  if (usepaw==1) then
   cond_string(1)='PAW' ; cond_values(1)=1
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'iscf',iscf,11,(/-3,-2,2,3,4,7,12,13,14,17,22/),0,0,iout)
  end if
! Mixing on density is only allowed for GS calculations
  if(optdriver>0)then
   cond_string(1)='optdriver' ; cond_values(1)=optdriver
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'iscf',iscf,1,(/9/),-1,9,iout)
  end if
! mixing on density is not allowed with some preconditioners
  if (iprctfvw /= 0) then
   cond_string(1)='iprctfvw' ; cond_values(1)=iprctfvw
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'iscf',iscf,1,(/9/),-1,9,iout)
  end if
! When pawoptmix=1 and nspden=4, iscf must be >=10
  if(pawoptmix/=0.and.nspden==4)then
   cond_string(1)='nspden'    ; cond_values(1)=nspden
   cond_string(2)='pawoptmix' ; cond_values(2)=pawoptmix
   call chkint(2,2,cond_string,cond_values,ierr,&
&   'iscf',iscf,1,(/10/),1,10,iout)
  end if


! istwfk
  if(response==1 .and. maxval( abs(istwfk(:)-1) ) >0)then
!  Force istwfk to be 1 for RF calculations
!  Other choices cannot be realized yet, because of the ddk perturbation.
   write(message,'(8a)' ) ch10,&
&   ' chkinp: ERROR -',ch10,&
&   '  When response==1, all the components of istwfk must be 1.',ch10,&
&   '  Not yet programmed for time-reversal symmetry.',ch10,&
&   '  Action : set istwfk to 1 for all k-points'
   call wrtout(iout,message,'COLL')
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
  if(nbandkss/=0 .and. maxval( abs(istwfk(:)-1) ) >0)then
   write(message,'(8a)' ) ch10,&
&   ' chkinp: ERROR -',ch10,&
&   '  When nbandkss/=0, all the components of istwfk must be 1.',ch10,&
&   '  Not yet programmed for time-reversal symmetry.',ch10,&
&   '  Action : set istwfk to 1 for all k-points'
   call wrtout(iout,message,'COLL')
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
  if(berryopt/=0 .and. maxval(istwfk(:))/=1)then
   write(message,'(8a)' ) ch10,&
&   ' chkinp: ERROR -',ch10,&
&   '  When berryopt/=0, all the components of istwfk must be 1.',ch10,&
&   '  Not yet programmed for time-reversal symmetry.',ch10,&
&   '  Action : set istwfk to 1 for all k-points'
   call wrtout(iout,message,'COLL')
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
  if ((wfoptalg==4.or.wfoptalg==14).and.maxval(istwfk(:)-2)>0) then
   write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&   ' chkinp : ERROR -',ch10,&
&   '  Only the gamma point can use time-reversal and wfoptalg=4 or 14',ch10,&
&   '  Action : put istwfk to 1 or remove k points with half integer coordinates ',ch10,&
&   '  Also contact ABINIT group to say that you need that option.'
   call wrtout(iout,message,'COLL')
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! ixc
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'ixc',ixc,29,(/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,20,21,22,23,26,27,30,31,32,33,34/),0,0,iout)
  if(iscf==-1)then
   cond_string(1)='iscf' ; cond_values(1)=-1
!  Make sure that ixc is 1, 7, 8, 20, 21 or 22
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'ixc',ixc,6,(/1,7,8,20,21,22/),0,0,iout)
  end if
  if(response==1)then
   cond_string(1)='response' ; cond_values(1)=1
!  Make sure that ixc is between 0 and 9, or 11, 12, 14, 15, or 23
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'ixc',ixc,15,(/0,1,2,3,4,5,6,7,8,9,11,12,14,15,23/),0,0,iout)
  end if
  if(nspden/=1)then
   cond_string(1)='nspden' ; cond_values(1)=nspden
!  Make sure that ixc is 0, 1 or the gga
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'ixc',ixc,15,(/0,1,7,8,9,11,12,13,14,15,16,17,23,26,27/),0,0,iout)
  end if

! kptnrm and kpt
! Coordinates components must be between -1 and 1.
  if(kptnrm<1.0-1.0d-10)then
   write(message, '(a,a,a,a,es22.14,a,a,a)' ) ch10,&
&   ' chkinp: ERROR -',ch10,&
&   '  The input variable kptnrm is',kptnrm,' while it must be >=1.0_dp.',&
&   ch10,'  Action : change the input variable kptnrm.'
   call wrtout(iout,message,'COLL')
   call wrtout(06,  message,'COLL')
   call leave_new('COLL')
  end if
  do ikpt=1,nkpt
   do mu=1,3
    if ( abs(kpt(mu,ikpt))> kptnrm*1.0000001_dp ) then
     write(message, '(a,a,a,a,i5,a,a,a,a,3es22.14,a,a,a,a)' ) ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '  For k point number',ikpt,'  the reduced coordinates',ch10,&
&     '  generated by the input variables kpt and kptnrm are',ch10,&
&     kpt(1,ikpt)/kptnrm,kpt(2,ikpt)/kptnrm,kpt(3,ikpt)/kptnrm,ch10,&
&     '  while they must be between -1.0_dp and 1.0_dp (included).',ch10,&
&     '  Action : check kpt and kptnrm in the input file.'
     call wrtout(iout,message,'COLL')
     call wrtout(06,  message,'COLL')
     call leave_new('COLL')
    end if
   end do
  end do

! jellslab
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'jellslab',jellslab,2,(/0,1/),0,0,iout)

  if (jellslab==1) then
!  slabwsrad must be positive
   cond_string(1)='jellslab' ; cond_values(1)=jellslab
   call chkdpr(1,0,cond_string,cond_values,ierr,&
&   'slabwsrad',slabwsrad,1,zero,iout)
!  slabzbeg must be positive
   call chkdpr(1,0,cond_string,cond_values,ierr,&
&   'slabzbeg',slabzbeg,1,zero,iout)
!  slabzend must be bigger than slabzbeg
   call chkdpr(1,0,cond_string,cond_values,ierr,&
&   'slabzend',slabzend,1,slabzbeg,iout)
!  rprimd(3,3) must be bigger than slabzend
   call chkdpr(1,0,cond_string,cond_values,ierr,&
&   'rprimd33',rprimd(3,3),1,slabzend,iout)
!  Third real space primitive translation has to be orthogonal to the other ones,
!  actually, for convenience it is useful that rprimd is something like:
!  a  b  0
!  c  d  0
!  0  0  e
   if(abs(rprimd(1,3))+abs(rprimd(2,3))+abs(rprimd(3,1))+abs(rprimd(3,2))>tol12) then
    write(message, '(a,a,a,a,a,a)' ) ch10,&
&    ' invars2 : ERROR - ',ch10,&
&    '  Third real space vector is not orthogonal to the other ones,',ch10,&
&    '  this is needed to use jellium'
    call wrtout(6,message,'COLL')
    ierr=ierr+1
   end if

!  Atoms have to be placed in the vacuum space
   do iatom=1,natom
    zatom=(xred(3,iatom)-anint(xred(3,iatom)-half+tol6))*rprimd(3,3)
    if(abs(zatom-slabzbeg)<tol8 .or. abs(zatom-slabzend)<tol8) then
     if(znucl(typat(iatom))>tol6) then
      write(message,'(4a,i5,a)') ch10,&
&      ' invars2 : WARNING -',ch10,&
&      ' atom number=',iatom,' lies precisely on the jellium edge !'
      call wrtout(6,message,'COLL')
     end if
     cycle
    end if
    if(zatom>slabzbeg .and. zatom<slabzend) then
     write(message, '(a,a,a,a,i5,a)' ) ch10,&
&     ' invars2 : ERROR - ',ch10,&
&     ' atom number=',iatom,' is inside the jellium slab.'
     call wrtout(6,message,'COLL')
     ierr=ierr+1
    end if
   end do
  end if

! kssform
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'kssform',kssform,4,(/0,1,2,3/),0,0,iout)

! localrdwf
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'localrdwf',localrdwf,2,(/0,1/),0,0,iout)
  if(mkmem==0)then
   cond_string(1)='mkmem' ; cond_values(1)=mkmem
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'localrdwf',localrdwf,1,(/1/),0,0,iout)
  end if
  if(mkqmem==0)then
   cond_string(1)='mkqmem' ; cond_values(1)=mkqmem
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'localrdwf',localrdwf,1,(/1/),0,0,iout)
  end if
  if(mk1mem==0)then
   cond_string(1)='mk1mem' ; cond_values(1)=mk1mem
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'localrdwf',localrdwf,1,(/1/),0,0,iout)
  end if

! mixalch
! For each type of atom, the sum of the psp components
! must be one.
  if(ntypalch>0)then
   do itypat=1,ntypalch
    sumalch=sum(mixalch(:,itypat))
    if(abs(sumalch-one)>tol10)then
     if(npspalch<=6)then
      write(message, '(2a,6es12.4)' ) ch10,&
&      ' chkinp : mixalch(:,itypat)=',mixalch(:,itypat)
     end if
     call wrtout(iout,message,'COLL')
     call wrtout(06,  message,'COLL')
     write(message, '(4a,i4,2a,f8.2,4a)' ) ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '  For the alchemical atom number',itypat,ch10,&
&     '  the sum of the pseudopotential coefficients is',sumalch,ch10,&
&     '  while it should be one.',ch10,&
&     '  Action : check the content of the input variable mixalch.'
     call wrtout(iout,message,'COLL')
     call wrtout(06,  message,'COLL')
     call leave_new('COLL')
    end if
   end do
  end if

! mffmem
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'mffmem',mffmem,2,(/0,1/),0,0,iout)

! natom
  if(prtgeo>0)then
   cond_string(1)='prtgeo' ; cond_values(1)=prtgeo
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'natom',natom,1,(/9999/),-1,9999,iout)
  end if

! nband
! Make sure all nband(nkpt) are > 0
  do isppol=1,nsppol
   do ikpt=1,nkpt
    if (nband(ikpt+(isppol-1)*nkpt)<=0) then
     cond_string(1)='ikpt' ; cond_values(1)=ikpt
     cond_string(2)='isppol' ; cond_values(2)=isppol
     call chkint(0,2,cond_string,cond_values,ierr,&
&     'nband',nband(ikpt+(isppol-1)*nkpt),1,(/1/),1,1,iout)
    end if
   end do
  end do

  if(mpi_enreg%nproc/=1 .and. nsppol==2 .and. dtsets(idtset)%usewvl == 0)then
   do ikpt=1,nkpt
    if (nband(ikpt)/=nband(ikpt+nkpt)) then
     write(message, '(8a,i4,a,2i5,a)' ) ch10,&
&     ' chkinp : in the parallel k point case, for each k point,',ch10,&
&     '  the number of bands in the spin up case must be equal to',ch10,&
&     '  the number of bands in the spin down case.',ch10,&
&     '  This is not the case for the k point number :',ikpt,&
&     '  The number of bands spin up and down are :',nband(ikpt),nband(ikpt+nkpt),&
&     '  Action : change nband, or use the sequential version of ABINIT.'
     call wrtout(iout,message,'COLL')
     call wrtout(06,  message,'COLL')
     call leave_new('COLL')
    end if
   end do
  end if

! nbdblock
! Must be larger or equal to 1
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'nbdblock',nbdblock,1,(/1/),1,1,iout)
! When wfoptalg==0, nbdblock must be 1
  if(mod(wfoptalg,10)==0)then
   cond_string(1)='wfoptalg' ; cond_values(1)=0
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'nbdblock',nbdblock,1,(/1/),0,0,iout)
  end if
! When wfoptalg==2, nbdblock must be 1
  if(wfoptalg==2)then
   cond_string(1)='wfoptalg' ; cond_values(1)=2
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'nbdblock',nbdblock,1,(/1/),0,0,iout)
  end if
! When wfoptalg==3, nbdblock must be 1, and iscf must be -2
  if(wfoptalg==3)then
   cond_string(1)='wfoptalg' ; cond_values(1)=3
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'nbdblock',nbdblock,1,(/1/),0,0,iout)
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'iscf',iscf,1,(/-2/),0,0,iout)
  end if
! When wfoptalg==4, nbdblock must be a divisor of nband
  if(mod(wfoptalg,10)==4)then
   do isppol=1,nsppol
    do ikpt=1,nkpt
     if(mod(nband(ikpt+(isppol-1)*nkpt),nbdblock)/=0) then
      write(message, '(8a)' ) ch10,&
&      ' chkinp: ERROR -',ch10,&
&      '  For the moment, when wfoptalg=4,',ch10,&
&      '  nband must be a multiple of nbdblock.',ch10,&
&      '  Action : check the value of the input variable nbdblock.'
      call wrtout(iout,message,'COLL')
      call wrtout(06,  message,'COLL')
      call leave_new('COLL')
     end if
    end do
   end do
  end if

! nberry
! must be between 0 and 20
  if(berryopt/=0)then
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'nberry',nberry,1,(/0/),1,0,iout)
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'nberry',nberry,1,(/20/),-1,20,iout)
   if(mpi_enreg%paral_compil==1)then
!   MPI Parallel case
    if ((nberry/=0).and.(berryopt > 0).and.(berryopt /= 4)) then
     write(message,'(a,a,a,a,a,a,a,a,i4,a,a,a)')ch10,&
&     ' chkinp : ERROR -',ch10,&
&     '  Berry phase calculation of polarisation with positive berryopt is not',ch10,&
&     '  allowed in the parallel version of ABINIT.',ch10,&
&     '  So, the value of nberry=',nberry,' is not allowed,',ch10,&
&     '  Action : change berryopt to negative values or change nberry, or use the sequential version.'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')
    end if
   end if
  end if

! nbandkss
! Must be larger or equal to -1
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'nbandkss',nbandkss,1,(/-1/),1,-1,iout)
! When nspinor==2, nbandkss must be 0 unless we use kssform=3
  if(nspinor==2.and.kssform/=3)then
   cond_string(1)='nspinor' ; cond_values(1)=2
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'nbandkss',nbandkss,1,(/0/),0,0,iout)
  end if
! When ionmov/=0
  if(ionmov/=0 .and. nbandkss/=0)then
   write(message,'(14a)') ch10,&
&   ' chkinp: WARNING -',ch10,&
&   '  Ions (or cell) are allowed to move (ionmov/=0),',ch10,&
&   '  and a _KSS file is requested (nbandkss/=0).',ch10,&
&   '  A _KSS file will be created at each geometry-optimisation step.',ch10,&
&   '  Note that this is time consuming !',ch10,&
&   '  Action : use datasets (one for geometry optimisation,',ch10,&
&   '           one for states output).'
   call wrtout(6,message,'COLL')
  end if

! ngeohist
! Must be larger or equal to 1
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'ngeohist',dtsets(idtset)%ngeohist,1,(/1/),1,1,iout)

! npwkss
! Must be larger or equal to -1
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'npwkss',npwkss,1,(/-1/),1,-1,iout)

! nfft and nfftdg
! Must have nfft<=nfftdg
  if (usepaw==1) then
   nfft  =ngfft(1)  *ngfft(2)  *ngfft(3)
   nfftdg=ngfftdg(1)*ngfftdg(2)*ngfftdg(3)
   cond_string(1)='nfft' ; cond_values(1)=nfft
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'nfftdg',nfftdg,1,(/0/),1,nfft,iout)
  end if

! nkpt
! Must be larger or equal to 1
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'nkpt',nkpt,1,(/1/),1,1,iout)
! If prtdos==2 or 3, must be larger or equal to 2
  if(prtdos==2 .or. prtdos==3)then
   cond_string(1)='prtdos' ; cond_values(1)=prtdos
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'nkpt',nkpt,1,(/2/),1,2,iout)
  end if
! Must be smaller than 50 if iscf=-2 (band structure)
! while prteig=0 and prtvol<2, except if kptopt>0
  if(iscf==-2 .and. prteig==0 .and. prtvol<2 .and. kptopt<=0)then
   cond_string(1)='iscf'   ; cond_values(1)=iscf
   cond_string(2)='prteig' ; cond_values(2)=prteig
   cond_string(3)='prtvol' ; cond_values(3)=prtvol
   call chkint(1,3,cond_string,cond_values,ierr,&
&   'nkpt',nkpt,1,(/50/),-1,50,iout)
  end if

! npband
! Must be larger or equal to 1
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'npband',npband,1,(/1/),1,1,iout)

! npfft
! Must be larger or equal to 1
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'npfft',npfft,1,(/1/),1,1,iout)
! If usepaw==1 and pawmixdg==0, npfft must be equal to 1
  if(usepaw==1 .and. pawmixdg==0)then
   cond_string(1)='usepaw  ' ; cond_values(1)=usepaw
   cond_string(2)='pawmixdg' ; cond_values(2)=pawmixdg
   call chkint(1,2,cond_string,cond_values,ierr,&
&   'npfft',npfft,1,(/1/),0,0,iout)
  end if

! npkpt
! Must be larger or equal to 1
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'npkpt',npkpt,1,(/1/),1,1,iout)

! nproj
! If there is more than one projector for some angular momentum
! channel of some pseudopotential
  do ilang=0,3
!  nprojmax(ilang)=maxval(pspheads(1:npsp)%nproj(ilang)) ! Likely problems with HP compiler
   nprojmax(ilang)=pspheads(1)%nproj(ilang)
   if(npsp>2)then
    do ii=2,npsp
     nprojmax(ilang)=max(pspheads(ii)%nproj(ilang),nprojmax(ilang))
    end do
   end if
  end do
  if(maxval(nprojmax(0:3))>1)then
   if(nbandkss/=0 .and. (kssform==1 .or. kssform==2))then
    ierr=ierr+1
   end if
  end if

! nqpt
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'nqpt',nqpt,2,(/0,1/),0,0,iout)

! nscforder
  call chkint(0, 0, cond_string, cond_values, ierr,&
&  'nscforder', dtsets(idtset)%nscforder, 10, &
&  (/ 8, 14, 16, 20, 24, 30, 40, 50, 60, 100 /), 0, 0, iout)

! nspden
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'nspden',nspden,3,(/1,2,4/),0,0,iout)
! When nsppol=2, nspden must be 2
  if(nsppol==2)then
   cond_string(1)='nsppol' ; cond_values(1)=2
!  Make sure that nspden is 2
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'nspden',nspden,1,(/2/),0,0,iout)
  end if
  if(nspden==2 .and. nsppol==1 .and. response==1)then
   write(message,'(16a)')ch10,&
&   ' chkinp : ERROR -',ch10,&
&   '  nspden==2 together with nsppol==1 is not allowed',ch10,&
&   '  for response function calculations.',ch10,&
&   '  For antiferromagnetic materials, use nspden==2 and nsppol=2.',ch10,&
&   '  In this case, Shubnikov symmetries will be used to decrease',ch10,&
&   '  the number of perturbations. In a future version, it will also be',ch10,&
&   '  used to decrease the number of spin components (to be coded).',ch10,&
&   '  Action : change nsppol to 1, or check nspden.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
  if(nspden==4.and.response==1)then
   write(message,'(6a)')ch10,&
&   ' chkinp : ERROR -',ch10,&
&   '  nspden==4 not allowed in response formalism.',ch10,&
&   '  Non collinear magnetism not yet implemented in perturbative treatment.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
! When iprcch<0 or 3, nspden must be 1 or 2
  if(iprcch<0.or.iprcch==3)then
   cond_string(1)='iprcch' ; cond_values(1)=iprcch
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'nspden',nspden,2,(/1,2/),0,0,iout)
  end if
! When ionmov=4 and iscf>10, nspden must be 1 or 2
  if(ionmov==4.and.iscf>10)then
   cond_string(1)='ionmov' ; cond_values(1)=ionmov
   cond_string(1)='iscf' ; cond_values(1)=iprcch
   call chkint(1,2,cond_string,cond_values,ierr,&
&   'nspden',nspden,2,(/1,2/),0,0,iout)
  end if
! When iprcel>20, nspden must be 1 or 2
  if(iprcel>20)then
   cond_string(1)='iprcel' ; cond_values(1)=iprcel
   call chkint(1,2,cond_string,cond_values,ierr,&
&   'nspden',nspden,2,(/1,2/),0,0,iout)
  end if

! nspinor
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'nspinor',nspinor,2,(/1,2/),0,0,iout)
! When nspden=2, nspinor must be 1
  if(nspden==2)then
   cond_string(1)='nspden' ; cond_values(1)=2
!  Make sure that nspinor is 1
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'nspinor',nspinor,1,(/1/),0,0,iout)
  end if
! When nspden=4, nspinor must be 2
  if(nspden==4)then
   cond_string(1)='nspden' ; cond_values(1)=4
!  Make sure that nspinor is 2
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'nspinor',nspinor,1,(/2/),0,0,iout)
  end if
! When iscf=-1, nspinor must be 1
  if(iscf==-1)then
   cond_string(1)='iscf' ; cond_values(1)=-1
!  Make sure that nsppol is 1
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'nspinor',nspinor,1,(/1/),0,0,iout)
  end if
! spin-orbit is not implemented for the strain perturbation
  if(rfstrs/=0)then
   cond_string(1)='rfstrs' ; cond_values(1)=rfstrs
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'nspinor',nspinor,1,(/1/),0,0,iout)
  end if

! nsppol
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'nsppol',nsppol,2,(/1,2/),0,0,iout)

! nsym
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'nsym',nsym,1,(/1/),1,1,iout)
  if(nspden==4)then
   cond_string(1)='nspden' ; cond_values(1)=4
!  Make sure that nsym is 1
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'nsym',nsym,1,(/1/),0,0,iout)
  end if

! ntypalch
  if (usepaw==1) then
   cond_string(1)='pspcod(atom_type)' ; cond_values(1)=7
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'ntypalch',ntypalch,1,(/0/),0,0,iout)
  end if

! occopt
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'occopt',occopt,8,(/0,1,2,3,4,5,6,7/),0,0,iout)
! When prtdos==1, occopt must be between 3 and 7
  if(prtdos==1)then
   cond_string(1)='prtdos' ; cond_values(1)=1
!  Make sure that occopt is 3,4,5,6, or 7
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'occopt',occopt,5,(/3,4,5,6,7/),0,0,iout)
  end if

! occ
! Do following tests only for occopt==0 or 2, when occupation numbers are needed
  if ((iscf>0.or.iscf==-1.or.iscf==-3) .and. (occopt==0 .or. occopt==2) ) then
!  make sure occupation numbers (occ(n)) were defined:
   sumocc=zero
   bantot=0
   do isppol=1,nsppol
    do ikpt=1,nkpt
     do iband=1,nband(ikpt+(isppol-1)*nkpt)
      bantot=bantot+1
      sumocc=sumocc+occ(bantot)
      if (occ(bantot)<zero) then
       write(message, '(a,a,a,a,2i6,a,e20.10,a,a,a)' )  ch10,&
&       ' chkinp: ERROR -',ch10,&
&       'iband,ikpt=',iband,ikpt,' has negative occ=',occ(bantot),' =>stop',&
&       ch10,'  Action : correct this occupation number in input file.'
       call wrtout(iout,message,'COLL')
       call wrtout(06,  message,'COLL')
       call leave_new('COLL')
      end if
     end do
    end do
   end do
   if (sumocc<=1.0d-8) then
    write(message, '(a,a,a,a,1p,e20.10,a,a,a)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '  Sum of occ=',sumocc, ' =>occ not defined => stop',ch10,&
&    '  Action : correct the array occ in input file.'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    ierr=ierr+1
   end if
  end if

! optcell
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'optcell',optcell,10,(/0,1,2,3,4,5,6,7,8,9/),0,0,iout)
! With berryopt=4, one must have optcell==0
  if(berryopt==4)then
   cond_string(1)='berryopt' ; cond_values(1)=berryopt
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'optcell',optcell,1,(/0/),0,0,iout)
  end if

! optdriver and PAW
  if(usepaw==1)then
   cond_string(1)='use_PAW' ; cond_values(1)=usepaw
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'optdriver',optdriver,4,(/0,1,3,4/),0,0,iout)
  end if

! optforces
! When ionmov>0, optforces must be >0
  if(ionmov>0)then
   cond_string(1)='ionmov' ; cond_values(1)=ionmov
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'optforces',optforces,2,(/1,2/),0,0,iout)
  end if
! When iscf=22, optforces must be 0 or 2
  if(iscf==22)then
   cond_string(1)='iscf' ; cond_values(1)=iscf
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'optforces',optforces,2,(/0,2/),0,0,iout)
  end if

! optstress
! When optcell>0, optstress must be >0
  if(optcell>0)then
   cond_string(1)='optcell' ; cond_values(1)=optcell
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'optstress',optstress,1,(/1/),0,0,iout)
  end if

! pawlcutd
  if (usepaw==1) then
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'pawlcutd',pawlcutd,1,(/0/),1,0,iout)
  end if

! pawlmix
  if (usepaw==1) then
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'pawlmix',pawlmix,1,(/0/),1,1,iout)
  end if

! pawmixdg
  if (usepaw==1) then
   if(ionmov==4)then
    cond_string(1)='ionmov' ; cond_values(1)=ionmov
    call chkint(1,1,cond_string,cond_values,ierr,&
&    'pawmixdg',pawmixdg,1,(/1/),0,0,iout)
   end if
   if(iscf==5.or.iscf==6.or.iscf==15.or.iscf==16)then
    cond_string(1)='iscf' ; cond_values(1)=iscf
    call chkint(1,1,cond_string,cond_values,ierr,&
&    'pawmixdg',pawmixdg,1,(/1/),0,0,iout)
   end if
  end if

! pawnhatxc
  if (usepaw==1) then
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'pawnhatxc',pawnhatxc,2,(/0,1/),0,0,iout)
  end if

! pawnzlm
  if (usepaw==1) then
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'pawnzlm',pawnzlm,2,(/0,1/),0,0,iout)
  end if

! pawoptmix
  if (usepaw==1) then
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'pawoptmix',pawoptmix,2,(/0,1/),0,0,iout)
  end if

! pawprtdos
  if (usepaw==1) then
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'pawprtdos',pawprtdos,3,(/0,1,2/),0,0,iout)
  end if

! pawprtvol
  if (usepaw==1) then
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'pawprtvol',pawprtvol,7,(/-3,-2,-1,0,1,2,3/),0,0,iout)
  end if

! pawspnorb
  if (usepaw==1) then
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'pawspnorb',pawspnorb,2,(/0,1/),0,0,iout)
   if (iscf<10) then
    cond_string(1)='iscf' ; cond_values(1)=iscf
    call chkint(1,1,cond_string,cond_values,ierr,&
&    'pawspnorb',pawspnorb,1,(/0/),0,0,iout)
   end if
   if (pawspnorb==1.and.(kptopt==1.or.kptopt==2)) then
    write(message, '(10a)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '  When spin-orbit coupling is activated (pawspnorb=1),',ch10,&
&    '  time-reversal symmetry is broken; k-points cannot',ch10,&
&    '  be generated using TR-symmetry.',ch10,&
&    '  Action: choose kptopt different from 1 or 2.'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    ierr=ierr+1
   end if
   if (pawspnorb==1.and.nspden==1) then
    write(message, '(10a)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '  When spin-orbit coupling is activated (pawspnorb=1),',ch10,&
&    '  nspden=1 may give incorrect results. This',ch10,&
&    '  feature has been temporary desactivated.',ch10,&
&    '  Action: choose nspden=4 in input file.'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    ierr=ierr+1
   end if
  end if

! pawstgylm
  if (usepaw==1) then
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'pawstgylm',pawstgylm,2,(/0,1/),0,0,iout)
  end if

! pawusecp
  if (usepaw==1) then
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'pawusecp',pawusecp,2,(/0,1/),0,0,iout)
   if (mkmem/=0)then
    cond_string(1)='mkmem' ; cond_values(1)=mkmem
    call chkint(1,1,cond_string,cond_values,ierr,&
&    'pawusecp',pawusecp,1,(/1/),0,0,iout)
   end if
   if (mk1mem/=0)then
    cond_string(1)='mk1mem' ; cond_values(1)=mk1mem
    call chkint(1,1,cond_string,cond_values,ierr,&
&    'pawusecp',pawusecp,1,(/1/),0,0,iout)
   end if
   if (mkqmem/=0)then
    cond_string(1)='mkqmem' ; cond_values(1)=mkqmem
    call chkint(1,1,cond_string,cond_values,ierr,&
&    'pawusecp',pawusecp,1,(/1/),0,0,iout)
   end if
  end if

! pawxcdev
  if (usepaw==1) then
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'pawxcdev',pawxcdev,3,(/0,1,2/),0,0,iout)
!  IF GGA must have pawxcdev>0
   if (ixc>9)then
    cond_string(1)='ixc' ; cond_values(1)=ixc
    call chkint(1,1,cond_string,cond_values,ierr,&
&    'pawxcdev',pawxcdev,2,(/1,2/),0,0,iout)
   end if
  end if

! prtdensph
  if (usepaw==1) then
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'prtdensph',prtdensph,2,(/0,1/),0,0,iout)
  end if

! prtdos
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'prtdos',prtdos,4,(/0,1,2,3/),0,0,iout)

! prtdosm
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'prtdosm',prtdosm,2,(/0,1/),0,0,iout)
  if(usepaw==1.and.pawprtdos==1)then
   cond_string(1)='pawprtdos' ; cond_values(1)=pawprtdos
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'prtdosm',prtdosm,1,(/0/),0,0,iout)
  end if

! prtfsurf
! only one shift allowed (gamma)
  if (prtfsurf /= 0) then
   if (abs(kptrlatt(1,2))+abs(kptrlatt(1,3))+abs(kptrlatt(2,3))+&
   abs(kptrlatt(2,1))+abs(kptrlatt(3,1))+abs(kptrlatt(3,2)) /= 0 ) then
    write(message,'(4a)') ch10,&
&    ' prtfsurf does not work with non-diagonal kptrlatt ', ch10, &
&    ' Action: set nshift 1 and shiftk 0 0 0'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if (nshiftk > 1) then
    write(message,'(4a)') ch10,&
&    ' prtfsurf does not work with multiple kpt shifts ', ch10, &
&    ' Action: set nshift 1 and shiftk 0 0 0'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if (sum(abs(shiftk)) /= 0) then
    write(message,'(4a)') ch10,&
&    ' prtfsurf does not work with non-zero kpt shift ', ch10, &
&    ' Action: set nshift 1 and shiftk 0 0 0'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
  end if

! prtstm
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'prtstm',prtstm,1,(/0/),1,0,iout)
  if(optdriver/=0)then
   cond_string(1)='optdriver' ; cond_values(1)=optdriver
   call chkint(0,1,cond_string,cond_values,ierr,&
&   'prtstm',prtstm,1,(/0/),0,0,iout)
  end if
  if(occopt/=7)then
   cond_string(1)='occopt' ; cond_values(1)=occopt
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'prtstm',prtstm,1,(/0/),0,0,iout)
  end if
  if(nstep/=1)then
   cond_string(1)='nstep' ; cond_values(1)=nstep
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'prtstm',prtstm,1,(/0/),0,0,iout)
  end if
  if(ionmov/=0)then
   cond_string(1)='ionmov' ; cond_values(1)=ionmov
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'prtstm',prtstm,1,(/0/),0,0,iout)
  end if
  if(tolwfr<tol6)then
   cond_string(1)='tolwfr' ; cond_values(1)=0
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'prtstm',prtstm,1,(/0/),0,0,iout)
  end if
  if(prtden/=0)then
   cond_string(1)='prtden' ; cond_values(1)=prtden
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'prtstm',prtstm,1,(/0/),0,0,iout)
  end if
  if(prtnabla>0)then
   cond_string(1)='prtnabla' ; cond_values(1)=prtnabla
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'prtstm',prtstm,1,(/0/),0,0,iout)
  end if
  if(prtvxc>0)then
   cond_string(1)='prtvxc' ; cond_values(1)=prtvxc
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'prtstm',prtstm,1,(/0/),0,0,iout)
  end if
  if(prtvha>0)then
   cond_string(1)='prtvha' ; cond_values(1)=prtvha
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'prtstm',prtstm,1,(/0/),0,0,iout)
  end if
  if(prtvhxc>0)then
   cond_string(1)='prtvhxc' ; cond_values(1)=prtvhxc
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'prtstm',prtstm,1,(/0/),0,0,iout)
  end if

! prtwant
#if !defined HAVE_WANNIER90
  if(prtwant==2) then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' chkinp: ERROR -',ch10,&
&   '         prtwant==2 is only relevant if wannier90 library is linked',ch10,&
&   '         Action: check compilation options'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
#endif

! prtwf
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'prtwf',prtwf,3,(/0,1,2/),0,0,iout)

! ratsph
! If PAW and (prtdos==3 or prtdensph==1), must be larger than PAW radius
  if(usepaw==1.and.(prtdos==3.or.prtdensph==1))then
   do itypat=1,ntypat
    if (pspheads(itypat)%pawheader%rpaw>ratsph(itypat)) then
     write(message, '(10a,i2,a,f15.12,3a)' ) ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '  Projected DOS/density is required in the framework of PAW !',ch10,&
&     '  The radius of spheres in which DOS/density has to be projected',ch10,&
&     '  must be larger or equal than the (max.) PAW radius !',ch10,&
&     '  Rpaw(atom_type ',itypat,')= ',pspheads(itypat)%pawheader%rpaw,' au',ch10,&
&     '  Action : modify value of ratsph in input file.'
     call wrtout(iout,message,'COLL')
     call wrtout(06,  message,'COLL')
     ierr=ierr+1
    end if
   end do
  end if

! rfatpol
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'rfatpol(1)',rfatpol(1),1,(/1/),1,1,iout)
  cond_string(1)='natom' ; cond_values(1)=natom
  call chkint(1,1,cond_string,cond_values,ierr,&
&  'rfatpol(2)',rfatpol(2),1,(/natom/),-1,natom,iout)

! rprimd
! With optcell beyond 4, one has constraints on rprimd.
  if(optcell==4 .or. optcell==7 )then
   cond_string(1)='optcell' ; cond_values(1)=4
   if(optcell==7)cond_values(1)=7
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'rprimd(1,2)',rprimd(1,2),0,0.0_dp,iout)
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'rprimd(1,3)',rprimd(1,3),0,0.0_dp,iout)
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'rprimd(2,1)',rprimd(2,1),0,0.0_dp,iout)
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'rprimd(3,1)',rprimd(3,1),0,0.0_dp,iout)
  else if(optcell==5 .or. optcell==8 )then
   cond_string(1)='optcell' ; cond_values(1)=5
   if(optcell==8)cond_values(1)=8
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'rprimd(2,1)',rprimd(2,1),0,0.0_dp,iout)
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'rprimd(2,3)',rprimd(2,3),0,0.0_dp,iout)
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'rprimd(1,2)',rprimd(1,2),0,0.0_dp,iout)
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'rprimd(3,2)',rprimd(3,2),0,0.0_dp,iout)
  else if(optcell==6 .or. optcell==9 )then
   cond_string(1)='optcell' ; cond_values(1)=6
   if(optcell==9)cond_values(1)=9
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'rprimd(3,1)',rprimd(3,1),0,0.0_dp,iout)
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'rprimd(3,2)',rprimd(3,2),0,0.0_dp,iout)
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'rprimd(1,3)',rprimd(1,3),0,0.0_dp,iout)
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'rprimd(2,3)',rprimd(2,3),0,0.0_dp,iout)
  end if

! so_psp
  do ipsp=1,npsp
!  Check that so_psp is between 0 and 3
   if ( so_psp(ipsp)<0 .or. so_psp(ipsp)>3 ) then
    write(message, '(a,a,a,a,i3,a,i3,a,a,a,a,a)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '  so_psp(',ipsp,' ) was input as',&
&    so_psp(ipsp),' .',ch10,&
&    '  Input value must be 0, 1, 2, or 3.',ch10,&
&    ' Action : modify value of so_psp (old name : so_typat) in input file.'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    ierr=ierr+1
   end if
!  If nspinor=1, the spin-orbit contribution cannot be taken into account
   if ( nspinor==1 .and. (so_psp(ipsp)==2 .or. so_psp(ipsp)==3) ) then
    write(message, '(a,a,a,a,i2,a,i3,a,a,a,a,a)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '  so_psp(',ipsp,') was input as',&
&    so_psp(ipsp),', with nspinor=1.',ch10,&
&    '  When nspinor=1, so_psp cannot be required to be 2 or 3.',ch10,&
&    '  Action : modify value of so_psp (old name : so_typat) or nspinor in input file.'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    ierr=ierr+1
   end if
!
!  TO BE ACTIVATED IN V5.5
!  Send a warning if the spin-orbit contribution cannot be included due to nspinor=1
!  if ( nspinor==1 .and. so_psp(ipsp)==1 .and. pspheads(ipsp)%pspso/=1 ) then
!  write(message, '(a,a,a,a,i2,a,i3,a,a,a,a,a)' ) ch10,&
!  &    ' chkinp: WARNING -',ch10,&
!  &    '  so_psp(',ipsp,') was input as',&
!  &     so_psp(ipsp),', with nspinor=1 and a pseudopotential containing spin-orbit contribution.',ch10,&
!  &    '  However, when nspinor=1, existing psp spin-orbit contribution is not taken into account.'
!  call wrtout(iout,message,'COLL')
!  call wrtout(06,  message,'COLL')
!  end if
!  END OF TO BE ACTIVATED
!
!  If nspden=4, so_psp must be 1
!  if ( so_psp(ipsp)/=1 .and. nspden==4 ) then
!  write(message, '(a,a,a,a,i2,a,i3,7a)' ) ch10,&
!  &    ' chkinp: ERROR -',ch10,&
!  &    '  so_psp(',ipsp,') was input as',&
!  &     so_psp(ipsp),', with nspden=4.',ch10,&
!  &    '  However, non-collinear magnetism is not yet implemented with spin-orbit.',ch10,&
!  &    '  so_psp must be 1 for each pseudopotential type.',ch10,&
!  &    '  Action : modify value of so_psp or nspden in input file.'
!  call wrtout(iout,message,'COLL')
!  call wrtout(06,  message,'COLL')
!  ierr=ierr+1
!  end if
  end do

! stmbias
  cond_string(1)='prtstm' ; cond_values(1)=prtstm
  if(prtstm/=0)then
!  If non-zero prtstm, stmbias cannot be zero : test is positive or zero
   if(stmbias>-tol10)then
!   Then, enforce positive
    call chkdpr(1,1,cond_string,cond_values,ierr,'stmbias',stmbias,1,2*tol10,iout)
   end if
  else
   call chkdpr(1,1,cond_string,cond_values,ierr,'stmbias',stmbias,0,zero,iout)
  end if

! symafm
  if(nsppol==1 .and. nspden==2)then
!  At least one of the symmetry operations must be antiferromagnetic
   if(minval(symafm(1:nsym))/=-1)then
    write(message, '(8a)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '  When nsppol==1 and nspden==2, at least one of the symmetry operations',ch10,&
&    '  must be anti-ferromagnetic (symafm=-1), in order to deduce the spin-down density',ch10,&
&    '  from the spin-up density.'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    write(message, '(7a)' ) &
&    '  However, it is observed that none of the symmetry operations is anti-ferromagnetic.',ch10,&
&    '  Action : Check the atomic positions, the input variables spinat, symrel, tnons, symafm.',ch10,&
&    '           In case your system is not antiferromagnetic (it might be ferrimagnetic ...),',ch10,&
&    '           you must use nsppol=2 with nspden=2 (the latter being the default when nsppol=2).'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    call leave_new('COLL')
   end if
  end if

! symrel and tnons
! Check the point group closure (should check the spatial group closure !!)
  call chkgrp(nsym,symafm,symrel)
! Check the orthogonality of the symmetry operations
! (lengths and absolute values of scalar products should be preserved)
  call chkorthsy(gprimd,iout,nsym,rmet,rprimd,symrel)

! symchi
  if (symchi/=0.and.symchi/=1.and.symchi/=2) then
   write(message, '(a,a,a,a,i3,a,a,a,a)' ) ch10,&
&   ' chkinp: ERROR -',ch10,&
&   '  symchi  was input as ',symchi,ch10,&
&   '  Input value must be 0, 1, or 2.',ch10,&
&   ' Action : modify value of symchi in input file.'
   call wrtout(iout,message,'COLL')
   call wrtout(06,  message,'COLL')
   ierr=ierr+1
  end if

! symsigma
  if (symsigma/=0.and.symsigma/=1.and.symsigma/=2) then
   write(message, '(a,a,a,a,i3,a,a,a,a)' ) ch10,&
&   ' chkinp: ERROR -',ch10,&
&   '  symsigma  was input as',symsigma,ch10,&
&   '  Input value must be 0, 1, or 2.',ch10,&
&   ' Action : modify value of symsigma in input file.'
   call wrtout(iout,message,'COLL')
   call wrtout(06,  message,'COLL')
   ierr=ierr+1
  end if

! MG now it is possible to perform a GW calculation with non-symmorphic operations if required by the user
! tnons
  if (symmorphi==0) then
   if(nbandkss/=0)then
    do isym=1,nsym
     if(sum(tnons(:,isym)**2)>tol6)then
      write(message, '(6a,i3,a,3f8.4,3a)' ) ch10,&
&      ' chkinp: ERROR -',ch10,&
&      '  When nbandkss/=0, all the components of tnons must be zero.',ch10,&
&      '  However, for the symmetry operation number ',isym,', tnons =',tnons(:,isym),'.',ch10,&
&      '  Action : use the symmetry finder (nsym=0) with symmorphi==0.'
      call wrtout(iout,message,'COLL')
      call wrtout(06,  message,'COLL')
      call leave_new('COLL')
     end if
    end do
   end if
   if(optdriver==3 .or. optdriver==4)then
    do isym=1,nsym
     if(sum(tnons(:,isym)**2)>tol6)then
      write(message, '(6a,i3,a,3f8.4,3a)' ) ch10,&
&      ' chkinp: ERROR -',ch10,&
&      '  When optdriver==3 or 4, all the components of tnons must be zero.',ch10,&
&      '  However, for the symmetry operation number ',isym,', tnons =',tnons(:,isym),'.',ch10,&
&      '  Action : use the symmetry finder (nsym=0) with symmorphi==0.'
      call wrtout(iout,message,'COLL')
      call wrtout(06,  message,'COLL')
      call leave_new('COLL')
     end if
    end do
   end if
  end if !of symmorphi

! toldff
  if(optforces/=1)then
   cond_string(1)='optforces' ; cond_values(1)=optforces
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'toldff',toldff,1,zero,iout)
   cond_string(1)='optforces' ; cond_values(1)=optforces
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'tolrff',tolrff,1,zero,iout)
  end if

! tolwfr
  if(iscf<0)then
   cond_string(1)='iscf' ; cond_values(1)=iscf
   call chkdpr(1,1,cond_string,cond_values,ierr,&
&   'tolwfr',tolwfr,1,0.0_dp,iout)
  end if

! tsmear
! Check that tsmear is non-zero positive for metallic occupation functions
  if(3<=occopt .and. occopt<=7)then
   if(tsmear<1.0d-12)then
    write(message, '(8a,i3,a,es16.6)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '  When occopt corresponds to metallic occupations, tsmear must',ch10,&
&    '  be a non-zero, positive number.',ch10,&
&    '  In the input file, occopt=',occopt,', while tsmear=',tsmear
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    call leave_new('COLL')
   end if
  end if

! usedmatpu
  if (usepaw==1.and.usepawu==1) then
!  abs(usedmatpu)<=nstep
   cond_string(1)='nstep' ; cond_values(1)=nstep
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'abs(usedmatpu)',abs(usedmatpu),1,(/nstep/),-1,nstep,iout)
!  lpawu must be constant or -1
   if (usedmatpu/=0) then
    do itypat=1,ntypat
     if (lpawu(itypat)/=-1.and.lpawu(itypat)/=maxval(lpawu(:))) then
      write(message, '(6a)' ) ch10,&
&      ' chkinp: ERROR -',ch10,&
&      '  When usedmatpu/=0 (use of an initial density matrix for LDA+U),',ch10,&
&      '  lpawu must be equal for all types of atoms on which +U is applied !'
      call wrtout(iout,message,'COLL')
      call wrtout(06,  message,'COLL')
      call leave_new('COLL')
     end if
    end do
   end if
  end if

! useexexch and lexexch
! Local exact-exchange and restrictions
  if(useexexch/=0)then
   cond_string(1)='useexexch' ; cond_values(1)=useexexch
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'useexexch',useexexch,1,(/1/),0,0,iout)
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'usepaw',usepaw,1,(/1/),0,0,iout)
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'pawxcdev',pawxcdev,2,(/1,2/),0,0,iout)
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'ixc',ixc,1,(/11/),0,0,iout)
   do itypat=1,ntypat
    cond_string(1)='lexexch' ; cond_values(1)=lexexch(itypat)
    call chkint(1,1,cond_string,cond_values,ierr,&
&    'lexexch',lexexch(itypat),5,(/-1,0,1,2,3/),0,0,iout)
   end do
  end if

! usepawu and lpawu
! PAW+U and restrictions
  if(usepawu/=0)then
   cond_string(1)='usepawu' ; cond_values(1)=usepawu
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'usepawu',usepawu,2,(/1,2/),0,0,iout)
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'usepaw',usepaw,1,(/1/),0,0,iout)
   do itypat=1,ntypat
    cond_string(1)='lpawu' ; cond_values(1)=lpawu(itypat)
    call chkint(1,1,cond_string,cond_values,ierr,&
&    'lpawu',lpawu(itypat),5,(/-1,0,1,2,3/),0,0,iout)
   end do
  end if

! useexexch AND usepawu
! Restriction when use together
  if(useexexch>0.and.usepawu>0)then
   do itypat=1,ntypat
    if (lpawu(itypat)/=lexexch(itypat).and.&
&    lpawu(itypat)/=-1.and.lexexch(itypat)/=-1) then
     write(message, '(8a,i2,3a)' ) ch10,&
&     ' chkinp: ERROR - ',ch10,&
&     '  When PAW+U (usepawu>0) and local exact-exchange (useexexch>0)',ch10,&
&     '  are selected together, they must apply on the same',ch10,&
&     '  angular momentum (lpawu/=lexexch forbidden, here for typat=',&
&     itypat,') !',ch10,'  Action: correct your input file.'
     call wrtout(iout,message,'COLL')
     call wrtout(06,  message,'COLL')
     call leave_new('COLL')
    end if
   end do
  end if

! useylm
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'useylm',useylm,2,(/0,1/),0,0,iout)
  if (usepaw==1) then
   cond_string(1)='pspcod(atom_type)' ; cond_values(1)=7
   call chkint(1,1,cond_string,cond_values,ierr,&
&   'useylm',useylm,1,(/1/),0,0,iout)
  end if

! wfoptalg
! Must be larger or equal to 0
  call chkint(0,0,cond_string,cond_values,ierr,&
&  'wfoptalg',wfoptalg,1,(/0/),1,0,iout)
! wfoptalg==0,1,10,11,4,14,5 if PAW
  if (usepaw==1) then
   cond_string(1)='usepaw' ; cond_values(1)=1
   call chkint(0,1,cond_string,cond_values,ierr,&
&   'wfoptalg',wfoptalg,7,(/0,1,10,11,4,14,5/),0,0,iout)
  end if
  if (fftalg/=400 .and. fftalg/=401) then   ! If fftalg/=400, cannot use wfoptalg=4 or 14 (lopbcg algo)
   cond_string(1)='fftalg' ; cond_values(1)=fftalg
   call chkint(0,1,cond_string,cond_values,ierr,&
&   'wfoptalg',wfoptalg,6,(/0,1,2,3,10,11/),0,0,iout)
  end if

! wtk
! Check that no k point weight is < 0:
  do ikpt=1,nkpt
   if (wtk(ikpt)< -tiny(0.0_dp) ) then
    write(message, '(a,a,a,a,i5,a,1p,e12.4,a,a,a)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '  At k point number',ikpt,'  wtk=',wtk(ikpt),' <0.',ch10,&
&    '  Action : check wtk in input file. Each wtk must be >=0.'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    call leave_new('COLL')
   end if
  end do

! xred
! Check that two atoms are not on top of each other
  if(natom>1)then
   allocate(frac(3,natom))
   do ia=1,natom
!   Map reduced coordinate xred(mu,ia) into [0,1)
    frac(1,ia)=xred(1,ia)-aint(xred(1,ia))+0.5_dp-sign(0.5_dp,xred(1,ia))
    frac(2,ia)=xred(2,ia)-aint(xred(2,ia))+0.5_dp-sign(0.5_dp,xred(2,ia))
    frac(3,ia)=xred(3,ia)-aint(xred(3,ia))+0.5_dp-sign(0.5_dp,xred(3,ia))
   end do
   do ia=1,natom-1
    do ib=ia+1,natom
     if( abs(frac(1,ia)-frac(1,ib))<1.0d-6 .and. &
&     abs(frac(2,ia)-frac(2,ib))<1.0d-6 .and. &
&     abs(frac(3,ia)-frac(3,ib))<1.0d-6         ) then
      write(message, '(a,a,a,a,i4,a,i4,a,a,a,a,a,a)' ) ch10,&
&      ' chkinp: ERROR - ',ch10,&
&      '  Atoms number',ia,' and',ib,' are located at the same point',&
&      ' of the unit cell',ch10,&
&      '  (periodic images are taken into account).',ch10,&
&      '  Action: change the coordinate of one of these atoms in the input file.'
      call wrtout(iout,message,'COLL')
      call wrtout(06,  message,'COLL')
      call leave_new('COLL')
     end if
    end do
   end do
   deallocate(frac)
  end if

! znucl
! Check that znucl and znuclpsp agree
  do ipsp=1,npsp
   if (abs(znucl(ipsp)-pspheads(ipsp)%znuclpsp)> tol12 ) then
    write(message, '(4a,i5,a,es12.4,a,a,es12.4,2a)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '  For pseudopotential ',ipsp,'  znucl from user input file= ',znucl(ipsp),ch10,&
&    '  while znucl from pseudopotential file=',pspheads(ipsp)%znuclpsp,ch10,&
&    '  Action : check znucl in input file, or check psp file. They must agree.'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    call leave_new('COLL')
   end if
  end do

! bandFFT
  if(dtsets(idtset)%paral_kgb==1) then
   if (mod(wfoptalg,10) /= 4) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of wfoptalg is found to be ',wfoptalg,ch10,&
&    '  This is not allowed in the case of band-FFT parallelization.',ch10,&
&    '  Action : put wfoptalg = 4 or 14 in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if (nspinor /= 1) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of nspinor is found to be ',nspinor,ch10,&
&    '  This is not allowed in the case of band-FFT parallelization.',ch10,&
&    '  Action : put nspinor = 1 in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if (nloalg /= 0) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of nloalg is found to be ',nloalg,ch10,&
&    '  (k+G) vectors have to be precomputed in the case of band-FFT parallelization.',ch10,&
&    '  Action : put nloalg < 10 in your input file (try nloalg 4)'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if (occopt == 2) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of occopt is found to be ',occopt,ch10,&
&    '  The number of bands have to remain constant in the case of band-FFT parallelization.',ch10,&
&    '  Action : put occopt /= 2 in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if(maxval(abs(istwfk(:)-1)) > 1)then
    write(message,'(8a)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '  One of the components of istwfk is not equal to 1.',ch10,&
&    '  Time-reversal symmetry is not yet programmed in the case of band-FFT parallelization.',ch10,&
&    '  Action : set istwfk to 1 for all k-points'
    call wrtout(iout,message,'COLL')
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if (mkmem == 0) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of mkmem is found to be ',mkmem,ch10,&
&    '  An out-of-core solution can''t be used in the case of band-FFT parallelization.',ch10,&
&    '  Action : put mkmem = nkpt in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
  end if

! WVL - wavelets checks and limitations
  if(dtsets(idtset)%usewvl == 1) then
   if (dtsets(idtset)%wvl_hgrid <= 0) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of wvl_hgrid is found to be ',dtsets(idtset)%wvl_hgrid,ch10,&
&    '  This value is mandatory and must be positive.',ch10,&
&    '  Action : put wvl_hgrid to a positive value in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if (dtsets(idtset)%iscf /= 2) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of iscf is found to be ',dtsets(idtset)%iscf,ch10,&
&    '  Only simple potential mixing is allowed with wavelets.',ch10,&
&    '  Action : put iscf = 2 in your input file.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if(dtsets(idtset)%icoulomb /= 1)then
    write(message,'(a,a,a,a,i3,a,a,a,a)' ) ch10,&
&    ' chkinp: ERROR -',ch10,&
&    '  The value of icoulomb is found to be ',dtsets(idtset)%icoulomb,ch10,&
&    '  The real space computation of hartree potential is mandatory with wavelets.',ch10,&
&    '  Action : put icoulomb = 1 in your input file.'
    call wrtout(iout,message,'COLL')
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if (dtsets(idtset)%nsym /= 1) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of nsym is found to be ',dtsets(idtset)%nsym,ch10,&
&    '  No symetry operations are allowed for isolated systems.',ch10,&
&    '  Action : put nsym = 1 in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if (dtsets(idtset)%optstress > 0) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of optstress is found to be ', dtsets(idtset)%optstress, ch10,&
&    '  There is no stress computation available with the wavelet code.',ch10,&
&    '  Action : put optstress = 0 in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if (usepaw == 1) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of usepaw is found to be ',usepaw,ch10,&
&    '  The wavelet computation is not allowed in the framework of PAW.',ch10,&
&    '  Action : use HGH or GTH pseudo-potentials'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if (dtsets(idtset)%nspden > 2) then
    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of nspden is found to be ', dtsets(idtset)%nspden, ch10, &
&    '  The wavelet computation is not allowed with non-colinear spin.',ch10,&
&    '  Action : put nspden = 1 or 2 in your input file'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if (dtsets(idtset)%nspden /= dtsets(idtset)%nsppol) then
    write(message,'(a,a,a,a,i3,a,a,i3,a,a)')ch10,&
&    ' chkinp : ERROR -',ch10,&
&    '  The value of nspden is found to be ', dtsets(idtset)%nspden, ch10, &
&    '  and the one of nsppol is found to be ', dtsets(idtset)%nsppol, ch10, &
&    '  In wavelet computation, nspden and nsppol must be equal.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
!  We check the consistency of occupation, empty bands are not allowed.
   if (dtsets(idtset)%nsppol == 2) then
    mband = dtsets(idtset)%nelect
   else
    mband = dtsets(idtset)%mband
   end if
   do ii = 1, mband, 1
    if (dtsets(idtset)%occ_orig(ii) < tol8) then
     write(message,'(a,a,a,a,f6.4,a,a,a,a,a,a)') ch10,&
&     ' chkinp : ERROR -',ch10,&
&     '  One value of occ is found to be ', dtsets(idtset)%occ_orig(ii), ch10, &
&     '  The wavelet computation is not allowed empty bands.',ch10,&
&     '  Action : use occopt = 1 for automatic band filling or', ch10, &
&     '  change occ value in your input file'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')
    end if
   end do
  end if


! If molecular dynamics or structural optimization is being done
! (ionmov>0), make sure not all atoms are fixed
! if (ionmov > 0) then
! if (natfix == natom) then
! write(message, '(a,a,a,a,i4,a,i5,a,a,i5,a,a,a,a,a,a)' ) ch10,&
! &   ' setup1: ERROR -',ch10,&
! &   '  ionmov is ',ionmov,' and number of fixed atoms is ',natfix,ch10,&
! &   '  while number of atoms natom is ',natom,'.',ch10,&
! &   '  Thus all atoms are fixed and option ionmov to move atoms',&
! &           ' is inconsistent.',ch10,&
! &   '  Action : change ionmov or natfix and iatfix in input file and resubmit.'
! call wrtout(06,message,'COLL')
! call leave_new('COLL')
! end if
! end if

! Should check that the symmetry operations are consistent with iatfixx,
! iatfixy and iatfixz (diagonal symmetry operations)

! Should check values of fftalg

! Should check values of nloalg

! rfasr=2 possible only when electric field response is computed.

! Must have nqpt=1 for rfphon=1

! ** Here ends the checking section **************************************

  deallocate(iatsph,istwfk,kpt,mixalch,shiftk)
  deallocate(nband,occ,ratsph,so_psp,symafm,symrel,tnons,typat,wtk,xred,znucl)
  deallocate(upawu,jpawu,lpawu,lexexch)

! End do loop on idtset (allocate statements are present)
 end do

 if (maxval(dtsets(:)%usewvl) > 0) then
  write(message,'(5A)') ch10,&
&  ' chkinp: COMMENT - comparison between wvl_hgrid and ecut',ch10,&
&  '  real-space mesh | eq. Ec around atoms | eq. Ec further from atoms'
  call wrtout(6,message,'COLL')
  wvl_hgrid = zero
  twvl = .false.
  do idtset=1,ndtset_alloc
!  Give an indication to the equivalent ecut corresponding to
!  given hgrid.
   if (dtsets(idtset)%usewvl == 1 .and. &
&   wvl_hgrid /= dtsets(idtset)%wvl_hgrid) then
    write(message,'(F11.3,A,F16.1,A,F16.1,A)') &
&    dtsets(idtset)%wvl_hgrid, " bohr  |", &
&    two * pi * pi / (dtsets(idtset)%wvl_hgrid ** 2), " Ht  | ", &
&    half * pi * pi / (dtsets(idtset)%wvl_hgrid ** 2), " Ht"
    call wrtout(6,message,'COLL')
    wvl_hgrid = dtsets(idtset)%wvl_hgrid
   end if
   twvl = twvl .or. (dtsets(idtset)%usewvl == 1 .and. &
&   dtsets(idtset)%accesswff /= 3)
  end do
  if (twvl) then
   write(message,'(8a)') ch10,&
&   ' chkinp: WARNING -',ch10,&
&   '  Restart files from wavelets in non ETSF format does not follow', ch10, &
&   '  the ABINIT standards.', ch10, &
&   '  Put accesswff to 3 to use ETSF retart files.'
   call wrtout(6,message,'COLL')
  end if
 end if

!If there was a problem, then stop.
 if(ierr/=0)then
  call leave_new('COLL')
 end if

!DEBUG
!write(6,*)' chkinp : exit '
!stop
!ENDDEBUG

end subroutine chkinp
!!***
