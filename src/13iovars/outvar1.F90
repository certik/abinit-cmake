!{\src2tex{textfont=tt}}
!!****f* ABINIT/outvar1
!! NAME
!! outvar1
!!
!! FUNCTION
!! Echo variables between acell and natom (by alphabetic order)
!! for the ABINIT code.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  choice= 1 if echo of preprocessed variables, 2 if echo after call driver
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  nqptdm=the number of q vectors provided by the user to calculate DM in GW
!!  iout=unit number for echoed output
!!  istatr=repetition rate for status file
!!  istatshft=shift of the repetition rate for status file
!!  jdtset_(0:ndtset_alloc)=actual index of the dataset (equal to dtsets(:)%jdtset)
!!  mxmband=maximum number of bands
!!  mxnatom=maximal value of input natom for all the datasets
!!  mxnatsph=maximal value of input natsph for all the datasets
!!  mxnatvshift=maximal value of input natsph for all the datasets
!!  mxnkptgw=maximal value of input nkptgw for all the datasets
!!  mxnkpt=maximal value of input nkpt for all the datasets
!!  mxnqptdm=maximal value of input nqptdm for all the datasets
!!  mxnsppol=maximal value of input nsppol for all the datasets
!!  mxnsym=maximum number of symmetries
!!  mxntypat=maximum number of type of atoms
!!  ndtset=number of datasets
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!      one data set. Use for most dimensioned arrays.
!!  npsp=number of pseudopotentials
!!  prtvol_glob= if 0, minimal output volume, if 1, no restriction.
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!  response= 1 if response variables must be output, 0 otherwise.
!!  response_(0:ndtset_alloc)= 1 if response variables must be output, 0 otherwise,
!!   for different datasets
!!  results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including evolving variables
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Note that this routine is called only by the processor me==0 .
!! In consequence, no use of message and wrtout routine.
!! The lines of code needed to output the defaults are preserved
!! (see last section of the routine, but are presently disabled)
!!
!!  Note that acell, occ, rprim, xred and vel might have been modified by the
!!  computation, so that their values if choice=1 or choice=2 will differ.
!!
!! PARENTS
!!      outvars
!!
!! CHILDREN
!!      prttagm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine outvar1 (choice,dtsets,iout,istatr,istatshft,&
& jdtset_,mxmband,mxnatom,mxnatsph,mxnatvshift,mxnkptgw,mxnkpt,mxnqptdm,mxnsppol,mxnsym,mxntypat,&
& ndtset,ndtset_alloc,npsp,prtvol_glob,pspheads,response,response_,results_out,usepaw)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_13iovars, except_this_one => outvar1
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,iout,istatr,istatshft,mxmband,mxnatom,mxnatsph,mxnatvshift
 integer,intent(in) :: mxnkpt,mxnkptgw,mxnqptdm,mxnsppol,mxnsym,mxntypat,ndtset
 integer,intent(in) :: ndtset_alloc,npsp,prtvol_glob,response,usepaw
!arrays
 integer,intent(in) :: jdtset_(0:ndtset_alloc),response_(ndtset_alloc)
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(pspheader_type),intent(in) :: pspheads(npsp)
 type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)

!Local variables-------------------------------
 character(len=*), parameter :: format01110 ="(1x,a9,1x,(t13,8i8) )"
 character(len=*), parameter :: format01150a="(1x,a9,a,1x,(t13,3es16.8))"
 character(len=*), parameter :: format01155 ="(1x,a9,1x,(t13,10i5))"
 character(len=*), parameter :: format01155a="(1x,a9,a,1x,(t13,10i5))"
 character(len=*), parameter :: format01160 ="(1x,a9,1x,(t13,3es18.10)) "
 character(len=*), parameter :: format01160a="(1x,a9,a,1x,(t13,3es18.10)) "
 character(len=*), parameter :: format01170 ="(1x,a9,a,1x,(t13,5f11.6)) "
!scalars
 integer,parameter :: nkpt_max=50
 integer :: allowed,first,iatom,iban,idtset,ii,ikpt,iscf,istatr_defo
 integer :: istatshft_defo,jdtset,kptopt,marr,mu
 integer :: multi_natom,multi_natfix,multi_natfixx
 integer :: multi_natfixy,multi_natfixz,multi_natsph,multi_natvshift,multi_nberry,multi_nkpt
 integer :: multi_nkptgw,multi_norb,multi_nqptdm,multi_nshiftk,multi_ntypalch
 integer :: multi_nsppol
 integer :: multi_ntypat,natfix,natfixx,natfixy,natfixz,natom,natsph,natvshift,nban
 integer :: nberry,ndtset_kptopt,nkpt,nkpt_eff,nkptgw,norb,npspalch,nqpt,nqptdm
 integer :: nshiftk,nsppol,nsym,ntypalch,ntypat,occopt,tnkpt
 real(dp) :: cpus,kpoint
 character(len=2) :: appen
 character(len=500) :: message
!arrays
 integer,allocatable :: iatfixio_(:,:),iatfixx_(:,:),iatfixy_(:,:)
 integer,allocatable :: iatfixz_(:,:),intarr(:,:),istwfk_2(:,:)
 integer,allocatable :: jdtset_kptopt(:),natfix_(:),natfixx_(:),natfixy_(:)
 integer,allocatable :: natfixz_(:)
 real(dp) :: acell(3)
 real(dp),allocatable :: dprarr(:,:)

! *************************************************************************

!DEBUG
!write(6,*)' outvar1 : enter '
!ENDDEBUG

!Must treat separately the translation of iatfix from the internal
!representation to the input/output representation
 allocate(natfix_(0:ndtset_alloc),iatfixio_(mxnatom,0:ndtset_alloc))
 allocate(natfixx_(0:ndtset_alloc),iatfixx_(mxnatom,0:ndtset_alloc))
 allocate(natfixy_(0:ndtset_alloc),iatfixy_(mxnatom,0:ndtset_alloc))
 allocate(natfixz_(0:ndtset_alloc),iatfixz_(mxnatom,0:ndtset_alloc))
 natfix_(0:ndtset_alloc)=0 ; iatfixio_(:,0:ndtset_alloc)=0
 natfixx_(0:ndtset_alloc)=0 ; iatfixx_(:,0:ndtset_alloc)=0
 natfixy_(0:ndtset_alloc)=0 ; iatfixy_(:,0:ndtset_alloc)=0
 natfixz_(0:ndtset_alloc)=0 ; iatfixz_(:,0:ndtset_alloc)=0
 do idtset=1,ndtset_alloc
! DEBUG
! write(6,*)' outvar1 : iatfix_ for idtset= ',idtset
! ENDDEBUG
  do iatom=1,dtsets(idtset)%natom
!  First look whether the atom is fixed along the three directions
   if( dtsets(idtset)%iatfix(1,iatom)+ &
&   dtsets(idtset)%iatfix(2,iatom)+ &
&   dtsets(idtset)%iatfix(3,iatom)   ==3 )then
    natfix_(idtset)=natfix_(idtset)+1
!   DEBUG
!   write(6,*)' outvar1: iatom,natfix_(idtset)=',iatom,natfix_(idtset)
!   ENDDEBUG
    iatfixio_(natfix_(idtset),idtset)=iatom
   else
!   Now examine each direction, one at a time
    if( dtsets(idtset)%iatfix(1,iatom) ==1)then
     natfixx_(idtset)=natfixx_(idtset)+1
     iatfixx_(natfixx_(idtset),idtset)=iatom
    end if
    if( dtsets(idtset)%iatfix(2,iatom) ==1)then
     natfixy_(idtset)=natfixy_(idtset)+1
     iatfixy_(natfixy_(idtset),idtset)=iatom
    end if
    if( dtsets(idtset)%iatfix(3,iatom) ==1)then
     natfixz_(idtset)=natfixz_(idtset)+1
     iatfixz_(natfixz_(idtset),idtset)=iatom
    end if
   end if
  end do
! DEBUG
! write(6,*)' natfix ...'
! write(6,*)natfix_(idtset),natfixx_(idtset),natfixy_(idtset),natfixz_(idtset)
! ENDDEBUG
 end do

!Maximal size of dprarr and intarr arrays
 marr=max(3*mxnatom,3*mxnkptgw,mxnkpt*mxnsppol*mxmband,3*mxnkpt,npsp,mxntypat,&
& 9*mxnsym,mxnatsph,mxnatvshift*mxnsppol*mxnatom)
 allocate(dprarr(marr,0:ndtset_alloc))
 allocate(intarr(marr,0:ndtset_alloc))

!Set up dimensions : determine whether these are different for different
!datasets.

 multi_nqptdm=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%nqptdm/=dtsets(idtset)%nqptdm)multi_nqptdm=1
  end do
 end if
 if(multi_nqptdm==0)nqptdm=dtsets(1)%nqptdm

 multi_natfix=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(natfix_(1)/=natfix_(idtset))multi_natfix=1
  end do
 end if
 if(multi_natfix==0)natfix=natfix_(1)

 multi_natfixx=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(natfixx_(1)/=natfixx_(idtset))multi_natfixx=1
  end do
 end if
 if(multi_natfixx==0)natfixx=natfixx_(1)

 multi_natfixy=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(natfixy_(1)/=natfixy_(idtset))multi_natfixy=1
  end do
 end if
 if(multi_natfixy==0)natfixy=natfixy_(1)

 multi_natfixz=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(natfixz_(1)/=natfixz_(idtset))multi_natfixz=1
  end do
 end if
 if(multi_natfixz==0)natfixz=natfixz_(1)

 multi_natom=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%natom/=dtsets(idtset)%natom)multi_natom=1
  end do
 end if
 if(multi_natom==0)natom=dtsets(1)%natom

 multi_natsph=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%natsph/=dtsets(idtset)%natsph)multi_natsph=1
  end do
 end if
 if(multi_natsph==0)natsph=dtsets(1)%natsph

 multi_natvshift=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%natvshift/=dtsets(idtset)%natvshift)multi_natvshift=1
  end do
 end if
 if(multi_natvshift==0)natvshift=dtsets(1)%natvshift


 multi_nberry=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%nberry/=dtsets(idtset)%nberry)multi_nberry=1
  end do
 end if
 if(multi_nberry==0)nberry=dtsets(1)%nberry

 multi_nkptgw=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%nkptgw/=dtsets(idtset)%nkptgw)multi_nkptgw=1
  end do
 end if
 if(multi_nkptgw==0)nkptgw=dtsets(1)%nkptgw

 multi_nkpt=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%nkpt/=dtsets(idtset)%nkpt)multi_nkpt=1
  end do
 end if
 if(multi_nkpt==0)nkpt=dtsets(1)%nkpt

 multi_norb=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%norb/=dtsets(idtset)%norb)multi_norb=1
  end do
 end if
 if(multi_norb==0)norb=dtsets(1)%norb

 multi_nshiftk=0
 if(sum((dtsets(1:ndtset_alloc)%kptopt)**2)/=0)then
  first=0
  do idtset=1,ndtset_alloc
   kptopt=dtsets(idtset)%kptopt
   if(kptopt>=1)then
    if(first==0)then
     first=1
     nshiftk=dtsets(idtset)%nshiftk
    else
     if(nshiftk/=dtsets(idtset)%nshiftk)multi_nshiftk=1
    end if
   end if
  end do
 end if

 multi_nsppol=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%nsppol/=dtsets(idtset)%nsppol)multi_nsppol=1
  end do
 end if
 if(multi_nsppol==0)nsppol=dtsets(1)%nsppol

 multi_ntypalch=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%ntypalch/=dtsets(idtset)%ntypalch)multi_ntypalch=1
  end do
 end if
 if(multi_ntypalch==0)ntypalch=dtsets(1)%ntypalch

 multi_ntypat=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%ntypat/=dtsets(idtset)%ntypat)multi_ntypat=1
  end do
 end if
 if(multi_ntypat==0)ntypat=dtsets(1)%ntypat

!Print each variable, one at a time

 intarr(1,:)=dtsets(:)%accesswff
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'accesswff','INT')

 do idtset=0,ndtset_alloc
  dprarr(1:3,idtset)=results_out(idtset)%acell(:)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,ndtset_alloc,'acell','LEN')

!DEBUG
!write(6,*)' outvar1 : before algalch'
!ENDDEBUG

!algalch
 if(multi_ntypalch==0)then
  do idtset=0,ndtset_alloc
   intarr(1:ntypalch,idtset)=dtsets(idtset)%algalch(1:ntypalch)
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypalch,ndtset_alloc,'algalch','INT')
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01155a)'algalch',appen,dtsets(idtset)%algalch(1:dtsets(idtset)%ntypalch)
  end do
 end if

!atvshift
 if(usepaw>0)then
  if(multi_natvshift==0 .and. multi_nsppol==0 .and. multi_natom==0)then
   if(natvshift/=0)then
    do idtset=0,ndtset_alloc
     dprarr(1:natvshift*nsppol*natom,idtset)=&
&     reshape(dtsets(idtset)%atvshift(1:natvshift,1:nsppol,1:natom),(/natvshift*nsppol*natom/) )
    end do
    call prttagm(dprarr,intarr,iout,jdtset_,-5,marr,natvshift*nsppol*natom,&
&    ndtset_alloc,'atvshift','DPR')
   end if
  else
   do idtset=1,ndtset_alloc
    if(dtsets(idtset)%natvshift/=0)then
     jdtset=jdtset_(idtset)
     if(jdtset<10)write(appen,'(i1)')jdtset
     if(jdtset>=10)write(appen,'(i2)')jdtset
     write(iout,format01170)'atvshift',appen,&
&     dtsets(idtset)%atvshift(1:dtsets(idtset)%natvshift,1:dtsets(idtset)%nsppol,1:dtsets(idtset)%natom)
    end if
   end do
  end if
 end if

 dprarr(1,:)=dtsets(:)%alpha
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'alpha','DPR')

!amu
 if(multi_ntypat==0)then
  do idtset=0,ndtset_alloc
   dprarr(1:ntypat,idtset)=dtsets(idtset)%amu(1:ntypat)
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,ndtset_alloc,'amu','DPR')
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01160a)'amu',appen,dtsets(idtset)%amu(1:dtsets(idtset)%ntypat)
  end do
 end if

 intarr(1,:)=dtsets(:)%awtr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'awtr','INT')

 intarr(1,:)=dtsets(:)%berryopt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'berryopt','INT')

 intarr(1,:)=dtsets(:)%bdberry(1)
 intarr(2,:)=dtsets(:)%bdberry(2)
 intarr(3,:)=dtsets(:)%bdberry(3)
 intarr(4,:)=dtsets(:)%bdberry(4)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,4,ndtset_alloc,'bdberry','INT')

!bdgw
 if(multi_nkptgw==0)then
! Required if for pathscale  to avoid failure with -ff_bounds_check
  if (nkptgw > 0) then
   do idtset=0,ndtset_alloc
    intarr(1:2*nkptgw,idtset)=&
&    reshape(dtsets(idtset)%bdgw(1:2,1:nkptgw),(/2*nkptgw/) )
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2*nkptgw,&
&   ndtset_alloc,'bdgw','INT')
  end if
 else
  do idtset=1,ndtset_alloc
   if(dtsets(idtset)%nkptgw>0)then
    jdtset=jdtset_(idtset)
    if(jdtset<10)write(appen,'(i1)')jdtset
    if(jdtset>=10)write(appen,'(i2)')jdtset
    write(iout,format01155a)'bdgw',appen,dtsets(idtset)%bdgw(1:2,1:dtsets(idtset)%nkptgw)
   end if
  end do
 end if


 dprarr(1,:)=dtsets(:)%boxcenter(1)
 dprarr(2,:)=dtsets(:)%boxcenter(2)
 dprarr(3,:)=dtsets(:)%boxcenter(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,ndtset_alloc,'boxcenter','DPR')

 dprarr(1,:)=dtsets(:)%boxcutmin
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'boxcutmin','DPR')

 if (usepaw==1) then

  dprarr(1,:)=dtsets(:)%bxctmindg
  call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'bxctmindg','DPR')

 end if

 intarr(1,:)=dtsets(:)%ceksph
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ceksph','INT')

 intarr(1,:)=dtsets(:)%chkexit
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'chkexit','INT')

 dprarr(1,:)=dtsets(:)%charge
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'charge','DPR')

 if(dtsets(1)%cpus>one)then
  cpus=dtsets(1)%cpus
  write(iout,'(1x,a9,1x,1p,t13,g10.2,t25,a)') 'cpus',cpus,'(seconds)'
  write(iout,'(1x,a9,1x,1p,t13,g10.2,t25,a)') 'cpum',cpus/60.0_dp,'(minutes)'
  write(iout,'(1x,a9,1x,1p,t13,g10.2,t25,a)') 'cpuh',cpus/3600.0_dp,'(hours)'
 end if

 do idtset=0, ndtset_alloc
  do ii = 1, ntypat
   dprarr(ii,idtset) = dtsets(idtset)%corecs(ii)
  end do ! end loop over ntypat
 end do ! end loop over datasets
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,ndtset_alloc,'corecs','DPR')

 dprarr(1,:)=dtsets(:)%dedlnn
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'dedlnn','ENE')

 intarr(1,:)=dtsets(:)%delayperm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'delayperm','INT')

!densty
 if(multi_ntypat==0)then
  do idtset=0,ndtset_alloc
!  Only one component of densty is used until now
   dprarr(1:ntypat,idtset)=dtsets(idtset)%densty(1:ntypat,1)
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,ndtset_alloc,'densty','DPR')
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01160a)'densty',appen,dtsets(idtset)%densty(1:dtsets(idtset)%ntypat,1)
  end do
 end if

 dprarr(1,:)=dtsets(:)%diecut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'diecut','ENE')

 dprarr(1,:)=dtsets(:)%diegap
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'diegap','ENE')

 dprarr(1,:)=dtsets(:)%dielam
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'dielam','DPR')

 dprarr(1,:)=dtsets(:)%dielng
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'dielng','LEN')

 dprarr(1,:)=dtsets(:)%diemac
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'diemac','DPR')

 dprarr(1,:)=dtsets(:)%diemix
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'diemix','DPR')

 dprarr(1,:)=dtsets(:)%dilatmx
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'dilatmx','DPR')

 dprarr(1,:)=dtsets(:)%dosdeltae
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'dosdeltae','ENE')

 do idtset=0,ndtset_alloc
  intarr(1:3,idtset)=dtsets(idtset)%dsifkpt(1:3)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,ndtset_alloc,'dsifkpt','INT')

 dprarr(1,:)=dtsets(:)%dtion
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'dtion','DPR')

 dprarr(1,:)=dtsets(:)%ecut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'ecut','ENE')

 dprarr(1,:)=dtsets(:)%ecuteps
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'ecuteps','ENE')

 dprarr(1,:)=dtsets(:)%ecutsigx
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'ecutsigx','ENE')

 dprarr(1,:)=dtsets(:)%ecutwfn
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'ecutwfn','ENE')

 dprarr(1,:)=dtsets(:)%ecutsm
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'ecutsm','ENE')

 dprarr(1,:)=dtsets(:)%effmass
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'effmass','DPR')

 dprarr(1,:)=dtsets(:)%efield(1)
 dprarr(2,:)=dtsets(:)%efield(2)
 dprarr(3,:)=dtsets(:)%efield(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,ndtset_alloc,'efield','DPR')

 intarr(1,:)=dtsets(:)%enunit
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'enunit','INT')

 dprarr(1,:)=dtsets(:)%eshift
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'eshift','ENE')

 dprarr(1,:)=dtsets(:)%exchmix
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'exchmix','DPR')
!etotal
 if(choice==2)then
  do idtset=1,ndtset_alloc
   iscf=dtsets(idtset)%iscf
   if(iscf>0 .or. iscf==-3)then
    if(ndtset>0)then
     jdtset=jdtset_(idtset)
     if(jdtset<10)write(appen,'(i1)')jdtset
     if(jdtset>=10)write(appen,'(i2)')jdtset
     write(iout,format01160a)'etotal',appen,results_out(idtset)%etotal
    else
     write(iout,format01160)'etotal',results_out(idtset)%etotal
    end if
   end if
  end do
 end if

 intarr(1,:)=dtsets(:)%exchn2n3d
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'exchn2n3d','INT')

 intarr(1,:)=dtsets(:)%ngfft(7)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'fftalg','INT')

 intarr(1,:)=dtsets(:)%ngfft(8)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'fftcache','INT')

 intarr(1,:)=dtsets(:)%fftgw
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'fftgw','INT')

 intarr(1,:)=dtsets(:)%fft_opt_lob
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'fft_opt_lob','INT')

!force
 if(choice==2)then
  do idtset=1,ndtset_alloc
   iscf=dtsets(idtset)%iscf
   if(iscf>0)then
    if(ndtset>0)then
     jdtset=jdtset_(idtset)
     if(jdtset<10)write(appen,'(i1)')jdtset
     if(jdtset>=10)write(appen,'(i2)')jdtset
     write(iout,format01160a)'fcart',appen,&
&     results_out(idtset)%fcart(:,1:dtsets(idtset)%natom)
    else
     write(iout,format01160)'fcart',&
&     results_out(idtset)%fcart(:,1:dtsets(idtset)%natom)
    end if
   end if
  end do
 end if

 dprarr(1,:)=dtsets(:)%fixmom
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'fixmom','DPR')

 dprarr(1,:)=dtsets(:)%rhoqpmix
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'rhoqpmix','DPR')

 dprarr(1,:)=dtsets(:)%freqremax
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'freqremax','ENE')

 dprarr(1,:)=dtsets(:)%freqspmax
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'freqspmax','ENE')

 dprarr(1,:)=dtsets(:)%freqsusin
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'freqsusin','DPR')

 dprarr(1,:)=dtsets(:)%freqsuslo
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'freqsuslo','DPR')

 dprarr(1,:)=dtsets(:)%friction
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'friction','DPR')

 intarr(1,:)=dtsets(:)%frzfermi
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'frzfermi','INT')

 intarr(1,:)=dtsets(:)%getcell
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getcell','INT')

 intarr(1,:)=dtsets(:)%getddk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getddk','INT')

 intarr(1,:)=dtsets(:)%getden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getden','INT')

 intarr(1,:)=dtsets(:)%getqps
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getqps','INT')

 intarr(1,:)=dtsets(:)%getscr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getscr','INT')

 intarr(1,:)=dtsets(:)%getsuscep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getsuscep','INT')

 intarr(1,:)=dtsets(:)%getkss
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getkss','INT')

 intarr(1,:)=dtsets(:)%getocc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getocc','INT')

 intarr(1,:)=dtsets(:)%getvel
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getvel','INT')

 intarr(1,:)=dtsets(:)%getwfk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getwfk','INT')

 intarr(1,:)=dtsets(:)%getwfq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getwfq','INT')

 intarr(1,:)=dtsets(:)%getxcart
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getxcart','INT')

 intarr(1,:)=dtsets(:)%getxred
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getxred','INT')

 intarr(1,:)=dtsets(:)%get1den
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'get1den','INT')

 intarr(1,:)=dtsets(:)%get1wf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'get1wf','INT')

 intarr(1,:)=dtsets(:)%gwcalctyp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gwcalctyp','INT')

!iatfix
 if(multi_natfix==0)then
  if(natfix/=0)then
   intarr(1:natfix,0:ndtset_alloc)=iatfixio_(1:natfix,0:ndtset_alloc)
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,natfix,&
&   ndtset_alloc,'iatfix','INT')
  end if
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01155a)'iatfix',appen,iatfixio_(1:natfix_(idtset),idtset)
  end do
 end if

!iatfixx
 if(multi_natfixx==0)then
  if(natfixx/=0)then
   intarr(1:natfixx,0:ndtset_alloc)=iatfixx_(1:natfixx,0:ndtset_alloc)
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,natfixx,&
&   ndtset_alloc,'iatfixx','INT')
  end if
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01155a)'iatfix',appen,iatfixx_(1:natfixx_(idtset),idtset)
  end do
 end if

!iatfixy
 if(multi_natfixy==0)then
  if(natfixy/=0)then
   intarr(1:natfixy,0:ndtset_alloc)=iatfixy_(1:natfixy,0:ndtset_alloc)
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,natfixy,&
&   ndtset_alloc,'iatfixy','INT')
  end if
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01155a)'iatfix',appen,iatfixy_(1:natfixy_(idtset),idtset)
  end do
 end if

!iatfixz
 if(multi_natfixz==0)then
  if(natfixz/=0)then
   intarr(1:natfixz,0:ndtset_alloc)=iatfixz_(1:natfixz,0:ndtset_alloc)
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,natfixz,&
&   ndtset_alloc,'iatfixz','INT')
  end if
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01155a)'iatfix',appen,iatfixz_(1:natfixz_(idtset),idtset)
  end do
 end if

!iatsph   need to be printed only if there is some occurence of prtdos==3
 do idtset=1,ndtset_alloc
  if(dtsets(idtset)%prtdos==3)then
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01155a)'iatsph',appen,dtsets(idtset)%iatsph(1:dtsets(idtset)%natsph)
  end if
 end do

 intarr(1,:)=dtsets(:)%iboxcut
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'iboxcut','INT')

 intarr(1,:)=dtsets(:)%icutcoul
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'icutcoul','INT')

 intarr(1,:)=dtsets(:)%icoulomb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'icoulomb','INT')

 intarr(1,:)=dtsets(:)%idyson
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'idyson','INT')

 intarr(1,:)=dtsets(:)%ieig2rf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ieig2rf','INT')

 intarr(1,:)=dtsets(:)%ikhxc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ikhxc','INT')

 intarr(1,:)=dtsets(:)%inclvkb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'inclvkb','INT')

 intarr(1,:)=dtsets(:)%intxc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'intxc','INT')

 intarr(1,:)=dtsets(:)%intexact
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'intexact','INT')

 intarr(1,:)=dtsets(:)%ionmov
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ionmov','INT')

 intarr(1,:)=dtsets(:)%iprcch
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'iprcch','INT')

 intarr(1,:)=dtsets(:)%iprcel
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'iprcel','INT')

 intarr(1,:)=dtsets(:)%iprctfvw
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'iprctfvw','INT')

 intarr(1,:)=dtsets(:)%iprcfc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'iprcfc','INT')

 intarr(1,:)=dtsets(:)%irdddk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdddk','INT')

 intarr(1,:)=dtsets(:)%irdkss
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdkss','INT')

 intarr(1,:)=dtsets(:)%irdqps
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdqps','INT')

 intarr(1,:)=dtsets(:)%irdscr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdscr','INT')

 intarr(1,:)=dtsets(:)%irdsuscep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdsuscep','INT')

 intarr(1,:)=dtsets(:)%irdwfk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdwfk','INT')

 intarr(1,:)=dtsets(:)%irdwfq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdwfq','INT')

 intarr(1,:)=dtsets(:)%ird1wf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ird1wf','INT')

 intarr(1,:)=dtsets(:)%iscf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'iscf','INT')

 intarr(1,:)=dtsets(:)%isecur
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'isecur','INT')

 istatr_defo=49
#if defined T3E
!The T3E has slow IO, so the rate for status files must be small by default
 istatr_defo=149
#endif
 if(istatr/=istatr_defo)write(iout,format01110) 'istatr',istatr

 istatshft_defo=1
 if(istatshft/=istatshft_defo)write(iout,format01110) 'istatshft',istatshft

!istwfk (must first restore the default istwf=0 for non-allowed k points)
 allocate(istwfk_2(mxnkpt,0:ndtset_alloc))
 do idtset=1,ndtset_alloc
  nqpt=dtsets(idtset)%nqpt
  do ikpt=1,dtsets(idtset)%nkpt
   allowed=1
   do ii=1,3
!   kpoint=dtsets(idtset)%kptns(ii)
    kpoint=dtsets(idtset)%kpt(ii,ikpt)/dtsets(idtset)%kptnrm
    if(nqpt/=0 .and. response_(idtset)==0)&
&    kpoint=kpoint+dtsets(idtset)%qptn(ii)
    if(abs(kpoint)>1.d-10 .and. abs(kpoint-0.5_dp)>1.e-10_dp )&
&    allowed=0
   end do
   if(allowed==0)then
    istwfk_2(ikpt,idtset)=0
   else
    istwfk_2(ikpt,idtset)=dtsets(idtset)%istwfk(ikpt)
   end if
  end do
 end do

!DEBUG
!write(6,*)' outvar1 '
!write(6,*)istwfk_2(:,:)
!ENDDEBUG

 if(multi_nkpt==0)then
! Might restrict the number of k points to be printed
  tnkpt=0
  nkpt_eff=nkpt
  if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
   nkpt_eff=nkpt_max
   tnkpt=1
  end if
  intarr(1:nkpt_eff,0)=0
  intarr(1:nkpt_eff,1:ndtset_alloc)=istwfk_2(1:nkpt_eff,1:ndtset_alloc)
  call prttagm(dprarr,intarr,iout,jdtset_,1,marr,nkpt_eff,&
&  ndtset_alloc,'istwfk','INT')
  if(tnkpt==1 .and. sum(istwfk_2(1:nkpt_eff,1:ndtset_alloc))/=0 ) &
&  write(iout,'(16x,a)' ) &
&  'outvar1 : prtvol=0, do not print more k-points.'

 else
  do idtset=1,ndtset_alloc
   tnkpt=0
   nkpt_eff=dtsets(idtset)%nkpt
   if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
    nkpt_eff=nkpt_max
    tnkpt=1
   end if
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   if(sum(istwfk_2(1:nkpt_eff,idtset))/=0)then
    write(iout,format01155a)'istwfk',appen,istwfk_2(1:nkpt_eff,idtset)
    if(tnkpt==1) write(iout,'(16x,a)' ) &
&    'outvar1 : prtvol=0, do not print more k-points.'
   end if
  end do
 end if
 deallocate(istwfk_2)

 intarr(1,:)=dtsets(:)%ixc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ixc','INT')

 intarr(1,:)=dtsets(:)%ixcpositron
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ixcpositron','INT')

 if (ndtset > 0) write(iout,format01155) 'jdtset',jdtset_(1:ndtset)

 intarr(1,:)=dtsets(:)%jellslab
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'jellslab','INT')

!kberry
 if(multi_nberry==0)then
  if(nberry/=0)then
   do idtset=0,ndtset_alloc
    intarr(1:3*nberry,idtset)=&
&    reshape( dtsets(idtset)%kberry(1:3,1:nberry), (/3*nberry/) )
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3*nberry,&
&   ndtset_alloc,'kberry','INT')
  end if
 else
  do idtset=1,ndtset_alloc
   if(dtsets(idtset)%nberry>0)then
    jdtset=jdtset_(idtset)
    if(jdtset<10)write(appen,'(i1)')jdtset
    if(jdtset>=10)write(appen,'(i2)')jdtset
    write(iout,format01155a)&
&    'kberry',appen,dtsets(idtset)%kberry(1:3,1:dtsets(idtset)%nberry)
   end if
  end do
 end if


!kpt
 if(multi_nkpt==0)then
! Might restrict the number of k points to be printed
  tnkpt=0
  nkpt_eff=nkpt
  if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
   nkpt_eff=nkpt_max
   tnkpt=1
  end if
  do idtset=0,ndtset_alloc
   dprarr(1:3*nkpt_eff,idtset)=&
&   reshape(dtsets(idtset)%kpt(1:3,1:nkpt_eff),(/3*nkpt_eff/) )
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3*nkpt_eff,&
&  ndtset_alloc,'kpt','DPR')
  if(tnkpt==1) write(iout,'(16x,a)' ) &
&  '       outvar1 : prtvol=0, do not print more k-points.'
 else
  do idtset=1,ndtset_alloc
   tnkpt=0
   nkpt_eff=dtsets(idtset)%nkpt
   if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
    nkpt_eff=nkpt_max
    tnkpt=1
   end if
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01150a)'kpt',appen,dtsets(idtset)%kpt(1:3,1:nkpt_eff)
   if(tnkpt==1) write(iout,'(16x,a)' ) &
&   'outvar1 : prtvol=0, do not print more k-points.'
  end do
 end if

!kptgw
 if(multi_nkptgw==0)then
! Required if for pathscale to avoid failure with -ff_bounds_check
  if (nkptgw > 0) then
   do idtset=0,ndtset_alloc
    dprarr(1:3*nkptgw,idtset)=&
&    reshape(dtsets(idtset)%kptgw(1:3,1:nkptgw),(/3*nkptgw/) )
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3*nkptgw,&
&   ndtset_alloc,'kptgw','DPR')
  end if
 else
  do idtset=1,ndtset_alloc
   if(dtsets(idtset)%nkptgw>0)then
    jdtset=jdtset_(idtset)
    if(jdtset<10)write(appen,'(i1)')jdtset
    if(jdtset>=10)write(appen,'(i2)')jdtset
    write(iout,format01150a)'kptgw',appen,dtsets(idtset)%kptgw(1:3,1:dtsets(idtset)%nkptgw)
   end if
  end do
 end if

 dprarr(1,:)=dtsets(:)%kptnrm
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'kptnrm','DPR')

 dprarr(1,:)=dtsets(:)%kptrlen
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'kptrlen','DPR')

 intarr(1,:)=dtsets(:)%kptopt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'kptopt','INT')

!kptrlatt
 if(sum((dtsets(1:ndtset_alloc)%kptopt)**2)/=0)then
  ndtset_kptopt=0
  intarr(1:9,0)=reshape( dtsets(0)%kptrlatt(:,:) , (/9/) )
  allocate(jdtset_kptopt(0:ndtset_alloc))
! Define the set of datasets for which kptopt>0
  do idtset=1,ndtset_alloc
   kptopt=dtsets(idtset)%kptopt
   if(kptopt>0)then
    ndtset_kptopt=ndtset_kptopt+1
    jdtset_kptopt(ndtset_kptopt)=jdtset_(idtset)
    intarr(1:9,ndtset_kptopt)=reshape( dtsets(idtset)%kptrlatt(:,:) , (/9/) )
   end if
  end do
  if(ndtset_kptopt>0)then
   call prttagm(dprarr,intarr,iout,jdtset_kptopt,3,marr,9,&
&   ndtset_kptopt,'kptrlatt','INT')
  end if
  deallocate(jdtset_kptopt)
 end if

 intarr(1,:)=dtsets(:)%kssform
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'kssform','INT')

 intarr(1,:)=dtsets(:)%ldgapp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ldgapp','INT')

 intarr(1,:)=dtsets(:)%localrdwf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'localrdwf','INT')

 intarr(1,:)=dtsets(:)%lofwrite
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'lofwrite','INT')

!ltypeorb
 if(multi_norb==0)then
  do idtset=0,ndtset_alloc
   intarr(1:norb,idtset)=dtsets(idtset)%ltypeorb(1:norb)
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,1,marr,norb,ndtset_alloc,'ltypeorb','INT')
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01155a)'ltypeorb',appen,dtsets(idtset)%ltypeorb(1:dtsets(idtset)%norb)
  end do
 end if


 dprarr(1,:)=dtsets(:)%mdftemp
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'mdftemp','DPR')

 dprarr(1,:)=dtsets(:)%mditemp
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'mditemp','DPR')

 dprarr(1,:)=dtsets(:)%mdwall
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'mdwall','LEN')

 intarr(1,:)=dtsets(:)%mffmem
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'mffmem','INT')

!mixalch
 if(multi_ntypalch==0)then
  do idtset=0,ndtset_alloc
   npspalch=dtsets(idtset)%npspalch
   dprarr(1:npspalch*ntypalch,idtset)=&
&   reshape(dtsets(idtset)%mixalch(1:npspalch,1:ntypalch),(/npspalch*ntypalch/))
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,1,marr,npspalch*ntypalch,ndtset_alloc,'mixalch','DPR')
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01160a)'mixalch',appen,&
&   dtsets(idtset)%mixalch(1:dtsets(idtset)%npspalch,1:dtsets(idtset)%ntypalch)
  end do
 end if

!DEBUG
!write(6,*)' outvar1 : after mixalch '
!ENDDEBUG


 intarr(1,:)=dtsets(:)%mkmem
 call prttagm(dprarr,intarr,iout,jdtset_,5,marr,1,ndtset_alloc,'mkmem','INT')

 if(response==1)then
  intarr(1,:)=dtsets(:)%mkqmem
  call prttagm(dprarr,intarr,iout,jdtset_,5,marr,1,ndtset_alloc,'mkqmem','INT')
 end if

 if(response==1)then
  intarr(1,:)=dtsets(:)%mk1mem
  call prttagm(dprarr,intarr,iout,jdtset_,5,marr,1,ndtset_alloc,'mk1mem','INT')
 end if

 intarr(1,:)=dtsets(:)%mqgrid
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'mqgrid','INT')
 if (usepaw==1) then
  intarr(1,:)=dtsets(:)%mqgriddg
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'mqgriddg','INT')
 end if

 intarr(1,:)=natfix_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'natfix','INT')

 intarr(1,:)=natfixx_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'natfixx','INT')

 intarr(1,:)=natfixy_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'natfixy','INT')

 intarr(1,:)=natfixz_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'natfixz','INT')

 intarr(1,:)=dtsets(0:ndtset_alloc)%natom
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'natom','INT')

 intarr(1,:)=dtsets(:)%nfreqim
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nfreqim','INT')

 intarr(1,:)=dtsets(:)%nfreqre
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nfreqre','INT')

 intarr(1,:)=dtsets(:)%nfreqsp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nfreqsp','INT')

 intarr(1,:)=dtsets(:)%npfft
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'npfft','INT')

 intarr(1,:)=dtsets(:)%npband
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'npband','INT')

 intarr(1,:)=dtsets(:)%npkpt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'npkpt','INT')

 intarr(1,:)=dtsets(:)%bandpp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'bandpp','INT')

 intarr(1,:)=dtsets(0:ndtset_alloc)%npulayit
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'npulayit','INT')

!nqptdm
 if(multi_nqptdm==0)then
! Required if for pathscale to avoid failure with -ff_bounds_check
  if (nqptdm > 0) then
   do idtset=0,ndtset_alloc
    dprarr(1:3*nqptdm,idtset)=&
&    reshape(dtsets(idtset)%qptdm(1:3,1:nqptdm),(/3*nqptdm/) )
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3*nqptdm,&
&   ndtset_alloc,'nqptdm','DPR')
  end if
 else
  do idtset=1,ndtset_alloc
   if(dtsets(idtset)%nqptdm>0)then
    jdtset=jdtset_(idtset)
    if(jdtset<10)write(appen,'(i1)')jdtset
    if(jdtset>=10)write(appen,'(i2)')jdtset
    write(iout,format01150a)'qptdm',appen,dtsets(idtset)%qptdm(1:3,1:dtsets(idtset)%nqptdm)
   end if
  end do
 end if


!

 deallocate(dprarr,intarr)
 deallocate(natfix_,iatfixio_)
 deallocate(natfixx_,iatfixx_)
 deallocate(natfixy_,iatfixy_)
 deallocate(natfixz_,iatfixz_)

!DEBUG
!write(6,*)' outvar1 : end of subroutine '
!if(.true.)stop
!ENDDEBUG

end subroutine outvar1
!!***
