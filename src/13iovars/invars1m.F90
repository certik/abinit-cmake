!{\src2tex{textfont=tt}}
!!****f* ABINIT/invars1m
!! NAME
!! invars1m
!!
!! FUNCTION
!! Initialisation phase : prepare the main input subroutine call by
!! reading all the NO MULTI variables, as well as the dimensions
!! needed for allocating the input arrays in abinit.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG, MKV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  iout=unit number of output file
!!  lenstr=actual length of string
!!  mpi_enreg=informations about MPI parallelization
!!  msym=default maximal number of symmetries
!!  mxnatom=maximal value of input natom for all the datasets
!!  ndtset= number of datasets to be read; if 0, no multi-dataset mode
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!               one data set.
!!  string*(*)=string of characters containing all input variables and data
!!  zion_max=maximal valence charge over all psps
!!
!! OUTPUT
!!  bravais_(11,0:ndtset_alloc)=characteristics of Bravais lattice (see symbrav.f)
!!  dmatpuflag=flag controlling the use of an initial density matrix in PAW+U (max. value over datasets)
!!  mband_upper_(0:ndtset_alloc)=list of mband_upper values
!!  mxlpawu=maximal value of input lpawu for all the datasets
!!  mxmband_upper=maximal value of input nband for all the datasets
!!  mxnatpawu=maximal value of number of atoms on which +U is applied for all the datasets
!!  mxnatsph=maximal value of input natsph for all the datasets
!!  mxnatvshift=maximal value of input natvshift for all the datasets
!!  mxncenter=maximal value of input ncenter for all the datasets
!!  mxnconeq=maximal value of input nconeq for all the datasets
!!  mxnkptgw=maximal value of input nkptgw for all the datasets
!!  mxnorb=maximal value of input norb for all the datasets
!!  mxnqptdm=maximal value of input nqptdm for all the datasets
!!  mxnkpt=maximal value of input nkpt for all the datasets
!!  mxnspinor=maximal value of input nspinor for all the datasets
!!  mxnsppol=maximal value of input nsppol for all the datasets
!!
!! SIDE EFFECTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are initialized here (see invars1.f for more details on the
!!   initialized records)
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      invars1,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine invars1m(bravais_,dmatpuflag,dtsets,iout,lenstr,mband_upper_,mpi_enreg,&
& msym,mxlpawu,mxmband_upper,mxnatom,mxnatpawu,mxnatsph,mxnatvshift,mxncenter,mxnconeq,&
& mxnkptgw,mxnkpt,mxnorb,mxnnos,mxnqptdm,mxnspinor,mxnsppol,mxnsym,&
& ndtset,ndtset_alloc,string,zion_max)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13iovars, except_this_one => invars1m
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,lenstr,msym,mxnatom,ndtset,ndtset_alloc
 integer,intent(out) :: dmatpuflag,mxlpawu,mxmband_upper,mxnatpawu,mxnatsph
 integer,intent(out) :: mxnatvshift,mxncenter,mxnconeq,mxnkpt,mxnkptgw,mxnnos,mxnorb
 integer,intent(out) :: mxnqptdm,mxnspinor,mxnsppol,mxnsym
 real(dp),intent(in) :: zion_max
 character(len=*),intent(inout) :: string
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(out) :: bravais_(11,0:ndtset_alloc)
 integer,intent(out) :: mband_upper_(0:ndtset_alloc)
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,iatom,idtset,ii,ipsp,isym,itypat,jdtset,lpawu,mband_upper
 integer :: natpawu,nkpt,tjdtset,tread
 character(len=500) :: message
!arrays
 integer :: bravais(11)
 integer,allocatable :: symafm_(:,:)
 integer,allocatable :: symrel_(:,:,:,:)
 real(dp),allocatable :: tnons_(:,:,:)
 integer,allocatable :: symafm(:),symrel(:,:,:)
 real(dp),allocatable :: tnons(:,:)

!******************************************************************

!Here, the default msym, and allocation of the three arrays that depend on msym.
 allocate(symrel_(3,3,msym,0:ndtset_alloc))
 allocate(symafm_(msym,0:ndtset_alloc))
 allocate(tnons_(3,msym,0:ndtset_alloc))

 allocate(symafm(msym),symrel(3,3,msym),tnons(3,msym))

!Set up default values (note that the default acell, amu
!mkmem, mkmem1,mkqmem, and nkpt must be overcome
 bravais_(:,0)=0
 do idtset=0,ndtset_alloc
  dtsets(idtset)%amu(:)=-one
  dtsets(idtset)%acell_orig(:)=zero
  dtsets(idtset)%berryopt=0
! Only one component of densty is used presently
  dtsets(idtset)%densty(:,:)=zero
  dtsets(idtset)%efield(:)=zero
  dtsets(idtset)%useexexch=0
  dtsets(idtset)%iatfix(:,:)=0
  dtsets(idtset)%jellslab=0
  dtsets(idtset)%kptopt=0
  dtsets(idtset)%lpawu(:)=-1
  dtsets(idtset)%lexexch(:)=-1
  dtsets(idtset)%nnos=0
  dtsets(idtset)%natsph=0
  dtsets(idtset)%natvshift=0
  dtsets(idtset)%nbdblock=1
  dtsets(idtset)%ncenter=0
  dtsets(idtset)%nconeq=0
  dtsets(idtset)%nkptgw=0
  dtsets(idtset)%npara=0
  dtsets(idtset)%npband=1
  dtsets(idtset)%npfft=1
  dtsets(idtset)%npkpt=1
  dtsets(idtset)%norb=0
  dtsets(idtset)%nqptdm=0
  dtsets(idtset)%optdriver=0
  dtsets(idtset)%parareel=0
  dtsets(idtset)%paral_kgb=0
  dtsets(idtset)%rprim_orig(:,:)=zero
  dtsets(idtset)%rprim_orig(1,1)=one
  dtsets(idtset)%rprim_orig(2,2)=one
  dtsets(idtset)%rprim_orig(3,3)=one
  dtsets(idtset)%slabzbeg=zero
  dtsets(idtset)%slabzend=zero
  dtsets(idtset)%spinat(:,:)=zero
  dtsets(idtset)%symmorphi=1
  dtsets(idtset)%typat(:)=0
  dtsets(idtset)%usedmatpu=0
  dtsets(idtset)%usepawu=0
  dtsets(idtset)%vel_orig(:,:)=zero
  dtsets(idtset)%xred_orig(:,:)=zero
  dtsets(idtset)%so_psp(:)=1
  dtsets(idtset)%w90nplot=0
  dtsets(idtset)%wfoptalg=0
 end do
 dtsets(:)%mkmem=-1
 dtsets(:)%mk1mem=-1
 dtsets(:)%mkqmem=-1
!natom is already initialized in invars0
 dtsets(0)%natom=-1
 dtsets(:)%natpawu=0
 dtsets(:)%nkpt=-1
 dtsets(:)%nspden=1
 dtsets(:)%nspinor=1
 dtsets(:)%nsppol=1
 dtsets(:)%nsym=0
 dtsets(:)%pawspnorb=0
 symafm_(:,0)=1
 symrel_(:,:,:,0)=0
 symrel_(1,1,:,0)=1 ; symrel_(2,2,:,0)=1 ; symrel_(3,3,:,0)=1
 tnons_(:,:,0)=0.0_dp

!Loop on datasets
 do idtset=1,ndtset_alloc

  jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0

  write(message, '(a,a,i6)') ch10,&
&  ' invars1m : enter jdtset=',jdtset
  call wrtout(6,message,'COLL')

! Input default values
  bravais(:)=bravais_(:,0)
  symafm(:)=symafm_(:,0)
  symrel(:,:,:)=symrel_(:,:,:,0)
  tnons(:,:)=tnons_(:,:,0)
  call invars1(bravais,dtsets(idtset),iout,jdtset,lenstr,&
&  mband_upper,mpi_enreg,msym,string,symafm,symrel,tnons,zion_max)

  bravais_(:,idtset)=bravais(:)
  mband_upper_ (idtset)=mband_upper
  symafm_(:,idtset)=symafm(:)
  symrel_(:,:,:,idtset)=symrel(:,:,:)
  tnons_(:,:,idtset)=tnons(:,:)

 end do

 mxmband_upper =maxval(mband_upper_ (1:ndtset_alloc))

 dmatpuflag=0;mxnatpawu=0;mxlpawu=0
 mxnatsph=dtsets(1)%natsph
 mxnatvshift=dtsets(1)%natvshift
 mxncenter=dtsets(1)%ncenter
 mxnconeq=dtsets(1)%nconeq
 mxnkptgw=dtsets(1)%nkptgw
 mxnkpt  =dtsets(1)%nkpt
 mxnnos  =dtsets(1)%nnos
 mxnorb  =dtsets(1)%norb
 mxnqptdm=dtsets(1)%nqptdm
 mxnspinor=dtsets(1)%nspinor
 mxnsppol=dtsets(1)%nsppol
 do ii=1,ndtset_alloc
  mxnatsph=max(dtsets(ii)%natsph,mxnatsph)
  mxnatvshift=max(dtsets(ii)%natvshift,mxnatvshift)
  mxncenter=max(dtsets(ii)%ncenter,mxncenter)
  mxnconeq=max(dtsets(ii)%nconeq,mxnconeq)
  mxnkptgw=max(dtsets(ii)%nkptgw,mxnkptgw)
  mxnkpt  =max(dtsets(ii)%nkpt,mxnkpt)
  mxnnos  =max(dtsets(ii)%nnos,mxnnos)
  mxnorb  =max(dtsets(ii)%norb,mxnorb)
  mxnqptdm=max(dtsets(ii)%nqptdm,mxnqptdm)
  mxnspinor=max(dtsets(ii)%nspinor,mxnspinor)
  mxnsppol=max(dtsets(ii)%nsppol,mxnsppol)
  if (dtsets(ii)%usepawu>0) then
   if (dtsets(ii)%usedmatpu/=0) dmatpuflag=1
   lpawu=maxval(dtsets(ii)%lpawu(:))
   mxlpawu=max(lpawu,mxlpawu)
   do iatom=1,dtsets(ii)%natom
    itypat=dtsets(ii)%typat(iatom)
    if (dtsets(ii)%lpawu(itypat)/=-1) dtsets(ii)%natpawu=dtsets(ii)%natpawu+1
   end do
   mxnatpawu=max(dtsets(ii)%natpawu,mxnatpawu)
  end if
 end do
 deallocate(symafm,symrel,tnons)

!mxnsym=maxval(dtsets(1:ndtset_alloc)%nsym) ! This might not work properly with HP compiler
 mxnsym=dtsets(1)%nsym
 do idtset=1,ndtset_alloc
  mxnsym=max(dtsets(idtset)%nsym,mxnsym)
 end do

 do idtset=0,ndtset_alloc
  allocate(dtsets(idtset)%atvshift(mxnatvshift,mxnsppol,mxnatom))
  allocate(dtsets(idtset)%bdgw(2,mxnkptgw))
  allocate(dtsets(idtset)%dmatpawu(2*mxlpawu+1,2*mxlpawu+1,max(mxnsppol,mxnspinor),mxnatpawu*dmatpuflag))
  allocate(dtsets(idtset)%kpt(3,mxnkpt))
  allocate(dtsets(idtset)%kptgw(3,mxnkptgw))
  allocate(dtsets(idtset)%kptns(3,mxnkpt))
  allocate(dtsets(idtset)%iatsph(mxnatsph))
  allocate(dtsets(idtset)%istwfk(mxnkpt))
  allocate(dtsets(idtset)%ltypeorb(mxnorb))
  allocate(dtsets(idtset)%nband(mxnkpt*mxnsppol))
  allocate(dtsets(idtset)%numorb(mxncenter))
  allocate(dtsets(idtset)%occ_orig(mxmband_upper*mxnkpt*mxnsppol))
  allocate(dtsets(idtset)%qmass(mxnnos))
  allocate(dtsets(idtset)%qptdm(3,mxnqptdm))
  allocate(dtsets(idtset)%rcoord(3,mxncenter))
  allocate(dtsets(idtset)%rtheta(3,mxnorb))
  allocate(dtsets(idtset)%symafm(mxnsym))
  allocate(dtsets(idtset)%symrel(3,3,mxnsym))
  allocate(dtsets(idtset)%tnons(3,mxnsym))
  allocate(dtsets(idtset)%wtatcon(3,mxnatom,mxnconeq))
  allocate(dtsets(idtset)%wtk(mxnkpt))
  allocate(dtsets(idtset)%w90lplot(dtsets(idtset)%w90nplot))
  dtsets(idtset)%symrel(:,:,:)=symrel_(:,:,1:mxnsym,idtset)
  dtsets(idtset)%symafm(:)    =symafm_(1:mxnsym,idtset)
  dtsets(idtset)%tnons (:,:)  =tnons_ (:,1:mxnsym,idtset)
 end do
 deallocate(symafm_,symrel_,tnons_)

!DEBUG
!write(6,*)' invars1m : exit'
!stop
!ENDDEBUG

end subroutine invars1m
!!***
