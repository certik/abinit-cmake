!{\src2tex{textfont=tt}}
!!****f* ABINIT/hdr_init
!! NAME
!! hdr_init
!!
!! FUNCTION
!! This subroutine initializes the header structured datatype
!! and most of its content from dtset and psps, and put default values for
!! evolving variables.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! bstruct <type(bandstructure_type)>=band structure information
!!  including Brillouin zone description
!! codvsn=code version
!! dtset <type(dataset_type)>=all input variables for this dataset
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! pertcase=index of the perturbation, or 0 if GS calculation
!! psps <type(pseudopotential_type)>=all the information about psps
!!
!! OUTPUT
!! hdr <type(hdr_type)>=the header, initialized, and for most part of
!!   it, contain its definite values, except for evolving variables
!!
!! PARENTS
!!      gstate,loper3,newsp,nonlinear,respfn
!!
!! CHILDREN
!!      date_and_time,leave_new,rhoij_alloc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine hdr_init(bstruct,codvsn,dtset,hdr,pawtab,pertcase,psps)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: pertcase
 character(len=6),intent(in) :: codvsn
 type(bandstructure_type),intent(in) :: bstruct
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(out) :: hdr
 type(pseudopotential_type),intent(in) :: psps
!arrays
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: bantot,date,iatom,itypat,natom,nkpt,npsp,nselect,nsppol,nsym,ntypat
 character(len=500) :: message
 character(len=8) :: date_time
!arrays
 integer,allocatable :: nlmn(:)

! *************************************************************************

!More checking would be needed ...
 if(dtset%ntypat/=psps%ntypat)then
  write(message, '(4a,2i6)' ) ch10,&
&  ' hdr_init : BUG -',ch10,&
&  '  dtset%ntypat and psps%ntypat differs. They are :',dtset%ntypat,psps%ntypat
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 if(dtset%npsp/=psps%npsp)then
  write(message, '(4a,2i6)' ) ch10,&
&  ' hdr_init : BUG -',ch10,&
&  '  dtset%npsp and psps%npsp differs. They are :',dtset%npsp,psps%npsp
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if


 call date_and_time(date_time)
 read(date_time,'(i8)')date

 natom  = dtset%natom
 npsp   = dtset%npsp
 nsym   = dtset%nsym
 ntypat = dtset%ntypat

 nkpt   = bstruct%nkpt
 nsppol = bstruct%nsppol
 bantot = bstruct%bantot

!DEBUG
!write(6,*)' hdr_init : before allocate'
!ENDDEBUG

!Allocate all components of hdr
 allocate(hdr%istwfk(nkpt))
 allocate(hdr%kptns(3,nkpt))
 allocate(hdr%lmn_size(npsp))
 allocate(hdr%nband(nkpt*nsppol))
 allocate(hdr%npwarr(nkpt)) ! Warning : npwarr here has only one dim
 allocate(hdr%occ(bantot))
 allocate(hdr%pspcod(npsp))
 allocate(hdr%pspdat(npsp))
 allocate(hdr%pspso(npsp))
 allocate(hdr%pspxc(npsp))
 allocate(hdr%so_psp(npsp))
 allocate(hdr%symafm(nsym))
 allocate(hdr%symrel(3,3,nsym))
 allocate(hdr%title(npsp))
 allocate(hdr%tnons(3,nsym))
 allocate(hdr%typat(natom))
 allocate(hdr%wtk(nkpt))
 allocate(hdr%xred(3,natom))
 allocate(hdr%zionpsp(npsp))
 allocate(hdr%znuclpsp(npsp))
 allocate(hdr%znucltypat(ntypat))
 if(psps%usepaw==1)then
  allocate(hdr%pawrhoij(natom))
  allocate(nlmn(ntypat))
  do itypat=1,ntypat
   nlmn(itypat)=pawtab(itypat)%lmn_size
  end do
  call rhoij_alloc(1,nlmn,dtset%nspden,dtset%nsppol,hdr%pawrhoij,dtset%typat)
  deallocate(nlmn)
 end if

!DEBUG
!write(6,*)' hdr_init : after allocate'
!ENDDEBUG

!Transfer data from dtset
 hdr%intxc    =dtset%intxc
 hdr%ixc      =dtset%ixc
 hdr%natom    =natom
 hdr%npsp     =npsp
 hdr%nspden   =dtset%nspden
 hdr%nspinor  =dtset%nspinor
 hdr%nsym     =nsym
 hdr%ntypat    =ntypat
 hdr%occopt   =dtset%occopt
 if(psps%usepaw==1)then
  hdr%ngfft(:) =dtset%ngfftdg(1:3)
 else if (dtset%usewvl == 1) then
  hdr%ngfft(:) = dtset%wvl_internal%dpSize(:)
 else
  hdr%ngfft(:) =dtset%ngfft(1:3)
 end if
 hdr%so_psp(:)=dtset%so_psp(:)
 hdr%symafm(1:min(size(dtset%symafm),size(hdr%symafm)))=dtset%symafm(1:min(size(dtset%symafm),size(hdr%symafm)))
 hdr%symrel(:,:,1:min(size(dtset%symrel,3),size(hdr%symrel,3))) =dtset%symrel(:,:,1:min(size(dtset%symrel,3),size(hdr%symrel,3)))
 hdr%typat(1:dtset%natom)  =dtset%typat(1:dtset%natom) ! PMA : in tests/v2/t11 size(dtset%typat) is bigger dtset%natom
 hdr%ecut     =dtset%ecut
 hdr%ecutsm   =dtset%ecutsm
 hdr%ecut_eff =dtset%ecut * (dtset%dilatmx)**2
 hdr%qptn(:)  =dtset%qptn(:)
 hdr%stmbias  =dtset%stmbias
 hdr%tnons(:,1:min(size(dtset%tnons,2),size(hdr%tnons,2)))    =dtset%tnons(:,1:min(size(dtset%tnons,2),size(hdr%tnons,2)))
 hdr%tphysel  =dtset%tphysel
 hdr%tsmear   =dtset%tsmear
 hdr%rprimd(:,:)=dtset%rprimd_orig(:,:)      ! Evolving data
 hdr%xred(:,1:dtset%natom)=dtset%xred_orig(:,1:dtset%natom)          ! Evolving data
!PMA : in tests/v2/t11 size(dtset%xred_orig,2) is bigger dtset%natom

!Transfer wavelets data.
 hdr%usewvl     = dtset%usewvl
!hdr%nwvlarr will be set later since the number
!of wavelets have not yet been computed.

!DEBUG
!write(6,*)' hdr_init : before transfer data from bstruct'
!ENDDEBUG

!Transfer data from bstruct
 hdr%bantot        =bantot
 hdr%nkpt          =nkpt
 hdr%nsppol        =nsppol
 hdr%istwfk(1:nkpt)=bstruct%istwfk(1:nkpt)
 hdr%nband(1:nkpt*nsppol) =bstruct%nband(1:nkpt*nsppol)
 hdr%npwarr(:)     =bstruct%npwarr(:)
 hdr%kptns(:,:)    =bstruct%kptns(:,:)
 hdr%wtk(:)        =bstruct%wtk(:)
 hdr%occ(:)        =bstruct%occ(:)           ! Evolving data

!DEBUG
!write(6,*)' hdr_init : before transfer data from psps'
!ENDDEBUG

!Transfer data from psps
 hdr%pspcod(:)=psps%pspcod(:)
 hdr%pspdat(:)=psps%pspdat(:)
 hdr%pspso (:)=psps%pspso (:)
 hdr%pspxc (:)=psps%pspxc (:)
 hdr%znuclpsp(:)=psps%znuclpsp(:)
 hdr%znucltypat(:)=psps%znucltypat(:)
 hdr%zionpsp(:)=psps%zionpsp(:)
 hdr%title (:)=psps%title (:)

!Transfer paw data
 hdr%usepaw=psps%usepaw
 if(psps%usepaw==1) then
  hdr%ecutdg   =dtset%pawecutdg
  hdr%lmn_size(1:npsp)=pawtab(1:npsp)%lmn_size
 else
  hdr%ecutdg=hdr%ecut
  hdr%lmn_size(:)=psps%lmnmax
 end if

!Initialize other known data
 hdr%codvsn   =codvsn
 hdr%date     =date
 hdr%headform =44              ! The present header format has been defined in v4.4
!!  headform=format of the header (44 if rdwr=2, while
!!   if rdwr=1, 22 for versions before 2.3, 23 for 2.3 until 3.3
!!   34 for 3.4, 40 for 4.0, 41 for 4.1 and later. Should not be 1,2,51,52,101 or 102
 hdr%pertcase =pertcase

!Default for other data  (all evolving data)
 hdr%etot     =1.0d20
 hdr%fermie   =1.0d20
 hdr%residm   =1.0d20

end subroutine hdr_init
!!***
