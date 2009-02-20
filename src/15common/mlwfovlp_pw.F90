!{\src2tex{textfont=tt}}
!!****f* ABINIT/mlwfovlp_pw
!! NAME
!! mlwfovlp_pw
!!
!! FUNCTION
!! Routine which computes PW part of overlap M_{mn}(k,b) 
!! for Wannier code (www.wannier.org f90 version).
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (BAmadon,FJollet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  g1(3,nkpt,nntot) = G vector shift which is necessary to obtain k1+b
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw.
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor) (see NOTES at beginning of scfcv)
!!  ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  ovikp(nkpt,nntot)= gives  nntot value of k2 (in the BZ) for each k1  (k2=k1+b mod(G))
!!  prtvol=control print volume and debugging output
!!
!! OUTPUT
!!  cm1(2,nkpt,nntot,mbandw,mbandw): overlap <u_(nk1)|u_(mk1+b)>.
!!  iwav(nsppol,nkpt,mbandw): shift for pw components in cg.
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      mlwfovlp
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine mlwfovlp_pw(cg,cm1,dtset,g1,iwav,kg,mband,mbandw,mkmem,mpsang,mpw,natom,&
& nfft,ngfft,nkpt,nntot,npwarr,nspden,nspinor,nsppol,ntypat,ovikp,prtvol)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mbandw,mkmem,mpsang,mpw,natom,nfft,nkpt,nntot
 integer,intent(in) :: nspden,nspinor,nsppol,ntypat,prtvol
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: g1(3,nkpt,nntot),kg(3,mpw*mkmem),ngfft(18),npwarr(nkpt)
 integer,intent(in) :: ovikp(nkpt,nntot)
 integer,intent(out) :: iwav(nsppol,nkpt,mband)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(out) :: cm1(2,nkpt,nntot,mband,mband)

!Local variables-------------------------------
!scalars
 integer :: iband,iband1,iband2,icgtemp,icpb,icspol,idum,ig,ig1,ig1b,ig2,ig2b
 integer :: ig3,ig3b,igk1,igk2,ii,ikg,ikpt,ikpt1,ikpt2,imntot,intot,ispden
 integer :: ispinor,isppol,n1,n2,n3,nband_k,npoint,npoint2,npw_k
 character(len=500) :: message
 character(len=fnlen) :: fildata
!arrays
 integer,allocatable :: icg(:,:),indpwk(:,:),invpwk(:,:),kg_k(:,:)

!************************************************************************

 write(message, '(a,a)' ) ch10,&
& '** mlwfovlp_pw : compute pw part of overlap'
 call wrtout(06,  message,'COLL')

!MJV 6/2008 : added error checking
 if (mkmem /= nkpt) then
  write(message, '(a,a,a,a)' ) ch10,&
&  '** mlwfovlp_pw : not all kpoints are in memory, or you are using this in parallel.',&
&  ch10,'Action: run a sequential job with mkmem = nkpt'
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!****************compute intermediate quantities  (index, shifts) ******
!------------compute index for g points--------------------------------
!ig is a plane waves which belongs to the sphere ecut for ikpt (they
!are npwarr(ikpt))
!npoint is the position in the grid of planes waves
!(they are nfft)
!indpwk is a application ig-> npoint
!invpwk is not an application (some npoint have no ig corresponding)
!cg are ordered with npw_k !
!----------------------------------------------------------------------
!------------compute index for g points--------------------------------
!----------------------------------------------------------------------
 write(message, '(a,a)' ) ch10,&
& '   first compute index for g-points'
 call wrtout(06,  message,'COLL')
 allocate(kg_k(3,mpw),indpwk(nkpt,mpw),invpwk(nkpt,nfft))
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 ikg=0
 invpwk=-1
 indpwk=0
 kg_k=0
 do ikpt=1,nkpt
  do npoint=1,nfft
   if(invpwk(ikpt,npoint)/=-1) then
    write(6,*) "error0 , invpwk is overwritten"
    write(6,*) ikpt,npoint
    stop
   end if
  end do
  npw_k=npwarr(ikpt)
! write(6,*) ikpt,npw_k,nfft
  kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
  do ig=1,npw_k
   if(ig.gt.mpw) then
    write(6,*)"error ig",ig,mpw
    stop
   end if
   if(indpwk(ikpt,ig)/=0) then
    write(6,*) "error, indpwk is overwritten"
    write(6,*) ikpt,ig,indpwk(ikpt,ig)
    stop
   end if
   ig1=modulo(kg_k(1,ig),n1)
   ig2=modulo(kg_k(2,ig),n2)
   ig3=modulo(kg_k(3,ig),n3)
   indpwk(ikpt,ig)=ig1+1+n1*(ig2+n2*ig3)
   npoint=indpwk(ikpt,ig)
   if(npoint.gt.nfft) then
    write(6,*)"error npoint"
    stop
   end if
!  write(6,*) ikpt,ig,npoint,invpwk(ikpt,npoint)
   if(invpwk(ikpt,npoint)/=-1) then
    write(6,*) "error, invpwk is overwritten"
    write(6,*) ikpt,ig,npoint,invpwk(ikpt,npoint)
    stop
   end if
   invpwk(ikpt,npoint)=ig
!  if(ikpt.eq.1) write(6,*) "ig npoint",ig, npoint
  end do
  ikg=ikg+npw_k
 end do
!write(6,*) "index for g points has been computed"

!write(6,*) "Computes shift for cg"
 write(message, '(a,a)' ) ch10,&
& '   mlwfovlp_pw : compute shifts for g-points '
 call wrtout(06,  message,'COLL')
!----------------------------------------------------------------------
!------------compute shifts for g points (icg,iwav)---------------------
!------------ (here mbandw is not used, because shifts are internal----
!variables of abinit)--------------------------------------------------
!----------------------------------------------------------------------
!write(6,*) mpw*nspinor*mband*mkmem*nsppol
 allocate(icg(nsppol,nkpt))
 icg=0
 icgtemp=0
 iwav=0
 do isppol=1,nsppol
  do ikpt=1,nkpt
   nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
!  write(6,*) ikpt+(isppol-1)*nkpt,nkpt
   npw_k=npwarr(ikpt)
   do iband=1,nband_k
    if(iband.gt.mband) then
     write(6,*)"error mband",iband,mband,nband_k
     stop
    end if
    iwav(isppol,ikpt,iband)= &
&    (iband-1)*npw_k*nspinor+icgtemp
!   write(6,*) "iwav", isppol,ikpt,iband,iwav(isppol,ikpt,iband)
   end do ! iband
   icgtemp=icgtemp+ npw_k*nspinor*nband_k
   icg(isppol,ikpt)=icgtemp
!  write(6,*) "icg", isppol,ikpt,icg(isppol,ikpt)
  end do  ! ikpt
 end do   ! isppol
!write(6,*) "shift for cg computed"

!----------------------------------------------------------------------
!------------test invpwk-----------------------------------------------
!----------------------------------------------------------------------
!write(6,*) "TEST INVPWK"
 ikpt=1
 isppol=1
 do ig=1,npwarr(ikpt)
  npoint=indpwk(ikpt,ig)
! write(6,*) "ig npoint    ",ig, npoint
! write(6,*) "ig npoint inv",invpwk(ikpt,npoint),npoint
 end do
 do ig3=1,n3
  do ig2=1,n2
   do ig1=1,n1
    npoint=ig1+(ig2-1)*n1+(ig3-1)*n2*n1
    ig=invpwk(ikpt,npoint)
!   if(ig/=-1)  write(6,*) "ig npoint",ig, npoint
   end do
  end do
 end do


!***********************************************************************
!**calculate overlap M_{mn}(k,b)=<\Psi_{k,m}|e^{-ibr}|\Psi_{k+b,n}>*****
!***********************************************************************
 write(message, '(a,a)' ) ch10,&
& '   mlwfovlp_pw : compute overlaps '
 call wrtout(06,  message,'COLL')
 write(message, '(a,a)' ) ch10,&
& "     nkpt  nntot  mbandw  mbandw"
 call wrtout(06,  message,'COLL')
 write(message, '(i6,2x,i6,2x,i6,2x,i6)' ) &
& nkpt,nntot,mbandw,mbandw
 call wrtout(06,  message,'COLL')
 cm1=zero
 write(message, '(a)' )  '  '
 call wrtout(06,  message,'COLL')
 do isppol=1,nsppol
  imntot=0
  do ikpt1=1,nkpt
   write(message, '(a,i6)' ) &
&   '   compute overlaps for k-point=',ikpt1
   call wrtout(06,  message,'COLL')
   do intot=1,nntot
    imntot=imntot+1
    ikpt2= ovikp(ikpt1,intot)
    do ig3=1,n3
     do ig2=1,n2
      do ig1=1,n1
!      write(6,*) isppol,ikpt1,iband1,iband2,intot
       npoint=ig1+(ig2-1)*n1+(ig3-1)*n2*n1
       if(npoint.gt.nfft) then
        write(6,*) "error npoint  b"
        stop
       end if
       ig1b=ig1+g1(1,ikpt1,intot)
       ig2b=ig2+g1(2,ikpt1,intot)
       ig3b=ig3+g1(3,ikpt1,intot)
!      write(6,*) ig1,ig2,ig3
!      write(6,*) ig1b,ig2b,ig3b
       if(ig1b.lt.1) ig1b=ig1b+n1
       if(ig2b.lt.1) ig2b=ig2b+n2
       if(ig3b.lt.1) ig3b=ig3b+n3
       if(ig1b.gt.n1) ig1b=ig1b-n1
       if(ig2b.gt.n2) ig2b=ig2b-n2
       if(ig3b.gt.n3) ig3b=ig3b-n3
       npoint2=ig1b+(ig2b-1)*n1+(ig3b-1)*n2*n1
       if(npoint2.gt.nfft) then
        write(6,*)"error npoint  c"
        stop
       end if
       igk1=invpwk(ikpt1,npoint)
       igk2=invpwk(ikpt2,npoint2)
       if(igk1/=-1.and.igk2/=-1) then
        do iband1=1,mbandw
         do iband2=1,mbandw
          cm1(1,ikpt1,intot,iband1,iband2)=cm1(1,ikpt1,intot,iband1,iband2)+ &
&          cg(1,igk1+iwav(isppol,ikpt1,iband1))*cg(1,igk2+iwav(isppol,ikpt2,iband2))&
&          + cg(2,igk1+iwav(isppol,ikpt1,iband1))*cg(2,igk2+iwav(isppol,ikpt2,iband2))
          cm1(2,ikpt1,intot,iband1,iband2)=cm1(2,ikpt1,intot,iband1,iband2)+ &
&          cg(1,igk1+iwav(isppol,ikpt1,iband1))*cg(2,igk2+iwav(isppol,ikpt2,iband2))&
&          - cg(2,igk1+iwav(isppol,ikpt1,iband1))*cg(1,igk2+iwav(isppol,ikpt2,iband2))
!         write(1111,*) ikpt1,intot,iband1,iband2
!         write(1111,*) iwav(isppol,ikpt1,iband1),iwav(isppol,ikpt2,iband2)
!         write(1111,*) ig3,ig2,ig1, cm1(1,ikpt1,intot,iband1,iband2), cm1(2,ikpt1,intot,iband1,iband2)
         end do ! iband2
        end do ! iband1
       end if
      end do
     end do
    end do
   end do ! intot
  end do ! ikpt1
 end do ! isppol
 deallocate(kg_k,indpwk,invpwk,icg)

!DEBUG
!write(6,*)' mlwfovlp_pw : exit'
!stop
!ENDDEBUG

 end subroutine mlwfovlp_pw
!!***
