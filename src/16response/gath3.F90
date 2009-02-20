!{\src2tex{textfont=tt}}
!!****f* ABINIT/gath3
!!
!! NAME
!! gath3
!!
!! FUNCTION
!! Sum (gather) the different parts of the 2nd-order matrix,
!! to get the matrix of second-order derivatives (d2matr)
!! Then, generates the dynamical matrix, not including the masses,
!! but the correct non-cartesian coordinates ( => d2cart)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! berryopt=option for berry phase treatment
!! blkflg(3,mpert,3,mpert)= ( 1 if the element of the dynamical
!!  matrix has been calculated ; 0 otherwise )
!! dyew(2,3,natom,3,natom)=Ewald part of the dyn.matrix
!! dyfrwf(3,3,natom)=frozen wf part of the dyn.matrix (except xc1)
!! dyfrx1(2,3,natom,3,natom)=xc core correction (1) part of the frozen-wf
!!  part of the dynamical matrix.
!! d2bbb(2,3,3,mpert,mband,mband*prtbbb)=band by band decomposition of some
!!       second order derivatives
!! d2nfr(2,3,mpert,3,mpert)=non-frozen wf part of the 2nd-order matr
!! eltcore(6,6)=core contribution to the elastic tensor
!! elteew(6+3*natom,6)=Ewald contribution to the elastic tsenor
!! eltfrhar(6,6)=hartree contribution to the elastic tensor
!! eltfrkin(6,6)=kinetic contribution to the elastic tensor
!! eltfrloc(6+3*natom,6)=local psp contribution to the elastic tensor
!! eltfrnl(6+3*natom,6)=non-local psp contribution to the elastic tensor
!! eltfrxc(6+3*natom,6)=exchange-correlation contribution to the elastic tensor
!! gprimd(3,3)=basis vector in the reciprocal space
!! mband=maximum number of bands
!! mpert =maximum number of ipert
!! natom=number of atoms in unit cell
!! ntypat=number of atom types
!! outd2=option for the output of the 2nd-order matrix :
!!  if outd2=1, non-stationary part
!!  if outd2=2, stationary part.
!! prtbbb=if 1, print the band-by-band decomposition, otherwise, prtbbb=0
!! rfasr= (0=> no acoustic sum rule [asr] imposed), (1=> asr is imposed,
!!  in the democratic way for the effective charges),
!! (2=> asr is imposed, in the aristocratic way for the effective
!!  charges)
!! rfdir(3)=define the directions for the perturbations
!! rfpert(mpert)=define the perturbations
!! rprimd(3,3)=dimensional primitive translations (bohr)
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume
!! zion(ntypat)=charge corresponding to the atom type
!!
!! OUTPUT
!! carflg(3,mpert,3,mpert)= ( 1 if the element of the cartesian
!!  2DTE matrix has been calculated correctly ; 0 otherwise )
!! d2cart(2,3,mpert,3,mpert)=
!!  dynamical matrix, effective charges, dielectric tensor,....
!!  all in cartesian coordinates
!! d2cart_bbb(2,3,3,mpert,mband,mband*prtbbb)=
!!  band by band decomposition of Born effective charges
!!  (calculated from phonon-type perturbation) in cartesian coordinates
!! d2matr(2,3,mpert,3,mpert)=2nd-order matrix (masses non included,
!!  no cartesian coordinates : simply second derivatives)
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      asria9,cart29,chneu9
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine gath3(berryopt,blkflg,carflg,dyew,dyfrwf,dyfrx1,&
& d2bbb,d2cart,d2cart_bbb,d2matr,d2nfr,&
& eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,&
& gprimd,mband,mpert,natom,ntypat,outd2,prtbbb,&
& rfasr,rfdir,rfpert,rprimd,typat,ucvol,zion)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_16response, except_this_one => gath3
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: berryopt,mband,mpert,natom,ntypat,outd2,prtbbb,rfasr
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: rfdir(3),rfpert(mpert),typat(natom)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert)
 integer,intent(out) :: carflg(3,mpert,3,mpert)
 real(dp),intent(in) :: d2bbb(2,3,3,mpert,mband,mband*prtbbb)
 real(dp),intent(in) :: d2nfr(2,3,mpert,3,mpert),dyew(2,3,natom,3,natom)
 real(dp),intent(in) :: dyfrwf(3,3,natom),dyfrx1(2,3,natom,3,natom)
 real(dp),intent(in) :: eltcore(6,6),elteew(6+3*natom,6),eltfrhar(6,6)
 real(dp),intent(in) :: eltfrkin(6,6),eltfrloc(6+3*natom,6)
 real(dp),intent(in) :: eltfrnl(6+3*natom,6),eltfrxc(6+3*natom,6),gprimd(3,3)
 real(dp),intent(in) :: rprimd(3,3),zion(ntypat)
 real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)
 real(dp),intent(out) :: d2cart_bbb(2,3,3,mpert,mband,mband*prtbbb)
 real(dp),intent(out) :: d2matr(2,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: chneut,i1,i2,iband,iblok,idir,idir1,idir2,ii,ipert,ipert1,ipert2,iw
 integer :: jj,mu,nblok,nsize,selectz
 real(dp) :: tmp
 character(len=500) :: message
!arrays
 integer :: flg1(3),flg2(3)
 real(dp) :: acg(3),sum(2),ter(3,3),vec1(3),vec2(3)
 real(dp),allocatable :: d2tmp(:,:,:,:,:),d2work(:,:,:,:,:),elfrtot(:,:)

! *********************************************************************

 write(6,*)' gath3 : enter '

!DEBUG
!write(6,*)' outd2,mpert =',outd2,mpert
!write(6,*)' blkflg(:,natom+2,:,natom+2)=',blkflg(:,natom+2,:,natom+2)
!ENDDEBUG

 if(outd2/=3)then

! Initialise the 2nd-derivative matrix
  d2matr(:,:,:,:,:)=0.0_dp

! Add the non-frozen-part, the
! Ewald part and the xc1 part of the frozen-wf part
  do ipert2=1,mpert
   do idir2=1,3
    do ipert1=1,mpert
     do idir1=1,3
      if( blkflg(idir1,ipert1,idir2,ipert2)==1 ) then
       do ii=1,2
        d2matr(ii,idir1,ipert1,idir2,ipert2)=&
&        d2nfr(ii,idir1,ipert1,idir2,ipert2)
        if(ipert1<=natom .and. ipert2<=natom)then
         d2matr(ii,idir1,ipert1,idir2,ipert2)=&
&         d2matr(ii,idir1,ipert1,idir2,ipert2)+&
&         dyew(ii,idir1,ipert1,idir2,ipert2)  +&
&         dyfrx1(ii,idir1,ipert1,idir2,ipert2)
        end if
       end do
      end if
     end do
    end do
   end do
  end do

! Add the frozen-wavefunction part
  do ipert2=1,natom
   do idir2=1,3
    do idir1=1,3
     if( blkflg(idir1,ipert2,idir2,ipert2)==1 ) then
      d2matr(1,idir1,ipert2,idir2,ipert2)=&
&      d2matr(1,idir1,ipert2,idir2,ipert2)&
&      +dyfrwf(idir1,idir2,ipert2)
     end if
    end do
   end do
  end do

! Section for strain perturbation
  if(rfpert(natom+3)==1 .or. rfpert(natom+4)==1) then
!  Make sure relevant columns of output are nulled
   d2matr(:,:,:,:,natom+3:natom+4)=0.0_dp
!  Accumulate all frozen parts of the elastic tensor
   allocate(elfrtot(6+3*natom,6))
   elfrtot(:,:)= elteew(:,:)+eltfrloc(:,:)+eltfrnl(:,:)+eltfrxc(:,:)

   elfrtot(1:6,1:6)=elfrtot(1:6,1:6)+&
&   eltcore(:,:)+eltfrhar(:,:)+eltfrkin(:,:)


   do ipert2=natom+3,natom+4
    do idir2=1,3
!    Internal strain components first
     do ipert1=1,natom
      do idir1=1,3
       if( blkflg(1,ipert1,idir2,ipert2)==1 ) then
        ii=idir1+6+3*(ipert1-1)
        jj=idir2+3*(ipert2-natom-3)
        d2matr(1,idir1,ipert1,idir2,ipert2)=&
&        d2nfr(1,idir1,ipert1,idir2,ipert2)+elfrtot(ii,jj)
       end if
      end do
     end do
!    Now, electric field - strain mixed derivative (piezoelectric tensor)
     ipert1=natom+2
     do idir1=1,3
      if( blkflg(idir1,ipert1,idir2,ipert2)==1 ) then
       d2matr(1,idir1,ipert1,idir2,ipert2)=&
&       d2nfr(1,idir1,ipert1,idir2,ipert2)
      end if
     end do
!    Now, strain-strain 2nd derivatives
     do ipert1=natom+3,natom+4
      do idir1=1,3
       if( blkflg(idir1,ipert1,idir2,ipert2)==1 ) then
        ii=idir1+3*(ipert1-natom-3)
        jj=idir2+3*(ipert2-natom-3)
        d2matr(1,idir1,ipert1,idir2,ipert2)=&
&        d2nfr(1,idir1,ipert1,idir2,ipert2)+elfrtot(ii,jj)
       end if
      end do
     end do
    end do
   end do
   deallocate(elfrtot)
  end if
! End section for strain perturbation

! The second-order matrix has been computed.

! Filter now components smaller in absolute value than 1.0d-20,
! for automatic testing reasons
  do ipert2=1,mpert
   do idir2=1,3
    do ipert1=1,mpert
     do idir1=1,3
      if( blkflg(idir1,ipert1,idir2,ipert2)==1 ) then
       do ii=1,2
        if( d2matr(ii,idir1,ipert1,idir2,ipert2)**2 < 1.0d-40)then
         d2matr(ii,idir1,ipert1,idir2,ipert2)=zero
        end if
       end do
      end if
     end do
    end do
   end do
  end do

! DEBUG
! write(6,*)' d2matr '
! ipert2=natom+2
! do idir2=1,3
! ipert1=natom+2
! do idir1=1,3
! write(6, '(4i4,2es16.6)' )idir1,ipert1,idir2,ipert2,&
! &       d2matr(1,idir1,ipert1,idir2,ipert2),&
! &       d2matr(2,idir1,ipert1,idir2,ipert2)
! end do
! end do
! ENDDEBUG

! Cartesian coordinates transformation
  iblok=1 ; nblok=1
  nsize=3*mpert*3*mpert

! In the case of finite electric field, the convention for the 
! direction of the electric field perturbation was NOT the usual convention ...
! So, there is a transformation to the usual conventions
! to be done first ...
  if(berryopt==4 .and. minval(abs(blkflg(:,natom+2,:,natom+2)))/=0)then
   if(minval(abs(blkflg(:,natom+2,:,natom+2)-1))/=0)then
    write(message,'(8a)')ch10,&
&    ' gath3 : ERROR -',ch10,&
&    '  In case of finite electric field, and electric field perturbation,',ch10,&
&    '  the three directions for the perturbations must be treated.',ch10,&
&    '  Action : set idir to 1 1 1, or forget about finite electric field.'
    call wrtout(std_out,message,'COLL')
    call leave_new('COLL')
   end if
   do ipert=1,mpert
    do idir=1,3
     do ii=1,2
      vec1(:)=d2matr(ii,idir,ipert,:,natom+2)
      flg1(:)=blkflg(idir,ipert,:,natom+2)
      call cart39(flg1,flg2,gprimd,1,1,rprimd,vec1,vec2)
      d2matr(ii,idir,ipert,:,natom+2)=vec2(:)*two_pi
      blkflg(idir,ipert,:,natom+2)=flg2(:)
     end do
    end do
   end do
   do ipert=1,mpert
    do idir=1,3
     do ii=1,2
      vec1(:)=d2matr(ii,:,natom+2,idir,ipert)
      flg1(:)=blkflg(:,natom+2,idir,ipert)
      call cart39(flg1,flg2,gprimd,1,1,rprimd,vec1,vec2)
      d2matr(ii,:,natom+2,idir,ipert)=vec2(:)*two_pi
      blkflg(:,natom+2,idir,ipert)=flg2(:)
     end do
    end do
   end do
!  Also to be done, a change of sign, that I do not understand (XG071110)
!  Perhaps due to d/dk replacing id/dk ? !
   d2matr(:,:,natom+2,:,natom+2)=-d2matr(:,:,natom+2,:,natom+2)
  end if

  call cart29(blkflg,d2matr,carflg,d2cart,&
&  gprimd,iblok,mpert,natom,nblok,nsize,ntypat,rprimd,typat,ucvol,zion)

! Band by band decomposition of the Born effective charges
  if(prtbbb==1)then
   allocate(d2work(2,3,mpert,3,mpert),d2tmp(2,3,mpert,3,mpert))
   tmp=0.0_dp
   do iband=1,mband
    d2work(:,:,:,:,:)=0.0_dp
    d2tmp(:,:,:,:,:)=0.0_dp
    d2work(:,:,natom+2,:,:) = d2bbb(:,:,:,:,iband,iband)
    call cart29(blkflg,d2work,carflg,d2tmp,&
&    gprimd,iblok,mpert,natom,nblok,nsize,ntypat,rprimd,typat,ucvol,zion)

!   Remove the ionic part
    do ipert1=1,natom
     do idir1=1,3
      d2tmp(1,idir1,natom+2,idir1,ipert1) = &
&      d2tmp(1,idir1,natom+2,idir1,ipert1) - zion(typat(ipert1))
     end do
    end do

    d2cart_bbb(:,:,:,:,iband,iband)=d2tmp(:,:,natom+2,:,:)
   end do
   deallocate(d2tmp,d2work)
  end if ! prtbbb==1

! 
! Now, the cartesian elements are ready for output
! carflg give the information on their correctness

! Imposition of the ASR on the analytical part of the DynMat
! Assume that if rfasr/=0, the whole cartesian matrix is correct
 end if

 if(rfasr/=0)then

  allocate(d2work(2,3,mpert,3,mpert))
  call asria9(rfasr,1,d2work,d2cart,mpert,natom)
  deallocate(d2work)

! Imposition of the ASR on the effective charges.
  if(rfpert(natom+2)==1)then
   chneut=rfasr
   selectz=0
   call chneu9(chneut,d2cart,mpert,natom,ntypat,selectz,typat,zion)
  end if

 end if

!Additional operations on cartesian strain derivatives
 if(rfpert(natom+3)==1 .or. rfpert(natom+4)==1) then
! Impose zero-net-force condition on internal strain tensor
  do ipert2=natom+3,natom+4
   do idir2=1,3
    vec1(:)=0.0_dp
    do ipert1=1,natom
     do idir1=1,3
      if(carflg(idir1,ipert1,idir2,ipert2)==1) then
       vec1(idir1)=vec1(idir1)+d2cart(1,idir1,ipert1,idir2,ipert2)
      end if
     end do
    end do
    vec1(:)=vec1(:)/dble(natom)
    do ipert1=1,natom
     do idir1=1,3
      if(carflg(idir1,ipert1,idir2,ipert2)==1) then
!      Note minus sign to convert gradients to forces
       d2cart(1,idir1,ipert1,idir2,ipert2)=&
&       -(d2cart(1,idir1,ipert1,idir2,ipert2)-vec1(idir1))
      end if
     end do
    end do
   end do
  end do
! Divide strain 2nd deriviative by ucvol to give elastic tensor
  do ipert2=natom+3,natom+4
   do idir2=1,3
    do ipert1=natom+3,natom+4
     do idir1=1,3
      if(carflg(idir1,ipert1,idir2,ipert2)==1) then
       d2cart(1,idir1,ipert1,idir2,ipert2)=&
&       d2cart(1,idir1,ipert1,idir2,ipert2)/ucvol
      end if
     end do
    end do
   end do
  end do
 end if

!calculate Born effective charges from electric field perturbation
!do ipert1=1,natom
!ter(:,:)=zero
!do idir1=1,3
!do ii=1,3
!do jj=1,3
!if(abs(gprimd(idir1,ii))>1.0d-10)then
!ter(idir1,ii)=ter(idir1,ii)+ d2nfr(1,idir1,natom+2,jj,ipert1)*gprimd(jj,ii)  
!endif
!enddo
!enddo
!add zion to bec
!ter(idir1,idir1)=ter(idir1,idir1)+zion(typat(ipert1))
!enddo
!d2cart(1,:,ipert1,:,natom+2)=ter(:,:)
!enddo
!carflg(:,1:natom,:,natom+2)=1

!Born effective charges from phonon perturbation
!do ipert1=1,natom
!ter(:,:)=zero
!do idir1=1,3
!do ii=1,3
!do jj=1,3
!if(abs(gprimd(idir1,ii))>1.0d-10)then
!ter(idir1,ii)=ter(idir1,ii)+ d2nfr(1,jj,ipert1,idir1,natom+2)*gprimd(jj,ii)
!endif
!enddo
!enddo
!! add zion to bec
!ter(idir1,idir1)=ter(idir1,idir1)+zion(typat(ipert1))
!enddo
!d2cart(1,:,natom+2,:,ipert1)=ter(:,:)
!enddo
!carflg(:,natom+2,:,1:natom)=1


!!Dielectric constant
!do ii=1,3
!do jj=1,3
!ter(ii,jj)=d2nfr(1,ii,natom+2,jj,natom+2)
!end do
!end do
!ter(:,:)=pi*four*ter(:,:)/ucvol
!
!do ii=1,3
!ter(ii,ii)=ter(ii,ii)+one
!end do
!d2cart(1,:,natom+2,:,natom+2)=ter(:,:)
!carflg(:,natom+2,:,1:natom+2)=1

!DEBUG
!Debugging, but needed to make the stuff work on the IBM Dirac ? !
!write(6,*)' d2cart '
!ipert2=natom+2
!do idir2=1,3
!ipert1=natom+2
!do idir1=1,3
!write(6, '(5i4,2d20.10)' )idir1,ipert1,idir2,ipert2,&
!&      carflg(idir1,ipert1,idir2,ipert2),&
!&      d2cart(1,idir1,ipert1,idir2,ipert2),&
!&      d2cart(2,idir1,ipert1,idir2,ipert2)
!end do
!end do
!ENDDEBUG

!DEBUG
!write(6,*)' gath3 : exit '
!ENDDEBUG

end subroutine gath3
!!***
