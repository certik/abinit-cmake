!{\src2tex{textfont=tt}}
!!****f* ABINIT/getshell
!! NAME
!! getshell
!!
!! FUNCTION
!! For each k-point, set up the shells of first neighbours and find
!! the weigths required for the finite difference expression
!! of Marzari and Vanderbilt (see PRB 56, 12847 (1997)).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gmet(3,3) = metric tensor of reciprocal space
!! kptopt = option for the generation of k points
!! kptrlatt = k-point lattice specification
!! kpt2(3,nkpt2) = reduced coordinates of the k-points in the
!!                 reduced part of the BZ (see below)
!! mkmem = number of k points which can fit in memory
!! mpi_enreg = informations about MPI parallelization
!! nkpt2 = number of k-points in the reduced BZ
!! nkpt3 = number of k-points in the full BZ
!! nshiftk = number of kpoint grid shifts
!! rmet(3,3) = metric tensor of real space
!! rprimd(3,3) = dimensional primitive translations (bohr)
!! shiftk = shift vectors for k point generation
!! wtk2 = weight assigned to each k point
!!
!! OUTPUT
!! kneigh(30,nkpt2) = for each k-point in the reduced part of the BZ
!!                    kneigh stores the index (ikpt) of the neighbouring
!!                    k-points
!! kptindex(2,nkpt3)
!!   kptindex(1,ikpt) = ikpt_rbz
!!     ikpt_rbz = index of the k-point in the reduced BZ
!!     ikpt = index of the k-point in the full BZ
!!   kptindex(2,ikpt) = 1: use time-reversal symmetry to transform the
!!                         wavefunction at ikpt_rbz to the wavefunction at ikpt
!!                      0: ikpt belongs already to the reduced BZ
!!                         (no transformation required)
!! kpt3(3,nkpt3) = reduced coordinates of the k-points in the full BZ
!! mvwtk(30,nkpt2) = weights required to evaluate the finite difference
!!                   formula of Marzari and Vanderbilt, computed for each
!!                   k-point in the reduced part of the BZ
!! mkmem_max = maximal number of k-points on each processor (MPI //)
!! nneigh = total number of neighbours required to evaluate the finite
!!          difference formula
!!
!! COMMENTS
!! The array kpt2 holds the reduced coordinates of the k-points in the
!! reduced part of the BZ. For example, in case time-reversal symmetry is
!! used (kptopt = 2) kpt2 samples half the BZ. Since some of the neighbours
!! of these k-points may lie outside the reduced BZ, getshell also needs the
!! coordinates of the k-points in the full BZ.
!! The coordinates of the k-points in the full BZ are stored in kpt3.
!! The weights mvwtk are computed for the k-points kpt2.
!!
!! In case no symmetry is used to reduce the number of k-points,
!! the arrays kpt2 and kpt3 are equal.
!!
!!
!! PARENTS
!!      nonlinear
!!
!! CHILDREN
!!      dgetrf,dgetri,getkgrid,leave_new,wrtout,xcomm_world,xmax_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine getshell(gmet,kneigh,kptindex,kptopt,kptrlatt,kpt2,&
& kpt3,mkmem,mkmem_max,mpi_enreg,mvwtk,&
& nkpt2,nkpt3,nneigh,nshiftk,rmet,rprimd,shiftk,wtk2)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13recipspace
 use interfaces_lib01hidempi
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kptopt,mkmem,nkpt2,nkpt3
 integer,intent(inout) :: nshiftk
 integer,intent(out) :: mkmem_max,nneigh
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(inout) :: kptrlatt(3,3)
 integer,intent(out) :: kneigh(30,nkpt2),kptindex(2,nkpt3)
 real(dp),intent(in) :: gmet(3,3),kpt2(3,nkpt2),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk),wtk2(nkpt2)
 real(dp),intent(out) :: kpt3(3,nkpt3),mvwtk(30,nkpt2)

!Local variables-------------------------------
!scalars
 integer :: bis,flag,ier,ii,ikpt,ikpt1,ikpt2,ikpt3,ineigh,info,irank,is1,ishell
 integer :: jj,kk,kptopt_used,mkmem_cp,ndiff,nkpt_computed,nshell,nsym1,orig
 integer :: spaceComm,wtkflg
 real(dp) :: dist_,dtm,kptrlen,last_dist,max_dist,resdm,s1
 character(len=500) :: message
!arrays
 integer :: dsifkpt(3),idiff(2,6),neigh(0:6,nkpt2),symafm_dummy(1),vacuum(3)
 integer,allocatable :: ipiv(:),symrel1(:,:,:)
 real(dp) :: dist(6),dk(3),dk_(3),mat(6,6),resdk(3),rvec(6),sgval(6)
 real(dp) :: shiftk_(3,8),work(30)
 real(dp),allocatable :: tnons1(:,:),wtk3(:)

!************************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'DGETRI' :: dgetri
!DEC$ ATTRIBUTES ALIAS:'DGETRF' :: dgetrf
#endif

!DEBUG
!write(6,*)'getshell : enter'
!write(6,*)'nkpt2 : ',nkpt2
!write(6,*)'nkpt3 : ',nkpt3
!stop
!ENDDEBUG

!In case of MPI //: compute maximum number of k-points per processor
 if (mpi_enreg%paral_compil_kpt == 1) then
! BEGIN TF_CHANGES
  call xcomm_world(mpi_enreg,spaceComm)
! END TF_CHANGES
  mkmem_cp=mkmem
  call xmax_mpi(mkmem_cp,mkmem_max,spaceComm,ier)
 else
  mkmem_max = mkmem
 end if

!------------- In case kptopt = 2 set up the whole k-point grid -------------

!kpt3(3,nkpt3) = reduced coordinates of k-points in the full BZ

 if (kptopt == 3) then

  allocate(wtk3(nkpt3))
  kpt3(:,:) = kpt2(:,:)
  wtk3(:) = wtk2(:)
  do ikpt = 1,nkpt3
   kptindex(1,ikpt) = ikpt
   kptindex(2,ikpt) = 0
  end do

 else if (kptopt == 2) then

  allocate(wtk3(nkpt3))
  dsifkpt(:) = 1 ; ii = 5 ; kptopt_used = 3
  symafm_dummy(1) = 1
  shiftk_(:,:) = 0._dp
  shiftk_(:,1:nshiftk) = shiftk(:,1:nshiftk)

  nsym1 = 1
  allocate(symrel1(3,3,nsym1),tnons1(3,nsym1))
  symrel1(:,:,1) = 0
  symrel1(1,1,1) = 1 ; symrel1(2,2,1) = 1 ; symrel1(3,3,1) = 1
  tnons1(:,:) = 0._dp
  vacuum(:) = 0

  call getkgrid(dsifkpt,ab_out,ii,kpt3,kptopt_used,kptrlatt,&
&  kptrlen,nsym1,nkpt3,nkpt_computed,nshiftk,nsym1,&
&  rprimd,shiftk_,symafm_dummy,symrel1,tnons1,&
&  vacuum,wtk3)

  if (nkpt_computed /= nkpt3) then
   write(message,'(a,a,a,a,i4,a,a,i4)') ch10,&
&   ' mv_3dte: BUG - ',ch10,&
&   ' The number of k-points in the whole BZ, nkpt_computed= ',nkpt_computed,&
&   ch10,&
&   ' is not twice the number of k-points in half the BZ, nkpt3=',nkpt3
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
   call leave_new('COLL')
  end if

  kptindex(:,:) = 0
  do ikpt3 = 1, nkpt3

   flag = 1
   do ikpt2 = 1, nkpt2

!   In case, the k-points differ only by one reciprocal lattice
!   vector, apply shift of one g-vector to kpt(:,ikpt3)

    dk_(:) = kpt3(:,ikpt3) - kpt2(:,ikpt2)
    dk(:) = dk_(:) - nint(dk_(:))
    if (dk(1)*dk(1) + dk(2)*dk(2) + dk(3)*dk(3) < tol10) then
     do ii = 1, 3
      if ((dk(ii)*dk(ii) < tol10).and.(dk_(ii)*dk_(ii) > tol10)) then
       kpt3(ii,ikpt3) = -1._dp*kpt3(ii,ikpt3)
      end if
     end do
    end if

    dk_(:) = kpt3(:,ikpt3) + kpt2(:,ikpt2)
    dk(:) = dk_(:) - nint(dk_(:))
    if (dk(1)*dk(1) + dk(2)*dk(2) + dk(3)*dk(3) < tol10) then
     do ii = 1, 3
      if ((dk(ii)*dk(ii) < tol10).and.(dk_(ii)*dk_(ii) > tol10)) then
       kpt3(ii,ikpt3) = -1._dp*kpt3(ii,ikpt3)
      end if
     end do
    end if

    dk(:) = kpt3(:,ikpt3) - kpt2(:,ikpt2)
    if (dk(1)*dk(1) + dk(2)*dk(2) + dk(3)*dk(3) < tol10) then
     kptindex(1,ikpt3) = ikpt2
     kptindex(2,ikpt3) = 0       ! no use of time-reversal symmetry
     flag = 0
     exit
    end if

    dk(:) = kpt3(:,ikpt3) + kpt2(:,ikpt2)
    if (dk(1)*dk(1) + dk(2)*dk(2) + dk(3)*dk(3) < tol10) then
     kptindex(1,ikpt3) = ikpt2
     kptindex(2,ikpt3) = 1       ! use time-reversal symmetry
     flag = 0
     exit
    end if

   end do     ! ikpt2

   if (flag == 1) then
    write(message,'(a,a,a,a,i4)') ch10,&
&    ' mv_3dte: BUG - ',ch10,&
&    ' Could not find a symmetric k-point for ikpt3=  ',&
&    ikpt3
    call wrtout(ab_out,message,'COLL')
    call wrtout(06,  message,'COLL')
    call leave_new('COLL')
   end if
  end do    ! ikpt3

 else

  write(message,'(a,a,a,a)') ch10,&
&  ' mv_3dte: ERROR - ',ch10,&
&  ' the only values for kptopt that are allowed are 2 and 3 '
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')

 end if   ! condition on kptopt


!--------- Compute the weights required for the Marzari-Vanderbilt ---------
!--------- finite difference formula ---------------------------------------


!Initialize distance between k-points
!The trace of gmet is an upper limit for its largest eigenvalue. Since the
!components of the distance vectors do not exceed 1, 3. * Tr[gmet] is
!an upper limit for the squared shell radius.
 dist_ = 0._dp
 do ii = 1,3
  dist_ = dist_ + gmet(ii,ii)
 end do
 max_dist = 3._dp * dist_

!Calculate an upper limit for the residuum
 resdm = rmet(1,1)*rmet(1,1) + rmet(2,2)*rmet(2,2) + rmet(3,3)*rmet(3,3)&
& + rmet(1,2)*rmet(1,2) + rmet(2,3)*rmet(2,3) + rmet(3,1)*rmet(3,1)

!Initialize shell loop 
 ishell = 0
 last_dist = 0._dp
 wtkflg = 0
 kneigh(:,:) = 0
 neigh(:,:) = 0

!Loop over shells until the residuum is zero
 do while ((wtkflg == 0).and.(resdm > tol8))
! Advance shell counter
  ishell = ishell + 1

! Initialize shell radius with upper limit
  dist(ishell) = max_dist

! Find the (squared) radius of the next shell
  do ikpt = 2,nkpt3
   dk(:) = kpt3(:,1) - kpt3(:,ikpt)
   dist_ = 0._dp
   do ii = 1,3
    do jj = 1,3
     dist_ = dist_ + dk(ii)*gmet(ii,jj)*dk(jj)
    end do
   end do

   if ((dist_ < dist(ishell)).and.(dist_ - last_dist > tol8)) then
    dist(ishell) = dist_
   end if
  end do

! DEBUG
! write(6,*)'ishell, dist = ',ishell,dist(ishell)
! ENDDEBUG
  last_dist = dist(ishell)

! For each k-point in halft the BZ get the shells of nearest neighbours.
! These neighbours can be out of the zone sampled by kpt2.
  do ikpt2 = 1, nkpt2              ! k-points in half the BZ
   orig = sum(neigh(0:ishell-1,ikpt2))
   nneigh = 0
   do ikpt3 = 1, nkpt3             ! whole k-point grid
    dk(:) = kpt3(:,ikpt3) - kpt2(:,ikpt2)
    dk_(:) = dk(:) - nint(dk(:))
    dist_ = 0._dp
    do ii = 1,3
     do jj = 1,3
      dist_ = dist_ + dk_(ii)*gmet(ii,jj)*dk_(jj)
     end do
    end do
    if (abs(dist_ - dist(ishell)) < tol8) then
     nneigh = nneigh + 1
     kneigh(orig+nneigh,ikpt2) = ikpt3
    end if
   end do
   neigh(ishell,ikpt2) = nneigh
  end do


! Check if the number of points in shell number ishell
! is the same for each k-point

  flag = 1
  do ikpt = 1,nkpt2
   if (neigh(ishell,ikpt) /= nneigh) flag = 0
  end do

  if (flag == 0) then
   write(message,'(a,a,a,a,i2,a,a)') ch10,&
&   ' getshell: BUG - ',ch10,&
&   ' The number of points in shell number',ishell,' is not the same',&
&   ' for each k-point.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
   call leave_new('COLL')
  end if

  if (nneigh == 0) then
   write(message,'(a,a,a,a)') ch10,&
&   ' getshell: BUG - ',ch10,&
&   ' Cannot find enough neighbor shells'
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
   wtkflg = 1
  end if
  
! Calculate the total number of neighbors 
  nneigh = sum(neigh(1:ishell,1))
! DEBUG
! write(6,*)'ishell = ',ishell,'nneigh = ',nneigh
! ENDDEBUG

! Find the weights needed to compute the finite difference expression
! of the ddk
! **********************************************************************

! mvwtk(:,:) = 0._dp

! The weights are calculated for ikpt=1. The results are copied later
  ikpt = 1

! Calculate the coefficients of the linear system to be solved
  mat(:,:) = 0._dp
  do is1 = 1, ishell
   orig = sum(neigh(0:is1-1,ikpt))
   bis = orig + neigh(is1,ikpt)
   do ineigh = orig+1, bis
    dk_(:) = kpt3(:,kneigh(ineigh,ikpt)) - kpt2(:,ikpt)
    dk(:) = dk_(:) - nint(dk_(:))
    mat(1,is1) = mat(1,is1) + dk(1)*dk(1)
    mat(2,is1) = mat(2,is1) + dk(2)*dk(2)
    mat(3,is1) = mat(3,is1) + dk(3)*dk(3)
    mat(4,is1) = mat(4,is1) + dk(1)*dk(2)
    mat(5,is1) = mat(5,is1) + dk(2)*dk(3)
    mat(6,is1) = mat(6,is1) + dk(3)*dk(1)
   end do
  end do

  rvec(1) = rmet(1,1)
  rvec(2) = rmet(2,2)
  rvec(3) = rmet(3,3)
  rvec(4) = rmet(1,2)
  rvec(5) = rmet(2,3)
  rvec(6) = rmet(3,1)
  
! DEBUG
! do ii = 1, 6
! write(6,*)mat(ii,1:ishell), ' : ', rvec(ii)
! end do
! ENDDEBUG

! Solve the linear least square problem 
  call dgelss(6,ishell,1,mat,6,rvec,6,sgval,tol8,irank,work,30,info)

  if( info /= 0 ) then
   write(message,'(a,a,a,a,a,a)') 'getshell : COMMENT -', ch10,&
&   ' Singular-value decomposition of the linear system determining the',ch10,&
&   ' weights failed (info /= 0).',ch10
   call wrtout(06,message,'COLL')
   write (*,*) 'info = ',info
   wtkflg = 1
  end if
  
! Check that the system has maximum rank
  if( irank == ishell ) then
!  System has full rank. Calculate the residuum
   s1 = resdm
   resdm = 0._dp
   do is1 = ishell + 1, 6
    resdm = resdm + rvec(is1) * rvec(is1)
   end do
   
   if( ishell == 6 .and. resdm > tol8 ) then
    write(message,'(a,a,a,a,a,a)') 'getshell : COMMENT -', ch10,&
&    ' Linear system determining the weights could not be solved',ch10,&
&    ' This should not happen.',ch10
    call wrtout(06,message,'COLL')
    wtkflg = 1
   end if
  else
!  The system is rank deficient
   ishell = ishell - 1
!  DEBUG
!  write(*,*) 'Shell not linear independent from previous shells. Skipped.'
!  ENDDEBUG
  end if
  
! DEBUG
! write(6,*) ishell, nneigh, irank, resdm
! ENDDEBUG

! end of loop over shells
 end do

!Copy weights
 ikpt=1
 do is1 = 1, ishell
  orig = sum(neigh(0:is1-1,ikpt))
  bis = orig + neigh(is1,ikpt)
  mvwtk(orig+1:bis,1) = rvec(is1)
 end do
 do ikpt = 2,nkpt2
  mvwtk(1:nneigh,ikpt) = mvwtk(1:nneigh,1)
 end do  ! ikpt

!Report weights
 write(6,*) 'Neighbors', neigh(1:ishell,1)
 write(6,*) 'Weights', rvec(1:ishell)
 write(6,*) mvwtk(1:nneigh,1) 

!Check the computed weights
 if (wtkflg == 0) then 
  do ikpt = 1, nkpt2
   do ii = 1,3
    do jj = 1,3
     s1 = 0._dp
     do ineigh = 1, nneigh
      dk_(:) = kpt3(:,kneigh(ineigh,ikpt)) - kpt2(:,ikpt)
      dk(:) = dk_(:) - nint(dk_(:))
      s1 = s1 + dk(ii)*dk(jj)*mvwtk(ineigh,ikpt)
     end do
     if (abs(s1 - rmet(ii,jj)) > tol6) wtkflg = 1
    end do
   end do
  end do

  if (wtkflg /= 0) then
   write(message,'(a,a,a,a)') ch10,&
&   ' getshell : BUG -',ch10,&
&   ' The calculated weights do not solve the linear system for all k-points.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
  end if
 end if

 if (wtkflg /= 0) then

  write(message,'(a,a,a,a)') ch10,&
&  ' getshell : BUG -',ch10,&
&  ' There is a problem with the finite difference expression of the ddk'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')

 else

  nshell = ishell

  write(message,'(a,a,a,a,a,a,a,i3,a,a,f16.7)') ch10,&
&  ' getshell : finite difference formula of Marzari and Vanderbilt',ch10,&
&  '            (see Marzari and Vanderbilt, PRB 56, 12847 (1997), Appendix B)',&
&  ch10,ch10,&
&  '            number of first neighbours  : ', neigh(1,1),ch10,&
&  '            weight : ',mvwtk(1,1)
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')

  if (nshell > 1) then
   is1 = neigh(1,1) + 1
   write(message,'(a,a,i3,a,a,f16.7)')ch10,&
&   '            number of second neighbours  : ', neigh(2,1),ch10,&
&   '            weight : ',mvwtk(is1,1)
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
  end if

  if (nshell > 2) then
   is1 = sum(neigh(1:2,1)) + 1
   write(message,'(a,a,i3,a,a,f16.7)')ch10,&
&   '            number of third neighbours  : ', neigh(3,1),ch10,&
&   '            weight : ',mvwtk(is1,1)
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
  end if

  if (nshell > 3) then
   is1 = sum(neigh(1:3,1)) + 1
   write(message,'(a,a,i3,a,a,f16.7)')ch10,&
&   '            number of fourth neighbours  : ', neigh(4,1),ch10,&
&   '            weight : ',mvwtk(is1,1)
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
  end if

  if (nshell > 4) then
   is1 = sum(neigh(1:4,1)) + 1
   write(message,'(a,a,i3,a,a,f16.7)')ch10,&
&   '            number of fifth neighbours  : ', neigh(5,1),ch10,&
&   '            weight : ',mvwtk(is1,1)
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
  end if

  if (nshell > 5) then
   is1 = sum(neigh(1:5,1)) + 1
   write(message,'(a,a,i3,a,a,f16.7)')ch10,&
&   '            number of sixth neighbours  : ', neigh(6,1),ch10,&
&   '            weight : ',mvwtk(is1,1)
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
  end if

 end if



!----------------------------------------------------------------------------

 deallocate(wtk3)


end subroutine getshell
!!***
