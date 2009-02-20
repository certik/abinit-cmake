!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcart2deloc
!! NAME
!! xcart2deloc
!!
!! FUNCTION
!!  calculate values of delocalized coordinates as a function of
!!  cartesian ones. First primitive internals, then B matrix, then F, then U
!!  then delocalized internals.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! angs= number of angles
!! bonds(2,2,nbond)=for a bond between iatom and jatom
!!              bonds(1,1,nbond) = iatom
!!              bonds(2,1,nbond) = icenter
!!              bonds(1,2,nbond) = jatom
!!              bonds(2,2,nbond) = irshift
!! carts(2,ncart)= index of total primitive internal, and atom (carts(2,:))
!! dihedrals(2,4,ndihed)=indexes to characterize dihedrals
!! dtset <type(dataset_type)>=all input variables for this dataset
!! nang(2,3,nang)=indexes to characterize angles
!! nbond=number of bonds
!! ncart=number of cartesian coordinates
!! ndihed= number of dihedrals
!! ninternal=nbond+nang+ndihed+ncart: number of internal coordinates
!! nrshift= dimension of rshift
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! rshift(3,nrshift)=shift in xred that must be done to find all neighbors of
!!                   a given atom within a given number of neighboring shells
!! xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! OUTPUT
!! bt_inv_matrix(3*(natom-1),3*natom)=Inverse of B^{T} matrix
!! deloc_int(3*(natom-1))=delocalized internal coordinates
!! prim_int(ninternal)=primitive internal coordinates
!!
!! SIDE EFFECTS
!! u_matrix(ninternal,3*(natom-1))=eigenvectors of BB^T matrix
!!
!! NOTES
!!
!! PARENTS
!!      deloc2xcart,delocint
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xcart2deloc(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&
& dtset,nrshift,rprimd,rshift,xcart,&
& bt_inv_matrix,u_matrix,deloc_int,prim_int)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_16geomoptim, except_this_one => xcart2deloc
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nang,nbond,ncart,ndihed,ninternal,nrshift
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,pointer :: angs(:,:,:),bonds(:,:,:),carts(:,:),dihedrals(:,:,:)
 real(dp),intent(in) :: rprimd(3,3),rshift(3,nrshift),xcart(3,dtset%natom)
 real(dp),intent(inout) :: u_matrix(ninternal,3*(dtset%natom-1))
 real(dp),intent(out) :: bt_inv_matrix(3*(dtset%natom-1),3*dtset%natom)
 real(dp),intent(out) :: deloc_int(3*(dtset%natom-1)),prim_int(ninternal)

!Local variables-------------------------------
!scalars
 integer :: i1,ibond
!arrays
 real(dp) :: b_matrix(ninternal,3*dtset%natom)
!no_abirules

! ******************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'DGEMV' :: dgemv
#endif

!DEBUG
!write (*,*) 'xcart2deloc : enter'
!do ibond=1,nbond
!do i1=1,2
!write (*,'(2I5)') bonds(:,i1,ibond)
!end do
!end do
!ENDDEBUG


 call calc_prim_int(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&
& dtset,nrshift,rprimd,rshift,xcart,prim_int)

 call calc_b_matrix(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&
& dtset,nrshift,rprimd,rshift,xcart,b_matrix)

 call calc_btinv_matrix(b_matrix,dtset,nang,nbond,ndihed,ncart,ninternal,&
& bt_inv_matrix,u_matrix)

!calculate value of delocalized internals


 call dgemv('T',ninternal,3*(dtset%natom-1),one,&
& u_matrix,ninternal,prim_int,1,zero,deloc_int,1)

!write (*,'(a,6E16.6)') 'xcart2deloc : deloc_int = ',  deloc_int

end subroutine xcart2deloc
!!***


!!****f* ABINIT/calc_btinv_matrix
!******************************************************************
!******************************************************************

!! NAME
!! calc_btinv_matrix
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine calc_btinv_matrix(b_matrix,dtset,nang,nbond,ndihed,ncart,ninternal,&
& bt_inv_matrix,u_matrix)

 use defs_basis
!
!note: bt_inv_matrix is inverse transpose of the delocalized
!coordinate B matrix. b_matrix is the primitive internal B matrix
!
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_16geomoptim, except_this_one => calc_btinv_matrix
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: ninternal,nang,nbond,ndihed,ncart
 real(dp),intent(in) :: b_matrix(ninternal,3*dtset%natom)
 real(dp),intent(out) :: bt_inv_matrix(3*(dtset%natom-1),3*dtset%natom)
 real(dp),intent(inout) :: u_matrix(ninternal,3*(dtset%natom-1))

!Local variables ------------------------------------
!scalars
 integer :: i1,i2,ii,info,iprim,lrwork,lwork,n1,n2,n3,n4,nb
!arrays
 real(dp) :: f_eigs(3*dtset%natom),f_matrix(3*dtset%natom,3*dtset%natom)
 real(dp) :: s_matrix(3*dtset%natom,3*dtset%natom)
 real(dp) :: s_red(3*dtset%natom,3*(dtset%natom-1))
 real(dp) :: u_matrix_old(ninternal,3*(dtset%natom-1))
 real(dp),allocatable :: rwork(:,:),tmp_fmat(:,:),tmp_zevec(:,:,:),tmpbmat(:,:)
 real(dp),allocatable :: work(:)

!******************************************************************
#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'DSYEV' :: dsyev
!DEC$ ATTRIBUTES ALIAS:'DGEMM' :: dgemm
#endif

!DEBUG
!write (*,*) 'T N 3*dtset%natom,3*dtset%natom,ninternal,one,',&
!&   ' b_matrix,ninternal,b_matrix,3*dtset%natom,zero,f_matrix,3*dtset%natom'
!write (*,*) 'T','N',3*dtset%natom,3*dtset%natom,ninternal,one,&
!&   'b_matrix',3*dtset%natom,'b_matrix',ninternal,zero,'f_matrix',3*dtset%natom
!ENDDEBUG


!DEBUG
!allocate (tmpbmat(nbond,3*dtset%natom))
!do iprim=1,nbond
!tmpbmat(iprim,:) = b_matrix(iprim,:)
!end do
!call dgemm('T','N',3*dtset%natom,3*dtset%natom,nbond,one,&
!&   tmpbmat,nbond,tmpbmat,nbond,zero,f_matrix,3*dtset%natom)
!write (*,*) 'F matrix = '
!do iprim=1,3*dtset%natom
!write (*,'(6E16.6)') f_matrix(:,iprim)
!end do
!deallocate (tmpbmat)
!allocate (tmpbmat(nbond+nang,3*dtset%natom))
!write (*,*) shape(b_matrix), shape(tmpbmat)
!do iprim=1,nbond+nang
!tmpbmat(iprim,:) = b_matrix(iprim,:)
!end do
!call dgemm('T','N',3*dtset%natom,3*dtset%natom,nbond+nang,one,&
!&   tmpbmat,nbond+nang,tmpbmat,nbond+nang,zero,f_matrix,3*dtset%natom)
!write (*,*) 'F matrix = '
!do iprim=1,3*dtset%natom
!write (*,'(6E16.6)') f_matrix(:,iprim)
!end do
!deallocate (tmpbmat)
!ENDDEBUG

!f matrix = B^{T} B
 call dgemm('T','N',3*dtset%natom,3*dtset%natom,ninternal,one,&
& b_matrix,ninternal,b_matrix,ninternal,zero,f_matrix,3*dtset%natom)
!DEBUG
!write (201,*) 'F matrix = '
!do iprim=1,3*dtset%natom
!write (201,'(6E16.6)') f_matrix(:,iprim)
!end do
!ENDDEBUG

!
!lwork = max(1,2*3*dtset%natom-1)
!lrwork = max(1,3*3*dtset%natom-2)
!
!allocate (work(lwork),rwork(2,lrwork))
!allocate (tmp_zevec(2,3*dtset%natom,3*dtset%natom))
!allocate (tmp_fmat(2,3*dtset%natom*(3*dtset%natom+1)/2))
!
!tmp_fmat(:,:) = zero
!ii=0
!do i1=1,3*dtset%natom
!do i2=i1,3*dtset%natom
!ii=ii+1
!tmp_fmat(1,ii) = f_matrix(i2,i1)
!write (*,'(a,2E16.6)') '  ', tmp_fmat(:,ii)
!end do
!write (*,*)
!end do
!
!call zhpev('V','L',3*dtset%natom,tmp_fmat,f_eigs,tmp_zevec,3*dtset%natom,work,rwork,info)
!
!
!s_matrix(:,:) = tmp_zevec(1,:,:)

 lwork = max(1,3*3*dtset%natom-1)
 allocate (work(lwork))
 s_matrix(:,:) = f_matrix(:,:)

 call dsyev('V','L',3*dtset%natom,s_matrix,3*dtset%natom,&
& f_eigs,work,lwork,info)

!write (*,*) 'info = ', info

!!DEBUG
!write (*,*) 'complex part? = '
!do iprim=1,3*dtset%natom
!write (*,'(6E16.6)') tmp_zevec(2,:,iprim)
!end do
!deallocate (tmp_zevec,tmp_fmat)
!
!write (201,*) 'S matrix = '
!do iprim=1,3*dtset%natom
!write (201,'(6E16.6)') s_matrix(:,iprim)
!end do
!write (201,'(a,6E16.6)')'f_eigs = '
!do ii=1,3*dtset%natom
!write(201,'(E16.6)',ADVANCE='NO') f_eigs(ii)
!end do
!write(201,*)
!ENDDEBUG

 if (abs(f_eigs(1)) + abs(f_eigs(2)) + abs(f_eigs(3)) > tol10 ) then
  write (*,*) 'Error: 3 lowest eigenvalues are not zero'
  write (*,*) '  internal coordinates do NOT span the full degrees of freedom !'
  write (*,'(6E16.6)') f_eigs
  stop
 end if
 if ( abs(f_eigs(4)) < tol10 ) then
  write (*,*) 'Error: fourth eigenvalue is zero'
  write (*,*) '  internal coordinates do NOT span the full degrees of freedom !'
  write (*,'(6E16.6)') f_eigs
  stop
 end if

!calculate U matrix from U = B * S_red * lambda^{-1/2}
 do ii=1,3*(dtset%natom-1)
  s_red(:,ii) = s_matrix(:,ii+3)/sqrt(f_eigs(ii+3))
 end do

 u_matrix_old(:,:) = u_matrix(:,:)

!DEBUG
!write (*,*) ' shapes ', shape(b_matrix)
!write (*,*) ' shapes ', shape(s_red)
!write (*,*) ' shapes ', shape(u_matrix)
!ENDDEBUG

 call dgemm('N','N',ninternal,3*(dtset%natom-1),3*dtset%natom,one,&
& b_matrix,ninternal,s_red,3*dtset%natom,zero,u_matrix,ninternal)


!DEBUG
!write (*,*) 'U^{T} matrix = '
!do iprim=1,ninternal
!write (*,'(6E16.6)') u_matrix(iprim,:)
!end do
!ENDDEBUG

!align eigenvectors, to preserve a form of continuity in convergences
!!!! eigenvalues are no longer in increasing order!!! but only s_red is reordered
!so that btinv is correct.
 call align_u_matrices(dtset,ninternal,u_matrix,u_matrix_old,s_matrix,f_eigs)

!calculate B_deloc^{-1} matrix for transformation of forces to deloc coord.
!(B^{T}_deloc)^{-1} = (B_deloc B^{T}_deloc)^{-1} B_deloc = lambda^{-3/2} S^{T} F
!= ( S lambda^{3/2} )^{T} F

!do ii=1,3*(dtset%natom-1)
!!s_red(:,ii) = s_matrix(:,ii+3)*sqrt(f_eigs(ii+3))
!s_red(:,ii) = s_matrix(:,ii+3)/f_eigs(ii+3)/sqrt(f_eigs(ii+3))
!end do
!
!
!call dgemm('T','N',3*(dtset%natom-1),3*dtset%natom,3*dtset%natom,one,&
!&   s_red,3*dtset%natom,f_matrix,3*dtset%natom,zero,&
!&   bt_inv_matrix,3*(dtset%natom-1))

!even better: B_deloc^{-1} = lambda^{-1/2} S^{T}
 do ii=1,3*(dtset%natom-1)
! s_red(:,ii) = s_matrix(:,ii+3)*sqrt(f_eigs(ii+3))
  bt_inv_matrix(ii,:) = s_matrix(:,ii+3)/sqrt(f_eigs(ii+3))
 end do


!DEBUG
!write (*,*) 'calc_btinv_matrix : BTinv matrix = '
!do iprim=1,3*dtset%natom
!write (*,'(6E16.6)') bt_inv_matrix(:,iprim)
!end do
!ENDDEBUG


end subroutine calc_btinv_matrix

!******************************************************************
!******************************************************************

 subroutine align_u_matrices(dtset,ninternal,u_matrix,u_matrix_old,s_matrix,f_eigs)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: ninternal
 type(dataset_type), intent(in) :: dtset
 real(dp), intent(in) :: u_matrix_old(ninternal,3*(dtset%natom-1))
 real(dp), intent(inout) :: u_matrix(ninternal,3*(dtset%natom-1))
 real(dp), intent(inout) :: s_matrix(3*dtset%natom,3*dtset%natom)
 real(dp), intent(inout) :: f_eigs(3*dtset%natom)

!Local variables ------------------------------
 integer :: ii, iint1, iint2, imax, lwork,info
 integer :: eigv_flag(3*(dtset%natom-1))
 integer :: eigv_ind(3*(dtset%natom-1))
 real(dp) :: ss, smax, tmpf(3*dtset%natom)
 real(dp) :: tmpu(ninternal,3*(dtset%natom-1)), tmps(3*dtset%natom,3*dtset%natom)
 real(dp) :: overlap(3*(dtset%natom-1),3*(dtset%natom-1))
 real(dp) :: a_matrix(3*(dtset%natom-1),3*(dtset%natom-1))
 real(dp), allocatable :: work(:)

!******************************************************************
!call dgemm('T','N',3*(dtset%natom-1),3*(dtset%natom-1),ninternal,one,&
!&   u_matrix,ninternal,u_matrix_old,ninternal,zero,&
!&   overlap,3*(dtset%natom-1))
!
!lwork = max(1,3*3*dtset%natom-1)
!allocate (work(lwork))
!a_matrix(:,:) = overlap(:,:)
!
!call dsyev('V','L',3*(dtset%natom-1),s_matrix,3*(dtset%natom-1),&
!& f_eigs,work,lwork,info)


 eigv_flag(:) = 0
 eigv_ind(:) = 0

!just permit a change in sign
 do iint1=1,3*(dtset%natom-1)
  ss = zero
  do ii=1,ninternal
   ss = ss + u_matrix_old(ii,iint1)*u_matrix(ii,iint1)
  end do
  if (ss < -tol12) then
   imax = -iint1
  else
   imax = iint1
  end if
  eigv_ind(iint1) = imax
  eigv_flag(abs(imax)) = 1
 end do

!permit a change in order or sign
!do iint1=1,3*(dtset%natom-1)
!smax = zero
!do iint2=1,3*(dtset%natom-1)
!if (eigv_flag(iint2) /= 0) cycle
!
!ss = zero
!do ii=1,ninternal
!ss = ss + u_matrix_old(ii,iint2)*u_matrix(ii,iint1)
!end do
!if (ss > smax) then
!imax = iint2
!smax = ss
!else if (-ss > smax) then
!imax = -iint2
!smax = abs(ss)
!end if
!
!end do
!eigv_ind(iint1) = imax
!eigv_flag(abs(imax)) = 1
!end do

!DEBUG
!write (*,*) 'align_u_matrices : eigv_flag = '
!write (*,'(8I4)') eigv_flag
!write (*,*) 'align_u_matrices : eigv_ind = '
!write (*,'(8I4)') eigv_ind
!ENDDEBUG

!DEBUG
!write (*,*) 'align_u_matrices : u_matrix, u_matrix_old before '
!do iint1=1,ninternal
!write (*,'(6E16.6)') u_matrix(iint1,:), u_matrix_old(iint1,:)
!end do
!ENDDEBUG

 tmpu(:,:) = u_matrix
 tmps(:,:) = s_matrix
 tmpf(:) = f_eigs
!exchange eigenvectors...
 do iint1=1,3*(dtset%natom-1)
  ss = one
  if (eigv_ind(iint1) < 0) ss = -one

  imax = abs(eigv_ind(iint1))

  tmpu(:,imax) = ss*u_matrix(:,iint1)

  tmps(:,imax+3) = ss*s_matrix(:,iint1+3)

  tmpf(imax+3) = f_eigs(iint1+3)
 end do

!DEBUG
!write (*,*) 'align_u_matrices : u_matrix, u_matrix_old after '
!do iint1=1,ninternal
!write (*,'(6E16.6)') tmpu(iint1,:), u_matrix_old(iint1,:)
!end do
!write (*,*) 'align_u_matrices : s_matrix, s_matrix_old after '
!do iint1=1,3*dtset%natom
!write (*,'(6E16.6)') tmps(iint1,4:3*dtset%natom),s_matrix(iint1,4:3*dtset%natom)
!end do
!ENDDEBUG

 u_matrix(:,:) = tmpu(:,:)
 s_matrix(:,:) = tmps(:,:)
 f_eigs(:) = tmpf(:)


end subroutine align_u_matrices


!******************************************************************
!******************************************************************

 subroutine cross_product(r1,r2,cp)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
 real(dp), intent(in) :: r1(3),r2(3)
 real(dp), intent(out) :: cp(3)

!Local variables ------------------------------

!******************************************************************
 cp(1) = r1(2)*r2(3)-r1(3)*r2(2)
 cp(2) = r1(3)*r2(1)-r1(1)*r2(3)
 cp(3) = r1(1)*r2(2)-r1(2)*r2(1)

end subroutine cross_product
!!***
