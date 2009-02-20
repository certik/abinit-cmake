!{\src2tex{textfont=tt}}
!!****f* ABINIT/deloc2xcart
!! NAME
!! deloc2xcart
!!
!! FUNCTION
!!  determine the cartesian coordinates which correspond to the given
!!  values of the delocalized coordinates. The relationship is non-linear,
!!  so use an iterative scheme, as in Baker JCP .105. 192 (1996)
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
!! deloc_int(3*(natom-1))=delocalized internal coordinates
!! dtset <type(dataset_type)>=all input variables for this dataset
!! nang(2,3,nang)=indexes to characterize angles
!! nbond=number of bonds
!! ncart=number of cartesian directions (used for constraints)
!! ndihed= number of dihedrals
!! ninternal=nbond+nang+ndihed+ncart: number of internal coordinates
!! nrshift= dimension of rshift
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! rshift(3,nrshift)=shift in xred that must be done to find all neighbors of
!!                   a given atom within a given number of neighboring shells
!!
!! OUTPUT
!! bt_inv_matrix(3*(natom-1),3*natom)=inverse of transpose of B matrix
!!
!! SIDE EFFECTS
!! u_matrix(ninternal,3*(natom-1))=eigenvectors of G = BB^T matrix
!! xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! NOTES
!!
!! PARENTS
!!      delocint
!!
!! CHILDREN
!!      dgemv,xcart2deloc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine deloc2xcart(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,nrshift,&
& dtset,rprimd,rshift,xcart,&
& deloc_int,btinv,u_matrix)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_16geomoptim, except_this_one => deloc2xcart
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nang,nbond,ncart,ndihed,ninternal,nrshift
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,pointer :: angs(:,:,:),bonds(:,:,:),carts(:,:),dihedrals(:,:,:)
 real(dp),intent(in) :: deloc_int(3*(dtset%natom-1)),rprimd(3,3)
 real(dp),intent(in) :: rshift(3,nrshift)
 real(dp),intent(inout) :: u_matrix(ninternal,3*(dtset%natom-1))
 real(dp),intent(inout) :: xcart(3,dtset%natom)
 real(dp),intent(out) :: btinv(3*(dtset%natom-1),3*dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: enough,i1,iatom,ibond,ideloc,ii,iint,iiter,iline,iprim,itest,niter
 integer :: nline
 real(dp) :: alpha,mindiff,mindiff1,mindiff2,mindiff3,minmix,minmix1,minmix2
 real(dp) :: minmix3,mix,s1,s2,tot_diff
!arrays
 real(dp) :: btinv_tmp(3*(dtset%natom-1),3*dtset%natom)
 real(dp) :: cgmat(3*dtset%natom,3*dtset%natom),cgmat_cgrad(3*dtset%natom)
 real(dp) :: cgrad(3*dtset%natom),cgrad_old(3*dtset%natom)
 real(dp) :: deloc_int_now(3*(dtset%natom-1)),prim_int(ninternal)
 real(dp) :: tmpxcart(3*dtset%natom),tmpxcart0(3*dtset%natom)
 real(dp) :: u_old(ninternal,3*(dtset%natom-1)),xdeloc_diff(3*(dtset%natom-1))
 real(dp),allocatable :: line_diff(:),line_mix(:)
!no_abirules

! ******************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'DGEMV' :: dgemv
#endif

!DEBUG
!write (*,*) 'deloc2xcart : enter'
!do ibond=1,nbond
!do i1=1,2
!write (*,'(2I5)') bonds(:,i1,ibond)
!end do
!end do
!ENDDEBUG

 niter = 200
 tmpxcart = reshape(xcart,(/3*dtset%natom/))
 nline = dtset%userib
 allocate (line_diff(nline+3),line_mix(nline+3))


 cgrad_old(:) = zero
 cgrad(:) = zero
 mix = dtset%userrb

 do iiter=1,niter
  if (iiter > 1) then
   tmpxcart0 = tmpxcart
   u_old(:,:) = u_matrix(:,:)
!  write (300,*) 'line minimization, iiter = ', iiter
   minmix1 = dtset%userrb
   mindiff1 = 1.0d10
   do iline=1,nline
!   do iline=-nline,nline
    mix = iline*dtset%userrb
    tmpxcart = tmpxcart0 + mix*cgrad

!   DEBUG
!   tmpxcart = tmpxcart0
!   tmpxcart(dtset%useric) = tmpxcart(dtset%useric)*(one + mix)
!   ENDDEBUG
    xcart = reshape(tmpxcart,(/3,dtset%natom/))
    u_matrix(:,:) = u_old(:,:)
    call xcart2deloc(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&
&    dtset,nrshift,rprimd,rshift,xcart,&
&    btinv_tmp,u_matrix,deloc_int_now,prim_int)
    tot_diff = sum(abs(deloc_int(:) - deloc_int_now(:)))
    write (300,'(a,I5,3E16.6)') '   ', iline,mix, tmpxcart(dtset%useric),tot_diff
    write (301,*) ' deloc_int_now   ', iline,mix, minmix1,deloc_int_now
    write (302,*) '  u_matrix  '
!   do ii=1,3*(dtset%natom-1)
!   write (400+ii,'(I5)',ADVANCE='NO') iline
!   do iint=1,ninternal
!   write (400+ii,'(E16.6)',ADVANCE='NO') u_matrix(iint,ii)
!   end do
!   write (400+ii,*)
!   end do
    if (tot_diff < mindiff1) then
     minmix1 = mix
     mindiff1 = tot_diff
    end if
   end do
   u_matrix(:,:) = u_old(:,:)
   tmpxcart = tmpxcart0
   mix = minmix1
  end if

! mix = dtset%userrb
  write (*,*) 'mix used is ', mix
  tmpxcart(:) = tmpxcart(:) + mix*cgrad(:)
  xcart = reshape(tmpxcart,(/3,dtset%natom/))
  call xcart2deloc(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&
&  dtset,nrshift,rprimd,rshift,xcart,&
&  btinv_tmp,u_matrix,deloc_int_now,prim_int)
! update the BT^{-1} matrix?
! if (iiter == 1) then
  btinv(:,:) = btinv_tmp(:,:)
! end if

! call xcart2deloc_fixb(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&
! & dtset,nrshift,rprimd,rshift,xcart,&
! & u_matrix,deloc_int_now,prim_int)

  xdeloc_diff(:) = deloc_int(:) - deloc_int_now(:)
  write (*,*) 'xdeloc_diff = '
  write (*,*) '     ', xdeloc_diff

  write (*,*) 'btinv = '
  do iatom=1,3*dtset%natom
   write (*,'(6E16.6)') btinv(:,iatom)
  end do

  tot_diff = sum(abs(xdeloc_diff))
  write (*,*) 'deloc2xcart : Iteration ', iiter, tot_diff
  write (*,*) '  tmpxcart = '
  write (*,'(3E16.6)') tmpxcart
  if (tot_diff < tol10) then
   write (*,*) 'deloc2xcart : converged '
   exit
  end if


  cgrad_old(:) = cgrad(:)

! gradient vector = btinv^{T} * xdeloc_diff
  call dgemv('T',3*(dtset%natom-1),3*dtset%natom,one,&
&  btinv,3*(dtset%natom-1),xdeloc_diff,1,zero,cgrad,1)
  write (*,*) 'deloc2xcart :cgrad'
  write (*,'(6E16.6)') cgrad

! ! calculate matrix cgmat: should be square, so we use Btinv * Btinv^{T}
! call dgemm('T','N',3*dtset%natom,3*dtset%natom,3*(dtset%natom-1),one,&
! &   btinv,3*(dtset%natom-1),btinv,3*(dtset%natom-1),zero,&
! &   cgmat,3*dtset%natom)
! call dgemv('N',3*dtset%natom,3*dtset%natom,one,&
! &  cgmat,3*dtset%natom,cgrad,1,zero,cgmat_cgrad,1)
! 
! write (*,*) 'deloc2xcart : cgmat,  cgmat_cgrad'
! write (*,'(6E16.6)')  cgmat, cgmat_cgrad

! call dgemv('T',3*(dtset%natom-1),3*dtset%natom,alpha,&
! &   btinv,3*(dtset%natom-1),xdeloc_diff,1,one,tmpxcart,1)
 end do
!end iiter do

 call xcart2deloc(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&
& dtset,nrshift,rprimd,rshift,xcart,&
& btinv,u_matrix,deloc_int_now,prim_int)
 write (*,*) 'deloc_int, after convergence = '
 write (*,*) '     ', deloc_int_now

 xdeloc_diff(:) = deloc_int(:) - deloc_int_now(:)
 write (*,*) 'xdeloc_diff = '
 write (*,*) '     ', xdeloc_diff
 tot_diff = sum(abs(xdeloc_diff))

 write (*,*) 'Primitive internal coordinate values:'
 do iprim=1,ninternal
  if (iprim <= nbond) then
   write (*,*) iprim, prim_int(iprim)
  else
   write (*,*) iprim, prim_int(iprim), prim_int(iprim)/pi*180.0_dp
  end if
 end do

 if (iiter == niter+1) then
  write (*,*) 'deloc2xcart : Error, xcart not converged in ', &
&  niter, 'iterations ', tot_diff
  stop
 end if

end subroutine deloc2xcart
!!***

!!****f* ABINIT/get_mix
!function get_mix(nline,iline,line_diff,line_mix)
!! NAME
!! get_mix
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

 function get_mix(x1,y1,x2,y2,x3,y3)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
 real(dp),intent(in) :: y1,y2,y3,x1,x2,x3
 real(dp) :: get_mix

!Local variables-------------------------------
!integer :: nline, iline
!real(dp) :: line_diff(nline+3),line_mix(nline+3)
!real(dp) :: get_mix
!scalars
 integer :: ii,imid
 real(dp) :: a,b,c,mix,sx,sx2,sy

! *************************************************************************

!make sure we dont reuse the same point
!imid=min(max(2,int(iline/2)),iline-2)
!y1=line_diff(iline-1)
!y2=line_diff(1)
!y3=line_diff(imid)
!x1=line_mix(iline-1)
!x2=line_mix(1)
!x3=line_mix(imid)


 a= (y3-y1)-(y2-y1)*(x3-x1)/(x2-x1)
 a= a / ( (x3**2-x1**2) - (x2**2-x1**2)*(x3-x1)/(x2-x1) )
 b = ( (y2-y1) - a*(x2**2-x1**2) ) / (x2-x1)

!sx = zero
!sx2 = zero
!do ii=1,iline-1
!sx = sx + line_mix(ii)
!sx2 = sx2 + line_mix(ii)**2
!sy = sy + line_diff(ii)
!end do
!y1=line_diff(1)
!y2=line_diff(iline-1)
!x1=line_mix(1)
!x2=line_mix(iline-1)
!
!a = y2-y1 - x2*(sy-y1*(iline-1))/sx
!a = a / (x2**2  - x2*sx2/sx)
!b = (sy-y1 - a*sx2) / sx

 write (*,*) 'get_mix : a,b = ', a, b
!y1*a is positive -> parabola goes the right way
!y1 should always be > 0
 if (y1*a > tol10) then
  mix = -b/two/a
 else
  write (*,*) 'got inverted parabola: use linear interpolation'
  mix = 999.0_dp
! stop
! mix = x1-y1*(x2-x1)/(y2-y1)
 end if

 if (mix < -tol12) then
  write (*,*) 'Error: negative mix'
! mix = half*(x1+x2)
  stop
 end if

 get_mix = mix

end function get_mix
!!***
