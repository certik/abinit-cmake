!{\src2tex{textfont=tt}}
!!****f* ABINIT/spline_paw_fncs
!! NAME
!! spline_paw_fncs
!!
!! FUNCTION
!! Compute radial PAW functions and their derivatives on a set of points in the PAW sphere.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (JJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! integer :: nnl : number of nl PAW basis functions in set
!! integer :: npts : number of points to perform fits on 
!! real(dp) :: points(npts) : SORTED vector of points to perform fits on
!! type(pawrad_type) :: pawrad : paw radial mesh data
!! type(pawtab_type) :: pawtab : paw wavefunctions around each type of atom
!!
!! OUTPUT
!! real(dp) :: phi(npts,nnl), dphi(npts,nnl), tphi(npts,nnl), dtphi(npts,nnl) : PAW functions
!!             phi, tphi and their radial derivatives evaluated at the input points by spline fits
!!
!! NOTES
!! The PAW basis functions are defined by $<r|\phi_i>=(u_i(r)/r)S_{lm}(\hat{r})$, evaluated on a radial
!! grid. This subroutine computes $u(r)$ and $d u(r)/dr$ on the set of points provided on input. Typically
!! the input points will be the the find grid points in the PAW sphere. They are presumed to be sorted
!! already on input.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine spline_paw_fncs(dphi,dtphi,nnl,npts,pawrad,pawtab,points,phi,tphi)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nnl,npts
!arrays
 real(dp),intent(in) :: points(npts)
 real(dp),intent(out) :: phi(npts,nnl),dphi(npts,nnl),tphi(npts,nnl),dtphi(npts,nnl)
 type(pawrad_type),intent(in) :: pawrad
 type(pawtab_type),intent(in) :: pawtab

!Local variables-------------------------------
!scalars
 integer :: inl,nder
 real(dp) :: ybcbeg, ybcend
!arrays
 real(dp),allocatable :: der(:),diag(:),ypp(:)

! ************************************************************************

!DEBUG
!write(*,*)' spline_paw_fncs : enter'
!ENDDEBUG

 nder = 1 ! compute first derivative in calls to nderiv_gen below
 allocate(der(pawrad%mesh_size),ypp(pawrad%mesh_size),diag(pawrad%mesh_size))
 do inl = 1, nnl

! spline phi onto points
  ypp(:) = zero; diag(:) = zero; ybcbeg = 0.0; ybcend = 0.0;
  call spline(pawrad%rad,pawtab%phi(:,inl),pawrad%mesh_size,ybcbeg,ybcend,ypp,diag)
  call splint(pawrad%mesh_size,pawrad%rad,pawtab%phi(:,inl),ypp,npts,points,phi(:,inl))

! next spline d phi/dr onto points
! need derivative of phi with respect to radius
  der(:) = zero
  call nderiv_gen(der,pawtab%phi(:,inl),nder,pawrad)
  ypp(:) = zero; diag(:) = zero; ybcbeg = 0.0; ybcend = 0.0;
  call spline(pawrad%rad,der,pawrad%mesh_size,ybcbeg,ybcend,ypp,diag)
  call splint(pawrad%mesh_size,pawrad%rad,der,ypp,npts,points,dphi(:,inl))

! next splint tphi onto points
  ypp(:) = zero; diag(:) = zero;
  call spline(pawrad%rad,pawtab%tphi(:,inl),pawrad%mesh_size,ybcbeg,ybcend,ypp,diag)
  call splint(pawrad%mesh_size,pawrad%rad,pawtab%tphi(:,inl),ypp,npts,points,tphi(:,inl))

! finally spline d tphi/dr onto points
! need derivative of tphi with respect to radius
  der(:) = zero
  call nderiv_gen(der,pawtab%tphi(:,inl),nder,pawrad)
  ypp(:) = zero; diag(:) = zero; ybcbeg = 0.0; ybcend = 0.0;
  call spline(pawrad%rad,der,pawrad%mesh_size,ybcbeg,ybcend,ypp,diag)
  call splint(pawrad%mesh_size,pawrad%rad,der,ypp,npts,points,dtphi(:,inl))

 end do ! end loop over nnl basis functions

 deallocate(der,ypp,diag)

!DEBUG
!write(6,*)' spline_paw_fncs : exit '
!stop
!ENDDEBUG

 end subroutine spline_paw_fncs
!!***
