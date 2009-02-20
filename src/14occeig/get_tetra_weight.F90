!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_tetra_weight
!! NAME
!! get_tetra_weight
!!
!! FUNCTION
!! calculate integration weights and their derivatives
!! from Blochl et al PRB 49 16223
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (MVer,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! eigen_in(nkpt)=eigenenergies for each k point
!! enemin=minimal energy for DOS
!! enemax=maximal energy for DOS
!! nene=number of energies for DOS
!! nkpt=number of irreducible kpoints
!! nkpt_fullbz=number of kpoints in full brillouin zone
!! max_occ=maximal occupation number (2 for nsppol=1, 1 for nsppol=2)
!! ntetra=number of tetrahedra
!! mtetra=maximum number of tetrahedra
!! tetra_full(4,2,mtetra)=for each irred tetrahedron, the list of k point vertices
!!   1 -> irred kpoint   2 -> fullkpt
!! tetra_mult(mtetra)=for each irred tetrahedron, its multiplicity
!! vv = volume of one tetrahedron in reciprocal space
!!
!! OUTPUT
!!  tweight(nkpt,nene) = integration weights for each irred kpoint from all adjacent tetrahedra
!!  dtweightde(nkpt,nene) = derivative of tweight wrt energy
!!
!! PARENTS
!!      elphon,mkphdos,tetrahedron
!!
!! CHILDREN
!!      sort_dp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine get_tetra_weight(eigen_in,enemin,enemax,&
&             max_occ,mtetra,nene,nkpt,ntetra,rcvol,tetra_full,&
&            tetra_mult,tweight,dtweightde,vv)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mtetra,nene,nkpt,ntetra
 real(dp),intent(in) :: enemax,enemin,max_occ,rcvol,vv
!arrays
 integer,intent(in) :: tetra_full(4,2,mtetra),tetra_mult(mtetra)
 real(dp),intent(in) :: eigen_in(nkpt)
 real(dp),intent(out) :: dtweightde(nkpt,nene),tweight(nkpt,nene)

!Local variables-------------------------------
!  needed for gaussian replacement of Dirac functions
!  the three coefficients of the DOS as quadratic form,
!    in the interval [eig(ikpt-1), eig(ikpt)]
!    for ikpt = 1 we add a point below eigen(1) which doesnt
!    contribute to the DOS in any tetrahedron
!scalars
 integer :: ieps,ikpt,itetra,nn1,nn2,nn3,nn4
 real(dp) :: cc,cc1,cc2,cc3,dcc1de,dcc2de,dcc3de,dccde,deltaene,eigmid,eps
 real(dp) :: epsilon21,epsilon31,epsilon32,epsilon41,epsilon42,epsilon43
 real(dp) :: gau_prefactor,gau_width,gau_width2,inv_epsilon21,inv_epsilon31
 real(dp) :: inv_epsilon32,inv_epsilon41,inv_epsilon42,inv_epsilon43
 real(dp) :: prefac_tetr,tmp,tmp1,tmp2,volconst,volconst_mult
!arrays
 integer :: ind_dum(4),ind_eigen_tmp(nkpt+1)
 real(dp) :: eigen_1tetra(4),eigen_tmp(nkpt+1)

! *********************************************************************

!DEBUG
!write (*,*) ' get_tetra_weight : enter'
!write (*,*) '   eigen_in ', eigen_in
!write (*,*) ' vv, rcvol, vv/rcvol = ', vv, rcvol, vv/rcvol
!write (*,*) '  nene, nkpt, ntetra ',  nene, nkpt, ntetra
!write (*,*) '  nene,nkpt,max_occ,ntetra,rcvol ', &
!&              nene,nkpt,max_occ,ntetra,rcvol
!ENDDEBUG

!DEBUG
!write (*,*) 'tetra_full = '
!write (*,*) 'tetra_full(:,:,2) = ', tetra_full(:,:,2)
!do itetra=1,ntetra
!write (*,'(4(I7,1x))' ) tetra_full(:,:,itetra)
!write (*,*)
!end do
!write (*,*) 'tetra_mult = '
!write (*,'(8(I7,1x))' ) (tetra_mult(itetra),itetra=1,ntetra)
!ENDDEBUG

 tweight(:,:) = zero
 dtweightde(:,:) = zero

 volconst = vv/rcvol/four
 if (nene <= 1) then
  write (*,*) 'get_tetra_weight : Error: nene must be at least 2'
  stop
 else
  deltaene = (enemax-enemin) / (nene-1)
 end if
!
!for each tetrahedron
!
 do itetra=1,ntetra
  volconst_mult = max_occ*volconst*float(tetra_mult(itetra))

! Here we need the original ordering to reference the correct irred kpoints
  eigen_1tetra(1) = eigen_in(tetra_full(1,1,itetra))
  eigen_1tetra(2) = eigen_in(tetra_full(2,1,itetra))
  eigen_1tetra(3) = eigen_in(tetra_full(3,1,itetra))
  eigen_1tetra(4) = eigen_in(tetra_full(4,1,itetra))
  ind_dum(1) = tetra_full(1,1,itetra)
  ind_dum(2) = tetra_full(2,1,itetra)
  ind_dum(3) = tetra_full(3,1,itetra)
  ind_dum(4) = tetra_full(4,1,itetra)

  call sort_dp(4,eigen_1tetra,ind_dum,tol14)

! DEBUG
! if (eigen_1tetra(1) >= -0.09 .and. eigen_1tetra(2) <= -0.04) then
! write (*,*) ' get_tetra_weight : eigen_1tetra(1:4)', eigen_1tetra(1:4)
! write (*,*) '  itetra, ind_dum = ', itetra, ind_dum
! end if
! ENDDEBUG

! 
! all notation are from Blochl PRB 49 16223 Appendix B
! 
  epsilon21 = eigen_1tetra(2)-eigen_1tetra(1)
  epsilon31 = eigen_1tetra(3)-eigen_1tetra(1)
  epsilon41 = eigen_1tetra(4)-eigen_1tetra(1)
  epsilon32 = eigen_1tetra(3)-eigen_1tetra(2)
  epsilon42 = eigen_1tetra(4)-eigen_1tetra(2)
  epsilon43 = eigen_1tetra(4)-eigen_1tetra(3)
  inv_epsilon21 = zero
  inv_epsilon31 = zero
  inv_epsilon41 = zero
  inv_epsilon32 = zero
  inv_epsilon42 = zero
  inv_epsilon43 = zero
  if (epsilon21 > tol6) then
   inv_epsilon21 = one / epsilon21
  end if
  if (epsilon31 > tol6) then
   inv_epsilon31 = one / epsilon31
  end if
  if (epsilon41 > tol6) then
   inv_epsilon41 = one / epsilon41
  end if
  if (epsilon32 > tol6) then
   inv_epsilon32 = one / epsilon32
  end if
  if (epsilon42 > tol6) then
   inv_epsilon42 = one / epsilon42
  end if
  if (epsilon43 > tol6) then
   inv_epsilon43 = one / epsilon43
  end if

  nn1 = int((eigen_1tetra(1)-enemin)/deltaene)+1
  nn2 = int((eigen_1tetra(2)-enemin)/deltaene)+1
  nn3 = int((eigen_1tetra(3)-enemin)/deltaene)+1
  nn4 = int((eigen_1tetra(4)-enemin)/deltaene)+1

! DEBUG
! write (*,*) 'deltaene = ', deltaene
! write (*,*) 'nene,enemin,enemax,eigen_1tetra(:) = ',&
! &              nene,enemin,enemax,eigen_1tetra(:)
! write (*,*) 'nn1,nn2,nn3,nn4 = ',&
! &              nn1,nn2,nn3,nn4
! ENDDEBUG

  nn1 = max(1,nn1)
  nn1 = min(nn1,nene)
  nn2 = max(1,nn2)
  nn2 = min(nn2,nene)
  nn3 = max(1,nn3)
  nn3 = min(nn3,nene)
  nn4 = max(1,nn4)
  nn4 = min(nn4,nene)

  eps = enemin+nn1*deltaene
! 
! interval enemin < eps < e1 nothing to do
! 
! 
! interval e1 < eps < e2
! 
  do ieps=nn1+1,nn2
   cc = volconst_mult*(eps-eigen_1tetra(1))*(eps-eigen_1tetra(1))*(eps-eigen_1tetra(1)) &
&   *inv_epsilon21*inv_epsilon31*inv_epsilon41
   tweight(ind_dum(1),ieps) = tweight(ind_dum(1),ieps) + &
&   cc*(four-(eps-eigen_1tetra(1))*(inv_epsilon21+inv_epsilon31+inv_epsilon41))
   tweight(ind_dum(2),ieps) = tweight(ind_dum(2),ieps) + &
&   cc*(eps-eigen_1tetra(1))*inv_epsilon21
   tweight(ind_dum(3),ieps) = tweight(ind_dum(3),ieps) + &
&   cc*(eps-eigen_1tetra(1))*inv_epsilon31
   tweight(ind_dum(4),ieps) = tweight(ind_dum(4),ieps) + &
&   cc*(eps-eigen_1tetra(1))*inv_epsilon41

   dccde = three*volconst_mult*(eps-eigen_1tetra(1))*(eps-eigen_1tetra(1)) &
&   *inv_epsilon21*inv_epsilon31*inv_epsilon41
   dtweightde(ind_dum(1),ieps) = dtweightde(ind_dum(1),ieps) + &
&   dccde*(four-(eps-eigen_1tetra(1))*(inv_epsilon21+inv_epsilon31+inv_epsilon41)) &
&   -cc*(inv_epsilon21+inv_epsilon31+inv_epsilon41)
   dtweightde(ind_dum(2),ieps) = dtweightde(ind_dum(2),ieps) + &
&   dccde*(eps-eigen_1tetra(1))*inv_epsilon21 + cc*inv_epsilon21
   dtweightde(ind_dum(3),ieps) = dtweightde(ind_dum(3),ieps) + &
&   dccde*(eps-eigen_1tetra(1))*inv_epsilon31 + cc*inv_epsilon31
   dtweightde(ind_dum(4),ieps) = dtweightde(ind_dum(4),ieps) + &
&   dccde*(eps-eigen_1tetra(1))*inv_epsilon41 + cc*inv_epsilon41

   eps = eps + deltaene
  end do
! DEBUG
! write (*,*) 'eps, enemin+nn2*deltaene ', eps, enemin+nn2*deltaene
! eps = enemin+nn2*deltaene
! ENDDEBUG
! 
! interval e2 < eps < e3
! 
  do ieps=nn2+1,nn3
   cc1 = volconst_mult*(eps-eigen_1tetra(1))*(eps-eigen_1tetra(1))&
&   *inv_epsilon31*inv_epsilon41
   cc2 = volconst_mult*(eps-eigen_1tetra(1))*(eps-eigen_1tetra(2))*(eigen_1tetra(3)-eps)&
&   *inv_epsilon41*inv_epsilon32*inv_epsilon31
   cc3 = volconst_mult*(eps-eigen_1tetra(2))*(eps-eigen_1tetra(2))*(eigen_1tetra(4)-eps)&
&   *inv_epsilon42*inv_epsilon32*inv_epsilon41
   tweight(ind_dum(1),ieps) = tweight(ind_dum(1),ieps) + &
&   cc1 + (cc1+cc2)*(eigen_1tetra(3)-eps)*inv_epsilon31 + &
&   (cc1+cc2+cc3)*(eigen_1tetra(4)-eps)*inv_epsilon41
   tweight(ind_dum(2),ieps) = tweight(ind_dum(2),ieps) + &
&   cc1+cc2+cc3+(cc2+cc3)*(eigen_1tetra(3)-eps)*inv_epsilon32 +&
&   cc3*(eigen_1tetra(4)-eps)*inv_epsilon42
   tweight(ind_dum(3),ieps) = tweight(ind_dum(3),ieps) + &
&   (cc1+cc2)*(eps-eigen_1tetra(1))*inv_epsilon31 + &
&   (cc2+cc3)*(eps-eigen_1tetra(2))*inv_epsilon32
   tweight(ind_dum(4),ieps) = tweight(ind_dum(4),ieps) + &
&   (cc1+cc2+cc3)*(eps-eigen_1tetra(1))*inv_epsilon41 + &
&   cc3*(eps-eigen_1tetra(2))*inv_epsilon42

   dcc1de = two*volconst_mult*(eps-eigen_1tetra(1))*inv_epsilon31*inv_epsilon41
   dcc2de = volconst_mult*inv_epsilon41*inv_epsilon32*inv_epsilon31*&
&   (-(eps-eigen_1tetra(1))*(eps-eigen_1tetra(2)) &
&   +(eps-eigen_1tetra(1))*(eigen_1tetra(3)-eps) &
&   +(eps-eigen_1tetra(2))*(eigen_1tetra(3)-eps))
   dcc3de = volconst_mult*inv_epsilon42*inv_epsilon32*inv_epsilon41*&
&   (two*(eps-eigen_1tetra(2))*(eigen_1tetra(4)-eps) &
&   -(eps-eigen_1tetra(2))*(eps-eigen_1tetra(2)))
   dtweightde(ind_dum(1),ieps) = dtweightde(ind_dum(1),ieps) + &
&   dcc1de+(dcc1de+dcc2de)*(eigen_1tetra(3)-eps)*inv_epsilon31-(cc1+cc2)*inv_epsilon31&
&   +(dcc1de+dcc2de+dcc3de)*(eigen_1tetra(4)-eps)*inv_epsilon41-(cc1+cc2+cc3)*inv_epsilon41
   dtweightde(ind_dum(2),ieps) = dtweightde(ind_dum(2),ieps) + &
&   dcc1de+dcc2de+dcc3de+(dcc2de+dcc3de)*(eigen_1tetra(3)-eps)*inv_epsilon32&
&   -(cc2+cc3)*inv_epsilon32 +dcc3de*(eigen_1tetra(4)-eps)*inv_epsilon42&
&   -cc3*inv_epsilon42
   dtweightde(ind_dum(3),ieps) = dtweightde(ind_dum(3),ieps) + &
&   (dcc1de+dcc2de)*(eps-eigen_1tetra(1))*inv_epsilon31 + &
&   (cc1+cc2)*inv_epsilon31+(dcc2de+dcc3de)*(eps-eigen_1tetra(2))*inv_epsilon32 + &
&   (cc2+cc3)*inv_epsilon32
   dtweightde(ind_dum(4),ieps) = dtweightde(ind_dum(4),ieps) + &
&   (dcc1de+dcc2de+dcc3de)*(eps-eigen_1tetra(1))*inv_epsilon41 + &
&   (cc1+cc2+cc3)*inv_epsilon41+dcc3de*(eps-eigen_1tetra(2))*inv_epsilon42 + &
&   cc3*inv_epsilon42
   eps = eps + deltaene
  end do
! DEBUG
! write (*,*) 'eps, enemin+nn3*deltaene ', eps, enemin+nn3*deltaene
! eps = enemin+nn3*deltaene
! ENDDEBUG
! 
! interval e3 < eps < e4
! 
  do ieps=nn3+1,nn4
   cc = volconst_mult*(eigen_1tetra(4)-eps)*(eigen_1tetra(4)-eps)*(eigen_1tetra(4)-eps)&
&   *inv_epsilon41*inv_epsilon42*inv_epsilon43
   tweight(ind_dum(1),ieps) = tweight(ind_dum(1),ieps) + &
&   volconst_mult - cc*(eigen_1tetra(4)-eps)*inv_epsilon41
   tweight(ind_dum(2),ieps) = tweight(ind_dum(2),ieps) + &
&   volconst_mult - cc*(eigen_1tetra(4)-eps)*inv_epsilon42
   tweight(ind_dum(3),ieps) = tweight(ind_dum(3),ieps) + &
&   volconst_mult - cc*(eigen_1tetra(4)-eps)*inv_epsilon43
   tweight(ind_dum(4),ieps) = tweight(ind_dum(4),ieps) + &
&   volconst_mult - cc*(four-(eigen_1tetra(4)-eps)*(inv_epsilon41+inv_epsilon42+inv_epsilon43))

   dccde = -three*volconst_mult*(eigen_1tetra(4)-eps)*(eigen_1tetra(4)-eps)&
&   *inv_epsilon41*inv_epsilon42*inv_epsilon43
   dtweightde(ind_dum(1),ieps) = dtweightde(ind_dum(1),ieps) &
&   -dccde*(eigen_1tetra(4)-eps)*inv_epsilon41 + cc*inv_epsilon41
   dtweightde(ind_dum(2),ieps) = dtweightde(ind_dum(2),ieps) &
&   -dccde*(eigen_1tetra(4)-eps)*inv_epsilon42 + cc*inv_epsilon42
   dtweightde(ind_dum(3),ieps) = dtweightde(ind_dum(3),ieps) &
&   -dccde*(eigen_1tetra(4)-eps)*inv_epsilon43 + cc*inv_epsilon43
   dtweightde(ind_dum(4),ieps) = dtweightde(ind_dum(4),ieps) &
&   -dccde*(four-(eigen_1tetra(4)-eps)*(inv_epsilon41+inv_epsilon42+inv_epsilon43)) &
&   -cc*(inv_epsilon41+inv_epsilon42+inv_epsilon43)
   eps = eps + deltaene
  end do
! 
! DEBUG
! write (*,*) 'eps, enemin+nn4*deltaene ', eps, enemin+nn4*deltaene
! eps = enemin+deltaene*nn4
! ENDDEBUG
! 
! 
! interval e4 < eps < enemax
! 
  do ieps=nn4+1,nene
   tweight(ind_dum(1),ieps) = tweight(ind_dum(1),ieps) + volconst_mult
   tweight(ind_dum(2),ieps) = tweight(ind_dum(2),ieps) + volconst_mult
   tweight(ind_dum(3),ieps) = tweight(ind_dum(3),ieps) + volconst_mult
   tweight(ind_dum(4),ieps) = tweight(ind_dum(4),ieps) + volconst_mult
!  dtweightde unchanged by this tetrahedron
   eps = eps + deltaene
  end do

! 
! if we have a fully degenerate tetrahedron,
! 1) the tweight is a Heaviside (step) function, which is correct above, but
! 2) the dtweightde should contain a Dirac function: add a Gaussian here
! 
  if (epsilon41 < tol6) then

!  to ensure the gaussian will integrate properly:
!  WARNING: this smearing could be problematic if too large
!  and doesnt integrate well if its too small
   gau_width = 10.0_dp*deltaene
!  gau_width = 0.001
   gau_width2 = 1.0 / gau_width / gau_width
   gau_prefactor = volconst_mult / gau_width / sqrt(pi)
!  DEBUG
!  write (*,*) 'get_tetra_weight : gau_width = ', gau_width
!  ENDDEBUG
!  
!  average position since bracket for epsilon41 is relatively large
   cc = (eigen_1tetra(1)+eigen_1tetra(2)+eigen_1tetra(3)+eigen_1tetra(4))/four
!  DEBUG
!  write (*,*) ' apply gaussian for Dirac delta: position = ',cc
!  ENDDEBUG
   eps = enemin
   do ieps=1,nene
    tmp = eps - cc
    dtweightde(ind_dum(4),ieps) = dtweightde(ind_dum(4),ieps) + &
&    gau_prefactor*exp(-tmp*tmp*gau_width2)
    eps = eps + deltaene
   end do
  end if
! end degenerate tetrahedron if


 end do
!end do itetra

!DEBUG
!write(6,*)' get_tetra_weight : exit '
!ENDDEBUG

end subroutine get_tetra_weight
!!***
