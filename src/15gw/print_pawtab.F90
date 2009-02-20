!{\src2tex{textfont=tt}}
!!****f* ABINIT/print_pawtab
!! NAME
!! print_pawtab
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Pawtab<pawtab_type> Only for PAW, TABulated data initialized at start
!!
!! OUTPUT
!!  Only writing  
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine print_pawtab(Pawtab,unitno,prtvol,mode_paral)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_13paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unitno
 character(len=4),intent(in),optional :: mode_paral
!arrays
 type(pawtab_type) :: Pawtab(:)

!Local variables-------------------------------
!scalars
 integer :: ierr,ityp,ntypat,unt,verb
 character(len=4) :: mode
 character(len=500) :: msg

! *************************************************************************

 verb=0      ; if (PRESENT(prtvol))     verb=prtvol                         
 unt=std_out ; if (PRESENT(unitno))     unt=unitno     
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 ntypat=SIZE(Pawtab(:))

 write(*,*) ' ==================================== '
 write(*,*) ' ==== Info on PAW TABulated data ==== '
 write(*,*) ' ==================================== '
 do ityp=1,ntypat 
  ! Print out integer values (dimensions)
  write(*,*)'                                 '
  write(*,*)'  ****************************** ' 
  write(*,*)'  **** Atom type ',ityp,' **** ' 
  write(*,*)'  ****************************** ' 
  write(*,*)'  Number of (n,l) elements ...................... ',Pawtab(ityp)%basis_size
  write(*,*)'  Number of (l,m,n) elements .................... ',Pawtab(ityp)%lmn_size
  write(*,*)'  Number of (i,j) elements (packed form) ........ ',Pawtab(ityp)%ij_size  
  write(*,*)'  Max L+1 leading to non-zero Gaunt ............. ',Pawtab(ityp)%l_size
  write(*,*)'  Max L+1 leading to non-zero Gaunt (pawlcutd) .. ',Pawtab(ityp)%lcut_size
  write(*,*)'  lmn2_size ..................................... ',Pawtab(ityp)%lmn2_size
  write(*,*)'  lmnmix_sz ..................................... ',Pawtab(ityp)%lmnmix_sz 
  write(*,*)'  Size of radial mesh ........................... ',Pawtab(ityp)%mesh_size 
  write(*,*)'  No of Q-points for tcorespl and tvalespl ...... ',Pawtab(ityp)%mqgrid 
  write(*,*)'  Radial shape function type .................... ',Pawtab(ityp)%shape_type
  write(*,*)'  shape_lambda .................................. ',Pawtab(ityp)%shape_lambda 
  write(*,*)'  Use pseudized core density .................... ',Pawtab(ityp)%usetcore
  write(*,*)'  Use pseudized valence density ................. ',Pawtab(ityp)%usetvale 
  write(*,*)'  Option for Vloc (1 Blochl, 2 Kresse) .......... ',Pawtab(ityp)%vlocopt  

  write(*,*)'  useexexch ',Pawtab(ityp)%useexexch
  write(*,*)'  usepawu ',Pawtab(ityp)%usepawu
  write(*,*)'  L on which local exact-exchange is applied .... ',Pawtab(ityp)%lexexch 
  write(*,*)'  L on which U is applied ....................... ',Pawtab(ityp)%lpawu 
  write(*,*)'  Number of (i,j) elements for PAW+U and EXX .... ',Pawtab(ityp)%ij_proj
  write(*,*)'  Number of projectors on which U or EXX acts ... ',Pawtab(ityp)%nproju 
  write(*,*)'  Has nabla = ',pawtab(ityp)%has_nabla
  !
  ! Real (real(dp)) scalars
  write(*,*)'  1/q d(tNcore(q))/dq for q=0 ',Pawtab(ityp)%dncdq0
  write(*,*)'  1/q d(tNvale(q))/dq for q=0 ',Pawtab(ityp)%dnvdq0
  write(*,*)'  XC energy for the core density ',Pawtab(ityp)%exccore
  write(*,*)'  Value of the J parameter [eV] ....................... ',Pawtab(ityp)%jpawu*Ha_eV
  write(*,*)'  Radius of the PAW sphere ............................ ',Pawtab(ityp)%rpaw
  write(*,*)'  Compensation charge radius (if > rshp, g(r)=0) ...... ',Pawtab(ityp)%rshp !(if r>rshp, g(r)=zero)
  write(*,*)'  Sigma parameter in gaussian shape function .......... ',Pawtab(ityp)%shape_sigma !(shape_type=2)
  write(*,*)'  Value of the U parameter [eV] ....................... ',Pawtab(ityp)%upawu*Ha_eV
 end do
 !
 ! Other huge arrays are not reported
 
end subroutine print_pawtab
!!***
