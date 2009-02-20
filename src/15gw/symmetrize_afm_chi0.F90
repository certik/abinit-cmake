!{\src2tex{textfont=tt}}
!!****f* ABINIT/symmetrize_afm_chi0
!! NAME
!! symmetrize_afm_chi0
!!
!! FUNCTION
!!  Reconstruct the (down, down) component of the irreducible polarizability
!!  starting from the (up,up) element in case of systems with AFM symmetries 
!!  (i.e nspden==2 and nsppol=1). Return the trace (up,up)+(down,down) of the 
!!  matrix as required by GW calculations.
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Cryst<Crystal_structure>= Information on symmetries and unit cell.
!!  Gsph<Gvectors_type>= The G-sphere used to descrive chi0.
!!  npwe=Number of G-vectors in chi0.
!!  nomega=number of frequencies.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! chi0(npwe,npwe,nomega)= In input the up-up component, in output the trace of chi0. 
!! The value of matrix elements that should be zero due to AFM symmetry properties are 
!! forced to be zero (see NOTES below).
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  For each set of paired FM-AFM symmetries, the down-down component of
!!  a generic response function in reciprocal space can be obtained according to:
!!
!!  chi^{down,down}_{G1,G2}(q) = chi^{up up}_{G1,G2}(q) e^{iS(G1-G2).(tnonsFM - tnonsAFM)}
!!
!! where S is the rotational part common to the FM-AFM pair, tnonsFM and tnonsAFM 
!! are the fractional translations associated to the ferromagnetic and antiferromagnetic
!! symmetry, respectively.
!! Note that, if for a given G1-G2 pair, the phase e^{iS(G1-G2).(tnonsFM - tnonsAFM) depends 
!! on the FM-AFM symmetry pair, then the corresponding matrix element of chi0 must be zero.
!! Actually this is manually enforced in the code because this property might not be 
!! perfectly satisfied due to round-off errors.
!!
!! TODO 
!!  It is possible to symmetrize chi0 without any the extra allocation for AFM_mat.
!!  More CPU demanding but safer in case of a large chi0 matrix.
!!
!! PARENTS
!!      cchi0,cchi0q0
!!
!! CHILDREN
!!      assert
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine symmetrize_afm_chi0(Cryst,Gsph,npwe,nomega,chi0)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : czero_gw, cone_gw
 use m_errors, only : assert

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwe,nomega
 type(Gvectors_type),intent(in) :: Gsph
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 complex(gwpc),intent(inout) :: chi0(npwe,npwe,nomega)

!Local variables ------------------------------
!scalars
 integer :: io,ig1,ig2,isymf,isyma,nsym,ipair,k0g,kg,npairs,nonzero
 complex(gwpc) :: phase,phase_old
 logical :: found
!arrays
 integer :: rotfm(3,3),rotafm(3,3),pairs2sym(2,Cryst%nsym/2)
 real(dp) :: tfm(3),tafm(3)
 complex(gwpc),allocatable :: AFM_mat(:)

!************************************************************************

 nsym = Cryst%nsym
 npairs = 0

 do isymf=1,nsym
  if (Cryst%symafm(isymf)==-1) CYCLE
  rotfm = Cryst%symrec(:,:,isymf)
  tfm = Cryst%tnons(:,isymf)
  found = .FALSE.

  do isyma=1,nsym
   if (Cryst%symafm(isyma)==1) CYCLE
   rotafm = Cryst%symrec(:,:,isyma)

   if (ALL(rotfm==rotafm)) then
    found=.TRUE.
    tafm = Cryst%tnons(:,isyma)
    npairs=npairs+1
    call assert((npairs<=nsym/2),"Wrong AFM group",__FILE__,__LINE__)
    pairs2sym(1,npairs)=isymf
    pairs2sym(2,npairs)=isyma
   end if
  end do !isyma

  call assert(found,"Corresponding AFM rotation not found",__FILE__,__LINE__)
 end do !isymf

 call assert((npairs==nsym/2),"Wrong AFM space group",__FILE__,__LINE__)

 allocate(AFM_mat(npwe*(npwe+1)/2))

 do ig2=1,npwe
  k0g=ig2*(ig2-1)/2 
  do ig1=1,ig2
   kg=k0g+ig1  
   nonzero=1

   do ipair=1,nsym/2
    isymf = pairs2sym(1,ipair)
    isyma = pairs2sym(2,ipair)
    phase = ( Gsph%phmSGt(ig1,isymf)*CONJG(Gsph%phmSGt(ig1,isyma)) ) * &
            ( Gsph%phmSGt(ig2,isymf)*CONJG(Gsph%phmSGt(ig2,isyma)) )

    if (ipair>1 .and. (ABS(phase_old-phase) > tol6)) then 
     nonzero = 0
     EXIT
    end if
    phase_old = phase
   end do !ipair

   AFM_mat(kg)=nonzero*(cone_gw + phase)
  end do !ig1
 end do !ig2


 do io=1,nomega
  
  ! Take care of diagonal.
  do ig1=1,npwe
   chi0(ig1,ig1,io)=2.0_gwp*chi0(ig1,ig1,io)
  end do

  ! Upper and lower triangle are treated differently:
  ! We take advantage of the fact the AFM_mat is hermitian to reduce memory.
  do ig2=2,npwe
   k0g=ig2*(ig2-1)/2 
   do ig1=1,ig2-1
    kg=k0g+ig1  
    chi0(ig1,ig2,io)=AFM_mat(kg)*chi0(ig1,ig2,io)
   end do
  end do

  do ig1=2,npwe
   k0g=ig1*(ig1-1)/2 
   do ig2=1,ig1-1
    kg=k0g+ig2 
    chi0(ig1,ig2,io)=CONJG(AFM_mat(kg))*chi0(ig1,ig2,io)
   end do
  end do

 end do !io 

 deallocate(AFM_mat)

end subroutine symmetrize_afm_chi0
!!***
