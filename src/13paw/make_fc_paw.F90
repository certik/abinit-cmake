!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_fc_paw
!! NAME
!! make_fc_paw
!!
!! FUNCTION
!! Compute the Fermi-contact term due to the PAW cores
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (JJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  natom=number of atoms in cell.
!!  ntypat=number of atom types
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  typat(natom)=type (integer) for each atom
!!
!! OUTPUT
!!  fc(natom), the Fermi-contact interaction at each site due to PAW

!! NOTES
!!
!! PARENTS
!!      calc_fc
!!
!! CHILDREN
!!      ratint
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine make_fc_paw(fc,natom,ntypat,pawrhoij,pawrad,pawtab,psps,typat)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(out) :: fc(natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iatom,itypat,ispden,jl,jlm,jlmn,jln,jm,j0lmn
 integer :: ilmn,il,km,iln,ilm,im,klmn,kln
 integer :: irhoij,nratint
 integer :: ig,imesh_size,mesh_size,fcgstep
 real(dp) :: gg,intg,yraterr,yratx
 real(dp) :: fcgh,fcgint,fcgerr,fcgint_old,fhatg
!arrays
 real(dp),allocatable :: ff(:),xratint(:),yratint(:)

! ************************************************************************

!DEBUG
!write(*,*)' make_fc_paw : enter'
!ENDDEBUG

 fc(:) = zero
 nratint = 10 ! number of points to use in Bulirsch-Stoer extrapolation
 allocate(xratint(nratint),yratint(nratint))

!loop over atoms in cell
 do iatom = 1, natom
  itypat = typat(iatom)
  mesh_size=pawrad(itypat)%mesh_size
  allocate(ff(mesh_size))

! loop over spin components
  do ispden=1,pawrhoij(iatom)%nspden

!  loop over basis elements for this atom
!  ----
   do jlmn=1,pawtab(itypat)%lmn_size
    jl= psps%indlmn(1,jlmn,itypat)  
    jm=psps%indlmn(2,jlmn,itypat)
    jlm = psps%indlmn(4,jlmn,itypat)
    jln=psps%indlmn(5,jlmn,itypat)
    j0lmn=jlmn*(jlmn-1)/2
    do ilmn=1,jlmn
     il= psps%indlmn(1,ilmn,itypat)
     im=psps%indlmn(2,ilmn,itypat)
     iln=psps%indlmn(5,ilmn,itypat)
     ilm = psps%indlmn(4,ilmn,itypat)
     klmn=j0lmn+ilmn
     kln = pawtab(itypat)%indklmn(2,klmn) ! need this for mesh selection below

     if(jl==il .and. jm==im) then ! this is the selection from the delta function on l=l' and m=m'
!     if (jl==0 .and. il==0) then ! select only s-states

!     Loop over non-zero elements of rhoij
      do irhoij=1,pawrhoij(iatom)%nrhoijsel
       if (klmn==pawrhoij(iatom)%rhoijselect(irhoij)) then ! rho_ij /= 0 for this klmn
        xratint(:)=pawrad(itypat)%rad(2:1+nratint)
        yratint(:)=(pawtab(itypat)%phiphj(2:1+nratint,kln)-pawtab(itypat)%tphitphj(2:1+nratint,kln))/&
&        (four_pi*pawrad(itypat)%rad(2:1+nratint)**2)
        call ratint(nratint,xratint,tol8,yratint,yraterr,yratx) ! this gives the value at r = 0 using
!       rational function extrapolation
        fc(iatom)=fc(iatom)+pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(irhoij,ispden)*yratx

!       gg = 0.000
!       fcgh = 0.001
!       fcgint = 0.000
!       fcgerr = 1.0
!       fcgstep = 0
!       do
!       gg = gg + fcgh
!       fcgstep = fcgstep + 1
!       do imesh_size=2,mesh_size
!       ff(imesh_size)=(pawtab(itypat)%phiphj(imesh_size,kln)-pawtab(itypat)%tphitphj(imesh_size,kln))&
!       &                        *sin(two_pi*gg*pawrad(itypat)%rad(imesh_size))/(two_pi*gg*pawrad(itypat)%rad(imesh_size)) 
!       end do
!       call deducer0(ff,mesh_size,pawrad(itypat))
!       call simp_gen(fhatg,ff,pawrad(itypat))
!       fcgint_old = fcgint
!       fcgint = fcgint + fcgh*four_pi*gg*gg*fhatg
!       fcgerr = abs(fcgint - fcgint_old)/abs(fcgint)
!       if (fcgstep > 1000000) then
!       exit
!       else if (fcgstep > 2000 .and. fcgerr < tol8) then
!       exit
!       end if
!       end do 
!       write(6,'(a,i8,a,f12.8,a,f12.8)')' fcgstep = ',fcgstep,' fcgint = ',fcgint,' fcgerr = ',fcgerr
!       fc(iatom)=fc(iatom)+pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(irhoij,ispden)*fcgint
       end if ! end selection on nonzero rho_ij
      end do ! end loop over nonzero rho_ij
     end if ! end selection on l=l' and m=m'
    end do ! end loop over ilmn
   end do ! end loop over jlmn
  end do ! end loop over spin densities
  deallocate(ff)
 end do     ! Loop on atoms

 deallocate(xratint,yratint)

!DEBUG
!write(6,*)' make_fc_paw : exit '
!stop
!ENDDEBUG

 end subroutine make_fc_paw
!!***
