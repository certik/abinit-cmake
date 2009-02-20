!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_fc_el
!! NAME
!! make_fc_el
!!
!! FUNCTION
!! compute the Fermi-contact term due to electron density
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (JJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! gcart(ngfft(1),ngfft(2),ngfft(3),3), G vectors in cartesian space
!! natom, number of atoms in unit cell
!! nfft,ngfft(18), number of FFT points and details of FFT
!! nhat(nfft,nspden), augmentation electron density $\hat{n}$
!! nspden, number of spin components
!! rhor(nfft,nspden), valence electron density, here $\tilde{n} + \hat{n}$
!! rprimd(3,3), conversion from crystal coordinates to cartesian coordinates
!!
!! OUTPUT
!! fc(natom), the Fermi-contact term at each atomic site due to rhor-nhat
!!
!! SIDE EFFECTS
!! xred(3,natom), location of atoms in crystal coordinates. It appears with intent(inout) here
!!                because routine xredxcart requires the ability to change xred, although we do not
!!                use that facility in the present call.
!!
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine make_fc_el(fc,gcart,natom,nfft,ngfft,nhat,nspden,paral_kgb,rhor,rprimd,xred)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12ffts
 use interfaces_12geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,paral_kgb
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gcart(ngfft(1),ngfft(2),ngfft(3),3),nhat(nfft,nspden)
 real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3)
 real(dp),intent(inout) :: xred(3,natom)
 real(dp),intent(out) :: fc(natom)

!Local variables-------------------------------
!scalars
 integer :: cplex,iatom,igfft1,igfft2,igfft3,ii,index,isign,ispden,jj
 integer :: tim_fourdp
 real(dp) :: cph,derivs,phase,sph,trace
 type(MPI_type) :: mpi_enreg
!arrays
 real(dp) :: gvec(3),ratom(3)
 real(dp),allocatable :: fc_c(:,:),fofg(:,:),fofr(:),xcart(:,:)
!no_abirules
!

! ************************************************************************

!DEBUG
!write(*,*)' make_efg_el : enter'
!ENDDEBUG

 allocate(fofg(2,nfft),fofr(nfft),fc_c(2,natom),xcart(3,natom))

 fc(:) = zero
 fc_c(:,:) = zero
 call xredxcart(natom,1,rprimd,xcart,xred) ! get atomic locations in cartesian coords

 tim_fourdp = 0 ! timing code, not using
 isign = -1 ! FT from R to G
 cplex = 1 ! fofr is real
!here we are only interested in the total electron particle density, which is rhor(:,1)-nhat(:,1)
!regardless of the value of nspden. This may change in the future depending on 
!developments with noncollinear magnetization and so forth. Such a change will
!require an additional loop over nspden.
!now construct electron particle density. Note that $\mathrm{rhor} = \tilde{n} + \hat{n}$, therefore we
!must subtract $\hat{n}$.
 fofr(:) = rhor(:,1) - nhat(:,1)
 call fourdp(cplex,fofg,fofr,isign,mpi_enreg,nfft,ngfft,paral_kgb,tim_fourdp) ! construct density in G space

 do igfft1 = 1, ngfft(1) ! summing over G space vector components...
  do igfft2 = 1, ngfft(2)
   do igfft3 = 1, ngfft(3)
    index = (igfft3-1)*ngfft(2)*ngfft(1) + (igfft2-1)*ngfft(1) + igfft1
    gvec(:) = gcart(igfft1,igfft2,igfft3,:) ! gvec is the current vector in G space
    do iatom = 1, natom ! sum over atoms in unit cell
     ratom(:) = xcart(:,iatom) ! extract location of atom iatom
     phase = two_pi*dot_product(gvec,ratom) ! argument of $e^{2\pi i G\cdot R}$
     cph = cos(phase)
     sph = sin(phase)
     fc_c(1,iatom) = fc_c(1,iatom) + fofg(1,index)*cph ! real part of fc
     fc_c(2,iatom) = fc_c(2,iatom) + fofg(2,index)*sph ! imaginary part, should be sum to zero if G grid dense enough
    end do ! end loop over atoms in cell
   end do ! end loop over ngfft(3)
  end do ! end loop over ngfft(2)
 end do ! end loop over ngfft(1)

 fc(:) = fc_c(1,:) ! extract real part only of final tensor

 deallocate(fofg,fofr,fc_c,xcart)

!DEBUG
!write(6,*)' make_fc_el : exit '
!stop
!ENDDEBUG

end subroutine make_fc_el
!!***
