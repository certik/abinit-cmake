!{\src2tex{textfont=tt}}
!!****f* ABINIT/dens_in_sph
!! NAME
!! dens_in_sph
!!
!! FUNCTION
!!   Calculate integrated density in sphere around each atom
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cg      = wavefunction coefficitents in recip space
!!  gmet    = metric in recip space
!!  istwfk  = storage mode for cg coefficients
!!  kg_k    = G vector indices
!!  natom   = number of atoms
!!  mpi_enreg=informations about MPI parallelization
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  npw_k   = number of plane waves for this kpoint
!!  ph3d    = phase factors for different atoms for all G vectors
!!  rmax(natom) = max radius to integrate to (in bohr)
!!
!! OUTPUT
!!  cmax    = integrated density for each atom for a rmax-radius sphere
!!
!! WARNING
!!  cg should not be modified by fourwf.
!!
!! PARENTS
!!      partial_dos_fractions,wffile
!!
!! CHILDREN
!!      dotprod_v,fftpac,fourdp,fourwf,ph1d3d,sphereboundary,sphericaldens
!!      sqnorm_g
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dens_in_sph(cmax,cg,gmet,istwfk,kg_k,natom,ngfft,mpi_enreg,npw_k,&
&                       paral_kgb,ph1d,rmax,ucvol)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12ffts
 use interfaces_12spacepar
 use interfaces_13nonlocal
 use interfaces_14occeig, except_this_one => dens_in_sph
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwfk,natom,npw_k,paral_kgb
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: kg_k(3,npw_k),ngfft(18)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: ph1d(2,(2*ngfft(1)+1+2*ngfft(2)+1+2*ngfft(3)+1)*natom)
 real(dp),intent(in) :: rmax(natom)
 real(dp),intent(inout) :: cg(2,npw_k)
 real(dp),intent(out) :: cmax(natom)

!Local variables -------------------------
!variables for sphereboundary
!variables for fourwf
!scalars
 integer :: cplex,i1,i2,i3,iatom,id1,id2,id3,ifft,ipw,mgfft,n1,n2,n3,n4,n5,n6
 integer :: nfft,nfftot,tim_fourwf=0
 real(dp) :: cmaxr,g1,g2,g3,norm,tmpr,weight
!arrays
 integer :: ngfft_here(18)
 integer,allocatable :: garr(:,:),gbound(:,:)
 real(dp) :: kpt(3)
 real(dp),allocatable :: denpot(:,:,:),fofgout(:,:),fofr(:,:,:,:),gnorm(:)
 real(dp),allocatable :: ph3d(:,:,:),phkxred(:,:),rhog(:,:),rhor(:)
 real(dp),allocatable :: sphrhog(:,:),vect(:,:)

! *********************************************************************

 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)
 n4=ngfft(4)
 n5=ngfft(5)
 n6=ngfft(6)
 nfftot = n1*n2*n3
 nfft=n1*n2*n3 
 ngfft_here(:) = ngfft(:)
!fourwf doesnt work with other options for mode 0 (fft G -> r)
 ngfft_here(7)=111
 ngfft_here(8)=256
 mgfft=maxval(ngfft_here(1:3))

 call sqnorm_g(norm,istwfk,mpi_enreg,npw_k,cg)

 if (abs(one-norm) > tol6) then
  write (6, '(a,f8.5)' ) ' dens_in_sph : this state is not normalized : norm=',norm
 end if

!-----------------------------------------------------------------
!inverse FFT of wavefunction to real space => density in real space
!-----------------------------------------------------------------
 allocate(gbound(2*mgfft+8,2))
 call sphereboundary(gbound,istwfk,kg_k,mgfft,npw_k)

 weight = one
 cplex=1
 allocate(denpot(cplex*n4,n5,n6))
 denpot(:,:,:)=zero
 allocate(fofgout(2,npw_k))
 allocate(fofr(2,n4,n5,n6))
 call fourwf(cplex,denpot,cg,fofgout,fofr,gbound,gbound, &
& istwfk,kg_k,kg_k,mgfft,mpi_enreg,1,ngfft_here,npw_k,&
& npw_k,n4,n5,n6,1,paral_kgb,tim_fourwf,weight)
 deallocate (fofgout,fofr)

 norm = sum(denpot(:,:,:))/nfftot
 if (abs(one-norm) > tol6) then
  write (6,'(a,f8.5)') ' dens_in_sph : this state is not normalized in real space : norm=',norm
 end if

!-----------------------------------------------------------------
!FFT of new density: we obtain n(G) in rhog(1,:)
!-----------------------------------------------------------------

!Change the packing of the reciprocal space density
 allocate (rhor(nfft))
 call fftpac(1,1,n1,n2,n3,n4,n5,n6,ngfft,rhor,denpot,1)

 allocate (rhog(2,nfft))
 call fourdp(1,rhog,rhor,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)

 deallocate (rhor,denpot,gbound)

 do ifft=1,nfft
  rhog(:,ifft) = rhog(:,ifft) / ucvol
 end do

!-----------------------------------------------------------------
!calculate norms of G vectors
!-----------------------------------------------------------------

 allocate(garr(3,nfft),gnorm(nfft))
 id3=ngfft(3)/2+2 ; id2=ngfft(2)/2+2 ; id1=ngfft(1)/2+2
 do i3=1,n3
  g3=i3-(i3/id3)*ngfft(3)-1
  do i2=1,n2
   g2=i2-(i2/id2)*ngfft(2)-1
   do i1=1,n1
    g1=i1-(i1/id1)*ngfft(1)-1
    ifft=i1+(i2-1)*n1+(i3-1)*n1*n2
    garr(1,ifft)=g1
    garr(2,ifft)=g2
    garr(3,ifft)=g3
    gnorm(ifft)=sqrt(gmet(1,1)*g1*g1 + &
&    two*gmet(2,1)*g2*g1 + &
&    two*gmet(3,1)*g3*g1 + &
&    gmet(2,2)*g2*g2 + &
&    gmet(3,2)*g3*g2 + &
&    gmet(3,3)*g3*g3)
   end do
  end do
 end do

!-----------------------------------------------------------------
!Compute the phases for the whole set of fft vectors
!-----------------------------------------------------------------

 allocate (ph3d(2,nfft,natom),phkxred(2,natom))
 kpt(:)=zero
 phkxred(1,:)=one
 phkxred(2,:)=zero
 call ph1d3d(1,natom,garr,kpt,natom,natom,nfft,ngfft(1),ngfft(2),ngfft(3),&
& phkxred,ph1d,ph3d)
 deallocate(phkxred)

!-----------------------------------------------------------------
!For each atom:
!call sphericaldens to calculate
!n(G) * 1/|G|^3  *  int_0^2*\pi*r_{max}*|G| 4 \pi y^2 j_0 (y) dy
!for all G vectors put into array sphrhog
!scalar product of phase factors with spherically convoluted density
!-----------------------------------------------------------------

 allocate(sphrhog(2,nfft))
 allocate(vect(2,nfft))

 do iatom=1,natom

  call sphericaldens(rhog,gnorm,nfft,rmax(iatom),sphrhog)

! Copy the phase factors, and take the compex conjugate, before evaluating
! the scalar product
  do ifft=1,nfft
   vect(1,ifft)= ph3d(1,ifft,iatom)
   vect(2,ifft)=-ph3d(2,ifft,iatom)
  end do
  cplex=2
  call dotprod_v(cplex,cmaxr,mpi_enreg,nfft,1,vect,sphrhog)
  cmax(iatom) = cmaxr

! DEBUG
! write (6, '(a,i4,a,es14.6,a,es12.6)' ) &
! &   ' dens_in_sph : At ', iatom, ' has ',cmaxr, &
! &      ' el.s in a sphere of rad ', rmax
! ENDDEBUG

 end do

 deallocate (rhog,gnorm)
 deallocate (garr,sphrhog,ph3d,vect)

end subroutine dens_in_sph
!!***
