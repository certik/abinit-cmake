!{\src2tex{textfont=tt}}
!!****f* ABINIT/mklocl_realspace
!!
!! NAME
!! mklocl_realspace
!!
!! FUNCTION
!! This method is equivalent to mklocl_recipspace except that
!! it uses real space pseudo-potentials. It is usefull for isolated
!! systems. Then the option 3 and 4 are not available for this
!! implementation.
!!
!! Optionally compute :
!!  option=1 : local ionic potential throughout unit cell
!!  option=2 : contribution of local ionic potential to E gradient wrt xred
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in unit cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms.
!!  option= (see above)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=electron density rho(G) (electrons/$\textrm{Bohr}^3$)
!!  rhor(nfft,nspden)=electron density in electrons/bohr**3.
!!    (needed if option==2 or if option==4)
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol=unit cell volume ($\textrm{Bohr}^3$).
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (if option==1) vpsp(nfft)=local crystal pseudopotential in real space.
!!  (if option==2) grtn(3,natom)=grads of Etot wrt tn. These gradients are in
!!                 reduced coordinates. Multiply them by rprimd to get
!!                 gradients in cartesian coordinates.
!!
!! SIDE EFFECTS
!!
!!
!! PARENTS
!!      mklocl
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mklocl_realspace(dtset, grtn, mgfft, mpi_enreg, natom, nattyp, nfft, ngfft, &
                          & nspden, ntypat, option, ph1d, psps, rhog, rhor, rmet, &
                          & rprimd, ucvol, vpsp, xred)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_14poisson
 use interfaces_15common, except_this_one => mklocl_realspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,natom,nfft,nspden,ntypat,option
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),rhog(2,nfft)
 real(dp),intent(in) :: rhor(nfft,nspden),rmet(3,3),rprimd(3,3)
 real(dp),intent(inout) :: xred(3,natom)
 real(dp),intent(out) :: grtn(3,natom),vpsp(nfft)

!Local variables-------------------------------
 character(len=1) :: geocode
  !testing variables
 !scalars
 integer,parameter :: nStep=2
 integer :: countParSeconde,i1,i2,i3,ia,ia1,ia2,igeo,ii,ind,itypat,ix,iy,iz,jj
 integer :: kk,me_fft,n1,n2,n3,n_interpol,nproc_fft,tpsStart,tpsStop
 real(dp),parameter :: min_rho_value=1.0d-12
 real(dp) :: aa,bb,cc,dd,delta,deltaV,dr,dr2div6,invdr,r,vol_interpol,x,y,z,hgx,hgy,hgz,entmp
 logical,parameter :: customRho=.false.,finiteDiff=.false.,testing=.false.
 logical :: doIt
 character(len=500) :: message
!arrays
 integer :: ngfft_interpol(18)
 real(dp) :: coord(3),coordXYZ(3),refValue(3),tsec(2)
 real(dp),allocatable :: coordCart_interpol(:,:),coordRed_interpol(:,:)
 real(dp),allocatable :: grad_sum(:,:),gridcart(:,:),gridred(:,:)
 real(dp),allocatable :: grtn_cart_interpol(:,:),grtn_diff(:,:)
 real(dp),allocatable :: rhog_interpol(:,:),rhog_testing(:,:),rhor_interpol(:)
 real(dp),allocatable :: rhor_testing(:),rhor_work(:),xcart(:,:),vhartr(:),gxyz(:,:)
 real(dp),pointer :: kernel(:)

! *************************************************************************

 if (dtset%icoulomb == 1) then
  geocode='F'
 else if (dtset%icoulomb == 2) then
  geocode='S'
 end if

!*************************************************************************

!Keep track of total time spent in mklocl
 if(option==2)then
  call timab(72,1,tsec)
 end if

 if (testing) then
  call system_clock(count_rate = countParSeconde)
  call system_clock(tpsStart, count_rate = countParSeconde)
 end if

 n1 = ngfft(1)
 n2 = ngfft(2)
 n3 = ngfft(3)
 me_fft = ngfft(11)
 nproc_fft = ngfft(10)

!!$ !this in principle is not anymore true
!!$ if (nspden /= 1) then
!!$
!!$  write(message, '(a,a,a,a)' ) ch10,&
!!$&  ' mklocl_realspace : ERROR - ',ch10,&
!!$&  '  real space computation is only possible without spin density'
!!$  call wrtout(06,message,'COLL')
!!$  call leave_new('COLL')
!!$ end if

!Store xcart for each atom
 allocate(xcart(3, natom))
 call xredxcart(natom, 1, rprimd, xcart, xred)
!Store cartesian coordinates for each grid points

 allocate(gridred(3, nfft))
 allocate(gridcart(3, nfft))
 ii = 0
 do i3 = 1, n3, 1
  coord(3) = real(i3 - 1, dp) / real(n3, dp)
  do i2 = 1, n2, 1
   coord(2) = real(i2 - 1, dp) / real(n2, dp)
   do i1 = 1, n1, 1
    ii = ii + 1
    coord(1) = real(i1 - 1, dp) / real(n1, dp)
    gridred(:, ii) = coord(:)
   end do
  end do
 end do
 call xredxcart(nfft, 1, rprimd, gridcart, gridred)
 deallocate(gridred)

!the branch with the HGH treatment of the PSP will presumably start here
!here we need to put the if statement for the PSP code =2,3,10 for GTH-HGH

!see whether all the PSP considered are of type GTH-HGH
 doIt=.true.
!doIt=.false.
 do ii=1,psps%npsp
  doIt=doIt .and.&
  (psps%pspcod(ii)==2 .or.psps%pspcod(ii)==3 .or. psps%pspcod(ii)==10)
 end do

 if (doIt) then

! WARNING: FOR THE MOMENT THIS ROUTINE DOES NOT WORK IN PARALLEL, DEPENDs ON THE COMMUNICATOR

! definition of the grid spacings as in the kernel routine
  hgx = rprimd(1,1)/(ngfft(1))
  hgy = rprimd(2,2)/(ngfft(2))
  hgz = rprimd(3,3)/(ngfft(3))

  call PSolver_kernel(dtset, 2, kernel, mpi_enreg, rprimd)

  if (option==1) then

!  !$       !calculate the center of mass of the atomic system such as to put the molecule in the middle of the simulation box
!  !$       xcm=0.d0
!  !$       ycm=0.d0
!  !$       zcm=0.d0
!  !$       do ii=1,dtset%natom
!  !$          xcm=xcm+xcart(1,ii)
!  !$          ycm=ycm+xcart(2,ii)
!  !$          zcm=zcm+xcart(3,ii)
!  !$       end do
!  !$       xcm=xcm/real(dtset%natom,dp)
!  !$       ycm=ycm/real(dtset%natom,dp)
!  !$       zcm=zcm/real(dtset%natom,dp)
!  !$
!  !$       xshift=0.5d0*real(n1-1,dp)*hgx-xcm
!  !$       yshift=0.5d0*real(n2-1,dp)*hgy-ycm
!  !$       zshift=0.5d0*real(n3-1,dp)*hgz-zcm
!  !$
!  !$       !do ii=1,dtset%natom
!  !$       !   xcart(1,ii)=xcart(1,ii)+xshift
!  !$       !   xcart(2,ii)=xcart(2,ii)+yshift
!  !$       !   xcart(3,ii)=xcart(3,ii)+zshift
!  !$       !end do

   call createIonicPotential_new(geocode,mpi_enreg%me, mpi_enreg%nproc, dtset%natom, &
&   dtset%ntypat, dtset%typat, psps%gth_params%psppar, &
&   int(psps%ziontypat), xcart,gridcart, hgx,hgy,hgz, real(0, dp), &
&   n1,n2,n3, kernel, vpsp,mpi_enreg%world_comm)

  else if (option ==2) then
   
!  the local forces with this formalism are calculated differently
   
!  Compute Hartree's potential from rhor.
   allocate(vhartr(nfft))
   call PSolver_hartree(dtset, entmp, mpi_enreg, rhor, rprimd, vhartr)


   allocate(gxyz(3, dtset%natom))
!  calculate local part of the forces grtn (inspired from BigDFT routine)
   call local_forces_new(geocode,mpi_enreg%me, mpi_enreg%nproc, dtset%ntypat, dtset%natom, &
&   dtset%typat, xcart, gridcart, psps%gth_params%psppar, &
&   int(psps%ziontypat), hgx,hgy,hgz, n1,n2,n3,&
&   rhor,vhartr, gxyz)
   deallocate(vhartr)

!  Forces should be in reduced coordinates.
   do ia = 1, dtset%natom, 1
    do igeo = 1, 3, 1
     grtn(igeo, ia) = - rprimd(1, igeo) * gxyz(1, ia) - &
&     rprimd(2, igeo) * gxyz(2, ia) - &
&     rprimd(3, igeo) * gxyz(3, ia)
    end do
   end do

!  Deallocate local variables
   deallocate(gxyz)



!  print *,'not implemented for the moment, wait a minute...'
!  stop
  end if


  deallocate(xcart)
  deallocate(gridcart)

! else statement for the non GTH-HGH PSP 
 else

! dr is the r step in the sampling psps%vlspl
  dr = psps%qgrid_vl(2)
  invdr = 1._dp / dr
  dr2div6 = dr * dr / 6._dp

  if (option == 1) then
!  Set 0 in vpsp before summing
   vpsp(:) = 0._dp
  else if (option == 2) then
!  Allocate array to store cartesian gradient computed with
!  an interpolation of rhor
   allocate(grtn_cart_interpol(3, natom))
   grtn_cart_interpol(:, :) = 0._dp

   n_interpol = nStep ** 3
   allocate(coordRed_interpol(3, nStep ** 3))
   allocate(coordCart_interpol(3, nStep ** 3))

   if (testing .and. customRho) then
!   Use a custom rho instead of the self-consistent one.
    allocate(rhor_testing(nfft))
    allocate(rhog_testing(2, nfft))
   end if

   allocate(rhor_interpol(nfft * n_interpol))
   allocate(rhor_work(nfft * n_interpol))
   allocate(rhog_interpol(2, nfft * n_interpol))

   if (testing .and. customRho) then
!   Testing only, changing rho with a centered gaussian
    do ii = 1, nfft, 1
!    using the position of the first atom as center.
     r = (gridcart(1, ii) - xcart(1, 1)) ** 2 + &
&     (gridcart(2, ii) - xcart(2, 1)) ** 2 + &
&     (gridcart(3, ii) - xcart(3, 1)) ** 2
     rhor_testing(ii) = exp(-r/4._dp)
    end do
!   Testing only, compute rhog_testing from rhor_testing
    call fourdp(1, rhog_testing, rhor_testing, -1, mpi_enreg, nfft, ngfft, dtset%paral_kgb,0)
   end if

!  Compute the interpolation of rho, using a fourrier transform
   rhog_interpol(:, :) = 0._dp
   ii = 0
   do i3 = 1, n3, 1
    if (i3 <= n3 / 2) then
     iz = i3
    else
     iz = n3 * nStep - n3 + i3
    end if
    do i2 = 1, n2, 1
     if (i2 <= n2 / 2) then
      iy = i2
     else
      iy = n2 * nStep - n2 + i2
     end if
     do i1 = 1, n1, 1
      ii = ii + 1
      if (i1 <= n1 / 2) then
       ix = i1
      else
       ix = n1 * nStep - n1 + i1
      end if
      jj = (iz - 1) * n2 * n1 * nStep ** 2 + (iy - 1) * n3 * nStep + ix
      if (testing .and. customRho) then
       rhog_interpol(:, jj) = rhog_testing(:, ii)
      else
       rhog_interpol(:, jj) = rhog(:, ii)
      end if
     end do
    end do
   end do

!  Compute the interpolation of rho from the Fourier transformation
   ngfft_interpol(:) = ngfft(:)
   ngfft_interpol(1:3) = (/ n1 * nStep, n2 * nStep, n3 * nStep /)
   ngfft_interpol(4:6) = (/ n1 * nStep + 1, n2 * nStep + 1, n3 * nStep /)
   call fourdp(1, rhog_interpol, rhor_work, 1, mpi_enreg, nfft * n_interpol, ngfft_interpol, dtset%paral_kgb,0)

!  Reorder rhor_interpol to be able to read it linearly
   jj = 0
   do i3 = 1, n3, 1
    do i2 = 1, n2, 1
     do i1 = 1, n1, 1
      do iz = 1, nStep, 1
       do iy = 1, nStep, 1
        do ix = 1, nStep, 1
         jj = jj + 1
         kk = ((i3 - 1) * nStep + iz - 1) ! z coordinate in the interpolated grid
         kk = kk * n1 * n2 * nStep ** 2
         kk = kk + ((i2 - 1) * nStep + iy - 1) * n1 * nStep ! adding y coordinate
         kk = kk + (i1 - 1) * nStep + ix ! adding x coordinate
         rhor_interpol(jj) = rhor_work(kk)
        end do
       end do
      end do
     end do
    end do
   end do
   deallocate(rhor_work)

!  Compute grid access in the interpolated volume
   ii = 0
   do iz = 1, nStep, 1
    z = real(iz - 1, dp) / real(nStep, dp)
    do iy = 1, nStep, 1
     y = real(iy - 1, dp) / real(nStep, dp)
     do ix = 1, nStep, 1
      x = real(ix - 1, dp) / real(nStep, dp)
      ii = ii + 1
      coordRed_interpol(:, ii) = (/ x, y, z /)
!     Assuming orthogonal box (should be change later)
      coordCart_interpol(:, ii) = (/ x * rprimd(1, 1) / real(n1, dp), &
&      y * rprimd(2, 2) / real(n2, dp), &
&      z * rprimd(3, 3) / real(n3, dp) /)
     end do
    end do
   end do

   vol_interpol = 1._dp / real(nStep, dp) ** 3
!  Compute the coordinates (integer) of each atom and deduce
!  the max extens of the integral summation.
!  !$  do ia = 1, natom, 1
!  !$   coordAtom(1, ia) = int(xred(1, ia) * n1) + 1
!  !$   coordAtom(2, ia) = int(xred(2, ia) * n2) + 1
!  !$   coordAtom(3, ia) = int(xred(3, ia) * n3) + 1
!  !$  end do
  end if

  if (testing .and. option == 2) then
   call system_clock(tpsStop, count_rate = countParSeconde)
   write(*,*) "Tps : ", real(tpsStop - tpsStart) / real(countParSeconde)
  end if

  ia1=1
  do itypat = 1, ntypat, 1
!  ia1,ia2 sets range of loop over atoms:
   ia2 = ia1 + nattyp(itypat) - 1

   do ii = 1, nfft, 1
    do ia = ia1, ia2, 1
     if (option == 1) then
!     Compute the potential
!     r is the distance between grid point and atom
      r = sqrt((gridcart(1, ii) - xcart(1, ia)) ** 2 + &
&      (gridcart(2, ii) - xcart(2, ia)) ** 2 + &
&      (gridcart(3, ii) - xcart(3, ia)) ** 2)

!     Coefficients needed to compute the spline.
      jj = int(r * invdr) + 1
      if (jj > psps%mqgrid_vl - 2) then
       write(message, '(a,a,a,a,a,a,i0,a,i0,a,a)' ) ch10,&
&       ' mklocl_realspace : ERROR - ',ch10,&
&       '  pseudo-potential local part sampling is not wide enough', ch10, &
&       '  want to access position ', jj, ' whereas mqgrid_vl = ', psps%mqgrid_vl, ch10, &
&       '  Action : no idea, contact developpers...'
       call wrtout(06,message,'COLL')
       call leave_new('COLL')
      end if
      delta = r - psps%qgrid_vl(jj)
      bb = delta * invdr
      aa = 1._dp - bb
      cc = aa * (aa ** 2 - 1._dp) * dr2div6
      dd = bb * (bb ** 2 - 1._dp) * dr2div6

!     compute V(r) from the spline, jj and jj + 1 is braketting r in
!     the sampling
      deltaV = aa * psps%vlspl(jj, 1, itypat) + bb * psps%vlspl(jj + 1, 1, itypat) + &
&      cc * psps%vlspl(jj, 2, itypat) + dd * psps%vlspl(jj + 1, 2, itypat)
!     Add on grid point ii the contribution of atom ia
      vpsp(ii) = vpsp(ii) + deltaV
     else if (option == 2) then
!     Compute the forces, as gradient of energy (V(r).rho(r))

!     Testing only - reference points
      if (.false.) then
!      r is the distance between grid point and atom
       r = sqrt((gridcart(1, ii) - xcart(1, ia)) ** 2 + &
&       (gridcart(2, ii) - xcart(2, ia)) ** 2 + &
&       (gridcart(3, ii) - xcart(3, ia)) ** 2)

!      Coefficients needed to compute the spline.
       jj = int(r * invdr) + 1
       delta = r - psps%qgrid_vl(jj)
       bb = delta * invdr
       aa = 1._dp - bb
       cc = aa * (aa ** 2 - 1._dp) * dr2div6
       dd = bb * (bb ** 2 - 1._dp) * dr2div6

!      When mesh position is on a node, forces are null.
       if (r /= 0._dp) then
!       This value deltaV is the first derivative of V(r) taken at r.
        deltaV = aa * psps%dvlspl(jj, 1, itypat) + bb * psps%dvlspl(jj + 1, 1, itypat) + &
&        cc * psps%dvlspl(jj, 2, itypat) + dd * psps%dvlspl(jj + 1, 2, itypat)
!       We multiply by rho(r) to have an energy.
        deltaV = deltaV * rhor(ii, 1) / r
        refValue(:) = - deltaV * (gridcart(:, ii) - xcart(:, ia))
        grtn_cart_interpol(:, ia) = grtn_cart_interpol(:, ia) + refValue(:)
       end if
      end if

!     Compute the interpolation for the point ii
      ind = (ii - 1) * n_interpol
      do kk = 1, n_interpol, 1
       ind = ind + 1

       if (rhor_interpol(ind) > min_rho_value) then
!       Assume orthogonal box...
        coordXYZ(1) = gridcart(1, ii) - xcart(1, ia) + coordCart_interpol(1, kk)
        coordXYZ(2) = gridcart(2, ii) - xcart(2, ia) + coordCart_interpol(2, kk)
        coordXYZ(3) = gridcart(3, ii) - xcart(3, ia) + coordCart_interpol(3, kk)
        r = coordXYZ(1) ** 2 + coordXYZ(2) ** 2 + coordXYZ(3) ** 2

        if (r /= 0._dp) then
         r = sqrt(r)
!        Coefficients needed to compute the spline.
         jj = int(r * invdr) + 1
         delta = r - psps%qgrid_vl(jj)
         bb = delta * invdr
         aa = 1._dp - bb
         cc = aa * (aa ** 2 - 1._dp) * dr2div6
         dd = bb * (bb ** 2 - 1._dp) * dr2div6
         deltaV = aa * psps%dvlspl(jj, 1, itypat) + &
&         bb * psps%dvlspl(jj + 1, 1, itypat) + &
&         cc * psps%dvlspl(jj, 2, itypat) + &
&         dd * psps%dvlspl(jj + 1, 2, itypat)
         deltaV = deltaV * rhor_interpol(ind) / r
         grtn_cart_interpol(1, ia) = grtn_cart_interpol(1, ia) - deltaV * coordXYZ(1)
         grtn_cart_interpol(2, ia) = grtn_cart_interpol(2, ia) - deltaV * coordXYZ(2)
         grtn_cart_interpol(3, ia) = grtn_cart_interpol(3, ia) - deltaV * coordXYZ(3)
!        do igeo = 1, 3, 1
!        grtn_cart_interpol(igeo, ia) = grtn_cart_interpol(igeo, ia) - deltaV * coordXYZ(igeo)
!        end do
        end if
       end if
      end do

!     =============
!     Testing only
!     =============
!     use of finite differences
      if (finiteDiff) then
       do igeo = 1, 3, 1
        coord(:) = 0._dp
        coord(igeo) = dr / 2.0_dp
        r = sqrt((gridcart(1, ii) - xcart(1, ia) + coord(1)) ** 2 + &
&        (gridcart(2, ii) - xcart(2, ia) + coord(2)) ** 2 + &
&        (gridcart(3, ii) - xcart(3, ia) + coord(3)) ** 2)

!       Coefficients needed to compute the spline.
        jj = int(r * invdr) + 1
        delta = r - psps%qgrid_vl(jj)
        bb = delta * invdr
        aa = 1._dp - bb
        cc = aa * (aa ** 2 - 1._dp) * dr2div6
        dd = bb * (bb ** 2 - 1._dp) * dr2div6

        deltaV = aa * psps%vlspl(jj, 1, itypat) + bb * psps%vlspl(jj + 1, 1, itypat) + &
&        cc * psps%vlspl(jj, 2, itypat) + dd * psps%vlspl(jj + 1, 2, itypat)


        coord(:) = 0._dp
        coord(igeo) = -dr / 2.0_dp
        r = sqrt((gridcart(1, ii) - xcart(1, ia) + coord(1)) ** 2 + &
&        (gridcart(2, ii) - xcart(2, ia) + coord(2)) ** 2 + &
&        (gridcart(3, ii) - xcart(3, ia) + coord(3)) ** 2)

!       Coefficients needed to compute the spline.
        jj = int(r * invdr) + 1
        delta = r - psps%qgrid_vl(jj)
        bb = delta * invdr
        aa = 1._dp - bb
        cc = aa * (aa ** 2 - 1._dp) * dr2div6
        dd = bb * (bb ** 2 - 1._dp) * dr2div6

        deltaV = deltaV - (aa * psps%vlspl(jj, 1, itypat) + &
&        bb * psps%vlspl(jj + 1, 1, itypat) + &
&        cc * psps%vlspl(jj, 2, itypat) + &
&        dd * psps%vlspl(jj + 1, 2, itypat))
        grtn_diff(igeo, ia) = grtn_diff(igeo, ia) - deltaV * rhor(ii, 1) / dr
       end do
      end if
!     =============
!     Testing only
!     =============

     end if
    end do
!   End loop over atoms of type itypat
   end do
!  End loop over real space grid points

   ia1 = ia2 + 1
  end do
! End loop over type of atoms

  deallocate(xcart)
  deallocate(gridcart)

  if(option==2)then
!  multiply the forces by the volume of a single box mesh.
   grtn_cart_interpol(:, :) = grtn_cart_interpol(:, :) * &
&   ucvol / real(n1 * n2 * n3, dp) * vol_interpol
!  Transform cartesian forces to reduce coordinates
   do ia = 1, natom, 1
    do igeo = 1, 3, 1
     grtn(igeo, ia) = rprimd(1, igeo) * grtn_cart_interpol(1, ia) + &
&     rprimd(2, igeo) * grtn_cart_interpol(2, ia) + &
&     rprimd(3, igeo) * grtn_cart_interpol(3, ia)
    end do
   end do
   deallocate(rhor_interpol)
   deallocate(rhog_interpol)
   deallocate(coordRed_interpol)
   deallocate(coordCart_interpol)
   if (testing .and. customRho) then
    deallocate(rhor_testing)
    deallocate(rhog_testing)
   end if

   call timab(72,2,tsec)

   if (testing) then
    call system_clock(tpsStop, count_rate = countParSeconde)
    write(*,*) "Tps : ", real(tpsStop - tpsStart) / real(countParSeconde)
    write(*,*) grtn_cart_interpol
    stop
   end if

  end if
 end if

end subroutine mklocl_realspace
!!***


subroutine createIonicPotential_new(geocode,iproc,nproc,nat,ntypes,iatype,psppar,nelpsp,rxyz,gridcart,&
  hxh,hyh,hzh,elecfield,n1i,n2i,n3i,pkernel,pot_ion,spaceworld)


#if defined HAVE_BIGDFT
use Poisson_Solver
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15common, except_this_one => createIonicPotential_new
 use interfaces_lib01hidempi
!End of the abilint section

implicit none

character(len=1), intent(in) :: geocode
integer, intent(in) :: iproc,nproc,ntypes,nat,n1i,n2i,n3i,spaceworld
real(kind=8), intent(in) :: hxh,hyh,hzh,elecfield
integer, dimension(nat), intent(in) :: iatype
integer, dimension(ntypes), intent(in) :: nelpsp
real(kind=8), dimension(3,n1i*n2i*n3i), intent(in) :: gridcart
real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
real(kind=8), dimension(3,nat), intent(in) :: rxyz
real(kind=8), dimension(*), intent(in) :: pkernel
real(kind=8), dimension(*), intent(inout) :: pot_ion
!local variables
logical :: perx,pery,perz,gox,goy,goz
integer :: iat,jat,i1,i2,i3,j1,j2,j3,isx,isy,isz,iex,iey,iez,ierr,ityp,jtyp,nspin
integer :: ind,i_all,i_stat,nloc,iloc
real(kind=8) :: pi,rholeaked,dist,rloc,charge,cutoff,x,y,z,r2,arg,xp,tt,rx,ry,rz
real(kind=8) :: tt_tot,rholeaked_tot
real(kind=8) :: ehart,eexcu,vexcu
real(kind=8), dimension(4) :: charges_mpi

pi=4.d0*atan(1.d0)
! Ionic charge (must be calculated for the PS active processes)
rholeaked=0.d0
! Ionic energy (can be calculated for all the processors)

!here we should insert the calculation of the ewald energy for the periodic BC case
!!$  eion=0.d0
!!$  do iat=1,nat
!!$     ityp=iatype(iat)
!!$     rx=rxyz(1,iat) 
!!$     ry=rxyz(2,iat)
!!$     rz=rxyz(3,iat)
!!$     !    ion-ion interaction
!!$     do jat=1,iat-1
!!$        dist=sqrt( (rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2 )
!!$        jtyp=iatype(jat)
!!$        eion=eion+real(nelpsp(jtyp)*nelpsp(ityp),kind=8)/dist
!!$     enddo
!!$  end do
!!$  if (iproc.eq.0) write(*,'(1x,a,1pe22.14)') 'ion-ion interaction energy',eion

!Creates charge density arising from the ionic PSP cores
!the n3pi dimension indicates the number of planes trated by each processor in the FFT parallelisation
!for a plane wave treatment this value depends on whether the direct space is divided in planes or not
!I don't know this variable, which in the future must be inserted at the place of n3pi (LG)
!if n3pi=0 this means that the processors doesn't calculate anything
!if (n3pi >0 ) then
 
!conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')
#if defined HAVE_BIGDFT

  !this initialise the array to zero, it will work only if bigdft library is enabled
  call razero(n1i*n2i*n3i,pot_ion)

  do iat=1,nat
     ityp=iatype(iat)
     rx=rxyz(1,iat)                         
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)

     rloc=psppar(0,0,ityp)
     charge=real(nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
     cutoff=10.d0*rloc

     isx=floor((rx-cutoff)/hxh)
     isy=floor((ry-cutoff)/hyh)
     isz=floor((rz-cutoff)/hzh)

     iex=ceiling((rx+cutoff)/hxh)
     iey=ceiling((ry+cutoff)/hyh)
     iez=ceiling((rz+cutoff)/hzh)


     do i3=isz,iez
        z=real(i3,kind=8)*hzh-rz
        call ind_positions(perz,i3,n3i,j3,goz) 
        do i2=isy,iey
           y=real(i2,kind=8)*hyh-ry
           call ind_positions(pery,i2,n2i,j2,goy)
           do i1=isx,iex
              x=real(i1,kind=8)*hxh-rx
              call ind_positions(perx,i1,n1i,j1,gox)
              r2=x**2+y**2+z**2
              if (goz  .and. goy  .and. gox ) then
                 ind=j1+(j2-1)*n1i+(j3-1)*n1i*n2i
                 r2=(gridcart(1,ind)-rx)**2+(gridcart(2,ind)-ry)**2+(gridcart(3,ind)-rz)**2
              end if
              arg=r2/rloc**2
              xp=exp(-.5d0*arg)
              if (goz  .and. goy  .and. gox ) then
                 ind=j1+(j2-1)*n1i+(j3-1)*n1i*n2i
                 pot_ion(ind)=pot_ion(ind)-xp*charge
              else 
                 rholeaked=rholeaked+xp*charge
              endif
           enddo
        enddo
     enddo

  enddo

!end if
! Check
tt=0.d0
do j3= 1,n3i
  do i2= 1,n2i
     do i1= 1,n1i
        ind=i1+(i2-1)*n1i+(j3-1)*n1i*n2i
        tt=tt+pot_ion(ind)
     enddo
  enddo
enddo

tt=tt*hxh*hyh*hzh
rholeaked=rholeaked*hxh*hyh*hzh

!print *,'test case input_rho_ion',iproc,i3start,i3end,n3i,2*n3+16,tt

 call xsum_mpi(tt,tt_tot,spaceworld,ierr)
 call xsum_mpi(rholeaked,rholeaked_tot,spaceworld,ierr)

if (iproc.eq.0) write(*,'(1x,a,f26.12,2x,1pe10.3)') &
    'total ionic charge, leaked charge ',tt_tot,rholeaked_tot

!here the value of the datacode must be kept fixed
!there can be some problems when running this stuff in parallel, if the ionic potential distribution does not agree with the
!plane distribution which is supposed to hold for the Poisson Solver
call PSolver(geocode,'D',iproc,nproc,n1i,n2i,n3i,0,hxh,hyh,hzh,&
    pot_ion,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.false.,1)

print *,'ehartree',ehart
!if (n3i > 0) then
  do iat=1,nat
     ityp=iatype(iat)

     rx=rxyz(1,iat)
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)

     if (iat==1) then
        print *,rx/hxh,ry/hyh,rz/hzh
        print *,rx,ry,rz
        print *,hxh,hyh,hzh
     end if

     !stop
     ! determine number of local terms
     nloc=0
     do iloc=1,4
        if (psppar(0,iloc,ityp).ne.0.d0) nloc=iloc
     enddo
     rloc=psppar(0,0,ityp)
     cutoff=10.d0*rloc

     isx=floor((rx-cutoff)/hxh)
     isy=floor((ry-cutoff)/hyh)
     isz=floor((rz-cutoff)/hzh)

     iex=ceiling((rx+cutoff)/hxh)
     iey=ceiling((ry+cutoff)/hyh)
     iez=ceiling((rz+cutoff)/hzh)

     if (nloc /= 0) then

        print *,'nloc=',nloc

        do i3=isz,iez
           z=real(i3,kind=8)*hzh-rz
           call ind_positions(perz,i3,n3i,j3,goz) 
           if (goz) then
              do i2=isy,iey
                 y=real(i2,kind=8)*hyh-ry
                 call ind_positions(pery,i2,n2i,j2,goy)
                 if (goy) then
                    do i1=isx,iex
                       x=real(i1,kind=8)*hxh-rx
                       call ind_positions(perx,i1,n1i,j1,gox)
                       if (gox) then
                          ind=j1+(j2-1)*n1i+(j3-1)*n1i*n2i
                          r2=(gridcart(1,ind)-rx)**2+(gridcart(2,ind)-ry)**2+(gridcart(3,ind)-rz)**2
                          !r2=x**2+y**2+z**2
                          arg=r2/rloc**2
                          xp=exp(-.5d0*arg)
                          tt=psppar(0,nloc,ityp)
                          do iloc=nloc-1,1,-1
                             tt=arg*tt+psppar(0,iloc,ityp)
                          enddo
                          ind=j1+(j2-1)*n1i+(j3-1)*n1i*n2i
                          pot_ion(ind)=pot_ion(ind)+xp*tt
                       end if
                    enddo
                 end if
              enddo
           end if
        end do

     end if

  enddo
#endif
!end if

!  !use rhopot to calculate the potential from a constant electric field along x direction
!  if (elecfield /= 0.d0) then
!     !constant electric field allowed only for free BC
!     if (geocode == 'F') then
!     if (iproc.eq.0) write(*,'(1x,a)') &
!          'The constant electric field is allowed only for Free BC'
!     stop
!     end if
!     if (iproc.eq.0) write(*,'(1x,a,1pe10.2)') &
!          'Adding constant electric field of intensity',elecfield,&
!          'Ha*Bohr'
!
!     if (n3pi > 0) then
!        do i3=1,n3pi
!           do i2= -14,2*n2+16
!              do i1= -14,2*n1+16
!                 ind=i1+15+(i2+14)*(2*n1+31)+(i3-1)*(2*n1+31)*(2*n2+31)
!                 pot_ion(ind)=pot_ion(ind)+0.5d0*elecfield*hxh*real(i1-n1,kind=8)
!              enddo
!           enddo
!        enddo
!     end if
!  end if


end subroutine createIonicPotential_new
!determine the index in which the potential must be inserted, following the BC
!determine also whether the index is inside or outside the box for free BC
subroutine ind_positions(periodic,i,n,j,go)


  implicit none

  logical, intent(in) :: periodic
  integer, intent(in) :: i,n
  logical, intent(out) :: go
  integer, intent(out) :: j

  if (periodic) then
     go=.true.
     j=modulo(i-1,n)+1
  else
     j=i
     if (i >= 1 .and. i <= n) then
        go=.true.
     else
        go=.false.
     end if
  end if

end subroutine ind_positions

subroutine ext_buffers(periodic,nl,nr)


  implicit none

  logical, intent(in) :: periodic
  integer, intent(out) :: nl,nr

  if (periodic) then
     nl=0
     nr=0
  else
     nl=14
     nr=15
  end if
end subroutine ext_buffers

subroutine local_forces_new(geocode,iproc,nproc,ntypes,nat,iatype,rxyz,gridcart,psppar,nelpsp,hxh,hyh,hzh,&
     n1,n2,n3,rho,pot,floc)
! Calculates the local forces acting on the atoms belonging to iproc


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15common, except_this_one => local_forces_new
!End of the abilint section

  implicit none

  !Arguments---------
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,ntypes,nat,n1,n2,n3
  real(kind=8), intent(in) :: hxh,hyh,hzh
  real(kind=8), dimension(3,n1*n2*n3), intent(in) :: gridcart
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(*), intent(in) :: rho,pot
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(kind=8), dimension(3,nat), intent(out) :: floc
  !Local variables---------
  logical :: perx,pery,perz,gox,goy,goz
  real(kind=8) :: pi,prefactor,cutoff,rloc,Vel,rhoel
  real(kind=8) :: fxerf,fyerf,fzerf,fxion,fyion,fzion,fxgau,fygau,fzgau,forceleaked,forceloc
  real(kind=8) :: rx,ry,rz,x,y,z,arg,r2,xp,dist,tt
  integer :: isx,isy,isz,iex,iey,iez,i1,i2,i3,j1,j2,j3,ind,iat,jat,ityp,jtyp,nloc,iloc,i_all,i_stat
  !array of coefficients of the derivative
  real(kind=8), dimension(4) :: cprime 

  pi=4.d0*atan(1.d0)

  if (iproc == 0) write(*,'(1x,a)',advance='no')'Calculate local forces...'
  forceleaked=0.d0

!conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  do iat=1,nat
     ityp=iatype(iat)
     !coordinates of the center
     rx=rxyz(1,iat) 
     ry=rxyz(2,iat) 
     rz=rxyz(3,iat)

     !inizialization of the forces
     !ion-electron term, error function part
     fxerf=0.d0
     fyerf=0.d0
     fzerf=0.d0
     !ion-electron term, gaussian part
     fxgau=0.d0
     fygau=0.d0
     fzgau=0.d0

     !building array of coefficients of the derivative of the gaussian part
     cprime(1)=2.d0*psppar(0,2,ityp)-psppar(0,1,ityp)
     cprime(2)=4.d0*psppar(0,3,ityp)-psppar(0,2,ityp)
     cprime(3)=6.d0*psppar(0,4,ityp)-psppar(0,3,ityp)
     cprime(4)=-psppar(0,4,ityp)

     ! determine number of local terms
     nloc=0
     do iloc=1,4
        if (psppar(0,iloc,ityp).ne.0.d0) nloc=iloc
     enddo

     !local part
     rloc=psppar(0,0,ityp)
     prefactor=real(nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**5)
     !maximum extension of the gaussian
     cutoff=10.d0*rloc
     isx=floor((rx-cutoff)/hxh)
     isy=floor((ry-cutoff)/hyh)
     isz=floor((rz-cutoff)/hzh)

     iex=ceiling((rx+cutoff)/hxh)
     iey=ceiling((ry+cutoff)/hyh)
     iez=ceiling((rz+cutoff)/hzh)

     !calculate the forces near the atom due to the error function part of the potential
     do i3=isz,iez
        z=real(i3,kind=8)*hzh-rz
        call ind_positions(perz,i3,n3,j3,goz) 
        do i2=isy,iey
           y=real(i2,kind=8)*hyh-ry
           call ind_positions(pery,i2,n2,j2,goy)
           do i1=isx,iex
              x=real(i1,kind=8)*hxh-rx
              call ind_positions(perx,i1,n1,j1,gox)
              r2=x**2+y**2+z**2
              if (goz  .and. goy  .and. gox ) then
                 ind=j1+(j2-1)*n1+(j3-1)*n1*n2
                 x=(gridcart(1,ind)-rx)
                 y=(gridcart(2,ind)-ry)
                 z=(gridcart(3,ind)-rz)
                 r2=x**2+y**2+z**2
              end if
              arg=r2/rloc**2
              xp=exp(-.5d0*arg)
              if (goz  .and. goy  .and. gox ) then
                 ind=j1+(j2-1)*n1+(j3-1)*n1*n2
                !gaussian part
                 if (nloc /= 0) then
                    tt=cprime(nloc)
                    do iloc=nloc-1,1,-1
                       tt=arg*tt+cprime(iloc)
                    enddo
                    rhoel=rho(ind)
                    forceloc=xp*tt*rhoel
                    fxgau=fxgau+forceloc*x
                    fygau=fygau+forceloc*y
                    fzgau=fzgau+forceloc*z
                 end if
                 !error function part
                 Vel=pot(ind)
                 fxerf=fxerf+xp*Vel*x
                 fyerf=fyerf+xp*Vel*y
                 fzerf=fzerf+xp*Vel*z
              else
                 forceleaked=forceleaked+xp*(1.d0+tt)
              endif
           end do
        end do
     end do

     !final result of the forces

     floc(1,iat)=(hxh*hyh*hzh*prefactor)*fxerf+(hxh*hyh*hzh/rloc**2)*fxgau
     floc(2,iat)=(hxh*hyh*hzh*prefactor)*fyerf+(hxh*hyh*hzh/rloc**2)*fygau
     floc(3,iat)=(hxh*hyh*hzh*prefactor)*fzerf+(hxh*hyh*hzh/rloc**2)*fzgau

  end do

  forceleaked=forceleaked*prefactor*hxh*hyh*hzh
  if (iproc.eq.0) write(*,'(a,1pe12.5)') 'done. Leaked force: ',forceleaked

end subroutine local_forces_new
