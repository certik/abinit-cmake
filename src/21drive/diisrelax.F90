!{\src2tex{textfont=tt}}
!!****f* ABINIT/diisRelax
!! NAME
!! diisRelax
!!
!! FUNCTION
!! Given a starting point xred that is a vector of length 3*natom
!! (reduced nuclei coordinates), and unit cell parameters (acell and rprim) this routine uses
!! the DIIS (direct inversion of the iterative subspace) to minize the
!! gradient (forces) on atoms. The preconditioning used to compute errors from
!! gradients is using an inversed hessian matrix obtained by a BFGS algorithm.
!! This method is known to converge to the nearest
!! point where gradients vanish. This is efficient to refine positions
!! around a saddle point for instance.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2006-2008 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cpus= cpu time limit in seconds
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs for the "coarse" grid (see NOTES below)
!!   | mkmem =number of k points which can fit in memory; set to 0 if use disk
!!   | mpw=maximum dimensioned size of npw.
!!   | natom=number of atoms in cell.
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   |      for the "coarse" grid (see NOTES below)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  ecore=core psp energy (part of total energy) (hartree)
!!  iapp=indicates the eventual suffix to be appended to the generic output root
!!         if 0 : no suffix to be appended (called directly from gstate)
!!         if positive : append "_TIM//iapp" (called from move or brdmin)
!!         if -1 : append "_TIM0" (called from brdmin)
!!         if -2, -3, -4, -5: append "_TIMA", ... ,"_TIMD", (called from move)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)= # atoms of each type.
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!       for the "fine" grid (see NOTES below)
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points and spins
!!
!! SIDE EFFECTS
!!  acell(3)=length scales of primitive translations (bohr)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=updated wavefunctions;  if mkmem>=nkpt,
!!         these are kept in a disk file.
!!  densymop_gs <type(dens_sym_operator_type)>=the density symmetrization
!!   operator (ground-state symmetries)
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations (see initberry.f)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  initialized= if 0 the initialization of the gstate run is not yet finished
!!  irrzon(nfft**(1-1/nsym),2,nspden/nsppol)=irreducible zone data
!!  occ(mband*nkpt*nsppol)=occupation number for each band (often 2) at each k point
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  phnons(2,nfft**(1-1/nsym),nspden/nsppol)=nonsymmorphic translation phases
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!   (should be made a pure output quantity)
!!  rhog(2,nfftf)=array for Fourier transform of electron density
!!  rhor(nfftf,nspden)=array for electron density in el./bohr**3
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  wffnew,wffnow=struct info for wf disk files.
!!  wvl <type(wvl_data)>=all wavelets data.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)= at input, previous reduced dimensionless atomic coordinates
!!                     at output, current xred is transferred to xred_old
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine diisRelax(acell, atindx, atindx1, cg, cpus, densymop_gs, dtefield, &
                   & dtfil, dtset, ecore, eigen, hdr, iapp, indsym, initialized, &
                   & irrzon, kg, mpi_enreg, nattyp, nfftf, npwarr, nspinor, occ, pawang, &
                   & pawfgr, pawrad, pawrhoij, pawtab, phnons, psps, pwind, pwind_alloc, pwnsfac, &
                   & resid, results_gs, rhog, rhor, rprimd, scf_history, symrec, &
                   & wffnew, wffnow, wvl, xred, xred_old, ylm, ylmgr)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_15common
 use interfaces_16geomoptim
 use interfaces_21drive, except_this_one => diisRelax
 use interfaces_linalg
!End of the abilint section

 implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: iapp,pwind_alloc
  integer, intent(inout) :: nfftf
  integer,intent(inout) :: initialized,nspinor
  real(dp),intent(in) :: cpus,ecore
  type(MPI_type),intent(inout) :: mpi_enreg
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  type(dens_sym_operator_type),intent(inout) :: densymop_gs
  type(efield_type),intent(inout) :: dtefield
  type(hdr_type),intent(inout) :: hdr
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(results_gs_type),intent(inout) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(wffile_type),intent(inout) :: wffnew,wffnow
  type(wvl_data), intent(inout) :: wvl
  !arrays
  integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
  integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
  !no_abirules
  integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
   !(nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise)
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat),npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
  integer, intent(inout) :: symrec(3,3,dtset%nsym)
  real(dp), intent(inout) :: acell(3),rprimd(3,3)
  real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
   !(nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), pointer :: rhog(:,:),rhor(:,:)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(inout) :: xred(3,dtset%natom),xred_old(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
  type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  !Local variables-------------------------------
  integer :: iStep, id, jd, currentId, j, info, iexit, rank, start
  integer, allocatable :: ipiv(:)
  real(dp) :: norm, ident(3, 3)
  character(len=32) :: statusOut
  character(len=500) :: message
  real(dp), parameter :: alphaPrecond = 0.005_dp
  real(dp), allocatable :: error(:, :, :)
  real(dp), allocatable :: savedCoord(:, :, :), savedGradient(:, :, :)
  real(dp), allocatable :: diisMatrix(:, :)
  real(dp), allocatable :: diisCoeff(:)
  real(dp), allocatable :: workArray(:)
  real(dp), allocatable :: workMatrix(:, :)
  real(dp), allocatable :: vel(:, :), xcart(:, :), gcart(:, :)
  real(dp), allocatable :: hessianInv(:, :)
!************************************************************************
#ifdef __VMS
!DEC$ ATTRIBUTES ALIAS:'DDOT' :: ddot
!DEC$ ATTRIBUTES ALIAS:'DCOPY' :: dcopy
!DEC$ ATTRIBUTES ALIAS:'DSYSV' :: dsysv
!DEC$ ATTRIBUTES ALIAS:'DGEMV' :: dgemv
#endif
!Allocate working arrays
!error: store the supposed error for each last dtset%ngeohist steps.
!it is required to compute the DIIS matrix.
!savedCoord: store the coordinates for each last dtset%ngeohist steps.
!it is used to compute the new positions.
!savedGradients: store cartesian forces.
!diisMatrix: store the matrix used to compute the coefficients.
!diisCoeff: store the coefficients computed from diisMatrix.
!workMatrix: Lapack work array.
!workArray: Lapack work array.
!ipiv: Lapack work array.
!xcart: store cartesian coordinates.
!hessianInv: store the inversed Hessian matrix (computed with a BFGS algorithm).
 allocate(error(3, dtset%natom, dtset%ngeohist))
 start = min(dtset%ngeohist - 1, 1)
 allocate(savedCoord(3, dtset%natom, start:dtset%ngeohist))
 allocate(savedGradient(3, dtset%natom, start:dtset%ngeohist))
 allocate(diisMatrix(dtset%ngeohist + 1, dtset%ngeohist + 1))
 allocate(workMatrix(dtset%ngeohist + 1, dtset%ngeohist + 1))
 allocate(diisCoeff(dtset%ngeohist + 1))
 allocate(ipiv(dtset%ngeohist + 1))
 allocate(xcart(3, dtset%natom))
 diisMatrix(:, :) = real(0, dp)
 allocate(hessianInv(3 * dtset%natom, 3 * dtset%natom))
!Transform xred to cartesian coordinates.
 call xredxcart(dtset%natom, 1, rprimd, xcart, xred)
!Initialise the Hessian matrix with identity.
 ident(:, :) = real(0, dp)
 ident(1, 1) = real(-1, dp)
 ident(2, 2) = real(-1, dp)
 ident(3, 3) = real(-1, dp)
 call hessinit(dtfil, dtset, hessianInv, ident, 3 * dtset%natom, real(0., dp))
!Begin DIIS relaxation.
 do iStep = 1, dtset%ntime, 1
  write(message, '(a,a,i3,a)' ) ch10, ' DIIS STEP NUMBER ', iStep,&
&  '  ---------------------------------------------------------'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
! Copy error and coordinates (2,dtset%ngeohist) -> (1,dtset%ngeohist-1)
! this is less efficient than a modulo, but simpler and not important
! for geometry optimisation.
  if (iStep > dtset%ngeohist) then
   do id = start + 1, dtset%ngeohist, 1
    savedGradient(:, :, id - 1) = savedGradient(:, :, id)
    savedCoord(:, :, id - 1)    = savedCoord(:, :, id)
!   do j = 2, id, 1
!   diisMatrix(j - 1, id - 1) = diisMatrix(j, id)
!   diisMatrix(id - 1, j - 1) = diisMatrix(id, j)
!   end do
   end do
  end if
  currentId = min(iStep, dtset%ngeohist)
! We store some informations : current positions.
  savedCoord(:, :, currentId) = xcart(:, :)
! Output coordinates (not velocities, prtvel = 0) and total energy
! MUST BE REMOVED
! call prtxvf(results_gs%fcart, dtset%iatfix, 06, dtset%natom, &
! & 0, vel, savedCoord(:, :, currentId))
! MUST BE REMOVED
! Compute the gradient for the current coordinates
  call xredxcart(dtset%natom, -1, rprimd, savedCoord(:, :, currentId), xred)
  call xredxcart(dtset%natom, -1, rprimd, savedCoord(:, :, max(1, currentId - 1)), xred_old)
  call scfcv(acell,atindx, atindx1, cg, cpus, densymop_gs, dtefield, dtfil, &
&  dtset, ecore, eigen, hdr, iapp, indsym, initialized,&
&  irrzon, kg, mpi_enreg, nattyp, nfftf, npwarr, nspinor, occ,&
&  pawang, pawfgr, pawrad, pawrhoij, pawtab,&
&  phnons, psps, pwind, pwind_alloc, pwnsfac, resid, results_gs, rhog, &
&  rhor, rprimd, scf_history, symrec, wffnew, wffnow, wvl, &
&  xred, xred_old, ylm, ylmgr)
! results_gs%fcart(:, :) = 0.d0
! results_gs%fcart(1, 2) = -40 * (xred(1, 2) - 0.5) ** 3 + 0.2 * (xred(1, 2) - 0.5)
! Output coordinates and forces (not velocities, prtvel = 0) and total energy
  call prtxvf(results_gs%fcart, dtset%iatfix, ab_out, dtset%natom, &
&  0, vel, savedCoord(:, :, currentId))
  call prtxvf(results_gs%fcart, dtset%iatfix, 06 , dtset%natom, &
&  0, vel, savedCoord(:, :, currentId))
! Check whether forces and stresses are below tolerance; if so, exit
! from the itime loop
  iexit = 0
  if (iStep == dtset%ntime) then
   iexit = 1
   statusOut = "Failed"
  else
   statusOut = "OK"
  end if
! Don't take care of stress, put optcell to 0.
  call fconv(results_gs%fcart, dtset%iatfix, iexit, iStep, dtset%natom, dtset%ntime,&
&  0, dtset%strfact, dtset%strtarget, results_gs%strten, dtset%tolmxf)
  if (iexit /= 0) then
   exit
  end if
! We store some informations : current forces.
  savedGradient(:, :, currentId) = results_gs%fcart(:, :)
! do id = 1, 3 * dtset%natom, 1
! write(0, "(24F4.1)") hessianInv(:, id)
! end do
! Precondition the error using the hessian matrix.
! In the quadratic approximation, we have:
! e = H^-1.g, where e is the error vectors, H the hessian
! and g the gradient.
  do id = 1, currentId, 1
   call DGEMV('N', 3 * dtset%natom, 3 * dtset%natom, real(1, dp), hessianInv, &
&   3 * dtset%natom, savedGradient(1, 1, id), 1, real(0, dp), error(1, 1, id), 1)
  end do
! write(0, *) "Hessianized error vectors"
! do id = 1, currentId, 1
! write(*,*) id, error(1, 2, id)
! write(0, "(24F4.1)") error(:, :, id)
! end do
! Create the DIIS matrix
  write(0,*) "DIIS matrix", currentId
  diisMatrix(currentId + 1, currentId + 1) = real(0, dp)
  do id = 1, currentId, 1
   do jd = id, currentId, 1
    diisMatrix(jd, id) = ddot(3 * dtset%natom, error(1, 1, id), 1, error(1, 1, jd), 1)
    diisMatrix(id, jd) = diisMatrix(jd, id)
   end do
   diisMatrix(id, currentId + 1) = real(-1, dp)
   diisMatrix(currentId + 1, id) = real(-1, dp)
   write(0,*) diisMatrix(1:currentId + 1, id)
  end do
  write(0,*) diisMatrix(1:currentId + 1, currentId + 1)
! Solve the system using Lapack
  diisCoeff(:) = real(0, dp)
  diisCoeff(currentId + 1) = real(-1, dp)
  write(0,*) "B vector", diisCoeff(1:currentId + 1)
  allocate(workArray((currentId + 1) ** 2))
  call DCOPY((currentId + 1) ** 2, diisMatrix(1:currentId + 1, 1:currentId + 1), 1, workMatrix, 1)
  	call DSYSV('L', currentId + 1, 1, workMatrix(1, 1), &
&  currentId + 1, ipiv, diisCoeff(1), currentId + 1, &
&  workArray, (currentId + 1) ** 2, info)
  deallocate(workArray)
  if (info /= 0) then
   write(*,*) "error solving DIIS matrix", info
   write(*,*) workMatrix(1:currentId + 1, 1:currentId + 1)
   exit
  end if
  write(0,*) "DIIS coeff :", diisCoeff(1:currentId + 1)
! Build the new coordinates, to do it, we compute a new error e,
! using the linear coefficients (temporary store it in error) applied
! on previous gradient: e=H^-1(sum_i c_i.g_i)
  xcart(:, :) = real(0, dp)
  error(:, :, currentId) = real(0, dp)
  do id = 1, currentId, 1
   xcart(:, :) = xcart(:, :) + savedCoord(:, :, id) * diisCoeff(id)
   error(:, :, currentId) = error(:, :, currentId) + savedGradient(:, :, id) * diisCoeff(id)
  end do
! write(0, *) "Linear gradients"
! write(*,*) error(1, 2, currentId)
! write(0, "(24F4.1)") error(:, :, currentId)
  call DGEMV('N', 3 * dtset%natom, 3 * dtset%natom, real(-1, dp), hessianInv, &
&  3 * dtset%natom, error(1, 1, currentId), 1, real(1, dp), xcart(1, 1), 1)
! if (iStep > 10) then
! stop
! end if
! We update the hessian matrix using a BFGS algorithm.
  if (iStep > 1) then
   call hessupdt(hessianInv, dtset%iatfix, dtset%natom, 3 * dtset%natom, &
&   reshape(savedCoord(:, :, currentId), (/ 3 * dtset%natom /)), &
&   reshape(savedCoord(:, :, currentId - 1), (/ 3 * dtset%natom /)), &
&   reshape(savedGradient(:, :, currentId), (/ 3 * dtset%natom /)), &
&   reshape(savedGradient(:, :, currentId - 1), (/ 3 * dtset%natom /)))
  end if
 end do
 call xredxcart(dtset%natom, -1, rprimd, xcart, xred)
!Free working arrays
 deallocate(error)
 deallocate(savedCoord)
 deallocate(savedGradient)
 deallocate(diisMatrix)
 deallocate(workMatrix)
 deallocate(diisCoeff)
 deallocate(ipiv)
 deallocate(xcart)
 deallocate(hessianInv)
!XML output of the status
 if (mpi_enreg%me == 0 .and. dtset%outputXML == 1) then
  write(ab_xml_out, "(A)") '    <geometryMinimisation type="diis">'
  write(ab_xml_out, "(A,A,A)") '      <status cvState="', trim(statusOut) , &
&  '" stop-criterion="tolmxf" />'
  write(ab_xml_out, "(A)") '    </geometryMinimisation>'
 end if
end subroutine diisRelax
!!***
