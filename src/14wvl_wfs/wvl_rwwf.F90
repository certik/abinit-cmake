!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_outwf
!! NAME
!! wvl_outwf
!!
!! FUNCTION
!! Simple wrapper around the read/write disk methods of BigDFT for wavefunctions.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=input variables.
!!  mpi_enreg=informations about MPI parallelization
!!  option= -2 for write with BigDFT format,
!!          -1 for reading wavelets coefficients with BigDFT format,
!!          2 for write,
!!          1 for read.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  wff <type(wffile_type)>=struct info for wavefunction
!!  wfs <type(wvl_wf_type)>=wavefunctions informations for wavelets.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  xred(3,natom)=reduced dimensionless atomic coordinates (in fact IN but here
!!                because of INOUT xredxcart() behavior).
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wvl_read(dtset, hdr0, hdr, mpi_enreg, option, rprimd, wff, wfs, xred)

  use defs_basis
  use defs_datatypes
  use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API
#endif
#if defined HAVE_ETSF_IO
  use etsf_io
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
!End of the abilint section

  implicit none

!Arguments -------------------------------
  !scalars
  integer, intent(in)                       :: option
  type(dataset_type), intent(in)            :: dtset
  type(hdr_type), intent(in)                :: hdr0
  type(hdr_type), intent(in)                :: hdr
  type(MPI_type), intent(in)                :: mpi_enreg
  type(wffile_type),intent(in)              :: wff
  type(wvl_wf_type), intent(inout)          :: wfs
  !arrays
  real(dp), intent(in)                      :: rprimd(3, 3)
  real(dp), intent(in)                      :: xred(3, dtset%natom)

!Local variables-------------------------------
  character(len = 500)  :: message
  logical               :: reformat
  integer               :: i, i1, i2, i3, iBand, iCoeff, iCoarse, iFine
  integer, allocatable  :: coeff_map(:,:,:)
  integer, allocatable, target  :: wvl_coord(:,:), wvl_ncoeff(:)
  real(dp), allocatable, target :: wvl_band(:)
  integer               :: coord(3), n_old(3)
  real(dp), allocatable :: xcart(:,:), psifscf(:,:,:), psigold(:,:,:,:,:,:)
  real(dp), allocatable :: xred_(:,:), xcart_old(:,:)
  real(dp)              :: center(3), center_old(3), hgrid_old
#if defined HAVE_ETSF_IO
  character(len=etsf_io_low_error_len)  :: errmess
  type(etsf_basisdata) :: basis_folder
  type(etsf_main)      :: main_folder
  type(etsf_electrons) :: electrons_folder
  type(etsf_geometry)  :: geometry_folder
  logical              :: lstat
  type(etsf_io_low_error) :: error
#endif

  if (abs(option) /= 1) then
    write(message,'(a,a,a,a,a,a,i5,a)') ch10,&
  &  ' wvl_read: BUG -', ch10, &
  &  '  Option argument is wrong,', ch10, &
  &  '  awaited values are -1 or  1 but option = ', option, '.'
    call wrtout(6,message,'PERS')
    call leave_new('PERS')
  end if

#if defined HAVE_BIGDFT
  ! Store xcart for each atom
  allocate(xcart(3, dtset%natom))
  allocate(xred_(3, dtset%natom))
  xred_ = xred
  call xredxcart(dtset%natom, 1, rprimd, xcart, xred_)

  write(message, '(a,a,a,a)' ) ch10,&
       &  ' wvl_read:  read wavefunctions from file.'
  call wrtout(6,message,'COLL')

  if (option > 0) then
     ! Read in the ABINIT way.
     if (wff%accesswff == 0 .or. (wff%accesswff == -1 .and. wff%master==wff%me)) then
#if defined MPI
        ! MPI is not taken into account yet.
        write(message, '(a,a,a,a)' ) ch10,&
             &  ' wvl_read:  no MPI implementation yet.'
        call wrtout(6,message,'COLL')
#else
3        allocate(psifscf(dtset%wvl_internal%dpSize(1), &
             & dtset%wvl_internal%dpSize(2), &
             & dtset%wvl_internal%nSize(3)))
        do iBand = 1, dtset%mband, 1
           call readonewave(wff%unwff, .false., iBand, mpi_enreg%me, &
                & dtset%wvl_internal%nSize(1), dtset%wvl_internal%nSize(2), &
                & dtset%wvl_internal%nSize(3), dtset%wvl_hgrid, &
                & xcart(:, 1), wfs%keys%nseg_c, wfs%keys%nseg_f, wfs%keys%nvctr_c, &
                & wfs%keys%nvctr_f, wfs%keys%keyg, wfs%keys%keyv, wfs%psi, &
                & wfs%eval(iBand), psifscf)
        end do
        deallocate(psifscf)
#endif
#if defined HAVE_ETSF_IO
     else if (wff%accesswff == 3) then
        ! Read a NetCDF file.
        ! coordinates_of_grid_points is used to store the geometric
        !  position of coefficients of wavelets i, as integer in
        ! dtset%wvl_internal%dpSize(:) dimensions.
        ! coefficients_of_wavefunctions is used to store the psi values for
        !  each wavelet.

        ! We check if we need reformating
        if (abs(hdr%rprimd(1,1) / hdr%ngfft(1) - &
             & hdr0%rprimd(1,1) / hdr0%ngfft(1)) < tol8 .and. &
             & maxval(abs(hdr%nwvlarr - hdr0%nwvlarr)) == 0 .and. & 
             & maxval(abs(hdr%ngfft   - hdr0%ngfft  )) == 0 ) then
           reformat = .false.
           write(message, '(a,a,a,a)' ) ch10,&
                &  ' wvl_read:  wavefunctions needs NO reformatting.'
           call wrtout(6,message,'COLL')           

           ! We allocate a temporary array to read the wavefunctions
           ! and reorder then into wfs%psi
           allocate(wvl_band(sum(hdr%nwvlarr)))
        else
           reformat = .true.
           write(message, '(a,a,a,a)' ) ch10,&
                &  ' wvl_read:  wavefunctions needs reformatting.'
           call wrtout(6,message,'COLL')

           ! Compute center of atoms
           center(:) = zero
           do i = 1, dtset%natom, 1
              center(:) = center(:) + xcart(:, i)
           end do
           center(:) = center(:) / real(dtset%natom, dp)
           ! Also for old system
           allocate(xcart_old(3, dtset%natom))
           xred_ = hdr0%xred
           call xredxcart(dtset%natom, 1, hdr0%rprimd, xcart_old, xred_)
           center_old(:) = zero
           do i = 1, dtset%natom, 1
              center_old(:) = center_old(:) + xcart_old(:, i)
           end do
           center_old(:) = center_old(:) / real(dtset%natom, dp)
           deallocate(xcart_old)

           ! We allocate a temporary array to read the wavefunctions
           ! and reorder then into wfs%psi
           allocate(wvl_band(sum(hdr0%nwvlarr)))
           ! Do the reformating by expanding old wavefunctions on the grid
           ! and interpolate it on the new grid using BigDFT reformatonewave().
           n_old = (hdr0%ngfft(1:3) - dtset%wvl_internal%buffer) / 2
           hgrid_old = hdr0%rprimd(1,1) / n_old(1)
           allocate(psigold(0:n_old(1), 2, 0:n_old(2), 2, 0:n_old(3), 2))
           allocate(psifscf(-7:dtset%wvl_internal%nSize(1) * 2 + 8, &
                & -7:dtset%wvl_internal%nSize(2) * 2 + 8, &
                & -7:dtset%wvl_internal%nSize(3) * 2 + 8))
           allocate(basis_folder%coordinates_of_basis_grid_points%data2D(3, &
                & hdr0%nwvlarr(1)))
        end if

        ! We read the basis set definition
        allocate(basis_folder%number_of_coefficients_per_grid_point%data1D(hdr0%nwvlarr(1)))
        call etsf_io_basisdata_get(wff%unwff, basis_folder, lstat, error)
        if (.not. lstat) then
           call etsf_io_low_error_to_str(errmess, error)
           write(message, "(A,A,A,A)") ch10, " wvl_read: ERROR -", ch10, &
                & errmess(1:min(475, len(errmess)))
           call wrtout(06, message, 'COLL')
        end if

        ! We just read the file from disk to memory, band per band.
        do iBand = 1, dtset%mband, 1
           main_folder%wfs_coeff__state_access = iBand
           main_folder%coefficients_of_wavefunctions%data1D => wvl_band
           call etsf_io_main_get(wff%unwff, main_folder, lstat, error)
           if (.not. lstat) then
              call etsf_io_low_error_to_str(errmess, error)
              write(message, "(A,A,A,A)") ch10, " wvl_read: ERROR -", ch10, &
                   & errmess(1:min(475, len(errmess)))
              call wrtout(06, message, 'COLL')
           end if

           ! Now we reorder
           if (reformat) then
              ! Expand the wf on the old grid.
              psigold = zero
              iCoeff = 1
              do i = 1, hdr0%nwvlarr(1), 1
                 coord = basis_folder%coordinates_of_basis_grid_points%data2D(:, i) / 2
                 psigold(coord(1), 1, coord(2), 1, coord(3), 1) = wvl_band(iCoeff)
                 iCoeff = iCoeff + 1
                 if (basis_folder%number_of_coefficients_per_grid_point%data1D(i) == 8) then
                    psigold(coord(1), 2, coord(2), 1, coord(3), 1) = wvl_band(iCoeff + 0)
                    psigold(coord(1), 1, coord(2), 2, coord(3), 1) = wvl_band(iCoeff + 1)
                    psigold(coord(1), 2, coord(2), 2, coord(3), 1) = wvl_band(iCoeff + 2)
                    psigold(coord(1), 1, coord(2), 1, coord(3), 2) = wvl_band(iCoeff + 3)
                    psigold(coord(1), 2, coord(2), 1, coord(3), 2) = wvl_band(iCoeff + 4)
                    psigold(coord(1), 1, coord(2), 2, coord(3), 2) = wvl_band(iCoeff + 5)
                    psigold(coord(1), 2, coord(2), 2, coord(3), 2) = wvl_band(iCoeff + 6)
                    iCoeff = iCoeff + 7
                 end if
              end do
              call reformatonewave(mpi_enreg%me, hgrid_old, n_old(1), n_old(2), &
                   & n_old(3), center_old, psigold, dtset%wvl_hgrid, &
                   & wfs%keys%nvctr_c, wfs%keys%nvctr_f, &
                   & dtset%wvl_internal%nSize(1), dtset%wvl_internal%nSize(2), &
                   & dtset%wvl_internal%nSize(3), center, wfs%keys%nseg_c, &
                   & wfs%keys%nseg_f, wfs%keys%keyg, wfs%keys%keyv, psifscf, &
                   & wfs%psi(:, iBand - mpi_enreg%me * wfs%mbandp))
           else
              ! No reformating, just reorder
              iCoeff = 1
              iFine = wfs%keys%nvctr_c + 1
              do iCoarse = 1, hdr0%nwvlarr(1), 1
                 wfs%psi(iCoarse, iBand - mpi_enreg%me * wfs%mbandp) = &
                      & wvl_band(iCoeff)
                 iCoeff  = iCoeff  + 1
                 if (basis_folder%number_of_coefficients_per_grid_point%data1D(iCoarse) == 8) then
                    wfs%psi(iFine:iFine + 6, iBand - mpi_enreg%me * wfs%mbandp) = &
                         & wvl_band(iCoeff:iCoeff + 6)
                    iFine  = iFine  + 7
                    iCoeff = iCoeff + 7
                 end if
              end do
           end if
        end do

        ! We read wfs%eval (TODO maybe removed later).
        electrons_folder%eigenvalues%data1D => wfs%eval
        call etsf_io_electrons_get(wff%unwff, electrons_folder, lstat, error)
        if (.not. lstat) then
           call etsf_io_low_error_to_str(errmess, error)
           write(message, "(A,A,A,A)") ch10, " wvl_read: ERROR -", ch10, &
                & errmess(1:min(475, len(errmess)))
           call wrtout(06, message, 'COLL')
        end if
        
        ! Deallocate temporary arrays
        deallocate(wvl_band)
        deallocate(basis_folder%number_of_coefficients_per_grid_point%data1D)
        if (reformat) then
           deallocate(basis_folder%coordinates_of_basis_grid_points%data2D)
           deallocate(psigold)
           deallocate(psifscf)
        end if
#endif

     else
        write(message,'(a,a,a,a,a,a,i5,a)') ch10,&
             &  ' wvl_read: BUG -', ch10, &
             &  '  wff%accesswff argument is wrong,', ch10, &
             &  '  awaited values are -1, 0 (or 3 if netcdf/etsf_io is available) but value = ', wff%accesswff, '.'
        call wrtout(6,message,'PERS')
        call leave_new('PERS')
     end if
  else
     call readmywaves(mpi_enreg%me, dtset%mband, wfs%mbandp, &
          & dtset%wvl_internal%nSize(1), dtset%wvl_internal%nSize(2), &
          & dtset%wvl_internal%nSize(3), dtset%wvl_hgrid, dtset%natom, &
          & xcart(:, 1), wfs%keys%nseg_c, wfs%keys%nseg_f, wfs%keys%nvctr_c, &
          & wfs%keys%nvctr_f, wfs%keys%keyg, wfs%keys%keyv, wfs%psi, wfs%eval)
  end if

  deallocate(xcart)
  deallocate(xred_)
#else
  write(message, '(a,a,a,a)' ) ch10,&
    &  ' wvl_outwf : BigDFT library is not compiled.', ch10, &
    &  '   Action, used the flag --enable-bigdft when configuring.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
#endif

end subroutine wvl_read
!!***

subroutine wvl_write(dtset, eigen, mpi_enreg, option, rprimd, wff, wfs, xred)

  use defs_basis
  use defs_datatypes
  use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API
#endif
#if defined HAVE_ETSF_IO
  use etsf_io
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
!End of the abilint section

  implicit none

!Arguments -------------------------------
  !scalars
  integer, intent(in)                       :: option
  type(dataset_type), intent(in)            :: dtset
  type(MPI_type), intent(in)                :: mpi_enreg
  type(wffile_type),intent(in)              :: wff
  type(wvl_wf_type), intent(in)             :: wfs
  !arrays
  real(dp), intent(in), target              :: eigen(dtset%mband)
  real(dp), intent(in)                      :: rprimd(3, 3)
  real(dp), intent(inout)                   :: xred(3, dtset%natom)

!Local variables-------------------------------
  character(len = 500)  :: message
  integer               :: iGrid, iCoeff, iCoarse, iFine
  integer               :: iorb, jj, j0, j1, i, ii, i0, i1, i2, i3
  integer               :: iseg, nseg, ipsi, npsi
  integer               :: buffer, iband
  integer, allocatable  :: coeff_map(:,:,:)
  integer, allocatable, target :: wvl_coord(:,:), wvl_ncoeff(:)
  real(dp), allocatable, target :: wvl_band(:)
  real(dp), allocatable :: xcart(:,:)
#if defined HAVE_ETSF_IO
  character(len=etsf_io_low_error_len)  :: errmess
  type(etsf_basisdata) :: basis_folder
  type(etsf_main)      :: main_folder
  type(etsf_electrons) :: electrons_folder
  logical              :: lstat
  type(etsf_io_low_error) :: error
#endif


  if (abs(option) /= 2) then
    write(message,'(a,a,a,a,a,a,i5,a)') ch10,&
  &  ' wvl_write: BUG -', ch10, &
  &  '  Option argument is wrong,', ch10, &
  &  '  awaited values are -2 or  2 but option = ', option, '.'
    call wrtout(6,message,'PERS')
    call leave_new('PERS')
  end if

#if defined HAVE_BIGDFT
  ! Store xcart for each atom
  allocate(xcart(3, dtset%natom))
  call xredxcart(dtset%natom, 1, rprimd, xcart, xred)

  write(message, '(a,a,a,a)' ) ch10,&
       &  ' wvl_write:  Write wavefunctions to file.'
  call wrtout(6,message,'COLL')

  if (option > 0) then
     ! Write in the ABINIT way.
     if (wff%accesswff == 0 .or. (wff%accesswff == -1 .and. wff%master==wff%me)) then
#if defined MPI
        ! MPI is not taken into account yet.
        write(message, '(a,a,a,a)' ) ch10,&
             &  ' wvl_write:  no MPI implementation yet.'
        call wrtout(6,message,'COLL')
#else
        iseg = wfs%keys%nseg_c
        nseg = wfs%keys%nseg_c + wfs%keys%nseg_f
        ipsi = wfs%keys%nvctr_c
        npsi = wfs%keys%nvctr_c + 7 * wfs%keys%nvctr_f
        do iorb = 1, dtset%mband
           call writeonewave(wff%unwff, .false., iorb, dtset%wvl_internal%nSize(1), &
                & dtset%wvl_internal%nSize(2), dtset%wvl_internal%nSize(3), &
                & dtset%wvl_hgrid, xcart(:, 1), wfs%keys%nseg_c, wfs%keys%nvctr_c, &
                & wfs%keys%keyg(:,1:iseg), wfs%keys%keyv(1:iseg), wfs%keys%nseg_f, &
                & wfs%keys%nvctr_f, wfs%keys%keyg(:, iseg + 1:nseg), &
                & wfs%keys%keyv(iseg + 1:nseg), wfs%psi(1:ipsi, iorb), &
                & wfs%psi(ipsi + 1:npsi, iorb), dtset%mband, wfs%eval)
        enddo
#endif
#if defined HAVE_ETSF_IO
     else if (wff%accesswff == 3) then
        ! Write a NetCDF file.
        ! coordinates_of_grid_points is used to store the geometric
        !  position of coefficients of wavelets i, as integer in
        ! dtset%wvl_internal%dpSize(:) dimensions.
        ! coefficients_of_wavefunctions is used to store the psi values for
        !  each wavelet.

        ! We write the basis set.
        ! =======================
        allocate(wvl_coord(3, wfs%keys%nvctr_c))
        allocate(wvl_ncoeff(wfs%keys%nvctr_c))
        ! Will store the grid index for a given geometric point
        allocate(coeff_map(dtset%wvl_internal%nSize(1), &
             & dtset%wvl_internal%nSize(2), dtset%wvl_internal%nSize(3)))
        ! coarse part
        iGrid = 0
        do iseg = 1, wfs%keys%nseg_c
           jj = wfs%keys%keyv(iseg)
           j0 = wfs%keys%keyg(1, iseg)
           j1 = wfs%keys%keyg(2, iseg)
           ii = j0 - 1
           i3 = ii / ((dtset%wvl_internal%nSize(1) + 1) * &
                & (dtset%wvl_internal%nSize(2) + 1))
           ii = ii - i3 * (dtset%wvl_internal%nSize(1) + 1) * &
                & (dtset%wvl_internal%nSize(2) + 1)
           i2 = ii / (dtset%wvl_internal%nSize(1) + 1)
           i0 = ii - i2 * (dtset%wvl_internal%nSize(1) + 1)
           i1 = i0 + j1 - j0
           do i = i0, i1
              iGrid = iGrid + 1
              coeff_map(i, i2, i3) = iGrid
              wvl_coord(:, iGrid) = (/ i * 2, i2 * 2, i3 * 2 /)
              wvl_ncoeff(iGrid) = 1
           enddo
        enddo
        ! fine part
        do iseg = 1, wfs%keys%nseg_f
           jj = wfs%keys%keyv(wfs%keys%nseg_c + iseg)
           j0 = wfs%keys%keyg(1, wfs%keys%nseg_c + iseg)
           j1 = wfs%keys%keyg(2, wfs%keys%nseg_c + iseg)
           ii = j0 - 1
           i3 = ii / ((dtset%wvl_internal%nSize(1) + 1) * &
                & (dtset%wvl_internal%nSize(2) + 1))
           ii = ii - i3 * (dtset%wvl_internal%nSize(1) + 1) * &
                & (dtset%wvl_internal%nSize(2) + 1)
           i2 = ii / (dtset%wvl_internal%nSize(1) + 1)
           i0 = ii - i2 * (dtset%wvl_internal%nSize(1) + 1)
           i1 = i0 + j1 - j0
           do i = i0, i1
              iGrid = coeff_map(i , i2 , i3)
              wvl_ncoeff(iGrid) = wvl_ncoeff(iGrid) + 7
           enddo
        enddo
        basis_folder%coordinates_of_basis_grid_points%data2D => wvl_coord
        basis_folder%number_of_coefficients_per_grid_point%data1D => wvl_ncoeff
        call etsf_io_basisdata_put(wff%unwff, basis_folder, lstat, error)
        if (.not. lstat) then
           call etsf_io_low_error_to_str(errmess, error)
           write(message, "(A,A,A,A)") ch10, " wvl_write: ERROR -", ch10, &
                & errmess(1:min(475, len(errmess)))
           call wrtout(06, message, 'COLL')
        end if

        ! We write the wavefuctions.
        ! ==========================
        ! We write the wavefunctions band per band.
        allocate(wvl_band(size(wfs%psi, 1)))
        ! We reorder the wfs%psi(:,iband) into wvl_band
        do iband = 1, dtset%mband, 1
           ! iCoeff is the index of the coefficient we are writing in wvl_band
           iCoeff  = 1
           ! iCoarse is the index of the coarse part in wfs%psi (=iGrid)
           iCoarse = 1
           ! iFine is the index of the fine part in wfs%psi
           iFine   = wfs%keys%nvctr_c + 1
           do iGrid = 1, wfs%keys%nvctr_c, 1
              wvl_band(iCoeff) = wfs%psi(iCoarse, iband - mpi_enreg%me * wfs%mbandp)
              iCoeff  = iCoeff + 1
              iCoarse = iCoarse + 1
              if (wvl_ncoeff(iGrid) == 8) then
                 wvl_band(iCoeff:iCoeff + 6) = &
                      & wfs%psi(iFine:iFine + 6, iband - mpi_enreg%me * wfs%mbandp)
                 iCoeff = iCoeff + 7
                 iFine  = iFine + 7
              end if
           end do
           main_folder%wfs_coeff__state_access = iband
           main_folder%coefficients_of_wavefunctions%data1D => wvl_band
           call etsf_io_main_put(wff%unwff, main_folder, lstat, error)
           if (.not. lstat) then
              call etsf_io_low_error_to_str(errmess, error)
              write(message, "(A,A,A,A)") ch10, " wvl_write: ERROR -", ch10, &
                   & errmess(1:min(475, len(errmess)))
              call wrtout(06, message, 'COLL')
           end if
        end do

        ! We write the electronic informations.
        ! =====================================
        electrons_folder%eigenvalues%data1D => eigen
        electrons_folder%eigenvalues__number_of_states = dtset%mband
        call etsf_io_electrons_put(wff%unwff, electrons_folder, lstat, error)
        if (.not. lstat) then
           call etsf_io_low_error_to_str(errmess, error)
           write(message, "(A,A,A,A)") ch10, " wvl_write: ERROR -", ch10, &
                & errmess(1:min(475, len(errmess)))
           call wrtout(06, message, 'COLL')
        end if

        ! We deallocate all arrays
        deallocate(wvl_band)
        deallocate(wvl_coord)
        deallocate(wvl_ncoeff)
        deallocate(coeff_map)
#endif
     else
        write(message,'(a,a,a,a,a,a,i5,a)') ch10,&
             &  ' wvl_write: BUG -', ch10, &
             &  '  wff%accesswff argument is wrong,', ch10, &
             &  '  awaited values are -1, 0 (or 3 if netcdf/etsf_io is available) but value = ', wff%accesswff, '.'
        call wrtout(6,message,'PERS')
        call leave_new('PERS')
     end if
  else
     ! Write in the BigDFT way.
     call  writemywaves(mpi_enreg%me, dtset%mband, wfs%mbandp, &
          & dtset%wvl_internal%nSize(1), dtset%wvl_internal%nSize(2), &
          & dtset%wvl_internal%nSize(3), dtset%wvl_hgrid,dtset%natom, &
          & xcart(:, 1), wfs%keys%nseg_c, wfs%keys%nseg_f, wfs%keys%nvctr_c, &
          & wfs%keys%nvctr_f, wfs%keys%keyg, wfs%keys%keyv, wfs%psi, wfs%eval)
  end if

  deallocate(xcart)
#else
  write(message, '(a,a,a,a)' ) ch10,&
    &  ' wvl_outwf : BigDFT library is not compiled.', ch10, &
    &  '   Action, used the flag --enable-bigdft when configuring.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
#endif

end subroutine wvl_write
