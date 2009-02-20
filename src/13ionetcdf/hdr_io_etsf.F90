!{\src2tex{textfont=tt}}
!!****f* ABINIT/hdr_io_etsf
!! NAME
!! hdr_io_etsf
!!
!! FUNCTION
!! This subroutine deals with the I/O of the hdr_type
!! structured variables (read/write/echo).
!! It handles variables according to the ETSF format, whenever
!! possible and uses new variables when not available in the ETSF
!! format.
!! According to the value of rdwr, it reads the header
!! of a file, writes it, or echo the value of the structured
!! variable to a file.
!! Note that, when reading, different records of hdr
!! are allocated here, according to the values of the
!! read variables. Records of hdr should be deallocated
!! correctly by a call to hdr_clean when hdr is not used anymore.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  rdwr= if 1, read the hdr structured variable from the header of the netCDF file,
!!        if 2, write the header to unformatted netCDF file
!!        if 3, echo part of the header to formatted file (records 1 and 2)
!!        if 4, echo the header to formatted file
!!        if 5, read the hdr without rewinding (unformatted), identical to 1 for netCDF
!!        if 6, read the hdr without rewinding (unformatted), identical to 2 for netCDF
!!  wff=
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  The following variables are both input or output :
!!  fform=kind of the array in the file
!!   if rdwr=1,5 : will be output ; if the reading fail, return fform=0
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!  hdr <type(hdr_type)>=the header structured variable
!!   if rdwr=1,5 : will be output
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!
!! NOTES
!!
!! PARENTS
!!      cut3d,inwffil,ioarr,outkss,outwf,rdkss,testlda
!!
!! CHILDREN
!!      etsf_io_basisdata_get,etsf_io_basisdata_put,etsf_io_dims_get
!!      etsf_io_electrons_get,etsf_io_electrons_put,etsf_io_geometry_get
!!      etsf_io_geometry_put,etsf_io_kpoints_get,etsf_io_kpoints_put
!!      etsf_io_low_error_to_str,etsf_io_low_read_dim,etsf_io_low_read_var
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,hdr_io_int,leave_new
!!      rhoij_alloc,strip,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine hdr_io_etsf(fform,hdr,rdwr,unitwff)

 use defs_basis
 use defs_datatypes
#if defined HAVE_ETSF_IO
 use etsf_io_low_level
 use etsf_io
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rdwr,unitwff
 integer,intent(inout) :: fform
 type(hdr_type),intent(inout) :: hdr

!Local variables-------------------------------
!no_abirules
#if defined HAVE_ETSF_IO
 type(etsf_dims) :: dims
 type(etsf_kpoints) :: kpoints
 type(etsf_basisdata) :: basisdata
 type(etsf_geometry) :: geometry
 type(etsf_electrons) :: electrons
 type(etsf_io_low_error) :: error_data
 logical :: lstat
 character(len = etsf_io_low_error_len) :: errmess
 character(len = 500)  :: message
 integer :: cplex, ilmn, irhoij, ispden, lmn2_size, nselect, rhoijdim1, value, nresolution
 integer :: headform, iatom, itypat
 character(len=etsf_charlen), target :: basis_set
 real(dp), target :: ecut, fermie
 real(dp), target :: rprimd(3, 3)
!temp variables
 real(dp), allocatable :: rhoij(:,:,:)
#endif

! *************************************************************************

!DEBUG
!write(6,*)' hdr_io_etsf : enter hdr_io_etsf'
!call flush(6)
!ENDDEBUG

#if defined HAVE_ETSF_IO
 write(message, '(A,I0)' ) ' hdr_io_etsf: accessing ABINIT specific data from unit ', unitwff
 call wrtout(std_out, message, 'COLL')

 if(rdwr==1 .or. rdwr==5)then
! We switch off from define mode.
  call etsf_io_low_set_write_mode(unitwff, lstat, error_data = error_data)
  if (.not.lstat) goto 1000

! In case the file is just ETSF valid, we ignore the missing variables
! and we use defualt values.
  hdr%codvsn   = "ETSF  "
  fform        = 1
  headform     = 56

! First, we read the declaration of code, fform ...
! We ignore errors, assuming that the file is at least ETSF valid.
  call etsf_io_low_read_var(unitwff, "codvsn", hdr%codvsn, 6, lstat, error_data = error_data)
  if (lstat) then
!  We pad the returned string with " " instead of "\0"
   call strip(hdr%codvsn)
  end if
  call etsf_io_low_read_var(unitwff, "fform", fform, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "headform", hdr%headform, lstat, error_data = error_data)

  if(headform <= 42)then
   write(message, '(4a,i3,3a,i8,3a)' ) ch10,&
&   ' hdr_io_etsf: ERROR -',ch10,&
&   '  headform is ',headform,', while it should be > 42.',ch10,&
&   '  Action : check the correctness of your file.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! Then, read dimensions handled by ETSF
  call etsf_io_dims_get(unitwff, dims, lstat, error_data)
  if (.not.lstat) goto 1000
! Copy dimensions to hdr structure
! FIXME: don't handle k_dependent = 1
  hdr%bantot   = dims%max_number_of_states * dims%number_of_kpoints * dims%number_of_spins
  hdr%natom    = dims%number_of_atoms
  hdr%nkpt     = dims%number_of_kpoints
  hdr%nspden   = dims%number_of_components
  hdr%nspinor  = dims%number_of_spinor_components
  hdr%nsppol   = dims%number_of_spins
  hdr%nsym     = dims%number_of_symmetry_operations
  hdr%ntypat   = dims%number_of_atom_species
  hdr%ngfft(1) = dims%number_of_grid_points_vector1
  hdr%ngfft(2) = dims%number_of_grid_points_vector2
  hdr%ngfft(3) = dims%number_of_grid_points_vector3

! We read other dimensions, not handled by ETSF format.
! In case the file is just ETSF valid, we ignore the missing dimensions
! and we use default values.
  hdr%npsp    = hdr%ntypat
  rhoijdim1   = 1
  hdr%usepaw  = 0
  hdr%usewvl  = 0
  nresolution = 0
  call etsf_io_low_read_dim(unitwff, "npsp", hdr%npsp, lstat, error_data = error_data)
  call etsf_io_low_read_dim(unitwff, "rhoijdim1", rhoijdim1, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "usepaw", hdr%usepaw, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "usewvl", hdr%usewvl, lstat, error_data = error_data)
  if (hdr%usewvl == 1) then
!  This value must be 2...
   call etsf_io_low_read_dim(unitwff, "number_of_wavelet_resolutions", &
&   nresolution, lstat, error_data = error_data)
   if (.not.lstat) goto 1000

!  We set the right ngfft, adding the padding space for wavelets.
   hdr%ngfft = hdr%ngfft + 31
  end if


! Allocate all parts of hdr that need to be
  allocate(hdr%istwfk(hdr%nkpt))
  allocate(hdr%lmn_size(hdr%npsp))
  allocate(hdr%nband(hdr%nkpt*hdr%nsppol))
  allocate(hdr%npwarr(hdr%nkpt)) ! Warning : npwarr here has only one dim
  allocate(hdr%pspcod(hdr%npsp))
  allocate(hdr%pspdat(hdr%npsp))
  allocate(hdr%pspso(hdr%npsp))
  allocate(hdr%pspxc(hdr%npsp))
  allocate(hdr%so_psp(hdr%npsp))
  allocate(hdr%symafm(hdr%nsym))
  allocate(hdr%symrel(3,3,hdr%nsym))
  allocate(hdr%typat(hdr%natom))
  allocate(hdr%kptns(3,hdr%nkpt))
  allocate(hdr%occ(hdr%bantot))
  allocate(hdr%tnons(3,hdr%nsym))
  allocate(hdr%wtk(hdr%nkpt))
  allocate(hdr%xred(3,hdr%natom))
  allocate(hdr%znuclpsp(hdr%npsp))
  allocate(hdr%znucltypat(hdr%ntypat))
  allocate(hdr%zionpsp(hdr%npsp))
  allocate(hdr%title(hdr%npsp))
  if(hdr%usepaw==1) allocate(hdr%pawrhoij(hdr%natom))

! We get then all variables included in ETSF
  if (hdr%usewvl == 0) then
   basisdata%kinetic_energy_cutoff => ecut
   basisdata%number_of_coefficients => hdr%npwarr
   call etsf_io_basisdata_get(unitwff, basisdata, lstat, error_data)
   if (.not.lstat) goto 1000
  else
   call etsf_io_low_read_var(unitwff, "number_of_wavelets", hdr%nwvlarr, &
&   lstat, error_data = error_data)
   if (.not.lstat) goto 1000
  end if
  electrons%fermi_energy => fermie
  electrons%number_of_states%data1D => hdr%nband
  electrons%occupations%data1D => hdr%occ
  call etsf_io_electrons_get(unitwff, electrons, lstat, error_data)
  if (.not.lstat) goto 1000
  geometry%primitive_vectors => rprimd
  geometry%reduced_symmetry_matrices => hdr%symrel
  geometry%atom_species => hdr%typat
  geometry%reduced_symmetry_translations => hdr%tnons
  geometry%reduced_atom_positions => hdr%xred
  geometry%atomic_numbers => hdr%znucltypat
  if (hdr%npsp == hdr%ntypat) then
   geometry%valence_charges => hdr%zionpsp
  end if
  call etsf_io_geometry_get(unitwff, geometry, lstat, error_data)
  if (.not.lstat) goto 1000
  kpoints%reduced_coordinates_of_kpoints => hdr%kptns
  kpoints%kpoint_weights => hdr%wtk
  call etsf_io_kpoints_get(unitwff, kpoints, lstat, error_data)
  if (.not.lstat) goto 1000
  hdr%fermie = fermie
  hdr%ecut   = ecut
  hdr%rprimd = rprimd
  hdr%znuclpsp(1:hdr%npsp) = hdr%znucltypat(1:hdr%npsp)

! We get all other variables
! In case the file is just ETSF valid, we ignore the missing variables
! and we use default values.
  hdr%date     = 0
  hdr%ecut_eff = hdr%ecut
  hdr%ecutsm   = real(0., dp)
  hdr%etot     = real(0., dp)
  hdr%intxc    = 0
  hdr%ixc      = 1
  hdr%occopt   = 1
  hdr%pertcase = 0
  hdr%qptn(:)  = 0
  hdr%residm   = real(0., dp)
  hdr%stmbias  = real(0., dp)
  hdr%tphysel  = real(0., dp)
  hdr%tsmear   = real(0., dp)
  hdr%ecutdg   = hdr%ecut
  call etsf_io_low_read_var(unitwff, "date", hdr%date, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "ecut_eff", hdr%ecut_eff, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "ecutsm", hdr%ecutsm, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "etot", hdr%etot, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "intxc", hdr%intxc, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "ixc", hdr%ixc, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "occopt", hdr%occopt, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "pertcase", hdr%pertcase, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "qptn", hdr%qptn, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "residm", hdr%residm, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "stmbias", hdr%stmbias, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "tphysel", hdr%tphysel, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "tsmear", hdr%tsmear, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "ecutdg", hdr%ecutdg, lstat, error_data = error_data)

! test for old wavefunction style
  if(hdr%ecutsm>tol6 .and. headform<44 .and. &
&  .not.(fform==51.or.fform==52.or.fform==101.or.fform==102)) then
   write(message, '(4a,es16.6,13a)' ) ch10,&
&   ' hdr_io_etsf : ERROR -',ch10,&
&   '  The value of ecutsm is',hdr%ecutsm, &
&   ', while the file has been produced prior to v4.4 .',ch10,&
&   '  The definition of the smearing function has changed,', &
&   ' so that you are not allowed',ch10,&
&   '  to restart from a old wavefunction file. By contrast,', &
&   ' you can restart from an old',ch10,&
&   '  potential or density file, and perform a self-consistent', &
&   ' cycle with a new ABINIT version.',ch10,&
&   '  Action : produce a density or potential file using the old', &
&   ' version of ABINIT, and restart from it.'
   call wrtout(std_out, message, 'COLL')
   call leave_new('COLL')
  end if

! Multidimensional variables.
! The case of istwfk is always 1, since ETSF don't use the time reversal symetry.
  hdr%istwfk(:)   = 1
  hdr%pspcod(:)   = 0
  hdr%pspdat(:)   = 0
  hdr%pspso(:)    = 0
  hdr%pspxc(:)    = 0
  hdr%title(:)    = ""
  hdr%so_psp(:) = 1
  hdr%symafm(:)   = 1
  hdr%lmn_size    = 1
  call etsf_io_low_read_var(unitwff, "pspcod", hdr%pspcod, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "pspdat", hdr%pspdat, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "pspso", hdr%pspso, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "pspxc", hdr%pspxc, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "so_psp", hdr%so_psp, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "symafm", hdr%symafm, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "title", hdr%title, 132, lstat, error_data = error_data)
  if (lstat) then
!  We pad the returned string with " " instead of "\0"
   do itypat = 1, size(hdr%title), 1
    call strip(hdr%title(itypat))
   end do
  end if
  call etsf_io_low_read_var(unitwff, "zionpsp", hdr%zionpsp, lstat, error_data = error_data)
  call etsf_io_low_read_var(unitwff, "znuclpsp", hdr%znuclpsp, lstat, error_data = error_data)

! Compared to 4.2, add lmn_size and
  if (headform>=44) then
   call etsf_io_low_read_var(unitwff, "lmn_size", hdr%lmn_size, lstat, error_data = error_data)
   ! DC: remove the PAW specific variables, they are not in accordance
   !     with 13io_mpi/hdr_io.F90. This part should be rewritten totally.
   if (hdr%usepaw==1) then
      write(message, '(12a)' ) ch10,&
           & ' hdr_io_etsf : ERROR -',ch10,&
           & '  The support for the internal variables of PAW are not yet',ch10,&
           & '  available with ETSF output. Restarting calculation from this',ch10,&
           & '  will not be possible.',ch10,&
           & '  Action : produce a density or potential file using the old',ch10,&
           & '  binary format of ABINIT, and restart from it.'
      call wrtout(std_out, message, 'COLL')
      call leave_new('COLL')
      do iatom=1,hdr%natom
         nullify(hdr%pawrhoij(iatom)%rhoijselect)
         nullify(hdr%pawrhoij(iatom)%rhoijp)
         hdr%pawrhoij(iatom)%ngrhoij = 0
         hdr%pawrhoij(iatom)%lmnmix_sz = 0
         hdr%pawrhoij(iatom)%use_rhoij_ = 0
         hdr%pawrhoij(iatom)%use_rhoijres = 0
      end do
!!$    call etsf_io_low_read_var(unitwff, "rhoijdim1", rhoijdim1, lstat, error_data = error_data)
!!$    allocate (rhoij(rhoijdim1,hdr%nspden,hdr%natom))
!!$    call etsf_io_low_read_var(unitwff, "rhoij", rhoij, lstat, error_data = error_data)
!!$    if (.not.lstat) goto 1000
!!$
!!$    cplex=1;if (rhoijdim1/=hdr%natom) cplex=2
!!$    call rhoij_alloc(cplex,hdr%lmn_size,hdr%nspden,hdr%nsppol,hdr%pawrhoij,hdr%typat)
!!$    do iatom=1,hdr%natom
!!$     itypat=hdr%typat(iatom)
!!$     lmn2_size=hdr%lmn_size(itypat)*(hdr%lmn_size(itypat)+1)/2
!!$     nselect=0
!!$     if (cplex==1) then
!!$      do ilmn=1,lmn2_size
!!$       if (any(abs(rhoij(ilmn,:,iatom))>tol10)) then
!!$        nselect=nselect+1
!!$        hdr%pawrhoij(iatom)%rhoijselect(nselect)=ilmn
!!$        do ispden=1,hdr%nspden
!!$         hdr%pawrhoij(iatom)%rhoijp(nselect,ispden)=rhoij(ilmn,ispden,iatom)
!!$        end do
!!$       end if
!!$      end do
!!$     else
!!$      do ilmn=1,lmn2_size
!!$       if (any(abs(rhoij(2*ilmn-1:2*ilmn,:,iatom))>tol10)) then
!!$        nselect=nselect+1
!!$        hdr%pawrhoij(iatom)%rhoijselect(nselect)=ilmn
!!$        hdr%pawrhoij(iatom)%rhoijp(2*nselect-1,ispden)=rhoij(2*ilmn-1,ispden,iatom)
!!$        hdr%pawrhoij(iatom)%rhoijp(2*nselect  ,ispden)=rhoij(2*ilmn  ,ispden,iatom)
!!$       end if
!!$      end do
!!$     end if
!!$     if (nselect<lmn2_size) then
!!$      hdr%pawrhoij(iatom)%rhoijselect(nselect+1:lmn2_size)=0
!!$      do ispden=1,hdr%nspden
!!$       hdr%pawrhoij(iatom)%rhoijp(cplex*nselect+1:cplex*lmn2_size,ispden)=zero
!!$      end do
!!$     end if
!!$     hdr%pawrhoij(iatom)%nrhoijsel=nselect
!!$    end do
!!$    deallocate(rhoij)
   end if
  end if

! BigDFT private variables.
! First implementation, we assume that the number of wavelet resolutions
! is 2. Latter, we may add this value to hdr.

  lstat = .true.




! -------------------------------------------------------------------------
! Writing the header of an unformatted file
! -------------------------------------------------------------------------
 else if(rdwr==2 .or. rdwr==6)then
! We switch to write mode.
  call etsf_io_low_set_write_mode(unitwff, lstat, error_data = error_data)
  if (.not.lstat) goto 1000

! Associate and write values to ETSF groups.
  if (hdr%usewvl == 0) then
!  Plane wave case.
   ecut = hdr%ecut
   write(basis_set, "(A)") "plane_waves"
   basisdata%basis_set => basis_set
   basisdata%kinetic_energy_cutoff => ecut
   basisdata%number_of_coefficients => hdr%npwarr
  else
!  Wavelet case.
   write(basis_set, "(A)") "daubechies_wavelets"
   basisdata%basis_set => basis_set
!  Required variable than should enter the standard.
   call etsf_io_low_write_var(unitwff, "number_of_wavelets", hdr%nwvlarr, &
&   lstat, error_data = error_data)
   if (.not.lstat) goto 1000
  end if
  call etsf_io_basisdata_put(unitwff, basisdata, lstat, error_data)
  if (.not.lstat) goto 1000
  fermie = hdr%fermie
  electrons%fermi_energy => fermie
  electrons%number_of_states%data1D => hdr%nband
  electrons%occupations%data1D => hdr%occ
  call etsf_io_electrons_put(unitwff, electrons, lstat, error_data)
  if (.not.lstat) goto 1000
  rprimd = hdr%rprimd
  geometry%primitive_vectors => rprimd
  geometry%reduced_symmetry_matrices => hdr%symrel
  geometry%atom_species => hdr%typat
  geometry%reduced_symmetry_translations => hdr%tnons
  geometry%reduced_atom_positions => hdr%xred
  geometry%atomic_numbers => hdr%znucltypat
  if (hdr%npsp == hdr%ntypat) then
   geometry%valence_charges => hdr%zionpsp
  end if
  call etsf_io_geometry_put(unitwff, geometry, lstat, error_data)
  if (.not.lstat) goto 1000
  kpoints%reduced_coordinates_of_kpoints => hdr%kptns
  kpoints%kpoint_weights => hdr%wtk
  call etsf_io_kpoints_put(unitwff, kpoints, lstat, error_data)
  if (.not.lstat) goto 1000

! goto 1000
! Write none-ETSF variables.
  call etsf_io_low_write_var(unitwff, "date", hdr%date, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "codvsn", hdr%codvsn, 6, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "ecut_eff", hdr%ecut_eff, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "ecutsm", hdr%ecutsm, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "etot", hdr%etot, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "headform", 44, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "fform", fform, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "intxc", hdr%intxc, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "ixc", hdr%ixc, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "occopt", hdr%occopt, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "pertcase", hdr%pertcase, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "residm", hdr%residm, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "stmbias", hdr%stmbias, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "tphysel", hdr%tphysel, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "tsmear", hdr%tsmear, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
! Version 44 add usepaw ecutdg
  call etsf_io_low_write_var(unitwff, "ecutdg", hdr%ecutdg, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "usepaw", hdr%usepaw, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
! Array variables.
  call etsf_io_low_write_var(unitwff, "pspcod", hdr%pspcod, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "pspdat", hdr%pspdat, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "pspso", hdr%pspso, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "pspxc", hdr%pspxc, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "qptn", hdr%qptn, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "so_psp", hdr%so_psp, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "symafm", hdr%symafm, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "title", hdr%title, 132, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_write_var(unitwff, "znuclpsp", hdr%znuclpsp, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  if (hdr%npsp /= hdr%ntypat) then
   call etsf_io_low_write_var(unitwff, "zionpsp", hdr%zionpsp, lstat, error_data = error_data)
   if (.not.lstat) goto 1000
  end if
! Version 44 add lmn_size and rhoij
  call etsf_io_low_write_var(unitwff, "lmn_size", hdr%lmn_size, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
   ! DC: remove the PAW specific variables, they are not in accordance
   !     with 13io_mpi/hdr_io.F90. This part should be rewritten totally.
  if (hdr%usepaw == 1) then
     write(message, '(12a)' ) ch10,&
          & ' hdr_io_etsf : WARNING -',ch10,&
          & '  The support for the internal variables of PAW are not yet',ch10,&
          & '  available with ETSF output. Restarting calculation from this',ch10,&
          & '  will not be possible.',ch10,&
          & '  Action : produce a density or potential file using the old',ch10,&
          & '  binary format of ABINIT, and restart from it.'
     call wrtout(std_out, message, 'COLL')
!!$   rhoijdim1 = maxval(hdr%lmn_size)
!!$   rhoijdim1 = hdr%pawrhoij(1)%cplex*rhoijdim1*(rhoijdim1+1)/2
!!$   call etsf_io_low_write_var(unitwff, "rhoijdim1", rhoijdim1, lstat, error_data = error_data)
!!$   allocate (rhoij(rhoijdim1,hdr%nspden,hdr%natom))
!!$   do iatom=1,hdr%natom
!!$    itypat=hdr%typat(iatom)
!!$    lmn2_size = hdr%lmn_size(itypat)*(hdr%lmn_size(itypat)+1)/2
!!$    cplex=hdr%pawrhoij(iatom)%cplex
!!$    do ispden=1,hdr%nspden
!!$     rhoij(1:cplex*lmn2_size,ispden,iatom)=zero
!!$     if (cplex==1) then
!!$      do irhoij=1,hdr%pawrhoij(iatom)%nrhoijsel
!!$       ilmn=hdr%pawrhoij(iatom)%rhoijselect(irhoij)
!!$       rhoij(ilmn,ispden,iatom)=hdr%pawrhoij(iatom)%rhoijp(irhoij,ispden)
!!$      end do
!!$     else
!!$      do irhoij=1,hdr%pawrhoij(iatom)%nrhoijsel
!!$       ilmn=hdr%pawrhoij(iatom)%rhoijselect(irhoij)
!!$       rhoij(2*ilmn-1,ispden,iatom)=hdr%pawrhoij(iatom)%rhoijp(2*irhoij-1,ispden)
!!$       rhoij(2*ilmn  ,ispden,iatom)=hdr%pawrhoij(iatom)%rhoijp(2*irhoij  ,ispden)
!!$      end do
!!$     end if
!!$    end do
!!$   end do
!!$   call etsf_io_low_write_var(unitwff, "rhoij", rhoij, lstat, error_data = error_data)
!!$   if (.not.lstat) goto 1000
!!$   deallocate (rhoij)
  end if
! BigDFT variables.
  call etsf_io_low_write_var(unitwff, "usewvl", hdr%usewvl, lstat, error_data = error_data)
  if (.not.lstat) goto 1000

 else if(rdwr==3 .or. rdwr==4)then

  call hdr_io_int(fform, hdr, rdwr, unitwff)
  lstat = .true.

 end if ! choice read/write/echo

 1000 continue
 if (.not. lstat) then
  call etsf_io_low_error_to_str(errmess, error_data)
  write(message, "(A,A,A,A)") ch10, " hdr_io_etsf: ERROR -", ch10, &
&  errmess(1:min(475, len(errmess)))
  call wrtout(std_out, message, 'COLL')
  call leave_new('COLL')
 end if

#endif


 return
 fform=0 ; return   ! This is to allow treatment of old epsm1 format

end subroutine hdr_io_etsf
!!***
