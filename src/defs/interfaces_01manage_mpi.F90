!!****m* ABINIT/interfaces_01manage_mpi
!! NAME
!! interfaces_01manage_mpi
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/01manage_mpi
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!! 
!!
!! SOURCE

module interfaces_01manage_mpi

 implicit none

interface
 subroutine clnmpi_band(nkpt,nsppol,mpi_enreg)
  use defs_datatypes
  implicit none
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(mpi_type), intent(inout) :: mpi_enreg
 end subroutine clnmpi_band
end interface

interface
 subroutine clnmpi_fft(dtset, mpi_enreg)
  use defs_datatypes
  implicit none
  type(dataset_type), intent(in) :: dtset
  type(mpi_type), intent(inout) :: mpi_enreg
 end subroutine clnmpi_fft
end interface

interface
 subroutine clnmpi_gs(dtset, mpi_enreg)
  use defs_datatypes
  implicit none
  type(dataset_type) :: dtset
  type(mpi_type) :: mpi_enreg
 end subroutine clnmpi_gs
end interface

interface
 subroutine distrb2(mband, nband, nkpt, nsppol, mpi_enreg)
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine distrb2
end interface

interface
 subroutine herald(code_name,code_version,iout)
  implicit none
  integer,intent(in) :: iout
  character(len=24),intent(in) :: code_name
  character(len=6),intent(in) :: code_version
 end subroutine herald
end interface

interface
 subroutine initmpi_band(mpi_enreg,nband,nkpt,nsppol)
  use defs_datatypes
  implicit none
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine initmpi_band
end interface

interface
 subroutine initmpi_fft(dtset,mpi_enreg)
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine initmpi_fft
end interface

interface
 subroutine initmpi_grid(dtset,mpi_enreg)
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine initmpi_grid
end interface

interface
 subroutine initmpi_gs(dtset,mpi_enreg)
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine initmpi_gs
end interface

interface
 subroutine initmpi_respfn(mpi_enreg, spaceComm)
  use defs_datatypes
  implicit none
  integer :: spaceComm
  type(mpi_type) :: mpi_enreg
 end subroutine initmpi_respfn
end interface

interface
 subroutine initmpi_seq(mpi_enreg)
  use defs_datatypes
  implicit none
  type(mpi_type),intent(out) :: mpi_enreg
 end subroutine initmpi_seq
end interface

interface
 subroutine leave_new(mode_paral)
  implicit none
  character(len=4),intent(in) :: mode_paral
 end subroutine leave_new
end interface

interface
 subroutine leave_test(mpi_enreg)
  use defs_datatypes
  implicit none
  type(mpi_type) :: mpi_enreg
 end subroutine leave_test
end interface

interface
 subroutine pre_scatter(array,array_allgather,n1,n2,n3,mpi_enreg,option)
  use defs_basis
  use defs_datatypes
  implicit none
  integer :: n1
  integer :: n2
  integer :: n3
  type(mpi_type) :: mpi_enreg
  character(*) :: option
  real(dp) :: array(n1,n2,n3/mpi_enreg%nproc_fft,1)
  real(dp) :: array_allgather(n1,n2,n3,1)
 end subroutine pre_scatter
end interface

interface
 SUBROUTINE build_grid_scalapack(grid,nbprocs, communicator)
  use defs_scalapack
  implicit none
  integer, intent(in) :: communicator
  integer,intent(in) :: nbprocs
  type(grid_scalapack),intent(out) :: grid
 end subroutine build_grid_scalapack
end interface

interface
 SUBROUTINE build_processor_scalapack(processor,grid,myproc, comm)
  use defs_scalapack
  implicit none
  integer,intent(in) :: comm
  integer,intent(in) :: myproc
  type(grid_scalapack),intent(in) :: grid
  type(processor_scalapack),intent(out) :: processor
 end subroutine build_processor_scalapack
end interface

interface
 SUBROUTINE init_scalapack(processor,communicator)
  use defs_scalapack
  implicit none
  integer, intent(in) :: communicator
  type(processor_scalapack),intent(out) :: processor
 end subroutine init_scalapack
end interface

interface
 SUBROUTINE end_scalapack(processor)
  use defs_scalapack
  implicit none
  type(processor_scalapack),intent(inout) :: processor
 end subroutine end_scalapack
end interface

interface
 SUBROUTINE init_matrix_scalapack(matrix,nbli_global,&  
  &  nbco_global,processor,tbloc)
  use defs_scalapack
  implicit none
  integer,intent(in) :: nbco_global
  integer,intent(in) :: nbli_global
  integer,intent(in),optional :: tbloc
  type(matrix_scalapack),intent(out) :: matrix
  type(processor_scalapack),intent(in),target :: processor
 end subroutine init_matrix_scalapack
end interface

interface
 SUBROUTINE matrix_copy(source,destination)
  use defs_scalapack
  implicit none
  type(matrix_scalapack),intent(out) :: destination
  type(matrix_scalapack),intent(in) :: source
 end subroutine matrix_copy
end interface

interface
 SUBROUTINE destruction_matrix_scalapack(matrix)
  use defs_scalapack
  implicit none
  type(matrix_scalapack),intent(inout) :: matrix
 end subroutine destruction_matrix_scalapack
end interface

interface
 FUNCTION matrix_get_local(matrix,i,j)
  use defs_basis
  use defs_scalapack
  implicit none
  integer, intent(in) :: i
  integer, intent(in) :: j
  type(matrix_scalapack),intent(in) :: matrix
  complex(dp) :: matrix_get_local
 end function matrix_get_local
end interface

interface
 SUBROUTINE matrix_set_local(matrix,i,j,value)
  use defs_basis
  use defs_scalapack
  implicit none
  integer, intent(in) :: i
  integer, intent(in) :: j
  type(matrix_scalapack),intent(out) :: matrix
  complex(dp) :: value
 end subroutine matrix_set_local
end interface

interface
 SUBROUTINE matrix_add_local(matrix,i,j,value)
  use defs_basis
  use defs_scalapack
  implicit none
  integer, intent(in) :: i
  integer, intent(in) :: j
  type(matrix_scalapack),intent(out) :: matrix
  complex(dp) :: value
 end subroutine matrix_add_local
end interface

interface
 FUNCTION idx_processor_is_local(matrix,i,j)
  use defs_scalapack
  implicit none
  integer, intent(in) :: i
  integer, intent(in) :: j
  logical :: idx_processor_is_local
  type(matrix_scalapack),intent(in) :: matrix
 end function idx_processor_is_local
end interface

interface
 FUNCTION idx_processor_concerned(matrix,idx,lico)
  use defs_scalapack
  implicit none
  integer, intent(in) :: idx
  integer :: idx_processor_concerned
  integer, intent(in) :: lico
  type(matrix_scalapack),intent(in) :: matrix
 end function idx_processor_concerned
end interface

interface
 SUBROUTINE idx_loc(matrix,i,j,iloc,jloc)
  use defs_scalapack
  implicit none
  integer, intent(in) :: i
  integer, intent(out) :: iloc
  integer, intent(in) :: j
  integer, intent(out) :: jloc
  type(matrix_scalapack),intent(in) :: matrix
 end subroutine idx_loc
end interface

interface
 FUNCTION glob_loc(matrix,idx,lico)
  use defs_scalapack
  implicit none
  integer :: glob_loc
  integer, intent(in) :: idx
  integer, intent(in) :: lico
  type(matrix_scalapack),intent(in) :: matrix
 end function glob_loc
end interface

interface
 SUBROUTINE idx_glob(matrix,iloc,jloc,i,j)
  use defs_scalapack
  implicit none
  integer, intent(out) :: i
  integer, intent(in) :: iloc
  integer, intent(out) :: j
  integer, intent(in) :: jloc
  type(matrix_scalapack),intent(in) :: matrix
 end subroutine idx_glob
end interface

interface
 FUNCTION loc_glob(matrix,proc,idx,lico)
  use defs_scalapack
  implicit none
  integer, intent(in) :: idx
  integer, intent(in) :: lico
  integer :: loc_glob
  type(matrix_scalapack),intent(in) :: matrix
  type(processor_scalapack),intent(in) :: proc
 end function loc_glob
end interface

interface
 FUNCTION matrix_get_global(matrix,i,j)
  use defs_basis
  use defs_scalapack
  implicit none
  integer, intent(in) :: i
  integer, intent(in) :: j
  type(matrix_scalapack),intent(in) :: matrix
  complex(dp) :: matrix_get_global
 end function matrix_get_global
end interface

interface
 SUBROUTINE matrix_set_global(matrix,i,j,value)
  use defs_basis
  use defs_scalapack
  implicit none
  integer, intent(in) :: i
  integer, intent(in) :: j
  type(matrix_scalapack),intent(inout) :: matrix
  complex(dp) :: value
 end subroutine matrix_set_global
end interface

interface
 SUBROUTINE matrix_add_global(matrix,i,j,value)
  use defs_basis
  use defs_scalapack
  implicit none
  integer, intent(in) :: i
  integer, intent(in) :: j
  type(matrix_scalapack),intent(inout) :: matrix
  complex(dp) :: value
 end subroutine matrix_add_global
end interface

interface
 SUBROUTINE matrix_check_global(matrix,reference)
  use defs_basis
  use defs_scalapack
  implicit none
  type(matrix_scalapack),intent(in) :: matrix
  complex(dp),dimension(:,:) :: reference
 end subroutine matrix_check_global
end interface

interface
 SUBROUTINE matrix_check_global_vector(matrix,reference,decli,nom)
  use defs_basis
  use defs_scalapack
  implicit none
  integer, optional :: decli
  type(matrix_scalapack),intent(in) :: matrix
  character(len=*),optional :: nom
  complex(dp),dimension(:) :: reference
 end subroutine matrix_check_global_vector
end interface

interface
 SUBROUTINE matrix_from_global(matrix,reference)
  use defs_basis
  use defs_scalapack
  implicit none
  type(matrix_scalapack),intent(inout) :: matrix
  real(dp),dimension(:) :: reference
 end subroutine matrix_from_global
end interface

interface
 SUBROUTINE matrix_to_reference(matrix,reference)
  use defs_basis
  use defs_scalapack
  implicit none
  type(matrix_scalapack),intent(in) :: matrix
  real(dp),dimension(:,:),intent(inout) :: reference
 end subroutine matrix_to_reference
end interface

interface
 SUBROUTINE matrix_storage_local_vector(matrix,array,decli,decco)
  use defs_basis
  use defs_scalapack
  implicit none
  integer, intent(in) :: decco
  integer, intent(in) :: decli
  type(matrix_scalapack),intent(inout) :: matrix
  complex(dp),dimension(:),intent(in) :: array
 end subroutine matrix_storage_local_vector
end interface

interface
 SUBROUTINE matrix_extract_vector(matrix,array,decli,decco,finli)
  use defs_basis
  use defs_scalapack
  implicit none
  integer, intent(in) :: decco
  integer, intent(in) :: decli
  integer, intent(in), optional :: finli
  type(matrix_scalapack),intent(in) :: matrix
  complex(dp),dimension(:),intent(out) :: array
 end subroutine matrix_extract_vector
end interface

interface
 SUBROUTINE matrix_pzgetrf(matrix)
  use defs_scalapack
  implicit none
  type(matrix_scalapack),intent(inout) :: matrix
 end subroutine matrix_pzgetrf
end interface

interface
 SUBROUTINE matrix_pzgetrs(matrix,vector)
  use defs_scalapack
  implicit none
  type(matrix_scalapack),intent(in) :: matrix
  type(matrix_scalapack),intent(inout) :: vector
 end subroutine matrix_pzgetrs
end interface

interface
 SUBROUTINE matrix_pzpotrf(matrix)
  use defs_scalapack
  implicit none
  type(matrix_scalapack),intent(inout) :: matrix
 end subroutine matrix_pzpotrf
end interface

interface
 SUBROUTINE matrix_pzpotrs(matrix,vector)
  use defs_scalapack
  implicit none
  type(matrix_scalapack),intent(in) :: matrix
  type(matrix_scalapack),intent(inout) :: vector
 end subroutine matrix_pzpotrs
end interface

interface
 SUBROUTINE matrix_pzgemm(matrix1,alpha,matrix2,beta,results)
  use defs_basis
  use defs_scalapack
  implicit none
  complex(dp), intent(in) :: alpha
  complex(dp), intent(in) :: beta
  type(matrix_scalapack),intent(in) :: matrix1
  type(matrix_scalapack),intent(in) :: matrix2
  type(matrix_scalapack),intent(inout) :: results
 end subroutine matrix_pzgemm
end interface

interface
 SUBROUTINE matrix_pzheevx(processor,matrix,results,eigen,communicator)
  use defs_scalapack
  implicit none
  integer,intent(in) :: communicator
  type(matrix_scalapack),intent(in) :: matrix
  type(processor_scalapack),intent(in) :: processor
  type(matrix_scalapack),intent(inout) :: results
  double precision,dimension(:),intent(inout) :: eigen
 end subroutine matrix_pzheevx
end interface

interface
 SUBROUTINE matrix_pzhegvx(processor,matrix1,matrix2,results,eigen,communicator)
  use defs_scalapack
  implicit none
  integer,intent(in) :: communicator
  type(matrix_scalapack),intent(in) :: matrix1
  type(matrix_scalapack),intent(in) :: matrix2
  type(processor_scalapack),intent(in) :: processor
  type(matrix_scalapack),intent(inout) :: results
  double precision,dimension(:),intent(inout) :: eigen
 end subroutine matrix_pzhegvx
end interface

interface
 SUBROUTINE NO_SCALAPACK
  implicit none
 end subroutine NO_SCALAPACK
end interface

interface
 subroutine split_work(ntasks,istart,istop,verbose)
  implicit none
  integer,intent(inout) :: istart
  integer,intent(inout) :: istop
  integer,intent(in) :: ntasks
  integer,optional,intent(in) :: verbose
 end subroutine split_work
end interface

interface
 subroutine split_work2(ntasks,nprocs,istart,istop,verbose)
  implicit none
  integer,intent(in) :: nprocs
  integer,intent(in) :: ntasks
  integer,optional,intent(in) :: verbose
  integer,intent(inout) :: istart(nprocs)
  integer,intent(inout) :: istop(nprocs)
 end subroutine split_work2
end interface

interface
 subroutine wrtout(unit,message,mode_paral)
  implicit none
  integer,intent(in) :: unit
  character(len=500),intent(inout) :: message
  character(len=4),intent(in) :: mode_paral
 end subroutine wrtout
end interface

end module interfaces_01manage_mpi
!!***
