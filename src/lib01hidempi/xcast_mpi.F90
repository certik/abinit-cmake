!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcast_mpi
!! NAME
!! xcast_mpi
!!
!! FUNCTION
!! This module contains functions that calls MPI routine,
!! if we compile the code using the MPI CPP flags.
!! xcast_mpi is the generic function.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (Rshaltaf,AR,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!--------------------------------------------------------------------

subroutine xcast_mpi_intv(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
  integer,intent(inout) :: xval
  integer,intent(in) :: spaceComm,master
  integer,intent(out)   :: ier
  ier=0
#if defined MPI
   call MPI_BCAST(xval,1,MPI_INTEGER,master,spaceComm,ier)
#endif
  end subroutine xcast_mpi_intv

!---------------------------------------------------------------

subroutine xcast_mpi_int1d(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(inout) :: xval(:)
 integer,intent(in) :: spaceComm,master
 integer,intent(out)   :: ier
!local
 integer :: n
 ier=0
 n=size(xval)
#if defined MPI 
                   call MPI_BCAST(xval,n,MPI_INTEGER,master,spaceComm,ier)
#endif
end subroutine xcast_mpi_int1d

!--------------------------------------------------------------------

subroutine xcast_mpi_int2d(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(inout) :: xval(:,:)
 integer,intent(in) :: spaceComm,master
 integer,intent(out)   :: ier

!Local variables-------------------
 integer :: n1,n2
 n1=size(xval,dim=1)
 n2=size(xval,dim=2)
 ier=0
#if defined MPI
                  call MPI_BCAST(xval,n1*n2,MPI_INTEGER,master,spaceComm,ier)
#endif
end subroutine xcast_mpi_int2d

!--------------------------------------------------------------------

subroutine xcast_mpi_int3d(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
  include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(inout) :: xval(:,:,:)
 integer,intent(in) :: spaceComm,master
 integer,intent(out)   :: ier

!Local variables-------------------
 integer :: n1,n2,n3
 n1=size(xval,dim=1)
 n2=size(xval,dim=2)
 n3=size(xval,dim=3)
 ier=0
#if defined MPI
                  call MPI_BCAST(xval,n1*n2*n3,MPI_INTEGER,master,spaceComm,ier)
#endif
end subroutine xcast_mpi_int3d

!-------------------------------------------------------------------------------------

subroutine xcast_mpi_dpv(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

  implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier
!Local variables-------------------
 ier=0
#if defined MPI
     call MPI_BCAST(xval,1,MPI_DOUBLE_PRECISION,master,spaceComm,ier)
#endif

  end subroutine xcast_mpi_dpv
!-------------------------------------------------------------------------------------------------
subroutine xcast_mpi_dp1d(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
    include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer::n
 n=size(xval,dim=1)
 ier=0
#if defined MPI
  call MPI_BCAST(xval,n,MPI_DOUBLE_PRECISION,master,spaceComm,ier)
#endif

end subroutine xcast_mpi_dp1d

!-------------------------------------------------------------------------------------------

subroutine xcast_mpi_dp2d(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
    include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer::n1,n2
 n1=size(xval,dim=1)
 n2=size(xval,dim=2) 
 ier=0
#if defined MPI
  call MPI_BCAST(xval,n1*n2,MPI_DOUBLE_PRECISION,master,spaceComm,ier)
#endif

end subroutine xcast_mpi_dp2d

!-------------------------------------------------------------------------------------------

subroutine xcast_mpi_dp3d(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
    include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer::n1,n2,n3
 n1=size(xval,dim=1)
 n2=size(xval,dim=2)
 n3=size(xval,dim=3)
 ier=0
#if defined MPI
  call MPI_BCAST(xval,n1*n2*n3,MPI_DOUBLE_PRECISION,master,spaceComm,ier)
#endif

end subroutine xcast_mpi_dp3d

!------------------------------------------------------------------------------------------

subroutine xcast_mpi_dp4d(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
    include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:,:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer::n1,n2,n3,n4
 n1=size(xval,dim=1)
 n2=size(xval,dim=2)
 n3=size(xval,dim=3)
 n4=size(xval,dim=4)
 ier=0
#if defined MPI
  call MPI_BCAST(xval,n1*n2*n3*n4,MPI_DOUBLE_PRECISION,master,spaceComm,ier)
#endif

end subroutine xcast_mpi_dp4d

!------------------------------------------------------------------------------------------

subroutine xcast_mpi_spv(xval,master,spaceComm,ier)

use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

implicit none

#if defined MPI && defined MPI1
  include 'mpif.h'
#endif

!Arguments-------------------------
 real,intent(inout) :: xval
 integer,intent(in) :: spaceComm,master
 integer,intent(out) :: ier

!Local variables-------------------
 
 ier=0
#if defined MPI 
 call MPI_BCAST(xval,1,MPI_REAL,master,spaceComm,ier)
#endif

end subroutine xcast_mpi_spv

!--------------------------------------------------------------------------------------------
subroutine xcast_mpi_sp1d(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
    include 'mpif.h'
#endif

!Arguments-------------------------
 real,intent(inout) :: xval(:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer::n
 n=size(xval,dim=1)
 ier=0
#if defined MPI 
  call MPI_BCAST(xval,n,MPI_REAL,master,spaceComm,ier)
#endif

end subroutine xcast_mpi_sp1d

!-------------------------------------------------------------------------------------------

subroutine xcast_mpi_sp2d(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real,intent(inout) :: xval(:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer::n1,n2
 n1=size(xval,dim=1)
 n2=size(xval,dim=2)
 ier=0
#if defined MPI
  call MPI_BCAST(xval,n1*n2,MPI_REAL,master,spaceComm,ier)
#endif

end subroutine xcast_mpi_sp2d
!--------------------------------------------------------------
subroutine xcast_mpi_sp3d(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
    include 'mpif.h'
#endif

!Arguments-------------------------
 real,intent(inout) :: xval(:,:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer::n1,n2,n3
 n1=size(xval,dim=1)
 n2=size(xval,dim=2)
 n3=size(xval,dim=3)
 ier=0
#if defined MPI
  call MPI_BCAST(xval,n1*n2*n3,MPI_REAL,master,spaceComm,ier)
#endif

end subroutine xcast_mpi_sp3d

!------------------------------------------------------------------------------------------

subroutine xcast_mpi_sp4d(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
    include 'mpif.h'
#endif

!Arguments-------------------------
 real,intent(inout) :: xval(:,:,:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer::n1,n2,n3,n4
 n1=size(xval,dim=1)
 n2=size(xval,dim=2)
 n3=size(xval,dim=3)
 n4=size(xval,dim=4)
 ier=0
#if defined MPI
  call MPI_BCAST(xval,n1*n2*n3*n4,MPI_REAL,master,spaceComm,ier)
#endif

end subroutine xcast_mpi_sp4d

!-----------------------------------------------------------------------------------------------

subroutine xcast_mpi_cplxv(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex,intent(inout) :: xval
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
 ier=0
#if defined MPI 
  call MPI_BCAST(xval,1,MPI_COMPLEX,master,spaceComm,ier)
#endif

 end subroutine xcast_mpi_cplxv

!-----------------------------------------------------------------------------------------------------

 subroutine xcast_mpi_cplx1d(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex,intent(inout) :: xval(:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer::n
 n=size(xval(:))
 ier=0
#if defined MPI 
  call MPI_BCAST(xval,n,MPI_COMPLEX,master,spaceComm,ier)
#endif

  end subroutine xcast_mpi_cplx1d
!----------------------------------------------------------------------------------------------------

subroutine xcast_mpi_cplx2d(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex,intent(inout) :: xval(:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
integer::n1,n2
n1=size(xval,dim=1)
n2=size(xval,dim=2)
 ier=0
#if defined MPI 
  call MPI_BCAST(xval,n1*n2,MPI_COMPLEX,master,spaceComm,ier)
#endif

 end subroutine xcast_mpi_cplx2d
!--------------------------------------------------------------------------------------------------
subroutine xcast_mpi_cplx3d(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex,intent(inout) :: xval(:,:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
integer::n1,n2,n3
n1=size(xval,dim=1)
n2=size(xval,dim=2)
n3=size(xval,dim=3)
 ier=0
#if defined MPI 
  call MPI_BCAST(xval,n1*n2*n3,MPI_COMPLEX,master,spaceComm,ier)
#endif

 end subroutine xcast_mpi_cplx3d




!--------------------------------------------------------------------------------------------------
subroutine xcast_mpi_cplx4d(xval,master,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex,intent(inout) :: xval(:,:,:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
integer::n1,n2,n3,n4
n1=size(xval,dim=1)
n2=size(xval,dim=2)
n3=size(xval,dim=3)
n4=size(xval,dim=4)
 ier=0
#if defined MPI
  call MPI_BCAST(xval,n1*n2*n3*n4,MPI_COMPLEX,master,spaceComm,ier)
#endif

 end subroutine xcast_mpi_cplx4d


!---------------------------------------------------------------------

subroutine xcast_mpi_dcv(xval,master,spaceComm,ier)

use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dp),intent(inout):: xval
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
  ier=0
#if defined MPI
  call MPI_BCAST(xval,1,MPI_DOUBLE_COMPLEX,master,spaceComm,ier)
#endif

end subroutine xcast_mpi_dcv
!-----------------------------------------------------------------------

subroutine xcast_mpi_dc1d(xval,master,spaceComm,ier)

use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dp),intent(inout):: xval(:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer::n
 n=size(xval(:))
  ier=0
#if defined MPI
  call MPI_BCAST(xval,n,MPI_DOUBLE_COMPLEX,master,spaceComm,ier)
#endif

end subroutine xcast_mpi_dc1d
!--------------------------------------------------------------
subroutine xcast_mpi_dc2d(xval,master,spaceComm,ier)

use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dp),intent(inout):: xval(:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
  integer::n1,n2
  n1=size(xval,dim=1)
  n2=size(xval,dim=2)
  ier=0
#if defined MPI
  call MPI_BCAST(xval,n1*n2,MPI_DOUBLE_COMPLEX,master,spaceComm,ier)
#endif

end subroutine xcast_mpi_dc2d
!-------------------------------------------------------------------
subroutine xcast_mpi_dc3d(xval,master,spaceComm,ier)

use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dp),intent(inout):: xval(:,:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out)   :: ier

!Local variables-------------------
  integer::n1,n2,n3
  n1=size(xval,dim=1)
  n2=size(xval,dim=2)
  n3=size(xval,dim=3)
  ier=0
#if defined MPI
  call MPI_BCAST(xval,n1*n2*n3,MPI_DOUBLE_COMPLEX,master,spaceComm,ier)
#endif

end subroutine xcast_mpi_dc3d
!!***
