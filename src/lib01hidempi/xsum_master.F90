!{\src2tex{textfont=tt}}
!!****f* ABINIT/xsum_master.F90
!! NAME
!! xsum_master.F90
!!
!! FUNCTION
!! This module contains functions that calls MPI routine,
!! if we compile the code using the MPI  CPP flags.
!! xsum_master is the generic function.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (AR,XG,MB)
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

subroutine xsum_master_dp1d(xval,master,spaceComm,ier)
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
 integer ,intent(in) :: master
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer :: n1

#if defined MPI
           real(dp) , allocatable :: xsum(:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
           if (nproc_space_comm /= 1) then
            n1 = size(xval,dim=1)
            !Accumulate xval on all proc. in spaceComm
            allocate(xsum(n1))
            call MPI_REDUCE(xval,xsum,n1,MPI_DOUBLE_PRECISION,&
            &  MPI_SUM,master,spaceComm,ier)
            xval (:) = xsum(:)
            deallocate(xsum)
           end if
#endif
end subroutine xsum_master_dp1d

!--------------------------------------------------------------------
subroutine xsum_master_dp2d(xval,master,spaceComm,ier)
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
 integer ,intent(in) :: master
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer :: n1,n2

#if defined MPI
           real(dp) , allocatable :: xsum(:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
           if (nproc_space_comm /= 1) then
            n1 = size(xval,dim=1)
            n2 = size(xval,dim=2)
            !Accumulate xval on all proc. in spaceComm
            allocate(xsum(n1,n2))
            call MPI_REDUCE(xval,xsum,n1*n2,MPI_DOUBLE_PRECISION,&
            &  MPI_SUM,master,spaceComm,ier)
            xval (:,:) = xsum(:,:)
            deallocate(xsum)
           end if
#endif
end subroutine xsum_master_dp2d

!--------------------------------------------------------------------

subroutine xsum_master_dp3d(xval,master,spaceComm,ier)
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
 integer ,intent(in) :: master
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer :: n1,n2,n3

#if defined MPI
           real(dp) , allocatable :: xsum(:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
           if (nproc_space_comm /= 1) then
            n1 = size(xval,dim=1)
            n2 = size(xval,dim=2)
            n3 = size(xval,dim=3)
            !Accumulate xval on all proc. in spaceComm
            allocate(xsum(n1,n2,n3))
            call MPI_REDUCE(xval,xsum,n1*n2*n3,MPI_DOUBLE_PRECISION,&
            &  MPI_SUM,master,spaceComm,ier)
            xval (:,:,:) = xsum(:,:,:)
            deallocate(xsum)
           end if
#endif
end subroutine xsum_master_dp3d

!--------------------------------------------------------------------

subroutine xsum_master_dp4d(xval,master,spaceComm,ier)
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
 integer ,intent(in) :: master
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer :: n1,n2,n3,n4

#if defined MPI
           real(dp) , allocatable :: xsum(:,:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
           if (nproc_space_comm /= 1) then
            n1 = size(xval,dim=1)
            n2 = size(xval,dim=2)
            n3 = size(xval,dim=3)
            n4 = size(xval,dim=4)
            !Accumulate xval on all proc. in spaceComm
            allocate(xsum(n1,n2,n3,n4))
            call MPI_REDUCE(xval,xsum,n1*n2*n3*n4,MPI_DOUBLE_PRECISION,&
            &  MPI_SUM,master,spaceComm,ier)
            xval (:,:,:,:) = xsum(:,:,:,:)
            deallocate(xsum)
           end if
#endif

end subroutine xsum_master_dp4d

!--------------------------------------------------------------------
subroutine xsum_master_dp5d(xval,master,spaceComm,ier)
 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif


 real(dp),intent(inout) :: xval(:,:,:,:,:)
 integer ,intent(in) :: master
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

#if defined MPI
           integer :: n1,n2,n3,n4,n5
           real(dp), allocatable :: xsum(:,:,:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
           if (nproc_space_comm /= 1) then
            n1 = size(xval,dim=1)
            n2 = size(xval,dim=2)
            n3 = size(xval,dim=3)
            n4 = size(xval,dim=4)
            n5 = size(xval,dim=5)
 !Accumulate xval on all proc. in spaceComm
            allocate(xsum(n1,n2,n3,n4,n5))
            call MPI_reduce(xval,xsum,n1*n2*n3*n4*n5,MPI_DOUBLE_PRECISION,&
            &  MPI_SUM,master,spaceComm,ier)
            xval (:,:,:,:,:) = xsum(:,:,:,:,:)
            deallocate(xsum)
           end if
#endif

end subroutine xsum_master_dp5d

!--------------------------------------------------------------------
subroutine xsum_master_dp6d(xval,master,spaceComm,ier)
 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif


 real(dp),intent(inout) :: xval(:,:,:,:,:,:)
 integer ,intent(in) :: master
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

#if defined MPI
           integer :: n1,n2,n3,n4,n5,n6
           real(dp), allocatable :: xsum(:,:,:,:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
           if (nproc_space_comm /= 1) then
            n1 = size(xval,dim=1)
            n2 = size(xval,dim=2)
            n3 = size(xval,dim=3)
            n4 = size(xval,dim=4)
            n5 = size(xval,dim=5)
            n6 = size(xval,dim=6)
 !Accumulate xval on all proc. in spaceComm
            allocate(xsum(n1,n2,n3,n4,n5,n6))
            call MPI_reduce(xval,xsum,n1*n2*n3*n4*n5*n6,MPI_DOUBLE_PRECISION,&
            &  MPI_SUM,master,spaceComm,ier)
            xval (:,:,:,:,:,:) = xsum(:,:,:,:,:,:)
            deallocate(xsum)
           end if
#endif

end subroutine xsum_master_dp6d

!--------------------------------------------------------------------

subroutine xsum_master_int4d(xval,master,spaceComm,ier)
 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

 integer ,intent(inout) :: xval(:,:,:,:)
 integer ,intent(in) :: master
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

#if defined MPI
           integer :: n1,n2,n3,n4
           integer, allocatable :: xsum(:,:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
           if (nproc_space_comm /= 1) then
            n1 = size(xval,dim=1)
            n2 = size(xval,dim=2)
            n3 = size(xval,dim=3)
            n4 = size(xval,dim=4)
!           Accumulate xval on all proc. in spaceComm
            allocate(xsum(n1,n2,n3,n4))
            call MPI_reduce(xval,xsum,n1*n2*n3*n4,MPI_INTEGER,&
            &  MPI_SUM,master,spaceComm,ier)
            xval (:,:,:,:) = xsum(:,:,:,:)
            deallocate(xsum)
           end if
#endif

end subroutine xsum_master_int4d

!--------------------------------------------------------------------

subroutine xsum_master_c2cplx(xval,master,spaceComm,ier)
 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex ,intent(inout) :: xval(:,:)
 integer ,intent(in) :: master
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer :: n1,n2

#if defined MPI
           complex , allocatable :: xsum(:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
           if (nproc_space_comm /= 1) then
            n1 = size(xval,dim=1)
            n2 = size(xval,dim=2)
            ! collect xval from processors on master in spaceComm
            allocate(xsum(n1,n2))
            call MPI_REDUCE(xval,xsum,n1*n2,MPI_COMPLEX,&
            &  MPI_SUM,master,spaceComm,ier)
            xval (:,:) = xsum(:,:)
            deallocate(xsum)
           end if
#endif
end subroutine xsum_master_c2cplx

!--------------------------------------------------------------------

subroutine xsum_master_c3cplx(xval,master,spaceComm,ier)
 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex ,intent(inout) :: xval(:,:,:)
 integer ,intent(in) :: master
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer :: n1,n2,n3

#if defined MPI
           complex , allocatable :: xsum(:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
           if (nproc_space_comm /= 1) then
            n1 = size(xval,dim=1)
            n2 = size(xval,dim=2)
            n3 = size(xval,dim=3)
            ! collect xval from processors on master in spaceComm
            allocate(xsum(n1,n2,n3))
            call MPI_REDUCE(xval,xsum,n1*n2*n3,MPI_COMPLEX,&
            &  MPI_SUM,master,spaceComm,ier)
            xval (:,:,:) = xsum(:,:,:)
            deallocate(xsum)
           end if
#endif
end subroutine xsum_master_c3cplx

!--------------------------------------------------------------------

subroutine xsum_master_c4cplx(xval,master,spaceComm,ier)
 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex ,intent(inout) :: xval(:,:,:,:)
 integer ,intent(in) :: master
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer :: n1,n2,n3,n4

#if defined MPI
           complex , allocatable :: xsum(:,:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
           if (nproc_space_comm /= 1) then
            n1 = size(xval,dim=1)
            n2 = size(xval,dim=2)
            n3 = size(xval,dim=3)
            n4 = size(xval,dim=4)
            ! collect xval from processors on master in spaceComm
            allocate(xsum(n1,n2,n3,n4))
            call MPI_REDUCE(xval,xsum,n1*n2*n3*n4,MPI_COMPLEX,&
            &  MPI_SUM,master,spaceComm,ier)
            xval (:,:,:,:) = xsum(:,:,:,:)
            deallocate(xsum)
           end if
#endif
end subroutine xsum_master_c4cplx

!--------------------------------------------------------------------


subroutine xsum_master_c2dpc(xval,master,spaceComm,ier)
 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc) ,intent(inout) :: xval(:,:)
 integer ,intent(in) :: master
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer :: n1,n2

#if defined MPI
           complex(dpc) , allocatable :: xsum(:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
           if (nproc_space_comm /= 1) then
            n1 = size(xval,dim=1)
            n2 = size(xval,dim=2)
            ! collect xval from processors on master in spaceComm
            allocate(xsum(n1,n2))
            call MPI_REDUCE(xval,xsum,n1*n2,MPI_DOUBLE_COMPLEX,&
            &  MPI_SUM,master,spaceComm,ier)
            xval (:,:) = xsum(:,:)
            deallocate(xsum)
           end if
#endif
end subroutine xsum_master_c2dpc

!--------------------------------------------------------------------

subroutine xsum_master_c3dpc(xval,master,spaceComm,ier)
 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc) ,intent(inout) :: xval(:,:,:)
 integer ,intent(in) :: master
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer :: n1,n2,n3

#if defined MPI
           complex(dpc) , allocatable :: xsum(:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
           if (nproc_space_comm /= 1) then
            n1 = size(xval,dim=1)
            n2 = size(xval,dim=2)
            n3 = size(xval,dim=3)
            ! collect xval from processors on master in spaceComm
            allocate(xsum(n1,n2,n3))
            call MPI_REDUCE(xval,xsum,n1*n2*n3,MPI_DOUBLE_COMPLEX,&
            &  MPI_SUM,master,spaceComm,ier)
            xval (:,:,:) = xsum(:,:,:)
            deallocate(xsum)
           end if
#endif
end subroutine xsum_master_c3dpc

!--------------------------------------------------------------------

subroutine xsum_master_c4dpc(xval,master,spaceComm,ier)
 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc) ,intent(inout) :: xval(:,:,:,:)
 integer ,intent(in) :: master
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer :: n1,n2,n3,n4

#if defined MPI
           complex(dpc) , allocatable :: xsum(:,:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
           if (nproc_space_comm /= 1) then
            n1 = size(xval,dim=1)
            n2 = size(xval,dim=2)
            n3 = size(xval,dim=3)
            n4 = size(xval,dim=4)
            ! collect xval from processors on master in spaceComm
            allocate(xsum(n1,n2,n3,n4))
            call MPI_REDUCE(xval,xsum,n1*n2*n3*n4,MPI_DOUBLE_COMPLEX,&
            &  MPI_SUM,master,spaceComm,ier)
            xval (:,:,:,:) = xsum(:,:,:,:)
            deallocate(xsum)
           end if
#endif
end subroutine xsum_master_c4dpc

!--------------------------------------------------------------------

!!***
