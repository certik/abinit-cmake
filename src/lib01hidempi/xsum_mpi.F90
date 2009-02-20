!{\src2tex{textfont=tt}}
!!****f* ABINIT/xsum_mpi
!! NAME
!! xsum_mpi
!!
!! FUNCTION
!! This module contains functions that calls MPI routine,
!! if we compile the code using the MPI CPP flags.
!! xsum_mpi is the generic function.
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

!--------------------------------------------------------------------

subroutine xsum_mpi_int(xval,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments----------------
 integer,intent(inout) :: xval(:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
 integer :: n1

#if defined MPI
           integer , allocatable :: xsum(:)
#endif

 ier=0
#if defined MPI
           if (spaceComm /= MPI_COMM_SELF) then
!          Accumulate xval on all proc. in spaceComm
            n1 = size(xval)
            allocate(xsum(n1))
            call MPI_ALLREDUCE(xval,xsum,n1,MPI_INTEGER,&
            &  MPI_SUM,spaceComm,ier)
            xval (:) = xsum(:)
            deallocate(xsum)
           end if
#endif
end subroutine xsum_mpi_int

!--------------------------------------------------------------------

subroutine xsum_mpi_intv(xval,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments----------------------
 integer,intent(inout) :: xval
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables----------------
#if defined MPI
           integer  :: xsum
#endif

 ier=0
#if defined MPI
           if (spaceComm /= MPI_COMM_SELF) then
           !Accumulate xval on all proc. in spaceComm
            call MPI_ALLREDUCE(xval,xsum,1,MPI_INTEGER,&
            &  MPI_SUM,spaceComm,ier)
            xval = xsum
           end if
#endif
end subroutine xsum_mpi_intv

!--------------------------------------------------------------------

subroutine xsum_mpi_intv2(xval,xsum,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments---------------------
 integer,intent(inout) :: xval,xsum
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables---------------

 ier=0
#if defined MPI
           if (spaceComm /=  MPI_COMM_SELF) then
           !Accumulate xval on all proc. in spaceComm
            call MPI_ALLREDUCE(xval,xsum,1,MPI_INTEGER,&
            &  MPI_SUM,spaceComm,ier)
           end if
#endif
end subroutine xsum_mpi_intv2

!--------------------------------------------------------------------

subroutine xsum_mpi_intn(xval,n1,spaceComm,ier)

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
 integer,intent(in)    :: n1
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
#if defined MPI
           integer , allocatable :: xsum(:)
           integer :: nproc_space_comm
#endif


 ier=0
#if defined MPI
           !Accumulate xval on all proc. in spaceComm
           if (spaceComm /=  MPI_COMM_SELF) then
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
             allocate(xsum(n1))
             call MPI_ALLREDUCE(xval,xsum,n1,MPI_INTEGER,&
             &  MPI_SUM,spaceComm,ier)
             xval (:) = xsum(:)
             deallocate(xsum)
            end if
           end if
#endif
end subroutine xsum_mpi_intn

!--------------------------------------------------------------------

subroutine xsum_mpi_int2t(xval,xsum,n1,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer ,intent(inout) :: xval(:),xsum(:)
 integer ,intent(in)    :: n1
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

 ier=0
#if defined MPI 
            if (spaceComm /=  MPI_COMM_SELF) then
             !Accumulate xval on all proc. in spaceComm
             call MPI_ALLREDUCE(xval,xsum,n1,MPI_INTEGER,&
             &  MPI_SUM,spaceComm,ier)
            end if
#endif
end subroutine xsum_mpi_int2t


!--------------------------------------------------------------------

subroutine xsum_mpi_int2d(xval,spaceComm,ier)

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
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
 integer ::   n1,n2

#if defined MPI
           integer , allocatable :: xsum(:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           if (spaceComm /=  MPI_COMM_SELF) then
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
             !Accumulate xval on all proc. in spaceComm
             n1 =size(xval,dim=1)
             n2 =size(xval,dim=2)
             allocate(xsum(n1,n2))
             call MPI_ALLREDUCE(xval,xsum,n1*n2,MPI_INTEGER,&
             &  MPI_SUM,spaceComm,ier)
             xval (:,:) = xsum(:,:)
             deallocate(xsum)
            end if
           end if
#endif
end subroutine xsum_mpi_int2d

!--------------------------------------------------------------------

subroutine xsum_mpi_int3d(xval,spaceComm,ier)

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
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
 integer ::   n1,n2,n3

#if defined MPI
           integer , allocatable :: xsum(:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           !Accumulate xval on all proc. in spaceComm
           if (spaceComm /=  MPI_COMM_SELF) then
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
             n1 =size(xval,dim=1)
             n2 =size(xval,dim=2)
             n3 =size(xval,dim=3)
             allocate(xsum(n1,n2,n3))
             call MPI_ALLREDUCE(xval,xsum,n1*n2*n3,MPI_INTEGER,&
             &  MPI_SUM,spaceComm,ier)
             xval (:,:,:) = xsum(:,:,:)
             deallocate(xsum)
            end if
           end if
#endif
end subroutine xsum_mpi_int3d

!--------------------------------------------------------------------

subroutine xsum_mpi_int4d(xval,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(inout) :: xval(:,:,:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
 integer ::   n1,n2,n3,n4

#if defined MPI
           integer , allocatable :: xsum(:,:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           !Accumulate xval on all proc. in spaceComm
           if (spaceComm /=  MPI_COMM_SELF) then
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
             n1 =size(xval,dim=1)
             n2 =size(xval,dim=2)
             n3 =size(xval,dim=3)
             n4 =size(xval,dim=4)
             allocate(xsum(n1,n2,n3,n4))
             call MPI_ALLREDUCE(xval,xsum,n1*n2*n3*n4,MPI_INTEGER,&
             &  MPI_SUM,spaceComm,ier)
             xval (:,:,:,:) = xsum(:,:,:,:)
             deallocate(xsum)
            end if
           end if
#endif
end subroutine xsum_mpi_int4d

!--------------------------------------------------------------------

subroutine xsum_mpi_dp(xval,spaceComm,ier)

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
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer ::   n1

#if defined MPI
           real(dp) , allocatable :: xsum(:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           if (spaceComm /=  MPI_COMM_SELF) then
           !Accumulate xval on all proc. in spaceComm
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
             n1 = size(xval)
             allocate(xsum(n1))
             call MPI_ALLREDUCE(xval,xsum,n1,MPI_DOUBLE_PRECISION,&
             &  MPI_SUM,spaceComm,ier)
             xval (:) = xsum(:)
             deallocate(xsum)
            end if
           end if
#endif

end subroutine xsum_mpi_dp

!--------------------------------------------------------------------

subroutine xsum_mpi_dpvt(xval,xsum,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in) :: xval
 real(dp),intent(out) :: xsum
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
#if defined MPI
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           !Accumulate xval on all proc. in spaceComm
           if (spaceComm /=  MPI_COMM_SELF) then
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
             call MPI_ALLREDUCE(xval,xsum,1,MPI_DOUBLE_PRECISION,&
             &  MPI_SUM,spaceComm,ier)
            end if
           end if
#endif
end subroutine xsum_mpi_dpvt

!--------------------------------------------------------------------

subroutine xsum_mpi_dpv(xval,spaceComm,ier)

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
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
#if defined MPI
           integer :: nproc_space_comm
           real(dp)  :: xsum
#endif

 ier=0
#if defined MPI
           !Accumulate xval on all proc. in spaceComm
           if (spaceComm /=  MPI_COMM_SELF) then
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
             call MPI_ALLREDUCE(xval,xsum,1,MPI_DOUBLE_PRECISION,&
             &  MPI_SUM,spaceComm,ier)
             xval  = xsum
            end if
           end if
#endif
end subroutine xsum_mpi_dpv

!--------------------------------------------------------------------

subroutine xsum_mpi_dpn(xval,n1,spaceComm,ier)

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
 integer ,intent(in)    :: n1
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
#if defined MPI
           integer :: nproc_space_comm
           real(dp) , allocatable :: xsum(:)
#endif

 ier=0
#if defined MPI
           if (spaceComm /=  MPI_COMM_SELF) then
            !Accumulate xval on all proc. in spaceComm
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
             allocate(xsum(n1))
             call MPI_ALLREDUCE(xval,xsum,n1,MPI_DOUBLE_PRECISION,&
             &  MPI_SUM,spaceComm,ier)
             xval (:) = xsum(:)
             deallocate(xsum)
            end if
           end if
#endif
end subroutine xsum_mpi_dpn

!--------------------------------------------------------------------

subroutine xsum_mpi_dp2d(xval,spaceComm,ier)

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
           if (spaceComm /=  MPI_COMM_SELF) then
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
             n1 = size(xval,dim=1)
             n2 = size(xval,dim=2)
             !Accumulate xval on all proc. in spaceComm
             allocate(xsum(n1,n2))
             call MPI_ALLREDUCE(xval,xsum,n1*n2,MPI_DOUBLE_PRECISION,&
             &  MPI_SUM,spaceComm,ier)
             xval (:,:) = xsum(:,:)
             deallocate(xsum)
            end if
           end if
#endif
end subroutine xsum_mpi_dp2d

!--------------------------------------------------------------------

subroutine xsum_mpi_dp3d(xval,spaceComm,ier)

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
           if (spaceComm /=  MPI_COMM_SELF) then
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
             n1 = size(xval,dim=1)
             n2 = size(xval,dim=2)
             n3 = size(xval,dim=3)
             !Accumulate xval on all proc. in spaceComm
             allocate(xsum(n1,n2,n3))
             call MPI_ALLREDUCE(xval,xsum,n1*n2*n3,MPI_DOUBLE_PRECISION,&
             &  MPI_SUM,spaceComm,ier)
             xval (:,:,:) = xsum(:,:,:)
             deallocate(xsum)
            end if
           end if
#endif
end subroutine xsum_mpi_dp3d

!--------------------------------------------------------------------

subroutine xsum_mpi_dp4d(xval,spaceComm,ier)

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
           if (spaceComm /=  MPI_COMM_SELF) then
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
             n1 = size(xval,dim=1)
             n2 = size(xval,dim=2)
             n3 = size(xval,dim=3)
             n4 = size(xval,dim=4)
             !Accumulate xval on all proc. in spaceComm
             allocate(xsum(n1,n2,n3,n4))
             call MPI_ALLREDUCE(xval,xsum,n1*n2*n3*n4,MPI_DOUBLE_PRECISION,&
             &  MPI_SUM,spaceComm,ier)
             xval (:,:,:,:) = xsum(:,:,:,:)
             deallocate(xsum)
            end if
           end if
#endif
end subroutine xsum_mpi_dp4d

!--------------------------------------------------------------------

subroutine xsum_mpi_dp5d(xval,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:,:,:,:)
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer :: n1,n2,n3,n4,n5

#if defined MPI
           real(dp) , allocatable :: xsum(:,:,:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           if (spaceComm /=  MPI_COMM_SELF) then
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
             n1 = size(xval,dim=1)
             n2 = size(xval,dim=2)
             n3 = size(xval,dim=3)
             n4 = size(xval,dim=4)
             n5 = size(xval,dim=5)
             !Accumulate xval on all proc. in spaceComm
             allocate(xsum(n1,n2,n3,n4,n5))
             call MPI_ALLREDUCE(xval,xsum,n1*n2*n3*n4*n5,MPI_DOUBLE_PRECISION,&
             &  MPI_SUM,spaceComm,ier)
             xval (:,:,:,:,:) = xsum(:,:,:,:,:)
             deallocate(xsum)
            end if
           end if
#endif
end subroutine xsum_mpi_dp5d

!--------------------------------------------------------------------

subroutine xsum_mpi_dp6d(xval,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:,:,:,:,:)
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
 integer :: n1,n2,n3,n4,n5,n6

#if defined MPI
           real(dp) , allocatable :: xsum(:,:,:,:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           if (spaceComm /=  MPI_COMM_SELF) then
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
             call MPI_ALLREDUCE(xval,xsum,n1*n2*n3*n4*n5*n6,MPI_DOUBLE_PRECISION,&
             &  MPI_SUM,spaceComm,ier)
             xval (:,:,:,:,:,:) = xsum(:,:,:,:,:,:)
             deallocate(xsum)
            end if
           end if
#endif
end subroutine xsum_mpi_dp6d

!--------------------------------------------------------------------

subroutine xsum_mpi_dp2t(xval,xsum,n1,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:),xsum(:)
 integer ,intent(in)    :: n1
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

 ier=0
#if defined MPI
            if (spaceComm /=  MPI_COMM_SELF) then
             !Accumulate xval on all proc. in spaceComm
             call MPI_ALLREDUCE(xval,xsum,n1,MPI_DOUBLE_PRECISION,&
             &  MPI_SUM,spaceComm,ier)
            end if
#endif
end subroutine xsum_mpi_dp2t

!--------------------------------------------------------------------

subroutine xsum_mpi_dp3d2t(xval,xsum,n1,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:,:),xsum(:,:,:)
 integer ,intent(in)    :: n1
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

 ier=0
#if defined MPI
           if (spaceComm /=  MPI_COMM_SELF) then
            !Accumulate xval on all proc. in spaceComm
            call MPI_ALLREDUCE(xval,xsum,n1,MPI_DOUBLE_PRECISION,&
            &  MPI_SUM,spaceComm,ier)
           end if
#endif
end subroutine xsum_mpi_dp3d2t

!--------------------------------------------------------------------

subroutine xsum_mpi_dp4d2t(xval,xsum,n1,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:,:,:),xsum(:,:,:,:)
 integer ,intent(in)    :: n1
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

 ier=0
#if defined MPI
           if (spaceComm /=  MPI_COMM_SELF) then
            !Accumulate xval on all proc. in spaceComm
            call MPI_ALLREDUCE(xval,xsum,n1,MPI_DOUBLE_PRECISION,&
            &  MPI_SUM,spaceComm,ier)
           end if
#endif
end subroutine xsum_mpi_dp4d2t



!--------------------------------------------------------------------
subroutine xsum_mpi_c2dc(xval,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc),intent(inout) :: xval(:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
 integer ::   n1,n2

#if defined MPI
           complex(dpc) , allocatable :: xsum(:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           if (spaceComm /= MPI_COMM_SELF) then
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
             !Accumulate xval on all proc. in spaceComm
             n1 =size(xval,dim=1)
             n2 =size(xval,dim=2)
             allocate(xsum(n1,n2))
             call MPI_ALLREDUCE(xval,xsum,n1*n2,MPI_DOUBLE_COMPLEX,&
             &  MPI_SUM,spaceComm,ier)
             xval (:,:) = xsum(:,:)
             deallocate(xsum)
            end if
           end if
#endif
end subroutine xsum_mpi_c2dc

!--------------------------------------------------------------------
subroutine xsum_mpi_c3dc(xval,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc),intent(inout) :: xval(:,:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
 integer ::   n1,n2,n3

#if defined MPI
           complex(dpc) , allocatable :: xsum(:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           if (spaceComm /= MPI_COMM_SELF) then
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
             !Accumulate xval on all proc. in spaceComm
             n1 =size(xval,dim=1)
             n2 =size(xval,dim=2)
             n3 =size(xval,dim=3)
             allocate(xsum(n1,n2,n3))
             call MPI_ALLREDUCE(xval,xsum,n1*n2*n3,MPI_DOUBLE_COMPLEX,&
             &  MPI_SUM,spaceComm,ier)
             xval (:,:,:) = xsum(:,:,:)
             deallocate(xsum)
            end if
           end if
#endif
end subroutine xsum_mpi_c3dc

!--------------------------------------------------------------------
subroutine xsum_mpi_c1cplx(xval,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

 complex,intent(inout) :: xval(:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

#if defined MPI
           integer ::   n1
           complex, allocatable :: xsum(:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
                if (spaceComm /= MPI_COMM_SELF) then
                 call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
                 if (nproc_space_comm /= 1) then
          !Accumulate xval on all proc. in spaceComm
                   n1 =size(xval,dim=1)
                   allocate(xsum(n1))
                   call MPI_ALLREDUCE(xval,xsum,n1,MPI_COMPLEX,&
          &  MPI_SUM,spaceComm,ier)
                  xval (:) = xsum(:)
                  deallocate(xsum)
                 end if
                end if
#endif
end subroutine xsum_mpi_c1cplx


subroutine xsum_mpi_c2cplx(xval,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

 complex,intent(inout) :: xval(:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

#if defined MPI
           integer ::   n1,n2
           complex, allocatable :: xsum(:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           if (spaceComm /= MPI_COMM_SELF) then
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
             !Accumulate xval on all proc. in spaceComm
              n1 =size(xval,dim=1)
              n2 =size(xval,dim=2)
              allocate(xsum(n1,n2))
              call MPI_ALLREDUCE(xval,xsum,n1*n2,MPI_COMPLEX,&
              &  MPI_SUM,spaceComm,ier)
              xval (:,:) = xsum(:,:)
              deallocate(xsum)
             end if
            end if
#endif
end subroutine xsum_mpi_c2cplx

!--------------------------------------------------------------------
subroutine xsum_mpi_c3cplx(xval,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

 complex,intent(inout) :: xval(:,:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

#if defined MPI
           integer ::   n1,n2,n3
           complex , allocatable :: xsum(:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           if (spaceComm /= MPI_COMM_SELF) then
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
            !Accumulate xval on all proc. in spaceComm
            n1 =size(xval,dim=1)
            n2 =size(xval,dim=2)
            n3 =size(xval,dim=3)
            allocate(xsum(n1,n2,n3))
            call MPI_ALLREDUCE(xval,xsum,n1*n2*n3,MPI_COMPLEX,&
            &  MPI_SUM,spaceComm,ier)
            xval (:,:,:) = xsum(:,:,:)
            deallocate(xsum)
           end if
          end if
#endif
end subroutine xsum_mpi_c3cplx
!--------------------------------------------------------------------
subroutine xsum_mpi_c4cplx(xval,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

 complex,intent(inout) :: xval(:,:,:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

#if defined MPI
           integer ::   n1,n2,n3,n4
           complex , allocatable :: xsum(:,:,:,:)
           integer :: nproc_space_comm
#endif

 ier=0
#if defined MPI
           if (spaceComm /= MPI_COMM_SELF) then
            call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
            if (nproc_space_comm /= 1) then
            !Accumulate xval on all proc. in spaceComm
            n1 =size(xval,dim=1)
            n2 =size(xval,dim=2)
            n3 =size(xval,dim=3)
            n4 =size(xval,dim=4)
            allocate(xsum(n1,n2,n3,n4))
            call MPI_ALLREDUCE(xval,xsum,n1*n2*n3*n4,MPI_COMPLEX,&
            &  MPI_SUM,spaceComm,ier)
            xval (:,:,:,:) = xsum(:,:,:,:)
            deallocate(xsum)
           end if
          end if
#endif
end subroutine xsum_mpi_c4cplx

!!***
