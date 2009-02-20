!{\src2tex{textfont=tt}}
!!****f* ABINIT/hessinit
!! NAME
!! hessinit
!!
!! FUNCTION
!! Initiliase an Hessian matrix, either from disk or using init_matrix.
!! The size ndim must be greater or equal than 3 * dtset%natom.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG, JCC).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  init_matrix(3,3)=matrix used for each atom (if iatfix = 0) for initialisation.
!!  ndim=size of the hessian and vectors
!!  ucvol=volume of the box (used when dtset%optcell is not null).
!!
!! OUTPUT
!!  hessin(ndim,ndim)=hessian matrix, initialised at output.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine hessinit(dtfil, dtset, hessin, init_matrix, ndim, ucvol)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndim
 real(dp),intent(in) :: ucvol
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
!arrays
 real(dp),intent(in) :: init_matrix(3,3)
 real(dp),intent(out) :: hessin(ndim,ndim)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=5
 integer :: hess_ok,iatom,idim,idir1,idir2,iexit,ii,ios,jj,ndim0
 real(dp) :: diag
 logical :: ex
 character(len=500) :: message
 character(len=fnlen) :: filhess
!arrays
 real(dp),allocatable :: bfgs(:),din(:),dout(:),hdelta(:)

! *********************************************************************

!Initialization of the inverse Hessian to the unit matrix
!Much better choices are possible--this simply corresponds to
!taking the first minimization step as the negative of the
!gradient, with the full length of the gradient vector as
!the step size.  Any spring type model would probably be a better
!starting guess.

 call status(0,dtfil%filstat,iexit,level,'hessian       ')
 
 if (ndim < 3 * dtset%natom) then
  write(message, '(a,a,a,a,a,a)' ) ch10,&
&  ' hessinit: ERROR -',ch10,&
&  '  the size of the given hessian matrix is to small.', ch10, &
&  '  This is an internal error, contact ABINIT developers.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!Special arrangement: if input hessian file exists, read data from there
 filhess=trim(dtfil%filnam_ds(3))//'_HES'
 inquire (file=filhess,iostat=ios,exist=ex)
 hess_ok=0
 if (ex) then
! Read inverse hessian data from file; format is
  open (unit=tmp_unit,file=filhess,form='formatted',status='old')
  read (tmp_unit,*)
  read (tmp_unit,*) ndim0
  if (ndim0/=ndim) then
!  Cannot read data because data file natom disagrees with current job
   write(message, '(8a,i10,a,i10,2a)' ) ch10,&
&   ' hessinit: WARNING -',ch10,&
&   '  Tried to read inverse hessian from file',trim(filhess),&
&   ' but',ch10,'  ndim of that file =',ndim0,&
&   ' , is not equal to input ndim =',ndim,ch10,&
&   ' => initialize inverse hessian with identity matrix.'
   call wrtout(06,message,'COLL')
   close(unit=tmp_unit)
  else
!  Read inverse hessian
   do jj=1,ndim
    read (tmp_unit,*)
    read (tmp_unit,*) (hessin(ii,jj),ii=1,ndim)
   end do
   close (unit=tmp_unit)
   write(message,*)' Inverse hessian has been input from input hessian file',&
&   trim(filhess)
   call wrtout(06,message,'COLL')
   hess_ok=1
  end if
 end if

!If hessin was not read, initialize inverse hessian with identity matrix
!in cartesian coordinates, which makes use of metric tensor gmet
!in reduced coordinates.
 if(hess_ok==0)then
  hessin(:,:)=zero
  do iatom=1,dtset%natom
   do idir1=1,3
    do idir2=1,3
!    Warning : implemented in reduced coordinates
     if ( dtset%iatfix(idir1,iatom) ==0 .and. dtset%iatfix(idir2,iatom) ==0 )then
      hessin(idir1+3*(iatom-1),idir2+3*(iatom-1))=init_matrix(idir1,idir2)
     end if
    end do
   end do
  end do
  if(dtset%optcell/=0)then
!  These values might lead to too large changes in some cases ...
   diag=dtset%strprecon*30.0_dp/ucvol
   if(dtset%optcell==1)diag=diag/three
   do idim=3*dtset%natom+1,ndim
    hessin(idim,idim)=diag
   end do
  end if
  write(message, '(a)' )' Inverse hessian has been initialized.'
  call wrtout(06,message,'COLL')
 end if

end subroutine hessinit
!!***
