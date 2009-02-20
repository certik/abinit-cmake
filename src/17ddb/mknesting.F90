!{\src2tex{textfont=tt}}
!!****f* ABINIT/mknesting
!! NAME
!! mknesting
!!
!! FUNCTION
!!  Calculate the nesting factor over the dense k-grid,
!!  interpolate the values along a given q path
!!  and write the data on file in the X-Y format or
!!  in the XCrysden format (XSF)
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nkpt = number of k points
!!  kpt(3,nkpt) = k points
!!  nkx, nky, nkz = number of k-point along each direction
!!  nband = number of bands to be considered in the calculation
!!  weight(nband,nkpt) =  integration weights for each k-point and band
!!  npoint_in = number of points requested along the trajectory
!!  nsegment_in = number of segments in the reciprocal space trajectory
!!  qpath_vertices_in = vertices of the reciprocal space trajectory
!!  elph_base_name = prefix of the output file
!!  gprimd(3,3) dimensional reciprocal lattice vectors
!!  prtnest = flags governing the format of the output file
!! OUTPUT
!!  only write to file
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      bfactor,interpol3d,leave_new,printxsf,wrtout,wstoconv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mknesting(nkpt,kpt,kptrlatt,nband,weight,nsegment_in,npoint_in,&
& qpath_vertices_in,elph_base_name,gprim,gprimd,rprim,brav,prtnest)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_17ddb, except_this_one => mknesting
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: brav,nband,nkpt,nsegment_in,prtnest
 character(len=fnlen),intent(in) :: elph_base_name
!arrays
 integer,intent(in) :: kptrlatt(3,3),npoint_in(nsegment_in)
 real(dp),intent(in) :: gprim(3,3),gprimd(3,3),kpt(3,nkpt)
 real(dp),intent(in) :: qpath_vertices_in(3,nsegment_in+1),rprim(3,3)
 real(dp),intent(in) :: weight(nband,nkpt)

!Local variables-------------------------------
!scalars
 integer :: ii,ikpt,indx,iost,ipoint,iseg,jkpt,kindex,maxrank,nkx,nky,nkz
 integer :: nresol,nsegment,option,unit_nest
 real(dp) :: kval,res
 character(len=500) :: message
 character(len=fnlen) :: fname
!arrays
 integer :: ngqpt(3),npoint(nsegment_in+1)
 integer,allocatable :: kptrank(:),ktable(:)
 real(dp) :: origin(3),qpath_vertices(3,nsegment_in+2),qpt(3),tmpkpt(3)
 real(dp),allocatable :: convkpt(:,:),nestfactor(:),nestordered(:)

! *************************************************************************

!DEBUG
!write(*,*)' mknesting : enter '
!ENDDEBUG

 if (     kptrlatt(1,2) /= 0 .or. kptrlatt(1,3) /= 0 .or. kptrlatt(2,1) /= 0       &
& .or. kptrlatt(2,3) /= 0 .or. kptrlatt(3,1) /= 0 .or. kptrlatt(3,2) /= 0 ) then
  write (message,'(7a)')ch10,' mknesting : WARNING-',ch10,                         &
&  ' kptrlatt should be diagonal in order to calculate the nesting factor,',ch10,&
&  ' skipping the nesting factor calculation ',ch10
  call wrtout(06,message,'COLL')
  call wrtout(ab_out,message,'COLL')
  return
 end if

 if (prtnest /= 1 .and. prtnest /= 2) then
  write(message,'(4a)')ch10,' mknesting : BUG-',ch10,&
&  ' prtnest should be 1 or 2'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 write(message,'(a,9i5)')' mknesting : kptrlatt = ',kptrlatt
 call wrtout(06,message,'COLL')

 nkx=kptrlatt(1,1)
 nky=kptrlatt(2,2)
 nkz=kptrlatt(3,3)

 if (nkpt /= nkx*nky*nkz) then
  write(message,'(5a)')ch10,' mknesting : ERROR-',ch10,&
&  ' Wrong input value for kptrlatt',ch10
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 allocate (nestfactor(nkpt))
 nestfactor(:)=zero

!NOTE: input weights are not normalised, the normalisation factor in introduced in bfactor
 call bfactor(nkpt,kpt,nkpt,kpt,weight,nband,nestfactor)

!================================================================================================
!use linear interpolation to plot the bfactor along the given q-path
!1) order the kpoints of the grid putting them in increasing x, then y, then z (FORTRAN convention)
!2) make table from input kpts to ordered kpts
!3) perform interpolation
!================================================================================================

 allocate (convkpt(3,nkpt))
 convkpt(:,:)=zero

!rank is used to order kpoints
 allocate (kptrank(nkpt))
 kptrank(:) = 0

 do ikpt=1,nkpt
  call wstoconv(kpt(1,ikpt),tmpkpt(1),res)
  call wstoconv(kpt(2,ikpt),tmpkpt(2),res)
  call wstoconv(kpt(3,ikpt),tmpkpt(3),res)
  convkpt(:,ikpt)=tmpkpt(:)
! DEBUG
! write(*,*)kpt(1,ikpt),tmpkpt(1)
! write(*,*)kpt(2,ikpt),tmpkpt(2)
! write(*,*)kpt(3,ikpt),tmpkpt(3)
! ENDDEBUG
  kptrank(ikpt) = 100000000.0_dp*(tmpkpt(3)+one) + &
&  100000.0_dp*(tmpkpt(2)+one) + &
&  100.0_dp*(tmpkpt(1)+one)
! DEBUG
! write(*,*)ikpt,tmpkpt,kptrank(ikpt)
! ENDDEBUG
 end do

 allocate (ktable(nkpt))
 ktable(:)=0

 kindex=nkpt
 do ikpt=1,nkpt
  maxrank=maxval(kptrank)
  findmax: do jkpt=1,nkpt
   if(kptrank(jkpt)==maxrank) then
    kptrank(jkpt)=-kptrank(jkpt)
    ktable(jkpt)=kindex
    kindex=kindex-1
    exit findmax
   end if
  end do findmax
 end do !ikpt
!DEBUG
!do ikpt=1,nkpt
!do jkpt=1,nkpt
!if (ktable(jkpt)==ikpt) write(*,'(3es16.8,i16)') convkpt(:,jkpt),kptrank(jkpt)
!end do
!end do
!stop
!ENDDEBUG

!fill the datagrid for the nesting factor using the Fortran convention and the conventional unit cell
!NOTE: the Fortran convention is a must if we want to plot the data
!in the BXSF format, useful for the linear interpolation since we use interpol3d.F90

 allocate (nestordered(nkpt))
 nestordered(:)=zero

 do ikpt=1,nkpt
  do jkpt=1,nkpt
   if (ktable(jkpt)==ikpt) then
    nestordered(ikpt)=nestfactor(jkpt)
!   write(77,'(4es16.8)')kpt(:,jkpt),nestordered(ikpt)
   end if
  end do
 end do

 deallocate (nestfactor)

!Get vertices along the input q-path
!add extra segment for last point: nvertices = nsegments + 1
 qpath_vertices(:,1:nsegment_in+1) = qpath_vertices_in(:,:)
 qpath_vertices(:,nsegment_in+2) = qpath_vertices_in(:,nsegment_in+1)
 npoint(1:nsegment_in) = npoint_in(:)
 npoint(nsegment_in+1) = 1
 nsegment=nsegment_in + 1

!open output file and write header
 unit_nest=110
 fname=trim(elph_base_name) // '_NEST'
 open (unit=unit_nest,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
  write (message,'(2a)')' mknesting : ERROR- opening file ',trim(fname)
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 write (unit_nest,'(a)')'#'
 write (unit_nest,'(a)')'# ABINIT package : Nesting factor file'
 write (unit_nest,'(a)')'#'
 write (unit_nest,'(a,i10,a)')'# Nesting factor calculated on ',sum(npoint),' Q-points'
 write (unit_nest,'(a)')'# Description of the Q-path :'
 write (unit_nest,'(a,i10)')'# Number of line segments = ',nsegment_in
 write (unit_nest,'(a)')'# Vertices of the Q-path and corresponding index = '
 indx=1
 do ii=1,nsegment_in+1
  write (unit_nest,'(a,3(E16.6,1x),i8)')'#  ',qpath_vertices_in(:,ii),indx
  indx=indx+npoint(ii)
 end do
 write (unit_nest,'(a)')'#'

!Get qpoint along the q-path from qpath_vertices(:,iseg)to the point
!just before qpath_vertices(:,iseg+1) and interpolate the nesting factor
 indx=1

 do iseg=1,nsegment
  nresol=npoint(iseg)
  do ipoint=0,npoint(iseg)-1
   qpt(:) = qpath_vertices(:,iseg)  &
&   +dble(ipoint)/dble(npoint(iseg))*(qpath_vertices(:,iseg+1)-qpath_vertices(:,iseg))
   call wstoconv(qpt(1),tmpkpt(1),res)
   call wstoconv(qpt(2),tmpkpt(2),res)
   call wstoconv(qpt(3),tmpkpt(3),res)

   call interpol3d(tmpkpt,nkx,nky,nkz,kval,nestordered)

   write(unit_nest,'(i5,18e16.5)')indx,kval
!  DEBUG
!  write (140,'(i6,1x,4(es16.8,1x))')indx,qpt(:),kval
!  ENDDEBUG
   indx = indx+1
  end do !end ipoint do
 end do !end iseg do
 close (unit_nest)

 if (prtnest==2) then !write also the nest factor in the XSF format
  fname=trim(elph_base_name) // '_NEST_XSF'
  open (unit=unit_nest,file=fname,status='unknown',iostat=iost)
  if (iost /= 0) then
   write (message,'(2a)')' mknesting : ERROR- opening file ',trim(fname)
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  origin(:)=zero
  option=0 !reciprocal space
  call printxsf(nkx,nky,nkz,nestordered,gprimd,origin,unit_nest,option)
  close (unit_nest)
 end if

 deallocate (nestordered)
 deallocate (convkpt)
 deallocate (ktable)
 deallocate (kptrank)

!DEBUG
!write(6,*)' mknesting : exit'
!stop
!ENDDEBUG

end subroutine mknesting
!!***
