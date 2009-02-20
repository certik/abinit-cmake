!{\src2tex{textfont=tt}}
!!****f* ABINIT/addout
!! NAME
!! addout
!!
!! FUNCTION
!! Output density and laplacian (see input variables denout and lapout)
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  aim_dtset=all input variables for aim
!!  also, uses the variables saved in the module "defs_aimprom"
!!
!! OUTPUT
!!  (print)
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet) : the use
!! of a module to transfer data should be avoided
!!
!! PARENTS
!!      drvaim
!!
!! CHILDREN
!!      vgh_rho
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine addout(aim_dtset)

 use defs_basis
 use defs_aimprom
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_14bader, except_this_one => addout
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(aim_dataset_type),intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: cod,dims,iat,ii,ipos,jj,nn,tgrd,volba
 real(dp) :: alfa,rho,rr,xx,yy
!arrays
 real(dp) :: grho(3),hrho(3,3),orig(3),vv(3)
 real(dp),allocatable :: dfld(:),lfld(:),nr(:),stp(:),uu(:,:)

!************************************************************************
 orig(:)=aim_dtset%vpts(:,1)
 if (aim_dtset%denout > 0) then
  dims=aim_dtset%denout
 elseif (aim_dtset%lapout > 0) then
  dims=aim_dtset%lapout
 end if

 select case (aim_dtset%dltyp)
  case (1)
   cod=1
  case (2)
   cod=2
   case default
   cod=0
 end select

 allocate(uu(3,dims),nr(dims),stp(dims))

 write(unto,*) 'grid:', aim_dtset%ngrid(1:dims)
 write(unto,*) 'kod :', cod
 tgrd=1
 do ii=1,dims
  tgrd=tgrd*aim_dtset%ngrid(ii)
  uu(:,ii)=aim_dtset%vpts(:,ii+1)-aim_dtset%vpts(:,1)
  nr(ii)=vnorm(uu(:,ii),0)
  stp(ii)=nr(ii)/(aim_dtset%ngrid(ii)-1)
  uu(:,ii)=uu(:,ii)/nr(ii)
 end do
 write(unto,*) 'tgrd :', tgrd
 do ii=1,dims
  write(unto,*) 'uu :', uu(1:3,ii)
 end do

 if (aim_dtset%denout > 0) then
  allocate(dfld(tgrd+1))
  dfld(:)=0._dp
 end if
 if (aim_dtset%lapout > 0) allocate(lfld(tgrd+1))

 select case (dims)
  case (1)
   nn=0
   do ii=0,aim_dtset%ngrid(1)-1
    nn=nn+1
    vv(:)=orig(:)+ii*stp(1)*uu(:,1)
    call vgh_rho(vv,rho,grho,hrho,rr,iat,ipos,cod)
    if (aim_dtset%denout > 0) dfld(nn)=rho
    if (aim_dtset%lapout > 0) lfld(nn)=hrho(1,1)+hrho(2,2)+hrho(3,3)
   end do
   if (aim_dtset%denout==1) then
    do ii=0,aim_dtset%ngrid(1)-1
     xx=ii*stp(1)
     write(untd,'(2E16.8)') xx, dfld(ii+1)
    end do
   end if
   if (aim_dtset%lapout==1) then
    do ii=0,aim_dtset%ngrid(1)-1
     xx=ii*stp(1)
     write(untl,'(2E16.8)') xx, lfld(ii+1)
    end do
   end if
  case (2)
   nn=0
   alfa=dot_product(uu(:,1),uu(:,2))
   alfa=acos(alfa)
   do ii=0,aim_dtset%ngrid(2)-1
    do jj=0,aim_dtset%ngrid(1)-1
     nn=nn+1
     vv(:)=orig(:)+jj*uu(:,2)*stp(2)+ii*stp(1)*uu(:,1)
     call vgh_rho(vv,rho,grho,hrho,rr,iat,ipos,cod)
     if (aim_dtset%denout > 0) dfld(nn)=rho
     if (aim_dtset%lapout > 0) lfld(nn)=hrho(1,1)+hrho(2,2)+hrho(3,3)
    end do
   end do
   write(unto,*) 'generace hotova', nn
   nn=0
   if (aim_dtset%denout==2) then
    do ii=0,aim_dtset%ngrid(2)-1
     do jj=0,aim_dtset%ngrid(1)-1
      nn=nn+1
      xx=jj*stp(1)+cos(alfa)*ii*stp(2)
      yy=sin(alfa)*ii*stp(2)
      write(untd,'(3E16.8)') xx, yy, dfld(nn)
     end do
     write(untd,*) ' '
    end do
   end if
   nn=0
   if (aim_dtset%lapout==2) then
    write(unto,*) 'lezes sem?'
    do ii=0,aim_dtset%ngrid(2)-1
     do jj=0,aim_dtset%ngrid(1)-1
      nn=nn+1
      xx=jj*stp(1)+cos(alfa)*ii*stp(2)
      yy=sin(alfa)*ii*stp(2)
      write(untl,'(3E16.8)') xx, yy, lfld(nn)
     end do
     write(untl,*) ' '
    end do
   end if
 end select
 deallocate(uu,stp,nr)
 if(aim_dtset%denout>0)deallocate(dfld)
 if(aim_dtset%lapout>0)deallocate(lfld)

end subroutine addout
!!***
