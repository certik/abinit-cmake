!{\src2tex{textfont=tt}}
!!****f* ABINIT/drvaim
!! NAME
!! drvaim
!!
!! FUNCTION
!! Main driver for the Bader analysis
!! it looks the values of the input variables
!! and calls corresponding procedures
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  aim_dtset
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  this routine acts primarily on the data contained in the aimprom module
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      aim
!!
!! CHILDREN
!!      addout,aim_follow,cpdrv,critics,graph,initaim,integrho,integvol,plint
!!      rsurf,surf,timein
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine drvaim(aim_dtset)

 use defs_basis
 use defs_parameters
 use defs_aimfields
 use defs_aimprom
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_14bader, except_this_one => drvaim
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(aim_dataset_type),intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: i1,i2,i3,iat,iatinit,ii,inxat,ipos,iposinit,ires,jii,jj,jri,nn
 integer :: npmax,npoints,nstep
 real(dp) :: dmax,dstlim,rr,ss,t1,t2,tf,wall
 logical :: debold,ortho,sfour,srch,sthree,stwo
 character(len=4) :: switch
!arrays
 real(dp) :: ev(3),grho(3),rbv(3),v1(3),v2(3),vt1(3),vt2(3),vv(3),xstart(3)
 real(dp) :: zz(3,3)

! *********************************************************************

!These input variables might be modified during what follows,
!so, they are copied outside of aim_dtset.
 inxat=aim_dtset%batom
 r0=aim_dtset%atrad
 h0=aim_dtset%folstp
 maxatdst=aim_dtset%maxatd
 maxcpdst=aim_dtset%maxcpd

 dstlim=maxcpdst

!Flags from the old version
!to be remove later
 deb=.false.
 stwo=.true.
 sthree=.true.
 sfour=.false.
 srch=.false.

 npmax=aim_npmaxin

!Main initialisation procedure -
!- it reads ABINIT density file and files
!with core densities and initialises the fields for
!spline interpolation

 call initaim(aim_dtset)


!CP SEARCHING

 if (aim_dtset%crit /= 0) then

  if (aim_dtset%crit==3) then
!  old version of the driver for searching CPs (original code)
   call critics(aim_dtset,inxat,stwo,sthree,sfour,dstlim)
  else
!  driver for searching CPs with Popellier algorithm
   call cpdrv(aim_dtset)
  end if

 end if

!
!BADER SURFACE CALCUL
!

 if (aim_dtset%isurf==1) then
! driver for determination of the Bader surface
  call surf(aim_dtset)
 end if

!
!CHARGE INTEGRATIOM
!

 if (aim_dtset%irho==1) then
  call integrho(aim_dtset)
 end if

!
!VOLUME INTEGRATION OF THE BADER ATOM
!

 if (aim_dtset%ivol==1) then
  call integvol()
 end if

!
!ONE RADIUS OF THE BADER SURFACE
!

 if (aim_dtset%irsur==1) then
  if (aim_dtset%isurf/=0) srch=.true.
  iat=aim_dtset%batom
  ss=r0
  call timein(t1,wall)
  call rsurf(aim_dtset,rr,grho,aim_dtset%th0,aim_dtset%phi0,ss,iat,npmax,srch)
  call timein(t2,wall)
  t2=t2-t1
  write(unts,'(2F12.8,F15.10)') aim_dtset%th0,aim_dtset%phi0,rr
  write(6,'(":RSUR ",2F12.8,2F15.10)') aim_dtset%th0,aim_dtset%phi0,rr,t2
 end if

!
!FOLLOW THE GRADIENT PATH FROM ONE POINT
!

 if (aim_dtset%foll==1) then
  iatinit=aim_dtset%batom
  iposinit=batcell
  if (aim_dtset%isurf/=0) srch=.true.
  debold=deb
  xstart(:)=aim_dtset%foldep(:)
  call timein(t1,wall)
  call aim_follow(aim_dtset,xstart,npmax,srch,iatinit,iposinit,iat,ipos,nstep)
  call timein(t2,wall)
  tf=t2-t1
  write(6,'(":TIME in aim_follow:", F12.4)') tf
 end if

 if (aim_dtset%plden == 1) then
! profile of the density integrated in plane xy
! belong the z-axes - not finished - cut3d better !
  call plint()
 end if

 if ((aim_dtset%denout > 0).or.(aim_dtset%lapout > 0)) then
! additional outputs of density and laplacian fields
! in the plane or line
  call addout(aim_dtset)
 end if

 if (aim_dtset%gpsurf == 1) then
! script for gnuplot - simple demonstration of the
! computed surface
  call graph(unts,untg)
 end if

!Deallocation of global variables allocated in initaim
!and declared in defs_aimfields.
 deallocate(dig1,dig2,dig3,llg1,llg2,llg3,cdig1,cdig2,cdig3)
 deallocate(ddx,ddy,ddz)
 deallocate(rrad,crho,sp2,sp3,sp4,corlim)
 deallocate(dvl)
 deallocate(ndat,rminl)
!Deallocation of global variables allocated in initaim
!and declared in defs_aimprom.
 deallocate(typat,xred,xatm)

end subroutine drvaim
!!***
