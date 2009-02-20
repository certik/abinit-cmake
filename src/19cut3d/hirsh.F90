!{\src2tex{textfont=tt}}
!!****f* ABINIT/hirsh
!! NAME
!! hirsh
!!
!! FUNCTION
!! Compute the Hirshfeld charges
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (XG,MVerstraete)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  grid_den(nrx,nry,nrz)= density on the grid
!!  natom = number of atoms in the unit cell
!!  nrx,nry,nrz= number of points in the grid for the three directions
!!  ntypat=number of types of atoms in unit cell.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  typat(natom)=type of each atom
!!  xcart(3,natom) = different positions of the atoms in the unit cell
!!  zion=(ntypat)gives the ionic charge for each type of atom
!!  znucl(ntypat)=gives the nuclear number for each type of atom
!!
!! OUTPUT
!!  write the Hirshfeld charge decomposition
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!      metric,spline,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine hirsh(grid_den,natom,nrx,nry,nrz,ntypat,rprimd,xcart,typat,zion,znucl)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12geometry
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrx,nry,nrz,ntypat
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: grid_den(nrx,nry,nrz),rprimd(3,3),zion(ntypat)
 real(dp),intent(in) :: znucl(ntypat)
 real(dp),intent(inout) :: xcart(3,natom)

!Local variables -------------------------
!scalars
 integer :: i1,i2,i3,iatom,icell,ierr,igrid,ii,inmax,inmin,ipoint,istep,itypat
 integer :: k1,k2,k3,mcells,mpoint,nfftot,ngoodpoints,npt
 real(dp) :: aa,bb,coeff1,coeff2,coeff3,den,factor,h_inv,hh,maxrad,minimal_den
 real(dp) :: param1,param2,rr,rr2,total_charge,total_weight,total_zion,ucvol,xx
 real(dp) :: yp1,ypn,yy
 character(len=fnlen) :: file_allelectron
!arrays
 integer :: highest(3),lowest(3)
 integer,allocatable :: ncells(:),npoint(:)
 real(dp) :: coordat(3),coordgrid(3),diff(3),gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp) :: vperp(3),width(3)
 real(dp),allocatable :: aeden(:,:),hcharge(:),hweight(:),local_den(:,:,:,:)
 real(dp),allocatable :: radii(:,:),step(:,:),sum_den(:,:,:),work(:)
 real(dp),allocatable :: xcartcells(:,:,:),xred(:,:),yder2(:)

! *********************************************************************

!1. Read the 1D all-electron atomic files
!Store the radii in radii(:,itypat), and the all-electron
!densities in aeden(:,itypat). The number of the last
!point with significant density is stored in npoint(itypat)

 minimal_den=tol6
 mpoint=4000
 allocate(npoint(ntypat),radii(4000,ntypat),aeden(4000,ntypat))
 do itypat=1,ntypat
  write(6, '(a)' )' Please, give the filename of the all-electron density file'
  write(6, '(a,es16.6)' )' for the first type of atom, with atomic number=',znucl(itypat)
  read(5, '(a)' )file_allelectron
  write(6,*)' The name you entered is : ',trim(file_allelectron),ch10
  open (unit=tmp_unit,iostat=ierr, file=trim(file_allelectron),form='formatted',status='old')
  if (ierr/=0) then
   write(6,*)'Problem opening the file => stop'
   stop
  else
   read(tmp_unit, *) param1, param2
   do ipoint=1,mpoint
!   Either the file is finished
    read(tmp_unit, *, end = 888) xx,yy
    radii(ipoint,itypat)=xx
    aeden(ipoint,itypat)=yy
!   Or the density is lower than the minimal significant value
    if(yy<minimal_den)exit
   end do
   888 continue
   npoint(itypat)=ipoint-1
   if(ipoint==mpoint)then
    write(6,*)' hirsh : mpoint is too low, increase its value to match ipoint.'
   end if
  end if
  close(tmp_unit)
 end do

!2. Compute the list of atoms that are sufficiently close to the cell

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 nfftot=nrx*nry*nrz
!DEBUG
!do k3=1,nrz
!do k2=1,nry
!do k1=1,nrx
!total_charge=total_charge+grid_den(k1,k2,k3)
!end do
!end do
!end do
!write(6,*)' total_charge=',total_charge*ucvol/dble(nfftot)
!ENDDEBUG

 allocate(xred(3,natom))
 call xredxcart(natom,-1,rprimd,xcart,xred)

!Compute the widths of the cell
!First width : perpendicular vector length
 vperp(:)=rprimd(:,1)-rprimd(:,2)*rmet(1,2)/rmet(2,2) &
& -rprimd(:,3)*rmet(1,3)/rmet(3,3)
 width(1)=sqrt(dot_product(vperp,vperp))
!Second width
 vperp(:)=rprimd(:,2)-rprimd(:,1)*rmet(2,1)/rmet(1,1) &
& -rprimd(:,3)*rmet(2,3)/rmet(3,3)
 width(2)=sqrt(dot_product(vperp,vperp))
!Third width
 vperp(:)=rprimd(:,3)-rprimd(:,1)*rmet(3,1)/rmet(1,1) &
& -rprimd(:,2)*rmet(3,2)/rmet(2,2)
 width(3)=sqrt(dot_product(vperp,vperp))

!Compute the number of cells that will make up the supercell
 allocate(ncells(natom))
 mcells=0
 do iatom=1,natom
  itypat=typat(iatom)
  maxrad=radii(npoint(itypat),itypat)
! Compute the lower and higher indices of the supercell
! for this atom
  do ii=1,3
   lowest(ii)=floor(xred(ii,iatom)-maxrad/width(ii))-1
   highest(ii)=ceiling(xred(ii,iatom)+maxrad/width(ii))+1
!  lowest(ii)=ceiling(-xred(ii,iatom)-maxrad/width(ii))
!  highest(ii)=floor(-xred(ii,iatom)+maxrad/width(ii)+1)
  end do
  ncells(iatom)=(highest(1)-lowest(1)+1)* &
&  (highest(2)-lowest(2)+1)* &
&  (highest(3)-lowest(3)+1)
! DEBUG
! write(6,*)' maxrad=',maxrad
! write(6,*)' lowest(:)=',lowest(:)
! write(6,*)' highest(:)=',highest(:)
! write(6,*)' ncells(iatom)=',ncells(iatom)
! ENDDEBUG
 end do
 mcells=maxval(ncells(:))

!Compute, for each atom, the set of image atoms
!in the whole supercell
 allocate(xcartcells(3,mcells,natom))
 do iatom=1,natom
  itypat=typat(iatom)
  maxrad=radii(npoint(itypat),itypat)
! Compute the lower and higher indices of the supercell
! for this atom

  do ii=1,3
   lowest(ii)=floor(xred(ii,iatom)-maxrad/width(ii))-1
   highest(ii)=ceiling(xred(ii,iatom)+maxrad/width(ii))+1
  end do
  icell=0
  do i1=lowest(1),highest(1)
   do i2=lowest(2),highest(2)
    do i3=lowest(3),highest(3)
     icell=icell+1
     xcartcells(:,icell,iatom)=xcart(:,iatom)+i1*rprimd(:,1)+i2*rprimd(:,2)+i3*rprimd(:,3)
    end do
   end do
  end do
 end do

!DEBUG
!write(6,*)' before all-electron pro-atom density'
!ENDDEBUG

!Compute, for each atom, the all-electron pro-atom density
!at each point in the primitive cell
 allocate(local_den(nrx,nry,nrz,natom),step(2,mpoint))
 allocate(work(mpoint),yder2(mpoint))
 coeff1=one/nrx
 coeff2=one/nry
 coeff3=one/nrz

 do iatom=1,natom
  itypat=typat(iatom)
  npt=npoint(itypat)
  maxrad=radii(npt,itypat)
  write(6,*)
  write(6, '(a,i4)' )' hirsh : accumulating density for atom ',iatom
! DEBUG
! write(6,*)' ncells(iatom)=',ncells(iatom)
! ENDDEBUG
  do istep=1,npt-1
   step(1,istep)=radii(istep+1,itypat) - radii(istep,itypat)
   step(2,istep)=one/step(1,istep)
  end do
! Approximate first derivative for small radii
  yp1=(aeden(2,itypat)-aeden(1,itypat))/(radii(2,itypat)-radii(1,itypat))
  ypn=zero
  call spline(radii(1:npt,itypat),aeden(1:npt,itypat),npt,yp1,ypn,yder2,work)

  local_den(:,:,:,iatom)=zero

! Big loop on the cells
  do icell=1,ncells(iatom)
!  DEBUG
!  write(6,*)' icell=',icell
!  ENDDEBUG

   coordat(:)=xcartcells(:,icell,iatom)

!  Big loop on the grid points
   do k3 = 1, nrz
    do k2 = 1, nry
     do k1 = 1, nrx
!     DEBUG
!     if(icell>=342)then
!     write(6,*)' k1,k2,k3=',k1,k2,k3
!     write(6,*)' coeff1,coeff2,coeff3=',coeff1,coeff2,coeff3
!     write(6,*)' rprimd=',rprimd
!     write(6,*)' coordat=',coordat
!     end if
!     ENDDEBUG

      coordgrid(:)=rprimd(:,1)*(k1-1)*coeff1+ &
&      rprimd(:,2)*(k2-1)*coeff2+ &
&      rprimd(:,3)*(k3-1)*coeff3
      diff(:)=coordgrid(:)-coordat(:)
      rr2=diff(1)**2+diff(2)**2+diff(3)**2

!     DEBUG
!     if(icell>=342)then
!     write(6,*)' rr2,maxrad',rr2,maxrad
!     end if
!     ENDDEBUG

      if(rr2<maxrad**2)then

       rr=sqrt(rr2)
!      Find the index of the radius by bissection
       if (rr < radii(1,itypat)) then
!       Linear extrapolation
        den=aeden(1,itypat)+(rr-radii(1,itypat))/(radii(2,itypat)-radii(1,itypat))&
&        *(aeden(2,itypat)-aeden(1,itypat))
       else
!       Use the spline interpolation
!       Find the index of the radius by bissection
        inmin=1
        inmax=npt
        igrid=1
        do
         if(inmax-inmin==1)exit
         igrid=(inmin+inmax)/2
         if(rr>=radii(igrid,itypat))then
          inmin=igrid
         else
          inmax=igrid
         end if
        end do
        igrid=inmin
!       DEBUG
!       write(6,*)' igrid',igrid
!       ENDDEBUG

        hh=step(1,igrid)
        h_inv=step(2,igrid)
        aa= (radii(igrid+1,itypat)-rr)*h_inv
        bb= (rr-radii(igrid,itypat))*h_inv
        den = aa*aeden(igrid,itypat) + bb*aeden(igrid+1,itypat)  &
&        +( (aa*aa*aa-aa)*yder2(igrid)         &
&        +(bb*bb*bb-bb)*yder2(igrid+1) ) *hh*hh*sixth
       end if ! Select small radius or spline

       local_den(k1,k2,k3,iatom)=local_den(k1,k2,k3,iatom)+den

      end if ! dist2<maxrad

     end do ! k1
    end do ! k2
   end do ! k3

  end do ! icell

 end do ! iatom

!Compute, the total all-electron density
!at each point in the primitive cell
 allocate(sum_den(nrx,nry,nrz))
 sum_den(:,:,:)=zero
 do iatom=1,natom
  sum_den(:,:,:)=sum_den(:,:,:)+local_den(:,:,:,iatom)
 end do

!DEBUG
!do k3=1,nrz
!do k2=1,nry
!do k1=1,nrx
!write(6, '(3i4,3es16.6)' )k1,k2,k3,local_den(k1,k2,k3,1:2),sum_den(k1,k2,k3)
!end do
!end do
!end do
!write(6,*)' hirsh : before accumulate the integral of the density'
!ENDDEBUG


!Accumulate the integral of the density, to get Hirshfeld charges
!There is a minus sign because the electron has a negative charge
 allocate(hcharge(natom),hweight(natom))
 ngoodpoints = 0
 hcharge(:)=zero
 hweight(:)=zero
 do k3=1,nrz
  do k2=1,nry
   do k1=1,nrx
!   Use minimal_den in order to avoid divide by zero
    if (abs(sum_den(k1,k2,k3)) > minimal_den) then
     ngoodpoints = ngoodpoints+1
     factor=grid_den(k1,k2,k3)/(sum_den(k1,k2,k3)+minimal_den)
     do iatom=1,natom
      hcharge(iatom)=hcharge(iatom)-local_den(k1,k2,k3,iatom)*factor
      hweight(iatom)=hweight(iatom)+local_den(k1,k2,k3,iatom)/(sum_den(k1,k2,k3)+minimal_den)
     end do
    end if
   end do
  end do
 end do
 hcharge(:)=hcharge(:)*ucvol/dble(nfftot)
 hweight(:)=hweight(:)/dble(nfftot)


!Check on the total charge
 total_zion=sum(zion(typat(1:natom)))
 total_charge=sum(hcharge(1:natom))
 total_weight=sum(hweight(1:natom))

!DEBUG
!write(6,*)' ngoodpoints = ', ngoodpoints, ' out of ', nfftot
!write(6,*)' total_weight=',total_weight
!write(6,*)' total_weight=',total_weight
!ENDDEBUG

!Output
 write(6,*)
 write(6,*)'    Atom       Zion       Hirshfeld Charge       Net charge '
 write(6,*)
 do iatom=1,natom
  write(6, '(i9,3es17.6)' )&
&  iatom,zion(typat(iatom)),hcharge(iatom),hcharge(iatom)+zion(typat(iatom))
 end do
 write(6,*)
!write(6,*)'    Total ',total_zion,total_charge,total_charge+total_zion
 write(6, '(a,3es17.6)')'    Total',total_zion,total_charge,total_charge+total_zion
 write(6,*)

 deallocate(aeden,hcharge,local_den)
 deallocate(ncells,npoint,radii,step,sum_den)
 deallocate(work,xcartcells,xred,yder2)

end subroutine hirsh
!!***
