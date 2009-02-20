!{\src2tex{textfont=tt}}
!!****f* ABINIT/thmeig
!! NAME
!! thmeig
!!
!! FUNCTION
!! This routine calculates the thermal corrections to the eigenvalues.
!! The output is this quantity for the input k point.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (PB, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!! 
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine thmeig(acell,amu,blkval2,eigvec,filnam,kpnt,mband,msize,natom,nkpt,nqpt,ntemper,ntypat,&
&                 phfreq,qphon,rprim,temperinc,tempermin,typat,wtq,xred)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,msize,natom,nkpt,nqpt,ntemper,ntypat
 real(dp),intent(in) :: temperinc,tempermin
 character(len=fnlen),intent(in) :: filnam
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: acell(3),amu(ntypat),blkval2(2,msize,mband,nkpt,nqpt)
 real(dp),intent(in) :: eigvec(2,3,natom,3*natom,nqpt),kpnt(3,nkpt,nqpt)
 real(dp),intent(in) :: phfreq(3*natom,nqpt),qphon(9,nqpt),rprim(3,3)
 real(dp),intent(in) :: wtq(3,nqpt),xred(3,natom)

!Local variables-------------------------------
!tolerance for non degenerated levels
!scalars
 integer :: gqpt,iatom1,iatom2,iband,idir1,idir2,ii,ikpt,imod,imod1,index
 integer :: index2,iqpt,isize,itemper,unitout
 real(dp),parameter :: qtol=2.0d-8
 real(dp) :: bosein,dot,fact1i,fact1r,fact2i,fact2r,facti,factr,im,qnrm,tmp
 real(dp) :: vec1i,vec1r,vec2i,vec2r,veci,vecr
 character(len=fnlen) :: outfile
!arrays
 real(dp) :: dedni(mband,nkpt),dednr(mband,nkpt),deigi(mband,nkpt)
 real(dp) :: deigr(mband,nkpt),dred(3),dwtermi(mband,nkpt),dwtermr(mband,nkpt)
 real(dp) :: multi(mband,nkpt),multr(mband,nkpt),norm(3),slope(2,mband,nkpt)
 real(dp) :: thmeigen(2,mband,nkpt),wghtq(nqpt),zeropoint(2,mband,nkpt)

! *********************************************************************

!wghtq(:)=wtq(1,:)
 wghtq(:)=2.0/nkpt
 slope(:,:,:) = zero
 zeropoint(:,:,:) = zero
 write(*,*)'wghtq ',wghtq
!unitout should be attributed in dtset to avoid conflicts
 unitout = 42
 outfile = trim(filnam)//"_TBS"

!open TBS file
 open (unit=unitout,file=outfile,form='formatted',status='unknown')

!Calculating the directions (usefull because derivatives are in reduced coordinates)
 do ii=1,3
  norm(ii) = acell(ii) * sqrt(rprim(1,ii)**2+rprim(2,ii)**2+rprim(3,ii)**2)
 end do

!Finding the Gamma point
 do iqpt=1,nqpt
  qnrm = qphon(1,iqpt)*qphon(1,iqpt)+qphon(2,iqpt)*qphon(2,iqpt)+qphon(3,iqpt)*qphon(3,iqpt)
  if(qnrm<qtol) gqpt=iqpt
 end do
 write(*,*)'thmeig: gqpt',gqpt
 write(unitout,'(a)')'thmeig: Thermal Eigenvalue corrections'
!Loop on temperatures
 do itemper= 1, ntemper
  tmp=tempermin+temperinc*float(itemper-1)
  thmeigen(:,:,:) = zero
! Sum on all phonon wavevectors and modes
  do iqpt=1,nqpt
   index2=0
   do imod=1,3*natom

!   Calculate the derivative
    deigr(:,:) = zero
    deigi(:,:) = zero
    dwtermr(:,:)=zero
    dwtermi(:,:)=zero
    dednr(:,:)=zero
    dedni(:,:)=zero
    index=0
    do iatom1=1,natom
     do idir1=1,3
      do iatom2=1,natom
       dred(:) = xred(:,iatom1) - xred(:,iatom2)
       dot = qphon(1,iqpt)*dred(1)+qphon(2,iqpt)*dred(2)+qphon(3,iqpt)*dred(3)
       index2=index2+1

!      Compute factor for SE term
       if(phfreq(imod,iqpt)<tol8)then
        factr = zero
        facti = zero
       else
        factr=cos(two_pi*dot) /sqrt(amu(typat(iatom1))*amu(typat(iatom2)))/phfreq(imod,iqpt)/amu_emass
        facti=sin(two_pi*dot) /sqrt(amu(typat(iatom1))*amu(typat(iatom2)))/phfreq(imod,iqpt)/amu_emass
       end if
       
       do idir2=1,3
        index = index+1

!       Compute products of polarization vectors
        vecr = eigvec(1,idir1,iatom1,imod,iqpt)*eigvec(1,idir2,iatom2,imod,iqpt)+&
&        eigvec(2,idir1,iatom1,imod,iqpt)*eigvec(2,idir2,iatom2,imod,iqpt)
        veci = eigvec(2,idir1,iatom1,imod,iqpt)*eigvec(1,idir2,iatom2,imod,iqpt)-&
&        eigvec(1,idir1,iatom1,imod,iqpt)*eigvec(2,idir2,iatom2,imod,iqpt)

        vec1r = eigvec(1,idir1,iatom1,imod,iqpt)*eigvec(1,idir2,iatom1,imod,iqpt)+&
&        eigvec(2,idir1,iatom1,imod,iqpt)*eigvec(2,idir2,iatom1,imod,iqpt)
        vec1i = eigvec(2,idir1,iatom1,imod,iqpt)*eigvec(1,idir2,iatom1,imod,iqpt)-&
&        eigvec(1,idir1,iatom1,imod,iqpt)*eigvec(2,idir2,iatom1,imod,iqpt)

        vec2r = eigvec(1,idir1,iatom2,imod,iqpt)*eigvec(1,idir2,iatom2,imod,iqpt)+&
&        eigvec(2,idir1,iatom2,imod,iqpt)*eigvec(2,idir2,iatom2,imod,iqpt)
        vec2i = eigvec(2,idir1,iatom2,imod,iqpt)*eigvec(1,idir2,iatom2,imod,iqpt)-&
&        eigvec(1,idir1,iatom2,imod,iqpt)*eigvec(2,idir2,iatom2,imod,iqpt)
        
!       Compute factor for DW term
        if(phfreq(imod,iqpt)<tol8)then
         fact2r = zero
         fact2i = zero
        else
         fact2r = wghtq(iqpt)*(vec1r/amu(typat(iatom1)) + vec2r/amu(typat(iatom2)))/phfreq(imod,iqpt)/&
&         amu_emass/2/norm(idir1)/norm(idir2)
         fact2i = wghtq(iqpt)*(vec1i/amu(typat(iatom1)) + vec2i/amu(typat(iatom2)))/phfreq(imod,iqpt)/&
&         amu_emass/2/norm(idir1)/norm(idir2)
        end if

        multr(:,:) =(blkval2(1,index,:,:,iqpt)*vecr - blkval2(2,index,:,:,iqpt)*veci)/(norm(idir1)*norm(idir2))
        multi(:,:) =(blkval2(1,index,:,:,iqpt)*veci + blkval2(2,index,:,:,iqpt)*vecr)/(norm(idir1)*norm(idir2))

!       Debye-Waller Term
        dwtermr(1:mband,1:nkpt) = dwtermr(1:mband,1:nkpt) + fact2r*blkval2(1,index,:,:,gqpt) - fact2i*blkval2(2,index,:,:,gqpt)
        dwtermi(1:mband,1:nkpt) = dwtermi(1:mband,1:nkpt) + fact2r*blkval2(2,index,:,:,gqpt) + fact2i*blkval2(1,index,:,:,gqpt)

!       Self-energy Term (Fan)
        deigr(1:mband,1:nkpt) = deigr(1:mband,1:nkpt) + wghtq(iqpt)*(factr*multr(1:mband,1:nkpt) - facti*multi(1:mband,1:nkpt))
        deigi(1:mband,1:nkpt) = deigi(1:mband,1:nkpt) + wghtq(iqpt)*(factr*multi(1:mband,1:nkpt) + facti*multr(1:mband,1:nkpt))

       end do !idir2
      end do !iatom2
     end do !idir1
    end do !iatom1
!   Eigenvalue derivative
    dednr(1:mband,1:nkpt) = deigr(1:mband,1:nkpt) - dwtermr(1:mband,1:nkpt)
    dedni(1:mband,1:nkpt) = deigi(1:mband,1:nkpt) - dwtermi(1:mband,1:nkpt)
    
!   Bose-Einstein distribution 
    if(phfreq(imod,iqpt)<tol6)then
     bosein = zero
    else
     bosein = one/(exp(phfreq(imod,iqpt)/(kb_HaK*tmp))-1) 
    end if

!   Calculate total
    thmeigen(1,1:mband,1:nkpt) = thmeigen(1,1:mband,1:nkpt) + dednr(1:mband,1:nkpt)*(bosein+half)
    thmeigen(2,1:mband,1:nkpt) = thmeigen(2,1:mband,1:nkpt) + dedni(1:mband,1:nkpt)*(bosein+half)

    if(itemper==1)then
!    Calculate slope of linear regime
     if(phfreq(imod,iqpt)<tol8)then
      slope(1,1:mband,1:nkpt) = slope(1,1:mband,1:nkpt) 
      slope(2,1:mband,1:nkpt) = slope(2,1:mband,1:nkpt) 
     else
      slope(1,1:mband,1:nkpt) = slope(1,1:mband,1:nkpt) + dednr(1:mband,1:nkpt)*(kb_HaK/phfreq(imod,iqpt))
      slope(2,1:mband,1:nkpt) = slope(2,1:mband,1:nkpt) + dedni(1:mband,1:nkpt)*(kb_HaK/phfreq(imod,iqpt))
     end if
!    Calculate zero-point renormalization
     zeropoint(1,1:mband,1:nkpt) = zeropoint(1,1:mband,1:nkpt) + dednr(1:mband,1:nkpt)*half
     zeropoint(2,1:mband,1:nkpt) = zeropoint(2,1:mband,1:nkpt) + dedni(1:mband,1:nkpt)*half
    end if
   end do ! imod
  end do !iqpt

! Write temperature independent results
  if(itemper==1)then
   write(unitout,'(a)')'Temperature independent results'
   do ikpt=1,nkpt
    write(unitout,'(a,3es16.8)')' Kpt :', kpnt(:,ikpt,1)
    do iband=1,mband
     write(unitout,'(4d22.14)') Ha_eV*slope(1,iband,ikpt),Ha_eV*slope(2,iband,ikpt),&
&     Ha_eV*zeropoint(1,iband,ikpt),Ha_eV*zeropoint(2,iband,ikpt)
    end do
   end do
   write(unitout,'(a)')'Temperature dependent results'
  end if
! Write result in file for each temperature
  write(unitout,'(a,es9.3,a)')'T :', tmp,' K'
  do ikpt=1,nkpt
   write(unitout,'(a,3es16.8)')' Kpt :', kpnt(:,ikpt,1)
   do iband=1,mband
    write(unitout,'(2d22.14)') Ha_eV*thmeigen(1,iband,ikpt), Ha_eV*thmeigen(2,iband,ikpt)
   end do
  end do
 end do !itemper


end subroutine thmeig
!!***

