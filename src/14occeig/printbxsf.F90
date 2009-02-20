!{\src2tex{textfont=tt}}
!!****f* ABINIT/printbxsf
!! NAME
!! printbxsf
!!
!! FUNCTION
!!  Print the Fermi surface in XCrysDen format
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2008 ABINIT group (MVerstraete,MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  eigen(mband,nkpt,nsppol) = eigenvalues in hartree
!!  ewind = energy window around the fermi level.
!!          if ewind /= 0 ==> a band is considered in the plot of FSurf
!!                            only if it is inside [ ef-ewind, ef+ewind ] for some k point
!!          if ewind == 0 ==> all bands will be keept in the _BXSF file
!!  fermie = Fermi energy (Hartree)
!!  gprimd(3,3) = dimensional primitive translations for reciprocal space (bohr^-1)
!!  kptrlatt(3,3) = reciprocal of lattice vectors for full kpoint grid
!!  mband = maximum number of bands
!!  nsppol = 1 for unpolarized, 2 for spin-polarized
!!  shiftk(3,nshiftk) =shift vector for k point grid
!!  fname = filename for the fortran file
!!
!! OUTPUT
!!  Only write
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      elphon,outscfcv
!!
!! CHILDREN
!!      canon9,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine printbxsf(eigen,ewind,fermie,gprimd,kptrlatt,mband,&
& nkptirred,kptirred,nsym,symrec,timrev,nsppol,shiftk,nshiftk,fname)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays,
!scalars
 integer,intent(in) :: mband,nkptirred,nshiftk,nsppol,nsym,timrev
 real(dp),intent(in) :: ewind,fermie
 character(len=fnlen) :: fname
!arrays
 integer,intent(in) :: kptrlatt(3,3),symrec(3,3,nsym)
 real(dp),intent(in) :: eigen(mband,nkptirred,nsppol),gprimd(3,3)
 real(dp),intent(in) :: kptirred(3,nkptirred),shiftk(3,nshiftk)

!Local variables-------------------------------
!scalars
 integer :: fform,found,iband,ierr,ik1,ik2,ik3,ikgrid,ikpt,ikpt1,indx,ios,iost
 integer :: isppol,isym,itim,maxband,minband,nk1,nk2,nk3,nkptfull,ubxsf
 real(dp) :: ene,mkval,res,ss,timsign
 character(len=500) :: message
!arrays
 integer,allocatable :: fulltoirred(:)
 real(dp) :: kconv(3),kpt(3),kptgrid(3),kptsym(3)

! *************************************************************************

!DEBUG
!write (message,'(3a,9i3)')' printbxsf : enter ',ch10,&
!& ' kptrlatt = ',kptrlatt(:,:)
!ENDDEBUG

!Error if klatt is no simple orthogonal lattice (in red space)
!for generalization to MP grids, need new version of XCrysDen

 if (     kptrlatt(1,2) /= 0 .or. kptrlatt(1,3) /= 0 .or. kptrlatt(2,1) /= 0       &
& .or. kptrlatt(2,3) /= 0 .or. kptrlatt(3,1) /= 0 .or. kptrlatt(3,2) /= 0 ) then
  write(message,'(5a)')' printbxsf : ERROR- ',ch10,     &
&  ' kptrlatt should be diagonal, for the FS calculation ',ch10,&
  ' action: use an orthogonal k-grid for the GS calculation '
  call wrtout(06,message,'COLL')
! call leave_new('COLL')
  return
 end if

 if (any(abs(shiftk(:,:)) > tol10 )) then
  write(message,'(5a)')' printbxsf : ERROR- ',ch10,        &
&  ' origin of the k grid should be (0,0,0) for the FS calculation ',ch10,&
  ' action: use an unshifted k-grid for the GS calculation '
  call wrtout(06,message,'COLL')
! call leave_new('COLL')
  return
 end if

!NOTE : Xcrysden uses aperiodical data-grid
 nk1 = kptrlatt(1,1)
 nk2=  kptrlatt(2,2)
 nk3=  kptrlatt(3,3)
 nkptfull=(nk1+1)*(nk2+1)*(nk3+1)

 allocate (fulltoirred(nkptfull),stat=ierr)
 if (ierr /= 0 ) then
  write (message,'(3a)')' printbxsf : ERROR- ',ch10,&
&  ' trying to allocate array fulltoirred '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!NOTE : Xcrysden uses the C-style for the Fermi Surface
 ikgrid=0
 do ik1=0,nk1
  do ik2=0,nk2
   do ik3=0,nk3

    ikgrid=ikgrid+1
    kptgrid(1)=1.*ik1/kptrlatt(1,1)
    kptgrid(2)=1.*ik2/kptrlatt(2,2)
    kptgrid(3)=1.*ik3/kptrlatt(3,3)
    call canon9(kptgrid(1),kpt(1),res)
    call canon9(kptgrid(2),kpt(2),res)
    call canon9(kptgrid(3),kpt(3),res)

!   find correspondence between the Xcrysden grid and the IBZ
!   NOTE1: this is the most heavy part of the code (doesnt scale well with nkptfull)
!   NOTE2: we can also use elph_ds%fulltoirred but I am too lazy to modify the code
    found=0
    irred: do ikpt1=1,nkptirred
     do isym=1,nsym
      do itim=0,timrev
       timsign = one-two*itim
       kptsym(:) = timsign*(symrec(:,1,isym)*kptirred(1,ikpt1) + &
&       symrec(:,2,isym)*kptirred(2,ikpt1) + &
&       symrec(:,3,isym)*kptirred(3,ikpt1))
       call canon9(kptsym(1),kconv(1),res)
       call canon9(kptsym(2),kconv(2),res)
       call canon9(kptsym(3),kconv(3),res)
!      is kconv equal to kpt?
       ss=  (kpt(1)-kconv(1))**2 +  &
&       (kpt(2)-kconv(2))**2 +  &
&       (kpt(3)-kconv(3))**2
       if (ss < tol6) then
        found=1
        fulltoirred(ikgrid)=ikpt1
        exit irred
       end if

      end do !itim
     end do !isym
    end do irred

    if (found /= 1) then
     write(message,'(3a,3es16.8,2a)')' printbxsf : ERROR- ',ch10,           &
&     ' kpt = ',kpt(:),ch10,                                                 &
&     ' has no symmetric among the irred kpoints used in the GS calculation '
     call wrtout(06,message,'COLL')
     call leave_new('COLL')
    end if

   end do !ik1
  end do !ik2
 end do !ik3

 if (abs(ewind) < tol12 ) then !keep all bands
  minband=1
  maxband=mband
 else
! if Fermi surface crosses band, include in plot of FSurf
! for all bands and 2 sppols (set as different bands)
  minband = mband
  maxband = 0
  ene=abs(ewind)
  do isppol=1,nsppol
   do iband=1,mband
    if(minval(eigen(iband,:,isppol))-fermie < -ene) then
     minband = iband
    end if
   end do
   do iband=mband,1,-1
    if(maxval(eigen(iband,:,isppol))-fermie > ene) then
     maxband = iband
    end if
   end do
  end do !end isppol

 end if !end if abs(energy_window)

 ubxsf=66
 open(unit=ubxsf,file=fname,status='unknown',form='formatted',iostat=iost)
 if (iost/=0) then
  write (message,'(2a)')' printbxsf : ERROR- opening file ',trim(fname)
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!write header
 write(ubxsf,*)' BEGIN_INFO'
 write(ubxsf,*)'   #'
 write(ubxsf,*)'   # this is a Band-XCRYSDEN-Structure-File for Visualization of Fermi Surface'
 write(ubxsf,*)'   # generated by the ABINIT package'
 write(ubxsf,*)'   #'
 write(ubxsf,*)'   #  bands between ',minband,' and ',maxband
 write(ubxsf,*)'   #'
 if (nsppol == 2 ) then
  write(ubxsf,*)'   # NOTE: the first band is relative to spin-up electrons,'
  write(ubxsf,*)'   # the second band to spin-down and so on .. '
  write(ubxsf,*)'   #'
 end if
 write(ubxsf,*)'   # Launch as: xcrysden --bxsf '
 write(ubxsf,*)'   #'
 write(ubxsf,'(a,es16.8)')'   Fermi Energy: ',fermie
 write(ubxsf,*)' END_INFO'
 write(ubxsf,*)' '
 write(ubxsf,*)' BEGIN_BLOCK_BANDGRID_3D'
 write(ubxsf,*)' band_energies'
 write(ubxsf,*)' BEGIN_BANDGRID_3D'

 write(ubxsf,*)' ',(maxband-minband+1)*nsppol
 write(ubxsf,*)' ',nk1+1,nk2+1,nk3+1
 write(ubxsf,*)' ',shiftk(:,1)
!NOTE : Angstrom units are used in the BXSF format
 write(ubxsf,*)' ',gprimd(:,1)/Bohr_Ang
 write(ubxsf,*)' ',gprimd(:,2)/Bohr_Ang
 write(ubxsf,*)' ',gprimd(:,3)/Bohr_Ang

!DEBUG
!write(*,*) 'ikpt, fulltoirred(ikpt)'
!do indx=1,nkptfull
!write(*,*) indx, fulltoirred(indx)
!end do
!ENDDEBUG

!print out data for all relevant bands and full kpt grid (redundant, yes)
!for each kpt in full zone, find equivalent irred kpt and print eigenval
 indx=0
 do iband=minband,maxband
  do isppol=1,nsppol
   write(ubxsf,*)' BAND: ',indx+minband
   write(ubxsf,'(7(es16.8))')(eigen(iband,fulltoirred(ikpt),isppol),&
&   ikpt=1,nkptfull)
   indx=indx+1
  end do
 end do

 write(ubxsf,*)'  END_BANDGRID_3D'
 write(ubxsf,*)' END_BLOCK_BANDGRID_3D'

 close (ubxsf)
 deallocate (fulltoirred)

!DEBUG
!write(std_out,*)' printbxsf : exit'
!stop
!ENDDEBUG

end subroutine printbxsf
!!***
