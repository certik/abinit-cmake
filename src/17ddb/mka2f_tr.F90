!{\src2tex{textfont=tt}}
!!****f* ABINIT/mka2f_tr
!!
!! NAME
!! mka2f_tr
!!
!! FUNCTION
!!  calculates the FS averaged Transport alpha^2F_tr alpha^2F_trout alpha^2F_trin functions
!!  calculates and outputs the associated electrical and thermal conductivities
!!  for the first task : copied from mka2F
!!
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (JPC)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYINGS=
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  elph_ds
!!    elph_ds%gkk2 = gkk2 matrix elements on full FS grid for each phonon mode
!!    elph_ds%nbranch = number of phonon branches = 3*natom
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%nFSkpt = number of kpts included in the FS integration
!!    elph_ds%phfrq = phonon frequencies
!!  FSintweight = integration weights on the FS
!!  FSirredtofull = indices of irreducible FS kpoints in full FSkpt array
!!  FSirredwtk = weights of irreducible kpoints
!!  FSkpt = coordinates of all FS kpoints
!!  gprim = reciprocal lattice vectors (maybe dimensioned...)
!!  gprimd = reciprocal lattice vectors (dimensionful)
!!  mustar = coulomb pseudopotential parameter
!!  n0 = DOS at the Fermi level calculated from the FSkpt integration weights
!!       eventually for 2 spin channels
!!  ucvol = Unit cell volume 
!!  natom = number of atoms
!!  nrpt = number of real-space points for FT interpolation
!!  phon_ds = datastructure with interatomic force constants to interpolate
!!     phonons
!!  rpt = coordinates of real-space points for FT interpolation
!!  wghatm = weights for real-space points for FT interpolation
!!
!! OUTPUT
!!  elph_ds
!!    
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      ftgam,inpphon,leave_new,simpson_int,wrtout,zgemm
!!
!! NOTES
!!   copied from ftiaf9.f
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mka2f_tr(elph_ds,FSkpt,gprim,n0,ucvol,natom,nrpt,phon_ds,rpt,wghatm,elph_tr_ds)

 use defs_basis
  use defs_datatypes
  use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_14occeig
 use interfaces_17ddb, except_this_one => mka2f_tr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrpt
 real(dp),intent(in) :: ucvol
 type(elph_tr_type) :: elph_tr_ds
 type(elph_type),intent(inout) :: elph_ds
 type(phon_type),intent(inout) :: phon_ds
!arrays
 real(dp),intent(in) :: gprim(3,3),n0(elph_ds%nsppol),rpt(3,nrpt)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 real(dp),intent(inout) :: FSkpt(3,elph_ds%nFSkpt)

!Local variables -------------------------
!x =w/(2kbT)
!scalars
 integer :: i1,iFSqpt,ibranch,iomega,iost,isppol,jbranch,nerr,unit_a2f_tr
 integer :: unit_lor,unit_rho,unit_therm
 real(dp) :: Temp,a2fprefactor,chgu,diagerr,firh,firhT,gaussfactor
 real(dp) :: gaussprefactor,gaussval,lambda_tr,lor0,lorentz,maxerr,maxx,omega
 real(dp) :: rho,tolexp,wtherm,xtr,xx
 character(len=500) :: message
 character(len=fnlen) :: fname
!arrays
 real(dp),parameter :: c0(2)=(/0.d0,0.d0/),c1(2)=(/1.d0,0.d0/)
 real(dp) :: displ(2,elph_ds%nbranch,elph_ds%nbranch),eigval(elph_ds%nbranch)
 real(dp) :: gam_now(2,elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: gammafact(elph_ds%nbranch),pheigval(elph_ds%nbranch)
 real(dp) :: pheigvec(2*elph_ds%nbranch*elph_ds%nbranch),phfrq(elph_ds%nbranch)
 real(dp) :: rho_T(1000),tmpa2f(elph_ds%na2f)
 real(dp) :: tmpgam1(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmpgam2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),allocatable :: integrho(:),matrx(:,:),tointegrho(:),zhpev1(:,:)
 real(dp),allocatable :: zhpev2(:)

! *********************************************************************
!calculate a2f_tr for frequencies between 0 and elph_ds%omega_max


#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zhpev
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zgemm
#endif

 write (*,*) 'mka2f_tr : enter '

 allocate (elph_tr_ds%a2f_1d_tr(elph_ds%na2f,elph_ds%nsppol))
 allocate (elph_tr_ds%a2f_1d_trin(elph_ds%na2f,elph_ds%nsppol))
 allocate (elph_tr_ds%a2f_1d_trout(elph_ds%na2f,elph_ds%nsppol))

!tolerance on gaussian being = 0
 tolexp = 1.d-100
 maxx = sqrt(-log(tolexp))

!! smearing width for deltas in the expression for alpha^2 F
!!  now set as input parameter in elphon
!elph_ds%a2fsmear = 0.00002_dp

!maximum value of frequency (a grid has to be chosen for the representation of alpha^2 F)
!WARNING! supposes this value has been set in mkelph_linwid.

!crc already set
!elph_ds%domega = (elph_ds%omega_max-elph_ds%omega_min)/(elph_ds%na2f-one)

 gaussprefactor = sqrt(piinv) / elph_ds%a2fsmear
 gaussfactor = one / elph_ds%a2fsmear


 allocate(matrx(2,(3*natom*(3*natom+1))/2))
 allocate(zhpev1(2,2*3*natom-1),zhpev2(3*3*natom-2))



 gammafact(:)=one
!ENDMG

 elph_tr_ds%a2f_1d_tr(:,:) = zero


 maxerr=0.
 nerr=0

 elph_tr_ds%a2f_1d_trin(:,:) = zero

 do isppol=1,elph_ds%nsppol

! loop over qpoint in full kpt grid (presumably dense)
  do iFSqpt=1,elph_ds%nFSkpt

!  DEBUG
!  write (*,*) 'mka2f_tr : iFSqpt = ', iFSqpt, ' / ', elph_ds%nFSkpt
!  ENDDEBUG

!  This reduced version of ftgkk supposes the kpoints have been integrated
!  in integrate_gamma. Do FT from real-space gamma grid to 1 qpt.
   call ftgam(wghatm,gam_now,elph_tr_ds%gamma_rpt_trin(:,:,isppol,:),gprim,natom,1,nrpt,0,rpt,&
&   FSkpt(:,iFSqpt))

!  MJV NOTE: this is done twice for both isppol. Unnecessary but cheap
   call inpphon(displ,pheigval,pheigvec,phfrq,phon_ds,FSkpt(:,iFSqpt))

!  Diagonalize gamma matrix at qpoint (complex matrix).

!  if doscalprod==0 we have to dot in the displacement vectors here
   if (elph_ds%doscalprod==0) then

    write(6,*)'doscalprod==0 in mka2f_tr is not coded: stop'
    stop
   else if (elph_ds%doscalprod == 1) then

    if (elph_ds%tkeepbands==0) then
     gammafact(:)= one
    else if (elph_ds%tkeepbands==1) then
!    gammafact(:)= two*pi*abs(phfrq(:))
     gammafact(:)= one
    else
     write (message,'(3a,i7)')' mka2f_tr : BUG- ',ch10,&
&     ' elph_ds%tkeepbands should be 0 or 1, while it is ',elph_ds%tkeepbands
     call wrtout(06,message,'COLL')
     call leave_new('COLL')
    end if

!   
!   NOTE: in these calls gam_now and pheigvec do not have the right rank, but blas usually does not care
!   
    call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, c1, gam_now, 3*natom,&
&    pheigvec, 3*natom, c0, tmpgam1, 3*natom)
    call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, c1, pheigvec, 3*natom,&
&    tmpgam1, 3*natom, c0, tmpgam2, 3*natom)
    diagerr = zero
    do ibranch=1,elph_ds%nbranch
     eigval(ibranch) = tmpgam2(1,ibranch,ibranch)
     do jbranch=1,ibranch-1
      diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))
     end do
     do jbranch=ibranch+1,elph_ds%nbranch
      diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))
     end do
    end do

    if (diagerr > tol12) then
     nerr=nerr+1
     maxerr=max(diagerr, maxerr)
!    write (*,*) 'mka2f_tr: error in diagonalization of gamma_tr with phon eigenvectors: ', diagerr
    end if

   else

    write (message,'(3a,i4)')' mka2f_tr: BUG-',ch10,&
&    ' Wrong value for elph_ds%doscalprod = ',elph_ds%doscalprod
    call wrtout(06,message,'COLL')
    call leave_new('COLL')

   end if
!  end doscalprod if



!  Add all contributions from the phonon modes at this qpoint to
!  a2f and the phonon dos.
   do ibranch=1,elph_ds%nbranch
!   write (*,*) 'mka2f_tr : ibranch,iFSqpt = ', ibranch,iFSqpt

!   if (phfrq(ibranch) < 0.00005) cycle
    if (abs(phfrq(ibranch)) < tol10) then
     a2fprefactor= zero
    else
     a2fprefactor = gammafact(ibranch)*eigval(ibranch)/(two_pi*abs(phfrq(ibranch))*n0(isppol))
    end if



    omega = elph_ds%omega_min
    tmpa2f(:) = zero
    do iomega=1,elph_ds%na2f
     xx = (omega-phfrq(ibranch))*gaussfactor
     gaussval = gaussprefactor*exp(-xx*xx)

     tmpa2f(iomega) = tmpa2f(iomega) + gaussval*a2fprefactor

     omega = omega+elph_ds%domega
    end do

    elph_tr_ds%a2f_1d_trin(:,isppol) = elph_tr_ds%a2f_1d_trin(:,isppol) + tmpa2f(:)

   end do
!  ! end ibranch do

  end do ! end iFSqpt do
  write (*,*) 'mka2f_tr: errors in diagonalization of gamma_tr with phon eigenvectors: ', nerr,maxerr
  
 end do ! end isppol

!write (*,*) 'a2f_1d = ', a2f_1d

!second 1 / elph_ds%nFSkpt factor for the integration weights
 elph_tr_ds%a2f_1d_trin(:,:) = elph_tr_ds%a2f_1d_trin(:,:) / elph_ds%nFSkpt
!write(6,*)'  elph_tr_ds%a2f_1d_trin', elph_tr_ds%a2f_1d_trin


 nerr=0
 maxerr=0.

!same for TROUT


 elph_tr_ds%a2f_1d_trout(:,:) = zero

 do isppol=1,elph_ds%nsppol

! loop over qpoint in full kpt grid (presumably dense)
  do iFSqpt=1,elph_ds%nFSkpt

!  DEBUG
!  write (*,*) 'mka2f_tr : iFSqpt = ', iFSqpt, ' / ', elph_ds%nFSkpt
!  ENDDEBUG

!  This reduced version of ftgkk supposes the kpoints have been integrated
!  in integrate_gamma. Do FT from real-space gamma grid to 1 qpt.
   call ftgam(wghatm,gam_now,elph_tr_ds%gamma_rpt_trout(:,:,isppol,:),gprim,natom,1,nrpt,0,rpt,&
&   FSkpt(:,iFSqpt))

!  NOTE MJV: this is now doubled for isppol AND trin/trout. Could all this
!  be merged?
   call inpphon(displ,pheigval,pheigvec,phfrq,phon_ds,FSkpt(:,iFSqpt))

!  Diagonalize gamma matrix at qpoint (complex matrix).

!  if doscalprod==0 we have to dot in the displacement vectors here
   if (elph_ds%doscalprod==0) then

    write(6,*)'not coded stop'
    stop
   else if (elph_ds%doscalprod == 1) then

    if (elph_ds%tkeepbands==0) then
     gammafact(:)= one
    else if (elph_ds%tkeepbands==1) then
!    gammafact(:)= two*pi*abs(phfrq(:))
     gammafact(:)= one
    else
     write (message,'(3a,i7)')' mka2f_tr : BUG- ',ch10,&
&     ' elph_ds%tkeepbands should be 0 or 1, while it is ',elph_ds%tkeepbands
     call wrtout(06,message,'COLL')
     call leave_new('COLL')
    end if

    call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, c1, gam_now, 3*natom,&
&    pheigvec, 3*natom, c0, tmpgam1, 3*natom)
    call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, c1, pheigvec, 3*natom,&
&    tmpgam1, 3*natom, c0, tmpgam2, 3*natom)

!   DEBUG
!   write (130,*) '# mka2f_tr : gamma diagonalized with phonon eigenvectors '
!   write (130,'(3(2(E20.10,1x)))') tmpgam2
!   ENDDEBUG

    diagerr = zero
    do ibranch=1,elph_ds%nbranch
     eigval(ibranch) = tmpgam2(1,ibranch,ibranch)
     do jbranch=1,ibranch-1
      diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))
     end do
     do jbranch=ibranch+1,elph_ds%nbranch
      diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))
     end do
    end do

    if (diagerr > tol12) then
     nerr=nerr+1
     maxerr=max(diagerr, maxerr)
!    write (*,*) 'mka2f_tr: error in diagonalization of gamma_trout with phon eigenvectors: ', diagerr
    end if

   else

    write (message,'(3a,i4)')' mka2f_tr: BUG-',ch10,&
&    ' Wrong value for elph_ds%doscalprod = ',elph_ds%doscalprod
    call wrtout(06,message,'COLL')
    call leave_new('COLL')

   end if
!  end doscalprod if



!  Add all contributions from the phonon modes at this qpoint to
!  a2f and the phonon dos.
   do ibranch=1,elph_ds%nbranch
!   write (*,*) 'mka2f_tr : ibranch,iFSqpt = ', ibranch,iFSqpt

!   if (phfrq(ibranch) < 0.00005) cycle
    if (abs(phfrq(ibranch)) < tol10) then
     a2fprefactor= zero
    else
     a2fprefactor = gammafact(ibranch)*eigval(ibranch)/(two_pi*abs(phfrq(ibranch))*n0(isppol))
    end if
!   ENDMG
!   if (a2fprefactor < tol12) cycle
!   write (*,*) ' a2fp,eigv,phf ',a2fprefactor,eigval(ibranch),phfrq(ibranch)

!   MG  Added gammafact to treat the doscalprod=1 case
!   Added avglambda and avgomglog to calculate lamda and omega_log using the sum over the kpt-grid.
!   If the k-grid is dense enough, these values should be better than the corresponding quantities
!   evaluated through the integration over omega that depends on the elph_ds%a2fsmear

!   ENDMG

!   omega = zero
    omega = elph_ds%omega_min
    tmpa2f(:) = zero
    do iomega=1,elph_ds%na2f
     xx = (omega-phfrq(ibranch))*gaussfactor
     gaussval = gaussprefactor*exp(-xx*xx)

     tmpa2f(iomega) = tmpa2f(iomega) + gaussval*a2fprefactor

     omega = omega+elph_ds%domega
    end do

    elph_tr_ds%a2f_1d_trout(:,isppol) = elph_tr_ds%a2f_1d_trout(:,isppol) + tmpa2f(:)

   end do
!  ! end ibranch do

  end do ! end iFSqpt do
 end do ! end isppol

!write (*,*) 'a2f_1d = ', a2f_1d

!second 1 / elph_ds%nFSkpt factor for the integration weights
 elph_tr_ds%a2f_1d_trout(:,:) = elph_tr_ds%a2f_1d_trout(:,:) / elph_ds%nFSkpt

 write (*,*) 'mka2f_tr: errors in diagonalization of gamma_tr with phon eigenvectors: ', nerr,maxerr

 elph_tr_ds%a2f_1d_tr(:,:) =   elph_tr_ds%a2f_1d_trout(:,:) -  elph_tr_ds%a2f_1d_trin(:,:)

 deallocate(matrx)
 deallocate(zhpev1)
 deallocate(zhpev2)

!write (*,'(6(E16.6,2x))') a2f_1d(:,:)

!output the elph_tr_ds%a2f_1d_tr
 unit_a2f_tr = 108
 fname = trim(elph_ds%elph_base_name) // '_A2F_TR'
 open (unit=unit_a2f_tr,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
  write (message,'(3a)')' mka2f_tr : ERROR- opening file ',trim(fname),' as new'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 write (*,*) ' a2f_TR function integrated over the FS'
 write (unit_a2f_tr,'(a)') '#'
 write (unit_a2f_tr,'(a)') '# ABINIT package : a2f_tr file'
 write (unit_a2f_tr,'(a)') '#'
 write (unit_a2f_tr,'(a)') '# a2f_tr function integrated over the FS. omega in a.u.'
 write (unit_a2f_tr,'(a,I10)') '#     number of kpoints integrated over : ', elph_ds%nFSkpt
 write (unit_a2f_tr,'(a,I10)') '#     number of energy points : ',elph_ds%na2f
 write (unit_a2f_tr,'(a,E16.6,a,E16.6,a)') '#       between omega_min = ', elph_ds%omega_min, &
& ' Ha and omega_max = ', elph_ds%omega_max, ' Ha'
 write (unit_a2f_tr,'(a,E16.6)') '#   and the smearing width for gaussians is ', elph_ds%a2fsmear
 write (unit_a2f_tr,'(a)') '#'

!done with header

 do isppol=1,elph_ds%nsppol
  write (unit_a2f_tr,'(a,E16.6)') '# The DOS at Fermi level is ', n0(isppol)
! omega = zero
  omega = elph_ds%omega_min
  do iomega=1,elph_ds%na2f
   write (*,*) omega, elph_tr_ds%a2f_1d_tr(iomega,isppol)
   write (unit_a2f_tr,'(4D16.6)') omega, elph_tr_ds%a2f_1d_tr(iomega,isppol),&
&   elph_tr_ds%a2f_1d_trout(iomega,isppol),&
&   elph_tr_ds%a2f_1d_trin(iomega,isppol)
   omega=omega+elph_ds%domega
  end do
  write (unit_a2f_tr,*)
 end do !isppol

 close (unit=unit_a2f_tr)

!calculation of transport properties
 allocate (integrho(elph_ds%na2f))  
 allocate (tointegrho(elph_ds%na2f))  

 unit_rho = 109
 fname = trim(elph_ds%elph_base_name) // '_RHO'
 open (unit=unit_rho,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
  write (message,'(3a)')' mka2f_tr : ERROR- opening file ',trim(fname),' as new'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
!print header to resistivity file
 write (unit_rho,*) '# Resistivity as a function of temperature.'
 write (unit_rho,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_rho,*) '#  '
 write (unit_rho,*) '#  Columns are: '
 write (unit_rho,*) '#  temperature[K]   rho[au]   rho [SI]        rho/temp [au]'
 write (unit_rho,*) '#  '

 unit_therm = 111
 fname = trim(elph_ds%elph_base_name) // '_WTH'
 open (unit=unit_therm,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
  write (message,'(3a)')' mka2f_tr : ERROR- opening file ',trim(fname),' as new'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!print header to thermal conductivity file
 write (unit_therm,*) '# Thermal conductivity/resistivity as a function of temperature.'
 write (unit_therm,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_therm,*) '#  '
 write (unit_therm,*) '#  Columns are: '
 write (unit_therm,*) '#  temperature[K]   themal rho[au]   thermal cond [au]   themal rho [SI]   thermal cond [SI]'
 write (unit_therm,*) '#  '

 unit_lor = 112
 fname = trim(elph_ds%elph_base_name) // '_LOR'
 open (unit=unit_lor,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
  write (message,'(3a)')' mka2f_tr : ERROR- opening file ',trim(fname),' as new'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!print header to lorentz file
 write (unit_lor,*) '# Lorentz number as a function of temperature.'
 write (unit_lor,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_lor,*) '#  '
 write (unit_lor,*) '#  Columns are: '
 write (unit_lor,*) '#  temperature[K]   Lorentz number[au]   Lorentz quantum = (pi*kb_HaK)**2/3'
 write (unit_lor,*) '#  '

 do isppol=1,elph_ds%nsppol
  write(888,*) '# tointegrho for isppol ', isppol
  omega = elph_ds%omega_min
  do iomega=1,elph_ds%na2f
   if(omega<=0) then
    omega=omega+elph_ds%domega
    cycle
   end if
   tointegrho(iomega)=two*elph_tr_ds%a2f_1d_tr(iomega,isppol)/omega
   write(888,*)omega, tointegrho(iomega)
   omega=omega+elph_ds%domega
  end do
  write(888,*)

  call simpson_int(elph_ds%na2f,elph_ds%domega,tointegrho,integrho)
  lambda_tr=integrho(elph_ds%na2f)
  write (message, '(a,i3,a,es16.6)' )&
&  '- mka2f_tr: TRANSPORT lambda for isppol ', isppol, ' =  ', lambda_tr
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')
 end do !end isppol do

!constant to change units of rho from au to SI
 chgu=2.173969*1.0d-7

 do isppol=1,elph_ds%nsppol

! prefactor for resistivity integral
  firh=6.*pi*ucvol*kb_HaK/(n0(isppol)*elph_tr_ds%FSelecveloc_sq(isppol))
  write(6,'(A,D19.7)')'factrho', firh*lambda_tr/2.

  write(unit_rho,*) '# Rho for isppol= ', isppol

  tointegrho(:)=0.
  do i1=1,1000  ! runs over termperature in K
   Temp=float(i1)
   firhT=firh*Temp
   omega = elph_ds%omega_min
   do iomega=1,elph_ds%na2f
    if(omega<=0) then
     omega=omega+elph_ds%domega
     cycle
    end if
    xtr=omega/(2*kb_HaK*Temp)
    tointegrho(iomega)=firhT*omega*elph_tr_ds%a2f_1d_tr(iomega,isppol)/(((2*Temp*kb_HaK)**2)*((exp(xtr)-exp(-xtr))/2)**2)
!   write(*,*)'omega tointegrho =',omega, tointegrho(iomega)
    omega=omega+elph_ds%domega
   end do

   call simpson_int(elph_ds%na2f,elph_ds%domega,tointegrho,integrho)
   rho=integrho(elph_ds%na2f)
   write(*,*)'TEMP RHO ',temp,rho
   write(unit_rho,'(4D17.10)')temp,rho,rho*chgu,rho/temp
   rho_T(i1)=rho
  end do ! temperature
  write(unit_rho,*)
 end do ! isppol

!-----------------------------


 do isppol=1,elph_ds%nsppol
! prefactor for integral of thermal conductivity
  firh=(18.*ucvol)/(pi*kb_HaK*n0(isppol)*elph_tr_ds%FSelecveloc_sq(isppol))

  write(unit_therm,*) '# Thermal resistivity for isppol= ', isppol
  write(unit_lor,*) '# Lorentz coefficient for isppol= ', isppol

! write(6,*)'ucvol',ucvol
! write(6,*)'pi',pi
! write(6,*)'kb',kb_HaK
! write(6,*)'n0',n0(isppol)
! write(6,*)'V2',elph_tr_ds%FSelecveloc_sq(isppol)
! write(6,'(A,D19.10)')'FACT1 WTH',firh
! write(6,'(A,D19.10)')'FACT2 WTH',firh*lambda_tr/2

  tointegrho(:)=0.
  do i1=1,1000

   Temp=float(i1)
!  write(*,*)'TEMP =',Temp
!  firhT=firh/Temp
   omega = elph_ds%omega_min
   do iomega=1,elph_ds%na2f
    if(omega<=0) then
     omega=omega+elph_ds%domega
     cycle
    end if
    xtr=omega/(2*kb_HaK*Temp)
!   factprt=xtr/(((exp(xtr)-exp(-xtr)/2))**2)
!   write(6,*)factprt
    tointegrho(iomega)=xtr**2/omega*&
&    ( elph_tr_ds%a2f_1d_tr(iomega,isppol)+&
&    4*xtr**2*elph_tr_ds%a2f_1d_trout(iomega,isppol)/pi**2+   &
&    2*xtr**2*elph_tr_ds%a2f_1d_trin(iomega,isppol)/pi**2)  &
&    /(((exp(xtr)-exp(-xtr))/2)**2)

!   write(*,*)'omega tointegrho =',omega, tointegrho(iomega)
    omega=omega+elph_ds%domega
   end do

   call simpson_int(elph_ds%na2f,elph_ds%domega,tointegrho,integrho)
!  write(6,*)'INTEGRALE WTH =',integrho(elph_ds%na2f)
   wtherm=integrho(elph_ds%na2f)*firh
   write(6,'(A,3D12.5)')'TEMP WTHERM ',temp,wtherm,1./wtherm

   write(unit_therm,'(5D12.5)')temp,wtherm,1./wtherm,wtherm/3.4057d9,1./(wtherm) *3.4057d9

   lorentz=rho_T(i1)/(wtherm*temp)
!  TODO: remove   MJV : this looks like a constant in a loop !
   lor0=(pi*kb_HaK)**2/3.
   write(unit_lor,*)temp,lorentz,lor0

  end do
  write(unit_therm,*)
  write(unit_lor,*)
 end do !end isppol do


 close (unit=unit_rho)
 close (unit=unit_therm)

 deallocate (integrho)
 deallocate (tointegrho)  
 deallocate (elph_tr_ds%a2f_1d_tr)
 deallocate (elph_tr_ds%a2f_1d_trin)
 deallocate (elph_tr_ds%a2f_1d_trout)
 
 
 deallocate (elph_tr_ds%gamma_qpt_trin)
 deallocate (elph_tr_ds%gamma_qpt_trout)
 deallocate (elph_tr_ds%gamma_rpt_trin)
 deallocate (elph_tr_ds%gamma_rpt_trout)
 write (*,*) ' mka2f_tr : end '


end subroutine mka2f_tr
!!***
