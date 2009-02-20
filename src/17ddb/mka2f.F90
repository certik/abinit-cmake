!{\src2tex{textfont=tt}}
!!****f* ABINIT/mka2f
!!
!! NAME
!! mka2f
!!
!! FUNCTION
!!  calculate the FS averaged alpha^2F function
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
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
!!  n0 = DOS at the Fermi level calculated from the FSkpt integration weights (event. 2 spin pol)
!!  natom = number of atoms
!!  nrpt = number of real-space points for FT interpolation
!!  phon_ds = datastructure with interatomic force constants to interpolate
!!     phonons
!!  rpt = coordinates of real-space points for FT interpolation
!!  wghatm = weights for real-space points for FT interpolation
!!
!! OUTPUT
!!  a2f_1d = 1D alpha
!!  dos_phon = density of states for phonons
!!  elph_ds
!!    elph_ds%a2f = alpha^2F function on full FS grid : omega,iFSkpt,iqpt
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

subroutine mka2f(a2f_1d,dos_phon,elph_ds,FSintweight,FSirredtofull,FSirredwtk,FSkpt,gprim,gprimd,mustar,&
&      n0,natom,nrpt,phon_ds,rpt,wghatm)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_14occeig
 use interfaces_17ddb, except_this_one => mka2f
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrpt
 real(dp),intent(in) :: mustar
 type(elph_type),intent(inout) :: elph_ds
 type(phon_type),intent(inout) :: phon_ds
!arrays
 integer,intent(in) :: FSirredtofull(elph_ds%nFSkptirred)
 real(dp),intent(in) :: FSintweight(elph_ds%nFSband,elph_ds%nFSkpt)
 real(dp),intent(in) :: FSirredwtk(elph_ds%nFSkptirred),gprim(3,3),gprimd(3,3)
 real(dp),intent(in) :: n0(elph_ds%nsppol),rpt(3,nrpt),wghatm(natom,natom,nrpt)
 real(dp),intent(inout) :: FSkpt(3,elph_ds%nFSkpt)
 real(dp),intent(out) :: a2f_1d(elph_ds%na2f),dos_phon(elph_ds%na2f)

!Local variables -------------------------
!scalars
 integer :: i1,i2,iFSqpt,iatom,ib1,ib2,ibranch,idir,ieqFSkpt1,ier,ii,iomega
 integer :: iost,ip,ipert1,ipert2,ipp,iqpt,irpt,isppol,jbranch,k1,kdir,mu,nu
 integer :: unit_a2f,unit_phdos
 real(dp) :: a2fprefactor,avgelphg,avglambda,avgomlog,diagerr,gaussfactor
 real(dp) :: gaussprefactor,gaussval,lambda_2,lambda_3,lambda_4,lambda_5
 real(dp) :: lambda_iso,lqn,maxx,omega,omegalog,omlog_qn,tc_macmill,tolexp,xx
 character(len=500) :: message
 character(len=fnlen) :: fname
!arrays
 real(dp),parameter :: c0(2)=(/0._dp,0._dp/),c1(2)=(/1._dp,0._dp/)
 real(dp) :: displ(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: eigval(elph_ds%nbranch),eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: gam_now(2,elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: gammafact(elph_ds%nbranch),imeigval(elph_ds%nbranch)
 real(dp) :: pheigval(elph_ds%nbranch)
 real(dp) :: pheigvec(2*elph_ds%nbranch*elph_ds%nbranch),phfrq(elph_ds%nbranch)
 real(dp) :: tmp_gkk2(elph_ds%nbranch,elph_ds%nFSkpt),tmpa2f(elph_ds%na2f)
 real(dp) :: tmpgam1(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmpgam2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmpphondos(elph_ds%na2f)
 real(dp),allocatable :: a2f1mom(:),a2f2mom(:),a2f3mom(:),a2f4mom(:)
 real(dp),allocatable :: a2f_1mom(:),a2f_1mom_int(:),a2flogmom(:)
 real(dp),allocatable :: a2flogmom_int(:),matrx(:,:),zhpev1(:,:),zhpev2(:)

! *********************************************************************
!calculate a2f for frequencies between 0 and elph_ds%omega_max

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zhpev
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zgemm
#endif

 write (*,*) 'mka2f : enter '

!tolerance on gaussian being = 0
 tolexp = 1.d-100
 maxx = sqrt(-log(tolexp))

!allocate (elph_ds%a2f(elph_ds%na2f,elph_ds%nFSkpt,elph_ds%nFSkpt))

!! smearing width for deltas in the expression for alpha^2 F
!!  now set as input parameter in elphon
!elph_ds%a2fsmear = 0.00002_dp

!maximum value of frequency (a grid has to be chosen for the representation of alpha^2 F)
!WARNING! supposes this value has been set in mkelph_linwid.
!elph_ds%omega_max = 0.002
 elph_ds%domega = (elph_ds%omega_max-elph_ds%omega_min)/(elph_ds%na2f-one)

 gaussprefactor = sqrt(piinv) / elph_ds%a2fsmear
 gaussfactor = one / elph_ds%a2fsmear

!DEBUG
!write (*,*) ' mka2f : FSirredwtk'
!write (*,'(6(E16.6))') FSirredwtk
!write (*,*) ' mka2f :sum(FSintweight) = ',sum(FSintweight)
!write (*,*) ' mka2f : omega_max,a2fsmear,na2f ', &
!& elph_ds%omega_max,elph_ds%a2fsmear,elph_ds%na2f
!write (*,'(6(E16.6))') gaussprefactor,gaussfactor,elph_ds%domega
!ENDDEBUG

 allocate(matrx(2,(3*natom*(3*natom+1))/2))
 allocate(zhpev1(2,2*3*natom-1),zhpev2(3*3*natom-2))

!
!output the a2f_1d header
!
 unit_a2f = 108
 fname = trim(elph_ds%elph_base_name) // '_A2F'
!
!only open the file for the first sppol
!
 open (unit=unit_a2f,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
  write (message,'(2a)')' mka2f : ERROR- opening file ',trim(fname)
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 write (*,*) ' a2f function integrated over the FS'
 write (unit_a2f,'(a)') '#'
 write (unit_a2f,'(a)') '# ABINIT package : a2f file'
 write (unit_a2f,'(a)') '#'
 write (unit_a2f,'(a)') '# a2f function integrated over the FS. omega in a.u.'
 write (unit_a2f,'(a,I10)') '#     number of kpoints integrated over : ', elph_ds%nFSkpt
 write (unit_a2f,'(a,I10)') '#     number of energy points : ',elph_ds%na2f
 write (unit_a2f,'(a,E16.6,a,E16.6,a)') '#       between omega_min = ', elph_ds%omega_min, &
& ' Ha and omega_max = ', elph_ds%omega_max, ' Ha'
 write (unit_a2f,'(a,E16.6)') '#   and the smearing width for gaussians is ', elph_ds%a2fsmear

!
!output the phonon DOS header
!
 unit_phdos = 109
 fname = trim(elph_ds%elph_base_name) // '_PDS'
 open (unit=unit_phdos,file=fname,status='replace',iostat=iost)
 if (iost /= 0) then
  write (message,'(3a)')' mka2f : ERROR- opening file ',trim(fname),' as new'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 
 write (*,*) 'phonon DOS integrated over the FS'
 write (unit_phdos,'(a)') '#'
 write (unit_phdos,'(a)') '# ABINIT package : phonon DOS file'
 write (unit_phdos,'(a)') '#'
 write (unit_phdos,'(a)') '# Phonon DOS integrated over the FS. omega in a.u. EXPERIMENTAL!!!'
 write (unit_phdos,'(a,I10)') '#     number of kpoints integrated over : ', elph_ds%nFSkpt
 write (unit_phdos,'(a,I10)') '#     number of energy points : ', elph_ds%na2f
 write (unit_phdos,'(a,E16.6,a,E16.6,a)') '#       between omega_min = ', elph_ds%omega_min, &
& ' Ha and omega_max = ', elph_ds%omega_max, ' Ha'
 write (unit_phdos,'(a,i4,a,E16.6)') '# The DOS at Fermi level for spin ', 1, ' is ', n0(1)
 if (elph_ds%nsppol==2) then
  write (unit_phdos,'(a,i4,a,E16.6)') '# The DOS at Fermi level for spin ', 2, ' is ', n0(2)
 end if
 write (unit_phdos,'(a,E16.6)') '#   and the smearing width for gaussians is ', elph_ds%a2fsmear
 write (unit_phdos,'(a)') '#'


!MG20060607 we need this factor to treat the doscalprod=1 case in the right way
 gammafact(:)=one
!ENDMG

!TODO add sppol loop here
 do isppol=1,elph_ds%nsppol

  write (*,*) '##############################################'
  write (*,*) 'mka2f : Treating spin polarization ', isppol
  write (*,*) '##############################################'

! Average of electron phonon coupling over the whole BZ
  avgelphg = zero
! MG20060607 Do the same for lambda and omega_log
  avglambda = zero
  avgomlog = zero

  a2f_1d(:) = zero
  dos_phon(:) = zero
! elph_ds%a2f(:,:,:) = zero

! loop over qpoint in full kpt grid (presumably dense)
  do iFSqpt=1,elph_ds%nFSkpt

!  DEBUG
!  write (*,*) 'mka2f : iFSqpt = ', iFSqpt, ' / ', elph_ds%nFSkpt
!  ENDDEBUG

!  This reduced version of ftgkk supposes the kpoints have been integrated
!  in integrate_gamma. Do FT from real-space gamma grid to 1 qpt.
!  TODO : call ftgam with only gamma_rpt for 1 spin channel
   call ftgam(wghatm,gam_now,elph_ds%gamma_rpt(:,:,isppol,:),gprim,natom,1,nrpt,0,rpt,&
&   FSkpt(:,iFSqpt))

   call inpphon(displ,pheigval,pheigvec,phfrq,phon_ds,FSkpt(:,iFSqpt))

!  Diagonalize gamma matrix at qpoint (complex matrix).

!  if doscalprod==0 we have to dot in the displacement vectors here
   if (elph_ds%doscalprod==0) then

    displ_red(:,:,:) = zero
    do jbranch=1,elph_ds%nbranch
     do iatom=1,natom
      do idir=1,3
       ibranch=idir+3*(iatom-1)
       do kdir=1,3
        k1 = kdir+3*(iatom-1)

        displ_red(1,ibranch,jbranch) = displ_red(1,ibranch,jbranch) + &
&        gprimd(kdir,idir)*displ(1,k1,jbranch)

        displ_red(2,ibranch,jbranch) = displ_red(2,ibranch,jbranch) + &
&        gprimd(kdir,idir)*displ(2,k1,jbranch)

       end do
      end do
     end do
    end do

!   DEBUG
!   write (*,*) 'mka2f : displ_red  = '
!   write (*,'(3(2E18.6,1x))') displ_red
!   write (97,'(9(2E14.6,x))')  displ(:,:,:)
!   write (96,'(9(2E14.6,x))')  displ_red(:,:,:)
!   write (98,'(9(2E14.6,x))')  gam_now(:,:)
!   ENDDEBUG

    eigval(:) = zero
    imeigval(:) = zero
    do jbranch=1,elph_ds%nbranch
     do ipert1=1,elph_ds%nbranch
      do ipert2=1,elph_ds%nbranch
       ipp=(ipert1-1)*elph_ds%nbranch + ipert2

!      calculate displ_red* gam_now* displ_red^{*T}
       eigval(jbranch) = eigval(jbranch)                                                  &
&       + displ_red(1,ipert2,jbranch)*gam_now(1,ipp)*displ_red(1,ipert1,jbranch)&
&       - displ_red(2,ipert2,jbranch)*gam_now(2,ipp)*displ_red(1,ipert1,jbranch)&
&       + displ_red(1,ipert2,jbranch)*gam_now(2,ipp)*displ_red(2,ipert1,jbranch)&
&       + displ_red(2,ipert2,jbranch)*gam_now(1,ipp)*displ_red(2,ipert1,jbranch)

!      could do away with this - should be 0
!      MG20060603   Fixed a minor bug present in version 5.1.3 (this quantity indeed should be zero and is not used!)
       imeigval(jbranch) = imeigval(jbranch)&
&       + displ_red(1,ipert2,jbranch)*gam_now(2,ipp)*displ_red(1,ipert1,jbranch)&
&       - displ_red(1,ipert2,jbranch)*gam_now(1,ipp)*displ_red(2,ipert1,jbranch)&
&       + displ_red(2,ipert2,jbranch)*gam_now(1,ipp)*displ_red(1,ipert1,jbranch)&
&       + displ_red(2,ipert2,jbranch)*gam_now(2,ipp)*displ_red(2,ipert1,jbranch)

      end do
     end do

     if (abs(imeigval(jbranch)) > tol8) then
      write (message,'(3a,i6,a,es16.8)')&
&      ' mka2f : WARNING-  imaginary values! ',ch10,&
&      ' branch = ',jbranch,' imeigval = ',imeigval(jbranch)
      call wrtout(06,message,'COLL')
     end if

    end do

!   if doscalprod==1 we have to diagonalize the matrix we interpolated.
   else if (elph_ds%doscalprod == 1) then

!   MG20060607 FIXME
!   I have modified nmsq_gam_sumfs and nmsq_gam to calculate the matrix
!   gkk_qpt using the same approach as in nmsq_pure_gkk, in this case gammafact is 1
!   this part can be removed after careful testing of the new implementation
    if (elph_ds%tkeepbands==0) then
     gammafact(:)= one
    else if (elph_ds%tkeepbands==1) then
!    gammafact(:)= two*pi*abs(phfrq(:))
     gammafact(:)= one
    else
     write (message,'(3a,i7)')' mka2f : BUG- ',ch10,&
&     ' elph_ds%tkeepbands should be 0 or 1, while it is ',elph_ds%tkeepbands
     call wrtout(06,message,'COLL')
     call leave_new('COLL')
    end if
!   ENDMG20060607

!   !  First method diagonalize the matrix explicitly
!   !  here copied from phfrq33
!   ier=0
!   ii=1
!   do i2=1,3*natom
!   do i1=1,i2
!   ipp=i1+(i2-1)*3*natom
!   if (ii > (3*natom*(3*natom+1))/2 ) then
!   write (*,*) 'PANIC', ii, i1, i2
!   stop
!   end if
!   matrx(1,ii)=gam_now(1,ipp)
!   matrx(2,ii)=gam_now(2,ipp)
!   ii=ii+1
!   end do
!   end do
!   !#if defined T3E
!   !     call CHPEV ('N','U',3*natom,matrx,eigval,eigvec,3*natom,zhpev1,&
!   !          &    zhpev2,ier)
!   !#else
!   call ZHPEV ('V','U',3*natom,matrx,eigval,eigvec,3*natom,zhpev1,&
!   &    zhpev2,ier)
!   !#endif
!   !
!   !   compare eigenvectors of phonon and gamma interpolations
!   !
!   write (130,*) '# mka2f : phonon eigenvectors at point ', FSkpt(:,iFSqpt)
!   write (130,'(3(2(E20.10,1x)))') pheigvec
!   write (130,*) '# mka2f : gamma eigenvectors '
!   write (130,'(3(2(E20.10,1x)))') eigvec

!   
!   Second method:
!   using phonon eigenvectors U, gamma_diag = U^{T*} gam_now U
!   
!   DEBUG
!   write (*,'(3(2(E20.10,1x)))') gam_now, pheigvec
!   write (*,'(3(2(E20.10,1x)))') CMPLX(one,zero), CMPLX(zero,zero)
!   write (*,*) 'N', 3*natom
!   ENDDEBUG

!   MJV NOTE : gam_now is being recast as a (3*natom)**2 matrix here
    call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, c1, gam_now, 3*natom,&
&    pheigvec, 3*natom, c0, tmpgam1, 3*natom)
    call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, c1, pheigvec, 3*natom,&
&    tmpgam1, 3*natom, c0, tmpgam2, 3*natom)

!   DEBUG
!   write (130,*) '# mka2f : gamma diagonalized with phonon eigenvectors '
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
     write (*,*) 'mka2f: error in diagonalization of gamma with phon eigenvectors: ', diagerr
    end if

   else

    write (message,'(3a,i4)')' mka2f: BUG-',ch10,&
&    ' Wrong value for elph_ds%doscalprod = ',elph_ds%doscalprod
    call wrtout(06,message,'COLL')
    call leave_new('COLL')

   end if
!  end doscalprod if

!  MG20060603MG
!  there was a bug in the calculation of the phonon DOS
!  since frequencies with small e-ph interaction were skipped inside the loop
!  In this new version all the frequencies (both positive and negative) are taken into account.
!  IDEA: it could be useful to calculate the PH-dos and the a2f
!  using several smearing values to perform a convergence study
!  Now the case doscalprod=1 is treated in the right way although it is not default anymore
!  FIXME to be checked
!  ENDMG

!  Add all contributions from the phonon modes at this qpoint to
!  a2f and the phonon dos.
   do ibranch=1,elph_ds%nbranch
!   write (*,*) 'mka2f : ibranch,iFSqpt = ', ibranch,iFSqpt

!   if (phfrq(ibranch) < 0.00005) cycle
    if (abs(phfrq(ibranch)) < tol10) then
     a2fprefactor= zero
     lqn=zero
     omlog_qn=zero
    else
     a2fprefactor = gammafact(ibranch)*eigval(ibranch)/(two_pi*abs(phfrq(ibranch))*n0(isppol))
     lqn= gammafact(ibranch)*eigval(ibranch)/(pi*phfrq(ibranch)**2*n0(isppol))
     omlog_qn =  lqn*log(abs(phfrq(ibranch)))
    end if
!   ENDMG
!   if (a2fprefactor < tol12) cycle
!   write (*,*) ' a2fp,eigv,phf ',a2fprefactor,eigval(ibranch),phfrq(ibranch)

!   Add contribution to average elphon coupling
!   MANY ISSUES WITH FINITE T SUMS. THIS IS DEFINITELY
!   NOT A CORRECT FORMULATION YET.

!   MG  Added gammafact to treat the doscalprod=1 case
!   Added avglambda and avgomglog to calculate lamda and omega_log using the sum over the kpt-grid.
!   If the k-grid is dense enough, these values should be better than the corresponding quantities
!   evaluated through the integration over omega that depends on the elph_ds%a2fsmear

    avgelphg = avgelphg + gammafact(ibranch)*eigval(ibranch)
    avglambda = avglambda + lqn
    avgomlog= avgomlog + omlog_qn
!   ENDMG

!   omega = zero
    omega = elph_ds%omega_min
    tmpa2f(:) = zero
    tmpphondos(:) = zero
    do iomega=1,elph_ds%na2f
     xx = (omega-phfrq(ibranch))*gaussfactor
     gaussval = gaussprefactor*exp(-xx*xx)

     tmpa2f(iomega) = tmpa2f(iomega) + gaussval*a2fprefactor
     tmpphondos(iomega) = tmpphondos(iomega) + gaussval

     omega = omega+elph_ds%domega
    end do

    a2f_1d(:) = a2f_1d(:) + tmpa2f(:)
    dos_phon(:) = dos_phon(:) + tmpphondos(:)

   end do
!  ! end ibranch do

  end do
! end iFSqpt do

! write (*,*) 'a2f_1d = ', a2f_1d

! second 1 / elph_ds%nFSkpt factor for the integration weights
  a2f_1d(:) = a2f_1d(:) / elph_ds%nFSkpt
  dos_phon(:) = dos_phon(:) / elph_ds%nFSkpt

! MG
  avglambda = avglambda/ elph_ds%nFSkpt
  avgomlog= avgomlog/ elph_ds%nFSkpt
  avgomlog = exp (avgomlog/avglambda)
  write(*,*) ' from mka2f: for sppol ', isppol
  write(*,*) ' lambda  = ',avglambda,' omega_log= ',avgomlog
! ENDMG

  write (*,'(a,I4,a,E16.6)') '# The DOS at Fermi level for spin ', isppol, ' is ', n0(isppol)
  write (unit_a2f,'(a,I4,a,E16.6)') '# The DOS at Fermi level for spin ', isppol, ' is ', n0(isppol)
  write (unit_a2f,'(a)') '#'
! write (*,'(6(E16.6,2x))') a2f_1d(:)

! omega = zero
  omega = elph_ds%omega_min
  do iomega=1,elph_ds%na2f
   write (*,*) omega, a2f_1d(iomega)
   write (unit_a2f,*) omega, a2f_1d(iomega)
   omega=omega+elph_ds%domega
  end do
  write (*,*)
  write (unit_a2f,*)

! output the phonon DOS, but only for the first sppol case
  if (isppol == 1) then
!  omega = zero
   write (*,'(a)') '# The phonon DOS:'
   omega = elph_ds%omega_min
   do iomega=1,elph_ds%na2f
    write (*,*) omega, dos_phon(iomega)
    write (unit_phdos,*) omega, dos_phon(iomega)
    omega=omega+elph_ds%domega
   end do
   write (*,*)
   
  end if

! Do isotropic calculation of lambda and output lambda, Tc(MacMillan)

! if (elph_ds%omega_min < -tol10) then
! write (*,*) ' mka2f: WARNING: omega_min < 0. Skipping calculation of Tc'
! else if (elph_ds%omega_min > tol8) then
! write (*,*) ' mka2f: ERROR: omega_min should be either negative or 0.0'
! stop
! else
  allocate (a2f_1mom(elph_ds%na2f),a2f_1mom_int(elph_ds%na2f))
  allocate (a2f1mom(elph_ds%na2f),a2f2mom(elph_ds%na2f))
  allocate (a2f3mom(elph_ds%na2f),a2f4mom(elph_ds%na2f))
  
! omega = elph_ds%domega
  omega = elph_ds%omega_min
  a2f_1mom(:) = zero
  do iomega=1,elph_ds%na2f
   if (abs(omega) > tol10) then
!   first inverse moment of alpha2F
    a2f_1mom(iomega) = two*a2f_1d(iomega)/abs(omega)
!   first positive moment of alpha2F
    a2f1mom(iomega) = two*a2f_1d(iomega)*abs(omega)
!   second positive moment of alpha2F
    a2f2mom(iomega) = a2f1mom(iomega)*abs(omega)
!   third positive moment of alpha2F
    a2f3mom(iomega) = a2f2mom(iomega)*abs(omega)
!   fourth positive moment of alpha2F
    a2f4mom(iomega) = a2f3mom(iomega)*abs(omega)
   end if
   omega=omega+elph_ds%domega
  end do
! 
! From Allen PRL 59 1460
! \lambda <\omega^n> = 2 \int_0^{\infty} d\omega [\alpha^2F / \omega] \omega^n
! 
  call simpson_int(elph_ds%na2f,elph_ds%domega,a2f_1mom,a2f_1mom_int)
  lambda_iso = a2f_1mom_int(elph_ds%na2f)
  call simpson_int(elph_ds%na2f,elph_ds%domega,a2f1mom,a2f_1mom_int)
  lambda_2 = a2f_1mom_int(elph_ds%na2f)
  call simpson_int(elph_ds%na2f,elph_ds%domega,a2f2mom,a2f_1mom_int)
  lambda_3 = a2f_1mom_int(elph_ds%na2f)
  call simpson_int(elph_ds%na2f,elph_ds%domega,a2f3mom,a2f_1mom_int)
  lambda_4 = a2f_1mom_int(elph_ds%na2f)
  call simpson_int(elph_ds%na2f,elph_ds%domega,a2f4mom,a2f_1mom_int)
  lambda_5 = a2f_1mom_int(elph_ds%na2f)

  deallocate (a2f_1mom,a2f_1mom_int)
  deallocate (a2f1mom,a2f2mom,a2f3mom,a2f4mom)

  write (*,*) 'mka2f: elphon coupling lambdas for isppol = ', isppol
  write (*,*) 'mka2f: isotropic lambda', lambda_iso

  write (*,*) 'mka2f: positive moments of alpha2F:'
  write (*,*) 'lambda <omega^2> = ', lambda_2
  write (*,*) 'lambda <omega^3> = ', lambda_3
  write (*,*) 'lambda <omega^4> = ', lambda_4
  write (*,*) 'lambda <omega^5> = ', lambda_5

! Get log moment of alpha^2F
  allocate (a2flogmom(elph_ds%na2f),a2flogmom_int(elph_ds%na2f))
  omega = elph_ds%omega_min
! omega = elph_ds%domega
  a2flogmom(:) = zero
  do iomega=1,elph_ds%na2f
   if (abs(omega) > tol10) then
    a2flogmom(iomega) = (two/lambda_iso)*a2f_1d(iomega)*log(abs(omega))/abs(omega)
   end if
   omega=omega+elph_ds%domega
  end do
  call simpson_int(elph_ds%na2f,elph_ds%domega,a2flogmom,a2flogmom_int)
  omegalog = exp(a2flogmom_int(elph_ds%na2f))

  deallocate (a2flogmom,a2flogmom_int)
  
  tc_macmill = omegalog/1.2_dp * exp((         -1.04_dp*(one+lambda_iso)        ) /&
&  (lambda_iso-mustar*(one+0.62_dp*lambda_iso)))

  if (elph_ds%nsppol > 1) then
   write (message, '(3a)' ) ch10,&
&   ' Warning : some of the following quantities should be integrated over spin', ch10
   call wrtout(6,message,'COLL')
   call wrtout(ab_out,message,'COLL')
  end if

  write (message, '(3a)' ) ch10,&
&  ' Superconductivity : isotropic evaluation of parameters from electron-phonon coupling.',ch10
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  write (message, '(a,es16.6)' )&
&  ' mka2f: isotropic lambda = ', lambda_iso
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  write (message, '(a,es16.6)' )&
&  ' mka2f: lambda <omega^2> = ', lambda_2
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  write (message, '(a,es16.6)' )&
&  ' mka2f: lambda <omega^3> = ', lambda_3
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  write (message, '(a,es16.6)' )&
&  ' mka2f: lambda <omega^4> = ', lambda_4
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  write (message, '(a,es16.6)' )&
&  ' mka2f: lambda <omega^5> = ', lambda_5
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  write (message, '(a,es16.6,a,es16.6,a)' )&
&  ' mka2f: omegalog  = ', omegalog, ' (Ha) ', omegalog/kb_HaK, ' (Kelvin) '
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  write (message, '(a,es16.6)' )&
&  ' mka2f: input mustar = ', mustar
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  write ( message, '(a,es16.6,a,es16.6,a)')&
&  ' mka2f: MacMillan Tc = ', tc_macmill, ' (Ha) ', tc_macmill/kb_HaK, ' (Kelvin) '
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')


 end do ! loop over sppol

 deallocate(matrx)
 deallocate(zhpev1)
 deallocate(zhpev2)

 close (unit=unit_a2f)
 close (unit=unit_phdos)

 write (*,*) ' mka2f : end '

end subroutine mka2f
!!***
