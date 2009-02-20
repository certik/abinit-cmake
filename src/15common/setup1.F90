!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup1
!!
!! NAME
!! setup1
!!
!! FUNCTION
!! Call near top of main routine to handle setup of various arrays,
!! filenames, checking of input data, etc.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  acell(3)=length scales (bohr)
!!  amu(ntypat)=atomic masses for each type of atom (amu)
!!  ecut_eff=effective energy cutoff (hartree) for planewave basis sphere
!!  ecutc_eff=- PAW only - effective energy cutoff (hartree) for the coarse grid
!!  iboxcut=1 if gsqcut has to define a sphere containing the whole FFT box
!!  intxc=0 for old, 1 for new xc quadrature
!!  ionmov=0 for no moving atoms; 1 for molecular dynamics; 2 for Broyden
!!  natom=number of atoms
!!  nband(nkpt*nsppol)=number of bands at each k point
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ngfftc(18)=contain all needed information about 3D FFT for the coarse grid
!!  nkpt=number of k points
!!  nqpt=1 if there is a q-wavevector
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetries
!!  ntypat=number of types of atoms
!!  qptn(3) see nqpt
!!  response=0 if called by gstate, =1 if called by respfn
!!  rprim(3,3)=dimensionless real space primitive translations
!!  typat(natom)=atom type for each atom
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!!  amass(natom)=atomic masses for each atom (in atomic units, where the electron mass is one)
!!  bantot=total number of bands at all k points
!!  gmet(3,3)=metric for reciprocal space inner products (bohr^-2)
!!  gprimd(3,3)=dimens. primitive translations for reciprocal space (bohr**-1)
!!  gsqcut_eff=Fourier cutoff on G^2 for "large sphere" of radius double
!!  gsqcutc_eff=(PAW) Fourier cutoff on G^2 for "large sphere" of radius double for the coarse FFT grid
!!   that of the basis sphere--appropriate for charge density rho(G),
!!   Hartree potential, and pseudopotentials, corresponding to ecut_eff
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol=unit cell volume (bohr^3)
!!
!! NOTES
!! SHOULD BE CLEANED !
!!
!! PARENTS
!!      gstate,nonlinear,respfn
!!
!! CHILDREN
!!      getcut,leave_new,metric,mkrdim,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine setup1(acell,amass,amu,bantot,&
&  ecut_eff,ecutc_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,iboxcut,intxc,ionmov,&
&  natom,nband,ngfft,ngfftc,nkpt,nqpt,nsppol,nsym,ntypat,&
&  qptn,response,rmet,rprim,rprimd,typat,ucvol,usepaw)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_13recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iboxcut,intxc,ionmov,natom,nkpt,nqpt,nsppol,nsym,ntypat
 integer,intent(in) :: response,usepaw
 integer,intent(out) :: bantot
 real(dp),intent(in) :: ecut_eff,ecutc_eff
 real(dp),intent(out) :: gsqcut_eff,gsqcutc_eff,ucvol
!arrays
 integer,intent(in) :: nband(nkpt*nsppol),ngfft(18),ngfftc(18),typat(natom)
 real(dp),intent(in) :: acell(3),amu(ntypat),qptn(3),rprim(3,3)
 real(dp),intent(out) :: amass(natom),gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),intent(out) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ia1,ia2,iatom,ii,ikpt,isppol,last,nkpt_me,nuniqu
 real(dp) :: boxcut,boxcutc,eps,norm
 character(len=500) :: message
!arrays
 real(dp) :: k0(3)

! ************************************************************************

!DEBUG
!write(6,*)' setup1 : enter '
!if(.true.)stop
!ENDDEBUG

!Compute bantot
 bantot=0
 do isppol=1,nsppol
  do ikpt=1,nkpt
   bantot=bantot+nband(ikpt+(isppol-1)*nkpt)
  end do
 end do

 if(nqpt>1.or.nqpt<0) then
  write(message, '(a,a,a,a,i12,a,a,a,a,a)' ) ch10,&
&  ' setup1 : ERROR -',ch10,&
&  '  nqpt =',nqpt,' is not allowed',ch10,&
&  '  (only 0 or 1 are allowed).',ch10,&
&  '  Action : correct your input file.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!Compute dimensional primitive translations rprimd
 call mkrdim(acell,rprim,rprimd)

!Obtain dimensional translations in reciprocal space gprimd,
!metrics and unit cell volume, from rprimd.
!Also output rprimd, gprimd and ucvol
 call metric(gmet,gprimd,ab_out,rmet,rprimd,ucvol)

!Assign masses to each atom (for MD)
 do iatom=1,natom
  amass(iatom)=amu_emass*amu(typat(iatom))
 end do

!Get boxcut for given acell, gmet, ngfft, and ecut_eff
!(center at 000 for groundstate, center at q for respfn):
!boxcut=ratio of basis sphere diameter to fft box side
 k0(:)=0.0_dp
 if(response==1 .and. nqpt==1)then
  k0(:)=qptn(:)
  write(message, '(a)' )&
&  ' setup1 : take into account q-point for computing boxcut.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(6,message,'COLL')
 end if
 if (usepaw==1) then
  write(message,'(2a)') ch10,' Coarse grid specifications (used for wave-functions):'
  call wrtout(ab_out,message,'COLL')
  call wrtout(6,message,'COLL')
  call getcut(boxcutc,ecutc_eff,gmet,gsqcutc_eff,iboxcut,ab_out,k0,ngfftc)
  write(message,'(2a)') ch10,' Fine grid specifications (used for densities):'
  call wrtout(ab_out,message,'COLL')
  call wrtout(6,message,'COLL')
  call getcut(boxcut,ecut_eff,gmet,gsqcut_eff,iboxcut,ab_out,k0,ngfft)
 else
  call getcut(boxcut,ecut_eff,gmet,gsqcut_eff,iboxcut,ab_out,k0,ngfft)
  gsqcutc_eff=gsqcut_eff
 end if

!Check that boxcut>=2 if intxc=1; otherwise intxc must be set=0
 if (boxcut<2.0_dp.and.intxc==1) then
  write(message, '(a,a,a,a,es12.4,a,a,a,a,a)' ) ch10,&
&  ' setup1: ERROR -',ch10,&
&  '  boxcut=',boxcut,' is < 2.0  => intxc must be 0;',ch10,&
&  '  Need larger ngfft to use intxc=1.',ch10,&
&  '  Action : you could increase ngfft, or decrease ecut, or put intxc=0.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!DEBUG
!write(6,*)' setup1 : exit '
!stop
!ENDDEBUG

end subroutine setup1
!!***
