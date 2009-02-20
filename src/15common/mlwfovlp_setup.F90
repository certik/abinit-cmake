!{\src2tex{textfont=tt}}
!!****f* ABINIT/mlwfovlp_setup
!! NAME
!! mlwfovlp_setup
!!
!! FUNCTION
!! Routine which creates table g1 and ovikp  necessary to compute
!! overlap for Wannier code (www.wannier.org f90 version).
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (BAmadon,FJollet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atom_symbols(natom)= table of symbol for each atom
!!                                          and each |p_lmn> non-local projector
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gprimd(3,3)        =dimensional reciprocal space primitive translations
!!  lwanniersetup= flag: only 1 is fully working.
!!  natom              =number of atoms in cell.
!!  mband=maximum number of bands
!!  mbandw=maximum number of bands
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpw=maximum dimensioned size of npw.
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nkpt=number of k points.
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  num_bands=number of bands actually used to construct the wannier function
!!  nwan= number of wannier fonctions (read in wannier90.win).
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  real_lattice(3,3)=dimensional primitive translations for real space
!!                 in format required by wannier90
!!  recip_lattice(3,3)=dimensional primitive translations for reciprocal space
!!                 in format required by wannier90
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  xcart(3,natom)=atomic coordinates in bohr
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  band_in(mband)   = band to take into account for wannier calculation
!!  g1(3,nkpt,nntot) = G vector shift which is necessary to obtain k1+b
!!                     from k2 in the case where k1+b does not belong to the 1st BZ.
!!  nband_inc = # of included bands
!!  nntot            = number of k-point neighbour
!!  ovikp(nkpt,nntot)= gives  nntot value of k2 (in the BZ) for each k1  (k2=k1+b mod(G))
!!  
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      mlwfovlp
!!
!! CHILDREN
!!      atmdata,leave_new,wannier_setup,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine mlwfovlp_setup(atom_symbols,band_in,dtset,eigen,gamma_only,&
& g1,gprimd,lwanniersetup,mband,mbandw,mkmem,mpw,natom,nattyp,nband_inc,nkpt,&
& nntot,nsppol,nspinor,ntypat,num_bands,num_nnmax,nwan,ovikp,&
& proj_l,proj_m,proj_radial,proj_site,proj_x,proj_z,proj_zona,&
& real_lattice,recip_lattice,rprimd,spinors,xcart,xred)

 use defs_basis
 use defs_datatypes
 use defs_wannier90


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments---------------------------
! scalars
!scalars
 integer,intent(in) :: lwanniersetup,mband,mkmem,mpw,natom,nkpt,nspinor,nsppol
 integer,intent(in) :: ntypat,num_nnmax
 integer,intent(out) :: mbandw,nband_inc,nntot,num_bands,nwan
 logical,intent(in) :: gamma_only,spinors
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: nattyp(ntypat)
 integer,intent(out) :: g1(3,nkpt,num_nnmax),ovikp(nkpt,num_nnmax)
 integer,intent(out) :: proj_l(mband),proj_m(mband),proj_radial(mband)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),gprimd(3,3),real_lattice(3,3)
 real(dp),intent(in) :: recip_lattice(3,3),rprimd(3,3),xred(3,natom)
 real(dp),intent(out) :: proj_site(3,mband),proj_x(3,mband),proj_z(3,mband)
 real(dp),intent(out) :: proj_zona(mband),xcart(3,natom)
 logical,intent(out) :: band_in(mband)
 character(len=3),intent(out) :: atom_symbols(natom)

!Local variables---------------------------
!scalars
 integer :: iatom,icb,ikpt,ikpt1,intot,itypat,jj
 real(dp) :: amu,rcov,znucl1
 character(len=2) :: symbol
 character(len=500) :: message
 character(len=9) :: seed_name
!arrays
 integer :: exclude_bands(mband),ngkpt(3)

! *************************************************************************

!^^^^^^^^^^^^^^^^read wannier90.nnkp^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 if(lwanniersetup==0) then 
  open(unit=111,file='wannier90.nnkp',form='formatted',status='old')
  read(111,*)
  read(111,*) nntot , mbandw, nwan
  write(message, '(a,a,i6,i6,i6)' )ch10,&
&  ' mlwfovlp_setup nntot,mbandw,nwan ', nntot,mbandw,nwan
  call wrtout(6,message,'COLL')
  if(mbandw.ne.mband) then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' mlwfovlp_setup : ERROR -',ch10,&
&   '  mbandw is not equal to mband ',ch10,&
&   '  Action : check wannier90.nnkp'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
   stop
  end if
  if(nwan>mbandw) then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' mlwfovlp_setup : ERROR -',ch10,&
&   '  nwan > mbandw ',ch10,&
&   '  Action : check wannier90.nnkp'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
  if(nwan==0) then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' mlwfovlp_setup : ERROR -',ch10,&
&   '  nwan = 0 ',ch10,&
&   '  Action : check wannier90.nnkp'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
  do ikpt=1,nkpt
   do intot=1,nntot
!   ikpt1: k point  (ikpt=ikpt1)
!   ovikp(intot,ikpt): neighbour number intot for ikpt
!   g1(1:3,intot,ikpt): non reciprocal space vector between the 2 k-points
    read(111,*)  &
&    ikpt1,ovikp(ikpt,intot),(g1(jj,ikpt,intot),jj=1,3)
    if(ikpt1.ne.ikpt) write(6,*) "warning: ikpt1 .ne ikpt : ?"
   end do
  end do
  close(111)
  write(6,*) "wannier90.nnkp has been read !"
  write(6,*) "exclude bands is not given in this case (not implemented) : STOP"
  stop
! ^^^^^^^^^^^^^^^^^^^^^^^ call wannier_setup begin^^^^^^^^^^^^^^^^^^^^^^^^
 else if (lwanniersetup==1) then
  num_bands=mband
! num_nnmax=12 !limit fixed for compact structure in wannier_setup.
  ovikp=0.d0
! "When nshiftk=1, kptrlatt is initialized as a diagonal (3x3) matrix, whose diagonal 
! elements are the three values ngkpt(1:3)"
  ngkpt(1)=dtset%kptrlatt(1,1)
  ngkpt(2)=dtset%kptrlatt(2,2) !  have to verif kptrlatt is diagonal
  ngkpt(3)=dtset%kptrlatt(3,3)
  seed_name="wannier90"
! if(psps%npsp.ne.psps%ntypat) then; write(6,*) "prb npsp"; stop; endif
  do iatom=1,natom
   itypat=dtset%typat(iatom)
   znucl1=dtset%znucl(itypat)
!  amu and rcov not initialized but are not useful here.
   call atmdata(amu,rcov,symbol,znucl1)
   symbol=trim(adjustl(symbol))
!  write(309,*) symbol
   atom_symbols(iatom)=symbol
   xcart(:,iatom)=rprimd(:,1)*xred(1,iatom)+&
&   rprimd(:,2)*xred(2,iatom)+&
&   rprimd(:,3)*xred(3,iatom)
  end do ! iatom
! write(6,*) xcart
! write(6,*) Bohr_Ang
! write(6,*) rprimd*Bohr_Ang
! write(6,*) gprimd/Bohr_Ang*two_pi
! write(6,*) seed_name
! write(6,*) ngkpt
! write(6,*) nkpt
! write(6,*) mband
! write(6,*) natom
! write(6,*) atom_symbols
  write(message, '(a,a)' )ch10,&
&  '** mlwfovlp_setup:  call wannier90 library subroutine wannier_setup'
  call wrtout(6,message,'COLL')
#if defined HAVE_WANNIER90

  call wannier_setup(seed_name,ngkpt,nkpt&                    !input
& ,real_lattice,recip_lattice,dtset%kpt&                      !input
& ,mband,natom,atom_symbols,xcart*Bohr_Ang&                   !input
& ,gamma_only,spinors&                                        !input
& ,nntot,ovikp,g1,num_bands,nwan&                             !output
& ,proj_site,proj_l,proj_m,proj_radial,proj_z&                !output
& ,proj_x,proj_zona,exclude_bands)                            !output
#endif
! write(6,*)  "1", nntot,nwan
! write(6,*)  "2",num_bands  ! states on which wannier functions are computed
! write(6,*)  "3", proj_site
! write(6,*)  "4",proj_l
! write(6,*)  "5",proj_m
! write(6,*)  "6", proj_radial
! write(6,*)  "7", proj_z
! write(6,*)  "8", proj_x
! write(6,*)  "9",proj_zona
! write(6,*)  "10",exclude_bands
! testdebug:  ovikp(1,1)=1
! ^^^^^^^^^^^^^^^^^^^^^^^ end ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  mbandw=mband
 end if  ! lwanniersetup
 band_in=.true.
 do icb=1,mband
  if(exclude_bands(icb).ne.0)  band_in(exclude_bands(icb))=.false.
 end do
 nband_inc=0
 do icb=1, mband
  if (band_in(icb)) then
   nband_inc=nband_inc+1
  end if
 end do
 if(mband.gt.num_bands) then
  write(message, '(a,a)' )ch10,&
&  '   Following bands are excluded from the calculation of wannier functions:'
  call wrtout(6,message,'COLL')
  write(message,fmt=10) exclude_bands(1:mband-num_bands)
  call wrtout(6,message,'COLL')
 end if
 write(message, '(a,i6,a)' )ch10,&
& nwan,' wannier functions will be computed (see wannier90.win)'
 call wrtout(6,message,'COLL')
!write(6,*) exclude_bands(icb),band_in(icb)
 10 format(300(2x,i5))
!^^^^^^^^^^^^^^^END OF READING
 write(message, '(a,i6,a)' )ch10,&
& num_bands,' bands will be used to extract wannier functions'
 call wrtout(6,message,'COLL')
 if(num_bands.lt.nwan) then
  write(message, '(a,a,a,a,a,a)' )ch10,&
&  ' mlwfovlp_setup : ERROR -',ch10,&
&  ' number of bands is lower than the number of wannier functions',ch10,&
&  ' Action : check input file and wannier90.win'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
  stop
 else if (num_bands==nwan) then
  write(message, '(a,a,a,a)' )ch10,&
&  '   Number of bands is equal to the number of wannier functions',ch10,&
&  '   Disentanglement will not be necessary'
  call wrtout(6,message,'COLL')
 else if  (num_bands.gt.nwan) then
  write(message, '(a,a,a,a)' )ch10,&
&  '   Number of bands is larger than the number of wannier functions',ch10,&
&  '   Disentanglement will be necessary'
  call wrtout(6,message,'COLL')
 end if
 write(message, '(2x,a,a,i3,1x,a)' )ch10,&
& '   Each k-point has', nntot,'neighbours'
 call wrtout(6,message,'COLL')


!DEBUG
!write(6,*)' mlwfovlp_setup : exit'
!stop
!ENDDEBUG

 end subroutine    mlwfovlp_setup
!!***

