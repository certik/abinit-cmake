!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkifc9
!!
!! NAME
!! mkifc9
!!
!! FUNCTION
!! This routine makes the interatomic force constants,
!! taking into account the dipole-dipole interaction,
!! and also perform some checks.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=length scales of cell (bohr)
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! anaddb_dtset= (derived datatype) contains all the input variables
!! blkflg(nsize,nblok)= flag of existence for each element of the DDB
!! blknrm(1,nblok)=norm of qpt providing normalization
!! blkqpt(1<ii<9,nblok)=q vector of a phonon mode (ii=1,2,3)
!! blktyp(nblok)=1 or 2 depending on non-stationary or stationary block
!!  3 for third order derivatives
!! blkval(2,3*mpert*3*mpert,nblok)= all the dynamical matrices
!! dielt(3,3)=dielectric tensor
!! gmet(3,3)=reciprocal space metric tensor in bohr**-2
!! gprim(3,3)=dimensionless reciprocal space primitive translations
!! indsym(4,msym*natom)=indirect indexing array for symmetries
!! iout=unit number for output of formatted data
!! mpert =maximum number of ipert
!! msym =maximum number of symmetries
!! natom=number of atoms in cell
!! nblok= number of bloks in the DDB
!! nsym=number of space group symmetries
!! ntypat=number of atom types
!! rmet(3,3)=metric tensor in real space (bohr^2)
!! rprim(3,3)= primitive translation vectors
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!! tcpui,twalli=initial values of cpu and wall clocktime
!! tnons(3,nsym)=fractional nonsymmorphic translations
!! typat(natom)=type integer for each atom in cell
!! ucvol=unit cell volume in bohr**3
!! xred(3,natom)=fractional dimensionless atomic coordinates
!! zeff(3,3,natom)=effective charge on each atom, versus electric
!!  field and atomic displacement
!!
!! OUTPUT
!! atmfrc(2,3,natom,3,natom,nrpt)
!!  = Interatomic Forces in real space
!! d2cart(2,3,mpert,3,mpert) work space for containing
!!  the dynamical matrix in cartesian coordinates
!! displ(2*3*natom*3*natom)= after phfrq3, contain
!!  the displacements of atoms in cartesian coordinates.
!! dyewq0(3,3,natom)=atomic self-interaction correction to the
!!  dynamical matrix. (only when dipdip=1)
!! eigval(3*natom)=contains the eigenvalues of the dynamical matrix
!! eigvec(2*3*natom*3*natom)= at the end, contain
!!  the eigenvectors of the dynamical matrix.
!! nrpt = number of vectors of the lattice in real space
!! phfrq(3*natom)=phonon frequencies (square root of the dynamical
!!  matrix eigenvalues, except if these are negative, and in this
!!  case, give minus the square root of the absolute value
!!  of the matrix eigenvalues). Hartree units.
!! rcan(3,natom) = Atomic position in canonical coordinates
!! rpt(3,nprt) = Canonical coordinates of the R points in the unit cell
!! trans(3,natom) = Atomic translations : xred = rcan + trans
!! wghatm(natom,natom,nrpt)
!!  = Weight associated to the couple of atoms
!!  and the R vector
!!
!! NOTES
!! 1. Should be executed by one processor only.
!! 2. Work variables are:
!!   d2cart,displ,eigval,eigvec,phfrq
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      asrif9,bigbx9,canat9,chkrp9,dymfz9,ewald9,ftiaf9,gtdyn9,hybrid9,mkrdim
!!      nanal9,phfrq3,prtph3,q0dy3,rsiaf9,smpbz,symdm9,timein,wght9,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mkifc9(acell,amu,anaddb_dtset,atmfrc,&
& blkflg,blknrm,blkqpt,blktyp,blkval,dielt,displ,dyewq0,&
& d2cart,eigval,eigvec,gmet,gprim,&
& indsym,iout,mpert,msym,natom,nblok,nrpt,&
& nsym,ntypat,phfrq,rcan,rmet,rprim,&
& rpt,symrec,symrel,tcpui,tnons,trans,twalli,typat,&
& ucvol,wghatm,xred,zeff)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_13recipspace
 use interfaces_16response
 use interfaces_17ddb, except_this_one => mkifc9
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iout,mpert,msym,natom,nblok,nsym,ntypat
 integer,intent(inout) :: nrpt
 real(dp),intent(in) :: tcpui,twalli,ucvol
 type(anaddb_dataset_type),intent(in) :: anaddb_dtset
!arrays
 integer,intent(in) :: blkflg(*),blktyp(*),indsym(4,nsym*natom)
 integer,intent(in) :: symrec(3,3,msym),symrel(3,3,msym),typat(natom)
 real(dp),intent(in) :: acell(3),amu(ntypat),blknrm(3,*),blkqpt(9,*)
 real(dp),intent(in) :: blkval(2,*),gmet(3,3),gprim(3,3),rmet(3,3),rprim(3,3)
 real(dp),intent(in) :: tnons(3,msym),xred(3,natom)
 real(dp),intent(inout) :: dielt(3,3),zeff(3,3,natom)
 real(dp),intent(out) :: atmfrc(2,3,natom,3,natom,nrpt)
 real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert),displ(2*3*natom*3*natom)
 real(dp),intent(out) :: dyewq0(3,3,natom),eigval(*),eigvec(*),phfrq(*)
 real(dp),intent(out) :: rcan(3,natom),rpt(3,nrpt),trans(3,natom)
 real(dp),intent(out) :: wghatm(natom,natom,nrpt)

!Local variables -------------------------
!scalars
 integer :: asr,brav,choice,dipdip,idir1,idir2,ii,ipert1,ipert2,iqpt,mqpt,mrpt
 integer :: nqpt,nqshft,option,plus,qtor,sumg0
 real(dp) :: qphnrm,tcpu,twall
 character(len=500) :: message
!arrays
 integer :: ngqpt(9),qptrlatt(3,3)
 integer,allocatable :: atifc(:)
 real(dp) :: qpt(3),rprimd(3,3)
 real(dp),allocatable :: dyew(:,:,:,:,:),dynmat(:,:,:,:,:,:),spqpt(:,:)

!******************************************************************

 asr=anaddb_dtset%asr
 brav=anaddb_dtset%brav
 dipdip=anaddb_dtset%dipdip
 ngqpt(1:3)=anaddb_dtset%ngqpt(1:3)
 nqshft=anaddb_dtset%nqshft

!DEBUG
!write(6,*)' mkifc9 : ngqpt(:)=',ngqpt(:)
!ENDDEBUG

 allocate(atifc(natom))
 atifc(:)=anaddb_dtset%atifc(:)

 call mkrdim(acell,rprim,rprimd)

 if(dipdip==1)then

  if(asr==1.or.asr==2)then
!  Calculation of the non-analytical part for q=0
   sumg0=0
   qpt(1)=0.0_dp
   qpt(2)=0.0_dp
   qpt(3)=0.0_dp
   allocate(dyew(2,3,natom,3,natom))
   call ewald9(acell,dielt,dyew,gmet,gprim,natom,&
&   qpt,rmet,rprim,sumg0,ucvol,xred,zeff)
   option=asr
   call q0dy3(natom,dyewq0,dyew,option)
   deallocate(dyew)
  else if (asr==0) then
   dyewq0(:,:,:)=0.0_dp
  end if
 end if

!Check if the rprim are coherent with the choice used in
!the interatomic forces generation
 write(6, '(a)' )' mkifc9 : enter chkrp9 '
 call chkrp9(brav,rprim)
 write(6, '(a)' )' mkifc9 : exit chkrp9 '

!Sample the Brillouin zone
 write(6, '(a)' )' mkifc9 : enter smpbz '
 option=1
 qptrlatt(:,:)=0
 qptrlatt(1,1)=ngqpt(1)
 qptrlatt(2,2)=ngqpt(2)
 qptrlatt(3,3)=ngqpt(3)
 mqpt=ngqpt(1)*ngqpt(2)*ngqpt(3)*nqshft
 if(brav==2)mqpt=mqpt/2
 if(brav==3)mqpt=mqpt/4
 allocate(spqpt(3,mqpt))
 call smpbz(brav,iout,qptrlatt,mqpt,nqpt,nqshft,option,anaddb_dtset%q1shft,spqpt)
 write(6, '(a)' )' mkifc9 : exit smpbz '

 allocate(dynmat(2,3,natom,3,natom,nqpt))

!Find symmetrical dynamical matrices
 write(6, '(a)' )' mkifc9 : enter symdm9 '
 call symdm9(acell,blkflg,blknrm,blkqpt,blktyp,blkval,&
& dynmat,gprim,indsym,mpert,natom,nblok,nqpt,nsym,anaddb_dtset%rfmeth,rprim,spqpt,&
& symrec,symrel,xred,tnons,ucvol)
 write(6, '(a)' )' mkifc9 : exit symdm9 '

 if(dipdip==1)then
! Take off the dipole-dipole part of the dynamical matrix
  write(6, '(a)' )' mkifc9 : will extract the dipole-dipole part,'
  write(6, '(a)' )' using ewald9, q0dy3 and nanal9 for every wavevector.'
  do iqpt=1,nqpt
   write(6, '(a,i4)' )'    wavevector number :',iqpt
   qpt(1)=spqpt(1,iqpt)
   qpt(2)=spqpt(2,iqpt)
   qpt(3)=spqpt(3,iqpt)
   sumg0=0

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-ewald9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(6,message,'COLL')
   allocate(dyew(2,3,natom,3,natom))
   call ewald9(acell,dielt,dyew,gmet,gprim,natom,&
&   qpt,rmet,rprim,sumg0,ucvol,xred,zeff)
   option=0

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-q0dy3 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(6,message,'COLL')
   call q0dy3(natom,dyewq0,dyew,option)
   plus=0

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-nanal9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(6,message,'COLL')
   call nanal9(dyew,dynmat,iqpt,natom,nqpt,plus)
   deallocate(dyew)
  end do
 end if

!Now, take care of the remaining part of the
!dynamical matrix

!Passage to canonical normalized coordinates
 write(6, '(a)' )' mkifc9 : enter canat9 '
 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-canat9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(6,message,'COLL')
 call canat9(brav,natom,rcan,rprim,trans,xred)
 write(6, '(a)' )' mkifc9 : exit canat9 '

!Multiply the dynamical matrix by a phase shift
 write(6, '(a)' )' mkifc9 : enter dymfz9 '
 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-dymfz9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(6,message,'COLL')
 option=1
 call dymfz9(dynmat,natom,nqpt,gprim,option,&
& rcan,spqpt,trans)
 write(6, '(a)' )' mkifc9 : exit dymfz9 '

!Create the Big Box of R vectors in real space
 write(6, '(a)' )' mkifc9 : enter bigbx9 '
 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-bigbx9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(6,message,'COLL')

 mrpt=nrpt ; choice=1
 call bigbx9(brav,choice,mrpt,ngqpt,nqshft,nrpt,rprim,rpt)
 write(6, '(a)' )' mkifc9 : exit bigbx9 '

!Weights associated to these R points and to atomic pairs
 write(6, '(a)' )' mkifc9 : enter wght9 '
 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-wght9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(6,message,'COLL')
 call wght9(brav,gprim,natom,ngqpt,nqpt,nqshft,&
& nrpt,anaddb_dtset%q1shft,rcan,rpt,wghatm)
 write(6, '(a)' )' mkifc9 : exit wght9 '

!Fourier transformation of the dynamical matrices
 qtor=1
 write(6, '(a)' )' mkifc9 : enter ftiaf9 '
 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-ftiaf9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(6,message,'COLL')

 call ftiaf9(atmfrc,dynmat,gprim,natom,nqpt,nrpt,qtor,rpt,spqpt,wghatm)
 write(6, '(a)' )' mkifc9 : exit ftiaf9 '

!Eventually impose Acoustic Sum Rule
!to the interatomic forces
 if(asr>0)then
  write(6, '(a)' )' mkifc9 : enter asrif9 '
  call asrif9(asr,atmfrc,natom,nrpt,rpt,wghatm)
  write(6, '(a)' )' mkifc9 : exit asrif9 '
 end if

!*** The interatomic forces have been calculated ! ***
 write(6,*)
 write(6, '(/,a,/,/,a)' )&
& ' mkifc9 : the interatomic forces have been obtained',&
& ' mkifc9 : analysis of real-space IFCs'
 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-hybrid9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(6,message,'COLL')

!Eventually artificially modify the IFC
 write(6, '(a)' )' mkifc9 : enter hybrid9'
 call  hybrid9(acell,asr,atmfrc,dielt,dipdip,dyew,dyewq0,&
& gmet,gprim,iout,natom,nrpt,rcan,rmet,&
& rprim,rpt,ucvol,wghatm,xred,zeff)
 write(6, '(a)' )' mkifc9 : exit hybrid9'


!Analysis of the real-space interatomic force constants
 if (anaddb_dtset%nsphere/=0 .or. anaddb_dtset%rifcsph>tol10 .or. anaddb_dtset%ifcout/=0) then

  call timein(tcpu,twall)
  write(message, '(a,f11.3,a,f11.3,a)' )&
&  '-rsiaf9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
  call wrtout(6,message,'COLL')

  call rsiaf9(acell,atifc,atmfrc,dielt,dipdip,dyewq0,&
&  gmet,gprim,anaddb_dtset%ifcana,anaddb_dtset%ifcout,iout,&
&  natom,nrpt,anaddb_dtset%nsphere,rcan,anaddb_dtset%rifcsph,rmet,rprim,rpt,&
&  tcpui,twalli,ucvol,wghatm,zeff)
  write(6, '(a)' )' mkifc9 : real space analysis of the IFC''s done '
 end if

!Eventually impose Acoustic Sum Rule
!to the interatomic forces
!(Note : here, after the analysis, in which the range
!of the short-range IFCs may have been changed
!That is why it is asked that asr be negative)
 if(asr<0)then
  asr=-asr

  call timein(tcpu,twall)
  write(message, '(a,f11.3,a,f11.3,a)' )&
&  '-asrif9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
  call wrtout(6,message,'COLL')

  write(6, '(a)' )' mkifc9 : enter asrif9 '
  call asrif9(asr,atmfrc,natom,nrpt,rpt,wghatm)
  write(6, '(a)' )' mkifc9 : exit asrif9 '
 end if

!Check that the starting values are well reproduced.
!(This is to be suppressed in a future version)
 write(6, '(a,a)' )' mkifc9 : now check that the starting values ',&
& ' are reproduced after the use of interatomic forces '
!call timein(tcpu,twall)
!write(6,1000)tcpu-tcpui,twall-twalli
 do iqpt=1,nqpt
  qpt(1)=spqpt(1,iqpt)
  qpt(2)=spqpt(2,iqpt)
  qpt(3)=spqpt(3,iqpt)
  qphnrm=1.0_dp

  call timein(tcpu,twall)
  write(message, '(a,f11.3,a,f11.3,a)' )&
&  '-gtdyn9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
  call wrtout(6,message,'COLL')

! The dynamical matrix d2cart is calculated here :
  call gtdyn9(acell,atmfrc,dielt,dipdip,&
&  dyewq0,d2cart,gmet,gprim,mpert,natom,&
&  nrpt,qphnrm,qpt,rcan,rmet,rprim,rpt,&
&  trans,ucvol,wghatm,xred,zeff)

! Calculation of the eigenvectors and eigenvalues
! of the dynamical matrix
  call phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,&
&  mpert,msym,natom,nsym,ntypat,phfrq,qphnrm,qpt,rprimd,anaddb_dtset%symdynmat,symrel,typat,ucvol,xred)

! Write the phonon frequencies (this is for checking
! purposes).
! Note : these phonon frequencies are not written on unit
! iout, only on unit 6.
  call prtph3(displ,0,anaddb_dtset%enunit,-1,-1,natom,phfrq,qphnrm,qpt)

 end do

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-end of mkifc9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(6,message,'COLL')

 deallocate(atifc,dynmat,spqpt)

end subroutine mkifc9
!!***
