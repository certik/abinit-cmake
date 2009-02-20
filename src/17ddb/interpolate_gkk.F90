!{\src2tex{textfont=tt}}
!!****f* ABINIT/interpolate_gkk
!!
!! NAME
!! interpolate_gkk
!!
!! FUNCTION
!! This routine interpolates the gkk matrices for all q vectors
!!  between points on the full FSkpt grid.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   acell(3)=length scales of cell (bohr)
!!   amu(ntypat)=mass of the atoms (atomic mass unit)
!!   atmfrc  = inter-atomic force constants from anaddb
!!   dielt(3,3) = dielectric tensor
!!   dipdip  =dipole dipole interaction flag
!!   dyewq0(3,3,natom)=atomic self-interaction correction to the
!!        dynamical matrix (only when dipdip=1)
!!   elph_ds = elphon datastructure with data and dimensions
!!   FSkptirred = coordinates of irreducible kpoints close to the FS
!!   FSkpt = coordinates of all kpoints close to the FS
!!   ftwghtgkk = weights for points in real space, in FT of gkk
!!   gmet(3,3) =metric in reciprocal space
!!   gprim(3,3) =dimensionless basis vectors of reciprocal space
!!   indsym = mapping of atoms btw themselves under symmetry
!!   mpert =maximum number of ipert
!!   msym =maximum number of symmetries
!!   natom=number of atoms in cell
!!   nrpt =number of real space points used to integrate IFC (for
!!        interpolation of dynamical matrices)
!!   nsym=number of space group symmetries
!!   ntypat = number of types of atoms
!!  phon_ds = datastructure with interatomic force constants to interpolate
!!     phonons
!!   rcan(3,natom) =canonical positions of atoms
!!   rmet(3,3)=metric tensor in real space (bohr^2)
!!   rprim(3,3)= primitive translation vectors (normalized)
!!   rprimd(3,3)= primitive translation vectors (dimensionful)
!!   rpt(3,nprt) =canonical positions of R points in the unit cell
!!   spqpt = coordinates of qpoints
!!   symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!!   trans(3,natom) = Atomic translations : xred = rcan + trans
!!   typat(natom)=type integer for each atom in cell
!!   ucvol=unit cell volume in bohr**3
!!   wghatm(natom,natom,nrpt) =Weight for the pair of atoms and the R vector
!!   xred(3,natom)=fractional dimensionless atomic coordinates
!!   zeff(3,3,natom) =effective charge on each atom, versus electric
!!        field and atomic displacement
!!
!! OUTPUT
!!   elph_ds = modified gkq
!!
!! NOTES
!!  inspired to some extent by epcouple.f from the DecAFT package by J. Kay Dewhurst
!!  most inputs taken from mkifc.f
!!  in anaddb set ifcflag 1 such that the IFC are calculated in atmfrc
!!    prior to calling elphon
!!
!! PARENTS
!!      get_all_gkk2
!!
!! CHILDREN
!!      canon9,chpev,ftgkk,inpphon,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine interpolate_gkk (acell,amu,atmfrc,dielt,dipdip,&
&      dyewq0,elph_ds,FSkptirred,FSkpt,ftwghtgkk,&
&      gmet,gprim,indsym,mpert,msym,natom,&
&      nrpt,nsym,ntypat,phon_ds,rcan,rmet,rprim,rprimd,rpt,spqpt,&
&      symrel,trans,typat,ucvol,wghatm,xred,zeff)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_17ddb, except_this_one => interpolate_gkk
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dipdip,mpert,msym,natom,nrpt,nsym,ntypat
 real(dp),intent(in) :: ucvol
 type(elph_type),intent(inout) :: elph_ds
 type(phon_type),intent(inout) :: phon_ds
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrel(3,3,nsym),typat(natom)
 real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt)
 real(dp),intent(in) :: FSkptirred(3,elph_ds%nFSkptirred),acell(3),amu(ntypat)
 real(dp),intent(in) :: atmfrc(2,3,natom,3,natom,nrpt),dielt(3,3)
 real(dp),intent(in) :: dyewq0(3,3,natom),ftwghtgkk(natom,nrpt),gmet(3,3)
 real(dp),intent(in) :: gprim(3,3),rcan(3,natom),rmet(3,3),rprim(3,3)
 real(dp),intent(in) :: rprimd(3,3),rpt(3,nrpt),spqpt(3,elph_ds%nqpt)
 real(dp),intent(in) :: trans(3,natom),wghatm(natom,natom,nrpt),xred(3,natom)
 real(dp),intent(in) :: zeff(3,3,natom)

!Local variables-------------------------------
  ! output variables for phfrq3
! variables for zhpev
! variables for phonon interpolation
!scalars
 integer :: dispindx,found,i1,i2,iFSkpt2,iFSqpt,iatom,ib1,ib2,idir,ier,ii,ikpt1
 integer :: ikpt2,iost,ipert1,ipert2,iqpt,iqpt1,iqpt2,isppol,nmaxint,qtor
 integer :: unit_gkkp
 real(dp) :: cpu,omegafactor,qphnrm,res,sumi,sumr,thisphfrq,wall
!arrays
 real(dp) :: displ(2,elph_ds%nbranch,elph_ds%nbranch),eigval(3*natom)
 real(dp) :: eigvec(3*3*natom*3*natom),pheigval(elph_ds%nbranch)
 real(dp) :: pheigvec(2*elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: phfrq_tmp(elph_ds%nbranch),qphon(3),qpt1(3),qpt2(3),redkpt(3)
 real(dp) :: tmpa(3),tmpai(3),tmpar(3),tmpkpt(3),tmpx(3),tmpxi(3),tmpxr(3)
 real(dp),allocatable :: gkk2_diag_tmp(:,:,:,:),gkk2_tmp(:,:,:,:,:,:,:)
 real(dp),allocatable :: gkk_tmp_full(:,:,:,:,:,:,:),matrx(:,:),zhpev1(:,:)
 real(dp),allocatable :: zhpev2(:)

! *************************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zhpev
#endif

!
!NOTE: mjv 18/5/2008 reverted to old style of ftgkk with all kpt done together.
!may want to modify this later to use the new cleaner format with 1 FT at a
!time.
!

 write (*,*) 'interpolate_gkk : enter'

 if (elph_ds%nsppol /= 1) then
  stop 'Error: interpolate_gkk not coded with nsppol>1 yet'
 end if
 isppol = 1


!------------------------------------------------------
!complete dynamical matrices for all qpts between points
!on full kpt grid (interpolation from IFC)
!------------------------------------------------------

!allocate (gkk_tmp(2,elph_ds%ngkkband,elph_ds%ngkkband,elph_ds%nbranch,elph_ds%nbranch,1,1))
!DEBUG
!allocate (gkk_tmp_full(2,elph_ds%ngkkband,elph_ds%ngkkband,elph_ds%nbranch,elph_ds%nFSband,elph_ds%nFSkpt))
!allocate (gkk_tmp_full(2,elph_ds%nbranch,elph_ds%nFSband,elph_ds%nFSband,elph_ds%nFSkpt))
!ENDDEBUG
 allocate (gkk2_tmp(2,elph_ds%ngkkband,elph_ds%ngkkband,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nFSkpt,1))
 allocate (gkk2_diag_tmp(elph_ds%ngkkband,elph_ds%ngkkband,elph_ds%nbranch,elph_ds%nFSkpt))

 allocate(zhpev1(2,2*3*natom-1),zhpev2(3*3*natom-2))
 allocate(matrx(2,(3*natom*(3*natom+1))/2))

 qphnrm = one
!in this part use the inverse Fourier transform to get 1 (arbitrary) qpt at a
!time
 ii = 0
 qtor = 0
 unit_gkkp = 150
 open (unit=unit_gkkp,file='gkkp_file_ascii',&
& form='formatted',&
& status='unknown',iostat=iost)
 if (iost /= 0) then
  write (*,*) 'interpolate_gkk : error opening gkkpfile as new'
  stop
 end if

!loop over all FS pairs.
!do ikpt1=1,elph_ds%nFSkptirred
!do iFSqpt=1,elph_ds%nFSkpt

!
!this should run through the sparse mesh of 2x2x2 kpoints
!
 do iFSqpt=1,elph_ds%nFSkpt
  res = 2.0_dp*(FSkpt(1,iFSqpt)+one)
  if (abs(res-int(res)) > tol10) cycle
  res = 2.0_dp*(FSkpt(2,iFSqpt)+one)
  if (abs(res-int(res)) > tol10) cycle
  res = 2.0_dp*(FSkpt(3,iFSqpt)+one)
  if (abs(res-int(res)) > tol10) cycle

! do ikpt1=1,1
! 
! NOTE: should be very easy to parallelize!
! 
! write (*,*) ' interpolate_gkk : ikpt1 = ',ikpt1, ' / ', elph_ds%nFSkptirred
  write (*,*) ' interpolate_gkk : ikpt1 = ',iFSqpt, ' / ', elph_ds%nFSkpt

! DEBUG
! write (*,*) ' interpolate_gkk : Warning debug version'
! cycle
! ENDDEBUG

  gkk2_tmp(:,:,:,:,:,:,:) = zero

! qphon = 1 - 2    ie.  1 = 2+qphon
  qphon(:) = FSkpt(:,iFSqpt)

! shouldnt be necessary here, but oh well
  call canon9(qphon(1),redkpt(1),res)
  call canon9(qphon(2),redkpt(2),res)
  call canon9(qphon(3),redkpt(3),res)

  qphon(:) = redkpt(:)
  redkpt(1) = qphon(1)*gprim(1,1)+qphon(2)*gprim(1,2)+qphon(3)*gprim(1,3)
  redkpt(2) = qphon(1)*gprim(2,1)+qphon(2)*gprim(2,2)+qphon(3)*gprim(2,3)
  redkpt(3) = qphon(1)*gprim(3,1)+qphon(2)*gprim(3,2)+qphon(3)*gprim(3,3)
  write (unit_gkkp,*) 'qp= ', redkpt

  call inpphon(displ,pheigval,pheigvec,phfrq_tmp,phon_ds,qphon)
  write (unit_gkkp,*) phfrq_tmp(:)*Ha_cmm1

  elph_ds%phfrq(:,iFSqpt) = phfrq_tmp
  ii = ii+1
! if (elph_ds%phfrqwrite == 1) then
! write (elph_ds%unitphfrq,REC=iFSqpt) phfrq_tmp
! if(ii > 0 .and. ii < 1000) write (*,'(a,i5,3E16.6,2x)') &
! &   ' wrote phfrq_tmp for time ', ii, phfrq_tmp
! end if

! phonon frequencies are in elph_ds%phfrq
! phonon eigenvectors are in eigvec
! real and imaginary parts
! phonon displacements = eigvec/sqrt(M_i) are in displ
! real and imaginary parts

! DEBUG
! test: uniform phonon frequency
! phfrq_tmp(:) = 0.0001_dp
! ENDDEBUG

! in this case we need to re calculate the gkk2
  if (elph_ds%gkk2exist == 1) cycle

! FT gamma matrices for all FSkpt points, and
! for qpoint = qphon(:) = FSkpt(iFSkpt)

  call ftgkk(wghatm,gkk2_tmp,elph_ds%gkk_rpt,elph_ds%gkqwrite,&
&  elph_ds%gkk_rptwrite,gprim,1,&
&  natom,elph_ds%nFSkpt,elph_ds%ngkkband,elph_ds%nFSkpt,1,nrpt,elph_ds%nsppol,&
&  qtor,rpt,qphon,elph_ds%unit_gkk_rpt,elph_ds%unitgkq)

! NOTE: Normally the eigenvectors of the gkk2_tmp should be the same as eigvec

! Diagonalize gamma matrices at qpoint (complex matrix) for all FSkpt.
! Copied from phfrq3
  do iFSkpt2=1,elph_ds%nFSkpt
   res = 8.0_dp*(FSkpt(1,iFSkpt2)+one)
   if (abs(res-int(res)) > tol10) cycle
   res = 8.0_dp*(FSkpt(2,iFSkpt2)+one)
   if (abs(res-int(res)) > tol10) cycle
   res = 8.0_dp*(FSkpt(3,iFSkpt2)+one)
   if (abs(res-int(res)) > tol10) cycle

   write (unit_gkkp,*) 'kp= ', FSkpt(:,iFSkpt2)

   do ib1=1,elph_ds%ngkkband
    do ib2=1,elph_ds%ngkkband
     ier=0
     ii=1
     do i2=1,3*natom
      do i1=1,i2
       matrx(1,ii)=gkk2_tmp(1,ib1,ib2,i1,i2,iFSkpt2,1)
       matrx(2,ii)=gkk2_tmp(2,ib1,ib2,i1,i2,iFSkpt2,1)
       ii=ii+1
      end do
     end do
#if defined T3E
     call CHPEV ('N','U',3*natom,matrx,eigval,eigvec,3*natom,zhpev1,&
&     zhpev2,ier)
#else
     call ZHPEV ('N','U',3*natom,matrx,eigval,eigvec,3*natom,zhpev1,&
&     zhpev2,ier)
#endif

     gkk2_diag_tmp(ib2,ib1,:,iFSkpt2) = eigval(:)
     do i1=1,3*natom
      write (unit_gkkp,*) elph_ds%minFSband-1+ib1,elph_ds%minFSband-1+ib2,i1,&
&      eigval(i1)
     end do
    end do
   end do
  end do

  if (elph_ds%gkk2write == 1) then
   write (*,*) 'WARNING COMMENTED WRITE TO BINARY FILE!!!'
!  write (elph_ds%unit_gkk2,REC=iFSqpt) gkk2_diag_tmp(:,:,:,:)
   write (*,'(a,i4,4(2E16.6,2x))') ' gkk2 loop ', &
&   iFSqpt,gkk2_diag_tmp(1,1,:,1:2),gkk2_diag_tmp(1,1,:,elph_ds%nFSkpt-1:elph_ds%nFSkpt)
!  &    ikpt1,gkk2_tmp(:,1,1,1,1,1:2),gkk2_tmp(:,1,1,elph_ds%nFSkpt-1:elph_ds%nFSkpt)
!  else if (elph_ds%gkk2write == 0 .and. elph_ds%gkk2exist == 0) then
  else if (elph_ds%gkk2write == 0) then
   elph_ds%gkk2(:,:,:,:,iFSqpt,isppol) = gkk2_diag_tmp(:,:,:,:)
!  elph_ds%gkk2(:,:,:,:,ikpt1) = gkk2_tmp
   write (*,*) ' interpolate_gkk : gkk2(b=1,b=1,:,kpt=1,iFSqpt) = '
   write (*,*) gkk2_diag_tmp(1,1,:,1)
  end if

 end do
!end do on iFSqpt

 deallocate(matrx,zhpev1,zhpev2)

!
!gkk2 have been written to file
!
 if (elph_ds%gkk2write == 1) then
  elph_ds%gkk2exist = 1
 end if

end subroutine interpolate_gkk
!----------------------------------------------------------
!!***
