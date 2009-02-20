!{\src2tex{textfont=tt}}
!!****f* ABINIT/read_gkk
!!
!! NAME
!! read_gkk
!!
!! FUNCTION
!! This routine reads in elphon matrix elements and completes them
!!  using the appropriate symmetries
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer, MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  amu = masses of atoms
!!  elph_ds = datastructure containing elphon matrix elements
!!  FSfullpqtofull = mapping of k+q to k
!!  FSfulltofull = correspondance between FS kpoints by syms
!!  FSintweight = integration weights on fermi surface for all kpts and bands
!!  FSkpt = fermi surface kpoints
!!  gprimd = reciprocal lattice vectors (dimensionful)
!!  indsym = map of atoms through symrel
!!  n1wf = number of 1WF files to be read and analyzed
!!  natom = number of atoms
!!  nband = number of bands per kpoint
!!  nsym = number of symmetries for full lattice
!!  ntypat = number of types of atoms
!!  phon_ds = phonon datastructure with real space interatomic force constants
!!  rprimd = real space primitive vectors of cell (trans to cart coord)
!!  spqpt = coordinates of full uniform qpt grid
!!  symrec = symmetry operations in reciprocal space
!!  symrel = symmetry operations in real space
!!  timrev = flag for using time reversal symmetry. No effect for the
!!     moment: info from symq
!!  tnons = translations associated to symrel
!!  typat = array of types of atoms
!!  ucvol = unit cell volume
!!  unitgkk = unit of GKK file for reading
!!
!! OUTPUT
!!  elph_ds = modified gkq
!!  gkk_qpt = el-ph matrix elements for irreducible qpoints and
!!    kpoints (as a function of the reduced symmetry for the qpoint)
!!  gkk_flag = flag array:
!!       -1 -> element is missing
!!        0 -> element is from symmetric qpt (Now done in complete_gkk)
!!        1 -> element is from symmetric pert
!!        2 -> element is kptsym of gkk file
!!        3 -> element was read from gkk file
!!  nqptirred = number of irreducible qpoints
!!  qptirred = irreducible qpoints found in gkk file
!!
!! NOTES
!!
!! PARENTS
!!      get_all_gkq
!!
!! CHILDREN
!!      canon9,completeperts,hdr_clean,hdr_io,inpphon,insy3,leave_new,mati3inv
!!      normsq_gkq,symq3,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine read_gkk(amu,elph_ds,FSfullpqtofull,FSfulltofull,FSintweight,FSkpt,    &
&                   gkk_flag,gprimd,indsym,n1wf,natom,nband,nqptirred,nsym,ntypat,&
&                   phon_ds,qptirred,rprimd,spqpt,symrec,symrel,timrev,tnons,typat,ucvol,unitgkk)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13io_mpi
 use interfaces_13recipspace
 use interfaces_14iowfdenpot
 use interfaces_16response
 use interfaces_17ddb, except_this_one => read_gkk
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1wf,natom,nband,nsym,ntypat,timrev,unitgkk
 integer,intent(out) :: nqptirred
 real(dp),intent(in) :: ucvol
 type(elph_type),intent(inout) :: elph_ds
 type(phon_type),intent(inout) :: phon_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
 integer,intent(in) :: FSfulltofull(2,nsym,elph_ds%nFSkpt),indsym(4,nsym,natom)
 integer,intent(in) :: symrec(3,3,nsym),symrel(3,3,nsym),typat(natom)
 integer,intent(out) :: gkk_flag(elph_ds%nbranch,elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol,elph_ds%nqpt)
 real(dp),intent(in) :: FSintweight(elph_ds%nFSband,elph_ds%nFSkpt)
 real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt),amu(ntypat),gprimd(3,3)
 real(dp),intent(in) :: rprimd(3,3),tnons(3,nsym)
 real(dp),intent(inout) :: spqpt(3,elph_ds%nqpt)
 real(dp),intent(out) :: qptirred(3,n1wf)

!Local variables-------------------------------
! 4.3.2004   gkq becomes a complex valued matrix and is summed over bands
!real(dp) :: gkk_qpt_tmp (2,elph_ds%nbranch,elph_ds%nFSband,elph_ds%nFSband,elph_ds%nFSkpt)
!real(dp) :: gkk_qpt_tmp2(2,elph_ds%nbranch,elph_ds%nFSband,elph_ds%nFSband,elph_ds%nFSkpt)
! local variables for insy3
!real(dp) :: ssym(3,3,nsym)
! output variables for inpphon
! real(dp) :: h1_mat_el(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
! real(dp) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
!
!scalars
 integer :: fform,gkqwrite,goodkpq,i1,i1wf,iFSkpt,iFSkptq,iatom,iatom1,iatom2
 integer :: ib,ib1,ib2,ibb,ibranch,idir,idir1,idir2,idir3,ierr,ii,ikpt1,index
 integer :: ipert,ipert1,ipert2,iqpt,iqptfull,isppol,isym,isym1,isymq1,itim
 integer :: itim1,jFSkpt,jatom,jbranch,jdir,jj,jpert,k1,kbranch,kdir,kk,ll,new
 integer :: nsym1,qtimrev,rdwr,symiFSkpt,sympertcase1,sympertcase2,syuse
 integer :: tdonecompl,test_flag,verify
 real(dp) :: eigentol,eigentol2,exparg,res,s1,s2,s3,ss,sumi,sumr,timsign
 character(len=4) :: qptstr
 character(len=500) :: message
 character(len=fnlen) :: xsffilnam
 type(hdr_type) :: hdr1
!arrays
 integer :: FSirrtok(3,elph_ds%nFSkpt),dummysymafm(nsym)
 integer :: irredpert(7,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nqpt)
 integer :: qptirredtofull(elph_ds%nqpt),symaf1(nsym),symq(4,2,nsym)
 integer :: symrc1(3,3,nsym),symrl1(3,3,nsym),symvdir(3)
 integer :: tmpflg(3,natom+2,3,natom+2),tmpflg_init(3,natom+2,3,natom+2)
 integer :: vdir(3)
 real(dp) :: cosarg(2,nsym),diffsymgkk(2,nsym),displ(2,3*natom,3*natom)
 real(dp) :: displ_red(2,3*natom,3*natom),eigval(3*natom)
 real(dp) :: eigvec(2,3*natom,3*natom),kpt(3),phfrq_tmp(3*natom),redkpt(3)
 real(dp) :: sinarg(2,nsym),tmp_mat(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_mat2(2,elph_ds%nbranch,elph_ds%nbranch),tmpa(3)
 real(dp) :: tmpval(2,3,natom+2,3,natom+2),tmpx(3),tnons1(3,nsym)
 real(dp) :: wf(elph_ds%nbranch)
 real(dp),allocatable :: accum_eigen1(:,:,:),eigen1(:,:,:),gkk_qpt_tmp(:,:,:,:)
 real(dp),allocatable :: h1_mat_el(:,:,:,:,:),h1_mat_el_sq(:,:,:,:,:)
 real(dp),allocatable :: qdata(:,:,:),qdata_tmp(:,:,:,:)

! *************************************************************************

!DEBUG
!write(*,*)' read_gkk : enter '
!write(*,*)' read_gkk : nband = ',nband
!write(*,*)' read_gkk : rprimd', rprimd
!write(*,*)' read_gkk : gprimd', gprimd
!write(*,*)' read_gkk : 1/sqrt(amu(typat(ipert))*amu_emass) = ',&
!&   one /sqrt(amu(typat(:))*amu_emass)
!write (127,*) 'read_gkk : FSfulltofull = '
!do iFSkpt=1,elph_ds%nFSkpt
!write (127,*) 'iFSkpt = ',iFSkpt
!write (127,'(4(2i5,4x))') (FSfulltofull(:,isym,iFSkpt),isym=1,nsym)
!end do
!ENDDEBUG

 allocate(h1_mat_el(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol),stat=ierr)
 if (ierr /= 0 ) then
  write (message,'(3a)')' read_gkk : ERROR- ',ch10,&
&  ' trying to allocate array h1_mat_el '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 h1_mat_el(:,:,:,:,:) = zero

 allocate(h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol),stat=ierr)
 if (ierr /= 0 ) then
  write (message,'(3a)')' read_gkk : ERROR- ',ch10,&
&  ' trying to allocate array h1_mat_el_sq '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!Tolerance for eigenvalues from 1WF files being non-zero
 eigentol = 1.0d-50
 eigentol2 = 2.0d-50

!MG array to store the e-ph quantities calculated over the input Q-grid
 allocate (qdata_tmp(elph_ds%nqpt,elph_ds%nbranch,elph_ds%nsppol,3))
 qdata_tmp(:,:,:,:)=zero

 nqptirred=0 !zero number of irred q-points found
 qptirred(:,:)=zero

 gkk_flag(:,:,:,:,:) = -1

 if (elph_ds%gkqwrite ==0) then
  elph_ds%gkk_qpt(:,:,:,:,:,:) = zero
 else if (elph_ds%gkqwrite == 1) then
  allocate(gkk_qpt_tmp(2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt),stat=ierr)
  if (ierr /= 0 ) then
   write (message,'(3a)')' read_gkk : ERROR- ',ch10,&
&   ' trying to allocate array gkk_qpt_tmp '
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if
  gkk_qpt_tmp(:,:,:,:) = zero
  do iqpt=1,elph_ds%nqpt
   write (*,*) ' initialized for ', iqpt
   write (elph_ds%unitgkq,REC=iqpt) gkk_qpt_tmp
  end do
  deallocate (gkk_qpt_tmp)
 else
  write (message,'(3a,i3)')&
&  ' read_gkk : BUG-',ch10,&
&  ' Wrong values for gkqwrite = ',elph_ds%gkqwrite
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if !gkqwrite

 allocate(eigen1(2,nband,nband),accum_eigen1(2,nband,nband))

 irredpert(:,:,:,:) = -999
 h1_mat_el(:,:,:,:,:) = zero

!===========================================================
!Loop over all files we have
!read in header for perturbation
!should check that all files are complete, have same header
!(taking into account the symmetries for the qpoint),
!represent the correct qpoints ...
!MG: this task should be performed in mrggkk
!===========================================================

 do i1wf=1,n1wf

  write (message,'(2a,i4,a,i4)')ch10,' read_gkk : reading 1WF header # ',i1wf,' /',n1wf
  call wrtout(06,message,'COLL')

! Could check for compatibility of natom, kpt grids, ecut, qpt with DDB grid...
! MG: Also this task should be done in mrggkk

  rdwr = 5 !read without rewinding
  call hdr_io(fform,hdr1,rdwr,unitgkk)
  if (fform == 0) then
   write (message,'(4a,i5,a)')ch10,' read_gkk : ERROR- :',ch10,&
&   ' 1WF header number ',i1wf,' was mis-read. fform == 0'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  write(message,'(a,i4)')' read_gkk : have read 1WF header #',i1wf
  call wrtout(06,message,'COLL')

! Find qpoint in grid
  new=1
  do iqptfull=1,elph_ds%nqpt
   kpt(:) = hdr1%qptn(:) - spqpt(:,iqptfull)
   call canon9(kpt(1),redkpt(1),res)
   call canon9(kpt(2),redkpt(2),res)
   call canon9(kpt(3),redkpt(3),res)
   ss=redkpt(1)**2+redkpt(2)**2+redkpt(3)**2
   if(ss < tol6) then
    new = 0
    exit !exit with iqtfull
   end if
  end do !iqptfull

  if (new == 1) then
!  Test should be at the end: dont care if there are additional
!  qpts in gkk file which are not on the main grid. Ignore them.
   write (message,'(4a,3es16.6,2a)')ch10,&
&   ' read_gkk : WARNING-  ',ch10,&
&   ' qpoint = ',hdr1%qptn(:),ch10,&
&   ' not found in the input q-grid. Ignoring this point '
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,message,'COLL')

   do isppol=1,hdr1%nsppol
    do ikpt1=1,hdr1%nkpt
     read(unitgkk) ((eigen1(:,ii,ib),ii=1,nband),ib=1,nband)
    end do
   end do
   cycle !cycle the loop on i1wf
  end if !end if (new ==1)

! Check whether other pieces of the DDB have used this qpt already
  new=1
  do iqpt=1,nqptirred
   kpt(:) = qptirred(:,iqpt) - hdr1%qptn(:)
   call canon9(kpt(1),redkpt(1),res)
   call canon9(kpt(2),redkpt(2),res)
   call canon9(kpt(3),redkpt(3),res)
   ss=redkpt(1)**2+redkpt(2)**2+redkpt(3)**2
   if(ss < tol6) then
    new=0
    exit  !MG We can use this information to avoid recalculating the dynamical matrix 
   end if !but we need to use a fixed format in GKK!
  end do !iqpt

  if (new==1) then  !we have a new valid irreducible qpoint, add it!
   nqptirred = nqptirred+1
   qptirred(:,nqptirred) = hdr1%qptn(:)
   iqpt = nqptirred
   tdonecompl = 0
   h1_mat_el(:,:,:,:,:) = zero
  end if

! now iqpt is the index of the present qpoint in the array qptirred
! and iqptfull is the index in the full spqpt array for future reference
  qptirredtofull(iqpt) = iqptfull

  write (message,'(a,i5,a,3es16.8)')&
&  ' read_gkk : full zone qpt number ',iqptfull,' is ',spqpt(:,iqptfull)
  call wrtout(06,message,'COLL')

! if this perturbation has already been filled (overcomplete gkk)
! check only 1st kpoint and spinpol, then check others
  verify = 0
  if (gkk_flag(hdr1%pertcase,hdr1%pertcase,1,1,iqptfull) /= -1) then
   do isppol=1,elph_ds%nsppol
    do iFSkpt=1,elph_ds%nFSkpt
     if (gkk_flag(hdr1%pertcase,hdr1%pertcase,iFSkpt,isppol,iqptfull) == -1) then
      write (message,'(4a)')ch10,&
&      ' read_gkk : ERROR-',ch10,&
&      ' partially filled perturbation '
      call wrtout(06,message,'COLL')
      write (*,*)hdr1%pertcase,iFSkpt,iqptfull
      call leave_new('COLL')
     end if
    end do ! iFSkpt
   end do ! isppol
   write(message,'(4a)')ch10,&
&   ' read_gkk : WARNING- ',ch10,&
&   ' gkk perturbation is already filled'
   call wrtout(06,message,'COLL')
   write(*,*)' hdr1%pertcase,iqptfull = ',hdr1%pertcase,iqptfull,&
&   gkk_flag(hdr1%pertcase,hdr1%pertcase,1,1,iqptfull)
   verify = 1
   write (125,*) '# matrix elements for symmetric perturbation'
!  Instead of reading eigen1 into void, verify == 1 checks
!  them later on wrt values in memory
!  do isppol=1,hdr1%nsppol
!  do ikpt1=1,hdr1%nkpt
!  read (unitgkk) ((eigen1(:,ii,ib),ii=1,nband),ib=1,nband)
!  end do
!  end do
!  cycle
  end if !gkk_flag

! Examine the symmetries of the q wavevector
! these will be used to complete the perturbations for other atoms and idir

! call symq3(nsym,qptirred(:,nqptirred),symq,symrec,qtimrev)
  call symq3(nsym,spqpt(:,iqptfull),symq,symrec,qtimrev)

! DEBUG
! write(*,*)'read_gkk : after symq3, timrev = ', qtimrev
! write(*,*)'read_gkk : symq = '
! write(*,'(3i7)') (isym,(symq(4,itim,isym),itim=1,2),isym=1,nsym)
! ENDDEBUG

! Determine dynamical matrix, phonon frequencies and displacement vector for qpoint
! MG: For each qpt we are calculating the dynamical matrix for each perturbation
! Possible solution: use a fixed format for the GKK file and a flag which
! calls inpphon only for the first perturbation of each q point
  
  write (message,'(2a)')ch10,' read_gkk : calling inpphon to calculate the dynamical matrix'
  call wrtout(6,message,'COLL')

  call inpphon(displ,eigval,eigvec,phfrq_tmp,phon_ds,spqpt(:,iqptfull))

! DEBUG
! write(*,*)'read_gkk : eigvec = '
! write(*,'(3(2e16.6,1x))')eigvec
! write(*,*)'read_gkk : displ = '
! write(*,'(3(2e16.6,1x))')displ
! ENDDEBUG

! Get displacement vectors for all branches in reduced coordinates
! used in scalar product with H(1)_atom,idir  matrix elements
! Calculate $displ_red = displ \cdot gprimd$ for each phonon branch

  displ_red(:,:,:) = zero
  do jbranch=1,elph_ds%nbranch
   do iatom=1,natom
    do idir=1,3
     ibranch=idir+3*(iatom-1)
     do kdir=1,3
      k1 = kdir+3*(iatom-1)
!     WARNING: could be non-transpose of rprimd matrix : to be checked.
!     23 june 2004: rprimd becomes gprimd
!     could be gprim and then multiply by acell...
!     Nope, checked and ok with gprimd 24 jun 2004
      displ_red(1,ibranch,jbranch) = displ_red(1,ibranch,jbranch) + &
&      gprimd(kdir,idir)*displ(1,k1,jbranch)

      displ_red(2,ibranch,jbranch) = displ_red(2,ibranch,jbranch) + &
&      gprimd(kdir,idir)*displ(2,k1,jbranch)

!     DEBUG
!     write(*,*)'ibranch,jbranch,idir,kdir,k1,rprimd(kdir,idir)'
!     write(*,*)ibranch,jbranch,idir,kdir,k1,rprimd(kdir,idir)
!     write(*,*)'displ(1,k1,jbranch),displ_red(1,ibranch,jbranch)'
!     write(*,*)displ(1,k1,jbranch),displ_red(1,ibranch,jbranch)
!     ENDDEBUG
     end do !kdir
    end do !idir
   end do !iatom

  end do !jbranch

  accum_eigen1(:,:,:) = zero

! prefactors for gk+q,n\prime;k,n matrix element
! COMMENT : in decaft there is a weird term in the mass factor, of M-zval(species)
! dont know why. Not needed to reproduce decaft results, though... 
! weight is squared in evaluation of
! gamma_{k,q,j} = 2 \pi omega_{q,j} sum_{nu,nu\prime} |g^{q,j}_{k+q,nu\prime; k,nu}|^2
! normally cancels with the 2 \pi omega_{q,j} factor in front of the sum...

  do jatom=1,natom
   do jdir=1,3
    jbranch=3*(jatom-1)+jdir

!   WARNING : the tolerance for which a phonon freq is 0 is a bit arbitrary.
!   Should be smaller.
    if (abs(phfrq_tmp(jbranch)) > tol10) then
!    It looks like the 1/sqrt(M) factor is already in displ from phfrq3
!    wf(jbranch) = one/sqrt(two*amu(typat(jatom))*amu_emass*phfrq_tmp(jbranch))
!    wf(jbranch) = one/sqrt(two*phfrq_tmp(jbranch))
!    wf(jbranch) = one/(two*phfrq_tmp(jbranch))
!    wf(jbranch) = one
!    wf(jbranch) = pi
!    Here add a factor of BZ_volume = 1/ucvol and a factor of pi for the linewidths
!    wf(jbranch) = pi/ucvol
!    Use same prefactor as decaft.
!    Is coherent with standard definition of g matrix elements and gamma linewidths
!    Not used in case doscalprod==0
     wf(jbranch) = one/sqrt(two*abs(phfrq_tmp(jbranch)))
    else
     wf(jbranch) = zero
    end if
   end do !jdir
  end do !jatom

! DEBUG
! write(*,'(a,3e16.6)') ' qpt ', spqpt(:,iqptfull)
! write(*,*) ' phfrq_tmp = '
! write(*,'(3e16.6)')  phfrq_tmp
! write(*,*) ' weights = '
! write(*,'(3e16.6)')  wf
! ENDDEBUG

! hdr1%pertcase = idir + (ipert-1)*3
! where ipert=iatom in the interesting cases
  idir = mod (hdr1%pertcase-1,3)+1
! vdir is the a,b, or c basis vector, for transformation under symops
  vdir(:) = 0
  vdir(idir) = 1
  ipert = int(dble(hdr1%pertcase-idir)/three)+1

  write (message,'(4a,i3,a,i3,a,i4,a)')ch10,&
&  ' read_gkk : calling insy3 to examine the symmetries of the full perturbation ',ch10,&
&  ' idir = ',idir,' ipert = ',ipert,' and Q point = ',iqptfull,ch10
  call wrtout(06,message,'COLL') 

! MG use this because insy3 writes to ab_out
! write (message,'(a,i3,a,i3,a,i4)')&
! & ' For perturbation : idir = ',idir,' ipert = ',ipert,' and Q point = ',iqptfull
! call wrtout(ab_out,message,'COLL')

! Examine the symmetries of the full perturbation these will be used to complete the kpoints
! DOESNT USE TIME REVERSAL IN insy3 except for gamma

  syuse = 0
  dummysymafm(:) = 1
  call insy3(gprimd,idir,indsym,ipert,natom,nsym,nsym1,2,dummysymafm,symaf1,&
&  symq,symrec,symrel,symrl1,syuse,tnons,tnons1)
  do isym1=1,nsym1
   call mati3inv(symrl1(:,:,isym1),symrc1(:,:,isym1))
  end do
! DEBUG
! write(*,*)' read_gkk : hdr1%pertcase ', hdr1%pertcase
! write(*,*)' read_gkk : nsym1 = ',nsym1
! write(*,*)' read_gkk : symrl1 = '
! write(*,'(3(3I4,3x))') symrl1(:,:,1:nsym1)
! write(*,*) ' read_gkk : symrc1 = '
! write(*,'(3(3I4,3x))') symrc1(:,:,1:nsym1)
! write(*,*) ' read_gkk : tnons1 = '
! write(*,'(3(3E16.6,3x))') tnons1(:,1:nsym1)
! write(*,*) ' read_gkk :spqpt  = '
! write(*,'(3E16.6,3x)') spqpt(:,iqptfull)
! ENDDEBUG

! !For present symmetries, pre-calculate matrices to transform
! !cartesian coordinates from one kpt to another
! !SHOULD NEED ONLY RED COORD symrl1 MATRICES
! ssym(:,:,:)=zero
! do isymq1=1,nsym1
! do ii=1,3
! do jj=1,3
! ssym(ii,jj,isymq1)=zero
! do kk=1,3
! do ll=1,3
! ssym(ii,jj,isymq1)=ssym(ii,jj,isymq1)+&
! &      rprimd(ii,kk)*symrl1(kk,ll,isymq1)*gprimd(jj,ll)
! end do
! end do
! end do
! end do
! end do

! MG06062006 NOTE: these quantities are not used any more
! For present symmetries, pre-calculate the phase factors for
! transforming matrix elements between kpoints
! NOT TESTED YET! MAY NEED MUCH MORE COMPLICATED SCHEME
! IF THE PHASE ENTERS IN THE WFK AS EXP(-TNON*(G+q))
! AND DOESNT FACTORIZE
! 17 May 2004 ... but this is used in pwscf as well!
  do isym1=1,nsym1
   exparg = tnons1(1,isym1)*spqpt(1,iqptfull)+ &
&   tnons1(2,isym1)*spqpt(2,iqptfull)+ &
&   tnons1(3,isym1)*spqpt(3,iqptfull)
   exparg = -two_pi*exparg
   cosarg(1,isym1) = cos(exparg)
   sinarg(1,isym1) = sin(exparg)
   cosarg(2,isym1) = cos(exparg)
   sinarg(2,isym1) =-sin(exparg)
  end do

  FSirrtok (:,:) = 0

! ========================================================
! Loop over irred kpts in file, and fill the default gkk
! ========================================================

! MG this array is not used anymore
  diffsymgkk(:,:) = zero

! MG NOTE : in the present implementation, if nsppol /=1 the code stops in rchkGSheader!
  do isppol=1,hdr1%nsppol !Loop over spins is trivial? Not tested.
   write (*,*) ' read_gkk : isppol = ', isppol 
   do ikpt1=1,hdr1%nkpt   !Loop over irred kpoints, WARNING  nkpt depends on qpoint and symmetry!

!   
!   this is the main read of the gkk matrix elements from the file (eigen1 arrays)
!   it has to be done exactly nsppol*nkpt times, and the FSkpt are completed
!   where appropriate in the loop below (normally succeeding only once for each kpt)
!   
    read(unitgkk) ((eigen1(:,ii,ib),ii=1,nband),ib=1,nband)

!   Check to see if kpoint is in FS set
!   WARNING! the kpoints in the file (kptns) could be ordered arbitrarily
    do iFSkpt=1,elph_ds%nFSkpt
     kpt(:) = hdr1%kptns(:,ikpt1)-FSkpt(:,iFSkpt)-spqpt(:,iqptfull)
     call canon9(kpt(1),redkpt(1),res)
     call canon9(kpt(2),redkpt(2),res)
     call canon9(kpt(3),redkpt(3),res)

     ss=redkpt(1)**2+redkpt(2)**2+redkpt(3)**2

!    this is not the point on the FS we are looking for
     if (ss > tol6) cycle

     if (verify == 1) then
      do ib1=1,elph_ds%nFSband
       do ib2=1,elph_ds%nFSband
        ibb = (ib1-1)*elph_ds%nFSband+ib2
        write (125,'(2(2E16.6,2x))') h1_mat_el(:,ibb,hdr1%pertcase,iFSkpt,isppol),&
&        eigen1(:,elph_ds%minFSband-1+ib2,elph_ds%minFSband-1+ib1)
       end do
      end do
!     DEBUGISSIMO
!     do isym=1,nsym
!     do itim=0,1
!     write (126,*) ' trying sym ', isym, iqptfull, sympertcase1
!     symiFSkpt=FSfulltofull(itim+1,isym,iFSkpt)
!     do ib1=1,elph_ds%nFSband
!     do ib2=1,elph_ds%nFSband
!     ibb = (ib1-1)*elph_ds%nFSband+ib2
!     diffsymgkk(itim+1,isym) = diffsymgkk(itim+1,isym) + &
!     &      abs(h1_mat_el(1,ibb,hdr1%pertcase,symiFSkpt,isppol)**2&
!     &      +h1_mat_el(2,ibb,hdr1%pertcase,symiFSkpt,isppol)**2&
!     &      -eigen1(1,elph_ds%minFSband-1+ib2,elph_ds%minFSband-1+ib1)**2&
!     &      -eigen1(2,elph_ds%minFSband-1+ib2,elph_ds%minFSband-1+ib1)**2)
!     write (126,'(2E16.6)') h1_mat_el(:,ibb,hdr1%pertcase,symiFSkpt,isppol)&
!     &      -eigen1(:,elph_ds%minFSband-1+ib2,elph_ds%minFSband-1+ib1)
!     end do
!     end do
!     end do
!     end do
!     ENDDEBUGISSIMO
     end if !verify

!    if this kpoint has already been filled (overcomplete gkk)
     if (gkk_flag(hdr1%pertcase,hdr1%pertcase,iFSkpt,isppol,iqptfull) /= -1) then
      write(*,*)' read_gkk warning : gkk element is already filled'
      write (*,*)' hdr1%pertcase,iFSkpt,isppol,iqptfull = ',&
&      hdr1%pertcase,iFSkpt,isppol,iqptfull,&
&      gkk_flag(hdr1%pertcase,hdr1%pertcase,iFSkpt,isppol,iqptfull)
      exit
     end if !gkk_flag

!    save this kpoint
     do ib1=1,elph_ds%nFSband
      do ib2=1,elph_ds%nFSband
       ibb = (ib1-1)*elph_ds%nFSband+ib2
       
!      17.05.04 corrected ib1ib2 order in eigen1 indices...

!      real
       res=eigen1(1,elph_ds%minFSband-1+ib2,elph_ds%minFSband-1+ib1)
!      if (abs(res*1.0d10) > eigentol2) then
       h1_mat_el(1,ibb,hdr1%pertcase,iFSkpt,isppol) = res
!      end if

!      imag
       res=eigen1(2,elph_ds%minFSband-1+ib2,elph_ds%minFSband-1+ib1)
!      if (abs(res*1.0d10) > eigentol2) then
       h1_mat_el(2,ibb,hdr1%pertcase,iFSkpt,isppol) = res
!      end if
      end do !ib2
     end do !ib1
     gkk_flag(hdr1%pertcase,hdr1%pertcase,iFSkpt,isppol,iqptfull) = 3

!    DEBUG
!    if (abs(spqpt(1,iqptfull))+abs(spqpt(2,iqptfull)-half)+&
!    &                     abs(spqpt(3,iqptfull)-half) < tol6) then
!    !write (*,'(a)') 'read_gkk : eigen1 = '
!    !write (*,'(3(2E18.6,1x))') ((eigen1(:,ii,jj),ii=1,3),jj=1,3)
!    accum_eigen1(:,:,:) = accum_eigen1(:,:,:) + eigen1(:,:,:)
!    end if
!    ENDDEBUG

!    ===============================================================
!    we now have contribution to g(k+q,k; \kappa,\alpha) from one
!    kpoint,and one perturbation,
!    NB: each perturbation will contribute to all the modes later!
!    find correspondence between this FSkpt and the others
!    provided sym conserves pert as well as qpoint: add to FSirrtok
!    
!    SHOULD ONLY DO THIS FOR THE SYMS NEEDED 
!    TO COMPLETE THE PERTURBATIONS!!!
!    ================================================================

     new=0
     do isym1=1,nsym1
      do itim1=0,qtimrev
       timsign=one-two*itim1
       kpt(:) = timsign*(symrc1(:,1,isym1)*FSkpt(1,iFSkpt)+&
&       symrc1(:,2,isym1)*FSkpt(2,iFSkpt)+&
&       symrc1(:,3,isym1)*FSkpt(3,iFSkpt))
       call canon9(kpt(1),redkpt(1),res)
       call canon9(kpt(2),redkpt(2),res)
       call canon9(kpt(3),redkpt(3),res)
       new=1
       do jFSkpt=1,elph_ds%nFSkpt
        ss=  (redkpt(1)-FSkpt(1,jFSkpt))**2&
&        +(redkpt(2)-FSkpt(2,jFSkpt))**2&
&        +(redkpt(3)-FSkpt(3,jFSkpt))**2
        if (ss < tol6) then
         new=0
         FSirrtok(1,jFSkpt) = iFSkpt
         FSirrtok(2,jFSkpt) = isym1
         FSirrtok(3,jFSkpt) = itim1
         exit !exit jFSkpt
        end if
       end do !jFSkpt
       if (new == 1) then
        write (message,'(4a,3es16.6,a,i5,a,i4,a)')ch10,&
&        ' read_gkk  ERROR- :',ch10,&
&        ' equivalent of kpt ',FSkpt(:,iFSkpt),' by sym ',isym1,' and itime ',itim1,' was not found'
        call wrtout(06,message,'COLL')
        call leave_new('COLL')
       end if
      end do !itim1
     end do !isim1

!    DEBUG
!    write(*,*)'iFSkpt = ', iFSkpt,' FSirrtok = '
!    write(*,'(3i12)')FSirrtok
!    ENDDEBUG

!    can we safely exit here? The ikpt has been identified to iFSkpt,
!    so no need to search the other FSkpt
     exit
    end do !iFSkpt
   end do !ikpt1

!  DEBUG
!  if (abs(spqpt(1,iqptfull)-half)+abs(spqpt(2,iqptfull)-half)+&
!  &abs(spqpt(3,iqptfull)) < tol6) then
!  write(*,'(a)') 'read_gkk : accum_eigen1 = '
!  write(*,'(3(2E18.6,1x))') ((accum_eigen1(:,ii,jj),ii=1,nband),jj=1,nband)
!  end if
!  ENDDEBUG
  end do !isppol

  if (verify == 1) then
   cycle
  end if

! check if irred kpoints found do reconstitute the FS kpts
  do iFSkpt=1,elph_ds%nFSkpt
   if (FSirrtok(1,iFSkpt) == 0) then
    write(message,'(4a,3es16.6,2a)')ch10,                        &
&    ' read_gkk : ERROR-',ch10,                                   &
&    ' kpt = ',FSkpt(:,iFSkpt),ch10,                              &
&    ' is not the symmetric of one of those found in the GKK file'
    call wrtout(06,message,'COLL')
    call leave_new('COLL')
   end if
  end do !iFSkpt

! ===============================================================
! We now have all irred kpoints : complete the others
! complete gkk for symmetric iFSkpt with sym1 which conserve
! the full perturbation+qpoint
! Not tested explicitly, but the results for Pb using reduced kpts look good
! should do same RF calculation with nsym=1 and check
! ===============================================================

  do symiFSkpt=1,elph_ds%nFSkpt
!  if the element has already been filled with another sym op, cycle
   if (gkk_flag(hdr1%pertcase,hdr1%pertcase,symiFSkpt,1,iqptfull) /= -1) cycle

   iFSkpt=FSirrtok(1,symiFSkpt)
   isym1 =FSirrtok(2,symiFSkpt)
   itim1 =FSirrtok(3,symiFSkpt)
   timsign = one-two*itim1

!  DEBUG
!  write(*,*)'New sym kpoint : iFSkpt,symiFSkpt,isym1,itim1,timsign '
!  write(*,*)iFSkpt,symiFSkpt,isym1,itim1,timsign
!  ENDDEBUG

!  copy kpt to symmetric kpoint
   do ibb=1,elph_ds%nFSband*elph_ds%nFSband
    h1_mat_el(1,ibb,hdr1%pertcase,symiFSkpt,:) = h1_mat_el(1,ibb,hdr1%pertcase,iFSkpt,:)
    h1_mat_el(2,ibb,hdr1%pertcase,symiFSkpt,:) = h1_mat_el(2,ibb,hdr1%pertcase,iFSkpt,:)
   end do
   gkk_flag(hdr1%pertcase,hdr1%pertcase,symiFSkpt,:,iqptfull) = 2
  end do !symiFSkpt

! normally at this point we have used all the gkk for all kpoints on the FS
! for the given irred perturbation: check
  do iFSkpt=1,elph_ds%nFSkpt
   if (gkk_flag(hdr1%pertcase,hdr1%pertcase,iFSkpt,1,iqptfull) == -1) then
    write (message,'(3a,i3,a,3es18.6,2a,i3,a,i3,a,3es18.6,a,a,i4,a,a)')&
&    ' read_gkk  ERROR- :',ch10,          &
&    ' For full qpt ', iqptfull,') ',spqpt(:,iqptfull),ch10,                                         &
&    ' the gkk element : pertcase = ',hdr1%pertcase,' ikpt = ',iFSkpt,' kpt = ',FSkpt(:,iFSkpt),ch10,&
&    ' and isppol ',1,ch10,&
&    ' was not found by symmetry operations on the irreducible kpoints given'
    call wrtout(06,message,'COLL')
    call leave_new('COLL')
   end if
   if (FSirrtok(1,iFSkpt) == 0) then
    write (message,'(3a)')' read_gkk : ERROR-',ch10,&
    ' One of the kpoints was not equivalent to an irred one found in the gkk file'
    call wrtout(06,message,'COLL')
    call leave_new('COLL')
   end if
  end do !iFSkpt

  write(message,'(a,i4)')' read_gkk : Done completing the kpoints for pert ',hdr1%pertcase
  call wrtout(06,message,'COLL')

  tmpflg(:,:,:,:) = 0

  do idir1=1,3
   do iatom1=1,natom
    ipert1 = (iatom1-1)*3+idir1
    do idir2=1,3
     do iatom2=1,natom
      ipert2 = (iatom2-1)*3+idir2
      if (gkk_flag(ipert1,ipert1,1,1,iqptfull) >= 0 .and. gkk_flag(ipert2,ipert2,1,1,iqptfull) >= 0) then
       tmpflg(idir1,iatom1,idir2,iatom2) = 1
      end if
     end do
    end do
   end do
  end do


! DEBUG
! write (*,*) 'gkk_flag(:,:,1,1,iqptfull) = ', gkk_flag(:,:,1,1,iqptfull)
! ENDDEBUG

! DEBUG
! call d2sym3 to check whether all perts can be reconstituted.
! NOTE: now called from within completeperts.F90
! tmpval(:,:,:,:,:) = zero  ! needs to be filled
! write (*,*) 'read_gkk : tmpflg before ', tmpflg
! call d2sym3(tmpflg,tmpval,indsym,natom+2,&
! &     natom,nsym,spqpt(:,iqptfull),symq,symrec,symrel,qtimrev)
! write (*,*) 'read_gkk : tmpflg after ', tmpflg
! ENDDEBUG

! ===============================================
! Full test: need all perturbations explicitly
! ===============================================

  test_flag = 0
  if (sum(tmpflg(:,1:natom,:,1:natom)) == 3*natom*3*natom .and. tdonecompl == 0) test_flag = 1
! write(*,*) 'read_gkk : test_flag = ', test_flag

! DEBUG
! other possibilities: try completion with symops...
! 
! test_flag = 1
! do idir1=1,3
! do iatom1=1,natom
! ipert1 = (iatom1-1)*3+idir1
! if (gkk_flag(ipert1,ipert1,1,1,iqptfull) < 0) test_flag = 0
! end do
! end do
! write (*,*) 'read_gkk : test_flag = ', test_flag
! 
! test_flag = 1
! do sympertcase1=1,elph_ds%nbranch
! do sympertcase2=1,elph_ds%nbranch
! !  if the perturbation is needed (irredpert < 0)
! !     then we should have it (gkk_flag >= 0)
! if (irredpert(1,sympertcase1,sympertcase2,iqptfull) < 0 .and. &
! &          (gkk_flag(sympertcase1,sympertcase1,1,1,iqptfull) < 0 .or. &
! &           gkk_flag(sympertcase2,sympertcase2,1,1,iqptfull) < 0)) then
! test_flag = 0
! exit
! end if
! end do
! end do
! ENDDEBUG

  write(*,*)'read_gkk: tdonecompl = ', tdonecompl

! de-activate completion of perts by symmetry for now.
! Must be called when all irreducible perturbations are in memory!!!!
  if (test_flag == 1 .and. tdonecompl == 0) then

!  write (*,*) ' read_gkk : enter fxgkkphase before completeperts'
!  call fxgkkphase(elph_ds,gkk_flag,h1_mat_el,iqptfull)

!  ========================================================================
!  Now use more general symops to complete the other equivalent
!  perturbations: the kpoints are also shuffled by these symops
!  afterwards h1_mat_el_sq contains gamma_\tau\alpha,\tau'\alpha'
!  in reduced coordinates
!  
!  \gamma_{\tau'\alpha',\tau\alpha} =
!  <psi_{k+q,ib2} | H(1)_{\tau'\alpha'} | psi_{k,ib1} >*
!  \cdot  <psi_{k+q,ib2} | H(1)_{\tau\alpha}   | psi_{k,ib1} >
!  
!  ========================================================================

   call completeperts(elph_ds,FSfulltofull,gkk_flag,h1_mat_el,h1_mat_el_sq,hdr1,indsym,&
&   iqptfull,irredpert,natom,nsym,spqpt,symq,symrec,symrel,qtimrev,tnons)
   tdonecompl = 1

  end if

! ==============================================================
! if we have all the perturbations for this qpoint, proceed
! with scalar product, norm squared, and add weight factors
! 
! SHOULD HAVE A TEST SO h1_mat_el IS NOT OVERWRITTEN
! BEFORE PREVIOUS QPOINT IS FINISHED!!!!!
! REMARK: there is just a factor of |rprimd(:,i)| between
! abinit and decaft matrix elements.
! ==============================================================

  test_flag = 1
  do isppol=1,elph_ds%nsppol
   do iFSkpt=1,elph_ds%nFSkpt
    do ibranch=1,elph_ds%nbranch
!    do jbranch=1,elph_ds%nbranch
     if (gkk_flag (ibranch,ibranch,iFSkpt,isppol,iqptfull) == -1) then
      test_flag = 0
      exit
     end if
!    end do
    end do
   end do
  end do

  if (test_flag /= 0) then

   write(message,'(a)')' read_gkk : enter normsq_gkq'
   call wrtout(06,message,'COLL')

!  DEBUG
!  write(100,*) 'matrix elements for band,band,kpt,qpt = ',iqptfull, ' both sppol'
!  do iFSkpt=1,elph_ds%nFSkpt
!  do ib1=1,elph_ds%nFSband
!  do ib2=1,elph_ds%nFSband
!  ibb=(ib1-1)*elph_ds%nFSband+ib2
!  write (100,'(4i5,3e10.2,6e16.6)')ib2,ib1,iFSkpt,&
!  &  iqptfull,FSkpt(:,iFSkpt),&
!  & h1_mat_el(:,ibb,:,iFSkpt,:)
!  end do
!  end do
!  end do
!  ENDDEBUG

!  ! NOTE: should have trivial effect if called _again_ here
!  !    checked that this does not change linewidths if all
!  !    perturbations are used explicitely (and all turned
!  !    by the same phase factors)
!  write (*,*) ' read_gkk : enter fxgkkphase before normsq'
!  call fxgkkphase(elph_ds,gkk_flag,h1_mat_el,iqptfull)

!  ! Optional symmetrization of the gkk before normsq. Supposes all
!  !    h1_mat_el are filled
!  !   BAD IDEA: dephasing between equivalent kpoints drives the
!  !       symmetrization to very small values
!  !   GOOD IDEA: once the phase is fixed with fxgkkphase the results
!  !       of symgamma are much more symmetrical, but the linewidths are
!  !       still much too small... do not know why.
!  call symgamma(elph_ds,FSfulltofull,h1_mat_el,&
!  &   indsym,iqptfull,natom,nsym,symq,symrec)


!  DEBUG
!  call completeperts(elph_ds,FSfulltofull,gkk_flag,h1_mat_el,hdr1,&
!  &      indsym,iqptfull,isppol,irredpert,natom,nsym,spqpt,symq,symrec,symrel,qtimrev,tnons)
!  ENDDEBUG

!  MG temporary array to save ph-linewidths before Fourier interpolation
   allocate (qdata(elph_ds%nbranch,elph_ds%nsppol,3))
   qdata(:,:,:)=zero

   call normsq_gkq(displ_red,eigvec,elph_ds,FSfullpqtofull,FSintweight,FSkpt,&
&   h1_mat_el_sq,iqptfull,phfrq_tmp,spqpt,wf,qdata)

   qdata_tmp(nqptirred,:,:,:)=qdata(:,:,:)
   deallocate (qdata)
!  h1_mat_el(:,:,:,:,:) = zero
  end if

  call hdr_clean(hdr1)

 end do !of i1wf 

!got all the gkk perturbations

 deallocate(eigen1)
 deallocate(h1_mat_el,h1_mat_el_sq)

!normally at this point we have the gkk for all kpoints on the FS
!for all the perturbations. Otherwise a 1WF file is missing.
!NOTE: still havent checked the qpoint grid completeness
 do iqpt=1,nqptirred
  iqptfull = qptirredtofull(iqpt)
  do isppol=1,elph_ds%nsppol
   do iFSkpt=1,elph_ds%nFSkpt
    do ipert=1,elph_ds%nbranch
     if (gkk_flag(ipert,ipert,iFSkpt,isppol,iqptfull) == -1) then
      write (message,'(3a,i5,1x,i5,1x,i5,1x,i5,a,a)')' read_gkk  : ERROR-',ch10,           &
&      ' gkk element',ipert,iFSkpt,isppol,iqptfull,' was not found by symmetry operations ',&
&      ' on the irreducible perturbations and qpoints given'
!     write (*,*) '  content = ', elph_ds%gkk_qpt(:,ipert,1,1,iFSkpt,iqptfull)
      call wrtout(06,message,'COLL')
      call leave_new('COLL')
     end if
    end do !ipert
   end do !iFSkpt
  end do !isppol
 end do !iqpt

 write(message,'(a)')'read_gkk : done completing the perturbations (and checked!)'
 call wrtout(06,message,'COLL')

!MG save phonon frequencies, ph-linewidths and lambda(q,n) values before Fourier interpolation
 allocate(elph_ds%qgrid_data(nqptirred,elph_ds%nbranch,elph_ds%nsppol,3))
 allocate(elph_ds%qirredtofull(nqptirred))

 do ii=1,nqptirred
  elph_ds%qirredtofull(ii)=qptirredtofull(ii)
  elph_ds%qgrid_data(ii,:,:,:)=qdata_tmp(ii,:,:,:)
 end do

 deallocate(qdata_tmp)

end subroutine read_gkk
!!***
