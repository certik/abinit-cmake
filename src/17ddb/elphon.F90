!{\src2tex{textfont=tt}}
!!****f* ABINIT/elphon
!!
!! NAME
!! elphon
!!
!! FUNCTION
!! This routine extracts the electron phonon coupling matrix
!! elements and calculates related properties - Tc, phonon linewidths...
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer, MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   a2fsmear = smearing for alpha2F function
!!   acell_in(3)= input length scales of cell (bohr)
!!   amu(ntypat)=mass of the atoms (atomic mass unit)
!!   atmfrc  = inter-atomic force constants from anaddb
!!   blkflg(3*mpert*3*mpert,nblok)= flag for existence of blok
!!   blkqpt(1<ii<9,nblok)=q vector of a phonon mode (ii=1,2,3)
!!   blknrm(1,nblok)=norm of qpt providing normalization
!!   blktyp =
!!   blkval(2,3*mpert*3*mpert,nblok)= all the dynamical matrices
!!   brav = type of Bravais lattice
!!   dielt(3,3) = dielectric tensor
!!   dipdip  =dipole dipole interaction flag
!!   dyewq0(3,3,natom)=atomic self-interaction correction to the
!!        dynamical matrix (only when dipdip=1)
!!   elphsmear = smearing width for gaussian integration
!!     or buffer in energy for calculations with tetrahedra (telphint=0)
!!   elph_fermie = input value of Fermi energy
!!     0 means use value from wfk file
!!   elph_base_name=base name for output files
!!   enunit = governs the units to be used for the output of the phonon frequencies and e-ph quantities
!!   gkk2exist= flag to read in gkk2 matrix elements
!!   gkk2write= flag to write out gkk2 matrix elements to disk
!!   gkk_rptexist= flag to read in real space gkk_rpt matrix elements
!!   gkk_rptwrite= flag to write out real space gkk_rpt matrix elements to disk
!!   gkqexist= flag to read in gkq matrix elements
!!   gkqwrite= flag to write out gkq matrix elements to disk
!!   phfrqexist= flag to read in phonon frequencies
!!   phfrqwrite= flag to write out phonon frequencies to disk
!!   gmet(3,3) =metric in reciprocal space (telphint=1)
!!   gprim(3,3) =dimensionless basis vectors of reciprocal space
!!   ifcflag = flag for IFC matrices in anaddb calling routine
!!     the IFCs are presumed to be known!
!!   ifltransport= flag for transport properties (no=0: yes=1 )
!!   indsym = mapping of atoms btw themselves under symmetry
!!   mpert =maximum number of ipert
!!   msym =maximum number of symmetries
!!   mustar = parameter for Coulombic pseudo-potential in McMillan T_c calculation
!!   natom=number of atoms in cell
!!   nblok= number of bloks in the DDB
!!   ngqpt(3)=integers defining the number of points in the qpt sampling
!!   nqpath=number of vertices in the path in reciprocal space, for band structure
!!      and phonon linewidth output
!!   nqshft= number of shift vectors for defining the sampling of q points
!!   nrpt =number of real space points used to integrate IFC (for
!!        interpolation of dynamical matrices)
!!   nsym=number of space group symmetries
!!   ntypat = number of types of atoms
!!   prtfsurf = integer flag for the output of the Fermi surface (XCrysden file format)
!!   prtnest = integer flag for the calculation of the nesting function
!!   qpath=vertices in the path in reciprocal space, for band structure
!!      and phonon linewidth output
!!   q1shft(3,4) =qpoint shifts considered
!!   rcan(3,natom) =canonical positions of atoms
!!   rmet(3,3)=metric tensor in real space (bohr^2)
!!   rprim_in(3,3)= input primitive translation vectors
!!   rpt(3,nprt) =canonical positions of R points in the unit cell
!!   symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!!   symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!!   telphint = flag for integration over the FS with 0=tetrahedra 1=gaussians
!!   tkeepbands = flag to keep gamma matrix dependence on electronic bands
!!   tnons(3,nsym)=fractional nonsymmorphic translations
!!   trans(3,natom) = Atomic translations : xred = rcan + trans
!!   typat(natom)=type integer for each atom in cell
!!   ucvol=unit cell volume in bohr**3
!!   unitgkk = fortran unit for file containing matrix elements, from mrggkk
!!   wghatm(natom,natom,nrpt) =Weight for the pair of atoms and the R vector
!!   xred(3,natom)=fractional dimensionless atomic coordinates
!!   zeff(3,3,natom) =effective charge on each atom, versus electric
!!        field and atomic displacement
!!
!! OUTPUT
!!
!! NOTES
!!  inspired to a large extent by epcouple.f from the DecAFT package by J. Kay Dewhurst
!!  most inputs taken from mkifc.f
!!  in anaddb ifcflag must be 1 such that the IFC are calculated in atmfrc prior to calling elphon
!!
!!  brav not taken into account propely in all of the code. (MG?)
!!
!!  could choose to make a full 3 dimensional kpt array (:,:,:). Easier for many operations
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      canon9,eliashberg_1d,get_all_gkk2,get_all_gkq,get_all_gkr,get_fs_kpts
!!      get_gkk_qpt_tr,get_tetra,get_tetra_weight,getkgrid,hdr_clean
!!      integrate_gamma,integrate_gamma_tr,leave_new,matr3inv,mka2f,mka2f_tr
!!      mka2fqgrid,mkfskgrid,mknesting,mkph_linwid,mkqptequiv,mkrdim
!!      order_fs_kpts,outelph,printbxsf,rchkgsheader,smpbz,symkpt,wrtout
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine elphon(a2fsmear,acell_in,amu,atmfrc,blkflg,blkqpt,blknrm,blktyp,blkval,&
& brav,ddkfilename,dielt,dipdip,dyewq0,elphsmear,elph_fermie,elph_base_name,enunit,ep_b_min,ep_b_max,&
& gkk2exist,gkk2write,gkk_rptexist,gkk_rptwrite,gkqexist,gkqwrite,&
& phfrqexist,phfrqwrite,prtfsurf,prtnest,gmet,gprim,&
& ifcflag,ifltransport,indsym,kptrlatt,mpert,mpi_enreg,msym,&
& mustar,natom,nblok,ngqpt,nqpath,nqshft,nrpt,nsym,ntypat,qpath,&
& q1shft,rcan,rmet,rprim_in,rpt,symrec,symrel,telphint,tkeepbands,&
& doscalprod,tnons,trans,typat,ucvol,unitgkk,wghatm,xred,zeff)

 use defs_basis
  use defs_datatypes
  use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13recipspace
 use interfaces_14iowfdenpot
 use interfaces_14occeig
 use interfaces_17ddb, except_this_one => elphon
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: brav,dipdip,doscalprod,enunit,ep_b_max,ep_b_min
 integer,intent(in) :: gkk2exist,gkk2write,gkk_rptexist,gkk_rptwrite,gkqexist
 integer,intent(in) :: gkqwrite,ifcflag,ifltransport,mpert,msym,natom,nblok
 integer,intent(in) :: nqpath,nrpt,nsym,ntypat,phfrqexist,phfrqwrite,prtfsurf
 integer,intent(in) :: prtnest,tkeepbands,unitgkk
 integer,intent(inout) :: nqshft
 real(dp),intent(in) :: a2fsmear,elph_fermie,elphsmear,mustar,ucvol
 character(len=fnlen) :: elph_base_name
 character(len=fnlen),intent(in) :: ddkfilename
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: blkflg(3*mpert*3*mpert,nblok),blktyp(nblok)
 integer,intent(in) :: indsym(4,nsym,natom),kptrlatt(3,3),ngqpt(3)
 integer,intent(in) :: symrec(3,3,nsym),symrel(3,3,nsym),typat(natom)
 real(dp),intent(in) :: acell_in(3),amu(ntypat),atmfrc(2,3,natom,3,natom,nrpt)
 real(dp),intent(in) :: blknrm(1,nblok),blkqpt(9,nblok)
 real(dp),intent(in) :: blkval(2,3*mpert*3*mpert,nblok),dielt(3,3)
 real(dp),intent(in) :: dyewq0(3,3,natom),gmet(3,3),gprim(3,3),q1shft(3,4)
 real(dp),intent(in) :: qpath(3,nqpath),rcan(3,natom),rmet(3,3),rprim_in(3,3)
 real(dp),intent(in) :: rpt(3,nrpt),tnons(3,nsym),trans(3,natom)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt),xred(3,natom),zeff(3,3,natom)

!Local variables-------------------------------
! integer,allocatable :: FStoqpt(:,:)
!scalars
 integer :: berryopt,dispindx,eivec,fform,iFSkpt,iFSkpt0,iFSqpt,iatom,ib1,ib2
 integer :: iband,ibranch,ieliash,ifermi,ii,ikpt,ikpt1,ikpt2,iost,iout,ipert1
 integer :: ipert2,iqpt,irpt,iscf,isppol,istat,isym,itim,jFSkpt,jFSqpt,jatom
 integer :: kFSkpt,maxFSband,minFSband,mqpt,mtetra,n1wf,nFSkpt,nFSqpt,nband
 integer :: neliash,nene,new,nfullkpt,ngkkband,nkpt,nomega,nptsym,nqpt1
 integer :: nqpt_computed,nqptirred,nsegment,ntetra,onegkk2size,onegkksize
 integer :: option,qtor,rdwr,telphint,timrev,unitfile
 real(dp) :: Tc,deltaene,elphsmear_tmp,enemax,enemin,factor,gaussfactor
 real(dp) :: gaussprefactor,intwei,max_occ,qphnrm,qptrlen,rcvol,realdp_ex
 real(dp) :: refdist,res,ss,sumi,sumr,tmpdist,vv,xx
 logical :: make_gkk2
 character(len=500) :: message
 character(len=fnlen) :: fname
 type(elph_tr_type) :: elph_tr_ds
 type(elph_type) :: elph_ds
 type(hdr_type) :: hdr
 type(phon_type) :: phon_ds
!arrays
 integer :: bravais(11),dsifkpt(3),ptsymrel(3,3,nsym),qptrlatt(3,3)
 integer :: rptrlatt(3,3),symafm(nsym),vacuum(3)
 integer,allocatable :: FSfullpqtofull(:,:),FSfulltofull(:,:,:)
 integer,allocatable :: FSfulltoirred(:,:),FSirredtoGS(:),FSirredtofull(:)
 integer,allocatable :: FSkptflag(:),FSkptirrank(:),gkk_flag(:,:,:,:,:)
 integer,allocatable :: indqpt1(:),npoint(:),qpttoqpt(:,:,:),tetra_full(:,:,:)
 integer,allocatable :: tetra_mult(:),tetra_wrap(:,:,:),tmpFSfulltofull(:,:,:)
 integer,allocatable :: tmpFSfulltoirred(:,:)
 real(dp) :: acell(3),displ(2*3*natom*3*natom),eigval(3*natom)
 real(dp) :: eigvec(3*3*natom*3*natom),ftwghtgkk(natom,nrpt),gprimd(3,3)
 real(dp) :: klatt(3,3),kpt(3),phfrq_tmp(3*natom),qpt(3),rlatt(3,3),rprim(3,3)
 real(dp) :: rprimd(3,3),shiftk(3),tmpkpt(3)
 real(dp),allocatable :: FSirredwtk(:),FSkpt(:,:),FSkptirred(:,:),FSqpt(:,:)
 real(dp),allocatable :: a2f_1d(:),d2cart(:,:,:,:,:),delta(:,:),dos_phon(:)
 real(dp),allocatable :: dtweightde(:,:),dummy_eig(:,:,:),dynmat(:,:,:,:)
 real(dp),allocatable :: eigenGS(:,:,:),gkk2_tmp(:,:,:,:),gkk_tmp(:,:,:,:,:,:)
 real(dp),allocatable :: gkk_tmp_full(:,:,:,:,:),qptirred(:,:),spqpt(:,:)
 real(dp),allocatable :: sprpt(:,:),tmpFSkpt(:,:),tmpFSqpt(:,:),tmp_eigen(:)
 real(dp),allocatable :: tmpshifts(:,:),tweight(:,:),wtq(:),wtq_folded(:)
 real(dp),allocatable :: zz(:,:)

! *************************************************************************

!DEBUG
!write (message,'(a)')' elphon : enter'
!call wrtout(06,message,'COLL')
!ENDDEBUG

!----------------------------------------------------------
!initializations
!----------------------------------------------------------

!number of phonon modes = 3 * natom
 elph_ds%nbranch = 3*natom

 elph_ds%tkeepbands = tkeepbands
 elph_ds%mustar=mustar

!maximum number of Matsubara frequencies.
!The precise number used depends on the value of Tc:
!they span $w_n = (2n+1) \pi T_c$  where $abs(w_n) < w_{cutoff}$
!ie $|n| < n_{cutoff} = ( \frac{w_{cutoff}}{\pi T_c} ) / 2$
 elph_ds%na2f = 400
 write (*,*) ' elphon : na2f = ', elph_ds%na2f
 elph_ds%a2fsmear = a2fsmear

!maximum number of iterations to converge Tc
 neliash = 10

!use time reversal symmetry always when possible for kpoint reduction,
!and suppose it has been used in WF generation
!not used for the moment: values are always taken from input files.
 timrev = 1


!==================================
!Initialization of some variables
!==================================

 neliash = 10                     !maximum number of iterations to converge Tc
 elph_ds%mustar=mustar            !input mustar
 elph_ds%nbranch = 3*natom        !number of phonon modes = 3 * natom
 elph_ds%tkeepbands = tkeepbands  !flag to sum over bands
 elph_ds%a2fsmear = a2fsmear      !smearing for Eliashber functions

 elph_ds%na2f = 400               !maximum number of Matsubara frequencies.
!The precise number used depends on the value of Tc:
!they span $w_n = (2n+1) \pi T_c$  where $abs(w_n) < w_{cutoff}$
!ie $|n| < n_{cutoff} = ( \frac{w_{cutoff}}{\pi T_c} ) / 2$
 
 timrev = 1                       !use time reversal symmetry always when possible for kpoint reduction,
!and suppose it has been used in WF generation
!not used for the moment: values are always taken from input files.

!save gkk data for full kpoints to file on disk
!or read them from file if it exists
!gkk2write == 0 gkk2exist == 0 makes the gkk2() be
!allocated and used in memory to store the data (faster)

 elph_ds%gkqwrite=gkqwrite         ; elph_ds%gkqexist=gkqexist
 elph_ds%gkk_rptwrite=gkk_rptwrite ; elph_ds%gkk_rptexist=gkk_rptexist
 elph_ds%gkk2write=gkk2write       ; elph_ds%gkk2exist=gkk2exist
 elph_ds%phfrqwrite=phfrqwrite     ; elph_ds%phfrqexist=phfrqexist

!This should never be turned off: symmetrization of elphon matrix elements
!in complete_gkk. See get_all_gkq
 elph_ds%tsymgkq=1

 elph_ds%elph_base_name = trim(elph_base_name)

!normalize input rprim and acell.
 do ii=1,3
  ss = sqrt(rprim_in(1,ii)**2+rprim_in(2,ii)**2+rprim_in(3,ii)**2)
  rprim(:,ii) = rprim_in(:,ii)/ss
  acell(ii) = acell_in(ii) * ss
 end do

!make dimension-ful rprimd and gprimd for transformation of derivatives to cartesian coordinates.
 call mkrdim(acell,rprim,rprimd)
 call matr3inv(rprimd,gprimd)

!===================
!Check some inputs
!===================

 if (nsym==1) then
  write (message,'(7a)')ch10,&
&  ' elphon: COMMENT- ',ch10,&
&  ' Symmetries are not used! ',ch10,&
&  ' Full matrix elements must be supplied for all perturbations and qpoints!',ch10
  call wrtout(06,message,'COLL')
  call wrtout(ab_out,message,'COLL')
  
  if (abs(tnons(1,1))+abs(tnons(2,1))+abs(tnons(3,1)) > tol10) then
   write (message,'(4a)')ch10,&
&   ' elphon : ERROR-',ch10,&
&   ' tnons should be (0,0,0) for unique symmetry Id'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if
 end if

 if (ifcflag/=1) then
  write(message,'(6a,i4)')ch10,&
&  ' elphon : ERROR-',ch10,&
&  ' ifcflag should be set to one,',ch10,&
&  ' the IFC matrices are supposed to exist.',ifcflag
  call wrtout(06,message,'COLL')
  call wrtout(ab_out,message,'COLL')
 end if

!=================================
!Set up the full grid of qpoints
!=================================

!qpoint lattice vectors (inverse, like kptrlatt)
 qptrlatt(:,:)=0
 qptrlatt(1,1)=ngqpt(1)
 qptrlatt(2,2)=ngqpt(2)
 qptrlatt(3,3)=ngqpt(3)

 if (nqshft /= 1) then
! MG Check this part, dont know if it works, never used shifted grids!
! try to reduce the qpoint grid to a single qshift, otherwise stop
! dummy args for call to getkgrid
! dsifkpt(1) = ngqpt(1)
! dsifkpt(2) = ngqpt(2)
! dsifkpt(3) = ngqpt(3)
  dsifkpt(:) = 1
  wtq(:) = one
  symafm(:) = 1
  vacuum(:) = 0
  iscf = 3
  
  mqpt = ngqpt(1)*ngqpt(2)*ngqpt(3)*nqshft
  allocate(spqpt(3,mqpt),wtq(mqpt),tmpshifts(3,8),stat=istat)
  if (istat /= 0) stop 'elphon: error in allocating spqpt(3,mqpt),wtq(mqpt),tmpshifts(3,8)'

  tmpshifts(:,:) = zero
  tmpshifts(:,1:4) = q1shft(:,:)

  iout=6

! call symbrav to fill bravais
  berryopt = 1
  call symbrav(berryopt,bravais,iout,nsym,nptsym,ptsymrel,rmet,rprimd)

  write (*,*) ' bravais ', bravais

! DEBUG
! ! call testkgrid to fill qptrlatt
! call testkgrid(bravais,iout,qptrlatt,qptrlen,&
! &   nsym,nqshft,nsym,0,rprimd,tmpshifts,symafm,symrel,tnons,vacuum)
! 
! write (*,*) 'qptrlatt before getkgrid = '
! write (*,*) qptrlatt
! ENDDEBUG

! call getkgrid to reduce qptrlatt to 1 shift 
  iscf = 3
! call getkgrid(dsifkpt,iout,iscf,spqpt,1,qptrlatt,qptrlen, &
! & nsym,mqpt,nqpt_computed,nqshft,nsym,rprimd,tmpshifts,symafm, &
! &                symrel,tnons,vacuum,wtq)

! just call with identity, to get full set of kpts in spqpt, but
! reduce qshfts
  
  call getkgrid(dsifkpt,iout,iscf,spqpt,3,qptrlatt,qptrlen, &
&  1,mqpt,nqpt_computed,nqshft,1,rprimd,tmpshifts,symafm, &
&  symrel,tnons,vacuum,wtq)
  deallocate (spqpt,wtq,tmpshifts)
  if (nqshft /= 1) then
   write (message,'(6a,i4)')ch10,&
&   ' elphon : ERROR- ',ch10,&
&   ' multiple qpt shifts not treated yet',ch10,&
&   ' -- should be possible ', nqshft
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if
 end if

 write(message,'(3a,9i3)')&
& ' elphon : enter smpbz ',ch10,&
& ' qptrlatt = ',qptrlatt 
 call wrtout(06,message,'COLL')

 option=1
!mqpt=ngqpt(1)*ngqpt(2)*ngqpt(3)*nqshft
 mqpt= qptrlatt(1,1)*qptrlatt(2,2)*qptrlatt(3,3) &
& +qptrlatt(1,2)*qptrlatt(2,3)*qptrlatt(3,1) &
& +qptrlatt(1,3)*qptrlatt(2,1)*qptrlatt(3,2) &
& -qptrlatt(1,2)*qptrlatt(2,1)*qptrlatt(3,3) &
& -qptrlatt(1,3)*qptrlatt(2,2)*qptrlatt(3,1) &
& -qptrlatt(1,1)*qptrlatt(2,3)*qptrlatt(3,2)

 allocate(spqpt(3,mqpt),sprpt(3,mqpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating spqpt ,sprpt'
 iout = 6
 call smpbz(brav,iout,qptrlatt,mqpt,elph_ds%nqpt,nqshft,option,q1shft,spqpt)

!save the q-grid for future reference
!MG In the future spqpt should be replaced by elph_ds%spqpt
 allocate(elph_ds%spqpt(3,elph_ds%nqpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating elph_ds%spqpt'

!reduce spqpt to correct zone
 do iqpt=1,elph_ds%nqpt
  call canon9(spqpt(1,iqpt),kpt(1),res)
  call canon9(spqpt(2,iqpt),kpt(2),res)
  call canon9(spqpt(3,iqpt),kpt(3),res)
  spqpt(:,iqpt) = kpt
  elph_ds%spqpt(:,iqpt)=kpt
 end do

!=================================================================
!Calculate weights, needed to estimate lambda using the weighted
!sum of the uninterpolated e-ph matrix elements
!=================================================================

 write (message,'(a)')' setqgrid : calling symkpt to find irred q points'
 call wrtout(6,message,'COLL')

 allocate(indqpt1(elph_ds%nqpt),wtq_folded(elph_ds%nqpt),wtq(elph_ds%nqpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating indqpt1 ,wtq_folded ,wtq'
 wtq(:) = one/dble(elph_ds%nqpt) !weights normalized to unity

 option=0 !do not write to ab_out
 call symkpt(gmet,indqpt1,elph_ds%spqpt,elph_ds%nqpt,nqpt1,nsym,option,symrec,&
& timrev,wtq,wtq_folded)

 write (message,'(2a,i6)')ch10,' Number of irreducibile q-points = ',nqpt1
 call wrtout(6,message,'COLL')
!call wrtout(ab_out,message,'COLL')
 elph_ds%nqptirred=nqpt1

 write (message,'(a)')' Irreducible q points with weights :'
 call wrtout(6,message,'COLL')
!call wrtout(ab_out,message,'COLL')

 do iqpt=1,elph_ds%nqpt
  if (wtq_folded(iqpt) /= zero) then
   write (message,'(1x,i4,a2,4es16.8)')iqpt,') ',elph_ds%spqpt(:,iqpt),wtq_folded(iqpt)
   call wrtout(6,message,'COLL')
!  call wrtout(ab_out,message,'COLL')
  end if
 end do

 write(message,'(a)')ch10
 call wrtout(6,message,'COLL')
!call wrtout(ab_out,message,'COLL')

!DEBUG
!write(*,*)'indqpt1',indqpt1
!write(*,*)'wtq_folded',wtq_folded
!write(*,*)'wtq',wtq
!ENDDEBUG

 allocate (elph_ds%wtq(elph_ds%nqpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating elph_ds%wtq'
 elph_ds%wtq(:)=wtq_folded(:)
!MEMO indqpt could be useful to test the qgrid read by abinit
 deallocate (indqpt1,wtq_folded,wtq)
!ENDMG20060531

!================================================================
!Summed weights for gkk FT. Is wghatm symmetric on the 2 atoms?
!present version does not use ftwghtgkk
!================================================================
 do iatom=1,natom
  do irpt=1,nrpt
   ftwghtgkk(iatom,irpt) = wghatm(iatom,iatom,irpt)
  end do
 end do

!DEBUG
!write (*,*)' elphon : wghatm = '
!do iatom=1,natom
!do jatom=1,natom
!write (*,'(4E16.5)') (wghatm(jatom,iatom,irpt),rpt(:,irpt),irpt=1,nrpt)
!end do
!end do
!ENDDEBUG

!====================================
!Read the GS header of the GKK file
!====================================
 
 write (message,'(2a)')ch10,' elphon : reading and checking the GS header of the GKK file'
 call wrtout (06,message,'COLL')

 call rchkGSheader(hdr,natom,nband,unitgkk)

 elph_ds%nsppol=hdr%nsppol

!MG NOTE that fermie has been removed and replaced by elph_ds%fermie
!Set elph_ds%fermie: either comes from anaddb input file or from wfk file
 elph_ds%fermie = hdr%fermie
 if (abs(elph_fermie) > tol10) then
  elph_ds%fermie = elph_fermie
  write(message,'(a,f12.6)')' Fermi level set by the user at :',elph_ds%fermie
  call wrtout(06,message,'COLL')
 end if

 do iatom=1,natom
  do irpt=1,nrpt
   ftwghtgkk(iatom,irpt) = wghatm(iatom,iatom,irpt)
  end do
 end do

!==================================================
!Read GS eigenvalues for each irreducible kpt and
!number of 1WF files contributing to the GKK file
!==================================================

!allocate arrays which depend on general variables
 allocate(eigenGS(nband,hdr%nkpt,elph_ds%nsppol),FSkptflag(hdr%nkpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating eigenGS FSkptflag'
!SPPOL
 do isppol=1,elph_ds%nsppol
  do ikpt=1,hdr%nkpt
   read(unitgkk) eigenGS(:,ikpt,isppol)
!  DEBUG
!  write(*,'(i4,2x,,10(e10.4,1x))')ikpt,eigenGS(:,ikpt,isppol)
!  ENDDEBUG
  end do
 end do

!read number of 1WF files contributing to the GKK file
 read(unitgkk) n1wf
 write(message,'(a,i4)')' elphon : number of perturbations in the gkk file = ',n1wf
 call wrtout(06,message,'COLL')

!====================================================================
!Setup of the kgrid :
!1) get bands near Ef
!2) order irred kpoints 
!3) reconstruct full kgrid from irred kpoints,
!4) make tables from IBZ to FULLZONE
!4) setup weights for integration (gaussian or tetrahedron method) 
!5) calculate DOS at Ef
!====================================================================

 elphsmear_tmp = elphsmear

 call get_fs_kpts(eigenGS,elph_ds,FSkptflag,elphsmear_tmp,hdr)

!DONETOHERE

 if (elph_ds%tkeepbands == 0) then !we are summing over bands
  elph_ds%ngkkband = 1
 else if (elph_ds%tkeepbands == 1) then
! use trivial integration weights since average over bands is done in normsq_gkk
! and keep the band dependency btw elph_ds%minFS and elph_ds%maxFS) 
  elph_ds%ngkkband = elph_ds%nFSband
 else
  write(message,'(4a,i4)')ch10,&
&  ' elphon : BUG- ',ch10,     &
&  ' tkeepbands must be 0 or 1 while it is : ',elph_ds%tkeepbands
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 minFSband=elph_ds%minFSband
 maxFSband=elph_ds%maxFSband

 write(message,'(a,i4,2x,i4)')' elphon : minFSband, maxFSband = ',minFSband,maxFSband
 call wrtout(06,message,'COLL')

!Check consistency of arguments, and that FSintweight will be well dimensioned
 if (elph_ds%tkeepbands == 1 .and. elph_ds%nFSband /= elph_ds%ngkkband) then
  write (*,*) 'elphon : Error: nFSband /= ngkkband'
  stop
 end if

 allocate(FSkptirred(3,elph_ds%nFSkptirred),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating FSkptirred'
 allocate(FSirredtoGS(elph_ds%nFSkptirred),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating FSirredtoGS'
 allocate(FSkptirrank(elph_ds%nFSkptirred),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating FSkptirrank'

 call order_fs_kpts(FSkptflag,FSkptirrank,FSkptirred,FSirredtoGS,&
& hdr,elph_ds%nFSkptirred)

!==========================================
!Set up full FSkpt grid from irred points
!==========================================

 allocate(tmpFSkpt(3,2*elph_ds%nFSkptirred*nsym),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating tmpFSkpt'
 allocate(tmpFSfulltoirred(3,2*elph_ds%nFSkptirred*nsym),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating tmpFSfulltoirred'
 allocate(tmpFSfulltofull(2,nsym,2*elph_ds%nFSkptirred*nsym),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating tmpFSfulltofull'
 allocate(FSirredwtk(elph_ds%nFSkptirred),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating FSirredwtk'
 allocate(FSirredtofull(elph_ds%nFSkptirred),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating FSirredtofull'

 call mkFSkgrid(FSirredwtk,FSirredtofull,tmpFSfulltoirred,&
& tmpFSfulltofull,tmpFSkpt,FSkptirred,&
& elph_ds%nFSkptirred,nFSkpt,nsym,symrec,timrev)

 elph_ds%nFSkpt = nFSkpt

!allocate FS arrays for kpoints only to size needed
 allocate(elph_ds%FSkpt(3,elph_ds%nFSkpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating elph_ds%FSkpt'
 allocate(FSkpt(3,elph_ds%nFSkpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating FSkpt'
 allocate(FSfulltofull(2,nsym,elph_ds%nFSkpt),FSfulltoirred(3,elph_ds%nFSkpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating FSfulltofull,FSfulltoirred'

!copy the kpoints and full --> irred FSkpt map
 do ikpt=1,elph_ds%nFSkpt
! MG FIXME 
! I dont like this part since we are duplicating FSkpt but I need elph_ds%FSkpt in my subroutines
! In the next version the setup of the kgrid could be done in the new subroutines setkgrid.F90
! and all these quantities should be saved in elph_ds%
  elph_ds%FSkpt(:,ikpt) = tmpFSkpt(:,ikpt)
! ENDMG
  FSkpt(:,ikpt) = tmpFSkpt(:,ikpt)
  FSfulltoirred(:,ikpt) = tmpFSfulltoirred(:,ikpt)
  do isym=1,nsym
   FSfulltofull(:,isym,ikpt) = tmpFSfulltofull(:,isym,ikpt)
  end do
 end do
 
 deallocate(tmpFSkpt,tmpFSfulltoirred,tmpFSfulltofull)

!in spinor or spin polarized case, orbitals have occupation <= 1 instead of 2
!max_occ = two
 max_occ = one
 if (hdr%nspinor > 1 .or. elph_ds%nsppol > 1) max_occ = half

!
!NOTE: should start encapsulation of intweight calculation HERE
!
!===================================
!Set up integration weights for FS
!telphint = 0 or 1
!1 = gaussians
!0 = tetrahedron method
!===================================

!SPPOL allocate for nsppol 
 allocate(elph_ds%FSintweight(elph_ds%nFSband,elph_ds%nFSkpt,elph_ds%nsppol),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating elph_ds%FSintweight'
 allocate(elph_ds%n0(elph_ds%nsppol),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating elph_ds%n0'

 if (telphint == 1) then

! ==============================================================
! Gaussian integration:
! Each FSkpt contributes a gaussian of integrated weight 1 
! for each band. The gaussian being centered at the Fermi level.
! ===============================================================

! gaussprefactor = sqrt(piinv)/elphsmear/elph_ds%nFSkpt
! took out factor 1/nkpt which intervenes only at integration time

! MJV 18/5/2008 does gaussprefactor need to contain max_occ?
  gaussprefactor = sqrt(piinv)/elphsmear
  gaussfactor = one/elphsmear
  elph_ds%FSintweight = zero
! SPPOL loop on isppol as well to get 2 sets of weights
  do isppol=1,elph_ds%nsppol
   do iFSkpt=1,elph_ds%nFSkpt
    ikpt = FSirredtoGS(FSfulltoirred(1,iFSkpt))
    do ib1=1,elph_ds%nFSband
     xx = gaussfactor*(eigenGS(minFSband-1+ib1,ikpt,isppol) - elph_ds%fermie)
!    DEBUG
!    write (*,*) 'xx = ',xx
!    ENDDEBUG
     elph_ds%FSintweight(ib1,iFSkpt,isppol) = exp(-xx*xx)*gaussprefactor
    end do
   end do
  end do
  

 else if (telphint == 0) then

! =========================
! Tetrahedron integration
! =========================

! Get tetrahedra, ie indexes of the full kpoints at their summits
  mtetra = 6 * elph_ds%nFSkpt
  ntetra = mtetra
  
! tetra_full(:,1,i) contains the irred kpt  number
! tetra_full(:,2,i) contains the fullbz_kpt number
! tetra_wrap(:,:,i) contains the a flag to wrap kpoints outside the IBZ (+-1)
! to get the irreducible tetrahedra, the number of equivalent
! tetrahedra is counted in tetra_mult and the inequivalent few
! (ntetra < mtetra) are packed into the beginning of tetra_full
! klatt here refer to the shortest kpt vectors, not qpt as below
! convert kptrlatt to double and invert

  rlatt(:,:) = kptrlatt(:,:)
  call matr3inv(rlatt,klatt)

  allocate (tetra_full(4,2,mtetra),tetra_wrap(3,4,mtetra),stat=istat)
  if (istat /= 0) stop 'elphon: error in allocating tetra_full,tetra_wrap'
  allocate(tetra_mult(mtetra),stat=istat)
  if (istat /= 0) stop 'elphon: error in allocating tetra_mult'

! call get_tetra (FSfulltoirred(1,:),gprimd,klatt,FSkpt,mtetra,&
! &       elph_ds%nFSkpt,ntetra,&
! &       tetra_full,tetra_mult,tetra_wrap,vv)
! 
! call get_tetra without using symmetries: trivial indkpt
! 
  call get_tetra(FSfulltofull(1,1,:),gprimd,klatt,FSkpt,mtetra,&
&  elph_ds%nFSkpt,ntetra,tetra_full,tetra_mult,tetra_wrap,vv)
! DEBUG
! write(*,*)' elphon : vv = ', vv
! write(*,*)' elphon : ntetra,mtetra = ', ntetra,mtetra
! write(*,*)' elphon : tetra_full = ', tetra_full
! write(*,*)' elphon : tetra_mult = ', tetra_mult
! ENDDEBUG

  rcvol = abs (gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
&  -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
&  +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))
! DEBUG
! write (*,*) ' elphon : rcvol = ', rcvol
! ENDDEBUG

! just do weights at FS
  nene = 100
! enemin = elph_ds%fermie
! enemax = elph_ds%fermie + elphsmear
  enemin = minval(eigenGS) - elphsmear
  enemax = maxval(eigenGS) + elphsmear
  deltaene = (enemax-enemin)/dble(nene-1)
  ifermi = int((elph_ds%fermie-enemin)/deltaene)+1
! redefine enemin enemax to be at rounded multiples of deltaene
  enemin = elph_ds%fermie - dble(ifermi)*deltaene
  enemax = elph_ds%fermie + dble(nene-ifermi-1)*deltaene

! trial : fix small window around elph_ds%fermie for tetrahedron weight calculation
  deltaene = 2*elphsmear/dble(nene-1)
  ifermi = int(nene/2)
  enemin = elph_ds%fermie - dble(ifermi-1)*deltaene
  enemax = enemin + dble(nene-1)*deltaene

! DEBUG
! write (*,*) 'elphon : deltaenes = ', deltaene,(enemax-enemin)/dble(nene-1)
! write (*,*) 'elphon : ifermi,nene,elph_ds%fermie,diff'
! write (*,*)  ifermi,nene,elph_ds%fermie,&
! &    elph_ds%fermie-(enemin+dble(ifermi)*deltaene)
! ENDDEBUG

  allocate(tmp_eigen(elph_ds%nFSkpt),stat=istat)
  if (istat /= 0) stop 'elphon: error in allocating elph_ds%nFSkpt'
  allocate(tweight(elph_ds%nFSkpt,nene),dtweightde(elph_ds%nFSkpt,nene),stat=istat)
  if (istat /= 0) stop 'elphon: error in allocating tweight,dtweightde'

  do iband = 1,elph_ds%nFSband
!  DEBUG
!  write (*,*) ' elphon : isppol, iband = ', isppol, iband
!  write (*,*) ' elphon : tetra_full(:,:,2) ', tetra_full(:,:,2)
!  ENDDEBUG

!  for each spin pol
   do isppol=1,elph_ds%nsppol
!   For this band get its contribution
    tmp_eigen(:) = zero
    do iFSkpt=1,elph_ds%nFSkpt
     ikpt = FSirredtoGS(FSfulltoirred(1,iFSkpt))
     tmp_eigen(iFSkpt) = eigenGS(minFSband+iband-1,ikpt,isppol)
    end do
!   DEBUG
!   write(113,*)'FS for band ', minFSband+iband-1, 'and sppol ',isppol
!   write(113,'(6E16.6)')  tmp_eigen
!   ENDDEBUG
!   calculate general integration weights at each irred kpoint as in Blochl et al PRB 49 16223
    call get_tetra_weight(tmp_eigen,enemin,enemax,max_occ,mtetra,nene,elph_ds%nFSkpt,&
&    ntetra,rcvol,tetra_full,tetra_mult,tweight,dtweightde,vv)
    elph_ds%FSintweight(iband,:,isppol) = dtweightde(:,ifermi)*elph_ds%nFSkpt
   end do

  end do
  deallocate(tmp_eigen,tweight,dtweightde)
  deallocate(tetra_full,tetra_wrap,tetra_mult)


 else if (telphint == 2) then ! range of bands occupied

! SPPOL eventually be able to specify bands for up and down separately
  elph_ds%FSintweight = zero
  do iFSkpt=1,elph_ds%nFSkpt
   ikpt = FSirredtoGS(FSfulltoirred(1,iFSkpt))
   do ib1=ep_b_min,ep_b_max
!   for the moment both spin channels same
    elph_ds%FSintweight(ib1,iFSkpt,:) = one
   end do
  end do
  
  write (*,*) ' elphon : DOS is calculated from states in bands ',ep_b_min,' to ',&
  ep_b_max
  
 else 
  write (message,'(4a)')ch10,&
&  ' elphon : ERROR-',ch10,&
&  'telphint should be 0 1 or 2 '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if !of if telphint

!SPPOL sum over spin channels to get total DOS
!channels decoupled => use separate values for DOS_up(Ef) resp down
 do isppol=1,elph_ds%nsppol
  elph_ds%n0(isppol) = sum(elph_ds%FSintweight(:,:,isppol))/elph_ds%nFSkpt
 end do

 if (elph_ds%nsppol == 1) then
  write (*,*) ' elphon : the estimated DOS(E_Fermi) = ', elph_ds%n0(1), ' states/Ha/spin '
  write (*,*) ' elphon : the total FS weight and # of kpoints = ',sum(elph_ds%FSintweight),elph_ds%nFSkpt
 else if (elph_ds%nsppol == 2) then
  write (*,*) ' elphon : the spin up   DOS(E_Fermi) = ', elph_ds%n0(1), ' states/Ha/spin '
  write (*,*) ' elphon : the spin down DOS(E_Fermi) = ', elph_ds%n0(2), ' states/Ha/spin '
  write (*,*) ' elphon : total DOS(E_Fermi) = ', elph_ds%n0(1)+elph_ds%n0(2), ' states/Ha '
  write (*,*) ' elphon : the spin up   FS weight and # of kpoints = ',&
&  sum(elph_ds%FSintweight(:,:,1)),elph_ds%nFSkpt
  write (*,*) ' elphon : the spin down FS weight and # of kpoints = ',&
&  sum(elph_ds%FSintweight(:,:,2)),elph_ds%nFSkpt
 else
  write (*,*) 'bad value for nsppol ', elph_ds%nsppol
  stop
 end if

 allocate(elph_ds%gkk_intweight(elph_ds%ngkkband,elph_ds%nFSkpt,elph_ds%nsppol),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating elph_ds%gkk_intweight'

 if (elph_ds%tkeepbands == 0) then
! use trivial integration weights  for single band,
! since average over bands is done in normsq_gkk
  elph_ds%gkk_intweight(1,:,:) = one

 else if (elph_ds%tkeepbands == 1) then
! use elph_ds%FSintweight since average over bands is not done in normsq_gkk
  elph_ds%gkk_intweight(:,:,:) = elph_ds%FSintweight(:,:,:)
 else
  write(message,'(4a,i4)')ch10,' elphon : BUG- ',ch10,            &
&  ' tkeepbands must be 0 or 1 while it is : ',elph_ds%tkeepbands
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!DEBUG
!write (109,*)'the integration weights on the FS are = '
!do isppol=1,elph_ds%nsppol
!do ib1=1,elph_ds%nFSband
!write (109,'(6E16.6)') elph_ds%FSintweight(ib1,:,isppol)
!write (109,*)
!end do
!end do
!ENDDEBUG

!
!NOTE: should end encapsulation of intweight calculation HERE
!


!Output of the Fermi Surface
 if (prtfsurf == 1) then
  fname=trim(elph_ds%elph_base_name) // '_BXSF'

! FIXME
! shiftk is not defined neither in the anaddb nor in the hdr data type
! an incorrect FS will be produced in case of a shifted k-grid used during the GS calculation
! check if we are using a unshifthed kgrid, obviously doesnt work in case
! of multiple shifts containg a zero translation but in this case prtbxsf should work
  shiftk=one
  do ii=1,hdr%nkpt
   if (all(hdr%kptns(:,ii) == zero)) shiftk=zero
  end do

! the input argument in printbxsf has shape (nband,nkpt,nspin), to be on the safe side we use:
  allocate (dummy_eig(nband,hdr%nkpt,elph_ds%nsppol),stat=istat)
  if (istat /= 0) stop 'elphon: error in allocating dummy_eig'
  dummy_eig(:,:,:)=eigenGS(:,:,:)

  call printbxsf(dummy_eig,real(0., dp),elph_ds%fermie,gprimd,kptrlatt,nband,hdr%nkpt,hdr%kptns,nsym,symrec,timrev,&
&  elph_ds%nsppol,shiftk,1,fname)

  deallocate (dummy_eig)
 end if !prtfsurf

!MG Do we need this?
!Get vectors of lattice in recip space
 rlatt(:,:) = qptrlatt(:,:)
 call matr3inv (rlatt,klatt)

 write(*,*) ' elphon : klatt, rlatt = '
 write(*,*) klatt
 write(*,*) rlatt
!END MG

!=========================================================    
!Get equivalence between a FSkpt pair and a qpt in spqpt
!only works if the qpt grid is complete (identical to
!the kpt one, with a basic shift of (0,0,0)
!FStoqpt(k1,k2) = q  => k1 = k2+q
!=========================================================    

!allocate(FStoqpt(elph_ds%nFSkpt,elph_ds%nFSkptirred),stat=istat)
!if (istat /= 0) stop 'elphon: error in allocating FStoqpt'
 allocate(FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating FSfullpqtofull'
 allocate(qpttoqpt(2,nsym,elph_ds%nqpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating qpttoqpt'

 write (message,'(2a)')ch10,' elphon : calling mkqptequiv to set up the FS qpoint set'
 call wrtout(06,message,'COLL')

 call mkqptequiv (FSfullpqtofull,FSkpt,FSkptirred,elph_ds%nFSkpt,&
& elph_ds%nFSkptirred,elph_ds%nqpt,nsym,qpttoqpt,spqpt,symrec)


!==========================================
!Set up dataset for phonon interpolations
!==========================================

 call setup_phon_ds(phon_ds,dipdip,mpert,nsym,natom,ntypat,nrpt,&
& ucvol,indsym,symrel,typat,acell,amu,atmfrc,dielt,dyewq0,gprim,gmet,&
& xred,zeff,rcan,rmet,rprim,rprimd,rpt,trans,wghatm)

!transfer ifltransport flag to structure
 elph_tr_ds%ifltransport=ifltransport
!transfer name of files file for ddk
 elph_tr_ds%ddkfilename=ddkfilename

!DEBUG
!mqpt = 1
!elph_ds%nqpt = 1
!spqpt(:,1) = (/0.5,0.5,0.0/)
!ENDDEBUG

!
!order fs qpts as well!
!
!call order_fs_kpts(FSkptflag,FSkptirrank,FSkptirred,FSirredtoGS,&
!& hdr,elph_ds%nFSkptirred)

!MG20060531 save the q-grid for future reference
!MEMO: In the future spqpt should be replaced by elph_ds%spqpt
 allocate (elph_ds%spqpt(3,elph_ds%nqpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating elph_ds%spqpt'
!ENDMG20060531

!reduce spqpt to correct zone
 do iqpt=1,elph_ds%nqpt
  call canon9(spqpt(1,iqpt),kpt(1),res)
  call canon9(spqpt(2,iqpt),kpt(2),res)
  call canon9(spqpt(3,iqpt),kpt(3),res)
  spqpt(:,iqpt) = kpt
  elph_ds%spqpt(:,iqpt)=kpt
 end do

!MG20060531
!==========================================================
!calculate weights, needed to estimate lambda using the weighted
!sum of the uninterpolated e-ph matrix elements
!==========================================================
 write (message,'(a)')' setqgrid : calling symkpt to find irred q points'
 call wrtout(6,message,'COLL')

 allocate (indqpt1(elph_ds%nqpt),wtq_folded(elph_ds%nqpt),wtq(elph_ds%nqpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating indqpt1,wtq_folded,wtq'
 wtq(:) = one/dble(elph_ds%nqpt) !weights normalized to unity

 option=0 !do not write to ab_out
 call symkpt(gmet,indqpt1,elph_ds%spqpt,elph_ds%nqpt,nqpt1,nsym,option,symrec,&
& timrev,wtq,wtq_folded)

 write (message,'(2a,i6)')ch10,' Number of irreducibile q-points = ',nqpt1
 call wrtout(6,message,'COLL')
!call wrtout(ab_out,message,'COLL')

 elph_ds%nqptirred=nqpt1
 write (message,'(a)')' Irreducible q points with weights :'
 call wrtout(6,message,'COLL')
!call wrtout(ab_out,message,'COLL')

 do ii=1,elph_ds%nqpt
  if (wtq_folded(ii) /= zero) then
   write (message,'(1x,i4,a2,4es16.8)')ii,') ',elph_ds%spqpt(:,ii),wtq_folded(ii)
   call wrtout(6,message,'COLL')
!  call wrtout(ab_out,message,'COLL')
  end if
 end do

 write(message,'(a)')ch10
 call wrtout(6,message,'COLL')
!call wrtout(ab_out,message,'COLL')

!DEBUG
!write(*,*)'indqpt1',indqpt1
!write(*,*)'wtq_folded',wtq_folded
!write(*,*)'wtq',wtq
!ENDDEBUG

 allocate (elph_ds%wtq(elph_ds%nqpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating elph_ds%wtq'
 elph_ds%wtq(:)=wtq_folded(:)
!MEMO indqpt could be useful to test the qgrid read by abinit
 deallocate (indqpt1,wtq_folded,wtq)
!ENDMG20060531

!Get vectors of lattice in recip space
!MJV: do we need this one? klatt already calculated and only needed for
!tetrahedron weights.
 rlatt(:,:) = qptrlatt(:,:)
 call matr3inv (rlatt,klatt)

 write(*,*) ' elphon : klatt, rlatt = '
 write(*,*) klatt
 write(*,*) rlatt

!DEBUG
!write (*,*) 'elphon : nbranch, nFSkpt = ',&
!&       elph_ds%nbranch, elph_ds%nFSkpt
!write (*,*) 'elphon : FSkpt(:,1:4) = '
!write (120,*) ' FSkpt(:,:) = '
!write (120,'(3F20.10)') FSkpt(:,:)
!write (*,*) 'elphon : FSirredwtk = 'anaddb_dtset%gkqwrite = 0
!write (*,'(6E16.6)') FSirredwtk
!write (*,*) 'elphon : FSfulltoirred = '
!write (*,'(2(3i5,2x))') FSfulltoirred
!write (*,*) 'elphon : FSfullpqtofull = '
!write (*,'(2(3i7,4x))') (((iFSkpt,iqpt,FSfullpqtofull(iFSkpt,iqpt)),iFSkpt=1,elph_ds%nFSkpt),iqpt=1,elph_ds%nqpt)
!write (*,*) 'elphon : FSirredtofull = '
!write (*,'(2(3i5,2x))') FSirredtofull
!write (122,*) 'elphon : FSfulltofull = '
!do iFSkpt=1,elph_ds%nFSkpt
!write (122,*) 'iFSkpt = ',iFSkpt
!write (122,'(4(2i5,4x))') (FSfulltofull(:,isym,iFSkpt),isym=1,nsym)
!end do
!write (*,*) 'elphon : qpttoqpt = '
!do iqpt=1,min(10,elph_ds%nqpt)
!write (*,'(10(i5,2x))') qpttoqpt(:,:,iqpt)
!end do
!ENDDEBUG

!===================================================
!Allocate all important arrays for FS integrations
!===================================================

!SPPOL allocate with additional dimension for sppol?
 allocate(gkk_flag(elph_ds%nbranch,elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol,elph_ds%nqpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating gkk_flag'
 allocate(elph_ds%phfrq(3*natom,elph_ds%nFSkpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating phfrq'

!Record sizes for matrices on disk: complex and real versions (for real and recip space resp!)
 onegkksize = 2*elph_ds%nbranch*elph_ds%nbranch*elph_ds%ngkkband*elph_ds%ngkkband*&
& elph_ds%nFSkpt*elph_ds%nsppol*kind(realdp_ex)
 elph_tr_ds%onegkksize=onegkksize

!MJV 18/5/2008 : looks like onegkk2size is identical to onegkksize - is it
!used?
 onegkk2size =  elph_ds%nbranch*elph_ds%nbranch*elph_ds%ngkkband*elph_ds%ngkkband*&
& elph_ds%nFSkpt*elph_ds%nsppol*kind(realdp_ex)

 write (message,'(a,3(a,i12,a))')ch10,               &
& ' elphon : kind(dp real)   = ',kind(realdp_ex),ch10,&
& '               onegkksize = ',onegkksize,ch10,     &
& '              onegkk2size = ',onegkk2size,ch10
 call wrtout(06,message,'COLL')

 allocate (qptirred(3,n1wf),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating qptirred'

 write (message,'(4a)')                         &
& ' elphon : preliminary setup completed ',ch10,&
& '          calling get_all_gkq to read in all the e-ph matrix elements',ch10
 call wrtout(06,message,'COLL')

!flag to do scalar product in gkq before interpolation:
!should also used in interpolate_gkk and mkph_linwid
 elph_ds%doscalprod=doscalprod
 if (elph_ds%doscalprod==0) then
  write (*,*) ' elphon: will NOT perform scalar product with phonon'
  write (*,*) '  displacement vectors in read_gkk. doscalprod==0'
 else if (elph_ds%doscalprod==1) then
  write (*,*) ' elphon: will perform scalar product with phonon'
  write (*,*) '  displacement vectors in read_gkk. doscalprod==1'
 else
  write (*,*) ' elphon: illegal value for doscalprod'
  stop
 end if

 call get_all_gkq (acell,amu,elph_ds,FSfullpqtofull,FSfulltofull,FSkpt,elph_ds%FSintweight,&
& gkk_flag,gprimd,indsym,mpert,natom,nband,nqptirred,nsym,ntypat,n1wf,onegkksize,  &
& phon_ds,qptirred,qpttoqpt,rprimd,spqpt,symrec,symrel,timrev,tnons,typat,ucvol,unitgkk,xred)

 if (elph_tr_ds%ifltransport==1 )then

  write (*,*) 'Transport calculation: allocate gkk_qpt_trin/out and call get_gkk_qpt_tr'

! crc 
  ngkkband=elph_ds%ngkkband
! TODO: do not allocate if arrays are paged to disk!
  if (elph_ds%gkqwrite == 0) then
!  allocate (elph_tr_ds%gkk_qpt_tr(2,ngkkband*ngkkband,elph_ds%nbranch*elph_ds%nbranch,&
!  &            elph_ds%nFSkpt,elph_ds%nsppol,elph_ds%nqpt),stat=istat)
!  if (istat /= 0) stop 'elphon: error in allocating '
   allocate (elph_tr_ds%gkk_qpt_trin(2,ngkkband*ngkkband,elph_ds%nbranch*elph_ds%nbranch,&
&   elph_ds%nFSkpt,elph_ds%nsppol,elph_ds%nqpt),stat=istat)
   if (istat /= 0) stop 'elphon: error in allocating elph_tr_ds%gkk_qpt_trin'
   allocate (elph_tr_ds%gkk_qpt_trout(2,ngkkband*ngkkband,elph_ds%nbranch*elph_ds%nbranch,&
&   elph_ds%nFSkpt,elph_ds%nsppol,elph_ds%nqpt),stat=istat)
   if (istat /= 0) stop 'elphon: error in allocating elph_tr_ds%gkk_qpt_trout'
!  crc
  end if
  
  call get_gkk_qpt_tr(elph_ds,mpi_enreg,nband,hdr,FSfulltoirred, FSirredtoGS,FSfullpqtofull,elph_tr_ds)

 end if

!MG20060531
!============================================================================
!Evaluate lambda and omega_log using the weighted sum over the irred q-points
!found in the GKK file. All the data we need are stored in elph_ds%qgrid_data
!============================================================================
 if (elph_ds%gkqexist == 0) then ! uninterpolated ph linwidths are stored in elph_ds
  write(message,'(a)')' elphon: calling outelph '
  call wrtout(06,message,'COLL')
  fname=trim(elph_ds%elph_base_name) // '_QPTS'
  call outelph(elph_ds,enunit,fname)
 end if
!ENDMG20060531

!========================================================
!Get FS averaged gamma matrices and Fourier transform to real space 
!========================================================
 
 write (message,'(3a)')ch10,' elphon : calling integrate_gamma ',ch10
 call wrtout(06,message,'COLL')

!NOTE: gprim (not gprimd) is used for all FT interpolations,
!to be consistent with the dimensions of the rpt, which come from anaddb.

 call integrate_gamma(elph_ds,FSfullpqtofull,gprim,elph_ds%n0,natom,nrpt,rpt,spqpt,wghatm)
 
 write(message,'(3a)')ch10,'elphon : out of integrate_gamma',ch10
 call wrtout(06,message,'COLL')

!============================================
!Get q path for Fourier interpolation where 
!ph linewidths are going to be interpolated
!TODO MJV 5/8/2007 could be made into a separate subroutine
!
!from HERE
!
!============================================

 nsegment = nqpath-1
 allocate (npoint(nsegment),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating npoint'

!get reference squared distance |q2-q1|^2 and calculate npoint for each segment
!MG In the previous version 5.1.3  there was a bug, since we were evaluating npoint as
!(tmpdist/refdist)^2 instead of the correct ratio (tmpdist/refdist)

 tmpkpt(:) = qpath(:,2)-qpath(:,1)
 refdist =  gmet(1,1)*tmpkpt(1)*tmpkpt(1) &
& +gmet(2,2)*tmpkpt(2)*tmpkpt(2) &
& +gmet(3,3)*tmpkpt(3)*tmpkpt(3) &
& +two*gmet(1,2)*tmpkpt(1)*tmpkpt(2) &
& +two*gmet(1,3)*tmpkpt(1)*tmpkpt(3) &
& +two*gmet(2,3)*tmpkpt(2)*tmpkpt(3)

 if (refdist <= tol10 ) then
  write (message,'(8a)')ch10,&
&  ' elphon : ERROR-',ch10,&
&  ' The first and second q point of the path are equal ',ch10,&
&  ' Impossible to calculate the reference distance. ',ch10,   &
&  ' Modify the q path in your input file '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 factor = 100.0_dp
 do ii=1,nqpath-1
  tmpkpt(:) = qpath(:,ii+1)-qpath(:,ii)
  tmpdist =   gmet(1,1)*tmpkpt(1)*tmpkpt(1) &
&  +gmet(2,2)*tmpkpt(2)*tmpkpt(2) &
&  +gmet(3,3)*tmpkpt(3)*tmpkpt(3) &
&  +two*gmet(1,2)*tmpkpt(1)*tmpkpt(2) &
&  +two*gmet(1,3)*tmpkpt(1)*tmpkpt(3) &
&  +two*gmet(2,3)*tmpkpt(2)*tmpkpt(3)
  npoint(ii) = nint(factor*(sqrt(tmpdist/refdist)))
 end do

!
!to HERE: end sub for qpath creation
!

!==========================================================
!calculate transport matrix elements, integrated over FS
!==========================================================
 if (elph_tr_ds%ifltransport==1 )then
  call integrate_gamma_tr(elph_ds,FSfullpqtofull,gprim,natom,nrpt,rpt,spqpt,wghatm,elph_tr_ds)
 end if

 write (*,*) ' elphon : npoint = ', npoint

!==============================================================
!Calculate phonon linewidths, interpolating on chosen qpoints
!==============================================================

 call mkph_linwid(elph_ds,elph_ds%FSintweight,FSfulltoirred,FSirredtofull,gmet,gprim,gprimd,&
& elph_ds%n0,natom,npoint,nrpt,nsegment,nsym,phon_ds,qpath,qpttoqpt,rpt,spqpt,wghatm)

!==============================================================
!the nesting factor calculation
!==============================================================
 if (prtnest==1 .or. prtnest==2) then
  write(message,'(a)')' elphon : calling mknesting to interpolate the nesting factor'
  call wrtout(06,message,'COLL')

  call mknesting(elph_ds%nFSkpt,FSkpt,kptrlatt,elph_ds%nFSband,elph_ds%FSintweight,    &
&  nsegment,npoint,qpath,elph_ds%elph_base_name,gprim,gprimd,rprim,brav,prtnest)
 end if
 
!======================================================
!Calculate alpha^2 F integrating over fine FSkpt grid
!======================================================

 allocate(a2f_1d(elph_ds%na2f),dos_phon(elph_ds%na2f),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating a2f_1d,dos_phon'
 
 write(message,'(a)')' elphon : calling mka2f '
 call wrtout(06,message,'COLL')

 call mka2f(a2f_1d,dos_phon,elph_ds,elph_ds%FSintweight,FSirredtofull,FSirredwtk,FSkpt,gprim,gprimd,mustar,&
& elph_ds%n0,natom,nrpt,phon_ds,rpt,wghatm)

!calculate transport spectral function and coefficients
 if (elph_tr_ds%ifltransport==1 )then
  call mka2f_tr(elph_ds,FSkpt,gprim,elph_ds%n0,ucvol,natom,nrpt,phon_ds,rpt,wghatm,elph_tr_ds)

! if transport gkk are kept in memory, deallocate now.
  if (elph_ds%gkqwrite == 0) then
   deallocate (elph_tr_ds%gkk_qpt_trin)
   deallocate (elph_tr_ds%gkk_qpt_trout)
  end if
 end if

!MG: evaluate a2F only using the input Q-grid (without using interpolated matrices)
!SCOPE: test the validity of the Fourier interpolation
!
!NOTE MJV: why limit to case gkqexist == 0 ?? Just need to implement the io?
!
 if (elph_ds%gkqexist == 0 )then
  write(message,'(a)')' elphon : calling mka2fQgrid'
  call wrtout(06,message,'COLL')
  fname=trim(elph_ds%elph_base_name) // '_A2F_QGRID'
  call mka2fQgrid(elph_ds,fname)
 end if
!ENDMG

!=============================================
!Eliashberg equation in 1-D (isotropic case)
!=============================================

 write (message,'(2a)')ch10,&
& ' Solving the 1-D Eliashberg equation (isotropic case)'
 call wrtout(06,message,'COLL')

 call eliashberg_1d(a2f_1d,elph_ds,mustar,elph_ds%n0,natom,nsym,phon_ds)

 deallocate (a2f_1d,dos_phon)

!MJV: 20070805 should exit here. None of the rest is tested or used yet to my knowledge

!========================================================================
!Now gkk contains the matrix elements of dH(1)/dxi i=1,2,3
!for kpoints on the FS but qpoints only in the given grid {Q}.
!
!1.) Need to complete the gkk elements for q and k\prime=k+q not 
!in the set of {k+Q} by Fourier interpolation on the Q.
!
!2.) Need to complete the dynamical matrices and phonon freqs for
!all q between points on the FS.
!
!3.) With the eigenvectors e_ph of the dyn mats, do the scalar product
!e_ph . gkk, which implies the gkk are turned to the eigenbasis of
!the phonons. Before the (non eigen-) modes are ordered
!atom1 xred1 atom1 xred2 atom1 xred3
!atom2 xred1 atom2 xred2 atom2 xred3 ...
!=======================================================================
!==========================================================
!FT of recip space gkk matrices to real space (gkk_rpt)
!NOTE: could be made into FFT, couldnt it? If shifts are
!used with a homogeneous grid
!==========================================================
 write (message,'(2a,i5)')ch10,&
& ' elphon : Fourier transform (q --> r) of the gkk matrices using nrpt = ',nrpt
 call wrtout(06,message,'COLL')

 call get_all_gkr(elph_ds,gprim,natom,nrpt,onegkksize,rpt,spqpt,wghatm)

!=========================================================
!complete gkk2 for all qpts between points
!on full kpt grid (interpolation from real space values)
!=========================================================

 make_gkk2=.false.
 
 if (.not. make_gkk2) then
  write(message,'(2a)')ch10,&
&  ' elphon : skipping full g(k,k") interpolation '
  call wrtout(06,message,'COLL')
 else

  write(message,'(2a)')ch10,&
&  ' elphon : Calling get_all_gkk2 to calculate gkk2 for q points over the full k grid'
  call wrtout(06,message,'COLL')
  
  call get_all_gkk2(acell,amu,atmfrc,dielt,dipdip,dyewq0,elph_ds,FSkptirred,FSkpt,   &
&  ftwghtgkk,gmet,gprim,indsym,mpert,msym,natom,nrpt,nsym,ntypat,   &
&  onegkksize,phon_ds,rcan,rmet,rprim,rprimd,rpt,spqpt,symrel,trans,&
&  typat,ucvol,wghatm,xred,zeff)
 end if

 allocate(zz(elph_ds%na2f,elph_ds%nFSkpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating zz'
 allocate(delta(elph_ds%na2f,elph_ds%nFSkpt),stat=istat)
 if (istat /= 0) stop 'elphon: error in allocating delta'
 
!=====================================================
!Here should be the anisotropic Eliashberg equations.
!=====================================================

!initialize delta function

!initialize T_c

!initialize delta

!iterate for calculation of T_c
 do ieliash=1,neliash

! ===========================
! calculate lambda function
! ===========================

! ========================================
! integrate lambda over FS -> Z function
! ========================================
  
! ========================================
! integrate delta*Z over FS -> new delta
! ========================================
  
! update T_c

 end do

!end iterate ieliash

!output T_c and related quantities


!clean and deallocate junk

 deallocate (elph_ds%wtq)
 deallocate(elph_ds%spqpt)

 call clean_phon_ds(phon_ds)

 call hdr_clean(hdr)

 end subroutine elphon
!!***
