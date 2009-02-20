!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkph_linwid
!!
!! NAME
!! mkph_linwid
!!
!! FUNCTION
!!  Calculate the phonon linewidths on a trajectory in q space
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  elph_ds = datastructure with phonon matrix elements
!!  FSintweight = integration weights on the FS
!!  FSfulltoirred = mapping of full FS kpts to irreducible ones
!!  FSirredtofull = indices of irreducible kpoints in the full ones
!!  gkk2 = gkk2 matrix elements on full FS grid for each phonon mode
!!  gmet = reciprocal space metric
!!  gprim = reciprocal lattice vectors
!!  gprimd = reciprocal-space lattice vectors (dimensionful)
!!  n0 = DOS at the Fermi level calculated from the FSkpt integration weights
!!  natom = number of atoms
!!  npoint_in = number of points requested along trajectory
!!  nrpt = number of real space points for FT interpolation
!!  nsegment_in = number of segments in reciprocal space trajectory
!!  nsym = number of symops
!!  phon_ds = datastructure with interatomic force constants
!!  qpath_vertices_in = vertices of reciprocal space trajectory
!!  qpttoqpt = mapping of qpoints under symops
!!  rpt = coordinates of real space points for FT interpolation
!!  spqpt = coordinates of qpoints
!!  wghatm = weights of pairs of atoms for FT interpolation
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      canon9,ftgam,inpphon,leave_new,wrtout,zgemm
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mkph_linwid(elph_ds,FSintweight,FSfulltoirred,FSirredtofull,&
& gmet,gprim,gprimd,n0,natom,npoint_in,nrpt,nsegment_in,nsym,phon_ds,&
& qpath_vertices_in,qpttoqpt,rpt,spqpt,wghatm)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_17ddb, except_this_one => mkph_linwid
!End of the abilint section

 implicit none

!Arguments ------------------------------------
  ! needed for phonon interpolation
!scalars
 integer,intent(in) :: natom,nrpt,nsegment_in,nsym
 type(elph_type),intent(inout) :: elph_ds
 type(phon_type),intent(inout) :: phon_ds
!arrays
 integer,intent(in) :: FSfulltoirred(3,elph_ds%nFSkpt)
 integer,intent(in) :: FSirredtofull(elph_ds%nFSkptirred)
 integer,intent(in) :: npoint_in(nsegment_in),qpttoqpt(2,nsym,elph_ds%nqpt)
 real(dp),intent(in) :: FSintweight(elph_ds%nFSband,elph_ds%nFSkpt,elph_ds%nsppol)
 real(dp),intent(in) :: gmet(3,3),gprim(3,3),gprimd(3,3),n0(elph_ds%nsppol)
 real(dp),intent(in) :: qpath_vertices_in(3,nsegment_in+1),rpt(3,nrpt)
 real(dp),intent(in) :: spqpt(3,elph_ds%nqpt),wghatm(natom,natom,nrpt)

!Local variables-------------------------------
  ! for diagonalization of gammma matrix
  ! output variables for gtdyn9+phfrq3
!scalars
 integer :: eivec=1,i1,i2,iFSkpt1,iFSkpt2,iFSkpt3,iatom,ib1,ib2,ibranch,idir
 integer :: ieqFSkpt1,ier,ii,indx,iost,ip,ipert1,ipert2,ipoint,ipp,iqpt
 integer :: iqptfull,irpt,iseg,isppol,jbranch,k1,kbranch,kdir,mu,nsegment,nu
 integer :: qtor,unit_bs,unit_lambda,unit_lwd
 real(dp) :: diagerr,gaussprefactor,gaussval,phnow,qphnrm=one,res,total_weight
 real(dp) :: weight
 character(len=500) :: message
 character(len=fnlen) :: fname
!arrays
 integer :: npoint(nsegment_in+1)
 real(dp),parameter :: c0(2)=(/0._dp,0._dp/),c1(2)=(/1._dp,0._dp/)
 real(dp) :: displ(2,3*natom,3*natom)
 real(dp) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: eigval(elph_ds%nbranch),eigvec(2*elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: elph_linwid_integ(elph_ds%nbranch,elph_ds%nqpt)
 real(dp) :: gam_now(2,elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: gam_rpt(2,elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: gammafact(elph_ds%nbranch),imeigval(elph_ds%nbranch)
 real(dp) :: lambda(elph_ds%nbranch),pheigval(3*natom)
 real(dp) :: pheigvec(2*3*natom*3*natom),phfrq_tmp(3*natom)
 real(dp) :: qpath_vertices(3,nsegment_in+2),qpt(3),redkpt(3)
 real(dp) :: tmp_gkk2(elph_ds%nbranch,elph_ds%nFSband,elph_ds%nFSband,elph_ds%nFSkpt)
 real(dp) :: tmpelph_linwid(elph_ds%nbranch)
 real(dp) :: tmpgam1(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmpgam2(2,elph_ds%nbranch,elph_ds%nbranch),tmpqpt(3)
 real(dp),allocatable :: matrx(:,:),zhpev1(:,:),zhpev2(:)

! *********************************************************************
!calculate elph_linwid for qpoints along path defined by qpath_vertices

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zhpev
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zgemm
#endif

 isppol=1

 write (*,*) 'mkph_linwid : enter '

!===================================================================
!Definition of the q path along which ph linwid will be interpolated
!===================================================================
!add extra segment for last point: nvertices = nsegments + 1
 qpath_vertices(:,1:nsegment_in+1) = qpath_vertices_in(:,:)
 qpath_vertices(:,nsegment_in+2) = qpath_vertices_in(:,nsegment_in+1)
 npoint(1:nsegment_in) = npoint_in(:)
 npoint(nsegment_in+1) = 1
 nsegment=nsegment_in + 1

!==========================================================
!Open _LWD file and write header
!==========================================================
 unit_lwd=108
 fname=trim(elph_ds%elph_base_name) // '_LWD'
 open (unit=unit_lwd,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
  write (message,'(2a)')' mkph_linwid : ERROR- opening file ',trim(fname)
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 write (unit_lwd,'(a)')'#'
 write (unit_lwd,'(a)')'# ABINIT package : Phonon linewidth file'
 write (unit_lwd,'(a)')'#'
 write (unit_lwd,'(a,i10,a)') '#  Phonon linewidths calculated on ',sum(npoint), &
& ' points along the qpath'
 write (unit_lwd,'(a)')'# Description of the Q-path :'
 write (unit_lwd, '(a,i10)')'# Number of line segments = ',nsegment_in
!write (unit_lwd, '(a)') '# starting and ending points = '
!write (unit_lwd, '(a,3(E16.6,1x))') "#  ", qpath_vertices_in(:,1)
!write (unit_lwd, '(a,3(E16.6,1x))') "#  ", qpath_vertices_in(:,nsegment_in+1)
 write (unit_lwd,'(a)')'# Vertices of the Q-path and corresponding index = '
 indx=1
 do ii=1,nsegment_in+1
  write (unit_lwd,'(a,3(e16.6,1x),i8)')'#  ',qpath_vertices_in(:,ii),indx
  indx=indx+npoint(ii)
 end do
 write (unit_lwd,'(a)')'#'

!==========================================================
!Open _BST file and write header
!==========================================================
 unit_bs=109
 fname=trim(elph_ds%elph_base_name) // '_BST'
 open (unit=unit_bs,file=fname,status='unknown')
 write (unit_bs, '(a)') '#'
 write (unit_bs, '(a)') '# ABINIT package : Phonon band structure file'
 write (unit_bs, '(a)') '#'
 write (unit_bs, '(a,I10,a)') '#  Phonon BS calculated on ', sum(npoint), &
& ' points along the qpath'
 write (unit_bs, '(a,I10)') '# Number of line segments = ', nsegment_in
!write (unit_bs, '(a)') '# starting and ending points = '
!write (unit_bs, '(a,3(E16.6,1x))') "#  ", qpath_vertices_in(:,1)
!write (unit_bs, '(a,3(E16.6,1x))') "#  ", qpath_vertices_in(:,nsegment_in+1)
!write (unit_bs, '(a)') '#'
 indx=1
 do ii=1,nsegment_in+1
  write (unit_bs,'(a,3(E16.6,1x),i8)')'#  ',qpath_vertices_in(:,ii),indx
  indx=indx+npoint(ii)
 end do
 write (unit_bs,'(a)')'#'

!MG20060606
!==========================================================
!open _LAMBDA file and write header
!contains \omega(q,n) and \lambda(q,n) and can be plotted using xmgrace
!==========================================================
 unit_lambda=110
 fname=trim(elph_ds%elph_base_name) // '_LAMBDA'
 open (unit=unit_lambda,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
  write (message,'(2a)')' mkph_linwid : ERROR- opening file ',trim(fname)
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 write (unit_lambda,'(a)')'#'
 write (unit_lambda,'(a)')'# ABINIT package : Lambda file'
 write (unit_lambda,'(a)')'#'
 write (unit_lambda,'(a,I10,a)')'#  Lambda(q,nu) calculated on ',sum(npoint),&
& ' Q-points'
 write (unit_lambda,'(a)')'# Description of the Q-path :'
 write (unit_lambda,'(a,I10)')'# Number of line segments = ',nsegment_in
 write (unit_lambda,'(a)')'# Vertices of the Q-path and corresponding index = '

 indx=1
 do ii=1,nsegment_in+1
  write (unit_lambda,'(a,3(E16.6,1x),i8)')'#  ',qpath_vertices_in(:,ii),indx
  indx=indx+npoint(ii)
 end do
 write (unit_lambda,'(a)')'#'
 write (unit_lambda,'(a)')'# index frequency lambda(q,n) frequency lambda(q,n) .... lambda_tot'
 write (unit_lambda,'(a)')'#'

!real space to q space
 qtor=0

!initialize the maximal phonon frequency
 elph_ds%omega_min = zero
 elph_ds%omega_max = zero

 write (*,*) ' mkph_linwid : shape(elph_ds%gamma_qpt) = ',&
& shape(elph_ds%gamma_qpt)


!
!Big do loop over spin polarizations
!could put in locally, so phonon stuff is not done twice...
!
 do isppol=1,elph_ds%nsppol

  indx=1

! DEBUG
! write (100,*) '# gamma_qpt on qpts actually calculated = '
! do iqpt=1,elph_ds%nqpt
! write (100,*) iqpt, spqpt(:,iqpt)
! write (100,'(9(2E14.6,1x))') elph_ds%gamma_qpt(:,:,isppol,iqpt)
! end do
! write (99,*) '# gamma_rpt on rpts actually calculated = '
! do irpt=1,nrpt
! write (99,*) irpt, rpt(:,irpt)
! write (99,'(9(2E14.6,1x))') elph_ds%gamma_rpt(:,:,isppol,irpt)
! end do
! ENDDEBUG

! Back FT real space integrated values to recip space
! DEBUG
! write (95,*) '# difference btw gamma_qpt and FTed from real space = '
! qtor = 0
! do iqpt=1,elph_ds%nqpt
! call ftgam(wghatm,gam_now,elph_ds%gamma_rpt,gprim,natom,1,nrpt,qtor,rpt,spqpt(:,iqpt))
! write (95,'(6E16.4,2x)') gam_now(:,:)-elph_ds%gamma_qpt(:,:,isppol,iqpt)
! end do
! write (94,*) '# difference btw gamma_rpt and back FTed from recip space = '
! qtor = 1
! do irpt=1,nrpt
! call ftgam(wghatm,elph_ds%gamma_qpt,gam_rpt,gprim,natom,elph_ds%nqpt,1,qtor,rpt(:,irpt),spqpt)
! write (94,'(6E16.4,2x)') gam_rpt(:,:)-elph_ds%gamma_rpt(:,:,isppol,irpt)
! elph_ds%gamma_rpt(:,:,isppol,irpt) = gam_rpt(:,:)
! end do
! write (93,*) '# explicit back then forward ftgam is done '
! write (93,*) '# elph_ds%gamma_qpt to gam_rpt back to gam_now '
! qtor = 0
! do iqpt=1,elph_ds%nqpt
! call ftgam(wghatm,gam_now,elph_ds%gamma_rpt,gprim,natom,1,nrpt,qtor,rpt,spqpt(:,iqpt))
! write (93,'(6E16.4,2x)') gam_now(:,:)-elph_ds%gamma_qpt(:,:,isppol,iqpt)
! end do
! ENDDEBUG


! Interpolation along specified path in q space

! DEBUG
! write (100,*) '#-------------------------------------------------------------'
! write (100,*) '#NOW INTERPOLATED VALUES OF gamma MATRICES'
! write (100,*) '#-------------------------------------------------------------'
! 
! write (96,*) '#displ_red vector at each qpoint :'
! write (97,*) '#displ vector at each qpoint :'
! write (98,*) '#bare gamma :'
! ENDDEBUG

! Output to the main output file
  write(message,'(a,a)')ch10,&
&  ' Output of the linewidths for the first point of each segment. Linewidths are given in Hartree.'
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  write (*,*) ' mkph_linwid : elph_ds%doscalprod = ', elph_ds%doscalprod

  qtor = 0
  do iseg=1,nsegment
   do ipoint=0,npoint(iseg)-1

!   Get qpoint along the path from qpath_vertices(:,iseg)
!   to the point just before qpath_vertices(:,iseg+1)
!   DEBUG
!   write (*,*) 'mkph_linwid : i and qpts ', iseg, ipoint, npoint(iseg)
!   ENDDEBUG

    qpt(:) = qpath_vertices(:,iseg) + dble(ipoint)/dble(npoint(iseg))*&
&    (qpath_vertices(:,iseg+1)-qpath_vertices(:,iseg))

!   DEBUG
!   write (140,'(3E16.6,1x)') qpt(:)
!   write (140,'(3(3E16.6,1x))') qpath_vertices(:,iseg), qpt(:), qpath_vertices(:,iseg+1)
!   ENDDEBUG

    call canon9(qpt(1),redkpt(1),res)
    call canon9(qpt(2),redkpt(2),res)
    call canon9(qpt(3),redkpt(3),res)
    qpt(:) = redkpt(:)

!   This reduced version of ftgkk supposes the kpoints have been integrated
!   in integrate_gamma. Do FT from real-space gamma grid to 1 qpt.
    call ftgam(wghatm,gam_now,elph_ds%gamma_rpt(:,:,isppol,:),gprim,natom,1,nrpt,qtor,rpt,qpt)

!   DEBUG
!   write (98,'(a)')  '# just after ftgam'
!   write (98,'(9(2E14.6,1x))')  gam_now(:,:)
!   ENDDEBUG

!   
!   get phonon freqs and eigenvectors anyway
!   
    call inpphon(displ,pheigval,pheigvec,phfrq_tmp,phon_ds,qpt)

!   
!   additional frequency factor for some cases
!   
    gammafact(:) = one
!   if the matrices contain only the bare gkk matrix elements:
!   MJV 27/2/2006 I think this is wrong, but tkeepbands not used
!   much yet anyway - for full band dependency of gkk
    if (elph_ds%tkeepbands == 1) then
!    gammafact(:) = two*pi*phfrq_tmp(:)
!    gammafact(:) = two*pi*abs(phfrq_tmp(:))
!    MG20060607 FIXME
!    I have modified nmsq_gam_sumfs and nmsq_gam to calculate the matrix
!    gkk_qpt using the same approach as in nmsq_pure_gkk, in this case gammafact is 1
!    this part can be removed after careful testing of the new implementation
     gammafact(:)=one
!    ENDMG
    end if


!   
!   if the matrices do not contain the scalar product
!   with the displ vectors yet do it now

    if (elph_ds%doscalprod == 0) then

!    DEBUG
!    write (*,*) 'phonon freqs ', phfrq_tmp
!    ENDDEBUG

     displ_red(:,:,:) = zero
     do jbranch=1,elph_ds%nbranch
      do iatom=1,natom
       do idir=1,3
        ibranch=idir+3*(iatom-1)
        do kdir=1,3
         k1 = kdir+3*(iatom-1)
         displ_red(1,ibranch,jbranch) = displ_red(1,ibranch,jbranch)&
&         + gprimd(kdir,idir)*displ(1,k1,jbranch)
         displ_red(2,ibranch,jbranch) = displ_red(2,ibranch,jbranch)&
&         + gprimd(kdir,idir)*displ(2,k1,jbranch)
        end do
       end do
      end do
     end do

!    DEBUG
!    write (97,'(9(2E14.6,1x))')  displ(:,:,:)
!    write (96,'(9(2E14.6,1x))')  displ_red(:,:,:)
!    write (98,'(9(2E14.6,1x))')  gam_now(:,:)
!    ENDDEBUG

     eigval(:) = zero
     imeigval(:) = zero

!    calculate displ_red* gam_now* displ_red^{*T}
     do jbranch=1,elph_ds%nbranch
      do ipert1=1,elph_ds%nbranch
       do ipert2=1,elph_ds%nbranch
        ipp = ipert2+(ipert1-1)*elph_ds%nbranch
        eigval(jbranch) = eigval(jbranch)&
&        + displ_red(1,ipert2,jbranch)*gam_now(1,ipp)*displ_red(1,ipert1,jbranch)&
&        - displ_red(2,ipert2,jbranch)*gam_now(2,ipp)*displ_red(1,ipert1,jbranch)&
&        + displ_red(1,ipert2,jbranch)*gam_now(2,ipp)*displ_red(2,ipert1,jbranch)&
&        + displ_red(2,ipert2,jbranch)*gam_now(1,ipp)*displ_red(2,ipert1,jbranch)

!       MG20060603   Fixed a minor bug present in version 5.1.3 (this quantity indeed should be zero and is not used!)
        imeigval(jbranch) = imeigval(jbranch)&
&        + displ_red(1,ipert2,jbranch)*gam_now(2,ipp)*displ_red(1,ipert1,jbranch)&
&        - displ_red(1,ipert2,jbranch)*gam_now(1,ipp)*displ_red(2,ipert1,jbranch)&
&        + displ_red(2,ipert2,jbranch)*gam_now(1,ipp)*displ_red(1,ipert1,jbranch)&
&        + displ_red(2,ipert2,jbranch)*gam_now(2,ipp)*displ_red(2,ipert1,jbranch)
!       ENDMG
       end do
      end do

      if (abs(imeigval(jbranch)) > tol8) then
       write (message,'(3a,i6,a,es16.8)')&
&       ' mkph_linwid : WARNING-  imaginary values! ',ch10,&
&       ' branch = ',jbranch,' imeigval = ',imeigval(jbranch)
       call wrtout(06,message,'COLL')
      end if

     end do

!    
!    if elph_ds%doscalprod is 1
!    
    else if (elph_ds%doscalprod == 1) then

!    DEBUG
!    write (100,'(9(2E14.6,1x))') gam_now(:,:)
!    ENDDEBUG

!    Diagonalize gamma matrix at qpoint (complex matrix).

!    !   Copied from phfrq33
!    !! Useless with new formulation and exact displ vector at qpoint :
!    !!    gam_now is already diagonal if doscalprod==0
!    ier=0
!    ii=1
!    allocate(matrx(2,(3*natom*(3*natom+1))/2))
!    do i2=1,3*natom
!    do i1=1,i2
!    ipp=i1+(i2-1)*3*natom
!    matrx(1,ii)=gam_now(1,ipp)
!    matrx(2,ii)=gam_now(2,ipp)
!    ii=ii+1
!    end do
!    end do
!    allocate(zhpev1(2,2*3*natom-1),zhpev2(3*3*natom-2))
!    #if defined T3E
!    call CHPEV ('V','U',3*natom,matrx,eigval,eigvec,3*natom,zhpev1,&
!    &    zhpev2,ier)
!    #else
!    call ZHPEV ('V','U',3*natom,matrx,eigval,eigvec,3*natom,zhpev1,&
!    &    zhpev2,ier)
!    #endif
!    deallocate(matrx,zhpev1,zhpev2)
!    
!    
!    using phonon eigenvectors U, gamma_diag = U^{T*} gam_now U
!    
!    DEBUG
!    write (*,'(3(2(E20.10,1x)))') gam_now, pheigvec
!    write (*,'(3(2(E20.10,1x)))') CMPLX(one,zero), CMPLX(zero,zero)
!    write (*,*) 'N', 3*natom
!    ENDDEBUG

!    MJV NOTE: gam_now is recast implicitly here to matrix 
     call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, c1, gam_now, 3*natom,&
&     pheigvec, 3*natom, c0, tmpgam1, 3*natom)
     call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, c1, pheigvec, 3*natom,&
&     tmpgam1, 3*natom, c0, tmpgam2, 3*natom)

!    DEBUG
!    write (130,*) '# mkph_linwid : gamma diagonalized with phonon eigenvectors '
!    write (130,'(3(2(E20.10,1x)))') tmpgam2
!    ENDDEBUG

     diagerr = zero
     do ibranch=1,elph_ds%nbranch

      eigval(ibranch) = tmpgam2(1,ibranch,ibranch)

      do jbranch=1,ibranch-1
       diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))+abs(tmpgam2(2,jbranch,ibranch))
      end do
      do jbranch=ibranch+1,elph_ds%nbranch
       diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))+abs(tmpgam2(2,jbranch,ibranch))
      end do
      diagerr = diagerr + abs(tmpgam2(2,ibranch,ibranch))
     end do

     if (diagerr > tol12) then
      write (*,*) 'mkph_linwid: error in diagonalization of gamma with phon eigenvectors: ', diagerr
     end if

    else

     write (message,'(3a,i4)')' mkph_linwid: BUG-',ch10,&
&     ' Wrong value for elph_ds%doscalprod = ',elph_ds%doscalprod
     call wrtout(06,message,'COLL')
     call leave_new('COLL')

    end if
!   end elph_ds%doscalprod if

!   ==========================================================
!   write data to files for each q point
!   ==========================================================

!   write (*,'(a,i5,i5,3E16.5,2x,6E16.6)',ADVANCE='NO') '---', iseg, ipoint, qpt, eigval
    write (unit_lwd,'(i5)', advance='no') indx
    write (unit_lwd,'(18E16.5)',advance='no') (eigval(ii)*gammafact(ii),ii=1,elph_ds%nbranch)
    write (unit_lwd,*)

!   only print phonon BS for isppol 1: independent of electron spins
    if (isppol==1) then
     write (unit_bs,'(i5)', advance='no') indx
     write (unit_bs,'(18E16.5)',advance='no') phfrq_tmp
     write (unit_bs,*)
    end if

!   MG20060606
    write (unit_lambda,'(i5)', advance='no') indx
    do ii=1,elph_ds%nbranch
     lambda(ii)=zero
     if (abs(phfrq_tmp(ii)) > tol10) lambda(ii)=eigval(ii)*gammafact(ii)/(pi*elph_ds%n0(isppol)*phfrq_tmp(ii)**2)
     write (unit_lambda,'(18es16.8)',advance='no')phfrq_tmp(ii),lambda(ii)
    end do
    write (unit_lambda,'(es16.8)',advance='no') sum(lambda)
    write (unit_lambda,*)
!   ENDMG

!   MG NOTE: I wrote a piece of code to output all these quantities using units
!   chosen by the user, maybe in version 5.2?
!   In this version the output of lambda(q,\nu) has been added

!   Output to the main output file, for first point in segment
    if(ipoint==0)then
     write(message,'(a,a,3es16.6,a,i4,a,a)')ch10,&
&     ' Q point =',qpt(:),'   isppol = ',isppol,ch10,&
&     ' Mode number    Frequency (Ha)  Linewidth (Ha)  Lambda(q,n)'
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,message,'COLL')
     do ii=1,elph_ds%nbranch
      write(message,'(i8,es20.6,2es16.6)' )ii,phfrq_tmp(ii),eigval(ii)*gammafact(ii),lambda(ii)
      call wrtout(6,message,'COLL')
      call wrtout(ab_out,message,'COLL')
     end do
    end if

!   find max/min phonon frequency along path chosen
!   presumed to be representative of full BZ to within 10 percent
    elph_ds%omega_min = min(elph_ds%omega_min,1.1_dp*phfrq_tmp(1))
    elph_ds%omega_max = max(elph_ds%omega_max,1.1_dp*phfrq_tmp(elph_ds%nbranch))

    indx = indx+1

   end do
!  end ipoint do
  end do
! end iseg do

! add blank lines to output files between sppol
  write(message,'(a)' ) ''
  call wrtout(unit_lwd,message,'COLL')
  call wrtout(unit_lambda,message,'COLL')
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

 end do
!end isppol do

 close (unit=unit_lwd)
 close (unit=unit_bs)
 close (unit=unit_lambda)

 write (*,*) ' elph_linwid : omega_min, omega_max = ', elph_ds%omega_min, elph_ds%omega_max

 write (*,*) ' elph_linwid : end '


end subroutine mkph_linwid
!!***
