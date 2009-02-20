!{\src2tex{textfont=tt}}
!!****f* ABINIT/initberry
!! NAME
!! initberry
!!
!! FUNCTION
!! Initialization of Berryphase calculation of the polarization, the
!! ddk and the response of an insulator to a homogenous electric field.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVeithen).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtfil <type(datafiles_type)> = variables related to files
!!  dtset <type(dataset_type)> = all input variables in this dataset
!!  gmet(3,3) = reciprocal space metric tensor in bohr**-2
!!  kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!!  mband = maximum number of bands
!!  mkmem = maximum number of k-points in core memory
!!  mpw = maximum number of plane waves
!!  nkpt = number of k points
!!  npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!!  nsppol = 1 for unpolarized, 2 for spin-polarized
!!  nsym = number of symmetry operations
!!  occ(mband*nkpt*nsppol) = occup number for each band at each k point
!!  rprimd(3,3) = dimensional primitive vectors
!!  symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!    reciprocal space primitive translations
!!
!! OUTPUT
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations
!!  pwind(pwind_alloc,2,3) = array used to compute the overlap matrix smat
!!                         between k-points k and k +- dk where dk is
!!                         parallel to the direction idir
!!    jpw = pwind(ipw,ifor,idir)
!!      * ipw = index of plane wave vector G for a given k-point k
!!      * ifor = 1: k + dk
!!               2: k - dk
!!      * idir = direction of the polarization/ddk calculation [dk(idir)
!!               is the only non-zero element of dk(:)]
!!      * jpw = index of plane wave vector G (+dG) at k +- dk
!!              where dG is a shift of one reciprocal lattice vector
!!              (required to close the strings of k-points using the
!!               periodic gauge condition)
!!    In case a G-vector of the basis sphere of plane waves at k
!!    does not belong to the basis sphere of plane waves at k+dk, jpw = 0.
!!   pwind_alloc = first dimension of pwind and pwnsfac
!!   pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!
!! SIDE EFFECTS
!!  mpi_enreg = informations about MPI parallelization
!!    kptdstrb(nproc,nneighbour,fmkmem_max*nsppol) : Array required
!!      by berryphase_new.f for MPI // over k-points. Defined
!!      for k-points in the fBZ
!!    kptdstrbi(nproc,nneighbour,mkmem_max*nsppol) : Same as kptdstrb
!!      but for k-points in the iBZ. Used by vtorho.f
!!           nproc = number of cpus
!!           nneighbour = number of neighbours for each k-point (= 6)
!!
!! TO DO
!!
!! NOTES
!! the interface of this routine is mandatory for src/21drive/gstate.F90
!! If you happen to modify it, please update src/defs/defs_berry.F90
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      getkgrid,kpgsph,leave_new,listkk,wrtout,xcomm_world,xmax_mpi
!!      xsum_mpi_int
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine initberry(dtefield,dtfil,dtset,gmet,kg,mband,mkmem,mpi_enreg,&
&              mpw,nkpt,npwarr,nsppol,nsym,occ,pwind,pwind_alloc,pwnsfac,&
&              rprimd,symrec)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13recipspace
 use interfaces_14wfs
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem,mpw,nkpt,nsppol,nsym
 integer,intent(out) :: pwind_alloc
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(efield_type),intent(out) :: dtefield
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt),symrec(3,3,nsym)
 integer,pointer :: pwind(:,:,:)
 real(dp),intent(in) :: gmet(3,3),occ(mband*nkpt*nsppol),rprimd(3,3)
 real(dp),pointer :: pwnsfac(:,:)

!Local variables-------------------------------
!scalars
 integer :: count,exchn2n3d,flag,flag_kpt,fnkpt_computed,iband,icg,idir,idum
 integer :: idum1,idumi,idumj,ierr,ifor,iforrot,ii,ikg,ikg1,ikpt,ikpt1,ikpt1f
 integer :: ikpt1i,ikpt2,ikpt_loc,ikptf,ikpti,ikstr,index,ineigh,ipw,ipwnsfac
 integer :: isppol,istr,istwf_k,isym,isym1,itrs,iunmark,jkstr,jpw,me_g0,mkmem_
 integer :: nband_k,nband_occ_k,nkstr,npw_k,npw_k1,spaceComm
 real(dp) :: ecut_eff,eg,eg_ev,phasedumi,phasedumr,rdum,rdum1
 character(len=500) :: message
!arrays
 integer :: dg(3),dg1(3),dk_(3),dsifkpt(3),iadum(3),iadum1(3),neigh(6)
 integer,allocatable :: aidum(:),aidum2(:),buffer(:),kg1_k(:,:),kpt_mark(:)
 real(dp) :: adum1(3),adum2(3),diffk(3),dk(3),dkrot(3),dum33(3,3),eg_dir(3)
 real(dp) :: kpt1(3),shiftk_(3,8)
 real(dp),allocatable :: dum_fwtk(:)

! *************************************************************************

!DEBUG
!write(*,*)'initberry: enter'
!call leave_new('COLL')
!ENDDEBUG

!----------------------------------------------------------------------------
!-------------------- Obtain k-point grid in the full BZ --------------------
!----------------------------------------------------------------------------

!Prepare first call to getkgrid (obtain number of k points in FBZ)
 dsifkpt(:) = 1
 dtefield%fnkpt = 0
 shiftk_(:,:) = 0._dp
 shiftk_(:,1:dtset%nshiftk) = dtset%shiftk(:,1:dtset%nshiftk)
 dk_(:) = 0

 allocate(dtefield%fkptns(3,dtefield%fnkpt),dum_fwtk(dtefield%fnkpt))

 call getkgrid(dsifkpt,ab_out,5,dtefield%fkptns,3,dtset%kptrlatt,&
& rdum,dtset%nsym,dtefield%fnkpt,fnkpt_computed,dtset%nshiftk,dtset%nsym,&
& rprimd,shiftk_,dtset%symafm,dtset%symrel,dtset%tnons,&
& dk_,dum_fwtk)

 dtefield%fnkpt = fnkpt_computed
 deallocate(dtefield%fkptns,dum_fwtk)

 allocate(dtefield%fkptns(3,dtefield%fnkpt),dum_fwtk(dtefield%fnkpt))

!second call to getkgrid is to obtain the full k point list
 call getkgrid(dsifkpt,ab_out,5,dtefield%fkptns,3,dtset%kptrlatt,&
& rdum,dtset%nsym,dtefield%fnkpt,fnkpt_computed,dtset%nshiftk,dtset%nsym,&
& rprimd,shiftk_,dtset%symafm,dtset%symrel,dtset%tnons,&
& dk_,dum_fwtk)

!call listkk to get mapping from FBZ to IBZ
 rdum=1.0d-5  ! cutoff distance to decide when two k points match
 allocate(dtefield%indkk_f2ibz(dtefield%fnkpt,6))

!ji: The following may need modification in the future
!**** no spin-polarization doubling ; allow use of time reversal symmetry ****

 call listkk(rdum,gmet,dtefield%indkk_f2ibz,dtset%kptns,dtefield%fkptns,nkpt,&
& dtefield%fnkpt,dtset%nsym,1,dtset%symafm,dtset%symrel,1)

!Construct i2fbz and f2ibz
 allocate(dtefield%i2fbz(nkpt))
 idum=0
 do ikpt=1,dtefield%fnkpt
  if (dtefield%indkk_f2ibz(ikpt,2)==1 .and. &
&  dtefield%indkk_f2ibz(ikpt,6) == 0 .and. &
&  maxval(abs(dtefield%indkk_f2ibz(ikpt,3:5))) == 0 ) then
   dtefield%i2fbz(dtefield%indkk_f2ibz(ikpt,1))=ikpt
   idum=idum+1
  end if
 end do
 if (idum/=nkpt)then
  write(message,'(a,a,a,a)')ch10,&
&  ' initberry: ERROR - ',ch10,&
&  '   Found wrong number of k-points in IBZ'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if


!------------------------------------------------------------------------------
!------------------- Compute variables related to MPI // ----------------------
!------------------------------------------------------------------------------

 if (mpi_enreg%paral_compil_kpt == 0) then  ! no MPI //

  dtefield%fmkmem = dtefield%fnkpt
  dtefield%fmkmem_max = dtefield%fnkpt
  dtefield%mkmem_max = nkpt

 else    ! MPI //

! Number of k-points in the FBZ for each cpu

  dtefield%fmkmem = 0
  do ikpt = 1, dtefield%fnkpt
   ikpti = dtefield%indkk_f2ibz(ikpt,1)
   nband_k = dtset%nband(ikpti)
   if (minval(abs(mpi_enreg%proc_distrb(ikpti,1:nband_k,1:nsppol) - &
&   mpi_enreg%me)) == 0) dtefield%fmkmem = dtefield%fmkmem + 1
  end do

! Maximum value of mkmem and fmkmem
! BEGIN TF_CHANGES
  call xcomm_world(mpi_enreg,spaceComm)
! END TF_CHANGES
  call xmax_mpi(dtefield%fmkmem,dtefield%fmkmem_max,spaceComm,ierr)

  mkmem_ = mkmem   ! I have to use the dummy variable mkmem_ because
! mkmem is declared as intent(in) while the first
! argument of xmax_mpi must be intent(inout)

  call xmax_mpi(mkmem_,dtefield%mkmem_max,spaceComm,ierr)

  allocate(mpi_enreg%kptdstrb(mpi_enreg%nproc,6,dtefield%fmkmem_max*nsppol*2))
  mpi_enreg%kptdstrb(:,:,:) = 0

  if (dtset%berryopt == 4) then
   allocate(mpi_enreg%kptdstrbi(mpi_enreg%nproc,6,dtefield%mkmem_max*nsppol*2))
   allocate(dtefield%cgqindex(3,6,nkpt*nsppol),dtefield%nneigh(nkpt))
   mpi_enreg%kptdstrbi(:,:,:) = 0
   dtefield%cgqindex(:,:,:) = 0 ; dtefield%nneigh(:) = 0
  end if

 end if

 pwind_alloc = mpw*dtefield%fmkmem_max
 allocate(pwind(pwind_alloc,2,3),pwnsfac(2,pwind_alloc))


!------------------------------------------------------------------------------
!---------------------- Compute efield_type variables -------------------------
!------------------------------------------------------------------------------

!Initialization of efield_type variables
 dtefield%efield_dot(:) = zero
 dtefield%dkvecs(:,:) = zero
 dtefield%maxnstr = 0    ; dtefield%maxnkstr  = 0
 dtefield%nstr(:) = 0    ; dtefield%nkstr(:) = 0
 allocate(dtefield%ikpt_dk(dtefield%fnkpt,2,3))
 allocate(dtefield%cgindex(nkpt,nsppol))
 allocate(dtefield%kgindex(nkpt))
 allocate(dtefield%fkgindex(dtefield%fnkpt))
 dtefield%ikpt_dk(:,:,:) = 0
 dtefield%cgindex(:,:) = 0
 dtefield%nband_occ = 0
 dtefield%kgindex(:) = 0
 dtefield%fkgindex(:) = 0

 if (dtset%berryopt == 4) then
  dtset%rfdir(1:3) = 1
 end if


!Compute spin degeneracy
 if (nsppol == 1) then
  dtefield%sdeg = two
 else if (nsppol == 2) then
  dtefield%sdeg = one
 end if

!Compute the number of occupied bands and check that
!it is the same for each k-point

 index = 0
 do isppol = 1, nsppol
  do ikpt = 1, nkpt

   nband_occ_k = 0
   nband_k = dtset%nband(ikpt + (isppol - 1)*nkpt)

   do iband = 1, nband_k
    index = index + 1
    if (abs(occ(index) - dtefield%sdeg) < tol8) nband_occ_k = nband_occ_k + 1
   end do

   if (dtset%berryopt == 4) then
    if (nband_k /= nband_occ_k) then
     write(message,'(a,a,a,a,a,a)')ch10,&
&     ' initberry: ERROR - ',ch10,&
&     '  In a finite electric field calculation, nband must be equal ',&
&     ch10,&
&     '  to the number of valence bands.'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')
    end if
   end if

   if ((ikpt > 1).or.(isppol > 1)) then
    if (dtefield%nband_occ /= nband_occ_k) then
     write(message,'(a,a,a,a)')ch10,&
&     ' initberry: ERROR - ',ch10,&
&     '   The number of valence bands is not the same for every k-point'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')
    end if
   else
    dtefield%nband_occ = nband_occ_k
   end if

  end do                ! close loop over ikpt
 end do                ! close loop over isppol


 if (dtset%berryopt == 4) then
  allocate(dtefield%smat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt*nsppol,2,3))

  dtefield%smat(:,:,:,:,:,:) = zero
 end if

 allocate(dtefield%sflag(dtefield%nband_occ,nkpt*nsppol,2,3))
 dtefield%sflag(:,:,:,:) = 0

!Compute the location of each wavefunction

 icg = 0
!ikg = 0
 do isppol = 1, nsppol
  do ikpt = 1, nkpt

   nband_k = dtset%nband(ikpt + (isppol-1)*nkpt)

   if (mpi_enreg%paral_compil_kpt == 1) then
    if (minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) - &
&    mpi_enreg%me)) /= 0) then
     cycle
    end if
   end if

   dtefield%cgindex(ikpt,isppol) = icg
   npw_k = npwarr(ikpt)
   icg = icg + dtset%nspinor*npw_k*nband_k

  end do
 end do

 ikg = 0
 do ikpt = 1, nkpt
  if (mpi_enreg%paral_compil_kpt == 1) then
   if ((minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,1) - &
&   mpi_enreg%me)) /= 0).and.&
&   (minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,nsppol) - &
&   mpi_enreg%me)) /= 0)) then
    cycle
   end if
  end if
  npw_k = npwarr(ikpt)
  dtefield%kgindex(ikpt) = ikg
  ikg = ikg + npw_k
 end do


!Compute the reciprocal lattice coordinates of the electric field
 if (dtset%berryopt == 4) then

  dtefield%efield_dot(1) = dot_product(dtset%efield(:),rprimd(:,1))
  dtefield%efield_dot(2) = dot_product(dtset%efield(:),rprimd(:,2))
  dtefield%efield_dot(3) = dot_product(dtset%efield(:),rprimd(:,3))

  write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
&  ' initberry: Reciprocal lattice coordinates of the electric field',ch10,&
&  '  efield_dot(1:3) = ',dtefield%efield_dot(1:3),ch10
  call wrtout(6,message,'COLL')

 end if

!------------------------------------------------------------------------------
!---------------------- Build the strings of k-points -------------------------
!------------------------------------------------------------------------------

 do idir = 1, 3

  if (dtset%rfdir(idir) == 1) then

!  Compute dk(:), the vector between a k-point and its nearest
!  neighbour along the direction idir

   dk(:) = zero
   dk(idir) = 1000._dp   ! initialize with a big number
   do ikpt = 2, dtefield%fnkpt
    diffk(:) = abs(dtefield%fkptns(:,ikpt) - dtefield%fkptns(:,1))
    if ((diffk(1) < dk(1)+tol8).and.(diffk(2) < dk(2)+tol8).and.&
&    (diffk(3) < dk(3)+tol8)) dk(:) = diffk(:)
   end do
   dtefield%dkvecs(:,idir) = dk(:)

!  For each k point, find k_prim such that k_prim= k + dk mod(G)
!  where G is a vector of the reciprocal lattice

   do ikpt = 1, dtefield%fnkpt

!   First: k + dk
    do ikpt1 = 1, dtefield%fnkpt
     diffk(:) = abs(dtefield%fkptns(:,ikpt1) - &
&     dtefield%fkptns(:,ikpt) - dk(:))
     if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
      dtefield%ikpt_dk(ikpt,1,idir) = ikpt1
      exit
     end if
    end do

!   Second: k - dk
    do ikpt1 = 1, dtefield%fnkpt
     diffk(:) = abs(dtefield%fkptns(:,ikpt1) - &
&     dtefield%fkptns(:,ikpt) + dk(:))
     if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
      dtefield%ikpt_dk(ikpt,2,idir) = ikpt1
      exit
     end if
    end do

   end do     ! ikpt

!  Find the string length, starting from k point 1
!  (all strings must have the same number of points)

   nkstr = 1
   ikpt1 = 1
   do ikpt = 1, dtefield%fnkpt
    ikpt1 = dtefield%ikpt_dk(ikpt1,1,idir)
    if (ikpt1 == 1) exit
    nkstr = nkstr + 1
   end do

!  Check that the string length is a divisor of nkpt
   if(mod(dtefield%fnkpt,nkstr) /= 0) then
    write(message,'(a,a,a,a,i5,a,i7)')ch10,&
&    ' berryphase: BUG -',ch10,&
&    '  The string length = ',nkstr,&
&    ', is not a divisor of fnkpt =',dtefield%fnkpt
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   dtefield%nkstr(idir) = nkstr
   dtefield%nstr(idir)  = dtefield%fnkpt/nkstr

  end if      ! dtset%rfdir(idir) == 1

  write(message,'(a,i1,a,i3,a,i3)')&
&  '  initberry: for direction ',idir,', nkstr = ',dtefield%nkstr(idir),&
&  ', nstr = ',dtefield%nstr(idir)
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

 end do     ! close loop over idir

 dtefield%maxnstr  = maxval(dtefield%nstr(:))
 dtefield%maxnkstr = maxval(dtefield%nkstr(:))
 allocate(dtefield%idxkstr(dtefield%maxnkstr,dtefield%maxnstr,3))
 dtefield%idxkstr(:,:,:) = 0


!Build the different strings

 allocate(kpt_mark(dtefield%fnkpt))
 do idir = 1, 3

  if (dtset%rfdir(idir) == 1) then

   iunmark = 1
   kpt_mark(:) = 0
   do istr = 1, dtefield%nstr(idir)

    do while(kpt_mark(iunmark) /= 0)
     iunmark = iunmark + 1
    end do
    dtefield%idxkstr(1,istr,idir) = iunmark
    kpt_mark(iunmark) = 1
    do ikstr = 2, dtefield%nkstr(idir)
     ikpt1 = dtefield%idxkstr(ikstr-1,istr,idir)
     ikpt2 = dtefield%ikpt_dk(ikpt1,1,idir)
     dtefield%idxkstr(ikstr,istr,idir) = ikpt2
     kpt_mark(ikpt2) = 1
    end do

   end do    ! istr

  end if         ! rfdir(idir) == 1

 end do           ! close loop over idir

 deallocate(kpt_mark)

!------------------------------------------------------------------------------
!------------ Build the array pwind that is needed to compute the -------------
!------------ overlap matrices at k +- dk                         -------------
!------------------------------------------------------------------------------

 ecut_eff = dtset%ecut*(dtset%dilatmx)**2
 exchn2n3d = 0 ; istwf_k = 1 ; ikg1 = 0
 pwind(:,:,:) = 0
 pwnsfac(1,:) = 1.0_dp
 pwnsfac(2,:) = 0.0_dp
 allocate(kg1_k(3,mpw))

 ipwnsfac = 0

 do idir = 1, 3

  if (dtset%rfdir(idir) == 1) then

   dk(:) = dtefield%dkvecs(:,idir)

   do ifor = 1, 2

    if (ifor == 2) dk(:) = -1._dp*dk(:)

!   Build pwind and kgindex
!   NOTE: The array kgindex is important for parallel execution.
!   In case nsppol = 2, it may happent that a particular processor
!   treats k-points at different spin polarizations.
!   In this case, it is not possible to address the elements of
!   pwind correctly without making use of the kgindex array.

    ikg = 0 ; ikpt_loc = 0 ; isppol = 1
    do ikpt = 1, dtefield%fnkpt

     ikpti = dtefield%indkk_f2ibz(ikpt,1)
     nband_k = dtset%nband(ikpti)
     ikpt1f = dtefield%ikpt_dk(ikpt,ifor,idir)
     ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)

     if (mpi_enreg%paral_compil_kpt == 1) then

      if ((minval(abs(mpi_enreg%proc_distrb(ikpti,1:nband_k,1) - &
&      mpi_enreg%me)) /= 0).and.&
&      (minval(abs(mpi_enreg%proc_distrb(ikpti,1:nband_k,nsppol) - &
&      mpi_enreg%me)) /= 0)) then
!      if (minval(abs(mpi_enreg%proc_distrb(ikpti,1:nband_k,1:nsppol) - &
!      &                mpi_enreg%me)) /= 0) then
       cycle
      end if

      ikpt_loc = ikpt_loc + 1

     end if

!    Build basis sphere of plane waves for the nearest neighbour of
!    the k-point (important for MPI //)

     kg1_k(:,:) = 0
     kpt1(:) = dtset%kptns(:,ikpt1i)
     call kpgsph(ecut_eff,exchn2n3d,gmet,ikg1,ikpt,istwf_k,kg1_k,kpt1,&
&     1,mpi_enreg,mpw,npw_k1)
     me_g0=mpi_enreg%me_g0


!    ji: fkgindex is defined here !
     dtefield%fkgindex(ikpt) = ikg

!    
!    Deal with symmetry transformations
!    

!    bra k-point k(b) and IBZ k-point kIBZ(b) related by
!    k(b) = alpha(b) S(b)^t kIBZ(b) + G(b)
!    where alpha(b), S(b) and G(b) are given by indkk_f2ibz
!    
!    For the ket k-point:
!    k(k) = alpha(k) S(k)^t kIBZ(k) + G(k) - GBZ(k)
!    where GBZ(k) takes k(k) to the BZ
!    

     isym  = dtefield%indkk_f2ibz(ikpt,2)
     isym1 = dtefield%indkk_f2ibz(ikpt1f,2)

!    Construct transformed G vector that enters the matching condition:
!    alpha(k) S(k)^{t,-1} ( -G(b) - GBZ(k) + G(k) )

     dg(:) = -dtefield%indkk_f2ibz(ikpt,3:5) &
&     -nint(-dtefield%fkptns(:,ikpt) - dk(:) - tol10 + &
&     dtefield%fkptns(:,ikpt1f)) &
&     +dtefield%indkk_f2ibz(ikpt1f,3:5)

     iadum(:)=0
     do idum=1,3
      iadum(:)=iadum(:)+ symrec(:,idum,isym1)*dg(idum)
     end do
     dg(:) = iadum(:)
     if ( dtefield%indkk_f2ibz(ikpt1f,6) == 1 ) dg(:) = -dg(:)

!    Construct S(k)^{t,-1} S(b)^{t}

     dum33(:,:) = zero
     do idumi=1,3
      do idumj=1,3
       dum33(:,idumi)=dum33(:,idumi)+ &
&       symrec(:,idumj,isym1)*dtset%symrel(idumi,idumj,isym)
      end do
     end do

!    Construct alpha(k) alpha(b)

     if (dtefield%indkk_f2ibz(ikpt,6) == dtefield%indkk_f2ibz(ikpt1f,6)) then
      itrs=0
     else
      itrs=1
     end if


     npw_k  = npwarr(ikpti)
!    npw_k1 = npwarr(ikpt1i)

!    loop over bra G vectors
     do ipw = 1, npw_k

!     NOTE: the bra G vector is taken for the sym-related IBZ k point,
!     not for the FBZ k point
      iadum(:) = kg(:,dtefield%kgindex(ikpti) + ipw)

!     Store non-symmorphic operation phase factor exp[i2\pi \alpha G \cdot t]

      if ( ipwnsfac == 0 ) then
       rdum=0.0_dp
       do idum=1,3
        rdum=rdum+dble(iadum(idum))*dtset%tnons(idum,isym)
       end do
       rdum=two_pi*rdum
       if ( dtefield%indkk_f2ibz(ikpt,6) == 1 ) rdum=-rdum
       pwnsfac(1,ikg+ipw) = cos(rdum)
       pwnsfac(2,ikg+ipw) = sin(rdum)
      end if

!     to determine r.l.v. matchings, we transformed the bra vector
!     Rotation
      iadum1(:)=0
      do idum1=1,3
       iadum1(:)=iadum1(:)+dum33(:,idum1)*iadum(idum1)
      end do
      iadum(:)=iadum1(:)
!     Time reversal
      if (itrs==1) iadum(:)=-iadum(:)
!     Translation
      iadum(:) = iadum(:) + dg(:)

      do jpw = 1, npw_k1
       iadum1(1:3) = kg1_k(1:3,jpw)
       if ( (iadum(1) == iadum1(1)).and. &
&       (iadum(2) == iadum1(2)).and. &
&       (iadum(3) == iadum1(3)) ) then
        pwind(ikg + ipw,ifor,idir) = jpw
!       write(6,'(a,2x,3i4,2x,i4)') 'Found !:',iadum1(:),jpw
        exit
       end if
      end do
     end do

     ikg  = ikg + npw_k

    end do    ! close loop over ikpt

    ipwnsfac = 1

   end do    ! close loop over ifor

  end if      ! rfdir(idir) == 1

 end do        ! close loop over idir

!Build mpi_enreg%kptdstrb
!array required to communicat the WFs between cpus in berryphase_new.f
!(MPI // over k-points)

 if (mpi_enreg%paral_compil_kpt == 1) then
  do idir = 1, 3
   if (dtset%rfdir(idir) == 1) then
    do ifor = 1, 2

     ikpt_loc = 0
     do isppol = 1, nsppol

      do ikpt = 1, dtefield%fnkpt

       ikpti = dtefield%indkk_f2ibz(ikpt,1)
       nband_k = dtset%nband(ikpti)
       ikpt1f = dtefield%ikpt_dk(ikpt,ifor,idir)
       ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)

       if (minval(abs(mpi_enreg%proc_distrb(ikpti,1:nband_k,isppol) - &
&       mpi_enreg%me)) /= 0) then
        cycle
       end if

       ikpt_loc = ikpt_loc + 1
       mpi_enreg%kptdstrb(mpi_enreg%me + 1,ifor+2*(idir-1),ikpt_loc) = &
&       ikpt1i + (isppol - 1)*nkpt

       mpi_enreg%kptdstrb(mpi_enreg%me+1,ifor+2*(idir-1),&
&       ikpt_loc+dtefield%fmkmem_max*nsppol) = &
&       ikpt1f + (isppol - 1)*dtefield%fnkpt

      end do   ! ikpt
     end do     ! isppol
    end do       ! ifor
   end if         ! dtset%rfdir(idir) == 1
  end do           ! idir
 end if             ! mpi_enreg%paral_compil_kpt == 1


!Build mpi_enreg%kptdstrbi
!(same as mpi_enreg%kptdstrb but for k-points in the iBZ),
!dtefield%cgqindex and dtefield%nneigh

 if ((dtset%berryopt == 4).and.(mpi_enreg%paral_compil_kpt == 1)) then

  ikpt_loc = 1
  do isppol = 1, nsppol
   do ikpt = 1, nkpt

    nband_k = dtset%nband(ikpt)
    ikptf = dtefield%i2fbz(ikpt)

    neigh(:) = 0 ; icg = 0 ; ikg = 0 ; flag_kpt = 0
    do idir = 1, 3

!    skip idir values for which efield_dot(idir) = 0
     if (abs(dtefield%efield_dot(idir)) < tol12) cycle

     do ifor = 1, 2

      ikpt1f = dtefield%ikpt_dk(ikptf,ifor,idir)
      ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)

      dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ikg
      ikg = ikg + npwarr(ikpt1i)

      flag = 0
      do ineigh = 1, (ifor+2*(idir-1))
       if (neigh(ineigh) == ikpt1i) then
        flag = 1
        dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ineigh
        dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
&        dtefield%cgqindex(2,ineigh,ikpt+(isppol-1)*nkpt)
        exit
       end if
      end do
      if (flag == 0) then
       neigh(ifor+2*(idir-1)) = ikpt1i
       dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
&       ifor+2*(idir-1)
       dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = icg
       if (isppol == 1) dtefield%nneigh(ikpt) = dtefield%nneigh(ikpt) + 1
       icg = icg + npwarr(ikpt1i)*dtset%nspinor*nband_k
      end if

      if (minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) - &
&      mpi_enreg%me)) /= 0) then
       cycle
      end if

      mpi_enreg%kptdstrbi(mpi_enreg%me + 1,ifor+2*(idir-1),&
&      ikpt_loc + dtefield%mkmem_max*nsppol) = &
&      ikpt1f + (isppol - 1)*dtefield%fnkpt

      flag_kpt = 1
      if (flag == 0) then
       mpi_enreg%kptdstrbi(mpi_enreg%me + 1,ifor+2*(idir-1),ikpt_loc) = &
&       ikpt1i + (isppol - 1)*nkpt
      end if

!     MVeithen: the if condition allows to avoid that the same wavefunction
!     is send several times to a particular cpu

     end do    ! ifor
    end do    ! idir

    if (flag_kpt == 1) ikpt_loc = ikpt_loc + 1

   end do    ! ikpt
  end do    ! isppol

 end if   ! berryopt and paral_compil_kpt

 if (mpi_enreg%paral_compil_kpt == 1) then

  count = mpi_enreg%nproc*6*dtefield%fmkmem_max*nsppol*2
  allocate(buffer(count))
  buffer(:) = reshape(mpi_enreg%kptdstrb(:,:,:),(/count/))
  call xsum_mpi_int(buffer,spaceComm,ierr)
  mpi_enreg%kptdstrb(:,:,:) = reshape(buffer(:),&
&  (/mpi_enreg%nproc,6,dtefield%fmkmem_max*nsppol*2/))
  deallocate(buffer)

  if (dtset%berryopt == 4) then
   count = mpi_enreg%nproc*6*dtefield%mkmem_max*nsppol*2
   allocate(buffer(count))
   buffer(:) = reshape(mpi_enreg%kptdstrbi(:,:,:),(/count/))
   call xsum_mpi_int(buffer,spaceComm,ierr)
   mpi_enreg%kptdstrbi(:,:,:) = reshape(buffer(:),&
&   (/mpi_enreg%nproc,6,dtefield%mkmem_max*nsppol*2/))
   deallocate(buffer)
  end if

 end if



!------------------------------------------------------------------------------
!------------------------ Estimate critical field -----------------------------
!------------------------------------------------------------------------------

!Compute the minimal value of the bandgap required to be below
!the critical field as defined by the relation
!| E_i*a_i | < E_g/n_i

 if (dtset%berryopt == 4) then

  do idir = 1, 3
   eg_dir(idir) = abs(dtefield%efield_dot(idir))*dtefield%nkstr(idir)
  end do
  eg = maxval(eg_dir)
  eg_ev = eg*Ha_eV

  write(message,'(a,a,a,a,a,a,a,a,f7.2,a,a)')ch10,&
&  ' initberry: COMMENT - ',ch10,&
&  '  As a rough estimate,',ch10,&
&  '  to be below the critical field, the bandgap of your system',ch10,&
&  '  should be larger than ',eg_ev,' eV.',ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(6,message,'COLL')

 end if

 deallocate(kg1_k)
 deallocate(dum_fwtk) !! by MM

end subroutine initberry
!!***
