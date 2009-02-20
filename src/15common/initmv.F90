!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmv
!! NAME
!! initmv
!!
!! FUNCTION
!! Initialize finite difference calculation of the ddk im mv_3dte.f
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (MVeithen)
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
!!  kneigh(30,nkpt2) = index of the neighbours of each k-point
!!  kptindex(2,nkpt3)= index of the k-points in the reduced BZ
!!                     related to a k-point in the full BZ
!!  kpt3(3,nkpt3) = reduced coordinates of k-points in the full BZ
!!  mband = maximum number of bands
!!  mkmem = number of k points which can fit in memory
!!  mkmem_max = maximal number of k-points on each processor (MPI //)
!!  mpi_enreg = informations about MPI parallelization
!!  mpw = maximum number of plane waves
!!  nband(nkpt*nsppol)=number of bands at each k point, for each polarization
!!  nkpt2 = number of k-points in the reduced BZ
!!  nkpt3 = number of k-points in the full BZ
!!  nneigh = total number of neighbours required to evaluate the finite
!!          difference formula
!!  npwarr(nkpt2)=number of planewaves at each k point
!!  nsppol = number of spin polarizations
!!  occ(mband*nkpt*nsppol) = occupation number for each band for each k
!!
!! OUTPUT
!! cgindex(nkpt2,nsppol) = for each k-point, cgindex tores the location
!!                         of the WF in the cg array
!! mpi_enreg%kptdstrb(nproc,nneigh,mkmem_max) = Array required for the MPI
!!         parallelization of the mv_3dte.f routine
!!        kptdstrb(mpi_enreg%me + 1,ineigh,ikpt_loc) = ikpt_rbz
!!        me = index of the current processor
!!        ineigh = index of a neighbour
!!        ikpt_loc = index of the iteration on ikpt on the current processor
!!        ikpt_rbz = index of a k-point in the reduced BZ
!! pwind(mpw,nneigh,mkmem) = array used to compute the overlap matrix smat
!!                           between k-points
!!                           (see initberry.f for more explanations)
!!
!! COMMENTS
!!
!! PARENTS
!!      nonlinear
!!
!! CHILDREN
!!      kpgio,leave_new,wrtout,xcomm_world,xsum_mpi_int
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine initmv(cgindex,dtfil,dtset,gmet,kg,kneigh,kptindex,&
&  kpt3,mband,mkmem,mkmem_max,mpi_enreg,mpw,nband,nkpt2,&
&  nkpt3,nneigh,npwarr,nsppol,occ,pwind)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13recipspace
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem,mkmem_max,mpw,nkpt2,nkpt3,nneigh,nsppol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),kneigh(30,nkpt2),kptindex(2,nkpt3)
 integer,intent(in) :: nband(nkpt2),npwarr(nkpt2)
 integer,intent(out) :: cgindex(nkpt2,nsppol),pwind(mpw,nneigh,mkmem)
 real(dp),intent(in) :: gmet(3,3),kpt3(3,nkpt3),occ(mband*nkpt2*nsppol)

!Local variables-------------------------------
!scalars
 integer :: count,flag,iband,icg,ierr,ii,ikg,ikg1,ikpt,ikpt2,ikpt_loc,ikpt_rbz
 integer :: index,ineigh,ipw,isppol,jpw,nband_k,nband_occ,nband_occ_k,npw_k
 integer :: npw_k1,orig,spaceComm
 real(dp) :: ecut_eff,sdeg
 character(len=500) :: message
 character(len=fnlen) :: tmpfil
!arrays
 integer :: dg(3)
 integer,allocatable :: buffer(:),kg1(:,:),kg1_k(:,:),npwar1(:),npwtot(:)
 real(dp) :: dk(3),dk_(3)
 real(dp),allocatable :: kpt1(:,:)

!************************************************************************

 if (mpi_enreg%paral_compil_kpt == 1) then
! BEGIN TF_CHANGES
  call xcomm_world(mpi_enreg,spaceComm)
! END TF_CHANGES
  mpi_enreg%kptdstrb(:,:,:) = 0
 end if

 ecut_eff = dtset%ecut*(dtset%dilatmx)**2
 allocate(kg1_k(3,mpw))
 allocate(kg1(3,mkmem*mpw),kpt1(3,nkpt2),npwar1(nkpt2),npwtot(nkpt2))
 kg1_k(:,:) = 0
 pwind(:,:,:) = 0
 cgindex(:,:) = 0
 tmpfil = trim(dtfil%filnam_ds(5))//'_KG1'

!Compute the number of occupied bands.
!Check that it is the same for every k-point and that
!nband(ikpt) is equal to this value

 if (nsppol == 1) then
  sdeg = two
 else if (nsppol == 2) then
  sdeg = one
 end if

 index = 0
 do isppol = 1, nsppol
  do ikpt = 1, nkpt2

   nband_occ_k = 0
   nband_k = nband(ikpt + (isppol - 1)*nkpt2)

   do iband = 1, nband_k
    index = index + 1
    if (abs(occ(index) - sdeg) < tol8) nband_occ_k = nband_occ_k + 1
   end do

   if (nband_k /= nband_occ_k) then
    write(message,'(a,a,a,a,a,a)')ch10,&
&    ' initmv: ERROR - ',ch10,&
&    '  In a non-linear response calculation, nband must be equal ',ch10,&
&    '  to the number of valence bands.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

   if ((ikpt > 1).or.(isppol > 1)) then
    if (nband_occ /= nband_occ_k) then
     write(message,'(a,a,a,a)')ch10,&
&     ' initmv: ERROR - ',ch10,&
&     '   The number of valence bands is not the same for every k-point'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')
    end if
   else
    nband_occ = nband_occ_k
   end if

  end do                ! close loop over ikpt
 end do                ! close loop over isppol

!Find the location of each wavefunction

 icg = 0
 do isppol = 1, nsppol
  do ikpt = 1, nkpt2

   nband_k = dtset%nband(ikpt)
   npw_k = npwarr(ikpt)

   if (mpi_enreg%paral_compil_kpt == 1) then
    if (minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) - &
&    mpi_enreg%me)) /= 0) then
     cycle
    end if
   end if

   cgindex(ikpt,isppol) = icg
   icg = icg + dtset%nspinor*npw_k*nband_k

  end do
 end do


!Build pwind and kptdstrb

 do ineigh = 1, nneigh

  do ikpt = 1, nkpt2
   ikpt2  = kneigh(ineigh,ikpt)
   ikpt_rbz = kptindex(1,ikpt2)   ! index of the k-point in the reduced BZ
   kpt1(:,ikpt) = dtset%kptns(:,ikpt_rbz)
  end do

! Set up the basis sphere of plane waves at kpt1
  kg1(:,:) = 0
  call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg1,tmpfil,&
&  kpt1,mkmem,dtset%nband,nkpt2,'PERS',mpi_enreg,mpw,&
&  npwar1,npwtot,dtset%nsppol,dtfil%unkg1)

  ikg = 0 ; ikg1 = 0 ; ikpt_loc = 0

  do ikpt = 1, nkpt2

   nband_k = dtset%nband(ikpt)
   ikpt2  = kneigh(ineigh,ikpt)
   ikpt_rbz = kptindex(1,ikpt2)   ! index of the k-point in the reduced BZ

   if (mpi_enreg%paral_compil_kpt == 1) then
    if (minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,1:dtset%nsppol) &
&    - mpi_enreg%me)) /= 0) cycle
   end if

   ikpt_loc = ikpt_loc + 1
   if (mpi_enreg%paral_compil_kpt == 1) then
    mpi_enreg%kptdstrb(mpi_enreg%me + 1,ineigh,ikpt_loc) = ikpt_rbz
   end if

   flag = 0
   npw_k = npwarr(ikpt)
   npw_k1 = npwarr(ikpt_rbz)
   dk_(:) = kpt3(:,ikpt2) - dtset%kptns(:,ikpt)
   dk(:)  = dk_(:) - nint(dk_(:))
   dg(:)  = nint(dk(:) - dk_(:))


   if (kptindex(2,ikpt2) == 0) then
    kg1_k(:,1:npw_k1) = kg1(:,ikg1+1:ikg1+npw_k1)
    if (dg(1)==0.and.dg(2)==0.and.dg(3)==0) flag = 1
   else
    kg1_k(:,1:npw_k1) = -1*kg1(:,ikg1+1:ikg1+npw_k1)
   end if

   orig = 1
   do ipw = 1, npw_k
    do jpw = orig, npw_k1

     if ((kg(1,ikg + ipw) == kg1_k(1,jpw) - dg(1)).and. &
&     (kg(2,ikg + ipw) == kg1_k(2,jpw) - dg(2)).and. &
&     (kg(3,ikg + ipw) == kg1_k(3,jpw) - dg(3)))  then

      pwind(ipw,ineigh,ikpt_loc) = jpw
      if (flag == 1)  orig = jpw + 1
      exit

     end if

    end do
   end do

   ikg = ikg + npw_k
   ikg1 = ikg1 + npw_k1

  end do     ! close loop over k-points

 end do    ! close loop over ineigh


 if (mpi_enreg%paral_compil_kpt == 1) then
  count = mpi_enreg%nproc*nneigh*mkmem_max
  allocate(buffer(count))
  buffer(:) = reshape(mpi_enreg%kptdstrb(:,:,:),(/count/))
  call xsum_mpi_int(buffer,spaceComm,ierr)
  mpi_enreg%kptdstrb(:,:,:) = reshape(buffer(:),&
&  (/mpi_enreg%nproc,nneigh,mkmem_max/))
  deallocate(buffer)
 end if


!----------------------------------------------------------------------------

 deallocate(kg1,kg1_k,kpt1,npwar1,npwtot)

end subroutine initmv
!!***
