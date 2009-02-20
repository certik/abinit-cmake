!{\src2tex{textfont=tt}}
!!****f* ABINIT/mv_3dte
!! NAME
!! mv_3dte
!!
!! FUNCTION
!! Compute the finite difference expression of the k-point derivative
!! using the PEAD formulation of the third-order energy
!! (see Nunes and Gonze PRB 63, 155107 (2001) Eq. 102)
!! and the finite difference formula of Marzari and Vanderbilt
!! (see Marzari and Vanderbilt, PRB 56, 12847 (1997), Appendix B)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol) = array for planewave
!!                                         coefficients of wavefunctions
!!  cgindex(nkpt2,nsppol) = for each k-point, cgindex stores the location
!!                          of the WF in the cg array
!!  cg1 = first-order wavefunction relative to the perturbations i1pert
!!  cg3 = first-order wavefunction relative to the perturbations i3pert
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  i1pert,i2pert,i3pert = type of perturbation that has to be computed
!!  i1dir,i3dir=directions of the corresponding perturbations
!!  kneigh(30,nkpt2) = index of the neighbours of each k-point
!!  kptindex(2,nkpt3)= index of the k-points in the reduced BZ
!!                     related to a k-point in the full BZ
!!  kpt3(3,nkpt3) = reduced coordinates of k-points in the full BZ
!!  mband = maximum number of bands
!!  mkmem = maximum number of k points which can fit in core memory
!!  mkmem_max = maximal number of k-points on each processor (MPI //)
!!  mk1mem = maximum number of k points for first-order WF
!!           which can fit in core memory
!!  mpert =maximum number of ipert
!!  mpi_enreg=MPI-parallelisation information
!!  mpw   = maximum number of planewaves in basis sphere (large number)
!!  mvwtk(30,nkpt) = weights to compute the finite difference ddk
!!  natom = number of atoms in unit cell
!!  nkpt2 = number of k-points in the reduced part of the BZ
!!          nkpt2 = nkpt/2 in case of time-reversal symmetry (kptopt = 2)
!!  nkpt3 = number of k-points in the full BZ
!!  nneigh = total number of neighbours required to evaluate the finite
!!           difference formula
!!  npwarr(nkpt) = array holding npw for each k point
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  pwind(mpw,nneigh,mkmem) = array used to compute the overlap matrix smat
!!                           between k-points
!!  rmet(3,3)=real space metric (bohr**2)
!!  tmpfil(15)=names for the temporary files based on dtfil%filnam_ds(5)
!!
!! OUTPUT
!!  d3_berry(2,3) = Berry-phase part of the third-order energy
!!
!! SIDE EFFECTS
!!  mpi_enreg=MPI-parallelisation information
!!
!! NOTES
!! For a given set of values of i1pert,i2pert,i3pert,i1dir and
!! i3dir, the routine computes the k-point derivatives for
!! 12dir = 1,2,3
!!
!! TODO
!!
!! PARENTS
!!      loop3dte
!!
!! CHILDREN
!!      mpi_recv,mpi_send,status,wrtout,xcomm_world,xsum_mpi_dp2d,dzgedi,dzgefa
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mv_3dte(cg,cgindex,cg1,cg3,dtset,dtfil,d3_berry,gmet,gprimd,&
&                   hdr,i1pert,i2pert,i3pert,i1dir,i3dir,kneigh,kptindex,&
&                   kpt3,mband,mkmem,mkmem_max,mk1mem,mpert,mpi_enreg,&
&                   mpw,mvwtk,natom,nkpt2,nkpt3,nneigh,npwarr,nspinor,&
&                   nsppol,occ,pwind,psps,rmet,tmpfil)

 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_lib00numeric
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!
!---  Arguments : integer scalars
 integer, intent(in) :: i1dir,i1pert,i2pert,i3dir,i3pert,mband,mk1mem
 integer, intent(in) :: mkmem,mkmem_max,mpert,mpw,natom
 integer, intent(in) :: nkpt2,nkpt3,nneigh,nspinor,nsppol
!
!---  Arguments : integer arrays
 integer, intent(in) :: cgindex(nkpt2,nsppol)
 integer, intent(in) :: kneigh(30,nkpt2),kptindex(2,nkpt3),npwarr(nkpt2)
 integer, intent(in) :: pwind(mpw,nneigh,mkmem)
!
!---  Arguments : real(dp) scalars
!
!---  Arguments : real(dp) arrays
 real(dp), intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp), intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp), intent(in) :: cg3(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp), intent(in) :: gmet(3,3),gprimd(3,3),kpt3(3,nkpt3)
 real(dp), intent(in) :: mvwtk(30,nkpt2),occ(mband*nkpt2*nsppol),rmet(3,3)
 real(dp), intent(out) :: d3_berry(2,3)
!
!---  Arguments : character variables
 character(len=fnlen), intent(in) :: tmpfil(15)
!
!---  Arguments : structured datatypes
 type(MPI_type), intent(inout) :: mpi_enreg
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(in) :: dtset
 type(hdr_type), intent(in) :: hdr
 type(pseudopotential_type), intent(in) :: psps

!Local variables-------------------------------
!
!---- Local variables : integer scalars
 integer :: bantot,count,counter,count1,flag,iband,icg,icg1
 integer :: ierr,iexit,ii,ikpt,ikpt_loc,ikpt2
 integer :: ikpt_rbz,ineigh,info,ipw,isppol,jband,jcg,jj,jkpt,job,jpw
 integer :: lband,ll,lpband,nband_k,npw_k,npw_k1,source,dest,tag
 integer :: spaceComm
 integer,parameter :: level=22
!
!---- Local variables : integer arrays
 integer,allocatable :: ipvt(:)
!
!---- Local variables : real(dp) scalars
 real(dp) :: alpha,c1,c2,c3,c4,c5,c6,dotnegi,dotnegr,dotposi,dotposr,dtm,mod_
!
!---- Local variables : real(dp) arrays
 real(dp) :: d3_aux(2,3),det(2,2),dk(3),dk_(3),dkcart(3),dummymat(3,3)
 real(dp) :: mod(0:10),z1(2),z2(2)
 real(dp),allocatable :: buffer(:,:),cgq(:,:),cg1q(:,:),cg3q(:,:),kpt(:,:)
 real(dp),allocatable :: qmat(:,:,:),s13mat(:,:,:),s1mat(:,:,:),s3mat(:,:,:)
 real(dp),allocatable :: smat(:,:,:),zgwork(:,:)
!
!---- Local variables : character variables
 character(len=500) :: message
!
!---- Local variables : structured datatypes
 type(dens_sym_operator_type) :: densymop_gs

#if defined MPI
             integer :: status1(MPI_STATUS_SIZE)
!BEGIN TF_CHANGES
             call xcomm_world(mpi_enreg,spaceComm)
!END TF_CHANGES
#endif


! ***********************************************************************

!DEBUG
!write(6,*)' mv_3dte : enter '
!stop
!ENDDEBUG

 call status(0,dtfil%filstat,iexit,level,'enter         ')

 write(message,'(8a)') ch10,&
& ' mv_3dte : finite difference expression of the k-point derivative',ch10,&
& '           is performed using the PEAD formulation of ',&
& 'the third-order energy',ch10,&
& '           (see Nunes and Gonze PRB 63, 155107 (2001) Eq. 102)',ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')


 isppol = 1
 bantot = 0
 ikpt_loc = 0
 d3_aux(:,:) = 0_dp

 allocate(s13mat(2,mband,mband),smat(2,mband,mband),s1mat(2,mband,mband))
 allocate(qmat(2,mband,mband),ipvt(mband),s3mat(2,mband,mband))
 allocate(zgwork(2,mband))


! Loop over k-points
! COMMENT: Every processor has to make mkmem_max iterations
!          even if mkmem < mkemem_max. This is due to the fact
!          that it still has to communicate its wavefunctions
!          to other processors even if it has no more overlap
!          matrices to compute.

 ikpt_loc = 0 ; ikpt = 0
 do while (ikpt_loc < mkmem_max)

  if (ikpt_loc < mkmem) ikpt = ikpt + 1

  if (mpi_enreg%paral_compil_kpt == 1) then
   if ((minval(abs(mpi_enreg%proc_distrb(ikpt,1:mband,1:dtset%nsppol) &
&      - mpi_enreg%me)) /= 0).and.(ikpt_loc < mkmem)) cycle
  end if

  ikpt_loc = ikpt_loc + 1
  npw_k = npwarr(ikpt)
  counter = 100*ikpt
  call status(counter,dtfil%filstat,iexit,level,'loop over k   ')

  ii = cgindex(ikpt,1)

! Loop on the  neighbours

  do ineigh = 1,nneigh

   s13mat(:,:,:) = zero
   smat(:,:,:) = zero
   s1mat(:,:,:) = zero
   s3mat(:,:,:) = zero
   qmat(:,:,:) = zero

   ikpt2  = kneigh(ineigh,ikpt)
   ikpt_rbz = kptindex(1,ikpt2)   ! index of the k-point in the reduced BZ
   jj = cgindex(ikpt_rbz,1)
   nband_k = dtset%nband(ikpt_rbz)
   npw_k1 = npwarr(ikpt_rbz)
   dk_(:) = kpt3(:,ikpt2) - dtset%kptns(:,ikpt)
   dk(:)  = dk_(:) - nint(dk_(:))

   count = nspinor*mband*npw_k1
   allocate(cgq(2,count))
   allocate(cg1q(2,count))
   allocate(cg3q(2,count))

#if defined MPI

          source = mpi_enreg%proc_distrb(ikpt_rbz,1,isppol)

          do jkpt = 1, nkpt2

           if ((jkpt == ikpt_rbz).and.(source /= mpi_enreg%me).and.&
&              (ikpt_loc <= mkmem)) then

            tag = jkpt

            allocate(buffer(2,3*count))
            call MPI_RECV(buffer,2*3*count,MPI_DOUBLE_PRECISION,&
                          source,tag,spaceComm,status1,ierr)

            cgq(:,1:count)  = buffer(:,1:count)
            cg1q(:,1:count) = buffer(:,count+1:2*count)
            cg3q(:,1:count) = buffer(:,2*count+1:3*count)
            deallocate(buffer)

           end if

!----------------------------------------------------------------------------
!--------------- Here: send the WF to all the cpus that need it -------------
!----------------------------------------------------------------------------

           do dest = 1, mpi_enreg%nproc

            if ((minval(abs(mpi_enreg%proc_distrb(jkpt,1:mband,1:dtset%nsppol) &
&               - mpi_enreg%me)) == 0).and.&
&               (mpi_enreg%kptdstrb(dest,ineigh,ikpt_loc) == jkpt)) then

             jcg = cgindex(jkpt,1)

             if (((dest-1) == mpi_enreg%me)) then

              cgq(:,1:count)  = cg(:,jcg+1:jcg+count)
              cg1q(:,1:count) = cg1(:,jcg+1:jcg+count)
              cg3q(:,1:count) = cg3(:,jcg+1:jcg+count)

             else

              tag = jkpt
              count1 = npwarr(jkpt)*mband*nspinor
              allocate(buffer(2,3*count1))
              buffer(:,1:count1)            = cg(:,jcg+1:jcg+count1)
              buffer(:,count1+1:2*count1)   = cg1(:,jcg+1:jcg+count1)
              buffer(:,2*count1+1:3*count1) = cg3(:,jcg+1:jcg+count1)

              call MPI_SEND(buffer,2*3*count1,MPI_DOUBLE_PRECISION,&
                           (dest-1),tag,spaceComm,status1,ierr)

              deallocate(buffer)

             end if

            end if

           end do          ! loop over dest

          end do          ! loop over jkpt

          if (ikpt_loc > mkmem) then
           deallocate(cgq,cg1q,cg3q)
           cycle
          end if

#else
!  no // over k-points

   cgq(:,1:count)  = cg(:,jj+1:jj+count)
   cg1q(:,1:count) = cg1(:,jj+1:jj+count)
   cg3q(:,1:count) = cg3(:,jj+1:jj+count)

#endif

! Compute overlap matrices

   if (kptindex(2,ikpt2) == 0) then  ! no time-reversal symmetry

    do ipw = 1, npw_k

     jpw = pwind(ipw,ineigh,ikpt_loc)
     if (jpw /= 0) then

      do iband = 1, nband_k
       do jband = 1, nband_k

        icg = ii + (iband-1)*npw_k + ipw
        jcg = (jband-1)*npw_k1 + jpw

        smat(1,iband,jband) = smat(1,iband,jband) + &
&         cg(1,icg)*cgq(1,jcg) + cg(2,icg)*cgq(2,jcg)
        smat(2,iband,jband) = smat(2,iband,jband) + &
&         cg(1,icg)*cgq(2,jcg) - cg(2,icg)*cgq(1,jcg)

        s13mat(1,iband,jband) = s13mat(1,iband,jband) + &
&         cg1(1,icg)*cg3q(1,jcg) + cg1(2,icg)*cg3q(2,jcg)
        s13mat(2,iband,jband) = s13mat(2,iband,jband) + &
&         cg1(1,icg)*cg3q(2,jcg) - cg1(2,icg)*cg3q(1,jcg)

        s1mat(1,iband,jband) = s1mat(1,iband,jband) + &
&         cg1(1,icg)*cgq(1,jcg) + cg1(2,icg)*cgq(2,jcg) + &
&         cg(1,icg)*cg1q(1,jcg) + cg(2,icg)*cg1q(2,jcg)
        s1mat(2,iband,jband) = s1mat(2,iband,jband) + &
&         cg1(1,icg)*cgq(2,jcg) - cg1(2,icg)*cgq(1,jcg) + &
&         cg(1,icg)*cg1q(2,jcg) - cg(2,icg)*cg1q(1,jcg)

        s3mat(1,iband,jband) = s3mat(1,iband,jband) + &
&         cg3(1,icg)*cgq(1,jcg) + cg3(2,icg)*cgq(2,jcg) + &
&         cg(1,icg)*cg3q(1,jcg) + cg(2,icg)*cg3q(2,jcg)
        s3mat(2,iband,jband) = s3mat(2,iband,jband) + &
&         cg3(1,icg)*cgq(2,jcg) - cg3(2,icg)*cgq(1,jcg) + &
&         cg(1,icg)*cg3q(2,jcg) - cg(2,icg)*cg3q(1,jcg)

       end do
      end do

     end if

    end do   ! ipw

   else                              ! use time-reversal symmetry

    do ipw = 1,npw_k

     jpw = pwind(ipw,ineigh,ikpt_loc)
     if (jpw /= 0) then

      do iband = 1, nband_k
       do jband = 1, nband_k

        icg = ii + (iband-1)*npw_k + ipw
        jcg = (jband-1)*npw_k1 + jpw

        smat(1,iband,jband) = smat(1,iband,jband) + &
&         cg(1,icg)*cgq(1,jcg) - cg(2,icg)*cgq(2,jcg)
        smat(2,iband,jband) = smat(2,iband,jband) - &
&         cg(1,icg)*cgq(2,jcg) - cg(2,icg)*cgq(1,jcg)

        s13mat(1,iband,jband) = s13mat(1,iband,jband) + &
&         cg1(1,icg)*cg3q(1,jcg) - cg1(2,icg)*cg3q(2,jcg)
        s13mat(2,iband,jband) = s13mat(2,iband,jband) - &
&         cg1(1,icg)*cg3q(2,jcg) - cg1(2,icg)*cg3q(1,jcg)

        s1mat(1,iband,jband) = s1mat(1,iband,jband) + &
&         cg1(1,icg)*cgq(1,jcg) - cg1(2,icg)*cgq(2,jcg) + &
&         cg(1,icg)*cg1q(1,jcg) - cg(2,icg)*cg1q(2,jcg)
        s1mat(2,iband,jband) = s1mat(2,iband,jband) - &
&         cg1(1,icg)*cgq(2,jcg) - cg1(2,icg)*cgq(1,jcg) - &
&         cg(1,icg)*cg1q(2,jcg) - cg(2,icg)*cg1q(1,jcg)

        s3mat(1,iband,jband) = s3mat(1,iband,jband) + &
&         cg3(1,icg)*cgq(1,jcg) - cg3(2,icg)*cgq(2,jcg) + &
&         cg(1,icg)*cg3q(1,jcg) - cg(2,icg)*cg3q(2,jcg)
        s3mat(2,iband,jband) = s3mat(2,iband,jband) - &
&         cg3(1,icg)*cgq(2,jcg) - cg3(2,icg)*cgq(1,jcg) - &
&         cg(1,icg)*cg3q(2,jcg) - cg(2,icg)*cg3q(1,jcg)

       end do
      end do

     end if

    end do   ! ipw

   end if

   deallocate(cgq,cg1q,cg3q)

!  Compute qmat, the inverse of smat

   job = 1  ! compute inverse only
   qmat(:,:,:) = smat(:,:,:)

   call dzgefa(qmat,dtset%nband(ikpt),dtset%nband(ikpt),ipvt,info)
   call dzgedi(qmat,dtset%nband(ikpt),dtset%nband(ikpt),ipvt,det,zgwork,job)

!DEBUG
!write(100,*)
!write(100,*)'ikpt = ',ikpt,'ineigh = ',ineigh
!do iband = 1,dtset%nband(ikpt)
!do jband = 1,dtset%nband(ikpt)
! c1 = 0_dp ; c2 = 0_dp
! do lband = 1,dtset%nband(ikpt)
!  c1 = c1 + smat(1,iband,lband)*qmat(1,lband,jband) - &
!&           smat(2,iband,lband)*qmat(2,lband,jband)
!  c2 = c2 + smat(1,iband,lband)*qmat(2,lband,jband) + &
!&           smat(2,iband,lband)*qmat(1,lband,jband)
! end do
! write(100,'(2(2x,i2),2(2x,f16.9))')iband,jband,&
!& c1,c2
!end do
!end do
!ENDDEBUG



!  Accumulate sum over bands

   dotposr = 0_dp ; dotposi = 0_dp
   dotnegr = 0_dp ; dotnegi = 0_dp
   do iband = 1, dtset%nband(ikpt)
    do jband = 1, dtset%nband(ikpt)

     dotposr = dotposr + &
&      s13mat(1,iband,jband)*qmat(1,jband,iband) - &
&      s13mat(2,iband,jband)*qmat(2,jband,iband)
     dotposi = dotposi + &
&      s13mat(1,iband,jband)*qmat(2,jband,iband) + &
&      s13mat(2,iband,jband)*qmat(1,jband,iband)


     do lband = 1, dtset%nband(ikpt)
      do lpband= 1, dtset%nband(ikpt)

       z1(1) = s1mat(1,iband,jband)*qmat(1,jband,lband) - &
&        s1mat(2,iband,jband)*qmat(2,jband,lband)
       z1(2) = s1mat(1,iband,jband)*qmat(2,jband,lband) + &
&        s1mat(2,iband,jband)*qmat(1,jband,lband)

       z2(1) = s3mat(1,lband,lpband)*qmat(1,lpband,iband) - &
&        s3mat(2,lband,lpband)*qmat(2,lpband,iband)
       z2(2) = s3mat(1,lband,lpband)*qmat(2,lpband,iband) + &
&        s3mat(2,lband,lpband)*qmat(1,lpband,iband)

       dotnegr = dotnegr + &
&        z1(1)*z2(1) - z1(2)*z2(2)
       dotnegi = dotnegi + &
&        z1(1)*z2(2) + z1(2)*z2(1)

      end do   ! lpband
     end do   ! lband

    end do   ! jband
   end do   ! iband

   d3_aux(1,:) = d3_aux(1,:) + &
&    dk(:)*mvwtk(ineigh,ikpt)*dtset%wtk(ikpt)*(2_dp*dotposr-dotnegr)
   d3_aux(2,:) = d3_aux(2,:) + &
&    dk(:)*mvwtk(ineigh,ikpt)*dtset%wtk(ikpt)*(2_dp*dotposi-dotnegi)

  end do        ! End loop over neighbours


  bantot = bantot + dtset%nband(ikpt)

 end do      ! End loop over k-points

!DEBUG
!write(*,*)'after loop over k-points'
!write(*,*)'me = ',mpi_enreg%me,'d3_aux = '
!write(*,*)d3_aux(:,1)
!write(*,*)d3_aux(:,2)
!write(*,*)d3_aux(:,3)
!call leave_test
!ENDDEBUG

 if (mpi_enreg%paral_compil_kpt == 1) then
  call xsum_mpi_dp2d(d3_aux,spaceComm,ierr)
 end if

!DEBUG
!write(*,*)'after xsum_mpi'
!call leave_test
!stop
!ENDDEBUG

 deallocate(s13mat,smat,s1mat)
 deallocate(qmat,ipvt,s3mat)
 deallocate(zgwork)

!  Take minus the imaginary part

 d3_berry(1,:) = -1_dp*d3_aux(2,:)
 d3_berry(2,:) = d3_aux(1,:)


 d3_berry(2,:) = 0_dp

!DEBUG
!write(100,*)'mv_3dte.f : d3_berry'
!write(100,*)'Perturbation',i1dir,i3dir
!write(100,*)
!write(100,*)'before transformation'
!write(100,*)'real part'
!write(100,'(3(2x,f20.9))')d3_berry(1,:)
!write(100,*)
!write(100,*)'imaginary part'
!write(100,'(3(2x,f20.9))')d3_berry(2,:)
!write(100,*)
!write(100,*)'after transformation'
!ENDDEBUG

!  Compute the projection on the basis vectors of
!  reciprocal space

 d3_aux(:,:) = 0_dp
 do ii = 1,3
  do jj = 1,3
   d3_aux(:,ii) = d3_aux(:,ii) + gmet(ii,jj)*d3_berry(:,jj)
  end do
 end do
 d3_berry(:,:) = d3_aux(:,:)

!Write out the berryphase part of the third order energy

 if (mpi_enreg%me == 0) then

  write(message,'(a,a,a)')ch10,&
&  ' Berryphase part of the third-order energy:',ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')

  if (i1pert < natom + 1) then
   write(message,'(a,i3,a,i3)')&
&   '            j1: Displacement of atom ',i1pert,&
&   ' along direction ',i1dir
  else if (i1pert == natom + 2) then
   write(message,'(a,i3)')&
&     '            j1: homogenous electric field along direction ',i1dir
  end if
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')

  write(message,'(a)')&
&  '            j2: k-point derivative along direction i2dir '
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')

  if (i3pert < natom + 1) then
   write(message,'(a,i3,a,i3,a)')&
&   '            j3: Displacement of atom ',i3pert,&
&   ' along direction ',i3dir,ch10
  else if (i3pert == natom + 2) then
   write(message,'(a,i3,a)')&
&     '            j3: homogenous electric field along direction ',i3dir,ch10
  end if
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')

  write(ab_out,'(5x,a5,8x,a9,5x,a14)')'i2dir','real part','imaginary part'
  write(6,'(5x,a5,8x,a9,5x,a14)')'i2dir','real part','imaginary part'
  do ii = 1,3
   write(ab_out,'(7x,i1,3x,f16.9,3x,f16.9)')ii,&
&   d3_berry(1,ii),d3_berry(2,ii)
   write(6,'(7x,i1,3x,f16.9,3x,f16.9)')ii,&
&   d3_berry(1,ii),d3_berry(2,ii)
  end do

 end if    ! mpi_enreg%me == 0

!DEBUG
!write(100,*)'real part'
!write(100,'(3(2x,f20.9))')d3_berry(1,:)
!write(100,*)
!write(100,*)'imaginary part'
!write(100,'(3(2x,f20.9))')d3_berry(2,:)
!ENDDEBUG



!close(dtfil%unwff1)
!close(dtfil%unwff2)

 call status(0,dtfil%filstat,iexit,level,'exit          ')

!DEBUG
!write(6,*)' mv_3dte : exit '
!ENDDEBUG

end subroutine mv_3dte
!!***
