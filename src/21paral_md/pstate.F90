!{\src2tex{textfont=tt}}
!!****f* ABINIT/pstate
!! NAME
!! pstate
!!
!! FUNCTION
!! Primary routine for conducting parareel/DFT calculations by CG minimization.
!! -- More or less a wrapper...
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, JYR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  acell
!!  codvsn=code version
!!  cpui=initial CPU time
!!  mband=maximum number of bands
!!  mgfft=maximum single fft dimension
!!  mkmem=maximum number of k points which can fit in core memory
!!  mpw=maximum number of planewaves in basis sphere (large number)
!!  natom=number of atoms in unit cell
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nkpt=number of k points
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=number of channels for spin-polarization (1 or 2)
!!  nsym=number of symmetry elements in space group
!!  walli=initial wall clock time
!!
!! OUTPUT
!!  npwtot(nkpt)=total number of plane waves at each k point
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!    forces and its components, the stress tensor) of a ground-state computation
!!
!! SIDE EFFECTS
!!  acell(3)=unit cell length scales (bohr)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  iexit=exit flag
!!  mpi_enreg=MPI-parallelisation information (some already initialized,
!!            some others to be initialized here)
!!  occ(mband*nkpt*nsppol)=occupation number for each band and k
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!    Before entering the first time in pstate, a significant part of
!!    psps has been initialized. All the remaining components of psps
!!    are to be initialized in the call to pspini .
!!  rprim(3,3)=dimensionless real space primitive translations
!!  vel(3,natom)=value of velocity
!!  xred(3,natom)=reduced atomic coordinates
!!
!! NOTES
!! * This routine drives several calls to gstate.f.
!! * Using restartxf in parallel when localrdwf==0 is not yet possible...
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      clnmpi_gs,gstate,initmpi_gs,mpi_bcast,mpi_comm_create,mpi_comm_group
!!      mpi_comm_rank,mpi_comm_size,mpi_group_incl,status,wrtout,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pstate(acell,codvsn,cpui,dtfil,dtset,iexit,mband,mgfft,mkmem,mpi_enreg,&
&                  mpw,natom,nfft,nkpt,npwtot,nspden,nspinor,nsppol,nsym,occ,&
&                  pawrad,pawtab,psps,results_gs,rprim,vel,walli,xred)

 use defs_basis
 use defs_datatypes
#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_21drive
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer, intent(in) :: mband,mgfft,mkmem,mpw,nfft
 integer, intent(inout) :: iexit,natom,nkpt,nspden,nspinor,nsppol,nsym
 real(dp), intent(in) :: cpui,walli
 character(len=6), intent(in) :: codvsn
 type(MPI_type), intent(inout) :: mpi_enreg
 type(datafiles_type), intent(inout) :: dtfil
 type(dataset_type), intent(inout) :: dtset
 type(pseudopotential_type), intent(inout) :: psps
 type(results_gs_type), intent(out) :: results_gs
 integer, intent(out) :: npwtot(nkpt)
 real(dp), intent(inout) :: acell(3),occ(mband*nkpt*nsppol),rprim(3,3),vel(3,natom)
 real(dp) :: xred(3,natom)
 type(pawrad_type), intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type), intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
 integer,parameter :: level=5
 integer :: ierr,ii,ipack,ipara,jpara,kpara,npack,npara,prtvol
 real(dp) :: total_time
 character(len=500) :: message
#if defined MPI 
           integer,allocatable :: ranks(:)
#endif
 type(pawang_type) :: pawang
 real(dp),allocatable :: s_vel(:,:,:),s_vel_tot(:,:,:),s_xred(:,:,:)
 real(dp),allocatable :: s_xred_tot(:,:,:),sc_vel(:,:,:),sc_xred(:,:,:)
 real(dp),allocatable :: u_vel(:,:,:,:),u_xred(:,:,:,:)

!***************************************************************************

 call status(0,dtfil%filstat,iexit,level,'enter pstate  ')

!Structured debugging if prtvol==-level
 prtvol=dtset%prtvol
 if(prtvol==-level)then
  write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' pstate : enter '
  call wrtout(06,message,'COLL')
 end if

 mpi_enreg%npara=dtset%npara
 npara=dtset%npara
 kpara=dtset%kpara
 npack=dtset%npack
 mpi_enreg%paralbd=0
 allocate(u_xred(3,natom,0:npara,0:1),u_vel(3,natom,0:npara,0:1))
 allocate(s_xred(3,natom,0:npara),s_vel(3,natom,0:npara))
 allocate(sc_xred(3,natom,0:npara),sc_vel(3,natom,0:npara))
#if defined MPI
           allocate(s_xred_tot(3,natom,0:npara),s_vel_tot(3,natom,0:npara))
#endif

!Starting point is accurate
 total_time=0.0_dp
 do ipack=0,npack
  u_xred(:,:,0,0)=xred(:,:)
  u_vel(:,:,0,0)=vel(:,:)
  total_time=dtset%dtion

!First "coarse" step
  write(message,'(a)')'@S0 LEGEND "coars1" '
  call wrtout(75, message, 'COLL')
  mpi_enreg%nproc_per_para=1
  mpi_enreg%nproc_fft=1
#if defined MPI
           call MPI_COMM_GROUP(MPI_COMM_WORLD,mpi_enreg%world_group,ierr)
           allocate(mpi_enreg%proc_distrb_para(0:mpi_enreg%npara-1,nkpt))
           allocate(mpi_enreg%kpt_comm_para(0:mpi_enreg%npara-1))
           allocate(mpi_enreg%kpt_group_para(0:mpi_enreg%npara-1))
           mpi_enreg%proc_distrb_para(0:mpi_enreg%npara-1,1:nkpt)=mpi_enreg%me
           allocate(ranks(1))
           do ipara=0,npara-1
            ranks(1)=mpi_enreg%me
            call MPI_GROUP_INCL(mpi_enreg%world_group,mpi_enreg%nproc_per_para, &
&                               ranks,mpi_enreg%kpt_group_para(ipara),ierr)
            call MPI_COMM_CREATE(MPI_COMM_WORLD,mpi_enreg%kpt_group_para(ipara), &
&                                mpi_enreg%kpt_comm_para(ipara),ierr)
           end do
           deallocate(ranks)
#endif
  do ipara=0,npara-1
   mpi_enreg%ipara=ipara
   total_time=(ipara+npara*ipack)*dtset%dtion
   mpi_enreg%master_group_para=mpi_enreg%me
   mpi_enreg%nproc_group_para=1
   mpi_enreg%num_group_para=1
   mpi_enreg%me_group_para=0
   if (mpi_enreg%me == 0) then
    write(77,*) 'first coarse in x,v',xred(1,1),vel(1,1),mpi_enreg%me
   end if
   dtset%tfkinfunc=1
   dtset%tfnewton=-0.4 !Thomas fermi
   call gstate(acell,codvsn,cpui,dtfil,dtset,&
&              iexit,mpi_enreg,&
&              npwtot,nspinor,occ,pawang,pawrad,pawtab,psps,&
&              results_gs,rprim,vel,walli,xred)
   if (mpi_enreg%me == 0) then
    write(77,*) 'fine gs%etotal', results_gs%etotal
    write(77,*) 'first coarse gs%etotal', results_gs%etotal
    write(75,*) total_time+dtset%dtion,xred(1,1)
   end if
   u_xred(:,:,ipara+1,0)=xred(:,:)
   u_vel(:,:,ipara+1,0)=vel(:,:)
  end do

!Parareal loop
  sc_xred(:,:,:)=zero
  sc_vel(:,:,:)=zero
  s_xred(:,:,:)=zero
  s_vel(:,:,:)=zero
  do jpara=0,kpara
   call initmpi_gs(dtset,mpi_enreg)
   write(message,'(a,i1,a,i2,a)') '@s',jpara+1,'  LEGEND "',jpara+1,'"'
   call wrtout(75, message, 'COLL')
   u_xred(:,:,0,1)=u_xred(:,:,0,0)
   u_vel(:,:,0,1)=u_vel(:,:,0,0)
!  Loop to parallelize
   mpi_enreg%jpara=jpara
   do ipara=0,npara-1
    mpi_enreg%ipara=ipara
    xred(:,:)=u_xred(:,:,ipara,0)
    vel(:,:)=u_vel(:,:,ipara,0)
#if defined MPI 
           if (minval(abs(mpi_enreg%proc_distrb_para(ipara,1:nkpt)-mpi_enreg%me))==0) then
            if ((mpi_enreg%nproc > mpi_enreg%nproc_per_para*mpi_enreg%npara) .and. &
&               (mpi_enreg%me >= mpi_enreg%nproc_per_para*mpi_enreg%npara)) then
             mpi_enreg%num_group_para=0
             mpi_enreg%me_group_para=-1
             mpi_enreg%nproc_group_para=-1
             mpi_enreg%master_group_para=-1
            else
             call MPI_COMM_RANK(mpi_enreg%kpt_comm_para(mpi_enreg%ipara), &
&                               mpi_enreg%me_group_para,ierr)
             if (mpi_enreg%me_group_para == 0)  then
              mpi_enreg%master_group_para=mpi_enreg%me
             else
              mpi_enreg%master_group_para=0
             end if
             call wrtout(mpi_enreg%master_group_para, message, 'INIT')
             call MPI_BCAST(mpi_enreg%master_group_para,1,MPI_INTEGER,&
&                           0,mpi_enreg%kpt_comm_para(mpi_enreg%ipara),ierr)
             call MPI_COMM_SIZE(mpi_enreg%kpt_comm_para(mpi_enreg%ipara), &
&                               mpi_enreg%nproc_group_para,ierr)
            end if
#endif
    dtset%tfkinfunc=0
    dtset%tfnewton=0.0
    call gstate(acell,codvsn,cpui,dtfil,dtset,&
&               iexit,mpi_enreg,&
&               npwtot,nspinor,occ,pawang,pawrad,pawtab,&
&               psps,results_gs,rprim,vel,walli,xred)
#if defined MPI
            else
             s_xred(:,:,ipara+1)=0.0_dp
             s_vel(:,:,ipara+1)=0.0_dp
             cycle
            end if
#endif
    s_xred(:,:,ipara+1)=xred(:,:)-u_xred(:,:,ipara+1,0)
    s_vel(:,:,ipara+1)=vel(:,:)-u_vel(:,:,ipara+1,0)
   end do

!Broadcast results
#if defined MPI 

           call xsum_mpi(s_xred,s_xred_tot,3*natom*(npara+1), MPI_COMM_WORLD,ierr)
           call xsum_mpi(s_vel,s_vel_tot  ,3*natom*(npara+1), MPI_COMM_WORLD,ierr)

           s_xred(:,:,:)=s_xred_tot(:,:,:)
           s_vel(:,:,:)=s_vel_tot(:,:,:)
           call wrtout(0, message, 'INIT')
#endif
   call clnmpi_gs(dtset,mpi_enreg)

#if defined MPI 
           mpi_enreg%nproc_per_para=1
           call MPI_COMM_GROUP(MPI_COMM_WORLD,mpi_enreg%world_group,ierr)
           allocate(mpi_enreg%proc_distrb_para(0:mpi_enreg%npara-1,nkpt))
           allocate(mpi_enreg%kpt_comm_para(0:mpi_enreg%npara-1))
           allocate(mpi_enreg%kpt_group_para(0:mpi_enreg%npara-1))
           mpi_enreg%proc_distrb_para(0:mpi_enreg%npara-1,1:nkpt)=mpi_enreg%me
           allocate(ranks(1))
           do ipara=0,npara-1
            ranks(1)=mpi_enreg%me
            call MPI_GROUP_INCL(mpi_enreg%world_group,mpi_enreg%nproc_per_para, &
&                               ranks,mpi_enreg%kpt_group_para(ipara),ierr)
            call MPI_COMM_CREATE(MPI_COMM_WORLD,mpi_enreg%kpt_group_para(ipara), &
&                                mpi_enreg%kpt_comm_para(ipara),ierr)
           end do
           deallocate(ranks)
#endif

!Coarse integration
   do ipara=0,npara-1
    mpi_enreg%ipara=ipara
    xred(:,:)=u_xred(:,:,ipara,1)
    vel(:,:)=u_vel(:,:,ipara,1)
    mpi_enreg%master_group_para=mpi_enreg%me
    mpi_enreg%nproc_group_para=1
    mpi_enreg%num_group_para=1
    mpi_enreg%me_group_para=0
    if (mpi_enreg%me == 0) then
     write(77,*)'coarse in, x,v',xred(1,1),vel(1,1)
    end if
    dtset%tfkinfunc=1
    dtset%tfnewton=-0.4
    call gstate(acell,codvsn,cpui,dtfil,dtset,&
&               iexit,mpi_enreg,&
&               npwtot,nspinor,occ,pawang,pawrad,pawtab,&
&               psps,results_gs,rprim,vel,walli,xred)

    if (mpi_enreg%me == 0) then
     write(77,*) 'coarse gs%etotal', results_gs%etotal
    end if
    sc_xred(:,:,ipara+1)=sc_xred(:,:,ipara+1)+s_xred(:,:,ipara+1)
    sc_vel(:,:,ipara+1)=sc_vel(:,:,ipara+1)+s_vel(:,:,ipara+1)
    u_xred(:,:,ipara+1,1)=xred(:,:)+sc_xred(:,:,ipara+1)
    u_vel(:,:,ipara+1,1)=vel(:,:)+sc_vel(:,:,ipara+1)
    if (mpi_enreg%me == 0) then
     write(75,*) (ipara+npara*ipack)*dtset%dtion,u_xred(1,1,ipara+1,1)
    end if
   end do
   u_xred(:,:,:,0)=u_xred(:,:,:,1)
   u_vel(:,:,:,0)=u_vel(:,:,:,1)
   call clnmpi_gs(dtset,mpi_enreg)

  end do
  xred(:,:)=u_xred(:,:,npara,1)
  vel(:,:)=u_vel(:,:,npara,1)
 end do
 deallocate(u_xred,u_vel)
 deallocate(s_xred,s_vel)
#if defined MPI 
           deallocate(s_xred_tot,s_vel_tot)
#endif
 deallocate(sc_xred,sc_vel)

end subroutine pstate
!!***
