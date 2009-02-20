!{\src2tex{textfont=tt}}
!!****f* ABINIT/suscep
!! NAME
!! suscep
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations of the polarisability
!! within the random phase approximation (RPA)
!!
!! COPYRIGHT
!! Copyright (C) 2000-2008 ABINIT group (XG,GMR,MF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  mband =maximum number of bands
!!  mgfft =maximum single fft dimension
!!  mpi_enreg=informations about MPI parallelization
!!  mpw   =maximum number of planewaves in basis sphere (large number)
!!  natom =number of atoms in unit cell
!!  nkpt  =number of k points
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=number of channels for spin-polarization (1 or 2)
!!  nsym=number of symmetry elements in space group
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (no direct output : results written)
!!
!! SIDE EFFECTS
!!  mkmem =maximum number of k points which can fit in core memory
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      distrb2,getfreqsus,getmpw,getng,hdr_clean,inwffil3,ioarr,kpgio,metric
!!      mkrdim,newocc,setsym,sphereboundary,status,suscep_dyn,suscep_stat,timab
!!      wffclose,wrtout,xcacfd,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine suscep(dtfil,dtset,iexit,&
& mband,mgfft,mkmem,mpi_enreg,mpw,natom,nfft,nkpt,&
& nspden,nspinor,nsppol,nsym,occ,xred,ncdims,ngfft)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_13recipspace
 use interfaces_14iowfdenpot
 use interfaces_14occeig
 use interfaces_17suscep, except_this_one => suscep
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iexit,mband,mgfft,mpw,natom,nfft,nkpt,nspden,nsppol,nsym
 integer,intent(inout) :: mkmem,nspinor
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(vardims_type),intent(inout) :: ncdims
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=3
 integer :: accessfil,dielop,fformr,ierr,ifreq,ii,ipw1,ipw2,isp,isp1,isp2
 integer :: lmax_diel,master,me,mgfftdiel,nband_mx,ncblocks,nfftdiel,nfreqsus
 integer :: npwdiel,prtvol,rdwr,rdwrpaw,susopt
 real(dp) :: diecut,dummy,ecutsus,entropy,etotal,fermie,residm,ucvol
 character(len=500) :: message
 character(len=fnlen) :: filkgs,filsustr,kgnam
 type(dens_sym_operator_type) :: densymop_diel
 type(hdr_type) :: hdr
 type(pawang_type) :: pawang
 type(wffile_type) :: wff1
!arrays
 integer :: ngfftdiel(18),npwarr_diel(1),npwtot_diel(1)
 integer,allocatable :: atindx1_dum(:),dimcprj(:),gbound_diel(:,:)
 integer,allocatable :: indsym(:,:,:),irrzondiel(:,:,:),kg(:,:),kg_diel(:,:)
 integer,allocatable :: nband_dum(:),npwarr(:),npwtot(:),symrec(:,:,:)
 real(dp) :: dielar(7),gmet(3,3),gprimd(3,3),kpt_diel(3),rmet(3,3),rprimd(3,3)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cg(:,:),dielinv(:,:,:,:,:),doccde(:),eigen(:),freq(:)
 real(dp),allocatable :: ph1ddiel(:,:),phnonsdiel(:,:,:),rhor(:,:)
 real(dp),allocatable :: susd_non_dyn(:,:,:,:),susmat(:,:,:,:,:)
 real(dp),allocatable :: susmat_dyn(:,:,:,:,:,:),wght_freq(:),ylmdiel(:,:)
 type(cprj_type),allocatable :: cprj(:,:)
 type(pawrhoij_type),allocatable :: rhoij_dum(:)
 type(pawtab_type),allocatable :: pawtab(:)

!***********************************************************************

 call timab(84,1,tsec)
 call status(0,dtfil%filstat,iexit,level,'enter         ')

 mpi_enreg%paralbd=0
 mpi_enreg%me_fft=0
 mpi_enreg%nproc_fft=1
 mpi_enreg%paral_fft=0
 mpi_enreg%paral_level=2

!
!If dtset%accesswff == 2 set all array outputs to netcdf format
!
 accessfil = 0
 if (dtset%accesswff == 2) then
  accessfil = 1
 end if
 if (dtset%accesswff == 3) then
  accessfil = 3
 end if

 master=0
!Init me
 call xme_init(mpi_enreg,me)

 if(mpi_enreg%paral_compil_kpt==1)then
  allocate(mpi_enreg%proc_distrb(nkpt,mband,nsppol))
  mpi_enreg%parareel=0
  call distrb2(mband, dtset%nband, nkpt, nsppol, mpi_enreg)
 end if

!Structured debugging if prtvol==-level
 prtvol=dtset%prtvol
 if(prtvol==-level)then
  write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,&
&  ' suscep : enter , debug mode '
  call wrtout(06,message,'COLL')
 end if

!Loop input variables
 nfreqsus=dtset%nfreqsus

 dielar(1)=dtset%diecut
 dielar(2)=dtset%dielng
 dielar(3)=dtset%diemac
 dielar(4)=dtset%diemix
 dielar(5)=dtset%diegap
 dielar(6)=dtset%dielam

!Unit numbers and file name for the _KGS file
 filkgs=trim(dtfil%filnam_ds(5))//'_KGS'
 filsustr=trim(dtfil%filnam_ds(5))//'_SUSTR'

!Impose mkmem=0 to read cg from disk
 mkmem=0

 call status(0,dtfil%filstat,iexit,level,'call inwffil3 ')

!Compute different geometric tensor, as well as ucvol, from rprimd
 call mkrdim(dtset%acell_orig,dtset%rprim_orig,rprimd)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Get diecut, and the fft grid to be used for the susceptibility computation
 diecut=abs(dielar(1))
 if( dielar(1)<0.0_dp )then
  ecutsus= dtset%ecut
 else
  ecutsus= ( sqrt( dtset%ecut) *0.5_dp + sqrt(diecut) *0.25_dp )**2
 end if

 ngfftdiel(1:3)=0 ; ngfftdiel(7)=101 ; ngfftdiel(8:18)=dtset%ngfft(8:18)
 call getng(dtset%boxcutmin,ecutsus,gmet,mpi_enreg%me_fft,mgfftdiel,nfftdiel,ngfftdiel,&
& mpi_enreg%nproc_fft,nsym,mpi_enreg%fft_option_lob,mpi_enreg%paral_fft,dtset%symrel)

!Compute the size of the dielectric matrix
 kpt_diel(1:3)=(/ 0.0_dp, 0.0_dp, 0.0_dp /)
 call getmpw(diecut,dtset%exchn2n3d,gmet,(/1/),kpt_diel,&
& mpi_enreg,npwdiel,1,ucvol)

!Now, performs allocation
 allocate(cg(2,mpw*nspinor*mband*mkmem*nsppol))
 allocate(eigen(mband*nkpt*nsppol))
 allocate(kg(3,mpw*mkmem))
 allocate(kg_diel(3,npwdiel))
 allocate(npwarr(nkpt),npwtot(nkpt))
 allocate(gbound_diel(2*mgfftdiel+8,2))
 allocate(irrzondiel(nfftdiel**(1-1/nsym),2,nspden/nsppol))
 allocate(phnonsdiel(2,nfftdiel**(1-1/nsym),nspden/nsppol))
 allocate(nband_dum(nsppol))

!Then, initialize and compute the values of different arrays
 call kpgio(dtset%ecut,dtset%exchn2n3d,gmet,dtset%istwfk,kg,filkgs,dtset%kptns,mkmem,&
& dtset%nband,nkpt,'PERS',mpi_enreg,mpw,npwarr,npwtot,nsppol,dtfil%unkg)
!This kpgio call for going from the suscep FFT grid to the diel sphere
!Note : kgnam is dummy, npwarr_diel is dummy, npwtot_diel is dummy, nband_dum is dummy
 nband_dum(:) = 1
 call kpgio(diecut,dtset%exchn2n3d,gmet,(/1/),kg_diel,kgnam,&
& kpt_diel,1,nband_dum,1,'COLL',mpi_enreg,npwdiel,npwarr_diel,npwtot_diel,&
& nsppol,tmp_unit)

 call sphereboundary(gbound_diel,1,kg_diel,mgfftdiel,npwdiel)

 allocate(indsym(4,nsym,natom),symrec(3,3,nsym))
 if (nsym>1) then
  call setsym(densymop_diel,indsym,irrzondiel,dtset%iscf,natom,&
&  nfftdiel,ngfftdiel,nspden,nsppol,nsym,phnonsdiel,&
&  dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)
 end if

!Read eigenvalues from the wavefunction file
!Also, initialize wff1 and hdr
 eigen(:)=0.0_dp
!mpi_enreg%paralbd=0
 call inwffil3(dtset,eigen,hdr,dtset%istwfk,mband,mpi_enreg,dtset%nband,&
& nkpt,npwarr,nsppol,prtvol,wff1,dtfil%unwff1,dtfil%fnamewffk)

!Compute new occupation numbers if needed
 allocate(doccde(mband*nkpt*nsppol))
 if(dtset%occopt>=3.and.dtset%occopt<=7) then
  call status(0,dtfil%filstat,iexit,level,'call newocc   ')
  call newocc(doccde,eigen,entropy,fermie,dtset%fixmom,mband,dtset%nband,&
&  dtset%nelect,nkpt,nspinor,nsppol,occ,dtset%occopt,prtvol,dtset%stmbias,&
&  dtset%tphysel,dtset%tsmear,dtset%wtk)
 end if

 dielop=2 ! Immediate computation of dielectric matrix
 if(nfreqsus==0) then

! Perform allocations
  allocate(dielinv(2,npwdiel,nspden,npwdiel,nspden))
  allocate(susmat(2,npwdiel,nspden,npwdiel,nspden))
  susmat(:,:,:,:,:)=0._dp

! Compute the static susceptibility matrix
  lmax_diel=0;allocate(atindx1_dum(dtset%natom))
  call suscep_stat(atindx1_dum,cg,cprj,densymop_diel,&
&  dielar,dielop,dimcprj,doccde,eigen,gbound_diel,gprimd,&
&  irrzondiel,dtset%istwfk,kg,kg_diel,lmax_diel,&
&  mband,mgfftdiel,mkmem,mpi_enreg,mpw,natom,dtset%nband,nfftdiel,ngfftdiel,&
&  nkpt,npwarr,npwdiel,nspden,nspinor,nsppol,nsym,dtset%ntypat,&
&  occ,dtset%occopt,0,pawang,pawtab,phnonsdiel,ph1ddiel,prtvol,&
&  susmat,dtset%symafm,dtset%symrel,dtset%tnons,dtset%typat,ucvol,&
&  dtfil%unkg,0,dtfil%unpaw,0,wff1,dtset%wtk,ylmdiel)
  deallocate(atindx1_dum)

! Print the susceptibility matrix
  do isp1=1,nspden
   do isp2=1,nspden
    write(6,'(5x,a,2i2)') 'Susceptibility matrix for spins=',isp1,isp2
    write(6,'(9x,a,13x,a,10x,a,10x,a)') "g","g'","real","imag"
    do ipw1=1,10
     do ipw2=ipw1,10
      write(6,'(2x,3i4,2x,3i4,2x,f12.8,2x,f12.8)') &
&      kg_diel(1:3,ipw1),kg_diel(1:3,ipw2),&
&      susmat(1,ipw1,isp1,ipw2,isp2),susmat(2,ipw1,isp1,ipw2,isp2)
     end do
    end do
   end do
  end do

! Perform deallocations
  deallocate(dielinv)
  deallocate(susmat)
  deallocate(nband_dum)

 else if(nfreqsus > 0) then

! Perform allocations
  allocate(freq(nfreqsus))
  allocate(rhor(nfft,nspden))
  allocate(wght_freq(nfreqsus))

! Perform initializations
  freq(:)=0.0_dp
  rhor(:,:)=0._dp
  wght_freq(:)=0.0_dp

! Read in the density, needed for ALDA kernel
  call status(0,dtfil%filstat,iexit,level,'call ioarr    ')
! Read rho(r) from a disk file
! Unit numbers and file name for the _KGS file
  rdwr=1;rdwrpaw=0
! Note : etotal is read here, and might serve in the tddft routine.
  fformr=52
  call ioarr(accessfil,rhor, dtset, etotal,fformr,dtfil%fildensin,hdr, mpi_enreg, &
&  nfft,rhoij_dum,rdwr,rdwrpaw,ngfft)
  call status(0,dtfil%filstat,iexit,level,'call fourdp   ')

! DEBUG
! leave in for MF to check the density
! dummy=0._dp
! do ii=1,nfft
! dummy=dummy+rhor(ii,1)
! end do
! write(6,*) '%suscep: nfft=',nfft
! write(6,*) '%suscep: dummy=',dummy*ucvol/dble(nfft)
! write(6,*) '%suscep: ucvol=',ucvol
! call flush(6)
! Compute up+down rho(G) by fft
! allocate(work(nfft))
! work(:)=rhor(:,1)
! call fourdp(1,rhog,work,-1,mpi_enreg,nfft,ngfft,0)
! deallocate(work)
! ENDDEBUG

! Create frequency grid and weights, 2 stands for preassigned grid
  call getfreqsus(freq,wght_freq,nfreqsus,dtset%optfreqsus,dtset%freqsuslo,dtset%freqsusin)

! DEBUG
! write(6,*)' suscep : after getfreqsus '
! call flush(6)
! ENDDEBUG

  call xcacfd(cg,densymop_diel,dielar,dielop,doccde,&
&  dtfil,dtset,eigen,filsustr,freq,gbound_diel,gmet,&
&  gprimd,irrzondiel,kg,kg_diel,mband,mgfftdiel,mkmem,&
&  mpi_enreg,mpw,nfft,nfftdiel,nfreqsus,dtset%ngfft,&
&  ngfftdiel,nkpt,npwarr,npwdiel,nspden,nspinor,nsppol,&
&  nsym,occ,phnonsdiel,rhor,rmet,rprimd,ucvol,wff1,wght_freq,ncdims)

! Perform deallocations
  deallocate(freq)
  deallocate(rhor)
  deallocate(wght_freq)

 else

! Leave intact for testing purposes
  nfreqsus=abs(nfreqsus)

! Perform allocations
  allocate(freq(nfreqsus))
  allocate(susd_non_dyn(2,npwdiel,nspden,nfreqsus))
  allocate(susmat_dyn(2,npwdiel,nspden,npwdiel,nspden,nfreqsus))

! Perform initializations
  freq(:)=0.0_dp
  susmat_dyn(:,:,:,:,:,:)=0.0_dp

! Create a linear frequency grid
  freq(1)=dtset%freqsuslo
  do ifreq=2,nfreqsus
   freq(ifreq)=freq(ifreq-1)+dtset%freqsusin
  end do

! DEBUG
! write(6,*) '%suscep: nband_mx=', nband_mx
! ENDDEBUG

! Compute the dynamical susceptibility matrices
  call suscep_dyn(cg,densymop_diel,dielar,dielop,doccde,dtset,&
&  eigen,freq,gbound_diel,&
&  irrzondiel,dtset%istwfk,kg,kg_diel,&
&  mband,mgfftdiel,mkmem,mpi_enreg,mpw,dtset%nband,nband_mx,nfftdiel,nfreqsus,&
&  ngfftdiel,nkpt,npwarr,npwdiel,nspden,nspinor,nsppol,nsym,&
&  occ,dtset%occopt,phnonsdiel,prtvol,&
&  susopt,susd_non_dyn,susmat_dyn,dtset%symafm,&
&  dtset%symrel,dtset%tnons,ucvol,dtfil%unkg,wff1,dtset%wtk)

! Print the dynamical susceptibility matrices
  do ifreq=1,nfreqsus
   write(6,'(/,2x,a,f12.8,a)') '---Susceptibility matrices for frequency=',freq(ifreq),'i'
   do isp1=1,nspden
    do isp2=1,nspden
     write(6,'(5x,a,2i2)') 'Susceptibility matrix for spins=',isp1,isp2
     write(6,'(9x,a,13x,a,10x,a,10x,a)') "g","g'","real","imag"
     do ipw1=1,10
      do ipw2=ipw1,10
       write(6,'(2x,3i4,2x,3i4,2x,f12.8,2x,f12.8)') &
&       kg_diel(1:3,ipw1),kg_diel(1:3,ipw2),&
&       susmat_dyn(1,ipw1,isp1,ipw2,isp2,ifreq),susmat_dyn(2,ipw1,isp1,ipw2,isp2,ifreq)
      end do
     end do
    end do
   end do
  end do

! Perform deallocations
  deallocate(freq)
  deallocate(susd_non_dyn)
  deallocate(susmat_dyn)

 end if !condition nfreqsus

!Performs deallocations
 deallocate(cg,doccde,eigen,gbound_diel)
 deallocate(indsym,irrzondiel)
 deallocate(kg_diel,kg,npwarr,npwtot,phnonsdiel,symrec)

 if(mkmem==0)then
! Sequential case
  if(mpi_enreg%paral_compil_kpt==0)then
!  Unit dtfil%unkg was opened in kpgio
   close (unit=dtfil%unkg,status='delete')

!  Parallel case
  else if(mpi_enreg%paral_compil_kpt==1)then

!  All procs close the file dtfil%unkg
   close(unit=dtfil%unkg)
   if(mpi_enreg%me==0)then
!   only proc 0 delete the file
    open(unit=dtfil%unkg,file=filkgs,form='unformatted',status='unknown')
    close(unit=dtfil%unkg,status='delete')
   end if

  end if
 end if

 call WffClose(wff1,ierr)

!Clean the header
 call hdr_clean(hdr)

 if(mpi_enreg%paral_compil_kpt==1)then
  deallocate(mpi_enreg%proc_distrb)
 end if

 write(message, '(a,a)' ) ch10,' suscep : exiting '
 call wrtout(06,message,'COLL')

 call status(0,dtfil%filstat,iexit,level,'exit          ')
 call timab(84,2,tsec)

end subroutine suscep
!!***
