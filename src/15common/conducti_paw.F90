!{\src2tex{textfont=tt}}
!!****f* ABINIT/conducti_paw
!! NAME
!! conducti_paw
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula for PAW formalism
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (VRecoules, PGhosh)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! NOTES
!!  bantot
!!  doccde(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy.
!!  dom=frequency range
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree).
!!  eigen11(2,nkpt,mband,mband,nsppol)=first-order eigenvalues (hartree)
!!  in direction x
!!  eigen12(2,nkpt,mband,mband,nsppol)=first-order eigenvalues (hartree)
!!  in direction y
!!  eigen13(2,nkpt,mband,mband,nsppol)=first-order eigenvalues (hartree)
!!  in direction z
!!  ecut=kinetic energy planewave cutoff (hartree).
!!  entropy= entropy associated with the smearing (adimensional)
!!  fermie= fermi energy (Hartree)
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{2}$).
!!  gmet_inv(3,3)=inverse of reciprocal space metric ($\textrm{bohr}^{2}$).
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1).
!!  kin11= Onsager kinetic coeficient=optical conductivity
!!  kin12= Onsager kinetic coeficient
!!  kin21= Onsager kinetic coeficient
!!  kin22= Onsager kinetic coeficient
!!  Kth=thermal conductivity
!!  mom=number of frequency for conductivity computation
!!  mband=maximum number of bands.
!!  natom = number of atoms in the unit cell.
!!  nband(nkpt*nsppol)=number of bands at each RF k point for each spin.
!!  nelect=number of electrons per unit cell
!!  nkpt=number of k points in the IBZ for this perturbation
!!  ngfft(3)=integer fft box dimensions.
!!  nspinor=number of spinorial components of the wavefunctions.
!!  nsppol=1 for unpolarized, 2 for spin-polarized.
!!  ntypat = number of atom types.
!!  occ(mband*nkpt*nsppol)=occupation number for each band and k.
!!  occopt==option for occupancies
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).sigx(mom,nphicor))
!!  rprimd(3,3)=real space primitive translations.
!!  of primitive translations.
!!  Sth=thermopower
!!  tsmear=smearing width (or temperature) in Hartree
!!  ucvol=unit cell volume in ($\textrm{bohr}^{3}$).
!!  wind=frequency windows for computations of sigma
!!  wtk(nkpt)=weight assigned to each k point.
!!  xspec=calculate spectro x using file *OPT2
!!  znucl(natom)=atomic number of atoms
!!  np_sum=noziere-pines sumrule
!!
!! PARENTS
!!
!! CHILDREN
!!      getnel,hdr_clean,hdr_io,jacobi,metric,wffclose,wffopen,wffreadeigk
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine conducti_paw(filnam,mpi_enreg)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_14iowfdenpot
 use interfaces_15common, except_this_one => conducti_paw
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 character(len=fnlen) :: filnam
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
!no_abirules
 integer :: accesswff,bantot,bd2tot_index,bdtot0_index,bdtot_index,dosdeltae
 integer :: fform0,formeig0,headform,iatom,iband,ierr,i,ii,jj,ikpt
 integer :: index_1,io_status,iom,isppol,jband,l1,l2,mband,mu,natom,nband1
 integer :: master,me,nrot,mom
 integer :: nband_k,nkpt,nlign,npw1,nrest,nspinor,nspinor1,nsppol,ntypat,nu
 integer :: occopt,rdwr,spaceComm,tim_rwwf,lnfm,nphicor,xspec=0
 integer,allocatable :: nband(:),ncor(:),lcor(:)
 real(dp) :: cond_d,deltae,diff_occ,ecut,entropy,etotal,fermie,maxocc
 real(dp) :: nelect,np_sum,np_sum_k1,np_sum_k2,omin,omax,dom,oml,residm,socc,socc_k,sig,swtk
 real(dp) :: tphysel,tsmear,ucvol,Tatm,del
 real(dp) :: gmet(3,3),gmet_inv(3,3),gprimd(3,3),gprimd_inv(3,3),rmet(3,3),rmet_inv(3,3),rprimd(3,3)
 real(dp),allocatable :: cond_nd(:,:,:),dhdk2_r(:,:,:,:),dhdk2_g(:,:)
 real(dp),allocatable ::doccde(:),doccde_k(:),trace(:)
 real(dp),allocatable :: eig0_k(:),eig0tmp(:),eig1_k(:,:,:,:),eigen0(:),eigen11(:)
 real(dp),allocatable :: eigen12(:),eigtmp(:),energy_cor(:)
 real(dp),allocatable :: eigen13(:),occ(:),occ_k(:),wtk(:),cond_tot(:),oml1(:)
 real(dp),allocatable :: kin11(:),kin12(:),kin21(:),kin22(:),sigx(:,:)
 real(dp),allocatable :: kin11_k(:),kin12_k(:),kin21_k(:),kin22_k(:),Kth(:),Stp(:)
 real(dp) :: z(3,3)
 real(dp) :: eig_cond(3)
!
 character(len=fnlen) :: filnam0,filnam1,filnam2,filnam3,filnam_gen,filnam_out
 type(hdr_type) :: hdr
 type(wffile_type) :: wff0,wff1,wff2,wff3
!
 type matrix_elts
  real(dp),pointer :: eigen11(:,:,:)
  real(dp),pointer :: eigen12(:,:,:)
  real(dp),pointer :: eigen13(:,:,:)
 end type matrix_elts
 type(matrix_elts), allocatable :: psinablapsi(:),psinablapsi2(:)

! *********************************************************************************
!BEGIN EXECUTABLE SECTION

!Read data file name
!write(6,'(a)')' Please, give the name of the data file ...'
!read(5, '(a)')filnam
!write(6,'(a)')' The name of the data file is :',filnam
 write(6,'(a)')' Please, give the name of the output file ...'
 read(5, '(a)')filnam_out
 write(6,'(a)')' The name of the output file is :',filnam_out
 lnfm=index(filnam_out,' ')-1
!Read data file
 open(15,file=filnam,form='formatted')
 rewind(15)
 read(15,*)
 read(15,'(a)')filnam_gen       ! generic name for the files
 filnam1=trim(filnam_gen)//'_OPT'
 if(xspec/=0) filnam2=trim(filnam_gen)//'_OPT2'
 filnam0=trim(filnam_gen)//'_WFK'
!read(15,'(a)')filnam1       ! optic file
!read(15,'(a)')filnam2       ! optic2 file
!read(15,'(a)')filnam0       ! ground-state data

!Open the Wavefunction and optic files
!These default values are typical of sequential use
 accesswff=0 ; spaceComm=abinit_comm_serial ; master=0 ; me=0
 call WffOpen(accesswff,spaceComm,filnam0,ierr,wff0,master,me,10)
 call WffOpen(accesswff,spaceComm,filnam1,ierr,wff1,master,me,11)
 if(xspec/=0)call WffOpen(accesswff,spaceComm,filnam2,ierr,wff1,master,me,12)

!Read the header from Ground state file
 rdwr=1
 call hdr_io(fform0,hdr,rdwr,wff0)

!Extract info from the header
 headform=hdr%headform
 bantot=hdr%bantot
 ecut=hdr%ecut_eff
 natom=hdr%natom
 nkpt=hdr%nkpt
 nspinor=hdr%nspinor
 nsppol=hdr%nsppol
 ntypat=hdr%ntypat
 occopt=hdr%occopt
 rprimd(:,:)=hdr%rprimd(:,:)
 allocate(nband(nkpt*nsppol),occ(bantot))
 fermie=hdr%fermie
 occ(1:bantot)=hdr%occ(1:bantot)
 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)

!Get mband, as the maximum value of nband(nkpt)
 mband=maxval(nband(:))

 write(6,*)
 write(6,'(a,3f10.5,a)' )' rprimd(bohr)      =',rprimd(1,1:3)
 write(6,'(a,3f10.5,a)' )'                    ',rprimd(2,1:3)
 write(6,'(a,3f10.5,a)' )'                    ',rprimd(3,1:3)
 write(6,'(a,i8)')       ' natom             =',natom
 write(6,'(a,2i8)')      ' nkpt,mband        =',nkpt,mband
 write(6, '(a, f10.5,a)' ) ' ecut              =',ecut,' Ha'
 write(6,'(a,f10.5,a,f10.5,a)' )' fermie            =',fermie,' Ha',fermie*Ha_eV,' eV'

!Prepare the reading of Wff files
 formeig0=0 ; tim_rwwf=0
 allocate(eigtmp(2*mband*mband),eig0tmp(mband))
!Read the eigenvalues of ground-state
 allocate(eigen0(mband*nkpt*nsppol))
 bdtot0_index=0 ; bdtot_index=0
 do isppol=1,nsppol
  do ikpt=1,nkpt
   nband1=nband(ikpt+(isppol-1)*nkpt)
   call WffReadEigK(eig0tmp,formeig0,headform,ikpt,isppol,mband,mpi_enreg,nband1,tim_rwwf,wff0)
   eigen0(1+bdtot0_index:nband1+bdtot0_index)=eig0tmp(1:nband1)
   bdtot0_index=bdtot0_index+nband1
  end do
 end do
 call WffClose(wff0,ierr)
 deallocate(eig0tmp)
!
 write(6,*)'Reading file optic'
!call hdr_skip(wff1,ierr)
 allocate(psinablapsi(nkpt))
 do isppol=1,nsppol
  do ikpt=1,nkpt
   nband1=nband(ikpt+(isppol-1)*nkpt)
   allocate(psinablapsi(ikpt)%eigen11(2,nband1,nband1))
   allocate(psinablapsi(ikpt)%eigen12(2,nband1,nband1))
   allocate(psinablapsi(ikpt)%eigen13(2,nband1,nband1))
   read(11)psinablapsi(ikpt)%eigen11
   read(11)psinablapsi(ikpt)%eigen12
   read(11)psinablapsi(ikpt)%eigen13
  end do
 end do
 write(6,*)'Reading file optic2 for X'
!call hdr_skip(wff1,ierr)
 if(xspec/=0) then
  read(12) nphicor
  allocate(ncor(nphicor),lcor(nphicor),energy_cor(nphicor))
  do i=1,nphicor
   read(12) ncor(i),lcor(i),energy_cor(i)
  end do
  
  allocate(psinablapsi2(nkpt))
  
  do isppol=1,nsppol
   do ikpt=1,nkpt
    nband1=nband(ikpt+(isppol-1)*nkpt)
    allocate(psinablapsi2(ikpt)%eigen11(2,nband1,nphicor))
    allocate(psinablapsi2(ikpt)%eigen12(2,nband1,nphicor))
    allocate(psinablapsi2(ikpt)%eigen13(2,nband1,nphicor))
    read(12)psinablapsi2(ikpt)%eigen11
    read(12)psinablapsi2(ikpt)%eigen12
    read(12)psinablapsi2(ikpt)%eigen13
   end do
  end do
 end if
 call WffClose(wff1,ierr)
!DEBUG
!do isppol=1,nsppol
!do ikpt=1,nkpt
!nband1=nband(ikpt+(isppol-1)*nkpt)
!do iband=1,nband1
!do jband=1,nphicor
!write(89,*)'real',psinablapsi2(ikpt)%eigen11(1,iband,jband)&
!&            ,psinablapsi2(ikpt)%eigen12(1,iband,jband),psinablapsi2(ikpt)%eigen13(1,iband,jband)
!write(89,*)'imag',psinablapsi2(ikpt)%eigen11(2,iband,jband)&
!&            ,psinablapsi2(ikpt)%eigen12(2,iband,jband),psinablapsi2(ikpt)%eigen13(2,iband,jband)
!enddo
!enddo
!enddo
!enddo
!ENDDEBUG
!---------------------------------------------------------------------------------
!gmet inversion
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!---------------------------------------------------------------------------------
!derivative of occupation wrt the energy.
 allocate(doccde(mband*nkpt*nsppol),wtk(nkpt))

 read(15,*)tsmear
 Tatm=tsmear*Ha_K
 write(6,'(a,f12.5,a,f12.5,a)') ' Temp              =',tsmear,' Ha ',Tatm,' Kelvin'    
!
 nlign=nkpt/6
 nrest=nkpt-6*nlign
 index_1=0
 do ii=1,nlign
  read(15,*)wtk(1+index_1:6+index_1)
  index_1=index_1+6
 end do
 if (nrest/=0) then
  read(15,*)wtk(6*nlign+1:nkpt)
 end if
!
 if (occopt==1) then
  write(6,'(a,i4)')  ' occopt            =',occopt
  doccde=0.0d0
 else
  tphysel=zero
  maxocc=two/(nsppol*nspinor)
  dosdeltae=zero
! call getnel(doccde,dosdeltae,eigen0,entropy,fermie,maxocc,mband,nband,&
! &  nelect,nkpt,nsppol,occ,occopt,1,tphysel,tsmear,11,wtk)
! DEBUG
! write(6,'(a,f10.5)')' getnel : nelect   =',nelect
! ENDDEBUG
 end if
!---------------------------------------------------------------------------------
!size of the frequency range
 read(15,*)dom,omin,omax,mom
 close(15)
 del=(omax-omin)/(mom-1)
 allocate(oml1(mom))
 do iom=1,mom
  oml1(iom)=omin+dble(iom-1)*del
 end do

 if(xspec/=0) then
  allocate(sigx(mom,nphicor))
 end if
 allocate(cond_nd(mom,3,3))
 allocate(trace(mom),cond_tot(mom))
 allocate(kin11(mom),kin12(mom),kin21(mom),kin22(mom))
 allocate(kin11_k(mom),kin12_k(mom),kin21_k(mom),kin22_k(mom))
 allocate(Kth(mom),Stp(mom))
 write(6,'(a,i8,3f10.5,a)')' npts,omin,omax,width      =',mom,omin,omax,dom,' Ha'


!---------------------------------------------------------------------------------


!
 
 kin11   = 0.0d0
 kin12   = 0.0d0
 kin21   = 0.0d0
 kin22   = 0.0d0
 np_sum  = 0.0d0
 socc    = 0.0d0

!
!LOOP OVER SPINS

 do isppol=1,nsppol
! 
  bdtot_index = 0
  bd2tot_index = 0
  deltae  = 0.0d0
! BIG FAT k POINT LOOP
! 
  do ikpt=1,nkpt
!  
   nband_k=nband(ikpt+(isppol-1)*nkpt)
!  
   allocate(eig0_k(nband_k),eig1_k(2,nband_k,nband_k,3))
   allocate(occ_k(nband_k),doccde_k(nband_k))
   allocate(dhdk2_r(nband_k,nband_k,3,3),dhdk2_g(nband_k,nband_k))

   cond_nd   = 0.0d0
   kin11_k   = 0.0d0
   kin12_k   = 0.0d0
   kin21_k   = 0.0d0
   kin22_k   = 0.0d0
   np_sum_k1 = 0.0d0
   np_sum_k2 = 0.0d0
   socc_k    = 0.0d0
   dhdk2_r   = 0.0d0
   dhdk2_g   = 0.0d0
!  
!  eigenvalue for k-point
   eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
!  first derivative eigenvalues for k-point
   do iband=1,nband_k
    do jband=1,nband_k
     eig1_k(1,iband,jband,1)=psinablapsi(ikpt)%eigen11(1,iband,jband)
     eig1_k(1,iband,jband,2)=psinablapsi(ikpt)%eigen12(1,iband,jband)
     eig1_k(1,iband,jband,3)=psinablapsi(ikpt)%eigen13(1,iband,jband)
     eig1_k(2,iband,jband,1)=psinablapsi(ikpt)%eigen11(2,iband,jband)
     eig1_k(2,iband,jband,2)=psinablapsi(ikpt)%eigen12(2,iband,jband)
     eig1_k(2,iband,jband,3)=psinablapsi(ikpt)%eigen13(2,iband,jband)
    end do
   end do
!  occupation numbers for k-point
   occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
!  derivative of occupation number for k-point
   doccde_k(:)=doccde(1+bdtot_index:nband_k+bdtot_index)

!  LOOP OVER BAND
   do iband=1,nband_k
    do jband=1,nband_k
!    
     do l1=1,3
      do l2=1,3
       dhdk2_r(iband,jband,l1,l2)=dhdk2_r(iband,jband,l1,l2)+(&
&       eig1_k(1,iband,jband,l1)*eig1_k(1,iband,jband,l2)&
&       +eig1_k(2,iband,jband,l1)*eig1_k(2,iband,jband,l2))
      end do
     end do

     do l1=1,3
      dhdk2_g(iband,jband)=dhdk2_g(iband,jband)+( &
&      eig1_k(1,iband,jband,l1)*eig1_k(1,iband,jband,l1) &
&      +eig1_k(2,iband,jband,l1)*eig1_k(2,iband,jband,l1))
     end do

     diff_occ = occ_k(iband)-occ_k(jband)
     if (dabs(diff_occ)>=tol8) then
!     
!     Conductivity for each omega
      omin = 0.0d0
      do iom=1,mom
       oml=oml1(iom)
       if (jband>iband) then
        sig= dhdk2_g(iband,jband)&
&        *(diff_occ)/oml*(dexp(-((eig0_k(jband)-eig0_k(iband)-oml)/dom)**2)&
&        -dexp(-((eig0_k(iband)-eig0_k(jband)-oml)/dom)**2))
        kin11_k(iom)=kin11_k(iom)+sig
        kin12_k(iom)=kin12_k(iom)-sig*(eig0_k(jband)-fermie)
        kin21_k(iom)=kin21_k(iom)-sig*(eig0_k(iband)-fermie)
        kin22_k(iom)=kin22_k(iom) + &
&        sig*(eig0_k(iband)-fermie)*(eig0_k(jband)-fermie)
       end if
       do l1=1,3
        do l2=1,3
         cond_nd(iom,l1,l2)=cond_nd(iom,l1,l2) +dhdk2_r(iband,jband,l1,l2)&
&         *(diff_occ)/oml*dexp(-((eig0_k(jband)-eig0_k(iband)-oml)/dom)**2)
        end do
       end do
      end do
      

!     
!     Sumrule start
      if (dabs(eig0_k(iband)-eig0_k(jband))>=tol10) then
       np_sum_k1=np_sum_k1 -dhdk2_g(iband,jband)&
&       *(diff_occ)/(eig0_k(iband)-eig0_k(jband))
      else
       np_sum_k2=np_sum_k2 - doccde_k(iband)*dhdk2_g(iband,jband)
      end if
!     

!     end loop over band
     end if
    end do
    socc_k=socc_k+occ_k(iband)
   end do
!  
   do iom=1,mom
    kin11(iom)=kin11(iom)+wtk(ikpt)*kin11_k(iom)
    kin12(iom)=kin12(iom)+wtk(ikpt)*kin12_k(iom)
    kin21(iom)=kin21(iom)+wtk(ikpt)*kin21_k(iom)
    kin22(iom)=kin22(iom)+wtk(ikpt)*kin22_k(iom)
    

   end do

   np_sum=np_sum + wtk(ikpt)*(np_sum_k1+np_sum_k2)
   socc=socc+wtk(ikpt)*socc_k
!  
!  validity limit
   deltae=deltae+(eig0_k(nband_k)-fermie)

   bd2tot_index=bd2tot_index+2*nband_k**2
   bdtot_index=bdtot_index+nband_k
   deallocate(eig0_k,eig1_k,occ_k,doccde_k,dhdk2_r,dhdk2_g)
!  End loop over k
  end do
! End loop over Spin
 end do

!spectro X ---------------------------------------------------------------------------------
!LOOP OVER SPINS
 if(xspec/=0) then
  sigx=0.
  do isppol=1,nsppol
   bdtot_index = 0
   bd2tot_index = 0
!  k POINT LOOP
   do ikpt=1,nkpt
    nband_k=nband(ikpt+(isppol-1)*nkpt)
    allocate(eig0_k(nband_k),eig1_k(2,nband_k,nband_k,3))
    allocate(occ_k(nband_k))
    allocate(dhdk2_g(nband_k,nphicor))
    dhdk2_g   = 0.0d0
!   
!   eigenvalue for k-point
    eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
!   first derivative eigenvalues for k-point
    do iband=1,nband_k
     do jband=1,nphicor
      eig1_k(1,iband,jband,1)=psinablapsi2(ikpt)%eigen11(1,iband,jband)
      eig1_k(1,iband,jband,2)=psinablapsi2(ikpt)%eigen12(1,iband,jband)
      eig1_k(1,iband,jband,3)=psinablapsi2(ikpt)%eigen13(1,iband,jband)
      eig1_k(2,iband,jband,1)=psinablapsi2(ikpt)%eigen11(2,iband,jband)
      eig1_k(2,iband,jband,2)=psinablapsi2(ikpt)%eigen12(2,iband,jband)
      eig1_k(2,iband,jband,3)=psinablapsi2(ikpt)%eigen13(2,iband,jband)
     end do
    end do

!   occupation numbers for k-point
    occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)

!   LOOP OVER BAND
    
    do iband=1,nband_k
     do jband=1,nphicor

      do l1=1,3
       dhdk2_g(iband,jband)=dhdk2_g(iband,jband)+( &
&       eig1_k(1,iband,jband,l1)*eig1_k(1,iband,jband,l1) &
&       +eig1_k(2,iband,jband,l1)*eig1_k(2,iband,jband,l1))
      end do

      diff_occ = 2.0-occ_k(iband)
      
!     if (dabs(diff_occ)>=tol8) then
!     
!     Conductivity for each omega
      omin = 0.0d0
      do iom=1,mom
       oml=abs(energy_cor(jband))+oml1(iom)
       sigx(iom,jband)=sigx(iom,jband)+ dhdk2_g(iband,jband) &
&       *(diff_occ)/oml*(-dexp(-((energy_cor(jband)-eig0_k(iband)-oml)/dom)**2)&
&       +dexp(-((eig0_k(iband)-energy_cor(jband)-oml)/dom)**2))         
      end do
     end do
    end do
    
    bd2tot_index=bd2tot_index+2*nband_k**2
    bdtot_index=bdtot_index+nband_k
!   end loop over k
   end do
!  end loop over spin
  end do

  do i=1,nphicor
   do iom=1,mom
    if(sigx(iom,i)<=tiny(1.e-30)) sigx(iom,i)=0.
   end do
  end do 

  sigx=sigx*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
  open(31,file=filnam_out(1:lnfm)//'_sigX',form='formatted')
  
  do iom=1,mom
   write(31,'(6(1x,e14.8))') (27.2*(abs(energy_cor(i))+oml1(iom)),sigx(iom,i),i=1,nphicor)
  end do
 end if
 
 write(6,'(a,3f10.5)')' sumrule           =',np_sum/socc/3.0d0,socc
 write(6,'(a,f10.5,a,f10.5,a)')&
& ' Emax-Efermi       =',deltae/dble(nkpt),' Ha',deltae/dble(nkpt)*Ha_eV,' eV'

 

 
 open(18,file=filnam_out(1:lnfm)//'_12',form='formatted')
 open(19,file=filnam_out(1:lnfm)//'_21',form='formatted')
 open(20,file=filnam_out(1:lnfm)//'_22',form='formatted')
 open(30,file=filnam_out(1:lnfm)//'_sig',form='formatted')
 write(30,'(a)')' # omega(ua) hbar*omega(eV)    cond(ua)             cond(ohm.cm)-1'
 open(41,file=filnam_out(1:lnfm)//'_Kth',form='formatted')
 write(41,'(a)')' #omega(ua) hbar*omega(eV)    thermal cond(ua)      Kth(W/m/K)'
 open(42,file=filnam_out(1:lnfm)//'_Stp',form='formatted')
 write(42,'(a)')' #omega(ua) hbar*omega(eV)    thermopower(ua)       Stp(microohm/K)'
 open(45,file=filnam_out(1:lnfm)//'.out',form='formatted')
 write(45,'(a)' )' Conducti output file:'
 write(45,'(a)' )' Contains all results produced by conducti utility'
 write(45,'(a)' )' '
 write(45,'(a)')' # omega(ua)       cond(ua)             thermal cond(ua)       thermopower(ua)'
!
!call isfile(filnam_out,'new')

!Keep this line : prevent silly (compiler ?) bug on HP 8000
 write(6,*)' conducti : after call isfile '
!
!Compute thermal conductivity and thermopower
 do iom=1,mom
  oml=oml1(iom)
  kin11(iom)=kin11(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
  kin21(iom)=kin21(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
  kin12(iom)=kin12(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
  kin22(iom)=kin22(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
  if (dabs(kin11(iom))<10.0d-20) kin11(iom)=0.0d0
  Kth(iom)=kin22(iom)
  Stp(iom)=zero
  if(kin11(iom)/=zero)  then
   Kth(iom)=Kth(iom)-(kin12(iom)*kin21(iom)/kin11(iom))
   Stp(iom)=kin12(iom)/(kin11(iom)*Tatm)
  end if
  if (dabs(Kth(iom))<10.0d-20) Kth(iom)=0.0d0
  if (dabs(Stp(iom))<10.0d-20) Stp(iom)=0.0d0
  write(18,*)oml,kin12(iom)
  write(19,*)oml,kin21(iom)
  write(20,*)oml,kin22(iom),kin22(iom)/Tatm*3.4057d9
  write(30,'(2f12.5,2es22.12)') oml,oml*Ha_eV,kin11(iom),kin11(iom)*Ohmcm
  write(41,'(2f12.5,2es22.12)') oml,oml*Ha_eV,Kth(iom),Kth(iom)*3.4057d9/Tatm
  write(42,'(2f12.5,2es22.12)') oml,oml*Ha_eV,Stp(iom),Stp(iom)*3.6753d-2
  write(45,'(1f12.5,3es22.12)') oml,kin11(iom),Kth(iom),Stp(iom)
 end do

!Calculate the imaginary part of the conductivity (principal value)
!+derived optical properties.

 call msig (kin11,mom,oml1)
 
 deallocate(nband,oml1)
 deallocate(occ)
 do ikpt=1,nkpt
  deallocate(psinablapsi(ikpt)%eigen11,psinablapsi(ikpt)%eigen12,psinablapsi(ikpt)%eigen13)
 end do
 deallocate (psinablapsi)
 deallocate(eigen0,doccde,wtk)
 deallocate(cond_nd)
 deallocate(kin11,kin22,kin12,kin21,kin11_k,kin22_k,kin12_k,kin21_k,Stp,Kth)

 close(17);close(18);close(19)
 close(20);close(30);close(40)
 close(41);close(42);close(45)

 call hdr_clean(hdr)

 end subroutine conducti_paw
!!***
