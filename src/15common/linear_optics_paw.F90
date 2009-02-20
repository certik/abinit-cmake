!{\src2tex{textfont=tt}}
!!****p* ABINIT/linear_optics_paw
!! NAME
!! linear_optics_paw
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! linear susceptiblity using matrix elements <-i Nabla> obtained from a
!! PAW ground state calculation. It uses formula 17 from Gadoc et al, 
!! Phys. Rev. B 73, 045112 (2006) together with a scissors correction. It uses
!! a Kramers-Kronig transform to compute the real part from the imaginary part, and
!! it will work on all types of unit cells. It outputs all tensor elements of
!! both the real and imaginary parts.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (VRecoules, PGhosh)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  filnam: base of file names to read data from
!!  mpi_enreg: mpi set up variable, not used in this code
!!
!! OUTPUT
!!  _real and _imag output files
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine linear_optics_paw(filnam,mpi_enreg)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_14iowfdenpot
 use interfaces_15gw
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 character(len=fnlen),intent(in) :: filnam
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
!no_abirules
 integer :: accesswff,bantot,bdtot0_index,bdtot_index,bd2tot_index,dic,djc,fform0,formeig0,headform
 integer :: iband,ierr,ii,ikpt,iom,iout,isppol,isym,jband,jj,jom,lnfm,master,me,mband
 integer :: method,mom,natom,nband1,nband_k,nkpt,nspinor,nsppol,nsym,occopt,only_check
 integer :: rdwr,spaceComm,tim_rwwf
 integer,allocatable :: nband(:),symrel(:,:,:)
 real(dp) :: del,diffw,dom,fij,gdelta,omin,omax,paijpbij(2),rcvol,sciss,wij,ucvol
 real(dp) :: diffwp, diffwm
 real(dp) :: e2rot(3,3),gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),rprimdinv(3,3),symd(3,3),symdinv(3,3)
 real(dp),allocatable :: e1(:,:,:),e2(:,:,:,:),epsilon_tot(:,:,:,:),eigen0(:),eigtmp(:),eig0_k(:)
 real(dp),allocatable :: eig1_k(:,:,:,:),eig0tmp(:),kpts(:,:),occ(:),occ_k(:),oml1(:),wtk(:)
 complex,allocatable :: eps_work(:)
!
 character(len=fnlen) :: filnam0,filnam1,filnam_gen,filnam_out
!
 type(hdr_type) :: hdr
 type(wffile_type) :: wff0,wff1
 type matrix_elts
  real(dp),pointer :: eigen11(:,:,:)
  real(dp),pointer :: eigen12(:,:,:)
  real(dp),pointer :: eigen13(:,:,:)
 end type matrix_elts
 type(matrix_elts), allocatable :: psinablapsi(:)

! *********************************************************************************
!BEGIN EXECUTABLE SECTION

 write(6,'(a)')' Give the name of the output file ...'
 read(5, '(a)') filnam_out
 write(6,'(a)')' The name of the output file is :',filnam_out
 lnfm=index(filnam_out,' ')-1

!Read data file
 open(15,file=filnam,form='formatted')
 rewind(15)
 read(15,*)
 read(15,'(a)')filnam_gen       ! generic name for the files
 filnam1=trim(filnam_gen)//'_OPT' ! nabla matrix elements file
 filnam0=trim(filnam_gen)//'_WFK' ! ground state data

!Open the Wavefunction and optic files
!These default values are typical of sequential use
 accesswff=0 ; spaceComm=abinit_comm_serial ; master=0 ; me=0
 call WffOpen(accesswff,spaceComm,filnam0,ierr,wff0,master,me,10)
 call WffOpen(accesswff,spaceComm,filnam1,ierr,wff1,master,me,11)

!Read the header from Ground state file
 rdwr=1
 call hdr_io(fform0,hdr,rdwr,wff0)

!Extract info from the header
 headform=hdr%headform
 bantot=hdr%bantot
 nkpt=hdr%nkpt
 allocate(kpts(3,nkpt),wtk(nkpt))
 kpts(:,:)=hdr%kptns(:,:) 
 wtk(:)=hdr%wtk(:)
 nspinor=hdr%nspinor
 nsppol=hdr%nsppol
 occopt=hdr%occopt
 rprimd(:,:)=hdr%rprimd(:,:)
 rprimdinv(:,:) = rprimd(:,:)
 call matrginv(rprimdinv,3,3) ! need the inverse of rprimd to symmetrize the tensors
 allocate(nband(nkpt*nsppol),occ(bantot))
 occ(1:bantot)=hdr%occ(1:bantot)
 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)
 nsym=hdr%nsym
 allocate(symrel(3,3,nsym))
 symrel(:,:,:)=hdr%symrel(:,:,:)

!Get mband, as the maximum value of nband(nkpt)
 mband=maxval(nband(:))

!get ucvol etc.
 iout = -1
 call metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

 write(6,*)
 write(6,'(a,3f10.5,a)' )' rprimd(bohr)      =',rprimd(1,1:3)
 write(6,'(a,3f10.5,a)' )'                    ',rprimd(2,1:3)
 write(6,'(a,3f10.5,a)' )'                    ',rprimd(3,1:3)
 write(6,*)
 write(6,'(a,3f10.5,a)' )' rprimdinv         =',rprimdinv(1,1:3)
 write(6,'(a,3f10.5,a)' )'                    ',rprimdinv(2,1:3)
 write(6,'(a,3f10.5,a)' )'                    ',rprimdinv(3,1:3)
 write(6,'(a,2i8)')      ' nkpt,mband        =',nkpt,mband

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
 call WffClose(wff1,ierr)

 read(15,*)sciss
 read(15,*)dom,omin,omax,mom
 close(15)
 allocate(oml1(mom),e1(3,3,mom),e2(2,3,3,mom))
 allocate(epsilon_tot(2,3,3,mom),eps_work(mom))
 del=(omax-omin)/(mom-1)
 do iom=1,mom
  oml1(iom)=omin+dble(iom-1)*del
 end do
 write(6,'(a,i8,4f10.5,a)')' npts,omin,omax,width,sciss      =',mom,omin,omax,dom,sciss,' Ha'

!loop over spin components
 do isppol=1,nsppol
  bdtot_index = 0
  bd2tot_index = 0
! loop over k points
  do ikpt=1,nkpt
!  
!  number of bands for this k point
   nband_k=nband(ikpt+(isppol-1)*nkpt)
   allocate(eig0_k(nband_k),eig1_k(2,nband_k,nband_k,3))
   allocate(occ_k(nband_k))
!  eigenvalues for this k-point
   eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
!  occupation numbers for this k-point
   occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
!  values of -i*nabla matrix elements for this k point
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
!  accumulate e2 for this k point, Eq. 17 from PRB 73, 045112 (2006)
   do iband = 1, nband_k
    do jband = 1, nband_k
     fij = occ_k(iband) - occ_k(jband) !occ number difference
     wij = eig0_k(iband) - eig0_k(jband) !energy difference
     if (abs(fij) > zero) then ! only consider states of differing occupation numbers
      do ii = 1, 3
       do jj = 1, 3
        paijpbij(1) = eig1_k(1,iband,jband,ii)*eig1_k(1,iband,jband,jj) + &
&        eig1_k(2,iband,jband,ii)*eig1_k(2,iband,jband,jj)
        paijpbij(2) = eig1_k(2,iband,jband,ii)*eig1_k(1,iband,jband,jj) - &
&        eig1_k(1,iband,jband,ii)*eig1_k(2,iband,jband,jj)
        do iom = 1, mom
!        original version
!        diffw = wij + sciss - oml1(iom) ! apply scissors term here
!        gdelta = exp(-diffw*diffw/(4.0*dom*dom))/(2.0*dom*sqrt(pi)) ! delta fnc resolved as Gaussian
!        e2(1,ii,jj,iom) = e2(1,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(1)*gdelta/(oml1(iom)*oml1(iom))
!        e2(2,ii,jj,iom) = e2(2,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(2)*gdelta/(oml1(iom)*oml1(iom))
         diffwm = wij - sciss + oml1(iom) ! apply scissors term here
         diffwp = wij + sciss - oml1(iom) ! apply scissors term here
         gdelta = exp(-diffwp*diffwp/(4.0*dom*dom))/(2.0*dom*sqrt(pi))
         e2(1,ii,jj,iom) = e2(1,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(1)*gdelta/(wij*wij)
         e2(2,ii,jj,iom) = e2(2,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(2)*gdelta/(wij*wij)
        end do ! end loop over spectral points
       end do ! end loop over jj = 1, 3
      end do ! end loop over ii = 1, 3
     end if ! end selection on fij /= 0
    end do ! end loop over jband
   end do ! end loop over iband

   deallocate(eig0_k,eig1_k,occ_k)
   bd2tot_index=bd2tot_index+2*nband_k**2
   bdtot_index=bdtot_index+nband_k
  end do ! end loop over k points
 end do ! end loop over spin polarizations

!here apply nsym symrel transformations to reconstruct full tensor from IBZ part
 epsilon_tot(:,:,:,:) = zero
 do isym = 1, nsym
  symd(:,:)=matmul(rprimd(:,:),matmul(symrel(:,:,isym),rprimdinv(:,:)))
  symdinv(:,:)=symd(:,:)
  call matrginv(symdinv,3,3)
  do iom = 1, mom
   e2rot(:,:)=matmul(symdinv(:,:),matmul(e2(1,:,:,iom),symd(:,:)))
   epsilon_tot(2,:,:,iom) = epsilon_tot(2,:,:,iom)+e2rot(:,:)/nsym
  end do
 end do

!generate e1 from e2 via KK transforma
 method=0 ! use naive integration ( = 1 for simpson)
 only_check=0 ! compute real part of eps in kk routine
 do ii = 1, 3
  do jj = 1, 3
   eps_work(:) = cmplx(0.0,epsilon_tot(2,ii,jj,:))
   call kramerskronig(mom,oml1,eps_work,method,only_check)
   epsilon_tot(1,ii,jj,:) = real(eps_work(:))
   if (ii /= jj) epsilon_tot(1,ii,jj,:) = epsilon_tot(1,ii,jj,:)- 1.0
  end do ! end loop over jj
 end do ! end loop over ii

 open(18,file=filnam_out(1:lnfm)//'_imag',form='formatted')
 open(19,file=filnam_out(1:lnfm)//'_real',form='formatted')

 write(18,'(a12,6a13)')' # Energy/Ha ','eps_2_xx','eps_2_yy','eps_2_zz',&
& 'eps_2_yz','eps_2_xz','eps_2_xy'
 write(19,'(a12,6a13)')' # Energy/Ha ','eps_1_xx','eps_1_yy','eps_1_zz',&
& 'eps_1_yz','eps_1_xz','eps_1_xy'

 do iom = 1, mom
  write(18,'(ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4)') oml1(iom),' ',&
&  epsilon_tot(2,1,1,iom),' ',epsilon_tot(2,2,2,iom),' ',epsilon_tot(2,3,3,iom),' ',&
&  epsilon_tot(2,2,3,iom),' ',epsilon_tot(2,1,3,iom),' ',epsilon_tot(2,1,2,iom)
  write(19,'(ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4)') oml1(iom),' ',&
&  epsilon_tot(1,1,1,iom),' ',epsilon_tot(1,2,2,iom),' ',epsilon_tot(1,3,3,iom),' ',&
&  epsilon_tot(1,2,3,iom),' ',epsilon_tot(1,1,3,iom),' ',epsilon_tot(1,1,2,iom)
 end do 

 close(18)
 close(19)



 deallocate(nband,oml1)
 deallocate(e2,e1)
 deallocate(occ)
 do ikpt=1,nkpt
  deallocate(psinablapsi(ikpt)%eigen11,psinablapsi(ikpt)%eigen12,psinablapsi(ikpt)%eigen13)
 end do
 deallocate (psinablapsi)
 deallocate(eigen0,wtk,kpts)


 call hdr_clean(hdr)

 end subroutine linear_optics_paw
!!***
