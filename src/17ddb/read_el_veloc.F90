!{\src2tex{textfont=tt}}
!!****f* ABINIT/read_el_veloc
!!
!! NAME
!! read_el_veloc
!!
!! FUNCTION
!! This routine reads the velocities of the electronic GS
!! for all kpts and bands
!! then maps them into the FS kpt states
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (JPCroc) based on conducti
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nkpttot,nbandtot
!!
!! OUTPUT
!! el_veloc(nkpttot,nbandtot,3)
!!  
!! PARENTS
!!      get_gkk_qpt_tr
!!
!! CHILDREN
!!      hdr_clean,hdr_io,hdr_skip,wffclose,wffopen,wffreadeigk
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine read_el_veloc(mpi_enreg,nbandtot,nkpttot,nsppol_in,elph_tr_ds)

 use defs_basis
  use defs_datatypes
  use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_13io_mpi
 use interfaces_14iowfdenpot
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: nbandtot,nkpttot,nsppol_in
 type(MPI_type) :: mpi_enreg
 type(elph_tr_type) :: elph_tr_ds

!Local variables-------------------------------
!scalars
 integer :: accesswff,bantot,bd2tot_index,bdtot0_index,bdtot_index,dosdeltae
 integer :: fform0,formeig,formeig0,headform,iatom,iband,ierr,ii,ikpt,index_1
 integer :: io_status,iom,isppol,jband,jj,l1,l2,lnfm,master,mband,me,mom,mu
 integer :: natom,nband1,nband_k,nkpt,nlign,npw1,nrest,nrot,nspinor,nspinor1
 integer :: nsppol,ntypat,nu,occopt,rdwr,spaceComm,tim_rwwf
 real(dp) :: Tatm,deltae,diff_occ,dom,ecut,entropy,etotal,fermie,maxocc,nelect
 real(dp) :: np_sum,np_sum_k1,np_sum_k2,omin,oml,residm,sig,socc,swtk,tphysel
 real(dp) :: tsmear,ucvol,usnelv,wind
 character(len=fnlen) :: filnam,filnam0,filnam1,filnam2,filnam3,filnam_out
 type(hdr_type) :: hdr
 type(wffile_type) :: wff0,wff1,wff2,wff3
!arrays
 integer,allocatable :: nband(:)
 real(dp) :: gmet(3,3),gmet_inv(3,3),gprimd(3,3),gprimd_inv(3,3),im_el_veloc(3)
 real(dp) :: rmet(3,3),rmet_inv(3,3),rprimd(3,3)
 real(dp),allocatable :: doccde(:),doccde_k(:),eig0_k(:),eig0tmp(:),eig1_k(:,:)
 real(dp),allocatable :: eigen0(:),eigen11(:),eigen12(:),eigen13(:),eigtmp(:)
 real(dp),allocatable :: occ(:),wtk(:)

! *********************************************************************************
!BEGIN EXECUTABLE SECTION

!Read data file name
!TODO: this should be standardized and read in anaddb always, not
!conditionally. Otherwise when new files are added to the anaddb files
!file...  Catastrophe!

 write(6,*)'enter read_el_veloc '

!write(6,'(a)')' Please give the name of the data file for transport properties...'
!read(5, '(a)')filnam
!write(6,'(a)')' The name of the data file is :',filnam

!Read data file
 open(15,file=trim(elph_tr_ds%ddkfilename),form='formatted')
 rewind(15)
 read(15,'(a)')filnam1       ! first ddk file
 read(15,'(a)')filnam2       ! second ddk file
 read(15,'(a)')filnam3       ! third ddk file
 read(15,'(a)')filnam0       ! ground-state data

!Open the Wavefunction files
!These default values are typical of sequential use
!filnam0 is the GS data file it should be duplicated to avoid nasty cross readings with the main anaddb
!filnam1,2,3 are the ddk files
 accesswff=0 ; spaceComm=abinit_comm_serial ; master=0 ; me=0
 call WffOpen(accesswff,spaceComm,filnam0,ierr,wff0,master,me,10)
 call WffOpen(accesswff,spaceComm,filnam1,ierr,wff1,master,me,11)
 call WffOpen(accesswff,spaceComm,filnam2,ierr,wff2,master,me,12)
 call WffOpen(accesswff,spaceComm,filnam3,ierr,wff3,master,me,13)

!Read the header from the first ddk file (might have been the GS file ?)
 rdwr=1
 call hdr_io(fform0,hdr,rdwr,wff1)

!Extract info from the header
 headform=hdr%headform
 bantot=hdr%bantot
 ecut=hdr%ecut_eff
 natom=hdr%natom
 nkpt=hdr%nkpt
 if(nkpt.ne.nkpttot) then
  write(6,*)'read_el_veloc ******** wrong number of kpoints',nkpt,nkpttot
  write(6,*)' istwfk = ', hdr%istwfk
  write(6,*)' kpt sets in ground state and ddk files must agree.'
  write(6,*)' Action: use kptopt 3 in all calculations,'
  write(6,*)' or code the completion here in read_el_veloc'
  stop
 end if
 nspinor=hdr%nspinor
 nsppol=hdr%nsppol
 if(nsppol.ne.nsppol_in) then
  write(6,*)'read_el_veloc ******** nsspol<>input nsppol'
  stop
 end if

 ntypat=hdr%ntypat
 occopt=hdr%occopt
 rprimd(:,:)=hdr%rprimd(:,:)
 allocate(nband(nkpt*nsppol),occ(bantot))
 fermie=hdr%fermie
 write (*,*) 'read_el_veloc : Warning using Fermi energy from GS,',&
& ' ignoring anaddb input'
 occ(1:bantot)=hdr%occ(1:bantot)
 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)

!Get mband, as the maximum value of nband(nkpt)
 mband=maxval(nband(:))

 write(6,*)
 write(6,*)'readings from read_kpt_veloc header'
!write(6,'(a,3f10.5,a)' )' rprimd(bohr)      =',rprimd(1,1:3)
!write(6,'(a,3f10.5,a)' )'                    ',rprimd(2,1:3)
!write(6,'(a,3f10.5,a)' )'                    ',rprimd(3,1:3)
 write(6,'(a,i8)')       ' natom             =',natom
 write(6,'(a,2i8)')      ' nkpt,mband        =',nkpt,mband
 write(6, '(a, f10.5,a)' ) ' ecut              =',ecut,' Ha'
 write(6,'(a,f10.5,a,f10.5,a)' )' fermie            =',fermie,' Ha',fermie*Ha_eV,' eV'

!Prepare the reading of ddk Wff files
 formeig0=0 ; formeig=1 ; tim_rwwf=0
 allocate(eigtmp(2*mband*mband),eig0tmp(mband))
 call hdr_skip(wff0,ierr)
 call hdr_skip(wff2,ierr)
 call hdr_skip(wff3,ierr)

!Read the eigenvalues of ground-state and ddk files
 allocate(eigen0(mband*nkpt*nsppol))
 allocate(eigen11(2*mband*mband*nkpt*nsppol))
 allocate(eigen12(2*mband*mband*nkpt*nsppol))
 allocate(eigen13(2*mband*mband*nkpt*nsppol))
 bdtot0_index=0 ; bdtot_index=0

 do isppol=1,nsppol
  do ikpt=1,nkpt
!  write(6,*)'ikpt bdtot0_index bdtot_index',ikpt, bdtot0_index, bdtot_index


   nband1=nband(ikpt+(isppol-1)*nkpt)
   call WffReadEigK(eig0tmp,formeig0,headform,ikpt,isppol,mband,mpi_enreg,nband1,tim_rwwf,wff0)
   eigen0(1+bdtot0_index:nband1+bdtot0_index)=eig0tmp(1:nband1)
   call WffReadEigK(eigtmp,formeig,headform,ikpt,isppol,mband,mpi_enreg,nband1,tim_rwwf,wff1)
   eigen11(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)
   call WffReadEigK(eigtmp,formeig,headform,ikpt,isppol,mband,mpi_enreg,nband1,tim_rwwf,wff2)
   eigen12(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)
   call WffReadEigK(eigtmp,formeig,headform,ikpt,isppol,mband,mpi_enreg,nband1,tim_rwwf,wff3)
   eigen13(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)
   bdtot0_index=bdtot0_index+nband1
   bdtot_index=bdtot_index+2*nband1**2
  end do
 end do ! end isppol

 call WffClose(wff0,ierr)
 call WffClose(wff1,ierr)
 call WffClose(wff2,ierr)
 call WffClose(wff3,ierr)

 deallocate(eigtmp,eig0tmp)

 allocate(eig0_k(nbandtot),eig1_k(2*nbandtot**2,3)) !
 bdtot_index = 0
 bd2tot_index = 0
 elph_tr_ds%el_veloc(:,:,:,:)=zero

 do isppol=1,nsppol
  im_el_veloc(:)=0.
  do ikpt=1,nkpt
!  write(10,*)'ikpt bdtot0_index bdtot_index',ikpt, bdtot0_index, bdtot_index
!  write(6,*)'ikpt',ikpt
   nband_k=nband(ikpt)

   if(nband_k.ne.nbandtot) then
    write(6,*)'read_el_veloc ******** change in nband number for ikpt ',ikpt,nband_k,nbandtot
    stop
   end if


!  


!  eigenvalue for k-point
   eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
!  first derivative eigenvalues for k-point
   eig1_k(:,1)=eigen11(1+bd2tot_index:2*nband_k**2+bd2tot_index)
   eig1_k(:,2)=eigen12(1+bd2tot_index:2*nband_k**2+bd2tot_index)
   eig1_k(:,3)=eigen13(1+bd2tot_index:2*nband_k**2+bd2tot_index)

!  
!  LOOP OVER BAND
!  write(6,*)'rprimd'
!  write(6,*)rprimd


   do iband=1,nband_k
    do l1=1,3
     do ii=1,3
      elph_tr_ds%el_veloc(ikpt,iband,l1,isppol)=elph_tr_ds%el_veloc(ikpt,iband,l1,isppol)+&
&      rprimd(l1,ii)*eig1_k(2*iband-1+(iband-1)*2*nband_k,ii)/two_pi
      im_el_veloc(l1)=im_el_veloc(l1)+&
&      rprimd(l1,ii)*eig1_k(2*iband+(iband-1)*2*nband_k,ii)/two_pi
     end do
    end do

!   crc veloc reelle et imaginaire ???
!   if(any(abs(im_el_veloc(:)).ge.1e-25)) then
!   write(10,*)
!   write(10,*)'ikpt,iband',ikpt,iband
!   write(10,'(I3,I3,2D13.5)')ikpt,iband,elph_tr_ds%el_veloc(ikpt,iband,1),im_el_veloc(1)
!   write(10,'(I3,I3,2D13.5)')ikpt,iband,elph_tr_ds%el_veloc(ikpt,iband,2),im_el_veloc(2)
!   write(10,'(I3,I3,2D13.5)')ikpt,iband,elph_tr_ds%el_veloc(ikpt,iband,3),im_el_veloc(3)
!   usnelv=1./sqrt((elph_tr_ds%el_veloc(ikpt,iband,1)**2+elph_tr_ds%el_veloc(ikpt,iband,2)**2+elph_tr_ds%el_veloc(ikpt,iband,3)**2))

!   end if
   end do
   bd2tot_index=bd2tot_index+2*nband_k**2
   bdtot_index=bdtot_index+nband_k
  end do
 end do ! end isppol


!end do


 deallocate(nband)
 deallocate(occ)
 deallocate(eigen11,eigen12,eigen13)
 deallocate(eigen0)
 write(6,*)'out of read_el_veloc '


 call hdr_clean(hdr)

end subroutine read_el_veloc
!!***
