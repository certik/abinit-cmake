!{\src2tex{textfont=tt}}
!!****f* ABINIT/wannier
!! NAME
!! wannier
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations of the wannier functions
!!
!! COPYRIGHT
!! Copyright (C) 2006-2008 ABINIT group (XG,JBattacharya)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  mband =maximum number of bands
!!  mpi_enreg=informations about MPI parallelization
!!  nkpt  =number of k points
!!  nsppol=number of channels for spin-polarization (1 or 2)
!!
!! OUTPUT
!!  (no direct output : results written)
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      distrb2,hdr_clean,hdr_io,hdr_io_netcdf,metric,status,wffclose,wffopen
!!      wrtout,xcomm_world,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wannier(dtfil,dtset,iexit,mband,mpi_enreg,nkpt,nsppol)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_14iowfdenpot
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iexit,mband,nkpt,nsppol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset

!Local variables-------------------------------
!scalars
 integer,parameter :: level=3
 integer :: accessfil,exchn2n3d,fform0,iatom,ierr,ii,insmet,master,me,natom,nr1
 integer :: nr2,nr3,nspden,ntypat,prtvol,rdwr,spaceworld
 real(dp) :: efermi,ucvol
 character(len=500) :: message
 character(len=fnlen) :: wffnm
 type(hdr_type) :: hdr
 type(wffile_type) :: wff1
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: xcart(:,:),xred(:,:)

!***********************************************************************

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 prtvol=dtset%prtvol
 if(prtvol==-level)then
  write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,&
&  ' wannier : enter , debug mode '
  call wrtout(06,message,'COLL')
 end if

!------------------------------------------------------------------------
!Init parallel computation ... if needed later ...
 mpi_enreg%paralbd=0
 mpi_enreg%me_fft=0
 mpi_enreg%nproc_fft=1
 mpi_enreg%paral_fft=0
 mpi_enreg%paral_level=2

!If dtset%accesswff == 2 set all array outputs to netcdf format
 accessfil = 0
 if (dtset%accesswff == 2) accessfil = 1
 if (dtset%accesswff == 3) accessfil = 3

 master=0
!Init me
 call xme_init(mpi_enreg,me)
!BEGIN TF_CHANGES
 call xcomm_world(mpi_enreg,spaceworld)
!END TF_CHANGES

 if(mpi_enreg%paral_compil_kpt==1)then
  allocate(mpi_enreg%proc_distrb(nkpt,mband,nsppol))
  mpi_enreg%parareel=0
  call distrb2(mband, dtset%nband, nkpt, nsppol, mpi_enreg)
 end if

!Parallelism has been initialized
!------------------------------------------------------------------------
!Initialize WF input files

 wffnm=trim(dtfil%fnamewffk)
!Open the file in which the input wfs are stored, also initialize wff1 handle
 call WffOpen(dtset%accesswff,spaceworld,wffnm,ierr,wff1,master,me,dtfil%unwff1)

 write(message, '(a,a)' )&
& ' wannier : will read wavefunctions from disk file ',trim(wffnm)
 call wrtout(ab_out,message,'COLL')
 write(message, '(a)' )&
& ' wannier : calculation will be mostly based on the data in the header of this file '
 call wrtout(ab_out,message,'COLL')

!Read header
 rdwr=1
 if (dtset%accesswff /= 2) then
  call hdr_io(fform0,hdr,rdwr,wff1)
#if defined HAVE_NETCDF
 else if (dtset%accesswff == 2) then
  call hdr_io_netcdf(fform0,hdr,rdwr,wff1)
 else if (dtset%accesswff == 2) then
  write (std_out,*) "FIXME: ETSF I/O support in wannier"
#endif
 end if

!WF input file has been initialized
!------------------------------------------------------------------------
!Echo part of the header
 rdwr=4
 call hdr_io(fform0,hdr,rdwr,ab_out)

!-------------------------------------------------------------------------
!Data initialized from the header
 natom=hdr%natom
 nr1=hdr%ngfft(1)
 nr2=hdr%ngfft(2)
 nr3=hdr%ngfft(3)
 nspden=hdr%nspden
 ntypat=hdr%ntypat
 rprimd(:,:)=hdr%rprimd(:,:)

!Need to know natom in order to allocate xcart
 allocate(xcart(3,natom),xred(3,natom))
 xred(:,:)=hdr%xred(:,:)
 do iatom=1,natom
  xcart(:,iatom)=rprimd(:,1)*xred(1,iatom)+&
&  rprimd(:,2)*xred(2,iatom)+&
&  rprimd(:,3)*xred(3,iatom)
 end do

!Echo the value of different input parameters
 write(*,*)'ECHO important input variables ...'
 write(*,*)
 write(*,*) ' Dimensional primitive vectors (ABINIT equivalent : rprimd):'
 write(*, '(3es16.6)' ) rprimd(1:3,1)
 write(*, '(3es16.6)' ) rprimd(1:3,2)
 write(*, '(3es16.6)' ) rprimd(1:3,3)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 if (abs(ucvol)<1.0d-12) then
  write(*,*)' At least two rprimd(,) vectors are collinear'
  write(*,*)' Please check the input rprim and acell, or rprimd.'
  write(*,*)' The program will stop'
  stop
 end if

 write(*, '(a,3i5)' ) '  Grid density (ABINIT equivalent : ngfft): ',nr1,nr2,nr3
 write(*,*) ' Number of atoms       :',natom
 write(*,*) ' Number of atomic types:',ntypat

 write(*,*)
 write(*,*) '  #    Atomic positions (cartesian coordinates - Bohr)'
 do iatom=1,natom
  write(6, '(i4,3es16.6)' )iatom,xcart(1:3,iatom)
 end do
 write(*,*)

 exchn2n3d=0
 efermi=hdr%fermie
 insmet=1
 if(hdr%occopt>=3)insmet=2

!XG060209 : Joydeep, I have prepared the call to localorb_S here. However,
!chr_inputfname should not be used anymore, since all the input variables
!are now defined in the usual ABINIT input file, as input variables ...
!Please, see all the new input variables defined in doc/input_variables/varwan.html,
!that you can also access from the Web.
!They are defined in the dtset structured datatype. One way to go
!would be to extract them here from dtset, as follows :
!ncenter=dtset%ncenter
!norb=dtset%norb
!...
!and them to pass all of them as arguments of localorb_S .
!Another thing that might be done would be to pass the full dtset as argument localorb_S
!and to extract them inside the routine.
!I leave this choice to you ...
!Note that efermi is already initialized, from hdr%fermie (see above),
!and insmet is initialized from occopt (see above)

!call localorb_S(chr_inputfname,hdr%ecut_eff,exchn2n3d,hdr%headform,hdr%istwfk,hdr%kptns,&
!& natom,hdr%nband,hdr%nkpt,hdr%npwarr,&
!& nr1,nr2,nr3,hdr%nspinor,hdr%nsppol,ntypat,rprimd,xcart,hdr%typat,hdr%znucltypat)

 write(*,*)" "
 write(*,*)" ###################################################### "
 write(*,*)" "
 write(*,*)" Localized orbital files fort.1**1 for spin up "
 write(*,*)"                     and fort.1**2 for spin dn written."
 write(*,*)" "
 write(*,*)" ###################################################### "


!-------------------------------------------------------------------------
!Cleaning

 deallocate(xcart,xred)

!Clean the file handle
 call WffClose(wff1,ierr)

!Clean the header
 call hdr_clean(hdr)

 if(mpi_enreg%paral_compil_kpt==1)then
  deallocate(mpi_enreg%proc_distrb)
 end if

 write(message, '(a,a)' ) ch10,' wannier : exiting '
 call wrtout(06,message,'COLL')

 call status(0,dtfil%filstat,iexit,level,'exit          ')

end subroutine wannier
!!***
