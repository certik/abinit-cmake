!{\src2tex{textfont=tt}}
!!****f* ABINIT/outscfcv
!! NAME
!! outscfcv
!!
!! FUNCTION
!! Output routine for the scfcv.F90 routine
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
!!  compch_fft=compensation charge, from FFT grid
!!  compch_sph=compensation charge, from sphere
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecrpj)=<p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                                 and each |p_lmn> non-local projector
!!  dimcprj(natom*usecprj)=array of dimensions of array cprj
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  etotal=total energy
!!  fermie= Fermi energy
!!  filapp character(len=fnlen)=generic output root name, with appendix
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kssform=govern the Kohn-Sham Structure file format
!!  mband=maximum number of bands
!!  mgfftc=maximum size of 1D FFTs for the PAW coarse grid
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw.
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nfft=(effective) number of FFT grid points (for this processor) (see NOTES at beginning of scfcv)
!!  ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!  nhat(nfft,nspden*usepaw)= compensation charge density  (PAW)
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetries in space group
!!  ntypat=number of types of atoms in unit cell.
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  occ(mband*nkpt*nsppol)=occupation number for each band (usually 2) for each k.
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr(natom) <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  ph1dc(2,3*(2*mgfftc+1)*natom)=one-dimensional structure factor information
!!            note:structure factors are given on the coarse grid for PAW
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(nfft,nspden)=total electron density in electrons/bohr**3, reciprocal space.
!!  rhor(nfft,nspden)=total electron density in electrons/bohr**3, real space.
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  ucvol=unit cell volume (bohr**3)
!!  usecprj=1 if cprj datastructure has been allocated
!!  usexcnhat= flag controling use of compensation density in the computation of Vxc
!!  vhartr(nfft)=Hartree potential
!!  vxc(nfft,nspden)=xc potential
!!  vxcavg=vxc average
!!  wffnow=information about wf disk file
!!  vtrial(nfft,nspden)=the trial potential
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)=real spherical harmonics for each G and k point
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      abi_etsf_electrons_put,abi_etsf_geo_put,bonds_lgth_angles,calc_cs
!!      calc_efg,calc_fc,calcdensph,ioarr,leave_new,mati3inv,mlwfovlp
!!      optics_paw,out1dm,outkss,outwant,partial_dos_fractions
!!      partial_dos_fractions_paw,pawmknhat,pawprt,printbxsf,prt_cml2,rhohxc
!!      tetrahedron,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine outscfcv(atindx,atindx1,cg,compch_fft,compch_sph,cprj,dimcprj,dtfil,dtset,ecut,eigen,etotal,&
& fermie,filapp,gmet,gprimd,gsqcut,hdr,kg,&
& kssform,mband,mgfftc,mkmem,mpi_enreg,mpsang,mpw,natom,&
& nattyp,nfft,ngfft,nhat,nkpt,npwarr,nspden,nspinor,nsppol,nsym,ntypat,n3xccc,occ,&
& pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,paw_ij,ph1dc,prtvol,psps,rhog,rhor,rmet,rprimd,&
& ucvol,usecprj,usexcnhat,wffnow,vhartr,vtrial,vxc,vxcavg,xccc3d,xred,ylm)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13ionetcdf
 use interfaces_13paw
 use interfaces_13xc
 use interfaces_13xml
 use interfaces_14iowfdenpot
 use interfaces_14occeig
 use interfaces_15common
 use interfaces_18seqpar
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kssform,mband,mgfftc,mkmem,mpsang,mpw,n3xccc,natom,nfft
 integer,intent(in) :: nkpt,nspden,nsppol,nsym,ntypat,prtvol,usecprj,usexcnhat
 integer,intent(inout) :: nspinor
 real(dp),intent(in) :: compch_fft,compch_sph,ecut,fermie,gsqcut,ucvol
 real(dp),intent(inout) :: etotal,vxcavg
 character(len=fnlen),intent(inout) :: filapp
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),dimcprj(natom*usecprj)
 integer,intent(in) :: kg(3,mpw*mkmem),nattyp(ntypat),ngfft(18),npwarr(nkpt)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1dc(2,3*(2*mgfftc+1)*natom)
 real(dp),intent(in) :: rmet(3,3),rprimd(3,3),vhartr(nfft),xccc3d(n3xccc)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(inout) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(inout) :: nhat(nfft,nspden*psps%usepaw),rhog(nfft,nspden)
 real(dp),intent(inout) :: rhor(nfft,nspden),vtrial(nfft,nspden)
 real(dp),intent(inout) :: vxc(nfft,nspden),xred(3,natom)
 type(cprj_type),intent(in) :: cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)
 type(paw_ij_type),intent(inout) :: paw_ij(natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: accessfil,coordn,fformr,fformv,i1,i2,i3,iat,ierr,ifft,ii,ii2,ii3
 integer :: ikpt,ilmn,ispden,isppol,isym,itypat,m_dos_flag,mbesslang
 integer :: ndosfraction,nfft_tmp,nhatgrdim,occopt,optxc,partial_dos_flag
 integer :: paw_dos_flag,prt1dm,prtcml,prtcs,prtden,prtdos,prtefg,prtfc,prtgeo
 integer :: prtnabla,prtpot,prtstm,prtvha,prtvhxc,prtvxc,rdwr,rdwrpaw,timrev
 real(dp) :: dum,enxc_dum
 character(len=500) :: message
 character(len=fnlen) :: bxsfname,fildata
!arrays
 integer :: ngfft_tmp(18)
 integer,allocatable :: npwarr1(:),symrec(:,:,:)
 real(dp) :: strsxc_dum(6),tsec(2)
 real(dp),allocatable :: dos_fractions(:,:,:,:),dos_fractions_m(:,:,:,:)
 real(dp),allocatable :: dos_fractions_paw1(:,:,:,:)
 real(dp),allocatable :: dos_fractions_pawt1(:,:,:,:),eigen2bxsf(:,:,:)
 real(dp),allocatable :: kxc_dum(:,:),nhatgr(:,:,:),rhocorval(:,:),rhor_tmp(:)
 real(dp),allocatable :: vhartr_tmp(:),vwork(:,:),vxc_dum(:,:)
 type(pawrhoij_type),allocatable :: pawrhoij_dum(:)

! *************************************************************************

!DEBUG
!write(6,*)' outscfcv : enter '
!ENDDEBUG

 if (usecprj==0.and.psps%usepaw==1.and. &
& (dtset%prtwant==2.or.dtset%prtnabla>0.or.dtset%prtdos==3.or.dtset%kssform==3)) then
  write (message,'(8a)')ch10,&
&  ' outscfcv : ERROR- ',ch10,&
&  ' cprj datastructure must be allocated',ch10,&
&  ' with options prtwant=2, prtnabla>0 or prtdos>3 !',ch10,&
&  ' Action: change pawusecp input keyword.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!begin wannier
 if (dtset%prtwant==2) then
  call mlwfovlp(cg,cprj,dtset,ecut,eigen,fermie,gprimd,gsqcut,kg,&
&  mband,mgfftc,mkmem,mpi_enreg,mpsang,mpw,natom,&
&  nattyp,nfft,ngfft,nkpt,npwarr,nspden,nspinor,nsppol,ntypat,&
&  pawang,pawrad,pawtab,prtvol,psps,rprimd,ucvol,xred)
 end if
!end wannier

!
!if accesswff == 2 then set all outputs to netcdf format
!if accesswff == 3 then set all outputs to ETSF format
!
 accessfil = 0
 if (dtset%accesswff == 2) accessfil = 1
 if (dtset%accesswff == 3) accessfil = 3
 if (dtset%accesswff == 1) accessfil = 4

 occopt=dtset%occopt;

 prtden=dtset%prtden ; prtpot=dtset%prtpot ; prtgeo=dtset%prtgeo
 prtcml=dtset%prtcml ; prtdos=dtset%prtdos ; prtstm=dtset%prtstm
 prt1dm=dtset%prt1dm ; prtvha=dtset%prtvha ; prtvhxc=dtset%prtvhxc
 prtvxc=dtset%prtvxc ; prtnabla=dtset%prtnabla; prtefg=dtset%prtefg
 prtcs=dtset%prtcs   ; prtfc=dtset%prtfc

!Warnings :
!- core charge is excluded from the charge density;
!- the potential is the INPUT vtrial.
 if(  mpi_enreg%paral_compil_kpt==0                         .or. &
& (mpi_enreg%me==0 .and. mpi_enreg%parareel == 0 .and. mpi_enreg%paral_compil_fft==0 ) .or. &
& (mpi_enreg%paral_compil_fft==1 .and. mpi_enreg%me_band==0 .and. mpi_enreg%me_kpt==0) .or. &
& (mpi_enreg%me_group_para==0 .and. mpi_enreg%parareel == 1)) then

  if (prtden/=0 .or. prtpot>0 .or. prtgeo>0 .or. prt1dm>0 .or. prtstm>0 .or. &
&  prtvha>0 .or. prtvhxc>0.or. prtvxc>0 .or. prtcml>0 .or. prtdos>=2 )then
  end if

! We output the density.
  if (prtden/=0) then
   rdwr=2 ; fformr=52 ; rdwrpaw=0
!  We create the file name.
   fildata=trim(filapp)//'_DEN'
!  We call ioarr with the adequat rhor values, depending on positron.
   if (dtset%positron==3) then
    allocate(rhocorval(nfft,nspden))
    do ifft=1,nfft
     rhocorval(ifft,1)=rhor(ifft,1)+xccc3d(ifft)
    end do
    if(nspden==2) then
     do ifft=1,nfft
      rhocorval(ifft,2)=rhor(ifft,2)+0.5_dp*xccc3d(ifft)
     end do
    end if
    call ioarr(accessfil,rhocorval, dtset, etotal,fformr,fildata,hdr, mpi_enreg, &
&    nfft,pawrhoij_dum,rdwr,rdwrpaw,ngfft)
    deallocate(rhocorval)
   else if (dtset%positron==1) then
    call ioarr(accessfil,rhor, dtset, etotal,fformr,fildata,hdr, mpi_enreg, &
&    nfft,pawrhoij_dum,rdwr,rdwrpaw,ngfft)
   else
    call ioarr(accessfil, rhor, dtset, etotal, fformr, fildata, hdr, mpi_enreg, &
&    nfft, pawrhoij_dum, rdwr, rdwrpaw, ngfft)
   end if
   if ( accessfil == 3 ) then
!   Complete the geometry informations with missing values from hdr_io().
    call abi_etsf_geo_put(dtset, fildata, psps, rprimd, xred)
!   Complete the electrons definition with missing values from hdr_io().
    call abi_etsf_electrons_put(dtset, fildata)
   end if
  end if

! We handle the output of wavefunctions.
  if (dtset%prtwf == 1) then
!  We create the file name.
   fildata = trim(filapp)//'_WFK'
!  In ETSF, some geometric informations are required for wave functions files.
   if ( accessfil == 3 ) then
!   Complete the geometry informations with missing values from hdr_io().
    call abi_etsf_geo_put(dtset, fildata, psps, rprimd, xred)
!   Complete the electrons definition with missing values from hdr_io().
    call abi_etsf_electrons_put(dtset, fildata)
   end if
  end if

  if (prtpot>0) then
   fildata=trim(filapp)//'_POT'
   rdwr=2 ; fformv=102 ; rdwrpaw=0
!  
!  MJV note: why is accessfil forced to 0???? This disables the writing of ETSF
!  format potentials!
!  
!  set to 1 for netcdf output
   accessfil = 0
   call ioarr(accessfil,vtrial, dtset, etotal,fformv,fildata,hdr, mpi_enreg, &
&   nfft,pawrhoij_dum,rdwr,rdwrpaw,ngfft)
   if ( accessfil == 3 ) then
!   Complete the geometry informations with missing values from hdr_io().
    call abi_etsf_geo_put(dtset, fildata, psps, rprimd, xred)
!   Complete the electrons definition with missing values from hdr_io().
    call abi_etsf_electrons_put(dtset, fildata)
   end if
  end if

  if (prtgeo>0) then
   coordn=prtgeo
   if ( accessfil == 3 ) then
    call abi_etsf_geo_put(dtset, filapp, psps, rprimd, xred)
   else
    call bonds_lgth_angles(coordn,filapp,natom,psps%ntypat,&
&    rprimd,dtset%typat,xred,dtset%znucl)
   end if
  end if

  if (prtcml>0) then
   call prt_cml2(filapp,natom,dtset%nsym,psps%ntypat,&
&   rprimd,dtset%spgroup,dtset%symrel,dtset%tnons,dtset%typat,xred,dtset%znucl)
  end if

  if (prtstm>0) then
   rdwr=2 ; fformr=52 ; rdwrpaw=0
!  set to 1 for netcdf output
   accessfil = 0
   fildata=trim(filapp)//'_STM'
   call ioarr(accessfil,rhor, dtset, etotal,fformr,fildata,hdr, mpi_enreg, &
&   nfft,pawrhoij_dum,rdwr,rdwrpaw,ngfft)
   if ( accessfil == 3 ) then
!   Complete the geometry informations with missing values from hdr_io().
    call abi_etsf_geo_put(dtset, fildata, psps, rprimd, xred)
!   Complete the electrons definition with missing values from hdr_io().
    call abi_etsf_electrons_put(dtset, fildata)
   end if
  end if

  if (prt1dm>0) then
   call out1dm(filapp,natom,nfft,ngfft,nspden,psps%ntypat,&
&   rhor,rprimd,dtset%typat,ucvol,vtrial,xred,dtset%znucl)
  end if

  if (prtvha>0) then
   fildata=trim(filapp)//'_VHA'
   rdwr=2 ; fformv=102 ; rdwrpaw=0
!  set to 1 for netcdf output
   allocate(vwork(nfft,nspden))
!  In the SCF part of the positron, we have to compute Vhartree acting on electrons
   if (dtset%positron==1.or.dtset%positron==2)  then
    allocate(vhartr_tmp(nfft),kxc_dum(nfft,0),vxc_dum(nfft,nspden))
    optxc = 1 ; nhatgrdim=0
    if (psps%usepaw==1.and.usexcnhat>0.and.dtset%xclevel==2) then
     nhatgrdim=1;allocate(nhatgr(nfft,nspden,3))
     call pawmknhat(dum,1,0,mpi_enreg,natom,nfft,ngfft,nhatgrdim,nspden,ntypat,&
&     dtset%paral_kgb,pawang,pawfgrtab,nhatgr,nhat,pawrhoij,pawtab,dtset%typat,ucvol)
    end if
    call rhohxc(dtset,enxc_dum,gsqcut,psps%usepaw,kxc_dum,mpi_enreg,nfft,ngfft,nhat,1,&
&    nhatgr,nhatgrdim,0,nspden,n3xccc,optxc,rhog,rhor,rprimd,strsxc_dum,usexcnhat,vhartr_tmp,vxc_dum,vxcavg,xccc3d)
    if (nhatgrdim>0) deallocate(nhatgr)
    do ispden=1,nspden
     vwork(:,ispden)=vhartr_tmp(:)
    end do
    deallocate(vhartr_tmp,kxc_dum,vxc_dum)
   else
    do ispden=1,nspden
     vwork(:,ispden)=vhartr(:)
    end do
   end if
   call ioarr(accessfil,vwork, dtset, etotal,fformv,fildata,hdr, mpi_enreg, &
&   nfft,pawrhoij_dum,rdwr,rdwrpaw,ngfft)
   deallocate(vwork)
   if ( accessfil == 3 ) then
!   Complete the geometry informations with missing values from hdr_io().
    call abi_etsf_geo_put(dtset, fildata, psps, rprimd, xred)
!   Complete the electrons definition with missing values from hdr_io().
    call abi_etsf_electrons_put(dtset, fildata)
   end if
  end if

  if (prtvhxc>0) then
   fildata=trim(filapp)//'_VHXC'
   rdwr=2 ; fformv=102 ; rdwrpaw=0
!  set to 1 for netcdf output
   allocate(vwork(nfft,nspden))
   do ispden=1,nspden
    vwork(:,ispden)=vhartr(:)+vxc(:,ispden)
   end do
   call ioarr(accessfil,vwork, dtset, etotal,fformv,fildata,hdr, mpi_enreg, &
&   nfft,pawrhoij_dum,rdwr,rdwrpaw,ngfft)
   deallocate(vwork)
   if ( accessfil == 3 ) then
!   Complete the geometry informations with missing values from hdr_io().
    call abi_etsf_geo_put(dtset, fildata, psps, rprimd, xred)
!   Complete the electrons definition with missing values from hdr_io().
    call abi_etsf_electrons_put(dtset, fildata)
   end if
  end if

  if (prtvxc>0) then
   fildata=trim(filapp)//'_VXC'
   rdwr=2 ; fformv=102 ; rdwrpaw=0
!  set to 1 for netcdf output
   call ioarr(accessfil,vxc, dtset, etotal,fformv,fildata,hdr, mpi_enreg, &
&   nfft,pawrhoij_dum,rdwr,rdwrpaw,ngfft)
   if ( accessfil == 3 ) then
!   Complete the geometry informations with missing values from hdr_io().
    call abi_etsf_geo_put(dtset, fildata, psps, rprimd, xred)
!   Complete the electrons definition with missing values from hdr_io().
    call abi_etsf_electrons_put(dtset, fildata)
   end if
  end if

 end if ! if master

!Generate DOS using the tetrahedron method
 if (prtdos>=2) then

  if(prtdos==2)partial_dos_flag = 0
  if(prtdos==3)partial_dos_flag = 1
  m_dos_flag=0
  if (partial_dos_flag==1) m_dos_flag=dtset%prtdosm
  paw_dos_flag=0
  if (psps%usepaw==1.and.partial_dos_flag==1.and.dtset%pawprtdos>=1) paw_dos_flag=1

! mjv : initialization is needed as mbesslang is used for allocation below
  mbesslang = 1
  if(partial_dos_flag==1)then
   mbesslang = 5
   ndosfraction=dtset%natsph*mbesslang
  else
   ndosfraction = 1
   mbesslang = 0
  end if

! For other types of partial DOSs, should use a pointer or something
! to be able to allocate dos_fractions inside partial_dos_fractions. XG20030506 : Mmmm... not sure !
  allocate(dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction))
  if (m_dos_flag==1) &
&  allocate(dos_fractions_m(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*mbesslang))
  if (psps%usepaw==1.and.partial_dos_flag==1) &
&  allocate(dos_fractions_paw1(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction),&
&  dos_fractions_pawt1(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction))
  if( partial_dos_flag==1)then
!  Generate fractions for partial DOSs if needed
!  partial_dos 1,2,3,4  give different decompositions
   if (psps%usepaw==0.or.dtset%pawprtdos/=2) then
    call partial_dos_fractions(cg,dos_fractions,dos_fractions_m,dtfil,dtset,hdr,mbesslang,mpi_enreg, &
&    m_dos_flag,ndosfraction,partial_dos_flag,wffnow)
   else
    dos_fractions=zero;if (m_dos_flag==1) dos_fractions_m=zero
   end if
   if (psps%usepaw==1) then
    call partial_dos_fractions_paw(atindx1,cprj,dimcprj,dos_fractions,dos_fractions_m,&
&    dos_fractions_paw1,dos_fractions_pawt1,dtfil,dtset,psps%indlmn,&
&    psps%lmnmax,mbesslang,mkmem,mpi_enreg,m_dos_flag,ndosfraction,&
&    paw_dos_flag,pawrad,pawtab,prtdos)
   end if
  else
   dos_fractions(:,:,:,1)=one
  end if

! Here, computation and output of DOS and partial DOS
  fildata=trim(filapp)//'_DOS'
  call tetrahedron (dos_fractions,dos_fractions_m,dos_fractions_paw1,dos_fractions_pawt1,&
&  dtset,fermie,eigen,fildata,mbesslang,m_dos_flag,ndosfraction,paw_dos_flag,rprimd)

  deallocate(dos_fractions)
  if (m_dos_flag==1) deallocate(dos_fractions_m)
  if (psps%usepaw==1.and.partial_dos_flag==1) deallocate(dos_fractions_paw1,dos_fractions_pawt1)

 end if ! prtdos > 1

!Output of integrated density inside atomic spheres
 if (dtset%prtdensph==1) call calcdensph(gmet,mpi_enreg,natom,nfft,ngfft,nspden,&
& ntypat,dtset%ratsph,rhor,rprimd,dtset%typat,ucvol,xred)

!If PAW, provide additional outputs
 if (psps%usepaw==1) then
! Output of compensation charge
  if (usexcnhat>0) then    !if (dtset%nstep>0.or.dtfil%ireadwf/=0)
   write(message, '(4a)' )ch10,' PAW TEST:',ch10,&
&   ' ==== Compensation charge inside spheres ============'
   if (compch_sph>-1.d4.and.compch_fft>-1.d4) &
&   write(message, '(3a)' ) trim(message),ch10,&
   ' The following values must be close to each other ...'
   if (compch_sph>-1.d4) write(message, '(3a,f22.15)' ) trim(message),ch10,&
&   ' Compensation charge over spherical meshes = ',compch_sph
   if (compch_fft>-1.d4) then
    if (pawfgr%usefinegrid==1) then
     write(message, '(3a,f22.15)' ) trim(message),ch10,&
&     ' Compensation charge over fine fft grid    = ',compch_fft
    else
     write(message, '(3a,f22.15)' ) trim(message),ch10,&
&     ' Compensation charge over fft grid         = ',compch_fft
    end if
   end if
   call wrtout(ab_out,message,'COLL')
   call wrtout(6,message,'COLL')
  end if
! Output of pseudopotential strength Dij and augmentation occupancies Rhoij
  call pawprt(psps%indlmn,dtset%enunit,psps%lmnmax,natom,ntypat,paw_ij,&
&  dtset%pawprtvol,pawrhoij,pawtab,dtset%typat)
 end if

!PAW + output for optical conductivity
 if (psps%usepaw==1.and.prtnabla>0) then
  fildata=trim(filapp)//'_OPT'
  call optics_paw(atindx1,cg,cprj,dimcprj,dtfil,dtset,ecut,fildata,gprimd,hdr,psps%indlmn,kg,psps%lmnmax,&
&  mband,mkmem,mpi_enreg,mpsang,mpw,natom,nattyp,nkpt,npwarr,nspinor,nsppol,pawrad,pawtab,wffnow)
 end if

!Optionally provide output for GW part of ABINIT
 if (dtset%nbandkss/=0) then
  call timab(233,1,tsec)
  call outkss(atindx,atindx1,dtfil,dtset,ecut,gmet,gprimd,hdr,kg,&
&  dtset%kssform,mband,mgfftc,mkmem,mpi_enreg,mpsang,mpw,natom,&
&  nattyp,nfft,nkpt,npwarr,nspinor,nspden,nsppol,nsym,psps%ntypat,occ,pawrad,pawtab,&
&  ph1dc,prtvol,psps,rmet,rprimd,ucvol,wffnow,vtrial,xred,ylm,cg,usecprj,cprj,eigen)
  call timab(233,2,tsec)
 end if

!Optionally provide output for WanT
 if (dtset%prtwant==1) then
  call outwant(dtfil,dtset,eigen,cg,kg,npwarr,mband,mpi_enreg,nkpt,nsppol,&
&  nspinor,mkmem,mpw,wffnow,dtset%prtwant)
 end if

!Optionally provide output for chemical shielding calculation
 if (prtcs > 0) then
  call calc_cs(dtset%corecs,natom,nspden,ntypat,occopt,pawang,pawrad,pawrhoij,pawtab,&
&  prtcs,dtset%typat,psps%usepaw)
 end if

!Optionally provide output for electric field gradient calculation
 if (prtefg > 0) then
  call calc_efg(gprimd,natom,nfft,ngfft,nhat,nspden,ntypat,dtset%paral_kgb,pawang,pawrad,pawrhoij,pawtab,&
&  dtset%ptcharge,prtefg,dtset%quadmom,rhor,rprimd,dtset%typat,ucvol,psps%usepaw,xred,psps%zionpsp)
 end if

!Optionally provide output for Fermi-contact term at nuclear positions
 if (prtfc > 0) then
  call calc_fc(gprimd,natom,nfft,ngfft,nhat,nspden,ntypat,dtset%paral_kgb,&
&  pawrad,pawrhoij,pawtab,psps,rhor,rprimd,dtset%typat,xred)
 end if

!Optionally provide Xcrysden output for the Fermi surface
 if (dtset%prtfsurf ==1) then
! if (nspinor ==2 ) then !in case of nspinor==2 stop
! write (message,'(4a)')' outscf : ERROR- ',ch10,&
! &   ' nspinor == 2 not yet tested',ch10
! call wrtout(06,message,'COLL')
! call leave_new('COLL')
! end if

  allocate (eigen2bxsf(mband,nkpt,nsppol),stat=ierr)
  if (ierr /= 0 ) then
   write (message,'(5a)')' outscfcv : ERROR- ',ch10,&
&   ' trying to allocate array eigen2bxsf ' ,ch10,&
&   ' skipping the output of the Fermi surface'
   call wrtout(06,message,'COLL')
  else
   ii=0
   do isppol=1,nsppol
    do ikpt=1,nkpt
!    DEBUG
!    write(100,*)isppol,ikpt,hdr%kptns(:,ikpt)
!    write(100,*)eigen(1+ii:mband+ii)
!    ENDDEBUG
     eigen2bxsf(:,ikpt,isppol)= eigen(1+ii:mband+ii)
     ii=ii+mband
    end do
   end do

!  Invert symrels => gives symrels for kpoints
   allocate(symrec(3,3,nsym))
   do isym=1,nsym
    call mati3inv (dtset%symrel(:,:,isym),symrec(:,:,isym))
   end do

   timrev=1 !includes the time inversion symmetry
   fildata=trim(filapp) // '_BXSF'

   call printbxsf(eigen2bxsf,real(0.,dp),fermie,gprimd,dtset%kptrlatt,mband,dtset%nkpt,hdr%kptns,&
&   nsym,symrec,timrev,hdr%nsppol,dtset%shiftk,dtset%nshiftk,fildata)
   deallocate(symrec,eigen2bxsf)

  end if
 end if !if prtfsurf

!DEBUG
!write(6,*)' outscfcv : exit'
!stop
!ENDDEBUG

end subroutine outscfcv
!!***
