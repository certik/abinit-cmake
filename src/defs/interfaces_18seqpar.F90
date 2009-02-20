!!****m* ABINIT/interfaces_18seqpar
!! NAME
!! interfaces_18seqpar
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/18seqpar
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!! 
!!
!! SOURCE

module interfaces_18seqpar

 implicit none

interface
 subroutine berryphase_new(cg,cprj,dtefield,dtfil,dtset,&  
  &  gprimd,hdr,kg,mband,&  
  &  mkmem,mpi_enreg,mpw,natom,nattyp,npwarr,nspinor,nsppol,ntypat,&  
  &  nkpt,option,pawang,pawrad,pawtab,pel,pelev,pion,pwind,pwind_alloc,pwnsfac,&  
  &  rprimd,typat,ucvol,unit_out,usecprj,usepaw,wffnow,xred,zion)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, intent(in) :: mband
  integer, intent(in) :: mkmem
  integer, intent(in) :: mpw
  integer, intent(in) :: natom
  integer, intent(in) :: nkpt
  integer, intent(in) :: nspinor
  integer, intent(in) :: nsppol
  integer, intent(in) :: ntypat
  integer, intent(in) :: option
  integer, intent(in) :: pwind_alloc
  integer, intent(in) :: unit_out
  integer, intent(in) :: usecprj
  integer, intent(in) :: usepaw
  type(efield_type), intent(inout) :: dtefield
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(in) :: dtset
  type(hdr_type), intent(inout) :: hdr
  type(mpi_type), intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp), intent(in) :: ucvol
  type(wffile_type), intent(inout) :: wffnow
  real(dp), intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  type(cprj_type),intent(in) :: cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)
  real(dp), intent(in) :: gprimd(3,3)
  integer, intent(in) :: kg(3,mpw*mkmem)
  integer, intent(in) :: nattyp(ntypat)
  integer, intent(in) :: npwarr(nkpt)
  type(pawrad_type), intent(in) :: pawrad(ntypat*usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)
  real(dp), intent(out) :: pel(3)
  real(dp), intent(out) :: pelev(3)
  real(dp), intent(out) :: pion(3)
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(in) :: rprimd(3,3)
  integer, intent(in) :: typat(natom)
  real(dp), intent(inout) :: xred(3,natom)
  real(dp), intent(in) :: zion(ntypat)
 end subroutine berryphase_new
end interface

interface
 subroutine cgwf(berryopt,cg,cgq,chkexit,cpus,dimffnl,dphase_k,dtefield,&  
  &  ffnl,filnam_ds1,filstat,&  
  &  gsc,gs_hamk,icg,igsc,ikg,ikpt,inonsc,&  
  &  isppol,kg_k,kinpw,lmnmax,matblk,mband,&  
  &  mcg,mcgq,mgfft,mgsc,mkgq,mkmem,mpi_enreg,mpsang,&  
  &  mpssoang,mpw,natom,nband,nbdblock,nkpt,nline,npw,npwarr,nspinor,&  
  &  nsppol,ntypat,nvloc,n4,n5,n6,ortalg,&  
  &  paral_kgb,ph3d,prtvol,pwind,pwind_alloc,pwnsfac,&  
  &  pwnsfacq,quit,resid,subham,subovl,subvnl,tolwfr,&  
  &  use_subovl,vlocal,wfoptalg,wtk,zshift)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: chkexit
  integer,intent(in) :: dimffnl
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: ikg
  integer,intent(in) :: ikpt
  integer,intent(in) :: inonsc
  integer,intent(in) :: isppol
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcgq
  integer,intent(in) :: mgfft
  integer,intent(in) :: mgsc
  integer,intent(in) :: mkgq
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: mpw
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: nbdblock
  integer,intent(in) :: nkpt
  integer,intent(in) :: nline
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: nvloc
  integer,intent(in) :: ortalg
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: pwind_alloc
  integer,intent(inout) :: quit
  integer,intent(in) :: use_subovl
  integer,intent(in) :: wfoptalg
  real(dp),intent(in) :: cpus
  type(efield_type),intent(inout) :: dtefield
  character(len=fnlen),intent(in) :: filnam_ds1
  character(len=fnlen),intent(in) :: filstat
  type(gs_hamiltonian_type),intent(in) :: gs_hamk
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: tolwfr
  real(dp),intent(in) :: wtk
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(in) :: cgq(2,mcgq)
  real(dp),intent(out) :: dphase_k(3)
  real(dp),intent(in) :: ffnl(npw,dimffnl,lmnmax,ntypat)
  real(dp),intent(inout) :: gsc(2,mgsc)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kinpw(npw)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(inout) :: ph3d(2,npw,matblk)
  integer,intent(in) :: pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(in) :: pwnsfacq(2,mkgq)
  real(dp),intent(out) :: resid(nband)
  real(dp),intent(out) :: subham(nband*(nband+1))
  real(dp),intent(out) :: subovl(nband*(nband+1)*use_subovl)
  real(dp),intent(out) :: subvnl(nband*(nband+1)*(1-gs_hamk%usepaw))
  real(dp),intent(inout) :: vlocal(n4,n5,n6,nvloc)
  real(dp),intent(in) :: zshift(nband)
 end subroutine cgwf
end interface

interface
 subroutine inwffil(ask_accurate,cg,dtset,ecut,ecut_eff,eigen,exchn2n3d,&  
  &  formeig,gmet,hdr,ireadwf,istwfk,kg,kptns,localrdwf,mband,&  
  &  mkmem,mpi_enreg,mpw,nband,ngfft,nkpt,npwarr,nspden,nspinor,&  
  &  nsppol,nsym,occ,optorth,psps,prtvol,rprimd,symafm,symrel,tnons,unkg,wff1,&  
  &  wffnow,unwff1,unwfnow,wffnm,wft1nm, wvl)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: ask_accurate
  integer, intent(in) :: exchn2n3d
  integer, intent(in) :: formeig
  integer, intent(in) :: ireadwf
  integer, intent(in) :: localrdwf
  integer, intent(in) :: mband
  integer, intent(in) :: mkmem
  integer, intent(in) :: mpw
  integer, intent(in) :: nkpt
  integer, intent(in) :: nspden
  integer, intent(inout) :: nspinor
  integer, intent(in) :: nsppol
  integer, intent(in) :: nsym
  integer, intent(in) :: optorth
  integer, intent(in) :: prtvol
  integer, intent(in) :: unkg
  integer, intent(in) :: unwff1
  integer, intent(in) :: unwfnow
  type(dataset_type), intent(in) :: dtset
  real(dp), intent(in) :: ecut
  real(dp), intent(in) :: ecut_eff
  type(hdr_type), intent(inout) :: hdr
  type(mpi_type), intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type), intent(out) :: wff1
  character(len=fnlen), intent(in) :: wffnm
  type(wffile_type), intent(out) :: wffnow
  character(len=fnlen), intent(in) :: wft1nm
  type(wvl_data), intent(inout) :: wvl
  integer, intent(in) :: ngfft(18)
  real(dp), intent(out) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp), intent(out) :: eigen((2*mband)**formeig*mband*nkpt*nsppol)
  real(dp), intent(in) :: gmet(3,3)
  integer, intent(in) :: istwfk(nkpt)
  integer, intent(in) :: kg(3,mpw*mkmem)
  real(dp), intent(in) :: kptns(3,nkpt)
  integer, intent(in) :: nband(nkpt*nsppol)
  integer, intent(in) :: npwarr(nkpt)
  real(dp), intent(inout) :: occ(mband*nkpt*nsppol)
  real(dp), intent(in) :: rprimd(3,3)
  integer, intent(in) :: symafm(nsym)
  integer, intent(in) :: symrel(3,3,nsym)
  real(dp), intent(in) :: tnons(3,nsym)
 end subroutine inwffil
end interface

interface
 subroutine iofn1(filnam,filstat,mpi_enreg)
  use defs_basis
  use defs_datatypes
  implicit none
  character(len=fnlen), intent(out) :: filstat
  type(mpi_type), intent(in) :: mpi_enreg
  character(len=fnlen), intent(out) :: filnam(5)
 end subroutine iofn1
end interface

interface
 subroutine iofn2(npsp,pspheads,mpi_enreg)
  use defs_datatypes
  implicit none
  integer,intent(in) :: npsp
  type(mpi_type), intent(inout) :: mpi_enreg
  type(pspheader_type),intent(out) :: pspheads(npsp)
 end subroutine iofn2
end interface

interface
 subroutine lobpcgIIwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&  
  &  kg_k,kinpw,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,mpssoang,natom,&  
  &  nband_k,nbdblock,npw_k,nspinor,ntypat,nvloc,n4,n5,n6,ph3d,prtvol,&  
  &  psps,resid_k,subham,subovl,subvnl,use_subovl,vlocal)
  use defs_basis
  use defs_datatypes
  implicit none
  integer :: dimffnl
  integer :: icg
  integer :: igsc
  integer :: lmnmax
  integer :: matblk
  integer :: mcg
  integer :: mgfft
  integer :: mgsc
  integer :: mpsang
  integer :: mpssoang
  integer :: n4
  integer :: n5
  integer :: n6
  integer :: natom
  integer :: nband_k
  integer :: nbdblock
  integer :: npw_k
  integer :: nspinor
  integer :: ntypat
  integer :: nvloc
  integer :: prtvol
  integer :: use_subovl
  type(datafiles_type) :: dtfil
  type(dataset_type) :: dtset
  type(gs_hamiltonian_type) :: gs_hamk
  type(mpi_type) :: mpi_enreg
  type(pseudopotential_type) :: psps
  real(dp) :: cg(2,mcg)
  real(dp) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  real(dp) :: gsc(2,mgsc)
  integer :: kg_k(3,npw_k)
  real(dp) :: kinpw(npw_k)
  real(dp) :: ph3d(2,npw_k,matblk)
  real(dp) :: resid_k(nband_k)
  real(dp) :: subham(nband_k*(nband_k+1))
  real(dp) :: subovl(nband_k*(nband_k+1)*use_subovl)
  real(dp) :: subvnl(nband_k*(nband_k+1)*(1-gs_hamk%usepaw))
  real(dp) :: vlocal(n4,n5,n6,nvloc)
 end subroutine lobpcgIIwf
end interface

interface
 subroutine lobpcgiiiwf(cg,dimffnl,dtfil,dtset,&  
  &  ffnl,gs_hamk,gsc,icg,igsc,iblock,&  
  &  kg_k,kinpw,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,&  
  &  mpssoang,natom,nbdblock,nband_k,npw_k,nspinor,ntypat,&  
  &  nvloc,n4,n5,n6,&  
  &  ph3d,prtvol,psps,resid_k,vlocal,&  
  &  subham,subvnl,subovl,blocksize,bblocksize,vectsize,pflag,frozen_count,&  
  &  blockvectorx,blockvectorbx,blockvectorax,blockvectory,blockvectorby,lambda,&  
  &  blockvectorp,blockvectorbp,blockvectorap&  
  &  )
  use defs_basis
  use defs_datatypes
  implicit none
  integer :: bblocksize
  integer :: blocksize
  integer :: dimffnl
  integer :: frozen_count
  integer :: iblock
  integer :: icg
  integer :: igsc
  integer :: lmnmax
  integer :: matblk
  integer :: mcg
  integer :: mgfft
  integer :: mgsc
  integer :: mpsang
  integer :: mpssoang
  integer :: n4
  integer :: n5
  integer :: n6
  integer :: natom
  integer :: nband_k
  integer :: nbdblock
  integer :: npw_k
  integer :: nspinor
  integer :: ntypat
  integer :: nvloc
  integer :: prtvol
  integer :: vectsize
  type(datafiles_type) :: dtfil
  type(dataset_type) :: dtset
  type(gs_hamiltonian_type) :: gs_hamk
  type(mpi_type) :: mpi_enreg
  type(pseudopotential_type) :: psps
  real(dp) :: blockvectorap(vectsize,blocksize)
  real(dp) :: blockvectorax(vectsize,blocksize)
  real(dp) :: blockvectorbp(vectsize,blocksize)
  real(dp) :: blockvectorbx(vectsize,blocksize)
  real(dp) :: blockvectorby(vectsize,bblocksize)
  real(dp) :: blockvectorp(vectsize,blocksize)
  real(dp) :: blockvectorx(vectsize,blocksize)
  real(dp) :: blockvectory(vectsize,bblocksize)
  real(dp) :: cg(2,mcg)
  real(dp) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  real(dp) :: gsc(2,mgsc)
  integer :: kg_k(3,npw_k)
  real(dp) :: kinpw(npw_k)
  real(dp) :: lambda(blocksize,blocksize)
  logical :: pflag(blocksize)
  real(dp) :: ph3d(2,npw_k,matblk)
  real(dp) :: resid_k(nband_k)
  real(dp) :: subham(nband_k*(nband_k+1))
  real(dp) :: subovl(nband_k*(nband_k+1))
  real(dp) :: subvnl(nband_k*(nband_k+1))
  real(dp) :: vlocal(n4,n5,n6,nvloc)
 end subroutine lobpcgiiiwf
end interface

interface
 subroutine lobpcgccIIwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&  
  &  kg_k,kinpw,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,mpssoang,natom,&  
  &  nband_k,nbdblock,npw_k,nspinor,ntypat,nvloc,n4,n5,n6,ph3d,prtvol,&  
  &  psps,resid_k,subham,subovl,subvnl,use_subovl,vlocal)
  use defs_basis
  use defs_datatypes
  implicit none
  integer :: dimffnl
  integer :: icg
  integer :: igsc
  integer :: lmnmax
  integer :: matblk
  integer :: mcg
  integer :: mgfft
  integer :: mgsc
  integer :: mpsang
  integer :: mpssoang
  integer :: n4
  integer :: n5
  integer :: n6
  integer :: natom
  integer :: nband_k
  integer :: nbdblock
  integer :: npw_k
  integer :: nspinor
  integer :: ntypat
  integer :: nvloc
  integer :: prtvol
  integer :: use_subovl
  type(datafiles_type) :: dtfil
  type(dataset_type) :: dtset
  type(gs_hamiltonian_type) :: gs_hamk
  type(mpi_type) :: mpi_enreg
  type(pseudopotential_type) :: psps
  real(dp) :: cg(2,mcg)
  real(dp) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  real(dp) :: gsc(2,mgsc)
  integer :: kg_k(3,npw_k)
  real(dp) :: kinpw(npw_k)
  real(dp) :: ph3d(2,npw_k,matblk)
  real(dp) :: resid_k(nband_k)
  real(dp) :: subham(nband_k*(nband_k+1))
  real(dp) :: subovl(nband_k*(nband_k+1)*use_subovl)
  real(dp) :: subvnl(nband_k*(nband_k+1)*(1-gs_hamk%usepaw))
  real(dp) :: vlocal(n4,n5,n6,nvloc)
 end subroutine lobpcgccIIwf
end interface

interface
 subroutine lobpcgcciiiwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&  
  &  kg_k,kinpw,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,&  
  &  mpssoang,natom,nbdblock,nband_k,npw_k,nspinor,ntypat,&  
  &  nvloc,n4,n5,n6,ph3d,prtvol,psps,resid_k,use_subovl,vlocal,&  
  &  subham,subvnl,subovl,blocksize,bblocksize,vectsize,pflag,&  
  &  blockvectorx,blockvectorbx,blockvectorax,blockvectory,blockvectorby,lambda,&  
  &  blockvectorp,blockvectorbp,blockvectorap&  
  &  )
  use defs_basis
  use defs_datatypes
  implicit none
  integer :: bblocksize
  integer :: blocksize
  integer :: dimffnl
  integer :: icg
  integer :: igsc
  integer :: lmnmax
  integer :: matblk
  integer :: mcg
  integer :: mgfft
  integer :: mgsc
  integer :: mpsang
  integer :: mpssoang
  integer :: n4
  integer :: n5
  integer :: n6
  integer :: natom
  integer :: nband_k
  integer :: nbdblock
  integer :: npw_k
  integer :: nspinor
  integer :: ntypat
  integer :: nvloc
  integer :: prtvol
  integer :: use_subovl
  integer :: vectsize
  type(datafiles_type) :: dtfil
  type(dataset_type) :: dtset
  type(gs_hamiltonian_type) :: gs_hamk
  real(dp) :: lambda
  type(mpi_type) :: mpi_enreg
  type(pseudopotential_type) :: psps
  complex(dp) :: blockvectorap(vectsize,blocksize)
  complex(dp) :: blockvectorax(vectsize,blocksize)
  complex(dp) :: blockvectorbp(vectsize,blocksize)
  complex(dp) :: blockvectorbx(vectsize,blocksize)
  complex(dp) :: blockvectorby(vectsize,bblocksize)
  complex(dp) :: blockvectorp(vectsize,blocksize)
  complex(dp) :: blockvectorx(vectsize,blocksize)
  complex(dp) :: blockvectory(vectsize,bblocksize)
  real(dp) :: cg(2,mcg)
  real(dp) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  real(dp) :: gsc(2,mgsc)
  integer :: kg_k(3,npw_k)
  real(dp) :: kinpw(npw_k)
  logical :: pflag(blocksize)
  real(dp) :: ph3d(2,npw_k,matblk)
  real(dp) :: resid_k(nband_k)
  real(dp) :: subham(nband_k*(nband_k+1))
  real(dp) :: subovl(nband_k*(nband_k+1))
  real(dp) :: subvnl(nband_k*(nband_k+1))
  real(dp) :: vlocal(n4,n5,n6,nvloc)
 end subroutine lobpcgcciiiwf
end interface

interface
 subroutine lobpcgccwf(cg,dimffnl,dtfil,dtset,ffnl,ffnl_gather,gs_hamk,gsc,icg,igsc,&  
  &  kg_k,kg_k_gather,kinpw,kinpw_gather,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,mpssoang,natom,&  
  &  nband_k,nbdblock,ndatarecv,npw_k,nspinor,ntypat,nvloc,n4,n5,n6,ph3d,ph3d_gather,prtvol,&  
  &  psps,resid_k,subham,subovl,totvnl,usebandfft,use_subovl,vlocal)
  use defs_basis
  use defs_datatypes
  implicit none
  integer :: dimffnl
  integer :: icg
  integer :: igsc
  integer :: lmnmax
  integer :: matblk
  integer :: mcg
  integer :: mgfft
  integer :: mgsc
  integer :: mpsang
  integer :: mpssoang
  integer :: n4
  integer :: n5
  integer :: n6
  integer :: natom
  integer :: nband_k
  integer :: nbdblock
  integer :: ndatarecv
  integer :: npw_k
  integer :: nspinor
  integer :: ntypat
  integer :: nvloc
  integer :: prtvol
  integer :: use_subovl
  integer :: usebandfft
  type(datafiles_type) :: dtfil
  type(dataset_type) :: dtset
  type(gs_hamiltonian_type) :: gs_hamk
  type(mpi_type) :: mpi_enreg
  type(pseudopotential_type) :: psps
  real(dp) :: cg(2,mcg)
  real(dp) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  real(dp), intent(in) :: ffnl_gather(ndatarecv,dimffnl,lmnmax,ntypat*usebandfft)
  real(dp) :: gsc(2,mgsc)
  integer :: kg_k(3,npw_k)
  integer, intent(in) :: kg_k_gather(3,ndatarecv*usebandfft)
  real(dp) :: kinpw(npw_k)
  real(dp), intent(in) :: kinpw_gather(ndatarecv*usebandfft)
  real(dp) :: ph3d(2,npw_k,matblk)
  real(dp), intent(in) :: ph3d_gather(2,ndatarecv,matblk*usebandfft)
  real(dp) :: resid_k(nband_k)
  real(dp) :: subham(nband_k*(nband_k+1))
  real(dp) :: subovl(nband_k*(nband_k+1)*use_subovl)
  real(dp) :: totvnl(2*nband_k*(1-gs_hamk%usepaw),nband_k*(1-gs_hamk%usepaw))
  real(dp) :: vlocal(n4,n5,n6,nvloc)
 end subroutine lobpcgccwf
end interface

interface
 subroutine lobpcgwf(cg,dimffnl,dtfil,dtset,ffnl,ffnl_gather,gs_hamk,gsc,icg,igsc,&  
  &  kg_k,kg_k_gather,kinpw,kinpw_gather,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,mpssoang,natom,&  
  &  nband_k,nbdblock,ndatarecv,npw_k,nspinor,ntypat,nvloc,n4,n5,n6,ph3d,ph3d_gather,prtvol,&  
  &  psps,resid_k,subham,subovl,totvnl,usebandfft,use_subovl,vlocal)
  use defs_basis
  use defs_datatypes
  implicit none
  integer :: dimffnl
  integer :: icg
  integer :: igsc
  integer :: lmnmax
  integer :: matblk
  integer :: mcg
  integer :: mgfft
  integer :: mgsc
  integer :: mpsang
  integer :: mpssoang
  integer :: n4
  integer :: n5
  integer :: n6
  integer :: natom
  integer :: nband_k
  integer :: nbdblock
  integer :: ndatarecv
  integer :: npw_k
  integer :: nspinor
  integer :: ntypat
  integer :: nvloc
  integer :: prtvol
  integer :: use_subovl
  integer :: usebandfft
  type(datafiles_type) :: dtfil
  type(dataset_type) :: dtset
  type(gs_hamiltonian_type) :: gs_hamk
  type(mpi_type) :: mpi_enreg
  type(pseudopotential_type) :: psps
  real(dp) :: cg(2,mcg)
  real(dp) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  real(dp), intent(in) :: ffnl_gather(ndatarecv,dimffnl,lmnmax,ntypat*usebandfft)
  real(dp) :: gsc(2,mgsc)
  integer :: kg_k(3,npw_k)
  integer, intent(in) :: kg_k_gather(3,ndatarecv*usebandfft)
  real(dp) :: kinpw(npw_k)
  real(dp), intent(in) :: kinpw_gather(ndatarecv*usebandfft)
  real(dp) :: ph3d(2,npw_k,matblk)
  real(dp), intent(in) :: ph3d_gather(2,ndatarecv,matblk*usebandfft)
  real(dp) :: resid_k(nband_k)
  real(dp) :: subham(nband_k*(nband_k+1))
  real(dp) :: subovl(nband_k*(nband_k+1)*use_subovl)
  real(dp) :: totvnl(nband_k*(1-gs_hamk%usepaw),nband_k*(1-gs_hamk%usepaw))
  real(dp) :: vlocal(n4,n5,n6,nvloc)
 end subroutine lobpcgwf
end interface

interface
 subroutine loper3(amass,atindx,atindx1,blkflg,codvsn,cpui,cpus,dimcprj,doccde,&  
  &  ddkfil,dtfil,dtset,dyew,dyfrlo,dyfrnl,dyfrx1,dyfrx2,d2bbb,d2lo,d2nl,&  
  &  eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,&  
  &  etotal,fermie,gsqcut_eff,iexit,indsym,kxc,&  
  &  mkmem,mkqmem,mk1mem,mpert,mpi_enreg,mpsang,nattyp,&  
  &  nfftf,nhat,nkpt,nkxc,nspden,nspinor,nsym,occ,&  
  &  paw_an,paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,&  
  &  pertsy,prtbbb,psps,rfpert,rhog,rhor,symq,symrec,timrev,tmpfil,&  
  &  usecprj,vtrial,vxcavg,walli,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, intent(out) :: iexit
  integer, intent(in) :: mk1mem
  integer, intent(in) :: mkmem
  integer, intent(in) :: mkqmem
  integer, intent(in) :: mpert
  integer, intent(in) :: mpsang
  integer, intent(in) :: nfftf
  integer, intent(in) :: nkpt
  integer, intent(in) :: nkxc
  integer, intent(in) :: nspden
  integer, intent(inout) :: nspinor
  integer, intent(in) :: nsym
  integer, intent(in) :: prtbbb
  integer, intent(in) :: timrev
  integer, intent(in) :: usecprj
  character(len=6), intent(in) :: codvsn
  real(dp), intent(in) :: cpui
  real(dp), intent(in) :: cpus
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(inout) :: dtset
  real(dp), intent(out) :: etotal
  real(dp), intent(inout) :: fermie
  real(dp), intent(in) :: gsqcut_eff
  type(mpi_type), intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type), intent(inout) :: psps
  real(dp), intent(in) :: vxcavg
  real(dp), intent(in) :: walli
  integer, intent(out) :: ddkfil(3)
  real(dp), intent(in) :: amass(dtset%natom)
  integer, intent(in) :: atindx(dtset%natom)
  integer, intent(in) :: atindx1(dtset%natom)
  integer, intent(out) :: blkflg(3,mpert,3,mpert)
  real(dp), intent(out) :: d2bbb(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb)
  real(dp), intent(out) :: d2lo(2,3,mpert,3,mpert)
  real(dp), intent(out) :: d2nl(2,3,mpert,3,mpert)
  integer, intent(in) :: dimcprj(dtset%natom)
  real(dp), intent(in) :: doccde(dtset%mband*nkpt*dtset%nsppol)
  real(dp), intent(in) :: dyew(2,3,dtset%natom,3,dtset%natom)
  real(dp), intent(in) :: dyfrlo(3,3,dtset%natom)
  real(dp), intent(in) :: dyfrnl(3,3,dtset%natom)
  real(dp), intent(in) :: dyfrx1(2,3,dtset%natom,3,dtset%natom)
  real(dp), intent(in) :: dyfrx2(3,3,dtset%natom)
  real(dp), intent(in) :: eltcore(6,6)
  real(dp), intent(in) :: elteew(6+3*dtset%natom,6)
  real(dp), intent(in) :: eltfrhar(6,6)
  real(dp), intent(in) :: eltfrkin(6,6)
  real(dp), intent(in) :: eltfrloc(6+3*dtset%natom,6)
  real(dp), intent(in) :: eltfrnl(6+3*dtset%natom,6)
  real(dp), intent(in) :: eltfrxc(6+3*dtset%natom,6)
  integer, intent(in) :: indsym(4,nsym,dtset%natom)
  real(dp), intent(in) :: kxc(nfftf,nkxc)
  integer, intent(in) :: nattyp(dtset%ntypat)
  real(dp), intent(in) :: nhat(nfftf,nspden*psps%usepaw)
  real(dp), intent(in) :: occ(dtset%mband*nkpt*dtset%nsppol)
  type(paw_an_type),intent(inout) :: paw_an(dtset%natom*psps%usepaw)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  integer, intent(in) :: pertsy(3,mpert)
  integer, intent(in) :: rfpert(mpert)
  real(dp), intent(in) :: rhog(2,nfftf)
  real(dp), intent(in) :: rhor(nfftf,nspden)
  integer, intent(in) :: symq(4,2,nsym)
  integer, intent(in) :: symrec(3,3,nsym)
  character(len=fnlen), intent(in) :: tmpfil(15)
  real(dp), intent(inout) :: vtrial(nfftf,nspden)
  real(dp), intent(inout) :: xred(3,dtset%natom)
 end subroutine loper3
end interface

interface
 subroutine mv_3dte(cg,cgindex,cg1,cg3,dtset,dtfil,d3_berry,gmet,gprimd,&  
  &  hdr,i1pert,i2pert,i3pert,i1dir,i3dir,kneigh,kptindex,&  
  &  kpt3,mband,mkmem,mkmem_max,mk1mem,mpert,mpi_enreg,&  
  &  mpw,mvwtk,natom,nkpt2,nkpt3,nneigh,npwarr,nspinor,&  
  &  nsppol,occ,pwind,psps,rmet,tmpfil)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, intent(in) :: i1dir
  integer, intent(in) :: i1pert
  integer, intent(in) :: i2pert
  integer, intent(in) :: i3dir
  integer, intent(in) :: i3pert
  integer, intent(in) :: mband
  integer, intent(in) :: mk1mem
  integer, intent(in) :: mkmem
  integer, intent(in) :: mkmem_max
  integer, intent(in) :: mpert
  integer, intent(in) :: mpw
  integer, intent(in) :: natom
  integer, intent(in) :: nkpt2
  integer, intent(in) :: nkpt3
  integer, intent(in) :: nneigh
  integer, intent(in) :: nspinor
  integer, intent(in) :: nsppol
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(in) :: dtset
  type(hdr_type), intent(in) :: hdr
  type(mpi_type), intent(inout) :: mpi_enreg
  type(pseudopotential_type), intent(in) :: psps
  real(dp), intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp), intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
  real(dp), intent(in) :: cg3(2,mpw*nspinor*mband*mk1mem*nsppol)
  integer, intent(in) :: cgindex(nkpt2,nsppol)
  real(dp), intent(out) :: d3_berry(2,3)
  real(dp), intent(in) :: gmet(3,3)
  real(dp), intent(in) :: gprimd(3,3)
  integer, intent(in) :: kneigh(30,nkpt2)
  real(dp), intent(in) :: kpt3(3,nkpt3)
  integer, intent(in) :: kptindex(2,nkpt3)
  real(dp), intent(in) :: mvwtk(30,nkpt2)
  integer, intent(in) :: npwarr(nkpt2)
  real(dp), intent(in) :: occ(mband*nkpt2*nsppol)
  integer, intent(in) :: pwind(mpw,nneigh,mkmem)
  real(dp), intent(in) :: rmet(3,3)
  character(len=fnlen), intent(in) :: tmpfil(15)
 end subroutine mv_3dte
end interface

interface
 subroutine outkss(atindx,atindx1,Dtfil,Dtset,ecut,gmet,gprimd,Hdr,kg,&  
  &  kssform,mband,mgfft,mkmem,MPI_enreg,mpsang,mpw,natom,&  
  &  nattyp,nfft,nkpt,npwarr,nspinor,nspden,nsppol,nsym,ntypat,occ,Pawrad,pawtab,&  
  &  ph1d,prtvol,Psps,rmet,rprimd,ucvol,Wffnow,vtrial,xred,ylm,cg,usecprj,Cprj,eigen)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: kssform
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(inout) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  integer,intent(in) :: usecprj
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(in) :: Dtset
  type(hdr_type),intent(inout) :: Hdr
  type(mpi_type),intent(inout) :: MPI_enreg
  type(pseudopotential_type),intent(in) :: Psps
  type(wffile_type),intent(inout) :: Wffnow
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ucvol
  type(cprj_type),intent(in) :: Cprj(natom,nspinor*mband*mkmem*nsppol*Psps%usepaw*usecprj)
  type(pawrad_type),intent(in) :: Pawrad(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in),target :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: vtrial(nfft,nspden)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*Psps%useylm)
 end subroutine outkss
end interface

interface
 subroutine outwf(cg,dtfil,dtset,eigen,filnam,hdr,kg,kptns,mband,mkmem,&  
  &  mpi_enreg,mpw,mxfh,natom,nband,nfft,ngfft,nkpt,npwarr,&  
  &  nqpt,nspinor,nsppol,nstep,nxfh,occ,resid,response,&  
  &  wffnow,wfs,xfhist)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: mband
  integer, intent(in) :: mkmem
  integer, intent(in) :: mpw
  integer, intent(in) :: mxfh
  integer, intent(in) :: natom
  integer, intent(in) :: nfft
  integer, intent(in) :: nkpt
  integer, intent(in) :: nqpt
  integer, intent(in) :: nspinor
  integer, intent(in) :: nsppol
  integer, intent(in) :: nstep
  integer, intent(in) :: nxfh
  integer, intent(in) :: response
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(in) :: dtset
  character(len=fnlen), intent(in) :: filnam
  type(hdr_type), intent(inout) :: hdr
  type(mpi_type), intent(inout) :: mpi_enreg
  type(wffile_type), intent(inout) :: wffnow
  type(wvl_wf_type), intent(in) :: wfs
  integer, intent(in) :: ngfft(18)
  real(dp), intent(inout) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp), intent(in) :: eigen((2*mband)**response*mband*nkpt*nsppol)
  integer, intent(in) :: kg(3,mpw*mkmem)
  real(dp), intent(in) :: kptns(3,nkpt)
  integer, intent(in) :: nband(nkpt*nsppol)
  integer, intent(in) :: npwarr(nkpt)
  real(dp), intent(in) :: occ(mband*nkpt*nsppol)
  real(dp), intent(in) :: resid(mband*nkpt*nsppol)
  real(dp), intent(in) :: xfhist(3,natom+4,2,mxfh)
 end subroutine outwf
end interface

interface
 subroutine prep_nonlop(atindx1,choice,cpopt,cprj_block,dimenl1,dimenl2,dimffnl,enl,enlout_block,&  
  &  ffnl_gather,gmet,gprimd,iblock,icall,idir,indlmn,istwf_k,kg_k_gather,&  
  &  kpg,kpt,lambdablock,lmnmax,matblk,&  
  &  blocksize,mgfft,mpi_enreg,mpsang,mpssoang,&  
  &  natom,nattyp,nbdblock,nband_k,dimtabs,ngfft,nkpg,nloalg,nnlout,npw_k,&  
  &  nspinor,ntypat,paral_kgb,paw_opt,phkxred,ph1d,ph3d_gather,prtvol,pspso,signs,sij,gsc,&  
  &  tim_nonlop,ucvol,useylm,cwavef,gvnlc)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: blocksize
  integer,intent(in) :: choice
  integer,intent(in) :: cpopt
  integer,intent(in) :: dimenl1
  integer,intent(in) :: dimenl2
  integer,intent(in) :: dimffnl
  integer,intent(in) :: dimtabs
  integer,intent(in) :: iblock
  integer,intent(in) :: icall
  integer,intent(in) :: idir
  integer,intent(in) :: istwf_k
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: natom
  integer,intent(in) :: nband_k
  integer,intent(in) :: nbdblock
  integer,intent(in) :: nkpg
  integer,intent(in) :: nnlout
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: paw_opt
  integer,intent(in) :: prtvol
  integer,intent(in) :: signs
  integer :: tim_nonlop
  integer,intent(in) :: useylm
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(5)
  integer,intent(in) :: atindx1(natom)
  type(cprj_type) :: cprj_block(natom,nspinor*blocksize*((cpopt+3)/3))
  real(dp),intent(inout) :: cwavef(2,npw_k*nspinor*blocksize)
  real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinor**2)
  real(dp),intent(out) :: enlout_block(nnlout*blocksize)
  real(dp),intent(in) :: ffnl_gather(dimtabs,dimffnl,lmnmax,ntypat)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: gsc(2,npw_k*nspinor*blocksize*(paw_opt/3))
  real(dp),intent(out) :: gvnlc(2,npw_k*nspinor*blocksize)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: kg_k_gather(3,dimtabs)
  real(dp),intent(in) :: kpg(npw_k,nkpg)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(in) :: lambdablock(blocksize)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(inout) :: ph3d_gather(2,dimtabs,matblk)
  real(dp),intent(in) :: phkxred(2,natom)
  integer,intent(in) :: pspso(ntypat)
  real(dp),intent(in) :: sij(dimenl1,ntypat*((paw_opt+1)/3))
 end subroutine prep_nonlop
end interface

interface
 subroutine respfn(codvsn,cpui,dtfil,dtset,etotal,iexit,&  
  &  mkmems,mpi_enreg,npwtot,&  
  &  nspinor,occ,pawang,pawrad,pawtab,psps,walli,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(inout) :: iexit
  integer,intent(inout) :: nspinor
  character(len=6),intent(in) :: codvsn
  real(dp),intent(in) :: cpui
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(out) :: etotal
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(inout) :: pawang
  type(pseudopotential_type),intent(inout) :: psps
  real(dp),intent(in) :: walli
  integer,intent(in) :: mkmems(3)
  integer,intent(inout) :: npwtot(dtset%nkpt)
  real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine respfn
end interface

interface
 subroutine subdiago(cg,filstat,eig_k,evec,gsc,icg,igsc,ikpt,inonsc,istwf_k,&  
  &  mcg,mgsc,mpi_enreg,nband_k,npw_k,nspinor,paral_kgb,prtvol,&  
  &  subham,subovl,use_subovl,usepaw)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: ikpt
  integer,intent(in) :: inonsc
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mcg
  integer,intent(in) :: mgsc
  integer,intent(in) :: nband_k
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: use_subovl
  integer,intent(in) :: usepaw
  character(len=fnlen),intent(in) :: filstat
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(out) :: eig_k(nband_k)
  real(dp),intent(out) :: evec(2*nband_k,nband_k)
  real(dp),intent(inout) :: gsc(2,mgsc)
  real(dp),intent(inout) :: subham(nband_k*(nband_k+1))
  real(dp),intent(inout) :: subovl(nband_k*(nband_k+1)*use_subovl)
 end subroutine subdiago
end interface

interface
 subroutine tddft(cg,dtfil,dtset,eigen,etotal,gmet,gprimd,gsqcut,&  
  &  kg,kxc,mband,mgfftdiel,mkmem,mpi_enreg,mpw,nfft,ngfftdiel,nkpt,nkxc,&  
  &  npwarr,nspinor,nsppol,occ,ucvol,wffnew)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, intent(in) :: mband
  integer, intent(in) :: mgfftdiel
  integer, intent(in) :: mkmem
  integer, intent(in) :: mpw
  integer, intent(in) :: nfft
  integer, intent(in) :: nkpt
  integer, intent(in) :: nkxc
  integer, intent(inout) :: nspinor
  integer, intent(in) :: nsppol
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(in) :: dtset
  real(dp), intent(in) :: etotal
  real(dp), intent(in) :: gsqcut
  type(mpi_type), intent(inout) :: mpi_enreg
  real(dp), intent(in) :: ucvol
  type(wffile_type), intent(inout) :: wffnew
  integer, intent(in) :: ngfftdiel(18)
  real(dp), intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp), intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp), intent(in) :: gmet(3,3)
  real(dp), intent(in) :: gprimd(3,3)
  integer, intent(in) :: kg(3,mpw*mkmem)
  real(dp), intent(in) :: kxc(nfft,nkxc)
  integer, intent(in) :: npwarr(nkpt)
  real(dp), intent(in) :: occ(mband*nkpt*nsppol)
 end subroutine tddft
end interface

interface
 subroutine vtorho(afford,atindx,atindx1,cg,compch_fft,cpus,dbl_nnsclo,&  
  &  densymop_diel,densymop_gs,dielop,dielstrt,dphase,dtefield,dtfil,dtset,&  
  &  eigen,energies,etotal,filapp,gbound_diel,&  
  &  gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&  
  &  istep,kg,kg_diel,kxc,lmax_diel,mgfftdiel,mpi_enreg,&  
  &  mpsang,natom,nattyp,nfftf,nfftdiel,ngfftdiel,nhat,nkxc,&  
  &  npwarr,npwdiel,nres2,nspinor,ntypat,nvresid,occ,optforces,&  
  &  optres,paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&  
  &  pwind,pwind_alloc,pwnsfac,resid,residm,rhog,rhor,&  
  &  rmet,rprimd,shiftvector,susmat,symrec,&  
  &  ucvol,wffnew,wffnow,val_min,val_max,vtrial,wvl,xred,ylm,ylmdiel)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: afford
  integer, intent(in) :: dbl_nnsclo
  integer, intent(in) :: dielop
  integer, intent(in) :: dielstrt
  integer, intent(in) :: istep
  integer, intent(in) :: lmax_diel
  integer, intent(in) :: mgfftdiel
  integer, intent(in) :: mpsang
  integer, intent(in) :: natom
  integer, intent(in) :: nfftdiel
  integer, intent(in) :: nfftf
  integer, intent(in) :: nkxc
  integer, intent(in) :: npwdiel
  integer, intent(inout) :: nspinor
  integer, intent(in) :: ntypat
  integer, intent(in) :: optforces
  integer, intent(in) :: optres
  integer, intent(in) :: pwind_alloc
  real(dp), intent(out) :: compch_fft
  real(dp), intent(in) :: cpus
  type(dens_sym_operator_type), intent(in) :: densymop_diel
  type(dens_sym_operator_type), intent(in) :: densymop_gs
  type(efield_type), intent(inout) :: dtefield
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(inout) :: dtset
  type(energies_type), intent(inout) :: energies
  real(dp), intent(in) :: etotal
  character(len=fnlen), intent(in) :: filapp
  real(dp), intent(in) :: gsqcut
  type(hdr_type), intent(inout) :: hdr
  type(mpi_type), intent(inout) :: mpi_enreg
  real(dp), intent(out) :: nres2
  type(pawang_type), intent(in) :: pawang
  type(pawfgr_type), intent(in) :: pawfgr
  type(pseudopotential_type), intent(in) :: psps
  real(dp), intent(out) :: residm
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: val_max
  real(dp), intent(in) :: val_min
  type(wffile_type), intent(inout) :: wffnew
  type(wffile_type), intent(inout) :: wffnow
  type(wvl_data), intent(inout) :: wvl
  integer, intent(in) :: ngfftdiel(18)
  integer, intent(in) :: atindx(natom)
  integer, intent(in) :: atindx1(natom)
  real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp), intent(out) :: dphase(3)
  real(dp), intent(out) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer, intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  real(dp), intent(in) :: gmet(3,3)
  real(dp), intent(in) :: gprimd(3,3)
  real(dp), intent(out) :: grnl(3*natom)
  integer, intent(in) :: indsym(4,dtset%nsym,natom)
  integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
  integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: kg_diel(3,npwdiel)
  real(dp), intent(inout) :: kxc(nfftf,nkxc)
  integer, intent(in) :: nattyp(ntypat)
  real(dp), intent(out) :: nhat(nfftf,dtset%nspden*psps%usepaw)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(out) :: nvresid(nfftf,dtset%nspden)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp), intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*natom)
  real(dp), intent(in) :: ph1ddiel(2,(3*(2*mgfftdiel+1)*natom)*psps%usepaw)
  real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
  real(dp), intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(inout) :: rhog(2,nfftf)
  real(dp), intent(inout) :: rhor(nfftf,dtset%nspden)
  real(dp), intent(in) :: rmet(3,3)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(in) :: shiftvector((dtset%mband+2)*dtset%nkpt)
  real(dp), intent(out) :: susmat(2,npwdiel*afford,dtset%nspden,npwdiel,dtset%nspden)
  integer, intent(in) :: symrec(3,3,dtset%nsym)
  real(dp), intent(inout) :: vtrial(nfftf,dtset%nspden)
  real(dp), intent(inout) :: xred(3,natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,mpsang*mpsang*psps%useylm)
  real(dp), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
 end subroutine vtorho
end interface

interface
 subroutine vtowfk(cg,cgq,cprj,cpus,dimcprj,dimffnl,dphase_k,dtefield,dtfil,dtset,&  
  &  eig_k,ek_k,enl_k,ffnl,ffnl_gather,grnl_k,gs_hamk,&  
  &  ibg,icg,ikg,ikpt,iscf,isppol,kg_k,kg_k_gather,kinpw,kinpw_gather,kpg_k,&  
  &  lmnmax,matblk,mcg,mcgq,mgfft,mkgq,mkmem,mpi_enreg,mpsang,&  
  &  mpssoang,mpw,natom,nband_k,ndatarecv,nkpg,nkpt,nnsclo_now,npw_k,npwarr,&  
  &  nspinor,ntypat,nvloc,n4,n5,n6,occ_k,optforces,ph3d,ph3d_gather,prtvol,psps,&  
  &  pwind,pwind_alloc,pwnsfac,pwnsfacq,resid_k,rhoaug,usebandfft,usecprj,vlocal,wtk,zshift)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, intent(in) :: dimffnl
  integer, intent(in) :: ibg
  integer, intent(in) :: icg
  integer, intent(in) :: ikg
  integer, intent(in) :: ikpt
  integer, intent(in) :: iscf
  integer, intent(in) :: isppol
  integer, intent(in) :: lmnmax
  integer, intent(in) :: matblk
  integer, intent(in) :: mcg
  integer, intent(in) :: mcgq
  integer, intent(in) :: mgfft
  integer, intent(in) :: mkgq
  integer, intent(in) :: mkmem
  integer, intent(in) :: mpsang
  integer, intent(in) :: mpssoang
  integer, intent(in) :: mpw
  integer, intent(in) :: n4
  integer, intent(in) :: n5
  integer, intent(in) :: n6
  integer, intent(in) :: natom
  integer, intent(in) :: nband_k
  integer, intent(in) :: ndatarecv
  integer, intent(in) :: nkpg
  integer, intent(in) :: nkpt
  integer, intent(in) :: nnsclo_now
  integer, intent(in) :: npw_k
  integer, intent(in) :: nspinor
  integer, intent(in) :: ntypat
  integer, intent(in) :: nvloc
  integer, intent(in) :: optforces
  integer, intent(in) :: prtvol
  integer, intent(in) :: pwind_alloc
  integer, intent(in) :: usebandfft
  integer, intent(in) :: usecprj
  real(dp), intent(in) :: cpus
  type(efield_type), intent(inout) :: dtefield
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(in) :: dtset
  type(gs_hamiltonian_type), intent(in) :: gs_hamk
  type(mpi_type), intent(inout) :: mpi_enreg
  type(pseudopotential_type), intent(in) :: psps
  real(dp), intent(in) :: wtk
  real(dp), intent(inout) :: cg(2,mcg)
  real(dp), intent(in) :: cgq(2,mcgq)
  type(cprj_type) :: cprj(natom,nspinor*dtset%mband*dtset%mkmem*dtset%nsppol*gs_hamk%usepaw)
  integer, intent(in) :: dimcprj(natom*gs_hamk%usepaw)
  real(dp), intent(out) :: dphase_k(3)
  real(dp), intent(out) :: eig_k(nband_k)
  real(dp), intent(out) :: ek_k(nband_k)
  real(dp), intent(out) :: enl_k(nband_k*(1-gs_hamk%usepaw))
  real(dp), intent(in) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  real(dp), intent(in) :: ffnl_gather(ndatarecv,dimffnl,lmnmax,ntypat*usebandfft)
  real(dp), intent(out) :: grnl_k(3*natom,nband_k*optforces)
  integer, intent(in) :: kg_k(3,npw_k)
  integer, intent(in) :: kg_k_gather(3,ndatarecv*usebandfft)
  real(dp), intent(in) :: kinpw(npw_k)
  real(dp), intent(in) :: kinpw_gather(ndatarecv*usebandfft)
  real(dp), intent(in) :: kpg_k(npw_k,nkpg)
  integer, intent(in) :: npwarr(nkpt)
  real(dp), intent(in) :: occ_k(nband_k)
  real(dp), intent(inout) :: ph3d(2,npw_k,matblk)
  real(dp), intent(inout) :: ph3d_gather(2,ndatarecv,matblk*usebandfft)
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(in) :: pwnsfacq(2,mkgq)
  real(dp), intent(out) :: resid_k(nband_k)
  real(dp), intent(inout) :: rhoaug(n4,n5,n6,nvloc)
  real(dp), intent(inout) :: vlocal(n4,n5,n6,nvloc)
  real(dp), intent(in) :: zshift(nband_k)
 end subroutine vtowfk
end interface

interface
 subroutine wfsinp(cg,cg_disk,ecut,ecut0,ecut_eff,eigen,exchn2n3d,&  
  &  formeig,gmet,gmet0,headform0,indkk,indkk0,istwfk,&  
  &  istwfk0,kptns,kptns0,localrdwf,mband,mban_dp_rd,&  
  &  mcg_disk,mkmem,mpi_enreg,mpw,mpw0,nband,nban_dp_rd,&  
  &  ngfft,nkassoc,nkpt,nkpt0,npwarr,npwarr0,nspinor,&  
  &  nspinor0,nsppol,nsppol0,nsym,occ,optorth,prtvol,restart,rprimd,&  
  &  sppoldbl,squeeze,symafm,symrel,tnons,wff1,wffnow)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, intent(in) :: exchn2n3d
  integer, intent(in) :: formeig
  integer, intent(in) :: headform0
  integer, intent(in) :: localrdwf
  integer, intent(in) :: mban_dp_rd
  integer, intent(in) :: mband
  integer, intent(in) :: mcg_disk
  integer, intent(in) :: mkmem
  integer, intent(in) :: mpw
  integer, intent(in) :: mpw0
  integer, intent(in) :: nkassoc
  integer, intent(in) :: nkpt
  integer, intent(in) :: nkpt0
  integer, intent(in) :: nspinor
  integer, intent(in) :: nspinor0
  integer, intent(in) :: nsppol
  integer, intent(in) :: nsppol0
  integer, intent(in) :: nsym
  integer, intent(in) :: optorth
  integer, intent(in) :: prtvol
  integer, intent(in) :: restart
  integer, intent(in) :: sppoldbl
  integer, intent(in) :: squeeze
  real(dp), intent(in) :: ecut
  real(dp), intent(in) :: ecut0
  real(dp), intent(in) :: ecut_eff
  type(mpi_type), intent(inout) :: mpi_enreg
  type(wffile_type), intent(inout) :: wff1
  type(wffile_type), intent(inout) :: wffnow
  integer, intent(in) :: ngfft(18)
  real(dp), intent(out) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp), intent(out) :: cg_disk(2,mcg_disk)
  real(dp), intent(out) :: eigen((2*mband)**formeig*mband*nkpt*nsppol)
  real(dp), intent(in) :: gmet(3,3)
  real(dp), intent(in) :: gmet0(3,3)
  integer, intent(in) :: indkk(nkpt*sppoldbl,6)
  integer, intent(in) :: indkk0(nkpt0,nkassoc)
  integer, intent(in) :: istwfk(nkpt)
  integer, intent(in) :: istwfk0(nkpt0)
  real(dp), intent(in) :: kptns(3,nkpt)
  real(dp), intent(in) :: kptns0(3,nkpt0)
  integer, intent(in) :: nban_dp_rd(nkpt0*nsppol0)
  integer, intent(in) :: nband(nkpt*nsppol)
  integer, intent(in) :: npwarr(nkpt)
  integer, intent(in) :: npwarr0(nkpt0)
  real(dp), intent(inout) :: occ(mband*nkpt*nsppol)
  real(dp), intent(in) :: rprimd(3,3)
  integer, intent(in) :: symafm(nsym)
  integer, intent(in) :: symrel(3,3,nsym)
  real(dp), intent(in) :: tnons(3,nsym)
 end subroutine wfsinp
end interface

interface
 subroutine wvl_vtorho(dtset, energies, istep, mpi_enreg,&  
  &  occ, proj, psps, residm, rhor, vtrial, wfs)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: istep
  type(dataset_type), intent(inout) :: dtset
  type(energies_type), intent(inout) :: energies
  type(mpi_type), intent(in) :: mpi_enreg
  type(wvl_projectors_type), intent(in) :: proj
  type(pseudopotential_type), intent(in) :: psps
  real(dp), intent(inout) :: residm
  type(wvl_wf_type), intent(inout) :: wfs
  real(dp), intent(in) :: occ(dtset%mband * dtset%nsppol)
  real(dp), intent(inout) :: rhor(dtset%nfft)
  real(dp), intent(in) :: vtrial(dtset%nfft * dtset%nspden)
 end subroutine wvl_vtorho
end interface

interface
 subroutine wvl_wfsinp_disk(dtset, hdr0, hdr, mpi_enreg, option,&  
  &  rprimd, wff, wfs, xred)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: option
  type(dataset_type), intent(in) :: dtset
  type(hdr_type), intent(in) :: hdr
  type(hdr_type), intent(in) :: hdr0
  type(mpi_type), intent(in) :: mpi_enreg
  type(wffile_type), intent(in) :: wff
  type(wvl_wf_type), intent(inout) :: wfs
  real(dp), intent(in) :: rprimd(3, 3)
  real(dp), intent(inout) :: xred(3, dtset%natom)
 end subroutine wvl_wfsinp_disk
end interface

interface
 subroutine wvl_wfsinp_reformat(acell, dtset, mpi_enreg, occ, psps,&  
  &  rprimd, wvl, xred, xred_old)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  type(dataset_type), intent(inout) :: dtset
  type(mpi_type), intent(inout) :: mpi_enreg
  type(pseudopotential_type), intent(in) :: psps
  type(wvl_data), intent(inout) :: wvl
  real(dp), intent(inout) :: acell(3)
  real(dp), intent(in) :: occ(dtset%mband * dtset%nkpt *dtset%nsppol)
  real(dp), intent(inout) :: rprimd(3,3)
  real(dp), intent(inout) :: xred(3, dtset%natom)
  real(dp), intent(inout) :: xred_old(3, dtset%natom)
 end subroutine wvl_wfsinp_reformat
end interface

interface
 subroutine wvl_wfsinp_scratch(dtset, mpi_enreg, psps, rprimd, wvl, xred)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  type(dataset_type), intent(in) :: dtset
  type(mpi_type), intent(in) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_data), intent(inout) :: wvl
  real(dp), intent(in) :: rprimd(3, 3)
  real(dp), intent(inout) :: xred(3, dtset%natom)
 end subroutine wvl_wfsinp_scratch
end interface

end module interfaces_18seqpar
!!***
