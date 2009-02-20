!!****m* ABINIT/interfaces_14iowfdenpot
!! NAME
!! interfaces_14iowfdenpot
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/14iowfdenpot
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

module interfaces_14iowfdenpot

 implicit none

interface
 subroutine WffReadEigK(eigen,formeig,headform,ikpt,isppol,mband,mpi_enreg,nband,tim_rwwf,wff)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: headform
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: nband
  integer,intent(in) :: tim_rwwf
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(out) :: eigen((2*mband)**formeig*mband)
 end subroutine WffReadEigK
end interface

interface
 subroutine WffReadSkipK(formeig,headform,ikpt,isppol,mpi_enreg,wff)
  use defs_datatypes
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: headform
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
 end subroutine WffReadSkipK
end interface

interface
 subroutine calcdensph(gmet,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,ratsph,rhor,rprimd,typat,ucvol,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: ratsph(ntypat)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine calcdensph
end interface

interface
 subroutine fappnd(filapp,filnam,iapp)
  use defs_basis
  implicit none
  integer,intent(in) :: iapp
  character(len=fnlen),intent(out) :: filapp
  character(len=fnlen),intent(in) :: filnam
 end subroutine fappnd
end interface

interface
 subroutine hdr_check(fform,fform0,hdr,hdr0,mode_paral,restart,restartpaw)
  use defs_datatypes
  implicit none
  integer,intent(in) :: fform
  integer,intent(in) :: fform0
  integer,intent(out) :: restart
  integer,intent(out) :: restartpaw
  type(hdr_type),intent(in) :: hdr
  type(hdr_type),intent(in) :: hdr0
  character(len=4),intent(in) :: mode_paral
 end subroutine hdr_check
end interface

interface
 subroutine hdr_clean(hdr)
  use defs_datatypes
  implicit none
  type(hdr_type),intent(inout) :: hdr
 end subroutine hdr_clean
end interface

interface
 subroutine hdr_init(bstruct,codvsn,dtset,hdr,pawtab,pertcase,psps)
  use defs_datatypes
  implicit none
  integer,intent(in) :: pertcase
  type(bandstructure_type),intent(in) :: bstruct
  character(len=6),intent(in) :: codvsn
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(out) :: hdr
  type(pseudopotential_type),intent(in) :: psps
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)
 end subroutine hdr_init
end interface

interface
 subroutine hdr_update(bantot,etot,fermie,hdr,natom,&  
  &  residm,rprimd,occ,pawrhoij,usepaw,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: bantot
  integer,intent(in) :: natom
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: etot
  real(dp),intent(in) :: fermie
  type(hdr_type),intent(out) :: hdr
  real(dp),intent(in) :: residm
  real(dp),intent(in) :: occ(bantot)
  type(pawrhoij_type),intent(in) :: pawrhoij(:)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine hdr_update
end interface

interface
 subroutine initwf(cg,eig_k,formeig,headform,icg,ikpt,ikptsp_old,&  
  &  isppol,mcg,mpi_enreg,&  
  &  nband_k,nkpt,npw,nspinor,occ_k,wff1)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: headform
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(inout) :: ikptsp_old
  integer,intent(in) :: isppol
  integer,intent(in) :: mcg
  integer,intent(in) :: nband_k
  integer,intent(in) :: nkpt
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff1
  real(dp),intent(out) :: cg(2,mcg)
  real(dp),intent(out) :: eig_k((2*nband_k)**formeig*nband_k)
  real(dp),intent(inout) :: occ_k(nband_k)
 end subroutine initwf
end interface

interface
 subroutine ioarr(accessfil,arr,dtset,etotal,fform,fildata,hdr,mpi_enreg,&  
  &  ncplxfft,pawrhoij,rdwr,rdwrpaw,ngfft)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: accessfil
  integer,intent(inout) :: fform
  integer,intent(in) :: ncplxfft
  integer,intent(in) :: rdwr
  integer,intent(in) :: rdwrpaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(inout) :: etotal
  character(len=fnlen),intent(in) :: fildata
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout),target :: arr(ncplxfft,dtset%nspden)
  type(pawrhoij_type),intent(inout) :: pawrhoij(hdr%natom*hdr%usepaw*rdwrpaw)
 end subroutine ioarr
end interface

interface
 subroutine mk_hdr_check_fmt(nelm,typfmt)
  implicit none
  integer,intent(in) :: nelm
  character(len=26),intent(out) :: typfmt
 end subroutine mk_hdr_check_fmt
end interface

interface
 subroutine out1dm(filapp,natom,nfft,ngfft,nspden,ntypat,&  
  &  rhor,rprimd,typat,ucvol,vtrial,xred,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  character(len=fnlen),intent(inout) :: filapp
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vtrial(nfft,nspden)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine out1dm
end interface

interface
 subroutine outwant(dtfil,dtset,eig,cg,kg,npwarr,mband,mpi_enreg,nkpt,nsppol,&  
  nspinor,mkmem,mpw,wff,prtwant)
  use defs_basis
  use defs_datatypes
  implicit none
  integer :: mband
  integer :: mkmem
  integer :: mpw
  integer :: nkpt
  integer :: nspinor
  integer :: nsppol
  integer :: prtwant
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
  real(dp) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp) :: eig(mband*nkpt*nsppol)
  integer :: kg(3,mpw*mkmem)
  integer :: npwarr(nkpt)
 end subroutine outwant
end interface

interface
 subroutine randac(debug,headform1,ikptsp_prev,ikpt,isppol,&  
  &  mband,nband,nkpt,nsppol,wffinp)
  use defs_datatypes
  implicit none
  integer,intent(in) :: debug
  integer,intent(in) :: headform1
  integer,intent(in) :: ikpt
  integer,intent(inout) :: ikptsp_prev
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(wffile_type),intent(inout) :: wffinp
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine randac
end interface

interface
 subroutine rdkss(Dtfil,Dtset,Pawtab,nsym_gw,nbndsA,nbvw,nkibz,npwvec,nspinor,nsppol,npwwfn,tit,symrec,gvec,&  
  &  en,occ,wf,Cprj_ibz,ntypat,natom,mpsang,tnons,vkbsign,vkb,vkbd,nel,MPI_enreg,my_minb,my_maxb,&  
  &  wf_val) ! Optional arguments 
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mpsang
  integer,intent(in) :: my_maxb
  integer,intent(in) :: my_minb
  integer,intent(in) :: natom
  integer,intent(in) :: nbndsA
  integer,intent(in) :: nbvw
  integer,intent(out) :: nel
  integer,intent(in) :: nkibz
  integer,intent(in) :: npwvec
  integer,intent(in) :: npwwfn
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym_gw
  integer,intent(in) :: ntypat
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(in) :: Dtset
  type(mpi_type),intent(in) :: MPI_enreg
  character(len=80),intent(out) :: tit(2)
  type(cprj_type),intent(out) :: Cprj_ibz(natom,nspinor*nbndsA*nkibz*nsppol*Dtset%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Dtset%ntypat*Dtset%usepaw)
  real(dp),intent(out) :: en(nkibz,nbndsA,nsppol)
  integer,intent(out) :: gvec(3,npwvec)
  real(dp),intent(out) :: occ(nkibz,nbndsA,nsppol)
  integer,intent(out) :: symrec(3,3,nsym_gw)
  real(dp),intent(out) :: tnons(3,nsym_gw)
  real(dp),intent(out) :: vkb(npwwfn,ntypat,mpsang,nkibz)
  real(dp),intent(out) :: vkbd(npwwfn,ntypat,mpsang,nkibz)
  real(dp),intent(out) :: vkbsign(mpsang,ntypat)
  complex(gwpc),intent(inout),target :: wf(npwwfn*nspinor,my_minb:my_maxb,nkibz,nsppol)
  complex(gwpc),intent(inout),optional,target :: wf_val(npwwfn*nspinor,nbvw,nkibz,nsppol)
 end subroutine rdkss
end interface

interface
 subroutine rdnpw(ikpt,isppol,nband_k,npw_k,nspinor,option,unitfile)
  implicit none
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(inout) :: nband_k
  integer,intent(inout) :: npw_k
  integer,intent(inout) :: nspinor
  integer,intent(in) :: option
  integer,intent(in) :: unitfile
 end subroutine rdnpw
end interface

interface
 subroutine read_wfrspa(state,dtfil,eigbnd,iband,isppol,imkmem,ndiel4,ndiel5,ndiel6,wfrspa_extract)
  use defs_basis
  use defs_datatypes
  implicit none
  integer :: iband
  integer,intent(in) :: imkmem
  integer :: isppol
  integer,intent(in) :: ndiel4
  integer,intent(in) :: ndiel5
  integer,intent(in) :: ndiel6
  integer,intent(in) :: state
  type(datafiles_type),intent(in) :: dtfil
  real(dp),intent(inout) :: eigbnd
  real(dp),intent(inout) :: wfrspa_extract(ndiel4,ndiel5,ndiel6)
 end subroutine read_wfrspa
end interface

interface
 subroutine testlda(Dtset,Dtfil,nsym,nbnds_kss,ng_kss,mpsang,gvec_p,Hdr,MPI_enreg,&  
  ibocc) ! Optional
  use defs_datatypes
  implicit none
  integer,intent(out) :: mpsang
  integer,intent(out) :: nbnds_kss
  integer,intent(out) :: ng_kss
  integer,intent(out) :: nsym
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(in) :: Dtset
  type(hdr_type),intent(out) :: Hdr
  type(mpi_type),intent(in) :: MPI_enreg
  integer,pointer :: gvec_p(:,:)
  integer,intent(inout),optional :: ibocc(Dtset%nsppol)
 end subroutine testlda
end interface

end module interfaces_14iowfdenpot
!!***
