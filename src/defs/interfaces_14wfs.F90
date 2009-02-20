!!****m* ABINIT/interfaces_14wfs
!! NAME
!! interfaces_14wfs
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/14wfs
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

module interfaces_14wfs

 implicit none

interface
 subroutine bestwfs(gcc_block,ghc_block,gscc_block,gscc_calc,&  
  &  gvnlc_block,gvnlc_calc,istwf_k,mpi_enreg,nbdblock,npw_k,nspinor,nvectin,nvectout,wfoptalg)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: gscc_calc
  integer,intent(in) :: gvnlc_calc
  integer,intent(in) :: istwf_k
  integer,intent(in) :: nbdblock
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: nvectin
  integer,intent(in) :: nvectout
  integer,intent(in) :: wfoptalg
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: gcc_block(2,npw_k*nspinor,nbdblock)
  real(dp),intent(inout) :: ghc_block(2,npw_k*nspinor,nbdblock)
  real(dp),intent(inout) :: gscc_block(2,npw_k*nspinor,nbdblock*gscc_calc)
  real(dp),intent(inout) :: gvnlc_block(2,npw_k*nspinor,nbdblock)
 end subroutine bestwfs
end interface

interface
 subroutine envlop(cg,ecut,gmet,icgmod,kg,kpoint,mcg,nband,npw,nspinor)
  use defs_basis
  implicit none
  integer,intent(in) :: icgmod
  integer,intent(in) :: mcg
  integer,intent(in) :: nband
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  real(dp),intent(in) :: ecut
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg(3,npw)
  real(dp),intent(in) :: kpoint(3)
 end subroutine envlop
end interface

interface
 subroutine fxphas(cg,gsc,icg,igsc,istwfk,mcg,mgsc,mpi_enreg,nband_k,npw_k,useoverlap)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: istwfk
  integer,intent(in) :: mcg
  integer,intent(in) :: mgsc
  integer,intent(in) :: nband_k
  integer,intent(in) :: npw_k
  integer,intent(in) :: useoverlap
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(inout) :: gsc(2,mgsc*useoverlap)
 end subroutine fxphas
end interface

interface
 subroutine getghc(cwavef,dimffnl,ffnl,filstat,ghc,gsc,gs_ham,&  
  &  gvnlc,kg_k,kinpw,lambda,lmnmax,&  
  &  matblk,mgfft,mpi_enreg,mpsang,mpssoang,&  
  &  natom,ndat,npw,nspinor,ntypat,nvloc,n4,n5,n6,&  
  &  paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,type_calc,vlocal)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: dimffnl
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: natom
  integer,intent(in) :: ndat
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: nvloc
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: sij_opt
  integer,intent(in) :: tim_getghc
  integer,intent(in) :: type_calc
  character(len=fnlen),intent(in) :: filstat
  type(gs_hamiltonian_type),intent(in) :: gs_ham
  real(dp) :: lambda
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cwavef(2,npw*nspinor*ndat)
  real(dp),intent(in) :: ffnl(npw,dimffnl,lmnmax,ntypat)
  real(dp),intent(inout) :: ghc(2,npw*nspinor*ndat)
  real(dp),intent(out) :: gsc(2,npw*nspinor*ndat*(sij_opt+1)/2)
  real(dp),intent(inout) :: gvnlc(2,npw*nspinor*ndat)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kinpw(npw)
  real(dp),intent(inout) :: ph3d(2,npw,matblk)
  real(dp),intent(inout) :: vlocal(n4,n5,n6,nvloc)
 end subroutine getghc
end interface

interface
 subroutine listkk(dksqmax,gmet,indkk,kptns1,kptns2,nkpt1,nkpt2,nsym,&  
  &  sppoldbl,symafm,symrel,timrev)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpt1
  integer,intent(in) :: nkpt2
  integer,intent(in) :: nsym
  integer,intent(in) :: sppoldbl
  integer,intent(in) :: timrev
  real(dp),intent(out) :: dksqmax
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(out) :: indkk(nkpt2*sppoldbl,6)
  real(dp),intent(in) :: kptns1(3,nkpt1)
  real(dp),intent(in) :: kptns2(3,nkpt2)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine listkk
end interface

interface
 subroutine precon(cg,eval,istwf_k,kinpw,mpi_enreg,npw,nspinor,optekin,pcon,vect)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: istwf_k
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: optekin
  real(dp),intent(in) :: eval
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(2,npw*nspinor)
  real(dp),intent(in) :: kinpw(npw)
  real(dp),intent(inout) :: pcon(npw)
  real(dp),intent(inout) :: vect(2,npw*nspinor)
 end subroutine precon
end interface

interface
 subroutine precon2(cg,eval,blocksize,istwf_k,kinpw,mpi_enreg,npw,nspinor,optekin,ghc,vect,vectsize)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: blocksize
  integer,intent(in) :: istwf_k
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: optekin
  integer,intent(in) :: vectsize
  type(mpi_type) :: mpi_enreg
  real(dp) :: cg(vectsize,blocksize)
  real(dp) :: eval(blocksize,blocksize)
  real(dp) :: ghc(vectsize,blocksize)
  real(dp) :: kinpw(npw)
  real(dp) :: vect(vectsize,blocksize)
 end subroutine precon2
end interface

interface
 subroutine prep_fourwf(rhoaug,blocksize,cwavef,wfraug,gs_hamk,istwf_k,iblock,icall,kg_k_gather,&  
  &  mgfft,mpi_enreg,nbdblock,nband_k,dimtabs,npw_k,n4,n5,n6,occ_k,paral_kgb,wtk)
  use defs_basis
  use defs_datatypes
  implicit none
  integer :: blocksize
  integer :: dimtabs
  integer :: iblock
  integer :: icall
  integer :: istwf_k
  integer :: mgfft
  integer :: n4
  integer :: n5
  integer :: n6
  integer :: nband_k
  integer :: nbdblock
  integer :: npw_k
  integer :: paral_kgb
  type(gs_hamiltonian_type) :: gs_hamk
  type(mpi_type) :: mpi_enreg
  real(dp) :: wtk
  real(dp) :: cwavef(2,npw_k*blocksize)
  integer :: kg_k_gather(3,dimtabs)
  real(dp) :: occ_k(nband_k)
  real(dp) :: rhoaug(n4,n5,n6)
  real(dp) :: wfraug(2,n4,n5,n6)
 end subroutine prep_fourwf
end interface

interface
 subroutine prep_getghc(cwavef,dimffnl,dtfil,ffnl_gather,gs_hamk,gvnlc,gwavef,swavef,iblock,icall,istwf_k,kg_k_gather,&  
  &  kinpw_gather,lambda,lmnmax,matblk,blocksize,mgfft,mpi_enreg,mpsang,mpssoang,natom,nbdblock,nband_k,dimtabs,npw_k,&  
  &  nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,ph3d_gather,prtvol,sij_opt,vlocal)
  use defs_basis
  use defs_datatypes
  implicit none
  integer :: blocksize
  integer :: dimffnl
  integer :: dimtabs
  integer :: iblock
  integer :: icall
  integer :: istwf_k
  integer :: lmnmax
  integer :: matblk
  integer :: mgfft
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
  integer :: paral_kgb
  integer :: prtvol
  integer :: sij_opt
  type(datafiles_type) :: dtfil
  type(gs_hamiltonian_type) :: gs_hamk
  real(dp) :: lambda
  type(mpi_type) :: mpi_enreg
  real(dp) :: cwavef(2,npw_k*nspinor*blocksize)
  real(dp) :: ffnl_gather(dimtabs,dimffnl,lmnmax,ntypat)
  real(dp) :: gvnlc(2,npw_k*nspinor*blocksize)
  real(dp) :: gwavef(2,npw_k*nspinor*blocksize)
  integer :: kg_k_gather(3,dimtabs)
  real(dp) :: kinpw_gather(dimtabs)
  real(dp) :: ph3d_gather(2,dimtabs,matblk)
  real(dp) :: swavef(2,npw_k*nspinor*blocksize)
  real(dp) :: vlocal(n4,n5,n6,nvloc)
 end subroutine prep_getghc
end interface

interface
 subroutine prep_index_wavef_bandpp(nproc_band,bandpp,&  
  nspinor,ndatarecv,&  
  recvcounts,rdispls,&  
  index_wavef_band)
  implicit none
  integer,intent(in) :: bandpp
  integer,intent(in) :: ndatarecv
  integer,intent(in) :: nproc_band
  integer,intent(in) :: nspinor
  integer,pointer :: index_wavef_band(:)
  integer,intent(in) :: rdispls(nproc_band)
  integer,intent(in) :: recvcounts(nproc_band)
 end subroutine prep_index_wavef_bandpp
end interface

interface
 subroutine prep_kg_sym_do(mpi_enreg,&  
  kg_k_gather,ndatarecv,&  
  kg_k_gather_sym,ndatarecv_tot,&  
  ndatasend_sym,idatarecv0,&  
  tab_proc,&  
  sendcounts_sym,sendcounts_sym_all,sdispls_sym,&  
  recvcounts_sym,recvcounts_sym_tot,rdispls_sym)
  use defs_datatypes
  implicit none
  integer,intent(out) :: idatarecv0
  integer,intent(in) :: ndatarecv
  integer,intent(out) :: ndatarecv_tot
  integer,intent(out) :: ndatasend_sym
  type(mpi_type),intent(in) :: mpi_enreg
  integer,pointer :: kg_k_gather_sym(:,:)
  integer,pointer :: rdispls_sym(:)
  integer,pointer :: recvcounts_sym(:)
  integer,pointer :: recvcounts_sym_tot(:)
  integer,pointer :: sdispls_sym(:)
  integer,pointer :: sendcounts_sym(:)
  integer,pointer :: sendcounts_sym_all(:)
  integer,pointer :: tab_proc(:)
  integer,intent(in) :: kg_k_gather(3,ndatarecv)
 end subroutine prep_kg_sym_do
end interface

interface
 subroutine prep_wavef_sym_do(mpi_enreg,nproc_band,bandpp,nspinor,&  
  ndatarecv,recvcounts,rdispls,&  
  ndatarecv_tot,ndatasend_sym,tab_proc,&  
  cwavef_alltoall,&  
  sendcounts_sym,sendcounts_sym_all,sdispls_sym,&  
  recvcounts_sym,recvcounts_sym_tot,rdispls_sym,&  
  ewavef_alltoall_sym,&  
  index_wavef_send)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: bandpp
  integer,intent(in) :: ndatarecv
  integer,intent(in) :: ndatarecv_tot
  integer,intent(in) :: ndatasend_sym
  integer,intent(in) :: nproc_band
  integer,intent(in) :: nspinor
  type(mpi_type),intent(in) :: mpi_enreg
  integer,pointer :: index_wavef_send(:)
  integer,pointer :: rdispls_sym(:)
  integer,pointer :: recvcounts_sym(:)
  integer,pointer :: recvcounts_sym_tot(:)
  integer,pointer :: sdispls_sym(:)
  integer,pointer :: sendcounts_sym(:)
  integer,pointer :: sendcounts_sym_all(:)
  integer,pointer :: tab_proc(:)
  real(dp),intent(inout) :: cwavef_alltoall(2,ndatarecv*nspinor*bandpp)
  real(dp),pointer :: ewavef_alltoall_sym(:,:)
  integer,intent(in) :: rdispls(nproc_band)
  integer,intent(in) :: recvcounts(nproc_band)
 end subroutine prep_wavef_sym_do
end interface

interface
 subroutine prep_wavef_sym_undo(mpi_enreg,nproc_band,bandpp,nspinor,&  
  ndatarecv,recvcounts,rdispls,&  
  ndatarecv_tot,ndatasend_sym,idatarecv0,tab_proc,&  
  gwavef_alltoall,&  
  sendcounts_sym,sendcounts_sym_all,sdispls_sym,&  
  recvcounts_sym,recvcounts_sym_tot,rdispls_sym,&  
  gwavef_alltoall_sym,&  
  index_wavef_send)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: bandpp
  integer,intent(in) :: idatarecv0
  integer,intent(in) :: ndatarecv
  integer,intent(in) :: ndatarecv_tot
  integer,intent(in) :: ndatasend_sym
  integer,intent(in) :: nproc_band
  integer,intent(in) :: nspinor
  type(mpi_type),intent(in) :: mpi_enreg
  integer,pointer :: index_wavef_send(:)
  integer,pointer :: rdispls_sym(:)
  integer,pointer :: recvcounts_sym(:)
  integer,pointer :: recvcounts_sym_tot(:)
  integer,pointer :: sdispls_sym(:)
  integer,pointer :: sendcounts_sym(:)
  integer,pointer :: sendcounts_sym_all(:)
  integer,pointer :: tab_proc(:)
  real(dp),intent(inout) :: gwavef_alltoall(2,ndatarecv*nspinor*bandpp)
  real(dp),pointer :: gwavef_alltoall_sym(:,:)
  integer,intent(in) :: rdispls(nproc_band)
  integer,intent(in) :: recvcounts(nproc_band)
 end subroutine prep_wavef_sym_undo
end interface

interface
 subroutine projbd(cg,direc,iband0,icg,iscg,istwf_k,mcg,mpi_enreg,mscg,nband,&  
  &  npw,nspinor,ortalg,printopt,scg,scprod,tim_projbd,useoverlap)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iband0
  integer,intent(in) :: icg
  integer,intent(in) :: iscg
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mcg
  integer,intent(in) :: mscg
  integer,intent(in) :: nband
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: ortalg
  integer,intent(in) :: printopt
  integer,intent(in) :: tim_projbd
  integer,intent(in) :: useoverlap
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(inout) :: direc(2,npw*nspinor)
  real(dp),intent(in) :: scg(2,mscg*useoverlap)
  real(dp),intent(out) :: scprod(2,nband)
 end subroutine projbd
end interface

interface
 subroutine pw_orthon(icg,igsc,istwf_k,mcg,mgsc,mpi_enreg,nelem,nvec,&  
  &  ortalgo,ovl_vecnm,useoverlap,vecnm)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mcg
  integer,intent(in) :: mgsc
  integer,intent(in) :: nelem
  integer,intent(in) :: nvec
  integer,intent(in) :: ortalgo
  integer,intent(in) :: useoverlap
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: ovl_vecnm(2,mgsc*useoverlap)
  real(dp),intent(inout) :: vecnm(2,mcg)
 end subroutine pw_orthon
end interface

interface
 subroutine sdirot(cg,evec,icg,mcg,ndim,num,npw)
  use defs_basis
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: mcg
  integer,intent(in) :: ndim
  integer,intent(in) :: npw
  integer,intent(in) :: num
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(in) :: evec(2*ndim,num)
 end subroutine sdirot
end interface

interface
 subroutine wfconv(ceksp2,cg1,cg2,debug,ecut1,ecut2,ecut2_eff,&  
  &  eig_k1,eig_k2,exchn2n3d,formeig,gmet1,gmet2,icg1,icg2,idum,&  
  &  ikpt1,ikpt10,ikpt2,indkk,inplace,isppol2,istwfk1,istwfk2,&  
  &  kg1,kg2,kptns1,kptns2,mband1,mband2,mcg1,mcg2,mgfft,mpi_enreg,mpw1,mpw2,nbd1,nbd2,&  
  &  ngfft,nkpt1,nkpt2,npw1,npw2,nspinor1,nspinor2,nsym,&  
  &  occ_k1,occ_k2,optorth,restart,rprimd2,sppoldbl,symafm,symrel,tnons)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ceksp2
  integer,intent(in) :: debug
  integer,intent(in) :: exchn2n3d
  integer,intent(in) :: formeig
  integer,intent(in) :: icg1
  integer,intent(in) :: icg2
  integer,intent(in) :: idum
  integer,intent(in) :: ikpt1
  integer,intent(inout) :: ikpt10
  integer,intent(in) :: ikpt2
  integer,intent(in) :: inplace
  integer,intent(in) :: isppol2
  integer,intent(in) :: mband1
  integer,intent(in) :: mband2
  integer,intent(in) :: mcg1
  integer,intent(in) :: mcg2
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpw1
  integer,intent(in) :: mpw2
  integer,intent(in) :: nbd1
  integer,intent(in) :: nbd2
  integer,intent(in) :: nkpt1
  integer,intent(in) :: nkpt2
  integer,intent(inout) :: npw1
  integer,intent(inout) :: npw2
  integer,intent(in) :: nspinor1
  integer,intent(in) :: nspinor2
  integer,intent(in) :: nsym
  integer,intent(in) :: optorth
  integer,intent(in) :: restart
  integer,intent(in) :: sppoldbl
  real(dp),intent(in) :: ecut1
  real(dp),intent(in) :: ecut2
  real(dp),intent(in) :: ecut2_eff
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: cg1(2,mcg1)
  real(dp),intent(inout) :: cg2(2,mcg2)
  real(dp),intent(inout) :: eig_k1(mband1*(2*mband1)**formeig)
  real(dp),intent(inout) :: eig_k2(mband2*(2*mband2)**formeig)
  real(dp),intent(in) :: gmet1(3,3)
  real(dp),intent(in) :: gmet2(3,3)
  integer,intent(in) :: indkk(nkpt2*sppoldbl,6)
  integer,intent(in) :: istwfk1(nkpt1)
  integer,intent(in) :: istwfk2(nkpt2)
  integer,intent(inout) :: kg1(3,mpw1)
  integer,intent(inout) :: kg2(3,mpw2)
  real(dp),intent(in) :: kptns1(3,nkpt1)
  real(dp),intent(in) :: kptns2(3,nkpt2)
  real(dp),intent(inout) :: occ_k1(mband1)
  real(dp),intent(inout) :: occ_k2(mband2)
  real(dp),intent(in) :: rprimd2(3,3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
 end subroutine wfconv
end interface

interface
 subroutine zprecon3(cg,eval,blocksize,istwf_k,kinpw,mpi_enreg,npw,nspinor,optekin,ghc,vect,vectsize)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: blocksize
  integer,intent(in) :: istwf_k
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: optekin
  integer,intent(in) :: vectsize
  type(mpi_type) :: mpi_enreg
  complex(dpc) :: cg(vectsize,blocksize)
  complex(dpc) :: eval(blocksize,blocksize)
  complex(dpc) :: ghc(vectsize,blocksize)
  real(dp) :: kinpw(npw)
  complex(dpc) :: vect(vectsize,blocksize)
 end subroutine zprecon3
end interface

end module interfaces_14wfs
!!***
