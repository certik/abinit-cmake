!!****m* ABINIT/interfaces_13paw
!! NAME
!! interfaces_13paw
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/13paw
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

module interfaces_13paw

 implicit none

interface
 subroutine chkpawovlp(natom,ntypat,pawovlp,pawtab,rmet,typat,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp) :: pawovlp
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: rmet(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine chkpawovlp
end interface

interface
 function clp(x)
  use defs_basis
  implicit none
  real(dp) :: clp
  real(dp),intent(in) :: x
 end function clp
end interface

interface
 function dbeta(cosbeta,ll,mp,mm)
  use defs_basis
  implicit none
  integer,intent(in) :: ll
  integer,intent(in) :: mm
  integer,intent(in) :: mp
  real(dp),intent(in) :: cosbeta
  real(dp) :: dbeta
 end function dbeta
end interface

interface
 function gaunt(ll,mm,l1,m1,l2,m2)
  use defs_basis
  implicit none
  integer,intent(in) :: l1
  integer,intent(in) :: l2
  integer,intent(in) :: ll
  integer,intent(in) :: m1
  integer,intent(in) :: m2
  integer,intent(in) :: mm
  real(dp) :: gaunt
 end function gaunt
end interface

interface
 subroutine gipaw_aug_fields(gipaw_aug,natom,nfft,ntypat,pawfgrtab,pawrad,pawrhoij,pawtab,psps,typat)
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  type(pseudopotential_type),intent(in) :: psps
  type(gipaw_type),intent(out) :: gipaw_aug(natom)
  type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine gipaw_aug_fields
end interface

interface
 subroutine gipaw_j_dia_aug(gipaw_aug,jdia,natom,nfft,ntypat,pawfgrtab,pawrhoij,pawtab,typat)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  type(gipaw_type),intent(in) :: gipaw_aug(natom)
  real(dp),intent(out) :: jdia(3,3,nfft)
  type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine gipaw_j_dia_aug
end interface

interface
 subroutine indgrid(coatofin,fintocoa,nfftc,nfftf,ngfftc,ngfftf)
  implicit none
  integer,intent(in) :: nfftc
  integer,intent(in) :: nfftf
  integer,intent(in) :: ngfftc(18)
  integer,intent(in) :: ngfftf(18)
  integer,intent(out) :: coatofin(nfftc)
  integer,intent(out) :: fintocoa(nfftf)
 end subroutine indgrid
end interface

interface
 subroutine initang(pawang)
  use defs_datatypes
  implicit none
  type(pawang_type),intent(inout) :: pawang
 end subroutine initang
end interface

interface
 subroutine initrhoij(indlmn,lexexch,lmnmax,lpawu,natom,nspden,nsppol,ntypat,pawrhoij,pawtab,spinat,typat)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: lmnmax
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: lexexch(ntypat)
  integer,intent(in) :: lpawu(ntypat)
  type(pawrhoij_type),intent(out) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: spinat(3,natom)
  integer,intent(in) :: typat(natom)
 end subroutine initrhoij
end interface

interface
 subroutine initylmr(mpsang,normchoice,npts,nrm,option,rr,ylmr,ylmr_gr)
  use defs_basis
  implicit none
  integer,intent(in) :: mpsang
  integer,intent(in) :: normchoice
  integer,intent(in) :: npts
  integer,intent(in) :: option
  real(dp),intent(in) :: nrm(npts)
  real(dp),intent(in) :: rr(3,npts)
  real(dp),intent(out) :: ylmr(mpsang*mpsang,npts)
  real(dp),intent(out) :: ylmr_gr(3*(option/2)+6*(option/3),mpsang*mpsang,npts)
 end subroutine initylmr
end interface

interface
 subroutine int_ang(ang_phipphj,mpsang)
  use defs_basis
  implicit none
  integer,intent(in) :: mpsang
  real(dp),intent(out) :: ang_phipphj(mpsang**2,mpsang**2,8)
 end subroutine int_ang
end interface

interface
 subroutine make_cs_dia(cs,natom,ntypat,pawang,pawrhoij,pawrad,pawtab,typat)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(out) :: cs(3,3,natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine make_cs_dia
end interface

interface
 subroutine make_efg_paw(efg,natom,ntypat,pawang,pawrhoij,pawrad,pawtab,typat)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(out) :: efg(3,3,natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine make_efg_paw
end interface

interface
 subroutine make_fc_paw(fc,natom,ntypat,pawrhoij,pawrad,pawtab,psps,typat)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(out) :: fc(natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine make_fc_paw
end interface

interface
 subroutine mkeuler(rot,cosbeta,cosalp,sinalp,cosgam,singam,isn)
  use defs_basis
  implicit none
  integer,intent(out) :: isn
  real(dp),intent(out) :: cosalp
  real(dp),intent(out) :: cosbeta
  real(dp),intent(out) :: cosgam
  real(dp),intent(out) :: sinalp
  real(dp),intent(out) :: singam
  real(dp),intent(in) :: rot(3,3)
 end subroutine mkeuler
end interface

interface
 subroutine nhatgrid(atindx1,gmet,mpi_enreg,natom,nattyp,nfft,ngfft,ntypat,&  
  &  optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,rprimd,typat,ucvol,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: optgr0
  integer,intent(in) :: optgr1
  integer,intent(in) :: optgr2
  integer,intent(in) :: optrad
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: nattyp(ntypat)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine nhatgrid
end interface

interface
 subroutine optics_paw(atindx1,cg,cprj,dimcprj,dtfil,dtset,ecut,fildata,gprimd,hdr,indlmn,kg,lmnmax,&  
  &  mband,mkmem,mpi_enreg,mpsang,mpw,natom,nattyp,nkpt,npwarr,nspinor,nsppol,&  
  &  pawrad,pawtab,wffnow)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: lmnmax
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: ecut
  character(len=fnlen),intent(in) :: fildata
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(inout) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  type(cprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
  integer,intent(in) :: dimcprj(natom)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indlmn(6,lmnmax,dtset%ntypat)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: nattyp(dtset%ntypat)
  integer,intent(in) :: npwarr(nkpt)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)
 end subroutine optics_paw
end interface

interface
 subroutine partial_dos_fractions_paw(atindx1,cprj,dimcprj,dos_fractions,dos_fractions_m,&  
  &  dos_fractions_paw1,dos_fractions_pawt1,&  
  &  dtfil,dtset,indlmn,lmnmax,mbesslang,mkmem,&  
  &  mpi_enreg,m_dos_flag,ndosfraction,paw_dos_flag,pawrad,pawtab,prtdos)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: lmnmax
  integer,intent(in) :: m_dos_flag
  integer,intent(in) :: mbesslang
  integer,intent(in) :: mkmem
  integer,intent(in) :: ndosfraction
  integer,intent(in) :: paw_dos_flag
  integer,intent(in) :: prtdos
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: atindx1(dtset%natom)
  type(cprj_type) :: cprj(dtset%natom,dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  integer,intent(in) :: dimcprj(dtset%natom)
  real(dp),intent(inout) :: dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction)
  real(dp),intent(inout) :: dos_fractions_m(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*mbesslang*m_dos_flag)
  real(dp),intent(out) :: dos_fractions_paw1(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*paw_dos_flag)
  real(dp),intent(out) :: dos_fractions_pawt1(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*paw_dos_flag)
  integer,intent(in) :: indlmn(6,lmnmax,dtset%ntypat)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)
 end subroutine partial_dos_fractions_paw
end interface

interface
 subroutine pawaccrhoij(atindx1,cplex,cwaveprj,cwaveprj1,dimpaw,ipert,isppol,natom,nspden,&  
  &  nspinor,nsppol,occ_k,option,pawrhoij,wtk_k)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: dimpaw
  integer,intent(in) :: ipert
  integer,intent(in) :: isppol
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: option
  real(dp) :: occ_k
  real(dp) :: wtk_k
  integer,intent(in) :: atindx1(natom)
  type(cprj_type),intent(in) :: cwaveprj(dimpaw,nspinor)
  type(cprj_type),intent(in) :: cwaveprj1(dimpaw,nspinor)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dimpaw)
 end subroutine pawaccrhoij
end interface

interface
 subroutine pawalloc(dtset,idtset,mpsang,mqgrid_vl,npsp,option,paw_size,paw_size_old,&  
  &  pawang,pawrad,pawtab,pspheads)
  use defs_datatypes
  implicit none
  integer,intent(in) :: idtset
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid_vl
  integer,intent(in) :: npsp
  integer,intent(in) :: option
  integer,intent(in) :: paw_size
  integer,intent(in) :: paw_size_old
  type(dataset_type),intent(in) :: dtset
  type(pawang_type),intent(inout) :: pawang
  type(pawrad_type),intent(inout) :: pawrad(paw_size)
  type(pawtab_type),intent(inout) :: pawtab(paw_size)
  type(pspheader_type),intent(in) :: pspheads(npsp)
 end subroutine pawalloc
end interface

interface
 subroutine pawdenpot(compch_sph,epaw,epawdc,ixc,natom,nspden,ntypat,nzlmopt,option,paw_an,&  
  &  paw_ij,pawang,pawprtvol,pawrad,pawrhoij,pawspnorb,pawtab,pawxcdev,typat,xclevel,znucl)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ixc
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: nzlmopt
  integer,intent(in) :: option
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: pawspnorb
  integer,intent(in) :: pawxcdev
  integer,intent(in) :: xclevel
  real(dp),intent(out) :: compch_sph
  real(dp),intent(out) :: epaw
  real(dp),intent(out) :: epawdc
  type(pawang_type),intent(in) :: pawang
  type(paw_an_type),intent(inout) :: paw_an(natom)
  type(paw_ij_type),intent(inout) :: paw_ij(natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
  real(dp) :: znucl(ntypat)
 end subroutine pawdenpot
end interface

interface
 subroutine pawdij(dtset,enunit,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,&  
  &  paw_an,paw_ij,pawang,pawfgrtab,pawprtvol,pawrad,&  
  &  pawspnorb,pawtab,pawxcdev,typat,ucvol,vtrial,vxc)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: enunit
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: pawspnorb
  integer,intent(in) :: pawxcdev
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  type(paw_an_type),intent(in) :: paw_an(natom)
  type(paw_ij_type),intent(inout) :: paw_ij(natom)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vtrial(nfft,nspden)
  real(dp),intent(in) :: vxc(nfft,nspden)
 end subroutine pawdij
end interface

interface
 subroutine pawgrnl(atindx1,dimnhat,dimvtrial,dyfrnl,grnl,mpi_enreg,natom,nattyp,nfft,ngfft,&  
  &  nhat,nlstr,nspden,nsym,ntypat,optgr,optgr2,optstr,pawang,pawfgrtab,&  
  &  pawrhoij,pawtab,rprimd,symrec,typat,vtrial)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: dimnhat
  integer,intent(in) :: dimvtrial
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: optgr
  integer,intent(in) :: optgr2
  integer,intent(in) :: optstr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(inout) :: dyfrnl(3,3,natom*optgr2)
  real(dp),intent(inout) :: grnl(3*natom*optgr)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(nfft,dimnhat)
  real(dp),intent(inout) :: nlstr(6*optstr)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vtrial(nfft,dimvtrial)
 end subroutine pawgrnl
end interface

interface
 subroutine pawgylm(gylm,gylmgr,gylmgr2,iatom,ifftsph,itypat,lm_size,nfgd,optgr0,optgr1,optgr2,&  
  &  pawtab,rfgd,rfgd_allocated)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iatom
  integer,intent(in) :: itypat
  integer,intent(in) :: lm_size
  integer,intent(in) :: nfgd
  integer,intent(in) :: optgr0
  integer,intent(in) :: optgr1
  integer,intent(in) :: optgr2
  integer,intent(in) :: rfgd_allocated
  type(pawtab_type),intent(in) :: pawtab
  real(dp),intent(out) :: gylm(nfgd,optgr0*lm_size)
  real(dp),intent(out) :: gylmgr(3,nfgd,optgr1*lm_size)
  real(dp),intent(out) :: gylmgr2(6,nfgd,optgr2*lm_size)
  integer,intent(in) :: ifftsph(nfgd)
  real(dp),intent(in) :: rfgd(3,nfgd)
 end subroutine pawgylm
end interface

interface
 subroutine pawgylmg(gprimd,gylmg,kg,kpg,kpt,lmax,nkpg,npw,ntypat,pawtab,ylm)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: lmax
  integer,intent(in) :: nkpg
  integer,intent(in) :: npw
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: gylmg(npw,lmax**2,ntypat)
  integer,intent(in) :: kg(3,npw)
  real(dp),intent(in) :: kpg(npw,nkpg)
  real(dp),intent(in) :: kpt(3)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ylm(npw,lmax**2)
 end subroutine pawgylmg
end interface

interface
 subroutine pawinit(ecutshp_eff,indlmn,lcutdens,lmix,lmnmax,mpsang,n1xccc,nphi,nsym,ntheta,ntypat,&  
  &  pawang,pawrad,pawspnorb,pawtab,pawxcdev)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: lcutdens
  integer,intent(in) :: lmix
  integer,intent(in) :: lmnmax
  integer,intent(in) :: mpsang
  integer,intent(inout) :: n1xccc
  integer,intent(in) :: nphi
  integer,intent(in) :: nsym
  integer,intent(in) :: ntheta
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawspnorb
  integer,intent(in) :: pawxcdev
  real(dp) :: ecutshp_eff
  type(pawang_type),intent(inout) :: pawang
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(inout) :: pawtab(ntypat)
 end subroutine pawinit
end interface

interface
 subroutine pawlsylm(pawang)
  use defs_datatypes
  implicit none
  type(pawang_type),intent(inout) :: pawang
 end subroutine pawlsylm
end interface

interface
 subroutine pawmknhat(compch_fft,ider,izero,mpi_enreg,natom,nfft,ngfft,nhatgrdim,nspden,&  
  &  ntypat,paral_kgb,pawang,pawfgrtab,pawgrnhat,pawnhat,pawrhoij,pawtab,typat,ucvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ider
  integer,intent(in) :: izero
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  real(dp),intent(out) :: compch_fft
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
  real(dp),intent(out) :: pawgrnhat(nfft,nspden,3*nhatgrdim)
  real(dp),intent(out) :: pawnhat(nfft,nspden)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine pawmknhat
end interface

interface
 subroutine pawmknhat3(cplex,idir,ipert,izero,mpi_enreg,natom,nfft,ngfft,nhat1dim,nrhoij1,nspden,&  
  &  ntypat,paral_kgb,pawang,pawfgrtab,pawnhat1,pawrhoij,pawrhoij1,pawtab,typat)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: izero
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nhat1dim
  integer,intent(in) :: nrhoij1
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: ngfft(18)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
  real(dp),intent(out) :: pawnhat1(cplex*nhat1dim,nspden)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawrhoij_type),intent(in) :: pawrhoij1(nrhoij1)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine pawmknhat3
end interface

interface
 subroutine pawmkrhoij(atindx1,cprj,dimcprj,istwfk,mband,mkmem,mpi_enreg,natom,nattyp,&  
  &  nband,nkpt,nspden,nspinor,nsppol,ntypat,occ,pawprtvol,pawrhoij,unpaw,wtk)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: unpaw
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: atindx1(natom)
  type(cprj_type),intent(in) :: cprj(natom,nspinor*mband*mkmem*nsppol)
  integer,intent(in) :: dimcprj(natom)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(pawrhoij_type),intent(inout) :: pawrhoij(natom)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine pawmkrhoij
end interface

interface
 subroutine pawmkrhoij3(atindx1,cplex,cprj,cprj1,dimcprj,dimpaw1,ipert,istwfk,mband,mkmem,mk1mem,&  
  &  mpi_enreg,natom,nattyp,nband,nkpt,nspden,nspinor,nsppol,ntypat,occ,pawprtvol,&  
  &  pawrhoij1,unpaw,unpaw1,usecprj,wtk)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: dimpaw1
  integer,intent(in) :: ipert
  integer,intent(in) :: mband
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: unpaw
  integer,intent(in) :: unpaw1
  integer,intent(in) :: usecprj
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: atindx1(natom)
  type(cprj_type),intent(in) :: cprj(dimpaw1,nspinor*mband*mkmem*nsppol*usecprj)
  type(cprj_type),intent(in) :: cprj1(dimpaw1,nspinor*mband*mk1mem*nsppol)
  integer,intent(in) :: dimcprj(natom)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(pawrhoij_type),intent(inout) :: pawrhoij1(dimpaw1)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine pawmkrhoij3
end interface

interface
 subroutine pawnabla_init(mpsang,lmnmax,ntypat,indlmn,pawrad,pawtab)
  use defs_datatypes
  implicit none
  integer,intent(in) :: lmnmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: ntypat
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(inout) :: pawtab(ntypat)
 end subroutine pawnabla_init
end interface

interface
 subroutine pawpolev(natom,nspden,ntypat,pawprtvol,pawrhoij,pawtab,pelev,typat)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawprtvol
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(out) :: pelev(3)
  integer,intent(in) :: typat(natom)
 end subroutine pawpolev
end interface

interface
 subroutine pawprt(indlmn,enunit,lmnmax,natom,ntypat,paw_ij,pawprtvol,pawrhoij,pawtab,typat)
  use defs_datatypes
  implicit none
  integer,intent(in) :: enunit
  integer,intent(in) :: lmnmax
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  type(paw_ij_type),intent(inout) :: paw_ij(natom)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine pawprt
end interface

interface
 subroutine pawpupot(ispden,nspden,paw_ij,pawprtvol,pawtab,vpawu,VUKS)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ispden
  integer,intent(in) :: nspden
  integer,intent(in) :: pawprtvol
  real(dp),intent(inout) :: VUKS
  type(paw_ij_type),intent(in) :: paw_ij
  type(pawtab_type),intent(in) :: pawtab
  real(dp),intent(inout) :: vpawu(pawtab%lpawu*2+1,pawtab%lpawu*2+1)
 end subroutine pawpupot
end interface

interface
 subroutine pawpuxinit(dmatpuopt,exchmix,jpawu,llexexch,llpawu,indlmn,lmnmax,ntypat,pawang,&  
  &  pawprtvol,pawrad,pawtab,upawu,useexexch,usepawu)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: dmatpuopt
  integer,intent(in) :: lmnmax
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: useexexch
  integer,intent(in) :: usepawu
  real(dp),intent(in) :: exchmix
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  real(dp),intent(in) :: jpawu(ntypat)
  integer,intent(in) :: llexexch(ntypat)
  integer,intent(in) :: llpawu(ntypat)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(inout) :: pawtab(ntypat)
  real(dp),intent(in) :: upawu(ntypat)
 end subroutine pawpuxinit
end interface

interface
 subroutine pawshpfun(ll,mesh,norm,pawtab,shapefunc)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ll
  type(pawrad_type),intent(in) :: mesh
  real(dp),intent(out) :: norm
  type(pawtab_type),intent(in) :: pawtab
  real(dp),intent(inout) :: shapefunc(mesh%mesh_size)
 end subroutine pawshpfun
end interface

interface
 subroutine pawsushat(atindx1,cprj_k,gbound_diel,gylmg_diel,iband1,iband2,istwf_k,kg_diel,&  
  &  lmax_diel,mgfftdiel,mpi_enreg,natom,nband,ndiel4,ndiel5,ndiel6,nfftdiel,&  
  &  ngfftdiel,npwdiel,nspinor,ntypat,optreal,paral_kgb,&  
  &  pawang,pawtab,ph3d_diel,typat,wfprod,wfraug)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iband1
  integer,intent(in) :: iband2
  integer,intent(in) :: istwf_k
  integer,intent(in) :: lmax_diel
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: ndiel4
  integer,intent(in) :: ndiel5
  integer,intent(in) :: ndiel6
  integer,intent(in) :: nfftdiel
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: optreal
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: ngfftdiel(18)
  integer,intent(in) :: atindx1(natom)
  type(cprj_type),intent(in) :: cprj_k(natom,nspinor*nband)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  real(dp),intent(in) :: gylmg_diel(npwdiel,lmax_diel**2,ntypat)
  integer,intent(in) :: kg_diel(3,npwdiel)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ph3d_diel(2,npwdiel,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: wfprod(2,npwdiel*(1-optreal))
  real(dp),intent(inout) :: wfraug(2,ndiel4,ndiel5,ndiel6*optreal)
 end subroutine pawsushat
end interface

interface
 subroutine pawuenergy(iatom,eldaumdc,eldaumdcdc,pawprtvol,pawtab,paw_ij,nspden)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iatom
  integer,intent(in) :: nspden
  integer,intent(in) :: pawprtvol
  real(dp),intent(inout) :: eldaumdc
  real(dp),intent(inout) :: eldaumdcdc
  type(paw_ij_type),intent(in) :: paw_ij
  type(pawtab_type),intent(in) :: pawtab
 end subroutine pawuenergy
end interface

interface
 subroutine pawxc(corexc,enxc,enxcdc,ixc,lm_size,lmselect,nhat,nspden,option,&  
  &  pawang,pawrad,rhor,usecore,usexcnhat,vxc,xclevel)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ixc
  integer,intent(in) :: lm_size
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: usecore
  integer,intent(in) :: usexcnhat
  integer,intent(in) :: xclevel
  real(dp),intent(out) :: enxc
  real(dp),intent(out) :: enxcdc
  type(pawang_type),intent(in) :: pawang
  type(pawrad_type),intent(in) :: pawrad
  real(dp),intent(in) :: corexc(pawrad%mesh_size)
  logical,intent(in) :: lmselect(lm_size)
  real(dp),intent(in) :: nhat(pawrad%mesh_size,lm_size,nspden)
  real(dp),intent(in) :: rhor(pawrad%mesh_size,lm_size,nspden)
  real(dp),intent(out) :: vxc(pawrad%mesh_size,pawang%angl_size,nspden)
 end subroutine pawxc
end interface

interface
 subroutine pawxcm(corexc,enxc,enxcdc,exexch,ixc,lm_size,lmselect,nhat,nspden,option,&  
  &  pawang,pawrad,pawxcdev,rhor,usecore,usexcnhat,vxc,xclevel)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: exexch
  integer,intent(in) :: ixc
  integer,intent(in) :: lm_size
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: pawxcdev
  integer,intent(in) :: usecore
  integer,intent(in) :: usexcnhat
  integer,intent(in) :: xclevel
  real(dp),intent(out) :: enxc
  real(dp),intent(out) :: enxcdc
  type(pawang_type),intent(in) :: pawang
  type(pawrad_type),intent(in) :: pawrad
  real(dp),intent(in) :: corexc(pawrad%mesh_size)
  logical,intent(in) :: lmselect(lm_size)
  real(dp),intent(in) :: nhat(pawrad%mesh_size,lm_size,nspden)
  real(dp),intent(in) :: rhor(pawrad%mesh_size,lm_size,nspden)
  real(dp),intent(out) :: vxc(pawrad%mesh_size,lm_size,nspden)
 end subroutine pawxcm
end interface

interface
 subroutine pawxcsph(exc,exexch,ixc,ndvxc,ngr2,ngrad,nrad,nspden,&  
  &  nspgrad,nvxcdgr,order,pawrad,rho_updn,vxc,xclevel)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: exexch
  integer,intent(in) :: ixc
  integer,intent(in) :: ndvxc
  integer,intent(in) :: ngr2
  integer,intent(in) :: ngrad
  integer,intent(in) :: nrad
  integer,intent(in) :: nspden
  integer,intent(in) :: nspgrad
  integer,intent(in) :: nvxcdgr
  integer,intent(in) :: order
  integer,intent(in) :: xclevel
  type(pawrad_type),intent(in) :: pawrad
  real(dp),intent(out) :: exc(nrad)
  real(dp),intent(in) :: rho_updn(nrad,nspden)
  real(dp),intent(out) :: vxc(nrad,nspden)
 end subroutine pawxcsph
end interface

interface
 subroutine pawxenergy(eexex,eexexdc,pawprtvol,pawrhoij,pawtab,nspden)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nspden
  integer,intent(in) :: pawprtvol
  real(dp),intent(inout) :: eexex
  real(dp),intent(inout) :: eexexdc
  type(pawrhoij_type),intent(in) :: pawrhoij
  type(pawtab_type),intent(in) :: pawtab
 end subroutine pawxenergy
end interface

interface
 subroutine pawxpot(nspden,pawprtvol,pawtab,paw_ij, pawrhoij)
  use defs_datatypes
  implicit none
  integer,intent(in) :: nspden
  integer,intent(in) :: pawprtvol
  type(paw_ij_type),intent(inout) :: paw_ij
  type(pawrhoij_type),intent(in) :: pawrhoij
  type(pawtab_type),intent(in) :: pawtab
 end subroutine pawxpot
end interface

interface
 function permutations(nn,kk)
  use defs_basis
  implicit none
  integer,intent(in) :: kk
  integer,intent(in) :: nn
  real(dp) :: permutations
 end function permutations
end interface

interface
 function phim(costheta,sintheta,mm)
  use defs_basis
  implicit none
  integer,intent(in) :: mm
  real(dp),intent(in) :: costheta
  real(dp) :: phim
  real(dp),intent(in) :: sintheta
 end function phim
end interface

interface
 subroutine realgaunt(l_max,ngnt,gntselect,realgnt)
  use defs_basis
  implicit none
  integer,intent(in) :: l_max
  integer,intent(out) :: ngnt
  integer,intent(out) :: gntselect((2*l_max-1)**2,l_max**2*(l_max**2+1)/2)
  real(dp),intent(out) :: realgnt((2*l_max-1)**2*(l_max)**4)
 end subroutine realgaunt
end interface

interface
 subroutine setnoccmmp(compute_dmat,dimdmat,dmatpawu,dmatudiag,impose_dmat,indsym,natom,natpawu,&  
  &  nspden,nspinor,nsppol,nsym,ntypat,paw_ij,pawang,pawprtvol,pawrhoij,pawtab,&  
  &  spinat,symafm,typat,useexexch,usepawu)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: compute_dmat
  integer,intent(in) :: dimdmat
  integer,intent(in) :: dmatudiag
  integer,intent(in) :: impose_dmat
  integer,intent(in) :: natom
  integer,intent(in) :: natpawu
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: useexexch
  integer,intent(in) :: usepawu
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: dmatpawu(dimdmat,dimdmat,nspinor*nsppol,natpawu*impose_dmat)
  integer,intent(in) :: indsym(4,nsym,natom)
  type(paw_ij_type),intent(inout) :: paw_ij(natom)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: spinat(3,natom)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: typat(natom)
 end subroutine setnoccmmp
end interface

interface
 subroutine setsymrhoij(gprimd,lmax,nsym,pawprtvol,rprimd,symafm,symrec,zarot)
  use defs_basis
  implicit none
  integer,intent(in) :: lmax
  integer,intent(in) :: nsym
  integer,intent(in) :: pawprtvol
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  real(dp),intent(out) :: zarot(2*lmax+1,2*lmax+1,lmax+1,nsym)
 end subroutine setsymrhoij
end interface

interface
 subroutine simple_j_dia(jdia,natom,nfft,pawfgrtab)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  real(dp),intent(out) :: jdia(3,3,nfft)
  type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)
 end subroutine simple_j_dia
end interface

interface
 subroutine smatrix_paw(cg,cgq,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,job,maxbd,mband,&  
  &  mcg_k,mcg_q,mcg1_k,minbd,mpw,nband_occ,npw_k1,npw_k2,nspinor,nsppol,&  
  &  pwind_k,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_k)
  use defs_basis
  implicit none
  integer,intent(in) :: ddkflag
  integer,intent(in) :: icg
  integer,intent(in) :: icg1
  integer,intent(in) :: itrs
  integer,intent(in) :: job
  integer,intent(in) :: maxbd
  integer,intent(in) :: mband
  integer,intent(in) :: mcg1_k
  integer,intent(in) :: mcg_k
  integer,intent(in) :: mcg_q
  integer,intent(in) :: minbd
  integer,intent(in) :: mpw
  integer,intent(in) :: nband_occ
  integer,intent(in) :: npw_k1
  integer,intent(in) :: npw_k2
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: shiftbd
  real(dp),intent(in) :: cg(2,mcg_k)
  real(dp),intent(out) :: cg1_k(2,mcg1_k)
  real(dp),intent(in) :: cgq(2,mcg_q)
  real(dp),intent(out) :: dtm_k(2)
  integer,intent(in) :: pwind_k(mpw)
  real(dp),intent(in) :: pwnsfac_k(4,mpw)
  integer,intent(inout) :: sflag_k(nband_occ)
  real(dp),intent(out) :: smat_inv(2,nband_occ,nband_occ)
  real(dp),intent(inout) :: smat_k(2,nband_occ,nband_occ)
 end subroutine smatrix_paw
end interface

interface
 subroutine smatrix_pawinit(cm2,cprj,ikpt1,ikpt2,&  
  &  g1,gprimd,kpt,mband,mbandw,mkmem,&  
  &  natom,nattyp,nband,&  
  &  nkpt,nspinor,nsppol,ntypat,pawang,pawrad,pawtab,rprimd,&  
  &  typat,usepaw,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ikpt1
  integer,intent(in) :: ikpt2
  integer,intent(in) :: mband
  integer,intent(in) :: mbandw
  integer,intent(in) :: mkmem
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  type(pawang_type),intent(in) :: pawang
  integer :: g1(3)
  real(dp),intent(inout) :: cm2(2,mbandw,mbandw)
  type(cprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: nband(nsppol*nkpt)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine smatrix_pawinit
end interface

interface
 subroutine spline_paw_fncs(dphi,dtphi,nnl,npts,pawrad,pawtab,points,phi,tphi)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nnl
  integer,intent(in) :: npts
  type(pawrad_type),intent(in) :: pawrad
  type(pawtab_type),intent(in) :: pawtab
  real(dp),intent(out) :: dphi(npts,nnl)
  real(dp),intent(out) :: dtphi(npts,nnl)
  real(dp),intent(out) :: phi(npts,nnl)
  real(dp),intent(in) :: points(npts)
  real(dp),intent(out) :: tphi(npts,nnl)
 end subroutine spline_paw_fncs
end interface

interface
 subroutine symdij(indlmn,indsym,lmnmax,natom,nsym,ntypat,&  
  &  paw_ij,pawang,pawprtvol,symafm,symrec,typat)
  use defs_datatypes
  implicit none
  integer,intent(in) :: lmnmax
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawprtvol
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: indsym(4,nsym,natom)
  type(paw_ij_type),intent(inout) :: paw_ij(natom)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
 end subroutine symdij
end interface

interface
 subroutine symrhoij(choice,indlmn,indsym,lmnmax,natom,nsym,ntypat,optrhoij,&  
  &  pawang,pawprtvol,pawrhoij,symafm,symrec,typat)
  use defs_datatypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: lmnmax
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: optrhoij
  integer,intent(in) :: pawprtvol
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: indsym(4,nsym,natom)
  type(pawrhoij_type),intent(inout) :: pawrhoij(natom)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
 end subroutine symrhoij
end interface

interface
 subroutine symrhoij3(indlmn,indsy1,ipert,lmnmax,natom,nrhoij1,nsym1,ntypat,&  
  &  pawang,pawprtvol,pawrhoij1,symaf1,symrc1,typat)
  use defs_datatypes
  implicit none
  integer,intent(in) :: ipert
  integer,intent(in) :: lmnmax
  integer,intent(in) :: natom
  integer,intent(in) :: nrhoij1
  integer,intent(in) :: nsym1
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawprtvol
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: indsy1(4,nsym1,natom)
  type(pawrhoij_type),intent(inout) :: pawrhoij1(nrhoij1)
  integer,intent(in) :: symaf1(nsym1)
  integer,intent(in) :: symrc1(3,3,nsym1)
  integer,intent(in) :: typat(natom)
 end subroutine symrhoij3
end interface

interface
 subroutine transgrid(cplex,mpi_enreg,nspden,optgrid,optin,optout,paral_kgb,pawfgr,rhog,rhogf,rhor,rhorf)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nspden
  integer,intent(in) :: optgrid
  integer,intent(in) :: optin
  integer,intent(in) :: optout
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawfgr_type),intent(in) :: pawfgr
  real(dp),intent(inout) :: rhog(2,pawfgr%nfftc)
  real(dp),intent(inout) :: rhogf(2,pawfgr%nfft)
  real(dp),intent(inout) :: rhor(cplex*pawfgr%nfftc,nspden)
  real(dp),intent(inout) :: rhorf(cplex*pawfgr%nfft,nspden)
 end subroutine transgrid
end interface

end module interfaces_13paw
!!***
