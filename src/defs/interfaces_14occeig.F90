!!****m* ABINIT/interfaces_14occeig
!! NAME
!! interfaces_14occeig
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/14occeig
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

module interfaces_14occeig

 implicit none

interface
 subroutine dens_in_sph(cmax,cg,gmet,istwfk,kg_k,natom,ngfft,mpi_enreg,npw_k,&  
  &  paral_kgb,ph1d,rmax,ucvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: istwfk
  integer,intent(in) :: natom
  integer,intent(in) :: npw_k
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: cg(2,npw_k)
  real(dp),intent(out) :: cmax(natom)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: ph1d(2,(2*ngfft(1)+1+2*ngfft(2)+1+2*ngfft(3)+1)*natom)
  real(dp),intent(in) :: rmax(natom)
 end subroutine dens_in_sph
end interface

interface
 subroutine dos_hdr_write(buffer,deltaene,dosdeltae,&  
  &  eigen,enemax,enemin,fermie,mband,nband,nene,&  
  &  nkpt,nsppol,occopt,prtdos,tphysel,tsmear,unitdos)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(out) :: nene
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: prtdos
  integer,intent(in) :: unitdos
  real(dp),intent(in) :: buffer
  real(dp),intent(out) :: deltaene
  real(dp),intent(in) :: dosdeltae
  real(dp),intent(out) :: enemax
  real(dp),intent(out) :: enemin
  real(dp),intent(in) :: fermie
  real(dp),intent(in) :: tphysel
  real(dp),intent(in) :: tsmear
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine dos_hdr_write
end interface

interface
 subroutine get_dos_1band (dos_fractions,enemin,enemax,&  
  &  integ_dos,nene,nkpt,ndosfraction,&  
  &  partial_dos,tweight,dtweightde)
  use defs_basis
  implicit none
  integer,intent(in) :: ndosfraction
  integer,intent(in) :: nene
  integer,intent(in) :: nkpt
  real(dp),intent(in) :: enemax
  real(dp),intent(in) :: enemin
  real(dp),intent(in) :: dos_fractions(nkpt,ndosfraction)
  real(dp),intent(in) :: dtweightde(nkpt,nene)
  real(dp),intent(out) :: integ_dos(nene,ndosfraction)
  real(dp),intent(out) :: partial_dos(nene,ndosfraction)
  real(dp),intent(in) :: tweight(nkpt,nene)
 end subroutine get_dos_1band
end interface

interface
 subroutine get_dos_1band_m (dos_fractions_m,enemin,enemax,&  
  &  integ_dos_m,nene,nkpt,ndosfraction_m,&  
  &  partial_dos_m,tweight,dtweightde)
  use defs_basis
  implicit none
  integer,intent(in) :: ndosfraction_m
  integer,intent(in) :: nene
  integer,intent(in) :: nkpt
  real(dp),intent(in) :: enemax
  real(dp),intent(in) :: enemin
  real(dp),intent(in) :: dos_fractions_m(nkpt,ndosfraction_m)
  real(dp),intent(in) :: dtweightde(nkpt,nene)
  real(dp),intent(out) :: integ_dos_m(nene,ndosfraction_m)
  real(dp),intent(out) :: partial_dos_m(nene,ndosfraction_m)
  real(dp),intent(in) :: tweight(nkpt,nene)
 end subroutine get_dos_1band_m
end interface

interface
 subroutine get_fsurf_1band(dtset,eigen_in,fermie,klatt,kpt_fullbz,&  
  &  mtetra,nfiner,nkpt_fullbz,ntetra,tetra_full,tetra_mult,tetra_wrap,tolfermi)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mtetra
  integer,intent(in) :: nfiner
  integer,intent(in) :: nkpt_fullbz
  integer,intent(in) :: ntetra
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: fermie
  real(dp),intent(in) :: tolfermi
  real(dp),intent(in) :: eigen_in(dtset%nkpt)
  real(dp),intent(in) :: klatt(3,3)
  real(dp),intent(in) :: kpt_fullbz(3,nkpt_fullbz)
  integer,intent(in) :: tetra_full(4,2,mtetra)
  integer,intent(in) :: tetra_mult(mtetra)
  integer,intent(in) :: tetra_wrap(3,4,mtetra)
 end subroutine get_fsurf_1band
end interface

interface
 subroutine get_tetra_weight(eigen_in,enemin,enemax,&  
  &  max_occ,mtetra,nene,nkpt,ntetra,rcvol,tetra_full,&  
  &  tetra_mult,tweight,dtweightde,vv)
  use defs_basis
  implicit none
  integer,intent(in) :: mtetra
  integer,intent(in) :: nene
  integer,intent(in) :: nkpt
  integer,intent(in) :: ntetra
  real(dp),intent(in) :: enemax
  real(dp),intent(in) :: enemin
  real(dp),intent(in) :: max_occ
  real(dp),intent(in) :: rcvol
  real(dp),intent(in) :: vv
  real(dp),intent(out) :: dtweightde(nkpt,nene)
  real(dp),intent(in) :: eigen_in(nkpt)
  integer,intent(in) :: tetra_full(4,2,mtetra)
  integer,intent(in) :: tetra_mult(mtetra)
  real(dp),intent(out) :: tweight(nkpt,nene)
 end subroutine get_tetra_weight
end interface

interface
 subroutine getnel(doccde,dosdeltae,eigen,entropy,fermie,maxocc,mband,nband,&  
  &  nelect,nkpt,nsppol,occ,occopt,option,tphysel,tsmear,unitdos,wtk)
  use defs_basis
  implicit none
  integer :: mband
  integer :: nkpt
  integer :: nsppol
  integer :: occopt
  integer :: option
  integer :: unitdos
  real(dp) :: dosdeltae
  real(dp),intent(out) :: entropy
  real(dp) :: fermie
  real(dp) :: maxocc
  real(dp),intent(out) :: nelect
  real(dp) :: tphysel
  real(dp) :: tsmear
  real(dp),intent(out) :: doccde(mband*nkpt*nsppol)
  real(dp) :: eigen(mband*nkpt*nsppol)
  integer :: nband(nkpt*nsppol)
  real(dp),intent(out) :: occ(mband*nkpt*nsppol)
  real(dp) :: wtk(nkpt)
 end subroutine getnel
end interface

interface
 subroutine init_bess_spl(mbess,bessargmax,bessint_delta,mlang,&  
  &  bess_spl,bess_spl_der,x_bess)
  use defs_basis
  implicit none
  integer,intent(in) :: mbess
  integer,intent(in) :: mlang
  real(dp),intent(in) :: bessargmax
  real(dp),intent(in) :: bessint_delta
  real(dp),intent(out) :: bess_spl(mbess,mlang)
  real(dp),intent(out) :: bess_spl_der(mbess,mlang)
  real(dp),intent(out) :: x_bess(mbess)
 end subroutine init_bess_spl
end interface

interface
 subroutine init_ylm_spl(mbessint,bessargmax,bessint_delta,mlang,spl_bessint)
  use defs_basis
  implicit none
  integer,intent(in) :: mbessint
  integer,intent(in) :: mlang
  real(dp),intent(in) :: bessargmax
  real(dp),intent(in) :: bessint_delta
  real(dp),intent(out) :: spl_bessint(mbessint,mlang)
 end subroutine init_ylm_spl
end interface

interface
 subroutine newocc(doccde,eigen,entropy,fermie,fixmom,mband,nband,&  
  &  nelect,nkpt,nspinor,nsppol,occ,occopt,prtvol,stmbias,tphysel,tsmear,wtk)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: prtvol
  real(dp),intent(out) :: entropy
  real(dp),intent(out) :: fermie
  real(dp),intent(in) :: fixmom
  real(dp),intent(in) :: nelect
  real(dp),intent(in) :: stmbias
  real(dp),intent(in) :: tphysel
  real(dp),intent(in) :: tsmear
  real(dp),intent(out) :: doccde(mband*nkpt*nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(out) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine newocc
end interface

interface
 subroutine occeig(doccde_k,doccde_kq,eig0_k,eig0_kq,nband_k,&  
  &  occopt,occ_k,occ_kq,rocceig)
  use defs_basis
  implicit none
  integer,intent(in) :: nband_k
  integer,intent(in) :: occopt
  real(dp),intent(in) :: doccde_k(nband_k)
  real(dp),intent(in) :: doccde_kq(nband_k)
  real(dp),intent(in) :: eig0_k(nband_k)
  real(dp),intent(in) :: eig0_kq(nband_k)
  real(dp),intent(in) :: occ_k(nband_k)
  real(dp),intent(in) :: occ_kq(nband_k)
  real(dp),intent(out) :: rocceig(nband_k,nband_k)
 end subroutine occeig
end interface

interface
 subroutine pareigocc(eigen,formeig,localrdwf,mpi_enreg,mband,nband,nkpt,nsppol,occ)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: localrdwf
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: eigen(mband*(2*mband)**formeig*nkpt*nsppol)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
 end subroutine pareigocc
end interface

interface
 subroutine partial_dos_fractions(cg,dos_fractions,dos_fractions_m,dtfil,&  
  &  dtset,hdr,mbesslang,mpi_enreg,m_dos_flag,ndosfraction,partial_dos,wffnow)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: m_dos_flag
  integer,intent(in) :: mbesslang
  integer,intent(in) :: ndosfraction
  integer,intent(in) :: partial_dos
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wffnow
  real(dp),intent(inout) :: cg(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp),intent(out) :: dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction)
  real(dp),intent(out) :: dos_fractions_m(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*mbesslang*m_dos_flag)
 end subroutine partial_dos_fractions
end interface

interface
 subroutine printbxsf(eigen,ewind,fermie,gprimd,kptrlatt,mband,&  
  &  nkptirred,kptirred,nsym,symrec,timrev,nsppol,shiftk,nshiftk,fname)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkptirred
  integer,intent(in) :: nshiftk
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  real(dp),intent(in) :: ewind
  real(dp),intent(in) :: fermie
  character(len=fnlen) :: fname
  integer,intent(in) :: kptrlatt(3,3)
  real(dp),intent(in) :: eigen(mband,nkptirred,nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: kptirred(3,nkptirred)
  real(dp),intent(in) :: shiftk(3,nshiftk)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine printbxsf
end interface

interface
 subroutine prtxcfermsurf(eigen,fermie,gprimd,klatt,indkpt,kpt_fullbz,&  
  &  mband,nkpt,nkpt_fullbz,nsppol,shiftk)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkpt_fullbz
  integer,intent(in) :: nsppol
  real(dp),intent(in) :: fermie
  real(dp),intent(in) :: eigen(mband,nkpt,nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indkpt(nkpt_fullbz)
  real(dp),intent(in) :: klatt(3,3)
  real(dp),intent(in) :: kpt_fullbz(3,nkpt_fullbz)
  real(dp),intent(in) :: shiftk(3)
 end subroutine prtxcfermsurf
end interface

interface
 subroutine recip_ylm (bessargmax,bess_fit,cgcband,iatsph,istwfk,kg,&  
  &  kpgnorm,nradint,nradintmax,mgfft,mlang,mpi_enreg,mpw,natom,natsph,ngfft,&  
  &  npw_k,ntypat,ph3d,prtsphere,rint,rmax,rprimd,sum_1ll_1atom,sum_1lm_1atom,typat,ucvol,ylm,znucl)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: istwfk
  integer,intent(in) :: mgfft
  integer,intent(in) :: mlang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: natsph
  integer,intent(in) :: npw_k
  integer,intent(in) :: nradintmax
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtsphere
  real(dp),intent(in) :: bessargmax
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: bess_fit(mpw,nradintmax,mlang)
  real(dp),intent(in) :: cgcband(2,npw_k)
  integer,intent(in) :: iatsph(natsph)
  integer,intent(in) :: kg(3,npw_k)
  real(dp),intent(in) :: kpgnorm(npw_k)
  integer,intent(in) :: nradint(natsph)
  real(dp),intent(in) :: ph3d(2,npw_k,natom)
  real(dp),intent(in) :: rint(nradintmax)
  real(dp),intent(in) :: rmax(natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: sum_1ll_1atom(mlang,natsph)
  real(dp),intent(out) :: sum_1lm_1atom(mlang*mlang,natsph)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: ylm(mpw,mlang*mlang)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine recip_ylm
end interface

interface
 subroutine simpson_int(nsimpson,simp_delta,simp_funct,simp_res)
  use defs_basis
  implicit none
  integer,intent(in) :: nsimpson
  real(dp),intent(in) :: simp_delta
  real(dp),intent(in) :: simp_funct(nsimpson)
  real(dp),intent(out) :: simp_res(nsimpson)
 end subroutine simpson_int
end interface

interface
 subroutine sphericaldens(fofg,gnorm,nfft,rmax,sphfofg)
  use defs_basis
  implicit none
  integer,intent(in) :: nfft
  real(dp),intent(in) :: rmax
  real(dp),intent(in) :: fofg(2,nfft)
  real(dp),intent(in) :: gnorm(nfft)
  real(dp),intent(out) :: sphfofg(2,nfft)
 end subroutine sphericaldens
end interface

interface
 subroutine tetrahedron (dos_fractions,dos_fractions_m,dos_fractions_paw1,dos_fractions_pawt1,&  
  &  dtset,fermie,eigen,fildata,mbesslang,m_dos_flag,ndosfraction,paw_dos_flag,rprimd)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: m_dos_flag
  integer,intent(in) :: mbesslang
  integer,intent(in) :: ndosfraction
  integer,intent(in) :: paw_dos_flag
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: fermie
  character(len=fnlen),intent(in) :: fildata
  real(dp),intent(in) :: dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction)
  real(dp),intent(in) :: dos_fractions_m(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*mbesslang*m_dos_flag)
  real(dp),intent(in) :: dos_fractions_paw1(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*paw_dos_flag)
  real(dp),intent(in) :: dos_fractions_pawt1(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*paw_dos_flag)
  real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine tetrahedron
end interface

end module interfaces_14occeig
!!***
