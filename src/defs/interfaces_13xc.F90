!!****m* ABINIT/interfaces_13xc
!! NAME
!! interfaces_13xc
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/13xc
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

module interfaces_13xc

 implicit none

interface
 subroutine calc_lifetime(ixcpositron,lifetime,nfft,nspden,positron,rhocore,rhore,rhototp,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: ixcpositron
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: positron
  real(dp),intent(out) :: lifetime
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: rhocore(nfft)
  real(dp),intent(in) :: rhore(nfft,nspden)
  real(dp),intent(in) :: rhototp(nfft)
 end subroutine calc_lifetime
end interface

interface
 subroutine calc_xc_ep(excapn,ixcpositron,positron,nfft,nspden,option,rhocore,rhor,rhore,rhorp,rhotote,rhototp,rsepts,rsppts,&  
  &  vhae,vhap,vpsp,vtrial,vxcapn)
  use defs_basis
  implicit none
  integer,intent(in) :: ixcpositron
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: positron
  real(dp),intent(out) :: excapn(nfft)
  real(dp),intent(in) :: rhocore(nfft)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rhore(nfft,nspden)
  real(dp),intent(in) :: rhorp(nfft,nspden)
  real(dp),intent(inout) :: rhotote(nfft)
  real(dp),intent(inout) :: rhototp(nfft)
  real(dp),intent(inout) :: rsepts(nfft)
  real(dp),intent(inout) :: rsppts(nfft)
  real(dp),intent(in) :: vhae(nfft,2)
  real(dp),intent(in) :: vhap(nfft,2)
  real(dp),intent(in) :: vpsp(nfft)
  real(dp),intent(inout) :: vtrial(nfft,nspden)
  real(dp),intent(inout) :: vxcapn(nfft)
 end subroutine calc_xc_ep
end interface

interface
 subroutine drivexc(exc,ixc,npts,nspden,order,rho_updn,vxc,ndvxc,ngr2,nvxcdgr,&  !Mandatory arguments
  &  dvxc,d2vxc,grho2_updn,vxcgr,exexch)    !Optional arguments
  use defs_basis
  implicit none
  integer,intent(in),optional :: exexch
  integer,intent(in) :: ixc
  integer,intent(in) :: ndvxc
  integer,intent(in) :: ngr2
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  integer,intent(in) :: nvxcdgr
  integer,intent(in) :: order
  real(dp),intent(out),optional :: d2vxc(npts)
  real(dp),intent(out),optional :: dvxc(npts,ndvxc)
  real(dp),intent(out) :: exc(npts)
  real(dp),intent(in),optional :: grho2_updn(npts,ngr2)
  real(dp),intent(in) :: rho_updn(npts,nspden)
  real(dp),intent(out) :: vxc(npts,nspden)
  real(dp),intent(out),optional :: vxcgr(npts,nvxcdgr)
 end subroutine drivexc
end interface

interface
 subroutine size_dvxc(ixc,ndvxc,ngr2,nspden,nvxcdgr,order)
  implicit none
  integer, intent(in) :: ixc
  integer, intent(out) :: ndvxc
  integer, intent(out) :: ngr2
  integer, intent(in) :: nspden
  integer, intent(out) :: nvxcdgr
  integer, intent(in) :: order
 end subroutine size_dvxc
end interface

interface
 subroutine hartre(cplex,gmet,gsqcut,izero,mpi_enreg,nfft,ngfft,paral_kgb,qphon,rhog,vhartr)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: izero
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(out) :: vhartr(cplex*nfft)
 end subroutine hartre
end interface

interface
 subroutine invcb(rhoarr,rspts,npts)
  use defs_basis
  implicit none
  integer,intent(in) :: npts
  real(dp),intent(in) :: rhoarr(npts)
  real(dp),intent(out) :: rspts(npts)
 end subroutine invcb
end interface

interface
 subroutine lifetime_bn(lifetime,npt,rhoer,rsepts,rhopr,rsppts,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  real(dp),intent(out) :: lifetime
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: rhoer(npt)
  real(dp),intent(in) :: rhopr(npt)
  real(dp),intent(in) :: rsepts(npt)
  real(dp),intent(in) :: rsppts(npt)
 end subroutine lifetime_bn
end interface

interface
 subroutine lifetime_psn(lifetime,npt,rhoer,rsepts,rhopr,rsppts,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  real(dp),intent(out) :: lifetime
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: rhoer(npt)
  real(dp),intent(in) :: rhopr(npt)
  real(dp),intent(in) :: rsepts(npt)
  real(dp),intent(in) :: rsppts(npt)
 end subroutine lifetime_psn
end interface

interface
 subroutine lifetime_rpa(lifetimeEq1,npt,rhoer,xccc3d,rsepts,rseptstot,rhopr,ucvol)
  use defs_basis
  implicit none
  integer :: npt
  real(dp) :: lifetimeEq1
  real(dp) :: ucvol
  real(dp) :: rhoer(npt)
  real(dp) :: rhopr(npt)
  real(dp) :: rsepts(npt)
  real(dp) :: rseptstot(npt)
  real(dp) :: xccc3d(npt)
 end subroutine lifetime_rpa
end interface

interface
 subroutine mkcore(corstr,dyfrx2,grxc,mpi_enreg,natom,nfft,nspden,ntypat,n1,n1xccc,&  
  &  n2,n3,option,rprimd,typat,ucvol,vxc,xcccrc,xccc1d,xccc3d,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: n1
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  type(mpi_type) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: corstr(6)
  real(dp),intent(out) :: dyfrx2(3,3,natom)
  real(dp),intent(out) :: grxc(3,natom)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vxc(nfft,nspden)
  real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat)
  real(dp),intent(inout) :: xccc3d(nfft)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mkcore
end interface

interface
 subroutine mkdenpos(iwarn,nfft,nspden,option,rhonow)
  use defs_basis
  implicit none
  integer,intent(inout) :: iwarn
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  real(dp),intent(inout) :: rhonow(nfft,nspden)
 end subroutine mkdenpos
end interface

interface
 subroutine mkvxc3(cplex,gmet,gsqcut,kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,n3xccc,option,&  
  &  paral_kgb,qphon,rhor1,rprimd,vxc1,xccc3d1)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhor1(cplex*nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vxc1(cplex*nfft,nspden)
  real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
 end subroutine mkvxc3
end interface

interface
 subroutine mkvxcgga3(cplex,gmet,gprimd,gsqcut,kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,&  
  &  paral_kgb,qphon,rhor1tmp,vxc1)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhor1tmp(cplex*nfft,2)
  real(dp),intent(out) :: vxc1(cplex*nfft,nspden)
 end subroutine mkvxcgga3
end interface

interface
 subroutine phase(ngfft,ph)
  use defs_basis
  implicit none
  integer,intent(in) :: ngfft
  real(dp),intent(out) :: ph(2*ngfft)
 end subroutine phase
end interface

interface
 subroutine rhohxc(dtset,enxc,gsqcut,izero,kxc,mpi_enreg,nfft,ngfft,&  
  &  nhat,nhatdim,nhatgr,nhatgrdim,nkxc,nspden,n3xccc,option,rhog,rhor,rprimd,&  
  &  strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d,k3xc)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: izero
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nhatdim
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: usexcnhat
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: enxc
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(out) :: vxcavg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out),optional :: k3xc(1:nfft)
  real(dp),intent(out) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: nhat(nfft,nspden*nhatdim)
  real(dp),intent(in) :: nhatgr(nfft,nspden,3*nhatgrdim)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: strsxc(6)
  real(dp),intent(out) :: vhartr(nfft)
  real(dp),intent(out) :: vxc(nfft,nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
 end subroutine rhohxc
end interface

interface
 subroutine xc_kernel(dtset,ixc,mpi_enreg,ngfft,nr,nsppol,rho,rprimd,igfft,npw,gmet,kernel,gvec,nq,qq)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ixc
  integer,intent(in) :: npw
  integer,intent(in) :: nq
  integer,intent(in) :: nr
  integer,intent(in) :: nsppol
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: gvec(3,npw)
  integer,intent(in) :: igfft(npw)
  complex,intent(out) :: kernel(npw,npw,nq)
  real(dp),intent(in) :: qq(3,nq)
  real(dp),intent(in) :: rho(nr,nsppol)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine xc_kernel
end interface

interface
 subroutine xcden (cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden,paral_kgb,qphon,rhor,rhonow)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: ishift
  integer,intent(in) :: nfft
  integer,intent(in) :: ngrad
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  type(mpi_type) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(out) :: rhonow(cplex*nfft,nspden,ngrad*ngrad)
  real(dp),intent(in) :: rhor(cplex*nfft,nspden)
 end subroutine xcden
end interface

interface
 subroutine xce_ap(exc,npt,rhoer,rsepts,rhopr,vxc)
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rhoer(npt)
  real(dp),intent(in) :: rhopr(npt)
  real(dp),intent(in) :: rsepts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xce_ap
end interface

interface
 subroutine xcepsn_tcdft(npt,rhoer,rsepts,rhopr,rsppts,vxc)
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  real(dp),intent(in) :: rhoer(npt)
  real(dp),intent(in) :: rhopr(npt)
  real(dp),intent(in) :: rsepts(npt)
  real(dp),intent(in) :: rsppts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xcepsn_tcdft
end interface

interface
 subroutine xchcth(dvxcdgr,exci,grho2_updn,ixc,npts,nspden,&  
  &  order,rho_updn,vxci)
  use defs_basis
  implicit none
  integer,intent(in) :: ixc
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  integer,intent(in) :: order
  real(dp),intent(out) :: dvxcdgr(npts,2)
  real(dp),intent(out) :: exci(npts)
  real(dp),intent(in) :: grho2_updn(npts,2*nspden-1)
  real(dp),intent(in) :: rho_updn(npts,nspden)
  real(dp),intent(out) :: vxci(npts,nspden)
 end subroutine xchcth
end interface

interface
 subroutine xchelu(exc,npt,order,rspts,vxc,dvxc)  ! dvxc is optional
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xchelu
end interface

interface
 subroutine xclb(grho2_updn,npts,nspden,rho_updn,vxci)
  use defs_basis
  implicit none
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  real(dp),intent(in) :: grho2_updn(npts,2*nspden-1)
  real(dp),intent(in) :: rho_updn(npts,nspden)
  real(dp),intent(inout) :: vxci(npts,nspden)
 end subroutine xclb
end interface

interface
 subroutine xcmult (dnexcdn,nfft,ngrad,nspden,nspgrad,rhonow)
  use defs_basis
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: ngrad
  integer,intent(in) :: nspden
  integer,intent(in) :: nspgrad
  real(dp),intent(in) :: dnexcdn(nfft,nspgrad)
  real(dp),intent(inout) :: rhonow(nfft,nspden,ngrad*ngrad)
 end subroutine xcmult
end interface

interface
 subroutine xcp_ap(exc,npt,rhoer,rsepts,vxc)
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rhoer(npt)
  real(dp),intent(in) :: rsepts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xcp_ap
end interface

interface
 subroutine xcpbe(exci,npts,nspden,option,order,rho_updn,vxci,ndvxci,ngr2,&  !Mandatory Arguments
  &  d2vxci,dvxcdgr,dvxci,exexch,grho2_updn)                          !Optional Arguments
  use defs_basis
  implicit none
  integer,intent(in),optional :: exexch
  integer,intent(in) :: ndvxci
  integer,intent(in) :: ngr2
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: order
  real(dp),intent(out),optional :: d2vxci(npts)
  real(dp),intent(out),optional :: dvxcdgr(npts,3)
  real(dp),intent(out),optional :: dvxci(npts,ndvxci)
  real(dp),intent(out) :: exci(npts)
  real(dp),intent(in),optional :: grho2_updn(npts,ngr2)
  real(dp),intent(in) :: rho_updn(npts,nspden)
  real(dp),intent(out) :: vxci(npts,nspden)
 end subroutine xcpbe
end interface

interface
 subroutine xcpot (cplex,dnexcdn,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden,&  
  &  nspgrad,paral_kgb,qphon,rhonow,vxc)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: ishift
  integer,intent(in) :: nfft
  integer,intent(in) :: ngrad
  integer,intent(in) :: nspden
  integer,intent(in) :: nspgrad
  integer,intent(in) :: paral_kgb
  type(mpi_type) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: dnexcdn(cplex*nfft,nspgrad)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhonow(cplex*nfft,nspden,ngrad*ngrad)
  real(dp),intent(out) :: vxc(cplex*nfft,nspden)
 end subroutine xcpot
end interface

interface
 subroutine xcppsn_tcdft(npt,rhoer,rsepts,rhopr,rsppts,vxc)
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  real(dp),intent(in) :: rhoer(npt)
  real(dp),intent(in) :: rhopr(npt)
  real(dp),intent(in) :: rsepts(npt)
  real(dp),intent(in) :: rsppts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xcppsn_tcdft
end interface

interface
 subroutine xcpzca(exc,npt,order,rhor,rspts,vxc,&  !Mandatory arguments
  &  dvxc)                            !Optional arguments
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rhor(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xcpzca
end interface

interface
 subroutine xcspol(exc,npts,nspden,order,rspts,vxc,zeta,ndvxc,&  !Mandatory arguments
  &  dvxc)                            !Optional arguments
  use defs_basis
  implicit none
  integer,intent(in) :: ndvxc
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npts,ndvxc)
  real(dp),intent(out) :: exc(npts)
  real(dp),intent(in) :: rspts(npts)
  real(dp),intent(out) :: vxc(npts,nspden)
  real(dp),intent(in) :: zeta(npts)
 end subroutine xcspol
end interface

interface
 subroutine xctetr(exc,npt,order,rhor,rspts,vxc,&  !Mandatory arguments
  &  d2vxc,dvxc)                    !Optional arguments
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: d2vxc(npt)
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rhor(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xctetr
end interface

interface
 subroutine xcwign(exc,npt,order,rhor,rspts,vxc,&  !Mandatory arguments
  &  dvxc)                           !Optional arguments
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rhor(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xcwign
end interface

interface
 subroutine xcxalp(exc,npt,order,rspts,vxc, dvxc)  ! dvxc is optional
  use defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xcxalp
end interface

end module interfaces_13xc
!!***
