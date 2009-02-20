!!****m* ABINIT/interfaces_13psp
!! NAME
!! interfaces_13psp
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/13psp
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

module interfaces_13psp

 implicit none

interface
 subroutine cc_derivatives(rad,ff,ff1,ff2,mmax,n1xccc,rchrg,xccc1d)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: n1xccc
  real(dp),intent(in) :: rchrg
  real(dp),intent(in) :: ff(mmax)
  real(dp),intent(in) :: ff1(mmax)
  real(dp),intent(in) :: ff2(mmax)
  real(dp),intent(in) :: rad(mmax)
  real(dp),intent(out) :: xccc1d(n1xccc,6)
 end subroutine cc_derivatives
end interface

interface
 subroutine der_int(ff,df,rr,dr,nlast,smf)
  use defs_basis
  implicit none
  integer,parameter :: nmax=2000
  integer,intent(in) :: nlast
  real(dp),intent(out) :: smf
  real(dp), intent(out) :: df(0:nmax)
  real(dp), intent(in) :: dr(0:nmax)
  real(dp), intent(in) :: ff(0:nmax)
  real(dp), intent(in) :: rr(0:nmax)
 end subroutine der_int
end interface

interface
 subroutine gg1cc(gg1cc_xx,xx)
  use defs_basis
  implicit none
  real(dp),intent(out) :: gg1cc_xx
  real(dp),intent(in) :: xx
 end subroutine gg1cc
end interface

interface
 subroutine gp1cc(gp1cc_xx,xx)
  use defs_basis
  implicit none
  real(dp),intent(out) :: gp1cc_xx
  real(dp),intent(in) :: xx
 end subroutine gp1cc
end interface

interface
 subroutine gpp1cc(gpp1cc_xx,xx)
  use defs_basis
  implicit none
  real(dp),intent(out) :: gpp1cc_xx
  real(dp),intent(in) :: xx
 end subroutine gpp1cc
end interface

interface
 subroutine psp10in(dtset, ekb, epsatm, ffspl, indlmn, ipsp, lmax,&  
  &  nproj, psps, pspso, vlspl, zion)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ipsp
  integer,intent(inout) :: lmax
  integer,intent(in) :: pspso
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: epsatm
  type(pseudopotential_type),intent(inout) :: psps
  real(dp),intent(in) :: zion
  real(dp),intent(out) :: ekb(psps%lnmax)
  real(dp),intent(out) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
  integer,intent(out) :: indlmn(6,psps%lmnmax)
  integer,intent(out) :: nproj(psps%mpssoang)
  real(dp),intent(out) :: vlspl(psps%mqgrid_ff,2)
 end subroutine psp10in
end interface

interface
 subroutine psp10nl(ekb,ffspl,hij,lmax,mproj,mpsang,mqgrid,nproj,qgrid,rr)
  use defs_basis
  implicit none
  integer,intent(in) :: lmax
  integer,intent(in) :: mproj
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  real(dp),intent(out) :: ekb(mpsang,mproj)
  real(dp),intent(out) :: ffspl(mqgrid,2,mpsang,mproj)
  real(dp),intent(in) :: hij(0:lmax,3,3)
  integer,intent(in) :: nproj(mpsang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rr(0:lmax)
 end subroutine psp10nl
end interface

interface
 subroutine psp1cc(fchrg,n1xccc,xccc1d)
  use defs_basis
  implicit none
  integer,intent(in) :: n1xccc
  real(dp),intent(in) :: fchrg
  real(dp),intent(out) :: xccc1d(n1xccc,6)
 end subroutine psp1cc
end interface

interface
 subroutine psp1in(dq,ekb,ekb1,ekb2,epsatm,epspsp,&  
  &  e990,e999,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&  
  &  mmax,mpsang,mqgrid,nproj,n1xccc,pspcod,&  
  &  qchrg,qgrid,rcpsp,rms,useylm,vlspl,xcccrc,xccc1d,&  
  &  zion,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: lloc
  integer,intent(in) :: lmax
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: pspcod
  integer,intent(in) :: useylm
  real(dp),intent(in) :: dq
  real(dp),intent(out) :: epsatm
  real(dp),intent(out) :: qchrg
  real(dp),intent(out) :: xcccrc
  real(dp),intent(in) :: zion
  real(dp),intent(in) :: znucl
  real(dp),intent(out) :: e990(mpsang)
  real(dp),intent(out) :: e999(mpsang)
  real(dp),intent(out) :: ekb(lnmax)
  real(dp),intent(out) :: ekb1(mpsang)
  real(dp),intent(out) :: ekb2(mpsang)
  real(dp),intent(out) :: epspsp(mpsang)
  real(dp),intent(out) :: ffspl(mqgrid,2,lnmax)
  integer,intent(out) :: indlmn(6,lmnmax)
  integer,intent(out) :: nproj(mpsang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(out) :: rcpsp(mpsang)
  real(dp),intent(out) :: rms(mpsang)
  real(dp),intent(out) :: vlspl(mqgrid,2)
  real(dp),intent(out) :: xccc1d(n1xccc,6)
 end subroutine psp1in
end interface

interface
 subroutine psp1lo(drad,epsatm,mmax,mqgrid,qgrid,q2vq,rad,&  
  &  vloc,wksincos,yp1,ypn,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: mqgrid
  real(dp),intent(out) :: epsatm
  real(dp),intent(out) :: yp1
  real(dp),intent(out) :: ypn
  real(dp),intent(in) :: zion
  real(dp),intent(in) :: drad(mmax)
  real(dp),intent(out) :: q2vq(mqgrid)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rad(mmax)
  real(dp),intent(in) :: vloc(mmax)
  real(dp),intent(inout) :: wksincos(mmax,2,2)
 end subroutine psp1lo
end interface

interface
 subroutine psp1nl(dr,ekb,ffspl,lloc,lmax,mmax,mpsang,mqgrid,&  
  &  qgrid,rad,vloc,vpspll,wfll,wksincos)
  use defs_basis
  implicit none
  integer,intent(in) :: lloc
  integer,intent(in) :: lmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: dr(mmax)
  real(dp),intent(out) :: ekb(mpsang)
  real(dp),intent(out) :: ffspl(mqgrid,2,mpsang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rad(mmax)
  real(dp),intent(in) :: vloc(mmax)
  real(dp),intent(in) :: vpspll(mmax,mpsang)
  real(dp),intent(in) :: wfll(mmax,mpsang)
  real(dp),intent(inout) :: wksincos(mmax,2,2)
 end subroutine psp1nl
end interface

interface
 subroutine psp2in(dtset,ekb,epsatm,ffspl,indlmn,ipsp,lmax,nproj,psps,vlspl,dvlspl,zion)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ipsp
  integer,intent(in) :: lmax
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: epsatm
  type(pseudopotential_type),intent(inout) :: psps
  real(dp),intent(in) :: zion
  real(dp),intent(out) :: dvlspl(psps%mqgrid_vl,2)
  real(dp),intent(out) :: ekb(psps%lnmax)
  real(dp),intent(out) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
  integer,intent(out) :: indlmn(6,psps%lmnmax)
  integer,intent(out) :: nproj(psps%mpsang)
  real(dp),intent(out) :: vlspl(psps%mqgrid_vl,2)
 end subroutine psp2in
end interface

interface
 subroutine psp2lo(cc1,cc2,cc3,cc4,dvloc,epsatm,mqgrid,qgrid,q2vq,&  
  &  rloc,vlspl_recipSpace,yp1,ypn,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: cc1
  real(dp),intent(in) :: cc2
  real(dp),intent(in) :: cc3
  real(dp),intent(in) :: cc4
  real(dp),intent(out) :: epsatm
  real(dp),intent(in) :: rloc
  logical,intent(in) :: vlspl_recipSpace
  real(dp),intent(out) :: yp1
  real(dp),intent(out) :: ypn
  real(dp),intent(in) :: zion
  real(dp),intent(out) :: dvloc(mqgrid)
  real(dp),intent(out) :: q2vq(mqgrid)
  real(dp),intent(in) :: qgrid(mqgrid)
 end subroutine psp2lo
end interface

interface
 subroutine psp2nl(ekb,ffspl,h1p,h1s,h2s,lnmax,mqgrid,qgrid,rrp,rrs)
  use defs_basis
  implicit none
  integer,intent(in) :: lnmax
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: h1p
  real(dp),intent(in) :: h1s
  real(dp),intent(in) :: h2s
  real(dp),intent(in) :: rrp
  real(dp),intent(in) :: rrs
  real(dp),intent(out) :: ekb(lnmax)
  real(dp),intent(out) :: ffspl(mqgrid,2,lnmax)
  real(dp),intent(in) :: qgrid(mqgrid)
 end subroutine psp2nl
end interface

interface
 subroutine psp2params_init(gth_params, npsp)
  use defs_datatypes
  implicit none
  integer,intent(in) :: npsp
  type(pseudopotential_gth_type),intent(out) :: gth_params
 end subroutine psp2params_init
end interface

interface
 subroutine psp2params_free(gth_params)
  use defs_datatypes
  implicit none
  type(pseudopotential_gth_type),intent(inout) :: gth_params
 end subroutine psp2params_free
end interface

interface
 subroutine psp3in(dtset, ekb, epsatm, ffspl, indlmn, ipsp, lmax,&  
  &  nproj, psps, pspso, vlspl, zion)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ipsp
  integer,intent(inout) :: lmax
  integer,intent(in) :: pspso
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: epsatm
  type(pseudopotential_type),intent(inout) :: psps
  real(dp),intent(in) :: zion
  real(dp),intent(out) :: ekb(psps%lnmax)
  real(dp),intent(out) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
  integer,intent(out) :: indlmn(6,psps%lmnmax)
  integer,intent(out) :: nproj(psps%mpssoang)
  real(dp),intent(out) :: vlspl(psps%mqgrid_ff,2)
 end subroutine psp3in
end interface

interface
 subroutine psp3nl(ekb,ffspl,h11s,h22s,h33s,h11p,h22p,h33p,h11d,h22d,&  
  &  h33d,h11f,mproj,mpsang,mqgrid,qgrid,rrd,rrf,rrp,rrs)
  use defs_basis
  implicit none
  integer,intent(in) :: mproj
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: h11d
  real(dp),intent(in) :: h11f
  real(dp),intent(in) :: h11p
  real(dp),intent(in) :: h11s
  real(dp),intent(in) :: h22d
  real(dp),intent(in) :: h22p
  real(dp),intent(in) :: h22s
  real(dp),intent(in) :: h33d
  real(dp),intent(in) :: h33p
  real(dp),intent(in) :: h33s
  real(dp),intent(in) :: rrd
  real(dp),intent(in) :: rrf
  real(dp),intent(in) :: rrp
  real(dp),intent(in) :: rrs
  real(dp),intent(out) :: ekb(mpsang,mproj)
  real(dp),intent(out) :: ffspl(mqgrid,2,mpsang,mproj)
  real(dp),intent(in) :: qgrid(mqgrid)
 end subroutine psp3nl
end interface

interface
 subroutine psp4cc(fchrg,n1xccc,xccc1d)
  use defs_basis
  implicit none
  integer,intent(in) :: n1xccc
  real(dp),intent(in) :: fchrg
  real(dp),intent(out) :: xccc1d(n1xccc,6)
 end subroutine psp4cc
end interface

interface
 subroutine psp5in(ekb,ekb1,ekb2,epsatm,epspsp,e990,e999,ffspl,indlmn,&  
  &  lloc,lmax,lmnmax,lnmax,mmax,mpsang,mpssoang,mqgrid,&  
  &  nproj,n1xccc,pspcod,pspso,qchrg,qgrid,rcpsp,rms,&  
  &  useylm,vlspl,xcccrc,xccc1d,zion,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: lloc
  integer,intent(in) :: lmax
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: pspcod
  integer,intent(in) :: pspso
  integer,intent(in) :: useylm
  real(dp),intent(out) :: epsatm
  real(dp),intent(out) :: qchrg
  real(dp),intent(out) :: xcccrc
  real(dp),intent(in) :: zion
  real(dp),intent(in) :: znucl
  real(dp),intent(out) :: e990(mpssoang)
  real(dp),intent(out) :: e999(mpssoang)
  real(dp),intent(out) :: ekb(lnmax)
  real(dp),intent(out) :: ekb1(mpssoang)
  real(dp),intent(out) :: ekb2(mpssoang)
  real(dp),intent(out) :: epspsp(mpssoang)
  real(dp),intent(out) :: ffspl(mqgrid,2,lnmax)
  integer,intent(out) :: indlmn(6,lmnmax)
  integer,intent(out) :: nproj(mpssoang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(out) :: rcpsp(mpssoang)
  real(dp),intent(out) :: rms(mpssoang)
  real(dp),intent(out) :: vlspl(mqgrid,2)
  real(dp),intent(out) :: xccc1d(n1xccc,6)
 end subroutine psp5in
end interface

interface
 subroutine psp5lo(al,epsatm,mmax,mqgrid,qgrid,q2vq,rad,&  
  &  vloc,yp1,ypn,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: al
  real(dp),intent(out) :: epsatm
  real(dp),intent(out) :: yp1
  real(dp),intent(out) :: ypn
  real(dp),intent(in) :: zion
  real(dp),intent(out) :: q2vq(mqgrid)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rad(mmax)
  real(dp),intent(in) :: vloc(mmax)
 end subroutine psp5lo
end interface

interface
 subroutine psp5nl(al,ekb,ffspl,lmax,mmax,mpsang,mqgrid,&  
  &  qgrid,rad,vloc,vpspll,wfll)
  use defs_basis
  implicit none
  integer,intent(in) :: lmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: al
  real(dp),intent(out) :: ekb(mpsang)
  real(dp),intent(out) :: ffspl(mqgrid,2,mpsang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rad(mmax)
  real(dp),intent(in) :: vloc(mmax)
  real(dp),intent(in) :: vpspll(mmax,mpsang)
  real(dp),intent(in) :: wfll(mmax,mpsang)
 end subroutine psp5nl
end interface

interface
 subroutine psp6cc(mmax,n1xccc,rchrg,xccc1d)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: n1xccc
  real(dp),intent(in) :: rchrg
  real(dp),intent(out) :: xccc1d(n1xccc,6)
 end subroutine psp6cc
end interface

interface
 subroutine psp6cc_drh(mmax,n1xccc,rchrg,xccc1d)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: n1xccc
  real(dp),intent(in) :: rchrg
  real(dp),intent(out) :: xccc1d(n1xccc,6)
 end subroutine psp6cc_drh
end interface

interface
 subroutine psp6ccpos(mmax,n1xccc,rchrg,xccc1d,vhtnzc,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: n1xccc
  real(dp),intent(in) :: rchrg
  real(dp),intent(in) :: znucl
  real(dp),intent(out) :: vhtnzc(mmax)
  real(dp),intent(out) :: xccc1d(n1xccc,6)
 end subroutine psp6ccpos
end interface

interface
 subroutine psp6in(ekb,epsatm,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&  
  &  mmax,mpsang,mqgrid,nproj,n1xccc,optnlxccc,positron,qchrg,qgrid,&  
  &  useylm,vlspl,xcccrc,xccc1d,zion,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: lloc
  integer,intent(in) :: lmax
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: optnlxccc
  integer,intent(in) :: positron
  integer,intent(in) :: useylm
  real(dp),intent(out) :: epsatm
  real(dp),intent(out) :: qchrg
  real(dp),intent(out) :: xcccrc
  real(dp),intent(in) :: zion
  real(dp),intent(in) :: znucl
  real(dp),intent(out) :: ekb(lnmax)
  real(dp),intent(out) :: ffspl(mqgrid,2,lnmax)
  integer,intent(out) :: indlmn(6,lmnmax)
  integer,intent(out) :: nproj(mpsang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(out) :: vlspl(mqgrid,2)
  real(dp),intent(out) :: xccc1d(n1xccc,6)
 end subroutine psp6in
end interface

interface
 subroutine psp7cc(core_mesh,n1xccc,rchrg,tncore,xccc1d)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: n1xccc
  type(pawrad_type),intent(in) :: core_mesh
  real(dp),intent(in) :: rchrg
  real(dp),intent(in) :: tncore(core_mesh%mesh_size)
  real(dp),intent(out) :: xccc1d(n1xccc,6)
 end subroutine psp7cc
end interface

interface
 subroutine psp7cg(dnqdq0,mqgrid,qgrid,nq,radmesh,nr,yp1,ypn)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mqgrid
  real(dp),intent(out) :: dnqdq0
  type(pawrad_type),intent(in) :: radmesh
  real(dp),intent(out) :: yp1
  real(dp),intent(out) :: ypn
  real(dp),intent(out) :: nq(mqgrid)
  real(dp),intent(in) :: nr(radmesh%mesh_size)
  real(dp),intent(in) :: qgrid(mqgrid)
 end subroutine psp7cg
end interface

interface
 subroutine psp7in(epsatm,ffspl,indlmn,ixc,lmax,lmnmax,lnmax,mmax,&  
  &  mqgrid_ff,mqgrid_vl,pawrad,pawtab,pawxcdev,&  
  &  pspso,qgrid_ff,qgrid_vl,vlspl,xcccrc,xclevel,zion)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ixc
  integer,intent(in) :: lmax
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mqgrid_ff
  integer,intent(in) :: mqgrid_vl
  integer,intent(in) :: pawxcdev
  integer,intent(in) :: pspso
  integer,intent(in) :: xclevel
  real(dp),intent(out) :: epsatm
  type(pawrad_type),intent(out) :: pawrad
  type(pawtab_type),intent(out) :: pawtab
  real(dp),intent(out) :: xcccrc
  real(dp),intent(in) :: zion
  real(dp),intent(out) :: ffspl(mqgrid_ff,2,lnmax)
  integer,intent(out) :: indlmn(6,lmnmax)
  real(dp),intent(in) :: qgrid_ff(mqgrid_ff)
  real(dp),intent(in) :: qgrid_vl(mqgrid_vl)
  real(dp),intent(out) :: vlspl(mqgrid_vl,2)
 end subroutine psp7in
end interface

interface
 subroutine bound_deriv(func,mesh,nn,yp1,ypn)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, intent(in) :: nn
  type(pawrad_type),intent(in) :: mesh
  real(dp), intent(out) :: yp1
  real(dp), intent(out) :: ypn
  real(dp), intent(in) :: func(nn)
 end subroutine bound_deriv
end interface

interface
 subroutine psp7lo(epsatm,mqgrid,qgrid,q2vq,radmesh,vloc,yp1,ypn,zion)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mqgrid
  real(dp),intent(out) :: epsatm
  type(pawrad_type),intent(in) :: radmesh
  real(dp),intent(out) :: yp1
  real(dp),intent(out) :: ypn
  real(dp),intent(in) :: zion
  real(dp),intent(out) :: q2vq(mqgrid)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: vloc(radmesh%mesh_size)
 end subroutine psp7lo
end interface

interface
 subroutine psp7nl(ffspl,indlmn,lmnmax,lnmax,mqgrid,qgrid,radmesh,wfll)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mqgrid
  type(pawrad_type),intent(in) :: radmesh
  real(dp),intent(out) :: ffspl(mqgrid,2,lnmax)
  integer,intent(in) :: indlmn(6,lmnmax)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: wfll(radmesh%mesh_size,lnmax)
 end subroutine psp7nl
end interface

interface
 subroutine psp8cc(mmax,n1xccc,rchrg,xccc1d)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: n1xccc
  real(dp),intent(in) :: rchrg
  real(dp),intent(out) :: xccc1d(n1xccc,6)
 end subroutine psp8cc
end interface

interface
 subroutine psp8in(ekb,epsatm,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&  
  &  mmax,mpsang,mqgrid,nproj,n1xccc,qchrg,qgrid,&  
  &  useylm,vlspl,xcccrc,xccc1d,zion,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: lloc
  integer,intent(in) :: lmax
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: useylm
  real(dp),intent(out) :: epsatm
  real(dp),intent(out) :: qchrg
  real(dp),intent(out) :: xcccrc
  real(dp),intent(in) :: zion
  real(dp),intent(in) :: znucl
  real(dp),intent(out) :: ekb(lnmax)
  real(dp),intent(out) :: ffspl(mqgrid,2,lnmax)
  integer,intent(out) :: indlmn(6,lmnmax)
  integer,intent(out) :: nproj(mpsang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(out) :: vlspl(mqgrid,2)
  real(dp),intent(out) :: xccc1d(n1xccc,6)
 end subroutine psp8in
end interface

interface
 subroutine psp8lo(amesh,epsatm,mmax,mqgrid,qgrid,q2vq,rad,vloc,yp1,ypn,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: amesh
  real(dp),intent(out) :: epsatm
  real(dp),intent(out) :: yp1
  real(dp),intent(out) :: ypn
  real(dp),intent(in) :: zion
  real(dp),intent(out) :: q2vq(mqgrid)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rad(mmax)
  real(dp),intent(in) :: vloc(mmax)
 end subroutine psp8lo
end interface

interface
 subroutine psp8nl(amesh,ffspl,lmax,lnmax,mmax,mpsang,mqgrid,nproj,qgrid,rad,vpspll)
  use defs_basis
  implicit none
  integer,intent(in) :: lmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: amesh
  real(dp),intent(out) :: ffspl(mqgrid,2,lnmax)
  integer,intent(in) :: nproj(mpsang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rad(mmax)
  real(dp),intent(in) :: vpspll(mmax,lnmax)
 end subroutine psp8nl
end interface

interface
 subroutine psp9in(filpsp,ekb,epsatm,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&  
  &  mmax,mpsang,mqgrid,nproj,n1xccc,optnlxccc,qchrg,qgrid,&  
  &  useylm,vlspl,xcccrc,xccc1d,zion,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: lloc
  integer,intent(in) :: lmax
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(out) :: mmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: optnlxccc
  integer,intent(in) :: useylm
  real(dp),intent(out) :: epsatm
  character(len=fnlen),intent(in) :: filpsp
  real(dp),intent(out) :: qchrg
  real(dp),intent(out) :: xcccrc
  real(dp),intent(in) :: zion
  real(dp),intent(in) :: znucl
  real(dp),intent(out) :: ekb(lnmax)
  real(dp),intent(out) :: ffspl(mqgrid,2,lnmax)
  integer,intent(out) :: indlmn(6,lmnmax)
  integer,intent(out) :: nproj(mpsang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(out) :: vlspl(mqgrid,2)
  real(dp),intent(out) :: xccc1d(n1xccc,6)
 end subroutine psp9in
end interface

interface
 subroutine pspatm(dq,dtset,ekb,epsatm,ffspl,indlmn,ipsp,pawrad,pawtab,&  
  &  psps,vlspl,dvlspl,xcccrc,xccc1d)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ipsp
  real(dp),intent(in) :: dq
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: epsatm
  type(pawrad_type),intent(out) :: pawrad
  type(pawtab_type),intent(out) :: pawtab
  type(pseudopotential_type),intent(inout) :: psps
  real(dp),intent(out) :: xcccrc
  real(dp),intent(out) :: dvlspl(psps%mqgrid_vl,2)
  real(dp),intent(out) :: ekb(psps%dimekb*(1-psps%usepaw))
  real(dp),intent(out) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
  integer,intent(out) :: indlmn(6,psps%lmnmax)
  real(dp),intent(out) :: vlspl(psps%mqgrid_vl,2)
  real(dp),intent(out) :: xccc1d(psps%n1xccc*(1-psps%usepaw),6)
 end subroutine pspatm
end interface

interface
 subroutine pspcor(ecore,epsatm,natom,ntypat,typat,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp),intent(out) :: ecore
  real(dp),intent(in) :: epsatm(ntypat)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine pspcor
end interface

interface
 subroutine pspini(dtset,ecore,gencond,gsqcut,gsqcutdg,level,&  
  &  pawrad,pawtab,psps,rprimd)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: gencond
  integer,intent(in) :: level
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: ecore
  real(dp),intent(in) :: gsqcut
  real(dp),intent(in) :: gsqcutdg
  type(pseudopotential_type), intent(inout) :: psps
  type(pawrad_type), intent(out) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type), intent(out) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine pspini
end interface

interface
 subroutine psxml2ab( psxml, znucl, zion, pspcod, pspxc, lmax, iwrite )
  use defs_basis
  use m_pseudo_types
  implicit none
  integer,intent(in) :: iwrite
  integer,intent(out) :: lmax
  integer,intent(out) :: pspcod
  integer,intent(out) :: pspxc
  type(pseudo_t),intent(in) :: psxml
  real(dp),intent(out) :: zion
  real(dp),intent(out) :: znucl
 end subroutine psxml2ab
end interface

interface
 subroutine radii_ps( vps, rofi, zval, nrval, lmxkb, nrgauss, rgauss, rgauss2)
  use defs_basis
  implicit none
  integer,intent(in) :: lmxkb
  integer,intent(out) :: nrgauss
  integer,intent(in) :: nrval
  real(dp),intent(out) :: rgauss
  real(dp),intent(out) :: rgauss2
  real(dp),intent(in) :: zval
  real(dp),intent(in) :: rofi(:)
  real(dp),intent(in) :: vps(:,0:)
 end subroutine radii_ps
end interface

interface
 subroutine sbf8(nm,xx,sb_out)
  use defs_basis
  implicit none
  integer,intent(in) :: nm
  real(dp),intent(in) :: xx
  real(dp),intent(out) :: sb_out(nm)
 end subroutine sbf8
end interface

interface
 subroutine shapebes(al,ql,ll,rc)
  use defs_basis
  implicit none
  integer :: ll
  real(dp) :: rc
  real(dp) :: al(2)
  real(dp) :: ql(2)
 end subroutine shapebes
end interface

interface
 subroutine sincos(iq,irmax,mmax,pspwk,rad,tpiq)
  use defs_basis
  implicit none
  integer,intent(in) :: iq
  integer,intent(in) :: irmax
  integer,intent(in) :: mmax
  real(dp),intent(in) :: tpiq
  real(dp),intent(inout) :: pspwk(mmax,2,2)
  real(dp),intent(in) :: rad(mmax)
 end subroutine sincos
end interface

interface
 subroutine smoothvlocal(  lmax, npts, scale, step, vlocal, vps, zval)
  use defs_basis
  implicit none
  integer, intent(in) :: lmax
  integer, intent(in) :: npts
  real(dp), intent(in) :: scale
  real(dp), intent(in) :: step
  real(dp), intent(in) :: zval
  real(dp), intent(out) :: vlocal(npts)
  real(dp), intent(in) :: vps(npts,0:lmax)
 end subroutine smoothvlocal
end interface

interface
 subroutine vlocal2( zval, nrval, a, rofi, drdi, s, vps, nrgauss,&  
  &  vlocal,nchloc,chlocal )
  use defs_basis
  implicit none
  integer,  intent(out) :: nchloc
  integer,  intent(inout) :: nrgauss
  integer,  intent(in) :: nrval
  real(dp), intent(in) :: a
  real(dp), intent(in) :: zval
  real(dp), intent(out) :: chlocal(:)
  real(dp), intent(in) :: drdi(:)
  real(dp), intent(in) :: rofi(:)
  real(dp), intent(in) :: s(:)
  real(dp), intent(out) :: vlocal(:)
  real(dp), intent(in) :: vps(:)
 end subroutine vlocal2
end interface

interface
 subroutine vlocal1( zval, nrval, a, rofi, drdi, s, rgauss, vlocal,&  
  &  nchloc, chlocal)
  use defs_basis
  implicit none
  integer,  intent(out) :: nchloc
  integer,  intent(in) :: nrval
  real(dp), intent(in) :: a
  real(dp), intent(inout) :: rgauss
  real(dp), intent(in) :: zval
  real(dp), intent(out) :: chlocal(:)
  real(dp), intent(in) :: drdi(:)
  real(dp), intent(in) :: rofi(:)
  real(dp), intent(in) :: s(:)
  real(dp), intent(out) :: vlocal(:)
 end subroutine vlocal1
end interface

interface
 subroutine vhrtre(R2RHO,V,R,DRDI,SRDRDI,NR,A)
  use defs_basis
  implicit none
  integer :: NR
  real(dp) :: A
  real(dp) :: DRDI(NR)
  real(dp) :: R(NR)
  real(dp) :: R2RHO(NR)
  real(dp) :: SRDRDI(NR)
  real(dp) :: V(NR)
 end subroutine vhrtre
end interface

interface
 subroutine solvbes(root,alpha,beta,ll,nq)
  use defs_basis
  implicit none
  integer :: ll
  integer :: nq
  real(dp) :: alpha
  real(dp) :: beta
  real(dp) :: root(nq)
 end subroutine solvbes
end interface

interface
 function vander(a,x) result(f)
  use defs_basis
  implicit none
  real(dp),intent(in) :: a
  real(dp) :: f
  real(dp),intent(in) :: x
 end function vander
end interface

end module interfaces_13psp
!!***
