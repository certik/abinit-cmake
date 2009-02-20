!!****m* ABINIT/interfaces_11util
!! NAME
!! interfaces_11util
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/11util
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

module interfaces_11util

 implicit none

interface
 subroutine acrossb(a,b,c)
  use defs_basis
  implicit none
  real(dp),intent(in) :: a(3)
  real(dp),intent(in) :: b(3)
  real(dp),intent(out) :: c(3)
 end subroutine acrossb
end interface

interface
 subroutine appdig(integ,string,strinn)
  implicit none
  integer,intent(in) :: integ
  character(len=*),intent(in) :: string
  character(len=*),intent(out) :: strinn
 end subroutine appdig
end interface

interface
 subroutine atmdata(amu,rcov,symbol,znucl)
  use defs_basis
  implicit none
  real(dp),intent(out) :: amu
  real(dp),intent(out) :: rcov
  character(len=2),intent(out) :: symbol
  real(dp),intent(in) :: znucl
 end subroutine atmdata
end interface

interface
 subroutine atmlength(densty,length,zion,znucl)
  use defs_basis
  implicit none
  real(dp),intent(in) :: densty
  real(dp),intent(out) :: length
  real(dp),intent(in) :: zion
  real(dp),intent(in) :: znucl
 end subroutine atmlength
end interface

interface
 subroutine besjm(arg,besjx,cosx,nn,nx,sinx,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: nn
  integer,intent(in) :: nx
  real(dp),intent(in) :: arg
  real(dp),intent(out) :: besjx(nx)
  real(dp),intent(in) :: cosx(nx)
  real(dp),intent(in) :: sinx(nx)
  real(dp),intent(in) :: xx(nx)
 end subroutine besjm
end interface

interface
 subroutine calc_psden(ff,mesh,nc,rc,step)
  use defs_basis
  implicit none
  integer,intent(in) :: mesh
  real(dp),intent(in) :: rc
  real(dp),intent(in) :: step
  real(dp),intent(out) :: ff(mesh)
  real(dp),intent(in) :: nc(mesh)
 end subroutine calc_psden
end interface

interface
 subroutine calc_psden_log(ff,mmax,nc,rc,rad)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  real(dp),intent(in) :: rc
  real(dp),intent(out) :: ff(mmax)
  real(dp),intent(in) :: nc(mmax)
  real(dp),intent(in) :: rad(mmax)
 end subroutine calc_psden_log
end interface

interface
 subroutine calc_vhtnzc(nc,rc,vhtnzc,mesh,rad,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: mesh
  real(dp),intent(inout) :: rc
  real(dp),intent(in) :: znucl
  real(dp),intent(in) :: nc(mesh)
  real(dp),intent(in) :: rad(mesh)
  real(dp),intent(out) :: vhtnzc(mesh)
 end subroutine calc_vhtnzc
end interface

interface
 subroutine canon9(num,red,shift)
  use defs_basis
  implicit none
  real(dp),intent(in) :: num
  real(dp),intent(out) :: red
  real(dp),intent(out) :: shift
 end subroutine canon9
end interface

interface
 subroutine chknm8(nmxpct,nmfond)
  implicit none
  character(len=9),intent(in) :: nmfond
  character(len=9),intent(in) :: nmxpct
 end subroutine chknm8
end interface

interface
 subroutine clsopn(wff)
  use defs_datatypes
  implicit none
  type(wffile_type),intent(inout) :: wff
 end subroutine clsopn
end interface

interface
 subroutine compmesh(mesh,r_for_intg)
  use defs_basis
  use defs_datatypes
  implicit none
  type(pawrad_type),intent(inout) :: mesh
  real(dp),intent(in) :: r_for_intg
 end subroutine compmesh
end interface

interface
 subroutine copymesh(mesh1,mesh2)
  use defs_datatypes
  implicit none
  type(pawrad_type),intent(in) :: mesh1
  type(pawrad_type),intent(out) :: mesh2
 end subroutine copymesh
end interface

interface
 subroutine ctrap(imax,ff,hh,ans)
  use defs_basis
  implicit none
  integer,intent(in) :: imax
  real(dp),intent(out) :: ans
  real(dp),intent(in) :: hh
  real(dp),intent(in) :: ff(imax)
 end subroutine ctrap
end interface

interface
 subroutine ctrap_gen(intg,func,radmesh)
  use defs_basis
  use defs_datatypes
  implicit none
  real(dp),intent(out) :: intg
  type(pawrad_type),intent(in) :: radmesh
  real(dp),intent(in) :: func(radmesh%int_meshsz)
 end subroutine ctrap_gen
end interface

interface
 subroutine deducer0(func,funcsz,radmesh)
  use defs_basis
  use defs_datatypes
  implicit none
  integer :: funcsz
  type(pawrad_type),intent(in) :: radmesh
  real(dp) :: func(funcsz)
 end subroutine deducer0
end interface

interface
 subroutine derf(derf_yy,yy)
  use defs_basis
  implicit none
  real(dp),intent(out) :: derf_yy
  real(dp),intent(in) :: yy
 end subroutine derf
end interface

interface
 subroutine derfc(derfc_yy,yy)
  use defs_basis
  implicit none
  real(dp),intent(out) :: derfc_yy
  real(dp),intent(in) :: yy
 end subroutine derfc
end interface

interface
 subroutine dtsetCopy(dtout, dtin)
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: dtin
  type(dataset_type),intent(out) :: dtout
 end subroutine dtsetCopy
end interface

interface
 subroutine  tells_sizes(chosen_size,name,index,default_size,actual_size)
  implicit none
  integer,intent(in) :: actual_size
  integer,intent(out) :: chosen_size
  integer,intent(in) :: default_size
  integer,intent(in) :: index
  character(len=12),intent(in) :: name
 end subroutine tells_sizes
end interface

interface
 subroutine dtsetFree(dtset)
  use defs_datatypes
  implicit none
  type(dataset_type),intent(inout) :: dtset
 end subroutine dtsetFree
end interface

interface
 subroutine energies_init(energies)
  use defs_datatypes
  implicit none
  type(energies_type),intent(out) :: energies
 end subroutine energies_init
end interface

interface
 function factorial(nn)
  use defs_basis
  implicit none
  integer,intent(in) :: nn
  real(dp) :: factorial
 end function factorial
end interface

interface
 subroutine hermit(chmin,chmout,ierr,ndim)
  use defs_basis
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: ndim
  real(dp),intent(in) :: chmin(ndim*ndim+ndim)
  real(dp),intent(out) :: chmout(ndim*ndim+ndim)
 end subroutine hermit
end interface

interface
 function ifromr(radmesh,rr)
  use defs_basis
  use defs_datatypes
  implicit none
  integer :: ifromr
  type(pawrad_type),intent(in) :: radmesh
  real(dp),intent(in) :: rr
 end function ifromr
end interface

interface
 subroutine int2char(iint,string)
  implicit none
  integer,intent(in) :: iint
  character(len=10),intent(out) :: string
 end subroutine int2char
end interface

interface
 subroutine int2char4(iint,string)
  implicit none
  integer,intent(in) :: iint
  character(len=4),intent(out) :: string
 end subroutine int2char4
end interface

interface
 subroutine interpol3d(r,nr1,nr2,nr3,denval,grid)
  use defs_basis
  implicit none
  integer,intent(in) :: nr1
  integer,intent(in) :: nr2
  integer,intent(in) :: nr3
  real(dp),intent(out) :: denval
  real(dp),intent(in) :: grid(nr1,nr2,nr3)
  real(dp),intent(in) :: r(3)
 end subroutine interpol3d
end interface

interface
 subroutine inupper(string)
  implicit none
  character(len=*),intent(inout) :: string
 end subroutine inupper
end interface

interface
 subroutine isfile(filnam,status)
  use defs_basis
  implicit none
  character(len=fnlen),intent(inout) :: filnam
  character(len=3),intent(in) :: status
 end subroutine isfile
end interface

interface
 subroutine jbessel(bes,besp,bespp,ll,order,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: ll
  integer,intent(in) :: order
  real(dp),intent(out) :: bes
  real(dp),intent(out) :: besp
  real(dp),intent(out) :: bespp
  real(dp),intent(in) :: xx
 end subroutine jbessel
end interface

interface
 subroutine lxyz(lp,mp,idir,ll,mm,lidir)
  use defs_basis
  implicit none
  integer,intent(in) :: idir
  integer,intent(in) :: ll
  integer,intent(in) :: lp
  integer,intent(in) :: mm
  integer,intent(in) :: mp
  complex(dpc),intent(out) :: lidir
 end subroutine lxyz
end interface

interface
 subroutine matcginv(a,lda,n)
  use defs_basis
  implicit none
  integer,intent(in) :: lda
  integer,intent(in) :: n
  complex(gwpc),intent(inout) :: a(lda,n)
 end subroutine matcginv
end interface

interface
 subroutine mati3inv(mm,mit)
  implicit none
  integer,intent(out) :: mit(3,3)
  integer,intent(in) :: mm(3,3)
 end subroutine mati3inv
end interface

interface
 subroutine matr3eigval(eigval,matr)
  use defs_basis
  implicit none
  real(dp),intent(out) :: eigval(3)
  real(dp),intent(in) :: matr(3,3)
 end subroutine matr3eigval
end interface

interface
 subroutine matr3inv(aa,ait)
  use defs_basis
  implicit none
  real(dp),intent(in) :: aa(3,3)
  real(dp),intent(out) :: ait(3,3)
 end subroutine matr3inv
end interface

interface
 subroutine matrginv(a,lda,n)
  use defs_basis
  implicit none
  integer,intent(in) :: lda
  integer,intent(in) :: n
  real(dp),intent(inout) :: a(lda,n)
 end subroutine matrginv
end interface

interface
 subroutine memerr(sub_name,array_name,nelements,kindp)
  implicit none
  integer,intent(in) :: nelements
  character(len=*),intent(in) :: array_name
  character(len=*),intent(in) :: kindp
  character(len=*),intent(in) :: sub_name
 end subroutine memerr
end interface

interface
 subroutine mkfilename(filnam,filnam_out,get,idtset,&  
  &  ird,jdtset_,ndtset,stringfil,stringvar,will_read)
  use defs_basis
  implicit none
  integer,intent(in) :: get
  integer,intent(in) :: idtset
  integer,intent(in) :: ird
  integer,intent(in) :: ndtset
  integer,intent(out) :: will_read
  character(len=fnlen),intent(out) :: filnam_out
  character(len=4),intent(in) :: stringfil
  character(len=9),intent(in) :: stringvar
  character(len=fnlen),intent(in) :: filnam(5)
  integer,intent(in) :: jdtset_(0:ndtset)
 end subroutine mkfilename
end interface

interface
 subroutine mkherm(array,ndim)
  use defs_basis
  implicit none
  integer,intent(in) :: ndim
  real(dp),intent(inout) :: array(2,ndim,ndim)
 end subroutine mkherm
end interface

interface
 subroutine mknormpath(nbounds,bounds,gmet,ndiv_small,ndiv,npt_tot,path)
  use defs_basis
  implicit none
  integer,intent(in) :: nbounds
  integer,intent(in) :: ndiv_small
  integer,intent(inout) :: npt_tot
  real(dp),intent(in) :: bounds(3,nbounds)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(inout) :: ndiv(nbounds-1)
  real(dp),intent(out),optional :: path(3,npt_tot)
 end subroutine mknormpath
end interface

interface
 subroutine mvrecord(ierr,nrec,unitfile)
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: nrec
  integer,intent(in) :: unitfile
 end subroutine mvrecord
end interface

interface
 subroutine nderiv(hh,yy,zz,ndim,norder)
  use defs_basis
  implicit none
  integer,intent(in) :: ndim
  integer,intent(in) :: norder
  real(dp),intent(in) :: hh
  real(dp),intent(in) :: yy(ndim)
  real(dp),intent(out) :: zz(ndim)
 end subroutine nderiv
end interface

interface
 subroutine nderiv_gen(der,func,nder,radmesh)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nder
  type(pawrad_type),intent(in) :: radmesh
  real(dp),intent(out) :: der(radmesh%mesh_size,nder)
  real(dp),intent(in) :: func(radmesh%mesh_size)
 end subroutine nderiv_gen
end interface

interface
 subroutine normev(evec,ndim,num)
  use defs_basis
  implicit none
  integer,intent(in) :: ndim
  integer,intent(in) :: num
  real(dp),intent(inout) :: evec(2*ndim,num)
 end subroutine normev
end interface

interface
 subroutine orthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg,sqgram,vectsize)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: blocksize
  integer,intent(in) :: vectsize
  type(mpi_type) :: mpi_enreg
  real(dp) :: blockvectorbx(vectsize,blocksize)
  real(dp) :: blockvectorx(vectsize,blocksize)
  real(dp) :: sqgram(blocksize,blocksize)
 end subroutine orthonormalize
end interface

interface
 function overlap_cmplx(wf1,wf2,usepaw,cprj1,cprj2,typat,pawtab,fact) result(cres)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: usepaw
  complex(dpc) :: cres
  real(dp),intent(in),optional :: fact
  integer,intent(in) :: typat(:)
  type(cprj_type),intent(in) :: cprj1(:)
  type(cprj_type),intent(in) :: cprj2(:)
  type(pawtab_type),intent(in) :: pawtab(:)
  complex(gwpc),intent(in) :: wf1(:)
  complex(gwpc),intent(in) :: wf2(:)
 end function overlap_cmplx
end interface

interface
 function overlap_real(wf1,wf2,usepaw,cprj1,cprj2,typat,pawtab,fact) result(res)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: usepaw
  real(dp),intent(in),optional :: fact
  integer,intent(in) :: typat(:)
  type(cprj_type),intent(in) :: cprj1(:)
  type(cprj_type),intent(in) :: cprj2(:)
  type(pawtab_type),intent(in) :: pawtab(:)
  real(dp) :: res(2)
  real(dp),intent(in) :: wf1(:,:)
  real(dp),intent(in) :: wf2(:,:)
 end function overlap_real
end interface

interface
 subroutine nullify_pawfgrtab(this)
  use defs_datatypes
  implicit none
  type(pawfgrtab_type),intent(inout) :: this(:)
 end subroutine nullify_pawfgrtab
end interface

interface
 subroutine destroy_pawfgrtab(this)
  use defs_datatypes
  implicit none
  type(pawfgrtab_type),intent(inout) :: this(:)
 end subroutine destroy_pawfgrtab
end interface

interface
 subroutine init_pawfgrtab(this,l_size_atm)
  use defs_datatypes
  implicit none
  integer,intent(in) :: l_size_atm(:)
  type(pawfgrtab_type),intent(inout) :: this(:)
 end subroutine init_pawfgrtab
end interface

interface
 subroutine print_pawfgrtab(this,unitno,prtvol,mode_paral)
  use defs_datatypes
  implicit none
  integer,intent(in),optional :: prtvol
  integer,intent(in),optional :: unitno
  character(len=4),intent(in),optional :: mode_paral
  type(pawfgrtab_type),intent(inout) :: this(:)
 end subroutine print_pawfgrtab
end interface

interface
 subroutine pl_deriv(mpsang,pl_d2,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: mpsang
  real(dp),intent(in) :: xx
  real(dp),intent(out) :: pl_d2(mpsang)
 end subroutine pl_deriv
end interface

interface
 subroutine plm_coeff(blm,mpsang,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: mpsang
  real(dp),intent(in) :: xx
  real(dp),intent(out) :: blm(5,mpsang*mpsang)
 end subroutine plm_coeff
end interface

interface
 subroutine plm_d2theta(mpsang,plm_d2t,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: mpsang
  real(dp),intent(in) :: xx
  real(dp),intent(out) :: plm_d2t(mpsang*mpsang)
 end subroutine plm_d2theta
end interface

interface
 function plm_dphi(ll,mm,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: ll
  integer,intent(in) :: mm
  real(dp) :: plm_dphi
  real(dp),intent(in) :: xx
 end function plm_dphi
end interface

interface
 function plm_dtheta(ll,mm,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: ll
  integer,intent(in) :: mm
  real(dp) :: plm_dtheta
  real(dp),intent(in) :: xx
 end function plm_dtheta
end interface

interface
 subroutine poisson(den,ll,qq,radmesh,rv)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ll
  real(dp),intent(out) :: qq
  type(pawrad_type),intent(in) :: radmesh
  real(dp),intent(in) :: den(radmesh%mesh_size)
  real(dp),intent(out) :: rv(radmesh%mesh_size)
 end subroutine poisson
end interface

interface
 subroutine print_ij(a_ij,adim,cplex,ndim,opt_io,opt_l,opt_l_index,opt_pack,opt_prtvol,pack2ij,test_value,unt,&  
  &  opt_sym,asym_ij)    !Optional arguments
  use defs_basis
  implicit none
  integer,intent(in) :: adim
  integer,intent(in) :: cplex
  integer,intent(in) :: ndim
  integer,intent(in) :: opt_io
  integer,intent(in) :: opt_l
  integer,intent(in) :: opt_pack
  integer,intent(in) :: opt_prtvol
  integer,intent(in),optional :: opt_sym
  integer,intent(in) :: unt
  real(dp),intent(in) :: test_value
  real(dp),intent(in) :: a_ij(cplex*adim)
  real(dp),intent(in),optional :: asym_ij(cplex*adim)
  integer,intent(in) :: opt_l_index(ndim*min(1+opt_l,1))
  integer,intent(in) :: pack2ij(adim*opt_pack)
 end subroutine print_ij
end interface

interface
 subroutine print_ngfft(ngfft,header,unitno,mode_paral,prtvol)
  implicit none
  integer,intent(in),optional :: prtvol
  integer,intent(in),optional :: unitno
  character(len=*),intent(in),optional :: header
  character(len=4),intent(in),optional :: mode_paral
  integer,intent(in) :: ngfft(18)
 end subroutine print_ngfft
end interface

interface
 subroutine prmat (mat, ni, nj, mi)
  use defs_basis
  implicit none
  integer,intent(in) :: mi
  integer,intent(in) :: ni
  integer,intent(in) :: nj
  real(dp),intent(in) :: mat(mi,nj)
 end subroutine prmat
end interface

interface
 subroutine ratint(npts,xin,xpt,yin,yerr,ypt)
  use defs_basis
  implicit none
  integer,intent(in) :: npts
  real(dp),intent(in) :: xpt
  real(dp),intent(out) :: yerr
  real(dp),intent(out) :: ypt
  real(dp),intent(in) :: xin(npts)
  real(dp),intent(in) :: yin(npts)
 end subroutine ratint
end interface

interface
 subroutine rhoij_alloc(cplex,nlmn,nspden,nsppol,pawrhoij,typat,&  ! Mandatory arguments
  &  ngrhoij,nlmnmix,use_rhoij_,use_rhoijres) ! Optional arguments
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in),optional :: ngrhoij
  integer,intent(in),optional :: nlmnmix
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in),optional :: use_rhoij_
  integer,intent(in),optional :: use_rhoijres
  integer,intent(in) :: nlmn(:)
  integer,intent(in) :: typat(:)
  type(pawrhoij_type),intent(inout) :: pawrhoij(:)
 end subroutine rhoij_alloc
end interface

interface
 subroutine rhoij_free(pawrhoij)
  use defs_datatypes
  implicit none
  type(pawrhoij_type),intent(inout) :: pawrhoij(:)
 end subroutine rhoij_free
end interface

interface
 subroutine rhoij_copy(pawrhoij_in,pawrhoij_out)
  use defs_datatypes
  implicit none
  type(pawrhoij_type),intent(in) :: pawrhoij_in(:)
  type(pawrhoij_type),intent(inout) :: pawrhoij_out(:)
 end subroutine rhoij_copy
end interface

interface
 subroutine rotmat(xaxis,zaxis,inversion_flag,umat)
  use defs_basis
  implicit none
  integer,intent(out) :: inversion_flag
  real(dp),intent(out) :: umat(3,3)
  real(dp),intent(in) :: xaxis(3)
  real(dp),intent(in) :: zaxis(3)
 end subroutine rotmat
end interface

interface
 subroutine simp_gen(intg,func,radmesh)
  use defs_basis
  use defs_datatypes
  implicit none
  real(dp),intent(out) :: intg
  type(pawrad_type),intent(in) :: radmesh
  real(dp),intent(in) :: func(radmesh%int_meshsz)
 end subroutine simp_gen
end interface

interface
 subroutine slxyzs(lp,mp,idir,ll,mm,sls_val)
  use defs_basis
  implicit none
  integer,intent(in) :: idir
  integer,intent(in) :: ll
  integer,intent(in) :: lp
  integer,intent(in) :: mm
  integer,intent(in) :: mp
  complex(dpc),intent(out) :: sls_val
 end subroutine slxyzs
end interface

interface
 subroutine status(counter,filstat,istatr,level,routine)
  use defs_basis
  implicit none
  integer,intent(in) :: counter
  integer,intent(in) :: istatr
  integer,intent(in) :: level
  character(len=fnlen),intent(in) :: filstat
  character(len=*),intent(in) :: routine
 end subroutine status
end interface

interface
 subroutine ylm_cmplx(lx,ylm,xx,yy,zz)
  use defs_basis
  implicit none
  integer,intent(in) :: lx
  real(dp),intent(in) :: xx
  real(dp),intent(in) :: yy
  real(dp),intent(in) :: zz
  complex(dpc),intent(out) :: ylm((lx+1)*(lx+1))
 end subroutine ylm_cmplx
end interface

interface
 function ylmc(il,im,kcart)
  use defs_basis
  implicit none
  integer,intent(in) :: il
  integer,intent(in) :: im
  complex(dpc) :: ylmc
  real(dp),intent(in) :: kcart(3)
 end function ylmc
end interface

interface
 subroutine ylmcd(il,im,kcart,dth,dphi)
  use defs_basis
  implicit none
  integer,intent(in) :: il
  integer,intent(in) :: im
  complex(dpc),intent(out) :: dphi
  complex(dpc),intent(out) :: dth
  real(dp),intent(in) :: kcart(3)
 end subroutine ylmcd
end interface

interface
 subroutine ys(lp,mp,ll,mm,ys_val)
  use defs_basis
  implicit none
  integer,intent(in) :: ll
  integer,intent(in) :: lp
  integer,intent(in) :: mm
  integer,intent(in) :: mp
  complex(dpc),intent(out) :: ys_val
 end subroutine ys
end interface

interface
 subroutine zorthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg,sqgram,vectsize)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: blocksize
  integer,intent(in) :: vectsize
  type(mpi_type) :: mpi_enreg
  complex(dpc) :: blockvectorbx(vectsize,blocksize)
  complex(dpc) :: blockvectorx(vectsize,blocksize)
  complex(dpc) :: sqgram(blocksize,blocksize)
 end subroutine zorthonormalize
end interface

end module interfaces_11util
!!***
