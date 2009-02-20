!!****m* ABINIT/interfaces_12spacepar
!! NAME
!! interfaces_12spacepar
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/12spacepar
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

module interfaces_12spacepar

 implicit none

interface
 subroutine dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw,option,vect1,vect2)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: istwf_k
  integer,intent(in) :: npw
  integer,intent(in) :: option
  real(dp),intent(out) :: doti
  real(dp),intent(out) :: dotr
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: vect1(2,npw)
  real(dp),intent(in) :: vect2(2,npw)
 end subroutine dotprod_g
end interface

interface
 subroutine dotprod_v(cplex,dotr,mpi_enreg,nfft,nspden,pot1,pot2)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  real(dp),intent(out) :: dotr
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: pot1(cplex*nfft,nspden)
  real(dp),intent(in) :: pot2(cplex*nfft,nspden)
 end subroutine dotprod_v
end interface

interface
 subroutine dotprod_vn(cplex,dens,dotr,doti,mpi_enreg,nfft,nfftot,nspden,option,pot,ucvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  real(dp),intent(out) :: doti
  real(dp),intent(out) :: dotr
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: dens(cplex*nfft,nspden)
  real(dp),intent(in) :: pot(cplex*nfft,nspden)
 end subroutine dotprod_vn
end interface

interface
 subroutine dotprodm_v(cplex,cpldot,dot,index1,index2,mpi_enreg,mult1,mult2,nfft,npot1,npot2,nspden,potarr1,potarr2)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cpldot
  integer,intent(in) :: cplex
  integer,intent(in) :: index1
  integer,intent(in) :: index2
  integer,intent(in) :: mult1
  integer,intent(in) :: mult2
  integer,intent(in) :: nfft
  integer,intent(in) :: npot1
  integer,intent(in) :: npot2
  integer,intent(in) :: nspden
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(out) :: dot(cpldot,mult1,mult2)
  real(dp),intent(in) :: potarr1(cplex*nfft,nspden,npot1)
  real(dp),intent(in) :: potarr2(cplex*nfft,nspden,npot2)
 end subroutine dotprodm_v
end interface

interface
 subroutine dotprodm_vn(cplex,cpldot,denarr,dot,id,ip,mpi_enreg,multd,multp,&  
  &  nden,nfft,nfftot,npot,nspden,potarr,ucvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cpldot
  integer,intent(in) :: cplex
  integer,intent(in) :: id
  integer,intent(in) :: ip
  integer,intent(in) :: multd
  integer,intent(in) :: multp
  integer,intent(in) :: nden
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: npot
  integer,intent(in) :: nspden
  type(mpi_type) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: denarr(cplex*nfft,nspden,nden)
  real(dp),intent(out) :: dot(cpldot,multp,multd)
  real(dp),intent(in) :: potarr(cplex*nfft,nspden,npot)
 end subroutine dotprodm_vn
end interface

interface
 subroutine matrixelmt_g(ai,ar,diag,istwf_k,mpi_enreg,needimag,npw,nspinor,vect1,vect2)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: istwf_k
  integer,intent(in) :: needimag
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  real(dp),intent(out) :: ai
  real(dp),intent(out) :: ar
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: diag(npw)
  real(dp),intent(in) :: vect1(2,npw*nspinor)
  real(dp),intent(in) :: vect2(2,npw*nspinor)
 end subroutine matrixelmt_g
end interface

interface
 subroutine mean_fftr(arraysp,meansp,mpi_enreg,nfft,nfftot,nspden)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: nspden
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: arraysp(nfft,nspden)
  real(dp),intent(out) :: meansp(nspden)
 end subroutine mean_fftr
end interface

interface
 subroutine meanvalue_g(ar,diag,filter,istwf_k,mpi_enreg,npw,nspinor,vect)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: filter
  integer,intent(in) :: istwf_k
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  real(dp),intent(out) :: ar
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: diag(npw)
  real(dp),intent(in) :: vect(2,npw*nspinor)
 end subroutine meanvalue_g
end interface

interface
 subroutine overlap_g(doti,dotr,mpw,npw_k1,npw_k2,pwind_k,vect1,vect2)
  use defs_basis
  implicit none
  integer,intent(in) :: mpw
  integer,intent(in) :: npw_k1
  integer,intent(in) :: npw_k2
  real(dp),intent(out) :: doti
  real(dp),intent(out) :: dotr
  integer,intent(in) :: pwind_k(mpw)
  real(dp),intent(in) :: vect1(1:2,0:mpw)
  real(dp),intent(in) :: vect2(1:2,0:mpw)
 end subroutine overlap_g
end interface

interface
 subroutine sqnorm_g(dotr,istwf_k,mpi_enreg,npwsp,vect)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: istwf_k
  integer,intent(in) :: npwsp
  real(dp),intent(out) :: dotr
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: vect(2,npwsp)
 end subroutine sqnorm_g
end interface

interface
 subroutine sqnorm_v(cplex,mpi_enreg,nfft,norm2,nspden,pot)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(out) :: norm2
  real(dp),intent(in) :: pot(cplex*nfft,nspden)
 end subroutine sqnorm_v
end interface

interface
 subroutine sqnormm_v(cplex,index,mpi_enreg,mult,nfft,norm2,npot,nspden,potarr)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: index
  integer,intent(in) :: mult
  integer,intent(in) :: nfft
  integer,intent(in) :: npot
  integer,intent(in) :: nspden
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(out) :: norm2(mult)
  real(dp),intent(in) :: potarr(cplex*nfft,nspden,npot)
 end subroutine sqnormm_v
end interface

end module interfaces_12spacepar
!!***
