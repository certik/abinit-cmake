!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_dyson
!! NAME
!! defs_dyson
!!
!! FUNCTION
!! to be described
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (..)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_dyson

  implicit none


 interface
  subroutine dyson_de(ikernel,kernel_diag,kernel_full,npwdiel,nspden,susmat)
   use defs_basis
   integer,intent(in) :: ikernel
   integer,intent(in) :: npwdiel
   integer,intent(in) :: nspden
   real(dp),pointer :: kernel_diag(:)
   real(dp),pointer :: kernel_full(:,:,:,:,:)
   real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
  end subroutine dyson_de
 end interface

 interface
  subroutine dyson_gl(ikernel,kernel_diag,kernel_full,npwdiel,nspden,&  
&  susd_LDG,susmat)
   use defs_basis
   integer,intent(in) :: ikernel
   integer,intent(in) :: npwdiel
   integer,intent(in) :: nspden
   real(dp),pointer :: kernel_diag(:)
   real(dp),pointer :: kernel_full(:,:,:,:,:)
   real(dp),intent(out) :: susd_LDG(npwdiel)
   real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
  end subroutine dyson_gl
 end interface

 interface
  subroutine dyson_ls(ikernel,kernel_diag,kernel_full,npwdiel,nspden,susmat)
   use defs_basis
   integer,intent(in) :: ikernel
   integer,intent(in) :: npwdiel
   integer,intent(in) :: nspden
   real(dp),pointer :: kernel_diag(:)
   real(dp),pointer :: kernel_full(:,:,:,:,:)
   real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
  end subroutine dyson_ls
 end interface

 interface
  subroutine dyson_sc(kernel_diag,npwdiel,nspden,susd_isc,susmat)
   use defs_basis
   integer,intent(in) :: npwdiel
   integer,intent(in) :: nspden
   real(dp),intent(in) :: kernel_diag(npwdiel)
   real(dp),intent(out) :: susd_isc(npwdiel)
   real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
  end subroutine dyson_sc
 end interface

 interface
  subroutine acfd_dyson(dtset,freq,gmet,gsq,idyson,ig_tiny,igsq_tiny,ikhxc,kg_diel,khxc,&  
&  ldgapp,mpi_enreg,ndyson,nfft,ngfft,npw_tiny,npwdiel,nspden,occopt,&  
&  option,rcut_coulomb,rhor,rhocut,rprimd,susd_data,suskxcrs,susmat,ucvol)
   use defs_basis
   use defs_datatypes
   integer,intent(in) :: idyson
   integer,intent(in) :: ikhxc
   integer,intent(in) :: ldgapp
   integer,intent(in) :: ndyson
   integer,intent(in) :: nfft
   integer,intent(in) :: npw_tiny
   integer,intent(in) :: npwdiel
   integer,intent(in) :: nspden
   integer,intent(in) :: occopt
   integer,intent(in) :: option
   integer,intent(in) :: suskxcrs
   type(dataset_type),intent(in) :: dtset
   real(dp),intent(in) :: freq
   type(MPI_type),intent(inout) :: mpi_enreg
   real(dp),intent(in) :: rcut_coulomb
   real(dp),intent(in) :: rhocut
   real(dp),intent(in) :: ucvol
   real(dp),intent(in) :: gmet(3,3)
   real(dp),intent(in) :: gsq(npwdiel)
   integer,intent(in) :: ig_tiny(npw_tiny,3)
   integer,intent(in) :: igsq_tiny(npw_tiny)
   integer,intent(in) :: kg_diel(3,npwdiel)
   real(dp),pointer :: khxc(:,:,:,:,:)
   integer,intent(in) :: ngfft(18)
   real(dp),intent(in) :: rhor(nfft,nspden)
   real(dp),intent(in) :: rprimd(3,3)
   real(dp),intent(out) :: susd_data(npwdiel,3)
   real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
  end subroutine acfd_dyson
 end interface

end module defs_dyson
!!***
