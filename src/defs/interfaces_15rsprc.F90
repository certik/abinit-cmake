!!****m* ABINIT/interfaces_15rsprc
!! NAME
!! interfaces_15rsprc
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/15rsprc
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

module interfaces_15rsprc

 implicit none

interface
 subroutine ladielmt(atindx,atindx1,cg,dielmat,dtfil,dtset,&  
  &  eigen,&  
  &  kg,kg_diel,mpi_enreg,&  
  &  nattyp,&  
  &  npwarr,npwdiel,nspinor,&  
  &  occ,&  
  &  ph1d,psps,rhor,rprimd,&  
  &  wffnow,xred,ylm,ladiel)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, intent(in) :: npwdiel
  integer,intent(inout) :: nspinor
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(in) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp),intent(in) :: dielmat(2,npwdiel,npwdiel)
  real(dp), intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(out) :: ladiel(dtset%nfft)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp), intent(in) :: rhor(dtset%nfft,dtset%nspden)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(in) :: xred(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine ladielmt
end interface

interface
 subroutine lavnl(atindx,atindx1,cg,dtfil,dtset,&  
  &  eigen,&  
  &  kg,lavnlr,mpi_enreg,&  
  &  nattyp,&  
  &  npwarr,nspinor,&  
  &  occ,&  
  &  ph1d,psps,rhor,rprimd,&  
  &  wffnow,xred,ylm)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(inout) :: nspinor
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(in) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp), intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  real(dp), intent(out) :: lavnlr(dtset%nfft,dtset%nspden)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp), intent(inout) :: rhor(dtset%nfft,dtset%nspden)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(in) :: xred(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine lavnl
end interface

interface
 subroutine moddiel_csrb(dielar,dtset,gprimd,mpi_enreg,rdiemac,rhor_in)
  use defs_basis
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: rdiemac(dtset%nfft,dtset%nspden)
  real(dp),intent(in) :: rhor_in(dtset%nfft,dtset%nspden)
 end subroutine moddiel_csrb
end interface

interface
 subroutine prcrskerker1(dtset,mpi_enreg,nfft,nspden,ngfft,dielar,etotal,gprimd,rprimd,vresid,vrespc,natom,xred,base)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  type(dataset_type),intent(in) :: dtset
  real(dp) :: etotal
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: base(nfft)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vresid(nfft,nspden)
  real(dp),intent(out) :: vrespc(nfft,nspden)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine prcrskerker1
end interface

interface
 subroutine prcrskerker2(dtset,nfft,nspden,ngfft,dielar,gprimd,rprimd,vresid,vrespc,natom,xred,mpi_enreg,ucvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vresid(nfft,nspden)
  real(dp),intent(out) :: vrespc(nfft,nspden)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine prcrskerker2
end interface

interface
 subroutine prctfvw1(atindx,atindx1,cg,deltae,densymop_gs,dtfil,dtset,eeig,&  
  &  efermi,eigen,ek,enl,etotal,fixmom,gsqcut,intxc,irrzon,ixc,&  
  &  kg,mband,mgfft,mkmem,mpi_enreg,mpsang,mpw,natom,nattyp,nfft,nfftf,ngfftf,&  
  &  nhat,nhatgr,nhatgrdim,&  
  &  nkpt,nkxc,npwarr,nspden,nspinor,nsppol,nsym,ntypat,n3xccc,occ,occopt,optene,optres,optxc,&  
  &  pawang,pawfgr,pawfgrtab,pawtab,&  
  &  phnons,ph1d,psps,resid,rhog,rhor,rprimd,&  
  &  usexcnhat,&  
  &  vin_old,vout_unmixed,vpsp,vtrial,&  
  &  wffnow,xccc3d,xred,ylm,inlavnl)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: intxc
  integer,intent(in) :: ixc
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftf
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(inout) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: optene
  integer,intent(in) :: optres
  integer,intent(in) :: optxc
  integer,intent(in) :: usexcnhat
  real(dp), intent(in) :: deltae
  type(dens_sym_operator_type),intent(in) :: densymop_gs
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: eeig
  real(dp), intent(in) :: efermi
  real(dp),intent(out) :: ek
  real(dp),intent(out) :: enl
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: fixmom
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  integer, intent(in) :: ngfftf(18)
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  real(dp), intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp), intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp), intent(in) :: inlavnl(nfftf,nspden)
  integer, intent(in) :: irrzon(nfft**(1-1/nsym),2,nspden/nsppol)
  integer, intent(in) :: kg(3,mpw*mkmem)
  integer, intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(nfftf,nspden*psps%usepaw)
  real(dp),intent(in) :: nhatgr(nfftf,nspden,3*nhatgrdim)
  integer, intent(in) :: npwarr(nkpt)
  real(dp), intent(in) :: occ(mband*nkpt*nsppol)
  type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)
  type(pawtab_type), intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp), intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp), intent(in) :: phnons(2,nfft**(1-1/nsym),nspden/nsppol)
  real(dp), intent(out) :: resid(mband*nkpt*nsppol)
  real(dp), intent(inout) :: rhog(2,nfftf)
  real(dp), intent(inout) :: rhor(nfftf,nspden)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(inout) :: vin_old(nfftf,nspden)
  real(dp), intent(inout) :: vout_unmixed(nfftf,nspden)
  real(dp), intent(inout) :: vpsp(nfftf)
  real(dp), intent(inout) :: vtrial(nfftf,nspden)
  real(dp),dimension(:),intent(inout) :: xccc3d(n3xccc)
  real(dp), intent(in) :: xred(3,natom)
  real(dp), intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 end subroutine prctfvw1
end interface

interface
 subroutine prctfvw2(atindx,atindx1,cg,deltae,densymop_gs,dtfil,dtset,eeig,&  
  &  efermi,eigen,ek,enl,etotal,fixmom,gsqcut,intxc,irrzon,ixc,&  
  &  kg,mband,mgfft,mkmem,mpi_enreg,mpsang,mpw,natom,nattyp,nfft,nfftf,ngfftf,&  
  &  nhat,nhatgr,nhatgrdim,nkpt,nkxc,&  
  &  npwarr,nspden,nspinor,nsppol,nsym,ntypat,n3xccc,occ,occopt,optene,optres,optxc,&  
  &  pawang,pawfgr,pawfgrtab,pawtab,&  
  &  phnons,ph1d,psps,resid,rhog,rhor,rprimd,&  
  &  usexcnhat,vin_old,vout_unmixed,vpsp,vtrial,&  
  &  wffnow,xccc3d,xred,ylm)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: intxc
  integer,intent(in) :: ixc
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftf
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(inout) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: optene
  integer,intent(in) :: optres
  integer,intent(in) :: optxc
  integer,intent(in) :: usexcnhat
  real(dp), intent(in) :: deltae
  type(dens_sym_operator_type),intent(in) :: densymop_gs
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: eeig
  real(dp), intent(in) :: efermi
  real(dp),intent(out) :: ek
  real(dp),intent(out) :: enl
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: fixmom
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  integer, intent(in) :: ngfftf(18)
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  real(dp), intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp), intent(in) :: eigen(mband*nkpt*nsppol)
  integer, intent(in) :: irrzon(nfft**(1-1/nsym),2,nspden/nsppol)
  integer, intent(in) :: kg(3,mpw*mkmem)
  integer, intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(nfft,nspden*psps%usepaw)
  real(dp),intent(in) :: nhatgr(nfft,nspden,3*nhatgrdim)
  integer, intent(in) :: npwarr(nkpt)
  real(dp), intent(in) :: occ(mband*nkpt*nsppol)
  type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)
  type(pawtab_type), intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp), intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp), intent(in) :: phnons(2,nfft**(1-1/nsym),nspden/nsppol)
  real(dp), intent(out) :: resid(mband*nkpt*nsppol)
  real(dp), intent(inout) :: rhog(2,nfftf)
  real(dp), intent(inout) :: rhor(nfftf,nspden)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(inout) :: vin_old(nfftf,nspden)
  real(dp), intent(inout) :: vout_unmixed(nfftf,nspden)
  real(dp), intent(inout) :: vpsp(nfftf)
  real(dp), intent(inout) :: vtrial(nfftf,nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp), intent(in) :: xred(3,natom)
  real(dp), intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 end subroutine prctfvw2
end interface

interface
 subroutine prctfw3(deltae,dtset,&  
  &  efermi,etotal,gsqcut,&  
  &  lavnlr,mpi_enreg,&  
  &  nhat,nhatgr,nhatgrdim,&  
  &  nkxc,n3xccc,optene,optxc,&  
  &  pawang,pawfgrtab,pawtab,&  
  &  psps,rhor_in,rprimd,&  
  &  usexcnhat,&  
  &  vpsp,vresid,vrespc,vtrial,&  
  &  xccc3d,xred,istep)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: istep
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: optene
  integer,intent(in) :: optxc
  integer,intent(in) :: usexcnhat
  real(dp), intent(in) :: deltae
  type(dataset_type),intent(in) :: dtset
  real(dp), intent(in) :: efermi
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  real(dp), intent(in) :: lavnlr(dtset%nfft,dtset%nspden)
  real(dp),intent(in) :: nhat(dtset%nfft,dtset%nspden*psps%usepaw)
  real(dp),intent(in) :: nhatgr(dtset%nfft,dtset%nspden,3*nhatgrdim)
  type(pawfgrtab_type),intent(in) :: pawfgrtab(dtset%natom)
  type(pawtab_type), intent(in) :: pawtab(dtset%ntypat*psps%usepaw)
  real(dp), intent(in) :: rhor_in(dtset%nfft,dtset%nspden)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(inout) :: vpsp(dtset%nfft)
  real(dp), intent(in) :: vresid(dtset%nfft,dtset%nspden)
  real(dp), intent(out) :: vrespc(dtset%nfft,dtset%nspden)
  real(dp), intent(in) :: vtrial(dtset%nfft,dtset%nspden)
  real(dp),dimension(:),intent(inout) :: xccc3d(n3xccc)
  real(dp), intent(in) :: xred(3,dtset%natom)
 end subroutine prctfw3
end interface

end module interfaces_15rsprc
!!***
