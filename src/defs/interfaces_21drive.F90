!!****m* ABINIT/interfaces_21drive
!! NAME
!! interfaces_21drive
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/21drive
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

module interfaces_21drive

 implicit none

interface
 subroutine afterscfloop(atindx,atindx1,cg,computed_forces,cprj,cpus,&  
  &  dimcprj,deltae,diffor,dtefield,dtfil,dtset,eigen,energies,etotal,&  
  &  favg,fcart,filapp,filfft,forold,fred,gresid,grewtn,grhf,&  
  &  grxc,gsqcut,hdr,indsym,&  
  &  istep,kg,kxc,maxfor,mgfftc,mgfftf,&  
  &  moved_atm_inside,mpi_enreg,&  
  &  n3xccc,nattyp,&  
  &  nfftf,ngfft,ngfftf,nhat,nkxc,npwarr,nvresid,&  
  &  occ,optres,optxc,paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,pel,pel_cg,&  
  &  ph1d,ph1df,pion,prtfor,psps,pwind,pwind_alloc,pwnsfac,res2,resid,residm,results_gs,&  
  &  rhocore,rhog,rhor,rhore,rhototp,&  
  &  rprimd,stress_needed,strsxc,strten,symrec,synlgr,tollist,usecprj,usexcnhat,&  
  &  vhartr,vpsp,vxc,vxcavg,wffnow,wvl,xccc3d,xred,xred_old,ylm,ylmgr)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(inout) :: computed_forces
  integer,intent(in) :: istep
  integer,intent(in) :: mgfftc
  integer,intent(in) :: mgfftf
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfftf
  integer,intent(in) :: nkxc
  integer,intent(in) :: optres
  integer,intent(in) :: optxc
  integer,intent(in) :: prtfor
  integer,intent(in) :: pwind_alloc
  integer,intent(in) :: stress_needed
  integer,intent(in) :: usecprj
  integer,intent(in) :: usexcnhat
  real(dp),intent(in) :: cpus
  real(dp),intent(in) :: deltae
  real(dp),intent(inout) :: diffor
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  type(energies_type),intent(inout) :: energies
  real(dp),intent(inout) :: etotal
  character(len=fnlen),intent(in) :: filapp
  character(len=fnlen),intent(in) :: filfft
  real(dp),intent(in) :: gsqcut
  type(hdr_type),intent(inout) :: hdr
  real(dp),intent(inout) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: res2
  real(dp),intent(in) :: residm
  type(results_gs_type),intent(inout) :: results_gs
  real(dp),intent(inout) :: vxcavg
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data),intent(inout) :: wvl
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: cg(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  type(cprj_type),intent(in) :: cprj(dtset%natom,dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol*usecprj)
  integer,intent(in) :: dimcprj(dtset%natom*usecprj)
  real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(out) :: favg(3)
  real(dp),intent(out) :: fcart(3,dtset%natom)
  real(dp),intent(inout) :: forold(3,dtset%natom)
  real(dp),intent(out) :: fred(3,dtset%natom)
  real(dp),intent(out) :: gresid(3,dtset%natom)
  real(dp),intent(in) :: grewtn(3,dtset%natom)
  real(dp),intent(out) :: grhf(3,dtset%natom)
  real(dp),intent(out) :: grxc(3,dtset%natom)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  real(dp),intent(out) :: kxc(nfftf,nkxc)
  integer,intent(in) :: nattyp(dtset%ntypat)
  real(dp),intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(inout) :: nvresid(nfftf,dtset%nspden)
  real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(dtset%natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%ntypat*psps%usepaw)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)
  real(dp),intent(inout) :: pel(3)
  real(dp),intent(in) :: pel_cg(3)
  real(dp),intent(inout) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(inout) :: ph1df(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(inout) :: pion(3)
  integer,intent(in) :: pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(in) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(in) :: rhocore(nfftf)
  real(dp),intent(inout) :: rhog(2,nfftf)
  real(dp),intent(inout) :: rhor(nfftf,dtset%nspden)
  real(dp),intent(in) :: rhore(nfftf,dtset%nspden)
  real(dp),intent(in) :: rhototp(nfftf)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: strsxc(6)
  real(dp),intent(out) :: strten(6)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(out) :: synlgr(3,dtset%natom)
  real(dp),intent(in) :: tollist(12)
  real(dp),intent(inout) :: vhartr(dtset%nfft)
  real(dp),intent(in) :: vpsp(nfftf)
  real(dp),intent(inout) :: vxc(nfftf,dtset%nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(out) :: xred_old(3,dtset%natom)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine afterscfloop
end interface

interface
 subroutine brdmin(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&  
  &  dtset,ecore,eigen,hdr,indsym,initialized,irrzon,&  
  &  kg,mpi_enreg,mxfh,&  
  &  nattyp,nfftf,npwarr,nspinor,nxfh,occ,&  
  &  pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,&  
  &  rhog,rhor,rprim,scf_history,symrec,wffnew,wffnow,vel,wvl,&  
  &  xfhist,xred,xred_old,ylm,ylmgr)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(inout) :: initialized
  integer,intent(in) :: mxfh
  integer,intent(inout) :: nfftf
  integer,intent(inout) :: nspinor
  integer,intent(inout) :: nxfh
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(dens_sym_operator_type),intent(inout) :: densymop_gs
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecore
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(results_gs_type),intent(out) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data),intent(inout) :: wvl
  real(dp), intent(inout) :: acell(3)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp), intent(out) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
  integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(inout) :: rprim(3,3)
  integer, intent(inout) :: symrec(3,3,dtset%nsym)
  real(dp), intent(in) :: vel(3,dtset%natom)
  real(dp), intent(inout) :: xfhist(3,dtset%natom+4,2,mxfh)
  real(dp), intent(inout) :: xred(3,dtset%natom)
  real(dp), intent(inout) :: xred_old(3,dtset%natom)
  real(dp), intent(inout) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(inout) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine brdmin
end interface

interface
 subroutine delocint(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&  
  &  dtset,ecore,eigen,hdr,indsym,initialized,irrzon,&  
  &  kg,mpi_enreg,mxfh,&  
  &  nattyp,nfftf,npwarr,nspinor,nxfh,occ,&  
  &  pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,&  
  &  rhog,rhor,rprim,scf_history,symrec,wffnew,wffnow,vel,wvl,&  
  &  xfhist,xred,xred_old,ylm,ylmgr)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(inout) :: initialized
  integer,intent(in) :: mxfh
  integer,intent(inout) :: nfftf
  integer,intent(inout) :: nspinor
  integer,intent(inout) :: nxfh
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(dens_sym_operator_type),intent(inout) :: densymop_gs
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecore
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(results_gs_type),intent(out) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data),intent(inout) :: wvl
  real(dp), intent(inout) :: acell(3)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp), intent(out) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
  integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(inout) :: rprim(3,3)
  integer, intent(inout) :: symrec(3,3,dtset%nsym)
  real(dp), intent(in) :: vel(3,dtset%natom)
  real(dp), intent(inout) :: xfhist(3,dtset%natom+4,2,mxfh)
  real(dp), intent(inout) :: xred(3,dtset%natom)
  real(dp), intent(inout) :: xred_old(3,dtset%natom)
  real(dp), intent(inout) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(inout) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine delocint
end interface

interface
 subroutine diisRelax(acell, atindx, atindx1, cg, cpus, densymop_gs, dtefield,&  
  &  dtfil, dtset, ecore, eigen, hdr, iapp, indsym, initialized,&  
  &  irrzon, kg, mpi_enreg, nattyp, nfftf, npwarr, nspinor, occ, pawang,&  
  &  pawfgr, pawrad, pawrhoij, pawtab, phnons, psps, pwind, pwind_alloc, pwnsfac,&  
  &  resid, results_gs, rhog, rhor, rprimd, scf_history, symrec,&  
  &  wffnew, wffnow, wvl, xred, xred_old, ylm, ylmgr)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: iapp
  integer,intent(inout) :: initialized
  integer, intent(inout) :: nfftf
  integer,intent(inout) :: nspinor
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(dens_sym_operator_type),intent(inout) :: densymop_gs
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecore
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(results_gs_type),intent(inout) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data), intent(inout) :: wvl
  real(dp), intent(inout) :: acell(3)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
  integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(inout) :: rprimd(3,3)
  integer, intent(inout) :: symrec(3,3,dtset%nsym)
  real(dp), intent(inout) :: xred(3,dtset%natom)
  real(dp), intent(inout) :: xred_old(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine diisRelax
end interface

interface
 subroutine driver(codvsn,cpui,dtfil,dtsets,filnam,filstat,&  
  &  mpi_enreg,ndtset,ndtset_alloc,npsp,pspheads,results_out,walli)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: npsp
  character(len=6),intent(in) :: codvsn
  real(dp),intent(in) :: cpui
  type(datafiles_type),intent(inout) :: dtfil
  character(len=fnlen),intent(in) :: filstat
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: walli
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  character(len=fnlen),intent(in) :: filnam(5)
  type(pspheader_type),intent(in) :: pspheads(npsp)
  type(results_out_type),intent(inout) :: results_out(0:ndtset_alloc)
 end subroutine driver
end interface

interface
 subroutine elpolariz(atindx1,cg,cprj,dtefield,dtfil,dtset,etotal,enefield,gprimd,hdr,&  
  &  kg,mband,mgfft,mkmem,mpi_enreg,mpw,natom,nattyp,nkpt,&  
  &  npwarr,nspinor,nsppol,ntypat,pawang,pawrad,pawtab,pel,pel_cg,pelev,pion,psps,pwind,pwind_alloc,&  
  &  pwnsfac,rprimd,ucvol,usecprj,wffnow,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: pwind_alloc
  integer,intent(in) :: usecprj
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(inout) :: enefield
  real(dp),intent(inout) :: etotal
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  type(cprj_type),intent(in) :: cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: npwarr(nkpt)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)
  real(dp),intent(inout) :: pel(3)
  real(dp),intent(in) :: pel_cg(3)
  real(dp),intent(inout) :: pelev(3)
  real(dp),intent(inout) :: pion(3)
  integer,intent(in) :: pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine elpolariz
end interface

interface
 subroutine gstate(acell,codvsn,cpui,dtfil,dtset,iexit,&  
  &  mpi_enreg,&  
  &  npwtot,nspinor,&  
  &  occ,pawang,pawrad,pawtab,psps,results_gs,rprim,vel,walli,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(inout) :: iexit
  integer,intent(inout) :: nspinor
  character(len=6),intent(in) :: codvsn
  real(dp),intent(in) :: cpui
  type(datafiles_type),intent(inout) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(inout) :: pawang
  type(pseudopotential_type),intent(inout) :: psps
  type(results_gs_type),intent(inout) :: results_gs
  real(dp),intent(in) :: walli
  real(dp),intent(inout) :: acell(3)
  integer,intent(out) :: npwtot(dtset%nkpt)
  real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(inout) :: rprim(3,3)
  real(dp),intent(inout) :: vel(3,dtset%natom)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine gstate
end interface

interface
 subroutine isotemp(amass,dtion,dtset,ekin,ktemp,mttk_vars,vel)
  use defs_basis
  use defs_datatypes
  implicit none
  real(dp),intent(in) :: dtion
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(out) :: ekin
  real(dp),intent(in) :: ktemp
  type(mttk_type) :: mttk_vars
  real(dp),intent(in) :: amass(dtset%natom)
  real(dp),intent(inout) :: vel(3,dtset%natom)
 end subroutine isotemp
end interface

interface
 subroutine isopress(amass,dtion,dtset,ekin,ktemp,strten,strtarget,ucvol,mttk_vars,vel,vlogv)
  use defs_basis
  use defs_datatypes
  implicit none
  real(dp),intent(in) :: dtion
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(out) :: ekin
  real(dp),intent(in) :: ktemp
  type(mttk_type) :: mttk_vars
  real(dp),intent(in) :: ucvol
  real(dp),intent(inout) :: vlogv
  real(dp),intent(in) :: amass(dtset%natom)
  real(dp),intent(inout) :: strtarget(6)
  real(dp),intent(inout) :: strten(6)
  real(dp),intent(inout) :: vel(3,dtset%natom)
 end subroutine isopress
end interface

interface
 subroutine isostress(amass,dtion,dtset,ekin,ktemp,strten,strtarget,ucvol,vel,mttk_vars)
  use defs_basis
  use defs_datatypes
  implicit none
  real(dp),intent(in) :: dtion
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(out) :: ekin
  real(dp),intent(in) :: ktemp
  type(mttk_type) :: mttk_vars
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: amass(dtset%natom)
  real(dp),intent(in) :: strtarget(6)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(inout) :: vel(3,dtset%natom)
 end subroutine isostress
end interface

interface
 subroutine loop3dte(blkflg,cg,cgindex,dtfil,dtset,d3lo,&  
  &  etotal,gmet,gprimd,gsqcut,gsqcut_eff,&  
  &  hdr,kg,kneigh,kptindex,kpt3,kxc,k3xc,mband,mgfft,mkmem,mkmem_max,mk1mem,&  
  &  mpert,mpi_enreg,mpw,mvwtk,natom,nfft,nkpt,nkpt3,nkxc,nneigh,nspinor,nsppol,&  
  &  npwarr,occ,psps,pwind,&  
  &  rfpert,rmet,rprimd,tmpfil,ucvol,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mkmem_max
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkpt3
  integer,intent(in) :: nkxc
  integer,intent(in) :: nneigh
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(inout) :: etotal
  real(dp),intent(in) :: gsqcut
  real(dp),intent(in) :: gsqcut_eff
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(out) :: blkflg(3,mpert,3,mpert,3,mpert)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  integer,intent(in) :: cgindex(nkpt,nsppol)
  real(dp),intent(out) :: d3lo(2,3,mpert,3,mpert,3,mpert)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: k3xc(nfft)
  integer,intent(in) :: kg(3,mk1mem*mpw)
  integer,intent(in) :: kneigh(30,nkpt)
  real(dp),intent(in) :: kpt3(3,nkpt3)
  integer,intent(in) :: kptindex(2,nkpt3)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: mvwtk(30,nkpt)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
  integer,intent(in) :: pwind(mpw,nneigh,mkmem)
  integer,intent(in) :: rfpert(3,mpert,3,mpert,3,mpert)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  character(len=fnlen),intent(in) :: tmpfil(15)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine loop3dte
end interface

interface
 subroutine moldyn(acell,amass,atindx,atindx1,cg,cpus,densymop_gs,&  
  &  dtefield,dtfil,dtset,ecore,eigen,hdr,indsym,initialized,&  
  &  irrzon,kg,mpi_enreg,mxfh,nattyp,nfftf,npwarr,nspinor,nxfh,occ,&  
  &  pawang,pawfgr,pawrad,pawrhoij,pawtab,&  
  &  phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprim,&  
  &  scf_history,symrec,wffnew,wffnow,vel,wvl,xfhist,xred,xred_old,ylm,ylmgr)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(inout) :: initialized
  integer,intent(in) :: mxfh
  integer,intent(inout) :: nfftf
  integer,intent(inout) :: nspinor
  integer,intent(inout) :: nxfh
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(dens_sym_operator_type),intent(inout) :: densymop_gs
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecore
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(results_gs_type),intent(out) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data),intent(inout) :: wvl
  real(dp),intent(inout) :: acell(3)
  real(dp),intent(in) :: amass(dtset%natom)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp),intent(out) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
  integer,intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: nattyp(psps%ntypat)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
  integer,intent(in) :: pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),pointer :: rhog(:,:)
  real(dp),pointer :: rhor(:,:)
  real(dp),intent(inout) :: rprim(3,3)
  integer,intent(inout) :: symrec(3,3,dtset%nsym)
  real(dp),intent(inout) :: vel(3,dtset%natom)
  real(dp),intent(inout) :: xfhist(3,dtset%natom+4,2,mxfh)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(inout) :: xred_old(3,dtset%natom)
  real(dp),intent(inout) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(inout) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine moldyn
end interface

interface
 subroutine move(acell,amass,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,&  
  &  ecore,eigen,hdr,indsym,initialized,irrzon,&  
  &  kg,mpi_enreg,&  
  &  nattyp,nfftf,npwarr,nspinor,occ,&  
  &  pawang,pawfgr,pawrad,pawrhoij,pawtab,&  
  &  phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,&  
  &  scf_history,symrec,wffnew,wffnow,vel,wvl,xred,xred_old,ylm,ylmgr)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(inout) :: initialized
  integer,intent(inout) :: nfftf
  integer,intent(inout) :: nspinor
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(dens_sym_operator_type),intent(inout) :: densymop_gs
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecore
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(results_gs_type),intent(out) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data),intent(inout) :: wvl
  real(dp),intent(inout) :: acell(3)
  real(dp), intent(in) :: amass(dtset%natom)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp), intent(out) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
  integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(inout) :: rprimd(3,3)
  integer, intent(inout) :: symrec(3,3,dtset%nsym)
  real(dp), intent(inout) :: vel(3,dtset%natom)
  real(dp), intent(inout) :: xred(3,dtset%natom)
  real(dp), intent(inout) :: xred_old(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine move
end interface

interface
 subroutine nonlinear(codvsn,dtfil,dtset,etotal,iexit,&  
  &  mband,mgfft,mkmem,mpi_enreg,mpw,natom,nfft,nkpt,npwtot,nspden,&  
  &  nspinor,nsppol,nsym,occ,pawrad,pawtab,psps,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(inout) :: natom
  integer,intent(in) :: nfft
  integer,intent(inout) :: nkpt
  integer,intent(inout) :: nspden
  integer,intent(inout) :: nspinor
  integer,intent(inout) :: nsppol
  integer,intent(inout) :: nsym
  character(len=6),intent(in) :: codvsn
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(inout) :: etotal
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(inout) :: psps
  integer,intent(out) :: npwtot(nkpt)
  real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
  type(pawrad_type),intent(inout) :: pawrad(psps%ntypat,psps%usepaw)
  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat,psps%usepaw)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine nonlinear
end interface

interface
 subroutine outscfcv(atindx,atindx1,cg,compch_fft,compch_sph,cprj,dimcprj,dtfil,dtset,ecut,eigen,etotal,&  
  &  fermie,filapp,gmet,gprimd,gsqcut,hdr,kg,&  
  &  kssform,mband,mgfftc,mkmem,mpi_enreg,mpsang,mpw,natom,&  
  &  nattyp,nfft,ngfft,nhat,nkpt,npwarr,nspden,nspinor,nsppol,nsym,ntypat,n3xccc,occ,&  
  &  pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,paw_ij,ph1dc,prtvol,psps,rhog,rhor,rmet,rprimd,&  
  &  ucvol,usecprj,usexcnhat,wffnow,vhartr,vtrial,vxc,vxcavg,xccc3d,xred,ylm)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: kssform
  integer,intent(in) :: mband
  integer,intent(in) :: mgfftc
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(inout) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  integer,intent(in) :: usecprj
  integer,intent(in) :: usexcnhat
  real(dp),intent(in) :: compch_fft
  real(dp),intent(in) :: compch_sph
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecut
  real(dp),intent(inout) :: etotal
  real(dp),intent(in) :: fermie
  character(len=fnlen),intent(inout) :: filapp
  real(dp),intent(in) :: gsqcut
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  real(dp),intent(inout) :: vxcavg
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(inout) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  type(cprj_type),intent(in) :: cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)
  integer,intent(in) :: dimcprj(natom*usecprj)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(inout) :: nhat(nfft,nspden*psps%usepaw)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(paw_ij_type),intent(inout) :: paw_ij(natom)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(in) :: ph1dc(2,3*(2*mgfftc+1)*natom)
  real(dp),intent(inout) :: rhog(nfft,nspden)
  real(dp),intent(inout) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(inout) :: vtrial(nfft,nspden)
  real(dp),intent(inout) :: vxc(nfft,nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 end subroutine outscfcv
end interface

interface
 subroutine scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&  
  &  dtset,ecore,eigen,hdr,iapp,indsym,initialized,&  
  &  irrzon,kg,mpi_enreg,nattyp,nfftf,npwarr,nspinor,occ,&  
  &  pawang,pawfgr,pawrad,pawrhoij,pawtab,&  
  &  phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,&  
  &  scf_history,symrec,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: iapp
  integer,intent(inout) :: initialized
  integer,intent(inout) :: nfftf
  integer,intent(inout) :: nspinor
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(dens_sym_operator_type),intent(inout) :: densymop_gs
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecore
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(inout) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(results_gs_type),intent(inout) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data),intent(inout) :: wvl
  real(dp), intent(inout) :: acell(3)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
  integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(inout) :: rprimd(3,3)
  integer, intent(inout) :: symrec(3,3,dtset%nsym)
  real(dp), intent(inout) :: xred(3,dtset%natom)
  real(dp), intent(inout) :: xred_old(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine scfcv
end interface

interface
 subroutine screening(acell,codvsn,Dtfil,Dtset,iexit,MPI_enreg,Pawang,Pawrad,Pawtab,Psps,rprim,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(inout) :: iexit
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(inout) :: Dtset
  type(mpi_type),intent(inout) :: MPI_enreg
  type(pawang_type),intent(inout) :: Pawang
  type(pseudopotential_type),intent(inout) :: Psps
  character(len=6),intent(in) :: codvsn
  type(pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Dtset%usepaw)
  type(pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Dtset%usepaw)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: xred(3,Dtset%natom)
 end subroutine screening
end interface

interface
 subroutine sigma(acell,codvsn,Dtfil,Dtset,iexit,MPI_enreg,Pawang,Pawrad,Pawtab,Psps,rprim,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(inout) :: iexit
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(inout) :: Dtset
  type(mpi_type),intent(inout) :: MPI_enreg
  type(pawang_type),intent(inout) :: Pawang
  type(pseudopotential_type),intent(inout) :: Psps
  character(len=6),intent(in) :: codvsn
  type(pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Psps%usepaw)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: xred(3,Dtset%natom)
 end subroutine sigma
end interface

interface
 subroutine testfi(etotal,filnam,filstat,fred,natom,strten,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  real(dp),intent(in) :: etotal
  character(len=fnlen),intent(in) :: filstat
  character(len=fnlen),intent(in) :: filnam(5)
  real(dp),intent(in) :: fred(3,natom)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine testfi
end interface

interface
 subroutine timana(mpi_enreg,natom,nband,ndtset,nfft,nkpt,npwtot,nsppol,timopt, papiopt)
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ndtset
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: papiopt
  integer,intent(in) :: timopt
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwtot(nkpt)
 end subroutine timana
end interface

interface
 subroutine wannier(dtfil,dtset,iexit,mband,mpi_enreg,nkpt,nsppol)
  use defs_datatypes
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine wannier
end interface

end module interfaces_21drive
!!***
