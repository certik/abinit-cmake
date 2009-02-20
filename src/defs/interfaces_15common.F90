!!****m* ABINIT/interfaces_15common
!! NAME
!! interfaces_15common
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/15common
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

module interfaces_15common

 implicit none

interface
 subroutine aprxdr(cplex,choice,dedv_mix,dedv_new,dedv_old,&  
  &  f_atm,f_fftgr,i_rhor2,i_vresid,moved_atm_inside,&  
  &  mpi_enreg,natom,nfft,nfftot,nspden,n_fftgr,rhor,ucvol,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: cplex
  integer,intent(in) :: i_rhor2
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: n_fftgr
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: nspden
  real(dp),intent(out) :: dedv_mix
  real(dp),intent(out) :: dedv_new
  real(dp),intent(out) :: dedv_old
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: i_vresid(3)
  real(dp),intent(in) :: f_atm(3,natom,n_fftgr)
  real(dp),intent(in) :: f_fftgr(cplex*nfft,nspden,n_fftgr)
  real(dp),intent(in) :: rhor(cplex*nfft,nspden)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine aprxdr
end interface

interface
 subroutine atm2fft(atindx1,atmrho,atmvloc,dyfrn,dyfrv,eei,gauss,gmet,gprimd,grn,grv,gsqcut,&  
  &  mgfft,mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,ntypat,&  
  &  optatm,optdyfr,optgr,optn,optn2,optstr,optv,paral_kgb,&  
  &  pawtab,ph1d,qgrid,qprtrb,rhog,strn,strv,ucvol,usepaw,vg,vprtrb,vspl)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: optatm
  integer,intent(in) :: optdyfr
  integer,intent(in) :: optgr
  integer,intent(in) :: optn
  integer,intent(in) :: optn2
  integer,intent(in) :: optstr
  integer,intent(in) :: optv
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: eei
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: qprtrb(3)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(out) :: atmrho(nfft*optn)
  real(dp),intent(out) :: atmvloc(nfft*optv)
  real(dp),intent(out) :: dyfrn(3,3,natom*optn*optdyfr)
  real(dp),intent(out) :: dyfrv(3,3,natom*optv*optdyfr)
  real(dp),intent(in) :: gauss(2,ntypat*(optn2/3))
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: grn(3,natom*optn*optgr)
  real(dp),intent(out) :: grv(3,natom*optv*optgr)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rhog(2,nfft*optv*max(optgr,optstr,optdyfr))
  real(dp),intent(out) :: strn(6*optn*optstr)
  real(dp),intent(out) :: strv(6*optv*optstr)
  real(dp),intent(in) :: vg(2,nfft*optn*max(optgr,optstr,optdyfr))
  real(dp),intent(in) :: vprtrb(2)
  real(dp),intent(in) :: vspl(mqgrid,2,ntypat*optv)
 end subroutine atm2fft
end interface

interface
 subroutine berryphase(atindx1,bdberry,cg,gprimd,istwfk,kberry,kg,kpt_,&  
  &  kptopt,kptrlatt,mband,&  
  &  mgfft,mkmem,mpi_enreg,mpw,natom,nattyp,nband,nberry,npwarr,nspinor,nsppol,ntypat,&  
  &  nkpt_,rprimd,typat,ucvol,unkg,wffnow,xred,zion)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: kptopt
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nberry
  integer,intent(in) :: nkpt_
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: unkg
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: bdberry(4)
  integer,intent(in) :: kptrlatt(3,3)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: gprimd(1:3,1:3)
  integer,intent(in) :: istwfk(nkpt_)
  integer,intent(in) :: kberry(3,nberry)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kpt_(3,nkpt_)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: nband(nkpt_*nsppol)
  integer,intent(in) :: npwarr(nkpt_)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine berryphase
end interface

interface
 subroutine brdene(etotal,etotal_prev,hessin,ndim,vin,vin_prev,vout,vout_prev)
  use defs_basis
  implicit none
  integer,intent(in) :: ndim
  real(dp),intent(in) :: etotal
  real(dp),intent(inout) :: etotal_prev
  real(dp),intent(in) :: hessin(ndim,ndim)
  real(dp),intent(inout) :: vin(ndim)
  real(dp),intent(inout) :: vin_prev(ndim)
  real(dp),intent(in) :: vout(ndim)
  real(dp),intent(inout) :: vout_prev(ndim)
 end subroutine brdene
end interface

interface
 subroutine bstruct_clean(bstruct)
  use defs_datatypes
  implicit none
  type(bandstructure_type),intent(inout) :: bstruct
 end subroutine bstruct_clean
end interface

interface
 subroutine bstruct_init(bantot,bstruct,doccde,eig,istwfk,kptns,&  
  &  nband,nkpt,npwarr,nsppol,occ,wtk)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: bantot
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(bandstructure_type),intent(out) :: bstruct
  real(dp),intent(in) :: doccde(bantot)
  real(dp),intent(in) :: eig(bantot)
  integer,intent(in) :: istwfk(nkpt)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(bantot)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine bstruct_init
end interface

interface
 subroutine calc_cs(corecs,natom,nspden,ntypat,occopt,pawang,pawrad,pawrhoij,pawtab,prtcs,typat,usepaw)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: prtcs
  integer,intent(in) :: usepaw
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: corecs(ntypat)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine calc_cs
end interface

interface
 subroutine calc_efg(gprimd,natom,nfft,ngfft,nhat,nspden,ntypat,paral_kgb,pawang,pawrad,pawrhoij,pawtab,&  
  &  ptcharge,prtefg,quadmom,rhor,rprimd,typat,ucvol,usepaw,xred,zion)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtefg
  integer,intent(in) :: usepaw
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: nhat(nfft,nspden*usepaw)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ptcharge(ntypat)
  real(dp),intent(in) :: quadmom(ntypat)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine calc_efg
end interface

interface
 subroutine calc_fc(gprimd,natom,nfft,ngfft,nhat,nspden,ntypat,paral_kgb,&  
  &  pawrad,pawrhoij,pawtab,psps,rhor,rprimd,typat,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  type(pseudopotential_type),intent(in) :: psps
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: nhat(nfft,nspden*psps%usepaw)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine calc_fc
end interface

interface
 subroutine clnup1(acell,dosdeltae,dtset,eigen,enunit,&  
  &  fermie,filnam,fred,iatfix,iscf,kptns,kptopt,mband,mkmem,mpi_enreg,mpw,&  
  &  natom,nband,nfft,ngfft,nkpt,nspden,nspinor,nsppol,nstep,occ,occopt,prtdos,&  
  &  prteig,prtfor,prtstm,prtvol,resid,rhor,rprimd,tphysel,tsmear,vxcavg,wtk,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: enunit
  integer,intent(in) :: iscf
  integer,intent(in) :: kptopt
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nstep
  integer,intent(in) :: occopt
  integer,intent(in) :: prtdos
  integer,intent(in) :: prteig
  integer,intent(in) :: prtfor
  integer,intent(in) :: prtstm
  integer,intent(in) :: prtvol
  real(dp),intent(in) :: dosdeltae
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: fermie
  character(len=fnlen),intent(in) :: filnam
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: tphysel
  real(dp),intent(in) :: tsmear
  real(dp),intent(in) :: vxcavg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: fred(3,natom)
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: resid(mband*nkpt*nsppol)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: wtk(nkpt)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine clnup1
end interface

interface
 subroutine clnup2(n1xccc,fred,gresid,grewtn,grxc,iscf,natom,prtfor,prtstr,prtvol,&  
  &  start,strten,synlgr,usepaw,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: iscf
  integer,intent(in) :: n1xccc
  integer,intent(in) :: natom
  integer,intent(in) :: prtfor
  integer,intent(in) :: prtstr
  integer,intent(in) :: prtvol
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: fred(3,natom)
  real(dp),intent(in) :: gresid(3,natom)
  real(dp),intent(in) :: grewtn(3,natom)
  real(dp),intent(in) :: grxc(3,natom)
  real(dp),intent(in) :: start(3,natom)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(in) :: synlgr(3,natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine clnup2
end interface

interface
 subroutine conducti_nc(filnam,mpi_enreg)
  use defs_basis
  use defs_datatypes
  implicit none
  character(len=fnlen) :: filnam
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine conducti_nc
end interface

interface
 subroutine conducti_paw(filnam,mpi_enreg)
  use defs_basis
  use defs_datatypes
  implicit none
  character(len=fnlen) :: filnam
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine conducti_paw
end interface

interface
 subroutine constrf(diffor,fcart,forold,fred,iatfix,ionmov,maxfor,natom,&  
  &  nconeq,prtvol,rprimd,wtatcon,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: ionmov
  integer,intent(in) :: natom
  integer,intent(in) :: nconeq
  integer,intent(in) :: prtvol
  real(dp),intent(out) :: diffor
  real(dp),intent(out) :: maxfor
  real(dp),intent(inout) :: fcart(3,natom)
  real(dp),intent(inout) :: forold(3,natom)
  real(dp),intent(out) :: fred(3,natom)
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: wtatcon(3,natom,nconeq)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine constrf
end interface

interface
 subroutine dielmt(dielar,dielinv,dielop,gmet,kg_diel,&  
  &  npwdiel,nspden,occopt,prtvol,susmat)
  use defs_basis
  implicit none
  integer,intent(in) :: dielop
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: occopt
  integer,intent(in) :: prtvol
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(out) :: dielinv(2,npwdiel,nspden,npwdiel,nspden)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 end subroutine dielmt
end interface

interface
 subroutine dielmt2(dielar,dielop,gmet,kg_diel,&  
  &  npwdiel,nspden,occopt,prtvol,susmat,&  
  &  dieldiag,dtset,ucvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: dielop
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: occopt
  integer,intent(in) :: prtvol
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(out) :: dieldiag(2,npwdiel,nspden)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 end subroutine dielmt2
end interface

interface
 subroutine dieltcel(dielar,dielinv,dielop,gmet,kg_diel,kxc,&  
  &  nfft,ngfft,nkxc,npwdiel,nspden,occopt,option,paral_kgb,prtvol,susmat)
  use defs_basis
  implicit none
  integer,intent(in) :: dielop
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: occopt
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(out) :: dielinv(2,npwdiel,nspden,npwdiel,nspden)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 end subroutine dieltcel
end interface

interface
 subroutine dsksta(ishm,nbandkss,npwkss,nkpt,nsym2)
  implicit none
  integer,intent(in) :: ishm
  integer,intent(in) :: nbandkss
  integer,intent(in) :: nkpt
  integer,intent(in) :: npwkss
  integer,intent(in) :: nsym2
 end subroutine dsksta
end interface

interface
 subroutine energy(atindx,atindx1,cg,compch_fft,densymop_gs,dtfil,dtset,energies,&  
  &  eigen,etotal,gsqcut,indsym,irrzon,kg,mpi_enreg,nattyp,nfftf,ngfft,ngfftf,nhat,&  
  &  nhatgr,nhatgrdim,npwarr,nspinor,n3xccc,occ,optene,paw_ij,pawang,pawfgr,&  
  &  pawfgrtab,pawrhoij,pawtab,phnons,ph1d,psps,resid,rhog,rhor,rprimd,strsxc,symrec,&  
  &  usexcnhat,vhartr,vtrial,vpsp,vxc,wffnow,wfs,xccc3d,xred,ylm)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfftf
  integer,intent(in) :: nhatgrdim
  integer,intent(inout) :: nspinor
  integer,intent(in) :: optene
  integer,intent(in) :: usexcnhat
  real(dp),intent(out) :: compch_fft
  type(dens_sym_operator_type),intent(in) :: densymop_gs
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(energies_type),intent(inout) :: energies
  real(dp),intent(out) :: etotal
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_wf_type),intent(in) :: wfs
  integer, intent(in) :: ngfft(18)
  integer, intent(in) :: ngfftf(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(in) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp), intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  integer :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
  integer :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  real(dp), intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw)
  real(dp),intent(in) :: nhatgr(nfftf,dtset%nspden,3*nhatgrdim)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(paw_ij_type), intent(in) :: paw_ij(dtset%natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(inout) :: rhog(2,nfftf)
  real(dp), intent(inout) :: rhor(nfftf,dtset%nspden)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(out) :: strsxc(6)
  integer, intent(in) :: symrec(3,3,dtset%nsym)
  real(dp), intent(out) :: vhartr(nfftf)
  real(dp), intent(in) :: vpsp(nfftf)
  real(dp), intent(out) :: vtrial(nfftf,dtset%nspden)
  real(dp), intent(out) :: vxc(nfftf,dtset%nspden)
  real(dp), intent(in) :: xccc3d(n3xccc)
  real(dp), intent(in) :: xred(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine energy
end interface

interface
 subroutine etheta(bcut,chc,detovc,detovd,dhc,dhd,efield_dot,e0,e1,&  
  &  hel,nkpt,nsppol,nstr,sdeg,theta)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  real(dp),intent(in) :: chc
  real(dp),intent(in) :: dhc
  real(dp),intent(in) :: dhd
  real(dp),intent(out) :: e0
  real(dp),intent(out) :: e1
  real(dp),intent(in) :: sdeg
  real(dp),intent(in) :: theta
  integer,intent(in) :: hel(2,3)
  integer,intent(in) :: nstr(3)
  real(dp),intent(in) :: bcut(2,3)
  real(dp),intent(in) :: detovc(2,2,3)
  real(dp),intent(in) :: detovd(2,2,3)
  real(dp),intent(in) :: efield_dot(3)
 end subroutine etheta
end interface

interface
 subroutine etotfor(atindx1,deltae,diffor,dtset,efield_dot,elast,energies,&  
  &  etotal,favg,fcart,forold,fred,gresid,grewtn,grhf,grnl,&  
  &  grxc,gsqcut,indsym,kxc,maxfor,mgfft,mpi_enreg,nattyp,&  
  &  nfft,ngfft,nhat,nkxc,ntypat,nvresid,n1xccc,n3xccc,optene,optforces,optres,&  
  &  pawang,pawfgrtab,pawrhoij,pawtab,pel,ph1d,pion,psps,rhog,rhor,rprimd,symrec,synlgr,&  
  &  ucvol,usepaw,usexcnhat,vhartr,vpsp,vxc,xccc3d,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: ntypat
  integer,intent(in) :: optene
  integer,intent(in) :: optforces
  integer,intent(in) :: optres
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  real(dp),intent(out) :: deltae
  real(dp),intent(out) :: diffor
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(inout) :: elast
  type(energies_type),intent(inout) :: energies
  real(dp),intent(out) :: etotal
  real(dp),intent(in) :: gsqcut
  real(dp),intent(out) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: efield_dot(3)
  real(dp),intent(out) :: favg(3)
  real(dp),intent(out) :: fcart(3,dtset%natom)
  real(dp),intent(inout) :: forold(3,dtset%natom)
  real(dp),intent(out) :: fred(3,dtset%natom)
  real(dp),intent(out) :: gresid(3,dtset%natom)
  real(dp),intent(in) :: grewtn(3,dtset%natom)
  real(dp),intent(out) :: grhf(3,dtset%natom)
  real(dp),intent(inout) :: grnl(3*dtset%natom)
  real(dp),intent(out) :: grxc(3,dtset%natom)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(inout) :: nhat(nfft,dtset%nspden*psps%usepaw)
  real(dp),intent(inout) :: nvresid(nfft,dtset%nspden)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(in) :: pel(3)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: pion(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(out) :: synlgr(3,dtset%natom)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(in) :: vpsp(nfft)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine etotfor
end interface

interface
 subroutine ewald(eew,gmet,grewtn,natom,ntypat,rmet,typat,ucvol,xred,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp),intent(out) :: eew
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(out) :: grewtn(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine ewald
end interface

interface
 subroutine ewald2(gmet,natom,ntypat,rmet,rprimd,stress,&  
  &  typat,ucvol,xred,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: stress(6)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine ewald2
end interface

interface
 subroutine extraprho(atindx1,dtset,gmet,gprimd,gsqcut,mgfft,mpi_enreg,mqgrid,nattyp,&  
  &  nfft,ngfft,ntypat,pawrhoij,pawtab,ph1d,qgrid,rhor,rprimd,scf_history,ucvol,&  
  &  usepaw,xred_new,xred_old,zion,znucl)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(scf_history_type),intent(inout) :: scf_history
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: nattyp(ntypat)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred_new(3,dtset%natom)
  real(dp),intent(inout) :: xred_old(3,dtset%natom)
  real(dp),intent(in) :: zion(ntypat)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine extraprho
end interface

interface
 subroutine fconv(fcart,iatfix,iexit,itime,natom,ntime,&  
  &  optcell,strfact,strtarget,strten,tolmxf)
  use defs_basis
  implicit none
  integer,intent(inout) :: iexit
  integer,intent(in) :: itime
  integer,intent(in) :: natom
  integer,intent(in) :: ntime
  integer,intent(in) :: optcell
  real(dp),intent(in) :: strfact
  real(dp),intent(in) :: tolmxf
  real(dp),intent(in) :: fcart(3,natom)
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: strtarget(6)
  real(dp),intent(in) :: strten(6)
 end subroutine fconv
end interface

interface
 subroutine filterpot(cplex,gmet,gsqcut,nfft,ngfft,nspden,paral_kgb,qphon,vpot)
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(inout) :: vpot(cplex*nfft,nspden)
 end subroutine filterpot
end interface

interface
 subroutine findmin(choice,dedv_1,dedv_2,dedv_predict,&  
  &  d2edv2_1,d2edv2_2,d2edv2_predict,&  
  &  etotal_1,etotal_2,etotal_predict,&  
  &  lambda_1,lambda_2,lambda_predict,status)
  use defs_basis
  implicit none
  integer,intent(in) :: choice
  integer,intent(out) :: status
  real(dp),intent(out) :: d2edv2_1
  real(dp),intent(out) :: d2edv2_2
  real(dp),intent(out) :: d2edv2_predict
  real(dp),intent(in) :: dedv_1
  real(dp),intent(in) :: dedv_2
  real(dp),intent(out) :: dedv_predict
  real(dp),intent(in) :: etotal_1
  real(dp),intent(in) :: etotal_2
  real(dp),intent(out) :: etotal_predict
  real(dp),intent(in) :: lambda_1
  real(dp),intent(in) :: lambda_2
  real(dp),intent(out) :: lambda_predict
 end subroutine findmin
end interface

interface
 subroutine fixsym(iatfix,indsym,natom,nsym)
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: iatfix(3,natom)
  integer,intent(in) :: indsym(4,nsym,natom)
 end subroutine fixsym
end interface

interface
 subroutine forces(atindx1,diffor,dtset,favg,fcart,forold,fred,gresid,grewtn,&  
  &  grhf,grnl,grxc,gsqcut,indsym,kxc,&  
  &  maxfor,mgfft,mpi_enreg,n1xccc,n3xccc,&  
  &  nattyp,nfft,ngfft,nkxc,ntypat,&  
  &  pawtab,ph1d,psps,rhog,rhor,rprimd,symrec,synlgr,&  
  &  vresid,vxc,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: ntypat
  real(dp),intent(out) :: diffor
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: gsqcut
  real(dp),intent(out) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(out) :: favg(3)
  real(dp),intent(inout) :: fcart(3,dtset%natom)
  real(dp),intent(inout) :: forold(3,dtset%natom)
  real(dp),intent(out) :: fred(3,dtset%natom)
  real(dp),intent(out) :: gresid(3,dtset%natom)
  real(dp),intent(in) :: grewtn(3,dtset%natom)
  real(dp),intent(out) :: grhf(3,dtset%natom)
  real(dp),intent(in) :: grnl(3*dtset%natom)
  real(dp),intent(out) :: grxc(3,dtset%natom)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(out) :: synlgr(3,dtset%natom)
  real(dp),intent(inout) :: vresid(nfft,dtset%nspden)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine forces
end interface

interface
 subroutine forstr(atindx1,cg,diffor,dtset,eigen,energies,favg,fcart,&  
  &  forold,fred,gresid,grewtn,grhf,grxc,gsqcut,indsym,&  
  &  kg,kxc,maxfor,mgfftf,mpi_enreg,n3xccc,nattyp,&  
  &  nfftf,ngfftf,nhat,nkxc,npwarr,nspinor,&  
  &  ntypat,nvresid,occ,optfor,optres,paw_ij,pawang,pawfgr,&  
  &  pawfgrtab,pawrhoij,pawtab,pel,ph1d,ph1df,pion,psps,rhog,rhor,rprimd,stress_needed,&  
  &  strsxc,strten,symrec,synlgr,ucvol,unkg,unylm,usexcnhat,vhartr,vpsp,&  
  &  vxc,wffnow,xccc3d,xred,ylm,ylmgr)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mgfftf
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfftf
  integer,intent(in) :: nkxc
  integer,intent(inout) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: optfor
  integer,intent(in) :: optres
  integer,intent(in) :: stress_needed
  integer,intent(in) :: unkg
  integer,intent(in) :: unylm
  integer,intent(in) :: usexcnhat
  real(dp),intent(inout) :: diffor
  type(dataset_type),intent(in) :: dtset
  type(energies_type),intent(in) :: energies
  real(dp),intent(in) :: gsqcut
  real(dp),intent(inout) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
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
  real(dp),intent(in) :: kxc(dtset%nfft,nkxc)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(inout) :: nvresid(nfftf,dtset%nspden)
  real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(dtset%natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(in) :: pel(3)
  real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*dtset%natom)
  real(dp),intent(in) :: pion(3)
  real(dp),intent(in) :: rhog(2,nfftf)
  real(dp),intent(inout) :: rhor(nfftf,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: strsxc(6)
  real(dp),intent(out) :: strten(6)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(out) :: synlgr(3,dtset%natom)
  real(dp),intent(in) :: vhartr(nfftf)
  real(dp),intent(in) :: vpsp(nfftf)
  real(dp),intent(in) :: vxc(nfftf,dtset%nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine forstr
end interface

interface
 subroutine forstrnps(atindx1,cg,ecut,ecutsm,effmass,eigen,&  
  &  grnl,indsym,istwfk,kg,kinstr,npsstr,kpt,mband,mgfft,mkmem,mpi_enreg,mpsang,&  
  &  mpw,natom,nattyp,nband,nfft,ngfft,nkpt,nloalg,npwarr,nspden,nspinor,nsppol,nsym,&  
  &  ntypat,occ,occopt,optfor,paw_ij,pawang,pawprtvol,pawtab,ph1d,psps,rprimd,&  
  &  stress_needed,symafm,symrec,typat,unkg,unylm,wffnow,wtk,xred,ylm,ylmgr)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(inout) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: optfor
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: stress_needed
  integer,intent(in) :: unkg
  integer,intent(in) :: unylm
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  real(dp),intent(in) :: effmass
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(5)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(out) :: grnl(3*natom)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(out) :: kinstr(6)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(out) :: npsstr(6)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(ntypat)
  real(dp),intent(in) :: wtk(nkpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang*psps%useylm)
 end subroutine forstrnps
end interface

interface
 subroutine fresid(dtset,gmet,gresid,gsqcut,mpi_enreg,nfft,ngfft,ntypat,option,&  
  &  pawtab,rhor,rprimd,ucvol,work,xred_new,xred_old,znucl)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(out) :: gresid(3,dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat*dtset%usepaw)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: work(nfft,dtset%nspden)
  real(dp),intent(in) :: xred_new(3,dtset%natom)
  real(dp),intent(in) :: xred_old(3,dtset%natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine fresid
end interface

interface
 subroutine fresidrsp(atindx1,dtset,gmet,gprimd,gresid,gsqcut,mgfft,mpi_enreg,mqgrid,nattyp,nfft,&  
  &  ngfft,ntypat,pawtab,ph1d,qgrid,ucvol,usepaw,vresid,zion,znucl)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: gresid(3,dtset%natom)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: vresid(nfft,dtset%nspden)
  real(dp),intent(in) :: zion(ntypat)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine fresidrsp
end interface

interface
 subroutine getgsc(atindx,atindx1,cg,dtfil,dtset,gmet,gprimd,gsc,&  
  &  kg,mband,mgfft,mkmem,mpi_enreg,mpsang,mpw,natom,&  
  &  nattyp,nkpt,npwarr,nspinor,nsppol,ntypat,option,pawtab,ph1d,psps,rmet,&  
  &  swfnow,tim_getgsc,ucvol,wffnew,wffnow,xred,ylm)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(inout) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  integer,intent(in) :: tim_getgsc
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: swfnow
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(inout) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: gsc(2,mpw*nspinor*mband*mkmem*nsppol*psps%usepaw)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: npwarr(nkpt)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang)
 end subroutine getgsc
end interface

interface
 subroutine gipaw_j_dia_bare(jdia,nfft,ngfft,nhat,nspden,rhor,rprimd)
  use defs_basis
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: jdia(3,3,nfft)
  real(dp),intent(in) :: nhat(nfft,nspden)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine gipaw_j_dia_bare
end interface

interface
 subroutine hartre1(cplex,gmet,gsqcut,nfft,ngfft,paral_kgb,qphon,rhog,vhartr,&  
  &  sum,rcut,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  real(dp),intent(in) :: rcut
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(out) :: sum(3)
  real(dp),intent(out) :: vhartr(cplex*nfft)
 end subroutine hartre1
end interface

interface
 subroutine initberry(dtefield,dtfil,dtset,gmet,kg,mband,mkmem,mpi_enreg,&  
  &  mpw,nkpt,npwarr,nsppol,nsym,occ,pwind,pwind_alloc,pwnsfac,&  
  &  rprimd,symrec)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(out) :: pwind_alloc
  type(efield_type),intent(out) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,pointer :: pwind(:,:,:)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),pointer :: pwnsfac(:,:)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine initberry
end interface

interface
 subroutine initmv(cgindex,dtfil,dtset,gmet,kg,kneigh,kptindex,&  
  &  kpt3,mband,mkmem,mkmem_max,mpi_enreg,mpw,nband,nkpt2,&  
  &  nkpt3,nneigh,npwarr,nsppol,occ,pwind)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mkmem_max
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt2
  integer,intent(in) :: nkpt3
  integer,intent(in) :: nneigh
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(out) :: cgindex(nkpt2,nsppol)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kneigh(30,nkpt2)
  real(dp),intent(in) :: kpt3(3,nkpt3)
  integer,intent(in) :: kptindex(2,nkpt3)
  integer,intent(in) :: nband(nkpt2)
  integer,intent(in) :: npwarr(nkpt2)
  real(dp),intent(in) :: occ(mband*nkpt2*nsppol)
  integer,intent(out) :: pwind(mpw,nneigh,mkmem)
 end subroutine initmv
end interface

interface
 subroutine initro(atindx,densty,gmet,gsqcut,izero,mgfft,mpi_enreg,mqgrid,natom,nattyp,&  
  &  nfft,ngfft,nspden,ntypat,paral_kgb,pawtab,ph1d,qgrid,rhog,rhor,spinat,ucvol,usepaw,zion,znucl)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: izero
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: gsqcut
  type(mpi_type) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(natom)
  real(dp),intent(in) :: densty(ntypat,4)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(out) :: rhog(2,nfft)
  real(dp),intent(out) :: rhor(nfft,nspden)
  real(dp),intent(in) :: spinat(3,natom)
  real(dp),intent(in) :: zion(ntypat)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine initro
end interface

interface
 subroutine ionion_realSpace(dtset, eew, grewtn, rprimd, xred, zion)
  use defs_basis
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: eew
  real(dp),intent(out) :: grewtn(3,dtset%natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(in) :: zion(dtset%ntypat)
 end subroutine ionion_realSpace
end interface

interface
 subroutine jellium(gmet,gsqcut,mpi_enreg,nfft,ngfft,nspden,&  
  &  option,paral_kgb,slabwsrad,rhog,rhor,rprimd,vjell,slabzstart,slabzend)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: slabwsrad
  real(dp),intent(in) :: slabzend
  real(dp),intent(in) :: slabzstart
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(inout) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,min(option,nspden))
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vjell(nfft)
 end subroutine jellium
end interface

interface
 subroutine jvec_to_B(cshield,gcart,jvec,natom,nfft,ngfft,paral_kgb,rprimd,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: cshield(3,3,natom)
  real(dp),intent(in) :: gcart(ngfft(1),ngfft(2),ngfft(3),3)
  real(dp),intent(in) :: jvec(3,3,nfft)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine jvec_to_B
end interface

interface
 subroutine linear_optics_paw(filnam,mpi_enreg)
  use defs_basis
  use defs_datatypes
  implicit none
  character(len=fnlen),intent(in) :: filnam
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine linear_optics_paw
end interface

interface
 subroutine linemin(bcut,chc,costh,detovc,detovd,dhc,dhd,dphase_aux1,&  
  &  efield_dot,iline,nkpt,nsppol,nstr,hel,phase_end,phase_init,sdeg,sinth,thetam)
  use defs_basis
  implicit none
  integer,intent(in) :: iline
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  real(dp),intent(in) :: chc
  real(dp),intent(out) :: costh
  real(dp),intent(in) :: dhc
  real(dp),intent(in) :: dhd
  real(dp),intent(in) :: sdeg
  real(dp),intent(out) :: sinth
  real(dp),intent(out) :: thetam
  integer,intent(out) :: hel(2,3)
  integer,intent(in) :: nstr(3)
  real(dp),intent(out) :: bcut(2,3)
  real(dp),intent(in) :: detovc(2,2,3)
  real(dp),intent(in) :: detovd(2,2,3)
  real(dp),intent(inout) :: dphase_aux1(3)
  real(dp),intent(in) :: efield_dot(3)
  real(dp),intent(out) :: phase_end(3)
  real(dp),intent(inout) :: phase_init(3)
 end subroutine linemin
end interface

interface
 subroutine linopt(nspin,omega,nkpt,wkpt,nsymcrys,symcrys,nstval,occv,evalv,efermi,pmat,&  
  v1,v2,nmesh,de,sc,brod,fnam)
  use defs_basis
  implicit none
  integer, intent(in) :: nkpt
  integer, intent(in) :: nmesh
  integer, intent(in) :: nspin
  integer, intent(in) :: nstval
  integer, intent(in) :: nsymcrys
  integer, intent(in) :: v1
  integer, intent(in) :: v2
  real(dp), intent(in) :: brod
  real(dp), intent(in) :: de
  real(dp), intent(in) :: efermi
  character(256), intent(in) :: fnam
  real(dp), intent(in) :: omega
  real(dp), intent(in) :: sc
  real(dp), intent(in) :: evalv(nstval,nspin,nkpt)
  real(dp), intent(in) :: occv(nstval,nspin,nkpt)
  complex(dpc), intent(in) :: pmat(nstval,nstval,nkpt,3,nspin)
  real(dp), intent(in) :: symcrys(3,3,nsymcrys)
  real(dp), intent(in) :: wkpt(nkpt)
 end subroutine linopt
end interface

interface
 subroutine make_efg_el(efg,gcart,natom,nfft,ngfft,nhat,nspden,paral_kgb,rhor,rprimd,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: efg(3,3,natom)
  real(dp),intent(in) :: gcart(ngfft(1),ngfft(2),ngfft(3),3)
  real(dp),intent(in) :: nhat(nfft,nspden)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine make_efg_el
end interface

interface
 subroutine make_efg_ion(efg,gcart,natom,ngfft,ntypat,rprimd,typat,ucvol,xred,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: efg(3,3,natom)
  real(dp),intent(in) :: gcart(ngfft(1),ngfft(2),ngfft(3),3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine make_efg_ion
end interface

interface
 subroutine make_fc_el(fc,gcart,natom,nfft,ngfft,nhat,nspden,paral_kgb,rhor,rprimd,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: fc(natom)
  real(dp),intent(in) :: gcart(ngfft(1),ngfft(2),ngfft(3),3)
  real(dp),intent(in) :: nhat(nfft,nspden)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine make_fc_el
end interface

interface
 subroutine memkss(mband,mgfft,mkmem,mpi_enreg,mproj,mpsang,mpw,natom,&  
  &  ngfft,nkpt,nspinor,nsym,ntypat)
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mproj
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
 end subroutine memkss
end interface

interface
 subroutine mklocl(dtset, dyfrlo,eei,gmet,gprimd,grtn,gsqcut,lpsstr,mgfft,&  
  &  mpi_enreg,natom,nattyp,nfft,ngfft,nspden,ntypat,option,ph1d,psps,qprtrb,&  
  &  rhog,rhor,rmet,rprimd,ucvol,vprtrb,vpsp,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: eei
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: qprtrb(3)
  real(dp),intent(out) :: dyfrlo(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: grtn(3,natom)
  real(dp),intent(out) :: lpsstr(6)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vprtrb(2)
  real(dp),intent(out) :: vpsp(nfft)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine mklocl
end interface

interface
 subroutine mklocl_realspace(dtset, grtn, mgfft, mpi_enreg, natom, nattyp, nfft, ngfft,&  
  &  nspden, ntypat, option, ph1d, psps, rhog, rhor, rmet,&  
  &  rprimd, ucvol, vpsp, xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: grtn(3,natom)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vpsp(nfft)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine mklocl_realspace
end interface

interface
 subroutine createIonicPotential_new(geocode,iproc,nproc,nat,ntypes,iatype,psppar,nelpsp,rxyz,gridcart,&  
  hxh,hyh,hzh,elecfield,n1i,n2i,n3i,pkernel,pot_ion,spaceworld)
  implicit none
  integer, intent(in) :: iproc
  integer, intent(in) :: n1i
  integer, intent(in) :: n2i
  integer, intent(in) :: n3i
  integer, intent(in) :: nat
  integer, intent(in) :: nproc
  integer, intent(in) :: ntypes
  integer, intent(in) :: spaceworld
  real(kind=8), intent(in) :: elecfield
  character(len=1), intent(in) :: geocode
  real(kind=8), intent(in) :: hxh
  real(kind=8), intent(in) :: hyh
  real(kind=8), intent(in) :: hzh
  real(kind=8), dimension(*), intent(in) :: pkernel
  real(kind=8), dimension(*), intent(inout) :: pot_ion
  real(kind=8), dimension(3,n1i*n2i*n3i), intent(in) :: gridcart
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
 end subroutine createIonicPotential_new
end interface

interface
 subroutine ind_positions(periodic,i,n,j,go)
  implicit none
  integer, intent(in) :: i
  integer, intent(out) :: j
  integer, intent(in) :: n
  logical, intent(out) :: go
  logical, intent(in) :: periodic
 end subroutine ind_positions
end interface

interface
 subroutine ext_buffers(periodic,nl,nr)
  implicit none
  integer, intent(out) :: nl
  integer, intent(out) :: nr
  logical, intent(in) :: periodic
 end subroutine ext_buffers
end interface

interface
 subroutine local_forces_new(geocode,iproc,nproc,ntypes,nat,iatype,rxyz,gridcart,psppar,nelpsp,hxh,hyh,hzh,&  
  n1,n2,n3,rho,pot,floc)
  implicit none
  integer, intent(in) :: iproc
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: nat
  integer, intent(in) :: nproc
  integer, intent(in) :: ntypes
  character(len=1), intent(in) :: geocode
  real(kind=8), intent(in) :: hxh
  real(kind=8), intent(in) :: hyh
  real(kind=8), intent(in) :: hzh
  real(kind=8), dimension(*), intent(in) :: pot
  real(kind=8), dimension(*), intent(in) :: rho
  real(kind=8), dimension(3,nat), intent(out) :: floc
  real(kind=8), dimension(3,n1*n2*n3), intent(in) :: gridcart
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
 end subroutine local_forces_new
end interface

interface
 subroutine mklocl_recipspace(dyfrlo,eei,gmet,gprimd,grtn,gsqcut,lpsstr,mgfft,&  
  &  mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,ntypat,option,paral_kgb,ph1d,qgrid,qprtrb,&  
  &  rhog,ucvol,vlspl,vprtrb,vpsp)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: eei
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: qprtrb(3)
  real(dp),intent(out) :: dyfrlo(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: grtn(3,natom)
  real(dp),intent(out) :: lpsstr(6)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: vlspl(mqgrid,2,ntypat)
  real(dp),intent(in) :: vprtrb(2)
  real(dp),intent(out) :: vpsp(nfft)
 end subroutine mklocl_recipspace
end interface

interface
 subroutine mklocl_wavelets(dtset, grtn, mpi_enreg, option, psps, rhor, rprimd,&  
  &  vpsp, xcart)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: option
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(inout) :: grtn(3,dtset%natom)
  real(dp),intent(in) :: rhor(dtset%nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: vpsp(dtset%nfft)
  real(dp),intent(inout) :: xcart(3,dtset%natom)
 end subroutine mklocl_wavelets
end interface

interface
 subroutine mkresi(cg,dimffnl,eig_k,ffnl,filstat,gs_hamk,icg,ikpt,isppol,kg_k,kinpw,lmnmax,&  
  &  matblk,mcg,mgfft,mpi_enreg,mpsang,mpssoang,natom,nband,npw,nspinor,ntypat,nvloc,n4,n5,n6,&  
  &  paral_kgb,ph3d,prtvol,resid_k,usepaw,vlocal)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: dimffnl
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: nvloc
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: usepaw
  character(len=fnlen),intent(in) :: filstat
  type(gs_hamiltonian_type),intent(in) :: gs_hamk
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(out) :: eig_k(nband)
  real(dp),intent(in) :: ffnl(npw,dimffnl,lmnmax,ntypat)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kinpw(npw)
  real(dp),intent(inout) :: ph3d(2,npw,matblk)
  real(dp),intent(out) :: resid_k(nband)
  real(dp),intent(inout) :: vlocal(n4,n5,n6,nvloc)
 end subroutine mkresi
end interface

interface
 subroutine mkrho(cg,densymop_gs,dtset,irrzon,kg,mpi_enreg,&  
  &  npwarr,nspinor,occ,phnons,rhog,rhor,tim_mkrho,ucvol,unkg,wffnow,wfs)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(inout) :: nspinor
  integer,intent(in) :: tim_mkrho
  integer,intent(in) :: unkg
  type(dens_sym_operator_type),intent(in) :: densymop_gs
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_wf_type),intent(in) :: wfs
  real(dp), intent(in) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(in) :: phnons(2,(dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3))**(1-1/dtset%nsym), &
  &         dtset%nspden/dtset%nsppol)
  real(dp), intent(out) :: rhog(2,dtset%nfft)
  real(dp), intent(out) :: rhor(dtset%nfft,dtset%nspden)
 end subroutine mkrho
end interface

interface
 subroutine mksubham(cg,ghc,ghc_block,gsc,gvnlc,gvnlc_block,iblock,icg,igsc,ikpt,isppol,istwf_k,&  
  &  isubh,isubo,mcg,mgsc,mpi_enreg,nband_k,nbdblock,npw_k,&  
  &  nspinor,subham,subovl,subvnl,usepaw,use_subovl,use_vnl,wfoptalg)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iblock
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: istwf_k
  integer,intent(inout) :: isubh
  integer,intent(inout) :: isubo
  integer,intent(in) :: mcg
  integer,intent(in) :: mgsc
  integer,intent(in) :: nband_k
  integer,intent(in) :: nbdblock
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: use_subovl
  integer,intent(in) :: use_vnl
  integer,intent(in) :: usepaw
  integer,intent(in) :: wfoptalg
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(inout) :: ghc(2,npw_k*nspinor)
  real(dp),intent(in) :: ghc_block(2,npw_k*nspinor,nbdblock)
  real(dp),intent(in) :: gsc(2,mgsc)
  real(dp),intent(inout) :: gvnlc(2,npw_k*nspinor)
  real(dp),intent(in) :: gvnlc_block(2,npw_k*nspinor,nbdblock*use_vnl)
  real(dp),intent(inout) :: subham(nband_k*(nband_k+1))
  real(dp),intent(inout) :: subovl(nband_k*(nband_k+1)*use_subovl)
  real(dp),intent(inout) :: subvnl(nband_k*(nband_k+1)*use_vnl)
 end subroutine mksubham
end interface

interface
 subroutine mlwfovlp(cg,cprj,dtset,ecut,eigen,fermie,gprimd,gsqcut,kg,&  
  &  mband,mgfftc,mkmem,mpi_enreg,mpsang,mpw,natom,&  
  &  nattyp,nfft,ngfft,nkpt,npwarr,nspden,nspinor,nsppol,ntypat,&  
  &  pawang,pawrad,pawtab,prtvol,psps,rprimd,ucvol,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mgfftc
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(inout) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: fermie
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer :: ngfft(18)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  type(cprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  integer :: kg(3,mpw*mkmem)
  integer :: nattyp(ntypat)
  integer :: npwarr(nkpt)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mlwfovlp
end interface

interface
 subroutine mlwfovlp_proj(A_matrix,band_in,cg,cprj,dtset,eigen,gprimd,kg,&  
  &  istwfk,iwav,lproj,mband,mkmem,mpi_enreg,mpw,natom,nattyp,&  
  &  nkpt,npwarr,nspden,nspinor,&  
  &  nsppol,ntypat,num_bands,nwan,pawtab,proj_l,proj_m,proj_radial,&  
  &  proj_site,proj_x,proj_z,proj_zona,prtvol,psps,ucvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: lproj
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: num_bands
  integer,intent(in) :: nwan
  integer,intent(in) :: prtvol
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp) :: ucvol
  complex(dpc),intent(out) :: A_matrix(num_bands,nwan,nkpt)
  logical,intent(in) :: band_in(mband)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  type(cprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: istwfk(nkpt)
  integer :: iwav(nsppol,nkpt,mband)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer :: nattyp(ntypat)
  integer,intent(in) :: npwarr(nkpt)
  type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  integer,intent(in) :: proj_l(mband)
  integer,intent(in) :: proj_m(mband)
  integer,intent(in) :: proj_radial(mband)
  real(dp),intent(in) :: proj_site(3,mband)
  real(dp),intent(in) :: proj_x(3,mband)
  real(dp),intent(in) :: proj_z(3,mband)
  real(dp),intent(in) :: proj_zona(mband)
 end subroutine mlwfovlp_proj
end interface

interface
 subroutine mlwfovlp_pw(cg,cm1,dtset,g1,iwav,kg,mband,mbandw,mkmem,mpsang,mpw,natom,&  
  &  nfft,ngfft,nkpt,nntot,npwarr,nspden,nspinor,nsppol,ntypat,ovikp,prtvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mbandw
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nntot
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  type(dataset_type),intent(in) :: dtset
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(out) :: cm1(2,nkpt,nntot,mband,mband)
  integer,intent(in) :: g1(3,nkpt,nntot)
  integer,intent(out) :: iwav(nsppol,nkpt,mband)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: ovikp(nkpt,nntot)
 end subroutine mlwfovlp_pw
end interface

interface
 subroutine mlwfovlp_radial(alpha,lmax,lmax2,radial,rvalue,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: lmax
  integer,intent(in) :: lmax2
  integer,intent(in) :: rvalue
  real(dp),intent(in) :: alpha
  real(dp),intent(in) :: xx
  real(dp),intent(out) :: radial(lmax2)
 end subroutine mlwfovlp_radial
end interface

interface
 subroutine mlwfovlp_setup(atom_symbols,band_in,dtset,eigen,gamma_only,&  
  &  g1,gprimd,lwanniersetup,mband,mbandw,mkmem,mpw,natom,nattyp,nband_inc,nkpt,&  
  &  nntot,nsppol,nspinor,ntypat,num_bands,num_nnmax,nwan,ovikp,&  
  &  proj_l,proj_m,proj_radial,proj_site,proj_x,proj_z,proj_zona,&  
  &  real_lattice,recip_lattice,rprimd,spinors,xcart,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: lwanniersetup
  integer,intent(in) :: mband
  integer,intent(out) :: mbandw
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(out) :: nband_inc
  integer,intent(in) :: nkpt
  integer,intent(out) :: nntot
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(out) :: num_bands
  integer,intent(in) :: num_nnmax
  integer,intent(out) :: nwan
  type(dataset_type),intent(in) :: dtset
  logical,intent(in) :: gamma_only
  logical,intent(in) :: spinors
  character(len=3),intent(out) :: atom_symbols(natom)
  logical,intent(out) :: band_in(mband)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(out) :: g1(3,nkpt,num_nnmax)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(out) :: ovikp(nkpt,num_nnmax)
  integer,intent(out) :: proj_l(mband)
  integer,intent(out) :: proj_m(mband)
  integer,intent(out) :: proj_radial(mband)
  real(dp),intent(out) :: proj_site(3,mband)
  real(dp),intent(out) :: proj_x(3,mband)
  real(dp),intent(out) :: proj_z(3,mband)
  real(dp),intent(out) :: proj_zona(mband)
  real(dp),intent(in) :: real_lattice(3,3)
  real(dp),intent(in) :: recip_lattice(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: xcart(3,natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mlwfovlp_setup
end interface

interface
 subroutine moddiel(cplex,dielar,mpi_enreg,nfft,ngfft,nspden,optreal,paral_kgb,qphon,rprimd,vresid,vrespc)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: optreal
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vresid(cplex*nfft,nspden)
  real(dp),intent(out) :: vrespc(cplex*nfft,nspden)
 end subroutine moddiel
end interface

interface
 subroutine msig (fcti,npti, xi)
  use defs_basis
  implicit none
  integer :: npti
  real(dp) :: fcti(npti)
  real(dp) :: xi(npti)
 end subroutine msig
end interface

interface
 subroutine newkpt(ceksp2,cg,debug,doorth,ecut1,ecut2,ecut2_eff,eigen,exchn2n3d,fill,&  
  &  formeig,gmet1,gmet2,headform1,indkk,iout,ireadwf,istwfk1,istwfk2,&  
  &  kg2,kptns1,kptns2,mband2,mcg,mkmem1,mkmem2,mpi_enreg,mpw1,mpw2,&  
  &  nband1,nband2,ngfft,nkpt1,nkpt2,npwarr1,npwarr2,nspinor1,nspinor2,&  
  &  nsppol1,nsppol2,nsym,occ,optorth,prtvol,restart,rprimd,sppoldbl,symafm,&  
  &  symrel,tnons,unkg2,wffinp,wffout)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ceksp2
  integer,intent(in) :: debug
  integer,intent(in) :: doorth
  integer,intent(in) :: exchn2n3d
  integer,intent(in) :: fill
  integer,intent(in) :: formeig
  integer,intent(in) :: headform1
  integer,intent(in) :: iout
  integer,intent(in) :: ireadwf
  integer,intent(in) :: mband2
  integer,intent(in) :: mcg
  integer,intent(in) :: mkmem1
  integer,intent(in) :: mkmem2
  integer,intent(in) :: mpw1
  integer,intent(in) :: mpw2
  integer,intent(in) :: nkpt1
  integer,intent(in) :: nkpt2
  integer,intent(in) :: nspinor1
  integer,intent(inout) :: nspinor2
  integer,intent(in) :: nsppol1
  integer,intent(in) :: nsppol2
  integer,intent(in) :: nsym
  integer,intent(in) :: optorth
  integer,intent(in) :: prtvol
  integer,intent(in) :: restart
  integer,intent(in) :: sppoldbl
  integer,intent(in) :: unkg2
  real(dp),intent(in) :: ecut1
  real(dp),intent(in) :: ecut2
  real(dp),intent(in) :: ecut2_eff
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wffinp
  type(wffile_type),intent(inout) :: wffout
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: cg(2,mcg)
  real(dp),intent(out) :: eigen(mband2*(2*mband2)**formeig*nkpt2*nsppol2)
  real(dp),intent(in) :: gmet1(3,3)
  real(dp),intent(in) :: gmet2(3,3)
  integer,intent(in) :: indkk(nkpt2*sppoldbl,6)
  integer,intent(in) :: istwfk1(nkpt1)
  integer,intent(in) :: istwfk2(nkpt2)
  integer,intent(in) :: kg2(3,mpw2*mkmem2)
  real(dp),intent(in) :: kptns1(3,nkpt1)
  real(dp),intent(in) :: kptns2(3,nkpt2)
  integer,intent(in) :: nband1(nkpt1*nsppol1)
  integer,intent(in) :: nband2(nkpt2*nsppol2)
  integer,intent(in) :: npwarr1(nkpt1)
  integer,intent(in) :: npwarr2(nkpt2)
  real(dp),intent(out) :: occ(mband2*nkpt2*nsppol2)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
 end subroutine newkpt
end interface

interface
 subroutine newrho(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,dtn_pc,dtset,etotal,fcart,ffttomix,filfft,&  
  &  f_atm,f_fftgr,f_paw,gmet,grhf,gsqcut,initialized,&  
  &  ispmix,istep,i_rhor,i_vresid,i_vrespc,i_vtrial,kg_diel,kxc,mgfft,mgfftdiel,mixtofft,&  
  &  moved_atm_inside,mpi_enreg,nattyp,nfft,nfftmix,ngfft,ngfftmix,nkxc,npawmix,npwdiel,&  
  &  nresid,ntypat,n_fftgr,n_index,n1xccc,pawrhoij,pawtab,&  
  &  ph1d,psps,rhog,rhor,rprimd,susmat,usepaw,vtrial,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: dbl_nnsclo
  integer,intent(in) :: dielstrt
  integer,intent(in) :: initialized
  integer,intent(in) :: ispmix
  integer,intent(in) :: istep
  integer,intent(in) :: mgfft
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n_fftgr
  integer,intent(in) :: n_index
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftmix
  integer,intent(in) :: nkxc
  integer,intent(in) :: npawmix
  integer,intent(in) :: npwdiel
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: etotal
  character(len=fnlen),intent(in) :: filfft
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftmix(18)
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(inout) :: dtn_pc(3,dtset%natom)
  real(dp),intent(inout) :: f_atm(3,dtset%natom,n_fftgr)
  real(dp),intent(inout) :: f_fftgr(ispmix*nfftmix,dtset%nspden,n_fftgr*dtset%mffmem)
  real(dp),intent(inout) :: f_paw(npawmix,n_fftgr*dtset%mffmem*usepaw)
  real(dp),intent(in) :: fcart(3,dtset%natom)
  integer,intent(in) :: ffttomix(nfft*(1-nfftmix/nfft))
  real(dp),intent(inout) :: gmet(3,3)
  real(dp),intent(in) :: grhf(3,dtset%natom)
  integer,intent(inout) :: i_rhor(n_index)
  integer,intent(inout) :: i_vresid(n_index)
  integer,intent(inout) :: i_vrespc(n_index)
  integer,intent(inout) :: i_vtrial(n_index)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  integer,intent(in) :: mixtofft(nfftmix*(1-nfftmix/nfft))
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(inout) :: nresid(nfft,dtset%nspden)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(out) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(in) :: vtrial(nfft,dtset%nspden)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine newrho
end interface

interface
 subroutine newvtr(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,&  
  &  dtn_pc,dtset,efermi,etotal,fcart,ffttomix,filfft,&  
  &  f_atm,f_fftgr,f_paw,gmet,grhf,gsqcut,&  
  &  initialized,ispmix,&  
  &  istep,i_rhor,i_vresid,i_vrespc,i_vtrial,&  
  &  kg_diel,kxc,mgfft,mgfftdiel,mixtofft,&  
  &  moved_atm_inside,mpi_enreg,nattyp,nfft,nfftmix,&  
  &  nhat,nhatgr,nhatgrdim,&  
  &  ngfft,ngfftmix,nkxc,npawmix,npwdiel,&  
  &  nstep,ntypat,n_fftgr,n_index,n1xccc,optres,optxc,&  
  &  pawrhoij,pawang,pawfgrtab,&  
  &  ph1d,&  
  &  psps,rhor,rprimd,susmat,usepaw,&  
  &  vhartr,vnew_mean,vpsp,vresid,&  
  &  vtrial,vxc,xred,&  
  &  atindx1,cg,deltae,densymop_gs,&  
  &  dtfil,eeig,eew,eigen,eii,ek,enl,entropy,epaw,epawdc,irrzon,kg,&  
  &  nfftf,&  
  &  ngfftf,npwarr,n3xccc,occ,optene,&  
  &  pawfgr,pawtab,phnons,&  
  &  resid,rhog,&  
  &  usexcnhat,&  
  &  wffnow,&  
  &  ylm,nspinor,xccc3d )
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: dbl_nnsclo
  integer,intent(in) :: dielstrt
  integer,intent(in) :: initialized
  integer,intent(in) :: ispmix
  integer,intent(in) :: istep
  integer,intent(in) :: mgfft
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: n_fftgr
  integer,intent(in) :: n_index
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftf
  integer,intent(in) :: nfftmix
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: npawmix
  integer,intent(in) :: npwdiel
  integer,intent(inout) :: nspinor
  integer,intent(in) :: nstep
  integer,intent(in) :: ntypat
  integer,intent(in) :: optene
  integer,intent(in) :: optres
  integer,intent(in) :: optxc
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  real(dp),intent(in) :: deltae
  type(dens_sym_operator_type),intent(in) :: densymop_gs
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: eeig
  real(dp),intent(in) :: eew
  real(dp),intent(in) :: efermi
  real(dp),intent(in) :: eii
  real(dp),intent(out) :: ek
  real(dp),intent(out) :: enl
  real(dp),intent(in) :: entropy
  real(dp),intent(in) :: epaw
  real(dp),intent(in) :: epawdc
  real(dp),intent(in) :: etotal
  character(len=fnlen),intent(in) :: filfft
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: ngfftmix(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(inout) :: dtn_pc(3,dtset%natom)
  real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(inout) :: f_atm(3,dtset%natom,n_fftgr)
  real(dp),intent(inout) :: f_fftgr(ispmix*nfftmix,dtset%nspden,n_fftgr*dtset%mffmem)
  real(dp),intent(inout) :: f_paw(npawmix,n_fftgr*dtset%mffmem*usepaw)
  real(dp),intent(in) :: fcart(3,dtset%natom)
  integer,intent(in) :: ffttomix(nfft*(1-nfftmix/nfft))
  real(dp),intent(inout) :: gmet(3,3)
  real(dp),intent(in) :: grhf(3,dtset%natom)
  integer,intent(inout) :: i_rhor(n_index)
  integer,intent(inout) :: i_vresid(n_index)
  integer,intent(inout) :: i_vrespc(n_index)
  integer,intent(inout) :: i_vtrial(n_index)
  integer,intent(in) :: irrzon(nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  integer,intent(in) :: mixtofft(nfftmix*(1-nfftmix/nfft))
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(nfftf,dtset%nspden*psps%usepaw)
  real(dp),intent(in) :: nhatgr(nfftf,dtset%nspden,3*nhatgrdim)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawfgrtab_type),intent(in) :: pawfgrtab(dtset%natom*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: phnons(2,nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
  real(dp),intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(inout) :: rhog(2,nfftf)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(in) :: vnew_mean(dtset%nspden)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(inout) :: vresid(nfft,dtset%nspden)
  real(dp),intent(inout) :: vtrial(nfft,dtset%nspden)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine newvtr
end interface

interface
 subroutine nlinopt(nspin,omega,nkpt,wkpt,nsymcrys,symcrys,nstval,evalv,efermi,&  
  pmat,v1,v2,v3,nmesh,de,sc,brod,tol,fnam)
  use defs_basis
  implicit none
  integer, intent(in) :: nkpt
  integer, intent(in) :: nmesh
  integer, intent(in) :: nspin
  integer, intent(in) :: nstval
  integer, intent(in) :: nsymcrys
  integer, intent(in) :: v1
  integer, intent(in) :: v2
  integer, intent(in) :: v3
  real(dp), intent(in) :: brod
  real(dp), intent(in) :: de
  real(dp), intent(in) :: efermi
  character(256), intent(in) :: fnam
  real(dp), intent(in) :: omega
  real(dp), intent(in) :: sc
  real(dp), intent(in) :: tol
  real(dp), intent(in) :: evalv(nstval,nspin,nkpt)
  complex(dpc), intent(inout) :: pmat(nstval,nstval,nkpt,3,nspin)
  real(dp), intent(in) :: symcrys(3,3,nsymcrys)
  real(dp), intent(in) :: wkpt(nkpt)
 end subroutine nlinopt
end interface

interface
 subroutine nres2vres(dtset,gsqcut,izero,kxc,mpi_enreg,nfft,ngfft,nhat,&  
  &  nkxc,nresid,n3xccc,optnc,optxc,pawang,pawfgrtab,pawrhoij,pawtab,&  
  &  rhog,rhor,rprimd,usepaw,usexcnhat,vresid,xccc3d)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: izero
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: optnc
  integer,intent(in) :: optxc
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(inout) :: nhat(nfft,dtset%nspden*usepaw)
  real(dp),intent(in) :: nresid(nfft,dtset%nspden)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%ntypat*usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vresid(nfft,dtset%nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
 end subroutine nres2vres
end interface

interface
 subroutine odamix(atindx,deltae,dtset,efield_dot,elast,energies,etotal,&  
  &  gsqcut,indsym,kxc,mgfft,mpi_enreg,nattyp,nfft,ngfft,nhat,&  
  &  nkxc,ntypat,nvresid,n1xccc,n3xccc,optene,optres,paw_ij,&  
  &  paw_an,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,pel,ph1d,&  
  &  pion,psps,rhog,rhor,rprimd,strsxc,symrec,ucvol,usepaw,&  
  &  usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,xccc3d,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: ntypat
  integer,intent(in) :: optene
  integer,intent(in) :: optres
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  real(dp),intent(out) :: deltae
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(inout) :: elast
  type(energies_type),intent(inout) :: energies
  real(dp),intent(out) :: etotal
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: vxcavg
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: efield_dot(3)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(inout) :: nhat(nfft,dtset%nspden*usepaw)
  real(dp),intent(inout) :: nvresid(nfft,dtset%nspden)
  type(paw_an_type),intent(inout) :: paw_an(dtset%natom)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: pel(3)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: pion(3)
  real(dp),intent(inout) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: strsxc(6)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(inout) :: vhartr(nfft)
  real(dp),intent(in) :: vpsp(nfft)
  real(dp),intent(inout) :: vtrial(nfft,dtset%nspden)
  real(dp),intent(inout) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine odamix
end interface

interface
 subroutine prcref(atindx,dielar,dielinv,&  
  &  dielstrt,dtn_pc,dtset,etotal,fcart,ffttomix,gmet,gsqcut,&  
  &  istep,kg_diel,kxc,&  
  &  mgfft,mgfftdiel,moved_atm_inside,mpi_enreg,&  
  &  nattyp,nfft,nfftprc,ngfft,ngfftprc,nkxc,npawmix,npwdiel,ntypat,n1xccc,&  
  &  optreal,optres,pawrhoij,pawtab,ph1d,psps,rhog,rhoijrespc,rhor,rprimd,&  
  &  susmat,vhartr,vpsp,vresid,vrespc,vxc,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: dielstrt
  integer,intent(in) :: istep
  integer,intent(in) :: mgfft
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: n1xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftprc
  integer,intent(in) :: nkxc
  integer,intent(in) :: npawmix
  integer,intent(in) :: npwdiel
  integer,intent(in) :: ntypat
  integer,intent(in) :: optreal
  integer,intent(in) :: optres
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftprc(18)
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(out) :: dtn_pc(3,dtset%natom)
  real(dp),intent(in) :: fcart(3,dtset%natom)
  integer,intent(in) :: ffttomix(nfft*(1-nfftprc/nfft))
  real(dp),intent(inout) :: gmet(3,3)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(ntypat)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(out) :: rhoijrespc(npawmix)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(in) :: vresid(nfftprc*optreal,dtset%nspden)
  real(dp),intent(out) :: vrespc(nfftprc*optreal,dtset%nspden)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine prcref
end interface

interface
 subroutine prcref_PMA(atindx,dielar,dielinv,&  
  &  dielstrt,dtn_pc,dtset,fcart,ffttomix,gmet,gsqcut,&  
  &  istep,kg_diel,kxc,lavnlr,&  
  &  mgfft,mgfftdiel,moved_atm_inside,mpi_enreg,&  
  &  nattyp,nfft,nfftprc,ngfft,ngfftprc,nkxc,npawmix,npwdiel,ntypat,n1xccc,&  
  &  optreal,optres,pawrhoij,ph1d,psps,rhog,rhoijrespc,rhor,rprimd,&  
  &  susmat,vhartr,vpsp,vresid,vrespc,vxc,xred,&  
  &  deltae,efermi,etotal,nfftf,nhat,nhatgr,nhatgrdim,optene,optxc,pawang,pawfgrtab,&  
  &  pawtab,use_lavnlr,usexcnhat,vtrial)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: dielstrt
  integer,intent(in) :: istep
  integer,intent(in) :: mgfft
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: n1xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftf
  integer,intent(in) :: nfftprc
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: npawmix
  integer,intent(in) :: npwdiel
  integer,intent(in) :: ntypat
  integer,intent(in) :: optene
  integer,intent(in) :: optreal
  integer,intent(in) :: optres
  integer,intent(in) :: optxc
  integer,intent(in) :: use_lavnlr
  integer,intent(in) :: usexcnhat
  real(dp),intent(in) :: deltae
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: efermi
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftprc(18)
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(out) :: dtn_pc(3,dtset%natom)
  real(dp),intent(in) :: fcart(3,dtset%natom)
  integer,intent(in) :: ffttomix(nfft*(1-nfftprc/nfft))
  real(dp),intent(inout) :: gmet(3,3)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: lavnlr(dtset%nfft,dtset%nspden*use_lavnlr)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(nfftf,dtset%nspden*psps%usepaw)
  real(dp),intent(in) :: nhatgr(nfftf,dtset%nspden,3*nhatgrdim)
  type(pawfgrtab_type),intent(in) :: pawfgrtab(dtset%natom*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(out) :: rhoijrespc(npawmix)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(in) :: vresid(nfftprc*optreal,dtset%nspden)
  real(dp),intent(out) :: vrespc(nfftprc*optreal,dtset%nspden)
  real(dp),intent(in) :: vtrial(dtset%nfft,dtset%nspden)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine prcref_PMA
end interface

interface
 subroutine prteigrs(eigen,enunit,fermie,filnam,iout,iscf,kptns,kptopt,mband,nband,&  
  &  nkpt,nnsclo_now,nsppol,occ,occopt,option,prteig,prtvol,resid,tolwfr,vxcavg,wtk)
  use defs_basis
  implicit none
  integer,intent(in) :: enunit
  integer,intent(in) :: iout
  integer,intent(in) :: iscf
  integer,intent(in) :: kptopt
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nnsclo_now
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: option
  integer,intent(in) :: prteig
  integer,intent(in) :: prtvol
  real(dp),intent(in) :: fermie
  character(len=fnlen),intent(in) :: filnam
  real(dp),intent(in) :: tolwfr
  real(dp),intent(in) :: vxcavg
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: resid(mband*nkpt*nsppol)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine prteigrs
end interface

interface
 subroutine prtene(dtset,energies,iout,usepaw)
  use defs_datatypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  type(energies_type),intent(in) :: energies
 end subroutine prtene
end interface

interface
 subroutine prtrhomxmn(iout,mpi_enreg,nfft,ngfft,nspden,option,rhor)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: rhor(nfft,nspden)
 end subroutine prtrhomxmn
end interface

interface
 subroutine prtxf(fred,iatfix,iout,iwfrc,natom,rprimd,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: iwfrc
  integer,intent(in) :: natom
  real(dp),intent(in) :: fred(3,natom)
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine prtxf
end interface

interface
 subroutine prtxvf(fcart,iatfix,iout,natom,prtvel,vel,xcart)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  integer,intent(in) :: prtvel
  real(dp),intent(in) :: fcart(3,natom)
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: vel(3,natom)
  real(dp),intent(in) :: xcart(3,natom)
 end subroutine prtxvf
end interface

interface
 subroutine rhophi(cx,phi,rho)
  use defs_basis
  implicit none
  real(dp),intent(out) :: phi
  real(dp),intent(out) :: rho
  real(dp),intent(in) :: cx(2)
 end subroutine rhophi
end interface

interface
 subroutine rhotov(dtset,energies,gsqcut,kxc,mpi_enreg,nfft,ngfft,&  
  &  nhat,nhatgr,nhatgrdim,nkxc,vresidnew,n3xccc,optene,optres,optxc,&  
  &  pawang,pawfgrtab,pawtab,rhog,rhor,rprimd,strsxc,ucvol,usepaw,usexcnhat,&  
  &  vhartr,vnew_mean,vpsp,vres_mean,vres2,vtrial,vxcavg,vxc,xccc3d)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: optene
  integer,intent(in) :: optres
  integer,intent(in) :: optxc
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  type(dataset_type),intent(in) :: dtset
  type(energies_type),intent(inout) :: energies
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: vres2
  real(dp),intent(out) :: vxcavg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: nhat(nfft,dtset%nspden*usepaw)
  real(dp),intent(in) :: nhatgr(nfft,dtset%nspden,3*nhatgrdim)
  type(pawfgrtab_type),intent(in) :: pawfgrtab(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: strsxc(6)
  real(dp),intent(inout) :: vhartr(nfft)
  real(dp),intent(out) :: vnew_mean(dtset%nspden)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(out) :: vres_mean(dtset%nspden)
  real(dp),intent(out) :: vresidnew(nfft,dtset%nspden)
  real(dp),intent(inout) :: vtrial(nfft,dtset%nspden)
  real(dp),intent(inout) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
 end subroutine rhotov
end interface

interface
 subroutine scfcge(cplex,dbl_nnsclo,dtn_pc,etotal,f_atm,&  
  &  f_fftgr,initialized,iscf,isecur,istep,&  
  &  i_rhor,i_vresid,i_vrespc,moved_atm_inside,mpi_enreg,&  
  &  natom,nfft,nfftot,nspden,n_fftgr,n_index,response,rhor,ucvol,vtrial,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(out) :: dbl_nnsclo
  integer,intent(in) :: initialized
  integer,intent(in) :: iscf
  integer,intent(in) :: isecur
  integer,intent(in) :: istep
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: n_fftgr
  integer,intent(in) :: n_index
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: nspden
  integer,intent(in) :: response
  real(dp),intent(in) :: etotal
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: dtn_pc(3,natom)
  real(dp),intent(inout) :: f_atm(3,natom,n_fftgr)
  real(dp),intent(inout) :: f_fftgr(cplex*nfft,nspden,n_fftgr)
  integer,intent(inout) :: i_rhor(n_index)
  integer,intent(inout) :: i_vresid(n_index)
  integer,intent(inout) :: i_vrespc(n_index)
  real(dp),intent(in) :: rhor(cplex*nfft,nspden)
  real(dp),intent(inout) :: vtrial(cplex*nfft,nspden)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine scfcge
end interface

interface
 subroutine scfeig(f_fftgr,istep,i_vresid1,i_vrespc1,nfft,nspden,n_fftgr,vtrial)
  use defs_basis
  implicit none
  integer,intent(in) :: i_vresid1
  integer,intent(in) :: i_vrespc1
  integer,intent(in) :: istep
  integer,intent(in) :: n_fftgr
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  real(dp),intent(inout) :: f_fftgr(nfft,nspden,n_fftgr)
  real(dp),intent(inout) :: vtrial(nfft,nspden)
 end subroutine scfeig
end interface

interface
 subroutine scfopt(cplex,dtn_pc,f_fftgr,f_paw,iscf,istep,i_vrespc,i_vtrial,move_atm,mpi_enreg,&  
  &  natom,nfft,npawmix,nspden,n_fftgr,n_index,pawoptmix,usepaw,vpaw,vtrial,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: iscf
  integer,intent(in) :: istep
  integer,intent(in) :: move_atm
  integer,intent(in) :: n_fftgr
  integer,intent(in) :: n_index
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: npawmix
  integer,intent(in) :: nspden
  integer,intent(in) :: pawoptmix
  integer,intent(in) :: usepaw
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: dtn_pc(3,natom)
  real(dp),intent(inout) :: f_fftgr(cplex*nfft,nspden,n_fftgr)
  real(dp),intent(inout) :: f_paw(npawmix,n_fftgr*usepaw)
  integer,intent(inout) :: i_vrespc(n_index)
  integer,intent(inout) :: i_vtrial(n_index)
  real(dp),intent(inout) :: vpaw(npawmix*usepaw)
  real(dp),intent(inout) :: vtrial(cplex*nfft,nspden)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine scfopt
end interface

interface
 subroutine scprqt(choice,cpus,deltae,diffor,dtset,&  
  &  eigen,etotal,favg,fcart,fermie,filapp,filnam1,initGS,&  
  &  iscf,istep,kpt,maxfor,moved_atm_inside,mpi_enreg,&  
  &  nband,nkpt,nstep,occ,optres,&  
  &  prtfor,quit,res2,resid,residm,response,tollist,usepaw,&  
  &  vxcavg,wtk,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: initGS
  integer,intent(in) :: iscf
  integer,intent(in) :: istep
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: nkpt
  integer,intent(in) :: nstep
  integer,intent(in) :: optres
  integer,intent(in) :: prtfor
  integer,intent(out) :: quit
  integer,intent(in) :: response
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: cpus
  real(dp),intent(in) :: deltae
  real(dp),intent(in) :: diffor
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: fermie
  character(len=fnlen),intent(in) :: filapp
  character(len=fnlen),intent(in) :: filnam1
  real(dp),intent(in) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: res2
  real(dp),intent(in) :: residm
  real(dp),intent(in) :: vxcavg
  real(dp),intent(in) :: eigen(dtset%mband*nkpt*dtset%nsppol)
  real(dp),intent(in) :: favg(3)
  real(dp),intent(in) :: fcart(3,dtset%natom)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: nband(nkpt*dtset%nsppol)
  real(dp),intent(in) :: occ(dtset%mband*nkpt*dtset%nsppol)
  real(dp),intent(in) :: resid(dtset%mband*nkpt*dtset%nsppol)
  real(dp),intent(in) :: tollist(12)
  real(dp),intent(in) :: wtk(nkpt)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine scprqt
end interface

interface
 subroutine setup1(acell,amass,amu,bantot,&  
  &  ecut_eff,ecutc_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,iboxcut,intxc,ionmov,&  
  &  natom,nband,ngfft,ngfftc,nkpt,nqpt,nsppol,nsym,ntypat,&  
  &  qptn,response,rmet,rprim,rprimd,typat,ucvol,usepaw)
  use defs_basis
  implicit none
  integer,intent(out) :: bantot
  integer,intent(in) :: iboxcut
  integer,intent(in) :: intxc
  integer,intent(in) :: ionmov
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nqpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: response
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: ecut_eff
  real(dp),intent(in) :: ecutc_eff
  real(dp),intent(out) :: gsqcut_eff
  real(dp),intent(out) :: gsqcutc_eff
  real(dp),intent(out) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftc(18)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(out) :: amass(natom)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(out) :: gmet(3,3)
  real(dp),intent(out) :: gprimd(3,3)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(in) :: qptn(3)
  real(dp),intent(out) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(out) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
 end subroutine setup1
end interface

interface
 subroutine setup2(dtset,epulay,iscf,&  
  &  npwtot,start,ucvol,wfs,xred)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: iscf
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: epulay
  real(dp),intent(in) :: ucvol
  type(wvl_wf_type),intent(in) :: wfs
  integer,intent(in) :: npwtot(dtset%nkpt)
  real(dp),intent(out) :: start(3,dtset%natom)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine setup2
end interface

interface
 subroutine setup_positron(dtset, fildensin,filstat,filvhain,mpi_enreg,n1,n2,n3,&  
  &  hdr,iexit,level,nfft,n1xccc,ntypat,&  
  &  xcccrc,xccc1d,etotal,rhocore,rhore,rhorp,rprimd,ucvol,&  
  &  vhae,vhap,xred,ngfft)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: level
  integer,intent(in) :: n1
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(inout) :: etotal
  character(len=fnlen),intent(in) :: fildensin
  character(len=fnlen),intent(in) :: filstat
  character(len=fnlen),intent(in) :: filvhain
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: rhocore(nfft)
  real(dp),intent(inout) :: rhore(nfft,2)
  real(dp),intent(inout) :: rhorp(nfft,2)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: vhae(nfft,2)
  real(dp),intent(inout) :: vhap(nfft,2)
  real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine setup_positron
end interface

interface
 subroutine setvtr(atindx1,dtset,energies,gmet,gprimd,&  
  &  grewtn,gsqcut,initialized,&  
  &  istep,kxc,mgfft,moved_atm_inside,moved_rhor,mpi_enreg,&  
  &  nattyp,nfft,ngfft,nhat,nhatgr,nhatgrdim,nkxc,ntypat,n1xccc,n3xccc,&  
  &  optene,pawtab,ph1d,psps,rhog,rhor,rmet,rprimd,strsxc,&  
  &  ucvol,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,&  
  &  xccc3d,xred,xred_old)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(inout) :: initialized
  integer,intent(in) :: istep
  integer,intent(in) :: mgfft
  integer,intent(inout) :: moved_atm_inside
  integer,intent(inout) :: moved_rhor
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: ntypat
  integer,intent(in) :: optene
  integer,intent(in) :: usexcnhat
  type(dataset_type),intent(in) :: dtset
  type(energies_type),intent(inout) :: energies
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: vxcavg
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: grewtn(3,dtset%natom)
  real(dp),intent(out) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(nfft,dtset%nspden*psps%usepaw)
  real(dp),intent(in) :: nhatgr(nfft,dtset%nspden,3*nhatgrdim)
  type(pawtab_type),intent(in) :: pawtab(ntypat*dtset%usepaw)
  real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: strsxc(6)
  real(dp),intent(inout) :: vhartr(nfft)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(inout) :: vtrial(nfft,dtset%nspden)
  real(dp),intent(inout) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(in) :: xred_old(3,dtset%natom)
 end subroutine setvtr
end interface

interface
 subroutine smatrix(cg,cgq,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,job,maxbd,mband,&  
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
 end subroutine smatrix
end interface

interface
 subroutine spin_current(atindx,atindx1,cg,dtfil,dtset,eigen,gmet,gprimd,hdr,kg,mpi_enreg,&  
  &  nattyp,nfftf,ph1d,psps,rhog,rhor,rmet,symrec,ucvol,wffnow,ylm,ylmgr)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nfftf
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(in) :: wffnow
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: cg(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: nattyp(dtset%ntypat)
  real(dp),intent(inout) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(in) :: rhog(2,nfftf)
  real(dp),intent(in) :: rhor(nfftf,dtset%nspden)
  real(dp),intent(in) :: rmet(3,3)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine spin_current
end interface

interface
 subroutine stress(atindx1,berryopt,dedlnn,eei,efield,ehart,eii,gsqcut,kinstr,&  
  &  mgfft,mpi_enreg,mqgrid,n1xccc,n3xccc,natom,nattyp,&  
  &  nfft,ngfft,nlstr,nspden,nsym,ntypat,paral_kgb,pawtab,pel,pion,ph1d,&  
  &  prtvol,qgrid,rhog,rprimd,strten,strsxc,symrec,typat,usepaw,vlspl,&  
  &  vxc,xccc1d,xccc3d,xcccrc,xred,zion)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: dedlnn
  real(dp),intent(in) :: eei
  real(dp),intent(in) :: ehart
  real(dp),intent(in) :: eii
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: efield(3)
  real(dp),intent(in) :: kinstr(6)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nlstr(6)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: pel(3)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: pion(3)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: strsxc(6)
  real(dp),intent(out) :: strten(6)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vlspl(mqgrid,2,ntypat)
  real(dp),intent(in) :: vxc(nfft,nspden)
  real(dp),intent(in) :: xccc1d(n1xccc*(1-usepaw),6,ntypat)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine stress
end interface

interface
 subroutine strhar(ehart,gprimd,gsqcut,harstr,mpi_enreg,nfft,ngfft,rhog,ucvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nfft
  real(dp),intent(in) :: ehart
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: harstr(6)
  real(dp),intent(in) :: rhog(2,nfft)
 end subroutine strhar
end interface

interface
 subroutine sygrad(fred,natom,dedt,nsym,symrec,indsym)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  real(dp),intent(in) :: dedt(3,natom)
  real(dp),intent(out) :: fred(3,natom)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine sygrad
end interface

interface
 subroutine symrhg(cplex,densymop,irrzon,mpi_enreg,nfft,nfftot,ngfft,nspden,nsppol,nsym,paral_kgb,&  
  &  phnons,rhog,rhor,symafm)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: paral_kgb
  type(dens_sym_operator_type),intent(in) :: densymop
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: irrzon(nfftot**(1-1/nsym),2,nspden/nsppol)
  real(dp),intent(in) :: phnons(2,nfftot**(1-1/nsym),nspden/nsppol)
  real(dp),intent(out) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(cplex*nfft,nspden)
  integer,intent(in) :: symafm(nsym)
 end subroutine symrhg
end interface

interface
 subroutine testsusmat(compute,dielop,dielstrt,dtset,istep)
  use defs_datatypes
  implicit none
  integer,intent(in) :: dielop
  integer,intent(in) :: dielstrt
  integer,intent(in) :: istep
  logical,intent(out) :: compute
  type(dataset_type),intent(in) :: dtset
 end subroutine testsusmat
end interface

interface
 subroutine uderiv(bdberry,cg,gprimd,hdr,istwfk,kberry,kg,kpt_,kptopt,kptrlatt,&  
  &  mband,mgfft,mkmem,mpi_enreg,mpw,natom,nband,nberry,npwarr,nspinor,nsppol,nkpt_,&  
  &  rprimd,unddk,unkg,wffnow,filnam)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: kptopt
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nberry
  integer,intent(in) :: nkpt_
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: unddk
  integer,intent(in) :: unkg
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: bdberry(4)
  integer,intent(in) :: kberry(3,20)
  integer,intent(in) :: kptrlatt(3,3)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  character(len=fnlen),intent(in) :: filnam(5)
  real(dp),intent(in) :: gprimd(1:3,1:3)
  integer,intent(in) :: istwfk(nkpt_)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kpt_(3,nkpt_)
  integer,intent(in) :: nband(nkpt_*nsppol)
  integer,intent(in) :: npwarr(nkpt_)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine uderiv
end interface

interface
 subroutine vtorhotf(densymop_gs,dtfil,dtset,ek,enl,entropy,fermie,grnl,&  
  &  irrzon,mpi_enreg,natom,nfft,nspden,nsppol,nsym,phnons,rhog,rhor,ucvol,vtrial)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  type(dens_sym_operator_type),intent(in) :: densymop_gs
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: ek
  real(dp),intent(out) :: enl
  real(dp),intent(out) :: entropy
  real(dp),intent(out) :: fermie
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: grnl(3*natom)
  integer,intent(in) :: irrzon((dtset%ngfft(1)*dtset%ngfft(1)*dtset%ngfft(1))**(1-1/nsym), &
  &         2,nspden/nsppol)
  real(dp),intent(in) :: phnons(2,(dtset%ngfft(1)*dtset%ngfft(1)*dtset%ngfft(1))**(1-1/nsym), &
  &         nspden/nsppol)
  real(dp),intent(inout) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,nspden)
  real(dp),intent(in) :: vtrial(nfft,nspden)
 end subroutine vtorhotf
end interface

interface
 function zfermim12(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermim12
 end function zfermim12
end interface

interface
 function zfermi12(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi12
 end function zfermi12
end interface

interface
 function zfermi1(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi1
 end function zfermi1
end interface

interface
 function zfermi32(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi32
 end function zfermi32
end interface

interface
 function zfermi2(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi2
 end function zfermi2
end interface

interface
 function zfermi52(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi52
 end function zfermi52
end interface

interface
 function zfermi3(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi3
 end function zfermi3
end interface

interface
 function ifermim12(ff)
  use defs_basis
  implicit none
  real(dp), intent(in) :: ff
  real(dp) :: ifermim12
 end function ifermim12
end interface

interface
 function ifermi12(ff)
  use defs_basis
  implicit none
  real(dp), intent(in) :: ff
  real(dp) :: ifermi12
 end function ifermi12
end interface

interface
 function ifermi32(ff)
  use defs_basis
  implicit none
  real(dp), intent(in) :: ff
  real(dp) :: ifermi32
 end function ifermi32
end interface

interface
 function ifermi52(ff)
  use defs_basis
  implicit none
  real(dp), intent(in) :: ff
  real(dp) :: ifermi52
 end function ifermi52
end interface

interface
 function fp12a1 (x)
  use defs_basis
  implicit none
  real(dp) :: fp12a1
  real(dp) :: x
 end function fp12a1
end interface

interface
 function fp32a1 (x)
  use defs_basis
  implicit none
  real(dp) :: fp32a1
  real(dp) :: x
 end function fp32a1
end interface

interface
 function xp12a1 (y)
  use defs_basis
  implicit none
  real(dp) :: xp12a1
  real(dp) :: y
 end function xp12a1
end interface

interface
 function fm12a1 (x)
  use defs_basis
  implicit none
  real(dp) :: fm12a1
  real(dp) :: x
 end function fm12a1
end interface

interface
 subroutine fm12a1t (cktf,rtnewt,tsmear,vtrial,rhor_middx,rhor_mid,&  
  &  nfft)
  use defs_basis
  implicit none
  integer :: nfft
  real(dp) :: cktf
  real(dp) :: rtnewt
  real(dp) :: tsmear
  real(dp) :: rhor_mid(nfft)
  real(dp) :: rhor_middx(nfft)
  real(dp) :: vtrial(nfft)
 end subroutine fm12a1t
end interface

interface
 subroutine waveformat(cg,cg_disk,cg_index,cg_new,dk,ii,ikpt,&  
  &  ikpt_,isgn,isppol,jj,jkpt,jkpt_,kg_kpt,kpt,kg_jl,maxband,mband,&  
  &  minband,mkmem,mpw,nkpt,nkpt_,npwarr,nsppol,nspinor,shift_g_2,tr)
  use defs_basis
  implicit none
  integer,intent(in) :: ii
  integer,intent(in) :: ikpt
  integer,intent(in) :: ikpt_
  integer,intent(in) :: isgn
  integer,intent(in) :: isppol
  integer,intent(in) :: jj
  integer,intent(in) :: jkpt
  integer,intent(in) :: jkpt_
  integer,intent(in) :: maxband
  integer,intent(in) :: mband
  integer,intent(in) :: minband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkpt_
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cg_disk(2,mpw*nspinor*mband,2)
  integer,intent(in) :: cg_index(mband,nkpt_,nsppol)
  real(dp),intent(out) :: cg_new(2,mpw,maxband)
  real(dp),intent(in) :: dk(3)
  integer,intent(in) :: kg_jl(3,mpw,2)
  integer,intent(in) :: kg_kpt(3,mpw*nspinor,nkpt_)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: npwarr(nkpt_)
  logical,intent(in) :: shift_g_2(nkpt,nkpt)
  real(dp),intent(in) :: tr(2)
 end subroutine waveformat
end interface

interface
 subroutine wvl_mkrho(dtset, mpi_enreg, occ, rhor, wfs)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(in) :: mpi_enreg
  type(wvl_wf_type),intent(in) :: wfs
  real(dp),intent(in) :: occ(dtset%mband)
  real(dp),intent(inout) :: rhor(dtset%nfft,dtset%nspden)
 end subroutine wvl_mkrho
end interface

interface
 subroutine wvl_newvtr(dtset, mpi_enreg, nele, offset, vhartr, vpsp, vtrial, vxc)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: nele
  integer,intent(out) :: offset
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(in) :: vhartr(dtset%nfft)
  real(dp),intent(in) :: vpsp(dtset%nfft)
  real(dp),intent(out) :: vtrial(dtset%nfft*dtset%nspden)
  real(dp),intent(in) :: vxc(dtset%nfft*dtset%nspden)
 end subroutine wvl_newvtr
end interface

end module interfaces_15common
!!***
