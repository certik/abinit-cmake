!!****m* ABINIT/interfaces_15gw
!! NAME
!! interfaces_15gw
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/15gw
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

module interfaces_15gw

 implicit none

interface
 subroutine assemblychi0_sym(ik_bz,nspinor,Ep,Ltg_q,green_w,npwepG0,rhotwg,Gsph_epsG0,chi0)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(in) :: npwepG0
  integer,intent(in) :: nspinor
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(little_group),intent(in) :: Ltg_q
  complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
  complex(dpc),intent(in) :: green_w(Ep%nomega)
  complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
 end subroutine assemblychi0_sym
end interface

interface
 subroutine mkrhotwg_sigma(ii,nspinor,npw,rhotwg,rhotwg_I)
  use defs_basis
  implicit none
  integer,intent(in) :: ii
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  complex(gwpc),intent(in) :: rhotwg(npw*nspinor**2)
  complex(gwpc),intent(inout) :: rhotwg_I(npw)
 end subroutine mkrhotwg_sigma
end interface

interface
 subroutine assemblychi0q0_sym(qpoint,ik_bz,isym_kbz,itim_kbz,gwcomp,nspinor,npwepG0,Ep,Cryst,Ltg_q,Gsph_epsG0,&  
  &  chi0,rhotwx,rhotwg,green_w,green_enhigh_w,deltaf_b1b2)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: gwcomp
  integer,intent(in) :: ik_bz
  integer,intent(in) :: isym_kbz
  integer,intent(in) :: itim_kbz
  integer,intent(in) :: npwepG0
  integer,intent(in) :: nspinor
  type(crystal_structure),intent(in) :: Cryst
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(little_group),intent(in) :: Ltg_q
  real(dp),intent(in) :: deltaf_b1b2
  complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
  complex(dpc),intent(in) :: green_enhigh_w(Ep%nomega)
  complex(dpc),intent(in) :: green_w(Ep%nomega)
  real(dp),intent(in) :: qpoint(3)
  complex(gwpc),intent(inout) :: rhotwg(npwepG0*nspinor**2)
  complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)
 end subroutine assemblychi0q0_sym
end interface

interface
 function q0limit(ii,qpoint,nspinor,rhotwx,b1,b2,b3)
  use defs_basis
  implicit none
  integer,intent(in) :: ii
  integer,intent(in) :: nspinor
  complex(gwpc) :: q0limit
  real(dp),intent(in) :: b1(3)
  real(dp),intent(in) :: b2(3)
  real(dp),intent(in) :: b3(3)
  real(dp),intent(in) :: qpoint(3)
  complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)
 end function q0limit
end interface

interface
 subroutine assemblychi0sf(ik_bz,nspinor,symchi,Ltg_q,npwepG0,npwe,rhotwg,Gsph_epsG0,&  
  &  factocc,my_wl,iomegal,wl,my_wr,iomegar,wr,nomegasf,chi0sf)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(in) :: iomegal
  integer,intent(in) :: iomegar
  integer,intent(in) :: my_wl
  integer,intent(in) :: my_wr
  integer,intent(in) :: nomegasf
  integer,intent(in) :: npwe
  integer,intent(in) :: npwepG0
  integer,intent(in) :: nspinor
  integer,intent(in) :: symchi
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(little_group),intent(in) :: Ltg_q
  real(dp),intent(in) :: factocc
  real(dp),intent(in) :: wl
  real(dp),intent(in) :: wr
  complex(gwpc),intent(inout) :: chi0sf(npwe,npwe,my_wl:my_wr)
  complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
 end subroutine assemblychi0sf
end interface

interface
 subroutine assemblychi0sfq0(qpoint,ikbz,isym_kbz,itim_kbz,nspinor,symchi,npwepG0,npwe,Cryst,Ltg_q,Gsph_epsG0,&  
  &  factocc,my_wl,iomegal,wl,my_wr,iomegar,wr,rhotwx,rhotwg,nomegasf,chi0sf)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ikbz
  integer,intent(in) :: iomegal
  integer,intent(in) :: iomegar
  integer,intent(in) :: isym_kbz
  integer,intent(in) :: itim_kbz
  integer,intent(in) :: my_wl
  integer,intent(in) :: my_wr
  integer,intent(in) :: nomegasf
  integer,intent(in) :: npwe
  integer,intent(in) :: npwepG0
  integer,intent(in) :: nspinor
  integer,intent(in) :: symchi
  type(crystal_structure),intent(in) :: Cryst
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(little_group),intent(in) :: Ltg_q
  real(dp),intent(in) :: factocc
  real(dp),intent(in) :: wl
  real(dp),intent(in) :: wr
  complex(gwpc),intent(inout) :: chi0sf(npwe,npwe,my_wl:my_wr)
  real(dp),intent(in) :: qpoint(3)
  complex(gwpc),intent(inout) :: rhotwg(npwepG0*nspinor**2)
  complex(gwpc),intent(in) :: rhotwx(3)
 end subroutine assemblychi0sfq0
end interface

interface
 subroutine init_Bands_Symmetries(BSym,only_trace,nspinor,nsppol,nbnds,kpt,nsym,symrec,tnons,EDIFF_TOL,ene_k,ierr)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: nbnds
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  type(bands_symmetries),intent(inout) :: BSym
  real(dp),intent(in) :: EDIFF_TOL
  logical,intent(in) :: only_trace
  real(dp),intent(in) :: ene_k(nbnds,nsppol)
  real(dp),intent(in) :: kpt(3)
  integer,intent(in) :: symrec(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
 end subroutine init_Bands_Symmetries
end interface

interface
 subroutine print_Bands_Symmetries(Bsym,unitno,prtvol,mode_paral)
  use defs_datatypes
  implicit none
  integer,intent(in),optional :: prtvol
  integer,intent(in),optional :: unitno
  type(bands_symmetries),intent(in) :: BSym
  character(len=4),intent(in),optional :: mode_paral
 end subroutine print_Bands_Symmetries
end interface

interface
 subroutine destroy_Bands_Symmetries(Bsym)
  use defs_datatypes
  implicit none
  type(bands_symmetries),intent(inout) :: BSym
 end subroutine destroy_Bands_Symmetries
end interface

interface
 subroutine destroy_Degenerate_Bands(Cplx)
  use defs_datatypes
  implicit none
  type(degenerate_bands),intent(inout) :: Cplx(:,:)
 end subroutine destroy_Degenerate_Bands
end interface

interface
 subroutine nullify_Bands_Symmetries(BSym)
  use defs_datatypes
  implicit none
  type(bands_symmetries),intent(inout) :: Bsym
 end subroutine nullify_Bands_Symmetries
end interface

interface
 subroutine nullify_Degenerate_Bands(Cplx)
  use defs_datatypes
  implicit none
  type(degenerate_bands),intent(inout) :: Cplx(:,:)
 end subroutine nullify_Degenerate_Bands
end interface

interface
 subroutine get_class(nsym,symrec,nclass,nelements,elements_idx)
  implicit none
  integer,intent(out) :: nclass
  integer,intent(in) :: nsym
  integer,intent(out) :: elements_idx(nsym,nsym)
  integer,intent(out) :: nelements(nsym)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine get_class
end interface

interface
 subroutine bz1(k,g,gmet)
  use defs_basis
  implicit none
  integer,intent(out) :: g(3)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(inout) :: k(3)
 end subroutine bz1
end interface

interface
 subroutine calc_coh(paral_kgb,nspinor,nsig_ab,nfftot,ngfft_gw,tim_fourdp,MPI_enreg,ktabi_k,ktabr_k,spinrot_k,&  
  &  wfr_jb,wfr_kb,npwx,npwc,gvec,epsm1q_o,vc_sqrt,i_sz,iqibz,same_band,sigcohme)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iqibz
  integer,intent(in) :: ktabi_k
  integer,intent(in) :: nfftot
  integer,intent(in) :: npwc
  integer,intent(in) :: npwx
  integer,intent(in) :: nsig_ab
  integer,intent(in) :: nspinor
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourdp
  type(mpi_type),intent(inout) :: MPI_enreg
  real(dp),intent(in) :: i_sz
  logical,intent(in) :: same_band
  integer,intent(in) :: ngfft_gw(18)
  complex(gwpc),intent(in) :: epsm1q_o(npwc,npwc)
  integer,intent(in) :: gvec(3,npwx)
  integer,intent(in) :: ktabr_k(nfftot)
  complex(gwpc),intent(out) :: sigcohme(nsig_ab)
  real(dp),intent(in) :: spinrot_k(4)
  complex(gwpc),intent(in) :: vc_sqrt(npwx)
  complex(gwpc),intent(in) :: wfr_jb(nfftot*nspinor)
  complex(gwpc),intent(in) :: wfr_kb(nfftot*nspinor)
 end subroutine calc_coh
end interface

interface
 subroutine calc_wfwfg(MPI_enreg,paral_kgb,tim_fourdp,ktabr_k,ktabi_k,nfftot,ngfft_gw,wfr_jb,wfr_kb,wfg2_jk)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ktabi_k
  integer,intent(in) :: nfftot
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourdp
  type(mpi_type),intent(inout) :: MPI_enreg
  integer,intent(in) :: ngfft_gw(18)
  integer,intent(in) :: ktabr_k(nfftot)
  complex(gwpc),intent(out) :: wfg2_jk(nfftot)
  complex(gwpc),intent(in) :: wfr_jb(nfftot)
  complex(gwpc),intent(in) :: wfr_kb(nfftot)
 end subroutine calc_wfwfg
end interface

interface
 subroutine calc_density(Wfs,Cryst,irottbf,nbnds,ngfftf,nfftf,igfftf,&  
  &  occ,rhor,Kmesh,MPI_enreg,tim_fourdp,use_MPI)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nbnds
  integer,intent(in) :: nfftf
  integer,intent(in) :: tim_fourdp
  type(crystal_structure),intent(in) :: Cryst
  type(bz_mesh_type),intent(in) :: Kmesh
  type(mpi_type),intent(inout) :: MPI_enreg
  type(wavefunctions_information),intent(inout) :: Wfs
  logical,intent(in) :: use_MPI
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: igfftf(Wfs%npwwfn)
  integer,intent(in) :: irottbf(nfftf,Cryst%nsym)
  real(dp),intent(in) :: occ(Kmesh%nibz,nbnds,Wfs%nsppol)
  real(dp),intent(out) :: rhor(nfftf,Wfs%nspden)
 end subroutine calc_density
end interface

interface
 subroutine test_charge(nfftf,nspden,rhor,ucvol,nhat,&  
  &  usepaw,usexcnhat,usefinegrid,compch_sph,compch_fft,omegaplasma)
  use defs_basis
  implicit none
  integer,intent(in) :: nfftf
  integer,intent(in) :: nspden
  integer,intent(in) :: usefinegrid
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  real(dp),intent(in) :: compch_fft
  real(dp),intent(in) :: compch_sph
  real(dp),intent(out) :: omegaplasma
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: nhat(nfftf,nspden*usepaw)
  real(dp),intent(in) :: rhor(nfftf,nspden)
 end subroutine test_charge
end interface

interface
 subroutine calc_ffm(epsm1,nq,npw,nomega,omega,gprimd,qq,ppmodel,gvec,nfmidm)
  use defs_basis
  implicit none
  integer,intent(in) :: nfmidm
  integer,intent(in) :: nomega
  integer,intent(in) :: npw
  integer,intent(in) :: nq
  integer,intent(in) :: ppmodel
  complex(gwpc),intent(in) :: epsm1(npw,npw,nomega,nq)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,npw)
  complex(gwpc),intent(in) :: omega(nomega)
  real(dp),intent(in) :: qq(3,nq)
 end subroutine calc_ffm
end interface

interface
 subroutine calc_rpa_functional(iq,Ep,Pvc,Qmesh,Dtfil,gmet,kxc,MPI_enreg,chi0)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iq
  type(datafiles_type),intent(in) :: Dtfil
  type(epsilonm1_parameters),intent(in) :: Ep
  type(mpi_type),intent(inout) :: MPI_enreg
  type(coulombian_type),intent(in) :: Pvc
  type(bz_mesh_type),intent(in) :: Qmesh
  complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
  real(dp),intent(in) :: gmet(3,3)
  complex(gwpc),intent(in) :: kxc(:,:)
 end subroutine calc_rpa_functional
end interface

interface
 subroutine calc_sig_noppm(npwc,npwx,nspinor,nomega,nomegae,nomegaer,nomegaei,rhotwgp,&  
  &  omega,epsm1q,omegame0i,theta_mu_minus_e0i,ket)
  use defs_basis
  implicit none
  integer,intent(in) :: nomega
  integer,intent(in) :: nomegae
  integer,intent(in) :: nomegaei
  integer,intent(in) :: nomegaer
  integer,intent(in) :: npwc
  integer,intent(in) :: npwx
  integer,intent(in) :: nspinor
  real(dp),intent(in) :: theta_mu_minus_e0i
  complex(gwpc),intent(in) :: epsm1q(npwc,npwc,nomegae)
  complex(gwpc),intent(inout) :: ket(npwc*nspinor,nomega)
  complex(gwpc),intent(in) :: omega(nomegae)
  real(dp),intent(in) :: omegame0i(nomega)
  complex(gwpc),intent(in) :: rhotwgp(npwx*nspinor)
 end subroutine calc_sig_noppm
end interface

interface
 subroutine calc_sig_ppm(PPm,nspinor,npwc,nomega,rhotwgp,botsq,otq,&  
  &  omegame0i,zcut,theta_mu_minus_e0i,eig,npwx,ket,sigcme)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc
  integer,intent(in) :: npwx
  integer,intent(in) :: nspinor
  type(ppmodel_type),intent(in) :: PPm
  real(dp),intent(in) :: theta_mu_minus_e0i
  real(dp),intent(in) :: zcut
  complex(gwpc),intent(in) :: botsq(npwc,PPm%dm2_botsq)
  complex(gwpc),intent(in) :: eig(PPm%dm_eig,PPm%dm_eig)
  complex(gwpc),intent(inout) :: ket(npwc*nspinor,nomega)
  real(dp),intent(in) :: omegame0i(nomega)
  real(dp),intent(in) :: otq(npwc,PPm%dm2_otq)
  complex(gwpc),intent(in) :: rhotwgp(npwx*nspinor)
  complex(gwpc),intent(out) :: sigcme(nomega)
 end subroutine calc_sig_ppm
end interface

interface
 subroutine calc_sig_ppm_comp(npwc,nomega,rhotwgp,botsq,otq,omegame0i_io,zcut,theta_mu_minus_e0i,ket,ppmodel,npwx,npwc1,npwc2)
  use defs_basis
  implicit none
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc
  integer,intent(in) :: npwc1
  integer,intent(in) :: npwc2
  integer,intent(in) :: npwx
  integer,intent(in) :: ppmodel
  real(dp),intent(in) :: omegame0i_io
  real(dp),intent(in) :: theta_mu_minus_e0i
  real(dp),intent(in) :: zcut
  complex(gwpc),intent(in) :: botsq(npwc,npwc1)
  complex(gwpc),intent(inout) :: ket(npwc,nomega)
  real(dp),intent(in) :: otq(npwc,npwc2)
  complex(gwpc),intent(in) :: rhotwgp(npwx)
 end subroutine calc_sig_ppm_comp
end interface

interface
 subroutine calc_vHxc_braket(Dtset,nkibz,nkcalc,kcalc2ibz,b1,b2,nbnds,nsig_ab,npwvec,gsqcutf_eff,&  
  &  nfftf_tot,nfftf,ngfftf,igfftf,Wf_info,vpsp,vhartr,vxc,Cprj,Pawtab,Paw_an,Pawang,&  
  &  Pawfgrtab,Pawrad,Paw_ij,MPI_enreg,Cryst,tim_fourdp,rhor,rhog,usexcnhat,nhat,nhatgr,nhatgrdim,&  
  &  vhartr_me,vxcval_me,vxc_me,vUpaw_me)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: b1
  integer,intent(in) :: b2
  integer,intent(in) :: nbnds
  integer,intent(in) :: nfftf
  integer,intent(in) :: nfftf_tot
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkcalc
  integer,intent(in) :: nkibz
  integer,intent(in) :: npwvec
  integer,intent(in) :: nsig_ab
  integer,intent(in) :: tim_fourdp
  integer,intent(in) :: usexcnhat
  type(crystal_structure),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(mpi_type),intent(inout) :: MPI_enreg
  type(pawang_type),intent(in) :: Pawang
  type(wavefunctions_information),intent(inout) :: Wf_info
  real(dp),intent(in) :: gsqcutf_eff
  integer,intent(in) :: ngfftf(18)
  type(cprj_type),intent(in) :: Cprj(Cryst%natom,Dtset%nspinor*nbnds*nkibz*Dtset%nsppol*Dtset%usepaw)
  type(paw_an_type),intent(in) :: Paw_an(Cryst%natom)
  type(paw_ij_type),intent(in) :: Paw_ij(Cryst%natom)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom)
  type(pawrad_type),intent(in) :: Pawrad(Cryst%ntypat)
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Dtset%usepaw)
  integer,intent(in) :: igfftf(npwvec)
  integer,intent(in) :: kcalc2ibz(nkcalc)
  real(dp),intent(in) :: nhat(nfftf,Dtset%nspden*Dtset%usepaw)
  real(dp),intent(in) :: nhatgr(nfftf,Dtset%nspden,3*nhatgrdim)
  real(dp),intent(in) :: rhog(2,nfftf)
  real(dp),intent(in) :: rhor(nfftf,Dtset%nspden)
  complex(dpc),intent(out) :: vUpaw_me(b1:b2,b1:b2,nkibz,Dtset%nsppol*nsig_ab)
  real(dp),intent(in) :: vhartr(nfftf)
  complex(dpc),intent(out) :: vhartr_me(b1:b2,b1:b2,nkibz,Dtset%nsppol*nsig_ab)
  real(dp),intent(in) :: vpsp(nfftf)
  real(dp),intent(in) :: vxc(nfftf,Dtset%nspden)
  complex(dpc),intent(out) :: vxc_me(b1:b2,b1:b2,nkibz,Dtset%nsppol*nsig_ab)
  complex(dpc),intent(out) :: vxcval_me(b1:b2,b1:b2,nkibz,Dtset%nsppol*nsig_ab)
 end subroutine calc_vHxc_braket
end interface

interface
 subroutine calc_wf_qp(MPI_enreg,nkibz,nbnds,nsize,nsppol,nspinor,&  
  &  m_lda_to_qp,my_minb,my_maxb,b1gw,b2gw,wf)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: b1gw
  integer,intent(in) :: b2gw
  integer,intent(in) :: my_maxb
  integer,intent(in) :: my_minb
  integer,intent(in) :: nbnds
  integer,intent(in) :: nkibz
  integer,intent(in) :: nsize
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(mpi_type),intent(inout) :: MPI_enreg
  complex(dpc),intent(in) :: m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)
  complex(gwpc),intent(inout) :: wf(nsize*nspinor,my_minb:my_maxb,nkibz,nsppol)
 end subroutine calc_wf_qp
end interface

interface
 subroutine calc_wf_qp_Wfval(MPI_enreg,nkibz,nbnds,nsize,nsppol,nspinor,&  
  &  m_lda_to_qp,my_minb,my_maxb,b1gw,b2gw,wf,nbvw,wfval)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: b1gw
  integer,intent(in) :: b2gw
  integer,intent(in) :: my_maxb
  integer,intent(in) :: my_minb
  integer,intent(in) :: nbnds
  integer,intent(in) :: nbvw
  integer,intent(in) :: nkibz
  integer,intent(in) :: nsize
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(mpi_type),intent(inout) :: MPI_enreg
  complex(dpc),intent(in) :: m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)
  complex(gwpc),intent(inout) :: wf(nsize,my_minb:my_maxb,nkibz,nsppol)
  complex(gwpc),intent(inout) :: wfval(nsize,nbvw,nkibz,nsppol)
 end subroutine calc_wf_qp_Wfval
end interface

interface
 subroutine update_cprj(natom,nkibz,nbnds,nsppol,nspinor,m_lda_to_qp,dimlmn,Cprj_ibz)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nbnds
  integer,intent(in) :: nkibz
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(cprj_type),intent(inout) :: Cprj_ibz(natom,nspinor*nbnds*nkibz*nsppol)
  integer,intent(in) :: dimlmn(natom)
  complex(dpc),intent(in) :: m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)
 end subroutine update_cprj
end interface

interface
 subroutine ccgradvnl(npwwfn,nkibz,gvec,gprimd,kibz,Cryst,mpsang,vkbsign,vkb,vkbd,gradvnl)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mpsang
  integer,intent(in) :: nkibz
  integer,intent(in) :: npwwfn
  type(crystal_structure),intent(in) :: Cryst
  real(dp),intent(in) :: gprimd(3,3)
  complex(gwpc),intent(out) :: gradvnl(3,npwwfn,npwwfn,nkibz)
  integer,intent(in) :: gvec(3,npwwfn)
  real(dp),intent(in) :: kibz(3,nkibz)
  real(dp),intent(in) :: vkb(npwwfn,Cryst%ntypat,mpsang,nkibz)
  real(dp),intent(in) :: vkbd(npwwfn,Cryst%ntypat,mpsang,nkibz)
  real(dp),intent(in) :: vkbsign(mpsang,Cryst%ntypat)
 end subroutine ccgradvnl
end interface

interface
 subroutine ccgradvnl_ylm(npwwfn,nkibz,Cryst,gvec,gprimd,kibz,&  
  &  mpsang,vkbsign,vkb,vkbd,l_fnl,l_fnld)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mpsang
  integer,intent(in) :: nkibz
  integer,intent(in) :: npwwfn
  type(crystal_structure),intent(in) :: Cryst
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,npwwfn)
  real(dp),intent(in) :: kibz(3,nkibz)
  complex(gwpc),intent(out) :: l_fnl(npwwfn,mpsang**2,Cryst%natom,nkibz)
  complex(gwpc),intent(out) :: l_fnld(3,npwwfn,mpsang**2,Cryst%natom,nkibz)
  real(dp),intent(in) :: vkb(npwwfn,Cryst%ntypat,mpsang,nkibz)
  real(dp),intent(in) :: vkbd(npwwfn,Cryst%ntypat,mpsang,nkibz)
  real(dp),intent(in) :: vkbsign(mpsang,Cryst%ntypat)
 end subroutine ccgradvnl_ylm
end interface

interface
 subroutine cchi0(Dtset,Cryst,Dtfil,qpoint,Ep,Psps,Kmesh,Gsph_epsG0,Gsph_wfn,gvec,Pawang,Pawtab,nbv,nbvw,occ,&  
  &  ngfft,igfft,nfftot,gwenergy,chi0,MPI_enreg,ktabr,Ltg_q,my_minb,my_maxb,rhox_spl,Cprj_bz,chi0sumrule,&  
  &  Wf,Wf_val)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: my_maxb
  integer,intent(in) :: my_minb
  integer,intent(in) :: nbvw
  integer,intent(in) :: nfftot
  type(crystal_structure),intent(in) :: Cryst
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(in) :: Dtset
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(gvectors_type),intent(in) :: Gsph_wfn
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_q
  type(mpi_type),intent(inout) :: MPI_enreg
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(wavefunctions_information),intent(inout),optional :: Wf
  type(wavefunctions_information),intent(inout),optional :: Wf_val
  integer,intent(in) :: ngfft(18)
  type(cprj_type),intent(in) :: Cprj_bz(Cryst%natom,Dtset%nspinor*Ep%nbnds*Kmesh%nbz*Ep%nsppol*Dtset%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
  complex(gwpc),intent(out) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
  real(dp),intent(out) :: chi0sumrule(Ep%npwe)
  integer,intent(in) :: gvec(3,Ep%npwvec)
  real(dp),intent(in) :: gwenergy(Kmesh%nibz,Ep%nbnds,Ep%nsppol)
  integer,intent(in),target :: igfft(Ep%npwvec,2*Ep%mG0(1)+1,2*Ep%mG0(2)+1,2*Ep%mG0(3)+1)
  integer,intent(in) :: ktabr(nfftot,Kmesh%nbz)
  integer,intent(in) :: nbv(Ep%nsppol)
  real(dp),intent(in) :: occ(Kmesh%nibz,Ep%nbnds,Ep%nsppol)
  real(dp),intent(in) :: qpoint(3)
  real(dp),intent(in) :: rhox_spl(Psps%mqgrid_ff,2,0:2*(Psps%mpsang-1), &
  &         Psps%lnmax*(Psps%lnmax+1)/2,Psps%ntypat*Dtset%usepaw)
 end subroutine cchi0
end interface

interface
 function mkG0w(omega,numerator,delta_ene,zcut,TOL_W0) result(green0w)
  use defs_basis
  implicit none
  real(dp),intent(in) :: TOL_W0
  real(dp),intent(in) :: delta_ene
  complex(dpc) :: green0w
  real(dp),intent(in) :: numerator
  complex(gwpc),intent(in) :: omega
  real(dp),intent(in) :: zcut
 end function mkG0w
end interface

interface
 subroutine cchi0q0(Dtset,Cryst,Dtfil,Ep,Psps,Kmesh,Gsph_epsG0,Gsph_wfn,gvec,Pawang,Pawtab,ktabr,&  
  &  nbv,nbvw,occ,ngfft,igfft,nfftot,energy,gwenergy,chi0,MPI_enreg,Ltg_q,&  
  &  my_minb,my_maxb,pawrhox_spl,Cprj_bz,vkbsign,vkb,vkbd,chi0sumrule,&  
  &  Wf,Wf_val) ! Optional
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: my_maxb
  integer,intent(in) :: my_minb
  integer,intent(in) :: nbvw
  integer,intent(in) :: nfftot
  type(crystal_structure),intent(in) :: Cryst
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(in) :: Dtset
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(gvectors_type),intent(in) :: Gsph_wfn
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_q
  type(mpi_type),intent(inout) :: MPI_enreg
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(wavefunctions_information),optional,intent(inout) :: Wf
  type(wavefunctions_information),optional,intent(inout) :: Wf_val
  integer,intent(in) :: ngfft(18)
  type(cprj_type),intent(in) :: Cprj_bz(Cryst%natom,Dtset%nspinor*Ep%nbnds*Kmesh%nbz*Ep%nsppol*Dtset%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
  complex(gwpc),intent(out) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
  real(dp),intent(out) :: chi0sumrule(Ep%npwe)
  real(dp),intent(in) :: energy(Kmesh%nibz,Ep%nbnds,Ep%nsppol)
  integer,intent(in) :: gvec(3,Ep%npwvec)
  real(dp),intent(in) :: gwenergy(Kmesh%nibz,Ep%nbnds,Ep%nsppol)
  integer,target,intent(in) :: igfft(Ep%npwvec,2*Ep%mG0(1)+1,2*Ep%mG0(2)+1,Ep%mG0(3)+1)
  integer,intent(in) :: ktabr(nfftot,Kmesh%nbz)
  integer,intent(in) :: nbv(Ep%nsppol)
  real(dp),intent(in) :: occ(Kmesh%nibz,Ep%nbnds,Ep%nsppol)
  real(dp),intent(in) :: pawrhox_spl(Psps%mqgrid_ff,2,0:2*(Psps%mpsang-1), &
  &         Psps%lnmax*(Psps%lnmax+1)/2,Psps%ntypat*Dtset%usepaw)
  real(dp),intent(in) :: vkb(Ep%npwwfn,Dtset%ntypat,Psps%mpsang,Kmesh%nibz)
  real(dp),intent(in) :: vkbd(Ep%npwwfn,Dtset%ntypat,Psps%mpsang,Kmesh%nibz)
  real(dp),intent(in) :: vkbsign(Psps%mpsang,Dtset%ntypat)
 end subroutine cchi0q0
end interface

interface
 subroutine paw_inabla(natom,typat,Pawtab,Cprj_b1,Cprj_b2,onsite)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  type(cprj_type),intent(in) :: Cprj_b1(natom)
  type(cprj_type),intent(in) :: Cprj_b2(natom)
  type(pawtab_type),intent(in) :: Pawtab(:)
  real(dp),intent(out) :: onsite(2,3)
  integer,intent(in) :: typat(natom)
 end subroutine paw_inabla
end interface

interface
 subroutine apply_gradvnl(npwwfn,wfg1,wfg2,gradvnl,res) 
  use defs_basis
  implicit none
  integer,intent(in) :: npwwfn
  complex(gwpc),intent(in) :: gradvnl(3,npwwfn,npwwfn)
  complex(gwpc),intent(out) :: res(3)
  complex(gwpc),intent(in) :: wfg1(npwwfn)
  complex(gwpc),intent(in) :: wfg2(npwwfn)
 end subroutine apply_gradvnl
end interface

interface
 subroutine apply_gradvnl_Ylm(npwwfn,wfg1,wfg2,natom,mpsang,fnl,fnld,res) 
  use defs_basis
  implicit none
  integer,intent(in) :: mpsang
  integer,intent(in) :: natom
  integer,intent(in) :: npwwfn
  complex(gwpc),intent(in) :: fnl(npwwfn,mpsang**2,natom)
  complex(gwpc),intent(in) :: fnld(3,npwwfn,mpsang**2,natom)
  complex(gwpc),intent(out) :: res(3)
  complex(gwpc),intent(in) :: wfg1(npwwfn)
  complex(gwpc),intent(in) :: wfg2(npwwfn)
 end subroutine apply_gradvnl_Ylm
end interface

interface
 subroutine cggfft(npwvec,ngfft1,ngfft2,ngfft3,gvec,igfft)
  implicit none
  integer,intent(in) :: ngfft1
  integer,intent(in) :: ngfft2
  integer,intent(in) :: ngfft3
  integer,intent(in) :: npwvec
  integer,intent(in) :: gvec(3,npwvec)
  integer,intent(out) :: igfft(npwvec,npwvec)
 end subroutine cggfft
end interface

interface
 subroutine cigfft(mG0,npwvec,ngfft,gvec,igfft)
  implicit none
  integer,intent(in) :: npwvec
  integer,intent(in) :: mG0(3)
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: gvec(3,npwvec)
  integer,intent(out) :: igfft(npwvec,2*mG0(1)+1,2*mG0(2)+1,2*mG0(3)+1)
 end subroutine cigfft
end interface

interface
 subroutine clcqpg(npwx,gvec,gprimd,qq,nq,qpg)
  use defs_basis
  implicit none
  integer,intent(in) :: npwx
  integer,intent(in) :: nq
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,npwx)
  real(dp),intent(out) :: qpg(npwx,nq)
  real(dp),intent(in) :: qq(3,nq)
 end subroutine clcqpg
end interface

interface
 subroutine cppm1par(npwc,nqibz,nomega,epsm1,omega,bigomegatwsq,omegatw,omegaplasma)
  use defs_basis
  implicit none
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc
  integer,intent(in) :: nqibz
  real(dp),intent(in) :: omegaplasma
  complex(gwpc),intent(inout) :: bigomegatwsq(npwc,npwc,nqibz)
  complex(gwpc),intent(in) :: epsm1(npwc,npwc,nomega,nqibz)
  complex(gwpc),intent(in) :: omega(nomega)
  complex(gwpc),intent(inout) :: omegatw(npwc,npwc,nqibz)
 end subroutine cppm1par
end interface

interface
 subroutine cppm2par(paral_kgb,npwc,nqiA,nomega,epsm1,bigomegatwsq,omegatw,&  
  &  ngfftf,gvec,gprimd,rhor,nfftf,Qmesh,gmet,&  
  &  iqiA) ! Optional
  use defs_basis
  use defs_datatypes
  implicit none
  integer,optional,intent(in) :: iqiA
  integer,intent(in) :: nfftf
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc
  integer,intent(in) :: nqiA
  integer,intent(in) :: paral_kgb
  type(bz_mesh_type),intent(in) :: Qmesh
  integer,intent(in) :: ngfftf(18)
  complex(gwpc),intent(inout) :: bigomegatwsq(npwc,npwc,nqiA)
  complex(gwpc),intent(in) :: epsm1(npwc,npwc,nomega,nqiA)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,npwc)
  complex(gwpc),intent(inout) :: omegatw(npwc,npwc,nqiA)
  real(dp),intent(inout) :: rhor(nfftf)
 end subroutine cppm2par
end interface

interface
 subroutine cppm3par(paral_kgb,npwc,nqiA,nomega,epsm1,bigomegatwsq,omegatw,&  
  &  ngfftf,gvec,gprimd,rho,nfftf,eigtot,Qmesh,&  
  &  iqiA) ! Optional
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in),optional :: iqiA
  integer,intent(in) :: nfftf
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc
  integer,intent(in) :: nqiA
  integer,intent(in) :: paral_kgb
  type(bz_mesh_type),intent(in) :: Qmesh
  integer,intent(in) :: ngfftf(18)
  complex(gwpc),intent(inout) :: bigomegatwsq(npwc,1,nqiA)
  complex(gwpc),intent(inout) :: eigtot(npwc,npwc,nqiA)
  complex(gwpc),intent(in) :: epsm1(npwc,npwc,nomega,nqiA)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,npwc)
  complex(gwpc),intent(inout) :: omegatw(npwc,1,nqiA)
  real(dp),intent(inout) :: rho(nfftf)
 end subroutine cppm3par
end interface

interface
 subroutine cppm4par(paral_kgb,npwc,nqiA,epsm1,nomega,bigomegatwsq,omegatw,ngfftf,gvec,gprimd,rho,nfftf,Qmesh,&  
  &  iqia) ! Optional
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in),optional :: iqiA
  integer,intent(in) :: nfftf
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc
  integer,intent(in) :: nqiA
  integer,intent(in) :: paral_kgb
  type(bz_mesh_type),intent(in) :: Qmesh
  integer,intent(in) :: ngfftf(18)
  complex(gwpc),intent(inout) :: bigomegatwsq(npwc,npwc,nqiA)
  complex(gwpc),intent(in) :: epsm1(npwc,npwc,nomega,nqiA)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,npwc)
  complex(gwpc),intent(inout) :: omegatw(npwc,1,nqiA)
  real(dp),intent(inout) :: rho(nfftf)
 end subroutine cppm4par
end interface

interface
 subroutine cqratio(npwc,gvec,nq,q,gmet,gprimd,qratio)
  use defs_basis
  implicit none
  integer,intent(in) :: npwc
  integer,intent(in) :: nq
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,npwc)
  real(dp),intent(in) :: q(3,nq)
  real(dp),intent(out) :: qratio(npwc,npwc,nq)
 end subroutine cqratio
end interface

interface
 subroutine init_Crystal_structure(Cryst,natom,ntypat,nsym,rprimd,typat,xred,timrev,use_antiferro,remove_inv,&  
  &  symrel,tnons,symafm) ! optional
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: timrev
  type(crystal_structure),intent(inout) :: Cryst
  logical,intent(in) :: remove_inv
  logical,intent(in) :: use_antiferro
  real(dp),intent(in) :: rprimd(3,3)
  integer,optional,intent(in) :: symafm(nsym)
  integer,optional,intent(in) :: symrel(3,3,nsym)
  real(dp),optional,intent(in) :: tnons(3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine init_Crystal_structure
end interface

interface
 subroutine init_Crystal_from_Hdr(Cryst,Hdr,timrev,remove_inv)
  use defs_datatypes
  implicit none
  integer,intent(in) :: timrev
  type(crystal_structure),intent(out) :: Cryst
  type(hdr_type),intent(in) :: Hdr
  logical,optional,intent(in) :: remove_inv
 end subroutine init_Crystal_from_Hdr
end interface

interface
 subroutine nullify_Crystal_structure(Cryst)
  use defs_datatypes
  implicit none
  type(crystal_structure),intent(inout) :: Cryst
 end subroutine nullify_Crystal_structure
end interface

interface
 subroutine destroy_Crystal_structure(Cryst)
  use defs_datatypes
  implicit none
  type(crystal_structure),intent(inout) :: Cryst
 end subroutine destroy_Crystal_structure
end interface

interface
 subroutine print_Crystal_structure(Cryst,unit,mode_paral,prtvol) 
  use defs_datatypes
  implicit none
  integer,optional,intent(in) :: prtvol
  integer,optional,intent(in) :: unit
  type(crystal_structure),intent(in) :: Cryst
  character(len=4),optional,intent(in) :: mode_paral
 end subroutine print_Crystal_structure
end interface

interface
 subroutine print_symmetries(nsym,symrel,tnons,symafm,unit,mode_paral)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  integer,optional,intent(in) :: unit
  character(len=4),optional,intent(in) :: mode_paral
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
 end subroutine print_symmetries
end interface

interface
 function is_equal_Crystal_structure(Cryst1,Cryst2) result(equal)
  use defs_datatypes
  implicit none
  integer :: equal
  type(crystal_structure),intent(in) :: Cryst1
  type(crystal_structure),intent(in) :: Cryst2
 end function is_equal_Crystal_structure
end interface

interface
 subroutine remove_inversion(nsym,symrel,tnons,nsym_out,symrel_out,tnons_out,pinv)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  integer,intent(out) :: nsym_out
  integer,intent(out) :: pinv
  integer,pointer :: symrel_out(:,:,:)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),pointer :: tnons_out(:,:)
 end subroutine remove_inversion
end interface

interface
 subroutine csigme(ikcalc,Dtset,Cryst,Dtfil,Sp,Sr,Er,Gsph_Max,Vcp,Kmesh,Qmesh,Ltg_k,PPm,&  
  &  Pawang,Pawtab,Psps,Cprj_bz,Wf_info,Wf_info_braket,MPI_enreg,pawrhox,gvec,ktabr,ngfft,igfft,&  
  &  nfftgw_tot,en_qp,occ_qp,efermi_qp,my_minb,my_maxb,ngfftf,nfftf,rhor)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ikcalc
  integer,intent(in) :: my_maxb
  integer,intent(in) :: my_minb
  integer,intent(in) :: nfftf
  integer,intent(in) :: nfftgw_tot
  type(crystal_structure),intent(in) :: Cryst
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(in) :: Dtset
  type(epsilonm1_results),intent(inout) :: Er
  type(gvectors_type),intent(in) :: Gsph_Max
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_k
  type(mpi_type),intent(inout) :: MPI_enreg
  type(ppmodel_type),intent(inout) :: PPm
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(bz_mesh_type),intent(in) :: Qmesh
  type(sigma_parameters),intent(in) :: Sp
  type(sigma_results),intent(inout) :: Sr
  type(coulombian_type),intent(in) :: Vcp
  type(wavefunctions_information),intent(inout) :: Wf_info
  type(wavefunctions_information),intent(inout) :: Wf_info_braket
  real(dp),intent(in) :: efermi_qp
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftf(18)
  type(cprj_type),intent(in) :: Cprj_bz(Cryst%natom,Dtset%nspinor*Sp%nbnds*Kmesh%nbz*Sp%nsppol*Dtset%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
  real(dp),intent(in) :: en_qp(Kmesh%nibz,Sp%nbnds,Sp%nsppol)
  integer,intent(in) :: gvec(3,Sp%npwvec)
  integer,intent(in),target :: igfft(Sp%npwvec,2*Sp%mG0(1)+1,2*Sp%mG0(2)+1,2*Sp%mG0(3)+1)
  integer,intent(in) :: ktabr(nfftgw_tot,Kmesh%nbz)
  real(dp),intent(in) :: occ_qp(Kmesh%nibz,Sp%nbnds,Sp%nsppol)
  real(dp),intent(in) :: pawrhox(2,Sp%npwx,Psps%lmnmax*(Psps%lmnmax+1)/2, &
  &         Cryst%natom,Qmesh%nbz*Dtset%usepaw)
  real(dp),intent(inout) :: rhor(nfftf,Dtset%nspden)
 end subroutine csigme
end interface

interface
 subroutine cutoff_cylinder(nq,qpt,ng,gvec,rcut,hcyl,pdir,&  
  &  boxcenter,rprimd,vccut,method,MPI_enreg)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: method
  integer,intent(in) :: ng
  integer,intent(in) :: nq
  type(mpi_type),intent(in) :: MPI_enreg
  real(dp),intent(in) :: hcyl
  real(dp),intent(in) :: rcut
  integer,intent(in) :: pdir(3)
  real(dp),intent(in) :: boxcenter(3)
  integer,intent(in) :: gvec(3,ng)
  real(dp),intent(in) :: qpt(3,nq)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vccut(ng,nq)
 end subroutine cutoff_cylinder
end interface

interface
 subroutine cutoff_m_elem(ep,kmesh,gvec,symrec,Wf,energy,gwenergy,z0,wdth,occ,direction,gprimd)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: direction
  type(wavefunctions_information),optional,intent(in) :: Wf
  type(epsilonm1_parameters),intent(in) :: ep
  type(bz_mesh_type),target,intent(in) :: kmesh
  real(dp),intent(in) :: wdth
  real(dp),intent(in) :: z0
  real(dp),intent(in) :: energy(ep%nkibz,ep%nbnds,ep%nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,ep%npwvec)
  real(dp),intent(in) :: gwenergy(ep%nkibz,ep%nbnds,ep%nsppol)
  real(dp),intent(in) :: occ(ep%nkibz,ep%nbnds,ep%nsppol)
  integer,intent(in) :: symrec(3,3,kmesh%nsym)
 end subroutine cutoff_m_elem
end interface

interface
 subroutine matrixelements(npwwfn,wfg1,wfg2,gvec,kpoint,res)
  use defs_basis
  implicit none
  integer,intent(in) :: npwwfn
  integer,intent(in) :: gvec(3,npwwfn)
  real(dp),intent(in) :: kpoint(3)
  complex(gwpc),intent(out) :: res(3)
  complex(gwpc),intent(in) :: wfg1(npwwfn)
  complex(gwpc),intent(in) :: wfg2(npwwfn)
 end subroutine matrixelements
end interface

interface
 subroutine matrixelements_cutoff(npwwfn,wfg1,wfg2,gvec,kpoint,z0,wdth,direction,res)
  use defs_basis
  implicit none
  integer,intent(in) :: direction
  integer,intent(in) :: npwwfn
  real(dp), intent(in) :: wdth
  real(dp), intent(in) :: z0
  integer,intent(in) :: gvec(3,npwwfn)
  real(dp),intent(in) :: kpoint(3)
  complex(gwpc),intent(out) :: res(3)
  complex(gwpc),intent(in) :: wfg1(npwwfn)
  complex(gwpc),intent(in) :: wfg2(npwwfn)
 end subroutine matrixelements_cutoff
end interface

interface
 subroutine cutoff_sphere(nq,qpt,ng,gvec,gmet,rcut,vc_cut)
  use defs_basis
  implicit none
  integer,intent(in) :: ng
  integer,intent(in) :: nq
  real(dp),intent(in) :: rcut
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: gvec(3,ng)
  real(dp),intent(in) :: qpt(3,nq)
  real(dp),intent(out) :: vc_cut(ng,nq)
 end subroutine cutoff_sphere
end interface

interface
 subroutine cutoff_surface(nq,qpt,ng,gvec,gprimd,gmet,rcut,boxcenter,pdir,alpha,vc_cut,method)
  use defs_basis
  implicit none
  integer,intent(in) :: method
  integer,intent(in) :: ng
  integer,intent(in) :: nq
  real(dp),intent(in) :: rcut
  integer,intent(in) :: pdir(3)
  real(dp),intent(in) :: alpha(3)
  real(dp),intent(in) :: boxcenter(3)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,ng)
  real(dp),intent(in) :: qpt(3,nq)
  real(dp),intent(out) :: vc_cut(ng,nq)
 end subroutine cutoff_surface
end interface

interface
 subroutine cvc(nq,iq,q,npwvec,gvec,gprimd,qplusg)
  use defs_basis
  implicit none
  integer,intent(in) :: iq
  integer,intent(in) :: npwvec
  integer,intent(in) :: nq
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,npwvec)
  real(dp),intent(in) :: q(3,nq)
  real(dp),intent(out) :: qplusg(npwvec)
 end subroutine cvc
end interface

interface
 subroutine check_zarot(nsym,symrec,timrev,npwvec,rprimd,gprimd,gvec,psps,pawang,grottb,grottbm1)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: npwvec
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: grottb(npwvec,timrev,nsym)
  integer,intent(in) :: grottbm1(npwvec,timrev,nsym)
  integer,intent(in) :: gvec(3,npwvec)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine check_zarot
end interface

interface
 subroutine get_TDij0(npw,natom,ntypat,typat,kpt,gmet,gvec,pawtab,cprj1,cprj2,cg1,cg2,res)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: npw
  integer,intent(in) :: ntypat
  complex(gwpc),intent(in) :: cg1(npw)
  complex(gwpc),intent(in) :: cg2(npw)
  type(cprj_type),intent(in) :: cprj1(natom)
  type(cprj_type),intent(in) :: cprj2(natom)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: gvec(3,npw)
  real(dp),intent(in) :: kpt(3)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(out) :: res(2)
  integer,intent(in) :: typat(natom)
 end subroutine get_TDij0
end interface

interface
 subroutine test_PAWH(natom,ntypat,typat,nspinor,Kmesh,gmet,gvec,nfftf,vpsps,usepaw,pawtab,cprj_ibz,&  
  &  b1gw,b2gw,sp,en,vxc,vhartr,Wfs,ngfftf,igfftf,mpi_enreg)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: b1gw
  integer,intent(in) :: b2gw
  integer,intent(in) :: natom
  integer,intent(in) :: nfftf
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  type(bz_mesh_type),intent(in) :: Kmesh
  type(wavefunctions_information),intent(in) :: Wfs
  type(mpi_type),intent(inout) :: mpi_enreg
  type(sigma_parameters),intent(in) :: sp
  integer,intent(in) :: ngfftf(18)
  type(cprj_type),intent(in) :: cprj_ibz(natom,nspinor*Kmesh%nibz*sp%nbnds*sp%nsppol)
  real(dp),intent(in) :: en(Kmesh%nibz,sp%nbnds,sp%nsppol)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: gvec(3,Wfs%npwwfn)
  integer,intent(in) :: igfftf(Wfs%npwwfn)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
  complex(gwpc),intent(in) :: vhartr(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,sp%nsppol)
  real(dp),intent(in) :: vpsps(nfftf)
  complex(gwpc),intent(in) :: vxc(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,sp%nsppol)
 end subroutine test_PAWH
end interface

interface
 subroutine hack_symmetries(str,nsym,timrev)
  implicit none
  integer,intent(inout) :: nsym
  integer,intent(inout) :: timrev
  character(len=*),intent(in) :: str
 end subroutine hack_symmetries
end interface

interface
 function difvxc(rho)
  use defs_basis
  implicit none
  real(dp) :: difvxc
  real(dp),intent(in) :: rho
 end function difvxc
end interface

interface
 function diffvx(rho)
  use defs_basis
  implicit none
  real(dp) :: diffvx
  real(dp),intent(in) :: rho
 end function diffvx
end interface

interface
 function diffvc(rho)
  use defs_basis
  implicit none
  real(dp) :: diffvc
  real(dp),intent(in) :: rho
 end function diffvc
end interface

interface
 function difrel(rho)
  use defs_basis
  implicit none
  real(dp) :: difrel
  real(dp),intent(in) :: rho
 end function difrel
end interface

interface
 function vxnr(rho)
  use defs_basis
  implicit none
  real(dp),intent(in) :: rho
  real(dp) :: vxnr
 end function vxnr
end interface

interface
 function vxcca(rho)
  use defs_basis
  implicit none
  real(dp),intent(in) :: rho
  real(dp) :: vxcca
 end function vxcca
end interface

interface
 function vxjas(rho)
  use defs_basis
  implicit none
  real(dp),intent(in) :: rho
  real(dp) :: vxjas
 end function vxjas
end interface

interface
 function vcjas(rho)
  use defs_basis
  implicit none
  real(dp),intent(in) :: rho
  real(dp) :: vcjas
 end function vcjas
end interface

interface
 function rel(rho)
  use defs_basis
  implicit none
  real(dp) :: rel
  real(dp),intent(in) :: rho
 end function rel
end interface

interface
 subroutine ckxcldar(nr,rho2,kxclda)
  use defs_basis
  implicit none
  integer,intent(in) :: nr
  complex,intent(out) :: kxclda(nr)
  real(dp),intent(in) :: rho2(nr)
 end subroutine ckxcldar
end interface

interface
 subroutine ckxcldag(ngfft1,ngfft1a,ngfft2,ngfft3,nr,paral_kgb,rho2,kxclda)
  use defs_basis
  implicit none
  integer,intent(in) :: ngfft1
  integer,intent(in) :: ngfft1a
  integer,intent(in) :: ngfft2
  integer,intent(in) :: ngfft3
  integer,intent(in) :: nr
  integer,intent(in) :: paral_kgb
  complex,intent(out) :: kxclda(nr)
  real(dp),intent(in) :: rho2(nr)
 end subroutine ckxcldag
end interface

interface
 subroutine dosym(op,iinv,k1,k2)
  use defs_basis
  implicit none
  integer,intent(in) :: iinv
  real(dp),intent(in) :: k1(3)
  real(dp),intent(out) :: k2(3)
  real(dp),intent(in) :: op(3,3)
 end subroutine dosym
end interface

interface
 subroutine dosymr(op,iinv,r1,ngfft,r2)
  use defs_basis
  implicit none
  integer,intent(in) :: iinv
  integer,intent(in) :: ngfft(3)
  real(dp),intent(in) :: op(3,3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(out) :: r2(3)
 end subroutine dosymr
end interface

interface
 function dotproductqrc(r,c,b1,b2,b3)
  use defs_basis
  implicit none
  complex(gwpc) :: dotproductqrc
  real(dp),intent(in) :: b1(3)
  real(dp),intent(in) :: b2(3)
  real(dp),intent(in) :: b3(3)
  complex(gwpc),intent(in) :: c(3)
  real(dp),intent(in) :: r(3)
 end function dotproductqrc
end interface

interface
 subroutine eps1_tc(dtset,mpi_enreg,ngfft,nr,rho,rprimd,igfft,&  
  &  npweps,gmet,gprimd,gvec,nomega,epsm1,nq,ucvol,qq,omega,dtfil,hdr,npwwfn,npwvec,nbnds)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nbnds
  integer,intent(in) :: nomega
  integer,intent(in) :: npweps
  integer,intent(in) :: npwvec
  integer,intent(in) :: npwwfn
  integer,intent(in) :: nq
  integer,intent(in) :: nr
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  complex(gwpc),intent(inout) :: epsm1(npweps,npweps,nomega,nq)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,npweps)
  integer,intent(in) :: igfft(npweps)
  complex(gwpc),intent(in) :: omega(nomega)
  real(dp),intent(in) :: qq(3,nq)
  real(dp),intent(in) :: rho(nr,dtset%nsppol)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine eps1_tc
end interface

interface
 subroutine fermi(hdr,nbnds,nkibz,fixmom,nsppol,wtk,en,occ,nel,nbv,fermie)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nbnds
  integer,intent(inout) :: nel
  integer,intent(in) :: nkibz
  integer,intent(in) :: nsppol
  type(hdr_type),intent(in) :: Hdr
  real(dp),intent(out) :: fermie
  real(dp),intent(in) :: fixmom
  real(dp),intent(in) :: en(nkibz,nbnds,nsppol)
  integer,intent(out) :: nbv(nsppol)
  real(dp),intent(inout) :: occ(nkibz,nbnds,nsppol)
  real(dp),intent(in) :: wtk(nkibz)
 end subroutine fermi
end interface

interface
 subroutine fft_onewfn(paral_kgb,nspinor,npwwfn,nfftot,wfg,wfr,igfft,ngfft,tim_fourdp,MPI_enreg)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nfftot
  integer,intent(in) :: npwwfn
  integer,intent(in) :: nspinor
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourdp
  type(mpi_type),intent(inout) :: MPI_enreg
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: igfft(npwwfn)
  complex(gwpc),intent(in) :: wfg(npwwfn*nspinor)
  complex(gwpc),intent(out) :: wfr(ngfft(1)*ngfft(2)*ngfft(3)*nspinor)
 end subroutine fft_onewfn
end interface

interface
 subroutine findggp(nsym,symrec,Gsphere,qpt,Gpairs_q)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nsym
  type(gpairs_type),intent(inout) :: Gpairs_q
  type(gvectors_type),intent(in) :: Gsphere
  real(dp),intent(in) :: qpt(3)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine findggp
end interface

interface
 subroutine init_Gpairs_type(Gpairs_q,qpt,Gsphere,Cryst)
  use defs_basis
  use defs_datatypes
  implicit none
  type(crystal_structure),intent(in) :: Cryst
  type(gpairs_type),intent(out) :: Gpairs_q
  type(gvectors_type),intent(in) :: Gsphere
  real(dp),intent(in) :: qpt(3)
 end subroutine init_Gpairs_type
end interface

interface
 subroutine destroy_Gpairs_type(Gpairs_q)
  use defs_datatypes
  implicit none
  type(gpairs_type),intent(inout) :: Gpairs_q
 end subroutine destroy_Gpairs_type
end interface

interface
 subroutine nullify_Gpairs_type(Gpairs_q)
  use defs_datatypes
  implicit none
  type(gpairs_type),intent(inout) :: Gpairs_q
 end subroutine nullify_Gpairs_type
end interface

interface
 subroutine DEBUG_Gpairs(Gpairs_q,Gsphere,nsym,symrec)
  use defs_datatypes
  implicit none
  integer,intent(in) :: nsym
  type(gpairs_type),intent(in) :: Gpairs_q
  type(gvectors_type),intent(in) :: Gsphere
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine DEBUG_Gpairs
end interface

interface
 subroutine findk(nkcalc,nkbz,xkcalc,kbz,kcalc,umklp_opt)
  use defs_basis
  implicit none
  integer,intent(in) :: nkbz
  integer,intent(in) :: nkcalc
  integer,intent(in) :: umklp_opt
  real(dp),intent(in) :: kbz(3,nkbz)
  integer,intent(out) :: kcalc(nkcalc)
  real(dp),intent(in) :: xkcalc(3,nkcalc)
 end subroutine findk
end interface

interface
 subroutine findnq(nkbz,kbz,nsym,symrec,nqibz,timrev)
  use defs_basis
  implicit none
  integer,intent(in) :: nkbz
  integer,intent(out) :: nqibz
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  real(dp),intent(in) :: kbz(3,nkbz)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine findnq
end interface

interface
 subroutine findq(nkbz,kbz,nsym,symrec,gprimd,nqibz,qibz,timrev,avoid_zero)
  use defs_basis
  implicit none
  integer,intent(in) :: nkbz
  integer,intent(in) :: nqibz
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  logical,intent(in) :: avoid_zero
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: kbz(3,nkbz)
  real(dp),intent(out) :: qibz(3,nqibz)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine findq
end interface

interface
 subroutine findqg0(iq,g0,kmkp,nqbz,qbz,mG0)
  use defs_basis
  implicit none
  integer,intent(out) :: iq
  integer,intent(in) :: nqbz
  integer,intent(out) :: g0(3)
  integer,intent(in) :: mG0(3)
  real(dp),intent(in) :: kmkp(3)
  real(dp),intent(in) :: qbz(3,nqbz)
 end subroutine findqg0
end interface

interface
 subroutine fkin(ff,kincontrib,omegame0,otw1,otw2,zcut)
  use defs_basis
  implicit none
  integer,intent(in) :: ff
  real(dp),intent(out) :: kincontrib
  real(dp),intent(in) :: omegame0
  real(dp),intent(in) :: otw1
  real(dp),intent(in) :: otw2
  real(dp),intent(in) :: zcut
 end subroutine fkin
end interface

interface
 subroutine fsumrule(nomega,omega,eps,omegaplasma,method)
  use defs_basis
  implicit none
  integer,intent(in) :: method
  integer,intent(in) :: nomega
  real(dp),intent(in) :: omegaplasma
  real(dp),intent(in) :: eps(nomega)
  real(dp),intent(in) :: omega(nomega)
 end subroutine fsumrule
end interface

interface
 subroutine get_Bands_Sym_GW(ntypat,natom,typat,nspinor,nsppol,nbnds,nkibz,Kmesh,Wfs,usepaw,pawtab,&  
  &  pawang,psps,dimlmn,cprj_ibz,mpi_enreg,ik_ibz,ene_k,nsym,symrec,tnons,indsym,only_trace,irottb,BSym,ierr,&  
  &  EDIFF_TOL) ! optional
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: natom
  integer,intent(in) :: nbnds
  integer,intent(in) :: nkibz
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  type(bands_symmetries),intent(out) :: BSym
  real(dp),intent(in),optional :: EDIFF_TOL
  type(bz_mesh_type),intent(in) :: Kmesh
  type(wavefunctions_information),intent(inout) :: Wfs
  type(mpi_type),intent(inout) :: mpi_enreg
  logical,intent(in) :: only_trace
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  type(cprj_type),intent(in) :: cprj_ibz(natom,nspinor*nbnds*nkibz*nsppol*usepaw)
  integer,intent(in) :: dimlmn(natom)
  real(dp),intent(in) :: ene_k(nbnds,nsppol)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: irottb(Wfs%nfft,nsym)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  integer,intent(in) :: symrec(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(in) :: typat(natom)
 end subroutine get_Bands_Sym_GW
end interface

interface
 subroutine rotate_cprj(ntypat,natom,nbnds,psps,pawang,pawtab,typat,isym,cprj_in,cprj_out)
  use defs_datatypes
  implicit none
  integer,intent(in) :: isym
  integer,intent(in) :: natom
  integer,intent(in) :: nbnds
  integer,intent(in) :: ntypat
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  type(cprj_type),intent(in) :: cprj_in(natom,nbnds)
  type(cprj_type),intent(out) :: cprj_out(natom,nbnds)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine rotate_cprj
end interface

interface
 subroutine get_ng0sh(nk1,kbz1,nk2,kbz2,nkfold,kfold,gmet,tolq0,mg0sh,opt_ng0)
  use defs_basis
  implicit none
  integer,intent(in) :: mg0sh
  integer,intent(in) :: nk1
  integer,intent(in) :: nk2
  integer,intent(in) :: nkfold
  real(dp),intent(in) :: tolq0
  integer,intent(out) :: opt_ng0(3)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: kbz1(3,nk1)
  real(dp),intent(in) :: kbz2(3,nk2)
  real(dp),intent(in) :: kfold(3,nkfold)
 end subroutine get_ng0sh
end interface

interface
 subroutine nullify_epsilonm1_parameters(Ep)
  use defs_datatypes
  implicit none
  type(epsilonm1_parameters),intent(inout) :: Ep
 end subroutine nullify_epsilonm1_parameters
end interface

interface
 subroutine destroy_epsilonm1_parameters(Ep)
  use defs_datatypes
  implicit none
  type(epsilonm1_parameters),intent(inout) :: Ep
 end subroutine destroy_epsilonm1_parameters
end interface

interface
 subroutine nullify_sigma_parameters(Sp)
  use defs_datatypes
  implicit none
  type(sigma_parameters),intent(inout) :: Sp
 end subroutine nullify_sigma_parameters
end interface

interface
 subroutine destroy_sigma_parameters(Sp)
  use defs_datatypes
  implicit none
  type(sigma_parameters),intent(inout) :: Sp
 end subroutine destroy_sigma_parameters
end interface

interface
 subroutine nullify_sigma_results(Sr)
  use defs_datatypes
  implicit none
  type(sigma_results),intent(inout) :: Sr
 end subroutine nullify_sigma_results
end interface

interface
 subroutine init_sigma_results(Sp,Kmesh,Sr)
  use defs_datatypes
  implicit none
  type(bz_mesh_type),intent(in) :: Kmesh
  type(sigma_parameters),intent(in) :: Sp
  type(sigma_results),intent(inout) :: Sr
 end subroutine init_sigma_results
end interface

interface
 subroutine destroy_sigma_results(Sr)
  use defs_datatypes
  implicit none
  type(sigma_results),intent(inout) :: Sr
 end subroutine destroy_sigma_results
end interface

interface
 subroutine make_transitions(chi0alg,nbnds,nbvw,nsppol,symchi,timrev,TOL_DELTA_OCC,zcut,&  
  &  max_rest,min_rest,my_max_rest,my_min_rest,kmesh,ltg_q,mpi_enreg,mG0,gwenergy,occ,qpoint)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: chi0alg
  integer,intent(in) :: nbnds
  integer,intent(in) :: nbvw
  integer,intent(in) :: nsppol
  integer,intent(in) :: symchi
  integer,intent(in) :: timrev
  real(dp),intent(in) :: TOL_DELTA_OCC
  type(bz_mesh_type),intent(in) :: kmesh
  type(little_group),intent(in) :: ltg_q
  real(dp),intent(out) :: max_rest
  real(dp),intent(out) :: min_rest
  type(mpi_type    ),intent(in) :: mpi_enreg
  real(dp),intent(out) :: my_max_rest
  real(dp),intent(out) :: my_min_rest
  real(dp),intent(in) :: zcut
  integer,intent(in) :: mG0(3)
  real(dp),intent(in) :: gwenergy(kmesh%nibz,nbnds,nsppol)
  real(dp),intent(in) :: occ(kmesh%nibz,nbnds,nsppol)
  real(dp),intent(in) :: qpoint(3)
 end subroutine make_transitions
end interface

interface
 subroutine nullify_transitions(self)
  use defs_datatypes
  implicit none
  type(transitions_type),intent(out) :: self
 end subroutine nullify_transitions
end interface

interface
 subroutine init_transitions(self,nkbz,nbnds,nbvw,nsppol,nomega,qpoint,ntrans)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nbnds
  integer,intent(in) :: nbvw
  integer,intent(in) :: nkbz
  integer,intent(in) :: nomega
  integer,intent(in) :: nsppol
  integer,optional,intent(in) :: ntrans
  type(transitions_type),intent(out) :: self
  real(dp),intent(in) :: qpoint(3)
 end subroutine init_transitions
end interface

interface
 subroutine destroy_transitions(self)
  use defs_datatypes
  implicit none
  type(transitions_type),intent(inout) :: self
 end subroutine destroy_transitions
end interface

interface
 subroutine copy_transitions(t_in,t_out)
  use defs_datatypes
  implicit none
  type(transitions_type),intent(in) :: t_in
  type(transitions_type),intent(out) :: t_out
 end subroutine copy_transitions
end interface

interface
 subroutine print_transitions(self,unit,prtvol)
  use defs_datatypes
  implicit none
  integer,optional,intent(in) :: prtvol
  integer,optional,intent(in) :: unit
  type(transitions_type),intent(in) :: self
 end subroutine print_transitions
end interface

interface
 subroutine get_my_extrema(self,my_min_res,my_max_res)
  use defs_basis
  use defs_datatypes
  implicit none
  real(dp),intent(out) :: my_max_res
  real(dp),intent(out) :: my_min_res
  type(transitions_type),intent(in) :: self
 end subroutine get_my_extrema
end interface

interface
 subroutine split_transitions(self,mpi_enreg)
  use defs_datatypes
  implicit none
  type(mpi_type),intent(in) :: mpi_enreg
  type(transitions_type),intent(inout) :: self
 end subroutine split_transitions
end interface

interface
 subroutine find_my_indeces(self,nomega,omega_mesh,my_w1,my_w2)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: my_w1
  integer,intent(out) :: my_w2
  integer,intent(in) :: nomega
  type(transitions_type),intent(in) :: self
  real(dp) :: omega_mesh(nomega)
 end subroutine find_my_indeces
end interface

interface
 subroutine gw2abi(nkibz,nbnds,nsppol,gw_array,abi_array)
  use defs_basis
  implicit none
  integer,intent(in) :: nbnds
  integer,intent(in) :: nkibz
  integer,intent(in) :: nsppol
  real(dp),intent(out) :: abi_array(nbnds*nkibz*nsppol)
  real(dp),intent(in) :: gw_array(nkibz,nbnds,nsppol)
 end subroutine gw2abi
end interface

interface
 subroutine abi2gw(nkibz,nbnds,nsppol,abi_array,gw_array)
  use defs_basis
  implicit none
  integer,intent(in) :: nbnds
  integer,intent(in) :: nkibz
  integer,intent(in) :: nsppol
  real(dp),intent(in) :: abi_array(nbnds*nkibz*nsppol)
  real(dp),intent(out) :: gw_array(nkibz,nbnds,nsppol)
 end subroutine abi2gw
end interface

interface
 subroutine completechi0_deltapart(ik_bz,qzero,symchi,npwe,npwvec,nomega,nspinor,&  
  &  nfftot,ngfft,gvec,igfft0,Gsph_wfn,Ltg_q,green_enhigh_w,wfwfg,chi0)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(in) :: nfftot
  integer,intent(in) :: nomega
  integer,intent(in) :: npwe
  integer,intent(in) :: npwvec
  integer,intent(in) :: nspinor
  integer,intent(in) :: symchi
  type(gvectors_type),intent(in) :: Gsph_wfn
  type(little_group),intent(in) :: Ltg_q
  logical,intent(in) :: qzero
  integer,intent(in) :: ngfft(18)
  complex(gwpc),intent(inout) :: chi0(npwe,npwe,nomega)
  complex(dpc),intent(in) :: green_enhigh_w(nomega)
  integer,intent(in) :: gvec(3,npwvec)
  integer,intent(in) :: igfft0(npwvec)
  complex(gwpc),intent(in) :: wfwfg(nfftot*nspinor**2)
 end subroutine completechi0_deltapart
end interface

interface
 subroutine output_chi0sumrule(qeq0,iq,npwe,omegaplasma,chi0sumrule,epsm1_w0,vc_sqrt)
  use defs_basis
  implicit none
  integer,intent(in) :: iq
  integer,intent(in) :: npwe
  real(dp),intent(in) :: omegaplasma
  logical,intent(in) :: qeq0
  real(dp),intent(inout) :: chi0sumrule(npwe)
  complex(gwpc),intent(in) :: epsm1_w0(npwe,npwe)
  complex(gwpc),intent(in) :: vc_sqrt(npwe)
 end subroutine output_chi0sumrule
end interface

interface
 subroutine accumulate_chi0sumrule(ik_bz,symchi,npwe,factor,delta_ene,&  
  &  Ltg_q,Gsph_epsG0,npwepG0,rhotwg,chi0sumrule)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(in) :: npwe
  integer,intent(in) :: npwepG0
  integer,intent(in) :: symchi
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(little_group),intent(in) :: Ltg_q
  real(dp),intent(in) :: delta_ene
  real(dp),intent(in) :: factor
  real(dp),intent(inout) :: chi0sumrule(npwe)
  complex(gwpc),intent(in) :: rhotwg(npwepG0)
 end subroutine accumulate_chi0sumrule
end interface

interface
 subroutine hdr_vs_dtset(Hdr,Dtset)
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: Dtset
  type(hdr_type) :: Hdr
 end subroutine hdr_vs_dtset
end interface

interface
 subroutine identk(kibz,nkibz,nkbzmx,nsym,timrev,symrec,symafm,use_antiferro,kbz,ktab,ktabi,ktabo,nkbz,wtk)
  use defs_basis
  implicit none
  integer,intent(out) :: nkbz
  integer,intent(in) :: nkbzmx
  integer,intent(in) :: nkibz
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  logical,intent(in) :: use_antiferro
  real(dp),intent(out) :: kbz(3,nkbzmx)
  real(dp),intent(in) :: kibz(3,nkibz)
  integer,intent(out) :: ktab(nkbzmx)
  integer,intent(out) :: ktabi(nkbzmx)
  integer,intent(out) :: ktabo(nkbzmx)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  real(dp),intent(out) :: wtk(nkibz)
 end subroutine identk
end interface

interface
 subroutine identq(qibz,nqibz,nqbzX,symrec,nsym,timrev,wtq,qbz,qtab,qtabi,qtabo,nqbz)
  use defs_basis
  implicit none
  integer,intent(out) :: nqbz
  integer,intent(in) :: nqbzX
  integer,intent(in) :: nqibz
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  real(dp),intent(out) :: qbz(3,nqbzX)
  real(dp),intent(in) :: qibz(3,nqibz)
  integer,intent(out) :: qtab(nqbzX)
  integer,intent(out) :: qtabi(nqbzX)
  integer,intent(out) :: qtabo(nqbzX)
  real(dp),intent(in) :: symrec(3,3,nsym)
  real(dp),intent(out) :: wtq(nqibz)
 end subroutine identq
end interface

interface
 subroutine wrppm(dtfil,hdr,ppmodel,npwc,npwc2,npwc3,gvec,nq,qpoint,bigomegatwsq,omegatw,eigvec)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: npwc
  integer,intent(in) :: npwc2
  integer,intent(in) :: npwc3
  integer,intent(in) :: nq
  integer,intent(in) :: ppmodel
  type(datafiles_type),intent(in) :: dtfil
  type(hdr_type),intent(inout) :: hdr
  complex(gwpc),intent(in) :: bigomegatwsq(npwc,npwc2,nq)
  complex(gwpc),intent(in) :: eigvec(npwc,npwc,nq)
  integer,intent(in) :: gvec(3,npwc)
  complex(gwpc),intent(in) :: omegatw(npwc,npwc3,nq)
  real(dp),intent(in) :: qpoint(3,nq)
 end subroutine wrppm
end interface

interface
 subroutine testppm(unitppm,ppmodel,npwc,npwc2,npwc3,nq,hdr)
  use defs_datatypes
  implicit none
  integer,intent(out) :: npwc
  integer,intent(out) :: npwc2
  integer,intent(out) :: npwc3
  integer,intent(out) :: nq
  integer,intent(out) :: ppmodel
  integer,intent(in) :: unitppm
  type(hdr_type),intent(out) :: hdr
 end subroutine testppm
end interface

interface
 subroutine rdppm(unitppm,ppmodel,npwc,npwc2,npwc3,gvec,nq,qpoint,bigomegatwsq,omegatw,eigvec)
  use defs_basis
  implicit none
  integer,intent(in) :: npwc
  integer,intent(in) :: npwc2
  integer,intent(in) :: npwc3
  integer,intent(in) :: nq
  integer,intent(in) :: ppmodel
  integer,intent(in) :: unitppm
  complex(gwpc),intent(out) :: bigomegatwsq(npwc,npwc2,nq)
  complex(gwpc),intent(out) :: eigvec(npwc,npwc,nq)
  integer,intent(out) :: gvec(3,npwc)
  complex(gwpc),intent(out) :: omegatw(npwc,npwc3,nq)
  real(dp),intent(out) :: qpoint(3,nq)
 end subroutine rdppm
end interface

interface
 subroutine mkppm(ppmodel,ig,igp,iq,nomega,omega,npwc,npwc2,npwc3,gvec,nq,q,bigomegatwsq,omegatw,eigvec,em1w)
  use defs_basis
  implicit none
  integer,intent(in) :: ig
  integer,intent(in) :: igp
  integer,intent(in) :: iq
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc
  integer,intent(in) :: npwc2
  integer,intent(in) :: npwc3
  integer,intent(in) :: nq
  integer,intent(in) :: ppmodel
  complex(gwpc),intent(in) :: bigomegatwsq(npwc,npwc2,nq)
  complex(gwpc),intent(in) :: eigvec(npwc,npwc,nq)
  complex(gwpc),intent(out) :: em1w(nomega)
  integer,intent(in) :: gvec(3,npwc)
  real(dp),intent(in) :: omega(nomega)
  complex(gwpc),intent(in) :: omegatw(npwc,npwc3,nq)
  real(dp),intent(in) :: q(3,nq)
 end subroutine mkppm
end interface

interface
 subroutine write_sigma_results_header(Sp,Er,Cryst,Kmesh,Qmesh)
  use defs_datatypes
  implicit none
  type(crystal_structure),intent(in) :: Cryst
  type(epsilonm1_results),intent(in) :: Er
  type(bz_mesh_type),intent(in) :: Kmesh
  type(bz_mesh_type),intent(in) :: Qmesh
  type(sigma_parameters),intent(in) :: Sp
 end subroutine write_sigma_results_header
end interface

interface
 subroutine write_sigma_results(Sp,Sr,ikcalc,ikibz,Kmesh,usepawu,en_lda)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ikcalc
  integer,intent(in) :: ikibz
  integer,intent(in) :: usepawu
  type(bz_mesh_type),intent(in) :: Kmesh
  type(sigma_parameters),intent(in) :: Sp
  type(sigma_results),intent(in) :: Sr
  real(dp),intent(in) :: en_lda(Kmesh%nibz,Sp%nbnds,Sp%nsppol)
 end subroutine write_sigma_results
end interface

interface
 subroutine rdgw(nkibz,nbnds,nbv,nsppol,kibz,gwenergy)
  use defs_basis
  implicit none
  integer,intent(in) :: nbnds
  integer,intent(in) :: nkibz
  integer,intent(in) :: nsppol
  real(dp),intent(out) :: gwenergy(nkibz,nbnds,nsppol)
  real(dp),intent(in) :: kibz(3,nkibz)
  integer,intent(in) :: nbv(nsppol)
 end subroutine rdgw
end interface

interface
 subroutine print_Sigma_perturbative(Sp,SigR,ik_ibz,iband,isp,usepawu,unit,prtvol,mode_paral)
  use defs_datatypes
  implicit none
  integer,intent(in) :: iband
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: isp
  integer,intent(in),optional :: prtvol
  integer,intent(in),optional :: unit
  integer,intent(in) :: usepawu
  type(sigma_results),intent(in) :: SigR
  type(sigma_parameters),intent(in) :: Sp
  character(len=4),intent(in),optional :: mode_paral
 end subroutine print_Sigma_perturbative
end interface

interface
 subroutine print_Sigma_SC(SigR,ik_ibz,iband,isp,SigP,Kmesh,en_lda,usepawu,unit,prtvol,mode_paral)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iband
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: isp
  integer,intent(in),optional :: prtvol
  integer,intent(in),optional :: unit
  integer,intent(in) :: usepawu
  type(bz_mesh_type),intent(in) :: Kmesh
  type(sigma_parameters),intent(in) :: SigP
  type(sigma_results),intent(in) :: SigR
  character(len=4),intent(in),optional :: mode_paral
  real(dp),intent(in) :: en_lda(Kmesh%nibz,SigP%nbnds,SigP%nsppol)
 end subroutine print_Sigma_SC
end interface

interface
 subroutine print_QP(nbnds,nkibz,nsppol,m_lda_to_qp,ene_qp,ib_start,ib_stop,unit,prtvol,mode_paral)
  use defs_basis
  implicit none
  integer,intent(in) :: ib_start
  integer,intent(in) :: ib_stop
  integer,intent(in) :: nbnds
  integer,intent(in) :: nkibz
  integer,intent(in) :: nsppol
  integer,intent(in),optional :: prtvol
  integer,intent(in),optional :: unit
  character(len=4),intent(in),optional :: mode_paral
  real(dp),intent(in) :: ene_qp(nkibz,nbnds,nsppol)
  complex(dpc),intent(in) :: m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)
 end subroutine print_QP
end interface

interface
 subroutine joint_dos(nqpt,qpt,mband,nkibz,kibz,nsppol,nband,Cryst,eigen,occ,method,step,broad)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: method
  integer,intent(in) :: nkibz
  integer,intent(in) :: nqpt
  integer,intent(in) :: nsppol
  type(crystal_structure),intent(in) :: Cryst
  real(dp),intent(in) :: broad
  real(dp),intent(in) :: step
  real(dp),intent(in) :: eigen(mband*nkibz*nsppol)
  real(dp),intent(in) :: kibz(3,nkibz)
  integer,intent(in) :: nband(nkibz*nsppol)
  real(dp),intent(in) :: occ(mband*nkibz*nsppol)
  real(dp),intent(in) :: qpt(3,nqpt)
 end subroutine joint_dos
end interface

interface
 subroutine kramerskronig(nomega,omega,eps,method,only_check)
  use defs_basis
  implicit none
  integer,intent(in) :: method
  integer,intent(in) :: nomega
  integer,intent(in) :: only_check
  complex,intent(inout) :: eps(nomega)
  real(dp),intent(in) :: omega(nomega)
 end subroutine kramerskronig
end interface

interface
 subroutine make_epsm1_driver(iq,npwe,nI,nJ,nomega,nomegaer,omega,&  
  &  approx_type,option_test,nqibz,qibz,Vcp,Dtfil,gmet,kxc,MPI_enreg,chi0)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: approx_type
  integer,intent(in) :: iq
  integer,intent(in) :: nI
  integer,intent(in) :: nJ
  integer,intent(in) :: nomega
  integer,intent(in) :: nomegaer
  integer,intent(in) :: npwe
  integer,intent(in) :: nqibz
  integer,intent(in) :: option_test
  type(datafiles_type),intent(in) :: Dtfil
  type(mpi_type),intent(in) :: MPI_enreg
  type(coulombian_type),intent(in) :: Vcp
  complex(gwpc),intent(inout) :: chi0(npwe*nI,npwe*nJ,nomega)
  real(dp),intent(in) :: gmet(3,3)
  complex(gwpc),intent(in) :: kxc(npwe*approx_type,npwe*approx_type)
  complex(gwpc),intent(in) :: omega(nomega)
  real(dp),intent(in) :: qibz(3,nqibz)
 end subroutine make_epsm1_driver
end interface

interface
 subroutine outeps(nge,nomega,omega,eps,Gsphere,Gpairs_q,title,unt,prtvol)
  use defs_datatypes
  implicit none
  integer,intent(in) :: nge
  integer,intent(in) :: nomega
  integer,intent(in) :: prtvol
  integer,intent(in) :: unt
  type(gpairs_type),intent(in) :: Gpairs_q
  type(gvectors_type),intent(in) :: Gsphere
  character(len=500),intent(in) :: title
  complex,intent(in) :: eps(nge,nge,nomega)
  complex,intent(in) :: omega(nomega)
 end subroutine outeps
end interface

interface
 subroutine paw_mkrhox(Cryst,pwff_spl,gmet,gvec,method,dim1_rhox,dim2_rhox,&  
  &  Psps,Pawang,Pawtab,qpt,npw,ylm_q,paw_rhox)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: dim1_rhox
  integer,intent(in) :: dim2_rhox
  integer,intent(in) :: method
  integer,intent(in) :: npw
  type(crystal_structure),intent(in) :: Cryst
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: gvec(3,npw)
  real(dp),intent(out) :: paw_rhox(2,npw,Psps%lmnmax*(Psps%lmnmax+1)/2,Cryst%natom)
  real(dp),intent(in) :: pwff_spl(Psps%mqgrid_ff,2,0:dim1_rhox,dim2_rhox,Cryst%ntypat)
  real(dp),intent(in) :: qpt(3)
  real(dp),intent(in) :: ylm_q(npw,(2*Psps%mpsang-1)**2)
 end subroutine paw_mkrhox
end interface

interface
 subroutine paw_mkrhox_spl(ntypat,Psps,Pawrad,Pawtab,pwff_spl,method,dim1_rhox,dim2_rhox)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: dim1_rhox
  integer,intent(in) :: dim2_rhox
  integer,intent(in) :: method
  integer,intent(in) :: ntypat
  type(pseudopotential_type), intent(in) :: Psps
  type(pawrad_type),intent(in) :: Pawrad(ntypat)
  type(pawtab_type),intent(in) :: Pawtab(ntypat)
  real(dp),intent(out) :: pwff_spl(Psps%mqgrid_ff,2,0:dim1_rhox,dim2_rhox,ntypat)
 end subroutine paw_mkrhox_spl
end interface

interface
 subroutine paw_rho_tw_g(npw,dim_rtwg,nspinor,natom,lmnmax,dimlmn,Cprj_k1,Cprj_k2,paw_rhox,rhotwg)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: dim_rtwg
  integer,intent(in) :: lmnmax
  integer,intent(in) :: natom
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  type(cprj_type),intent(in) :: Cprj_k1(natom,nspinor)
  type(cprj_type),intent(in) :: Cprj_k2(natom,nspinor)
  integer,intent(in) :: dimlmn(natom)
  real(dp),intent(in) :: paw_rhox(2,npw,lmnmax*(lmnmax+1)/2,natom)
  complex(gwpc),intent(inout) :: rhotwg(npw*dim_rtwg)
 end subroutine paw_rho_tw_g
end interface

interface
 subroutine paw_symcprj(Pawtab,Cryst,nspinor,mband,nband,&  
  &  nsppol,Psps,Kmesh,Cprj_ibz,Pawang,dimlmn,Cprj_bz)
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(crystal_structure),intent(in) :: Cryst
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(cprj_type),intent(out) :: Cprj_bz(Cryst%natom,nspinor*mband*Kmesh%nbz*nsppol)
  type(cprj_type),intent(in) :: Cprj_ibz(Cryst%natom,nspinor*mband*Kmesh%nibz*nsppol)
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat)
  integer,intent(in) :: dimlmn(Cryst%natom)
  integer,intent(in) :: nband(Kmesh%nibz*nsppol)
 end subroutine paw_symcprj
end interface

interface
 subroutine destroy_paw_ij(Paw_ij)
  use defs_datatypes
  implicit none
  type(paw_ij_type),intent(inout) :: Paw_ij(:)
 end subroutine destroy_paw_ij
end interface

interface
 subroutine init_paw_ij(Paw_ij,cplex,cplex_dij,nspinor,nsppol,nspden,pawspnorb,natom,ntypat,typat,Pawtab,&  
  &  has_dijhartree,has_dijhat,has_dijxc,has_dijxc_val,has_dijso,has_dijU) ! Optional
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: cplex_dij
  integer,optional,intent(in) :: has_dijU
  integer,optional,intent(in) :: has_dijhartree
  integer,optional,intent(in) :: has_dijhat
  integer,optional,intent(in) :: has_dijso
  integer,optional,intent(in) :: has_dijxc
  integer,optional,intent(in) :: has_dijxc_val
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawspnorb
  type(paw_ij_type),intent(inout) :: Paw_ij(natom)
  type(pawtab_type),intent(in) :: Pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine init_paw_ij
end interface

interface
 subroutine destroy_Paw_an(Paw_an)
  use defs_datatypes
  implicit none
  type(paw_an_type),intent(inout) :: Paw_an(:)
 end subroutine destroy_Paw_an
end interface

interface
 subroutine init_pawfgr(Dtset,k0,gmet,Pawfgr,mgfftf,nfftf,ecut_eff,ecutdg_eff,&  
  &  gsqcutc_eff,gsqcutf_eff,ngfftc,ngfftf)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: mgfftf
  integer,intent(out) :: nfftf
  type(dataset_type),intent(in) :: Dtset
  type(pawfgr_type),intent(out) :: Pawfgr
  real(dp),intent(out) :: ecut_eff
  real(dp),intent(out) :: ecutdg_eff
  real(dp),intent(out) :: gsqcutc_eff
  real(dp),intent(out) :: gsqcutf_eff
  integer,intent(out) :: ngfftc(18)
  integer,intent(out) :: ngfftf(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: k0(3)
 end subroutine init_pawfgr
end interface

interface
 subroutine nullify_paw_ij(Paw_ij)
  use defs_datatypes
  implicit none
  type(paw_ij_type),intent(inout) :: Paw_ij(:)
 end subroutine nullify_paw_ij
end interface

interface
 subroutine nullify_paw_an(Paw_an)
  use defs_datatypes
  implicit none
  type(paw_an_type),intent(inout) :: Paw_an(:)
 end subroutine nullify_paw_an
end interface

interface
 subroutine init_paw_an(natom,ntypat,nspden,cplex,pawxcdev,pawspnorb,typat,Pawang,Pawtab,Paw_an,&  
  &  has_vxcval) ! Optional
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,optional,intent(in) :: has_vxcval
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawspnorb
  integer,intent(in) :: pawxcdev
  type(pawang_type),intent(in) :: Pawang
  type(paw_an_type),intent(inout) :: Paw_an(:)
  type(pawtab_type),intent(in) :: Pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine init_paw_an
end interface

interface
 subroutine pawr(iat,Pawtab,Pawrad,Pawang,Psps,natom,ntypat,typat,xcart,lmn2_size,rc_onsite)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iat
  integer,intent(in) :: lmn2_size
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(pawrad_type),intent(in) :: Pawrad(ntypat)
  type(pawtab_type),intent(in) :: Pawtab(ntypat)
  real(dp),intent(inout) :: rc_onsite(3,lmn2_size)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xcart(3,natom)
 end subroutine pawr
end interface

interface
 subroutine pclock(itimpt)
  implicit none
  integer,intent(in) :: itimpt
 end subroutine pclock
end interface

interface
 subroutine PPmodel_symmetrize(PPm,Gsph,Qmesh,iq_bz,botsq,otq,eig) 
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iq_bz
  type(gvectors_type),intent(in) :: Gsph
  type(ppmodel_type),intent(in) :: PPm
  type(bz_mesh_type),intent(in) :: Qmesh
  complex(gwpc),intent(out) :: botsq(PPm%npwc,PPm%dm2_botsq)
  complex(gwpc),intent(out) :: eig(PPm%dm_eig,PPm%dm_eig)
  real(dp),intent(out) :: otq(PPm%npwc,PPm%dm2_otq)
 end subroutine PPmodel_symmetrize
end interface

interface
 subroutine nullify_PPmodel(PPm)
  use defs_datatypes
  implicit none
  type(ppmodel_type),intent(inout) :: PPm
 end subroutine nullify_PPmodel
end interface

interface
 subroutine destroy_PPmodel(PPm)
  use defs_datatypes
  implicit none
  type(ppmodel_type),intent(inout) :: PPm
 end subroutine destroy_PPmodel
end interface

interface
 subroutine init_PPmodel(PPm,Qmesh,ppmodel,npwc,mqmem,drude_plsmf)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mqmem
  integer,intent(in) :: npwc
  integer,intent(in) :: ppmodel
  type(ppmodel_type),intent(out) :: PPm
  type(bz_mesh_type),intent(in) :: Qmesh
  real(dp),intent(in) :: drude_plsmf
 end subroutine init_PPmodel
end interface

interface
 subroutine print_paw_ij(Paw_ij,pawprtvol)
  use defs_datatypes
  implicit none
  integer,intent(in) :: pawprtvol
  type(paw_ij_type),intent(in) :: Paw_ij(:)
 end subroutine print_paw_ij
end interface

interface
 subroutine print_pawtab(Pawtab,unitno,prtvol,mode_paral)
  use defs_datatypes
  implicit none
  integer,intent(in),optional :: prtvol
  integer,intent(in),optional :: unitno
  character(len=4),intent(in),optional :: mode_paral
  type(pawtab_type) :: Pawtab(:)
 end subroutine print_pawtab
end interface

interface
 subroutine print_psps(psps,unit,prtvol,mode_paral)
  use defs_datatypes
  implicit none
  integer,intent(in),optional :: prtvol
  integer,intent(in),optional :: unit
  character(len=4),intent(in),optional :: mode_paral
  type(pseudopotential_type),intent(in) :: psps
 end subroutine print_psps
end interface

interface
 subroutine plot_psps(psps,root_filename)
  use defs_basis
  use defs_datatypes
  implicit none
  type(pseudopotential_type),intent(in) :: psps
  character(len=fnlen),intent(in),optional :: root_filename
 end subroutine plot_psps
end interface

interface
 subroutine q0fit(nq,q,gvec,nomega,omega,npwvec,chi0,qcut,metal,&  
  &  nop,op,ninv,gprimd)
  use defs_basis
  implicit none
  integer,intent(in) :: ninv
  integer,intent(in) :: nomega
  integer,intent(in) :: nop
  integer,intent(in) :: npwvec
  integer,intent(in) :: nq
  logical,intent(in) :: metal
  real(dp),intent(in) :: qcut
  complex(gwpc),intent(inout) :: chi0(nq,npwvec,npwvec,nomega)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,npwvec)
  complex(gwpc),intent(in) :: omega(nomega)
  real(dp),intent(in) :: op(3,3,nop)
  real(dp),intent(in) :: q(3,nq)
 end subroutine q0fit
end interface

interface
 subroutine rdqps(gwcalctyp,Dtfil,Kmesh,nbnds,nspden,nsppol,nscf,nfftot,en,nbsc,en_qp,m_lda_to_qp,rho_p)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: gwcalctyp
  integer,intent(in) :: nbnds
  integer,intent(out) :: nbsc
  integer,intent(in) :: nfftot
  integer,intent(out) :: nscf
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: Dtfil
  type(bz_mesh_type),intent(in) :: Kmesh
  real(dp),intent(in) :: en(Kmesh%nibz,nbnds,nsppol)
  real(dp),intent(out) :: en_qp(Kmesh%nibz,nbnds,nsppol)
  complex(dpc),intent(out) :: m_lda_to_qp(nbnds,nbnds,Kmesh%nibz,nsppol)
  real(dp),intent(inout) :: rho_p(nfftot,nspden)
 end subroutine rdqps
end interface

interface
 subroutine rdscr(optfil,unt,fname,npweA,nqibz,nqibzA,nomegaA,qibz,omegaA,gmet,&  
  &  epsm1,MPI_enreg,localrdwf,nqlwlA,verbose,qlwl,uwing,lwing,&  
  &  iqiA) ! Optional
  use defs_basis
  use defs_datatypes
  implicit none
  integer,optional,intent(in) :: iqiA
  integer,intent(in) :: localrdwf
  integer,intent(in) :: nomegaA
  integer,intent(in) :: npweA
  integer,intent(in) :: nqibz
  integer,intent(in) :: nqibzA
  integer,intent(in) :: nqlwlA
  integer,intent(in) :: optfil
  integer,intent(in) :: unt
  type(mpi_type),intent(in) :: MPI_enreg
  character(len=fnlen),intent(in) :: fname
  logical,intent(in) :: verbose
  complex(gwpc),intent(inout) :: epsm1(npweA,npweA,nomegaA,nqibzA)
  real(dp),intent(in) :: gmet(3,3)
  complex(gwpc),intent(inout) :: lwing(npweA,nomegaA,nqlwlA*optfil)
  complex(gwpc),intent(inout) :: omegaA(nomegaA)
  real(dp),intent(out) :: qibz(3,nqibz)
  real(dp),intent(inout) :: qlwl(3,nqlwlA*optfil)
  complex(gwpc),intent(inout) :: uwing(npweA,nomegaA,nqlwlA*optfil)
 end subroutine rdscr
end interface

interface
 subroutine nullify_epsilonm1_results(Er)
  use defs_datatypes
  implicit none
  type(epsilonm1_results),intent(inout) :: Er
 end subroutine nullify_epsilonm1_results
end interface

interface
 subroutine destroy_epsilonm1_results(Er)
  use defs_datatypes
  implicit none
  type(epsilonm1_results),intent(inout) :: Er
 end subroutine destroy_epsilonm1_results
end interface

interface
 subroutine print_epsilonm1_results(Er,unit,prtvol,mode_paral)
  use defs_datatypes
  implicit none
  integer,optional,intent(in) :: prtvol
  integer,optional,intent(in) :: unit
  type(epsilonm1_results),intent(in) :: Er
  character(len=4),optional,intent(in) :: mode_paral
 end subroutine print_epsilonm1_results
end interface

interface
 subroutine Epsm1_symmetrizer(iq_bz,nomega,npwc,Er,Gsph,Qmesh,remove_exchange,epsm1_qbz) 
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iq_bz
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc
  type(epsilonm1_results),intent(in) :: Er
  type(gvectors_type),intent(in) :: Gsph
  type(bz_mesh_type),intent(in) :: Qmesh
  logical,intent(in) :: remove_exchange
  complex(gwpc),intent(out) :: epsm1_qbz(npwc,npwc,nomega)
 end subroutine Epsm1_symmetrizer
end interface

interface
 subroutine init_Er_from_file(Er,optfil,unt,fname,mqmem,accesswff,localrdwf,MPI_enreg)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: accesswff
  integer,intent(in) :: localrdwf
  integer,intent(in) :: mqmem
  integer,intent(in) :: optfil
  integer,intent(in) :: unt
  type(epsilonm1_results),intent(inout) :: Er
  type(mpi_type),intent(in) :: MPI_enreg
  character(len=fnlen),intent(in) :: fname
 end subroutine init_Er_from_file
end interface

interface
 subroutine mkdump_Er(Er,Dtset,Dtfil,unt_dump,fname_dump,accesswff,localrdwf,gmet,MPI_enreg)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: accesswff
  integer,intent(in) :: localrdwf
  integer,intent(in) :: unt_dump
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(in) :: Dtset
  type(epsilonm1_results),intent(inout) :: Er
  type(mpi_type),intent(in) :: MPI_enreg
  character(len=fnlen),intent(in) :: fname_dump
  real(dp),intent(in) :: gmet(3,3)
 end subroutine mkdump_Er
end interface

interface
 subroutine get_epsm1(Er,Vcp,Dtfil,approx_type,option_test,accesswff,localrdwf,kxcg,gmet,MPI_enreg,iqibzA)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: accesswff
  integer,intent(in) :: approx_type
  integer,optional,intent(in) :: iqibzA
  integer,intent(in) :: localrdwf
  integer,intent(in) :: option_test
  type(datafiles_type),intent(in) :: Dtfil
  type(epsilonm1_results),intent(inout) :: Er
  type(mpi_type),intent(in) :: MPI_enreg
  type(coulombian_type),intent(in) :: Vcp
  real(dp),intent(in) :: gmet(3,3)
  complex(gwpc),intent(in) :: kxcg(:,:)
 end subroutine get_epsm1
end interface

interface
 subroutine rho_tw_g(paral_kgb,nspinor,npwvec,nr,ngfft,map2sphere,igfftg0,&  
  &  wfn1,i1,ktabr1,ktabp1,spinrot1,&  
  &  wfn2,i2,ktabr2,ktabp2,spinrot2,&  
  &  dim_rtwg,rhotwg,tim_fourdp,MPI_enreg)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: dim_rtwg
  integer,intent(in) :: i1
  integer,intent(in) :: i2
  integer,intent(in) :: map2sphere
  integer,intent(in) :: npwvec
  integer,intent(in) :: nr
  integer,intent(in) :: nspinor
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourdp
  type(mpi_type),intent(inout) :: MPI_enreg
  complex(gwpc),intent(in) :: ktabp1
  complex(gwpc),intent(in) :: ktabp2
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: igfftg0(npwvec*map2sphere)
  integer,intent(in) :: ktabr1(nr)
  integer,intent(in) :: ktabr2(nr)
  complex(gwpc),intent(out) :: rhotwg(npwvec*dim_rtwg)
  real(dp),intent(in) :: spinrot1(4)
  real(dp),intent(in) :: spinrot2(4)
  complex(gwpc),intent(in) :: wfn1(nr*nspinor)
  complex(gwpc),intent(in) :: wfn2(nr*nspinor)
 end subroutine rho_tw_g
end interface

interface
 subroutine setmesh(gmet,gvec,ngfft,npwvec,npwsigx,npwwfn,nfftot,method,mG0,Cryst,enforce_sym)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: enforce_sym
  integer,intent(in) :: method
  integer,intent(out) :: nfftot
  integer,intent(in) :: npwsigx
  integer,intent(in) :: npwvec
  integer,intent(in) :: npwwfn
  type(crystal_structure),intent(in) :: Cryst
  integer,intent(in) :: mG0(3)
  integer,intent(inout) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: gvec(3,npwvec)
 end subroutine setmesh
end interface

interface
 subroutine select_divisor_mesh(n,nsym,tnons)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  integer,intent(inout) :: n(3)
  real(dp),intent(in) :: tnons(3,nsym)
 end subroutine select_divisor_mesh
end interface

interface
 integer function divisor_tns(dd)
 use defs_basis
 implicit none
 real(dp),intent(in) :: dd
end function divisor_tns
end interface

interface
 integer function mcm(ii,jj)
 implicit none
 integer,intent(in) :: ii
 integer,intent(in) :: jj
end function mcm
end interface

interface
 logical function check_rot_fft(nsym,symrel,nr1,nr2,nr3)
 implicit none
 integer,intent(in) :: nr1
 integer,intent(in) :: nr2
 integer,intent(in) :: nr3
 integer,intent(in) :: nsym
 integer,intent(in) :: symrel(3,3,nsym)
end function check_rot_fft
end interface

interface
 subroutine setshells(ecut,npw,nsh,nsym,gmet,gprimd,symrel,tag,ucvol)
  use defs_basis
  implicit none
  integer,intent(inout) :: npw
  integer,intent(inout) :: nsh
  integer,intent(in) :: nsym
  real(dp),intent(inout) :: ecut
  character(len=*),intent(in) :: tag
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine setshells
end interface

interface
 subroutine setup_FFT_rotation(nsym,symrec,tnons,nfft,ngfft,irottb)
  use defs_basis
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nsym
  integer,intent(in) :: ngfft(18)
  integer,intent(out) :: irottb(nfft,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
 end subroutine setup_FFT_rotation
end interface

interface
 subroutine FFT_rotations(Cryst,ngfft,irottb)
  use defs_datatypes
  implicit none
  type(crystal_structure),intent(in) :: Cryst
  integer,intent(in) :: ngfft(18)
  integer,intent(out) :: irottb(ngfft(1)*ngfft(2)*ngfft(3),Cryst%nsym)
 end subroutine FFT_rotations
end interface

interface
 subroutine setup_G_rotation(only_one_kpt,nsym,symrec,timrev,npw,gvec,g2sh,nsh,shlim,grottb,grottbm1)
  implicit none
  integer,intent(in) :: npw
  integer,intent(in) :: nsh
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  logical,intent(in) :: only_one_kpt
  integer,intent(in) :: g2sh(npw)
  integer,intent(inout) :: grottb(npw,timrev,nsym)
  integer,intent(inout) :: grottbm1(npw,timrev,nsym)
  integer,intent(in) :: gvec(3,npw)
  integer,intent(in) :: shlim(nsh+1)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine setup_G_rotation
end interface

interface
 subroutine init_Gvectors_type(only_one_kpt,Gsph,Cryst,ng,gvec,gmet,gprimd)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ng
  type(crystal_structure),intent(in) :: Cryst
  type(gvectors_type),intent(out) :: Gsph
  logical,intent(in) :: only_one_kpt
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,ng)
 end subroutine init_Gvectors_type
end interface

interface
 subroutine print_Gvectors(Gsph,unit,prtvol,mode_paral)
  use defs_datatypes
  implicit none
  integer,intent(in),optional :: prtvol
  integer,intent(in),optional :: unit
  type(gvectors_type),intent(in) :: Gsph
  character(len=4),intent(in),optional :: mode_paral
 end subroutine print_Gvectors
end interface

interface
 subroutine destroy_Gvectors(Gsph)
  use defs_datatypes
  implicit none
  type(gvectors_type),intent(inout) :: Gsph
 end subroutine destroy_Gvectors
end interface

interface
 subroutine nullify_Gvectors(Gsph)
  use defs_datatypes
  implicit none
  type(gvectors_type),intent(inout) :: Gsph
 end subroutine nullify_Gvectors
end interface

interface
 subroutine setup_coulombian(Dtset,Gsph,Qmesh,Kmesh,ng,rprimd,ngfft,MPI_enreg,Vcp)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ng
  type(dataset_type),intent(in) :: Dtset
  type(gvectors_type),intent(in) :: Gsph
  type(bz_mesh_type),intent(in) :: Kmesh
  type(mpi_type),intent(inout) :: MPI_enreg
  type(bz_mesh_type),intent(in) :: Qmesh
  type(coulombian_type),intent(out) :: Vcp
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine setup_coulombian
end interface

interface
 subroutine plot_Vc(Vcp,Qmesh,Gsph,ng,vc,MPI_enreg)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ng
  type(gvectors_type),intent(in) :: Gsph
  type(mpi_type),intent(inout) :: MPI_enreg
  type(bz_mesh_type),intent(in) :: Qmesh
  type(coulombian_type),intent(in) :: Vcp
  real(dp),intent(in) :: vc(ng,Qmesh%nibz)
 end subroutine plot_Vc
end interface

interface
 subroutine print_coulombian(Vcp,unit,prtvol,mode_paral)
  use defs_datatypes
  implicit none
  integer,intent(in),optional :: prtvol
  integer,intent(in),optional :: unit
  type(coulombian_type),intent(in) :: Vcp
  character(len=4),intent(in),optional :: mode_paral
 end subroutine print_coulombian
end interface

interface
 subroutine cutoff_table(Vcp)
  use defs_datatypes
  implicit none
  type(coulombian_type),intent(inout) :: Vcp
 end subroutine cutoff_table
end interface

interface
 subroutine cutoff_density(ngfft,nspden,nsppol,Vcp,rhor,MPI_enreg)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  type(mpi_type),intent(in) :: MPI_enreg
  type(coulombian_type),intent(in) :: Vcp
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: rhor(ngfft(1)*ngfft(2)*ngfft(3),nspden)
 end subroutine cutoff_density
end interface

interface
 subroutine destroy_Coulombian(Vcp) 
  use defs_datatypes
  implicit none
  type(coulombian_type),intent(inout) :: Vcp
 end subroutine destroy_Coulombian
end interface

interface
 subroutine setup_Kmesh(nkibz,kibz,Cryst,Kmesh,prtvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nkibz
  integer,intent(in) :: prtvol
  type(crystal_structure),intent(in) :: Cryst
  type(bz_mesh_type),intent(inout) :: Kmesh
  real(dp),intent(in) :: kibz(3,nkibz)
 end subroutine setup_Kmesh
end interface

interface
 subroutine destroy_bz_mesh_type(Kmesh)
  use defs_datatypes
  implicit none
  type(bz_mesh_type),intent(inout) :: Kmesh
 end subroutine destroy_bz_mesh_type
end interface

interface
 subroutine print_BZ_mesh(Kmesh,unit,prtvol,mode_paral)
  use defs_datatypes
  implicit none
  integer,intent(in),optional :: prtvol
  integer,intent(in),optional :: unit
  type(bz_mesh_type),intent(in) :: Kmesh
  character(len=4),intent(in),optional :: mode_paral
 end subroutine print_BZ_mesh
end interface

interface
 subroutine setup_k_rotation(nsym,symrec,timrev,nbz,kbz,krottb,krottbm1)
  use defs_basis
  implicit none
  integer,intent(in) :: nbz
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  real(dp),intent(in) :: kbz(3,nbz)
  integer,intent(inout) :: krottb(nbz,timrev,nsym)
  integer,intent(inout) :: krottbm1(nbz,timrev,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine setup_k_rotation
end interface

interface
 subroutine get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym,itim,ph_mkbzt)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(out) :: ik_ibz
  integer,intent(out) :: isym
  integer,intent(out) :: itim
  type(bz_mesh_type),intent(in) :: Kmesh
  complex(gwpc),intent(out),optional :: ph_mkbzt
  real(dp),intent(out) :: kbz(3)
 end subroutine get_BZ_item
end interface

interface
 subroutine get_IBZ_item(Kmesh,ik_ibz,kibz,wtk)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ik_ibz
  type(bz_mesh_type),intent(in) :: Kmesh
  real(dp),intent(out) :: wtk
  real(dp),intent(out) :: kibz(3)
 end subroutine get_IBZ_item
end interface

interface
 subroutine get_BZ_diff(Kmesh,k1,k2,idiff_bz,G0,nfound)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: idiff_bz
  integer,intent(out) :: nfound
  type(bz_mesh_type),intent(in) :: Kmesh
  integer,intent(out) :: G0(3)
  real(dp),intent(in) :: k1(3)
  real(dp),intent(in) :: k2(3)
 end subroutine get_BZ_diff
end interface

interface
 logical function is_samek(k1,k2,G0)
 use defs_basis
 implicit none
 integer,intent(out) :: G0(3)
 real(dp),intent(in) :: k1(3)
 real(dp),intent(in) :: k2(3)
end function is_samek
end interface

interface
 logical function has_BZ_item(Kmesh,item,ikbz,G0)
 use defs_basis
 use defs_datatypes
 implicit none
 integer,intent(out) :: ikbz
 type(bz_mesh_type),intent(in) :: Kmesh
 integer,intent(out) :: G0(3)
 real(dp),intent(in) :: item(3)
end function has_BZ_item
end interface

interface
 logical function has_IBZ_item(Kmesh,item,ikibz,G0)
 use defs_basis
 use defs_datatypes
 implicit none
 integer,intent(out) :: ikibz
 type(bz_mesh_type),intent(in) :: Kmesh
 integer,intent(out) :: G0(3)
 real(dp),intent(in) :: item(3)
end function has_IBZ_item
end interface

interface
 subroutine setup_little_group(ext_pt,Kmesh,gmet,Cryst,npwvec,gvec,npwe,use_umklp,prtvol,Ltg)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: npwe
  integer,intent(in) :: npwvec
  integer,intent(in) :: prtvol
  integer,intent(in) :: use_umklp
  type(crystal_structure),intent(in) :: Cryst
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(inout) :: Ltg
  real(dp),intent(in) :: ext_pt(3)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: gvec(3,npwvec)
 end subroutine setup_little_group
end interface

interface
 subroutine nullify_little_group(Ltg)
  use defs_datatypes
  implicit none
  type(little_group),intent(inout) :: Ltg
 end subroutine nullify_little_group
end interface

interface
 subroutine destroy_little_group(Ltg)
  use defs_datatypes
  implicit none
  type(little_group),intent(inout) :: Ltg
 end subroutine destroy_little_group
end interface

interface
 subroutine print_little_group(Ltg,unit,prtvol,mode_paral)
  use defs_datatypes
  implicit none
  integer,intent(in),optional :: prtvol
  integer,intent(in),optional :: unit
  type(little_group),intent(in) :: Ltg
  character(len=4),intent(in),optional :: mode_paral
 end subroutine print_little_group
end interface

interface
 subroutine setup_ppmodel(PPm,paral_kgb,Qmesh,Sp,Er,MPI_enreg,&  
  &  nfftf,gvec,ngfftf,gmet,gprimd,rhor_tot,&  
  &  epsm1q,iqiA) !Optional
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in),optional :: iqiA
  integer,intent(in) :: nfftf
  integer,intent(in) :: paral_kgb
  type(epsilonm1_results),intent(inout) :: Er
  type(mpi_type),intent(inout) :: MPI_enreg
  type(ppmodel_type),intent(inout) :: PPm
  type(bz_mesh_type),intent(in) :: Qmesh
  type(sigma_parameters),intent(in) :: Sp
  integer,intent(in) :: ngfftf(18)
  complex(gwpc),intent(in),optional :: epsm1q(Sp%npwc,Sp%npwc,Er%nomega)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,Sp%npwc)
  real(dp),intent(inout) :: rhor_tot(nfftf)
 end subroutine setup_ppmodel
end interface

interface
 subroutine setup_Qmesh(nqibz,Cryst,prtvol,qibz,Qmesh)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nqibz
  integer,intent(in) :: prtvol
  type(crystal_structure),intent(in) :: Cryst
  type(bz_mesh_type),intent(inout) :: Qmesh
  real(dp),intent(in) :: qibz(3,nqibz)
 end subroutine setup_Qmesh
end interface

interface
 subroutine find_Qmesh(Cryst,gprimd,Kmesh,nqsm,qsmall,Qmesh,prtvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nqsm
  integer,intent(in) :: prtvol
  type(crystal_structure),intent(in) :: Cryst
  type(bz_mesh_type),intent(in) :: Kmesh
  type(bz_mesh_type),intent(out) :: Qmesh
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: qsmall(3,nqsm)
 end subroutine find_Qmesh
end interface

interface
 subroutine setup_sigma(acell,rprim,ngfftf,Dtset,Dtfil,MPI_enreg,mpsang_kss,ngfft_gw,Hdr_kss,&  
  &  Cryst,Kmesh,Qmesh,Gsph_Max,Vcp,Er,Sp)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: mpsang_kss
  type(crystal_structure),intent(out) :: Cryst
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(inout) :: Dtset
  type(epsilonm1_results),intent(inout) :: Er
  type(gvectors_type),intent(out) :: Gsph_Max
  type(hdr_type),intent(out) :: Hdr_kss
  type(bz_mesh_type),intent(out) :: Kmesh
  type(mpi_type),intent(inout) :: MPI_enreg
  type(bz_mesh_type),intent(out) :: Qmesh
  type(sigma_parameters),intent(inout) :: Sp
  type(coulombian_type),intent(out) :: Vcp
  integer,intent(out) :: ngfft_gw(18)
  integer,intent(in) :: ngfftf(18)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
 end subroutine setup_sigma
end interface

interface
 subroutine sizefft(m,n)
  implicit none
  integer,intent(in) :: m
  integer,intent(out) :: n
 end subroutine sizefft
end interface

interface
 subroutine solve_Dyson(ikcalc,nomega_sigc,Sp,Kmesh,sigxme_tmp,sigcme_tmp,en_qp,Sr)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ikcalc
  integer,intent(in) :: nomega_sigc
  type(bz_mesh_type),intent(in) :: Kmesh
  type(sigma_parameters),intent(in) :: Sp
  type(sigma_results),intent(inout) :: Sr
  real(dp),intent(in) :: en_qp(Kmesh%nibz,Sp%nbnds,Sp%nsppol)
  complex(gwpc),intent(in) :: sigcme_tmp(nomega_sigc,Sp%minbnd(ikcalc):Sp%maxbnd(ikcalc), &
  &         Sp%minbnd(ikcalc):Sp%maxbnd(ikcalc),Sp%nsppol*Sp%nsig_ab)
  complex(gwpc),intent(in) :: sigxme_tmp(Sp%minbnd(ikcalc):Sp%maxbnd(ikcalc), &
  &         Sp%minbnd(ikcalc):Sp%maxbnd(ikcalc),Sp%nsppol*Sp%nsig_ab)
 end subroutine solve_Dyson
end interface

interface
 subroutine approxdelta(nomegasf,omegasf,egwdiff_re,smear,iomegal,iomegar,wl,wr,spmeth)
  use defs_basis
  implicit none
  integer,intent(out) :: iomegal
  integer,intent(out) :: iomegar
  integer,intent(in) :: nomegasf
  integer,intent(in) :: spmeth
  real(dp),intent(in) :: egwdiff_re
  real(dp),intent(in) :: smear
  real(dp),intent(out) :: wl
  real(dp),intent(out) :: wr
  real(dp),intent(in) :: omegasf(nomegasf)
 end subroutine approxdelta
end interface

interface
 subroutine calc_kkweight(ne,omegae,nsp,omegasp,delta,omegamax,kkw)
  use defs_basis
  implicit none
  integer,intent(in) :: ne
  integer,intent(in) :: nsp
  real(dp),intent(in) :: delta
  real(dp),intent(in) :: omegamax
  complex(dpc),intent(out) :: kkw(ne,nsp)
  complex(gwpc),intent(in) :: omegae(ne)
  real(dp),intent(in) :: omegasp(nsp)
 end subroutine calc_kkweight
end interface

interface
 subroutine setup_spectral(nomega,omega,nomegasf,omegasf,max_rest,min_rest,my_max_rest,my_min_rest,&  
  &  method,zcut,omegaplasma,my_wl,my_wr,kkweight)
  use defs_basis
  implicit none
  integer,intent(in) :: method
  integer,intent(out) :: my_wl
  integer,intent(out) :: my_wr
  integer,intent(in) :: nomega
  integer,intent(in) :: nomegasf
  real(dp),intent(in) :: max_rest
  real(dp),intent(in) :: min_rest
  real(dp),intent(in) :: my_max_rest
  real(dp),intent(in) :: my_min_rest
  real(dp),intent(in) :: omegaplasma
  real(dp),intent(in) :: zcut
  complex(dpc),intent(out) :: kkweight(nomega,nomegasf)
  complex(gwpc),intent(in) :: omega(nomega)
  real(dp),intent(out) :: omegasf(nomegasf)
 end subroutine setup_spectral
end interface

interface
 subroutine split_sigc(sp,sr,npwc1,npwc2,jb,isppol,io,ioe0j,theta_mu_minus_e0i,&  
  &  rhotwg,omegame0i,otq,botsq,sigccoh,sigcsex)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: io
  integer,intent(in) :: ioe0j
  integer,intent(in) :: isppol
  integer,intent(in) :: jb
  integer,intent(in) :: npwc1
  integer,intent(in) :: npwc2
  type(sigma_parameters),intent(in) :: sp
  type(sigma_results),intent(in) :: sr
  real(dp),intent(in) :: theta_mu_minus_e0i
  complex(gwpc) :: botsq(sp%npwc,npwc1)
  real(dp),intent(in) :: omegame0i(sp%nomegasr+sp%nomegasrd)
  real(dp),intent(in) :: otq(sp%npwc,npwc2)
  complex(gwpc),intent(in) :: rhotwg(sp%npwx)
  complex(gwpc),intent(inout) :: sigccoh(sp%nbnds,sp%nsppol)
  complex(gwpc),intent(inout) :: sigcsex(sp%nbnds,sp%nsppol)
 end subroutine split_sigc
end interface

interface
 subroutine symf12(npw,nomega,omega,eps,Gsphere,Gpairs_q,nsym,symrec,phgt,mode)
  use defs_datatypes
  implicit none
  integer,intent(in) :: mode
  integer,intent(in) :: nomega
  integer,intent(in) :: npw
  integer,intent(in) :: nsym
  type(gpairs_type),intent(in) :: Gpairs_q
  type(gvectors_type),intent(in) :: Gsphere
  complex,intent(inout) :: eps(npw,npw,nomega)
  complex,intent(in) :: omega(nomega)
  complex,intent(in) :: phgt(npw,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine symf12
end interface

interface
 subroutine symmetrize_afm_chi0(Cryst,Gsph,npwe,nomega,chi0)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nomega
  integer,intent(in) :: npwe
  type(crystal_structure),intent(in) :: Cryst
  type(gvectors_type),intent(in) :: Gsph
  complex(gwpc),intent(inout) :: chi0(npwe,npwe,nomega)
 end subroutine symmetrize_afm_chi0
end interface

interface
 subroutine testscr(optfil,unt,fname,nqibzR,nqlwlR,nomegaR,npweR,npwwfn_used,nbnds_used,title,&  
  &  fform,MPI_enreg,localrdwf,qibz_p,qlwl_p,omega_p,gvec_p,Hdr_scr)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: fform
  integer,intent(in) :: localrdwf
  integer,intent(out) :: nbnds_used
  integer,intent(out) :: nomegaR
  integer,intent(out) :: npweR
  integer,intent(out) :: npwwfn_used
  integer,intent(out) :: nqibzR
  integer,intent(out) :: nqlwlR
  integer,intent(in) :: optfil
  integer,intent(in) :: unt
  type(hdr_type),intent(out) :: Hdr_scr
  type(mpi_type),intent(in) :: MPI_enreg
  character(len=fnlen),intent(in) :: fname
  integer,pointer :: gvec_p(:,:)
  character(len=80),intent(out) :: title(2)
  complex(dpc),pointer :: omega_p(:)
  real(dp),pointer :: qibz_p(:,:)
  real(dp),pointer :: qlwl_p(:,:)
 end subroutine testscr
end interface

interface
 subroutine init_wf_info_1(Wf_info,gwmem,paral_kgb,npwwfn,my_minb,my_maxb,nk,nsppol,nspden,nspinor)
  use defs_datatypes
  implicit none
  integer,intent(in) :: gwmem
  integer,intent(in) :: my_maxb
  integer,intent(in) :: my_minb
  integer,intent(in) :: nk
  integer,intent(in) :: npwwfn
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: paral_kgb
  type(wavefunctions_information),intent(inout) :: Wf_info
 end subroutine init_wf_info_1
end interface

interface
 subroutine init_wf_info_2(Wf_info,igfft0,ngfft)
  use defs_datatypes
  implicit none
  type(wavefunctions_information),intent(inout) :: Wf_info
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: igfft0(Wf_info%npwwfn)
 end subroutine init_wf_info_2
end interface

interface
 subroutine destroy_wf_info(Wf_info)
  use defs_datatypes
  implicit none
  type(wavefunctions_information),intent(inout) :: Wf_info
 end subroutine destroy_wf_info
end interface

interface
 subroutine reinit_wf_info(Wf_info)
  use defs_datatypes
  implicit none
  type(wavefunctions_information),intent(inout) :: Wf_info
 end subroutine reinit_wf_info
end interface

interface
 subroutine get_wfr(Wf_info,MPI_enreg,ib,ik,is,wfr)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ib
  integer,intent(in) :: ik
  integer,intent(in) :: is
  type(mpi_type),intent(inout) :: MPI_enreg
  type(wavefunctions_information),intent(inout) :: Wf_info
  complex(gwpc),intent(out) :: wfr(Wf_info%nfft*Wf_info%nspinor)
 end subroutine get_wfr
end interface

interface
 subroutine duplicate_wf_info(MPI_enreg,Wf_info,Wf_info_braket,kcalc,Kmesh)
  use defs_datatypes
  implicit none
  type(bz_mesh_type),intent(in) :: Kmesh
  type(mpi_type),intent(in) :: MPI_enreg
  type(wavefunctions_information),intent(inout) :: Wf_info
  type(wavefunctions_information),intent(inout) :: Wf_info_braket
  integer,intent(in) :: kcalc(Wf_info_braket%nk)
 end subroutine duplicate_wf_info
end interface

interface
 subroutine nullify_wf_info(Wf_info)
  use defs_datatypes
  implicit none
  type(wavefunctions_information),intent(inout) :: Wf_info
 end subroutine nullify_wf_info
end interface

interface
 subroutine print_wavefunctions_information(Wf_info,unitno,prtvol,mode_paral)
  use defs_datatypes
  implicit none
  integer,optional,intent(in) :: prtvol
  integer,optional,intent(in) :: unitno
  type(wavefunctions_information),intent(in) :: Wf_info
  character(len=4),optional,intent(in) :: mode_paral
 end subroutine print_wavefunctions_information
end interface

interface
 subroutine rotate_wfg(Wfs,Kmesh,iband,ik_bz,isppol,grottbm1,phmGt,wfg_rot)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iband
  integer,intent(in) :: ik_bz
  integer,intent(in) :: isppol
  type(bz_mesh_type),intent(in) :: Kmesh
  type(wavefunctions_information),intent(in) :: Wfs
  integer,target,intent(in) :: grottbm1(:,:,:)
  complex(gwpc),intent(in) :: phmGt(:,:)
  complex(gwpc),intent(out) :: wfg_rot(Wfs%npwwfn)
 end subroutine rotate_wfg
end interface

interface
 subroutine rotate_wfr(Wfs,Kmesh,iband,ik_bz,isppol,irottb,MPI_enreg,ur_rot)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iband
  integer,intent(in) :: ik_bz
  integer,intent(in) :: isppol
  type(bz_mesh_type),intent(in) :: Kmesh
  type(mpi_type),intent(inout) :: MPI_enreg
  type(wavefunctions_information),intent(inout) :: Wfs
  integer,intent(in),target :: irottb(Wfs%nfft,Kmesh%nsym)
  complex(gwpc),intent(out) :: ur_rot(Wfs%nfft)
 end subroutine rotate_wfr
end interface

interface
 subroutine wrqps(Dtfil,Sp,Kmesh,nspden,nscf,nfftot,Sr,m_lda_to_qp,rho_qp)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nfftot
  integer,intent(in) :: nscf
  integer,intent(in) :: nspden
  type(datafiles_type),intent(in) :: Dtfil
  type(bz_mesh_type),intent(in) :: Kmesh
  type(sigma_parameters),intent(in) :: Sp
  type(sigma_results),intent(in) :: Sr
  complex(dpc),intent(in) :: m_lda_to_qp(Sp%nbnds,Sp%nbnds,Kmesh%nibz,Sp%nsppol)
  real(dp),intent(in) :: rho_qp(nfftot,nspden)
 end subroutine wrqps
end interface

interface
 subroutine wrscr(iq,optfil,unt,fname,Hdr,Dtset,npwe,npwwfn_used,nbnds_used,&  
  &  nqibz,nqlwl,nomega,qibz,omega,gvec,gmet,epsm1,title,&  
  &  qlwl,uwing,lwing)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iq
  integer,intent(in) :: nbnds_used
  integer,intent(in) :: nomega
  integer,intent(in) :: npwe
  integer,intent(in) :: npwwfn_used
  integer,intent(in) :: nqibz
  integer,intent(in) :: nqlwl
  integer,intent(in) :: optfil
  integer,intent(in) :: unt
  type(dataset_type),intent(in) :: Dtset
  type(hdr_type),intent(inout) :: Hdr
  character(len=fnlen),intent(in) :: fname
  character(len=80),intent(in) :: title(2)
  complex(gwpc),intent(in) :: epsm1(npwe,npwe,nomega)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: gvec(3,npwe)
  complex(gwpc),intent(in) :: lwing(npwe,nomega,nqlwl*optfil)
  complex(gwpc),intent(in) :: omega(nomega)
  real(dp),intent(in) :: qibz(3,nqibz)
  real(dp),intent(in) :: qlwl(3,nqlwl*optfil)
  complex(gwpc),intent(in) :: uwing(npwe,nomega,nqlwl*optfil)
 end subroutine wrscr
end interface

end module interfaces_15gw
!!***
