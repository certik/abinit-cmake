!{\src2tex{textfont=tt}}
!!****f* ABINIT/sigma
!! NAME
!! sigma
!!
!! FUNCTION
!! Calculate the matrix elements self-energy operator.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MT, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! codvsn=code version
!! Dtfil<type(datafiles_type)>=variables related to files
!! Dtset<type(dataset_type)>=all input variables for this dataset
!! iexit=exit flag
!! MPI_enreg=information about MPI parallelization
!! Pawang<type(pawang_type)>=paw angular mesh and related data
!! Pawrad(ntypat*usepaw)<type(pawrad_type)>=paw radial mesh and related data
!! Pawtab(ntypat*usepaw)<type(pawtab_type)>=paw tabulated starting data
!! Psps<type(pseudopotential_type)>=variables related to pseudopotentials
!!   Before entering the first time in sigma, a significant part of Psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,ntypat,n1xccc,usepaw,useylm, 
!!   and the arrays dimensioned to npsp. All the remaining components of Psps are to be initialized in 
!!   the call to pspini. The next time the code enters screening, Psps might be identical to the
!!   one of the previous Dtset, in which case, no reinitialisation is scheduled in pspini.F90.
!! rprim(3,3)=dimensionless real space primitive translations
!! xred(3,natom) = reduced atomic coordinates
!!
!! OUTPUT
!!  Output is written on the main abinit output file. Some results are stored in external files
!!
!! PARENTS
!!      driver
!!
!! NOTES
!!
!! ON THE USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut) for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ... are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg) for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ... Total density, potentials, ... are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used. It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf) are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! CHILDREN
!!      abi_etsf_electrons_put,abi_etsf_geo_put,assert,calc_density,calc_ffm
!!      calc_vhxc_braket,calc_wf_qp,chkpawovlp,cigfft,cprj_alloc,cprj_free
!!      csigme,cutoff_density,destroy_bz_mesh_type,destroy_coulombian
!!      destroy_crystal_structure,destroy_epsilonm1_parameters
!!      destroy_epsilonm1_results,destroy_gvectors_type,destroy_little_group
!!      destroy_paw_an_type,destroy_paw_ij_type,destroy_pawfgrtab
!!      destroy_ppmodel,destroy_sigma_parameters,destroy_sigma_results
!!      destroy_wf_info,duplicate_wf_info,energies_init,eps1_tc,fappnd,fermi
!!      findk,findqg0,flush_unit,fourdp,get_ng0sh,getph,gw2abi,hdr_clean
!!      hdr_vs_dtset,init_crystal_from_hdr,init_gvectors_type,init_pawfgr
!!      init_pawfgrtab,init_ppmodel,init_sigma_results,init_wf_info_1
!!      init_wf_info_2,initmpi_seq,initylmg,int2char4,ioarr,kgindex,memerr
!!      metric,mkrdim,nhatgrid,nullify_crystal_structure
!!      nullify_epsilonm1_parameters,nullify_epsilonm1_results
!!      nullify_little_group,nullify_pawfgrtab,nullify_ppmodel
!!      nullify_sigma_parameters,nullify_sigma_results,nullify_wf_info
!!      paw_mkrhox,paw_mkrhox_spl,paw_symcprj,pawdenpot,pawdij,pawinit
!!      pawmknhat,pawmkrhoij,pawprt,pawpuxinit,pclock,print_crystal_structure
!!      print_ngfft,print_pawtab,print_psps,prtrhomxmn,pspini,rdgw,rdkss,rdqps
!!      rdscr,reinit_wf_info,rhoij_alloc,rhoij_copy,rhoij_free,setmesh
!!      setshells,setsymrhoij,setup_coulombian,setup_fft_rotation,setup_kmesh
!!      setup_little_group,setup_ppmodel,setup_qmesh,setvtr,split_work2,status
!!      symdij,symrhoij,test_charge,testlda,testscr,timab,update_cprj
!!      write_sigma_results,write_sigma_results_header,wrqps,wrtout,xcomm_init
!!      xmaster_init,xme_init,xproc_max
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine sigma(acell,codvsn,Dtfil,Dtset,iexit,MPI_enreg,Pawang,Pawrad,Pawtab,Psps,rprim,xred)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : GW_TOLQ0, unt_gw, unt_sig, unt_sgr, unt_sgm
 use m_errors, only : assert
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_13ionetcdf
 use interfaces_13nonlocal
 use interfaces_13paw
 use interfaces_13psp
 use interfaces_13recipspace
 use interfaces_14iowfdenpot
 use interfaces_15common
 use interfaces_15gw
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: iexit
 character(len=6),intent(in) :: codvsn
 type(Datafiles_type),intent(in) :: Dtfil
 type(Dataset_type),intent(inout) :: Dtset
 type(MPI_type),intent(inout) :: MPI_enreg
 type(Pawang_type),intent(inout) :: Pawang
 type(Pseudopotential_type),intent(inout) :: Psps
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3),xred(3,Dtset%natom)
 type(Pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Psps%usepaw)
 type(Pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Psps%usepaw)

!Local variables-------------------------------
 character(len=50),parameter :: FILE__='sigma.F90'
!scalars
 integer,parameter :: level=25,tim_fourdp=5
 integer,save :: nsym_old=-1
 integer :: accessfil,b1gw,b2gw,bantot,choice,cplex,cplex_dij,dim1_rhox,dim2_rhox
 integer :: enforce_sym,fformr,has_dijso,has_dijU,iab,iapp,iat,ib,ib1,ib2,ider
 integer :: idx_dwn,idx_up,ierr,ifft,ig,ii,ik,ik_bz,ik_ibz,ikbs_proc_rank
 integer :: ikcalc,ikibz,ikq,ilmn,indx,initialized,io,iq,iqm,irank,is,is_idx
 integer :: ispden,isppol,istat,istep,isym,itypat,izero,jj,jlmn,k0lmn,klmn
 integer :: lm_size,lmn2_size,lpawumax,master,mband_,method,mgfft,mgfftf
 integer :: mgfftgw,mkmem_,mod10,moved_atm_inside,moved_rhor,mpsang,mpw_
 integer :: mqmem__,my_maxb,my_minb,my_nbnds,n3xccc,nG01d,nG02d,nG03d,nbcw
 integer :: nbsc,nbvw,nbw,ndij,nel=0,nfftc_tot,nfftf,nfftf_tot,nfftgw
 integer :: nfftgw_tot,nhatgrdim,nkxc,nprocs,nscf
 integer :: ntasks,nzlmopt,optder,optene,optgr0,optgr1,optgr2,option
 integer :: optrad,optrhoij,psp_gencond,rank,rdwr,rdwrpaw,restart,restartpaw
 integer :: rhoxsp_method,shift,spaceComm,two_lmaxp1,umklp_opt
 integer :: use_umklp,usexcnhat
 real(dp) :: band_ene,boxcut,boxcutc,compch_fft,compch_sph,diecut_eff_dum
 real(dp) :: drude_plsmf,dummy,ecore,ecut_eff,ecutdg_eff
 real(dp) :: efermi,efermi_qp,ehartree,exchange_energy,g1,g2,g3,gsq,gsqcutc_eff
 real(dp) :: gsqcutf_eff,gsqcutf_eff_tmp,norm,oldefermi,sij,ucvol,vxcavg
 real(dp) :: vxcavg_qp
 logical,parameter :: prior_to_56=.TRUE.,update_energies=.FALSE.
 logical :: ltest
 character(len=4) :: tag
 character(len=500) :: msg
 character(len=fnlen) :: filapp,fname
 type(Bands_Symmetries) :: BSym
 type(Bandstructure_type) :: Bstruct
 type(BZ_mesh_type) :: Kmesh,Qmesh
 type(Coulombian_type) :: Vcp
 type(Crystal_structure) :: Cryst
 type(Energies_type) :: KS_energies,QP_energies
 type(Epsilonm1_results) :: Er
 type(Gvectors_type) :: Gsph_Max
 type(Hdr_type) :: Hdr_kss,Hdr_sigma
 type(MPI_type) :: MPI_enreg_seq
 type(PPmodel_type) :: PPm
 type(Pawfgr_type) :: Pawfgr
 type(Sigma_parameters) :: Sp
 type(Sigma_results) :: Sr
 type(Wavefunctions_information) :: Wf_info,Wf_info_braket
!arrays
 integer :: g0(3),ibocc(Dtset%nsppol),ngfft_gw(18),ngfftc(18)
 integer :: ngfftf(18)
 integer,allocatable :: dimlmn(:),gvec(:,:),igfftf(:),irottb(:,:),irottbf(:,:)
 integer,allocatable :: istart(:),istop(:),istwfk_(:),ktabr(:,:)
 integer,allocatable :: l_size_atm(:),nband_(:),nbv(:),nlmn(:),npwarr(:)
 integer,allocatable :: qg(:,:)
 integer,allocatable,target :: igfft(:,:,:,:)
 integer,pointer :: igfft0(:),trec(:,:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),k0(3),kgwmk(3),rmet(3,3),rprimd(3,3)
 real(dp) :: strsxc(6),tsec(2)
 real(dp),allocatable :: abiocc_qp(:),doccde(:),eigen(:),grewtn(:,:)
 real(dp),allocatable :: ks_energy(:,:,:),ks_occ(:,:,:),kxc(:,:),kxc_qp(:,:)
 real(dp),allocatable :: nhat(:,:),nhat_qp(:,:),nhatgr(:,:,:),nhatgr_qp(:,:,:)
 real(dp),allocatable :: occfact(:),pawrhox(:,:,:,:,:),pawrhox_spl(:,:,:,:,:)
 real(dp),allocatable :: ph1d(:,:),ph1df(:,:),qp_energy(:,:,:)
 real(dp),allocatable :: qp_occ(:,:,:),qtmp(:,:),rhog(:,:),rhog_qp(:,:)
 real(dp),allocatable :: rhor(:,:),rhor_p(:,:),rhor_qp(:,:),sr_gwenergy(:,:,:)
 real(dp),allocatable :: vhartr_qp_vtr(:),vhartr_vtr(:),vkb(:,:,:,:)
 real(dp),allocatable :: vkbd(:,:,:,:),vkbsign(:,:),vpsp(:),vtrial(:,:)
 real(dp),allocatable :: vtrial_qp(:,:),vxc_qp_vtr(:,:),vxc_vtr(:,:),xccc3d(:)
 real(dp),allocatable :: xred_dummy(:,:),ylm_q(:,:),ylmgr_q(:,:,:)
 real(dp),pointer :: ttns(:,:)
 complex(dpc) :: ovlp(2)
 complex(dpc),allocatable :: ctmp(:,:),hbare(:,:,:,:),hlda(:,:,:,:),htmp(:,:,:,:)
 complex(dpc),allocatable :: m_lda_to_qp(:,:,:,:),uks2qp(:,:),vUpaw(:,:,:,:),vUpaw_qp(:,:,:,:)
 complex(dpc),allocatable :: vhartr(:,:,:,:),vxc(:,:,:,:),vxc_qp(:,:,:,:)
 complex(dpc),allocatable :: vxc_val(:,:,:,:),vxcval_qp(:,:,:,:)
 complex(gwpc),allocatable :: kxcg(:,:)
 complex(gwpc),pointer :: cgdwn(:),cgup(:) 
 logical,allocatable :: mask(:)
 character(len=80) :: title(2)
 character(len=fnlen) :: tmpfil(7)
 type(Cprj_type),allocatable :: Cprj_bz(:,:),Cprj_ibz(:,:)
 type(Little_group),allocatable :: Ltg_k(:)
 type(Paw_an_type),allocatable :: KS_paw_an(:),QP_paw_an(:)
 type(Paw_ij_type),allocatable :: KS_paw_ij(:),QP_paw_ij(:)
 type(Pawfgrtab_type),allocatable :: Pawfgrtab(:)
 type(Pawrhoij_type),allocatable :: Pawrhoij(:),Pawrhoij_dum(:),QP_pawrhoij(:)

!************************************************************************

#ifdef FC_PGI6
 write(msg,'(6a)')ch10,&
& ' sigma.F90: COMMENT - ',ch10,&
& ' Due to a bug in PGI v6, the compilation of sigma.F90 has been skipped ',ch10,&
& ' To perform GW calculations use a more recent version of PGI, alternatively you might try a different compiler ' 
 call wrtout(std_out,msg,'COLL') 
 call wrtout(ab_out,msg,'COLL') 
 call leave_new('COLL')
#else

#if defined DEBUG_MODE
 write(msg,'(a)')' sigma : enter '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 call pclock(0)
 call timab(401,1,tsec) ! overall time
 call timab(402,1,tsec) ! sigma(1)

 write(msg,'(8a)')&
& ' SIGMA: Calculation of the GW corrections ',ch10,ch10,&
& ' Based on a program developped by R.W. Godby, V. Olevano, G. Onida, and L. Reining.',ch10,&
& ' Incorporated in ABINIT by V. Olevano, G.-M. Rignanese, and M. Torrent.',ch10
 call wrtout(std_out,msg,'COLL') 
 call wrtout(ab_out,msg,'COLL')

#if defined HAVE_GW_DPC
 write(msg,'(a,i2)')'.Using double precision arithmetic ; gwpc = ',gwpc
 ltest=(gwpc==8) 
 write(msg,'(6a)')ch10,&
& ' Number of bytes for double precision complex /=8 ',ch10,&
& ' Cannot continue due to kind mismatch in BLAS library ',ch10,&
& ' Some BLAS interfaces are not generated by abilint '
 call assert(ltest,msg,__FILE__,__LINE__)
#else 
 write(msg,'(a,i2)')'.Using single precision arithmetic ; gwpc = ',gwpc
#endif
 call wrtout(std_out,msg,'COLL') 
 !call wrtout(ab_out,msg,'COLL') 
 !
 ! === Initialize MPI related variables and parallelization level ===
 ! gwpara 0 => sequential run
 !        1 => parallelism is over k-points (memory is not parallelized)
 !        2 => parallelism is over bands: bands are divided among processors,
 !             but each proc has all the states where GW corrections are required
 call xcomm_init  (MPI_enreg,spaceComm)  
 call xme_init    (MPI_enreg,rank     )          
 call xmaster_init(MPI_enreg,master   )  
 call xproc_max   (nprocs,ierr)

 if (nprocs==1) Dtset%gwpara=0  
 MPI_enreg%gwpara      = Dtset%gwpara           
 MPI_enreg%parareel    = 0  
 MPI_enreg%paralbd     = 0
 MPI_enreg%paral_level = 2   ! This means k-points but it is not used
 MPI_enreg%me_fft      = 0 
 MPI_enreg%nproc_fft   = 1
 ! * Fake MPI_type for sequential part
 call initmpi_seq(MPI_enreg_seq) ; MPI_enreg_seq%nproc_fft=1 ; MPI_enreg_seq%me_fft=0
!
!Create names for the temporary files based on Dtfil%filnam_ds(5) by appending adequate string.
!'_WF1' -> Dtfil%unwft1
!'_WF2' -> Dtfil%unwft2
!'_KG'  -> Dtfil%unkg
!'_DUM' -> tmp_unit (real dummy name)
!'_YLM' -> Dtfil%unylm
!'_PAW' -> Dtfil%unpaw
!
!Concerning the IO mode:
!localrdwf==1 ==> every processor has access to files (default)
!localrdwf==0 ==> only master has access to files
!
!if accesswff == 2 then set all outputs to netcdf format
!if accesswff == 3 then set all outputs to ETSF format
!
 accessfil=0
 if (Dtset%accesswff==2) accessfil=1
 if (Dtset%accesswff==3) accessfil=3
 if (Dtset%accesswff==1) accessfil=4

!Prepare the name of the auxiliary files DEN, DOS, EIG...
 iapp = 0
 call fappnd(filapp,dtfil%filnam_ds(4),iapp)

 tmpfil(1)=TRIM(Dtfil%filnam_ds(5))//'_WF1'
 tmpfil(2)=TRIM(Dtfil%filnam_ds(5))//'_WF2'
 tmpfil(3)=TRIM(Dtfil%filnam_ds(5))//'_KG'
 tmpfil(4)=TRIM(Dtfil%filnam_ds(5))//'_DUM'
 tmpfil(6)=TRIM(Dtfil%filnam_ds(5))//'_YLM'
 tmpfil(7)=TRIM(Dtfil%filnam_ds(5))//'_PAW'
!
!* Parallel case: the index of the processor must be appended
 if (MPI_enreg%paral_compil_kpt==1) then
  call int2char4(MPI_enreg%me,tag) 
  jj=1 ; if (MPI_enreg%paral_compil_mpio==1 .and. Dtset%accesswff==1) jj=3
  do ii=jj,7 
   tmpfil(ii)=TRIM(tmpfil(ii))//'_P-'//tag 
  end do
 end if
!
!=== Some variables need to be initialized/nullify at start ===

!
!Initialize the band structure datatype.
!£ this is for the next version, should fix problems with nkpt and symmorphy
!£ bantot=0
!£ do isppol=1,Dtset%nsppol
!£  do ikibz=1,Dtset%nkpt
!£   if (Dtset%nband(1)/=Dtset%nband(ikibz+(isppol-1)*Dtset%nkpt)) stop "nband must be costant"
!£   bantot=bantot+Dtset%nband(ikibz+(isppol-1)*Dtset%nkpt)
!£  end do
!£ end do
!£ allocate(doccde(bantot),eigen(bantot),npwarr(Dtset%nkpt),occfact(bantot))
!£ doccde(:)=zero ; eigen(:)=zero ; npwarr(:)=0 ; occfact(:)=zero 
!£ !occfact(:)=Dtset%occ_orig(1:Dtset%mband*Dtset%nkpt*Dtset%nsppol)
!£
!£ call bstruct_init(bantot,Bstruct,doccde,eigen,Dtset%istwfk,Dtset%kptns,&
!£& Dtset%nband,Dtset%nkpt,npwarr,Dtset%nsppol,occfact,Dtset%wtk) !should check wtk wrt Kmesh%wt
!£ deallocate(doccde,eigen,npwarr,occfact)

 call energies_init(KS_energies)
 usexcnhat=0

!=== Dimensional primitive translations rprimd (from input), gprimd, metrics and unit cell volume ===
 call mkrdim(acell,rprim,rprimd)  
 call metric(gmet,gprimd,ab_out,rmet,rprimd,ucvol) 
!
!=== Define FFT grid(s) sizes === 
!* BE CAREFUL!, these meshes are used to evaluate the density and the matrix elements of v_Hxc  
!They are defined according to ecut and ecutdg. See also NOTES in the comments at the beginning of this file.
!* The FFT mesh used for the oscillator matrix elements is defined in setmesh.F90 and is usually much smaller
!TODO should add a comment ecutwfn<ecut
 k0(:)=zero
 call init_pawfgr(Dtset,k0,gmet,Pawfgr,mgfftf,nfftf,ecut_eff,ecutdg_eff,gsqcutc_eff,gsqcutf_eff,ngfftc,ngfftf)
 nfftc_tot=ngfftc(1)*ngfftc(2)*ngfftc(3)
 nfftf_tot=ngfftf(1)*ngfftf(2)*ngfftf(3)
 call print_ngfft(ngfftf,'FFT mesh for density and matrix elements of V_Hxc')

 ! === Set up parameters of the calculation and basic data structures ===
 call nullify_sigma_parameters(Sp)
 call nullify_epsilonm1_results(Er)

 call setup_sigma(acell,rprim,ngfftf,Dtset,Dtfil,MPI_enreg,mpsang,ngfft_gw,Hdr_kss,&
& Cryst,Kmesh,Qmesh,Gsph_Max,Vcp,Er,Sp)

 mod10=MOD(Sp%gwcalctyp,10)
 b1gw=Sp%minbdgw
 b2gw=Sp%maxbdgw
 nfftgw_tot=ngfft_gw(1)*ngfft_gw(2)*ngfft_gw(3)
 nfftgw=nfftgw_tot 
 mgfftgw=MAXVAL(ngfft_gw(1:3))

 if (ANY(ABS(Cryst%xred-xred)>tol6)) STOP 'BUG xred'
 if (ANY(ABS(rprimd-Cryst%rprimd)>tol6)) STOP 'BUG rprimd'
 if (Psps%ntypat/=Cryst%ntypat) STOP 'BUG ntypat'
 !
 ! === Read pseudopotential files ===
 call status(0,Dtfil%filstat,iexit,level,'call pspini   ')
 call pspini(Dtset,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,level,Pawrad,Pawtab,Psps,Cryst%rprimd)
!if (psp_gencond==1) 
 call print_psps(Psps,std_out,0,'COLL')
 ltest=(Psps%mpsang==mpsang) 
 call assert(ltest,'Psps%mpsang/=mpsang_kss',__FILE__,__LINE__)
 KS_energies%e_corepsp=ecore/Cryst%ucvol
!
!TRYING TO RECREATE AN "ABINIT ENVIRONMENT"
 allocate(nband_(Kmesh%nibz*Sp%nsppol)) ; nband_(:)=Sp%nbnds
 allocate(istwfk_(Kmesh%nibz)) ; istwfk_(:)=1
!
!Initialize MPI_enreg%proc_distrb according to gwpara. In case of parallelism 
!on k-points (gwpara==1) proc_distrb is redefined on the fly for each q-point.
 my_minb=1 ; my_maxb=Sp%nbnds ; my_nbnds=my_maxb-my_minb+1

 allocate(MPI_enreg%proc_distrb(Kmesh%nibz,Sp%nbnds,Sp%nsppol)) ; MPI_enreg%proc_distrb(:,:,:)=rank
 if (MPI_enreg%gwpara==2) then 
! === Setup distrb in case of band parallelism ===
  write(msg,'(2a)')ch10,' sigma : loop over bands done in parallel '
  call wrtout(std_out,msg,'PERS')
  allocate(istart(nprocs),istop(nprocs))
  MPI_enreg%proc_distrb(:,:,:)=-999
  call split_work2(Sp%nbnds,nprocs,istart,istop)
  my_minb=istart(rank+1) ; my_maxb=istop(rank+1)
  if (my_minb>my_maxb) then 
   write(msg,'(6a,2(i6,a),a)')ch10,&
&   ' sigma : ERROR - ',ch10,&
&   ' One or more processors has zero number of bands ',ch10,&
&   ' my_minb = ',my_minb,' my_maxb = ',my_maxb,ch10,&
&   ' This is a waste, decrease the number of processors '
   call wrtout(std_out,msg,'PERS') 
   call leave_new('COLL')
  end if 
  do irank=0,nprocs-1
   MPI_enreg%proc_distrb(:,istart(irank+1):istop(irank+1),:)=irank
  end do 
  deallocate(istart,istop)
  my_nbnds=my_maxb-my_minb+1
! * Announce the treatment of bands by each proc
  do irank=0,nprocs-1
   if (irank==rank) then 
    write(msg,'(4(a,i4))')&
&    ' treating ',my_nbnds,' bands from ',my_minb,' up to ',my_maxb,' by node ',irank
    call wrtout(std_out,msg,'PERS')
   end if 
  end do 
 end if 
!
!=== Allocate basic arrays and KS electronic structure ===
!TODO change ordering of indeces in ks_energy and ks_occ
 allocate(gvec(3,Sp%npwvec))
 allocate(ks_energy(Kmesh%nibz,Sp%nbnds,Sp%nsppol)) ;  ks_energy(:,:,:)=zero 
 allocate(ks_occ   (Kmesh%nibz,Sp%nbnds,Sp%nsppol)) ;  ks_occ(:,:,:)=zero 

 call nullify_wf_info(Wf_info)
 call init_wf_info_1(Wf_info,Dtset%gwmem,Dtset%paral_kgb,Sp%npwwfn,&
& my_minb,my_maxb,Kmesh%nibz,Sp%nsppol,Dtset%nspden,Dtset%nspinor)
 !
 ! ================================
 ! ==== Read KS band structure ====
 ! ================================ 
 call timab(403,1,tsec) ! rdkss
 allocate(trec(3,3,Cryst%nsym),ttns(3,Cryst%nsym))
 allocate(vkbsign(mpsang,Cryst%ntypat))
 allocate(vkb(Sp%npwwfn,Cryst%ntypat,mpsang,Kmesh%nibz),STAT=istat)
 if (istat/=0) call memerr(FILE__,'vkb',Sp%npwwfn*Cryst%ntypat*mpsang*Kmesh%nibz,'spc')
 allocate(vkbd(Sp%npwwfn,Cryst%ntypat,mpsang,Kmesh%nibz),STAT=istat)
 if (istat/=0) call memerr(FILE__,'vkbd',Sp%npwwfn*Cryst%ntypat*mpsang*Kmesh%nibz,'spc')

 allocate(Cprj_ibz(Dtset%natom,Dtset%nspinor*Kmesh%nibz*Sp%nbnds*Sp%nsppol*Dtset%usepaw)) 

 if (Dtset%usepaw==1) then 
  allocate(dimlmn(Cryst%natom),nlmn(Cryst%ntypat))
  do iat=1,Cryst%natom
   dimlmn(iat)=Pawtab(Cryst%typat(iat))%lmn_size
  end do
  do itypat=1,Cryst%ntypat
   nlmn(itypat)=Pawtab(itypat)%lmn_size
  end do
  call cprj_alloc(Cprj_ibz,0,dimlmn)
  allocate(Pawrhoij(Cryst%natom)) 
  call rhoij_alloc(1,nlmn,Dtset%nspden,Dtset%nsppol,Pawrhoij,Cryst%typat)
 end if

 nbvw=0 ! Not used 
 call rdkss(Dtfil,Dtset,Pawtab,Cryst%nsym,Sp%nbnds,nbvw,Kmesh%nibz,Sp%npwvec,Dtset%nspinor,Sp%nsppol,Sp%npwwfn,title,trec,&
& gvec,ks_energy,ks_occ,Wf_info%wfg,Cprj_ibz,Cryst%ntypat,Cryst%natom,mpsang,ttns,vkbsign,vkb,vkbd,nel,MPI_enreg,&
& my_minb,my_maxb)
 deallocate(vkbsign,vkb,vkbd)

 ltest=ALL(gvec(:,:)==Gsph_Max%gvec(:,1:Sp%npwvec))
 deallocate(gvec)
 call assert(ltest,'gvec/=Gsph_Max%gvec',__FILE__,__LINE__)
 ltest=(ALL(istwfk_==Hdr_kss%istwfk)) 
 call assert(ltest,'istwfk_ /= hdr%istwfk',__FILE__,__LINE__)
 ltest=ALL(trec==Cryst%symrec)      
 call assert(ltest,'BUG in trec',__FILE__,__LINE__)
 ltest=ALL((ttns-Cryst%tnons)<tol6) 
 call assert(ltest,'BUG in ttns',__FILE__,__LINE__)  
 deallocate(trec,ttns)
 call timab(403,2,tsec) ! rdkss
 ! 
 ! ============================
 ! ==== PAW initialization ====
 ! ============================
  if (Dtset%usepaw==1) then 

   call chkpawovlp(Cryst%natom,Cryst%ntypat,Dtset%pawovlp,Pawtab,Cryst%rmet,Cryst%typat,Cryst%xred)
!  
!  if (psp_gencond==1) then
   call timab(553,1,tsec)
!  1- Initialize values for several PAW arrays 
!  TODO this part should be done at the very beginning, otherwise test on orthonormalization in rdkss fails!
   diecut_eff_dum=ABS(Dtset%diecut)*Dtset%dilatmx**2
   call pawinit(diecut_eff_dum,Psps%indlmn,Dtset%pawlcutd,Dtset%pawlmix,Psps%lmnmax,Psps%mpsang,Psps%n1xccc,&
&   Dtset%pawnphi,Cryst%nsym,Dtset%pawntheta,Cryst%ntypat,Pawang,Pawrad,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev)
   call print_pawtab(Pawtab)
   call timab(553,2,tsec)
!  end if
!  FIXME Note that here nsym_old comes from KSS.
!  if (psp_gencond==1) then !.or. nsym_old/=nsym) then
   call setsymrhoij(gprimd,Pawang%l_max-1,Cryst%nsym,Dtset%pawprtvol,&
&   Cryst%rprimd,Dtset%symafm,Cryst%symrec,Pawang%zarot)
   nsym_old=Cryst%nsym
!  end if
!  * Initialize and compute data for LDA+U
   if (Dtset%usepawu>0.or.Dtset%useexexch>0) then
    call pawpuxinit(Dtset%dmatpuopt,Dtset%exchmix,Dtset%jpawu,Dtset%lexexch,Dtset%lpawu,&
&    Psps%indlmn,Psps%lmnmax,Cryst%ntypat,Pawang,Dtset%pawprtvol,Pawrad,Pawtab,Dtset%upawu,&
&    Dtset%useexexch,Dtset%usepawu)
   end if
!  3-
!  if (mkmem_==0) then
!  open(Dtfil%unpaw,file=tmpfil(7),form='unformatted',status='unknown')
!  rewind(unit=Dtfil%unpaw)
!  end if

!  === Get Pawrhoij from the header of the KSS file === 
   call rhoij_copy(Hdr_kss%pawrhoij,Pawrhoij)
!  === Re-symmetrize symrhoij ===
!  choice=1 ; optrhoij=1 
!  call symrhoij(choice,Psps%indlmn,Cryst%indsym,Psps%lmnmax,Cryst%natom,Cryst%nsym,Cryst%ntypat,&
!  &   optrhoij,Pawang,Dtset%pawprtvol,Pawrhoij,Dtset%symafm,Cryst%symrec,Cryst%typat)
!  
!  === Evaluate form factor of radial part of phi.phj-tphi.tphj === 
   rhoxsp_method=2 
   if (rhoxsp_method==1) then ! Arnaud-Alouani 
    dim1_rhox=2*(Psps%mpsang-1)  
    dim2_rhox=Psps%lnmax*(Psps%lnmax+1)/2  
   else if (rhoxsp_method==2) then ! Shiskin-Kresse
    dim1_rhox=MAXVAL(Pawtab(:)%l_size)**2
    dim2_rhox=MAXVAL(Pawtab(:)%lmn2_size) 
   end if
   allocate(pawrhox_spl(Psps%mqgrid_ff,2,0:dim1_rhox,dim2_rhox,Cryst%ntypat))  

   call paw_mkrhox_spl(Cryst%ntypat,Psps,Pawrad,Pawtab,pawrhox_spl,rhoxsp_method,dim1_rhox,dim2_rhox)
!  
!  === Variables/arrays related to the fine FFT grid ===
   allocate(nhat(nfftf,Dtset%nspden)) ; nhat(:,:)=zero
   allocate(Pawfgrtab(Cryst%natom),l_size_atm(Cryst%natom)) 
   do iat=1,Cryst%natom 
    l_size_atm(iat)=Pawtab(Cryst%typat(iat))%l_size 
   end do
   call nullify_pawfgrtab(Pawfgrtab) 
   call init_pawfgrtab(Pawfgrtab,l_size_atm) 
   deallocate(l_size_atm)
   compch_fft=greatest_real
   usexcnhat=MAXVAL(Pawtab(:)%vlocopt)
!  * 0 if Vloc in atomic data is Vbare    (Blochl s formulation)
!  * 1 if Vloc in atomic data is VH(tnzc) (Kresse s formulation)
   write(*,*)' sigma : using usexcnhat = ',usexcnhat 
!  
!  === Identify parts of the rectangular grid where the density has to be calculated ===
   call status(0,Dtfil%filstat,iexit,level,'call nhatgrid ')

   optgr0=Dtset%pawstgylm ; optgr1=0 ; optgr2=0 ; optrad=1-Dtset%pawstgylm
   if (Dtset%xclevel==2.and.usexcnhat>0) optgr1=Dtset%pawstgylm 

   call nhatgrid(Cryst%atindx1,gmet,MPI_enreg_seq,Cryst%natom,Cryst%nattyp,nfftf,ngfftf,Cryst%ntypat,&
&   optgr0,optgr1,optgr2,optrad,Pawfgrtab,Pawtab,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)
  end if ! End of PAW Initialization
! 
! Create Sigma header. 
! TODO this is for the next version, problems with symmorphy and k-points
! £  Bstruct%npwarr=Sp%npwvec  ! only now npwarr(:) has been read from KSS as well as occ
! £  call hdr_init(Bstruct,codvsn,Dtset,Hdr_sigma,Pawtab,0,Psps)
! £  allocate(occfact(bantot)) ; jj=0
! £  do isppol=1,Dtset%nsppol
! £   do ik=1,Dtset%nkpt
! £    do ib=1,Hdr_kss%nband(ik+Dtset%nkpt*(isppol-1))
! £     ii=(ik-1)*Hdr_kss%nband(ik)+(isppol-1)*Dtset%nkpt*Hdr_kss%nband(ik)+ib  !assuming nband(k) costant
! £     if (ib<=Sp%nbnds) then 
! £      jj=jj+1
! £      occfact(jj)=Hdr_kss%occ(ii)
! £     end if
! £    end do
! £   end do
! £  end do
! £  call hdr_update(bantot,1.0d20,1.0d20,Hdr_sigma,Cryst%natom,1.0d20,Dtset%rprimd_orig,&
! £  occfact,Pawrhoij,Dtset%usepaw,Dtset%xred_orig)
! £  call bstruct_clean(Bstruct) ; deallocate(occfact)
! £  call hdr_check(fform0,fform_kss,Hdr_sigma,Hdr_kss,'COLL',restart,restartpaw)
 call pclock(10)
 call pclock(20)
!
!=== Get KS occupation numbers and nbv(1:nsppol), valence band index ===
!* fixmom is passed to fermi.F90 to fix the problem with newocc in case of magnetic metals
 allocate(nbv(Sp%nsppol))
 call fermi(Hdr_kss,Sp%nbnds,Kmesh%nibz,Dtset%fixmom,Sp%nsppol,Kmesh%wt,ks_energy,ks_occ,nel,nbv,efermi)
 call pclock(30)
!
!=== In case of symmetrization, find little group ===
 allocate(Ltg_k(Sp%nkcalc))
 do ikcalc=1,Sp%nkcalc 
  call nullify_little_group(Ltg_k(ikcalc)) 
 end do

 if (Sp%symsigma/=0) then
  do ikcalc=1,Sp%nkcalc
   use_umklp=1 !££ 1 for 5.6 switch on umklapp but I should add an automatic test
   call setup_little_group(Sp%xkcalc(:,ikcalc),Qmesh,gmet,Cryst,Sp%npwvec,Gsph_Max%gvec,&
&   0,use_umklp,Dtset%prtvol,Ltg_k(ikcalc))
  end do
 end if

 allocate(Cprj_bz(Cryst%natom,Dtset%nspinor*Sp%nbnds*Kmesh%nbz*Sp%nsppol*Dtset%usepaw))

 if (Dtset%usepaw==1) then 
! === Symmetrize projected wave functions <Proj_i|Cnk> in the Full Brillouin zone ===
  call cprj_alloc(Cprj_bz,0,dimlmn)

! TODO add rotation in spinor space
  call paw_symcprj(Pawtab,Cryst,Dtset%nspinor,Sp%nbnds,nband_,&
&  Dtset%nsppol,Psps,Kmesh,Cprj_ibz,Pawang,dimlmn,Cprj_bz)
! 
! === Set up REAL Ylm(q+G) up to 2*l_max for each q-point ===
! * Note that Sp%npwx is always >= Sp%npwc
  mqmem__=Qmesh%nbz ; two_lmaxp1=2*Psps%mpsang-1 ; mpw_=Sp%npwx ; optder=0 
  allocate(ylm_q(mpw_*mqmem__,two_lmaxp1**2*Dtset%usepaw))
  allocate(ylmgr_q(mpw_*mqmem__,3+6*(optder/2),two_lmaxp1**2*Dtset%usepaw))
  allocate(npwarr(Qmesh%nbz),qg(3,mpw_*mqmem__)) ; npwarr(:)=mpw_
  allocate(qtmp(3,mqmem__))
  do iq=1,mqmem__
!  Copy the mesh but setting the small q equal to 0
   indx=(iq-1)*mpw_ ; qg(:,indx+1:indx+mpw_)=Gsph_Max%gvec(:,1:mpw_)
   qtmp(:,iq)=Qmesh%bz(:,iq) ; if (normv(qtmp(:,iq),gmet,'G')<GW_TOLQ0) qtmp(:,iq)=zero 
  end do

  call status(0,Dtfil%filstat,iexit,level,'call initylmg ')
! === Get real spherical harmonics in G space ===
! * Dtset%nband and Dtset%nsppol are not used in sequential mode.
  call initylmg(gprimd,qg,qtmp,mqmem__,MPI_enreg_seq,two_lmaxp1,mpw_,Dtset%nband,Qmesh%nbz,&
  npwarr,Dtset%nsppol,optder,Cryst%rprimd,Dtfil%unkg,Dtfil%unylm,ylm_q,ylmgr_q)
  deallocate(npwarr,qg)
! 
! === Evaluate oscillator matrix elements, pawrhox, btw partial waves ===
  allocate(pawrhox(2,Sp%npwx,Psps%lmnmax*(Psps%lmnmax+1)/2,Cryst%natom,Qmesh%nbz*Dtset%usepaw))  
  do iq=1,Qmesh%nbz
   indx=(iq-1)*Sp%npwx 
   call paw_mkrhox(Cryst,pawrhox_spl,gmet,Gsph_Max%gvec,rhoxsp_method,dim1_rhox,dim2_rhox,&
&   Psps,Pawang,Pawtab,qtmp(:,iq),Sp%npwx,ylm_q(indx+1:indx+Sp%npwx,:),pawrhox(:,:,:,:,iq)) 
  end do
  deallocate(qtmp,ylm_q,ylmgr_q)
 end if
 !
 ! === Set up tables for FFT ===
 ! * Sp%mG0 gives, for each reduced direction, the max G0 component to account for umklapp
 ! * igfft0 contains the FFT index of G used in fourdp.
 nG01d=2*Sp%mG0(1)+1 ; nG02d=2*Sp%mG0(2)+1 ; nG03d=2*Sp%mG0(3)+1
 allocate(igfft(Sp%npwvec,nG01d,nG02d,nG03d))

 call cigfft(Sp%mG0,Sp%npwvec,ngfft_gw,Gsph_Max%gvec,igfft)
 igfft0 => igfft(:,Sp%mG0(1)+1,Sp%mG0(2)+1,Sp%mG0(3)+1)
 call timab(402,2,tsec) ! sigma(1)
!
!FIXME This part might be CPU consuming in isolated systems
!Should take into account possible FFT //. Try to encapsulate EVERYTHING into Wf_info
 call pclock(40)
 call timab(404,1,tsec) ! sigma(2)
 allocate(irottb(nfftgw,Cryst%nsym)) ; irottb(:,:)=0 

 call setup_FFT_rotation(Cryst%nsym,Cryst%symrec,Cryst%tnons,nfftgw,ngfft_gw,irottb)
!£££ call FFT_rotations(Cryst,ngfft_gw,irottb)
 call pclock(50)
!
!$u(R^{-1}(r-\tau))$ where S=\transpose R^{-1} and k_BZ = S k_IBZ
!irottb is the FFT index of $R^{-1} (r-\tau)$.
 allocate(ktabr(nfftgw_tot,Kmesh%nbz))
 do ik_bz=1,Kmesh%nbz
  isym=Kmesh%tabo(ik_bz)
  do ifft=1,nfftgw_tot
   ktabr(ifft,ik_bz)=irottb(ifft,isym)
  end do
 end do

 write(msg,'(a)')' sigma: using augmented FFT mesh for density '
 call wrtout(std_out,msg,'COLL')
 call print_ngfft(ngfftf)
!£££
!
!=== Compute structure factor phases and large sphere cut-off ===
!WARNING cannot use Dtset%mgfft, this has to be checked better
!mgfft=MAXVAL(ngfftc(:))
!allocate(ph1d(2,3*(2*mgfft+1)*Cryst%natom),ph1df(2,3*(2*mgfftf+1)*Cryst%natom))
 write(std_out,*)' CHECK ',Dtset%mgfft,mgfftf
 if (Dtset%mgfft/=mgfftf) write(std_out,*)"WARNING Dtset%mgfftf /= mgfftf"
 allocate(ph1d(2,3*(2*Dtset%mgfft+1)*Cryst%natom),ph1df(2,3*(2*mgfftf+1)*Cryst%natom))
 call status(0,Dtfil%filstat,iexit,level,'call getph    ')
 call getph(Cryst%atindx,Cryst%natom,ngfftc(1),ngfftc(2),ngfftc(3),ph1d,Cryst%xred)
 if (Psps%usepaw==1.and.Pawfgr%usefinegrid==1) then
  call getph(Cryst%atindx,Cryst%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,Cryst%xred)
 else
  ph1df(:,:)=ph1d(:,:)
 end if
!
!=== Tables for the dense FFT mesh used for rhor ===
 allocate(igfftf(Sp%npwvec),mask(Sp%npwvec)) ; mask(:)=.TRUE.
 call kgindex(igfftf,Gsph_Max%gvec,mask,MPI_enreg_seq,ngfftf,Sp%npwvec)
 deallocate(mask)

 allocate(irottbf(nfftf,Cryst%nsym))
 call setup_FFT_rotation(Cryst%nsym,Cryst%symrec,Cryst%tnons,nfftf,ngfftf,irottbf)
!£££ call FFT_rotations(Cryst,ngfftf,irottbf)
 call pclock(60)

 call init_wf_info_2(Wf_info,igfft0,ngfft_gw)
!
!=== Initialize the object Wf_info_braket (wavefunctions for GW corrections) ===
 call nullify_wf_info(Wf_info_braket)
 call init_wf_info_1(Wf_info_braket,Dtset%gwmem,Dtset%paral_kgb,Sp%npwwfn,b1gw,b2gw,Sp%nkcalc,&
& Sp%nsppol,Dtset%nspden,Dtset%nspinor)

 call init_wf_info_2(Wf_info_braket,igfft0,ngfft_gw)
 call duplicate_wf_info(MPI_enreg,Wf_info,Wf_info_braket,Sp%kcalc,Kmesh)

#if defined DEBUG_MODE
!allocate(eigen(mband_*Kmesh%nibz*Sp%nsppol)) 
!call gw2abi(Kmesh%nibz,mband_,Sp%nsppol,ks_energy,eigen)
!call joint_dos(Qmesh%nibz,Qmesh%ibz,mband_,Kmesh%nibz,Kmesh%ibz,Dtset%nsppol,Dtset%nband,&
!& Cryst%timrev,Cryst%nsym,Cryst%symrec,Cryst%tnons,eigen,Hdr_kss%occ,1,0.1/Ha_eV,0.3/Ha_eV)
!deallocate(eigen)
!if (Dtset%usepaw==1) then
!call check_zarot(Cryst%nsym,Cryst%symrec,Cryst%timrev,Sp%npwvec,Cryst%rprimd,gprimd,Gsph_Max%gvec,Psps,
!& Pawang,Gsph_Max%rottb,Gsph_Max%rottbm1)
!end if
!do ik=1,Kmesh%nibz
!call nullify_Bands_Symmetries(BSym)
!call get_Bands_Sym_GW(Cryst%ntypat,Cryst%natom,Cryst%typat,Dtset%nspinor,Sp%nsppol,Sp%nbnds,Kmesh%nibz,Kmesh,Wf_info,&
!& Dtset%usepaw,Pawtab,Pawang,Psps,dimlmn,Cprj_ibz,MPI_enreg,ik,ks_energy(ik,:,:),Cryst%nsym,Cryst%symrec,Cryst%tnons,&
!& Cryst%indsym,.FALSE.,irottb,BSym,ierr)
!call print_Bands_Symmetries(Bsym,unitno=77)
!call destroy_Bands_Symmetries(Bsym)
!end do
#endif
!
!========================
!=== COMPUTE DENSITY ====
!========================
!
!=== Evaluate plane wave part of rhor ===
!* If gwpara=2 do the calculation in parallel inside calc_density
 allocate(rhor(nfftf,Dtset%nspden))
 call calc_density(Wf_info,Cryst,irottbf,Sp%nbnds,ngfftf,nfftf,igfftf,&
& ks_occ,rhor,Kmesh,MPI_enreg,tim_fourdp,.TRUE.)

!TODO this has to be done in a better way, moreover wont work for PAW
 call cutoff_density(ngfftf,Dtset%nspden,Dtset%nsppol,Vcp,rhor,MPI_enreg)
!
!=== Only for PAW: treat onsite contributions ===
 if (Dtset%usepaw==1) then 
  
! * Compute compensation density to be added to the density coming from smooth WFs
  nhatgrdim=0 ; if (Dtset%xclevel==2) nhatgrdim=usexcnhat*Dtset%pawnhatxc ; ider=2*nhatgrdim ; izero=0
  if (nhatgrdim>0) allocate(nhatgr(nfftf,Dtset%nspden,3))

  call pawmknhat(compch_fft,ider,izero,MPI_enreg_seq,Cryst%natom,nfftf,ngfftf,nhatgrdim,Dtset%nspden,&
&  Cryst%ntypat,Dtset%paral_kgb,Pawang,Pawfgrtab,nhatgr,nhat,Pawrhoij,Pawtab,Cryst%typat,Cryst%ucvol)

! === Initialize Variables/arrays related to the PAW spheres ===
! * Evaluate onsite energies, potentials, densities.
! * Initialize also "lmselect" (index of non-zero LM-moments of densities).
! TODO call init_paw_ij in scfcv and respfn, fix small issues
  allocate(KS_paw_ij(Cryst%natom)) 
  call nullify_paw_ij(KS_paw_ij)
  cplex=1 ; cplex_dij=Dtset%nspinor ; has_dijso=Dtset%pawspnorb ; has_dijU=Dtset%usepawu

  call init_paw_ij(KS_paw_ij,cplex,cplex_dij,Dtset%nspinor,Dtset%nsppol,&
&  Dtset%nspden,Dtset%pawspnorb,Cryst%natom,Cryst%ntypat,Cryst%typat,Pawtab,&
&  has_dijhartree=1,has_dijhat=0,has_dijxc=1,has_dijxc_val=1,has_dijso=has_dijso,has_dijU=has_dijU)

  allocate(KS_paw_an(Cryst%natom))
  call nullify_paw_an(KS_paw_an)
  call init_paw_an(Cryst%natom,Cryst%ntypat,Dtset%nspden,cplex,Dtset%pawxcdev,&
&  Dtset%pawspnorb,Cryst%typat,Pawang,Pawtab,KS_paw_an,has_vxcval=1)

  call status(0,Dtfil%filstat,iexit,level,'call pawdenpot')
! 
! === Calculate v_xc on the radial mesh with and without core charge ===
  nzlmopt=-1 ; option=0 ; compch_sph=greatest_real 
  call pawdenpot(compch_sph,KS_energies%e_paw,KS_energies%e_pawdc,Dtset%ixc,Cryst%natom,Dtset%nspden,&
&  Cryst%ntypat,nzlmopt,option,KS_paw_an,KS_paw_ij,Pawang,Dtset%pawprtvol,Pawrad,Pawrhoij,Dtset%pawspnorb,&
&  Pawtab,Dtset%pawxcdev,Cryst%typat,Dtset%xclevel,Psps%znuclpsp)

 end if !PAW

 call test_charge(nfftf,Dtset%nspden,rhor,Cryst%ucvol,nhat,Dtset%usepaw,&
& usexcnhat,Pawfgr%usefinegrid,compch_sph,compch_fft,drude_plsmf)
!
!=== Add compensation charge on FFT mesh for PAW then get rho(G) ===
 if (Dtset%usepaw==1) rhor(:,:)=rhor(:,:)+nhat(:,:) 
 allocate(rhog(2,nfftf)) 
 call fourdp(1,rhog,rhor(:,1),-1,MPI_enreg,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp)
 call prtrhomxmn(std_out,MPI_enreg,nfftf,ngfftf,Dtset%nspden,1,rhor)
!
!The following steps have been gathered in the setvtr routine:
!- get Ewald energy and Ewald forces
!- compute local ionic pseudopotential vpsp
!- eventually compute 3D core electron density xccc3d
!- eventually compute vxc and vhartr
!- set up vtrial
!**** NOTE THAT Vxc CONTAINS THE CORE DENSITY CONTRIBUTION ****
!
 call status(0,Dtfil%filstat,iexit,level,'call setvtr   ')

 allocate(grewtn(3,Cryst%natom),xred_dummy(3,Cryst%natom)) ; xred_dummy(:,:)=xred(:,:)
 nkxc=0 
 if (Dtset%nspden==1) nkxc=2 
 if (Dtset%nspden>=2) nkxc=3 ! check GGA and spinor that is messy !!!
 allocate(kxc(nfftf,nkxc))

 n3xccc=0 ; if (Psps%n1xccc/=0) n3xccc=nfftf
 allocate(xccc3d(n3xccc),vhartr_vtr(nfftf),vtrial(nfftf,Dtset%nspden),vpsp(nfftf),vxc_vtr(nfftf,Dtset%nspden))

 if (prior_to_56.and.(Dtset%usepaw/=1).and..TRUE.) then 
! FIXME old implementation with small cutoff.  now rho(G) is on the doubled grid 
  gsqcutf_eff_tmp=-one
  do ig=1,Sp%npwvec
   g1=real(Gsph_Max%gvec(1,ig)) ; g2=real(Gsph_Max%gvec(2,ig)) ; g3=real(Gsph_Max%gvec(3,ig))
   gsq=      gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
&       two*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
   gsqcutf_eff_tmp=max(gsqcutf_eff_tmp,gsq)
!  call getcut(boxcut,Sp%ecutwfn,gmet,gsqcutf_eff_tmp,Dtset%iboxcut,std_out,k0,ngfftf)
  end do
 else 
! TODO for v5.6  this is the Correct cutoff in G space, but the matrix elements of v_H 
! in tests v87 and t88 (SCGW) will change a bit
! (we recall getcut only if prior_to_56 ngfft==ngfft_gw) 
! call getcut(boxcut,Sp%ecutwfn,gmet,gsqcutf_eff_tmp,Dtset%iboxcut,std_out,k0,ngfftf)
  gsqcutf_eff_tmp=gsqcutf_eff
 end if

 optene=4 ; moved_atm_inside=0 ; moved_rhor=0 ; initialized=1 ; istep=1
 call setvtr(Cryst%atindx1,Dtset,KS_energies,gmet,gprimd,grewtn,gsqcutf_eff_tmp,initialized,istep,kxc,mgfftf,&
& moved_atm_inside,moved_rhor,MPI_enreg_seq,Cryst%nattyp,nfftf,ngfftf,nhat,nhatgr,nhatgrdim,nkxc,Cryst%ntypat,&
& Psps%n1xccc,n3xccc,optene,Pawtab,ph1df,Psps,rhog,rhor,Cryst%rmet,Cryst%rprimd,strsxc,Cryst%ucvol,usexcnhat,&
& vhartr_vtr,vpsp,vtrial,vxc_vtr,vxcavg,xccc3d,xred_dummy,Cryst%xred) 
!FIXME here xred is INOUT due to ionion_realSpace and xredcart, why?

!£££ call print_energies(KS_energies,1,std_out)

 if (Dtset%usepaw==1) then 
! === Dij computation (unsymmetrized quantities!) ===
  call status(0,Dtfil%filstat,iexit,level,'call pawdij   ')

  call pawdij(Dtset,Dtset%enunit,MPI_enreg_seq,Cryst%natom,nfftf,ngfftf,Dtset%nspden,Cryst%ntypat,&
&  KS_paw_an,KS_paw_ij,Pawang,Pawfgrtab,Dtset%pawprtvol,Pawrad,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,&
&  Cryst%typat,Cryst%ucvol,vtrial,vxc_vtr)

  call status(0,Dtfil%filstat,iexit,level,'call symdij   ')

  call symdij(Psps%indlmn,Cryst%indsym,Psps%lmnmax,Cryst%natom,Cryst%nsym,Cryst%ntypat,KS_paw_ij,Pawang,&
&  Dtset%pawprtvol,Dtset%symafm,Cryst%symrec,Cryst%typat)

! === Output of pseudopotential strength Dij and augmentation occupancies Rhoij ===
  call pawprt(Psps%indlmn,Dtset%enunit,Psps%lmnmax,Cryst%natom,Cryst%ntypat,KS_paw_ij,&
&  Dtset%pawprtvol,Pawrhoij,Pawtab,Cryst%typat)

 end if !PAW
!
!=== Calculate Vxc(b1,b2,k,s)=<b1,k,s|v_{xc}|b2,k,s>  for all the states included in GW ===
!* vxc_val is calculated without non linear core correction while vxc contains NLCC
 allocate(vxc    (b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab)) 
 allocate(vxc_val(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab)) 
 allocate(vhartr (b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab)) 
 allocate(vUpaw  (b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab)) 

!FIXME here Im not consistent as vxc_vtr in PAW contains the core contribution
!TODO the interface has to be cleaned but after the merge with Marc 
 call calc_vHxc_braket(Dtset,Kmesh%nibz,Sp%nkcalc,Kmesh%tab(Sp%kcalc(:)),b1gw,b2gw,Sp%nbnds,&
& Sp%nsig_ab,Sp%npwvec,gsqcutf_eff,nfftf_tot,nfftf,ngfftf,igfftf,Wf_info_braket,vpsp,vhartr_vtr,&
& vxc_vtr,Cprj_ibz,Pawtab,KS_paw_an,Pawang,Pawfgrtab,Pawrad,KS_paw_ij,MPI_enreg,Cryst,tim_fourdp,&
& rhor,rhog,usexcnhat,nhat,nhatgr,nhatgrdim,vhartr,vxc_val,vxc,vUpaw)

#if defined DEBUG_MODE
!if (Dtset%usepaw==1) then 
!call test_PAWH(Cryst%natom,Cryst%ntypat,Cryst%typat,Dtset%nspinor,Kmesh,gmet,Gsph_Max%gvec,nfftf,vpsp,&
!&  Dtset%usepaw,Pawtab,Cprj_ibz,b1gw,b2gw,Sp,ks_energy,vxc,vhartr,Wf_info,ngfftf,igfftf,MPI_enreg)
!end if
#endif

!FB: Please Matteo, do not break this coding!
!When gwcalctyp>10, the order of the bands can be interexchanged after
!the diagonalization. Therefore, we have to correctly assign the matrix
!elements to the corresponding bands and we cannot skip the following even though
!it looks unuseful.
 if (Dtset%gwcalctyp>=10) then 
! === Self-consistency both on eigenvalues and eigenfunctions ===
! Calculate matrix elements of the Hartree potential for GW states
  write(msg,'(2a)')ch10,' *************** KS Energies *******************'
  call wrtout(std_out,msg,'COLL')
 end if
 call pclock(70)

 if (Dtset%gwcalctyp<10) then 
! === One-shot GW: copy KS density, energies, ks_occ numbers and Ef and pass them to csigme ===
! Keep in mind however that these are KS quantities!
  allocate(qp_energy (Kmesh%nibz,Sp%nbnds,Sp%nsppol))
  allocate(qp_occ    (Kmesh%nibz,Sp%nbnds,Sp%nsppol))
  qp_occ(:,:,:)=ks_occ(:,:,:) ! This arrays will be passed to cisgme and hartrham
  allocate(rhor_qp(nfftf,Dtset%nspden))
! MG FIXME this is only to pass the automatic tests !, should be removed
! after the partial merge
  call calc_density(Wf_info,Cryst,irottbf,Sp%nbnds,ngfftf,nfftf,igfftf,&
&  qp_occ,rhor_qp,Kmesh,MPI_enreg,tim_fourdp,.TRUE.)

  call test_charge(nfftf,Dtset%nspden,rhor,Cryst%ucvol,nhat,Dtset%usepaw,&
&  usexcnhat,Pawfgr%usefinegrid,compch_sph,compch_fft,drude_plsmf)
! 
! These quantities will be passed to hartrham and csigme.
  qp_energy(:,:,:)=ks_energy(:,:,:) 
  rhor_qp(:,:)=rhor(:,:) 
  efermi_qp=efermi
  deallocate(ks_occ)
 else 
! 
! === Self-consistent GW. Read QP wavefunctions of the previous step (FBruneval) ===
! 1>>> Unitary transformation 2>>> QP energies 3>>> QP density 
  allocate(m_lda_to_qp(Sp%nbnds,Sp%nbnds,Kmesh%nibz,Sp%nsppol)) ; m_lda_to_qp(:,:,:,:)=czero
  allocate(qp_energy(Kmesh%nibz,Sp%nbnds,Sp%nsppol))
  allocate(rhor_p(nfftf,Dtset%nspden)) ! Previous density for mixing problem if case of G //
  rhor_p(:,:)=rhor(:,:) 
  
  call rdqps(Sp%gwcalctyp,Dtfil,Kmesh,Sp%nbnds,Dtset%nspden,Sp%nsppol,&
&  nscf,nfftf,ks_energy,nbsc,qp_energy,m_lda_to_qp,rhor_p)
! 
! === Compute QP wfg as linear combination of KS states ===
! * Wf_info%wfg is modified after calc_wf_qp
! * WARNING the first dimension of MPI_enreg MUST be Kmesh%nibz 
! TODO here we should use nbsc instead of nbnds

  call calc_wf_qp(MPI_enreg,Kmesh%nibz,Sp%nbnds,Wf_info%npwwfn,Sp%nsppol,Dtset%nspinor,&
&  m_lda_to_qp,my_minb,my_maxb,b1gw,b2gw,Wf_info%wfg)
! 
! === Reinit the storage mode of Wf_info and Wf_info_braket as wfs have just changed ===
! * Update also the wavefunctions for GW corrections on each processor
  call reinit_wf_info(Wf_info) 
  call reinit_wf_info(Wf_info_braket)
  call duplicate_wf_info(MPI_enreg,Wf_info,Wf_info_braket,Sp%kcalc,Kmesh)

  if (Dtset%usepaw==1) then
   call update_cprj(Cryst%natom,Kmesh%nibz,Sp%nbnds,Sp%nsppol,Dtset%nspinor,m_lda_to_qp,dimlmn,Cprj_ibz)
  end if
! 
! === Compute QP occupation numbers ===
  write(msg,'(3a)')ch10,' sigma : calculating QP occupation numbers ',ch10
  call wrtout(std_out,msg,'COLL')

  allocate(qp_occ(Kmesh%nibz,Sp%nbnds,Sp%nsppol)) 
  qp_occ(:,:,:)=ks_occ(:,:,:)
  call fermi(Hdr_kss,Sp%nbnds,Kmesh%nibz,Dtset%fixmom,Sp%nsppol,Kmesh%wt,qp_energy,qp_occ,nel,nbv,efermi_qp)
! 
! === Compute QP density using updated wfg ===
! * If gwpara==2 do the calculation in parallel
  allocate(rhor_qp(nfftf,Dtset%nspden))
  call calc_density(Wf_info,Cryst,irottbf,Sp%nbnds,ngfftf,nfftf,igfftf,&
&  qp_occ,rhor_qp,Kmesh,MPI_enreg,tim_fourdp,.TRUE.)

! === Additional allocations/calculations for SCGW+PAW ===
  if (Dtset%usepaw==1) then 

   call energies_init(QP_energies) ; QP_energies%e_corepsp=ecore/Cryst%ucvol
!  === Calculate new rhoij_qp from updated Cprj_ibz, note that use_rhoij_==1 ===
   allocate(QP_pawrhoij(Cryst%natom)) 
   call rhoij_alloc(1,nlmn,Dtset%nspden,Dtset%nsppol,QP_pawrhoij,Cryst%typat,use_rhoij_=1,use_rhoijres=1)
!  
!  Here cprj are unsorted, see  ctocprj.F90
   mband_=Sp%nbnds ; mkmem_=Kmesh%nibz 
   allocate(abiocc_qp(mband_*Kmesh%nibz*Dtset%nsppol))
   call gw2abi(Kmesh%nibz,mband_,Dtset%nsppol,qp_occ,abiocc_qp)

   call pawmkrhoij(Cryst%atindx1,Cprj_ibz,dimlmn,istwfk_,mband_,mkmem_,MPI_enreg_seq,Cryst%natom,&
&   Cryst%nattyp,nband_,Kmesh%nibz,Dtset%nspden,Dtset%nspinor,Dtset%nsppol,Cryst%ntypat,abiocc_qp,&
&   Dtset%pawprtvol,QP_pawrhoij,Dtfil%unpaw,Kmesh%wt)
   deallocate(abiocc_qp)

   choice=1 ; optrhoij=1 
   call symrhoij(choice,Psps%indlmn,Cryst%indsym,Psps%lmnmax,Cryst%natom,Cryst%nsym,Cryst%ntypat,&
&   optrhoij,Pawang,Dtset%pawprtvol,QP_pawrhoij,Dtset%symafm,Cryst%symrec,Cryst%typat)

   do iat=1,Cryst%natom 
    QP_pawrhoij(iat)%use_rhoij_=0 ; deallocate(QP_pawrhoij(iat)%rhoij_) 
   end do

   allocate(nhat_qp(nfftf,Dtset%nspden)) ; nhat_qp(:,:)=zero
   nhatgrdim=0 ; if (Dtset%xclevel==2) nhatgrdim=usexcnhat ; ider=2*nhatgrdim ; izero=0
   if (nhatgrdim>0) allocate(nhatgr_qp(nfftf,Dtset%nspden,3))

   call pawmknhat(compch_fft,ider,izero,MPI_enreg_seq,Cryst%natom,nfftf,ngfftf,nhatgrdim,Dtset%nspden,&
&   Cryst%ntypat,Dtset%paral_kgb,Pawang,Pawfgrtab,nhatgr_qp,nhat_qp,QP_pawrhoij,Pawtab,Cryst%typat,Cryst%ucvol)

   ! === Variables/arrays related to the PAW spheres for the QP Hamiltonian ===
   ! TODO call init_paw_ij in scfcv and respfn, fix small issues
   cplex=1 ; cplex_dij=Dtset%nspinor
   allocate(QP_paw_ij(Cryst%natom)) 
   call nullify_paw_ij(QP_paw_ij)
   call init_paw_ij(QP_paw_ij,cplex,cplex_dij,Dtset%nspinor,Dtset%nsppol,&
&   Dtset%nspden,Dtset%pawspnorb,Cryst%natom,Cryst%ntypat,Cryst%typat,Pawtab,&
&   has_dijhartree=1,has_dijhat=0,has_dijxc=0,has_dijxc_val=0,has_dijso=0,has_dijU=0)

   allocate(QP_paw_an(Cryst%natom))
   call nullify_paw_an(QP_paw_an)
   call init_paw_an(Cryst%natom,Cryst%ntypat,Dtset%nspden,cplex,Dtset%pawxcdev,&
&   Dtset%pawspnorb,Cryst%typat,Pawang,Pawtab,QP_paw_an,has_vxcval=1)
!  
!  === Evaluate on-site" energies, potentials, densities using QP density ===
!  Initialize also "lmselect" (index of non-zero LM-moments of densities).
   call status(0,Dtfil%filstat,iexit,level,'call pawdenpot')

   nzlmopt=-1 ; option=0 ; compch_sph=greatest_real
   call pawdenpot(compch_sph,QP_energies%e_paw,QP_energies%e_pawdc,Dtset%ixc,Cryst%natom,Dtset%nspden,&
&   Cryst%ntypat,nzlmopt,option,QP_paw_an,QP_paw_ij,Pawang,Dtset%pawprtvol,Pawrad,QP_pawrhoij,Dtset%pawspnorb,&
&   Pawtab,Dtset%pawxcdev,Cryst%typat,Dtset%xclevel,Psps%znuclpsp)

!  === Re-symmetrize PAW Cprj_bz in the full BZ ===
!  TODO add rotation in spinor space
   call paw_symcprj(Pawtab,Cryst,Dtset%nspinor,Sp%nbnds,nband_,&
&   Dtset%nsppol,Psps,Kmesh,Cprj_ibz,Pawang,dimlmn,Cprj_bz)
  end if

  call test_charge(nfftf,Dtset%nspden,rhor_qp,Cryst%ucvol,nhat_qp,Dtset%usepaw,&
&  usexcnhat,Pawfgr%usefinegrid,compch_sph,compch_fft,drude_plsmf) 

  if (Dtset%usepaw==1) rhor_qp(:,:)=rhor_qp(:,:)+nhat_qp(:,:) 
  call prtrhomxmn(std_out,MPI_enreg,nfftf,ngfftf,Dtset%nspden,1,rhor_qp)
! 
! === Output QP density ===
  if (rank==master.and.Dtset%prtden/=0) then
   rdwr=2 ; fformr=52 ; rdwrpaw=0
   fname=TRIM(filapp)//'_QP_DEN'
   call ioarr(accessfil,rhor_qp,Dtset,dummy,fformr,fname,Hdr_kss,MPI_enreg,&
&   nfftf,Pawrhoij_dum,rdwr,rdwrpaw,ngfftf)
   if (accessfil==3) then
!   Complete the geometry and electronic information with missing values from hdr_io().
    call abi_etsf_geo_put(Dtset,fname,Psps,Cryst%rprimd,Cryst%xred)
    call abi_etsf_electrons_put(Dtset,fname)
   end if
  end if
! 
! ==== Simple mixing of the densities to damp oscillations in the Hartree potential (FB061217) ===
! TODO Implement similar trick for PAW+GW, nhat is missing
  write(msg,'(2a,f5.3,a)')ch10,' sigma: mixing QP densities using rhoqpmix = ',Dtset%rhoqpmix,ch10
  call wrtout(std_out,msg,'COLL')
  rhor_qp(:,:)=rhor_p(:,:)+Dtset%rhoqpmix*(rhor_qp(:,:)-rhor_p(:,:))

  allocate(rhog_qp(2,nfftf)) 
  call fourdp(1,rhog_qp,rhor_qp(:,1),-1,MPI_enreg,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp)

! TODO this should be reported in the main output
  write(msg,'(a,f9.4)')' *** From QP density: number of electrons *** ',rhog(1,1)*Cryst%ucvol
  call wrtout(std_out,msg,'COLL')
! Calculate the band energy
  band_ene=zero
  do is=1,Sp%nsppol
   do ik=1,Kmesh%nibz
    do ib=1,Sp%nbnds
     band_ene=band_ene+Kmesh%wt(ik)*qp_energy(ik,ib,is)*qp_occ(ik,ib,is)
    end do 
   end do 
  end do
  write(msg,'(2a,es17.8)')ch10,' QP-Band Energy    [Ha] = ',band_ene
  call wrtout(std_out,msg,'COLL')

  call status(0,Dtfil%filstat,iexit,level,'call setvtr   ')
  nkxc=0 
  if (Dtset%nspden==1) nkxc=2 
  if (Dtset%nspden>=2) nkxc=3 !check GGA and spinor that is messy !!!
  allocate(kxc_qp(nfftf,nkxc))
! 
! **** NOTE THAT Vxc CONTAINS CORE DENSITY CONTRIBUTION ****
  n3xccc=0 ; if (Psps%n1xccc/=0) n3xccc=nfftf
  allocate(vhartr_qp_vtr(nfftf),vtrial_qp(nfftf,Dtset%nspden),vxc_qp_vtr(nfftf,Dtset%nspden))

  if (prior_to_56.and.(Dtset%usepaw/=1).and..TRUE.) then 
!  FIXME this is the old implementation with small cutoff: now rho(G) is on the augmented grid 
   gsqcutf_eff_tmp=-one
   do ig=1,Sp%npwvec
    g1=real(Gsph_Max%gvec(1,ig)) ; g2=real(Gsph_Max%gvec(2,ig)) ; g3=real(Gsph_Max%gvec(3,ig))
    gsq=      gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
&    two*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
    gsqcutf_eff_tmp=max(gsqcutf_eff_tmp,gsq)
!   call getcut(boxcut,Sp%ecutwfn,gmet,gsqcutf_eff_tmp,Dtset%iboxcut,std_out,k0,ngfftf)
   end do
  else 
!  TODO for v5.6  this is the Correct cutoff in G space, but the matrix elements of v_H 
!  in tests v87 and t88 (SCGW) will change a bit
!  (we recall getcut only if prior_to_56 and ngfft==ngfft_gw) 
!  call getcut(boxcut,Sp%ecutwfn,gmet,gsqcutf_eff_tmp,Dtset%iboxcut,std_out,k0,ngfftf)
   gsqcutf_eff_tmp=gsqcutf_eff
  end if

  optene=4 ; moved_atm_inside=0 ; moved_rhor=0 ; initialized=1 ; istep=1
  call setvtr(Cryst%atindx1,Dtset,QP_energies,gmet,gprimd,grewtn,gsqcutf_eff_tmp,initialized,istep,kxc_qp,mgfftf,&
&  moved_atm_inside,moved_rhor,MPI_enreg_seq,Cryst%nattyp,nfftf,ngfftf,nhat_qp,nhatgr,nhatgrdim,nkxc,Cryst%ntypat,&
&  Psps%n1xccc,n3xccc,optene,Pawtab,ph1df,Psps,rhog_qp,rhor_qp,Cryst%rmet,Cryst%rprimd,strsxc,Cryst%ucvol,usexcnhat,&
&  vhartr_qp_vtr,vpsp,vtrial_qp,vxc_qp_vtr,vxcavg_qp,xccc3d,xred_dummy,Cryst%xred)
! FIXME here xred is INOUT due to ionion_realSpace and xredcart, why?

  ehartree=half*SUM(rhor_qp(:,1)*vhartr_qp_vtr(:))/DBLE(nfftf)*Cryst%ucvol
  write(msg,'(2a,es17.8)')ch10,' QP-Hartree Energy [Ha] = ',ehartree
  call wrtout(std_out,msg,'COLL')

  deallocate(kxc_qp,rhor_p)
! Since plasmonpole model 2-3-4 depend on the Fourier components of the density
! in case of self-consistency we might calculate here the ppm coefficients using rhor_qp
 end if
 deallocate(irottb,irottbf)
!
!=== Setup the Hhartree Hamiltonian := T + v_{loc} + v_{nl} + v_{H}^s (FBruneval) ===
!1) Calculate the KS hamiltonian hlda(b1,b1,k,s)= <b1,k,s| H_s | b1,k,s>
 allocate(hlda(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab)) 
 hlda(:,:,:,:)=czero

 if (Dtset%nspinor==1) then
  do is=1,Sp%nsppol
   do ik=1,Kmesh%nibz
    do ib=b1gw,b2gw
     hlda(ib,ib,ik,is)=ks_energy(ik,ib,is)
    end do 
   end do 
  end do
 else 
! === Spinorial case === 
! * Note that here vxc contains the contribution of the core.
! * Scale ovlp if orthonormalization is not satisfied as npwwfn might be < npwvec.
! Here we fill only the entries relative to IBZ, because Wf_info_braket is needed 
! TODO add spin-orbit case
! FIXME this part has to be Rationalized
  do is=1,Sp%nsppol
   do ikcalc=1,Sp%nkcalc 
    ikibz=Kmesh%tab(Sp%kcalc(ikcalc)) ! Irred k-point for GW
    do ib=b1gw,b2gw
     cgup  => Wf_info_braket%wfg(1:Wf_info_braket%npwwfn,ib,ikcalc,is)
     cgdwn => Wf_info_braket%wfg(Wf_info_braket%npwwfn+1:2*Wf_info_braket%npwwfn,ib,ikcalc,is)
     shift=Sp%nspinor*Sp%nbnds*(ikibz-1)
     idx_up =shift+(2*ib-1) ; idx_dwn=idx_up+1
     ovlp(1) = overlap_cmplx(cgup ,cgup ,Dtset%usepaw,Cprj_ibz(:,idx_up ),Cprj_ibz(:,idx_up ),Cryst%typat,Pawtab)
     ovlp(2) = overlap_cmplx(cgdwn,cgdwn,Dtset%usepaw,Cprj_ibz(:,idx_dwn),Cprj_ibz(:,idx_dwn),Cryst%typat,Pawtab)
!    write(77,*)ovlp(1),ovlp(2)
     norm=REAL(ovlp(1)+ovlp(2))
     ovlp(1)=REAL(ovlp(1)/norm)
     ovlp(2)=REAL(ovlp(2)/norm)
!    ovlp(2) = cone-ovlp(1)
     hlda(ib,ib,ikibz,1)=ks_energy(ikibz,ib,1)*ovlp(1)-vxc(ib,ib,ikibz,3)
     hlda(ib,ib,ikibz,2)=ks_energy(ikibz,ib,1)*ovlp(2)-vxc(ib,ib,ikibz,4)
     hlda(ib,ib,ikibz,3)=vxc(ib,ib,ikibz,3)
     hlda(ib,ib,ikibz,4)=vxc(ib,ib,ikibz,4)
    end do
   end do
  end do
 end if
!
!=== Initialize Sigma results ===
!TODO it is better if we use ragged arrays indexed by the k-point
 call nullify_sigma_results(Sr)
 call init_sigma_results(Sp,Kmesh,Sr)
!
!2) Calculate hhartree(b1,b2,k,s)= <b1,k,s|T+v_{loc}+v_{nl}+v_{H}|b2,k,s>
!* The representation depends if we are updating the wfs or not.
!* vUpaw is zero unless we are using LDA+U as starting point, see calc_vHxc_braket
!* Note that vH matrix elements are calculated using the true uncutted interaction.

 if (Dtset%gwcalctyp<10) then
! === Self-consistent only on energies, use the KS representation ===
  Sr%hhartree(:,:,:,:)=hlda(:,:,:,:)-vxc_val(:,:,:,:)-vUpaw(:,:,:,:)
 else
! === Self-consistent both on energies and wavefunctions ===
! * Get bare Hamiltonian  $H_{bare}= T+v_{loc}+ v_{nl}$ in the KS representation
  allocate(hbare(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab))
  hbare(:,:,:,:)=hlda(:,:,:,:)-vhartr(:,:,:,:)-vxc_val(:,:,:,:)-vUpaw(:,:,:,:)

! * Change basis from KS to QP, hbare is overwritten: A_{QP} = U^\dagger A_{KS} U 
  allocate(htmp(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab))
  allocate(ctmp(b1gw:b2gw,b1gw:b2gw),uks2qp(b1gw:b2gw,b1gw:b2gw))
  htmp(:,:,:,:)=hbare(:,:,:,:) ; hbare(:,:,:,:)=czero

  do is=1,Sp%nsppol
   do ik=1,Kmesh%nibz
    uks2qp(:,:) = m_lda_to_qp(b1gw:b2gw,b1gw:b2gw,ik,is)

    do iab=1,Sp%nsig_ab
     is_idx=is ; if (Sp%nsig_ab>1) is_idx=iab
     ctmp(:,:)=MATMUL(htmp(:,:,ik,is_idx),uks2qp(:,:))
     hbare(:,:,ik,is_idx)=MATMUL(TRANSPOSE(CONJG(uks2qp(:,:))),ctmp(:,:))
    end do

   end do !ik
  end do !is 
  deallocate(htmp,ctmp,uks2qp)

! * Calculate the QP Hartree potential ===
  write(msg,'(2a)')ch10,' *************** QP Energies *******************'
  call wrtout(std_out,msg,'COLL')

  allocate(vxc_qp   (b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab)) 
  allocate(vxcval_qp(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab)) 
  allocate(vUpaw_qp (b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab)) 

! FIXME here Im not consistent as vxc_vtr in PAW contains the core contribution, vpsp can be removed
! TODO the interface has to be cleaned but after the merge with Marc 
  call calc_vHxc_braket(Dtset,Kmesh%nibz,Sp%nkcalc,Kmesh%tab(Sp%kcalc(:)),b1gw,b2gw,Sp%nbnds,&
&  Sp%nsig_ab,Sp%npwvec,gsqcutf_eff,nfftf_tot,nfftf,ngfftf,igfftf,Wf_info_braket,vpsp,vhartr_qp_vtr,&
&  vxc_qp_vtr,Cprj_ibz,Pawtab,KS_paw_an,Pawang,Pawfgrtab,Pawrad,QP_paw_ij,MPI_enreg,Cryst,tim_fourdp,&
&  rhor_qp,rhog_qp,usexcnhat,nhat_qp,nhatgr,nhatgrdim,vhartr,vxcval_qp,vxc_qp,vUpaw_qp)

  deallocate(rhog_qp,vxc_qp,vhartr_qp_vtr,vxc_qp_vtr,vUpaw_qp)

  Sr%hhartree(:,:,:,:)=hbare(:,:,:,:)+vhartr(:,:,:,:)
 end if ! gwcalctyp<10
!
!=== Free some memory ===
 if (allocated(vhartr)) deallocate(vhartr)
 if (allocated(hbare )) deallocate(hbare )
 if (allocated(hlda  )) deallocate(hlda  )
!
!=== Prepare the storage of QP amplitudes and energies ===
!* Initialize with KS wavefunctions and energies.
 Sr%eigvec_qp(:,:,:,:)=czero 
 Sr%en_qp_diago(:,:,:)=zero
 do ib=1,Sp%nbnds
  Sr%en_qp_diago(ib,:,:)=ks_energy(:,ib,:)
  Sr%eigvec_qp(ib,ib,:,:)=cone
 end do 
!
!=== Store <n,k,s|V_xc[n_val]|n,k,s> and <n,k,s|V_U|n,k,s> ===
!* Note that we store the matrix elements of V_xc in the KS basis set, not in the QP basis set
!* Matrix elements of V_U are zero unless we are using LDA+U as starting point 
 do ib=b1gw,b2gw
  Sr%vxcme(ib,:,:)=vxc_val(ib,ib,:,:)
  Sr%vUme (ib,:,:)=vUpaw  (ib,ib,:,:)
 end do
 if (Dtset%usepaw==0) then
  ltest=(ALL(ABS(Sr%vUme)<tol6))
  call assert(ltest,'Sr%vUme differs from zero',__FILE__,__LINE__)
 end if
 deallocate(vxc,vxc_val,vUpaw)
 call pclock(80)
!
!=== Initial guess for GW energies ===
!* Save also energies of the previous iteration.
 do is=1,Sp%nsppol
  do ik=1,Kmesh%nibz
   do ib=1,Sp%nbnds
    Sr%e0 (ib,ik,is)=qp_energy(ik,ib,is)
    Sr%egw(ib,ik,is)=qp_energy(ik,ib,is) 
   end do
!  TODO add a check at the beginning on nbv
   Sr%e0gap(ik,is)=Sr%e0(nbv(is)+1,ik,is)-Sr%e0(nbv(is),ik,is)
  end do
 end do
!
!=== If required apply a scissor operator or update the energies ===
 Sp%soenergy=Dtset%soenergy
 if (Sp%soenergy>0.1d-4) then
  write(msg,'(6a,f10.5,a)')ch10,&
&  ' sigma : performing a first self-consistency',ch10,&
&  '  update of the energies in G by a scissor operator',ch10, &
&  '  applying a scissor operator of [eV] ',Sp%soenergy*Ha_eV,ch10
  call wrtout(std_out,msg,'COLL') 
  do is=1,Sp%nsppol
   if (Sp%nbnds>=nbv(is)+1) then 
    Sr%egw(nbv(is)+1:Sp%nbnds,:,is)=Sr%egw(nbv(is)+1:Sp%nbnds,:,is)+Sp%soenergy
!   RS patch for scissor operator
    if (Sp%nbnds>=nbv(is)+1) qp_energy(:,nbv(is)+1:Sp%nbnds,is)=qp_energy(:,nbv(is)+1:Sp%nbnds,is)+Sp%soenergy
   end if
  end do
 else if (update_energies) then
! TODO This should be done in rdqps.
  allocate(sr_gwenergy(Sp%nbnds,Kmesh%nibz,Sp%nsppol))
  write(msg,'(4a)')ch10,&
&  ' sigma : performing a first self-consistency',ch10,&
&  '  update of the energies in G by a previous GW calculation'
  call wrtout(std_out,msg,'COLL') 
  call wrtout(ab_out,msg,'COLL')

  call rdgw(Kmesh%nibz,Sp%nbnds,nbv,Sp%nsppol,Kmesh%ibz,sr_gwenergy)

  do ib=1,Sp%nbnds
   do ik=1,Kmesh%nibz
    do is=1,Sp%nsppol
     Sr%egw(ib,ik,is)=sr_gwenergy(ik,ib,is)
    end do
   end do
  end do
  deallocate(sr_gwenergy)
! TODO here we should call fermi to recalculate the new fermi level
 end if
 call pclock(90)
!
!In case of AC refer all the energies wrt to the fermi level
!Take care because results from ppmodel cannot be used for AC 
!FIXME check ks_energy or qp_energy (in case of SCGW?)
 if (mod10==1) then 
! All these quantities will be passed to csigme
! if I skipped the self-consistent part then here I have to use fermi
  qp_energy  = qp_energy -efermi_qp
  Sr%egw = Sr%egw-efermi_qp
  Sr%e0  = Sr%e0 -efermi_qp
  oldefermi=efermi_qp 
  efermi=zero
 end if 

 ! === Retrieve e^-1 either from the _SCR or the _SUSC file ===
 !  * If Er%mqmem==0, we allocate and read a single q-slice inside csigme. 
 ! TODO Er%nomega should be initialized so that only the frequencies really needed are stored in memory
 ! TODO TDDFT not yet operative
 if (Er%mqmem/=0) then 
  allocate(kxcg(1,1))
  call get_epsm1(Er,Vcp,Dtfil,0,0,Dtset%accesswff,Dtset%localrdwf,kxcg,gmet,MPI_enreg)
  deallocate(kxcg)
 else 
  write(msg,'(3a)')ch10,&
&  ' sigma : gwmem is smaller than 10,',&
&  ' allocating single slice of screening (slower but less memory)' 
  call wrtout(std_out,msg,'COLL')
 end if 

 !call q0fit(Qmesh%nibz,Qmesh%ibz,Gsph_Max%gvec,Er%nomega,Er%omega,Er%npwe,Er%epsm1,qcut,metal,&
 !& Cryst%nsym,REAL(Cryst%symrec,dp),Cryst%timrev,gprimd)

 !Added by Rshaltaf for the vertex correction inclusion
 !We should discuss on this subroutine, it requires too much memory and it is obsolete
 !Should also check the dimensions declared in input, now I dont have time
 !For sure it would not work if we swith to prior_to_56=.TRUE.
 if (Dtset%gwgamma==1) then
  call eps1_tc(Dtset,MPI_enreg,ngfft_gw,nfftgw_tot,rhor,Cryst%rprimd,igfft(1:Sp%npwc,Sp%mG0(1)+1,Sp%mG0(2)+1,Sp%mG0(3)+1),&
&  Sp%npwc,gmet,gprimd,Gsph_Max%gvec,Er%nomega,Er%epsm1,Qmesh%nibz,Cryst%ucvol,Qmesh%ibz,Er%omega,Dtfil,Hdr_kss,&
&  Er%npwwfn_used,Sp%npwvec,Er%nbnds_used)
 end if
 ! If nfreqmidm == -1000, the default, the frequency moment will not be calculated
 if (Dtset%nfreqmidm/=-1000) then
  if (Dtset%nfreqmidm<0) then
   write(msg,'(a,i4,a,a)')&
&   'Calculating the ',Dtset%nfreqmidm,'th frequency moment of imaginary part of inverse DM',ch10
  else if (Dtset%nfreqmidm>0) then
   write(msg,'(a,i4,a,a)')&
&   'Calculating the ',Dtset%nfreqmidm,'th frequency moment of imaginary part of DM',ch10
  else
   write(msg,'(a,i4,a,a)')&
&   'Calculating the ',Dtset%nfreqmidm,&
&   'th frequency moment of imaginary part of full polarizibility',ch10
  end if
  call wrtout(std_out,msg,'COLL')
  call calc_ffm(Er%epsm1,Qmesh%nibz,Sp%npwc,Er%nomega,Er%omega,gprimd,Qmesh%ibz,Sp%ppmodel,Gsph_Max%gvec,Dtset%nfreqmidm)
 end if

 ! === Calculate plasmonpole model parameters ===
 ! TODO In case of PAW rhor contains only tn + nhat, PPmodels requiring n(G) not tested
 ! TODO Maybe its better if we use mqmem as input variable
 ! TODO should pass rhor_qp but drude_plsmf of KS density!
 if (Dtset%usepaw==1.and.(Sp%ppmodel/=1.and.Sp%ppmodel/=0)) STOP 'PAW + this ppmodel not tested'
 call nullify_PPmodel(PPm)
 call init_PPmodel(PPm,Qmesh,Sp%ppmodel,Sp%npwc,Er%mqmem,drude_plsmf)

 if (Er%mqmem/=0) then
  call setup_ppmodel(PPm,Dtset%paral_kgb,Qmesh,Sp,Er,MPI_enreg,nfftf,Gsph_Max%gvec,ngfftf,&
&  gmet,gprimd,rhor(:,1))
 end if
 !
 ! === Write header and open files to output the final results ===
 if (rank==master) then
  call write_sigma_results_header(Sp,Er,Cryst,Kmesh,Qmesh)
  fname=TRIM(Dtfil%filnam_ds(4))//'_GW'  ; open(unt_gw,file=fname,status='unknown',form='formatted')
  write(unt_gw,*)Sp%nkcalc,Sp%nsppol
  fname=TRIM(Dtfil%filnam_ds(4))//'_SIG' ; open(unt_sig,file=fname,status='unknown',form='formatted')
  fname=TRIM(Dtfil%filnam_ds(4))//'_SGR' ; open(unt_sgr,file=fname,status='unknown',form='formatted')
  if (mod10==1) then 
   ! Open File for Matsubara axis 
   fname=TRIM(Dtfil%filnam_ds(4))//'_SGM' ; open(unt_sgm,file=fname,status='unknown',form='formatted')
  end if 
 end if 
 call timab(404,2,tsec) ! sigma(2)
 call pclock(100)
!
!=================================================================
!=== Calculate self-energy and output results for each k-point ===
!=================================================================
!
!TODO here we have what I call the stupid input file BUG (or stupid man BUG):
!if one calculates the GW corrections in the same k-point twice the results
!for the correlation part of sigma are different, since in the second
!calculation the starting point is updated and is different from the LDA value
!this can be useful to update the energy but I have to add a check in chkinp
!
 call timab(405,1,tsec) ! sigma(csigme)

 do ikcalc=1,Sp%nkcalc
  write(msg,'(2a,i5)')ch10,&
&  ' Calculating GW corrections for k-point number : ',ikcalc
  call wrtout(std_out,msg,'COLL')

  ikibz=Kmesh%tab(Sp%kcalc(ikcalc)) ! Irred k-point for GW

  if (MPI_enreg%gwpara==1) then
!  === Parallelization over k-points, redefine the distribution of k-points ===
!  NOTE that Kmesh%nbz is used as first dimension of proc_distrb. This implies that 
!  proc_distrb *MUST* be redefined if it is passed to routines employing the IBZ indexing.
   deallocate(MPI_enreg%proc_distrb) 
   allocate(MPI_enreg%proc_distrb(Kmesh%nbz,Sp%nbnds,Sp%nsppol)) 
!  If nprocs>ntasks, proc_distrb==-999 for rank>ntasks-1 and no harm should be done
   MPI_enreg%proc_distrb=-999 
   allocate(istart(nprocs),istop(nprocs))
   if (Sp%symsigma==0) then 
!   * No symmetries, divide the full BZ among procs
    ntasks=Kmesh%nbz 
    call split_work2(ntasks,nprocs,istart,istop)
    do irank=0,nprocs-1
     ii=istart(irank+1) ; jj=istop(irank+1)
     MPI_enreg%proc_distrb(ii:jj,:,:)=irank
    end do 
   else if (Sp%symsigma/=0) then  
!   * Divide the IBZ_q among procs. Distrb might be not so efficient for particular qs
!   Here proc_distrb is -999 for all the k-points not in the IBZ_q 
    ntasks=SUM(Ltg_k(ikcalc)%ibzq(:)) 
    call split_work2(ntasks,nprocs,istart,istop)  
!   Identify q and G0 where q+G0=k_j-k_i
    do irank=0,nprocs-1
     do ik=1,Kmesh%nbz
      kgwmk(:)= Sp%xkcalc(:,ikcalc)-Kmesh%bz(:,ik) ! Warn xkcalc must be inside the BZ
      do iq=istart(irank+1),istop(irank+1)
       call findqg0(iqm,g0,kgwmk,Qmesh%nbz,Qmesh%bz,Sp%mG0) 
       if (Ltg_k(ikcalc)%bz2ibz(iqm)==iq) MPI_enreg%proc_distrb(ik,:,:)=irank
      end do
     end do 
    end do
   end if 
   deallocate(istart,istop)
!  === Announce the treatment of k-points by each proc ===
   do ik=1,Kmesh%nbz 
    do is=1,Sp%nsppol
     if (MPI_enreg%proc_distrb(ik,Sp%nbnds,is)==rank) then 
      write(msg,'(3(a,i4))')'P sigma : treating k-point ',ik,' and spin ',is,' by node ',rank
      call wrtout(std_out,msg,'PERS')
     end if
    end do 
   end do 
  end if ! gwpara==1
! 
! === Call csigme to calculate matrix elements of Sigma ===
  call status(ikcalc,Dtfil%filstat,iexit,level,'call csigme   ')

  call csigme(ikcalc,Dtset,Cryst,Dtfil,Sp,Sr,Er,Gsph_Max,Vcp,Kmesh,Qmesh,Ltg_k(ikcalc),PPm,&
&  Pawang,Pawtab,Psps,Cprj_bz,Wf_info,Wf_info_braket,MPI_enreg,pawrhox,Gsph_Max%gvec,ktabr,ngfft_gw,igfft,&
&  nfftgw_tot,qp_energy,qp_occ,efermi_qp,my_minb,my_maxb,ngfftf,nfftf,rhor)
! 
! === Calculate direct gap for each spin and print out final results ===
  do is=1,Sp%nsppol
   if (Sp%maxbnd(ikcalc)>=nbv(is)+1) then
    Sr%egwgap (ikibz,is)=  Sr%egw(nbv(is)+1,ikibz,is) -  Sr%egw(nbv(is),ikibz,is)
    Sr%degwgap(ikibz,is)= Sr%degw(nbv(is)+1,ikibz,is) - Sr%degw(nbv(is),ikibz,is)
   else
!   The "gap" cannot be computed
    Sr%e0gap  (ikibz,is)=zero
    Sr%egwgap (ikibz,is)=zero
    Sr%degwgap(ikibz,is)=zero
   end if
  end do
  if (rank==master) then 
   call write_sigma_results(Sp,Sr,ikcalc,ikibz,Kmesh,Dtset%usepawu,ks_energy)
  end if
  call pclock(100+ikcalc)
 end do !ikcalc
 call timab(405,2,tsec) ! sigma(csigme)
!
!=== If self-consistent recalculate new occupation numbers as well as the Fermi level ===
 if (Sp%gwcalctyp>=10) then 
  do ib=1,Sp%nbnds
   qp_energy(:,ib,:)=Sr%en_qp_diago(ib,:,:)
  end do
  call fermi(Hdr_kss,Sp%nbnds,Kmesh%nibz,Dtset%fixmom,Sp%nsppol,Kmesh%wt,qp_energy,qp_occ,nel,nbv,efermi_qp)
  write(msg,'(a,3x,2(es16.6,a))')' New Fermi energy : ',efermi_qp,' Ha ,',efermi_qp*Ha_eV,' eV'
  call wrtout(std_out,msg,'COLL') 
  call wrtout(ab_out,msg,'COLL')
 end if
!
!=== If all k-points and all occupied bands are calculated, output EXX ===
 if (rank==master.and.Sp%nkcalc==Kmesh%nibz.and.ALL(Sp%minbnd(:)==1).and.ALL(Sp%maxbnd(:)>=MAXVAL(nbv(:)))) then
  exchange_energy=zero
  do is=1,Sp%nsppol
   do ik=1,Kmesh%nibz
    do ib=b1gw,b2gw
     if (Sp%nsig_ab==1) then
      exchange_energy = exchange_energy + half*qp_occ(ik,ib,is)*Kmesh%wt(ik)*Sr%sigxme(ib,ik,is)
     else
      exchange_energy = exchange_energy + half*qp_occ(ik,ib,is)*Kmesh%wt(ik)*SUM(Sr%sigxme(ib,ik,:))
     end if
    end do
   end do
  end do
  write(msg,'(a,2(es16.6,a))')' New Exchange energy : ',exchange_energy,' Ha ,',exchange_energy*Ha_eV,' eV'
  call wrtout(std_out,msg,'COLL') 
  call wrtout(ab_out,msg,'COLL')
 end if
!
!=== Write SCF data in case of self-consistent calculation === 
!* Save Sr%en_qp_diago, Sr%eigvec_qp and m_lda_to_qp in the _QPS file.
!* Note that rhor_qp is the QP density read from the input _QPS file.
!TODO 1) Write information on the FFT grid 
!2) Add abinit header for PAW
 if (rank==master.and.Sp%gwcalctyp>=10) then 
  call wrqps(Dtfil,Sp,Kmesh,Dtset%nspden,nscf,nfftf,Sr,m_lda_to_qp,rhor_qp)
 end if 
!
!----------------------------- END OF THE CALCULATION ------------------------
!
!=== Close Files === 
 if (rank==master) then
  close(unt_gw )
  close(unt_sig)
  close(unt_sgr)
  if (mod10==1) close(unt_sgm) 
 end if

 call status(0,Dtfil%filstat,iexit,level,'deallocate    ')

 deallocate(igfft,igfftf,ktabr)
 deallocate(ks_energy,qp_occ,qp_energy,rhor_qp,rhor)
 deallocate(nbv,nband_,istwfk_,ph1d,ph1df)
 deallocate(vhartr_vtr,vtrial,vpsp,vxc_vtr,rhog)
 deallocate(kxc,xccc3d,grewtn,xred_dummy)
 deallocate(Pawfgr%fintocoa,Pawfgr%coatofin)

 if (allocated(ks_occ     )) deallocate(ks_occ) 
 if (allocated(m_lda_to_qp)) deallocate(m_lda_to_qp)
 if (allocated(vtrial_qp  )) deallocate(vtrial_qp) 

 if (associated(MPI_enreg%proc_distrb)) deallocate(MPI_enreg%proc_distrb)
!
!=== Destroy the dinamic arrays in the local data structures ===
!* Optional deallocation for PAW
 if (Dtset%usepaw==1) then 
  deallocate(dimlmn,nlmn,pawrhox_spl,pawrhox)
  deallocate(nhat) ; if (nhatgrdim>0) deallocate(nhatgr)
  call cprj_free(Cprj_ibz)          ; deallocate(Cprj_ibz)
  call cprj_free(Cprj_bz )          ; deallocate(Cprj_bz )
  call rhoij_free(Pawrhoij)         ; deallocate(Pawrhoij)
  call destroy_pawfgrtab(Pawfgrtab) ; deallocate(Pawfgrtab)
  call destroy_paw_ij(KS_paw_ij)    ; deallocate(KS_paw_ij)
  call destroy_paw_an(KS_paw_an)    ; deallocate(KS_paw_an)
  if (Dtset%gwcalctyp>=10) then 
   deallocate(nhat_qp) ; if (nhatgrdim>0) deallocate(nhatgr_qp)
   call rhoij_free(QP_pawrhoij)   ; deallocate(QP_pawrhoij)
   call destroy_paw_ij(QP_paw_ij) ; deallocate(QP_paw_ij)
   call destroy_paw_an(QP_paw_an) ; deallocate(QP_paw_an)
  end if
 end if

 call destroy_wf_info(Wf_info)
 call destroy_wf_info(Wf_info_braket)
 do ikcalc=1,Sp%nkcalc
  call destroy_little_group(Ltg_k(ikcalc))
 end do
 call destroy_BZ_mesh_type(Kmesh)
 call destroy_BZ_mesh_type(Qmesh)
 call destroy_Gvectors(Gsph_Max)
 call destroy_Coulombian(Vcp) 
 call destroy_Crystal_structure(Cryst)
 call destroy_Sigma_results(Sr)
 call destroy_Sigma_parameters(Sp)
 call destroy_Epsilonm1_results(Er)
 call destroy_PPmodel(PPm)
 !TODO avoid preprocessor in testscr skipping hdr_comm inside the procedure if !MPI
!£ call hdr_clean(Hdr_sigma)
 call hdr_clean(Hdr_kss)

 call timab(401,2,tsec)
 call pclock(9999)
 call status(0,Dtfil%filstat,iexit,level,'exit          ')

#if defined DEBUG_MODE
 write(msg,'(a)')' sigma : ended'
 call wrtout(std_out,msg,'PERS') 
 call flush_unit(std_out)
#endif

#endif
!Workaround for buggy PGI6

end subroutine sigma
!!***
