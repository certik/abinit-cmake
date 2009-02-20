!{\src2tex{textfont=tt}}
!!****f* ABINIT/screening
!! NAME
!! screening
!!
!! FUNCTION
!! Calculate screening and dielectric functions
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (GMR, VO, LR, RWG, MT, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! codvsn=code version
!! Dtfil <type(datafiles_type)>=variables related to files
!! iexit= exit flag
!! Pawang <type(pawang_type)>=paw angular mesh and related data
!! Pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!! Pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  Before entering the first time in screening, a significant part of Psps has been initialized:
!!  the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid, ntypat,n1xccc,usepaw,useylm, 
!!  and the arrays dimensioned to npsp. All the remaining components of Psps are to be initialized in 
!!  the call to pspini. The next time the code enters screening, Psps might be identical to the
!!  one of the previous Dtset, in which case, no reinitialisation is scheduled in pspini.F90.
!! rprim(3,3)=dimensionless real space primitive translations
!! xred(3,natom) = reduced atomic coordinates
!!
!! OUTPUT
!! Output is written on the main output file.
!! The symmetrical inverse dielectric matrix is stored in the _SCR file
!!
!! SIDE EFFECTS
!!  Dtset<type(dataset_type)>=all input variables for this dataset
!!
!! NOTES
!! USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut) for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ... are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg) for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...Total density, potentials, ... are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used. It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf) are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! PARENTS
!!      driver,drivergw
!!
!! CHILDREN
!!      calc_wf_qp,ccgradvnl,cchi0,cchi0q0,cigfft,ckxcldag,cvc,density,distrb2
!!      fermi,fftwfn,findnq,findq,findshells,hdr_clean,hermitianize,identk
!!      lattice,leave_new,matcginv,matrginv,metric,mkrdim,pclock,printcm,printv
!!      rdgw,rdkss,rdlda,rdldaabinit,rdqps,setmesh,setshells,surot,testlda
!!      timab,wrscr,wrtout
!!
!! SOURCE
#include "abi_errors.h" 

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine screening(acell,codvsn,Dtfil,Dtset,iexit,MPI_enreg,Pawang,Pawrad,Pawtab,Psps,rprim,xred)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : GW_TOLQ0, GW_TOLQ, czero_gw
 use m_errors
 use m_numeric_tools, only : print_arr
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_13nonlocal
 use interfaces_13paw
 use interfaces_13psp
 use interfaces_13xc
 use interfaces_14iowfdenpot
 use interfaces_15common
 use interfaces_15gw
 use interfaces_lib00numeric
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
 type(Pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Dtset%usepaw)
 type(Pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Dtset%usepaw)

!Local variables ------------------------------
 character(len=50),parameter :: FILE__='screening.F90'
 character(len=4) :: ctype='RPA ',tag
!scalars
 integer,parameter :: NOMEGAGAUSS=30,NOMEGAREAL=201,level=23,tim_fourdp=4
 integer,save :: nsym_old=-1
 integer :: choice,cplex,cplex_dij,dim1_rhox,dim2_rhox,dim_wings,enforce_sym,i3,iat,ib,ider,ierr
 integer :: ifft,ig,ii,ij_size,ik,ik_bz,ik_ibz,ikq,ilmn,initialized,io,iomega
 integer :: iq,iqp,irank,is,ispden,isppol,istat,isym,itypat,itype,izero,jj
 integer :: label,lm_size,lmn2_size,master,method,mg0sh,mgfftf,mgfftgw
 integer :: mkmem_,mod10,moved_atm_inside,moved_rhor,mpsang,my_maxb,my_minb
 integer :: my_nbnds,nG01d,nG02d,nG03d,nbcw,nbnds_kss,nbsc,nbvw,nel=0
 integer :: nfftc_tot,nfftf,nfftf_tot,nfftgw,nfftgw_tot,ng_kss,nhatgrdim,nprocs
 integer :: npwe_pG0,nqptdm,nqsm,nscf,nsheps_pG0,nsym_kss,ntasks,nzlmopt,optfil
 integer :: optgr0,optgr1,optgr2,option,approx_type,option_test,optrad,optrhoij,psp_gencond,rank
 integer :: rhoxsp_method,spaceComm,timrev,unem1ggp=886,unt_scr,unt_susc,use_umklp,usexcnhat
 real(dp),parameter :: OMEGAERMAX=100.0/Ha_eV
 real(dp) :: boxcut,boxcutc,compch_fft,compch_sph,diecut_eff_dum,domegareal,e0
 real(dp) :: ecore,ecut_eff,ecutdg_eff,ecuteps_pG0,efermi,epaw,epawdc
 real(dp) :: gsqcutc_eff,gsqcutf_eff,omegaplasma,ucvol
 logical,parameter :: prior_to_56=.TRUE.  ! For the time being do not touch
 logical :: found,ltest,only_one_kpt,qeq0,read_occupied,remove_inv
 logical :: update_energies=.FALSE.
 character(len=10) :: string
 character(len=500) :: msg
 character(len=80) :: bar
 character(len=fnlen) :: filnam,fname_scr,fname_susc
 type(BZ_mesh_type) :: Kmesh,Qmesh
 type(Coulombian_type) :: Vcp
 type(Crystal_structure) :: Cryst
 type(Epsilonm1_parameters) :: Ep
 type(Gpairs_type) :: Gpairs_q
 type(Gvectors_type) :: Gsph_epsG0,Gsph_wfn
 type(Hdr_type) :: Hdr_kss
 type(MPI_type) :: MPI_enreg_seq
 type(Pawfgr_type) :: Pawfgr
 type(Wavefunctions_information) :: Wf,Wf_val
!arrays
 integer :: ibocc(Dtset%nsppol),ng0sh_opt(3),ngfft_gw(18),ngfftc(18),ngfftf(18)
 integer,pointer :: gvec_p(:,:)
 integer,allocatable :: dimlmn(:),gvec(:,:),igfftf(:),irottb(:,:),irottbf(:,:)
 integer,allocatable :: istart(:),istop(:),ktabr(:,:),l_size_atm(:),nband_(:)
 integer,allocatable :: nbv(:),nlmn(:)
 integer,allocatable,target :: igfft(:,:,:,:)
 integer,pointer :: igfft0(:),trec(:,:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),k0(3),qtmp(3),rm1t(3),rmet(3,3),rprimd(3,3)
 real(dp) :: tsec(2)
 real(dp),allocatable :: chi0sumrule(:),en_qp(:,:,:),gw_energy(:,:,:),kibz(:,:)
 real(dp),allocatable :: ks_energy(:,:,:),nhat(:,:),nhatgr(:,:,:),occ(:,:,:)
 real(dp),allocatable :: occ_val(:,:,:),pawrhox_spl(:,:,:,:,:),qptdm(:,:)
 real(dp),allocatable :: qsmall(:,:),qsmall_dum(:,:),rhog(:,:),rhor(:,:),rhor_p(:,:)
 real(dp),allocatable :: vkb(:,:,:,:),vkbd(:,:,:,:),vkbsign(:,:),z(:),zw(:)
 real(dp),pointer :: ttns(:,:)
 complex(gwpc),allocatable :: chitmp(:,:),kxc(:,:)
 complex(gwpc),allocatable :: lwing_dum(:,:,:),uwing_dum(:,:,:)
 complex(gwpc),allocatable :: lwing(:,:,:),uwing(:,:,:)
 complex(dpc),allocatable :: m_lda_to_qp(:,:,:,:)
 complex(gwpc),allocatable,target :: chi0(:,:,:)
 complex(gwpc),pointer :: epsm1(:,:,:),vc_sqrt(:)
 logical,allocatable :: mask(:)
 character(len=80) :: title(2)
 character(len=fnlen) :: tmpfil(7)
 type(Cprj_type),allocatable :: Cprj_bz(:,:),Cprj_ibz(:,:)
 type(Little_group),allocatable :: Ltg_q(:)
 type(Paw_an_type),allocatable :: Paw_an(:)
 type(Paw_ij_type),allocatable :: Paw_ij(:)
 type(Pawfgrtab_type),allocatable :: Pawfgrtab(:)
 type(Pawrhoij_type),allocatable :: Pawrhoij(:)

#undef HAVE_GW_CUTOFF 
#if defined HAVE_GW_CUTOFF
 ! /** Variables added for cutoffed matrix elements **/
 real(dp) :: width, z0
 integer :: direction 
 !integer, allocatable :: igjell(:,:)
#endif

!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' screening : enter '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 call pclock(0)        
 call timab(301,1,tsec) ! overall time
 call timab(302,1,tsec) ! screening(1)

 write(msg,'(7a)')&
& ' SCREENING: Calculation of the susceptibility and dielectric matrices ',ch10,ch10,&
& ' Based on a program developped by R.W. Godby, V. Olevano, G. Onida, and L. Reining.',ch10,&
& ' Incorporated in ABINIT by V. Olevano, G.-M. Rignanese, and M. Torrent.',ch10
 call wrtout(ab_out,msg,'COLL') 
 call wrtout(std_out,msg,'COLL')

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
!=== Initialize MPI related quantities and parallelization level ===
!gwpara 0 --> sequential run
!1 --> parallelism over k-points 
!2 --> parallelism over bands: each proc has both fully and partially 
!occupied states while conduction bands are divided 
!
 call xcomm_init  (MPI_enreg,spaceComm)  
 call xmaster_init(MPI_enreg,master   ) 
 call xme_init    (MPI_enreg,rank     )          
 call xproc_max(nprocs,ierr)

 if (nprocs==1) Dtset%gwpara=0  
 MPI_enreg%gwpara     = Dtset%gwpara  
 MPI_enreg%parareel   = 0  
 MPI_enreg%paralbd    = 0
 MPI_enreg%paral_level= 2 ! This means k-points but it is not used
 MPI_enreg%me_fft     = 0 
 MPI_enreg%nproc_fft  = 1
!* Fake MPI_type for sequential part
 call initmpi_seq(MPI_enreg_seq) ; MPI_enreg_seq%nproc_fft=1 ; MPI_enreg_seq%me_fft=0
!
 ltest=(Dtset%mkmem/=0)
 call assert(ltest,'Option mkmem=0 not implemented',__FILE__,__LINE__)
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
!localrdwf==1 --> every processor has access to files (default)
!localrdwf==0 --> only master has access to files
!
 tmpfil(1)=TRIM(Dtfil%filnam_ds(5))//'_WF1'
 tmpfil(2)=TRIM(Dtfil%filnam_ds(5))//'_WF2'
 tmpfil(3)=TRIM(Dtfil%filnam_ds(5))//'_KG'
 tmpfil(4)=TRIM(Dtfil%filnam_ds(5))//'_DUM'
 tmpfil(6)=TRIM(Dtfil%filnam_ds(5))//'_YLM'
 tmpfil(7)=TRIM(Dtfil%filnam_ds(5))//'_PAW'
!* Parallel case: the index of the processor must be appended
 if (MPI_enreg%paral_compil_kpt==1) then
  call int2char4(MPI_enreg%me,tag) 
  jj=1 ; if (MPI_enreg%paral_compil_mpio==1 .and. Dtset%accesswff==1) jj=3
  do ii=jj,7 
   tmpfil(ii)=TRIM(tmpfil(ii))//'_P-'//tag 
  end do
 end if
!
!=== Nullify the pointers in the data types ===
 call nullify_epsilonm1_parameters(Ep) 
 call nullify_Gpairs_type(Gpairs_q)
 call nullify_Crystal_structure(Cryst)
!
!=== Get dimensional primitive translations rprimd (from input file), gprimd, metrics.. ===
 call mkrdim(acell,rprim,rprimd) 
 call metric(gmet,gprimd,ab_out,rmet,rprimd,ucvol)
!
!Define FFT grid(s) sizes (be careful!, these quantities refer to the ground-state part of the code. 
!The FFT mesh used for GW is defined in setmesh.F90. See also NOTES in the comments at the beginning of this file.
!FIXME At the moment the mesh is defined according to ecut but ecutwfn should be used
 k0(:)=zero
 call init_pawfgr(Dtset,k0,gmet,Pawfgr,mgfftf,nfftf,ecut_eff,ecutdg_eff,gsqcutc_eff,gsqcutf_eff,ngfftc,ngfftf)

 nfftc_tot=ngfftc(1)*ngfftc(2)*ngfftc(3)
 nfftf_tot=ngfftf(1)*ngfftf(2)*ngfftf(3)
!
!=== Open and read pseudopotential files ===
 call status(0,Dtfil%filstat,iexit,level,'call pspini   ')
 call pspini(Dtset,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,level,Pawrad,Pawtab,Psps,rprimd)
 if (psp_gencond==1) call print_psps(Psps,std_out,0,'COLL')
!
!=== Set up basic parameters of the calculation ===
 Ep%gwcalctyp=Dtset%gwcalctyp
 Ep%plasmon_pole_model   =.TRUE.  ; Ep%static             =.FALSE. 
 Ep%analytic_continuation=.FALSE. ; Ep%contour_deformation=.FALSE.

 mod10=MOD(Ep%gwcalctyp,10)
 if (mod10/=0.and.mod10/=8)            Ep%plasmon_pole_model=.FALSE.
 if (mod10==1)                         Ep%analytic_continuation=.TRUE.
 if (mod10==2.or.mod10==9)             Ep%contour_deformation=.TRUE.
 if (mod10==5.or.mod10==6.or.mod10==7) Ep%static=.TRUE.

 Ep%rpa  =.TRUE.  ; Ep%testparticle=.TRUE. 
 Ep%tddft=.FALSE. ; Ep%testelectron=.FALSE.
 Ep%nbnds  =Dtset%nband(1) 
 Ep%symchi =Dtset%symchi
 Ep%inclvkb=Dtset%inclvkb ; if (Dtset%usepaw/=0) Ep%inclvkb=0
 Ep%zcut   =Dtset%zcut   

 write(msg,'(2a,i4,2a,f10.6,a)')ch10,&
& ' GW calculation type              = ',Ep%gwcalctyp,ch10,&
& ' zcut to avoid poles in chi0 [eV] = ',Ep%zcut*Ha_eV,ch10
 call wrtout(std_out,msg,'COLL')
!
!=== Define consistently npw, nsh, and ecut for wavefunctions and dielectric matrix ===
 call setshells(Dtset%ecutwfn,Dtset%npwwfn,Dtset%nshwfn,Dtset%nsym,gmet,gprimd,Dtset%symrel,'wfn',ucvol)
 call setshells(Dtset%ecuteps,Dtset%npweps,Dtset%nsheps,Dtset%nsym,gmet,gprimd,Dtset%symrel,'eps',ucvol)

 Ep%npwe=Dtset%npweps  ; Ep%npwwfn=Dtset%npwwfn ; Ep%npwvec=MAX(Ep%npwe,Ep%npwwfn)
!
!=== Read parameters from KSS and verifify them === 
!* Get also ibocc (Max occ band index for each spin). 
!* For metals use a threshold after which bands are considered empty
!TODO in case of SCGW vale and conduction has to be recalculated to avoid errors 
!if a metal becomes semiconductor or viceversa.
 Ep%awtr=Dtset%awtr 
 call testlda(Dtset,Dtfil,nsym_kss,nbnds_kss,ng_kss,mpsang,gvec_p,Hdr_kss,MPI_enreg,ibocc)
 deallocate(gvec_p)

!Copy important dimension from header
 Ep%nsppol=Hdr_kss%nsppol
 Ep%nkibz =Hdr_kss%nkpt

 call hdr_vs_dtset(Hdr_kss,Dtset) 
 remove_inv=(nsym_kss/=Hdr_kss%nsym) 

 timrev = 2 ! This information is not reported in the header
!1 --> do not use time-reversal symmetry 
!2 --> take advantage of time-reversal symmetry
 call init_Crystal_from_Hdr(Cryst,Hdr_kss,timrev,remove_inv)
 call print_Crystal_structure(Cryst,mode_paral='PERS')

 if (ANY(ABS(Cryst%xred-xred)>tol6)) STOP "BUG xred"
 if (ANY(ABS(rprimd-Cryst%rprimd)>tol6)) STOP "BUG rprimd"

 allocate(kibz(3,Ep%nkibz))
 kibz(:,:)=Hdr_kss%kptns(:,:)
!END COPY

 if (Ep%npwvec>ng_kss) then
  Ep%npwvec=ng_kss
  if (Ep%npwwfn> ng_kss) Ep%npwwfn=ng_kss
  if (Ep%npwe  > ng_kss) Ep%npwe  =ng_kss
  write(msg,'(7a,3(a,i6,a))')ch10,&
&  ' screening: WARNING - ',ch10,&
&  '  Number of G-vectors found less then required. ',ch10,&
&  '  Calculation will proceed with ',ch10,&
&  '   npwvec = ',Ep%npwvec,ch10,&
&  '   npweps = ',Ep%npwe  ,ch10,&
&  '   npwwfn = ',Ep%npwwfn,ch10
  call wrtout(std_out,msg,'COLL')
 end if
 if (Ep%nbnds>nbnds_kss) then
  Ep%nbnds=nbnds_kss
  write(msg,'(6a,i4,a)')ch10,&
&  ' screening: WARNING -',ch10,&
&  '  Number of bands found less then required. ',ch10,&
&  '  Calculation will proceed with nbnds = ',nbnds_kss,ch10
  call wrtout(std_out,msg,'COLL')
 end if
!
!RECREATE AN "ABINIT ENVIRONMENT"
 allocate(nband_(Ep%nkibz*Ep%nsppol)) ; nband_(:)=Ep%nbnds
!
!If read_occupied is true, valence and partially occupied are stored in a different array 
!This method is mandatory in case of spectral method and/or gwpara==2
!If awtr==1 we evaluate the Adler-Wiser expression taking advantage of time-reversal (speed-up~2)
 read_occupied=.FALSE. ; nbvw=0
 if (Cryst%timrev==2.and.(Ep%awtr==1.or.MPI_enreg%gwpara==2.or.Dtset%spmeth/=0)) then 
  read_occupied=.TRUE. ; nbvw=MAXVAL(ibocc) ; nbcw=Ep%nbnds-nbvw
! These are tha indeces used to allocate valence and conduction states
! nbvw = Maximum number of fully/partially occupied states considering the spin 
! nbcw = Maximum number of unoccupied states considering the spin 
! TODO Here for semiconducting systems we have to be sure that each processor has all the 
! states considered in the SCGW, moreover nbsc<nbvw
  write(msg,'(4a,i5,2a,i5,a)')ch10,&
  ' screening : taking advantage of time-reversal symmetry ',ch10,&
&  ' maximum band index for partially occupied states nbvw = ',nbvw,ch10,&
&  ' remaining bands to be divided among processors   nbcw = ',nbcw,ch10
  call wrtout(std_out,msg,'PERS')
 end if 

 my_minb=1 ; my_maxb=Ep%nbnds ; my_nbnds=my_maxb-my_minb+1
 allocate(MPI_enreg%proc_distrb(Ep%nkibz,Ep%nbnds,Ep%nsppol))  
 MPI_enreg%proc_distrb(:,:,:)=rank

 if (MPI_enreg%gwpara==2) then 
  write(msg,'(2a)')ch10,&
&  ' loop over bands done in parallel (assuming time reversal!)'
  call wrtout(std_out,msg,'COLL')
! write(*,*)'rank ',rank,' has ',Ep%nbnds,Ep%nkibz,Ep%nsppol,nbvw,nbcw
  MPI_enreg%proc_distrb(:,:,:)=-999
  allocate(istart(nprocs),istop(nprocs))
! Divide conduction states or divide all states in the case of completeness trick
! MG FIXME this part has to be checked, conflict due to merge with FB
  if (Dtset%gwcomp==0) then
!  * Each proc has fully and partially occupied states 
   MPI_enreg%proc_distrb(:,1:nbvw,:)=rank
   call split_work2(nbcw,nprocs,istart,istop)
   my_minb=nbvw+istart(rank+1) ; my_maxb=nbvw+istop(rank+1)
   do irank=0,nprocs-1
    MPI_enreg%proc_distrb(:,nbvw+istart(irank+1):nbvw+istop(irank+1),:)=irank
   end do 
  else
   call split_work2(Ep%nbnds,nprocs,istart,istop)
   my_minb=istart(rank+1) ; my_maxb=istop(rank+1)
   do irank=0,nprocs-1
    MPI_enreg%proc_distrb(:,istart(irank+1):istop(irank+1),:)=irank
   end do 
  end if
  my_nbnds=my_maxb-my_minb+1 
  ltest=(my_nbnds>=1) 
  write(msg,'(4a,2(i4,a),a)')ch10,&
&  ' One or more processors has zero number of bands ',ch10,&
&  ' my_minb = ',my_minb,' my_maxb = ',my_maxb,ch10,&
&  ' This is a waste, decrease the number of processors '
  call assert(ltest,msg,__FILE__,__LINE__)
  deallocate(istart,istop)
! === Announce the treatment of bands by each node ===
  do irank=0,nprocs-1
   if (irank==rank) then 
    write(msg,'(4(a,i4))')' treating ',my_nbnds,&
&    ' bands from ',my_minb,' up to ',my_maxb,' by node ',irank
    call wrtout(std_out,msg,'PERS')
   end if 
  end do 
 end if 
!
!=== Allocate basic arrays and electronic structure ===
 allocate(gvec(3,Ep%npwvec))
 allocate(gw_energy(Ep%nkibz,Ep%nbnds,Ep%nsppol)) ; gw_energy(:,:,:)=zero 
 allocate(ks_energy(Ep%nkibz,Ep%nbnds,Ep%nsppol)) ; ks_energy(:,:,:)=zero 
 allocate(     occ(Ep%nkibz,Ep%nbnds,Ep%nsppol))  ;       occ(:,:,:)=zero

 call nullify_wf_info(Wf)
 call init_wf_info_1(Wf,Dtset%gwmem,Dtset%paral_kgb,Ep%npwwfn,&
& my_minb,my_maxb,Ep%nkibz,Dtset%nsppol,Dtset%nspden,Dtset%nspinor)
!
!=== Allocate Kleynmann-Bylander form factors, derivatives and sign ===
 allocate(vkb (Ep%npwwfn,Cryst%ntypat,mpsang,Ep%nkibz),STAT=istat)
 if (istat/=0) call memerr(FILE__,'vkb',Ep%npwwfn*Cryst%ntypat*mpsang*Ep%nkibz,'spc')
 allocate(vkbd(Ep%npwwfn,Cryst%ntypat,mpsang,Ep%nkibz),STAT=istat)
 if (istat/=0) call memerr(FILE__,'vkbd',Ep%npwwfn*Cryst%ntypat*mpsang*Ep%nkibz,'spc')
 allocate(vkbsign(mpsang,Cryst%ntypat))
!
!=================================================================
!=== Read KS band structure and information about the material ===
!=================================================================
! * accesswf defines Fortran-IO or ETSF-IO

 allocate(Cprj_ibz(Cryst%natom,Dtset%nspinor*Ep%nbnds*Ep%nkibz*Ep%nsppol*Dtset%usepaw))

 if (Dtset%usepaw==1) then 
  allocate(dimlmn(Cryst%natom)) 
  do iat=1,Cryst%natom
   dimlmn(iat)=Pawtab(Cryst%typat(iat))%lmn_size
  end do
  call cprj_alloc(Cprj_ibz,0,dimlmn)

  allocate(Pawrhoij(Cryst%natom),nlmn(Cryst%ntypat))
  do itypat=1,Cryst%ntypat
   nlmn(itypat)=Pawtab(itypat)%lmn_size
  end do
  call rhoij_alloc(1,nlmn,Dtset%nspden,Dtset%nsppol,Pawrhoij,Cryst%typat)
  deallocate(nlmn)
 end if

 allocate(trec(3,3,Cryst%nsym),ttns(3,Cryst%nsym))

 if (read_occupied) then 
! * Store occupied and empty states in two different arrays to use faster equation based 
! on time reversal symmetry or spectral method. Note that if nsppol==2 than wfg_val 
! might also contain unoccupied states, this case is treated inside cchi0 and cchi0q0
! TODO gwpara==2 should be the default method in parallel executions
  call nullify_wf_info(Wf_val)
  call init_wf_info_1(Wf_val,Dtset%gwmem,Dtset%paral_kgb,Ep%npwwfn,&
&  1,nbvw,Ep%nkibz,Ep%nsppol,Dtset%nspden,Dtset%nspinor)

  call rdkss(Dtfil,Dtset,Pawtab,Cryst%nsym,Ep%nbnds,nbvw,Ep%nkibz,Ep%npwvec,Dtset%nspinor,Ep%nsppol,Ep%npwwfn,title,trec,&
&  gvec,ks_energy,occ,Wf%wfg,Cprj_ibz,Cryst%ntypat,Cryst%natom,mpsang,ttns,vkbsign,vkb,vkbd,nel,MPI_enreg,my_minb,my_maxb,& 
&  wf_val=Wf_val%wfg) ! Optional argument
 else 
! * Old implementation: store all the wavefunctions in a unique arrays
! * It can be used if the system does not present time reversal symmetry 
  call rdkss(Dtfil,Dtset,Pawtab,Cryst%nsym,Ep%nbnds,nbvw,Ep%nkibz,Ep%npwvec,Dtset%nspinor,Ep%nsppol,Ep%npwwfn,title,trec,&
&  gvec,ks_energy,occ,Wf%wfg,Cprj_ibz,Cryst%ntypat,Cryst%natom,mpsang,ttns,vkbsign,vkb,vkbd,nel,MPI_enreg,my_minb,my_maxb) 
 end if 

 ltest=ALL(trec==Cryst%symrec)      
 call assert(ltest,'BUG in trec',__FILE__,__LINE__)
 ltest=ALL((ttns-Cryst%tnons)<tol6) 
 call assert(ltest,'BUG in ttns',__FILE__,__LINE__)  
 if (Cryst%nsym/=Dtset%nsym .and. Dtset%usepaw==1) stop "nsym /= Dtset%nsym, check pawinit and symrhoij"
 deallocate(trec,ttns)

 ltest=(ALL(ABS(Cryst%rprimd-rprimd)<tol6))
 call assert(ltest,'Unit cell mismatch in Cryst%',__FILE__,__LINE__)
! 
! ==========================
! === PAW initialization ===
! ==========================
 if (Dtset%usepaw==1) then 
  call chkpawovlp(Cryst%natom,Cryst%ntypat,Dtset%pawovlp,Pawtab,Cryst%rmet,Cryst%typat,Cryst%xred)
!  
! === Initialize values for several basic arrays ===
! if (psp_gencond==1) then
  call timab(553,1,tsec)
! TODO this part should be done at the very beginning, otherwise test on orthonormalization in rdkss fails!
! Check pawxcdev>2 since gaunt coefficients are allocated with different size 
  diecut_eff_dum=ABS(Dtset%diecut)*Dtset%dilatmx**2
  call pawinit(diecut_eff_dum,Psps%indlmn,Dtset%pawlcutd,Dtset%pawlmix,Psps%lmnmax,Psps%mpsang,Psps%n1xccc,&
&  Dtset%pawnphi,Cryst%nsym,Dtset%pawntheta,Cryst%ntypat,Pawang,Pawrad,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev)
  call timab(553,2,tsec)
! end if
! === Evaluate <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j> for long wavelength limit ===
  call pawnabla_init(Psps%mpsang,Psps%lmnmax,Cryst%ntypat,Psps%indlmn,Pawrad,Pawtab)
! FIXME see above. Note that here nsym_old comes from KSS.
! if (psp_gencond==1.or.nsym_old/=Cryst%nsym) then
  call setsymrhoij(gprimd,Pawang%l_max-1,Cryst%nsym,Dtset%pawprtvol,rprimd,Dtset%symafm,Cryst%symrec,Pawang%zarot)
  nsym_old=Cryst%nsym
! end if
! === Initialize and compute data for LDA+U ===
  if (Dtset%usepawu>0.or.Dtset%useexexch>0) then
   call pawpuxinit(Dtset%dmatpuopt,Dtset%exchmix,Dtset%jpawu,Dtset%lexexch,Dtset%lpawu,&
&   Psps%indlmn,Psps%lmnmax,Cryst%ntypat,Pawang,Dtset%pawprtvol,Pawrad,Pawtab,Dtset%upawu,&
&   Dtset%useexexch,Dtset%usepawu)
  end if
! === Eventually open temporary storage file ===
! FIXME also mkmem_ not yet defined
! if (mkmem_==0) then
!  open(Dtfil%unpaw,file=tmpfil(7),form='unformatted',status='unknown')
!  rewind(unit=Dtfil%unpaw)
! end if

! === Get Pawrhoij from the header of the KSS file ===
  call rhoij_copy(Hdr_kss%Pawrhoij,Pawrhoij)

! === Re-symmetrize symrhoij ===
! choice=1 ; optrhoij=1 
! call symrhoij(choice,Psps%indlmn,Cryst%indsym,Psps%lmnmax,Cryst%natom,Cryst%nsym,Cryst%ntypat,&
! & optrhoij,Pawang,Dtset%pawprtvol,Pawrhoij,Dtset%symafm,Cryst%symrec,Cryst%typat)
!  
! === Evaluate form factors for radial part of phi.phj-tphi.tphj ===
  rhoxsp_method=2 
  if (rhoxsp_method==1) then ! Arnaud-Alouani
   dim2_rhox=Psps%lnmax*(Psps%lnmax+1)/2 
   dim1_rhox=2*(Psps%mpsang-1) 
  else if (rhoxsp_method==2) then  ! Shiskin-Kresse
   dim1_rhox=MAXVAL(Pawtab(:)%l_size)**2
   dim2_rhox=MAXVAL(pawtab(:)%lmn2_size) 
  end if
  allocate(pawrhox_spl(Psps%mqgrid_ff,2,0:dim1_rhox,dim2_rhox,Cryst%ntypat))  
  call paw_mkrhox_spl(Cryst%ntypat,Psps,Pawrad,Pawtab,pawrhox_spl,rhoxsp_method,dim1_rhox,dim2_rhox)
!  
! === Variables/arrays related to the fine FFT grid ===
  allocate(nhat(nfftf,Dtset%nspden)) 
  nhat(:,:)=zero
  allocate(Pawfgrtab(Cryst%natom),l_size_atm(Cryst%natom)) 
  do iat=1,Cryst%natom 
   l_size_atm(iat)=Pawtab(Cryst%typat(iat))%l_size 
  end do
  call nullify_pawfgrtab(Pawfgrtab) 
  call    init_pawfgrtab(Pawfgrtab,l_size_atm) 
  deallocate(l_size_atm)
  compch_fft=greatest_real
  usexcnhat=MAXVAL(Pawtab(:)%vlocopt)
! * 0 --> Vloc in atomic data is Vbare    (Blochl s formulation)
! * 1 --> Vloc in atomic data is VH(tnzc) (Kresse s formulation)
  write(*,*)' screening : using usexcnhat = ',usexcnhat 
!  
! === Identify parts of the rectangular grid where the density has to be calculated ===
  call status(0,Dtfil%filstat,iexit,level,'call nhatgrid ')
  optgr0=Dtset%pawstgylm ; optgr1=0 ; optgr2=0 ; optrad=1-Dtset%pawstgylm
  if (Dtset%xclevel==2.and.usexcnhat>0) optgr1=Dtset%pawstgylm 

  call nhatgrid(Cryst%atindx1,gmet,MPI_enreg_seq,Cryst%natom,Cryst%nattyp,nfftf,ngfftf,Cryst%ntypat,&
&  optgr0,optgr1,optgr2,optrad,Pawfgrtab,Pawtab,rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)
 else
!  
! FB: We need to allocate these arrays anyway, since they are passed to subroutines
  allocate(pawrhox_spl(1,1,1,1,0))  
  allocate(nhat(1,0)) 
 end if !End of PAW initialization
 call pclock(1)
!
!=== Create structure describing k-point sampling ===
!TODO Kmesh%bz should be [-half,half[ but this modification will be painful!
 call setup_Kmesh(Ep%nkibz,kibz,Cryst,Kmesh,Dtset%prtvol)
!
!=== Find Q-mesh, and do setup for long wavelength limit ===
!For the moment hard-coded values but they should be passed through Dtset%
 nqsm=1 ; allocate(qsmall(3,nqsm)) ; qsmall(1,1)=0.000010 ; qsmall(2,1)=0.000020 ; qsmall(3,1)=0.000030
 call find_Qmesh(Cryst,gprimd,Kmesh,nqsm,qsmall,Qmesh,Dtset%prtvol)
 deallocate(qsmall)
 Ep%nqlwl=nqsm ; allocate(Ep%qlwl(3,Ep%nqlwl)) ; Ep%qlwl(:,:)=Qmesh%small  !? duplication
!
!=== Find optimal value for G-sphere enlargment due to oscillator matrix elements ===
 mg0sh=25
 call get_ng0sh(Kmesh%nbz,Kmesh%bz,Qmesh%nibz,Qmesh%ibz,Kmesh%nbz,Kmesh%bz,gmet,GW_TOLQ0,mg0sh,ng0sh_opt)
 Ep%mG0(:)=ng0sh_opt(:)
 
 Ep%npwepG0=Ep%npwe
 allocate(Ltg_q(Qmesh%nibz))
 do iq=1,Qmesh%nibz 
  call nullify_little_group(Ltg_q(iq)) 
 end do

 if (Ep%symchi/=0) then
! === In case of symmetrization, find little group ===
! * In long-wavelength limit we consider a small but finite q. However the oscillators are 
! evaluated setting q==0. Thus it is possible to take advantage of symmetries also when q --> 0.
  do iq=1,Qmesh%nibz
!  TODO Switch on use_umklp
   qtmp(:)=Qmesh%ibz(:,iq) ; if (normv(qtmp,gmet,'G')<GW_TOLQ0) qtmp(:)=zero ; use_umklp=0
   call setup_little_group(qtmp,Kmesh,gmet,Cryst,Ep%npwvec,gvec,Ep%npwe,&
&   use_umklp,Dtset%prtvol,Ltg_q(iq))
  end do

  ecuteps_pG0=MAXVAL(Ltg_q(:)%max_kin_gmG0)+tol6 ; npwe_pG0=0 ; nsheps_pG0=0
  write(*,*)' Due to umklapp processes : ecuteps_pg0=',ecuteps_pG0
  call setshells(ecuteps_pG0,npwe_pG0,nsheps_pG0,Cryst%nsym,gmet,gprimd,Cryst%symrel,'eps_pG0',Cryst%ucvol)
  Ep%npwepG0=npwe_pG0
 end if

 ltest=(Ep%npwepG0<=Ep%npwvec)
 write(msg,'(4a,i5,a,i5)')ch10,&
& ' npwepG0 > npwvec, decrease npweps or increase npwwfn. ',ch10,&
& ' npwepG0 = ',Ep%npwepG0,' npwvec = ',Ep%npwvec
 call assert(ltest,msg,__FILE__,__LINE__)
 call pclock(2)
!
!=== Setup of the FFT mesh used for oscillator matrix elements ===
!ngfft_gw(7:18) is the same as Dtset%ngfft(7:18), initialized before entering screening. 
!Here we just redefine ngfft_gw(1:6) according to the following options :
!
!method==0 ==> FFT grid read from __fft.in__ (only for debugging purpose)
!method==1 ==> normal FFT grid 
!method==2 ==> slightly augmented FFT grid to calculate exactly rho_tw_g (see setmesh.F90)
!method==3 ==> doubled FFT grid, to treat exactly the convolution defining the density,
!useful in sigma if ppmodel=[2,3,4] since rho(G-Gp) or to calculate matrix elements of v_Hxc.
!
!enforce_sym==1 ==> enforce a direct space FFT mesh compatible with all symmetries operation
!enforce_sym==0 ==> Find the smallest FFT grid compatibile with the library, do not care about symmetries
!
 ngfft_gw(1:18)=Dtset%ngfft(1:18) ; method=2
 if (Dtset%fftgw==00 .or. Dtset%fftgw==01) method=0
 if (Dtset%fftgw==10 .or. Dtset%fftgw==11) method=1
 if (Dtset%fftgw==20 .or. Dtset%fftgw==21) method=2
 if (Dtset%fftgw==30 .or. Dtset%fftgw==31) method=3
 enforce_sym=MOD(Dtset%fftgw,10) 

!call setmesh(gmet,gvec,Dtset%ngfft,Ep%npwvec,Ep%npwe,Ep%npwwfn,nfftgw_tot,method,Ep%mG0,Cryst,enforce_sym)
 call setmesh(gmet,gvec,ngfft_gw,Ep%npwvec,Ep%npwepG0,Ep%npwwfn,nfftgw_tot,method,Ep%mG0,Cryst,enforce_sym)

 nfftgw=nfftgw_tot ; mgfftgw=MAXVAL(ngfft_gw(1:3))
!
!=== Calculate igfft table for FFT and Ep%mG0 defines the number of shells to be added ===
!Ep%mG0 gives, for each reduced direction, the max G0 component to account for umklapp processes
!igfft0 contains the FFT index of the G"s used by fourdp.
 nG01d=2*Ep%mG0(1)+1 ; nG02d=2*Ep%mG0(2)+1 ; nG03d=2*Ep%mG0(3)+1
 allocate(igfft(Ep%npwvec,nG01d,nG02d,nG03d))

 call cigfft(Ep%mG0,Ep%npwvec,ngfft_gw,gvec,igfft)
 igfft0 => igfft(:,Ep%mG0(1)+1,Ep%mG0(2)+1,Ep%mG0(3)+1) 
 call pclock(3)
!
!=== Set up table indicating rotations of r-points (i.e $R^{-1}(r-\tau)$) and G vectors ===
!This part might be CPU consuming in isolated systems.
 allocate(irottb(nfftgw,Cryst%nsym))
 call setup_FFT_rotation(Cryst%nsym,Cryst%symrec,Cryst%tnons,nfftgw,ngfft_gw,irottb)
!£££ call FFT_rotations(Cryst,ngfft,irottb)

!=== Create structure describing the G-sphere used for chi0/espilon and Wfns ===
!* The cutoff is >= ecuteps to allow for umklapp
!MG, Fabien, I modified setup_G_rotation to speed up the loops, 
!sincerely I dont like only_one_kpt because some entries in Gsphere are not correctly filled!
 only_one_kpt=(Kmesh%nbz==1) 
 call nullify_Gvectors(Gsph_epsG0)
 call init_Gvectors_type(only_one_kpt,Gsph_epsG0,Cryst,Ep%npwepG0,gvec,gmet,gprimd)
 call nullify_Gvectors(Gsph_wfn)
 call init_Gvectors_type(only_one_kpt,Gsph_wfn  ,Cryst,Ep%npwvec ,gvec,gmet,gprimd)
!TODO deallocate gvec, everything should be passed through objects
!
!£££ TODO For the moment fine mesh is equal to GW mesh, just to check PPs implementation
 if (prior_to_56 .and. Dtset%usepaw/=1) then
  ngfftf(:)=ngfft_gw(:)
  nfftf_tot=ngfftf(1)*ngfftf(2)*ngfftf(3)
  nfftf=nfftf_tot 
  mgfftf=mgfftgw
  write(msg,'(a)')'screening: WARNING - enforcing OLD GW mesh for density '
 else 
  write(msg,'(a)')'screening: using augmented FFT mesh '
 end if
 call wrtout(std_out,msg,'COLL') 
 call print_ngfft(ngfftf)
!£££
!
!=== Setup of the coulombian interaction, a cutoff can be employed ===
 call setup_coulombian(Dtset,Gsph_epsG0,Qmesh,Kmesh,Ep%npwe,rprimd,ngfftf,MPI_enreg,Vcp)
 call pclock(4)
!
!=== Tables for the dense FFT mesh used for rhor ===
 allocate(igfftf(Ep%npwvec),mask(Ep%npwvec)) ; mask(:)=.TRUE.
 call kgindex(igfftf,gvec,mask,MPI_enreg_seq,ngfftf,Ep%npwvec)
 deallocate(mask)

 allocate(irottbf(nfftf,Cryst%nsym))
 call setup_FFT_rotation(Cryst%nsym,Cryst%symrec,Cryst%tnons,nfftf,ngfftf,irottbf)
!£££ call FFT_rotations(Cryst,ngfftf,irottbf)
!
!=== FFT index of R^-1 (r-t) to symmetrize u_Sk  ===
 allocate(ktabr(nfftgw,Kmesh%nbz))
 do ik_bz=1,Kmesh%nbz
  isym=Kmesh%tabo(ik_bz)
  do ifft=1,nfftgw
   ktabr(ifft,ik_bz)=irottb(ifft,isym)
  end do
 end do

 allocate(Cprj_bz(Cryst%natom,Dtset%nspinor*Ep%nbnds*Kmesh%nbz*Ep%nsppol*Dtset%usepaw))

 if (Dtset%usepaw==1) then 
! === Symmetrize <Proj_i|Cnk> to get them in the Full Brillouin zone ===
  call cprj_alloc(Cprj_bz,0,dimlmn)

! TODO add rotation in spinor space
  call paw_symcprj(Pawtab,Cryst,Dtset%nspinor,Ep%nbnds,nband_,&
&  Dtset%nsppol,Psps,Kmesh,Cprj_ibz,Pawang,dimlmn,Cprj_bz)
 end if

 deallocate(kibz)
 call timab(302,2,tsec) ! screening(1)
!
!=== Save FFT related information in Wf% and Wf_val% ===
 call timab(303,1,tsec) ! screening(fftwfn 
 call init_wf_info_2(Wf,igfft0,ngfft_gw)
 if (read_occupied) call init_wf_info_2(Wf_val,igfft0,ngfft_gw)
 call timab(303,2,tsec)
 call pclock(5)

 call timab(304,1,tsec) ! KS => QP; [wfrg]
 if (Ep%gwcalctyp>=10) then 
! 
! === Self-consistent GW: read QP wavefunctions of the previous step ===
! * Transformation matrices, QP energies and density (FBruneval).
  allocate(rhor_p(nfftf,Dtset%nspden),en_qp(Kmesh%nibz,Ep%nbnds,Ep%nsppol))
  allocate(m_lda_to_qp(Ep%nbnds,Ep%nbnds,Kmesh%nibz,Ep%nsppol))
! TODO add the abinit header to fix the problem with the FFT grid
!      reduce memory required by m_lda_to_qp!
  call rdqps(Ep%gwcalctyp,Dtfil,Kmesh,Ep%nbnds,Dtset%nspden,Ep%nsppol,&
&  nscf,nfftf,ks_energy,nbsc,en_qp,m_lda_to_qp,rhor_p)

  ks_energy(:,:,:)=en_qp(:,:,:) ; deallocate(rhor_p)
! 
! === Skip tranformation if iteration is number 0, otherwise update only the wfg treated with GW ===
! *** TODO rewrite this part, using nbsc should be faster ***
  if (nscf/=0) then 

   if (.not.read_occupied) then 
!   All bands on each processor
    call calc_wf_qp(MPI_enreg,Kmesh%nibz,Ep%nbnds,Wf%npwwfn,Ep%nsppol,Dtset%nspinor,&
&    m_lda_to_qp,my_minb,my_maxb,1,nbsc,Wf%wfg)
    call reinit_wf_info(Wf)
   else  
!   It may be better if we add a pointer inside Wf to store the valence states!
    call calc_wf_qp_Wfval(MPI_enreg,Kmesh%nibz,Ep%nbnds,Wf%npwwfn,Ep%nsppol,Dtset%nspinor,&
&    m_lda_to_qp,my_minb,my_maxb,1,nbsc,Wf%wfg,nbvw,Wf_val%wfg)
    call reinit_wf_info(Wf) 
    call reinit_wf_info(Wf_val) 
   end if 
   if (Dtset%usepaw==1) then 
!   * Update and re-symmetrize PAW cprj in the full BZ
    call update_cprj(Cryst%natom,Kmesh%nibz,Ep%nbnds,Ep%nsppol,Dtset%nspinor,m_lda_to_qp,dimlmn,Cprj_ibz)

!   TODO add rotation in spinor space
    call paw_symcprj(Pawtab,Cryst,Dtset%nspinor,Ep%nbnds,nband_,&
&    Dtset%nsppol,Psps,Kmesh,Cprj_ibz,Pawang,dimlmn,Cprj_bz)
   end if

  end if ! nscf/=0
  deallocate(m_lda_to_qp,en_qp)
 end if !gwcalctyp >= 10
 call timab(304,2,tsec) 
 call pclock(6)
!
!=== Get occupation numbers for metals or nbv (valence band index for each spin) ===
!* fixmom is passed to fermi.F90 to fix the problem with newocc in case of magnetic metals
 call timab(305,1,tsec) ! screening(2)
 allocate(nbv(Ep%nsppol))
 call fermi(Hdr_kss,Ep%nbnds,Kmesh%nibz,Dtset%fixmom,Ep%nsppol,Kmesh%wt,ks_energy,occ,nel,nbv,efermi)
 gw_energy(:,:,:)=ks_energy(:,:,:)
!
!=== Update of the eigenvalues from external file ===
!TODO this part should be treated in a clearer way, maybe introducing an input variable, 
!moreover what happens if the user asks for both a scissor and an update of the energies? 
 inquire(file='in.gw',exist=update_energies)
 Ep%soenergy=Dtset%soenergy
 if (Ep%soenergy>0.1d-4) then
  write(msg,'(5a,f7.3,a)')&
&  ' screening : performing a first self-consistency',ch10,&
&  ' update of the energies in W by a scissor operator',ch10,&
&  ' applying a scissor operator of [eV] : ',Ep%soenergy*Ha_eV,ch10
  call wrtout(std_out,msg,'COLL') 
  call wrtout(ab_out,msg,'COLL')
  do is=1,Ep%nsppol
   do ib=nbv(is)+1,Ep%nbnds
    gw_energy(:,ib,is)=ks_energy(:,ib,is)+Ep%soenergy
   end do
  end do
 else if (update_energies) then
  write(msg,'(4a)')&
&  ' screening : performing a first self-consistency',ch10,&
&  ' update of the energies in W by a previous GW calculation',ch10
  call wrtout(std_out,msg,'COLL') 
  call wrtout(ab_out,msg,'COLL')
! TODO this sub should be cleaned and modified a bit, for ex add spin and header
  call rdgw(Kmesh%nibz,Ep%nbnds,nbv,Ep%nsppol,Kmesh%ibz,gw_energy)
 end if
!
#if defined HAVE_GW_CUTOFF
!qui bisogna mettere un check per controllare che z0 e d siano 0<= z0, d<= 1
!MG It wont work in case of gwpara==2
 if(.true.) then
  z0 = dtset%userra
  width  = dtset%userrb
  direction = dtset%useria
  call cutoff_m_elem(Ep,Kmesh,gvec,Cryst%symrec,Wf,ks_energy,gw_energy,z0,width,occ,direction,gprimd)

  return
 end if
#endif
!
!========================
!=== COMPUTE DENSITY ====
!========================
!
!Evaluate "smooth" part (complete charge in case of NC pseudos)
 allocate(rhor(nfftf,Dtset%nspden))
 if (read_occupied) then 
! === Valence states are stored on each proc ===
! * note use_MPI=.FALSE. ===
  call calc_density(Wf_val,Cryst,irottbf,Ep%nbnds,ngfftf,nfftf,igfftf,&
&  occ,rhor,Kmesh,MPI_enreg,tim_fourdp,.FALSE.)
 else 
! === Each proc has the full set of WFR ===
  call calc_density(Wf,Cryst,irottbf,Ep%nbnds,ngfftf,nfftf,igfftf,&
&  occ,rhor,Kmesh,MPI_enreg,tim_fourdp,.FALSE.)
 end if

!TODO this has to be done in a better way, moreover wont work for PAW
 call cutoff_density(ngfftf,Dtset%nspden,Dtset%nsppol,Vcp,rhor,MPI_enreg)
!
!=== Only for PAW: treat onsite contributions ===
 if (Dtset%usepaw==1) then 

! === Compensation charge to be added to the density coming from smooth WFs ===
  nhatgrdim=0 ; if (Dtset%xclevel==2) nhatgrdim=usexcnhat*Dtset%pawnhatxc ; ider=2*nhatgrdim ; izero=0
  if (nhatgrdim>0) allocate(nhatgr(nfftf,Dtset%nspden,3))

  call pawmknhat(compch_fft,ider,izero,MPI_enreg_seq,Cryst%natom,nfftf,ngfftf,nhatgrdim,Dtset%nspden,&
&  Cryst%ntypat,Dtset%paral_kgb,Pawang,Pawfgrtab,nhatgr,nhat,Pawrhoij,Pawtab,Cryst%typat,Cryst%ucvol)
! 
! === Evaluate onsite energies, potentials, densities ===
! * Initialize Variables/arrays related to the PAW spheres.
! * Initialize also lmselect (index of non-zero LM-moments of densities).
! TODO call init_paw_ij in scfcv and respfn, fix small issues

  cplex=1 ; cplex_dij=Dtset%nspinor 
  allocate(Paw_ij(Cryst%natom)) 
  call nullify_paw_ij(Paw_ij)
  call init_paw_ij(Paw_ij,cplex,cplex_dij,Dtset%nspinor,Dtset%nsppol,&
&  Dtset%nspden,Dtset%pawspnorb,Cryst%natom,Cryst%ntypat,Cryst%typat,Pawtab,has_dijhartree=1)

  allocate(Paw_an(Cryst%natom))
  call nullify_paw_an(Paw_an)
  call init_paw_an(Cryst%natom,Cryst%ntypat,Dtset%nspden,cplex,Dtset%pawxcdev,&
&  Dtset%pawspnorb,Cryst%typat,Pawang,Pawtab,Paw_an,has_vxcval=0)

  call status(0,Dtfil%filstat,iexit,level,'call pawdenpot')

  nzlmopt=-1 ; option=0 ; compch_sph=greatest_real
  call pawdenpot(compch_sph,epaw,epawdc,Dtset%ixc,Cryst%natom,Dtset%nspden,&
&  Cryst%ntypat,nzlmopt,option,Paw_an,Paw_ij,Pawang,Dtset%pawprtvol,Pawrad,Pawrhoij,Dtset%pawspnorb,&
&  Pawtab,Dtset%pawxcdev,Cryst%typat,Dtset%xclevel,Psps%znuclpsp)
! 
! * Dij are not calculated since they are not use (call to setvtrial is needed)
 end if

 call test_charge(nfftf,Dtset%nspden,rhor,ucvol,nhat,Dtset%usepaw,&
& usexcnhat,Pawfgr%usefinegrid,compch_sph,compch_fft,omegaplasma)

!=== Add compensation charge on FFT mesh for PAW then get rho(G) ===
 if (Dtset%usepaw==1) rhor(:,:)=rhor(:,:)+nhat(:,:) 
 call prtrhomxmn(std_out,MPI_enreg,nfftf,ngfftf,Dtset%nspden,1,rhor)
 allocate(rhog(2,nfftf)) 
 call fourdp(1,rhog,rhor(:,1),-1,MPI_enreg,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp)
 call pclock(7)
!
!=== Define the frequency mesh for epsilon according to the method used ===
 Ep%nomegaei=1 
 Ep%nomegaer=1 ; if (Ep%static) Ep%nomegaer=0

 if (Ep%plasmon_pole_model.and.Dtset%nfreqre==1.and.Dtset%nfreqim==0) then
! === For ppmodels 2,3,4, only omega=0 is needed ===
  Ep%nomegaer=1 ; Ep%nomegaei=0
  write(msg,'(7a)')ch10,&
&  ' The inverse dielectric matrix will be calculated on zero frequency only',ch10,&
&  ' please note that the calculated epsilon^-1 cannot be used ',ch10,&
&  ' to calculate QP corrections using plasmonpole model 1',ch10
  call wrtout(std_out,msg,'COLL') 
  call wrtout(ab_out,msg,'COLL')
 end if 
!
!=== Max number of omega along the imaginary axis ===
 if (Ep%analytic_continuation.or.Ep%contour_deformation) then
  Ep%nomegaei=Dtset%nfreqim
  if (Ep%nomegaei==0) then 
   Ep%nomegaei=NOMEGAGAUSS
   write(msg,'(4a,i5)')ch10,&
&   ' screening : WARNING - ',ch10,&
&   ' number of imaginary frequencies set to default = ',NOMEGAGAUSS
   call wrtout(std_out,msg,'COLL')
  end if
 end if
!=== Range and total number of real frequencies ===
 if (Ep%contour_deformation) then
  Ep%nomegaer=Dtset%nfreqre ; Ep%omegaermax=Dtset%freqremax
  if (Ep%nomegaer==0) then 
   Ep%nomegaer=NOMEGAREAL
   write(msg,'(4a,i5)')ch10,&
&   ' screening : WARNING - ',ch10,&
&   ' number of real frequencies set to default = ',NOMEGAREAL
   call wrtout(std_out,msg,'COLL')
  end if 
  if (Ep%omegaermax<0.1d-4) then 
   Ep%omegaermax=OMEGAERMAX
   write(msg,'(4a,f8.4)')ch10,&
&   ' screening : WARNING - ',ch10,&
&   ' Max real frequencies set to default [Ha] = ',OMEGAERMAX
   call wrtout(std_out,msg,'COLL')
  end if 
  domegareal=Ep%omegaermax/(Ep%nomegaer-1)
 end if
!
!=== Calculate frequency mesh. First omega is always zero without broadening ===
!FIXME what about metals, I guess we should add eta, this means that we have to 
!know if the system is metallic, for example using occopt
 Ep%nomega=Ep%nomegaer+Ep%nomegaei 
 allocate(Ep%omega(Ep%nomega))
 Ep%omega(1)=czero_gw 
 do io=2,Ep%nomegaer
  Ep%omega(io)=CMPLX((io-1)*domegareal,0.0)
 end do
 if (Ep%plasmon_pole_model .and. Ep%nomega==2) then
  e0=Dtset%ppmfrq ; if (e0<0.1d-4) e0=omegaplasma
  Ep%omega(2)=CMPLX(0.0,e0)
 end if

 if (Ep%analytic_continuation) then
! === For AC use Gauss-Legendre quadrature method === 
! * To calculate \int_0^\infty dx f(x) we calculate \int_0^1 dz f(1/z - 1)/z^2.
! * Note that the grid is not log as required by CD, so we cannot use the same SCR file 
  allocate(z(Ep%nomegaei),zw(Ep%nomegaei)) ! knots and weights for AC 
  call coeffs_gausslegint(zero,one,z,zw,Ep%nomegaei)
  do io=1,Ep%nomegaei
   Ep%omega(Ep%nomegaer+io)=CMPLX(0.0,1/z(io)-1)
  end do
  deallocate(z,zw)
 else if (Ep%contour_deformation) then
  e0=Dtset%ppmfrq ; if (e0<0.1d-4) e0=omegaplasma
  do io=1,Ep%nomegaei
   Ep%omega(Ep%nomegaer+io)=CMPLX(0.0,e0/3.*(EXP(2./(Ep%nomegaei+1)*LOG(4.)*io)-1.))
  end do
 end if

 Ep%spmeth=Dtset%spmeth ; Ep%nomegasf=Dtset%nomegasf ; Ep%spsmear=Dtset%spbroad
 if (Ep%spmeth/=0) then 
  write(msg,'(2a,i3,2a,i8)')ch10,&
&  ' screening : using spectral method = ',Ep%spmeth,ch10,&
&  ' Number of frequencies for imaginary part= ',Ep%nomegasf
  call wrtout(std_out,msg,'COLL')
  if (Ep%spmeth==2) then 
   write(msg,'(a,f8.5,a)')' gaussian broadening = ',Ep%spsmear*Ha_eV,' [eV]'
   call wrtout(std_out,msg,'COLL')
  end if  
 end if
!
!== Print frequency mesh for chio ==
 write(msg,'(2a)')ch10,' calculating chi0 at frequencies [eV] :'
 call wrtout(std_out,msg,'COLL') 
 call wrtout(ab_out,msg,'COLL')
 do io=1,Ep%nomega
  write(msg,'(i3,2es16.6)')io,Ep%omega(io)*Ha_eV
  call wrtout(std_out,msg,'COLL') 
  call wrtout(ab_out,msg,'COLL')
 end do
 call pclock(8)

 allocate(kxc(0,0))
 if (Ep%tddft) then
! MG TODO this has to be rewritten
  allocate(chitmp(Ep%npwe,Ep%npwe),STAT=istat) ; if (istat/=0) call memerr(FILE__,'chitmp',Ep%npwe**2,'gwpc')
  ltest=(Ep%nsppol==1)
  call assert(ltest,'TDDFT CHI and nsppol==2 not implemented',__FILE__,__LINE__)
! 
! Calculate exchange-correlation kernel Kxc in G-space, put it into Kxc(1+igx+igy*...)
! TODO this has to be checked. xc_kernel is messy, we could call setvtr.F90 
  deallocate(kxc)
  allocate(kxc(Ep%npwe,Ep%npwe),STAT=istat) ; if (istat/=0) call memerr(FILE__,'kxc',Ep%npwe**2,'spc')
! call ckxcldag(Dtset%paral_kgb,ngfft1,ngfft1a,ngfft2,ngfft3,nfftgw_tot,rhor(:,1),kxc)
  call xc_kernel(Dtset,Dtset%ixc,MPI_enreg,ngfft_gw,nfftgw_tot,Ep%nsppol,rhor,rprimd,igfft0,&
&  Ep%npwe,gmet,kxc,gvec,Qmesh%nibz,Qmesh%ibz)

  do ig=1,Ep%npwe
!  chitmp(ig,1)=kxc(igfft(ig,Ep%mG0(1)+1,Ep%mG0(2)+1,Ep%mG0(3)+1))
   chitmp(ig,1)=kxc(ig,1)
  end do
  write(msg,'(a)')' kxc(g)'
  call wrtout(std_out,msg,'COLL')
  call print_arr(chitmp(:,1))
  deallocate(chitmp)
 end if !of if (Ep%tddft)
!
!=== Symmetrized dielectric matrix and array for chi0sumrule check ===
 Ep%nI=1 ; Ep%nJ=1 
 if (Dtset%nspinor==2) then 
  if (Dtset%usepaw==1.and.Dtset%pawspnorb>0) then
   Ep%nI=1 ; Ep%nJ=4 
  end if
! For spin-spin interaction
! Ep%nI=4 ; Ep%nJ=4 
  if (Ep%npwepG0/=Ep%npwe) STOP "If spinor npwepG0 must be == npwe"
  if (Ep%symchi/=0) STOP "symchi/0 and spinor not available"
 end if
 write(*,*)"using Ep%nI, Ep%nJ",Ep%nI,Ep%nJ

 allocate(chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega),stat=istat)
 allocate(chi0sumrule(Ep%npwe))
 if (istat/=0) then
  call memerr(FILE__,'chi0',Ep%npwe**2*Ep%nomega,'spc')
 end if
!
!=== Enable calculations of chi0 on one or more selected q-points ===
 if (Dtset%nqptdm>0) then
  write(msg,'(6a)')ch10,&
&  ' Dielectric matrix will be calculated on some ',ch10,&
&  ' selected q points provided by the user through the input variables ',ch10,&
&  ' nqptdm and qptdm'
  call wrtout(std_out,msg,'COLL')  
  call wrtout(ab_out,msg,'COLL')
  ltest=(Dtset%nqptdm<=Qmesh%nibz)
  write(msg,'(a)')'nqptdm should not exceed the number of q points in the IBZ'
  call assert(ltest,msg,__FILE__,__LINE__)
  allocate(qptdm(3,Dtset%nqptdm))
! Size of 2-dimension of Dtset%qptdm might be larger
  qptdm(:,:)=Dtset%qptdm(:,1:Dtset%nqptdm) 
! * Check if the provided q-points are correct
  do iq=1,Dtset%nqptdm
   found=.FALSE.
   do iqp=1,Qmesh%nibz
    qtmp(:)=qptdm(:,iq)-Qmesh%ibz(:,iqp) 
    found=(normv(qtmp,gmet,'G')<GW_TOLQ) ; if (found) EXIT
   end do
   write(msg,'(a)')'One or more qptdm points do not satisfy q=k1-k2'
   call assert(found,msg,__FILE__,__LINE__)
  end do 
 end if 

!if (Dtset%prtvol==10) then 
!=== Open file to write independent matrix elements of \epsilon^-1 ===
 if (rank==master) then 
  filnam=TRIM(Dtfil%filnam_ds(4))//'_EM1'
  call isfile(filnam,'new')
  open(unem1ggp,file=filnam,status='unknown',form='formatted')
 end if 
!end if 
 call pclock(9)
 call timab(305,2,tsec) ! screening(2)
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END OF INITIALIZATION PART <<<<<<<<<<<<<<<<<<<<<<<<<<<
!
 write(msg,'(84a)')ch10,('-',ii=1,80),ch10,&
& ' screening : beginning of the main loop',ch10
 call wrtout(std_out,msg,'COLL')
!call print_wavefunctions_information(Wf)
!call print_wavefunctions_information(Wf_val)
!
!=== Loop over q-points. Calculate \epsilon^{-1} and save on disc ===
 do iq=1,Qmesh%nibz
  call timab(306,1,tsec) !loop1
! 
! * Selective q-points calculation
  found=.FALSE. ; label=iq
  if (Dtset%nqptdm>0) then
   do iqp=1,Dtset%nqptdm
    qtmp(:)=Qmesh%ibz(:,iq)-Dtset%qptdm(:,iqp) ; found=(normv(qtmp,gmet,'G')<GW_TOLQ)
    if (found) then 
     label=iqp ; EXIT !iqp
    end if 
   end do
   if (.not.found) CYCLE !iq
  end if

  bar=REPEAT('-',80)
  write(msg,'(4a,1x,a,i2,a,f9.6,2(",",f9.6),3a)')ch10,ch10,bar,ch10,&
&  ' q-point number ',label,'        q = (',(Qmesh%ibz(ii,iq),ii=1,3),') [r.l.u.]',ch10,bar
  call wrtout(std_out,msg,'COLL') 
  call wrtout(ab_out,msg,'COLL')
  qeq0=(normv(Qmesh%ibz(:,iq),gmet,'G')<GW_TOLQ0)
!  === If q==0 allocate wings for each q-direction ===
  dim_wings=0 ; if (qeq0) dim_wings=Ep%nqlwl
  allocate(lwing(Ep%npwe,Ep%nomega,dim_wings))
  allocate(uwing(Ep%npwe,Ep%nomega,dim_wings))
  lwing=zero
  uwing=zero

! #if defined DEBUG_MODE
! if (Dtset%prtvol==10) then
! === Find the independent set of G-Gp pairs for this q-point. ===
! Useful to write the independent matrix elements of epsilon^-1 or to speed-up the KK transform 
! In the long wavelength limit we set q==0, because we still can use symmetryes for the Body.
  qtmp(:)=Qmesh%ibz(:,iq) ; if (normv(qtmp,gmet,'G')<GW_TOLQ0) qtmp(:)=zero
  call init_Gpairs_type(Gpairs_q,qtmp,Gsph_epsG0,Cryst)
! end if 
! #endif
  call pclock(10*iq)

  if (MPI_enreg%gwpara==1) then
!  === Parallelization over k-points, redefine the distribution of k-points ===
!  IMPORTANT: Note that Kmesh%nbz is used in the first dimension of proc_distrb. This implies that 
!  proc_distrb *MUST* be redefined if it is passed to routines employing the IBZ indexing.
   deallocate(MPI_enreg%proc_distrb) ; allocate(MPI_enreg%proc_distrb(Kmesh%nbz,Ep%nbnds,Ep%nsppol)) 
!  If nprocs>ntasks, proc_distrb==-999 for rank>ntasks-1 and no harm should be done
   MPI_enreg%proc_distrb(:,:,:)=-999
   allocate(istart(nprocs),istop(nprocs))

   if (Ep%symchi==0) then  
!   === No symmetries, divide the full BZ among procs ===
    ntasks=Kmesh%nbz 
    call split_work2(ntasks,nprocs,istart,istop)
    do irank=0,nprocs-1
     ii=istart(irank+1) ; jj=istop(irank+1)
     MPI_enreg%proc_distrb(ii:jj,:,:)=irank
    end do 
   else if (Ep%symchi/=0) then  
!   === Divide IBZ_q among procs, Distrb might be not so efficient for particular qs ===
!   Here proc_distrb is -999 for all the k-points not in the IBZ_q 
    ntasks=SUM(Ltg_q(iq)%ibzq(:)) 
    call split_work2(ntasks,nprocs,istart,istop)
    do irank=0,nprocs-1
     do ik=1,Kmesh%nbz
      do ikq=istart(irank+1),istop(irank+1)
       if (Ltg_q(iq)%bz2ibz(ik)==ikq) MPI_enreg%proc_distrb(ik,:,:)=irank
      end do 
     end do
    end do 
   end if 
   deallocate(istart,istop)
!  == Announce the treatment of k-points by each proc ==
   do ik=1,Kmesh%nbz 
    do is=1,Ep%nsppol
     if (MPI_enreg%proc_distrb(ik,Ep%nbnds,is)==rank) then 
      write(msg,'(3(a,i4))')'P : treating k-point ',ik,' and spin ',is,' by node ',rank
      call wrtout(std_out,msg,'PERS')
     end if 
    end do 
   end do
  end if !gwpara==1
! 
! === Special treatment of the long wavelenght limit ===
  if (qeq0) then 
   call timab(307,1,tsec)
   call status(iq,Dtfil%filstat,iexit,level,'call cchi0q0  ')

   if (read_occupied) then 
!   === Adler-Wiser with time-reversal ===
    call cchi0q0(Dtset,Cryst,Dtfil,Ep,Psps,Kmesh,Gsph_epsG0,Gsph_wfn,gvec,Pawang,Pawtab,ktabr,&
&    nbv,nbvw,occ,ngfft_gw,igfft,nfftgw,ks_energy,gw_energy,chi0,MPI_enreg,Ltg_q(iq),&
&    my_minb,my_maxb,pawrhox_spl,Cprj_bz,vkbsign,vkb,vkbd,chi0sumrule,&
&    Wf=Wf,Wf_val=Wf_val) ! Optional arguments
   else 
#ifndef FC_PGI6
!   === Adler-Wiser without time-reversal ===
    call cchi0q0(Dtset,Cryst,Dtfil,Ep,Psps,Kmesh,Gsph_epsG0,Gsph_wfn,gvec,Pawang,Pawtab,ktabr,&
&    nbv,nbvw,occ,ngfft_gw,igfft,nfftgw,ks_energy,gw_energy,chi0,MPI_enreg,Ltg_q(iq),&
&    my_minb,my_maxb,pawrhox_spl,Cprj_bz,vkbsign,vkb,vkbd,chi0sumrule,&
&    Wf=Wf) ! Optional arguments
#else 
    write(msg,'(a)')' calculation not implemented under PGI6, please use awtr=1, or change compiler!'
    call wrtout(std_out,msg,'COLL') 
    call leave_new('COLL')
#endif
   end if 
!  
!  === Free KB form factors === 
!  FIXME here we have a memory leak if nqptdm
   deallocate(vkb,vkbd,vkbsign)
   call timab(307,2,tsec)
  else 
!  === Calculate cchi0 for q/=0 ===
   call timab(308,1,tsec) ! screening(cchi0)
   call status(iq,Dtfil%filstat,iexit,level,'call cchi0    ')

   if (read_occupied) then
!   === Adler-Wiser with time reversal ===
    call cchi0(Dtset,Cryst,Dtfil,Qmesh%ibz(:,iq),Ep,Psps,Kmesh,Gsph_epsG0,Gsph_wfn,gvec,Pawang,Pawtab,nbv,nbvw,occ,&
&    ngfft_gw,igfft,nfftgw,gw_energy,chi0,MPI_enreg,ktabr,Ltg_q(iq),my_minb,my_maxb,pawrhox_spl,Cprj_bz,chi0sumrule,&
&    Wf=Wf,Wf_val=Wf_val) ! Optional arguments
   else 
#ifndef FC_PGI6
!   === Adler-Wiser without time-reversal ===
    call cchi0(Dtset,Cryst,Dtfil,Qmesh%ibz(:,iq),Ep,Psps,Kmesh,Gsph_epsG0,Gsph_wfn,gvec,Pawang,Pawtab,nbv,nbvw,occ,&
&    ngfft_gw,igfft,nfftgw,gw_energy,chi0,MPI_enreg,ktabr,Ltg_q(iq),my_minb,my_maxb,pawrhox_spl,Cprj_bz,chi0sumrule,&
&    Wf=Wf) ! Optional arguments
#else 
    write(msg,'(a)')' calculation not implemented under PGI6, please use awtr=1, or change compiler!'
    call wrtout(std_out,msg,'COLL') 
    call leave_new('COLL')
#endif
   end if
   call timab(308,2,tsec)
  end if ! cchi0 or cchi0q0
! 
! === Print chi0(q,G,Gp,omega), then calculate epsilon and epsilon^-1 for this q ===
! * Only master works but this part could be parallelized over frequencies
  call pclock(10*iq+1)

  do iomega=1,Ep%nomega
   write(msg,'(1x,a,i4,a,2f9.4,a)')&
&   ' chi0(G,G'') at the ',iomega,' th omega',Ep%omega(iomega)*Ha_eV,' [eV]'
   call wrtout(std_out,msg,'COLL') 
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(1x,a,i3,a,i4,a)')' chi0(q =',iq, ', omega =',iomega,', G,G'')'
   if (Dtset%nqptdm>0) write(msg,'(a,i3,a,i4,a)')'  chi0(q=',iqp,', omega=',iomega,', G,G'')'
   call wrtout(std_out,msg,'COLL') 
!  call wrtout(ab_out,msg,'COLL')
   call print_arr(chi0(:,:,iomega))
   call print_arr(chi0(:,:,iomega),max_r=2,unit=ab_out)
  end do

! Divide by the volume 
! TODO this should be done in cchi0, but I have to update all the test, sigh!
  chi0(:,:,:)=chi0(:,:,:)/ucvol
! 
! === Write chi0 on _SUSC file ===
  if (rank==master.and..TRUE.) then
   title(1)='CHI0 file: chi0'
   if (qeq0) then
    string='0' ; if (Ep%inclvkb/=0) call int2char(Ep%inclvkb,string)
    title(1)=title(1)(1:21)//', calculated using inclvkb = '//string
   end if 
   optfil=1 ! Write _SUSC file.
   unt_susc=Dtfil%unchi0 ; fname_susc=TRIM(Dtfil%filnam_ds(4))//'_SUS'
   call wrscr(iq,optfil,unt_susc,fname_susc,Hdr_kss,Dtset,Ep%npwe,Ep%npwwfn,Ep%nbnds,&
&   Qmesh%nibz,Ep%nqlwl,Ep%nomega,Qmesh%ibz,Ep%omega,gvec,gmet,chi0,title,&
&   Ep%qlwl,uwing,lwing)
  end if
! 
! Quick and dirty coding of the RPA functional
! development stage!
  if(.false.) call calc_rpa_functional(iq,Ep,Vcp,Qmesh,Dtfil,gmet,kxc,MPI_enreg,chi0)
! 
! =======================================================
! === Calculate \tilde\epsilon^{-1} overwriting chi0 ====
! =======================================================
  ! TODO change all these logical flags to integer
  if (Ep%rpa  ) approx_type=0 
  if (Ep%tddft) approx_type=1
  if (Ep%testparticle) option_test=0
  if (Ep%testelectron) option_test=1
  call make_epsm1_driver(iq,Ep%npwe,Ep%nI,Ep%nJ,Ep%nomega,Ep%nomegaer,Ep%omega,&
&  approx_type,option_test,Qmesh%nibz,Qmesh%ibz,Vcp,Dtfil,gmet,kxc,MPI_enreg,chi0)

  epsm1   => chi0
  vc_sqrt => Vcp%vc_sqrt(:,iq)  ! Contains vc^{1/2}(q,G), complex-valued due to a possible cutoff
! 
! === Output the sum rule evaluation ===
  chi0sumrule(:)=chi0sumrule(:)/ucvol
  call output_chi0sumrule(qeq0,iq,Ep%npwe,omegaplasma,chi0sumrule,epsm1(:,:,1),vc_sqrt)

! === Write heads and wings on the main output ===
  if (qeq0) then 
   write(msg,'(1x,2a)')' Heads and wings of the symmetrical epsilon^-1(G,G'') ',ch10
   call wrtout(ab_out,msg,'COLL')
   do iomega=1,Ep%nomega
    write(msg,'(2x,a,i4,a,2f9.4,a)')&
&    ' Upper and lower wings at the ',iomega,' th omega',Ep%omega(iomega)*Ha_eV,' [eV]'
    call wrtout(ab_out,msg,'COLL')
    call print_arr(epsm1(1,:,iomega),max_r=9,unit=ab_out)
    call print_arr(epsm1(:,1,iomega),max_r=9,unit=ab_out)
    write(msg,'(a)')ch10 
    call wrtout(ab_out,msg,'COLL')
   end do
  end if

  if (rank==master) then 
   call timab(310,1,tsec) ! wrscr
!  if (Dtset%prtvol==10) then 
!  #if defined DEBUG_MODE
!  Independent matrix elements of \tilde epsilon^-1 (the later only if prtvol=10)
   write(msg,'(a,3(f10.6),a)')' Symmetrical epsilon^-1(G,G'') at q = ( ',Qmesh%ibz(:,iq),' ) [r.l.u.]'
   call outeps(Ep%npwe,Ep%nomega,Ep%omega,epsm1,Gsph_epsG0,Gpairs_q,msg,unem1ggp,Dtset%prtvol) 
!  #endif 
!  end if 
!  
!  === Write symmetrical epsilon^-1 on file ===
   title(1)='SCR file: epsilon^-1'
   if (qeq0) then
    string='0' ; if (Ep%inclvkb/=0) call int2char(Ep%inclvkb,string)
    title(1)=title(1)(1:21)//', calculated using inclvkb = '//string
   end if 
   if (Ep%testparticle) title(2)='TESTPARTICLE' 
   if (Ep%testelectron) title(2)='TESTELECTRON'
   ctype='RPA' ; if (Ep%tddft) ctype='TDDFT'
   title(2)(14:17)=ctype !this has to be modified

   optfil=0 ! Write \tilde\epsilon^-1
   unt_scr=Dtfil%unscr ; fname_scr=TRIM(Dtfil%filnam_ds(4))//'_SCR'
   allocate(lwing_dum(Ep%npwe,Ep%nomega,optfil),uwing_dum(Ep%npwe,Ep%nomega,optfil),qsmall_dum(3,optfil))

   call wrscr(iq,optfil,unt_scr,fname_scr,Hdr_kss,Dtset,Ep%npwe,Ep%npwwfn,Ep%nbnds,&
&   Qmesh%nibz,0,Ep%nomega,Qmesh%ibz,Ep%omega,gvec,gmet,epsm1,title,&
&   qsmall_dum,lwing_dum,uwing_dum)
   deallocate(lwing_dum,uwing_dum,qsmall_dum,STAT=istat)

  end if ! Only master TODO this part can be parallelized but I have to use xsum_mpi in cchi0 and cchi0q0   

  deallocate(lwing,uwing,STAT=istat)

  call pclock(10*iq+9)
  call timab(310,2,tsec)
  call timab(306,2,tsec)
 end do ! Loop over q-points
!
!----------------------------- END OF THE CALCULATION ------------------------
!
!=== Close Files ===
!if (Dtset%prtvol==10 .and. rank==0) close(unem1ggp)
 if (rank==master) close(unem1ggp)
!
!=== Deallocate memory ===
 call status(0,Dtfil%filstat,iexit,level,'deallocate    ')
 deallocate(chi0,kxc)
 deallocate(ks_energy,gw_energy)
 deallocate(nband_,occ,rhor,rhog,nbv,ktabr)
 deallocate(gvec,igfft,irottb,igfftf,irottbf)
 deallocate(Pawfgr%fintocoa,Pawfgr%coatofin)
 if (Dtset%nqptdm>0) deallocate(qptdm)

 !if (allocated(kxc   )) deallocate(kxc)
 if (allocated(dimlmn)) deallocate(dimlmn )
 if (associated(MPI_enreg%proc_distrb)) deallocate(MPI_enreg%proc_distrb)
!
!=== Destroy the dinamic arrays in the local data structures ===
!
!* Optional deallocation for PAW
 if (Dtset%usepaw==1) then 
  deallocate(nhat,pawrhox_spl) ; if (nhatgrdim>0) deallocate(nhatgr)
  call cprj_free(Cprj_ibz)          ; deallocate(Cprj_ibz)
  call cprj_free(Cprj_bz )          ; deallocate(Cprj_bz )
  call rhoij_free(Pawrhoij)         ; deallocate(Pawrhoij)
  call destroy_pawfgrtab(Pawfgrtab) ; deallocate(Pawfgrtab)
  call destroy_paw_ij(Paw_ij)       ; deallocate(Paw_ij)
  call destroy_paw_an(Paw_an)       ; deallocate(Paw_an)
 end if

 call destroy_wf_info(Wf)
 if (read_occupied) call destroy_wf_info(Wf_val) 
 do iq=1,Qmesh%nibz
  call destroy_little_group(Ltg_q(iq))
 end do
 call destroy_BZ_mesh_type(Kmesh)
 call destroy_BZ_mesh_type(Qmesh)
 call destroy_Crystal_structure(Cryst)
 call destroy_Gvectors(Gsph_epsG0)
 call destroy_Gvectors(Gsph_wfn)
 call destroy_Gpairs_type(Gpairs_q)
 call destroy_Coulombian(Vcp) 
 call destroy_Epsilonm1_parameters(Ep)
 call hdr_clean(Hdr_kss)

 call pclock(9999)
 call timab(301,2,tsec)
 call status(0,Dtfil%filstat,iexit,level,'exit          ')

#if defined DEBUG_MODE
 write(msg,'(a)')' screening ended ' 
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

!ABI_WARNING("HELLO")
!ABI_COMMENT("HELLO")
!ABI_ERROR("HELLO")
!ABI_BUG("HELLO")
!MEM_CHECK(ii)
!IO_CHECK(msg,ii,ii)

end subroutine screening
!!***
