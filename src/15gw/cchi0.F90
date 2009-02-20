!{\src2tex{textfont=tt}}
!!****f* ABINIT/cchi0
!! NAME
!! cchi0
!!
!! FUNCTION
!! Main calculation of the independent-particle susceptibility chi0 for qpoint!=0
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! Dtset <type(dataset_type)>=all input variables in this dataset
!! Dtfil <type(datafiles_type)>=variables related to files
!! k_mesh <type(bz_mesh_type)>= datatype gathering parameters related to the k-point sampling 
!!    %nibz=number of k-points in the IBZ
!!    %nbz=number of k-points in the BZ
!!    %bz(3,nbz)=reduced coordinates for k-points in the full Brillouin zone
!!    %ibz(3,nibz)=reduced coordinates for k-points in the irreducible wedge 
!!    %tab(nbz)=mapping between a kpt in the BZ (array bz) and the irred point in the array ibz
!!    %tabi(nbz)= -1 if inversion is needed to obtain this particular kpt in the BZ, 1 means identity
!!    %tabo(nbz)= for each point in the BZ, the index of the symmetry operation S in reciprocal 
!!      space which rotates k_IBZ onto \pm k_BZ (depending on tabi)
!!    %tabp(nbz)= For each k_BZ, it gives the phase factors associated to non-symmorphic operations, i.e 
!!      e^{-i 2 \pi k_IBZ \cdot R{^-1}t} == e{-i 2\pi k_BZ cdot t} where : 
!!      \transpose R{-1}=S and (S k_IBZ) = \pm k_BZ (depending on ktabi) 
!!    %tabr(nfftot,nbz) For each point r on the real mesh and for each k-point in the BZ, tabr 
!!      gives the index of (R^-1 (r-t)) in the FFT array where R=\transpose S^{-1} and k_BZ=S k_IBZ. 
!!      t is the fractional translation associated to R
!! Ep<type(epsilonm1_parameters_type)>= Parameters related to the calculation of the inverse dielectric matrix.
!!    %nbnds=number of bands summed over
!!    %npwe=number of planewaves for the irreducible polarizability X^0_GGp
!!    %npwvec=maximum number of G vectors (between Ep%npwe and Ep%npwwfn) 
!!     used to define the dimension of some arrays e.g igfft
!!    %npwwfn=number of planewaves for wavefunctions (input variable, might be modified in screening)
!!    %nsppol=1 for unpolarized, 2 for spin-polarized
!!    %nomega=total number of frequencies in X^0 (both real and imaginary)
!!    %nomegasf=number of real frequencies used to sample the imaginary part of X^0 (spectral method)
!!    %spmeth=1 if we use the spectral method, 0 for standard Adler-Wiser expression 
!!    %spsmear=gaussian broadening used to approximate the delta distribution
!!    %zcut=small imaginary shift to avoid poles in X^0
!! gwenergy(Kmesh%nibz,Ep%nbnds,Ep%nsppol)=GW energies, for self-consistency purposes
!! Gsph_wfn<gvectors_type data type> The G-sphere used to describe the wavefunctions.
!! Gsph_epsG0<gvectors_type data type> The G-sphere used to describe chi0/eps. (including umklapp G0 vectors)
!!    %ng=number of G vectors for chi0
!!    %rottbm1(ng,2,nsym)=contains the index (IS^{-1}) G 
!!    %phmGt(Ep%npwe,nsym)=phase factors e^{-iG \cdot t} needed to symmetrize oscillator matrix elements and epsilon 
!!    %gprimd(3,3)=dimensional reciprocal space primitive translations (b^-1)
!! gvec(3,Ep%npwvec)= G vectors in reduced coordiinates
!! gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!! igfft(Ep%npwvec,2*Ep%mG0(1)+1,2*Ep%mG0(2)+1,Ep%mG0(3)+1)=index of G-G0 in the FFT grid for 
!!  each G0 vector (see the cigfft.F90 routine)
!! MPI_enreg=datatype gathering information on parallelism.
!!    %gwpara= 0 no parallelism, 1 parallelism is over k-points, 2 bands are spread btw processors
!! nbvw=number of bands in the arrays wfrv
!! nbv=number of valence bands
!! ngfft(18)= array containing all the information for FFT (see input variable)
!! Cryst<Crystal_structure>= data type gathering info on symmetries and unit cell
!!    %natom=number of atoms
!!    %nsym=number of symmetries 
!!    %xred(3,natom)=reduced coordinated of atoms
!!    %typat(natom)=type of each atom
!!    %rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!    %timrev= 2 if time reversal can be used, 1 otherwise
!! nfftot=number of points in FFT grid
!! occ(Kmesh%nibz,Ep%nbnds,Ep%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! qpoint(3)=reciprocal space coordinates of the q wavevector
!! Ltg_q= little group datatype FIXME to be better describe (also cleaned)
!! Pawang<pawang_type> angular mesh discretization and related data:
!! Pawtab(Psps%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! Wf%wfr(nfftot,my_minb:my_maxb,Kmesh%nibz,Ep%nsppol) = (optional) wavefunctions in real space 
!! wfrv(Ep%npwwfn,nbvw,Kmesh%nibz,Ep%nsppol)= (optional) array containing fully and partially occupied states in t space
!! Cprj_bz(natom,Dtset%nspinor*Ep%nbnds*Kmesh%nbz*Ep%nsppol) <type(cprj_type)>=
!!  projected input wave functions <Proj_i|Cnk> with all NL projectors for each k-point in the full Brillouin zone
!!
!! OUTPUT
!!  chi0(Ep%npwe,Ep%npwe,Ep%nomega)=independent-particle susceptibility matrix at wavevector qpoint and
!!   each frequeny defined by Ep%omega and Ep%nomega
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!      assemblychi0,leave_new,rho_tw_g,timab,wrtout,xcomm_init
!!      xmaster_init,xme_init,xsum_master
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cchi0(Dtset,Cryst,Dtfil,qpoint,Ep,Psps,Kmesh,Gsph_epsG0,Gsph_wfn,gvec,Pawang,Pawtab,nbv,nbvw,occ,&
& ngfft,igfft,nfftot,gwenergy,chi0,MPI_enreg,ktabr,Ltg_q,my_minb,my_maxb,rhox_spl,Cprj_bz,chi0sumrule,&
& Wf,Wf_val)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : GW_TOL_DOCC, GW_TOL_W0, czero_gw
 use m_errors, only : assert
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13nonlocal
 use interfaces_13recipspace
 use interfaces_15gw, except_this_one => cchi0
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_maxb,my_minb,nbvw,nfftot
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Crystal_structure),intent(in) :: Cryst
 type(Datafiles_type),intent(in) :: Dtfil
 type(Dataset_type),intent(in) :: Dtset
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(Gvectors_type),intent(in) :: Gsph_epsG0,Gsph_wfn
 type(Little_group),intent(in) :: Ltg_q
 type(MPI_type),intent(inout) :: MPI_enreg
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
 type(Wavefunctions_information),intent(inout),optional :: Wf,Wf_val
!arrays
 integer,intent(in) :: gvec(3,Ep%npwvec),ktabr(nfftot,Kmesh%nbz),nbv(Ep%nsppol)
 integer,intent(in) :: ngfft(18)
 integer,intent(in),target :: igfft(Ep%npwvec,2*Ep%mG0(1)+1,2*Ep%mG0(2)+1,2*Ep%mG0(3)+1)
 real(dp),intent(in) :: gwenergy(Kmesh%nibz,Ep%nbnds,Ep%nsppol)
 real(dp),intent(in) :: occ(Kmesh%nibz,Ep%nbnds,Ep%nsppol),qpoint(3)
 real(dp),intent(in) :: rhox_spl(Psps%mqgrid_ff,2,0:2*(Psps%mpsang-1),Psps%lnmax*(Psps%lnmax+1)/2,Psps%ntypat*Dtset%usepaw)
 real(dp),intent(out) :: chi0sumrule(Ep%npwe)
 complex(gwpc),intent(out) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
 type(Cprj_type),intent(in) :: Cprj_bz(Cryst%natom,Dtset%nspinor*Ep%nbnds*Kmesh%nbz*Ep%nsppol*Dtset%usepaw)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)

!Local variables ------------------------------
 character(len=50),parameter :: FILE__='screening.F90'
 !type(transitions_type) :: trans_q
!scalars
 integer,parameter :: level=24,tim_fourdp=1
 integer :: dim1_rhox,dim2_rhox,dim_rtwg,i1,i2,iat,ib,ib1,ib2,ibsp,ierr,ig,ig01
 integer :: ig02,ig03,ig1,ig2,igp,ii,ik_bz,ik_ibz,ikmq_bz,ikmq_ibz,indx_k_bz
 integer :: indx_kmq_bz,io,iomegal,iomegar,iosf,ir,is,ispinor,istat,isym,isym_k
 integer :: isym_kmq,itim,itim_k,itim_kmq,k0g,kg,master,mqmem_,my_wl,my_wr
 integer :: natom,nfound,nkpt_summed,nprocs,nqpt_,nspinor,optder,rank
 integer :: rhoxsp_method,shift,spaceComm,spad,two_lmaxp1
 real(dp) :: deltaeGW_b1kmq_b2k,deltaeGW_enhigh_b2k,deltaf_b1kmq_b2k,domega
 real(dp) :: e_b1_kmq,enhigh,f_b1_kmq,factor,max_rest,min_rest,my_max_rest
 real(dp) :: my_min_rest,numerator,omegal,omegar,sgn,spin_fact,weight,wl,wr
 complex(dpc) :: ct
 complex(gwpc) :: ph_mkmqt,ph_mkt
 logical :: ltest,qzero,time_reversal
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: G0(3),wtk_ltg(Kmesh%nbz)
 integer,allocatable :: dimlmn(:),npwarr(:),tabr_k(:),tabr_kmq(:)
 integer,pointer :: igfft0(:),igfftG0(:)
 real(dp) :: gmet(3,3),gprimd(3,3),kbz(3),kmq_bz(3),spinrot_k(4),spinrot_kmq(4)
 real(dp) :: tsec(2)
 real(dp),allocatable :: omegasf(:),pawrhox(:,:,:,:),qptns(:,:),ylm_q(:,:)
 real(dp),allocatable :: ylmgr_q(:,:,:)
 complex(dpc),allocatable :: green_enhigh_w(:),green_w(:),kkweight(:,:)
 complex(gwpc),allocatable :: chi0sf(:,:,:),rhotwg(:),wfr1(:)
 complex(gwpc),allocatable :: wfr2(:),wfwfg(:)
 type(Cprj_type),allocatable :: Cprj_k(:,:),Cprj_kmq(:,:)

!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(2a)')ch10,' cchi0 : enter '
 call wrtout(std_out,msg,'COLL')
 call flush_unit(std_out)
#endif

 call timab(331,1,tsec) ! cchi0
 !
 ! === Check optional arguments === 
 if (PRESENT(Wf_val)) then 
  ltest=(PRESENT(Wf_val).and.PRESENT(Wf))
  call assert(ltest,'Both Wf_val% and Wf% must be present',__FILE__,__LINE__)
 end if
 !
 ! === Initialize MPI related quantities ===
 call xcomm_init  (MPI_enreg,spaceComm) 
 call xmaster_init(MPI_enreg,master   )  
 call xme_init    (MPI_enreg,rank     )          
 call xproc_max(nprocs,ierr)
 !
 ! == Copy some values ===
 nspinor = Dtset%nspinor
 natom   = Cryst%natom
 gprimd(:,:)=Gsph_epsG0%gprimd(:,:) 
 gmet(:,:)  =Gsph_epsG0%gmet(:,:)

 dim_rtwg=1
 if (nspinor==2) then
  !can reduce size depending on Ep%nI and Ep%nj
  dim_rtwg=4
 end if
 !
 ! === Initialize the completeness correction  === 
 if (Dtset%gwcomp==1) then
  enhigh=MAXVAL(gwenergy(:,Ep%nbnds,:)) + Dtset%gwencomp
  allocate(wfwfg(nfftot*nspinor**2))
  allocate(green_enhigh_w(Ep%nomega))
  write(msg,'(a,f8.2,a)')' use the completeness correction with the energy ',enhigh*Ha_eV,' [eV]'
  call wrtout(std_out,msg,'COLL')
 end if
 !
 ! === Setup weights (2 for spin unpolarized sistem, 1 for polarized) ===
 ! * spin_fact is used to normalize the occupation factors to one.
 ! * Consider also the AFM case.
 SELECT CASE (Ep%nsppol)
 CASE (1)
  weight=two/Kmesh%nbz ; spin_fact=half     
  if (Dtset%nspden==2) then
   weight=one/Kmesh%nbz ; spin_fact=half     
  end if
  if (Dtset%nspinor==2) then
   weight=one/Kmesh%nbz ; spin_fact=one     
  end if
 CASE (2) 
  weight=one/Kmesh%nbz ; spin_fact=one
 CASE DEFAULT 
  call assert(.FALSE.,'Wrong value of nsppol',__FILE__,__LINE__)
 END SELECT

 ! === Weight for points in the IBZ_q ===
 wtk_ltg(:)=1
 if (Ep%symchi==1) then
  do ik_bz=1,Ltg_q%nbz
   wtk_ltg(ik_bz)=0
   if (Ltg_q%ibzq(ik_bz)/=1) CYCLE ! Only k in IBZ_q
   wtk_ltg(ik_bz)=SUM(Ltg_q%wtksym(:,:,ik_bz))
  end do
 end if

 time_reversal=.FALSE.
 if (PRESENT(Wf_val)) then 
  time_reversal=.TRUE.
  ! Use faster algorithm based on time reversal (presently only in parallel)
  ! Note that care has to be used in case of metals and/or spin since wfrv can contain unoccupied states 
  write(msg,'(4a)')ch10,&
&  ' Valence states are in memory',ch10,&
&  ' Using faster implementation based on time-reversal symmetry' 
  call wrtout(std_out,msg,'COLL')
 end if 

 if (rank==master) then 
  write(std_out,*)' use symmetries    = ',Ep%symchi
  write(std_out,*)' use time reversal = ',time_reversal
  write(std_out,*)' spectral method   = ',Ep%spmeth
 end if 

 if (Dtset%usepaw==1) then 
  ! === Set up REAL Ylm(q+G) up to 2*l_max for this q-point ===
  mqmem_=1 ; nqpt_=1 ; optder=0
  allocate(npwarr(nqpt_),qptns(3,nqpt_)) 
  npwarr(:)=Ep%npwepG0 ; qptns(:,1)=qpoint(:)
  two_lmaxp1=2*Psps%mpsang-1
  allocate(ylm_q(Ep%npwepG0*mqmem_,two_lmaxp1**2))
  allocate(ylmgr_q(Ep%npwepG0*mqmem_,3+6*(optder/2),two_lmaxp1**2))

  call status(0,Dtfil%filstat,0,level,'call initylmg ')
  call initmpi_seq(MPI_enreg_seq) 
  !
  ! * Dtset%nband and Dtset%nsppol are not used in sequential mode.
  call initylmg(gprimd,gvec,qptns,mqmem_,MPI_enreg_seq,two_lmaxp1,Ep%npwepG0,Dtset%nband,nqpt_,&
   npwarr,Dtset%nsppol,optder,Cryst%rprimd,Dtfil%unkg,Dtfil%unylm,ylm_q,ylmgr_q)
  !
  ! * Evaluate oscillator matrix elements btw partial waves.
  allocate(pawrhox(2,Ep%npwepG0,Psps%lmnmax*(Psps%lmnmax+1)/2,natom))  
  rhoxsp_method=2
  if (rhoxsp_method==1) then !Alouani-Method
   dim1_rhox=2*(Psps%mpsang-1) 
   dim2_rhox=Psps%lnmax*(Psps%lnmax+1)/2  
  else if (rhoxsp_method==2) then  ! Shiskin-Kresse
   dim1_rhox=MAXVAL(Pawtab(:)%l_size)**2
   dim2_rhox=MAXVAL(Pawtab(:)%lmn2_size) 
  end if
  call paw_mkrhox(Cryst,rhox_spl,gmet,gvec,rhoxsp_method,dim1_rhox,dim2_rhox,&
&  Psps,Pawang,Pawtab,qpoint,Ep%npwepG0,ylm_q,pawrhox)

  allocate(dimlmn(natom)) 
  do iat=1,natom
   dimlmn(iat)=Pawtab(Cryst%typat(iat))%lmn_size
  end do
  allocate(Cprj_k  (natom,nspinor*Ep%nbnds)) ; call cprj_alloc(Cprj_k,  0,dimlmn)
  allocate(Cprj_kmq(natom,nspinor*Ep%nbnds)) ; call cprj_alloc(Cprj_kmq,0,dimlmn)
 end if

 allocate(rhotwg(Ep%npwepG0*nspinor**2))
 allocate(tabr_k(nfftot),tabr_kmq(nfftot))
 allocate(wfr1(Wf%nfftot*nspinor),wfr2(Wf%nfftot*nspinor))

 SELECT CASE (Ep%spmeth) 
 CASE (0)
  write(msg,'(a)')' Calculating chi0(q,omega,G,G")   ' ; call wrtout(std_out,msg,'COLL')
  allocate(green_w(Ep%nomega))
 CASE (1,2)
  write(msg,'(a)')' Calculating Im chi0(q,omega,G,G")' ; call wrtout(std_out,msg,'COLL')
  !
  ! === Find Max and min resonant transitions for this q, report also values for this proc ===
  call make_transitions(1,Ep%nbnds,nbvw,Dtset%nsppol,Ep%symchi,Cryst%timrev,GW_TOL_DOCC,Ep%zcut,&
&  max_rest,min_rest,my_max_rest,my_min_rest,Kmesh,Ltg_q,MPI_enreg,Ep%mG0,gwenergy,occ,qpoint)

  if (MPI_enreg%gwpara==1) STOP "make_transitions is buggy" !there is a problem in make_transitions due to MPI_enreg
  !
  ! === Calculate frequency dependent weights for Hilbert transform ===
  allocate(omegasf(Ep%nomegasf),kkweight(Ep%nomega,Ep%nomegasf))
  !my_wl=1 ; my_wr=Ep%nomegasf
  call setup_spectral(Ep%nomega,Ep%omega,Ep%nomegasf,omegasf,max_rest,min_rest,my_max_rest,my_min_rest,&
&  0,Ep%zcut,zero,my_wl,my_wr,kkweight)

  if (MPI_enreg%gwpara==2.and..not.(PRESENT(Wf_val))) then 
   write(msg,'(a)')' cchi0 : BUG : valence wfs in r space are not in memory'
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end if 
  write(msg,'(2a,2i5,a,i4,2a)')ch10,&
&  ' cchi0 : allocating chi0sf using my_wl and my_wr = ',my_wl,my_wr,' ( ',my_wr-my_wl+1,' )',ch10
  call wrtout(std_out,msg,'PERS')
  allocate(chi0sf(Ep%npwe,Ep%npwe,my_wl:my_wr),STAT=istat)
  if (istat/=0) call memerr(FILE__,'chi0sf',Ep%npwe**2*(my_wr-my_wl+1)*Ep%nomegasf,'spc')
  chi0sf(:,:,:)=czero_gw
 CASE DEFAULT
  call assert(.FALSE.,'Wrong value for spmeth', __FILE__,__LINE__)
 END SELECT

 nkpt_summed=Kmesh%nbz 
 if (Ep%symchi==1) then
  nkpt_summed=Ltg_q%nibz_ltg 
  call print_little_group(Ltg_q,std_out,Dtset%prtvol,'COLL')
 end if
 write(msg,'(a,i6,a)')' Calculation status : ',nkpt_summed,' to be completed '
 call wrtout(std_out,msg,'COLL')
 !
 ! === Loop over k-points in BZ ===
 chi0sumrule(:)=zero
 chi0(:,:,:)=czero_gw 

 do ik_bz=1,Kmesh%nbz

  if (Ep%symchi==1) then  
   if (Ltg_q%ibzq(ik_bz)/=1) CYCLE  ! Only IBZ_q 
  end if 

  call status(ik_bz,Dtfil%filstat,0,level,'loop ikpt     ')
  !
  ! * Get ik_ibz, non-symmorphic phase, ph_mkt, and symmetries from ik_bz.
  call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k,ph_mkt)

  ! * Get index of k-q in the BZ, stop if not found as the weight=one/nkbz is not correct.
  call get_BZ_diff(Kmesh,kbz,qpoint,ikmq_bz,G0,nfound) ; if (nfound==0) call leave_new('COLL') 

  ! * Get ikmq_ibz, non-symmorphic phase, ph_mkmqt, and symmetries from ikmq_bz.
  call get_BZ_item(Kmesh,ikmq_bz,kmq_bz,ikmq_ibz,isym_kmq,itim_kmq,ph_mkmqt)

  ! * Copy tables for rotated FFT points
  tabr_k(:)  =ktabr(:,  ik_bz)
  spinrot_k(:)=Cryst%spinrot(:,isym_k) 

  tabr_kmq(:)=ktabr(:,ikmq_bz)
  spinrot_kmq(:)=Cryst%spinrot(:,isym_kmq) 

  ig01=G0(1)+Ep%mG0(1)+1
  ig02=G0(2)+Ep%mG0(2)+1
  ig03=G0(3)+Ep%mG0(3)+1
  igfftG0 => igfft(:,ig01,ig02,ig03)
  igfft0  => igfft(:,Ep%mG0(1)+1,Ep%mG0(2)+1,Ep%mG0(3)+1)
  ! 
  ! === Loop on spin to calculate trace $\chi_{up,up}+\chi_{down,down}$ ===
  ! * Only $\chi_{up,up} for AFM.
  do is=1,Ep%nsppol 
   !  
   ! Parallelization over k-points in the full BZ
   ! Note that spin para is not active but the check should be done here
   ! array proc_distrb does not depend on band index (there is a check in screening)
   if (MPI_enreg%gwpara==1) then 
    if ((MPI_enreg%proc_distrb(ik_bz,1,is))/=rank) cycle
    if (ANY(MPI_enreg%proc_distrb(ik_bz,:,is)==-999)) then 
     write(*,*)' cchi0 : wrong distrb' ; call leave_new('COLL')
    end if 
   end if
   write(msg,'(2(a,i4),a,i2,a,i3)')&
&   ' ik = ',ik_bz,' / ',Kmesh%nbz,' is = ',is,' done by processor ',rank
   call wrtout(std_out,msg,'PERS')

   if (Dtset%usepaw==1) then 
    ! === Load cprj for k and k-q ===
    ! * Do not take care of umklapp G0 in k-q as the phase is included ===
    shift=nspinor*Ep%nbnds*Kmesh%nbz*(is-1)
    indx_k_bz  =nspinor*Ep%nbnds*(  ik_bz-1)+shift
    indx_kmq_bz=nspinor*Ep%nbnds*(ikmq_bz-1)+shift
    ibsp=0
    do ib=1,Ep%nbnds
     do ispinor=1,nspinor
      ibsp=ibsp+1
      do iat=1,natom
       Cprj_k  (iat,ibsp)%cp(:,:)=Cprj_bz(iat,indx_k_bz  +ibsp)%cp(:,:)
       Cprj_kmq(iat,ibsp)%cp(:,:)=Cprj_bz(iat,indx_kmq_bz+ibsp)%cp(:,:)
      end do
     end do
    end do
   end if
   !  
   ! === Loop over "conduction" states ===
   do ib1=1,Ep%nbnds 
    if (MPI_enreg%gwpara==2) then
     if (MPI_enreg%proc_distrb(ik_ibz,ib1,is)/=rank) CYCLE
    end if

    e_b1_kmq=gwenergy(ikmq_ibz,ib1,is)
    f_b1_kmq=     occ(ikmq_ibz,ib1,is)
    !   
    ! === Loop over "valence" states ===
    do ib2=1,Ep%nbnds 
     !MG FIXME this has been commented by Fabien 
     if (MPI_enreg%gwpara==2.and.Dtset%gwcomp==0) then
      if (MPI_enreg%proc_distrb(ik_ibz,ib2,is)/=rank) CYCLE
     end if

     deltaf_b1kmq_b2k=spin_fact*(f_b1_kmq-occ(ik_ibz,ib2,is))

     if (Dtset%gwcomp==0) then
      ! === Skip negligible transitions ===
      if (ABS(deltaf_b1kmq_b2k) < GW_TOL_DOCC) CYCLE 
     else
      ! when the completeness correction is used,
      ! we need to also consider transitions with vanishing deltaf
      if (occ(ik_ibz,ib2,is) < GW_TOL_DOCC) CYCLE  
     end if
     deltaeGW_b1kmq_b2k=e_b1_kmq-gwenergy(ik_ibz,ib2,is)

     SELECT CASE (Ep%spmeth)

     CASE (0)
      ! === Standard Adler-Wiser expression ===
      ! * Add the small imaginary of the Time-Ordered RF only for non-zero real omega
      ! FIXME What about metals?
      SELECT CASE (time_reversal)
       CASE (.TRUE.)
        if (Dtset%gwcomp==0) then ! cannot be completely skipped in case of completeness correction
         if (ib1<ib2) CYCLE ! Here we GAIN a factor ~2 
        end if
        do io=1,Ep%nomega
         green_w(io) = mkG0w(Ep%omega(io),deltaf_b1kmq_b2k,deltaeGW_b1kmq_b2k,Ep%zcut,GW_TOL_W0)

         if (Dtset%gwcomp==1) then 
          ! * Calculate the completeness correction
          numerator= -spin_fact*occ(ik_ibz,ib2,is)
          deltaeGW_enhigh_b2k=enhigh-gwenergy(ik_ibz,ib2,is)
          ! Completeness correction is NOT valid for real frequencies
          if (REAL(Ep%omega(io))<GW_TOL_W0) then
           green_enhigh_w(io) = mkG0w(Ep%omega(io),numerator,deltaeGW_enhigh_b2k,Ep%zcut,GW_TOL_W0)
          else
           green_enhigh_w(io) = czero_gw
          end if
          if (deltaf_b1kmq_b2k<0.d0) then
           green_w(io)= green_w(io) - green_enhigh_w(io)
          else 
           ! Disregard green_w, since it is already accounted for through the time-reversal
           green_w(io)=             - green_enhigh_w(io)    
          end if
         end if !gwcomp==1
        end do !io

        ! === Add the "delta part", symmetrization is done inside the routine ===
        ! MG Be careful here as time-reversal itim is not the same as tabi 
        ! TODO Add PAW, doesnt work for spinor
        if (Dtset%gwcomp==1.and.ib1==ib2) then
         call get_wfr(Wf_val,MPI_enreg,ib2,ik_ibz,is,wfr2)
         call calc_wfwfg(MPI_enreg,Dtset%paral_kgb,tim_fourdp,tabr_k,Kmesh%tabi(ik_bz),&
&                        nfftot,ngfft,wfr2,wfr2,wfwfg)
         qzero=.FALSE.
         call completechi0_deltapart(ik_bz,qzero,Ep%symchi,Ep%npwe,Ep%npwvec,Ep%nomega,nspinor,&
&         nfftot,ngfft,gvec,igfft0,Gsph_wfn,Ltg_q,green_enhigh_w,wfwfg,chi0)
        end if

       CASE (.FALSE.)
        ! === Old slow implementation, sum over all possible trasitions ===
        sgn=deltaeGW_b1kmq_b2k/ABS(deltaeGW_b1kmq_b2k) 
        do io=1,Ep%nomega
         if (REAL(Ep%omega(io))>GW_TOL_W0) then
          green_w(io)= deltaf_b1kmq_b2k/(Ep%omega(io)+deltaeGW_b1kmq_b2k-(0.,1.)*sgn*Ep%zcut)
         else
          green_w(io)= deltaf_b1kmq_b2k/(Ep%omega(io)+deltaeGW_b1kmq_b2k)
         end if
        end do
       END SELECT

     CASE (1,2) 
      ! === Spectral method, WARNING here Im assuming time-reversal ===
      if (deltaeGW_b1kmq_b2k<0) CYCLE
      call approxdelta(Ep%nomegasf,omegasf,deltaeGW_b1kmq_b2k,Ep%spsmear,iomegal,iomegar,wl,wr,Ep%spmeth) 
     END SELECT
     !    
     ! === Form rho-twiddle(r)=u^*_{b1,kmq_bz}(r) u_{b2,kbz}(r) and its FFT transform ===
     SELECT CASE (time_reversal)
     CASE (.TRUE.)
      if (ib1>nbvw) then 
       if (ib2<=nbvw) then 
        call get_wfr(Wf,    MPI_enreg,ib1,ikmq_ibz,is,wfr1)
        call get_wfr(Wf_val,MPI_enreg,ib2,  ik_ibz,is,wfr2)
       else ! unlikely
        call get_wfr(Wf,MPI_enreg,ib1,ikmq_ibz,is,wfr1)
        call get_wfr(Wf,MPI_enreg,ib2,  ik_ibz,is,wfr2)
       end if
      else
       if (ib2>nbvw) STOP "BUG"
       if (rank/=master.and.Dtset%gwcomp==0) CYCLE ! Only master for "valence-valence" transitions
       call get_wfr(Wf_val,MPI_enreg,ib1,ikmq_ibz,is,wfr1)
       call get_wfr(Wf_val,MPI_enreg,ib2,  ik_ibz,is,wfr2)
      end if
     CASE (.FALSE.)
      call get_wfr(Wf,MPI_enreg,ib1,ikmq_ibz,is,wfr1)
      call get_wfr(Wf,MPI_enreg,ib2,  ik_ibz,is,wfr2)
     END SELECT

     call rho_tw_g(Dtset%paral_kgb,nspinor,Ep%npwepG0,nfftot,ngfft,1,igfftG0,&
&     wfr1,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,&
&     wfr2,itim_k  ,tabr_k  ,ph_mkt  ,spinrot_k,&
&     dim_rtwg,rhotwg,tim_fourdp,MPI_enreg)

     ! === For PAW add on-site contribution, projectors are already in BZ ===
     if (Dtset%usepaw==1) then 
! FIXME Find a clever way to deal with spinors
      spad=(nspinor-1)
      i1=ib1 ; if (nspinor==2) i1=(2*ib1-1)
      i2=ib2 ; if (nspinor==2) i2=(2*ib2-1)
      call paw_rho_tw_g(Ep%npwepG0,dim_rtwg,nspinor,natom,Psps%lmnmax,dimlmn,&
&      Cprj_kmq(:,i1:i1+spad),Cprj_k(:,i2:i2+spad),pawrhox,rhotwg)
     end if

     SELECT CASE (Ep%spmeth)
     CASE (0)
      ! === Adler-Wiser ===
      call assemblychi0_sym(ik_bz,nspinor,Ep,Ltg_q,green_w,Ep%npwepG0,rhotwg,Gsph_epsG0,chi0)
     CASE (1,2)
      ! === Spectral method ===
      ! FIXME Does not work with spinor
      call assemblychi0sf(ik_bz,nspinor,Ep%symchi,Ltg_q,Ep%npwepG0,Ep%npwe,rhotwg,Gsph_epsG0,&
&      deltaf_b1kmq_b2k,my_wl,iomegal,wl,my_wr,iomegar,wr,Ep%nomegasf,chi0sf)
     CASE DEFAULT
      call assert(.FALSE.,'Wrong value for spmeth', __FILE__,__LINE__)
     END SELECT 
     !
     ! Accumulating the sum rule on chi0 Eq. (5.284) in G. D. Mahan Many-Particle Physics 3rd edition.
     ! FIXME Does not work with spinor
     !chi0sumrule(:)=chi0sumrule(:) + spin_fact*occ(ik_ibz,ib2,is)*deltaeGW_b1kmq_b2k*ABS(rhotwg(:))**2
     factor=spin_fact*occ(ik_ibz,ib2,is)
     call accumulate_chi0sumrule(ik_bz,Ep%symchi,Ep%npwe,factor,deltaeGW_b1kmq_b2k,&
&     Ltg_q,Gsph_epsG0,Ep%npwepG0,rhotwg,chi0sumrule)
     !
     ! Include also the completeness correction in the sum rule
     if (Dtset%gwcomp==1) then 
      !chi0sumrule(:)=chi0sumrule(:) - spin_fact*occ(ik_ibz,ib2,is)*deltaeGW_enhigh_b2k*ABS(rhotwg(:))**2
      factor=-spin_fact*occ(ik_ibz,ib2,is)
      call accumulate_chi0sumrule(ik_bz,Ep%symchi,Ep%npwe,factor,deltaeGW_enhigh_b2k,&
&      Ltg_q,Gsph_epsG0,Ep%npwepG0,rhotwg,chi0sumrule)
      !if (ib1==Ep%nbnds) chi0sumrule(:)=chi0sumrule(:) + spin_fact*occ(ik_ibz,ib2,is)*deltaeGW_enhigh_b2k
      if (ib1==Ep%nbnds) chi0sumrule(:)=chi0sumrule(:) + wtk_ltg(ik_bz)*spin_fact*occ(ik_ibz,ib2,is)*deltaeGW_enhigh_b2k
     end if

    end do !ib2
   end do !ib1
  end do !is
 end do !ik_bz
 call leave_test(MPI_enreg)
 !
 ! === After big loop over transitions, now MPI ===
 ! * Master took care of the contribution in case of metallic|spin polarized systems.
 SELECT CASE (Ep%spmeth) 
 CASE (0)
  ! === Adler-Wiser: Sum contributions from each proc ===
  ! * Looping on frequencies to avoid problems with the size of the MPI packet
  do io=1,Ep%nomega
   call xsum_mpi(chi0(:,:,io),spaceComm,ierr)
  end do
  call leave_test(MPI_enreg)
  chi0(:,:,:)=weight*chi0(:,:,:)
 CASE (1,2)
  ! === Spectral method: perform Hilbert transform ===
  write(msg,'(2a,i3,a)')ch10,&
&  ' Performing Hilbert transform using method ',Ep%spmeth,' It might take a while ...'
  call wrtout(std_out,msg,'COLL')
  !
  ! First coding, The loop over ig1, ig2 could be optimised taking into account symmetries 
  do ig1=1,Ep%npwe
   do ig2=1,Ep%npwe
    do io=1,Ep%nomega 
     ct=czero
     do iosf=my_wl,my_wr
      ct=ct+kkweight(io,iosf)*chi0sf(ig1,ig2,iosf)
     end do
     chi0(ig1,ig2,io)=ct
    end do 
   end do
  end do 
  ! === Sum contributions from each proc === 
  ! * Looping on frequencies to avoid problems with the size of the MPI packet
  do io=1,Ep%nomega
   call xsum_master(chi0(:,:,io),master,spaceComm,ierr)
  end do 
  call leave_test(MPI_enreg)
  chi0(:,:,:)=weight*chi0(:,:,:)
 CASE DEFAULT 
  call assert(.FALSE.,'Wrong value for spmeth', __FILE__,__LINE__)
 END SELECT
 !
 ! === MPI for the sum rule ===
 call xsum_mpi(chi0sumrule(:),spaceComm,ierr)
 chi0sumrule(:)=chi0sumrule(:)*pi*weight   ! the pi factor comes from Im[1/(x-ieta)] = pi delta(x)
 !
 ! **********************************************
 ! **** Now master has chi0(q,G,Gp,Ep%omega) ****
 ! **********************************************
 ! Impose Hermiticity (valid only for zero or purely imaginary frequencies)
 ! MG what about metals, where we have poles around zero?
 do io=1,Ep%nomega
  if (ABS(REAL(Ep%omega(io)))<0.00001) then
   do ig2=1,Ep%npwe
    do ig1=1,ig2-1
     chi0(ig2,ig1,io)=CONJG(chi0(ig1,ig2,io))
    end do
   end do
  end if
 end do 

 if (Cryst%use_antiferro) then 
  call symmetrize_afm_chi0(Cryst,Gsph_epsG0,Ep%npwe,Ep%nomega,chi0)
 end if
 !
 ! === Deallocate memory ===
 deallocate(rhotwg,tabr_k,tabr_kmq)
 deallocate(wfr1,wfr2)

 if (Dtset%usepaw==1) then 
  deallocate(npwarr,qptns,dimlmn)
  deallocate(ylm_q,ylmgr_q,pawrhox)
  call cprj_free(Cprj_k  ) ; deallocate(Cprj_k  )
  call cprj_free(Cprj_kmq) ; deallocate(Cprj_kmq)
 end if

 if (allocated(green_enhigh_w)) deallocate(green_enhigh_w)
 if (allocated(wfwfg   )) deallocate(wfwfg   )
 if (allocated(kkweight)) deallocate(kkweight)
 if (allocated(omegasf )) deallocate(omegasf )
 if (allocated(green_w )) deallocate(green_w )  
 if (allocated(chi0sf  )) deallocate(chi0sf  )
 !call destroy_transitions(trans_q)

 call status(0,Dtfil%filstat,0,level,'exit          ')
 call timab(331,2,tsec)

#if defined DEBUG_MODE
 write(msg,'(a)')' cchi0 : exit '
 call wrtout(std_out,msg,'PERS') 
 call flush_unit(std_out)
#endif

end subroutine cchi0
!!***

!!****if* ABINIT/mkG0w
!! NAME
!! mkG0w
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function mkG0w(omega,numerator,delta_ene,zcut,TOL_W0) result(green0w)

 use defs_basis
 use m_gwdefs, only : j_dpc

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: TOL_W0,delta_ene,numerator,zcut
 complex(dpc) :: green0w
 complex(gwpc),intent(in) :: omega

!Local variables ------------------------------
!scalars
 real(dp) :: sgn

!************************************************************************

 if(delta_ene**2>tol14)then

  sgn=delta_ene/ABS(delta_ene)

  if (REAL(omega)>TOL_W0) then
   green0w =  numerator / (omega + delta_ene - j_dpc*sgn*zcut)&
&            -numerator / (omega - delta_ene + j_dpc*sgn*zcut)
  else
   green0w =  numerator / (omega + delta_ene)&
&            -numerator / (omega - delta_ene)
  endif

 else
   green0w = zero
 end if

end function mkG0w
!!***
