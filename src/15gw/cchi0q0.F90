!{\src2tex{textfont=tt}}
!!****f* ABINIT/cchi0q0
!! NAME
!! cchi0q0
!!
!! FUNCTION
!! Calculate chi0 in the limit q-->0
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  Dtset <type(dataset_type)>=all input variables in this dataset
!!  Dtfil <type(datafiles_type)>=variables related to files
!!  Ep= datatype gathering differening parameters related to the calculation of the inverse dielectric matrix
!!  energy(Kmesh%nibz,Ep%nbnds,Ep%nsppol)=KS energies
!!  Gsph_wfn<gvectors_data_type>: Info on the G-sphere used for the wavefunctions.
!!  Gsph_epsG0<gvectors_data_type>: Info on the G-sphere used to describe chi0/espilon (including umklapp)
!!    %ng=number of G vectors
!!    %rottbm1(ng,2,nsym)=contains the index (IS^{-1}) G  in the array gvec
!!    %phmGt(ng,nsym)=phase factor e^{-iG.\tau} needed to symmetrize oscillator matrix elements and chi0
!!    %gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!!    %gprimd(3,3)=dimensional reciprocal space primitive translations (b^-1)
!!  gwenergy(Kmesh%nibz,Ep%nbnds,Ep%nsppol)=GW energies, for self-consistency purposes
!!  igfft(Ep%npwvec,2*Ep%mG0(1)+1,Ep%mG0(2)+1,Ep%mG0(3)+1)= index of G-G0 planewaves for 
!!   each G0 vectors (see cigfft.F90 routine)
!!  Ep%inclvkb=flag to include (or not) the grad of Vkb
!!  Ltg_q= little group datatype
!!  MPI_enreg= datatype gathering information on parallelism, variables used 
!!    %gwpara 
!!     0 no parallelism 
!!     1 parallelism is over k-points
!!     2 bands are spread btw processors
!!  nbvw=number of bands in the (optional) arrays wfrv,wfgv
!!  Kmesh<bz_mesh_type> The k-point mesh 
!!   %kbz(3,nbz)=k-point coordinates, full Brillouin zone
!!   %tab(nbz)= table giving for each k-point in the BZ (kBZ), the corresponding 
!!   irreducible point (kIBZ), where kBZ= (IS) kIBZ and I is either the inversion or the identity
!!   %tabi(nbzx)= for each point in the BZ defines whether inversion  has to be 
!!   considered in the relation kBZ=(IS) kIBZ (1 => only S; -1 => -S)  
!!   %tabo(nbzx)= the symmetry operation S that takes kIBZ to each kBZ
!!   %tabp(nbzx)= phase factor associated to tnons e^{-i 2 \pi k\cdot R{^-1}t}
!!  ktabr(nfftot,Kmesh%nbz) index of R^-(r-t) in the FFT array, where k_BZ = (IS) k_IBZ and S = \transpose R^{-1} 
!!  Ep%nbnds=number of bands
!!  nbv=number of valence bands
!!  ngfft(18)= array containing all the information for FFT (see input variable)
!!  Ep%nomega=number of frequencies
!!  Cryst<Crystal_structure>= data type gathering info on symmetries and unit cell
!!   %natom=number of atoms 
!!   %nsym=number of symmetry operations
!!   %symrec(3,3,nsym)=symmetry operations in reciprocal space
!!   %typat(natom)=type of each atom
!!   %xred(3,natom)=reduced coordinated of atoms
!!   %rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!   %timrev=2 if time-reversal symmetry can be used, 1 otherwise
!!  Ep%npwe=number of planewaves for sigma exchange (input variable)
!!  Ep%npwvec=dimension of igfft
!!  Ep%npwwfn=number of planewaves for wavefunctions (input variable)
!!  nfftot=number of points of FFT grid
!!  Ep%nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(Kmesh%nibz,Ep%nbnds,Ep%nsppol)=occupation numbers, for each k point in IBZ, and each band
!!  Ep%omega(Ep%nomega)=frequencies
!!  Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!     %mpsang=1+maximum angular momentum for nonlocal pseudopotential
!!  Pawang<pawang_type> angular mesh discretization and related data:
!!  Wf%wfg(Ep%npwwfn,my_minb:my_maxb,Kmesh%nibz,Ep%nsppol)= (optional) 
!!   Wfs in real space, for each band treated by this processor
!!  Wf%wfr(nfftot,my_minb:my_maxb,Kmesh%nibz,Ep%nsppol)= (optional) 
!!   wfs in G space for each band treated by this processor
!!  Wf_val%wfg(Ep%npwwfn,nbvw,Kmesh%nbz,Ep%nsppol)= (optional) 
!!   array containing fully and partially occupied states in G space
!!  Wf_val%wfr(nfftot,nbvw,Kmesh%nibz,Ep%nsppol) = (optional) array containing unoccupied states in real space
!!  Cprj_bz(natom,Dtset%nspinor*Ep%nbnds*Kmesh%nbz*Ep%nsppol*Psps%usepaw) <type(Cprj_type)>=
!!  projected input wave functions <Proj_i|Cnk> with all NL projectors for each k-point in the full Brillouin zone
!!
!! OUTPUT
!!  chi0(Ep%npwe,Ep%npwe,Ep%nomega)=independent-particle susceptibility matrix for wavevector qq,
!!   and frequencies defined by Ep%omega
!!
!! NOTES
!!  The terms "head", "wings" and "body" of chi(G,Gp) refer to
!!  G=Gp=0, either G or Gp=0, and neither=0 respectively
!!
!! TODO
!!  Check npwepG0 before Switching on umklapp 
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cchi0q0(Dtset,Cryst,Dtfil,Ep,Psps,Kmesh,Gsph_epsG0,Gsph_wfn,gvec,Pawang,Pawtab,ktabr,&
& nbv,nbvw,occ,ngfft,igfft,nfftot,energy,gwenergy,chi0,MPI_enreg,Ltg_q,&
& my_minb,my_maxb,pawrhox_spl,Cprj_bz,vkbsign,vkb,vkbd,chi0sumrule,&
& Wf,Wf_val) ! Optional

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : GW_TOL_DOCC, GW_TOL_W0, czero_gw
 use m_errors, only : assert
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13nonlocal
 use interfaces_13recipspace
 use interfaces_15gw, except_this_one => cchi0q0
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_maxb,my_minb,nbvw,nfftot
 type(Crystal_structure),intent(in) :: Cryst
 type(Dataset_type),intent(in) :: Dtset
 type(Datafiles_type),intent(in) :: Dtfil
 type(MPI_type),intent(inout) :: MPI_enreg
 type(Little_group),intent(in) :: Ltg_q
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Gvectors_type),intent(in) :: Gsph_epsG0,Gsph_wfn
 type(Pseudopotential_type),intent(in) :: Psps
 type(Pawang_type),intent(in) :: Pawang
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
 type(Cprj_type),intent(in) :: Cprj_bz(Cryst%natom,Dtset%nspinor*Ep%nbnds*Kmesh%nbz*Ep%nsppol*Dtset%usepaw)
 type(Wavefunctions_information),optional,intent(inout) :: Wf,Wf_val
!arrays
 integer,intent(in) :: gvec(3,Ep%npwvec) 
 integer,target,intent(in) :: igfft(Ep%npwvec,2*Ep%mG0(1)+1,2*Ep%mG0(2)+1,Ep%mG0(3)+1) 
 integer,intent(in) :: ktabr(nfftot,Kmesh%nbz) 
 integer,intent(in) :: ngfft(18)
 integer,intent(in) :: nbv(Ep%nsppol)
 real(dp),intent(in) ::   energy(Kmesh%nibz,Ep%nbnds,Ep%nsppol)
 real(dp),intent(in) :: gwenergy(Kmesh%nibz,Ep%nbnds,Ep%nsppol)
 real(dp),intent(in) ::      occ(Kmesh%nibz,Ep%nbnds,Ep%nsppol)
 real(dp),intent(in) :: pawrhox_spl(Psps%mqgrid_ff,2,0:2*(Psps%mpsang-1),Psps%lnmax*(Psps%lnmax+1)/2,Psps%ntypat*Dtset%usepaw)
 real(dp),intent(in) :: vkb (Ep%npwwfn,Dtset%ntypat,Psps%mpsang,Kmesh%nibz) 
 real(dp),intent(in) :: vkbd(Ep%npwwfn,Dtset%ntypat,Psps%mpsang,Kmesh%nibz),vkbsign(Psps%mpsang,Dtset%ntypat)
 real(dp),intent(out) :: chi0sumrule(Ep%npwe)
 complex(gwpc),intent(out) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)

!Local variables ------------------------------
!scalars
 integer,parameter :: level=24,tim_fourdp=1,enough=10
 integer :: counter,rhoxsp_method,dim1_rhox,dim2_rhox,nspinor,ispinor,ibsp,iab,nab
 type(MPI_type) :: MPI_enreg_seq
 integer :: gwpara,ib,ib1,ib2,ig1,ig2,ig,igp,itim,ik,ik_bz,ik_ibz,io,isym,ir,is,istat
 integer :: nkpt_summed,spad,dim_rtwg
 integer :: iat,ilm,ii,i1,i2,two_lmaxp1
 integer :: mqmem_,nqpt_,optder,shift,indx_k_bz
 integer :: master,nprocs,spaceComm,ierr,rank,my_wl,my_wr
 integer :: iomegal,iomegar,iosf,method,natom,spad1,spad2
 real(dp) :: domega,spin_fact,deltaf_b1b2,weight,factor
 real(dp) :: max_rest,min_rest,my_max_rest,my_min_rest
 real(dp) :: enhigh,deltaeGW_enhigh_b2
 real(dp) :: omegal,omegar,wl,wr,numerator,deltaeGW_b1b2,sgn
 complex(dpc) :: deltaeKS_b1b2
 complex(gwpc) :: ct
 logical :: time_reversal,qzero,ltest
 character(len=500) :: msg
 character(len=50),parameter :: FILE__='cchi0q0.F90'
!arrays
 integer :: wtk_ltg(Kmesh%nbz)
 integer,allocatable :: dimlmn(:),npwarr(:)
 integer :: spinorwf_pad(2,4) 
 integer,allocatable :: tabr_k(:)
 integer,pointer :: igfft0(:)
 real(dp) :: kbz(3),qcart2red(3,3),qpoint(3)
 real(dp) :: gmet(3,3),gprimd(3,3),igradcart_paw(2,3),igradred_paw(2,3)
 real(dp) :: spinrot_kbz(4)
 real(dp),allocatable :: pawrhox(:,:,:,:)
 real(dp),allocatable :: qptns(:,:),omegasf(:)
 real(dp),allocatable :: ylm_q(:,:),ylmgr_q(:,:,:)
 complex(gwpc) :: rhotwx(3,Dtset%nspinor**2)
 complex(gwpc),allocatable :: rhotwg(:)
 complex(dpc),allocatable :: green_w(:),green_enhigh_w(:)
 complex(gwpc) :: ph_mkt,dum(3)
 complex(gwpc),allocatable :: chi0sf(:,:,:)
 complex(dpc),allocatable :: kkweight(:,:)
 complex(gwpc),allocatable :: wfr1(:),wfr2(:) !what about pointers?
 complex(gwpc),allocatable :: wfwfg(:)
 complex(gwpc),pointer :: wfg1(:),wfg2(:),ug1(:),ug2(:)
 complex(gwpc),allocatable :: fnl(:,:,:,:),fnld(:,:,:,:,:),gradvnl(:,:,:,:)
 type(Cprj_type),allocatable :: Cprj_k(:,:)
!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(2a)')ch10,' cchi0q0 : enter'
 call wrtout(std_out,msg,'PERS')
 call flush_unit(std_out)
#endif
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
 gwpara=MPI_enreg%gwpara
 !
 ! == Copy some values ===
 natom      = Cryst%natom
 nspinor    = Dtset%nspinor

 gmet(:,:)  =Gsph_epsG0%gmet(:,:) 
 gprimd(:,:)=Gsph_epsG0%gprimd(:,:)
 qcart2red = TRANSPOSE(Cryst%rprimd)

 dim_rtwg=1
 if (nspinor==2) then
  !can reduce size depending on Ep%nI and Ep%nj
  dim_rtwg=4
 end if
 spinorwf_pad(:,:)=RESHAPE((/0,0,Wf%npwwfn,Wf%npwwfn,0,Wf%npwwfn,Wf%npwwfn,0/),(/2,4/))
 nab=Wf%nspinor**2
 !
 ! === If required get derivatives of nonlocal operator ===
 if (Dtset%usepaw==0) then 
  SELECT CASE (Ep%inclvkb)
  CASE (0)
   ! === Add only <nk|\nabla|mk> ===
   write(msg,'(4a)')ch10,&
&   ' cchi0q0 : WARNING - ',ch10,&
&   ' neglecting term <n,k|[Vnl,iqr]|m,k> '
   call wrtout(std_out,msg,'COLL')
  CASE (1)
   ! === Legendre polynomials (CPU and mem ~npwwfn^2) ===
   ! * gradvnl = (grad_K+grad_Kp) Vnl(K,Kp) in reciprocal lattice units.
   allocate(gradvnl(3,Ep%npwwfn,Ep%npwwfn,Kmesh%nibz),STAT=istat) 
   if (istat/=0) call memerr(FILE__,'gradvnl',3*Ep%npwwfn**2*Kmesh%nibz,'gwpc')

   call ccgradvnl(Ep%npwwfn,Kmesh%nibz,gvec,gprimd,Kmesh%ibz,Cryst,Psps%mpsang,vkbsign,vkb,vkbd,gradvnl)
  CASE (2) 
   ! === Complex spherical harmonics (CPU and mem \propto npwwfn) ===
   allocate(fnl(Ep%npwwfn,Psps%mpsang**2,natom,Kmesh%nibz),STAT=istat)
   if (istat/=0) call memerr(FILE__,'fnl',Ep%npwwfn*Psps%mpsang**2*natom*Kmesh%nibz,'gwpc')
   allocate(fnld(3,Ep%npwwfn,Psps%mpsang**2,natom,Kmesh%nibz),stat=istat)
   if (istat/=0) call memerr(FILE__,'fnld',3*Ep%npwwfn*Psps%mpsang**2*natom*Kmesh%nibz,'gwpc')
 
   call ccgradvnl_ylm(Ep%npwwfn,Kmesh%nibz,Cryst,gvec,gprimd,Kmesh%ibz,Psps%mpsang,vkbsign,vkb,vkbd,fnl,fnld)
  CASE DEFAULT
   call assert(.FALSE.,'Wrong value of inclvkb',__FILE__,__LINE__)
  END SELECT
 else
  ! === Preliminary computation only in case of PAW+LDA+U ===
  if (Dtset%usepawu/=0) then 
   ! TODO Add LDA+U STUFF to be generalized to nspden==4 and AFM
  end if
 end if
 !
 ! === Initialize the completeness correction === 
 allocate(green_enhigh_w(Ep%nomega)) ; green_enhigh_w=czero
 if (Dtset%gwcomp==1) then
  enhigh=MAXVAL(gwenergy(:,Ep%nbnds,:))+Dtset%gwencomp
  allocate(wfwfg(nfftot*nspinor**2))
  write(msg,'(a,f8.2,a)')' use the completeness correction with the energy ',enhigh*Ha_eV,' [eV] '
  call wrtout(std_out,msg,'COLL')
 end if
 !
 ! === Setup weight (2 for spin unpolarized systems, 1 for polarized) ===
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

 qpoint(:)=Ep%qlwl(:,1) !awful but I have to completely rewrite this part
 if (rank==master) then
  write(*,*)' number of small points ',Ep%nqlwl ; if (Ep%nqlwl/=1) STOP
  write(*,*)' using small q = ',qpoint
 end if
 !
 ! === FFT index of G to be passed to rho_tw_g ===
 igfft0 => igfft(:,Ep%mG0(1)+1,Ep%mG0(2)+1,Ep%mG0(3)+1)

 time_reversal=.FALSE.
 if (PRESENT(Wf_val)) then 
  time_reversal=.TRUE.
  ! Use faster algorithm based on time reversal (presently only in paralell)
  ! Note that special care has to be used in case of metals and/or spin 
  ! since wfrv could also contain unoccupied states 
  write(msg,'(3a)')&
&  ' valence states in real and reciprocal space are in memory ',ch10,&
&  ' using faster equation based on time reversal symmetry '
  call wrtout(std_out,msg,'COLL')
 end if 

 if (rank==master) then 
  write(std_out,*)' symmetrization flag = ',Ep%symchi
  write(std_out,*)' use time reversal   = ',time_reversal
  write(std_out,*)' spectral method     = ',Ep%spmeth
 end if 

 if (Dtset%usepaw==1) then 
  ! === Set up REAL Ylm(q+G) for this q-point ===
  mqmem_=1 ; nqpt_=1 ; optder=0
  allocate(npwarr(nqpt_),qptns(3,nqpt_)) ; npwarr(:)=Ep%npwepG0
  !*** FIXME still do not know wheter use q==0 or q/=0; check npwepG0, look at assemblychi0 
  qptns=zero !;  qptns(:,1)=qpoint(:)
  two_lmaxp1=2*Psps%mpsang-1
  allocate(ylm_q(Ep%npwepG0*mqmem_,two_lmaxp1**2))
  allocate(ylmgr_q(Ep%npwepG0*mqmem_,3+6*(optder/2),two_lmaxp1**2))
  call status(0,Dtfil%filstat,0,level,'call initylmg ')
  call initmpi_seq(MPI_enreg_seq)  
  ! Note: Dtset%nband and Dtset%nsppol are not used in sequential mode
  call initylmg(gprimd,gvec,qptns,mqmem_,MPI_enreg_seq,two_lmaxp1,Ep%npwepG0,Dtset%nband,nqpt_,&
   npwarr,Dtset%nsppol,optder,Cryst%rprimd,Dtfil%unkg,Dtfil%unylm,ylm_q,ylmgr_q)
  !
  ! === Evaluate oscillator matrix elements btw partial waves ===
  allocate(pawrhox(2,Ep%npwepG0,Psps%lmnmax*(Psps%lmnmax+1)/2,Cryst%natom))  
  rhoxsp_method=2
  dim1_rhox=2*(Psps%mpsang-1) !(0:2*(Psps%mpsang-1)%lnmax+1)/2  !Alouani-Method
  dim2_rhox=Psps%lnmax*(Psps%lnmax+1)/2  !Alouani-Method
  if (rhoxsp_method==2)  then 
   dim1_rhox=MAXVAL(Pawtab(:)%l_size)**2
   dim2_rhox=MAXVAL(Pawtab(:)%lmn2_size) !lmn2_size_max
  end if
  call paw_mkrhox(Cryst,pawrhox_spl,gmet,gvec,rhoxsp_method,dim1_rhox,dim2_rhox,&
&  Psps,Pawang,Pawtab,qptns,Ep%npwepG0,ylm_q,pawrhox)

  allocate(dimlmn(Cryst%natom)) 
  do iat=1,Cryst%natom
   dimlmn(iat)=Pawtab(Cryst%typat(iat))%lmn_size
  end do
  allocate(Cprj_k(Cryst%natom,nspinor*Ep%nbnds)) ; call cprj_alloc(Cprj_k,0,dimlmn)
 end if

 allocate(rhotwg(Ep%npwe*nspinor**2),tabr_k(nfftot))
 allocate(wfr1(Wf%nfftot*nspinor),wfr2(Wf%nfftot*nspinor))

 SELECT CASE (Ep%spmeth) 
 CASE (0)
  allocate(green_w(Ep%nomega))
  write(msg,'(2a)')ch10,' Calculating chi0(q=(0,0,0),omega,G,G")'    ; call wrtout(std_out,msg,'COLL')
 CASE (1,2) 
  write(msg,'(2a)')ch10,' Calculating Im chi0(q=(0,0,0),omega,G,G")' ; call wrtout(std_out,msg,'COLL')
  !
  ! === Find max and min resonant transitions for this q, report values for this processor ===
  call make_transitions(1,Ep%nbnds,nbvw,Dtset%nsppol,Ep%symchi,Cryst%timrev,GW_TOL_DOCC,Ep%zcut,&
&  max_rest,min_rest,my_max_rest,my_min_rest,Kmesh,Ltg_q,MPI_enreg,Ep%mG0,gwenergy,occ,(/zero,zero,zero/))
  if (gwpara==1) STOP " make_transitions is buggy" !there is a problem in make_transitions due to MPI_enreg
  !
  ! === Calculate frequency dependent weights for Kramers Kronig transform ===
  allocate(omegasf(Ep%nomegasf),kkweight(Ep%nomega,Ep%nomegasf))
  !my_wl=1 ; my_wr=Ep%nomegasf
  call setup_spectral(Ep%nomega,Ep%omega,Ep%nomegasf,omegasf,max_rest,min_rest,my_max_rest,my_min_rest,&
&  0,Ep%zcut,zero,my_wl,my_wr,kkweight)

  if (gwpara==2.and..not.PRESENT(Wf_val)) then 
   write(msg,'(a)')' valence wfs in r-space are not in memory, not coded yet'
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end if 
  write(msg,'(2a,2i5,a,i5,2a)')ch10,&
&  ' cchi0q0: allocating chi0sf using my_wl and my_wr = ',my_wl,my_wr,' ( ',my_wr-my_wl+1,' )',ch10
  call wrtout(std_out,msg,'PERS')
  allocate(chi0sf(Ep%npwe,Ep%npwe,my_wl:my_wr),stat=istat)
  if (istat/=0) call memerr(FILE__,'chi0sf',Ep%npwe**2*(my_wr-my_wl+1)*Ep%nomegasf,'spc')
  chi0sf=czero_gw
 CASE DEFAULT
  call assert(.FALSE.,'Wrong value of spmeth',__FILE__,__LINE__)
 END SELECT

 nkpt_summed=Kmesh%nbz 
 if (Ep%symchi/=0) then 
  nkpt_summed=Ltg_q%nibz_ltg 
  call print_little_group(Ltg_q,std_out,Dtset%prtvol,'COLL')
 end if
 write(msg,'(a,i6,a)')' Calculation status ( ',nkpt_summed,' to be completed) :'
 call wrtout(std_out,msg,'COLL')
 !
 ! === Loop over k-points in BZ ===
 chi0(:,:,:)=czero_gw
 chi0sumrule(:)=zero

 do ik_bz=1,Kmesh%nbz 

  if (Ep%symchi==1) then 
   if (Ltg_q%ibzq(ik_bz)/=1) CYCLE ! Only IBZ_q 
  end if

  call status(ik_bz,Dtfil%filstat,0,level,'loop ikpt     ')
  !
  ! === Get ik_ibz, non-symmorphic phase and symmetries from ik_bz === 
  call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym,itim,ph_mkt)
  tabr_k=ktabr(:,ik_bz) ! Table for rotated FFT points
  spinrot_kbz(:)=Cryst%spinrot(:,isym)
  !
  ! === Loop on spin to calculate $\chi_{up,up} + \chi_{down,down}$ ===
  do is=1,Ep%nsppol
   !
   ! Parallelization over k-points in the full BZ
   ! Note that spin para is not active but the check should be done here
   ! array proc_distrb does not depend on band index (there is a check in screening)
   if (gwpara==1.and.nprocs>1) then
    if (MPI_enreg%proc_distrb(ik_bz,1,is)/=rank) CYCLE
!DEBUG for the time being do not comment
    if (ANY(MPI_enreg%proc_distrb(ik_bz,:,is)==-999)) then 
     write(*,*)'cchi0q0 : wrong distrb' ; call leave_new('COLL')
    end if 
!ENDDEBUG
   end if
   write(msg,'(2(a,i4),a,i2,a,i3)')&
&    ' ik = ',ik_bz,' / ',Kmesh%nbz,' is = ',is,' done by processor ',rank
   call wrtout(std_out,msg,'PERS') 

   if (Dtset%usepaw==1) then 
    ! === Load cprj for this k-point ===
    shift=nspinor*Ep%nbnds*Kmesh%nbz*(is-1)
    indx_k_bz=nspinor*Ep%nbnds*(ik_bz-1)+shift
    ibsp=0
    do ib=1,Ep%nbnds
     do ispinor=1,nspinor
      ibsp=ibsp+1
      do iat=1,Cryst%natom
       Cprj_k(iat,ibsp)%cp(:,:)=Cprj_bz(iat,indx_k_bz+ibsp)%cp(:,:)
      end do
     end do
    end do
   end if
   !
   ! /***********************************************************************/
   !
   ! Conventions: 1) a symmetry in real space acts as R_t f(r) = f(R^-1(r-t))
   !              2) S=\transpose R^-1 
   !              3) kbz=S kibz 
   !
   !  The wavefunctions for the k-point in the BZ are (assuming nondegenerate states):
   !
   !  u(G,b, Sk) = u ( S^-1G,b,k)* e^{-i(Sk+G)*t) 
   !  u(G,b,-Sk) = u*(-S^-1G,b,k)* e^{ i(Sk-G)*t) 
   !
   !  u(r,b, Sk) = u (R^-1(r-t),b,k) e^{-iSk*t} 
   !  u(r,b,-Sk) = u*(R^-1(r-t),b,k) e^{ iSK*t} 
   !
   !  The gradient of Vnl(K,Kp) for the k-point in the BZ should be:
   !   gradvnl(SG,SGp,Sk)=S gradvnl(G,Gp,kibz)
   ! FIXME should check the expression in case of non zero tnons but Im lazy now 
   !
   ! /***********************************************************************/
   !
   ! ==== Loop over "conduction" states ===
   do ib1=1,Ep%nbnds 
    if (gwpara==2) then
     if (MPI_enreg%proc_distrb(ik_ibz,ib1,is)/=rank) CYCLE
    end if     
    if(ik_bz==1) then
     write(msg,'(2(a,i4),a,i2,a,i3)')&
&    ' ib1 = ',ib1,' / ',Ep%nbnds,' is = ',is,' done by processor ',rank
     call wrtout(std_out,msg,'PERS') 
    endif
    !   
    ! === Loop over "valence" states ===
    do ib2=1,Ep%nbnds 
     !MG FIXME this has been commented by Fabien 
     !Note that completeness uses Wf_val% of ib2
     if (MPI_enreg%gwpara==2.and.Dtset%gwcomp==0) then
      if (MPI_enreg%proc_distrb(ik_ibz,ib2,is)/=rank) CYCLE
     end if

     deltaf_b1b2=spin_fact*(occ(ik_ibz,ib1,is)-occ(ik_ibz,ib2,is))
     !
     ! === Skip negligible transitions ===
     if (Dtset%gwcomp==0) then
      if (ABS(deltaf_b1b2) < GW_TOL_DOCC) CYCLE 
     else
      ! when the completeness trick is used,
      ! we need to also consider transitions with vanishing deltaf
      if (occ(ik_ibz,ib2,is) < GW_TOL_DOCC) CYCLE  
     end if
     deltaeKS_b1b2=   energy(ik_ibz,ib1,is)-  energy(ik_ibz,ib2,is)
     deltaeGW_b1b2= gwenergy(ik_ibz,ib1,is)-gwenergy(ik_ibz,ib2,is)

     SELECT CASE (Ep%spmeth)
     CASE (0) 
      ! === Adler-Wiser expression === 
      ! Add small imaginary of the Time-Ordered resp function but only for non-zero real omega
      ! FIXME What about metals?
      SELECT CASE (time_reversal)
      CASE (.TRUE.)
       if (Dtset%gwcomp==0) then ! cannot be completely skipped in case of completeness correction
        if (ib1<ib2) CYCLE ! Here we GAIN a factor ~2 
       end if

       do io=1,Ep%nomega 
        if(abs(deltaeGW_b1b2)>GW_TOL_W0) green_w(io) = mkG0w(Ep%omega(io),deltaf_b1b2,deltaeGW_b1b2,Ep%zcut,GW_TOL_W0)

        if (Dtset%gwcomp==1) then
         ! * Calculate the completeness correction
         numerator= -spin_fact*occ(ik_ibz,ib2,is)
         deltaeGW_enhigh_b2=enhigh-gwenergy(ik_ibz,ib2,is)
         ! Completeness correction is NOT valid for real frequencies
         if (REAL(Ep%omega(io))<GW_TOL_W0) then
          green_enhigh_w(io) = mkG0w(Ep%omega(io),numerator,deltaeGW_enhigh_b2,Ep%zcut,GW_TOL_W0)
         else
          green_enhigh_w(io) = czero_gw
         endif
         if (deltaf_b1b2<0.d0) then
          green_w(io)= green_w(io) - green_enhigh_w(io)
         else 
          ! Disregard green_w, since it is already accounted for through the time-reversal
          green_w(io)=             - green_enhigh_w(io)    
         end if
        end if !gwcomp==1

       end do !io

       ! === Add the "delta part", symmetrization is done inside the routine ===
       ! MG Be careful here as time-reversal itim is not the same as tabi 
       ! TODO Add PAW
       if (Dtset%gwcomp==1.and.ib1==ib2) then
        call get_wfr(Wf_val,MPI_enreg,ib2,ik_ibz,is,wfr2)
        call calc_wfwfg(MPI_enreg,Dtset%paral_kgb,tim_fourdp,tabr_k(:),Kmesh%tabi(ik_bz),&
&                       nfftot,ngfft,wfr2,wfr2,wfwfg)
        qzero=.TRUE.
        call completechi0_deltapart(ik_bz,qzero,Ep%symchi,Ep%npwe,Ep%npwvec,Ep%nomega,nspinor,&
&        nfftot,ngfft,gvec,igfft0,Gsph_wfn,Ltg_q,green_enhigh_w,wfwfg,chi0)
       end if

      CASE (.FALSE.)
       ! === Adler-Wiser without time-reversal ===
       sgn=deltaeGW_b1b2/ABS(deltaeGW_b1b2) 
       do io=1,Ep%nomega
        if (REAL(Ep%omega(io))>GW_TOL_W0) then
         green_w(io)= deltaf_b1b2 / (Ep%omega(io)+deltaeGW_b1b2-(0.,1.)*sgn*Ep%zcut)
        else
         green_w(io)= deltaf_b1b2 / (Ep%omega(io)+deltaeGW_b1b2)
        end if
       end do
      END SELECT
     CASE (1,2)
      ! === Spectral method, here time-reversal is always assumed ===
      if (deltaeGW_b1b2<0) CYCLE
      call approxdelta(Ep%nomegasf,omegasf,deltaeGW_b1b2,Ep%spsmear,iomegal,iomegar,wl,wr,Ep%spmeth)
     END SELECT 
     !
     ! === Form rho-twiddle(r)=u^*_{b1,k}(r) u_{b2,k}(r) and its FFT transform ===
     ! * Form also rho-twiddle(b1,b2,q,G=0) using small vector q and k.p perturbation theory
     SELECT CASE (time_reversal)
     CASE (.TRUE.)
      if (ib1>nbvw) then
       if (ib2<=nbvw) then 
        wfg1 => Wf%wfg(:,ib1,ik_ibz,is) ; wfg2 => Wf_val%wfg(:,ib2,ik_ibz,is)
        call get_wfr(Wf,    MPI_enreg,ib1,ik_ibz,is,wfr1)
        call get_wfr(Wf_val,MPI_enreg,ib2,ik_ibz,is,wfr2)
       else ! Unlikely
        wfg1 => Wf%wfg(:,ib1,ik_ibz,is) ; wfg2 => Wf%wfg(:,ib2,ik_ibz,is)
        call get_wfr(Wf,MPI_enreg,ib1,ik_ibz,is,wfr1)
        call get_wfr(Wf,MPI_enreg,ib2,ik_ibz,is,wfr2)
       end if
      else
       if (ib2>nbvw) STOP "BUG"
       if (rank/=master.and.Dtset%gwcomp==0) CYCLE ! Only master for "valence-valence" transitions
       wfg1 => Wf_val%wfg(:,ib1,ik_ibz,is) ; wfg2 => Wf_val%wfg(:,ib2,ik_ibz,is)
       call get_wfr(Wf_val,MPI_enreg,ib1,ik_ibz,is,wfr1)
       call get_wfr(Wf_val,MPI_enreg,ib2,ik_ibz,is,wfr2)
      end if
     CASE (.FALSE.)
      wfg1 => Wf%wfg(:,ib1,ik_ibz,is) ; wfg2 => Wf%wfg(:,ib2,ik_ibz,is)
      call get_wfr(Wf,MPI_enreg,ib1,ik_ibz,is,wfr1)
      call get_wfr(Wf,MPI_enreg,ib2,ik_ibz,is,wfr2)
     END SELECT 

     call rho_tw_g(Dtset%paral_kgb,nspinor,Ep%npwe,nfftot,ngfft,1,igfft0,&
&     wfr1,itim,tabr_k,ph_mkt,spinrot_kbz,&
&     wfr2,itim,tabr_k,ph_mkt,spinrot_kbz,&
&     dim_rtwg,rhotwg,tim_fourdp,MPI_enreg)
     !
     ! === Add PAW onsite contribution, projectors are already in BZ ===
     if (Dtset%usepaw==1) then 
! FIXME Find a clever way to deal with spinors
      spad=(nspinor-1)
      i1=ib1 ; if (nspinor==2) i1=(2*ib1-1)
      i2=ib2 ; if (nspinor==2) i2=(2*ib2-1)
      call paw_rho_tw_g(Ep%npwe,dim_rtwg,nspinor,Cryst%natom,Psps%lmnmax,dimlmn,&
&      Cprj_k(:,i1:i1+spad),Cprj_k(:,i2:i2+spad),pawrhox,rhotwg)
     end if
     !
     ! === -iq.<c,k|\nabla|v,k> is always calculated ===
     ! * nonsymmorphic phases cancel each other.
     ! * For nspinor==2, we neglect off-diagonal elements which require rv12 and rv21
     !
     ! 1) Plane wave contribution for -i\nabla
     rhotwx(:,:)=czero_gw
     do iab=1,nab
      spad1 = spinorwf_pad(1,iab)
      spad2 = spinorwf_pad(2,iab)
      ug1 => wfg1(spad1+1:spad1+Wf%npwwfn)
      ug2 => wfg2(spad2+1:spad2+Wf%npwwfn)
      do ig=1,Wf%npwwfn
       ct=CONJG(ug1(ig))*ug2(ig)
       rhotwx(:,iab)=rhotwx(:,iab)+gvec(:,ig)*ct
      end do
     end do
     ! 2) Onsite contribution for PAW
     ! TODO check this part, mind reduced coordinates
     ! FIXME Find a clever way to deal with spinors
     if (Dtset%usepaw==1) then 
      ! This should be a function to allow inling but check g95
      spad=(nspinor-1)
      i1=ib1 ; if (nspinor==2) i1=(2*ib1-1)
      i2=ib2 ; if (nspinor==2) i2=(2*ib2-1)
      call paw_inabla(Cryst%natom,Cryst%typat,Pawtab,Cprj_k(:,i1),Cprj_k(:,i2),igradcart_paw) 
      igradred_paw=igradcart_paw
      !igradred_paw(1,:)=MATMUL(qcart2red,igradcart_paw(1,:))
      !igradred_paw(2,:)=MATMUL(qcart2red,igradcart_paw(2,:))
      ! TODO Add spinorial case
      rhotwx(:,1)=rhotwx(:,1)+CMPLX(igradred_paw(1,:),igradred_paw(2,:),gwpc)
     end if
     !
     ! === For PPS add term <c,k|[Vnl,iqr]|v,k> only if required ===
     ! * Two different algorithms are coded, the second one is much faster
     ! * Must multiply by ucvol because yet gradvnl contain 1/ucvol
     SELECT CASE (Ep%inclvkb)
     CASE (1)
      ! === Legendre Polynomials. CPU and MEM ~ npwwfn**2 ===
      ! FIXME check the phase in case of simmorphic, 
      ! TODO Add spinorial case
      call apply_gradvnl(Ep%npwwfn,wfg1,wfg2,gradvnl(:,:,:,ik_ibz),dum) !this should be a function
      rhotwx(:,1)=rhotwx(:,1)+dum(:)
     CASE (2)
      ! === Complex spherical harmonics (much faster!) ===
      ! FIXME check the phase in case of simmorphic ; check better timereversal case
      ! TODO Add spinorial case
      call apply_gradvnl_Ylm(Wf%npwwfn,wfg1,wfg2,Cryst%natom,Psps%mpsang,fnl(:,:,:,ik_ibz),fnld(:,:,:,:,ik_ibz),dum) 
      rhotwx(:,1)=rhotwx(:,1)+dum(:)
     END SELECT

     SELECT CASE (Ep%spmeth) 
     CASE (0) 
      ! ==== Adler-Wiser expression, to be consistent here we use the KS eigenvalues (?) ===
      if(abs(deltaeKS_b1b2)>GW_TOL_W0) then
       rhotwx(:,:)=-rhotwx(:,:)/deltaeKS_b1b2
      else
       rhotwx(:,:)=czero_gw
      endif

      call assemblychi0q0_sym(qpoint,ik_bz,isym,itim,Dtset%gwcomp,nspinor,Ep%npwepG0,Ep,&
&      Cryst,Ltg_q,Gsph_epsG0,chi0,rhotwx(:,1),rhotwg,green_w,green_enhigh_w,deltaf_b1b2) 

     CASE (1,2)
      ! === Spectral method, to be consistent here we use the KS eigenvalues ===
      rhotwx(:,:)=-rhotwx(:,:)/deltaeKS_b1b2
      call assemblychi0sfq0(qpoint,ik_bz,isym,itim,nspinor,Ep%symchi,Ep%npwepG0,Ep%npwe,Cryst,Ltg_q,Gsph_epsG0,&
&      deltaf_b1b2,my_wl,iomegal,wl,my_wr,iomegar,wr,rhotwx(:,1),rhotwg,Ep%nomegasf,chi0sf)

     CASE DEFAULT
      call assert(.FALSE.,'Wrong value of spmeth',__FILE__,__LINE__)
     END SELECT
     !
     ! Accumulating the sum rule on chi0 Eq. (5.284) in G. D. Mahan Many-Particle Physics 3rd edition
     !chi0sumrule(:)=chi0sumrule(:) + spin_fact * occ(ik_ibz,ib2,is) * deltaeGW_b1b2 * abs(rhotwg(:))**2
     factor=spin_fact*occ(ik_ibz,ib2,is)
     call accumulate_chi0sumrule(ik_bz,Ep%symchi,Ep%npwe,factor,deltaeGW_b1b2,&
&     Ltg_q,Gsph_epsG0,Ep%npwepG0,rhotwg,chi0sumrule)
     !
     ! Include  also the completeness correction in the sum rule
     if (Dtset%gwcomp==1) then
      !chi0sumrule(:)=chi0sumrule(:) - spin_fact * occ(ik_ibz,ib2,is) * deltaeGW_enhigh_b2 * abs(rhotwg(:))**2
      factor=-spin_fact*occ(ik_ibz,ib2,is)
      call accumulate_chi0sumrule(ik_bz,Ep%symchi,Ep%npwe,factor,deltaeGW_enhigh_b2,&
&      Ltg_q,Gsph_epsG0,Ep%npwepG0,rhotwg,chi0sumrule)
      !if (ib1==Ep%nbnds) chi0sumrule(:)=chi0sumrule(:) + spin_fact * occ(ik_ibz,ib2,is) * deltaeGW_enhigh_b2
      if (ib1==Ep%nbnds) chi0sumrule(:)=chi0sumrule(:) + wtk_ltg(ik_bz)*spin_fact*occ(ik_ibz,ib2,is)*deltaeGW_enhigh_b2
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
  ! First coding, the loop over ig igp could be optimised taking into account symmetries 
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
  ! == Sum contributions from each proc ===
  ! * Looping on frequencies to avoid problems with the size of the MPI packet
  do io=1,Ep%nomega
   call xsum_mpi(chi0(:,:,io),spaceComm,ierr)
  end do 
  call leave_test(MPI_enreg)
  chi0(:,:,:)=weight*chi0(:,:,:)
 CASE DEFAULT
  call assert(.FALSE.,'Wrong value for spmeth', __FILE__,__LINE__)
 END SELECT
 !
 ! Apply MPI for the sum rule
 call xsum_mpi(chi0sumrule(:),spaceComm,ierr)
 chi0sumrule(:)=chi0sumrule(:) * pi * weight   ! the pi factor comes from Im[1/(x-ieta)] = pi delta(x)
 ! /*** Now master has chi0(q,G,Gp,Ep%omega) ***/

 if (Cryst%use_antiferro) then 
  call symmetrize_afm_chi0(Cryst,Gsph_epsG0,Ep%npwe,Ep%nomega,chi0)
 end if
 !
 ! === Deallocate memory ===
 deallocate(rhotwg,tabr_k)
 deallocate(wfr1,wfr2)

 if (Dtset%usepaw==1) then 
  deallocate(npwarr,qptns,dimlmn)
  deallocate(ylm_q,ylmgr_q,pawrhox)
  call cprj_free(Cprj_k) ; deallocate(Cprj_k)
 end if
 if (Ep%inclvkb==1) deallocate(gradvnl)
 if (Ep%inclvkb==2) deallocate(fnl,fnld)

 if (allocated(green_enhigh_w)) deallocate(green_enhigh_w)
 if (allocated(wfwfg   ))       deallocate(wfwfg   )
 if (allocated(kkweight))       deallocate(kkweight)
 if (allocated(omegasf ))       deallocate(omegasf )
 if (allocated(green_w ))       deallocate(green_w )
 if (allocated(chi0sf  ))       deallocate(chi0sf  )

 call status(0,Dtfil%filstat,0,level,'exit          ')
 call pclock(180)

#if defined DEBUG_MODE
 write(msg,'(a)')' cchi0q0 : exit '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine cchi0q0
!!***

!TODO  The following routines should be redefined as functions to allow inlining. 
!      Unfortunately g95 crashes if inabla et al are defined as functions

!! NAME
!! paw_inabla
!!
!! FUNCTION
!!  Calculate the PAW onsite contribution to the matrix elements of the  i\nabla operator.
!!
!! INPUTS
!!  natom=number of atoms in unit cell
!!  Pawtab(ntypat)=Only for PAW, TABulated data initialized at start
!!    %lmn_size Number of (l,m,n) elements for the paw basis
!!    %nabla_ij(3,lmn_size,lmn_size)) Onsite contribution
!!     <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j> for each type
!!
!! OUTPUT
!!  onsite(2,3)=Onsite contribution to  $i<wfg1|\nabla|wfg2>$
!!
subroutine paw_inabla(natom,typat,Pawtab,Cprj_b1,Cprj_b2,onsite)
    
 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
 real(dp),intent(out) :: onsite(2,3)
!arrays
 integer,intent(in) :: typat(natom)
 type(Pawtab_type),intent(in) :: Pawtab(:)
 type(Cprj_type),intent(in) :: Cprj_b1(natom),Cprj_b2(natom)

!Local variables-------------------------------
 integer :: ig,iatom,itypat,lmn_size,ilmn,jlmn
 real(dp) :: re_p,im_p
!arrays
 real(dp),pointer :: nabla_ij(:,:,:)
! *************************************************************************
 
  onsite(:,:)=zero
  do iatom=1,natom 
   itypat=typat(iatom)
   lmn_size=Pawtab(itypat)%lmn_size
   nabla_ij => Pawtab(itypat)%nabla_ij(:,:,:) 
   !
   !=== Unpacked loop over lmn channels ====
   do jlmn=1,lmn_size
    do ilmn=1,lmn_size

     re_p =  Cprj_b1(iatom)%cp(1,ilmn)*Cprj_b2(iatom)%cp(1,jlmn) &
&           +Cprj_b1(iatom)%cp(2,ilmn)*Cprj_b2(iatom)%cp(2,jlmn) 
 
     im_p =  Cprj_b1(iatom)%cp(1,ilmn)*Cprj_b2(iatom)%cp(2,jlmn) &
&           -Cprj_b1(iatom)%cp(2,ilmn)*Cprj_b2(iatom)%cp(1,jlmn)

     onsite(1,1)=onsite(1,1) - im_p*nabla_ij(1,ilmn,jlmn)
     onsite(1,2)=onsite(1,2) - im_p*nabla_ij(2,ilmn,jlmn)
     onsite(1,3)=onsite(1,3) - im_p*nabla_ij(3,ilmn,jlmn)

     onsite(2,1)=onsite(2,1) + re_p*nabla_ij(1,ilmn,jlmn)
     onsite(2,2)=onsite(2,2) + re_p*nabla_ij(2,ilmn,jlmn)
     onsite(2,3)=onsite(2,3) + re_p*nabla_ij(3,ilmn,jlmn)

    end do !ilmn
   end do !jlmn
  end do !iatom

end subroutine paw_inabla
!!***

!! NAME
!! apply_gradvnl
!!
!! FUNCTION
!!  Simple function to calculate matrix elements of the gradient of the non-localc operator 
!!  when Legendre polynomial are employed. Wavefunctions are supposed to be complex.
!!
!! INPUTS
!!  gradvnl(3,npwwfn,npwwfn)= the gradient at this k-point
!!  npwwfn=number of G vectors for wavefunctions
!!  wfg1(npwwfn),wfg1(npwwfn)= bra and ket in reciprocal space
!!
!! OUTPUT
!!  res(3)= sum_{G_1,G_2} c(G_1)^\dagger c(G_2) V^{nl}_{G_1,G_2}
!!
subroutine apply_gradvnl(npwwfn,wfg1,wfg2,gradvnl,res) 
    
 use defs_basis
 use m_gwdefs, only : czero_gw

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: npwwfn
 complex(gwpc),intent(in) :: wfg1(npwwfn),wfg2(npwwfn)
 complex(gwpc),intent(in) :: gradvnl(3,npwwfn,npwwfn)
 complex(gwpc),intent(out) :: res(3)

!Local variables-------------------------------
 integer :: ig1,ig2 
 complex(gwpc) :: ct
! *************************************************************************
 
 res(:)=czero_gw
 do ig1=1,npwwfn
  do ig2=1,npwwfn
   ct=CONJG(wfg1(ig1))*wfg2(ig2)
   res(:)=res(:)+ct*gradvnl(:,ig1,ig2)
  end do
 end do

end subroutine apply_gradvnl
!!***


!! NAME
!! apply_gradvnl_Ylm
!!
!! FUNCTION
!!  Simple function to calculate matrix elements of the gradient of the non-local operator 
!!  when Legendre polynomial are employed. Wavefunctions are supposed to be complex.
!!
!! INPUTS
!!  gradvnl(3,npwwfn,npwwfn)= the gradient at this k-point
!!  npwwfn=number of G vectors for wavefunctions
!!  wfg1(npwwfn),wfg1(npwwfn)= bra and ket in reciprocal space
!!  fnl(npwwfn,mpsang**2,natom)= Kleynmann-Bylander form factor for this k-point
!!  fnld(3,npwwfn,mpsang**2,natom)= Derivative of the KB form factor for this k-point
!!
!! OUTPUT
!!  
!!
subroutine apply_gradvnl_Ylm(npwwfn,wfg1,wfg2,natom,mpsang,fnl,fnld,res) 
    
 use defs_basis
 use m_gwdefs, only : czero_gw

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: npwwfn,natom,mpsang
 complex(gwpc),intent(in) :: wfg1(npwwfn),wfg2(npwwfn)
 complex(gwpc),intent(in) :: fnl(npwwfn,mpsang**2,natom) 
 complex(gwpc),intent(in) :: fnld(3,npwwfn,mpsang**2,natom)
 complex(gwpc),intent(out) :: res(3)

!Local variables-------------------------------
 integer :: iat,ig,ilm
 complex(gwpc) :: cta1,cta4
 complex(gwpc) :: cta2(3),cta3(3)
! *************************************************************************
 
 res(:)=czero_gw

 do iat=1,natom 
  do ilm=1,mpsang*mpsang 
   cta1=czero_gw ; cta2(:)=czero_gw
   cta4=czero_gw ; cta3(:)=czero_gw
   ! === Here we take advantage of the property Y_(l-m)= (-i)^m Y_lm^* ===
   do ig=1,npwwfn
    cta1   = cta1    + wfg1(ig) * fnl(ig,ilm,iat)
    cta2(:)= cta2(:) + wfg2(ig) * fnld(:,ig,ilm,iat)
    cta3(:)= cta3(:) + wfg1(ig) * fnld(:,ig,ilm,iat)
    cta4   = cta4    + wfg2(ig) * fnl(ig,ilm,iat)
   end do 
   res(:)= res(:) +CONJG(cta1)*cta2(:) +CONJG(cta3(:))*cta4
  end do !ilm
 end do !iat

end subroutine apply_gradvnl_Ylm
!!***
