!{\src2tex{textfont=tt}}
!!****f* ABINIT/csigme
!! NAME
!! csigme
!!
!! FUNCTION
!! Calculate diagonal and off-diagonal matrix elements of the self-energy operator
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (FB, GMR, VO, LR, RWG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! Cprj_bz(natom,Dtset%nspinor*Sp%nbnds*Kmesh%nbz*Sp%nsppol) <type(cprj_type)>=
!!  projected input wave functions <Proj_i|Cnk> with all NL projectors for each k-point in the full Brillouin zone
!! Dtfil=filenames (only Dtfil%filsrc is used if the screening must be read from file)
!! Dtset <type(dataset_type)>=all input variables in this dataset
!!    %accesswff
!!    %localrdwf
!!    %paral_kgb 
!!    %nspinor=Number of spinorial components.
!!    %usepaw=1 for PAW, 0 for NC pseudopotentials.
!!    %gwcomp=If 1 use an extrapolar approximation to accelerate convergence.
!!    %gwencomp=Extrapolar energy.
!! efermi_qp=KS or QP Fermi energy
!! en_qp(Kmesh%nibz,Sp%nbnds,Sp%nsppol)=KS or QP energies for k-points, bands and spin
!! Er <Epsilonm1_results> (see the definition of this structured datatype)
!!    %mqmem=if 0 use out-of-core method in which a single q-slice of espilon is read inside the loop over k
!!    %nomega_i=Number of imaginary frequencies.
!!    %nomega_r=Number of real frequencies.
!!    %nomega=Total number of frequencies.
!! Gsph_Max<Gvectors_type>= info on biggest G-sphere
!!    %nsym=number of symmetry operations
!!    %rottb(Sp%npwvec,timrev,nsym)=index of (IS) G where I is the identity or the inversion
!!      operation and G is one of the npwvec vectors in reciprocal space
!!    %timrev=2 if time-reversal symmetry is used, 1 otherwise
!! gvec(3,Sp%npwvec)=integer coordinates of each plane wave in reciprocal space
!! igfft(Sp%npwvec,2*Sp%mG0(1)+1,2*Sp%mG0(2)+1,Sp%mG0(3)+1)=index of G-G0 in the FFT grid 
!! ikcalc=index in the array Sp%kcalc of the k-point where GW corrections are calculated 
!! Ltg_k datatype containing information on the little group
!! ktabr(nfftgw_tot,Kmesh%nbz)= index of R^{-1}(r-t) in the FFT array where kBZ=(IS) kIBZ
!! Kmesh <BZ_mesh_type> 
!!    %nbz=Number of points in the BZ
!!    %nibz=Number of points in IBZ
!!    %kibz(3,nibz)=k-point coordinates, irreducible Brillouin zone
!!    %kbz(3,nbz)=k-point coordinates, full Brillouin zone
!!    %ktab(nbz)= table giving for each k-point in the BZ (kBZ), the corresponding
!!    %ktabi(nbz)= for each k-point in the BZ defines whether inversion has to be considered 
!!    %ktabp(nbz)= phase factor associated to tnons
!! MPI_enreg= datatype gathering information on parallelism, variables used 
!!    %gwpara= 0 no parallelism 
!!             1 parallelism is over k-points                                     
!!             2 bands are spread btw processors
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nfftgw_tot=number of points of the FFT grid for GW wavefunctions
!! occ_qp(Kmesh%nibz,Sp%nbnds,Sp%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!! Vcp <Coulombian_type datatype> containing information on the cutoff technique
!!    %vc_sqrt(npwx,nqibz)= square-root of the coulombian potential for q-points in the IBZ
!! Pawang<pawang_type> angular mesh discretization and related data:
!! Pawtab(Psps%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! Qmesh <bz_mesh_type> : datatype gathering information of the q-mesh used
!!    %ibz=q points where $\tilde\epsilon^{-1}$ has been computed
!!    %bz(3,nqbz)=coordinates of all q-points in BZ
!!    %qtab(nqbz)= table giving for each q-point in the BZ (qBZ), the corresponding
!!      irreducible point (qIBZ), where qBZ= (IS) qIBZ and I is either the inversion or the identity operation
!!    %qtabi(nqbz)= for each q-point in the BZ defines whether the inversion has to be considered 
!!    %qtabo(nqbz)= the symmetry operation S that takes qIBZ to each qBZ
!! Sp <sigma_parameters> (see the definition of this structured datatype)
!!    %npwvec= Max betwee npweps and npwwfn used to dimension arrays
!! Cryst<Crystal_structure>=Info on unit cell and symmetries
!!    %natom=number of atoms in unit cell
!!    %ucvol=unit cell volume
!!    %nsym=number of symmetry operations
!!    %typat(natom)=type of each atom
!!  much slower but it requires less memory
!! pawrhox(2,Sp%npx,Psps%lmnmax*(Psps%lmnmax+1)/2,natom,Qmesh%nbz*Dtset%usepaw): array containing 
!!   $<phj/r|e^{-i(q+G)}|phi/r>-<tphj/r|e^{-i(q+G)}|tphi/r>$ in packed form.
!! PPm<PPmodel_type>= Datatype gathering information on the Plasmonpole technique (see also PPmodel_symmetrize).
!!    %model=type of ppmodel      
!!    %npwc=number of G-vectors for the correlation part.
!!    %dm2_otq =size of second dimension of otq array
!!    %dm2_bots=size of second dimension of botsq arrays
!!    %dm_eig  =size of second dimension of eig arrays
!!
!! OUTPUT
!!  Sr=sigma_results (see the definition of this structured datatype)
!!
!! NOTES
!!  The treatment of the divergence of Gygi+Baldereschi (PRB 1986) is included.
!!  The calculation of energy derivative is based on finite elements.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      calc_coh,calc_sig_noppm,calc_sig_ppm,cgemv,cggfft,diago_hamilt,distrb2
!!      findqg0,mpi_barrier,rho_tw_g,timab,xcomm_init,xmaster_init
!!      xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine csigme(ikcalc,Dtset,Cryst,Dtfil,Sp,Sr,Er,Gsph_Max,Vcp,Kmesh,Qmesh,Ltg_k,PPm,&
& Pawang,Pawtab,Psps,Cprj_bz,Wf_info,Wf_info_braket,MPI_enreg,pawrhox,gvec,ktabr,ngfft,igfft,&
& nfftgw_tot,en_qp,occ_qp,efermi_qp,my_minb,my_maxb,ngfftf,nfftf,rhor)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : czero_gw, cone_gw, j_gw
 use m_numeric_tools, only : linfit, pade, dpade, newrap_step
 use m_errors, only : assert
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13nonlocal
 use interfaces_15gw, except_this_one => csigme
 use interfaces_lib00numeric
 use interfaces_lib01hidempi
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikcalc,my_maxb,my_minb,nfftgw_tot,nfftf
 real(dp),intent(in) :: efermi_qp
 type(Crystal_structure),intent(in) :: Cryst
 type(MPI_type),intent(inout) :: MPI_enreg
 type(BZ_mesh_type),intent(in) :: Kmesh,Qmesh
 type(Coulombian_type),intent(in) :: Vcp
 type(Datafiles_type),intent(in) :: Dtfil
 type(Dataset_type),intent(in) :: Dtset
 type(Epsilonm1_results),intent(inout) :: Er
 type(Gvectors_type),intent(in) :: Gsph_Max
 type(Little_group),intent(in) :: Ltg_k
 type(PPmodel_type),intent(inout) :: PPm
 type(Pseudopotential_type),intent(in) :: Psps
 type(Pawang_type),intent(in) :: Pawang
 type(Sigma_parameters),intent(in) :: Sp
 type(Sigma_results),intent(inout) :: Sr
 type(Wavefunctions_information),intent(inout) :: Wf_info,Wf_info_braket
!arrays
 integer,intent(in) :: gvec(3,Sp%npwvec)
 integer,intent(in) :: ngfft(18),ngfftf(18)
 integer,intent(in),target :: igfft(Sp%npwvec,2*Sp%mG0(1)+1,2*Sp%mG0(2)+1,2*Sp%mG0(3)+1)
 integer,intent(in) ::  ktabr(nfftgw_tot,Kmesh%nbz)
 real(dp),intent(in) :: en_qp(Kmesh%nibz,Sp%nbnds,Sp%nsppol)
 real(dp),intent(in) :: occ_qp(Kmesh%nibz,Sp%nbnds,Sp%nsppol)
 real(dp),intent(in) :: pawrhox(2,Sp%npwx,Psps%lmnmax*(Psps%lmnmax+1)/2,Cryst%natom,Qmesh%nbz*Dtset%usepaw) 
 real(dp),intent(inout) :: rhor(nfftf,Dtset%nspden)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
 type(Cprj_type),intent(in) :: Cprj_bz(Cryst%natom,Dtset%nspinor*Sp%nbnds*Kmesh%nbz*Sp%nsppol*Dtset%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: level=26,enough=10,tim_fourdp=2
 character(len=50),parameter :: FILE__='csigme.F90'
 integer :: iab,counter,gwpara,iat,ib,ib1,ib2,ierr,ig,ig01,ig02,ig03,iggp,igp,ii,iik,itim_q,i2
 integer :: ikbz,ikibz,io,ioe0j,iop,isym_q,iq_bz,iq_ibz,is,isppol,istat,isym,jb,is_idx
 integer :: jik,jj,jkbz,jkibz,kb,master,maxbnd,rank,minbnd,nspinor
 integer :: nomega,nprocs,nq_summed,ispinor,ibsp,dimcprj_gw
 integer :: spad,spadx,spadc,spadx1,spadx2,spadc1,spadc2
 integer :: spaceComm,sumcxtab,wtqm,wtqp,mod10,nomega_sigc
 integer :: shift,indx_ki,base_kj,isym_kgw,isym_ki,natom
 real(dp) :: e0i,fact_sp,theta_mu_minus_e0i,tol_empty,z2,en_high,norm
 complex(dpc) :: ct,omegame0i2_ac,omegame0i_ac
 complex(gwpc) :: sigxme,ph_mkgwt,ph_mkt
 complex(dpc) :: scprod
 logical :: average_real=.TRUE.,cohsex,ltest
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: g0(3)
 integer :: spinor_padx(2,4),spinor_padc(2,4) 
 integer,allocatable :: cxtab(:,:,:),grottb(:,:,:),dimlmn(:)
 integer,pointer :: igfftg0(:),qtab(:),qtabi(:),qtabo(:)
 real(dp) :: ki(3),kj(3),kjmki(3),omegap(Er%nomega_i),omegap2(Er%nomega_i)
 real(dp) :: tsec(2)
 real(dp) :: z(Er%nomega_i),zw(Er%nomega_i)
 real(dp) :: gmet(3,3),gprimd(3,3)
 real(dp) :: spinrot_kbz(4),spinrot_kgw(4)
 real(dp),allocatable :: e0pde(:,:,:),omegame0i(:),otq(:,:)
 real(dp),allocatable :: otq_transp(:,:)
 complex(gwpc) :: sigcohme(Sp%nsig_ab)
 complex(dpc) :: ovlp(2)
 complex(gwpc),allocatable :: vc_sqrt_qbz(:),rhotwg(:),rhotwgp(:),sigsex(:)
 complex(gwpc),allocatable :: botsq(:,:),botsq_conjg_transp(:,:),dummy_cme(:,:,:,:)
 complex(gwpc),allocatable :: dummy_xme(:,:,:),eig(:,:),epsm1cqwz2(:,:,:)
 complex(gwpc),allocatable :: epsm1_qbz(:,:,:),epsm1q_trcc(:,:,:),epsm1_comp(:,:)
 complex(gwpc),allocatable :: integr(:,:,:),ket(:,:),ket1(:,:),ket2(:,:)
 complex(gwpc),allocatable :: kxcg(:,:)
 complex(gwpc),allocatable :: rhotwg_ki(:,:),sigccoh(:,:)
 complex(gwpc),allocatable :: sigcme2(:,:),sigcme_3(:),sigcme_ac(:,:)
 complex(gwpc),allocatable :: sigcme_new(:),sigcme_tmp(:,:,:,:),sigcsex(:,:)
 complex(gwpc),allocatable :: sigctmp(:,:),sigxme_tmp(:,:,:)
 complex(gwpc),allocatable :: wfr_jb(:,:,:),wfr_tmp(:)
 complex(gwpc),pointer :: cgup_jb(:),cgup_kb(:),cgup_tmp(:)
 complex(gwpc),pointer :: cgdwn_jb(:),cgdwn_kb(:),cgdwn_tmp(:)
 type(Cprj_type),allocatable :: Cprj_kj(:,:),Cprj_ki(:,:)
!
 complex(dpc),external :: ZDOTC
!************************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'CDOTC' :: cdotc
!DEC$ ATTRIBUTES ALIAS:'CGEMV' :: cgemv
#endif

#if defined DEBUG_MODE
 write(msg,'(a)')' cisgme : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif
 !
 ! === Initial check ===
 ltest=(Sr%nomega==Sp%nomegasr) 
 call assert(ltest,'Sr%nomega/=Sp%nomegasr',__FILE__,__LINE__)
 ltest=(Sp%npwvec==Gsph_Max%ng) 
 call assert(ltest,'mismatch in ng',__FILE__,__LINE__)
 ltest=(Sp%npwc==PPm%npwc)
 call assert(ltest,'mismatch in npwc',__FILE__,__LINE__)
 ltest=(Sp%npwc==Er%npwe)
 call assert(ltest,'Sp%npwc/=Er%npwe',__FILE__,__LINE__)
 if (Sp%nspinor==2) then
  ltest=(Sp%symsigma==0)
  call assert(ltest,'symsigma and nspinor=2 not implemented',__FILE__,__LINE__)
 end if
 ltest=(Wf_info%nfftot==nfftgw_tot)
 call assert(ltest,'nfftot/nfftgw_tot',__FILE__,__LINE__)
 !
 ! === Start clocks ===
 call timab(421,1,tsec) ! csigme (tot)
 call timab(422,1,tsec) ! csigme(1)
 !
 ! === Initialize MPI variables ===
 call xcomm_init  (MPI_enreg,spaceComm)
 call xme_init    (MPI_enreg,rank     )         
 call xmaster_init(MPI_enreg,master   ) 
 call xproc_max(nprocs,ierr)
 ! * Fake MPI type for sequential part
 call initmpi_seq(MPI_enreg_seq) ; MPI_enreg_seq%nproc_fft=1 ; MPI_enreg_seq%me_fft=0
 gwpara=MPI_enreg%gwpara

 ! === Copy some values ===
 nspinor = Dtset%nspinor
 natom   = Cryst%natom

 gmet(:,:)  =Gsph_Max%gmet(:,:) 
 gprimd(:,:)=Gsph_Max%gprimd(:,:)

 spinor_padx(:,:)=RESHAPE((/0,0,Sp%npwx,Sp%npwx,0,Sp%npwx,Sp%npwx,0/),(/2,4/))
 spinor_padc(:,:)=RESHAPE((/0,0,Sp%npwc,Sp%npwc,0,Sp%npwc,Sp%npwc,0/),(/2,4/))

 if (Dtset%gwcomp==1) then
  en_high=MAXVAL(en_qp(:,Sp%nbnds,:)) + Dtset%gwencomp
  write(msg,'(6a,e10.4,a)')ch10,&
&  ' Use an extrapolar approximation to accelerate convergence',ch10,&
&  ' with respect to the number of bands included',ch10,&
&  ' with extrapolar energy: ',en_high*Ha_eV,' [eV]'
  call wrtout(std_out,msg,'COLL')
 end if
 ! === Out-of-core solution for epsilon ===
 if (Er%mqmem==0) then 
  write(msg,'(3a)')ch10,&
&  ' csigme : reading ',&
&  ' a single slice from file (slower but less memory required)'
  call wrtout(std_out,msg,'COLL')
 end if 
 !
 ! min and Max band indeces for GW corrections (for this k-point)
 minbnd=Sp%minbnd(ikcalc) ; ib1=minbnd 
 maxbnd=Sp%maxbnd(ikcalc) ; ib2=maxbnd

 allocate(e0pde(Sp%nomegasrd,minbnd:maxbnd,Sp%nsppol)) ; e0pde(:,:,:)=zero
 allocate(rhotwg_ki(Sp%npwx*nspinor,minbnd:maxbnd)) ; rhotwg_ki(:,:)=czero_gw
 allocate(rhotwg   (Sp%npwx*nspinor))
 allocate(rhotwgp  (Sp%npwx*nspinor))
 allocate(vc_sqrt_qbz(Sp%npwx))
 !
 ! === Normalization of theta_mu_minus_e0i ===
 ! * If nsppol==2, occ_qp $\in [0,1]$
 SELECT CASE (Sp%nsppol)
 CASE (1) 
  fact_sp=half ; tol_empty=0.01  ! below this value the state is assumed empty
  if (Sp%nspinor==2) then
   fact_sp=one ; tol_empty=0.005  ! below this value the state is assumed empty
  end if
 CASE (2)
  fact_sp=one  ; tol_empty=0.005 ! to be consistent and obtain similar results if a metallic
 CASE DEFAULT                    ! spin unpolarized system is treated using nsppol==2
  call assert(.FALSE.,'Wrong nsppol',__FILE__,__LINE__)
 END SELECT
 !
 ! Data related to the q-point sampling
 qtab  => Qmesh%tab (1:Qmesh%nbz)
 qtabi => Qmesh%tabi(1:Qmesh%nbz)
 qtabo => Qmesh%tabo(1:Qmesh%nbz)

 ! Here I need two G spheres one for epsilon, the other one for vc, this part has to be cleaned
 allocate(grottb(Sp%npwvec,Gsph_Max%timrev,Gsph_Max%nsym))
 grottb(:,:,:)=Gsph_Max%rottb(:,:,:)
 !
 ! === Print type of calculation ===
 mod10=MOD(Sp%gwcalctyp,10)
 if (mod10==0) write(msg,'(a)')' standard GW with PPM'
 if (mod10==1) write(msg,'(a)')' standard GW without PPM (analytical continuation)'
 if (mod10==2) write(msg,'(a)')' standard GW without PPM (contour deformation)'
 if (mod10==5) write(msg,'(a)')' Hartree-Fock calculation'
 if (mod10==6) write(msg,'(a)')' Screened Exchange calculation'
 if (mod10==7) write(msg,'(a)')' COHSEX calculation'
 if (mod10==8) write(msg,'(a)')' model GW with PPM'
 if (mod10==9) write(msg,'(a)')' model GW without PPM'
 call wrtout(std_out,msg,'COLL')

 if (Sp%gwcalctyp<10) then
  write(msg,'(a)')' Perturbative Calculation'
  if (Sp%gwcalctyp==1) write(msg,'(a)')' Newton Raphson method '
 else if (Sp%gwcalctyp<20) then
  write(msg,'(a)')' Self-Consistent on Energies only'
 else
  write(msg,'(a)')' Self-Consistent on Energies and Wavefunctions'
 end if
 call wrtout(std_out,msg,'COLL')
 !
 ! === Set up logical flags for Sigma calculation ===
 if (mod10==1.and.Sp%gwcalctyp/=1) then 
  write(msg,'(a)')' not yet implemented'
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 cohsex=.FALSE. 
 if (MOD(Sp%splitsigc,2)==1) then 
  ! Decomposition of Sigma_c into Coulomb-hole and screened-exchange
  ! FIXME this part is broken, should be rewritten!
  cohsex=.TRUE. 
  allocate(sigccoh(Sp%nbnds,Sp%nsppol*Sp%nsig_ab)) ; sigccoh(:,:)=zero 
  allocate(sigcsex(Sp%nbnds,Sp%nsppol*Sp%nsig_ab)) ; sigcsex(:,:)=zero 
 end if
 !
 ! === Index of the GW point in the BZ array, its image in IBZ and time-reversal ===
 jkbz=Sp%kcalc(ikcalc) 
 call get_BZ_item(Kmesh,jkbz,kj,jkibz,isym_kgw,jik,ph_mkgwt)
 spinrot_kgw(:)=Cryst%spinrot(:,isym_kgw)
 !
 write(msg,'(2a,3f8.3,2a,2(i3,a))')ch10,&
& ' Calculating <nk|Sigma|nk> at k = ',kj(:),ch10,&
& ' bands n = from ',ib1,' to ',ib2,ch10
 call wrtout(std_out,msg,'COLL')
 !
 ! Diagonal elements of $\Sigma_x$ and $\Sigma_c(\omega)$
 Sr%sigxme(ib1:ib2,jkibz,:)  =zero
 Sr%sigcme(ib1:ib2,jkibz,:,:)=czero
 !
 ! === Load wavefunctions for GW corrections ===
 allocate(wfr_tmp(Wf_info%nfftot*nspinor))
 allocate(wfr_jb(nfftgw_tot*nspinor,ib1:ib2,Sp%nsppol))
 wfr_jb(:,:,:)=czero_gw
 do is=1,Sp%nsppol
  do jb=ib1,ib2
   call get_wfr(Wf_info_braket,MPI_enreg,jb,ikcalc,is,wfr_tmp)
   wfr_jb(:,jb,is)=wfr_tmp(:) 
  end do
 end do 
 !
 ! === Additional allocations for PAW ===
 if (Dtset%usepaw==1) then 
  allocate(dimlmn(natom)) 
  do iat=1,natom
   dimlmn(iat)=Pawtab(Cryst%typat(iat))%lmn_size
  end do
  allocate(Cprj_ki(natom,nspinor)) ; call cprj_alloc(Cprj_ki,0,dimlmn)
  ! * Load cprj for GW states, note the different shape
  dimcprj_gw=nspinor*(ib2-ib1+1)*Dtset%nsppol
  allocate(Cprj_kj(natom,ib1:ib1+dimcprj_gw-1))
  call cprj_alloc(Cprj_kj,0,dimlmn)
  ibsp=ib1-1
  do is=1,Dtset%nsppol
   shift  =nspinor*Sp%nbnds*Kmesh%nbz*(is-1)
   base_kj=nspinor*Sp%nbnds*(jkbz-1)+shift
   do ib=ib1,ib2
    do ispinor=1,nspinor
     ibsp=ibsp+1
     do iat=1,natom
      Cprj_kj(iat,ibsp)%cp(:,:)=Cprj_bz(iat,base_kj+ibsp)%cp(:,:)
     end do
    end do
   end do
  end do
 end if
 !
 ! === Setup frequencies around the KS\QP eigenvalues to compute Sigma derivatives (notice the spin) ===
 ! TODO it is better using an odd Sp%nomegasrd so that the KS\QP eigenvalue is in the middle
 ioe0j=Sp%nomegasrd/2+1
 do is=1,Sp%nsppol
  do jb=ib1,ib2
   do io=1,Sp%nomegasrd
    e0pde(io,jb,is)=Sr%egw(jb,jkibz,is)+Sp%deltae*(io-ioe0j)
    Sr%omegasrd(jb,jkibz,io,is)=e0pde(io,jb,is)
   end do
  end do
 end do
 !
 ! === Initialize quantities needed by AC ===
 if (mod10==1) then
  ! * Calculate Gauss-Legendre quadrature knots and weights
  call coeffs_gausslegint(zero,one,z,zw,Er%nomega_i)
  ! First frequencies are always real 
  do io=1,Er%nomega_i
   if (ABS(AIMAG(one*Er%omega(Er%nomega_r+io))-(one/z(io)-one)) > 0.0001) then 
    write(msg,'(6a)')ch10,&
&    ' csigme : ERROR - ',ch10,&
&    '  Frequencies in the SCR file are not compatible with the analytic continuation.',ch10,&
&    '  Verify the frequencies in the SCR file'
    call wrtout(std_out,msg,'COLL') 
    write(*,*)AIMAG(Er%omega(Er%nomega_r+io)),(one/z(io)-one) 
    call leave_new('COLL')
   end if 
  end do
  ! === To calculate \int_0^\infty domegap f(omegap), we calculate \int_0^1 dz f(1/z-1)/z^2 ===
  omegap(:)=one/z(:)-one 
  omegap2(:)=omegap(:)*omegap(:)
  allocate(epsm1cqwz2(Sp%npwc,Sp%npwc,Er%nomega_i),STAT=istat)
  allocate(integr(Sp%npwc,Sp%npwc,Sp%nomegasi),STAT=istat)
  if (istat/=0) call memerr(FILE__,'integr',Sp%npwc**2*Sp%nomegasi,'spc')
 end if ! of AC
 !
 ! Calculate total number of frequencies and allocate related arrays
 ! NOTE sigcme2 is used to accumulate the diagonal matrix elements over k-points and 
 ! GW bands, used only in case of ppmodel 3 and 4 (TODO save memory)
 nomega=Sp%nomegasr+Sp%nomegasrd
 allocate(sigcme2(nomega,ib1:ib2),sigcme_3(nomega))
 sigcme2(:,:)=czero_gw ; sigcme_3(:)=czero_gw
 !
 ! Allocate arrays used to accumulate the matrix elements of \Sigma_c over 
 ! k-points and bands. Note that for AC requires only the imaginary frequencies
 nomega_sigc=Sp%nomegasr+Sp%nomegasrd
 if (mod10==1) nomega_sigc=Sp%nomegasi
 allocate(sigctmp(nomega_sigc,Sp%nsig_ab))
 allocate(ket(Sp%npwc*nspinor,nomega_sigc))
 if (mod10==6.or.mod10==7) allocate(sigsex(Sp%npwc))

 allocate(sigxme_tmp(            ib1:ib2,ib1:ib2,Sp%nsppol*Sp%nsig_ab))
 allocate(sigcme_tmp(nomega_sigc,ib1:ib2,ib1:ib2,Sp%nsppol*Sp%nsig_ab))
 sigxme_tmp(:,:,:)  =czero_gw 
 sigcme_tmp(:,:,:,:)=czero_gw

 !FIXME This quantities are only used for model GW if I am not wrong
 allocate(ket1(Sp%npwc*nspinor,nomega))
 allocate(ket2(Sp%npwc*nspinor,nomega))
 allocate(omegame0i(nomega))
 !
 ! ********************************************************
 ! /*** On the symmetrization of Sigma matrix elements ***/
 ! ********************************************************
 !  If  Sk = k+G0 then  M_G(k, Sq)= e^{-i( Sq+G).t} M_{ S^-1(G}   (k,q)
 !  If -Sk = k+G0 then  M_G(k,-Sq)= e^{-i(-Sq+G).t} M_{-S^-1(G)}^*(k,q)
 !
 ! Notice the absence of G0 in the expression. Moreover, when we sum over the little group, it turns out 
 ! that there is a cancellation of the phase factor associated to the non-symmorphic operations due to a 
 ! similar term coming from the symmetrization of \epsilon^{-1}. Mind however that the nonsymmorphic phase
 ! has to be considered when epsilon^-1 is reconstructed starting from the q-points in the IBZ. 
 !
 ! In case of Sigma calculations, the unitary transformation relating wavefunctions
 ! at symmetric k-points should be taken into account during the symmetrization 
 ! of the oscillator matrix elements. In case of G_oW_o and GW_o calculations, however,
 ! it is possible to make an invariant by just including all the degenerate states and 
 ! averaging the final results over the degenerate subset. Here we divide the states 
 ! where the QP energies are required into complexes. Note however that this approach is not 
 ! based on group theory, and it might lead to spurious results in case of accidental degeneracies.
 !
 nq_summed=Kmesh%nbz
 if (Sp%symsigma/=0) then
  call print_little_group(Ltg_k,std_out,Dtset%prtvol,'COLL')
  nq_summed=SUM(Ltg_k%ibzq(:))
  ! 
  ! === Find number of complexes and number of bands in each complex ===
  allocate(cxtab(maxbnd-minbnd+1,Sp%nsppol,maxbnd-minbnd+1))
  cxtab(:,:,:)=0
  do isppol=1,Sp%nsppol
   do ib=1,maxbnd-minbnd+1
    do jb=1,maxbnd-minbnd+1
     ! The tolerance is a little bit arbitrary (0.001 eV)
     ! It could be reduced, in particular in case of nearly accidental degeneracies
     if (ABS(en_qp(jkibz,ib-1+minbnd,isppol)-en_qp(jkibz,jb-1+minbnd,isppol))<0.001/Ha_ev) then
      cxtab(ib,isppol,jb)=1
     end if
    end do
   end do
  end do 
  if (ANY(cxtab(:,:,:)/=0)) then 
   write(msg,'(2a,3f8.3,a)')ch10,&
&   ' found degenerate states at k-point = ',kj(:),ch10
   call wrtout(std_out,msg,'COLL')
   do isppol=1,Sp%nsppol
    do ib=1,maxbnd-minbnd+1
     do jb=ib+1,maxbnd-minbnd+1
      if (cxtab(ib,isppol,jb)==1) then 
       write(msg,'(a,i2,a,i4,a,i4)')' (spin ',isppol,')',ib-1+minbnd,' <====> ',jb-1+minbnd
       call wrtout(std_out,msg,'COLL')
       ! If two states do not belong to the same complex => matrix elements of v_xc differ
!       if (ABS(Sr%vxcme(ib,jkibz,isppol)-Sr%vxcme(jb,jkibz,isppol))>ABS(tol6*Sr%vxcme(jb,jkibz,isppol))) then 
!        write(msg,'(11a)')ch10,&
!&       ' csigme : COMMENT -',ch10,&
!&       ' It seems that an accidental degeneracy occurs at this k-point ',ch10,&
!&       ' In this case, using symsigma=1 might lead to spurious results since the algorithm ',ch10,&
!&       ' will treat all these states as degenerate, and won''t be able to remove the degeneracy. ',ch10,&
!&       ' In order to avoid this deficiency, please, run the calculation using symsigma=0',ch10
!       call wrtout(std_out,msg,'COLL') !; call wrtout(ab_out,msg,'COLL') !; call leave_new('COLL')
!       end if
      end if 
     end do 
    end do 
   end do 
  end if 
 end if !symsigma

 write(msg,'(2a,i6,a)')ch10,&
& ' calculation status ( ',nq_summed,' to be completed):'
 call wrtout(std_out,msg,'COLL') 

 if (PPm%model/=0) then
  ! * There might be zero-sized arrays and some compiler might crash
  ! FIXME NOTE here otq is defined as real but omegatw is complex, this has to be investigated
  allocate(botsq(PPm%npwc,PPm%dm2_botsq),STAT=istat)
  allocate(otq(PPm%npwc,PPm%dm2_otq),STAT=istat)
  allocate(eig(PPm%dm_eig,PPm%dm_eig),STAT=istat) 
 end if

 ! Here we have a problem in case of CD since epsm1q might be huge
 ! TODO if single q (ex molecule) dont allocate epsm1q, avoid waste of memory 
 if (mod10==1.or.mod10==2.or.mod10==6.or.mod10==7.or.mod10==9) then
  allocate(epsm1_qbz(Sp%npwc,Sp%npwc,Er%nomega),STAT=istat) 
  if (istat/=0) call memerr(FILE__,'epsm1_qbz',Sp%npwc*2*Er%nomega,'spc')
 end if
 if (mod10==9) then ! MG should check whether I need also epsm1q
  allocate(epsm1q_trcc(Sp%npwc,Sp%npwc,Er%nomega),STAT=istat)
  if (istat/=0) call memerr(FILE__,'epsm1q_trcc',Sp%npwc*2*Er%nomega,'spc')
 end if 
 call timab(422,2,tsec)

 call timab(499,1,tsec) 
 counter=0 

 ! === Loop over k_i in BZ ===
 do ikbz=1,Kmesh%nbz
  call timab(423,1,tsec) ! csigme (initq)
  ! 
  ! === Parallelization over k-points and spin ===
  ! * For the spin there is another check in the inner loop
  if (gwpara==1) then
   if (MINVAL(ABS(MPI_enreg%proc_distrb(ikbz,:,:)-rank))/=0) CYCLE
   if (ANY(MPI_enreg%proc_distrb(ikbz,:,:)==-999)) stop '-999'
  end if
  ! 
  ! * Find irreducible k-point
  call get_BZ_item(Kmesh,ikbz,ki,ikibz,isym_ki,iik,ph_mkt)
  spinrot_kbz(:)=Cryst%spinrot(:,isym_ki)

  ! * Identify q and G0 where q+G0=k_GW-k_i
  kjmki(:)=kj(:)-ki(:)
  call findqg0(iq_bz,g0,kjmki,Qmesh%nbz,Qmesh%bz,Sp%mG0) 

  wtqp=1 ; wtqm=0
  if (Sp%symsigma/=0) then
   ! * Sum only q"s in IBZ_k. In this case matrix elements are weighted 
   !   according to wtqp and wtqm. wtqm is for time-reversal. 
   if (Ltg_k%ibzq(iq_bz)/=1) CYCLE 
   wtqp=0 ; wtqm=0
   do isym=1,Ltg_k%nsym_sg
    wtqp=wtqp+Ltg_k%wtksym(1,isym,iq_bz)
    wtqm=wtqm+Ltg_k%wtksym(2,isym,iq_bz)
   end do
  end if

  counter=counter+1
  call status(ikbz,Dtfil%filstat,0,level,'loop ikpt     ')
  !if (counter<=enough) then 
   write(msg,'(2(a,i4)a,i3)')' csigme : ik = ',ikbz,' / ',Kmesh%nbz,' done by processor ',rank
   call wrtout(std_out,msg,'PERS')
  !end if
  ! 
  ! === Find irred q-point, and define the G_o shift for the FFT  ===
  ! Here there is a problem with the small q
  iq_ibz=qtab(iq_bz) ; isym_q=qtabo(iq_bz) ; itim_q=(3-qtabi(iq_bz))/2
  ig01=g0(1)+Sp%mG0(1)+1
  ig02=g0(2)+Sp%mG0(2)+1
  ig03=g0(3)+Sp%mG0(3)+1
  igfftg0 => igfft(:,ig01,ig02,ig03)

  if (Er%mqmem==0) then  
   ! === Read a q-slice of epsilon^{-1}|chi0 in Er%epsm1(:,:,:,1) (much slower but less memory) ===
   allocate(kxcg(1,1)) !this should be input
   call get_epsm1(Er,Vcp,Dtfil,0,0,Dtset%accesswff,Dtset%localrdwf,kxcg,gmet,MPI_enreg,iqibzA=iq_ibz)
   deallocate(kxcg)
   call setup_ppmodel(PPm,Dtset%paral_kgb,Qmesh,Sp,Er,MPI_enreg,nfftf,gvec,ngfftf,&
&   gmet,gprimd,rhor(:,1),Er%epsm1(:,:,:,1),iqiA=iq_ibz)
  end if 

   ! === Symmetrize PPM parameters and epsm1 (q_IBZ --> q_BZ) ===
   ! * NOTE: We are not considering umklapp with G0/=0. In this case, 
   !   indeed the equation is different since we have to use G-G0. 
   !   A check, however, is performed in sigma.
   ! * If gwcomp==1 and mod10=1,2,9, one needs both to set up botsq and epsm1_q

  if (mod10==0.or.mod10==8.or.Dtset%gwcomp==1) then
   ! FIXME still I have to solve possible issue with npwc!
   call PPmodel_symmetrize(PPm,Gsph_Max,Qmesh,iq_bz,botsq,otq,eig) 
  end if

  if (mod10==1.or.mod10==2.or.mod10==9) then
   ! === Numerical integration or model GW with contour deformation or Analytic Continuation ===
   !  TODO In case of AC we should symmetrize only the imaginary frequencies
    call Epsm1_symmetrizer(iq_bz,Er%nomega,Sp%npwc,Er,Gsph_Max,Qmesh,.TRUE.,epsm1_qbz) 

   if (mod10==1) then 
    ! === For AC prepare the first term: \Sum w_i 1/z_i^2 f(1/z_i - 1) ===
    ! * Note that the first frequencies are always real, skip them.
    ! * Obviously memory is not optimized.
    do iop=1,Er%nomega_i
     z2=z(iop)*z(iop)
     epsm1cqwz2(:,:,iop)= zw(iop)*epsm1_qbz(:,:,Er%nomega_r+iop)/z2
    end do
   end if 

   if (mod10==9) then
    ! === For model GW we need transpose(conjg(epsm1_qbz)) ===
    do io=1,Er%nomega
     epsm1q_trcc(:,:,io)=TRANSPOSE(CONJG(epsm1_qbz(:,:,io)))
    end do
   end if

  else if (mod10==6.or.mod10==7) then
   ! === SEX or COHSEX ===
   ! * Only omega==0 is needed 
   call Epsm1_symmetrizer(iq_bz,1,Sp%npwc,Er,Gsph_Max,Qmesh,.TRUE.,epsm1_qbz) 
  end if !gwcalctyp
  ! 
  ! === Get Fourier components of the Coulombian interaction in the BZ ===
  ! * In 3D systems, neglecting umklapp,  vc(Sq,sG)=vc(q,G)=4pi/|q+G|
  ! * The same relation holds for 0-D systems, but not in 1-D or 2D systems. It depends on S.
  do ig=1,Sp%npwx
   vc_sqrt_qbz(grottb(ig,itim_q,isym_q))=Vcp%vc_sqrt(ig,iq_ibz)
  end do
  call timab(423,2,tsec) ! csigme (initq)
  ! 
  ! === Loop over spin and bands ===
  ! * Bands are summed up, spin is external parameter.
  do isppol=1,Sp%nsppol
   do ib=1,Sp%nbnds
    !   
    ! === Parallelism over spin ===
    ! * This processor has this k-point but what about spin?
    if (gwpara==1) then
     if (MPI_enreg%proc_distrb(ikbz,ib,isppol)/=rank) CYCLE
     if (ANY(MPI_enreg%proc_distrb(ikbz,:,:)==-999)) STOP '-999'
    end if
    if (gwpara==2) then ! Note ikibz
     if (MPI_enreg%proc_distrb(ikibz,ib,isppol)/=rank) CYCLE
    end if
    !
    ! === Skip empty state ib for HF, SEX, and COHSEX ===
    if (occ_qp(ikibz,ib,isppol)<tol_empty.and.(mod10==5.or.mod10==6.or.mod10==7)) CYCLE

    call get_wfr(Wf_info,MPI_enreg,ib,ikibz,isppol,wfr_tmp)

    if (Dtset%usepaw==1) then 
     ! === Load cprj for point ki, this spin or spinor and *THIS* band ===
     ! ??? MG I could avoid doing this but I have to exchange spin and bands ???
     ! For sure there is a better way to do this!
     shift  =nspinor*Sp%nbnds*Kmesh%nbz*(isppol-1)
     indx_ki=nspinor*Sp%nbnds*(ikbz-1)+shift
     do ispinor=1,nspinor
      ibsp=ib+(ispinor-1)
      do iat=1,natom
       Cprj_ki(iat,ispinor)%cp(:,:)=Cprj_bz(iat,indx_ki+ibsp)%cp(:,:)
      end do
     end do
    end if

    if (mod10==1) then
     integr(:,:,:)=czero_gw
     do io=1,Sp%nomegasi
      ! === Calculate integral over omegap with Gauss-Legendre quadrature ===
      ! * -1/pi \int domegap epsm1c*(omega-e0i) / ( (omega-e0i)^2 + omegap^2)
      ! * Note that energies are calculated wrt the Fermi level.
      omegame0i_ac  = Sp%omegasi(io)-en_qp(ikibz,ib,isppol)
      omegame0i2_ac = omegame0i_ac*omegame0i_ac
      do iop=1,Er%nomega_i
       do iggp=0,Sp%npwc*Sp%npwc-1
        ig=iggp/Sp%npwc+1
        igp= iggp-(ig-1)*Sp%npwc+1
        ! \int domegap epsm1c/((omega-e0i)^2 + omegap^2)
        integr(ig,igp,io)= integr(ig,igp,io) + epsm1cqwz2(ig,igp,iop)/(omegame0i2_ac + omegap2(iop))
       end do
      end do
      integr(:,:,io)=integr(:,:,io)*omegame0i_ac
     end do
     integr(:,:,:)=-integr(:,:,:)*piinv
    end if

    do jb=ib1,ib2
     !    
     ! === Get all <k-q,ib,s|e^{-i(q+G).r}|s,jb,k>, at once ===
     call rho_tw_g(Dtset%paral_kgb,nspinor,Sp%npwx,nfftgw_tot,ngfft,1,igfftg0,&
&     wfr_tmp(:)         ,iik,ktabr(:,ikbz),ph_mkt  ,spinrot_kbz,  & 
&     wfr_jb(:,jb,isppol),jik,ktabr(:,jkbz),ph_mkgwt,spinrot_kgw,& 
&     nspinor,rhotwg_ki(:,jb),tim_fourdp,MPI_enreg)

     if (Dtset%usepaw==1) then 
      ! === Add on-site contribution, projectors are already in BZ ===
      !TODO Recheck this!
! FIXME Find a clever way to deal with spinors
      shift=nspinor*(ib2-ib1+1)*(isppol-1) 
      i2=jb+shift ; if (nspinor==2) i2=(2*jb-1)
      spad=(nspinor-1)
      call paw_rho_tw_g(Sp%npwx,nspinor,nspinor,natom,Psps%lmnmax,dimlmn,&
&      Cprj_ki(:,:),Cprj_kj(:,i2:i2+spad),pawrhox(:,:,:,:,iq_bz),rhotwg_ki(:,jb))
     end if
     !
     ! === Multiply by the square root of the coulombian ===
     ! * In 3-D systems, the factor sqrt(4pi) is included)
     do ii=1,nspinor
      spad=(ii-1)*Sp%npwx
      rhotwg_ki(spad+1:spad+Sp%npwx,jb)=rhotwg_ki(spad+1:spad+Sp%npwx,jb)*vc_sqrt_qbz(1:Sp%npwx)
     end do
     !    
     ! === Treat analytically the case q --> 0 ===
     ! * The oscillator is evaluated at q=O as it is considered constant in the small cube around Gamma
     !   while the colulombian term is integrated.
     ! * In the scalar case we have nonzero contribution only if ib==jb
     ! * For nspinor==2 evalute <ib,up|jb,up> and <ib,dwn|jb,dwn>, 
     !   impose orthonormalization since npwwfn might be < npwvec.
     if (ikbz==jkbz) then
      if (nspinor==1) then 
       rhotwg_ki(1,jb)=czero_gw
       if (ib==jb) rhotwg_ki(1,jb)=CMPLX(SQRT(Vcp%i_sz),0.0_gwp)
      else
       cgup_tmp  => Wf_info%wfg(1:Wf_info%npwwfn, ib,ikibz,isppol)
       cgdwn_tmp => Wf_info%wfg(Wf_info%npwwfn+1:2*Wf_info%npwwfn,ib,ikibz,isppol)
       cgup_jb   => Wf_info_braket%wfg(1:Wf_info_braket%npwwfn,jb,ikcalc,isppol)
       cgdwn_jb  => Wf_info_braket%wfg(Wf_info_braket%npwwfn+1:2*Wf_info_braket%npwwfn,jb,ikcalc,isppol)
! FIXME Find a clever way to deal with spinors
      !TODO Recheck this!
       !shift=nspinor*(ib2-ib1+1)*(isppol-1) 
       !i2=jb+shift ; if (nspinor==2) 
       i2=(2*jb-1)
       ovlp(1) = overlap_cmplx( cgup_tmp, cgup_jb,Dtset%usepaw,Cprj_ki(:,1),Cprj_kj(:,i2  ),Cryst%typat,Pawtab)
       ovlp(2) = overlap_cmplx(cgdwn_tmp,cgdwn_jb,Dtset%usepaw,Cprj_ki(:,2),Cprj_kj(:,i2+1),Cryst%typat,Pawtab)
       !ovlp(2) = -ovlp(1)
       !if (ib==jb) ovlp(2)=cone_gw-ovlp(1)
       if (ib==jb) then 
        norm=REAL(ovlp(1)+ovlp(2))
        ovlp(1)=REAL(ovlp(1)/norm)
        ovlp(2)=REAL(ovlp(2)/norm)
       else 
        scprod=ovlp(1)+ovlp(2)
        ovlp(1)=ovlp(1)-scprod*half
        ovlp(2)=ovlp(2)-scprod*half
       end if
       rhotwg_ki(1        ,jb)=CMPLX(SQRT(Vcp%i_sz),0.0_gwp)*ovlp(1)
       rhotwg_ki(Sp%npwx+1,jb)=CMPLX(SQRT(Vcp%i_sz),0.0_gwp)*ovlp(2)
      end if
     end if
     ! Got all matrix elements from minbnd up to maxbnd
    end do !jb 

    theta_mu_minus_e0i=fact_sp*occ_qp(ikibz,ib,isppol)

    ! Starting point to evaluate the derivative of Sigma and the Spectral function
    e0i=en_qp(ikibz,ib,isppol)

    ! Frequencies for the spectral function, e0i=en_qp(ikibz,ib,isppol)
    ! FIXME the interval is not centered on eoi ! WHY?
    if (Sp%nomegasr>0) omegame0i(1:Sr%nomega)=REAL(Sp%omegasf(1:Sr%nomega))-e0i

    do kb=ib1,ib2 
    
     call timab(424,1,tsec) ! csigme (SigX)  FIXME this is wrong
     !    
     ! Get frequencies $\omega$-\epsilon_in$ to evaluate  $d\Sigma/dE$, note the spin
     ! subtract e_KS since we have stored e_KS+ Delta \omega in Sr%omegasrd, not required for AC 
     do io=Sp%nomegasr+1,nomega
      omegame0i(io)=REAL(Sr%omegasrd(kb,jkibz,io-Sp%nomegasr,isppol))-e0i
     end do 
     !    
     ! === Get the ket \Sigma|\phi_{k,kb}> according to the method ===
     rhotwgp(:)=rhotwg_ki(:,kb)
     ket (:,:)=czero_gw 
     ket1(:,:)=czero_gw
     ket2(:,:)=czero_gw

     SELECT CASE (mod10)

     CASE (0)
      ! === GW WITH Plasmon-Pole Model ===
      ! * Note that ppmodel 3 or 4 work only in case of standard perturbative approach!
      !   Moreover, for ppmodel 3 and 4, spinorial case is not allowed
      call calc_sig_ppm(PPm,nspinor,Sp%npwc,nomega,rhotwgp,botsq,otq,&
&      omegame0i,Sp%zcut,theta_mu_minus_e0i,eig,Sp%npwx,ket,sigcme_3)

      if (PPm%model==3.or.PPm%model==4) then
       sigcme2(:,kb)=sigcme2(:,kb) + (wtqp+wtqm)*REAL(sigcme_3(:)) + (wtqp-wtqm)*j_gw*AIMAG(sigcme_3(:))
      end if

     CASE (1)
      ! === GW with Analytical continuation ===
      ! * Evaluate \sum_Gp integr_GGp(omegasi) rhotw_Gp === 
      ! TODO this part can be optimized
      do io=1,Sp%nomegasi
       do ispinor=1,Sp%nspinor
        spadx=(ispinor-1)*Sp%npwx
        spadc=(ispinor-1)*Sp%npwc

        do ig=1,Sp%npwc
         ct=czero
         do igp=1,Sp%npwc
          ct=ct+integr(ig,igp,io)*rhotwgp(igp+spadx)
         end do
         ket(ig+spadc,io)=ct 
        end do !ig

       end do !ispinor
      end do !io

     CASE (2)
      ! === GW calculation with contour deformation ===
      call calc_sig_noppm(Sp%npwc,Sp%npwx,Sp%nspinor,nomega,Er%nomega,Er%nomega_r,Er%nomega_i,rhotwgp,&
&      Er%omega,epsm1_qbz,omegame0i,theta_mu_minus_e0i,ket)

     CASE (6,7)
      ! === SEX or COHSEX, here only SEX part ===
      ! TODO add check on theta_mu_minus_e0i
      do ispinor=1,nspinor
       spadx=(ispinor-1)*Sp%npwx
       spadc=(ispinor-1)*Sp%npwc
#if defined HAVE_GW_DPC
       call ZGEMV('N',Sp%npwc,Sp%npwc,cone_gw,epsm1_qbz(:,:,1),Sp%npwc,rhotwgp(1+spadx),1,czero_gw,sigsex,1)
#else
       call CGEMV('N',Sp%npwc,Sp%npwc,cone_gw,epsm1_qbz(:,:,1),Sp%npwc,rhotwgp(1+spadx),1,czero_gw,sigsex,1)
#endif
       sigsex(:)= -theta_mu_minus_e0i*sigsex(:)
       ! here nomega==1 as SEX is energy independent.
       do io=1,nomega 
        ket(spadc+1:spadc+Sp%npwc,io)=sigsex(:)
       end do
      end do

     CASE (8)
      ! === MODEL GW calculation WITH PPm ===
      ! TODO Spinor not tested.
      allocate(sigcme_new(nomega))
      ! * Calculate \Sigma(E_k) |k> to obtain <j|\Sigma(E_k)|k>
      call calc_sig_ppm(PPm,nspinor,Sp%npwc,nomega,rhotwgp(1:Sp%npwc),botsq,otq,&
&      omegame0i,Sp%zcut,theta_mu_minus_e0i,eig,Sp%npwx,ket1,sigcme_new)

      if (Sp%gwcalctyp==28) then
       if (PPm%model/=1.and.PPm%model/=2) then ! This is needed to have Sp%npwc=PPm%dm2_botsq=PPm%dm2_otq
        write(msg,'(6a)')ch10,&
&        ' csigme : ERROR -',ch10,&
&        ' For the time being, gwcalctyp=28 cannot be used with ppmodel=3 or ppmodel=4.',ch10,&
&        ' Use another Plasmon Pole Model when gwcalctyp=28.'
        call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
       end if
       allocate(botsq_conjg_transp(PPm%dm2_botsq,Sp%npwc))
       botsq_conjg_transp=TRANSPOSE(botsq) ! Keep these two lines separated, otherwise gfortran messes up
       botsq_conjg_transp=CONJG(botsq_conjg_transp)
       allocate(otq_transp(PPm%dm2_otq,PPm%npwc),STAT=istat)
       otq_transp=TRANSPOSE(otq)

       call calc_sig_ppm(PPm,nspinor,Sp%npwc,nomega,rhotwgp(1:Sp%npwc),botsq_conjg_transp,otq_transp,&
&       omegame0i,Sp%zcut,theta_mu_minus_e0i,eig,Sp%npwx,ket2,sigcme_3)

       deallocate(botsq_conjg_transp)
       deallocate(otq_transp)
       ket(:,:)=(ket1(:,:)+ket2(:,:))*0.5
      else
       ket(:,:)=ket1(:,:)
      end if
      deallocate(sigcme_new)

     CASE (9)
      ! === MODEL GW calculation numerical integration ===
      ! * Calculate \Sigma(E_k)|k> to obtain <j|\Sigma(E_k)|k>
      call calc_sig_noppm(Sp%npwc,Sp%npwx,Sp%nspinor,nomega,Er%nomega,Er%nomega_r,Er%nomega_i,rhotwgp,&
&      Er%omega,epsm1_qbz,omegame0i,theta_mu_minus_e0i,ket1)
      if (Sp%gwcalctyp==29) then
       ! * Calculate <k|\Sigma(E_k) to obtain <k|\Sigma(E_k)|j>^*
       call calc_sig_noppm(Sp%npwc,Sp%npwx,Sp%nspinor,nomega,Er%nomega,Er%nomega_r,Er%nomega_i,rhotwgp,&
&       Er%omega,epsm1q_trcc,omegame0i,theta_mu_minus_e0i,ket2)
       ket(:,:)= (ket1(:,:)+ket2(:,:))*0.5
      else
       ket(:,:)=ket1(:,:)
      end if

     CASE DEFAULT 
      continue
     END SELECT ! Sp%gwcalctyp

     if (Dtset%gwcomp==1) then
      ! TODO spinor not implemented
      call calc_sig_ppm_comp(Sp%npwc,nomega,rhotwgp,botsq,otq,dble(Sr%egw(kb,jkibz,isppol)-en_high),Sp%zcut,&
&       theta_mu_minus_e0i,ket,PPm%model,Sp%npwx,PPm%dm2_botsq,PPm%dm2_otq)
     end if

     ! === Loop on the Left band index ===
     ! TODO Use hermitianity for HF, SEX, COHSEX to reduce the loop to: do jb=minbnd,kb
     do jb=ib1,ib2  

      ! * If not self-consistent, only diagonal elements
      if (Sp%gwcalctyp<20 .and. jb/=kb) CYCLE
      rhotwg(:)=rhotwg_ki(:,jb)

      ! === Calculate bare exchange <\phi_j|\Sigma_x|\phi_k> ===
      ! * Do the scalar product only if ib is occupied.

      ! TODO this is for 5.6.2 some of the sc-gw tests are slightly affected.
      !if (theta_mu_minus_e0i/fact_sp >= tol_empty) then 
       do iab=1,Sp%nsig_ab
        spadx1=spinor_padx(1,iab)
        spadx2=spinor_padx(2,iab)
#if defined HAVE_GW_DPC
        sigxme=-ZDOTC(Sp%npwx,rhotwg(spadx1+1),1,rhotwgp(spadx2+1),1)*theta_mu_minus_e0i
#else
        sigxme=-CDOTC(Sp%npwx,rhotwg(spadx1+1),1,rhotwgp(spadx2+1),1)*theta_mu_minus_e0i
#endif
        ! === Accumulate and, in case, symmetrize Sigma_x ===
        ! * -wtqm comes from time-reversal (exchange of band indeces)
        is_idx=isppol ; if (nspinor==2) is_idx=iab
        sigxme_tmp(jb,kb,is_idx) = sigxme_tmp(jb,kb,is_idx) + &
&        (wtqp+wtqm)*REAL(sigxme) + (wtqp-wtqm)*j_gw*AIMAG(sigxme)
       end do
      !end if
      call timab(424,2,tsec) ! csigme (SigX) !FIXME this is wrong

      ! === Calculate <\phi_j|\Sigma_c|\phi_k> ===
      ! * Different freqs according to method (AC or Perturbative), see nomega_sigc.
      call timab(425,1,tsec) ! csigme (SigC)

      do iab=1,Sp%nsig_ab
       spadx1=spinor_padx(1,iab)
       spadc2=spinor_padc(2,iab)
       do io=1,nomega_sigc
#if defined HAVE_GW_DPC
        sigctmp(io,iab)=ZDOTC(Sp%npwc,rhotwg(spadx1+1),1,ket(spadc2+1,io),1)
#else
        sigctmp(io,iab)=CDOTC(Sp%npwc,rhotwg(spadx1+1),1,ket(spadc2+1,io),1)
#endif
       end do
      end do

      if (ib==1) then
       ! === Evaluate Static COH ===
       ! * Calculated once as it does not depend on the index ib summed over.
       ! FIXME PAW not yet implemented, see calc_wfwfg 
       ! mapping onto FFT box in paw_rho_tw_g is needed.
       if (mod10==7) then
        call calc_coh(Dtset%paral_kgb,Sp%nspinor,Sp%nsig_ab,nfftgw_tot,ngfft,tim_fourdp,MPI_enreg,&
&        jik,ktabr(:,jkibz),spinrot_kgw,wfr_jb(:,jb,isppol),wfr_jb(:,kb,isppol),Sp%npwx,Sp%npwc,gvec,&
&        epsm1_qbz(:,:,1),vc_sqrt_qbz,Vcp%i_sz,iq_ibz,(jb==kb),sigcohme)
        do io=1,nomega_sigc ! Should be 1
         sigctmp(io,:) = sigctmp(io,:)+sigcohme(:)
        end do
        !write(std_out,'(a,(es14.6))')' === sigcohme ',REAL(SUM(sigcohme(:)))/(Cryst%ucvol*Kmesh%nbz)*Ha_eV
        !&' === sigcohme ',(REAL(sigcohme(iab))/(Cryst%ucvol*Kmesh%nbz)*Ha_eV,iab=1,MIN(2,Sp%nsig_ab))
       end if

       ! === The static contribution from completeness relation is calculated once ===
       ! FIXME spinor should be tested
       if (Dtset%gwcomp==1) then
        allocate(epsm1_comp(Sp%npwc,Sp%npwc))
        epsm1_comp(:,:) = botsq(:,:) / ( otq(:,:) * ( (Sr%egw(kb,jkibz,isppol)-en_high)-otq(:,:) ) )
        call calc_coh(Dtset%paral_kgb,Sp%nspinor,Sp%nsig_ab,nfftgw_tot,ngfft,tim_fourdp,MPI_enreg,&
&        jik,ktabr(:,jkibz),spinrot_kgw,wfr_jb(:,jb,isppol),wfr_jb(:,kb,isppol),Sp%npwx,Sp%npwc,gvec,&
&        epsm1_comp,vc_sqrt_qbz,Vcp%i_sz,iq_ibz,(jb==kb),sigcohme)
        deallocate(epsm1_comp)
        do io=1,nomega_sigc
         sigctmp(io,:) = sigctmp(io,:)+sigcohme(:)
        end do
       end if
      end if
      !
      ! === Accumulate and, in case, symmetrize matrix elements of Sigma_c ===
      do iab=1,Sp%nsig_ab
       is_idx=isppol ; if (nspinor==2) is_idx=iab
       sigcme_tmp(:,jb,kb,is_idx)=sigcme_tmp(:,jb,kb,is_idx) + &
&       (wtqp+wtqm)*REAL(sigctmp(:,iab)) + (wtqp-wtqm)*j_gw*AIMAG(sigctmp(:,iab))
      end do

      ! Decomposition Sigma_c into Coulomb-hole and screened-exchange (GMR)
      ! TODO Commented as this routine is broken should be reinstated
      !if (cohsex) call split_sigc(Sp,Sr,PPm%dm2_botsq,PPm%dm2_otq,jb,isppol,io,ioe0j,theta_mu_minus_e0i,&
      !& rhotwg,omegame0i,otq,botsq,sigccoh,sigcsex) ! MG WARNING io now is meaningless 

      call timab(425,2,tsec) ! csigme(SigC)
     end do !jb used to calculate matrix elements of $\Sigma$

     ! shaltaf (030406): this has to be done in a clean way later
     ! FIXME doesn work with spinor
     if (mod10==0.and.(PPm%model==3.or.PPm%model==4)) then
       sigcme_tmp(:,kb,kb,isppol)= sigcme2(:,kb)
     end if

    end do !kb to calculate matrix elements of $\Sigma$
   end do !ib
  end do !end loop over spin
  !££if (counter=enough) then 
   !££call timab(499,4,tsec)
   !££write(msh,'(a)')' estimated CPU  time [s] : ',my_nkpt*tsec(2)/enough
   !££write(msh,'(a)')' estimated WALL time [s] : ',my_nkpt*tsec(1)/enough 
  !££end if
 end do !ikbz
 !
 ! === Got all diagonal (off-diagonal) matrix elements for this k ===
 write(msg,'(a)')' csigme : finished the k point loop'
 call wrtout(std_out,msg,'COLL')
 deallocate(sigcme2,sigcme_3)
 !
 ! === Master sums up contributions from all the CPUs ===
 call xsum_master(sigcme_tmp(:,:,:,:),master,spaceComm,ierr)
 call xsum_master(sigxme_tmp(:,:,:)  ,master,spaceComm,ierr)
 if (rank/=master) RETURN
 !
 ! === Multiply by constants ===
 ! * For 3D systems sqrt(4pi) is included in vc_sqrt_qbz ===
 sigxme_tmp(:,:,:)  = (one/(Cryst%ucvol*Kmesh%nbz))*sigxme_tmp(:,:,:)
 sigcme_tmp(:,:,:,:)= (one/(Cryst%ucvol*Kmesh%nbz))*sigcme_tmp(:,:,:,:)

 call timab(426,1,tsec) ! b-loop
 !
 ! === If we have summed over the IBZ_q now we have to average over complexes ===
 ! * Presently only diagonal terms are considered
 ! * TODO QP-SCGW required a more involved approach, there is a check in sigma
 ! FIXME doesn work with spinor
 if (Sp%symsigma/=0) then
  allocate(dummy_xme(-minbnd+maxbnd+1,-minbnd+maxbnd+1,Sp%nsppol))
  if (mod10==1) then 
   ! FIXME here there is a problem in case of AC with symmetries
   allocate(dummy_cme(Sp%nomegasi,-minbnd+maxbnd+1,-minbnd+maxbnd+1,Sp%nsppol))
  else 
   allocate(dummy_cme(nomega,-minbnd+maxbnd+1,-minbnd+maxbnd+1,Sp%nsppol))
  end if 
  dummy_cme=czero_gw ; dummy_xme=czero_gw
  ! === Average over degenerate diagonal elements ===
  do is=1,Sp%nsppol
   do ib=1,maxbnd-minbnd+1
    sumcxtab=0
    do jb=1,maxbnd-minbnd+1
     if (cxtab(ib,is,jb)==1) then
      dummy_xme(ib,ib,is)=dummy_xme(ib,ib,is)+sigxme_tmp(minbnd-1+jb,minbnd-1+jb,is)
      ! NOTE: frequencies should be equal, another good reason,
      ! to use a strict criterion for the tollerance on eigenvalues.
      dummy_cme(:,ib,ib,is)=dummy_cme(:,ib,ib,is)+sigcme_tmp(:,minbnd-1+jb,minbnd-1+jb,is)
     end if
     sumcxtab=sumcxtab+cxtab(ib,is,jb)
    end do
    dummy_xme  (ib,ib,is)=dummy_xme  (ib,ib,is)/sumcxtab
    dummy_cme(:,ib,ib,is)=dummy_cme(:,ib,ib,is)/sumcxtab
   end do
  end do

  ! === Copy values ===
  do is=1,Sp%nsppol
   do ib=1,maxbnd-minbnd+1
    do jb=1,maxbnd-minbnd+1
     ! if (ib/=jb) CYCLE
     sigxme_tmp(ib-1+minbnd,jb-1+minbnd,is)=dummy_xme(ib,jb,is)
     ! this is to check another scheme in case of AC
     if (mod10==1.and.average_real) CYCLE
     sigcme_tmp(:,ib-1+minbnd,jb-1+minbnd,is)=dummy_cme(:,ib,jb,is)
    end do
   end do
  end do
 end if

 ! =====================================================
 ! ==== Solve Dyson equation storing results in Sr% ====
 ! =====================================================
 ! * Use perturbative approach or AC to find QP corrections.
 ! * If qp-GW, diagonalize also H0+Sigma in the KS basis set to get the 
 !   new QP amplitudes and energies (Sr%eigvec_qp and Sr%en_qp_diago.
 ! TODO AC with spinor not implemented yet. 
 ! TODO Diagonalization of Sigma+hhartre with AC is wrong.

 call solve_Dyson(ikcalc,nomega_sigc,Sp,Kmesh,sigxme_tmp,sigcme_tmp,en_qp,Sr)

 if (cohsex) then
  ! 4pi is included in vc_sqrt_qbz
  sigccoh(:,:)=(one/(Cryst%ucvol*Kmesh%nbz))*sigccoh(:,:)
  sigcsex(:,:)=(one/(Cryst%ucvol*Kmesh%nbz))*sigcsex(:,:)
  write(std_out,'(/a)')' COH and SEX decomposition of Sig_c(E0) (eV):'
  write(std_out,'(a5,3a10)')' band','     Sig_c','   Sig_coh','   Sig_sex'
  do is=1,Sp%nsppol
   do jb=minbnd,maxbnd
    write(std_out,'(i5,3f10.5," (Re)")')jb,REAL(Sr%sigcmee0(jb,jkibz,is))*Ha_eV,&
&    REAL(sigccoh(jb,is))*Ha_eV,REAL(sigcsex(jb,is))*Ha_eV
    write(std_out,'(5x,3f10.5," (Im)")')   AIMAG(Sr%sigcmee0(jb,jkibz,is))*Ha_eV,&
&    AIMAG(sigccoh(jb,is))*Ha_eV,AIMAG(sigcsex(jb,is))*Ha_eV
   end do
  end do
 end if
 !
 ! === Deallocate memory ===
 if (Dtset%usepaw==1) then
  deallocate(dimlmn)
  call cprj_free(Cprj_ki) ; deallocate(Cprj_ki)
  call cprj_free(Cprj_kj) ; deallocate(Cprj_kj)
 end if

 deallocate(wfr_jb,wfr_tmp)
 deallocate(ket,ket1,ket2)
 deallocate(rhotwg_ki,rhotwg,rhotwgp,vc_sqrt_qbz)
 deallocate(e0pde,omegame0i,sigcme_tmp,sigctmp)
 deallocate(grottb)

 if (mod10==6.or.mod10==7) deallocate(sigsex)

 if (allocated(sigccoh    ))  deallocate(sigccoh)
 if (allocated(sigcsex    ))  deallocate(sigcsex)
 if (allocated(epsm1_qbz  ))  deallocate(epsm1_qbz)
 if (allocated(epsm1q_trcc))  deallocate(epsm1q_trcc)
 if (allocated(cxtab      ))  deallocate(cxtab)
 if (allocated(dummy_xme  ))  deallocate(dummy_xme)
 if (allocated(dummy_cme  ))  deallocate(dummy_cme)
 if (allocated(sigcme_ac  ))  deallocate(sigcme_ac)
 if (allocated(epsm1cqwz2 ))  deallocate(epsm1cqwz2)
 if (allocated(integr     ))  deallocate(integr)
 if (allocated(botsq      ))  deallocate(botsq)
 if (allocated(otq        ))  deallocate(otq) 
 if (allocated(eig        ))  deallocate(eig)

 call status(0,Dtfil%filstat,0,level,'exit          ')
 call timab(426,2,tsec)
 call timab(421,2,tsec)

#if defined DEBUG_MODE
 write(msg,'(a)')' csigme : exit '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine csigme
!!***
