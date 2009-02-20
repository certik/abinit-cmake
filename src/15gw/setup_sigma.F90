!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_sigma
!! NAME
!! setup_sigma
!!
!! FUNCTION
!!  Initialize the data type containing parameters for a sigma calculation.
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! Dtset<type(dataset_type)>=all input variables for this dataset
!! Dtfil<type(datafiles_type)>=variables related to files
!! rprim(3,3)=dimensionless real space primitive translations
!! MPI_enreg=information about MPI parallelization
!! ngfft(18)=information on the (fine) FFT grid used for the density.
!!
!! OUTPUT
!! Sp<Sigma_parameters> 
!! Kmesh <BZ_mesh_type> 
!! Qmesh <BZ_mesh_type> 
!! Cryst<Crystal_structure>=Info on unit cell and symmetries
!! Gsph_Max<Gvectors_type>=Info on G-sphere
!! Hdr_kss
!! Vcp<Coulombian_type>= datatype gathering information on the coulombian interaction and the cutoff technique.
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine setup_sigma(acell,rprim,ngfftf,Dtset,Dtfil,MPI_enreg,mpsang_kss,ngfft_gw,Hdr_kss,&
& Cryst,Kmesh,Qmesh,Gsph_Max,Vcp,Er,Sp)
    
 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_14iowfdenpot
 use interfaces_15gw, except_this_one => setup_sigma
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: mpsang_kss
 type(Datafiles_type),intent(in) :: Dtfil
 type(MPI_type),intent(inout) :: MPI_enreg
 type(Dataset_type),intent(inout) :: Dtset
 type(Sigma_parameters),intent(inout) :: Sp
 type(Epsilonm1_results),intent(inout) :: Er
 type(BZ_mesh_type),intent(out) :: Kmesh,Qmesh
 type(Crystal_structure),intent(out) :: Cryst
 type(Gvectors_type),intent(out) :: Gsph_Max
 type(Hdr_type),intent(out) :: Hdr_kss 
 type(Coulombian_type),intent(out) :: Vcp
!arrays
 integer,intent(in) :: ngfftf(18)
 integer,intent(out) :: ngfft_gw(18)
 real(dp),intent(in) :: acell(3),rprim(3,3)

!Local variables-------------------------------
!scalars
 integer :: enforce_sym,io,method,mg0sh,mod10,mqmem,nbnds_kss
 integer :: nfftgw_tot,ng_kss,nsym_kss
 integer :: optfil,timrev,umklp_opt,unt
 real(dp),parameter :: OMEGASIMIN=0.01d0
 real(dp) :: domegas,domegasi,ucvol
 logical :: ltest,only_one_kpt,remove_inv
 character(len=500) :: msg                   
 character(len=fnlen) :: fname

!arrays
 integer :: ng0sh_opt(3)
 integer,pointer :: gvec_p(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 
! *************************************************************************
 
#if defined DEBUG_MODE
 write(msg,'(a)')' setup_sigma : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif
 !
 ! === Basic parameters ===
 Sp%ppmodel    = Dtset%ppmodel
 Sp%splitsigc  = Dtset%splitsigc 
 Sp%gwcalctyp  = Dtset%gwcalctyp
 Sp%nbnds      = Dtset%nband(1) 
 Sp%symsigma   = Dtset%symsigma
 Sp%zcut       = Dtset%zcut

 write(msg,'(2a,i2,2a,f10.6,a)')ch10,&
& ' GW calculation type          = ',Sp%gwcalctyp,ch10,&
& ' zcut for avoiding poles [eV] = ',Sp%zcut*Ha_eV,ch10
 call wrtout(std_out,msg,'COLL')
 !
 ! === For HF, SEX or COHSEX use Hybertsen-Louie PPM (only $\omega=0$) ===
 ! * Use fake screening for HF.
 mod10=MOD(Sp%gwcalctyp,10)
 if (mod10==5.or.mod10==6.or.mod10==7) Sp%ppmodel=2
 if (mod10<5.and.MOD(Sp%gwcalctyp,1)/=1) then
  ! * One shot GW (PPM or contour deformation).
  Sp%nomegasrd  =Dtset%nomegasrd 
  Sp%omegasrdmax=Dtset%omegasrdmax
  Sp%deltae     =(2*Sp%omegasrdmax)/(Sp%nomegasrd-1)
 else
  ! * For AC no need to evaluate derivative by finite differences.
  Sp%nomegasrd  =1 
  Sp%omegasrdmax=zero 
  Sp%deltae     =zero
 end if

 write(msg,'(4a,i4,2(2a,f10.6),a)')ch10,&
& ' Parameter used to calculate Sigma derivatives : ',ch10,&
& '  number of points     = ',Sp%nomegasrd,ch10,&
& '  frequency range [eV] = ',Sp%omegasrdmax*Ha_eV,ch10,&
& '  frequency step  [eV] = ',Sp%deltae*Ha_eV,ch10
 call wrtout(std_out,msg,'COLL') 
 !
 !=== For analytic continuation define the number of imaginary frequencies for Sigma ===
 ! * Tests show than more than 12 freqs in the Pade approximant worsen the results!
 if (mod10==1) then 
  Sp%nomegasi  =Dtset%nomegasi
  Sp%omegasimax=Dtset%omegasimax 
  Sp%omegasimin=OMEGASIMIN 
  write(msg,'(4a,i3,2(2a,f8.3),a)')ch10,&
&  ' Parameters for analytic continuation : ',ch10,&
&  '  number of imaginary frequencies for sigma =  ',Sp%nomegasi,ch10,&
&  '  min frequency for sigma on imag axis [eV] =  ',Sp%omegasimin*Ha_eV,ch10,&
&  '  max frequency for sigma on imag axis [eV] =  ',Sp%omegasimax*Ha_eV,ch10
  call wrtout(std_out,msg,'COLL')

  allocate(Sp%omegasi(Sp%nomegasi))
  ! * Logarithmic mesh along the imaginary axis
  !domegasi=(Sp%omegasimax/Sp%omegasimin)**(one/(Sp%nomegasi-1))
  !Sp%omegasi(1)=czero ; ldi=domegasi
  !do io=2,Sp%nomegasi
  ! omega(io)=cmplx(0,ldi*Sp%omegasimin)
  ! Sp%omegasi(io)=ldi*domegasi
  !end do
  ! * Linear mesh along the imaginary axis
  domegasi=Sp%omegasimax/(Sp%nomegasi-1)
  do io=1,Sp%nomegasi
   Sp%omegasi(io)=CMPLX(zero,(io-1)*domegasi)
  end do
  write(msg,'(4a)')ch10,&
&  ' setup_sigma : calculating Sigma(iw)',&
&  ' at imaginary frequencies [eV] (Fermi Level set to 0) ',ch10
  call wrtout(std_out,msg,'COLL') 
  call wrtout(ab_out,msg,'COLL')
  do io=1,Sp%nomegasi
   write(msg,'(2(f10.3,2x))')Sp%omegasi(io)*Ha_eV
   call wrtout(std_out,msg,'COLL') 
   call wrtout(ab_out,msg,'COLL')
  end do

  ltest=(Sp%omegasimax>0.1d-4.and.Sp%nomegasi>0)
  call assert(ltest,'Wrong value of omegasimax or nomegasi',__FILE__,__LINE__)
  if (Sp%gwcalctyp/=1) then 
   ! For AC, only one shot GW is allowed
   write(msg,'(2a)')ch10,' SC-GW with Analytic continuation is not implemented' 
   call wrtout(std_out,msg,'COLL')  
   call leave_new('COLL')
  end if 
 end if 

 if (Sp%symsigma/=0.and.Sp%gwcalctyp>=20) then
  write(msg,'(2a)')ch10,' SC-GW with symmetries is not available' 
  call wrtout(std_out,msg,'COLL') 
  call leave_new('COLL')
 end if
 !
 !=== Setup parameters for Spectral function ===
 Sp%nomegasr  =Dtset%nfreqsp 
 Sp%omegasrmax=Dtset%freqspmax
 if (Sp%nomegasr>0) then
  domegas=2*Sp%omegasrmax/(Sp%nomegasr-1)
  allocate(Sp%omegasf(Sp%nomegasr))
  do io=1,Sp%nomegasr
   Sp%omegasf(io)=-Sp%omegasrmax+domegas*(io-1)
  end do
  write(msg,'(4a,i8,2(2a,f8.3),a)')ch10,&
&  ' Parameters for the calculation of the spectral function : ',ch10,&
&  '  Number of points    = ',Sp%nomegasr,ch10,&
&  '  Max frequency  [eV] = ',Sp%omegasrmax*Ha_eV,ch10,&
&  '  Frequency step [eV] = ',domegas*Ha_eV,ch10
  call wrtout(std_out,msg,'COLL')
 else
  !In indefo all these quantities are set to zero
  !Sp%nomegasr=1 
  !allocate(Sp%omegasf(Sp%nomegasr))
  !Sp%omegasf(1)=0
 end if
 !
 ! === Check input ===
 if (Sp%ppmodel==3.or.Sp%ppmodel==4) then 
  if (Sp%gwcalctyp>=10) then 
   write(msg,'(2a,i3,a)')ch10,&
&   ' sigma : the ppmodel chosen and gwcalctyp ',Dtset%gwcalctyp,' are not yet compatible'
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
  end if
  if (Sp%nspinor==2) then 
   write(msg,'(2a,i3,a)')ch10,&
&   ' sigma : the ppmodel chosen and nspinor ',Sp%nspinor,' are not yet compatible'
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
  end if
 end if 
 !
 ! === Dimensional primitive translations rprimd (from input), gprimd, metrics and unit cell volume ===
 call mkrdim(acell,rprim,rprimd)  
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol) 
 !
 ! === Define consistently npw, nsh, and ecut for wavefunctions and Sigma_X ===
 call setshells(Dtset%ecutwfn, Dtset%npwwfn, Dtset%nshwfn, Dtset%nsym,gmet,gprimd,Dtset%symrel,'wfn',ucvol)
 call setshells(Dtset%ecutsigx,Dtset%npwsigx,Dtset%nshsigx,Dtset%nsym,gmet,gprimd,Dtset%symrel,'mat',ucvol)
 Sp%npwwfn=Dtset%npwwfn 
 Sp%npwx  =Dtset%npwsigx 
 Sp%npwvec=MAX(Sp%npwwfn,Sp%npwx)
 !
 ! === Read parameters of the KSS, verifify them and retrieve all G-vectors ===
 call testlda(Dtset,Dtfil,nsym_kss,nbnds_kss,ng_kss,mpsang_kss,gvec_p,Hdr_kss,MPI_enreg)

 if (Sp%npwvec>ng_kss) then
  Sp%npwvec=ng_kss 
  if (Sp%npwwfn>ng_kss) Sp%npwwfn=ng_kss 
  if (Sp%npwx  >ng_kss) Sp%npwx  =ng_kss
  write(msg,'(5a,3(a,i8,a))')ch10,&
&  ' setup_sigma: WARNING - ',ch10,&
&  ' Number of G-vectors found less than required',ch10,&
&  '  calculation will proceed with npwvec  = ',Sp%npwvec,ch10,&
&  '  calculation will proceed with npwsigx = ',Sp%npwx,ch10,&
&  '  calculation will proceed with npwwfn  = ',Sp%npwwfn,ch10
  call wrtout(std_out,msg,'COLL')
 end if

 if (Sp%nbnds>nbnds_kss) then
  Sp%nbnds      =nbnds_kss 
  Dtset%nband(:)=nbnds_kss
  write(msg,'(6a,i4,a)')ch10,&
&  ' setup_sigma: WARNING - ',ch10,&
&  ' Number of bands found less then required',ch10,&
&  ' calculation will proceed with nbnds = ',nbnds_kss,ch10
  call wrtout(std_out,msg,'COLL')
 end if

 ! === Get important dimensions from the KSS header ===
 call hdr_vs_dtset(Hdr_kss,Dtset) 
 Sp%nsppol =Hdr_kss%nsppol
 Sp%nspinor=Hdr_kss%nspinor 
 Sp%nsig_ab=Hdr_kss%nspinor**2  !FIXME Is it useful calculating only diagonal terms?

 ! === Create crystal_structure data type ===
 remove_inv=(nsym_kss/=Hdr_kss%nsym) 
 timrev=  2 ! This information is not reported in the header
            !1 => do not use time-reversal symmetry 
            !2 => take advantage of time-reversal symmetry

 call nullify_crystal_structure(Cryst)
 call init_Crystal_from_Hdr(Cryst,Hdr_kss,timrev,remove_inv)
 !call hack_symmetries('setup_sigma',Cryst%nsym,Cryst%timrev)
 call print_Crystal_structure(Cryst)

 !==== Set up of the k-points and tables in the whole BZ ===
 call setup_Kmesh(Hdr_kss%nkpt,Hdr_kss%kptns,Cryst,Kmesh,Dtset%prtvol)

 !=== Setup of k-points and bands for the GW corrections ===
 ! * maxbdgw and minbdgw are the Max and min band index for GW corrections over k-points. 
 !   They are used to dimension wfr_gw and calculate the matrix elements.
 if (Dtset%nkptgw==0) then
  ! * All k-points in the IBZ and all bands up to 2*nbv.
  Sp%nkcalc=Kmesh%nibz
  allocate(Sp%xkcalc(3,Sp%nkcalc),Sp%kcalc(Sp%nkcalc))
  allocate(Sp%minbnd(Sp%nkcalc),Sp%maxbnd(Sp%nkcalc))
  Sp%xkcalc(:,:)=Kmesh%ibz(:,:) ; Sp%minbnd(:)=1 ; Sp%minbdgw=1 
  STOP 'Fix nbv in setup_sigma'
#if 0
  !FIXME fix problem with nbv, maybe it is better if I use nele since it is constant 
  ! during the iterations unlike nbv
  Sp%maxbnd(:)=2*MAXVAL(nbv) ; if (2*MAXVAL(nbv)>Sp%nbnds) Sp%maxbnd(:)=Sp%nbnds
  Sp%maxbdgw=MAXVAL(Sp%maxbnd(:)) 
#endif
 else
  ! * Treat only the k-points and bands specified in the input file.
  Sp%nkcalc=Dtset%nkptgw
  allocate(Sp%xkcalc(3,Sp%nkcalc),Sp%kcalc(Sp%nkcalc))
  allocate(Sp%minbnd(Sp%nkcalc),Sp%maxbnd(Sp%nkcalc))
  Sp%xkcalc(:,:)=Dtset%kptgw(:,:)
  Sp%minbnd(:)=Dtset%bdgw(1,:) ; Sp%minbdgw=MINVAL(Sp%minbnd) 
  Sp%maxbnd(:)=Dtset%bdgw(2,:) ; Sp%maxbdgw=MAXVAL(Sp%maxbnd) 
  if (ANY(Sp%maxbnd(:)>Sp%nbnds)) then
   write(msg,'(8a)')ch10,&
&   ' setup_sigma : ERROR - ',ch10,&
&   ' At least one band where the GW corrections are required ',ch10,&
&   ' exceeds the number of treated bands.',ch10,&
&   ' ACtion : Increase the number of bands in the input file (or in the KSS file) '
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
  end if
 end if
 !
 !=== Check if the k-points are in the BZ ===
 !FB Honestly the code is not able to treat k-points, which are not in the IBZ.
 !This extension should require to change the code in different places.
 !Therefore, one should by now prevent the user from calculating sigma for a k-point not in the IBZ.
 umklp_opt=0
 call findk(Sp%nkcalc,Kmesh%nbz,Sp%xkcalc,Kmesh%bz,Sp%kcalc,umklp_opt)
 !call findk(Sp%nkcalc,Kmesh%nibz,Sp%xkcalc,Kmesh%ibz,Sp%kcalc,umklp_opt)

 ! === Read external file and initialize basic dimension of Er% ===
 ! TODO this has to be rewritten, its better if we use mqmem as input variable
 mqmem=0 ; if (Dtset%gwmem/10==1) mqmem=1
 optfil=0 ! Read epsilon^-1

 unt  =Dtfil%unscr
 fname=Dtfil%filscr 
 if (optfil==1) then 
  unt  =Dtfil%unchi0
  fname=Dtfil%filchi0 
 end if

 call init_Er_from_file(Er,optfil,unt,fname,mqmem,Dtset%accesswff,Dtset%localrdwf,MPI_enreg)

 call print_epsilonm1_results(Er,unit=std_out)
 call hdr_vs_dtset(Er%Hdr,Dtset) 

 write(msg,'(2(a,i4),a)')' nomegaer= ',Er%nomega_r,' nomegaei= ',Er%nomega_i,ch10
 call wrtout(std_out,msg,'COLL')

 ! === If required use a matrix for $\Sigma_c$ which is smaller than that stored on file ===
 ! * By default the entire matrix is read and used,
 if (Dtset%npweps>0.or.Dtset%ecuteps>0.or.Dtset%nsheps>0) then
  ! * Define consistently npweps, nsheps, and ecuteps for \Sigma_c according the input
  call setshells(Dtset%ecuteps,Dtset%npweps,Dtset%nsheps,Dtset%nsym,gmet,gprimd,Dtset%symrel,'eps',ucvol)
  if (Dtset%npweps>Er%npwe) then
   write(msg,'(4a,i8,2a,i8)')ch10,&
&   ' setup_sigma : WARNING -',ch10,&
&   '  Number of G-vectors found is less than the value required = ',Dtset%npweps,ch10,&
&   '  Calculation will proceed with Max available npwe = ',Er%npwe
   call wrtout(std_out,msg,'COLL') 
   call wrtout(ab_out,msg,'COLL')
  else 
   ! Redefine the no. of G"s for W.
   Er%npwe=Dtset%npweps
  end if 
 end if
 
 Sp%npwc=Er%npwe          
 if (Sp%npwc>Sp%npwx) then 
  Sp%npwc=Sp%npwx
  write(msg,'(4a)')ch10,&
&  ' setup_sigma : COMMENT -',ch10,&
&  '  Found npw_correlation > npw_exchange, Imposing npwc=npwx '
  call wrtout(std_out,msg,'COLL')
  ! There is a good reason for doing so, see csigme.F90 and the size of the arrays 
  ! rhotwgp and rhotwgp: we need to define a max size and we opt for Sp%npwx.
 end if 
 Er%npwe=Sp%npwc

 ! === Setup of q-mesh in the whole BZ ===
 ! * Stop if a nonzero umklapp is needed to reconstruct the BZ. In this case, indeed, 
 !   epsilon^-1(Sq) should be symmetrized in csigme using a different expression.
 call setup_Qmesh(Er%nqibz,Cryst,Dtset%prtvol,Er%qibz,Qmesh)

 ! === Find optimal value for G-sphere enlargment due to oscillator matrix elements ===
 ! * Here I have to be sure that Qmesh%bz is always inside the BZ, not always true size bz is buggy
 ! * -one is used because we loop over all the possibile differences, unlike screening
 mg0sh=5
 call get_ng0sh(Sp%nkcalc,Sp%xkcalc,Kmesh%nbz,Kmesh%bz,Qmesh%nbz,Qmesh%bz,gmet,-one,mg0sh,ng0sh_opt)
 Sp%mG0(:)=ng0sh_opt(:)

 ! === Make biggest G-sphere of Sp%npwvec vectors ===
 only_one_kpt=(Kmesh%nbz==1)
 call nullify_Gvectors(Gsph_Max)
 call init_Gvectors_type(only_one_kpt,Gsph_Max,Cryst,Sp%npwvec,gvec_p,gmet,gprimd)
 deallocate(gvec_p)

 ! === Get Fourier components of the Coulombian for all q-points in the IBZ ===
 ! * If required, use a cutoff in the interaction 
 ! * Pcv%vc_sqrt contains Vc^{-1/2}
 ! * Setup also the analytical calculation of the q->0 component
 Qmesh%nsmall=0 !FIXME this has to be done in a clever way
 ! FIXME recheck ngfftf since I got different charge outside the cutoff region
 call setup_coulombian(Dtset,Gsph_Max,Qmesh,Kmesh,Sp%npwx,Cryst%rprimd,ngfftf,MPI_enreg,Vcp)

 ! === Setup of the FFT mesh for the oscilator strengths === 
 ! * ngfft_gw(7:18)==Dtset%ngfft(7:18) which is initialized before entering screening.
 ! * Here we redefine ngfft_gw(1:6) according to the following options :
 !
 ! method==0 --> FFT grid read from fft.in (debugging purpose)
 ! method==1 --> Normal FFT mesh
 ! method==2 --> Slightly augmented FFT grid to calculate exactly rho_tw_g (see setmesh.F90)
 ! method==3 --> Doubled FFT grid, same as the the FFT for the density,
 !
 ! enforce_sym==1 ==> Enforce a FFT mesh compatible with all the symmetry operation and FFT library
 ! enforce_sym==0 ==> Find the smallest FFT grid compatbile with the library, do not care about symmetries

 ngfft_gw(1:18)=Dtset%ngfft(1:18) 
 method=2
 if (Dtset%fftgw==00 .or. Dtset%fftgw==01) method=0
 if (Dtset%fftgw==10 .or. Dtset%fftgw==11) method=1
 if (Dtset%fftgw==20 .or. Dtset%fftgw==21) method=2
 if (Dtset%fftgw==30 .or. Dtset%fftgw==31) method=3
 enforce_sym=MOD(Dtset%fftgw,10) 

 call setmesh(gmet,Gsph_Max%gvec,ngfft_gw,Sp%npwvec,MAX(Sp%npwx,Sp%npwc),Sp%npwwfn,&
& nfftgw_tot,method,Sp%mG0,Cryst,enforce_sym)

 call print_ngfft(ngfft_gw,'FFT mesh for oscillator strengths')

#if defined DEBUG_MODE
 write(msg,'(a)')' setup_sigma : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine setup_sigma
!!***
