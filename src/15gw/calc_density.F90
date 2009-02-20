!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_density
!! NAME
!! calc_density
!!
!! FUNCTION
!! Calculate the charge density rhor on the fine FFT grid in real space.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  igfftf(Wfs%npwwfn)=index of G in the fine FFT grid
!!  irottbf(nfftf,nsym)= irottbf(r,R) gives the index of $R^{-1}(r-\tau)$ in the fine FFT array 
!!  where R is one of the nsym symmetry operation in real space and \tau is the associated fract. translation. 
!!  MPI_enreg= datatype gathering information on parallelism.
!!     %gwpara= if 2 bands are spread btw processors
!!  nbnds=number of bands.
!!  ngfftf(18)=array containing all the information for the "fine" FFT 
!!  Cryst<Crystal_structure> Info on the crystalline structure
!!     %nsym=number of symmetry operations.
!!     %ucvol=unit cell volume.
!!  nfftf=total number of points on the fine FFT grid (for this processor)
!!  occ(%nibz,nbnds,%nsppol)= occupation numbers.
!!  Kmesh<bz_mesh_type>= Info on the k-sampling
!!     %nibz=number of irreducible k-point
!!     %nbz=number of k-points in the full Brillouin zone.
!!     %wt(nibz)=irreducible k-points weights.
!!     %timrev=2 if time-reversal symmetry is used, 1 otherwise.
!!  tim_fourdp=timing code of fourdp (can be set to 0 if not attributed)
!!  use_MPI=only used in case of gwpara==2. If .FALSE. do not communicate,
!!   since all processors have valence states (only used in screening)
!!  Wfs<wavefunctions_information)=datatype gathering info on wavefunctions
!!     %nspinor=number of spinorial components
!!     %nsppol=1 for unpolarized, 2 for spin-polarized
!!     %nspden=number of spin-density components
!!
!! OUTPUT
!!  rhor(nfftf,%nspden)=the density in the real space on the fine FFT grid. 
!!   If nsppol==2 total charge in first half, spin-up component in second half.
!!
!! NOTES 
!! In case of PAW calculations:
!!    All computations are done on the fine FFT grid.
!!    All variables (nfftf,ngfftf,mgfftf) refer to this fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!    Developpers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid.
!! In case of norm-conserving calculations: 
!!    The mesh is the usual augmented FFT grid to treat correctly the convolution.
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine calc_density(Wfs,Cryst,irottbf,nbnds,ngfftf,nfftf,igfftf,&
& occ,rhor,Kmesh,MPI_enreg,tim_fourdp,use_MPI)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : j_gw
 use m_errors, only : assert
 use m_io_tools, only : flush_unit, get_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13recipspace
 use interfaces_15common
 use interfaces_15gw, except_this_one => calc_density
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbnds,nfftf,tim_fourdp
 logical,intent(in) :: use_MPI
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Crystal_structure),intent(in) :: Cryst
 type(MPI_type),intent(inout) :: MPI_enreg
 type(Wavefunctions_information),intent(inout) :: Wfs
!arrays
 integer,intent(in) :: igfftf(Wfs%npwwfn),irottbf(nfftf,Cryst%nsym),ngfftf(18)
 real(dp),intent(in) :: occ(Kmesh%nibz,nbnds,Wfs%nsppol)
 real(dp),intent(out) :: rhor(nfftf,Wfs%nspden)

!Local variables ------------------------------
!scalars
 integer :: cplex,ib,ier,ierr,ik,iop,ir,is,ispinor,master,n1,n2,n3,nfftotf
 integer :: npwwfn,nspden,nspinor,nsppol,nsym,rank,spaceComm,unt
 real(dp) :: ucvol
 logical :: DEBUG,ltest,use_fineFFT
 character(len=100) :: frmt
 character(len=500) :: msg
 character(len=fnlen) :: filnam
 type(Dens_sym_operator_type) :: densymop
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer,allocatable :: irrzon(:,:,:)
 real(dp),allocatable :: phnons(:,:,:),rhog(:,:),rhor_down(:),rhor_mx(:)
 real(dp),allocatable :: rhor_my(:)
 complex(gwpc),allocatable :: wfr_x(:),wfr_y(:)
 complex(gwpc),allocatable,target :: wfr(:)
 complex(gwpc),pointer :: cwavef1(:),cwavef2(:)

!*************************************************************************

#if defined DEBUG_MODE
 write(msg,'(2a)')ch10,&
& ' calc_density: calculating charge density'
 call wrtout(std_out,msg,'COLL')
 call flush_unit(std_out)
#endif
 !
 ! === Initialize some MPI related variables ===
 call xcomm_init  (MPI_enreg,spaceComm)  
 call xme_init    (MPI_enreg,rank     )          
 call xmaster_init(MPI_enreg,master   )

 nspinor = Wfs%nspinor
 nspden  = Wfs%nspden
 nsppol  = Wfs%nsppol
 npwwfn  = Wfs%npwwfn

 nsym    = Cryst%nsym
 ucvol   = Cryst%ucvol

 use_fineFFT=.FALSE.
 if (ANY(Wfs%ngfft(1:3)/=ngfftf(1:3))) then 
  use_fineFFT=.TRUE.
  write(msg,'(a,3(i3,1x))')' calc_density : using fine FFT grid ',ngfftf(1:3)
  call wrtout(std_out,msg,'COLL')
 end if
 !
 ! === Calculate IBZ contribution to the charge density ===
 allocate(wfr(nfftf*nspinor))

 if (nspinor==2) then
  allocate(wfr_x(nfftf),wfr_y(nfftf))
  if (nspden==4) then
   allocate(rhor_down(nfftf),rhor_mx(nfftf),rhor_my(nfftf))
   rhor_down(:)=zero 
   rhor_mx  (:)=zero 
   rhor_my  (:)=zero 
  else 
   !TODO
   call assert(.FALSE.,'nspden and nspinor=1 not implemeted yet') 
  end if
 end if

 ! === Get unsymmetrized density ===
 rhor(:,:)=zero 
 do is=1,nsppol
  do ik=1,Kmesh%nibz
   do ib=1,nbnds
    !
    ! * Skip if occupation is less than tol8 or if band does not belong to rank
    if (MPI_enreg%gwpara==2.and.use_MPI) then
     if (MPI_enreg%proc_distrb(ik,ib,is)/=rank) CYCLE
    end if
    if (ABS(occ(ik,ib,is))<tol8) CYCLE
    !
    ! === Get wavefunction in real space ===
    if (use_fineFFT) then 
     call fft_onewfn(Wfs%paral_kgb,nspinor,npwwfn,nfftf,Wfs%wfg(:,ib,ik,is),wfr,&
&     igfftf,ngfftf,tim_fourdp,MPI_enreg) 
    else
     call get_wfr(Wfs,MPI_enreg,ib,ik,is,wfr)
    end if

    cwavef1 => wfr(1:nfftf) 
    do ir=1,nfftf
     rhor(ir,is)=rhor(ir,is)+occ(ik,ib,is)*CONJG(cwavef1(ir))*cwavef1(ir)*Kmesh%wt(ik)/ucvol
    end do

    if (nspinor==2) then 
     cwavef2 => wfr(1+nfftf:2*nfftf)
     ! $(\Psi^{1}+\Psi^{2})$
     wfr_x(:)=cwavef1(:)+cwavef2(:)
     ! $(\Psi^{1}-i\Psi^{2})$
     wfr_y(:)=cwavef1(:)-j_gw*cwavef2(:)
     do ir=1,nfftf
      rhor_down(ir)=rhor_down(ir)+occ(ik,ib,is)*CONJG(cwavef2(ir))*cwavef2(ir)*Kmesh%wt(ik)/ucvol
      rhor_mx  (ir)=rhor_mx  (ir)+occ(ik,ib,is)*CONJG(wfr_x  (ir))*wfr_x  (ir)*Kmesh%wt(ik)/ucvol
      rhor_my  (ir)=rhor_my  (ir)+occ(ik,ib,is)*CONJG(wfr_y  (ir))*wfr_y  (ir)*Kmesh%wt(ik)/ucvol
     end do
    end if

   end do !ib
  end do !ik
 end do !is

 if (nspden==4) then
  !rhor(:,1) now contains rho_up
  rhor(:,2)=rhor_mx(:)
  rhor(:,3)=rhor_my(:)
  rhor(:,4)=rhor_down(:)
  !HACK to get collinear case. 
  !rhor(:,2)=two*rhor(:,1)
  !rhor(:,3)=two*rhor(:,1)
  !rhor(:,4)=    rhor(:,1)
  !write(*,*)' TEST DIAG RHO     ',SUM(rhor(:,1))*ucvol/nfftf,SUM(rhor(:,4))*ucvol/nfftf
  !write(*,*)' TEST OFF-DIAG RHO ',SUM(rhor(:,2))*ucvol/nfftf,SUM(rhor(:,3))*ucvol/nfftf
 end if
 if (use_MPI.and.MPI_enreg%gwpara==2) call xsum_mpi(rhor,spaceComm,ierr)
 !
 ! === Symmetrization in G-space implementing also the AFM case ===
 n1=ngfftf(1) 
 n2=ngfftf(2) 
 n3=ngfftf(3) 
 nfftotf=n1*n2*n3

 allocate(irrzon(nfftotf**(1-1/nsym),2,nspden/nsppol))
 allocate(phnons(2,nfftotf,nspden/nsppol))

 if (nsym/=1) then
  call irrzg(densymop,irrzon,nspden,nsppol,nsym,n1,n2,n3,&
&  phnons,Cryst%symafm,Cryst%symrel,Cryst%tnons)
 end if

 ! * Fake MPI_type for sequential part
 call initmpi_seq(MPI_enreg_seq) ; MPI_enreg_seq%nproc_fft=1 ; MPI_enreg_seq%me_fft=0

 cplex=1 
 allocate(rhog(2,cplex*nfftf)) !this might be output

 call symrhg(cplex,densymop,irrzon,MPI_enreg_seq,nfftf,nfftotf,ngfftf,nspden,nsppol,&
& nsym,Wfs%paral_kgb,phnons,rhog,rhor,Cryst%symafm)

 deallocate(rhog)
 deallocate(phnons,irrzon)

 write(msg,'(a,f9.4)')&
& ' planewave contribution to nelect: ',SUM(rhor(:,1))*ucvol/nfftf 
 call wrtout(std_out,msg,'COLL') 

 if (nspden==4) then
 write(msg,'(a,3f9.4)')&
&  ' mx, my, mz: ',SUM(rhor(:,2))*ucvol/nfftf,SUM(rhor(:,3))*ucvol/nfftf,SUM(rhor(:,4))*ucvol/nfftf 
  call wrtout(std_out,msg,'COLL') 
 end if

 DEBUG=.FALSE.
 if (DEBUG.and.rank==master) then 
  filnam='__rho__.dat' ; call isfile(filnam,'new')
  unt=get_unit() ; open(unit=unt,file=filnam)
  write(frmt,*)'(2x,',Wfs%nspden,'(1x,f8.3))'
  do ir=1,nfftf
   write(unt,frmt)(rhor(ir,:))
  end do 
  close(unt)
 end if 

 deallocate(wfr)

 if (nspinor==2) then
  deallocate(wfr_x,wfr_y)
  if (nspden==4) deallocate(rhor_down,rhor_mx,rhor_my)
 end if

#if defined DEBUG_MODE
 write(msg,'(a)')' calc_density : ended'
 call wrtout(std_out,msg,'COLL')
 call flush_unit(std_out)
#endif

end subroutine calc_density
!!***


!!****f* ABINIT/test_charge
!! NAME
!! test_charge
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
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine test_charge(nfftf,nspden,rhor,ucvol,nhat,&
& usepaw,usexcnhat,usefinegrid,compch_sph,compch_fft,omegaplasma)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,nspden,usefinegrid,usepaw,usexcnhat
 real(dp),intent(in) :: compch_fft,compch_sph,ucvol
 real(dp),intent(out) :: omegaplasma
!arrays
 real(dp),intent(in) :: nhat(nfftf,nspden*usepaw),rhor(nfftf,nspden)

!Local variables ------------------------------
!scalars
 real(dp) :: nel,nel_fft,nel_pw,nel_sph,rhoav,rs
 character(len=500) :: msg

!*************************************************************************

 ! === For PAW output of compensation charges ===
 if (usepaw==1.and.usexcnhat>0) then ! TODO I still dont understand this if!  
  write(msg,'(4a)')ch10,' PAW TEST:',ch10,&
&  ' ==== Compensation charge inside spheres ============'
  if (compch_sph<greatest_real.and.compch_fft<greatest_real) &
&  write(msg,'(3a)')TRIM(msg),ch10,' The following values must be close...'
  if (compch_sph<greatest_real) write(msg,'(3a,f22.15)')TRIM(msg),ch10,&
&  ' Compensation charge over spherical meshes = ',compch_sph
  if (compch_fft<greatest_real) then
   if (usefinegrid==1) then
    write(msg,'(3a,f22.15)')TRIM(msg),ch10,&
&    ' Compensation charge over fine fft grid    = ',compch_fft
   else
    write(msg,'(3a,f22.15)')TRIM(msg),ch10,&
&    ' Compensation charge over fft grid         = ',compch_fft
   end if
  end if
  call wrtout(ab_out,msg,'COLL') ; call wrtout(std_out,msg,'COLL')
  write(msg,'(a)')ch10
  call wrtout(ab_out,msg,'COLL') ; call wrtout(std_out,msg,'COLL')
 end if

 nel_pw=SUM(rhor(:,1))*ucvol/nfftf ; nel=nel_pw
 if (usepaw==1) then 
  nel_sph=nel_pw+compch_sph
  nel_fft=nel_pw+compch_fft
  nel=nel_sph
 end if
 rhoav=nel/ucvol ; rs=(three/(four_pi*rhoav))**third
 if (usepaw/=1) then 
  write(msg,'(a,f9.4)')' total number of electrons per unit cell = ',nel
 else 
  write(msg,'(2(a,f9.4),a)')&
&  ' total number of electrons per unit cell = ',nel_sph,' (Spherical mesh), ',nel_fft,' (FFT mesh)'
 end if
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,f9.6)')' average of density, n = ',rhoav
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,f9.4)')' r_s = ',rs
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 omegaplasma=SQRT(four_pi*rhoav)
 write(msg,'(a,f9.4,2a)')' omega_plasma = ',omegaplasma*Ha_eV,' [eV]',ch10
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

end subroutine test_charge
!!***
