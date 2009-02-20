!{\src2tex{textfont=tt}}
!!****f* ABINIT/spin_current
!! NAME
!! spin_current
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (Mver)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=inverse of atindx
!!  cg(2,mpw*dtset%nspinor*mband*mkmem*nsppol)=wavefunctions
!!                  (may be read from disk instead of input)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gmet = reciprocal space metric
!!  gprimd = dimensionful reciprocal space vectors
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  kg(3,mpw*mkmem)=reduced (integer) coordinates of G vecs in basis sphere
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(dtset%ntypat)=number of atoms of each type
!!  nfftf = fft grid dimensions for fine grid
!!  ph1d = phase factors in 1 radial dimension
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum
!!  rhog(2,nfftf)=Fourier transform of total electron density (including compensation density in PAW)
!!  rhor(nfftf,nspden)=total electron density (including compensation density in PAW)
!!  rmet = real space metric tensor
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  ucvol = unit cell volume
!!  wffnow=unit number for current wf disk file
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!   only output to file
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

subroutine spin_current(atindx,atindx1,cg,dtfil,dtset,eigen,gmet,gprimd,hdr,kg,mpi_enreg,&
  &   nattyp,nfftf,ph1d,psps,rhog,rhor,rmet,symrec,ucvol,wffnow,ylm,ylmgr)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_13nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(in) :: wffnow
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),nattyp(dtset%ntypat)
 integer,intent(in) :: symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: cg(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
 real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol),gmet(3,3)
 real(dp),intent(in) :: gprimd(3,3),rhog(2,nfftf),rhor(nfftf,dtset%nspden)
 real(dp),intent(in) :: rmet(3,3)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)

!Local variables-------------------------------
 ! variables for ph3d mkffnl and company
 ! dummy variables for nonlop
 ! real variables for nonlop
!scalars
 integer :: choice,cplex,cpopt_dummy,dimenl1,dimenl2,dimffnl,fft_option,i1,i1p
 integer :: i2,i2p,i3,i3p,ia,iatom,iband,icartdir,icg,ider,idir_dummy,ig,igp
 integer :: ikg,ikpt,iocc,iost,irealsp,irealsp_p,ispindir,ispinor,ispinorp
 integer :: matblk,npw,only_SO,paw_opt_dummy,signs,spcur_unit
 real(dp) :: arg,lambda_dummy,prefact_nk
 character(len=500) :: message
 character(len=fnlen) :: filnam
!arrays
 integer,allocatable :: gbound(:,:),kg_k(:,:)
 real(dp),allocatable :: density_matrix(:,:,:,:,:),dpsidr(:,:,:,:,:,:)
 real(dp),allocatable :: dummy_denpot(:,:,:),dummy_fofgout(:,:),enlout_dummy(:)
 real(dp),allocatable :: ffnl(:,:,:,:),gpsi(:,:,:,:),kgcart(:,:),kpg_dummy(:,:)
 real(dp),allocatable :: ph3d(:,:,:),phkxred(:,:),position_op(:,:,:,:)
 real(dp),allocatable :: psi(:,:,:),psi_r(:,:,:,:,:),sij_dummy(:,:)
 real(dp),allocatable :: spincurrent(:,:,:,:,:),svectout_dummy(:,:),vectin(:,:)
 real(dp),allocatable :: vectin_ft(:,:),vectout(:,:),vectout_ft(:,:,:,:)
 real(dp),allocatable :: vso_realrecip(:,:,:,:,:),vso_realspace(:,:,:,:,:)
 character :: spin_symbol(3)
 type(cprj_type),allocatable :: cprjin_dummy(:)

! *************************************************************************
!source

 write (*,*) ' Entering subroutine spin_current '
 write (*,*) ' dtset%ngfft = ', dtset%ngfft
 write (*,*) ' hdr%istwfk = ', hdr%istwfk

!===================== init and checks ============================  
!check if nspinor is 2
 if (dtset%nspinor /= 2) then
  write(message, '(a,a,a,a,i6,a,i6,a,a)' ) ch10,&
&  ' spin_current : ERROR -',ch10,&
&  '  nspinor must be 2, but it is ',dtset%nspinor,ch10
  call wrtout(06,message,'PERS')
  call leave_new('PERS')
 end if

 if (dtset%nsppol /= 1) then
  write(message, '(a,a,a,a,i6,a)' ) ch10,&
&  ' spin_current : ERROR -',ch10,&
&  ' nsppol must be 1 but it is ',dtset%nsppol,ch10
  call wrtout(06,message,'PERS')
  call leave_new('PERS')
 end if

 if (dtset%mkmem /= dtset%nkpt) then
  write(message, '(a,a,a,a,i6,a,i6,a,a)' ) ch10,&
&  ' spin_current : ERROR -',ch10,&
&  ' mkmem =  ',dtset%mkmem,' must be equal to nkpt ',dtset%nkpt,&
&  ch10,' keep all kpt in memory'
  call wrtout(06,message,'PERS')
  call leave_new('PERS')
 end if

 if (dtset%usepaw /= 0) then
  write(message, '(a,a,a,a,i6,a,a,a)' ) ch10,&
&  ' spin_current : ERROR -',ch10,&
&  ' usepaw =  ',dtset%usepaw,' must be equal to 0 ',&
&  ch10,' Not functional for PAW case yet.'
  call wrtout(06,message,'PERS')
  call leave_new('PERS')
 end if

 cplex=2
 fft_option = 0 ! just do direct fft
 spin_symbol = (/'x','y','z'/)

 
 write (*,*) ' psps%mpsang,psps%mpssoang ', psps%mpsang,psps%mpssoang

!variables for nonlop
 choice = 1 ! NL energy contribution, not derivatives
 signs = 2 ! get function of G instead of contracted KS matrix element
!only_SO 1 gets the full SO potential  (V_SO L.S) (G,s,G',s')
!only_SO 2 gets a partial SO potential (V_SO   S) (G,s,G',s') then FT wrt G,G'
 only_SO = 2

 cpopt_dummy = -1
 idir_dummy = 0 ! should not be used
 lambda_dummy = zero
 paw_opt_dummy=0

!dimensions for ffnl and nonlop
 dimenl1 = psps%dimekb
 dimenl2 = dtset%ntypat
 dimffnl=1
 matblk=dtset%natom

!allocate stuff for nonlop that does not depend on npw/kpt
 allocate (cprjin_dummy(dtset%natom*((cpopt_dummy+3)/3)))
 allocate (sij_dummy(dimenl1,dtset%ntypat*((paw_opt_dummy+1)/3)))
 allocate (enlout_dummy(1))

!======================= main code ================================  
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!first get normal contribution to current, as psi tau dpsidr + dpsidr tau psi
!where tau are 1/2 the pauli matrices
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!init plane wave coeff counter
 icg = 0
!init plane wave counter
 ikg = 0
!init occupation/band counter
 iocc = 1
 
!rspace point, cartesian direction, spin pol=x,y,z
 allocate (spincurrent(dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),3,3)) 
 spincurrent = zero

 allocate(dummy_denpot(cplex*dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6)))
 allocate(dpsidr(2,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),dtset%nspinor,3))
 allocate(psi_r(2,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),dtset%nspinor))

 allocate(gbound(2*dtset%mgfft+8,2))

 allocate (density_matrix(2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%nspinor,&
& dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%nspinor))
 density_matrix= zero

!loop over kpoints
 do ikpt=1,dtset%nkpt


! number of plane waves for this kpt
  npw = hdr%npwarr(ikpt)

! allocate arrays dep on number of pw
  allocate (kg_k(3,npw))
  allocate (gpsi(2,npw,dtset%nspinor,3)) ! cmplx,ng,nspinor,cartesian direction
  allocate (psi(2,npw,dtset%nspinor))
  allocate (kgcart(3,npw))
  allocate(dummy_fofgout(2,dtset%mpw))

! get cartesian coordinates of k+G vectors around this kpoint
  do ig=1,npw
   kgcart(:,ig) = matmul(gprimd(:,:),dtset%kpt(:,ikpt)+kg(:,ikg+ig))
   kg_k (:,ig) = kg(:,ikg+ig)
  end do

! get gbound
  call sphereboundary(gbound,dtset%istwfk(ikpt),kg_k,dtset%mgfft,npw)

! loop over bands
  do iband=1,dtset%nband(ikpt)

!  prefactor for sum over bands and kpoints
   prefact_nk = hdr%occ(iocc) * dtset%wtk(ikpt)

!  initialize this wf
   gpsi=zero
   psi=zero
   psi(:,1:npw,1) = cg(:,icg+1:icg+npw)

!  multiply psi by G
   do ig=1,npw
    gpsi(:,ig,:,1) = kgcart(1,ig)*psi(:,ig,:) 
    gpsi(:,ig,:,2) = kgcart(2,ig)*psi(:,ig,:) 
    gpsi(:,ig,:,3) = kgcart(3,ig)*psi(:,ig,:) 
   end do

!  loop over spinorial components
   do ispinor=1,dtset%nspinor
!   FT Gpsi_x to real space
    call fourwf(cplex,dummy_denpot,gpsi(:,:,ispinor,1),dummy_fofgout,&
&    dpsidr(:,:,:,:,ispinor,1),gbound,gbound,&
&    hdr%istwfk(ikpt),kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw,&
&    npw,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),&
&    fft_option,dtset%paral_kgb,0,one)

!   FT Gpsi_y to real space
    call fourwf(cplex,dummy_denpot,gpsi(:,:,ispinor,2),dummy_fofgout,&
&    dpsidr(:,:,:,:,ispinor,2),gbound,gbound,&
&    hdr%istwfk(ikpt),kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw,&
&    npw,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),&
&    fft_option,dtset%paral_kgb,0,one)

!   FT Gpsi_z to real space
    call fourwf(cplex,dummy_denpot,gpsi(:,:,ispinor,3),dummy_fofgout,&
&    dpsidr(:,:,:,:,ispinor,3),gbound,gbound,&
&    hdr%istwfk(ikpt),kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw,&
&    npw,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),&
&    fft_option,dtset%paral_kgb,0,one)

!   FT psi to real space
    call fourwf(cplex,dummy_denpot,psi(:,:,ispinor),dummy_fofgout,&
&    psi_r(:,:,:,:,ispinor),gbound,gbound,&
&    hdr%istwfk(ikpt),kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw,&
&    npw,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),&
&    fft_option,dtset%paral_kgb,0,one)

   end do ! ispinor

!  get 3 pauli matrix contributions to the current: x,y,z, cart dir, spin dir
   do icartdir=1,3

!   x pauli spin matrix 
    spincurrent(:,:,:,icartdir,1) =  spincurrent(:,:,:,icartdir,1) + prefact_nk * &
!   Re(psi_r(up)^* dpsidr(down))
&    real(psi_r(1,:,:,:,1)*dpsidr(1,:,:,:,2,icartdir)  &
&    + psi_r(2,:,:,:,1)*dpsidr(2,:,:,:,2,icartdir)  &
!   Re(psi_r(down)^* dpsidr(up))
&    + psi_r(1,:,:,:,2)*dpsidr(1,:,:,:,1,icartdir)  &
&    + psi_r(2,:,:,:,2)*dpsidr(2,:,:,:,1,icartdir))
!   y pauli spin matrix
    spincurrent(:,:,:,icartdir,2) =  spincurrent(:,:,:,icartdir,2) + prefact_nk * &
!   Re(-i psi_r(up)^* dpsidr(down))
&    real(psi_r(1,:,:,:,1)*dpsidr(2,:,:,:,2,icartdir)  &
&    - psi_r(2,:,:,:,1)*dpsidr(1,:,:,:,2,icartdir)  &
!   Re(i psi_r(down)^* dpsidr(up))
&    - psi_r(1,:,:,:,2)*dpsidr(2,:,:,:,1,icartdir)  &
&    + psi_r(2,:,:,:,2)*dpsidr(1,:,:,:,1,icartdir))
!   z pauli spin matrix
    spincurrent(:,:,:,icartdir,3) =  spincurrent(:,:,:,icartdir,3) + prefact_nk * &
!   Re(psi_r(up)^* dpsidr(up))
&    real(psi_r(1,:,:,:,1)*dpsidr(1,:,:,:,1,icartdir)  &
&    - psi_r(2,:,:,:,1)*dpsidr(2,:,:,:,1,icartdir)  &
!   Re(-psi_r(down)^* dpsidr(down))
&    - psi_r(1,:,:,:,2)*dpsidr(1,:,:,:,2,icartdir)  &
&    + psi_r(2,:,:,:,2)*dpsidr(2,:,:,:,2,icartdir))
   end do ! end icartdir

!  
!  accumulate non local density matrix in real space
!  NOTE: if we are only using the local part of the current, this becomes the
!  density! (much lighter to calculate)
!  
   do ispinor=1,dtset%nspinor
    do i3=1,dtset%ngfft(3)
     do i2=1,dtset%ngfft(2)
      do i1=1,dtset%ngfft(1)
       irealsp = i1 + (i2-1)*dtset%ngfft(1) + (i3-1)*dtset%ngfft(2)*dtset%ngfft(1)
       
       do ispinorp=1,dtset%nspinor
        do i3p=1,dtset%ngfft(3)
         do i2p=1,dtset%ngfft(2)
          do i1p=1,dtset%ngfft(1)
           irealsp_p = i1p + (i2p-1)*dtset%ngfft(1) + (i3p-1)*dtset%ngfft(2)*dtset%ngfft(1)
           
           density_matrix(1,irealsp,ispinor,irealsp_p,ispinorp) = &
&           density_matrix(1,irealsp,ispinor,irealsp_p,ispinorp) + &
&           prefact_nk * (psi_r(1,i1,i2,i3,ispinor)*psi_r(1,i1p,i2p,i3p,ispinorp)&
&           -  psi_r(2,i1,i2,i3,ispinor)*psi_r(2,i1p,i2p,i3p,ispinorp))
           density_matrix(2,irealsp,ispinor,irealsp_p,ispinorp) = &
&           density_matrix(2,irealsp,ispinor,irealsp_p,ispinorp) + &
&           prefact_nk * (psi_r(1,i1,i2,i3,ispinor)*psi_r(2,i1p,i2p,i3p,ispinorp)&
&           +  psi_r(2,i1,i2,i3,ispinor)*psi_r(1,i1p,i2p,i3p,ispinorp))
          end do
         end do
        end do
       end do !end ispinorp do
       
      end do
     end do
    end do
   end do !end ispinor do

!  update pw counter
   icg=icg+npw
   iocc=iocc+1
  end do ! iband

  ikg=ikg+npw

! deallocate arrays dep on npw for this kpoint
  deallocate (kg_k)
  deallocate (gpsi,psi,kgcart,dummy_fofgout)

 end do ! ikpt

!prefactor for contribution to spin current
!prefactor is 1/2 * 1/2 * 2 Re(.):
!1/2 from the formula for the current
!1/2 from the use of the normalized Pauli matrices
!2 from the complex conjugate part 
!total = 1/2
 spincurrent = half * spincurrent

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!add electric field term to current. Non local term in case of pseudopotential SO
!present theory is that it is equal to (r V_SO(r,r') + V_SO(r,r')r')
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!choose which kpt we will use to get V_SO (closest to Gamma probably best)
 ikg=0
 do ikpt=1,dtset%nkpt
  if ( sum(abs(dtset%kpt(:,ikpt))) < tol10) exit
  ikg=ikg+hdr%npwarr(ikpt)
 end do
 write (*,*) 'Found Gamma to be ikpt ', ikpt, dtset%kpt(:,ikpt)
 write (*,*) ' ikg = ', ikg

 npw = hdr%npwarr(ikpt)

 allocate (kg_k(3,npw))
 kg_k = kg(:,ikg+1:ikg+npw)


!rebuild phkxred
 allocate (phkxred(2,dtset%natom))
 do ia=1,dtset%natom
  iatom=atindx(ia)
  arg=two_pi*(dtset%kpt(1,ikpt)*hdr%xred(1,ia)&
&  +dtset%kpt(2,ikpt)*hdr%xred(2,ia)&
&  +dtset%kpt(3,ikpt)*hdr%xred(3,ia))
  phkxred(1,iatom)=cos(arg)
  phkxred(2,iatom)=sin(arg)
 end do

!rebuild ph3d
 allocate(ph3d(2,npw,matblk))
 call ph1d3d(1,dtset%natom,kg_k,dtset%kpt(:,ikpt),matblk,dtset%natom,npw,&
& dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),&
& phkxred,ph1d,ph3d)


!rebuild ffnl
 ider=0
 allocate(ffnl(npw,dimffnl,psps%lmnmax,dtset%ntypat))
 allocate (kpg_dummy(npw,0))
 call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
& gmet,gprimd,ider,idir_dummy,psps%indlmn,kg_k,&
& kpg_dummy,dtset%kpt(:,ikpt),psps%lmnmax,&
& psps%lnmax,psps%mpsang,psps%mqgrid_ff,0,&
& npw,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
& psps%usepaw,psps%useylm,ylm,ylmgr)

!get gbound
 call sphereboundary(gbound,dtset%istwfk(ikpt),kg_k,dtset%mgfft,npw)

!allocations for nonlop
 allocate (vectin (2,dtset%nspinor*npw))
 allocate (vectout(2,dtset%nspinor*npw))

 allocate (svectout_dummy(2,dtset%nspinor*npw*(paw_opt_dummy/3)))
 allocate (vectin_ft(2,npw))
 allocate (vectout_ft(2,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6)))
 allocate (vso_realrecip(2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%nspinor,&
& npw,dtset%nspinor))
 allocate (vso_realspace(2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%nspinor,&
& dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%nspinor))

!for each spinorial component   
 do ispinorp=1,dtset%nspinor
! for each planewave G',
  do igp=1,npw
!  make wavefunction with only that component
!  probably to be changed: loop over ks states and call nonlop with them
!  eventually premultiplying with r????
!  
!  Aaaaah maybe not: want full spatial
!  dependency and nonlop gives you a projected quantity summed over the G of the
!  KS state
!  
!  This is actually a barbaric way of extracting the so potential 1 GG' pair
!  at a time
   vectin = zero
   vectin(1,(ispinorp-1)*npw+igp) = one
   
!  and call nonlop -> get <G|V_SO|G'> for all G
!  added flag to not calculate scalar relativistic term, only SO
   call nonlop(atindx1,choice,cpopt_dummy,cprjin_dummy,dimenl1,dimenl2,dimffnl,dimffnl,&
&   psps%ekb,enlout_dummy,ffnl,ffnl,gmet,gprimd,idir_dummy,psps%indlmn,dtset%istwfk(ikpt),&
&   kg_k,kg_k,kpg_dummy,kpg_dummy,dtset%kpt(:,ikpt),dtset%kpt(:,ikpt),&
&   lambda_dummy,psps%lmnmax,matblk,dtset%mgfft,&
&   mpi_enreg,psps%mpsang,psps%mpssoang,dtset%natom,nattyp,dtset%ngfft,0,0,dtset%nloalg,&
&   1,npw,npw,dtset%nspinor,dtset%ntypat,only_SO,paw_opt_dummy,phkxred,&
&   phkxred,ph1d,ph3d,ph3d,psps%pspso,signs,sij_dummy,svectout_dummy,&
&   0,ucvol,psps%useylm,vectin,vectout)
   

!  FT wrt G, one spinorial component of vectout at a time
   do ispinor=1,dtset%nspinor
    vectin_ft = vectout(:,(ispinor-1)*npw+1:(ispinor)*npw)

    call fourwf(cplex,dummy_denpot,vectin_ft,dummy_fofgout,&
&    vectout_ft,gbound,gbound,&
&    hdr%istwfk(ikpt),kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw,&
&    npw,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),&
&    fft_option,dtset%paral_kgb,0,one)

    vso_realrecip(:,:,ispinor,igp,ispinorp)=&
&    reshape(vectout_ft(:,1:dtset%ngfft(1),1:dtset%ngfft(2),1:dtset%ngfft(3)),&
&    (/2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)/))
   end do  ! ispinor
  end do  ! igp

 end do ! ispinorp

 deallocate (kpg_dummy,svectout_dummy)

!FT wrt G'
 do ispinor=1,dtset%nspinor
  do irealsp=1,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)
   do ispinorp=1,dtset%nspinor
    vectin_ft = vso_realrecip(:,irealsp,ispinor,1:npw,ispinorp)

    call fourwf(cplex,dummy_denpot,vectin_ft,dummy_fofgout,&
&    vectout_ft,gbound,gbound,&
&    hdr%istwfk(ikpt),kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw,&
&    npw,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),&
&    fft_option,dtset%paral_kgb,0,one)

    vso_realspace(:,irealsp,ispinor,:,ispinorp) = &
&    reshape(vectout_ft(:,1:dtset%ngfft(1),1:dtset%ngfft(2),1:dtset%ngfft(3)),&
&    (/2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)/))
   end do 
  end do 
 end do
 deallocate (vso_realrecip)

!DEBUG check symmetric quality of vso_realspace
!do ispinor=1,dtset%nspinor
!do irealsp=1,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)
!do ispinorp=1,dtset%nspinor
!do irealsp_p=1,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)
!write (1666,'(2E16.10)') vso_realspace(:,irealsp,ispinor,irealsp_p,ispinorp) &
!&                      - vso_realspace(:,irealsp_p,ispinorp,irealsp,ispinor)
!end do 
!end do
!end do 
!end do
!ENDDEBUG

!make array of positions for all points on grid
 allocate (position_op(3,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3)))
 do i3=1,dtset%ngfft(3)
  do i2=1,dtset%ngfft(2)
   do i1=1,dtset%ngfft(1)
    position_op(:,i1,i2,i3) = matmul(hdr%rprimd,(/i1-one,i2-one,i3-one/))&
&    /(/dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3)/)
   end do
  end do 
 end do

!anticommutator of VSO with position operator

!multiply by density matrix

!add to spin current

 deallocate (kg_k,vectin, vectout) ! these do depend on npw
 deallocate (ffnl,phkxred,ph3d)  ! these do depend on npw


 deallocate (cprjin_dummy,sij_dummy,enlout_dummy) ! nonlop dummies indep of npw
 deallocate (dummy_denpot,dpsidr,psi_r,gbound) ! dummies for fourwf indep of npw

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!output SO potential (non local) for each pair of real space points
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
 filnam=trim(dtfil%filnam_ds(4))//"_VSO_rrp"
 spcur_unit=200
 open (file=filnam,unit=spcur_unit,status='unknown',iostat=iost)
 if (iost /= 0) then
  write (message,'(2a)')' spin_current: ERROR- opening file ',trim(filnam)
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!print header
 write (spcur_unit,'(a)') &
& '#  SO operator (nonlocal) as a function of real space point rprime, for fixed r'
 write (spcur_unit,'(a,3(I5,1x))') &
& '#  fft grid is ', dtset%ngfft(1), dtset%ngfft(2),   dtset%ngfft(3)
 write (spcur_unit,'(a)') &
& '#  cart xprime * cart yprime * cart zprime ***   V_SO '
!
!NOTE: have chosen actual dims of grid (n123) instead of fft box, for which n45
!may be different - forced to be odd for FT
!
 i1=1
 i2=1
 i3=1
 write (spcur_unit,'(a,3(E12.5,1x))') &
& '# position of first r point for V_SO(r,rprime): ', &
& position_op(:,i1,i2,i3)
!look at a given spinorial component of V_SO matrix:
 ispinor=1
 ispinorp=2

!do i3=1,dtset%ngfft(3)
!do i2=1,dtset%ngfft(2)
!do i1=1,dtset%ngfft(1)

 irealsp = i1 + (i2-1)*dtset%ngfft(1) + (i3-1)*dtset%ngfft(2)*dtset%ngfft(1)
 do i3p=1,dtset%ngfft(3)
  do i2p=1,dtset%ngfft(2)
   do i1p=1,dtset%ngfft(1)
    irealsp_p = i1p + (i2p-1)*dtset%ngfft(1) + (i3p-1)*dtset%ngfft(2)*dtset%ngfft(1)
!   write (spcur_unit,'(3(E12.5,1x),3x,3(E12.5,1x),3x,2(E20.10,1x))')&
    write (spcur_unit,'(3(E12.5,1x),3x,2(E20.10,1x))')&
&    position_op(:,i1p,i2p,i3p), &
&    vso_realspace(:,irealsp,ispinor,irealsp_p,ispinorp)
   end do
  end do
 end do

!end do
!end do
!end do

 close (spcur_unit)


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!output 3 components of current for each real space point
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
 do ispindir=1,3
  filnam=trim(dtfil%filnam_ds(4))//"_SPCUR_"//spin_symbol(ispindir)
  spcur_unit=200
  open (file=filnam,unit=spcur_unit,status='unknown',iostat=iost)
  if (iost /= 0) then
   write (message,'(2a)')' spin_current: ERROR- opening file ',trim(filnam)
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

! print header
  write (spcur_unit,'(a)') '#'
  write (spcur_unit,'(a)') '#  spin current density, for all real space points'
  write (spcur_unit,'(a,3(I5,1x))') '#  fft grid is ', dtset%ngfft(1), dtset%ngfft(2),   dtset%ngfft(3)
  write (spcur_unit,'(a,a,a)') '# ', spin_symbol(ispindir), '-spin current, as a vector'
  write (spcur_unit,'(a,a)') '#  cart x     *  cart y    *  cart z    ***',&
&  ' x component of j   *  y component of j  * z component of j   '
! 
! NOTE: have chosen actual dims of grid (n123) instead of fft box, for which n45
! may be different
! 
  do i3=1,dtset%ngfft(3)
   do i2=1,dtset%ngfft(2)
    do i1=1,dtset%ngfft(1)
     write (spcur_unit,'(3(E12.5,1x),3x,3(E20.10,1x))')&
&     position_op(:,i1,i2,i3), &
&     spincurrent(i1,i2,i3,1:3,ispindir)
    end do
   end do
  end do

  close (spcur_unit)

 end do ! end ispindir

 deallocate (vso_realspace)
 deallocate (spincurrent)

 write (*,*) ' Exiting subroutine spin_current '

end subroutine spin_current


!!***
