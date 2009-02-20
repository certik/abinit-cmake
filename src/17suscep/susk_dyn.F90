!{\src2tex{textfont=tt}}
!!****f* ABINIT/susk_dyn
!! NAME
!! susk_dyn
!!
!! FUNCTION
!! Compute the contribution of one k point to the imaginary-frequency
!! susceptibility matrix from input wavefunctions, band occupations,
!! and k point wts. Include the usual sum-over-state terms, but also the
!! corrections due to the change of the Fermi level in the metallic
!! case, as well as implicit sum over higher lying conduction
!! states, thanks to the closure relation (referred to as an extrapolation).
!! Compared to the routine suskmm, there is no particular attention
!! to the use of the memory, so the code is simpler.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MF, XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  bdtot_index=index for the number of the band
!!  cg(2,mcg)=wfs in G space
!!  doccde(mband*nkpt*nsppol)=derivative of occupancies wrt
!!           the energy for each band and k point
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  extrap: if==1, the closure relation (an extrapolation) must be used
!!  freq(nfreq)=array for frequencies (hartree)
!!  gbound(2*mgfftdiel+8,2)=G sphere boundary for going from WF sphere to
!!      medium size FFT grid
!!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for going from medium size
!!      FFT grid to small sphere.
!!  icg=index for cg
!!  ikpt=number of the k point
!!  isp=number of the current spin
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kg_k(3,npw)=coordinates of planewaves in basis sphere.
!!  mband=maximum number of bands
!!  mcg=dimension of cg
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of
!!     the dielectric matrix
!!  mkmem=maximum number of k points in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum allowed value for npw
!!  nband_k=number of bands at this k point for that spin polarization
!!  ndiel4,ndiel5,ndiel6= FFT dimensions, modified to avoid cache trashing
!!  nfftdiel=number of fft grid points for the computation of the diel matrix
!!  nfreq=size of frequency grid
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  npwdiel=third and fifth dimension of the susmat array.
!!  npw_k=number of plane waves at this k point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(mband*nkpt*nsppol)=
!!          occupation numbers for each band (usually 2.0) at each k point
!!  occopt=option for occupancies
!!  occ_deavg(mband)=factor for extrapolation (occup. divided by an energy gap)
!!  occ_freq(2,mband,nfreq)=array holding average weights for extrapolation
!!  susopt=option for susceptibility matrix
!!        =0 static case (not used)
!!        =1 dynamical case
!!        =2 susmat_dyn holds square modulus of the one-particle density matrix
!!  ucvol=unit cell volume (Bohr**3)
!!  wtk(nkpt)=k point weights (they sum to 1.0)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! These quantities are accumulated in this routine:
!! drhode(2,npwdiel,nsppol)=weighted density, needed to compute the
!!   effect of change of fermi energy
!!  rhoextrap_dyn(2,ndiel4,ndiel5,ndiel6,ifreq)=density-like array,
!!   needed for extrapolation procedure.
!!  sumdocc=sum of weighted occupation numbers, needed to compute the
!!   effect of change of fermi energy
!!  susmat_dyn(2,npwdiel,nsppol,npwdiel,nsppol,nfreq)=the frequency
!!   dependent susceptibility matrix in reciprocal space
!!
!! NOTES
!! The array susmat_dyn can be very large, so that it is not advisable
!! to treat more than one frequency at time.
!!
!! WARNINGS
!! a - Argument occopt>=3 not allowed, since contributions due to Fermi
!!     level changes are not implemented.
!! b - For susopt==2 the square of the one-particle density matrix is
!!     computed (and later used to get the exchange energy). The present
!!     implementation is meant for closed shell cases only.
!!
!! TODO
!! a - resolve Warnings
!! b - this routine extends the susk routine, and might be merged with it
!! c - time and try to optimize wfprod in real space for large fft meshes,
!!     possibly save time by precomputing the wavefunction products for
!!     all bands, might require to put/take them to/from disk
!! e - why not save memory by making susmat_dyn single precision?
!!
!! PARENTS
!!      suscep_dyn,suscep_kxc_dyn
!!
!! CHILDREN
!!      cgerc,cher,fourwf,timab,zgerc,zher
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine susk_dyn(bdtot_index,cg,doccde,drhode,dtset,eigen,extrap,freq,&
&  gbound,gbound_diel,icg,ikpt,isp,istwfk,kg_diel,kg_k,&
&  mband,mcg,mgfftdiel,mkmem,mpi_enreg,mpw,&
&  nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,nfreq,ngfftdiel,nkpt,&
&  npwdiel,npw_k,nspden,nspinor,nsppol,occ,occopt,&
&  occ_deavg,occ_freq,rhoextrap_dyn,sumdocc,susmat_dyn,susopt,ucvol,wtk)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_12ffts
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!This type is defined in defs_mpi
!scalars
 integer,intent(in) :: bdtot_index,extrap,icg,ikpt,isp,mband,mcg,mgfftdiel
 integer,intent(in) :: mkmem,mpw,nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,nfreq
 integer,intent(in) :: nkpt,npw_k,npwdiel,nspden,nspinor,nsppol,occopt,susopt
 real(dp),intent(in) :: ucvol
 real(dp),intent(inout) :: sumdocc
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: gbound(2*mgfftdiel+8,2),gbound_diel(2*mgfftdiel+8,2)
 integer,intent(in) :: istwfk(nkpt),kg_diel(3,npwdiel),kg_k(3,npw_k)
 integer,intent(in) :: ngfftdiel(18)
 real(dp),intent(in) :: cg(2,mcg),doccde(mband*nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),freq(nfreq)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),occ_deavg(mband)
 real(dp),intent(in) :: occ_freq(2,mband,nfreq),wtk(nkpt)
 real(dp),intent(inout) :: drhode(2,npwdiel,nsppol)
 real(dp),intent(inout) :: rhoextrap_dyn(2,ndiel4,ndiel5,ndiel6,nfreq)
 real(dp),intent(inout) :: susmat_dyn(2,npwdiel,nsppol,npwdiel,nsppol,nfreq)

!Local variables-------------------------------
! real(dp), allocatable :: cg_disk(:,:)
!scalars
 integer :: i1,i2,i3,iband,ibd1,ibd2,ifreq,ipw,ipw1,ipw2,isp1,isp2,istwf_k
 integer :: ndiel1,ndiel2,ndiel3,testocc,tim_fourwf
 real(dp) :: ai,ar,eigdiff,norm,normr,occdiff,tolocc,weight,wght1,wght2
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwavef(:,:),dummy(:,:),rhoaug(:,:,:),wfprod(:,:)
 real(dp),allocatable :: wfraug(:,:,:,:),wfrspa(:,:,:,:,:),wght_dyn(:,:)

! *************************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZGERC' :: zgerc
#endif

!DEBUG
!write(6,*)'%susk_dyn: ENTER '
!write(6,*)' extrap:', extrap
!write(6,*)' nfreq :', nfreq
!write(6,*)' mpw   :', mpw
!write(6,*)' npw_k :', npw_k
!write(6,*)' occopt:', occopt
!write(6,*)' susopt:', susopt
!if(.true.)stop
!ENDDEBUG

 if(occopt>=3) then
  write(6,*) ' %susk_dyn: occopt>=3 not implemented: stop ' !MF
  stop
 end if

 call timab(87,1,tsec)

 ndiel1=ngfftdiel(1) ; ndiel2=ngfftdiel(2) ; ndiel3=ngfftdiel(3)
 istwf_k=istwfk(1)

 testocc=1
 weight=0._dp

!DEBUG
!write(6,*)' susk : set testocc to 0 '
!testocc=0
!write(6,*)' susk : set extrap to 0 '
!extrap=0
!ENDDEBUG

!DEBUG
!write(6,*) ' istwf_k=',istwf_k
!ENDDEBUG

 allocate(cwavef(2,mpw),dummy(2,1))
 allocate(rhoaug(ndiel4,ndiel5,ndiel6),wfraug(2,ndiel4,ndiel5,ndiel6))
!MF allocate for used bands 'nband_k' only
!allocate(wfprod(2,npwdiel),wfrspa(2,ndiel4,ndiel5,ndiel6,mband))
 allocate(wfprod(2,npwdiel),wfrspa(2,ndiel4,ndiel5,ndiel6,nband_k))
 allocate(wght_dyn(2,nfreq))

 rhoaug(:,:,:)=0._dp
 wfraug(:,:,:,:)=0._dp
 wfprod(:,:)=0._dp
 wfrspa(:,:,:,:,:)=0._dp
 wght_dyn(:,:)=0.0_dp

!Loop over bands to fft and store Fourier transform of wavefunction

!DEBUG
!if(extrap==1)then
!write(6,*)' %susk_dyn: extrapolation: weights for diagonal term:'
!write(6,*)'      row1: ifreq, wght_dyn | row2: weight (==static case)'
!end if
!ENDDEBUG

 do iband=1,nband_k
! Obtain Fourier transform in fft box
  cwavef(:,1:npw_k)=cg(:,1+(iband-1)*npw_k+icg:iband*npw_k+icg)

! DEBUG
! write(6,*)' susk_dyn : will call fourwf for band',iband
! ENDDEBUG

  tim_fourwf=10
  call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&  istwf_k,kg_k,kg_k,mgfftdiel,mpi_enreg,1,ngfftdiel,npw_k,1,ndiel4,ndiel5,ndiel6,&
&  0,dtset%paral_kgb,tim_fourwf,weight)
  wfrspa(:,:,:,:,iband)=wfraug(:,:,:,:)

  if( (occopt>=3 .and. testocc==1) .or. extrap==1 .or. susopt==2)then
!  In the case of metallic occupation, or if the extrapolation
!  over higher bands is included, must compute the
!  Fourier transform of the density of each band, then
!  generate the part of the susceptibility matrix due
!  varying occupation numbers.

!  MF Not needed as occopt>=3 not implemented
!  weight=-2.0_dp*occ_deavg(iband)*wtk(ikpt)/ucvol

!  DEBUG
!  write(6,*)' susk : debug one band contribution '
!  weight=0.0_dp
!  if(iband==1)weight=-2.0_dp*occ_deavg(iband)*wtk(ikpt)/ucvol
!  ENDDEBUG

!  Case where the band density does not need to be accumulated
   if(extrap==0 .or. susopt==2) then
    do i3=1,ndiel3
     do i2=1,ndiel2
      do i1=1,ndiel1
       wfraug(1,i1,i2,i3)=wfraug(1,i1,i2,i3)**2+wfraug(2,i1,i2,i3)**2
       wfraug(2,i1,i2,i3)=0.0_dp
      end do
     end do
    end do

!   Accumulate density in real space
   else if(extrap==1) then
    do ifreq=1,nfreq
     wght_dyn(:,ifreq)=2.0_dp*occ_freq(:,iband,ifreq)*wtk(ikpt)/ucvol

!    DEBUG
!    write(6,*) ' %susk_dyn: show on-diagonal weights:'
!    write(6,*) '  bands: iband    :',iband
!    write(6,*) '  occ_freq(1)(2)  :',occ_freq(:,iband,ifreq)
!    write(6,*) '  wght_dyn(1)(2)  :',wght_dyn(:,ifreq)
!    ENDDEBUG

!    MF  use that rhoextrap is real
     do i3=1,ndiel3
      do i2=1,ndiel2
       do i1=1,ndiel1
        wfraug(1,i1,i2,i3)=wfraug(1,i1,i2,i3)**2+wfraug(2,i1,i2,i3)**2
        rhoextrap_dyn(1,i1,i2,i3,ifreq)=rhoextrap_dyn(1,i1,i2,i3,ifreq) &
&        +wght_dyn(1,ifreq)*wfraug(1,i1,i2,i3)
        wfraug(2,i1,i2,i3)=0.0_dp
       end do
      end do
     end do
!    End frequency loop
    end do

!   End condition extrap or susopt
   end if

!  Tentative call to rhohxc for ACFD-LDA
!  output : enxc,kxc,strsxc,vhartr,vxc,vxcavg,k3xc
!  input  : option=3 (Warning : new option, kxc without Hartree !)
!  ixc (to take it from the parent routine)
!  rhog(:,:)=zero
!  rhodiel(ndiel1,ndiel2,ndiel3) to be transferred to a single index array
!  rprimd to take from the parent
!  xccc3d(:,:,:)=zero
!  
!  call rhohxc(dtset,enxc,kxc,mpi_enreg,nfftdiel,ngfftdiel,&
!  &   1,nspden,0,3,rhog,rhor,rprimd,strsxc,&
!  & vhartr,vxc,vxcavg,xccc3d,k3xc)

!  DEBUG
!  if(iband==1)then
!  write(6,*)' wfraug ='
!  do i3=1,ndiel3,4
!  write(6,*)1,1,i3,wfraug(1,1,1,i3)
!  end do
!  end if
!  ENDDEBUG

!  Performs the Fourier Transform of the density of the band,
!  and store it in wfprod
   tim_fourwf=11
   call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&   istwf_k,kg_diel,kg_diel,&
&   mgfftdiel,mpi_enreg,1,ngfftdiel,1,npwdiel,ndiel4,ndiel5,ndiel6,3,dtset%paral_kgb,tim_fourwf,weight)

!  For density matrix related case
   if(susopt==2) then
    ifreq=1
    weight=2._dp*wtk(ikpt)/ucvol
    do ipw2=1,npwdiel
     do ipw1=ipw2,npwdiel
      susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)+&
&      weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
      susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)+&
&      weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
     end do
    end do

!   For susceptibility matrix case, perform now the summation of terms related
!   to direct change of eigenvalues or extrapolation over higher bands
   else
!   Loop over frequencies
    do ifreq=1,nfreq

!    MF  Not needed as occopt>=3 not implemented
!    wght1=0.0_dp ; wght2=0.0_dp
!    if(occopt>=3 .and. testocc==1) wght1=doccde(iband+bdtot_index)*wtk(ikpt)/ucvol
!    if(extrap==1) wght2=2.0_dp*occ_deavg(iband)*wtk(ikpt)/ucvol
!    weight=wght1+wght2
!    MF

!    DEBUG
!    write(6,*)' susk : debug '
!    if(extrap==1 .and. iband==1) &
!    &    wght2=2.0_dp*occ_deavg(iband)*wtk(ikpt)/ucvol
!    ENDDEBUG

     if(extrap==1) wght_dyn(:,ifreq)=-wght_dyn(:,ifreq)

!    Update susmat_dyn explicitely-------------------------
     if(.false.) then
      do ipw2=1,npwdiel
!      Only fills lower half of the matrix (here, the susceptibility matrix)
!      Note that wfprod of the first index must behave like a density,
!      so that it is used as generated by fourwf, while wfprod of the
!      second index will be implicitely used to make a scalar product
!      with a potential change, meaning that its complex conjugate must be
!      used. This explains the following signs...

!      DEBUG
!      write(6, '(a,i3,2es14.6)' )&
!      &      ' ipw,wfprod',ipw2,wfprod(1,ipw2),wfprod(2,ipw2)
!      ENDDEBUG

!      MF    fill full matrix explicite to check time reversal symmetry in suscep_dyn
       do ipw1=1,npwdiel
        ar=(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
        ai=(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
        susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)&
&        +wght_dyn(1,ifreq)*ar-wght_dyn(2,ifreq)*ai
        susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)&
&        +wght_dyn(2,ifreq)*ar+wght_dyn(1,ifreq)*ai
       end do
      end do

!     MF  Update susmat_dyn by blas rank-1 updating--------------
!     Note: can use hermitian rank-1 updating if any of these apply
!     wfcts are real (istwf_k > 1)
!     system has inversion symmetry
!     frequency is zero
     else

!     MF  general rank-1 update by blas
#if defined T3E
      call CGERC(npwdiel,npwdiel,wght_dyn(1,ifreq),&
&      wfprod(1,1),1,&
&      wfprod(1,1),1,&
&      susmat_dyn(1,1,isp,1,isp,ifreq),npwdiel)
#else
      call ZGERC(npwdiel,npwdiel,wght_dyn(1,ifreq),&
&      wfprod(1,1),1,&
&      wfprod(1,1),1,&
&      susmat_dyn(1,1,isp,1,isp,ifreq),npwdiel)
#endif

!     End condition for updating susmat_dyn
     end if

!    MF  Not needed as occopt>=3 not implemented
!    if( occopt>=3 .and. testocc==1) then
!    Accumulate product of band densities by their doccde, for the
!    computation of the effect of change of Fermi level.
!    do ipw=1,npwdiel
!    drhode(1,ipw,isp)=drhode(1,ipw,isp)+wfprod(1,ipw)*wght1
!    drhode(2,ipw,isp)=drhode(2,ipw,isp)+wfprod(2,ipw)*wght1
!    end do
!    Also accumulate weighted sum of doccde
!    sumdocc=sumdocc+wght1
!    end if
!    MF

!    End frequency loop
    end do
!   End condition susopt
   end if
!  End condition of metallic occupancies or extrapolation or susopt
  end if
! End loop on iband
 end do

 call timab(87,2,tsec)

!DEBUG
!write(6,*)' %susk_dyn: END get real space wavefunctions'
!write(6,*)' susk_dyn : stop '
!stop
!ENDDEBUG

!--Wavefunctions have been generated in real space--------------------------

 call timab(88,1,tsec)

!Compute product of wavefunctions for different bands
 tolocc=1.0d-3
 if(nband_k>1)then
  do ibd1=1,nband_k-1
   do ibd2=ibd1+1,nband_k
!   If the occupation numbers are sufficiently different, or
!   if extrapolation is used and the corresponding factor is not zero,
!   then there is a contribution
    occdiff=occ(ibd1+bdtot_index)-occ(ibd2+bdtot_index)
    if( abs(occdiff)>tolocc      .or. &
&    ( extrap==1 .and.            &
&    ( abs(occ_deavg(ibd1)) + abs(occ_deavg(ibd2)) ) >tolocc ) &
&    .or. susopt==2 ) then

     eigdiff=eigen(ibd1+bdtot_index)-eigen(ibd2+bdtot_index)

!    DEBUG
!    write(6,*)' susk : contribution from bands',ibd1,ibd2
!    write(6,*)'   occ diff =',occdiff
!    write(6,*)'   eig diff =',eigdiff
!    ENDDEBUG

!    Store the contribution in wfraug
     do i3=1,ndiel3
      do i2=1,ndiel2
       do i1=1,ndiel1
        wfraug(1,i1,i2,i3)=wfrspa(1,i1,i2,i3,ibd1)*wfrspa(1,i1,i2,i3,ibd2)&
&        +wfrspa(2,i1,i2,i3,ibd1)*wfrspa(2,i1,i2,i3,ibd2)
        wfraug(2,i1,i2,i3)=wfrspa(2,i1,i2,i3,ibd1)*wfrspa(1,i1,i2,i3,ibd2)&
&        -wfrspa(1,i1,i2,i3,ibd1)*wfrspa(2,i1,i2,i3,ibd2)
       end do
      end do
     end do

!    DEBUG
!    norm=0.0_dp ; normr=0.0_dp
!    do i3=1,ndiel3
!    do i2=1,ndiel2
!    do i1=1,ndiel1
!    norm=norm+wfraug(1,i1,i2,i3)**2+wfraug(2,i1,i2,i3)**2
!    normr=normr+wfraug(1,i1,i2,i3)**2
!    end do
!    end do
!    end do
!    write(6,*)' norm in real space =',norm/dble(nfftdiel)
!    write(6,*)' norm of real part  =',normr/dble(nfftdiel)
!    ENDDEBUG

!    Performs the Fourier Transform of the product, and store it in wfprod
     tim_fourwf=11
     call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&     istwf_k,kg_diel,kg_diel,&
&     mgfftdiel,mpi_enreg,1,ngfftdiel,1,npwdiel,ndiel4,ndiel5,ndiel6,3,dtset%paral_kgb,tim_fourwf,weight)

!    DEBUG
!    write(6,*)' weight =',weight
!    norm=0.0_dp
!    do ipw=1,npwdiel
!    norm=norm+wfprod(1,ipw)**2+wfprod(2,ipw)**2
!    end do
!    write(6,*)' norm in reciprocal space  =',norm
!    ENDDEBUG

!    For density matrix related case
     if(susopt==2) then

!     DEBUG
!     write(6,*) ' %susk_dyn: off-diagonal part ibd1, ibd2=',ibd1,ibd2
!     ENDDEBUG

      ifreq=1
      weight=4._dp*wtk(ikpt)/ucvol
      do ipw2=1,npwdiel
       do ipw1=ipw2,npwdiel
        susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)+&
&        weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
        susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)+&
&        weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
       end do
      end do

!     For susceptibility matrix case
     else
!     Loop over frequencies
      do ifreq=1,nfreq

!      Step 1 - determine the weight for each frequency
       wght_dyn(:,ifreq)=0.0_dp
       if(abs(occdiff)>tolocc) then
        wght1=occdiff/max( 1.e-20_dp,freq(ifreq)**2+eigdiff**2 ) * 2.0_dp*wtk(ikpt)/ucvol
!       MF note: eigdiff is *minus* the difference (e_unoccupied - e_occupied)
        wght_dyn(1,ifreq)=eigdiff*wght1
        wght_dyn(2,ifreq)=-freq(ifreq)*wght1
       else
        wght1=0.0_dp
       end if

!      DEBUG
!      write(6,*) ' %susk_dyn: show off-diagonal weights:'
!      write(6,*) '  ibd1,ibd2            :',ibd1,ibd2
!      write(6,*) '  occdiff,eigdiff      :',occdiff,eigdiff
!      write(6,*) '  wght_dyn(1)(2)       :',wght_dyn(:,ifreq)
!      ENDDEBUG

       if(extrap==1)then
        wght_dyn(1,ifreq)=wght_dyn(1,ifreq)- 2.0_dp*wtk(ikpt)/ucvol &
&        *(occ_freq(1,ibd1,ifreq)+occ_freq(1,ibd2,ifreq))
        wght_dyn(2,ifreq)=wght_dyn(2,ifreq)- 2.0_dp*wtk(ikpt)/ucvol &
&        *(occ_freq(2,ibd1,ifreq)+occ_freq(2,ibd2,ifreq))
        weight=(occ_deavg(ibd1)+occ_deavg(ibd2)) * 2.0_dp*wtk(ikpt)/ucvol !MF!old

!       DEBUG
!       write(6,*) '  ibd1: occ_freq(1)(2):',occ_freq(:,ibd1,ifreq)*2.0*wtk(ikpt)/ucvol
!       write(6,*) '  ibd2: occ_freq(1)(2):',occ_freq(:,ibd2,ifreq)*2.0*wtk(ikpt)/ucvol
!       write(6,*) '  wght_dyn(1)(2)      :',wght_dyn(:,ifreq), 'done'
!       ENDDEBUG

       end if

!      Step 2 - sum contribution for each frequency
!      Update susmat_dyn explicitely-------------------------
       if(.false.) then
        do ipw2=1,npwdiel
!        Only fills lower half of the matrix (here, the susceptibility matrix)
!        Note that wfprod of the first index must behave like a density,
!        so that it is used as generated by fourwf, while wfprod of the
!        second index will be implicitely used to make a scalar product
!        with a potential change, meaning that its complex conjugate must be
!        used. This explains the following signs...

!        MF      fill full matrix explicite to check time reversal symmetry in suscep_dyn
         do ipw1=1,npwdiel
          ar=(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
          ai=(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
          susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)&
&          +wght_dyn(1,ifreq)*ar-wght_dyn(2,ifreq)*ai
          susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)&
&          +wght_dyn(2,ifreq)*ar+wght_dyn(1,ifreq)*ai
         end do
        end do

!       Update susmat_dyn by blas rank-1 updating--------------
       else

!       MF  general rank-1 update by blas
#if defined T3E
        call CGERC(npwdiel,npwdiel,wght_dyn(1,ifreq),&
&        wfprod(1,1),1,&
&        wfprod(1,1),1,&
&        susmat_dyn(1,1,isp,1,isp,ifreq),npwdiel)
#else
        call ZGERC(npwdiel,npwdiel,wght_dyn(1,ifreq),&
&        wfprod(1,1),1,&
&        wfprod(1,1),1,&
&        susmat_dyn(1,1,isp,1,isp,ifreq),npwdiel)
#endif

!       End condition for updating susmat_dyn
       end if

!      End frequency loop
      end do
!     End condition susopt
     end if

!    MF begin time reversed part-------------------------------------------
     if(.false.) then
!     Store the contribution in wfraug
      do i3=1,ndiel3
       do i2=1,ndiel2
        do i1=1,ndiel1
         wfraug(1,i1,i2,i3)=wfrspa(1,i1,i2,i3,ibd1)*wfrspa(1,i1,i2,i3,ibd2)&
&         +wfrspa(2,i1,i2,i3,ibd1)*wfrspa(2,i1,i2,i3,ibd2)
         wfraug(2,i1,i2,i3)=-wfrspa(2,i1,i2,i3,ibd1)*wfrspa(1,i1,i2,i3,ibd2)&
&         +wfrspa(1,i1,i2,i3,ibd1)*wfrspa(2,i1,i2,i3,ibd2)
        end do
       end do
      end do

!     DEBUG
!     norm=0.0_dp ; normr=0.0_dp
!     do i3=1,ndiel3
!     do i2=1,ndiel2
!     do i1=1,ndiel1
!     norm=norm+wfraug(1,i1,i2,i3)**2+wfraug(2,i1,i2,i3)**2
!     normr=normr+wfraug(1,i1,i2,i3)**2
!     end do
!     end do
!     end do
!     write(6,*)' norm in real space =',norm/dble(nfftdiel)
!     write(6,*)' norm of real part  =',normr/dble(nfftdiel)
!     ENDDEBUG

!     Performs the Fourier Transform of the product, and store it in wfprod
      tim_fourwf=11
      call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&      istwf_k,kg_diel,kg_diel,&
&      mgfftdiel,mpi_enreg,1,ngfftdiel,1,npwdiel,ndiel4,ndiel5,ndiel6,3,dtset%paral_kgb,tim_fourwf,weight)

!     For density matrix related case
      if(susopt==2) then

!      DEBUG
!      write(6,*) ' %susk_dyn: off-diagonal part ibd1, ibd2=',ibd1,ibd2
!      ENDDEBUG

       ifreq=1
       weight=4._dp*wtk(ikpt)/ucvol
       do ipw2=1,npwdiel
        do ipw1=ipw2,npwdiel
         susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)+&
&         weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
         susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)+&
&         weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
        end do
       end do

!      For susceptibility matrix case
      else
!      Loop over frequencies
       do ifreq=1,nfreq

!       Step 1 - determine the weight for each frequency
        wght_dyn(:,ifreq)=0.0_dp
        eigdiff=0._dp
        if(abs(occdiff)>tolocc) then
         wght1=occdiff/max( 1.e-20_dp,freq(ifreq)**2+eigdiff**2 ) * 2.0_dp*wtk(ikpt)/ucvol
!        MF note: eigdiff is *minus* the difference (e_unoccupied - e_occupied)
         wght_dyn(1,ifreq)=eigdiff*wght1
         wght_dyn(2,ifreq)=freq(ifreq)*wght1
        else
         wght1=0.0_dp
        end if

!       DEBUG
!       write(6,*) ' %susk_dyn: show off-diagonal weights:'
!       write(6,*) '  ibd1,ibd2            :',ibd1,ibd2
!       write(6,*) '  occdiff,eigdiff      :',occdiff,eigdiff
!       write(6,*) '  wght_dyn(1)(2)       :',wght_dyn(:,ifreq)
!       ENDDEBUG

        if(extrap==1)then
         wght_dyn(1,ifreq)=wght_dyn(1,ifreq)- 2.0_dp*wtk(ikpt)/ucvol &
&         *(occ_freq(1,ibd1,ifreq)+occ_freq(1,ibd2,ifreq))
         wght_dyn(2,ifreq)=wght_dyn(2,ifreq)- 2.0_dp*wtk(ikpt)/ucvol &
&         *(occ_freq(2,ibd1,ifreq)+occ_freq(2,ibd2,ifreq))
         weight=(occ_deavg(ibd1)+occ_deavg(ibd2)) * 2.0_dp*wtk(ikpt)/ucvol !MF!old

!        DEBUG
!        write(6,*) '  ibd1: occ_freq(1)(2):',occ_freq(:,ibd1,ifreq)*2.0*wtk(ikpt)/ucvol
!        write(6,*) '  ibd2: occ_freq(1)(2):',occ_freq(:,ibd2,ifreq)*2.0*wtk(ikpt)/ucvol
!        write(6,*) '  wght_dyn(1)(2)      :',wght_dyn(:,ifreq), 'done'
!        ENDDEBUG

        end if

!       Step 2 - sum contribution for each frequency
!       Update susmat_dyn explicitely-------------------------
        if(.true.) then
         do ipw2=1,npwdiel
!         Only fills lower half of the matrix (here, the susceptibility matrix)
!         Note that wfprod of the first index must behave like a density,
!         so that it is used as generated by fourwf, while wfprod of the
!         second index will be implicitely used to make a scalar product
!         with a potential change, meaning that its complex conjugate must be
!         used. This explains the following signs...

!         MF      fill complete matrix
          do ipw1=1,npwdiel
           ar=(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
           ai=(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
           susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)&
&           +wght_dyn(1,ifreq)*ar-wght_dyn(2,ifreq)*ai
           susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)&
&           +wght_dyn(2,ifreq)*ar+wght_dyn(1,ifreq)*ai
          end do
         end do

!        Update susmat_dyn by blas rank-1 updating--------------
        else

#if defined T3E
         call CHER('l',npwdiel,wght_dyn(1,ifreq),&
&         wfprod(1,1),1,&
&         susmat_dyn(1,1,isp,1,isp,ifreq),npwdiel)
#else
         call ZHER('l',npwdiel,wght_dyn(1,ifreq),&
&         wfprod(1,1),1,&
&         susmat_dyn(1,1,isp,1,isp,ifreq),npwdiel)
#endif

!        End condition for updating susmat_dyn
        end if

!       End frequency loop
       end do
!      End condition susopt
      end if

!     End condition false
     end if
!    MF end time reversed part---------------------------------------------

!    End condition of different occupation numbers or extrapolation
    end if
!   End internal loop over bands
   end do
!  End external loop over bands
  end do
! End condition of having more than one band
 end if

!DEBUG
!write(6,*)' %susk_dyn: ! susmat_dyn(): diagonal elements'
!do ifreq=1,nfreq
!write(6,*) 'ifreq:',ifreq
!do ipw1=1,npwdiel
!write(6,*) ipw1,susmat_dyn(1,ipw1,1,ipw1,1,ifreq),susmat_dyn(2,ipw1,1,ipw1,1,ifreq)
!end do
!end do
!write(6,*)' %susk_dyn : DONE & EXIT'; write(6,*)
!ENDDEBUG

 deallocate(cwavef)
 deallocate(dummy)
 deallocate(rhoaug)
 deallocate(wfprod)
 deallocate(wfraug)
 deallocate(wfrspa)
 deallocate(wght_dyn)

 call timab(88,2,tsec)

end subroutine susk_dyn
!!***
