!{\src2tex{textfont=tt}}
!!****f* ABINIT/suskmm
!! NAME
!! suskmm
!!
!! FUNCTION
!! Compute the contribution of one k point to the susceptibility matrix
!! from input wavefunctions, band occupations, and k point wts.
!! Include the usual sum-over-state terms, but also the
!! corrections due to the change of the Fermi level in the metallic
!! case, as well as implicit sum over higher lying conduction
!! states, thanks to the closure relation (referred to as an extrapolation).
!!
!! This routine is similar to susk, but use blocking on wavefunctions
!! to decrease memory requirements, at the expense of CPU time.
!!
!! NOTES
!! There is still room for optimization !!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  bdtot_index=index for the number of the band
!!  cg(2,mcg)=wf in G space
!!  cprj_k(natom,nspinor*nband_k)= wave functions projected with non-local projectors:
!!                                 cprj_k=<p_i|Cnk> where p_i is a non-local projector.
!!  doccde(mband*nkpt*nsppol)=derivative of occupancies wrt
!!           the energy for each band and k point
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  extrap: if==1, the closure relation (an extrapolation) must be used
!!  gbound(2*mgfftdiel+8,2)=G sphere boundary for going from WF sphere to
!!      medium size FFT grid
!!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for going from medium size
!!      FFT grid to small sphere.
!!  gylmg_diel(npwdiel,lmax_diel**2,ntypat*usepaw)= -PAW only- Fourier transform of g_l(r).Y_ml(r) shape functions
!!                                                   for dielectric matrix
!!  icg=index for cg
!!  ikpt=number of the k point
!!  isp=number of the current spin
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kg_k(3,npw)=coordinates of planewaves in basis sphere.
!!  lmax_diel=1+max. value of l angular momentum used for dielectric matrix
!!  mband=maximum number of bands
!!  mcg=dimension of cg
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of
!!     the dielectric matrix
!!  mkmem=maximum number of k points in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum allowed value for npw
!!  natom=number of atoms in cell
!!  nband_k=number of bands at this k point for that spin polarization
!!  ndiel4,ndiel5,ndiel6= FFT dimensions, modified to avoid cache trashing
!!  nfftdiel=number of fft grid points for the computation of the diel matrix
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  npwdiel=third and fifth dimension of the susmat array.
!!  npw_k=number of plane waves at this k point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  occ(mband*nkpt*nsppol)=
!!          occupation numbers for each band (usually 2.0) at each k point
!!  occopt=option for occupancies
!!  occ_deavg(mband)=factor for extrapolation (occup. divided by an energy gap)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph3d_diel(2,npwdiel,natom*usepaw)=3-dim structure factors, for each atom and plane wave, for dielectric matrix
!!  typat(natom)=type (integer) for each atom
!!  ucvol=unit cell volume (Bohr**3)
!!  usepaw=flag for PAW
!!  wtk(nkpt)=k point weights (they sum to 1.0)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  drhode(2,npwdiel,nsppol)=weighted density, needed to compute the
!!   effect of change of fermi energy
!!  rhoextrap(ndiel4,ndiel5,ndiel6)=density-like array, needed for the
!!   extrapolation procedure.
!!  sumdocc=sum of weighted occupation numbers, needed to compute the
!!   effect of change of fermi energy
!!  susmat(2,npwdiel,nsppol,npwdiel,nsppol)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!
!! PARENTS
!!      suscep_stat
!!
!! CHILDREN
!!      fourwf,leave_new,pawsushat,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine suskmm(atindx1,bdtot_index,cg,cprj_k,doccde,drhode,eigen,extrap,gbound,&
&  gbound_diel,gylmg_diel,icg,ikpt,isp,istwfk,kg_diel,kg_k,&
&  lmax_diel,mband,mcg,mgfftdiel,mkmem,mpi_enreg,mpw,&
&  natom,nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,ngfftdiel,nkpt,&
&  npwdiel,npw_k,nspden,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,paral_kgb,&
&  pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&
&  susmat,typat,ucvol,usepaw,wtk)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_13paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!This type is defined in defs_mpi
!scalars
 integer,intent(in) :: bdtot_index,extrap,icg,ikpt,isp,lmax_diel,mband,mcg
 integer,intent(in) :: mgfftdiel,mkmem,mpw,natom,nband_k,ndiel4,ndiel5,ndiel6
 integer,intent(in) :: nfftdiel,nkpt,npw_k,npwdiel,nspden,nspinor,nsppol,ntypat
 integer,intent(in) :: occopt,paral_kgb,usepaw
 real(dp),intent(in) :: ucvol
 real(dp),intent(inout) :: sumdocc
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: atindx1(natom),gbound(2*mgfftdiel+8,2)
 integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2),istwfk(nkpt)
 integer,intent(in) :: kg_diel(3,npwdiel),kg_k(3,npw_k),ngfftdiel(18)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: cg(2,mcg),doccde(mband*nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: gylmg_diel(npwdiel,lmax_diel**2,ntypat*usepaw)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),occ_deavg(mband)
 real(dp),intent(in) :: ph3d_diel(2,npwdiel,natom*usepaw),wtk(nkpt)
 real(dp),intent(inout) :: drhode(2,npwdiel,nsppol)
 real(dp),intent(inout) :: rhoextrap(ndiel4,ndiel5,ndiel6)
 real(dp),intent(inout) :: susmat(2,npwdiel,nsppol,npwdiel,nsppol)
 type(cprj_type) :: cprj_k(natom,nspinor*nband_k*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
! real(dp), allocatable :: cg_disk(:,:)
!scalars
 integer :: i1,i2,i3,iband,iband_shift,iband_shift2,ibd1,ibd2,ibdshft1,ibdshft2
 integer :: iblk1,iblk2,ipw,ipw1,ipw2,isp1,isp2,istwf_k,mblk,nblk,nbnd_current
 integer :: nbnd_in_blk,nbnd_in_blk1,ndiel1,ndiel2,ndiel3,testocc,tim_fourwf
 real(dp) :: ai,ar,eigdiff,norm,normr,occdiff,tolocc,weight,wght1,wght2
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwavef(:,:),dummy(:,:),rhoaug(:,:,:),wfprod(:,:)
 real(dp),allocatable :: wfraug(:,:,:,:),wfrspa1(:,:,:,:,:),wfrspa2(:,:,:,:,:)

! *************************************************************************

!DEBUG
!write(6,*)' suskmm : enter '
!if(.true.)stop
!ENDDEBUG

 ndiel1=ngfftdiel(1) ; ndiel2=ngfftdiel(2) ; ndiel3=ngfftdiel(3)
 istwf_k=istwfk(1)

 testocc=1

 allocate(cwavef(2,mpw),dummy(2,1))
 allocate(rhoaug(ndiel4,ndiel5,ndiel6),wfraug(2,ndiel4,ndiel5,ndiel6))
 allocate(wfprod(2,npwdiel))

!Prepare the blocking : compute the number of blocks,
!the number of bands in each normal block,
!and the number in the first one, usually smaller.

!Consider that if the number of bands is large, there are at most 8 blocks
 if(nband_k>=48)then
  mblk=8
  nbnd_in_blk=(nband_k-1)/mblk+1
! If the number of bands is medium, place 6 bands per block
 else if(nband_k>=12)then
  nbnd_in_blk=6
! Otherwise, must have at least 2 blocks
 else if(nband_k>=2)then
  mblk=2
  nbnd_in_blk=(nband_k-1)/mblk+1
 else
  write(message, '(a,a,a,a,a,a,i2,a,a,a)') ch10,&
&  ' suskmm : ERROR -',ch10,&
&  '  The number of bands must be larger or equal to 2, in suskmm.',ch10,&
&  '  It is equal to ',nband_k,'.',ch10,&
&  '  Action : choose another preconditioner.'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if

!Compute the effective number of blocks, and the number of bands in
!the first block.
 nblk=(nband_k-1)/nbnd_in_blk+1
 nbnd_in_blk1=nband_k-(nblk-1)*nbnd_in_blk

!DEBUG
!write(6,*)' suskmm : nband_k,nblk,nbnd_in_blk,nbnd_in_blk1 '
!write(6,*)nband_k,nblk,nbnd_in_blk,nbnd_in_blk1
!stop
!ENDDEBUG

!wfrspa1 will contain the wavefunctions of the slow sampling (iblk1)
 allocate(wfrspa1(2,ndiel4,ndiel5,ndiel6,nbnd_in_blk))
!wfrspa2 will contain the wavefunctions of the rapid sampling (iblk2)
 allocate(wfrspa2(2,ndiel4,ndiel5,ndiel6,nbnd_in_blk))

!First loop over blocks
 do iblk1=1,nblk

  call timab(87,1,tsec)

! Initialisation
  if(iblk1==1)then

   nbnd_current=nbnd_in_blk1
   iband_shift=0
!  Loop over bands to fft and store Fourier transform of wavefunction
   do iband=1,nbnd_current
!   Obtain Fourier transform in fft box
    cwavef(:,1:npw_k)=cg(:,1+(iband-1)*npw_k+icg:iband*npw_k+icg)
    tim_fourwf=12
    call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&    istwf_k,kg_k,kg_k,mgfftdiel,mpi_enreg,1,ngfftdiel,npw_k,1,ndiel4,ndiel5,ndiel6,&
&    0,paral_kgb,tim_fourwf,weight)
    wfrspa1(:,:,:,:,iband)=wfraug(:,:,:,:)
   end do

  else

!  The Fourier transform of wavefunctions have already been obtained
   nbnd_current=nbnd_in_blk
   iband_shift=nbnd_in_blk1+(iblk1-2)*nbnd_in_blk

  end if

! Loop over bands of this block, to generate band-diagonal
! contributions to sumdocc, drhode, rhoextrap, and susmat.

! DEBUG
! write(6,*)' suskmm : 1'
! ENDDEBUG

  do iband=1,nbnd_current

   if( (occopt>=3 .and. testocc==1) .or. extrap==1 )then
!   In the case of metallic occupation, or if the extrapolation
!   over higher bands is included, must compute the
!   Fourier transform of the density of each band, then
!   generate the part of the susceptibility matrix due
!   varying occupation numbers.
    weight=-two*occ_deavg(iband+iband_shift)*wtk(ikpt)/ucvol
    do i3=1,ndiel3
     do i2=1,ndiel2
      do i1=1,ndiel1
       wfraug(1,i1,i2,i3)=wfrspa1(1,i1,i2,i3,iband)**2&
&       +wfrspa1(2,i1,i2,i3,iband)**2
       wfraug(2,i1,i2,i3)=zero
      end do
     end do
!    If extrapolation, accumulate density in real space
     if(extrap==1.and.usepaw==0)then
      do i2=1,ndiel2
       do i1=1,ndiel1
        rhoextrap(i1,i2,i3)=rhoextrap(i1,i2,i3)+weight*wfraug(1,i1,i2,i3)
       end do
      end do
     end if
    end do

!   In case of PAW, add compensation charge contribution
    if (usepaw==1.and.extrap==1) then
     call pawsushat(atindx1,cprj_k,gbound_diel,gylmg_diel,iband,iband,istwf_k,kg_diel,&
&     lmax_diel,mgfftdiel,mpi_enreg,natom,nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,&
&     ngfftdiel,npwdiel,nspinor,ntypat,1,paral_kgb,&
&     pawang,pawtab,ph3d_diel,typat,dummy,wfraug)
     rhoextrap(:,:,:)=rhoextrap(:,:,:)+weight*wfraug(1,:,:,:)
    end if

!   Performs the Fourier Transform of the density of the band,
!   and store it in wfprod
    tim_fourwf=13
    call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&    istwf_k,kg_diel,kg_diel,&
&    mgfftdiel,mpi_enreg,1,ngfftdiel,1,npwdiel,ndiel4,ndiel5,ndiel6,3,paral_kgb,tim_fourwf,weight)
!   In case of PAW, add compensation charge contribution if not already done
    if (usepaw==1.and.extrap==0) then
     call pawsushat(atindx1,cprj_k,gbound_diel,gylmg_diel,iband,iband,istwf_k,kg_diel,&
&     lmax_diel,mgfftdiel,mpi_enreg,natom,nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,&
&     ngfftdiel,npwdiel,nspinor,ntypat,0,paral_kgb,&
&     pawang,pawtab,ph3d_diel,typat,wfprod,dummy)
    end if

!   Perform now the summation of terms related to direct change of eigenvalues
!   or extrapolation over higher bands
    wght1=zero ; wght2=zero
    if(occopt>=3 .and. testocc==1)then
     wght1=doccde(iband+iband_shift+bdtot_index)*wtk(ikpt)/ucvol
    end if

    if(extrap==1) wght2=two*occ_deavg(iband+iband_shift)*wtk(ikpt)/ucvol

    weight=wght1+wght2
    do ipw2=1,npwdiel
!    Only fills lower half of the matrix (here, the susceptibility matrix)
!    Note that wfprod of the first index must behave like a density,
!    so that it is used as generated by fourwf, while wfprod of the
!    second index will be implicitely used to make a scalar product
!    with a potential change, meaning that its complex conjugate must be
!    used. This explains the following signs...
     do ipw1=ipw2,npwdiel
      susmat(1,ipw1,isp,ipw2,isp)=susmat(1,ipw1,isp,ipw2,isp)+&
&      weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
      susmat(2,ipw1,isp,ipw2,isp)=susmat(2,ipw1,isp,ipw2,isp)+&
&      weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
     end do
    end do
    if( occopt>=3 .and. testocc==1) then
!    Accumulate product of band densities by their doccde, for the
!    computation of the effect of change of Fermi level.
     do ipw=1,npwdiel
      drhode(1,ipw,isp)=drhode(1,ipw,isp)+wfprod(1,ipw)*wght1
      drhode(2,ipw,isp)=drhode(2,ipw,isp)+wfprod(2,ipw)*wght1
     end do
!    Also accumulate weighted sum of doccde
     sumdocc=sumdocc+wght1
    end if

!   End condition of metallic occupancies or extrapolation
   end if

!  End loop on iband
  end do

! DEBUG
! write(6,*)' suskmm : 2'
! ENDDEBUG

  call timab(87,2,tsec)

! -- Compute now off-band-diagonal terms ------------------------------------

  call timab(88,1,tsec)

  tolocc=1.0d-3

! Compute product of wavefunctions for different bands, inside the block
  if(nbnd_current>1)then
   do ibd1=1,nbnd_current-1
    ibdshft1=ibd1+iband_shift
    do ibd2=ibd1+1,nbnd_current
     ibdshft2=ibd2+iband_shift

!    If the occupation numbers are sufficiently different, or
!    if extrapolation is used and the corresponding factor is not zero,
!    then there is a contribution
     occdiff=occ(ibdshft1+bdtot_index)-occ(ibdshft2+bdtot_index)
     if( abs(occdiff)>tolocc      .or. &
&     ( extrap==1 .and.            &
&     ( abs(occ_deavg(ibdshft1)) + abs(occ_deavg(ibdshft2)) ) >tolocc ) &
&     ) then

      eigdiff=eigen(ibdshft1+bdtot_index) - eigen(ibdshft2+bdtot_index)

!     Store the contribution in wfraug
      do i3=1,ndiel3
       do i2=1,ndiel2
        do i1=1,ndiel1
         wfraug(1,i1,i2,i3)=wfrspa1(1,i1,i2,i3,ibd1)*wfrspa1(1,i1,i2,i3,ibd2)&
&         +wfrspa1(2,i1,i2,i3,ibd1)*wfrspa1(2,i1,i2,i3,ibd2)
         wfraug(2,i1,i2,i3)=wfrspa1(2,i1,i2,i3,ibd1)*wfrspa1(1,i1,i2,i3,ibd2)&
&         -wfrspa1(1,i1,i2,i3,ibd1)*wfrspa1(2,i1,i2,i3,ibd2)
        end do
       end do
      end do

!     Performs the Fourier Transform of the product, and store it in wfprod
      tim_fourwf=13
      call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&      istwf_k,kg_diel,kg_diel,&
&      mgfftdiel,mpi_enreg,1,ngfftdiel,1,npwdiel,ndiel4,ndiel5,ndiel6,3,paral_kgb,tim_fourwf,weight)

!     In case of PAW, add compensation charge contribution
      if (usepaw==1) then
       call pawsushat(atindx1,cprj_k,gbound_diel,gylmg_diel,ibd1,ibd2,istwf_k,kg_diel,&
&       lmax_diel,mgfftdiel,mpi_enreg,natom,nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,&
&       ngfftdiel,npwdiel,nspinor,ntypat,0,paral_kgb,&
&       pawang,pawtab,ph3d_diel,typat,wfprod,dummy)
      end if

!     Perform now the summation
      wght1=zero ; wght2=zero
      if(abs(occdiff)>tolocc)wght1= occdiff/eigdiff * two*wtk(ikpt)/ucvol
      if(extrap==1)then
       wght2=(occ_deavg(ibdshft1)+occ_deavg(ibdshft2)) * two*wtk(ikpt)/ucvol
      end if
      weight=wght1+wght2

      do ipw2=1,npwdiel
!      Only fills lower half of the matrix (here, the susceptibility matrix)
!      Note that wfprod of the first index must behave like a density,
!      so that it is used as generated by fourwf, while wfprod of the
!      second index will be implicitely used to make a scalar product
!      with a potential change, meaning that its complex conjugate must be
!      used. This explains the following signs...
       do ipw1=ipw2,npwdiel
        susmat(1,ipw1,isp,ipw2,isp)=susmat(1,ipw1,isp,ipw2,isp)+&
&        weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
        susmat(2,ipw1,isp,ipw2,isp)=susmat(2,ipw1,isp,ipw2,isp)+&
&        weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
       end do
      end do

!     End condition of different occupation numbers or extrapolation
     end if
!    End internal loop over bands
    end do
!   End external loop over bands
   end do
!  End condition of having more than one band
  end if

! DEBUG
! write(6,*)' suskmm : 3'
! ENDDEBUG

! Loop on secondary block, with fast varying index, in decreasing order.
  if(iblk1/=nblk)then
   do iblk2=nblk,iblk1+1,-1
    iband_shift2=nbnd_in_blk1+(iblk2-2)*nbnd_in_blk

!   Loop over bands to fft and store Fourier transform of wavefunction
    iband_shift2=nbnd_in_blk1+(iblk2-2)*nbnd_in_blk
    do iband=1,nbnd_in_blk
!    Obtain Fourier transform in fft box
     cwavef(:,1:npw_k)= &
&     cg(:,1+(iband+iband_shift2-1)*npw_k+icg:(iband+iband_shift2)*npw_k+icg)
     tim_fourwf=12
     call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&     istwf_k,kg_k,kg_k,mgfftdiel,mpi_enreg,1,ngfftdiel,npw_k,1,&
&     ndiel4,ndiel5,ndiel6,0,paral_kgb,tim_fourwf,weight)
     wfrspa2(:,:,:,:,iband)=wfraug(:,:,:,:)
    end do

    do ibd1=1,nbnd_current
     ibdshft1=ibd1+iband_shift
     do ibd2=1,nbnd_in_blk
      ibdshft2=ibd2+iband_shift2

!     If the occupation numbers are sufficiently different, or
!     if extrapolation is used and the corresponding factor is not zero,
!     then there is a contribution
      occdiff=occ(ibdshft1+bdtot_index)-occ(ibdshft2+bdtot_index)
      if( abs(occdiff)>tolocc      .or. &
&      ( extrap==1 .and.            &
&      ( abs(occ_deavg(ibdshft1)) + abs(occ_deavg(ibdshft2)) ) >tolocc ) &
&      ) then

       eigdiff=eigen(ibdshft1+bdtot_index) - eigen(ibdshft2+bdtot_index)

!      Store the contribution in wfraug
       do i3=1,ndiel3
        do i2=1,ndiel2
         do i1=1,ndiel1
          wfraug(1,i1,i2,i3)=wfrspa1(1,i1,i2,i3,ibd1)*wfrspa2(1,i1,i2,i3,ibd2)&
&          +wfrspa1(2,i1,i2,i3,ibd1)*wfrspa2(2,i1,i2,i3,ibd2)
          wfraug(2,i1,i2,i3)=wfrspa1(2,i1,i2,i3,ibd1)*wfrspa2(1,i1,i2,i3,ibd2)&
&          -wfrspa1(1,i1,i2,i3,ibd1)*wfrspa2(2,i1,i2,i3,ibd2)
         end do
        end do
       end do

!      Performs the Fourier Transform of the product, and store it in wfprod
       tim_fourwf=13
       call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&       istwf_k,kg_diel,kg_diel,&
&       mgfftdiel,mpi_enreg,1,ngfftdiel,1,npwdiel,ndiel4,ndiel5,ndiel6,3,paral_kgb,tim_fourwf,weight)

!      In case of PAW, add compensation charge contribution
       if (usepaw==1) then
        call pawsushat(atindx1,cprj_k,gbound_diel,gylmg_diel,ibd1,ibdshft2,istwf_k,kg_diel,&
&        lmax_diel,mgfftdiel,mpi_enreg,natom,nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,&
&        ngfftdiel,npwdiel,nspinor,ntypat,0,paral_kgb,&
&        pawang,pawtab,ph3d_diel,typat,wfprod,dummy)
       end if

!      Perform now the summation
       wght1=zero ; wght2=zero
       if(abs(occdiff)>tolocc)wght1= occdiff/eigdiff * two*wtk(ikpt)/ucvol
       if(extrap==1)then
        wght2=(occ_deavg(ibdshft1)+occ_deavg(ibdshft2)) * two*wtk(ikpt)/ucvol
       end if
       weight=wght1+wght2

       do ipw2=1,npwdiel
!       Only fills lower half of the matrix (here, the susceptibility matrix)
!       Note that wfprod of the first index must behave like a density,
!       so that it is used as generated by fourwf, while wfprod of the
!       second index will be implicitely used to make a scalar product
!       with a potential change, meaning that its complex conjugate must be
!       used. This explains the following signs...
        do ipw1=ipw2,npwdiel
         susmat(1,ipw1,isp,ipw2,isp)=susmat(1,ipw1,isp,ipw2,isp)+&
&         weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
         susmat(2,ipw1,isp,ipw2,isp)=susmat(2,ipw1,isp,ipw2,isp)+&
&         weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
        end do
       end do

!      End condition of different occupation numbers or extrapolation
      end if
!     End internal loop over bands
     end do
!    End external loop over bands
    end do
!   End loop on bloks
   end do

!  Finish the loop on blok with iblk2=iblk1+1, so can use the
!  FFTd wavefunctions for the next iblk1.
   do iband=1,nbnd_in_blk
    wfrspa1(:,:,:,:,iband)=wfrspa2(:,:,:,:,iband)
   end do

!  End condition of iblk1/=nblk
  end if

  call timab(88,2,tsec)

! End loop on iblk1
 end do

!DEBUG
!write(6,*)' suskmm : exit '
!do ipw1=1,npwdiel
!write(6,*)ipw1,susmat(1,ipw1,1,ipw1,1),susmat(2,ipw1,1,ipw1,1)
!end do
!write(6,*)' suskmm : end of susmat '
!stop
!ENDDEBUG

 deallocate(cwavef,dummy,rhoaug,wfprod,wfraug,wfrspa1,wfrspa2)

end subroutine suskmm
!!***
