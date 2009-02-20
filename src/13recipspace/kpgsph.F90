!{\src2tex{textfont=tt}}
!!****f* ABINIT/kpgsph
!! NAME
!! kpgsph
!!
!! FUNCTION
!! Use reciprocal space metric gmet(3,3) to set up the list
!! of G vectors inside a sphere out to $ (1/2)*(2*\pi*(k+G))^2=ecut $.
!! If mkmem=0 and mpw=0, then only count the number of planewaves
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ecut=planewave kinetic energy cutoff (hartrees)
!!  exchn2n3d=if 1, n2 and n3 are exchanged
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  ikg=shift to be given to the location of the output data in the array kg
!!  ikpt=number of the k-point
!!  istwf_k=option parameter that describes the storage of wfs
!!  kpt(3)=reduced coords of k point (in terms of recip latt vecs)
!!  mkmem =maximum number of k points which can fit in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum number of planewaves as dimensioned in calling routine
!!
!! OUTPUT
!!  kg(3,mpw*mkmem)=dimensionless coords of resulting G vecs (integer)
!!  mpi_enreg%me_g0=if 1, the plane wave G(0 0 0) is in the set of plane waves (and is the first)
!!  npw=resulting number of planewaves inside ecut centered at kpt
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  Must take into account the time-reversal symmetry
!!  when istwf_k is not 1.
!!
!! PARENTS
!!      getmpw,initberry,kpgio,newsp,setshells,wfconv
!!
!! CHILDREN
!!      leave_new,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine kpgsph(ecut,exchn2n3d,gmet,ikg,ikpt,istwf_k,kg,kpt,mkmem,mpi_enreg,mpw,npw)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: exchn2n3d,ikg,ikpt,istwf_k,mkmem,mpw
 integer,intent(out) :: npw
 real(dp),intent(in) :: ecut
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(out) :: kg(3,mpw*mkmem)
 real(dp),intent(in) :: gmet(3,3),kpt(3)

!Local variables-------------------------------
!For the bandFFT case
!scalars
 integer :: i1,i2,ig,ig1p,ig1pmax,ig2,ig2p,ig2pmax,ig2pmin,ig3,ig3p,ig3pmax
 integer :: ig3pmin,ii,in,ind,n2,np_band,npw_remain,npw_split
 real(dp) :: gmet11,gmet_trace,gmin,gs,gs_fact,gs_part,gscut,v1,v2,v3,xx
 logical :: ipw_ok
 character(len=500) :: message
!arrays
 integer :: ngrid(3),nmax(3),nmin(3)
 integer,allocatable :: array_ipw(:),array_ipw_1_n2(:),ig1arr(:),ig2arr(:)
 integer,allocatable :: ig3arr(:),kg_small(:,:)
 real(dp) :: kmax(3),minor(3),numer(3),tsec(2)
 real(dp),allocatable :: kg1arr(:),kg2arr(:),kg3arr(:)

! *************************************************************************

!DEBUG
!write(6,*)' kpgsph : enter , exchn2n3d=',exchn2n3d
!write(6,*)' istwf_k ',istwf_k
!write(6,*)' gmet ',gmet
!write(6,*)' kpt ',kpt
!write(6,*)' ecut ',ecut
!write(6,*)' ikgi ',ikg
!write(6,*)' me_fft ',mpi_enreg%me_fft
!write(6,*)' mkmem ',mkmem
!write(6,*)' mpw ',mpw
!write(6,*)' nproc_fft',mpi_enreg%nproc_fft
!ENDDEBUG

 call timab(23,1,tsec)
 if(istwf_k<1 .or. istwf_k>9)then
  write(message, '(a,a,a,a,a,a,i4,a)' ) ch10,&
&  ' kpgsph : BUG -',ch10,&
&  '  The variable istwf_k must be between 1 and 9, while',ch10,&
&  '  the argument of the routine istwf_k =',istwf_k,'.'
  call wrtout(06,message,'PERS')
  call leave_new('PERS')
 end if

 if(ikg+mpw>mkmem*mpw)then
  write(message, '(a,a,a,a,a,a,a,a,i8,a,i8,a,i8,a,a,a,a)' ) ch10,&
&  ' kpgsph : BUG -',ch10,&
&  '  The variables ikg, mkmem, and mpw  must satisfy ikg<=(mkmem-1)*mpw,',ch10,&
&  '  while the arguments of the routine are',ch10,&
&  '  ikg =',ikg,', mkmem =',mkmem,', mpw =',mpw,ch10,&
&  '  Probable cause: Known error in invars1 for parallel spin-polarized case.',ch10,&
&  '  Temporary solution: Change the number of parallel processes.'
  call wrtout(06,message,'PERS')
  call leave_new('PERS')
 end if

!PATCH kpgsph init npband car mpw<=0
 if (mpw<=0) then
  np_band=0
 end if

 if(mpw > 0) then
  np_band=1
  if(mpi_enreg%mode_para=='b') then
   np_band=mpi_enreg%sizecart(2)
  end if
  allocate(kg_small(3,mpw*np_band))
 end if

!A larger array, that will be split on the correct processor
!G**2 cutoff, gscut=Ecut/2 /Pi^2

 gscut=0.5_dp*ecut*piinv**2

!In reduced coordinates, determine maximal value of k+G and G
!for each direction

 minor(1)=gmet(2,2)*gmet(3,3)-gmet(2,3)**2
 numer(1)=gmet(1,2)**2*gmet(3,3)-2.0_dp*gmet(1,2)*gmet(1,3)*gmet(2,3) &
& +gmet(1,3)**2*gmet(2,2)
 minor(2)=gmet(1,1)*gmet(3,3)-gmet(1,3)**2
 numer(2)=gmet(2,3)**2*gmet(1,1)-2.0_dp*gmet(1,2)*gmet(1,3)*gmet(2,3) &
& +gmet(2,1)**2*gmet(3,3)
 minor(3)=gmet(2,2)*gmet(1,1)-gmet(1,2)**2
 numer(3)=gmet(3,2)**2*gmet(1,1)-2.0_dp*gmet(1,2)*gmet(1,3)*gmet(2,3) &
& +gmet(1,3)**2*gmet(2,2)

!Take the trace of the gmet tensor as dimensional reference
 gmet_trace=gmet(1,1)+gmet(2,2)+gmet(3,3)

 do ii=1,3
  xx=gmet(ii,ii)*minor(ii)-numer(ii)
  if(xx<tol10*gmet_trace**3 .or. minor(ii)<tol10*gmet_trace**2)then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' kpgsph : BUG -',ch10,&
&   '  The metric tensor seem incorrect'
   call wrtout(06,message,'PERS')
   call leave_new('PERS')
  end if
  kmax(ii)=sqrt(gscut*minor(ii)/xx)
  nmax(ii)=floor(kmax(ii)-kpt(ii)+tol10)
  nmin(ii)=ceiling(-kmax(ii)-kpt(ii)-tol10)
  ngrid(ii)=nmax(ii)-nmin(ii)+1
 end do
!perform looping over fft box grid of size ngfft(1)*ngfft(2)*ngfft(3):
 ig=0;ind=0
 in=0
 gmet11=gmet(1,1)

!Set up standard search sequence for grid points, in standard storage mode :
!0 1 2 3 ... nmax nmin ... -1
!If the mode is not standard, then some part of the FFT grid must be selected
!
 allocate(ig1arr(ngrid(1)),ig2arr(ngrid(2)),ig3arr(ngrid(3)))
 allocate(kg1arr(ngrid(1)),kg2arr(ngrid(2)),kg3arr(ngrid(3)))

 do ig1p=1,ngrid(1)
  ig1arr(ig1p)=ig1p-1
  if(ig1p-1>nmax(1))ig1arr(ig1p)=ig1p-ngrid(1)-1
  kg1arr(ig1p)=kpt(1)+dble(ig1arr(ig1p))
 end do
!For the second direction, the number of points might depend on istwf_k
 ig2pmax=ngrid(2)
 if(istwf_k>=2 .and. exchn2n3d==0)ig2pmax=nmax(2)+1
 allocate(array_ipw(-ig2pmax:ig2pmax))
 array_ipw(:)=0
 do ig2p=1,ig2pmax
  ig2arr(ig2p)=ig2p-1
  if(ig2p-1>nmax(2))ig2arr(ig2p)=ig2p-ngrid(2)-1
  kg2arr(ig2p)=kpt(2)+dble(ig2arr(ig2p))
 end do
!For the third direction, the number of points might depend on istwf_k
 ig3pmax=ngrid(3)
 if(istwf_k>=2 .and. exchn2n3d==1)ig3pmax=nmax(3)+1

 do ig3p=1,ig3pmax
  ig3arr(ig3p)=ig3p-1
  if(ig3p-1>nmax(3))ig3arr(ig3p)=ig3p-ngrid(3)-1
  kg3arr(ig3p)=kpt(3)+dble(ig3arr(ig3p))
 end do

!Performs loop on all grid points.
 if(exchn2n3d==0)then
  mpi_enreg%me_g0=0

  do ig3p=1,ngrid(3)
   ig3=ig3arr(ig3p)
   v3=kg3arr(ig3p)
   ig2pmin=1
   if( istwf_k>=2 .and. istwf_k<=5 .and. ig3<0)then
    ig2pmin=2
   end if
!  ig2pmax was initialized previously
   do ig2p=ig2pmin,ig2pmax
    ig2=ig2arr(ig2p)
    ipw_ok=(istwf_k==1 .and. mpi_enreg%me_fft==modulo(ig2,mpi_enreg%nproc_fft) .or. &
&    (istwf_k>=2 .and. mpi_enreg%me_fft==modulo(ig2,mpi_enreg%nproc_fft)))
!   PAY ATTENTION : old if was ipw_ok=(me_fft==modulo(ig2-1,nproc_fft))
!   change due to //isation : the proc 0 must have me_g0=1
    if (ipw_ok) then
!    old if   if (me_fft==modulo(ig2+1,mpi_enreg%nproc_fft)) then
     if (ig2==0) mpi_enreg%me_g0=1
     v2=kg2arr(ig2p)
     gs_part=gmet(2,2)*v2*v2+gmet(3,3)*v3*v3+2.0_dp*gmet(2,3)*v2*v3
     gs_fact=2.0_dp*(gmet(1,2)*v2+gmet(3,1)*v3)
     ig1pmax=ngrid(1)
     if( (istwf_k==2.or.istwf_k==3) .and. ig3p==1 .and. ig2p==1)ig1pmax=nmax(1)+1
     do ig1p=1,ig1pmax
      v1=kg1arr(ig1p)
      gmin=gs_part+v1*(gs_fact+v1*gmet11)
!     If inside sphere:
      if (gmin<=gscut) then
       ig=ig+1  ! inside sphere
       if (ig<=mpw*np_band) then
!       Keep coords of pw:
        kg_small(1,ig)=ig1arr(ig1p)
        kg_small(2,ig)=ig2
        kg_small(3,ig)=ig3
       end if
       array_ipw(ig2)=array_ipw(ig2)+1
      end if
!     End loop on ig1p
     end do
    end if
!   End loop on ig2p
   end do
!  End loop on ig3p
  end do

! Add for future use of ind_fft_planes
! not completely tested
  if (mpi_enreg%paral_compil_fft==1) then
   n2 =size(mpi_enreg%ind_fft_planes,dim=2)
!  creation array array_ipw_1_n2
   if(istwf_k ==1 )then
    allocate(array_ipw_1_n2(n2))
   else
    allocate(array_ipw_1_n2(n2))
   end if
   array_ipw_1_n2(:)=0
   do i2=0,ig2pmax
    array_ipw_1_n2(i2+1)=array_ipw(i2)
   end do
   do i2=-1,-ig2pmax,-1
    array_ipw_1_n2(i2+n2+1)=array_ipw(i2)
   end do
   if(istwf_k>=2 .and. istwf_k<=5)then
    do i2=1,ig2pmax
     array_ipw_1_n2(n2-i2+1)=array_ipw(i2)
    end do
    if(istwf_k==2) then
     array_ipw_1_n2(1)=array_ipw_1_n2(1)*2-1
    else
     array_ipw_1_n2(1)=array_ipw_1_n2(1)*2
    end if

   else
    do i2=1,ig2pmax
     array_ipw_1_n2(n2+1-i2)=array_ipw(i2)
    end do
   end if
   mpi_enreg%nplanes_fft(ikpt)=0
   mpi_enreg%ind_fft_planes(ikpt,:)=0
   do i2=1,n2
    if (mpi_enreg%me_fft==modulo(i2+1,mpi_enreg%nproc_fft)) then
     mpi_enreg%nplanes_fft(ikpt)=mpi_enreg%nplanes_fft(ikpt)+1
     mpi_enreg%ind_fft_planes(ikpt,mpi_enreg%nplanes_fft(ikpt))=i2
    end if
   end do
   deallocate(array_ipw_1_n2)
  end if

 else

! ig2pmax was initialized previously
  mpi_enreg%me_g0=0
  do ig2p=1,ngrid(2)
   ig2=ig2arr(ig2p)
!  MPIWF Here, one select the set of planes in case of FFT parallelism
   ipw_ok=(istwf_k==1 .and. mpi_enreg%me_fft==modulo(ig2,mpi_enreg%nproc_fft) .or. &
&   (istwf_k>=2 .and. mpi_enreg%me_fft==modulo(ig2,mpi_enreg%nproc_fft)))
!  PAY ATTENTION : old if was ipw_ok=(((ig2p-1)/((ig2pmax+nproc_fft-1)/nproc_fft))==me_fft)
!  change due to //isation : the proc 0 must have me_g0=1
   if (ipw_ok) then
    if(ig2==0 .and. istwf_k>=2) mpi_enreg%me_g0=1
    v2=kg2arr(ig2p)
    ig3pmin=1
    if( (istwf_k==2 .or. istwf_k==3 .or. istwf_k==6 .or. istwf_k==7) .and. ig2<0)then
     ig3pmin=2
    end if
    do ig3p=ig3pmin,ig3pmax
     ig3=ig3arr(ig3p)
     v3=kg3arr(ig3p)
     gs_part=gmet(2,2)*v2*v2+gmet(3,3)*v3*v3+2.0_dp*gmet(2,3)*v2*v3
     gs_fact=2.0_dp*(gmet(1,2)*v2+gmet(3,1)*v3)
     ig1pmax=ngrid(1)
     if( (istwf_k==2.or.istwf_k==3) .and. ig3p==1 .and. ig2p==1)ig1pmax=nmax(1)+1
     do ig1p=1,ig1pmax
      v1=kg1arr(ig1p)
      gmin=gs_part+v1*(gs_fact+v1*gmet11)
!     If inside sphere:
      if (gmin<=gscut) then
       ig=ig+1  ! inside sphere
!      Make sure not to overrun array, or simply do not store if mpw=0
       if (ig<=mpw*np_band) then
!       Keep coords of pw:
        kg_small(1,ig)=ig1arr(ig1p)
        kg_small(2,ig)=ig2
        kg_small(3,ig)=ig3
       end if
      end if
!     End loop on ig1p
     end do
!    End loop on ig3p
    end do
   end if ! if the ig2 plane is to be treated by this processor
!  End loop on ig2p
  end do

 end if ! exchn2n3d==0 or ==1

!Total number of G vectors at this k point is assigned: npw
 npw=ig

!If npw exceeds array dimension mpw, while non-zero mkmem, call a halt
!if (npw>mpw .and. mkmem/=0) then
!write(message, '(a,a,a,a,i10,a,i10,a)' )ch10,&
!&   ' kpgsph : BUG -',ch10,&
!&   '  npw=',npw,' > mpw=',mpw,'.'
!call wrtout(06,message,'PERS')
!call leave_new('PERS')
!end if
 deallocate(ig1arr,ig2arr,ig3arr,kg1arr,kg2arr,kg3arr)

!bandFFT
 if(mpi_enreg%mode_para=='b') then
  npw_split=npw
  npw=npw/mpi_enreg%sizecart(2)
  if(mpi_enreg%coords(2) /= mpi_enreg%sizecart(2)-1) then
   if(mpw > 0) then ! This is for the case when we only compute npw and put mpw=0
    kg_small(:,1:npw)=kg_small(:,mpi_enreg%coords(2)*npw+1:(mpi_enreg%coords(2)+1)*npw)
   end if
  else
   npw_remain=modulo(npw_split,mpi_enreg%sizecart(2))
   if(mpw > 0) then ! This is for the case when we only compute npw and put mpw=0
    kg_small(:,1:npw+npw_remain)=kg_small(:,mpi_enreg%coords(2)*npw+1:(mpi_enreg%coords(2)+1)*npw+npw_remain)
   end if
   npw=npw+npw_remain
  end if
 end if

 if(mpw > 0) then
  do i1=1,npw
   kg(:,i1+ikg)=kg_small(:,i1)
  end do
! DEBUG
! write(6,*) 'In the loop npw lt mpw'
! ENDDEBUG
  deallocate(kg_small)
 end if

 deallocate(array_ipw)
!Take care of the me_g0 flag
 if(mpi_enreg%mode_para=='b') then
  if(mpi_enreg%coords(2)==0.and.mpi_enreg%me_g0==1) then
!  In this case, the processors had the 0 G vector before the new distribution, and still keeps it
   mpi_enreg%me_g0=1
  else
!  All other cases
   mpi_enreg%me_g0=0
  end if
 end if
 call timab(23,2,tsec)

!DEBUG
!write(6,*)' kpgsph : exit with npw=',npw
!do ig=1,min(npw,mpw)
!write(6, '(4i6)' )ig,kg(1:3,ig+ikg)
!end do
!if(npw<=mpw)stop
!ENDDEBUG

end subroutine kpgsph
!!***
