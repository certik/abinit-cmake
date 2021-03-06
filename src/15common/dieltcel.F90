!{\src2tex{textfont=tt}}
!!****f* ABINIT/dieltcel
!! NAME
!! dieltcel
!!
!! FUNCTION
!! Compute either test charge or electronic dielectric matrices
!! from susceptibility matrix
!! Diagonalize it, then invert it.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, LSI)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dielar(7)=input parameters for dielectric matrix:
!!              diecut,dielng,diemac,diemix,diegap,dielam.
!!  dielop=option for this routine (in development)
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kxc(nfft,nkxc)=exchange-correlation kernel,
!!       needed if the electronic dielectric matrix is computed
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  npwdiel=size of the dielinv and susmat arrays.
!!  nspden=number of spin-density components
!!  occopt=option for occupancies
!!  option=1 for Test Charge dielectric matrix, 2 for electronic dielectric matrix
!!  prtvol=control print volume and debugging output
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!
!! OUTPUT
!!  dielinv(2,npwdiel,nspden,npwdiel,nspden)=inverse of the (non-hermitian)
!!      TC dielectric matrix in reciprocal space.
!!
!! NOTES
!! Output (not cleaned)
!! !!! Spin behaviour is not obvious !!!
!! Will not work in the spin-polarized, metallic case.
!!
!! PARENTS
!!      prcref
!!
!! CHILDREN
!!      chpev,fourdp,leave_new,timab,wrtout,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dieltcel(dielar,dielinv,dielop,gmet,kg_diel,kxc,&
&  nfft,ngfft,nkxc,npwdiel,nspden,occopt,option,paral_kgb,prtvol,susmat)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dielop,nfft,nkxc,npwdiel,nspden,occopt,option,paral_kgb
 integer,intent(in) :: prtvol
!arrays
 integer,intent(in) :: kg_diel(3,npwdiel),ngfft(18)
 real(dp),intent(in) :: dielar(7),gmet(3,3),kxc(nfft,nkxc)
 real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 real(dp),intent(out) :: dielinv(2,npwdiel,nspden,npwdiel,nspden)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,ieig,ier,ifft,ii,index,ipw,ipw0,ipw1,ipw2,ipw3,ispden,j1
 integer :: j2,j3,jj,k1,k2,k3,n1,n2,n3
 real(dp) :: ai,ai2,ar,ar2,eiginv,elementi,elementr,gred1,gred2,gred3,gsquar,si
 real(dp) :: sr,tpisq
 character(len=500) :: message
 type(MPI_type) :: mpi_enreg
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: dielmat(:,:,:,:,:),eig_msusinvsqr(:),eig_msussqr(:)
 real(dp),allocatable :: eig_sus(:),eig_sym(:),invsqrsus(:,:,:,:,:)
 real(dp),allocatable :: khxc(:,:,:,:,:),kxcg(:,:),sqrsus(:,:,:,:,:),sush(:)
 real(dp),allocatable :: susvec(:,:,:),symdielmat(:,:,:,:,:),symh(:)
 real(dp),allocatable :: symvec(:,:,:,:,:),wkxc(:),work(:,:,:,:,:)
 real(dp),allocatable :: work2(:,:,:,:,:),zhpev1(:,:),zhpev2(:)
!no_abirules
!DEBUG
!Used to moderate divergence effect near rho=0 of the Kxc
!this limit value is truly empirical (exprmt on small Sr cell).
 real(dp) :: kxc_min=-200.0
!ENDDEBUG

! *************************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zhpev
#endif

!DEBUG
!write(6,*)' dieltcel : enter '
!if(.true.)stop
!ENDDEBUG

 call timab(96,1,tsec)

!tpisq is (2 Pi) **2:
 tpisq=(two_pi)**2

 if(nspden/=1 .and. (occopt>=3 .and. occopt<=7) )then
  write(message, '(a,a,a,a,a,a)' ) ch10,&
&  ' dieltcel : BUG -',ch10,&
&  '  In the present version of the code, one cannot produce',ch10,&
&  '  the dielectric matrix in the metallic, spin-polarized case.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if



!-Diagonalize the susceptibility matrix

 allocate(sush(npwdiel*(npwdiel+1)),susvec(2,npwdiel,npwdiel))
 allocate(eig_msusinvsqr(npwdiel),eig_msussqr(npwdiel),eig_sus(npwdiel))
 allocate(zhpev1(2,2*npwdiel-1),zhpev2(3*npwdiel-2))
 allocate(work(2,npwdiel,nspden,npwdiel,nspden))
 allocate(work2(2,npwdiel,nspden,npwdiel,nspden))
 allocate(sqrsus(2,npwdiel,nspden,npwdiel,nspden))
 allocate(invsqrsus(2,npwdiel,nspden,npwdiel,nspden))

!At some time, should take care of different spin channels
 do ispden=1,nspden

  if(nspden/=1)then
   write(6,*)' dieltcel : stop, nspden/=1'
   stop
  end if

! Store the susceptibility matrix in proper mode before calling zhpev
  index=1
  do ii=1,npwdiel
   do jj=1,ii
    sush(index  )=susmat(1,jj,1,ii,1)
    sush(index+1)=susmat(2,jj,1,ii,1)
    index=index+2
   end do
  end do

  ier=0
#if defined T3E
  call CHPEV ('V','U',npwdiel,sush,eig_sus,susvec,npwdiel,zhpev1,&
&  zhpev2,ier)
#else
  call ZHPEV ('V','U',npwdiel,sush,eig_sus,susvec,npwdiel,zhpev1,&
&  zhpev2,ier)
#endif

! DEBUG
! write(6,*)' dieltcel : print eigenvalues of the susceptibility matrix'
! do ii=1,npwdiel
! write(6, '(i5,es16.6)' )ii,eig_sus(ii)
! end do
! ENDDEBUG

  do ii=1,npwdiel
   if(-eig_sus(ii)>1.d-12)then
    eig_msussqr(ii)=sqrt(-eig_sus(ii))
    eig_msusinvsqr(ii)=1._dp/eig_msussqr(ii)
   else if(-eig_sus(ii)< -1.d-12)then
    write(message, '(a,a,a,a)' ) ch10,&
&    ' dieltcel : BUG -',ch10,&
&    '  Found positive eigenvalue of susceptibility matrix.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   else
!   Set the eigenvalue corresponding to a constant potential change to 1,
!   while it will be set to zero in Khx.
    eig_msussqr(ii)=1._dp
    eig_msusinvsqr(ii)=1._dp
   end if
  end do

! Compute square root of minus susceptibility matrix
! and inverse square root of minus susceptibility matrix
  do ii=1,npwdiel
   work(:,:,1,ii,1)=susvec(:,:,ii)*eig_msussqr(ii)
   work2(:,:,1,ii,1)=susvec(:,:,ii)*eig_msusinvsqr(ii)
  end do
  do ipw2=1,npwdiel
   do ipw1=ipw2,npwdiel
    ar=0._dp ; ai=0._dp ; ar2=0._dp ; ai2=0._dp
    do ii=1,npwdiel
     sr=susvec(1,ipw2,ii) ; si=susvec(2,ipw2,ii)
     ar =ar  +work(1,ipw1,1,ii,1)*sr  +work(2,ipw1,1,ii,1)*si
     ai =ai  +work(2,ipw1,1,ii,1)*sr  -work(1,ipw1,1,ii,1)*si
     ar2=ar2 +work2(1,ipw1,1,ii,1)*sr +work2(2,ipw1,1,ii,1)*si
     ai2=ai2 +work2(2,ipw1,1,ii,1)*sr -work2(1,ipw1,1,ii,1)*si
    end do
    sqrsus(1,ipw1,1,ipw2,1)=ar
    sqrsus(2,ipw1,1,ipw2,1)=ai
    invsqrsus(1,ipw1,1,ipw2,1)=ar2
    invsqrsus(2,ipw1,1,ipw2,1)=ai2
    if(ipw1/=ipw2)then
     sqrsus(1,ipw2,1,ipw1,1)=ar
     sqrsus(2,ipw2,1,ipw1,1)=-ai
     invsqrsus(1,ipw2,1,ipw1,1)=ar2
     invsqrsus(2,ipw2,1,ipw1,1)=-ai2
    end if
   end do
  end do

! DEBUG
! Checks whether sqrsus and invsqrsus are inverse of each other.
! do ipw1=1,npwdiel
! do ipw2=1,npwdiel
! elementr=0.0_dp
! elementi=0.0_dp
! do ipw3=1,npwdiel
! elementr=elementr+sqrsus(1,ipw1,1,ipw3,1)*invsqrsus(1,ipw3,1,ipw2,1)&
! &                    -sqrsus(2,ipw1,1,ipw3,1)*invsqrsus(2,ipw3,1,ipw2,1)
! elementi=elementi+sqrsus(1,ipw1,1,ipw3,1)*invsqrsus(2,ipw3,1,ipw2,1)&
! &                    +sqrsus(2,ipw1,1,ipw3,1)*invsqrsus(1,ipw3,1,ipw2,1)
! end do
! if(elementr**2+elementi**2 > 1.0d-12)then
! if( ipw1 /= ipw2 .or. &
! &        ( abs(elementr-1.0_dp)>1.0d-6 .or. abs(elementi)>1.0d-6 ))then
! write(6,*)' dieltcel : sqrsus and invsqrsus are not (pseudo)',&
! &        'inverse of each other'
! write(6,*)' ipw1, ipw2 =',ipw1,ipw2
! write(6,*)' elementr,elementi=',elementr,elementi
! stop
! end if
! end if
! end do
! end do
! ENDDEBUG

! End loop over spins
 end do

 deallocate(eig_msusinvsqr,eig_msussqr,eig_sus,sush,susvec)

!-Compute the Hxc kernel

 allocate(khxc(2,npwdiel,nspden,npwdiel,nspden))
 allocate(symdielmat(2,npwdiel,nspden,npwdiel,nspden))

 khxc(:,:,:,:,:)=0.0_dp

!Compute Hartree kernel
 do ipw1=1,npwdiel
  gred1=dble(kg_diel(1,ipw1))
  gred2=dble(kg_diel(2,ipw1))
  gred3=dble(kg_diel(3,ipw1))
  gsquar=tpisq*(gmet(1,1)*gred1**2+gmet(2,2)*gred2**2+gmet(3,3)*gred3**2 &
&  +2.0_dp*( (gmet(1,2)*gred2+gmet(1,3)*gred3)* gred1 +      &
&  gmet(2,3)*gred2*gred3)                        )
! Distinguish G=0 from other elements
  if(gsquar>1.0d-12)then
   khxc(1,ipw1,1,ipw1,1)= 4.0_dp*pi/gsquar
  else
!  G=0
   ipw0=ipw1
  end if
 end do

!Eventually add the xc part
 if(option>=2)then

  allocate(wkxc(nfft),kxcg(2,nfft))
  wkxc(:)=kxc(:,1)
! DEBUG
! write(0,*) kxc(:,1)
! Used to moderate divergenc effect near rho=0 of the Kxc (see above).
! wkxc(:)=merge(kxc(:,1), kxc_min, kxc(:,1) > kxc_min)
! ENDDEBUG
  mpi_enreg%me_fft=0
  mpi_enreg%nproc_fft=1
  call fourdp(1,kxcg,wkxc,-1,mpi_enreg,nfft,ngfft,paral_kgb,0) ! trsfrm R to G

! Compute difference in G vectors
  n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
  do ipw2=1,npwdiel
   if(ipw2/=ipw0)then

    j1=kg_diel(1,ipw2) ; j2=kg_diel(2,ipw2) ; j3=kg_diel(3,ipw2)
!   Fills diagonal
    khxc(1,ipw2,1,ipw2,1)=khxc(1,ipw2,1,ipw2,1)+kxcg(1,1)
    khxc(2,ipw2,1,ipw2,1)=khxc(2,ipw2,1,ipw2,1)+kxcg(2,1)

    if(ipw2/=npwdiel)then
!    Fills off-diagonal part of the matrix, except G=0
     do ipw1=ipw2+1,npwdiel
      if(ipw1/=ipw0)then
       i1=kg_diel(1,ipw1) ; i2=kg_diel(2,ipw1) ; i3=kg_diel(3,ipw1)
!      Use of two mod calls handles both i1-j1>=ndiel1 AND i1-j1<0
       k1=mod(n1+mod(i1-j1,n1),n1)
       k2=mod(n2+mod(i2-j2,n2),n2)
       k3=mod(n3+mod(i3-j3,n3),n3)
       ifft=k1+1+n1*(k2+n2*k3)
!      The signs of imaginary contributions have been checked
       khxc(1,ipw1,1,ipw2,1)=kxcg(1,ifft)
       khxc(2,ipw1,1,ipw2,1)=kxcg(2,ifft)
       khxc(1,ipw2,1,ipw1,1)=kxcg(1,ifft)
       khxc(2,ipw2,1,ipw1,1)=-kxcg(2,ifft)
      end if
     end do
    end if

   end if
  end do

  deallocate(wkxc,kxcg)

! Endif option 2
 end if

!Now, get the symmetric dielectric matrix
!Premultiplication by square root of minus susceptibility matrix
 do ipw2=1,npwdiel
  do ipw1=1,npwdiel
   ar=0._dp ; ai=0._dp
   do ii=1,npwdiel
    ar=ar+sqrsus(1,ipw1,1,ii,1)*khxc(1,ii,1,ipw2,1) &
&    -sqrsus(2,ipw1,1,ii,1)*khxc(2,ii,1,ipw2,1)
    ai=ai+sqrsus(2,ipw1,1,ii,1)*khxc(1,ii,1,ipw2,1) &
&    +sqrsus(1,ipw1,1,ii,1)*khxc(2,ii,1,ipw2,1)
   end do
   work(1,ipw1,1,ipw2,1)=ar
   work(2,ipw1,1,ipw2,1)=ai
  end do
 end do
!Postmultiplication by square root of minus susceptibility matrix
 do ipw2=1,npwdiel
! do ipw1=ipw2,npwdiel
  do ipw1=1,npwdiel
   ar=0._dp ; ai=0._dp
   do ii=1,npwdiel
    ar=ar+work(1,ipw1,1,ii,1)*sqrsus(1,ii,1,ipw2,1) &
&    -work(2,ipw1,1,ii,1)*sqrsus(2,ii,1,ipw2,1)
    ai=ai+work(2,ipw1,1,ii,1)*sqrsus(1,ii,1,ipw2,1) &
&    +work(1,ipw1,1,ii,1)*sqrsus(2,ii,1,ipw2,1)
   end do
   symdielmat(1,ipw1,1,ipw2,1)=ar
   symdielmat(2,ipw1,1,ipw2,1)=ai
!  if(ipw1/=ipw2)then
!  symdielmat(1,ipw2,1,ipw1,1)=ar
!  symdielmat(2,ipw2,1,ipw1,1)=-ai
!  end if
  end do
! Add the unity matrix
  symdielmat(1,ipw2,1,ipw2,1)=1._dp+symdielmat(1,ipw2,1,ipw2,1)
 end do

 deallocate(khxc)

 allocate(symh(npwdiel*(npwdiel+1)),symvec(2,npwdiel,nspden,npwdiel,nspden))
 allocate(eig_sym(npwdiel))

!Store the symmetrized dielectric matrix in proper mode before calling zhpev
 index=1
 do ii=1,npwdiel
  do jj=1,ii
   symh(index  )=symdielmat(1,jj,1,ii,1)
   symh(index+1)=symdielmat(2,jj,1,ii,1)
   index=index+2
  end do
 end do

 ier=0
#if defined T3E
 call CHPEV ('V','U',npwdiel,symh,eig_sym,symvec,npwdiel,zhpev1,&
& zhpev2,ier)
#else
 call ZHPEV ('V','U',npwdiel,symh,eig_sym,symvec,npwdiel,zhpev1,&
& zhpev2,ier)
#endif

 if(prtvol>=10)then
  write(message, '(a,a,a,5es12.4)' )ch10,&
&  ' Five largest eigenvalues of the symmetrized dielectric matrix:',&
&  ch10,eig_sym(npwdiel:npwdiel-4:-1)
  call wrtout(ab_out,message,'COLL')
 end if

 write(message, '(a,a)' )ch10,&
& ' dieltcel : 15 largest eigenvalues of the symmetrized dielectric matrix'
 call wrtout(6,message,'COLL')
 write(message, '(a,5es12.5)' )'  1-5  :',eig_sym(npwdiel:npwdiel-4:-1)
 call wrtout(6,message,'COLL')
 write(message, '(a,5es12.5)' )'  6-10 :',eig_sym(npwdiel-5:npwdiel-9:-1)
 call wrtout(6,message,'COLL')
 write(message, '(a,5es12.5)' )'  11-15:',eig_sym(npwdiel-10:npwdiel-14:-1)
 call wrtout(6,message,'COLL')
 write(message, '(a,a)' )ch10,&
& ' dieltcel : 5 smallest eigenvalues of the symmetrized dielectric matrix'
 call wrtout(6,message,'COLL')
 write(message, '(a,5es12.5)' )'  1-5  :',eig_sym(1:5)
 call wrtout(6,message,'COLL')

!Invert the hermitian dielectric matrix,
 work(:,:,:,:,:)=0.0_dp
 do ieig=1,npwdiel
  eiginv=1.0_dp/eig_sym(ieig)
  do ipw2=1,npwdiel
!  do ipw1=ipw2,npwdiel
   do ipw1=1,npwdiel
    work(1,ipw1,1,ipw2,1)=work(1,ipw1,1,ipw2,1)+&
&    (symvec(1,ipw1,1,ieig,1)*symvec(1,ipw2,1,ieig,1)+ &
&    symvec(2,ipw1,1,ieig,1)*symvec(2,ipw2,1,ieig,1) ) * eiginv
    work(2,ipw1,1,ipw2,1)=work(2,ipw1,1,ipw2,1)+&
&    (symvec(2,ipw1,1,ieig,1)*symvec(1,ipw2,1,ieig,1)- &
&    symvec(1,ipw1,1,ieig,1)*symvec(2,ipw2,1,ieig,1) ) * eiginv
   end do
  end do
 end do
!if(npwdiel>1)then
!do ipw2=2,npwdiel
!do ipw1=1,ipw2-1
!work(1,ipw1,1,ipw2,1)= work(1,ipw2,1,ipw1,1)
!work(2,ipw1,1,ipw2,1)=-work(2,ipw2,1,ipw1,1)
!end do
!end do
!end if

 deallocate(eig_sym,symh,symvec)

!DEBUG
!Checks whether the inverse of the symmetric dielectric matrix
!has been correctly generated
!do ipw1=1,npwdiel
!do ipw2=1,npwdiel
!elementr=0.0_dp
!elementi=0.0_dp
!do ipw3=1,npwdiel
!elementr=elementr+work(1,ipw1,1,ipw3,1)*symdielmat(1,ipw3,1,ipw2,1)&
!&                    -work(2,ipw1,1,ipw3,1)*symdielmat(2,ipw3,1,ipw2,1)
!elementi=elementi+work(1,ipw1,1,ipw3,1)*symdielmat(2,ipw3,1,ipw2,1)&
!&                    +work(2,ipw1,1,ipw3,1)*symdielmat(1,ipw3,1,ipw2,1)
!end do
!if(elementr**2+elementi**2 > 1.0d-12)then
!if( ipw1 /= ipw2 .or. &
!&        ( abs(elementr-1.0_dp)>1.0d-6 .or. abs(elementi)>1.0d-6 ))then
!write(6,*)' dieltcel : the inversion procedure is not correct '
!write(6,*)' ipw1, ipw2 =',ipw1,ipw2
!write(6,*)' elementr,elementi=',elementr,elementi
!stop
!end if
!end if
!end do
!end do
!write(6,*)'dieltcel : matrix has been inverted successfully '
!ENDDEBUG

 deallocate(symdielmat)

!Then get the inverse of the asymmetric
!dielectric matrix, as required for the preconditioning.
!Premultiplication by square root of minus susceptibility matrix
 do ipw2=1,npwdiel
  do ipw1=1,npwdiel
   ar=0._dp ; ai=0._dp
   do ii=1,npwdiel
    ar=ar+invsqrsus(1,ipw1,1,ii,1)*work(1,ii,1,ipw2,1) &
&    -invsqrsus(2,ipw1,1,ii,1)*work(2,ii,1,ipw2,1)
    ai=ai+invsqrsus(2,ipw1,1,ii,1)*work(1,ii,1,ipw2,1) &
&    +invsqrsus(1,ipw1,1,ii,1)*work(2,ii,1,ipw2,1)
   end do
   work2(1,ipw1,1,ipw2,1)=ar
   work2(2,ipw1,1,ipw2,1)=ai
  end do
 end do
!Postmultiplication by square root of minus susceptibility matrix
 do ipw2=1,npwdiel
  do ipw1=1,npwdiel
   ar=0._dp ; ai=0._dp
   do ii=1,npwdiel
    ar=ar+work2(1,ipw1,1,ii,1)*sqrsus(1,ii,1,ipw2,1) &
&    -work2(2,ipw1,1,ii,1)*sqrsus(2,ii,1,ipw2,1)
    ai=ai+work2(2,ipw1,1,ii,1)*sqrsus(1,ii,1,ipw2,1) &
&    +work2(1,ipw1,1,ii,1)*sqrsus(2,ii,1,ipw2,1)
   end do
   dielinv(1,ipw1,1,ipw2,1)=ar
   dielinv(2,ipw1,1,ipw2,1)=ai
  end do
 end do

 deallocate(invsqrsus,sqrsus,work,work2)
 deallocate(zhpev1,zhpev2)

 call timab(96,2,tsec)

!DEBUG
!write(6,*)' dieltcel : exit, will stop'
!stop
!ENDDEBUG

end subroutine dieltcel
!!***
