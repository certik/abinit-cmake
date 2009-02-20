!{\src2tex{textfont=tt}}
!!****f* ABINIT/dielmt2
!! NAME
!! dielmt2
!!
!!
!! FUNCTION
!! Compute dielectric matrix from susceptibility matrix
!! !!!DO NOT !Diagonalize it, or invert it.
!! ouput the diagonal part of the matrix

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
!!  npwdiel=size of the dielinv and susmat arrays.
!!  nspden=number of spin-density components
!!  occopt=option for occupancies
!!  prtvol=control print volume and debugging output
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!
!! OUTPUT
!!  dielinv(2,npwdiel,nspden,npwdiel,nspden)=inverse of the (non-hermitian)
!!      TC dielectric matrix in reciprocal space.
!!
!! NOTES
!! Warning : will not work in the spin-polarized, metallic case.
!! Output (not cleaned)
!! !!! Spin behaviour is not obvious !!!
!!
!! TODO
!! Write equation below (hermitian matrix)
!!
!! PARENTS
!!      prcref
!!
!! CHILDREN
!!      chpev,leave_new,timab,wrtout,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dielmt2(dielar,dielop,gmet,kg_diel,&
&  npwdiel,nspden,occopt,prtvol,susmat,&
& dieldiag,dtset,ucvol)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_14iowfdenpot
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 !real(dp),intent(out) :: dielinv(2,npwdiel,nspden,npwdiel,nspden)
!scalars
 integer,intent(in) :: dielop,npwdiel,nspden,occopt,prtvol
 real(dp),intent(in) :: ucvol
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: kg_diel(3,npwdiel)
 real(dp),intent(in) :: dielar(7),gmet(3,3)
 real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 real(dp),intent(out) :: dieldiag(2,npwdiel,nspden)

!Local variables-------------------------------
!scalars
 integer :: ieig,ier,ii,index,ipw,ipw1,ipw2,ipw3,isp,jj,npwsp
 real(dp) :: eiginv,elementi,elementr,gfact,gfactinv,gred1,gred2,gred3,gsquar
 real(dp) :: tpisq
 character(len=500) :: message
 character(len=fnlen) :: filapp
!arrays
 real(dp) :: buffer(dtset%nfft,dtset%nspden)
 real(dp) :: dielinv(2,npwdiel,nspden,npwdiel,nspden),tsec(2)
 real(dp) :: xred_cp(3,size(dtset%xred_orig,2))
 real(dp),allocatable :: dielh(:),dielmat(:,:,:,:,:),dielvec(:,:,:)
 real(dp),allocatable :: dielwk(:,:,:),eig_diel(:),zhpev1(:,:),zhpev2(:)

! *************************************************************************
 xred_cp=dtset%xred_orig
#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zhpev
#endif

!DEBUG
!write(6,*)' dielmt2 : enter '
!if(.true.)stop
!ENDDEBUG

!tpisq is (2 Pi) **2:
 tpisq=(two_pi)**2

 if(nspden==4)then
  write(6,*)' dielmt2 : does not work yet for nspden=4'
  stop
 end if

 call timab(90,1,tsec)

 if(nspden==2 .and. (occopt>=3 .and. occopt<=7) )then
  write(message, '(a,a,a,a,a,a)' ) ch10,&
&  ' dielmt2 : BUG -',ch10,&
&  '  In the present version of the code, one cannot produce',ch10,&
&  '  the dielectric matrix in the metallic, spin-polarized case.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!-Compute now the hermitian dielectric matrix------------------------------

!for the spin-unpolarized case, 1 - 4pi (1/G) chi0(G,Gp) (1/Gp)
!for the spin-polarized case,
!( 1  0 ) - 4pi  1/2  ( 1/G  1/G )   ( chi0 up      0    )   ( 1/Gp 1/Gp )
!( 0  1 )             ( 1/G  1/G )   (   0     chi0 down )   ( 1/Gp 1/Gp )
!which is equal to
!( 1  0 ) - 4pi  1/2 (1/G  0 ) (chi0 up+down   chi0 up+down) (1/Gp 0  )
!( 0  1 )            ( 0  1/G) (chi0 up+down   chi0 up+down) ( 0  1/Gp)
!So, if spin-polarized, sum spin up and spin down contributions

 allocate(dielmat(2,npwdiel,nspden,npwdiel,nspden))

 if(nspden/=1)then
  do ipw2=1,npwdiel
   do ipw1=1,npwdiel
    dielmat(1,ipw1,1,ipw2,1)=susmat(1,ipw1,1,ipw2,1)+susmat(1,ipw1,2,ipw2,2)
    dielmat(2,ipw1,1,ipw2,1)=susmat(2,ipw1,1,ipw2,1)+susmat(2,ipw1,2,ipw2,2)
   end do
  end do
 else
  do ipw2=1,npwdiel
   do ipw1=1,npwdiel
    dielmat(1,ipw1,1,ipw2,1)=susmat(1,ipw1,1,ipw2,1)
    dielmat(2,ipw1,1,ipw2,1)=susmat(2,ipw1,1,ipw2,1)
   end do
  end do
 end if

!Compute 1/G factors and include them in the dielectric matrix
 do ipw1=1,npwdiel
  gred1=dble(kg_diel(1,ipw1))
  gred2=dble(kg_diel(2,ipw1))
  gred3=dble(kg_diel(3,ipw1))
  gsquar=tpisq*(gmet(1,1)*gred1**2+gmet(2,2)*gred2**2+gmet(3,3)*gred3**2 &
&  +2.0d0*( (gmet(1,2)*gred2+gmet(1,3)*gred3)* gred1 +      &
&  gmet(2,3)*gred2*gred3)                        )
! Distinguish G=0 from other elements
  if(gsquar>1.0d-12)then
!  $ gfact=\sqrt (4.0d0 \pi/gsquar/dble(nspden))$
   gfact=sqrt(4.0d0*pi/gsquar)
   do ipw2=1,npwdiel
!   Must multiply both rows and columns, and also changes the sign
    dielmat(1,ipw2,1,ipw1,1)=-dielmat(1,ipw2,1,ipw1,1)*gfact
    dielmat(2,ipw2,1,ipw1,1)=-dielmat(2,ipw2,1,ipw1,1)*gfact
    dielmat(1,ipw1,1,ipw2,1)= dielmat(1,ipw1,1,ipw2,1)*gfact
    dielmat(2,ipw1,1,ipw2,1)= dielmat(2,ipw1,1,ipw2,1)*gfact
   end do
  else
!  Zero the G=0 elements, head and wings
   do ipw2=1,npwdiel
    dielmat(1,ipw2,1,ipw1,1)=0.0d0
    dielmat(2,ipw2,1,ipw1,1)=0.0d0
    dielmat(1,ipw1,1,ipw2,1)=0.0d0
    dielmat(2,ipw1,1,ipw2,1)=0.0d0
   end do
  end if
 end do

!Complete the matrix in the spin-polarized case
 if(nspden/=1)then
  do ipw1=1,npwdiel
   do ipw2=1,npwdiel
    dielmat(1,ipw1,1,ipw2,2)=dielmat(1,ipw1,1,ipw2,1)
    dielmat(2,ipw1,1,ipw2,2)=dielmat(2,ipw1,1,ipw2,1)
    dielmat(1,ipw1,2,ipw2,1)=dielmat(1,ipw1,1,ipw2,1)
    dielmat(2,ipw1,2,ipw2,1)=dielmat(2,ipw1,1,ipw2,1)
    dielmat(1,ipw1,2,ipw2,2)=dielmat(1,ipw1,1,ipw2,1)
    dielmat(2,ipw1,2,ipw2,2)=dielmat(2,ipw1,1,ipw2,1)
   end do
  end do
 end if

!DEBUG
!write(6,*)' dielmt2 : make dielmat equal to identity matrix '
!do ipw1=1,npwdiel
!do ipw2=1,npwdiel
!dielmat(1,ipw1,1,ipw2,1)=0.0d0
!dielmat(2,ipw1,1,ipw2,1)=0.0d0
!end do
!end do
!ENDDEBUG
 dieldiag(:,:,:)=zero
!Add the diagonal part
 do isp=1,nspden
  do ipw=1,npwdiel
   dielmat(1,ipw,isp,ipw,isp)=1.0d0+dielmat(1,ipw,isp,ipw,isp)
   dieldiag(1,ipw,isp)=dielmat(1,ipw,isp,ipw,isp)
  end do
 end do


 write(filapp,*)  'rdielmt'
 call out1dm(filapp,dtset%natom,dtset%nfft,dtset%ngfft,dtset%nspden,dtset%ntypat,&
& buffer,&
& dtset%rprimd_orig,dtset%typat,ucvol,&
& buffer,&
& xred_cp,dtset%znucl)


!!$ !-The hermitian dielectric matrix is computed ------------------------------
!!$ !-Now, diagonalize it ------------------------------------------------------
!!$
!!$ !In RPA, everything is projected on the spin-symmetrized
!!$ !space. This was coded here (for the time being).
!!$
!!$ !Diagonalize the hermitian dielectric matrix
!!$
!!$ !npwsp=npwdiel*nspden
!!$ npwsp=npwdiel
!!$
!!$ allocate(dielh(npwsp*(npwsp+1)))
!!$ allocate(dielvec(2,npwsp,npwsp))
!!$ allocate(eig_diel(npwsp))
!!$ allocate(zhpev1(2,2*npwsp-1),zhpev2(3*npwsp-2))
!!$ ier=0
!!$ !Store the dielectric matrix in proper mode before calling zhpev
!!$ index=1
!!$ do ii=1,npwdiel
!!$  do jj=1,ii
!!$   dielh(index  )=dielmat(1,jj,1,ii,1)
!!$   dielh(index+1)=dielmat(2,jj,1,ii,1)
!!$   index=index+2
!!$  end do
!!$ end do
!!$ !If spin-polarized and non RPA, need to store other parts of the matrix
!!$ !if(nspden/=1)then
!!$ ! do ii=1,npwdiel
!!$ !  Here, spin-flip contribution
!!$ !  do jj=1,npwdiel
!!$ !   dielh(index  )=dielmat(1,jj,1,ii,2)
!!$ !   dielh(index+1)=dielmat(2,jj,1,ii,2)
!!$ !   index=index+2
!!$ !  end do
!!$ !  Here spin down-spin down upper matrix
!!$ !  do jj=1,ii
!!$ !   dielh(index  )=dielmat(1,jj,2,ii,2)
!!$ !   dielh(index+1)=dielmat(2,jj,2,ii,2)
!!$ !   index=index+2
!!$ !  end do
!!$ ! end do
!!$ !end if
!!$
!!$#if defined T3E
!!$ call CHPEV ('V','U',npwsp,dielh,eig_diel,dielvec,npwdiel,zhpev1,&
!!$&   zhpev2,ier)
!!$#else
!!$ call ZHPEV ('V','U',npwsp,dielh,eig_diel,dielvec,npwdiel,zhpev1,&
!!$&   zhpev2,ier)
!!$#endif
!!$ deallocate(zhpev1,zhpev2)
!!$
!!$ if(prtvol>=10)then
!!$  write(message, '(a,a,a,5es12.4)' )ch10,&
!!$&  ' Five largest eigenvalues of the hermitian RPA dielectric matrix:',&
!!$&  ch10,eig_diel(npwdiel:npwdiel-4:-1)
!!$  call wrtout(ab_out,message,'COLL')
!!$ end if
!!$
!!$ write(message, '(a,a)' )ch10,&
!!$& ' dielmt2 : 15 largest eigenvalues of the hermitian RPA dielectric matrix'
!!$ call wrtout(6,message,'COLL')
!!$ write(message, '(a,5es12.5)' )'  1-5  :',eig_diel(npwdiel:npwdiel-4:-1)
!!$ call wrtout(6,message,'COLL')
!!$ write(message, '(a,5es12.5)' )'  6-10 :',eig_diel(npwdiel-5:npwdiel-9:-1)
!!$ call wrtout(6,message,'COLL')
!!$ write(message, '(a,5es12.5)' )'  11-15:',eig_diel(npwdiel-10:npwdiel-14:-1)
!!$ call wrtout(6,message,'COLL')
!!$ write(message, '(a,a)' )ch10,&
!!$& ' dielmt2 : 5 smallest eigenvalues of the hermitian RPA dielectric matrix'
!!$ call wrtout(6,message,'COLL')
!!$ write(message, '(a,5es12.5)' )'  1-5  :',eig_diel(1:5)
!!$ call wrtout(6,message,'COLL')
!!$
!!$ !Invert the hermitian dielectric matrix,
!!$ dielinv(:,:,:,:,:)=0.0d0
!!$ do ieig=1,npwdiel
!!$  eiginv=1.0d0/eig_diel(ieig)
!!$  do ipw2=1,npwdiel
!!$   do ipw1=1,npwdiel
!!$    dielinv(1,ipw1,1,ipw2,1)=dielinv(1,ipw1,1,ipw2,1)+&
!!$&          (dielvec(1,ipw1,ieig)*dielvec(1,ipw2,ieig)+ &
!!$&           dielvec(2,ipw1,ieig)*dielvec(2,ipw2,ieig) ) * eiginv
!!$    dielinv(2,ipw1,1,ipw2,1)=dielinv(2,ipw1,1,ipw2,1)+&
!!$&          (dielvec(2,ipw1,ieig)*dielvec(1,ipw2,ieig)- &
!!$&           dielvec(1,ipw1,ieig)*dielvec(2,ipw2,ieig) ) * eiginv
!!$   end do
!!$  end do
!!$ end do
!!$
!!$ deallocate(dielh,dielvec,eig_diel)
!!$
!!$ !DEBUG
!!$ !Checks whether the inverse of the hermitian dielectric matrix
!!$ !has been correctly generated
!!$ ! do ipw1=1,npwdiel
!!$ !  do ipw2=1,npwdiel
!!$ !   elementr=0.0d0
!!$ !   elementi=0.0d0
!!$ !   do ipw3=1,npwdiel
!!$ !    elementr=elementr+dielinv(1,ipw1,1,ipw3,1)*dielmat(1,ipw3,1,ipw2,1)&
!!$ !&                    -dielinv(2,ipw1,1,ipw3,1)*dielmat(2,ipw3,1,ipw2,1)
!!$ !    elementi=elementi+dielinv(1,ipw1,1,ipw3,1)*dielmat(2,ipw3,1,ipw2,1)&
!!$ !&                    +dielinv(2,ipw1,1,ipw3,1)*dielmat(1,ipw3,1,ipw2,1)
!!$ !   end do
!!$ !   if(elementr**2+elementi**2 > 1.0d-12)then
!!$ !    if( ipw1 /= ipw2 .or. &
!!$ !&        ( abs(elementr-1.0d0)>1.0d-6 .or. abs(elementi)>1.0d-6 ))then
!!$ !     write(6,*)' dielmt2 : the inversion procedure is not correct '
!!$ !     write(6,*)' ipw1, ipw2 =',ipw1,ipw2
!!$ !     write(6,*)' elementr,elementi=',elementr,elementi
!!$ !     stop
!!$ !    end if
!!$ !   end if
!!$ !  end do
!!$ ! end do
!!$ ! write(6,*)'dielmt2 : matrix has been inverted successfully '
!!$ !ENDDEBUG
!!$
!!$ !Then get the inverse of the asymmetric
!!$ !dielectric matrix, as required for the preconditioning.
!!$
!!$ !Inverse of the dielectric matrix : ( 1 - 4pi (1/G^2) chi0(G,Gp) )^(-1)
!!$ !In dielinv there is now (1 - 4pi (1/G) chi0(G,Gp) (1/Gp) )^(-1)
!!$ !So, evaluate dielinv_after(G,Gp) =
!!$ !                  (4pi/G^2)^(1/2) dielinv_before(G,Gp) (4pi/Gp^2)^(-1/2)
!!$ !In RPA, can focus on the spin-averaged quantities
!!$ do ipw1=1,npwdiel
!!$  gred1=dble(kg_diel(1,ipw1))
!!$  gred2=dble(kg_diel(2,ipw1))
!!$  gred3=dble(kg_diel(3,ipw1))
!!$  gsquar=tpisq*(gmet(1,1)*gred1**2+gmet(2,2)*gred2**2+gmet(3,3)*gred3**2 &
!!$&               +2.0d0*( (gmet(1,2)*gred2+gmet(1,3)*gred3)* gred1 +      &
!!$&                         gmet(2,3)*gred2*gred3)                        )
!!$ ! Distinguish G=0 from other elements
!!$  if(gsquar>1.0d-12)then
!!$   gfact=sqrt(4.0d0*pi/gsquar)
!!$   gfactinv=1.0d0/gfact
!!$   do ipw2=1,npwdiel
!!$ !   Must multiply both rows and columns
!!$    dielinv(1,ipw2,1,ipw1,1)=dielinv(1,ipw2,1,ipw1,1)*gfactinv
!!$    dielinv(2,ipw2,1,ipw1,1)=dielinv(2,ipw2,1,ipw1,1)*gfactinv
!!$    dielinv(1,ipw1,1,ipw2,1)=dielinv(1,ipw1,1,ipw2,1)*gfact
!!$    dielinv(2,ipw1,1,ipw2,1)=dielinv(2,ipw1,1,ipw2,1)*gfact
!!$   end do
!!$  else
!!$ !  Zero the G=0 elements, head and wings
!!$   do ipw2=1,npwdiel
!!$    dielinv(1,ipw2,1,ipw1,1)=0.0d0
!!$    dielinv(2,ipw2,1,ipw1,1)=0.0d0
!!$    dielinv(1,ipw1,1,ipw2,1)=0.0d0
!!$    dielinv(2,ipw1,1,ipw2,1)=0.0d0
!!$   end do
!!$  end if
!!$ end do

 deallocate(dielmat)

 call timab(90,2,tsec)

end subroutine dielmt2

!!***
