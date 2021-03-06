!{\src2tex{textfont=tt}}
!!****f* ABINIT/elast9
!!
!! NAME
!! elast9
!!
!! FUNCTION
!! Get the elastic and compliance tensors, both clamped ion and relaxed ion,
!! under the fixed electric field boundary condition; in which realxed ion
!! tensors can generate two output tensors one is conventional, the other
!! considers the sress correction.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XW)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! anaddb_dtset= (derived datatype) contains all the input variables
!! blkval(2,3,mpert,3,mpert,nblok)=
!!   second derivatives of total energy with respect to electric fields
!!   atom displacements,strain,...... all in cartesian coordinates
!! iblok= bolk number in DDB file
!! iblok_stress= blok number which contain stress tensor
!! instrain=force response internal strain tensor
!! iout=out file number
!! mpert=maximum number of ipert
!! natom=number of atoms in unit cell
!! nblok=number of total bloks in DDB file
!! ucvol=unit cell volume
!!
!! OUTPUT
!! elast=relaxed-ion elastic tensor(without stress correction) (6*6) in Voigt notation
!!
!! NOTES
!! The elastic (compliance) tensors calculated here are under boundary conditions of
!! fixed Electric Field, different from those in piezo9.F90 which are under fixed
!! Displacement Field and incorporate piezoelectric corrections.
!!
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      matrginv,wrtout,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine elast9(anaddb_dtset,blkval,elast,iblok,iblok_stress,instrain,iout,mpert,&
&            natom,nblok,ucvol)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments -------------------------------------------
!scalars
 integer,intent(in) :: iblok,iblok_stress,iout,mpert,natom,nblok
 real(dp),intent(in) :: ucvol
 type(anaddb_dataset_type),intent(in) :: anaddb_dtset
!arrays
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok),instrain(3*natom,6)
 real(dp),intent(out) :: elast(6,6)

!Local variables------------------------------------
!scalars
 integer :: ier,ii1,ii2,ipert1,ipert2,ivarA,ivarB
 character(len=500) :: message
!arrays
 real(dp) :: Amatr(3*natom-3,3*natom-3),Apmatr(3*natom,3*natom)
 real(dp) :: Bmatr(2,((3*natom-3)*(3*natom-2))/2)
 real(dp) :: Bpmatr(2,(3*natom*(3*natom+1))/2),Cmatr(3*natom-3,3*natom-3)
 real(dp) :: Cpmatr(3*natom,3*natom),Nmatr(3*natom,3*natom),compl_clamped(6,6)
 real(dp) :: compl_relaxed(6,6),compl_stress(6,6),eigval(3*natom-3)
 real(dp) :: eigvalp(3*natom),eigvec(2,3*natom-3,3*natom-3)
 real(dp) :: eigvecp(2,3*natom,3*natom),elast_clamped(6,6),elast_relaxed(6,6)
 real(dp) :: elast_stress(6,6),kmatrix(3*natom,3*natom),new1(6,3*natom)
 real(dp) :: new2(6,6),stress(6),zhpev1(2,2*3*natom-4),zhpev1p(2,2*3*natom-1)
 real(dp) :: zhpev2(3*3*natom-5),zhpev2p(3*3*natom-2)

!***************************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zhpev
#endif

!extration of the elastic constants from the blkvals

 do ivarA=1,6
  do ivarB=1,6
!  because the elastic constant is 6*6,
!  so we should judge if the idir is larger than 3
!  or not
   if(ivarA>3) then
    ii1=ivarA-3
    ipert1=natom+4  !for the shear modulus
   else if(ivarA<=3) then
    ii1=ivarA
    ipert1=natom+3  !for the diagonal part
   end if
   if(ivarB>3) then
    ii2=ivarB-3
    ipert2=natom+4  !for the shear modulus
   else if(ivarB<=3) then
    ii2=ivarB
    ipert2=natom+3  !for the diagonal part
   end if
   elast(ivarA,ivarB)=blkval(1,ii1,ipert1,ii2,ipert2,iblok)
  end do
 end do

!then consider the volume,cause the unit above is in
!Hartree, in fact the elastic constant should be in
!the units of presure, the energy/volume
!And then transform the unit to si unit using GPa
!from Hartree/Bohr^3

 do ivarA=1,6
  do ivarB=1,6
   elast(ivarA,ivarB)=(elast(ivarA,ivarB)/ucvol)*HaBohr3_GPa
  end do
 end do

!then should consider the two situations:clamped and relaxed
!ions respectively,give the initial value of elast_clamped
 compl_clamped(:,:)=elast(:,:)

!then do the matrix mulplication of instrain*K*instrain to get the
!correction of the relaxed ion quantities

 if(anaddb_dtset%elaflag==2 .or. anaddb_dtset%elaflag==3&
& .or. anaddb_dtset%elaflag==4 .or. anaddb_dtset%elaflag==5)then
! extracting force matrix at gamma
  do ipert1=1,natom
   do ii1=1,3
    ivarA=ii1+3*(ipert1-1)
    do ipert2=1,natom
     do ii2=1,3
      ivarB=ii2+3*(ipert2-1)
      kmatrix(ivarA,ivarB)=blkval(1,ii1,ipert1,ii2,ipert2,iblok)
     end do
    end do
   end do
  end do

! DEBUG
! write(6,'(/,a,/)')'the k matrix before inverse'
! do ii1=1,3*natom
! write(6,'(6es16.3)')kmatrix(ii1,1),kmatrix(ii1,2),kmatrix(ii1,3),&
! & kmatrix(ii1,4),kmatrix(ii1,5),kmatrix(ii1,6)
! end do
! ENDDEBUG

! according to formula, invert the kmatrix(3natom,3natom)
  Apmatr(:,:)=kmatrix(:,:)

  Nmatr(:,:)=0.0_dp
  do ivarA=1,3*natom
   do ivarB=1,3*natom
    if (mod(ivarA,3)==0 .and. mod(ivarB,3)==0)then
     Nmatr(ivarA,ivarB)=1
    end if
    if (mod(ivarA,3)==1 .and. mod(ivarB,3)==1)then
     Nmatr(ivarA,ivarB)=1
    end if
    if (mod(ivarA,3)==2 .and. mod(ivarB,3)==2)then
     Nmatr(ivarA,ivarB)=1
    end if
   end do
  end do

! DEBUG
! write(6,'(/,a,/)')'the direct in verse of the Kmatrix'
! do ivarA=1,3*natom
! write(6,'(/)')
! do ivarB=1,3*natom
! write(6,'(es16.6)')kmatrix(ivarA,ivarB)
! end do
! end do
! ENDDEBUG


! strarting the pseudoinervering processes
! then get the eigenvectors of the big matrix,give values to matrixBp
  ii1=1
  do ivarA=1,3*natom
   do ivarB=1,ivarA
    Bpmatr(1,ii1)=Nmatr(ivarB,ivarA)
    ii1=ii1+1
   end do
  end do
  Bpmatr(2,:)=0.0_dp  !the imaginary part of the force matrix
! then call the subroutines CHPEV and ZHPEV to get the eigenvectors
  call ZHPEV ('V','U',3*natom,Bpmatr,eigvalp,eigvecp,3*natom,&
&  zhpev1p,zhpev2p,ier)

! DEBUG
! the rigenval and eigenvec
! write(6,'(/,a,/)')'the eigenvalues and eigenvectors'
! do ivarA=1,3*natom
! write(6,'(/)')
! write(6,'(es16.6)')eigvalp(ivarA)
! end do
! do ivarA=1,3*natom
! write(6,'(/)')
! do ivarB=1,3*natom
! write(6,'(es16.6)')eigvecp(1,ivarB,ivarA)
! end do
! end do
! ENDDEBUG

! the do the muplication to get the reduced matrix,in two steps
  Cpmatr(:,:)=0.0_dp
  do ivarA=1,3*natom
   do ivarB=1,3*natom
    do ii1=1,3*natom
     Cpmatr(ivarA,ivarB)=Cpmatr(ivarA,ivarB)+eigvecp(1,ii1,ivarA)*&
&     Apmatr(ii1,ivarB)
    end do
   end do
  end do

  Apmatr(:,:)=0.0_dp
  do ivarA=1,3*natom
   do ivarB=1,3*natom
    do ii1=1,3*natom
     Apmatr(ivarA,ivarB)=Apmatr(ivarA,ivarB)+Cpmatr(ivarA,ii1)*&
&     eigvecp(1,ii1,ivarB)
    end do
   end do
  end do

! DEBUG
! the blok diagonal parts
! write(6,'(/,a,/)')'Apmatr'
! do ivarA=1,3*natom
! write(6,'(/)')
! do ivarB=1,3*natom
! write(6,'(es16.6)')Apmatr(ivarA,ivarB)
! end do
! end do
! ENDDEBUG


! check the last three eigenvalues whether too large
  ivarB=0
  do ivarA=3*natom-2,3*natom
   if (ABS(Apmatr(ivarA,ivarA))>tol6)then
    ivarB=1
   end if
  end do
  if(ivarB==1)then
   write(message,'(a,a,a,a,a,a,a,a,a,a,3es16.6)')ch10,&
&   ' elast9 : WARNING - :',ch10,&
&   '  Acoustic sum rule violation met : the eigenvalues of accoustic mode',ch10,&
&   '  are too large at Gamma point',ch10,&
&   '  increase cutoff energy or k-points sampling.',ch10,&
&   '  The three eigenvalues are:',Apmatr(3*natom-2,3*natom-2),&
&   Apmatr(3*natom-1,natom-1),Apmatr(3*natom,3*natom)
   call wrtout(06, message, 'COLL')
   call wrtout(iout,message,'COLL')
  end if
! then give the value of reduced matrix form Apmatr to Amatr
  do ivarA=1,3*natom-3
   do ivarB=1,3*natom-3
    Amatr(ivarA,ivarB)=Apmatr(ivarA,ivarB)
   end do
  end do

! now the reduced matrix is in the Amatr, the convert it
! first give the give the value of Bmatr from Amatr
  ii1=1
  do ivarA=1,3*natom-3
   do ivarB=1,ivarA
    Bmatr(1,ii1)=Amatr(ivarB,ivarA)
    ii1=ii1+1
   end do
  end do
  Bmatr(2,:)=0.0_dp
! then call the subroutines CHPEV and ZHPEV to get the eigenvectors and the
! eigenvalues
  call ZHPEV ('V','U',3*natom-3,Bmatr,eigval,eigvec,3*natom-3,&
&  zhpev1,zhpev2,ier)

! check the unstable phonon modes, if the first is negative then print
! warning message
  if(eigval(1)<-1.0*tol8)then
   write(message,'(a,a,a,a,a,a)') ch10,&
&   ' elast9 : WARNING - :',ch10,&
&   '  Unstable eigenvalue detected in force constant matrix at Gamma point.',ch10,&
&   '  The system under calculation is physically unstable.'
   call wrtout(06, message, 'COLL')
   call wrtout(iout,message,'COLL')
  end if

! the do the matrix muplication to get pseudoinverse inverse matrix
  Cmatr(:,:)=0.0_dp
  Amatr(:,:)=0.0_dp
  do ivarA=1,3*natom-3
   Cmatr(ivarA,ivarA)=1.0_dp/eigval(ivarA)
  end do
  do ivarA=1,3*natom-3
   do ivarB=1,3*natom-3
    do ii1=1,3*natom-3
     Amatr(ivarA,ivarB)=Amatr(ivarA,ivarB)+eigvec(1,ivarA,ii1)*&
&     Cmatr(ii1,ivarB)
    end do
   end do
  end do

! then the second mulplication
  Cmatr(:,:)=0.0_dp
  do ivarA=1,3*natom-3
   do ivarB=1,3*natom-3
    do ii1=1,3*natom-3
     Cmatr(ivarA,ivarB)=Cmatr(ivarA,ivarB)+&
&     Amatr(ivarA,ii1)*eigvec(1,ivarB,ii1)
    end do
   end do
  end do

! DEBUG
! write(6,'(/,a,/)')'the pseudo inverse of the force matrix'
! do ivarA=1,3*natom
! write(6,'(/)')
! do ivarB=1,3*natom
! write(6,'(es16.6)')Cmatr(ivarA,ivarB)
! end do
! end do
! ENDDEBUG

! so now the inverse of the reduced matrix is in the matrixC
! now do another mulplication to get the pseudoinverse of the original
  Cpmatr(:,:)=0.0_dp
  Apmatr(:,:)=0.0_dp
  do ivarA=1,3*natom-3
   do ivarB=1,3*natom-3
    Cpmatr(ivarA,ivarB)=Cmatr(ivarA,ivarB)
   end do
  end do

! times the eigvecp
  do ivarA=1,3*natom
   do ivarB=1,3*natom
    do ii1=1,3*natom
     Apmatr(ivarA,ivarB)=Apmatr(ivarA,ivarB)+eigvecp(1,ivarA,ii1)*&
&     Cpmatr(ii1,ivarB)
    end do
   end do
  end do
  Cpmatr(:,:)=0.0_dp
  do ivarA=1,3*natom
   do ivarB=1,3*natom
    do ii1=1,3*natom
     Cpmatr(ivarA,ivarB)=Cpmatr(ivarA,ivarB)+&
&     Apmatr(ivarA,ii1)*eigvecp(1,ivarB,ii1)
    end do
   end do
  end do

! now the inverse in in Cpmatr
  kmatrix(:,:)=Cpmatr(:,:)
! transfer the inverse of k-matrix back to the k matrix
! so now the inverse of k matrix is in the kmatrix
! ending the part for pseudoinversing the K matrix
! then do the first matrix mulplication
  new1(:,:)=0.0_dp
  do ii1=1,6
   do ii2=1,3*natom
    do ivarA=1,3*natom
     new1(ii1,ii2)=new1(ii1,ii2)+instrain(ivarA,ii1)*&
&     kmatrix(ivarA,ii2)
    end do
   end do
  end do
! then do the second matrix mulplication, and change the value of kmatrix
  new2(:,:)=0.0_dp
  do ii1=1,6
   do ii2=1,6
    do ivarB=1,3*natom
     new2(ii1,ii2)=new2(ii1,ii2)+new1(ii1,ivarB)*&
&     instrain(ivarB,ii2)
    end do
   end do
  end do
! then finish the matrix mupl., consider the unit cellvolume
! and the unit change next step
  do ivarA=1,6
   do ivarB=1,6
    new2(ivarA,ivarB)=(new2(ivarA,ivarB)/ucvol)*HaBohr3_GPa
   end do
  end do
! then the relaxed one should be the previous one minus the new2 element
  do ivarA=1,6
   do ivarB=1,6
    elast_clamped(ivarA,ivarB)=elast(ivarA,ivarB)
    elast_relaxed(ivarA,ivarB)=elast_clamped(ivarA,ivarB)-&
&    new2(ivarA,ivarB)
   end do
  end do
 end if
!the above end if end if for elaflag=2 or elafalg=3 or elafalg=4,
!or elafalg=5 in line 125

!DEBUG
!write(6,'(/,a,/)')'debug the unit cell volume'
!write(6,'(2es16.6)')ucvol,HaBohr3_GPa
!ENDDEBUG

!then give the initial value of the compl_relaxed(6,6)
 compl_relaxed(:,:)=elast_relaxed(:,:)

!*******************************************************************
 if(anaddb_dtset%elaflag==1.or. anaddb_dtset%elaflag==3)then
! print out the clamped-ion elastic constants to output file
  write(message,'(3a)')ch10,&
&  ' Elastic Tensor (clamped ion) (unit:10^2GP):',ch10
  call wrtout(06,message,'COLL')
  do ivarA=1,6
   write(6,'(6f12.7)')elast(ivarA,1)/100.00_dp,elast(ivarA,2)/100.00_dp,&
&   elast(ivarA,3)/100.00_dp,elast(ivarA,4)/100.00_dp,&
&   elast(ivarA,5)/100.00_dp,elast(ivarA,6)/100.00_dp
  end do

  call wrtout(iout,message,'COLL')
  do ivarA=1,6
   write(iout,'(6f12.7)')elast(ivarA,1)/100.00_dp,elast(ivarA,2)/100.00_dp,&
&   elast(ivarA,3)/100.00_dp,elast(ivarA,4)/100.00_dp,&
&   elast(ivarA,5)/100.00_dp,elast(ivarA,6)/100.00_dp
  end do
 end if

 if(anaddb_dtset%elaflag==2.or.anaddb_dtset%elaflag==3&
& .or. anaddb_dtset%elaflag==4.or. anaddb_dtset%elaflag==5)then
  if(anaddb_dtset%instrflag==0)then
   write(message,'(a,a,a,a,a,a,a,a)' )ch10,&
&   ' WARNING: in order to get the elastic  tensor(relaxed ion), ',ch10,&
&   '  one needs information about internal strain ',ch10,&
&   '  one should set  instrflag==1;',ch10,&
&   '  otherwise the program will continue but give wrong values.'
   call wrtout(6,message,'coll')
   call wrtout(iout,message,'coll')
  end if

  write(message,'(5a)')ch10,&
&  ' Elastic Tensor (relaxed ion) (unit:10^2GP):',ch10,&
&  '  (at fixed electric field boundary condition)',ch10
  call wrtout(06,message,'COLL')
  do ivarA=1,6
   write(6,'(6f12.7)')elast_relaxed(ivarA,1)/100.00_dp,&
&   elast_relaxed(ivarA,2)/100.00_dp,elast_relaxed(ivarA,3)/100.00_dp,&
&   elast_relaxed(ivarA,4)/100.00_dp,elast_relaxed(ivarA,5)/100.00_dp,&
&   elast_relaxed(ivarA,6)/100.00_dp
  end do

  call wrtout(iout,message,'COLL')
  do ivarA=1,6
   write(iout,'(6f12.7)')elast_relaxed(ivarA,1)/100.00_dp,&
&   elast_relaxed(ivarA,2)/100.00_dp,elast_relaxed(ivarA,3)/100.00_dp,&
&   elast_relaxed(ivarA,4)/100.00_dp,elast_relaxed(ivarA,5)/100.00_dp,&
&   elast_relaxed(ivarA,6)/100.00_dp
  end do
 end if

!then print the corresponding compliances

 if(anaddb_dtset%elaflag==1.or.anaddb_dtset%elaflag==3)then
! compl(:,:)=elast_clamped(:,:) !convert the elastic tensor
  call matrginv(compl_clamped,6,6)
  write(message,'(a,a,a)')ch10,&
&  ' Compliance Tensor (clamped ion) (unit: 10^-2GP^-1):',ch10
  call wrtout(06,message,'COLL')

  do ivarB=1,6
   write(6,'(6f12.7)')compl_clamped(ivarB,1)*100.00_dp,&
&   compl_clamped(ivarB,2)*100.00_dp,&
&   compl_clamped(ivarB,3)*100.00_dp,compl_clamped(ivarB,4)*100.00_dp,&
&   compl_clamped(ivarB,5)*100.00_dp,&
&   compl_clamped(ivarB,6)*100.00_dp
  end do

  call wrtout(iout,message,'COLL')

  do ivarB=1,6
   write(iout,'(6f12.7)')compl_clamped(ivarB,1)*100.00_dp,&
&   compl_clamped(ivarB,2)*100.00_dp,&
&   compl_clamped(ivarB,3)*100.00_dp,compl_clamped(ivarB,4)*100.00_dp,&
&   compl_clamped(ivarB,5)*100.00_dp,&
&   compl_clamped(ivarB,6)*100.00_dp
  end do
 end if

 if(anaddb_dtset%elaflag==2.or.anaddb_dtset%elaflag==3&
& .or. anaddb_dtset%elaflag==4 .or. anaddb_dtset%elaflag==5)then
! compl(:,:)=elast_relaxed(:,:)
  call matrginv(compl_relaxed,6,6)
  if(anaddb_dtset%instrflag==0)then
   write(message,'(a,a,a,a,a,a,a,a)' )ch10,&
&   ' WARNING: in order to get the compliance tensor(relaxed ion), ',ch10,&
&   '  one needs information about internal strain ',ch10,&
&   '  one should set  instrflag==1;',ch10,&
&   '  otherwise the program will continue but give wrong values.'
   call wrtout(6,message,'coll')
   call wrtout(iout,message,'coll')
  end if
  write(message,'(5a)')ch10,&
&  ' Compliance Tensor (relaxed ion)  (unit: 10^-2GP^-1):',ch10,&
&  '  (at fixed electric field boundary condition)',ch10
  call wrtout(06,message,'COLL')

  do ivarB=1,6
   write(6,'(6f12.7)')compl_relaxed(ivarB,1)*100.00_dp,&
&   compl_relaxed(ivarB,2)*100.00_dp,&
&   compl_relaxed(ivarB,3)*100.00_dp,compl_relaxed(ivarB,4)*100.00_dp,&
&   compl_relaxed(ivarB,5)*100.00_dp,&
&   compl_relaxed(ivarB,6)*100.00_dp
  end do
  call wrtout(iout,message,'COLL')

  do ivarB=1,6
   write(iout,'(6f12.7)')compl_relaxed(ivarB,1)*100.00,&
&   compl_relaxed(ivarB,2)*100.00_dp,&
&   compl_relaxed(ivarB,3)*100.00_dp,compl_relaxed(ivarB,4)*100.00_dp,&
&   compl_relaxed(ivarB,5)*100.00_dp,&
&   compl_relaxed(ivarB,6)*100.00_dp
  end do

 end if

!DEBUG
!write(6,'(/,6es16.6,/)')ucvol
!ENDDEBUG

!befor the end , make sure the tensor elast(6,6)
!will have the relaxed ion values
 elast(:,:)=elast_relaxed(:,:)

!begin the part of computing stress corrected elastic tensors
 if(anaddb_dtset%elaflag==5)then

! DEBUG
! check the iblok number of first derivative of energy
! write(6,'(/,a,/)')'iblok number at 8:00Pm'
! write(6, '(i)')iblok_stress
! write(6, '(a,f12.7)')'the total energy', blkval(1,1,1)
! write(6,*)'',blkval(1,:,:,:,:,iblok_stress)
! write(6,*)'',blkval(1,:,7,1,1,iblok_stress)
! ENDDEBUG

! firts give the corect stress values
! diagonal parts
  stress(1)=blkval(1,1,natom+3,1,1,iblok_stress)
  stress(2)=blkval(1,2,natom+3,1,1,iblok_stress)
  stress(3)=blkval(1,3,natom+3,1,1,iblok_stress)
! the shear parts
  stress(4)=blkval(1,1,natom+4,1,1,iblok_stress)
  stress(5)=blkval(1,2,natom+4,1,1,iblok_stress)
  stress(6)=blkval(1,3,natom+4,1,1,iblok_stress)
! then convert the unit from atomic to the GPa unit
  do ivarA=1,6
   stress(ivarA)=stress(ivarA)*HaBohr3_GPa
  end do
! give the initial values of elast_stress tensor
  elast_stress(:,:)=elast_relaxed(:,:)
! notice that only the first three rows need to be corrected
  do ivarA=1,3
   do ivarB=1,6
    elast_stress(ivarA,ivarB)=elast_stress(ivarA,ivarB)-stress(ivarB)
   end do
  end do
! then compute the values of compliance tensor with stress correction
  compl_stress(:,:)=elast_stress(:,:)
  call matrginv(compl_stress,6,6)
! then print out the results of stress corrected elastic and compliance tensors
  if(anaddb_dtset%instrflag==0)then
   write(message,'(a,a,a,a,a,a,a,a)' )ch10,&
&   ' WARNING: in order to get the elastic tensor (relaxed ion with stress correction), ',ch10,&
&   '  one needs information about internal strain ',ch10,&
&   '  one should set  instrflag==1;',ch10,&
&   '  otherwise the program will continue but give wrong values.'
   call wrtout(6,message,'coll')
   call wrtout(iout,message,'coll')
  end if
  write(message,'(5a)')ch10,&
&  ' Elastic Tensor (relaxed ion with stress corrected) (unit:10^2GP)',ch10,&
&  '  (at fixed electric field boundary condition)',ch10
  call wrtout(06,message,'COLL')
  call wrtout(iout,message,'COLL')
  do ivarA=1,6
   write(6,'(6f12.7)')elast_stress(ivarA,1)/100.00_dp,elast_stress(ivarA,2)/100.00_dp,&
&   elast_stress(ivarA,3)/100.00_dp,elast_stress(ivarA,4)/100.00_dp,&
&   elast_stress(ivarA,5)/100.00_dp,elast_stress(ivarA,6)/100.00_dp
  end do
  do ivarA=1,6
   write(iout,'(6f12.7)')elast_stress(ivarA,1)/100.00_dp,elast_stress(ivarA,2)/100.00_dp,&
&   elast_stress(ivarA,3)/100.00_dp,elast_stress(ivarA,4)/100.00_dp,&
&   elast_stress(ivarA,5)/100.00_dp,elast_stress(ivarA,6)/100.00_dp
  end do

! then the complinace tensors with stress correction
  write(message,'(5a)')ch10,&
&  ' Compliance Tensor (relaxed ion with stress correction) (unit: 10^-2(GP)^-1):',ch10,&
&  '  (at fixed electric field boundary condition)',ch10
  call wrtout(06,message,'COLL')
  call wrtout(iout,message,'COLL')
  do ivarB=1,6
   write(6,'(6f12.7)')compl_stress(ivarB,1)*100.00_dp,&
&   compl_stress(ivarB,2)*100.00_dp,&
&   compl_stress(ivarB,3)*100.00_dp,compl_stress(ivarB,4)*100.00_dp,&
&   compl_stress(ivarB,5)*100.00_dp,&
&   compl_stress(ivarB,6)*100.00_dp
  end do
  do ivarB=1,6
   write(iout,'(6f12.7)')compl_stress(ivarB,1)*100.00_dp,&
&   compl_stress(ivarB,2)*100.00_dp,&
&   compl_stress(ivarB,3)*100.00_dp,compl_stress(ivarB,4)*100.00_dp,&
&   compl_stress(ivarB,5)*100.00_dp,&
&   compl_stress(ivarB,6)*100.00_dp
  end do
 end if
!end the if 510th line
!end the part of stress corrected elastic and compliance tensors

end subroutine elast9

!!***
