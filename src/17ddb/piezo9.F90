!{\src2tex{textfont=tt}}
!!****f* ABINIT/piezo9
!!
!! NAME
!! piezo9
!!
!! FUNCTION
!! Get the piezoelectric tensor (e-tensor), both clamped ion and relaxed ion;
!! Compute physical(relaxed ion) piezoeletric (d, g, h) tensors;
!! Compute relaxed ion and free stress dielectric tensor;
!! Compute relaxed ion elatic and compliance tensors under fixed
!! displacement field boundary conditions.
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
!! dielt_rlx=relaxed ion dielectric tensor
!! iblok= bolk number in DDB file contains 2 derivative of energy
!! instrain=force response internal strain tensor
!! iout=out file number
!! mpert=maximum number of ipert
!! natom=number of atoms in unit cell
!! nblok=number of total bloks in DDB file
!! ucvol=unit cell volume
!!
!! OUTPUT
!! piezo = piezoelectric tensor
!!
!! NOTES
!! The elastic (compliance) tensors calculated here are under fixed D-field boundary
!! condition, which include piezoelectric corrections to the elastic (compliance)
!! tensors calculated in elast9.F90 whose boundary condition is fixed E-field.
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

subroutine piezo9(anaddb_dtset,blkval,dielt_rlx,elast,iblok,instrain,iout,mpert,&
&            natom,nblok,piezo,ucvol)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments-------------------------------------------
!scalars
 integer,intent(in) :: iblok,iout,mpert,natom,nblok
 real(dp),intent(in) :: ucvol
 type(anaddb_dataset_type),intent(in) :: anaddb_dtset
!arrays
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok),dielt_rlx(3,3)
 real(dp),intent(in) :: elast(6,6),instrain(3*natom,6)
 real(dp),intent(out) :: piezo(6,3)

!Local variables---------------------------------------
!scalars
 integer :: idir1,idir2,ier,ii1,ii2,ipert1,ipert2,ivarA,ivarB
 character(len=500) :: message
!arrays
 real(dp) :: Amatr(3*natom-3,3*natom-3),Apmatr(3*natom,3*natom)
 real(dp) :: Bmatr(2,((3*natom-3)*(3*natom-2))/2)
 real(dp) :: Bpmatr(2,(3*natom*(3*natom+1))/2),Cmatr(3*natom-3,3*natom-3)
 real(dp) :: Cpmatr(3*natom,3*natom),Nmatr(3*natom,3*natom),beta_tensor(3,3)
 real(dp) :: compliance(6,6),compliance_dis(6,6)
 real(dp) :: d2cart_relaxed(2,3,mpert,3,mpert,nblok),d_tensor(6,3)
 real(dp) :: dielt_stress(3,3),eigval(3*natom-3),eigvalp(3*natom)
 real(dp) :: eigvec(2,3*natom-3,3*natom-3),eigvecp(2,3*natom,3*natom)
 real(dp) :: elast_dis(6,6),g_tensor(3,6),h_tensor(3,6)
 real(dp) :: kmatrix(3*natom,3*natom),new1(6,3*natom),piezo_clamped(6,3)
 real(dp) :: piezo_correction(6,3),piezo_relaxed(6,3),zhpev1(2,2*3*natom-4)
 real(dp) :: zhpev1p(2,2*3*natom-1),zhpev2(3*3*natom-5),zhpev2p(3*3*natom-2)
 real(dp) :: zstar1(3,3*natom),zstar2(3*natom,3)

!****************************************************************

!DEBUG
!write(6,*)' piezo9 : enter '
!ENDDEBUG

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zhpev
#endif

!extracion of the clamped ion piezoelectric constants from blkvals

!the six strain perturbations
 do ivarA=1,6
! the three E-field perturbations
  do ivarB=1,3
!  judge if the ivarA>3 or not
   if(ivarA>3) then
    idir1=ivarA-3
    ipert1=natom+4
!   for the shear part of the strain
   else if(ivarA<=3) then
    idir1=ivarA
    ipert1=natom+3
!   for the diagonal part of strain
   end if
   idir2=ivarB
   ipert2=natom+2 !for the E-field perturbation only
   piezo(ivarA,ivarB)=blkval(1,idir2,ipert2,idir1,ipert1,iblok)
  end do
 end do

!consider the volume nad the -Qe before the piezo
!according to the (30) in notes, the units are tranformed from
!atomic units to the SI units

 do ivarA=1,6
  do ivarB=1,3
   piezo(ivarA,ivarB)=piezo(ivarA,ivarB)*AmuBohr2_Cm2
!  now it is in the SI unit
  end do
 end do

!DEBUG
!write(06,'es16.6')ucvol
!ENDDEBUG

!give the values of d2cart_relaxed as the same as blkval
!and also give the initial values of piezo_clamped and piezo_relaxed

 d2cart_relaxed(:,:,:,:,:,:)=blkval(:,:,:,:,:,:)
 piezo_clamped(:,:)=piezo(:,:)

!********************************************************************
!print the main results of the piezoelectric constants
 if(anaddb_dtset%piezoflag==1.or.anaddb_dtset%piezoflag==3&
& .or. anaddb_dtset%piezoflag==7)then
  write(message,'(3a)')ch10,&
&  ' Proper piezoelectric constants (clamped ion) (unit:c/m^2)',ch10
  call wrtout(06,message,'COLL')

  do ivarA=1,6
   write(6,'(3f16.8)')piezo_clamped(ivarA,1),piezo_clamped(ivarA,2),&
&   piezo_clamped(ivarA,3)
  end do

  call wrtout(iout,message,'COLL')
  do ivarA=1,6
   write(iout,'(3f16.8)')piezo_clamped(ivarA,1),piezo_clamped(ivarA,2),&
&   piezo_clamped(ivarA,3)
  end do
 end if

!the next is the calculation of the relaxed ion piezoelectric constants
!first extract the K(force constant) matrix

!if (piezoflag==2 .or. anaddb_dtset%piezoflag==3)then
!extracting force matrix at gamma
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

 Apmatr(:,:)=kmatrix(:,:)

!DEBUG
!kmatrix values
!write(6,'(/,a,/)')'the force constant matrix'
!do ivarA=1,3*natom
!write(6,'(/)')
!do ivarB=1,3*natom
!write(6,'(es16.6)')kmatrix(ivarB,ivarA)
!end do
!end do
!ENDDEBUG

 Nmatr(:,:)=zero

 do ivarA=1,3*natom
  do ivarB=1,3*natom
   if (mod(ivarA,3)==0 .and. mod(ivarB,3)==0)then
    Nmatr(ivarA,ivarB)=one
   end if
   if (mod(ivarA,3)==1 .and. mod(ivarB,3)==1)then
    Nmatr(ivarA,ivarB)=one
   end if
   if (mod(ivarA,3)==2 .and. mod(ivarB,3)==2)then
    Nmatr(ivarA,ivarB)=one
   end if
  end do
 end do

!DEBUG
!do ivarA=1,3*natom
!write(6,'(/)')
!do ivarB=1,3*natom
!write(6,'(es16.6)')Nmatr(ivarB,ivarA)
!end do
!end do
!ENDDEBUG

!starting the pseudoinversing processes
!then get the eigenvectors of the big matrix,give values to matrixBp
 ii1=1
 do ivarA=1,3*natom
  do ivarB=1,ivarA
   Bpmatr(1,ii1)=Nmatr(ivarB,ivarA)
   ii1=ii1+1
  end do
 end do
!the imaginary part of the force matrix
 Bpmatr(2,:)=0.0_dp
!then call the subroutines CHPEV and ZHPEV to get the eigenvectors
 call ZHPEV ('V','U',3*natom,Bpmatr,eigvalp,eigvecp,3*natom,&
& zhpev1p,zhpev2p,ier)

!DEBUG
!the eigenval and eigenvec
!write(6,'(/,a,/)')'the eigenvalues and eigenvectors'
!do ivarA=1,3*natom
!write(6,'(/)')
!write(6,'(es16.6)')eigvalp(ivarA)
!end do
!do ivarA=1,3*natom
!write(6,'(/)')
!do ivarB=1,3*natom
!write(6,'(es16.6)')eigvecp(1,ivarB,ivarA)
!end do
!end do
!ENDDEBUG

!then do the muplication to get the reduced matrix,in two steps
!After this the force constant matrix is decouple in two bloks,
!accoustic and optical ones
 Cpmatr(:,:)=0.0_dp
 do ivarA=1,3*natom
  do ivarB=1,3*natom
   do ii1=1,3*natom
    Cpmatr(ivarA,ivarB)=Cpmatr(ivarA,ivarB)+eigvecp(1,ii1,ivarA)*&
&    Apmatr(ii1,ivarB)
   end do
  end do
 end do

 Apmatr(:,:)=0.0_dp
 do ivarA=1,3*natom
  do ivarB=1,3*natom
   do ii1=1,3*natom
    Apmatr(ivarA,ivarB)=Apmatr(ivarA,ivarB)+Cpmatr(ivarA,ii1)*&
&    eigvecp(1,ii1,ivarB)
   end do
  end do
 end do

!DEBUG
!the blok diago
!write(6,'(/,a,/)')'Apmatr'
!do ivarA=1,3*natom
!write(6,'(/)')
!do ivarB=1,3*natom
!write(6,'(es16.6)')Apmatr(ivarA,ivarB)
!end do
!end do
!ENDDEBUG

!check the last three eigenvalues whether too large or not
 ivarB=0
 do ivarA=3*natom-2,3*natom
  if (ABS(Apmatr(ivarA,ivarA))>tol6)then
   ivarB=1
  end if
 end do
 if(ivarB==1)then
  write(message,'(a,a,a,a,a,a,a,a,a,a,3es16.6)')ch10,&
&  ' piezo9 : WARNING -',ch10,&
&  '  Acoustic sum rule violation met : the eigenvalues of accoustic mode',ch10,&
&  '  are too large at Gamma point',ch10,&
&  '  Increase cutoff energy or k-points sampling.',ch10,&
&  '  The three eigenvalues are:',Apmatr(3*natom-2,3*natom-2),&
&  Apmatr(3*natom-1,natom-1),Apmatr(3*natom,3*natom)
  call wrtout(06, message, 'COLL')
  call wrtout(iout,message,'COLL')
 end if

!give the value of reduced matrix form Apmatr to Amatr
 do ivarA=1,3*natom-3
  do ivarB=1,3*natom-3
   Amatr(ivarA,ivarB)=Apmatr(ivarA,ivarB)
  end do
 end do
!now the reduced matrix is in the matrixA, the convert it
!first give the give the value of matixB from matrixA
 ii1=1
 do ivarA=1,3*natom-3
  do ivarB=1,ivarA
   Bmatr(1,ii1)=Amatr(ivarB,ivarA)
   ii1=ii1+1
  end do
 end do
 Bmatr(2,:)=0.0_dp

!call the subroutines CHPEV and ZHPEV to get the eigenvectors and the
!eigenvalues
 call ZHPEV ('V','U',3*natom-3,Bmatr,eigval,eigvec,3*natom-3,&
& zhpev1,zhpev2,ier)

!check the unstable phonon modes, if the first is negative then print
!warning message
 if(eigval(1)<-1.0*tol8)then
  write(message,'(a,a,a,a,a,a)') ch10,&
&  ' piezo9 : WARNING -',ch10,&
&  '  Unstable eigenvalue detected in force constant matrix at Gamma point',ch10,&
&  '  The system under calculation is physically unstable.'
  call wrtout(06, message, 'COLL')
  call wrtout(iout,message,'COLL')
 end if

!do the matrix multiplication to get pseudoinverse inverse matrix
 Cmatr(:,:)=0.0_dp
 Amatr(:,:)=0.0_dp
 do ivarA=1,3*natom-3
  Cmatr(ivarA,ivarA)=1.0_dp/eigval(ivarA)
 end do

 do ivarA=1,3*natom-3
  do ivarB=1,3*natom-3
   do ii1=1,3*natom-3
    Amatr(ivarA,ivarB)=Amatr(ivarA,ivarB)+eigvec(1,ivarA,ii1)*&
&    Cmatr(ii1,ivarB)
   end do
  end do
 end do

!the second mulplication
 Cmatr(:,:)=0.0_dp
 do ivarA=1,3*natom-3
  do ivarB=1,3*natom-3
   do ii1=1,3*natom-3
    Cmatr(ivarA,ivarB)=Cmatr(ivarA,ivarB)+&
&    Amatr(ivarA,ii1)*eigvec(1,ivarB,ii1)
   end do
  end do
 end do

!DEBUG
!write(6,'(/,a,/)')'the pseudo inverse of the force matrix'
!do ivarA=1,3*natom
!write(6,'(/)')
!do ivarB=1,3*natom
!write(6,'(es16.6)')Cmatr(ivarA,ivarB)
!end do
!end do
!ENDDEBUG

!so now the inverse of the reduced matrix is in the matrixC
!now do another mulplication to get the pseudoinverse of the original
 Cpmatr(:,:)=0.0_dp
 Apmatr(:,:)=0.0_dp
 do ivarA=1,3*natom-3
  do ivarB=1,3*natom-3
   Cpmatr(ivarA,ivarB)=Cmatr(ivarA,ivarB)
  end do
 end do

!now times the eigvecp
 do ivarA=1,3*natom
  do ivarB=1,3*natom
   do ii1=1,3*natom
    Apmatr(ivarA,ivarB)=Apmatr(ivarA,ivarB)+eigvecp(1,ivarA,ii1)*&
&    Cpmatr(ii1,ivarB)
   end do
  end do
 end do

 Cpmatr(:,:)=0.0_dp
 do ivarA=1,3*natom
  do ivarB=1,3*natom
   do ii1=1,3*natom
    Cpmatr(ivarA,ivarB)=Cpmatr(ivarA,ivarB)+&
&    Apmatr(ivarA,ii1)*eigvecp(1,ivarB,ii1)
   end do
  end do
 end do

!now the inverse in in Cpmatr
 kmatrix(:,:)=Cpmatr(:,:)
!transfer the inverse of k-matrix back to the k matrix
!so now the inverse of k matrix is in the kmatrix
!ending the part for pseudoinversing the K matrix

!we still need the z-star matrix
 do idir1=1,3
  d2cart_relaxed(1,idir1,natom+2,idir1,natom+2,iblok)=&
&  d2cart_relaxed(1,idir1,natom+2,idir1,natom+2,iblok)-1.0_dp
 end do

 do idir1=1,3
  do idir2=1,3
   do ii1=1,2
    d2cart_relaxed(ii1,idir1,natom+2,idir2,natom+2,iblok)=&
&    d2cart_relaxed(ii1,idir1,natom+2,idir2,natom+2,iblok)/four_pi
   end do
  end do
 end do

 do ivarA=1,3
  idir1=ivarA
  ipert1=natom+2
  do ipert2=1,natom
   do idir2=1,3
    ivarB=idir2+3*(ipert2-1)
    zstar1(ivarA,ivarB)=d2cart_relaxed(1,idir1,ipert1,&
&    idir2,ipert2,iblok)
   end do
  end do
 end do

!then get the inverse of the zstar1 for zstar2(3*natom,3)
 do ivarA=1,3*natom
  do ivarB=1,3
   zstar2(ivarA,ivarB)=zstar1(ivarB,ivarA)
  end do
 end do
!the the matrix I need for the multiplication is in kmatrix and
!zstar2
!the first matrix mulplication
 new1(:,:)=0.0_dp
 do ii1=1,6
  do ii2=1,3*natom
   do ivarA=1,3*natom
    new1(ii1,ii2)=new1(ii1,ii2)+instrain(ivarA,ii1)*&
&    kmatrix(ivarA,ii2)
   end do
  end do
 end do

!do the second matrix mulplication
 piezo_correction(:,:)=0.0_dp
 do ii1=1,6
  do ii2=1,3
   do ivarA=1,3*natom
    piezo_correction(ii1,ii2)=piezo_correction(ii1,ii2)+&
&    new1(ii1,ivarA)*&
&    zstar2(ivarA,ii2)
   end do
  end do
 end do

!then consider the volume and the change the unit form atomic to SI
 do ii1=1,6
  do ii2=1,3
   piezo_correction(ii1,ii2)= (piezo_correction(ii1,ii2)&
&   /ucvol)*AmuBohr2_Cm2
   piezo_relaxed(ii1,ii2)=piezo_clamped(ii1,ii2)+&
&   piezo_correction(ii1,ii2)
  end do
 end do
!end the calculation of piezoelectric constants

!then print out the relaxed ion piezoelectric constants

 if(anaddb_dtset%piezoflag==2.or.anaddb_dtset%piezoflag==3&
& .or. anaddb_dtset%piezoflag==7)then
  if(anaddb_dtset%instrflag==0)then
   write(message,'(a,a,a,a,a,a,a,a)' )ch10,&
&   ' WARNING: in order to get the piezoelectric tensor (relaxed ion), ',ch10,&
&   '  one needs information about internal strain ',ch10,&
&   '  one should set  instrflag==1;',ch10,&
&   '  otherwise the program will continue but will give wrong values.'
   call wrtout(6,message,'coll')
   call wrtout(iout,message,'coll')
  end if
  write(message,'(3a)')ch10,&
&  ' Proper piezoelectric constants (relaxed ion) (unit:c/m^2)',ch10
  call wrtout(06,message,'COLL')
  do ivarA=1,6
   write(6,'(3f16.8)')piezo_relaxed(ivarA,1),piezo_relaxed(ivarA,2),&
&   piezo_relaxed(ivarA,3)
  end do

  call wrtout(iout,message,'COLL')
  do ivarA=1,6
   write(iout,'(3f16.8)')piezo_relaxed(ivarA,1),piezo_relaxed(ivarA,2),&
&   piezo_relaxed(ivarA,3)
  end do
 end if

!DEBUG
!check the values of the relaxed ion dielectric tensor
!write(message,'(a,a,a,a)')ch10,'debug the dielt tensor values ',&
!&  '(unit:c/m^2)',ch10
!call wrtout(06,message,'COLL')
!do ivarA=1,3
!write(6, '(3f16.8)')dielt_rlx(ivarA,1),dielt_rlx(ivarA,2),dielt_rlx(ivarA,3)
!end do
!END DEBUG

!DEBUG
!print the relaxed ion elast tensor
!write(message,'(a,a,a,a)')ch10,' debugElastic Tensor(relaxed ion)',&
!&  '(unit:10^2GP,VOIGT notation):',ch10
!call wrtout(06,message,'COLL')
!do ivarA=1,6
!write(6,'(6f12.7)')elast(ivarA,1)/100.00_dp,elast(ivarA,2)/100.00_dp,&
!&   elast(ivarA,3)/100.00_dp,elast(ivarA,4)/100.00_dp,&
!&   elast(ivarA,5)/100.00_dp,elast(ivarA,6)/100.00_dp
!end do
!ENDDEBUG

!Start to compute the piezoelectric d tensors
!first initialize the d_tensor values
!first make sure the elastic tensor is not zero
 if(anaddb_dtset%elaflag>1)then
  d_tensor(:,:)=0
! then get the relaxed ion compliance tensor
  compliance(:,:)=elast(:,:)
! DEBUG
! write(6,*)' call matrginv 1 '
! ENDDEBUG
  call matrginv(compliance,6,6)
  do ivarA=1,6
   do ivarB=1,3
    do ii1=1,6
     d_tensor(ivarA,ivarB)=d_tensor(ivarA,ivarB)+compliance(ivarA,ii1)*&
&     piezo_relaxed(ii1,ivarB)
    end do
   end do
  end do
 end if
!then convert in to the right unit pc/N
 do ivarA=1,6
  do ivarB=1,3
   d_tensor(ivarA,ivarB)=1000*d_tensor(ivarA,ivarB)
  end do
 end do
!then print out the results of d tensor in log and output files
 if(anaddb_dtset%piezoflag==4 .or. anaddb_dtset%piezoflag==7)then
  if(anaddb_dtset%instrflag==0 .or. anaddb_dtset%elaflag==0&
&  .or. anaddb_dtset%elaflag==1)then
   write(message,'(12a)' )ch10,&
&   ' WARNING:in order to get the piezoelectric d tensor(relaxed ion),', ch10,&
&   ' one needs the elastic tensor(relaxed ion) and piezoelectric e tensor',ch10,&
&   ' the latter needs the information of internal strain;',ch10,&
&   ' please check that both instrflag and elaflag are set to correct numbers',ch10,&
&   ' (elaflag= 2,3,4, or 5; instrflag=1)',ch10,&
&   ' otherwise the program  will continue but give wrong values.'
   call wrtout(06,message,'COLL')
   call wrtout(iout,message,'COLL')
  end if
  write(message,'(3a)')ch10,&
&  ' Piezoelectric d tensor (relaxed ion) (unit:pc/N)',ch10
  call wrtout(06,message,'COLL')
  do ivarA=1,6
   write(06,'(3f16.8)')d_tensor(ivarA,1),d_tensor(ivarA,2),d_tensor(ivarA,3)
  end do
  call wrtout(iout,message,'COLL')
  do ivarA=1,6
   write(iout,'(3f16.8)')d_tensor(ivarA,1),d_tensor(ivarA,2),d_tensor(ivarA,3)
  end do
 end if
!end the part of piezoelectric d tensor (relaxed ion only).

!then start to compute the piezoelectric g tensor
!according to the equation, we first need to know the information
!of the free-stress dielectric tensor
!first make sure dielt_rlx exits, so we do not invert zero matrix
 if(anaddb_dtset%dieflag>0)then
  do ivarA=1,3
   do ivarB=1,3
    dielt_stress(ivarA,ivarB)=zero
    do ii1=1,6
     dielt_stress(ivarA,ivarB)=dielt_stress(ivarA,ivarB)+&
&     piezo_relaxed(ii1,ivarA)*d_tensor(ii1,ivarB)
    end do
   end do
  end do

! then combine the relaxed ion(fixed strain) dielectric
! tensor and also restore the unit
  do ivarA=1,3
   do ivarB=1,3
    dielt_stress(ivarA,ivarB)=dielt_rlx(ivarA,ivarB)+&
&    dielt_stress(ivarA,ivarB)/(8.854187817)
   end do
  end do

! DEBUG
! write(message,'(a,a,a,a)')ch10,'debug the free stress dielectric tensor ',&
! &  '(unit:pc/N)',ch10
! call wrtout(06,message,'COLL')
! do ivarA=1,3
! write(06,'(3f16.8)')dielt_stress(ivarA,1),dielt_stress(ivarA,2),&
! &    dielt_stress(ivarA,3)
! end do
! ENDDEBUG

! then get the g tensor
  beta_tensor(:,:)=0
  beta_tensor(:,:)=dielt_stress(:,:)
! DEBUG
! write(6,*)' call matrginv 2 '
! ENDDEBUG

  call matrginv(beta_tensor,3,3)
  do ivarA=1,3
   do ivarB=1,6
    g_tensor(ivarA,ivarB)=zero
    do ii1=1,3
     g_tensor(ivarA,ivarB)=g_tensor(ivarA,ivarB)+beta_tensor(ivarA,ii1)*&
&     d_tensor(ivarB,ii1)
    end do
   end do
  end do
! then restore the unit to be m^2/C
  do ivarA=1,3
   do ivarB=1,6
    g_tensor(ivarA,ivarB)=g_tensor(ivarA,ivarB)/(8.854187817)
   end do
  end do
 end if
!then print out the final results of the g tensors(relaxed ion)
 if(anaddb_dtset%piezoflag==5 .or. anaddb_dtset%piezoflag==7)then
  if(anaddb_dtset%instrflag==0 .or. anaddb_dtset%elaflag==0&
&  .or.  anaddb_dtset%elaflag==1 .or.  anaddb_dtset%elaflag==0&
&  .or. anaddb_dtset%dieflag==2 .or. anaddb_dtset%dieflag==1)then
   write(message,'(a,a,a,a,a,a,a,a)' )ch10,&
&   ' WARNING:in order to get the piezoelectric g tensor(relaxed ion),',ch10,&
&   ' need internal strain, dielectric(relaxed-ion) and elastic(realxed ion)',ch10,&
&   ' please set instrflag==1, elaflag==2,3,4 or 5, dieflag==3 or 4',ch10,&
&   ' otherwise the program will still continue but give wrong values.'
   call wrtout(6,message,'coll')
   call wrtout(iout,message,'coll')
  end if

  write(message,'(3a)')ch10,&
&  ' Piezoelectric g tensor (relaxed ion) (unit:m^2/c)',ch10
  call wrtout(06,message,'COLL')
  do ivarA=1,6
   write(06,'(3f16.8)')g_tensor(1,ivarA),g_tensor(2,ivarA),g_tensor(3,ivarA)
  end do
  call wrtout(iout,message,'COLL')
  do ivarA=1,6
   write(iout,'(3f16.8)')g_tensor(1,ivarA),g_tensor(2,ivarA),g_tensor(3,ivarA)
  end do
 end if
!end the part of piezoelectric g tensor (relaxed ion only).

!then start the part for computation of h tensor
!first make sure the dielt_rlx is not zero in the memory
 if(anaddb_dtset%dieflag>0)then
  beta_tensor(:,:)=0
  beta_tensor(:,:)=dielt_rlx(:,:)
! DEBUG
! write(6,*)' call matrginv 3, dielt_rlx(:,:)= ',dielt_rlx(:,:)
! ENDDEBUG

  call matrginv(beta_tensor,3,3)
  do ivarA=1,3
   do ivarB=1,6
    h_tensor(ivarA,ivarB)=zero
    do ii1=1,3
     h_tensor(ivarA,ivarB)=h_tensor(ivarA,ivarB)+beta_tensor(ivarA,ii1)*&
&     piezo_relaxed(ivarB,ii1)
    end do
   end do
  end do
! then restore the unit to be N/c
  do ivarA=1,3
   do ivarB=1,6
    h_tensor(ivarA,ivarB)=1000.0*(h_tensor(ivarA,ivarB)/(8.854187817))
   end do
  end do
 end if
!then print out the final results of h tensors
 if(anaddb_dtset%piezoflag==6 .or. anaddb_dtset%piezoflag==7)then
  if(anaddb_dtset%instrflag==0 .or. anaddb_dtset%dieflag==1 .or. &
&  anaddb_dtset%dieflag==2)then
   write(message,'(a,a,a,a,a,a,a,a)' )ch10,&
&   ' WARNING: in order to get the h tensor, ',ch10,&
&   ' one needs information about internal strain and dielectric(relaxed ion)',ch10,&
&   ' one should set dieflag==3 or 4 and instrflag==1;',ch10,&
&   ' otherwise the program will continue but give wrong values.'
   call wrtout(6,message,'coll')
   call wrtout(iout,message,'coll')
  end if
  write(message,'(3a)')ch10,&
&  ' Piezoelectric h tensor (relaxed ion) (unit:GN/c)',ch10
  call wrtout(06,message,'COLL')
  do ivarA=1,6
   write(06,'(3f16.8)')h_tensor(1,ivarA),h_tensor(2,ivarA),h_tensor(3,ivarA)
  end do
  call wrtout(iout,message,'COLL')
  do ivarA=1,6
   write(iout,'(3f16.8)')h_tensor(1,ivarA),h_tensor(2,ivarA),h_tensor(3,ivarA)
  end do
 end if
!end the part of piezoelectric h tensor (relaxed ion only).

!start to compute the free stress dielectric tensor
 if(anaddb_dtset%dieflag==4)then
  write(message, '(a,a)')ch10,'************************************************'
  call wrtout(6,message,'coll')
  call wrtout(iout,message,'coll')
  if(anaddb_dtset%instrflag==0 .or. anaddb_dtset%elaflag==0 .or. &
&  anaddb_dtset%elaflag==1)then
   write(message,'(a,a,a,a,a,a,a,a)' )ch10,&
&   ' WARNING: in order to get the free stress dielectric tensor,',ch10,&
&   ' one needs internal strain and elastic (relaxed ion)', ch10,&
&   ' we need set elaflag==2,3,4 or 5 and instrflag==1.',ch10,&
&   ' otherwise the program may continue but give wrong and nonsense values.'
   call wrtout(6,message,'coll')
   call wrtout(iout,message,'coll')
  end if
  write(message,'(a,a,a)')ch10,&
&  ' Free stress dielectric tensor (dimensionless)',ch10
  call wrtout(06,message,'COLL')
  do ivarA=1,3
   write(06,'(3f16.8)')dielt_stress(ivarA,1),dielt_stress(ivarA,2),&
&   dielt_stress(ivarA,3)
  end do
  call wrtout(iout,message,'COLL')
  do ivarA=1,3
   write(iout,'(3f16.8)')dielt_stress(ivarA,1),dielt_stress(ivarA,2),&
&   dielt_stress(ivarA,3)
  end do
 end if
!end the part of printing out the free stress dielectric tensor

!then begin the computation of the fixed displacement
!elastic and compliance tensor(relaxed ion)
 do ivarA=1,6
  do ivarB=1,6
   elast_dis(ivarA,ivarB)=zero
   do ii1=1,3
    elast_dis(ivarA,ivarB)=elast_dis(ivarA,ivarB)+&
&    h_tensor(ii1,ivarA)*piezo_relaxed(ivarB,ii1)
   end do
  end do
 end do
!then should add the relaxed ion fixed E-field values
 do ivarA=1,6
  do ivarB=1,6
   elast_dis(ivarA,ivarB)=elast_dis(ivarA,ivarB)+&
&   elast(ivarA,ivarB)
  end do
 end do
!then print out the fixed displacement elastic tensor
 if(anaddb_dtset%elaflag==4)then
  if(anaddb_dtset%instrflag==0 .or. anaddb_dtset%dieflag==1 .or. &
&  anaddb_dtset%dieflag==2 .or. anaddb_dtset%dieflag==0)then
   write(message,'(a,a,a,a,a,a,a,a)' )ch10,&
&   ' WARNING: in order to get the elatic(fixed D field) tensor, ',ch10,&
&   ' one needs information about internal strain and dielectric(relaxed ion)',ch10,&
&   ' one should set dieflag==3 or 4 and instrflag==1;',ch10,&
&   ' otherwise the program will continue but give wrong values.'
   call wrtout(6,message,'coll')
   call wrtout(iout,message,'coll')
  end if
  write(message, '(a,a)')ch10,'************************************************'
  call wrtout(6,message,'coll')
  call wrtout(iout,message,'coll')
  write(message,'(5a)')ch10,&
&  ' Elastic Tensor (relaxed ion) (unit:10^2GP)',ch10,&
&  '  (at fixed displacement field boundary condition)',ch10
  call wrtout(06,message,'COLL')
  do ivarA=1,6
   write(06,'(6f12.7)')elast_dis(ivarA,1)/100.00_dp,elast_dis(ivarA,2)/100.00_dp,&
&   elast_dis(ivarA,3)/100.00_dp,elast_dis(ivarA,4)/100.00_dp,&
&   elast_dis(ivarA,5)/100.00_dp,elast_dis(ivarA,6)/100.00_dp
  end do
  call wrtout(iout,message,'COLL')
  do ivarA=1,6
   write(iout,'(6f12.7)')elast_dis(ivarA,1)/100.00_dp,elast_dis(ivarA,2)/100.00_dp,&
&   elast_dis(ivarA,3)/100.00_dp,elast_dis(ivarA,4)/100.00_dp,&
&   elast_dis(ivarA,5)/100.00_dp,elast_dis(ivarA,6)/100.00_dp
  end do
! then invert the above to get the corresponding compliance tensor
  compliance_dis(:,:)=0
  compliance_dis(:,:)=elast_dis(:,:)
! DEBUG
! write(6,*)' call matrginv 4 '
! ENDDEBUG

  call matrginv(compliance_dis,6,6)
! then print out the compliance tensor at fixed displacement field
  write(message,'(5a)')ch10,&
&  ' Compliance  Tensor (relaxed ion) (unit: 10^-2(GP)^-1)',ch10,&
&  '  (at fixed displacement field boundary condition)',ch10
  call wrtout(06,message,'COLL')
  do ivarB=1,6
   write(6,'(6f12.7)')compliance_dis(ivarB,1)*100.00_dp,&
&   compliance_dis(ivarB,2)*100.00_dp,&
&   compliance_dis(ivarB,3)*100.00_dp,compliance_dis(ivarB,4)*100.00_dp,&
&   compliance_dis(ivarB,5)*100.00_dp,&
&   compliance_dis(ivarB,6)*100.00_dp
  end do
  call wrtout(iout,message,'COLL')

  do ivarB=1,6
   write(iout,'(6f12.7)')compliance_dis(ivarB,1)*100.00,&
&   compliance_dis(ivarB,2)*100.00_dp,&
&   compliance_dis(ivarB,3)*100.00_dp,compliance_dis(ivarB,4)*100.00_dp,&
&   compliance_dis(ivarB,5)*100.00_dp,&
&   compliance_dis(ivarB,6)*100.00_dp
  end do
 end if
!end if the elaflag==4 for the fixed didplacement field elastic tensor
!end the part for computation of elastic at fixed displacement field

!DEBUG
!write(6,*)' piezo9 : exit '
!ENDDEBUG

end subroutine piezo9

!!***
