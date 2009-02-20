!{\src2tex{textfont=tt}}
!!****f* ABINIT/getghc
!!
!! NAME
!! getghc
!!
!! FUNCTION
!! Compute <G|H|C> for input vector |C> expressed in reciprocal space
!! Result is put in array ghc.
!! <G|Vnonlocal|C> is also returned in gvnlc.
!! if required, <G|S|C> is returned in gsc (S=overlap - PAW only)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, LSI, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cwavef(2,npw*nspinor*ndat)=planewave coefficients of wavefunction.
!! dimffnl=second dimension of ffnl (1+number of derivatives)
!! ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!! filstat=name of the status file
!! gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian to be applied
!! kg_k(3,npw)=G vec coordinates wrt recip lattice transl.
!! kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!! lambda=factor to be used when computing <G|H-lambda.S|C> - only for sij_opt=-1
!!        Typically lambda is the eigenvalue (or its guess)
!! lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!       =if useylm=0, max number of (l,n)   comp. over all type of psps
!! matblk=dimension of the array ph3d
!! mgfft=maximum size for 1D FFTs
!! mpi_enreg=informations about MPI parallelization
!! mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!! mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!! natom=number of atoms in unit cell.
!! ndat=number of FFT to do in //
!! npw=number of planewaves in basis for given k point.
!! nspinor=number of spinorial components of the wavefunctions
!! ntypat=number of types of atoms in cell.
!! nvloc=final dimension of vlocal (usually 1, but 4 for non-collinear)
!! n4,n5,n6 used for dimensionning of vlocal
!! ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!! prtvol=control print volume and debugging output
!! sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H|C> have to be computed
!!    (S=overlap)       if  1, matrix elements <G|S|C> have to be computed in gsc in addition to ghc
!!                      if -1, matrix elements <G|H-lambda.S|C> have to be computed in ghc (gsc not used)
!! tim_getghc=timing code of the calling subroutine(can be set to 0 if not attributed)
!! type_calc= option governing which part of Hamitonian is to be applied:
!             0: whole Hamiltonian
!!            1: local part only
!!            2: non-local+kinetic only (added to the exixting Hamiltonian)
!! vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
!!
!! OUTPUT
!!  ghc(2,npw*nspinor*ndat)=matrix elements <G|H|C> (if sij_opt>=0)
!!                                       or <G|H-lambda.S|C> (if sij_opt=-1)
!!  gvnlc(2,npw*nspinor*ndat)=matrix elements <G|Vnonlocal|C> (if sij_opt>=0)
!!                                         or <G|Vnonlocal-lambda.S|C> (if sij_opt=-1)
!!    (sometimes desired for computing nonlocal part of total energy, but can be ignored).
!!  if (sij_opt=1)
!!    gsc(2,npw*nspinor*ndat)=matrix elements <G|S|C> (S=overlap).
!!
!! PARENTS
!!      cgwf,cgwf3,lobpcgIIwf,lobpcgccIIwf,lobpcgccwf,lobpcgwf,mkresi,prep_getghc
!!
!! CHILDREN
!!      fourwf,leave_new,nonlop,status,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine getghc(cwavef,dimffnl,ffnl,filstat,ghc,gsc,gs_ham,&
&  gvnlc,kg_k,kinpw,lambda,lmnmax,&
&  matblk,mgfft,mpi_enreg,mpsang,mpssoang,&
&  natom,ndat,npw,nspinor,ntypat,nvloc,n4,n5,n6,&
&  paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,type_calc,vlocal)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_13nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!This type is defined in defs_mpi
!scalars
 integer,intent(in) :: dimffnl,lmnmax,matblk,mgfft,mpsang,mpssoang,n4,n5,n6
 integer,intent(in) :: natom,ndat,npw,nspinor,ntypat,nvloc,paral_kgb,prtvol
 integer,intent(in) :: sij_opt,tim_getghc,type_calc
 real(dp) :: lambda
 character(len=fnlen),intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
 type(gs_hamiltonian_type),intent(in) :: gs_ham
!arrays
 integer,intent(in) :: kg_k(3,npw)
 real(dp),intent(in) :: ffnl(npw,dimffnl,lmnmax,ntypat),kinpw(npw)
 real(dp),intent(inout) :: cwavef(2,npw*nspinor*ndat),ghc(2,npw*nspinor*ndat)
 real(dp),intent(inout) :: gvnlc(2,npw*nspinor*ndat),ph3d(2,npw,matblk)
 real(dp),intent(inout) :: vlocal(n4,n5,n6,nvloc)
 real(dp),intent(out) :: gsc(2,npw*nspinor*ndat*(sij_opt+1)/2)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,level=10,re=1
 integer :: choice,cplex,cpopt,i1,i2,i3,idat,idir,iexit,ig,igspinor,ii,ipw,ir
 integer :: ispinor,nkpg,nnlout,paw_opt,signs,tim_fourwf,tim_nonlop
 real(dp) :: dum,enlk_direct,ghcim,ghcre,weight
 character(len=500) :: message
!arrays
 real(dp) :: enlout(1),nonlop_dum(1,1),tsec(2)
 real(dp),allocatable :: cwavef_nonlop(:,:),cwavef_sp(:,:),ghc1(:,:),ghc2(:,:)
 real(dp),allocatable :: ghc3(:,:),ghc4(:,:),ghc_sp(:,:),gsc_nonlop(:,:)
 real(dp),allocatable :: gvnlc_nonlop(:,:),kpg_dum(:,:),vlocal_tmp(:,:,:)
 real(dp),allocatable :: work(:,:,:,:)
 type(cprj_type) :: cwaveprj(1,1)

! *********************************************************************

!Keep track of total time spent in getghc:
 call timab(200+tim_getghc,1,tsec)

!DEBUG
!write(6,*)' enter getghc '
!write(6,*)' getghc : cwavef(:,1)=',cwavef(:,1)
!stop
!ENDDEBUG

 if(prtvol<0)then
  call status(0,filstat,iexit,level,'enter         ')
 end if

!Structured debugging if prtvol==-level
 if(prtvol==-level)then
  write(message,'(80a,a,a)') ('=',ii=1,80),ch10,&
&  ' getghc : enter, debugging '
  call wrtout(06,message,'PERS')
 end if

 if ((type_calc==0).or.(type_calc==1)) then

  allocate(work(2,n4,n5,n6*ndat))

! Apply the local potential to the wavefunction
! Start from wavefunction in reciprocal space cwavef
! End with function ghc in reciprocal space also.
  weight=1.0_dp

! DEBUG
! write(6,*)' getghc : will call fourwf '
! write(6,*)' vlocal='
! do i1=1,n4
! do i2=1,n5
! do i3=1,n6
! write(6,*)i1,i2,i3,vlocal(i1,i2,i3,1)
! end do
! end do
! end do
! write(6,*)' cwavef='
! do ii=1,npw
! write(6,*)ii,cwavef(1,ii),cwavef(2,ii)
! end do
! ENDDEBUG

! Application of the local potential
  tim_fourwf=1

! Treat scalar local potentials
  if(nvloc == 1) then

   call fourwf(1,vlocal,cwavef,ghc,work,gs_ham%gbound,gs_ham%gbound,&
&   gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&   npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight)
   if(nspinor==2)then
    allocate(cwavef_sp(2,npw),ghc_sp(2,npw))
!   $OMP PARALLEL DO PRIVATE(ipw) &
!   $OMP&SHARED(cwavef,cwavef_sp,npw)
    do ipw=1,npw
     cwavef_sp(1,ipw)=cwavef(1,ipw+npw)
     cwavef_sp(2,ipw)=cwavef(2,ipw+npw)
    end do
!   $OMP END PARALLEL DO
    call fourwf(1,vlocal,cwavef_sp,ghc_sp,work,gs_ham%gbound,gs_ham%gbound,&
&    gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&    npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight)
!   $OMP PARALLEL DO PRIVATE(ipw) &
!   $OMP&SHARED(ghc,ghc_sp,npw)
    do ipw=1,npw
     ghc(1,ipw+npw)=ghc_sp(1,ipw)
     ghc(2,ipw+npw)=ghc_sp(2,ipw)
    end do
!   $OMP END PARALLEL DO
    deallocate(cwavef_sp,ghc_sp)
   end if
!  Treat non-collinear local potentials
  else if(nvloc==4) then

   allocate(cwavef_sp(2,npw))
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(cwavef,cwavef_sp,npw)
   do ipw=1,npw
    cwavef_sp(1,ipw)=cwavef(1,ipw+npw)
    cwavef_sp(2,ipw)=cwavef(2,ipw+npw)
   end do
!  $OMP END PARALLEL DO
   allocate(ghc1(2,npw),ghc2(2,npw),ghc3(2,npw),ghc4(2,npw))

   allocate(vlocal_tmp(n4,n5,n6))
!  v11*phi1=ghc1
   vlocal_tmp(:,:,:)=vlocal(:,:,:,1)
   call fourwf(1,vlocal_tmp,cwavef,ghc1,work,gs_ham%gbound,gs_ham%gbound,&
&   gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&   npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight)
!  v22*phi2=ghc2
   vlocal_tmp(:,:,:)=vlocal(:,:,:,2)
   call fourwf(1,vlocal_tmp,cwavef_sp,ghc2,work,gs_ham%gbound,gs_ham%gbound,&
&   gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&   npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight)
   deallocate(vlocal_tmp)

   cplex=2
   allocate(vlocal_tmp(cplex*n4,n5,n6))
!  (re(v12)-im(v12))*phi1=ghc3
   do i3=1,n6
    do i2=1,n5
     do i1=1,n4
      vlocal_tmp(2*i1-1,i2,i3)= vlocal(i1,i2,i3,3)
      vlocal_tmp(2*i1  ,i2,i3)=-vlocal(i1,i2,i3,4)
     end do
    end do
   end do
   call fourwf(cplex,vlocal_tmp,cwavef,ghc3,work,gs_ham%gbound,gs_ham%gbound,&
&   gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&   npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight)
!  (re(v12)+im(v12))*phi2=ghc4
   do i3=1,n6
    do i2=1,n5
     do i1=1,n4
      vlocal_tmp(2*i1,i2,i3)=-vlocal_tmp(2*i1,i2,i3)
     end do
    end do
   end do
   call fourwf(cplex,vlocal_tmp,cwavef_sp,ghc4,work,gs_ham%gbound,gs_ham%gbound,&
&   gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&   npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight)
!  Construct ghc from pieces
!  (v11,v22,Re(v12)+iIm(v12);Re(v12)-iIm(v12))(psi1;psi2): matrix product
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(ghc,ghc1,ghc2,ghc3,ghc4,npw)
   do ipw=1,npw
    ghc(1,ipw)=ghc1(1,ipw)+ghc4(1,ipw)
    ghc(2,ipw)=ghc1(2,ipw)+ghc4(2,ipw)
    ghc(1,npw+ipw)=ghc3(1,ipw)+ghc2(1,ipw)
    ghc(2,npw+ipw)=ghc3(2,ipw)+ghc2(2,ipw)
   end do
!  $OMP END PARALLEL DO

   deallocate(ghc1,ghc2,ghc3,ghc4,cwavef_sp,vlocal_tmp)

  end if ! nvloc==1

  deallocate(work)

  if(prtvol<0)then
   call status(0,filstat,iexit,level,'call nonlop   ')
  end if

 end if


 if ((type_calc==0).or.(type_calc==2)) then

! DEBUG
! write(6,*)' getghc: before nonlop'
! ENDDEBUG

  signs=2 ; choice=1 ; nnlout=1 ; idir=0 ; tim_nonlop=1 ; nkpg=0 ; cpopt=-1
  paw_opt=gs_ham%usepaw ; if (sij_opt/=0) paw_opt=sij_opt+3
  if(ndat==1)then
   if(gs_ham%usepaw==0)then
    call nonlop(gs_ham%atindx1,choice,cpopt,cwaveprj,gs_ham%dimekb1,gs_ham%dimekb2,dimffnl,dimffnl,&
&    gs_ham%ekb,enlout,ffnl,ffnl,gs_ham%gmet,gs_ham%gprimd,idir,&
&    gs_ham%indlmn,gs_ham%istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_ham%kpoint,gs_ham%kpoint,&
&    lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_ham%nattyp,&
&    gs_ham%ngfft,nkpg,nkpg,gs_ham%nloalg,nnlout,npw,npw,nspinor,ntypat,0,paw_opt,&
&    gs_ham%phkxred,gs_ham%phkxred,gs_ham%ph1d,ph3d,ph3d,gs_ham%pspso,&
&    signs,nonlop_dum,nonlop_dum,tim_nonlop,gs_ham%ucvol,gs_ham%useylm,cwavef,gvnlc)
   else
    call nonlop(gs_ham%atindx1,choice,cpopt,cwaveprj,gs_ham%dimekb1,gs_ham%dimekb2,dimffnl,dimffnl,&
&    gs_ham%ekb,enlout,ffnl,ffnl,gs_ham%gmet,gs_ham%gprimd,idir,&
&    gs_ham%indlmn,gs_ham%istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_ham%kpoint,gs_ham%kpoint,&
&    lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_ham%nattyp,&
&    gs_ham%ngfft,nkpg,nkpg,gs_ham%nloalg,nnlout,npw,npw,nspinor,ntypat,0,paw_opt,&
&    gs_ham%phkxred,gs_ham%phkxred,gs_ham%ph1d,ph3d,ph3d,gs_ham%pspso,&
&    signs,gs_ham%sij,gsc,tim_nonlop,gs_ham%ucvol,gs_ham%useylm,cwavef,gvnlc)
   end if
  else
   allocate(cwavef_nonlop(2,npw*nspinor),gvnlc_nonlop(2,npw*nspinor))
   if (sij_opt==1) allocate(gsc_nonlop(2,npw*nspinor))
   do idat=1,ndat
    cwavef_nonlop(re,1:npw*nspinor)=cwavef(re,1+npw*nspinor*(idat-1):npw*nspinor*idat)
    cwavef_nonlop(im,1:npw*nspinor)=cwavef(im,1+npw*nspinor*(idat-1):npw*nspinor*idat)
    if(gs_ham%usepaw==0)then
     call nonlop(gs_ham%atindx1,choice,cpopt,cwaveprj,gs_ham%dimekb1,gs_ham%dimekb2,dimffnl,dimffnl,&
&     gs_ham%ekb,enlout,ffnl,ffnl,gs_ham%gmet,gs_ham%gprimd,idir,&
&     gs_ham%indlmn,gs_ham%istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_ham%kpoint,gs_ham%kpoint,&
&     lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_ham%nattyp,&
&     gs_ham%ngfft,nkpg,nkpg,gs_ham%nloalg,nnlout,npw,npw,nspinor,ntypat,0,paw_opt,&
&     gs_ham%phkxred,gs_ham%phkxred,gs_ham%ph1d,ph3d,ph3d,gs_ham%pspso,&
&     signs,nonlop_dum,nonlop_dum,tim_nonlop,gs_ham%ucvol,gs_ham%useylm,cwavef_nonlop,gvnlc_nonlop)
    else
     call nonlop(gs_ham%atindx1,choice,cpopt,cwaveprj,gs_ham%dimekb1,gs_ham%dimekb2,dimffnl,dimffnl,&
&     gs_ham%ekb,enlout,ffnl,ffnl,gs_ham%gmet,gs_ham%gprimd,idir,&
&     gs_ham%indlmn,gs_ham%istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_ham%kpoint,gs_ham%kpoint,&
&     lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_ham%nattyp,&
&     gs_ham%ngfft,nkpg,nkpg,gs_ham%nloalg,nnlout,npw,npw,nspinor,ntypat,0,paw_opt,&
&     gs_ham%phkxred,gs_ham%phkxred,gs_ham%ph1d,ph3d,ph3d,gs_ham%pspso,&
&     signs,gs_ham%sij,gsc_nonlop,tim_nonlop,gs_ham%ucvol,gs_ham%useylm,cwavef_nonlop,gvnlc_nonlop)
    end if
    gvnlc(re,1+npw*nspinor*(idat-1):npw*nspinor*idat)=gvnlc_nonlop(re,1:npw*nspinor)
    gvnlc(im,1+npw*nspinor*(idat-1):npw*nspinor*idat)=gvnlc_nonlop(im,1:npw*nspinor)
    if (sij_opt==1) then
     gsc(re,1+npw*nspinor*(idat-1):npw*nspinor*idat)=gsc_nonlop(re,1:npw*nspinor)
     gsc(im,1+npw*nspinor*(idat-1):npw*nspinor*idat)=gsc_nonlop(im,1:npw*nspinor)
    end if
   end do
   deallocate(cwavef_nonlop,gvnlc_nonlop)
   if (sij_opt==1) deallocate(gsc_nonlop)
  end if
! DEBUG
! write(6,*)' getghc: after nonlop'
! ENDDEBUG

  if(prtvol<0)then
   call status(0,filstat,iexit,level,'assemble      ')
  end if

  if(prtvol==-level)then
   write(message,'(a)')' getghc : components of ghc '
   call wrtout(06,message,'PERS')
   write(message,'(a,a)')' ig re/im     ghc     ',&
&   '   kinpw         cwavef      glocc        gvnlc '
   call wrtout(06,message,'PERS')
  end if

! Assemble modified kinetic, local and nonlocal contributions
! to <G|H|C(n,k)>. Take also into account build-in debugging.
  if(prtvol/=-level)then
   do idat=1,ndat
    do ispinor=1,nspinor
!    $OMP PARALLEL DO PRIVATE(idat,ig,igspinor) &
!    $OMP&SHARED(cwavef,ghc,gvnlc,kinpw,ndat,npw)
     do ig=1,npw
      igspinor=ig+npw*(ispinor-1)+npw*nspinor*(idat-1)
      if(kinpw(ig)<huge(zero)*1.d-11)then
       ghc(re,igspinor)=kinpw(ig)*cwavef(re,igspinor)+ghc(re,igspinor)+gvnlc(re,igspinor)
       ghc(im,igspinor)=kinpw(ig)*cwavef(im,igspinor)+ghc(im,igspinor)+gvnlc(im,igspinor)
      else
       ghc(re,igspinor)=zero
       ghc(im,igspinor)=zero
       if (sij_opt==1) then
        gsc(re,igspinor)=zero
        gsc(im,igspinor)=zero
       end if
      end if
     end do ! ig
    end do ! ispinor
!   $OMP END PARALLEL DO
   end do ! idat
  else
!  Here, debugging section
   do idat=1,ndat
    do ispinor=1,nspinor
!    $OMP PARALLEL DO PRIVATE(ghcre,ghcim,idat,ig,igspinor) &
!    $OMP&SHARED(cwavef,ghc,gvnlc,kinpw,ndat,npw)
     do ig=1,npw
      igspinor=ig+npw*(ispinor-1)+npw*nspinor*(idat-1)
      if(kinpw(ig)<huge(zero)*1.d-11)then
       ghcre=kinpw(ig)*cwavef(re,igspinor)+ghc(re,igspinor)+gvnlc(re,igspinor)
       ghcim=kinpw(ig)*cwavef(im,igspinor)+ghc(im,igspinor)+gvnlc(im,igspinor)
      else
       ghcre=zero
       ghcim=zero
       if (sij_opt==1) then
        gsc(re,igspinor)=zero
        gsc(im,igspinor)=zero
       end if
      end if
      if(ig<=3 .or. mod(ig,50)==0 )then
!      DEBUG
!      if(ig<=25 .or. mod(ig,20)==0 )then
!      ENDDEBUG
       write(message,'(i3,a,5es13.6)')igspinor,'  1  ',ghcre,&
&       kinpw(ig),cwavef(re,igspinor),ghc(re,igspinor),gvnlc(re,igspinor)
       call wrtout(06,message,'PERS')
       write(message,'(a,5es13.6)')'     2  ',ghcim,&
&       kinpw(ig),cwavef(im,igspinor),ghc(im,igspinor),gvnlc(im,igspinor)
       call wrtout(06,message,'PERS')
      end if
      ghc(re,igspinor)=ghcre
      ghc(im,igspinor)=ghcim
     end do ! ig
    end do ! ispinor
!   $OMP END PARALLEL DO
   end do ! idat
  end if

! Structured debugging : if prtvol=-level, stop here.
  if(prtvol==-level)then
   write(message,'(a,a,a,a,i2,a)') ch10,&
&   ' getghc : exit ',&
&   ch10,'  prtvol=-',level,', debugging mode => stop '
   call wrtout(06,message,'PERS')
   call status(0,filstat,iexit,level,'debug => stop ')
   call leave_new('PERS')
  end if

  if(prtvol<0)then
   call status(0,filstat,iexit,level,'exit          ')
  end if

 end if

 call timab(200+tim_getghc,2,tsec)

!DEBUG
!write(6,*)' getghc : exit'
!ENDDEBUG

end subroutine getghc
!!***
