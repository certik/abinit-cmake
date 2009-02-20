!!{\src2tex{textfont=tt}}
!!****m* ABINIT/ftfvw2
!! NAME
!! ftfvw2
!!
!! FUNCTION
!! Module providing the ability to compute the
!! penalty function and its first derivative associated
!! with some density residuals and the TFvW2 variational charge mixing
!!
!! when trying to solve problem of the form AX = B where A and B are known operators
!! we can minimize the penalty function XAX - 2XB to find the X which fulfilled AX = B
!! the penalty function, pf_,  is  XAX - 2XB and
!! dpf_ is the derivative of the penalty function with respect to X
!! in this function the opA include a laplacian which must be calculated inside the pf and dpf
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!
!!
!! OUTPUT
!!
!!
!! SIDE EFFECTS
!!
!! WARNINGS
!! This is a module.
!! It is made of two functions and one init subroutine
!!
!! NOTES
!!
!! PARENTS
!! prctfw2
!!
!! CHILDREN
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif
module ftfvw2

  use defs_basis
  use defs_datatypes
  use interfaces_11util        ! THIS IS MANDATORY TO CALL dotproduct
  use interfaces_12spacepar
  use interfaces_13recipspace  ! THIS IS MANDATORY TO CALL LAPLACIAN

  implicit none
  !! common variables copied from input
  integer,private                  :: nfft,nspden,ngfft(18)
  real(dp), allocatable,private    :: opA(:,:),opB(:,:),phi(:,:),buffer1(:,:),buffer2(:,:),sqrtrhor(:,:)
  real(dp), private                :: gprimd(3,3),gmet(3,3),gsqcut,qphon(3),Z,ucvol
  type(dataset_type),private       :: dtset
  type(MPI_type),private           :: mpi_enreg
  !! common variables computed
  logical,private                  :: ok=.false.
contains
  !!-------------------------------------------------------------------!!
  !! initialisation subroutine
  !!-------------------------------------------------------------------!!
  !! Copy every variables required for the energy calculation
  !! Allocate the required memory
  !!-------------------------------------------------------------------!!
  subroutine ftfvw2__init(dtset_in,nfft_in,ngfft_in,nspden_in,gprimd_in,gmet_in,gsqcut_in,mpi_enreg_in,&
       & sqrtrhor_in,ucvol_in,Z_in,opA_in,opB_in )


    implicit none
    !Arguments ------------------------------------
    integer,intent(in)  :: nfft_in,ngfft_in(18),nspden_in
    real(dp),intent(in) :: opA_in(nfft_in,nspden_in),opB_in(nfft_in,nspden_in)
    real(dp),intent(in) :: sqrtrhor_in(nfft_in,nspden_in)
    real(dp),intent(in) :: gprimd_in(3,3),gmet_in(3,3)
    real(dp),intent(in) :: gsqcut_in,Z_in,ucvol_in
    type(dataset_type)  ,intent(in)    :: dtset_in
    type(MPI_type),intent(in) :: mpi_enreg_in

    !!allocation and data transfer
    !!Thought it would have been more logical to use the privates intrinsic of the module as
    !!input variables it seems that it is not possible...
    if(.not.ok) then
       dtset = dtset_in
       nspden=nspden_in
       ngfft=ngfft_in
       nfft=nfft_in
       allocate(opA(size(opA_in,1),size(opA_in,2)))
       allocate(opB(size(opB_in,1),size(opB_in,2)))
       allocate(buffer1(size(opB_in,1),size(opB_in,2)))
       allocate(buffer2(size(opB_in,1),size(opB_in,2)))
       allocate(sqrtrhor(size(opB_in,1),size(opB_in,2)))
       opA=opA_in
       opB=opB_in
       sqrtrhor=sqrtrhor_in
       gprimd=gprimd_in
       gmet=gmet_in
       gsqcut=gsqcut_in
       mpi_enreg=mpi_enreg_in
       qphon=zero
       ucvol=ucvol_in
       Z=Z_in
       ok = .true.
    end if
  end subroutine ftfvw2__init

  !!-------------------------------------------------------------------!!
  !! routine to change opB
  !!-------------------------------------------------------------------!!
  !! Copy the new opB into memory after basic check
  !!-------------------------------------------------------------------!!
  subroutine ftfvw2__change_opB(nv1,nv2,newopB)


    integer,intent(in) :: nv1,nv2
    real(dp),intent(in) :: newopB(nv1,nv2)
    if(ok) then
       if(size(newopB,1)==size(opB,1) .and. size(newopB,2)==size(opB,2)) opB=newopB
    end if
  end subroutine ftfvw2__change_opB
  !!-------------------------------------------------------------------!!
  !! ending subroutine
  !!-------------------------------------------------------------------!!
  !! deallocate memory areas
  !!-------------------------------------------------------------------!!
  subroutine ftfvw2__end()
    use defs_basis
    use defs_datatypes

    implicit none
    if(ok) then
       !! set ok to false which prevent using the pf and dpf
       ok = .false.
       !! free memory
       deallocate(opA,opB,buffer1,buffer2,sqrtrhor)
    end if
  end subroutine ftfvw2__end

  !!-------------------------------------------------------------------!!
  !! affectation subroutine
  !!-------------------------------------------------------------------!!
  !! do the required renormalisation when providing a new value for
  !! the density after application of the gradient
  !!-------------------------------------------------------------------!!
  subroutine ftfvw2__newphi2(nv1,nv2,x, grad, vrespc)


    implicit none
    !Arguments ------------------------------------
    integer,intent(in)     :: nv1,nv2
    real(dp),intent(in)    :: x
    real(dp),intent(inout) ::grad(nv1,nv2)
    real(dp),intent(inout) ::vrespc(nv1,nv2)
    grad(:,:)=x*grad(:,:)
    vrespc(:,:)=vrespc(:,:)+grad(:,:)
  end subroutine ftfvw2__newphi2

  !!-------------------------------------------------------------------!!
  !! penalty function associated with the preconditionned residuals    !!
  !!-------------------------------------------------------------------!!
  function ftfvw2__pf(nv1,nv2,phi)



!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12ffts
 use interfaces_13recipspace
 use interfaces_13xc
 use interfaces_lib01cg
!End of the abilint section

    implicit none
    !Arguments ------------------------------------
    integer,intent(in)  :: nv1,nv2
    real(dp),intent(in) ::phi(nv1,nv2)
    real(dp)::ftfvw2__pf
    !Local variables-------------------------------
    real(dp):: phig(2,size(phi,1)),vhphi(size(phi,1))
    integer :: cplex,izero,isign,tim_fourdp
    integer :: ifft,ispden
    real(dp),parameter   :: identity(4)=(/1.0_dp,1.0_dp,0.0_dp,0.0_dp/)

    ! Added by YP [HAVE_FORTRAN_INTERFACES]

    if(ok) then
       buffer1=phi
       tim_fourdp=0
       cplex=1
       isign=-1
       buffer2=phi*sqrtrhor
       call fourdp(cplex,phig,buffer2(:,1),isign,mpi_enreg,nfft,ngfft,dtset%paral_kgb,tim_fourdp)
       call hartre(cplex,gmet,gsqcut,izero,mpi_enreg,nfft,ngfft,dtset%paral_kgb,qphon,phig,vhphi)


       call laplacian(gprimd,mpi_enreg,nfft,nspden,ngfft,dtset%paral_kgb,rdfuncr=buffer1,laplacerdfuncr=buffer2)

       !if(dotproduct(opA(:,:)*phi(:,:)-half*buffer2(:,:),phi) .lt. zero )then
       !   write(0,*) 'dotprod is negative',dotproduct(opA(:,:)*phi(:,:)-half*buffer2(:,:),phi)
       !end if

       do ispden=1,nspden
          buffer2(:,ispden)=((opA(:,ispden))*phi(:,ispden)&
               & -half*buffer2(:,ispden))&
               & -two*opB(:,ispden) &
               & +two*vhphi(:)*identity(ispden)*sqrtrhor(:,ispden)
       end do



       !do ifft=1,nfft
       !   do ispden=1,nspden
       !      if(phi(ifft,ispden) < zero) buffer2(ifft,ispden)=buffer2(ifft,ispden)+phi(ifft,ispden)**2*1000._dp
       !   end do
       !end do
       buffer1=zero
       buffer1(:,1)=vhphi(:)
       ftfvw2__pf=dotproduct(nv1,nv2,phi,buffer2)!+half*dotproduct(phi*sqrtrhor,buffer1)

    else
       ftfvw2__pf=zero
    end if
  end function ftfvw2__pf


  !!-------------------------------------------------------------------!!
  !! derivative of the penalty function
  !!-------------------------------------------------------------------!!
  function ftfvw2__dpf(nv1,nv2,phi)



!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12ffts
 use interfaces_12spacepar
 use interfaces_13recipspace
 use interfaces_13xc
!End of the abilint section

    implicit none
    !Arguments ------------------------------------
    integer,intent(in)  :: nv1,nv2
    real(dp),intent(in) :: phi(nv1,nv2)
    real(dp) :: ftfvw2__dpf(nv1,nv2)
    !Local variables-------------------------------
    real(dp):: phig(2,size(phi,1)),vhphi(size(phi,1)),dotr,doti
    integer :: cplex,izero,isign,tim_fourdp
    integer :: ifft,ispden,option
    real(dp),parameter   :: identity(4)=(/1.0_dp,1.0_dp,0.0_dp,0.0_dp/)
    if(ok) then
       buffer1=phi
       tim_fourdp=0
       cplex=1
       isign=-1
       buffer2=phi*sqrtrhor
       call fourdp(cplex,phig,buffer2(:,1),isign,mpi_enreg,nfft,ngfft,dtset%paral_kgb,tim_fourdp)
       call hartre(cplex,gmet,gsqcut,izero,mpi_enreg,nfft,ngfft,dtset%paral_kgb,qphon,phig,vhphi)

       call laplacian(gprimd,mpi_enreg,nfft,nspden,ngfft,dtset%paral_kgb,rdfuncr=buffer1,laplacerdfuncr=buffer2)
       do ispden=1,nspden
          ftfvw2__dpf(:,ispden)= (two*(opA(:,ispden)*phi(:,ispden)-opB(:,ispden))&
               & - buffer2(:,ispden)) &
               & + four*sqrtrhor(:,ispden)*vhphi(:)*identity(ispden)


       end do

       !! We keep only the part of phi which is orthogonal to sqrtrhor
       !! In order to do that we must use a normed version of sqrtrhor
       !write(0,*) 'avant le projection'
       option=1
       call dotprod_vn(cplex,& !complex density/pot
            &ftfvw2__dpf(:,:),&          !the density
            &dotr,&  !resulting dorproduct integrated over r
            &doti,&          !imaginary part of the integral
            &mpi_enreg,&     !
            &nfft,&          !number of localy(cpu) attributed grid point
            &nfft,&        !real total number of grid point
            &nspden,&        !nspden
            &option,&        !1=compute only the real part 2=compute also the imaginary part
            &sqrtrhor(:,:),&          !the potential
            &ucvol)          !cell volume
      ftfvw2__dpf=ftfvw2__dpf-sqrtrhor*dotr/Z
       !write(0,*) 'apres la projection',dotr,Z




    else
       ftfvw2__dpf = zero
    end if
  end function ftfvw2__dpf


end module ftfvw2
!!***
