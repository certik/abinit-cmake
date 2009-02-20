!{\src2tex{textfont=tt}}
!!****f* ABINIT/rhoij_alloc
!! NAME
!! rhoij_alloc
!!
!! FUNCTION
!! This module contains functions used to manipulate
!! variables of structured datatype pawrhoij_type.
!! pawrhoij_type variables are rhoij occupancies
!! matrix used within PAW formalism
!!
!! COPYRIGHT
!! Copyright (C) 2007-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine rhoij_alloc(cplex,nlmn,nspden,nsppol,pawrhoij,typat,&      ! Mandatory arguments
&                      ngrhoij,nlmnmix,use_rhoij_,use_rhoijres) ! Optional arguments

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nspden,nsppol
 integer,intent(in),optional :: ngrhoij,nlmnmix,use_rhoij_,use_rhoijres
!arrays
 integer,intent(in) :: nlmn(:),typat(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer :: irhoij,itypat,lmn2_size,nn1,nn2,nrhoij

! *************************************************************************

 nrhoij=size(pawrhoij);nn1=size(nlmn);nn2=size(typat)
 if (nrhoij/=nn2.or.maxval(typat)>nn1) stop "Error in rhoij_alloc: wrong sizes ! "

 do irhoij=1,nrhoij

  itypat=typat(irhoij)
  lmn2_size=nlmn(itypat)*(nlmn(itypat)+1)/2

! Scalars initializations
  pawrhoij(irhoij)%cplex=cplex
  pawrhoij(irhoij)%lmn_size=nlmn(itypat)
  pawrhoij(irhoij)%lmn2_size=lmn2_size
  pawrhoij(irhoij)%nspden=nspden
  pawrhoij(irhoij)%nsppol=nsppol
  pawrhoij(irhoij)%nrhoijsel=0
  pawrhoij(irhoij)%lmnmix_sz=0
  pawrhoij(irhoij)%ngrhoij=0
  pawrhoij(irhoij)%use_rhoij_=0
  pawrhoij(irhoij)%use_rhoijres=0

! Mandatory pointers allocations
  allocate(pawrhoij(irhoij)%rhoijselect(lmn2_size))
  allocate(pawrhoij(irhoij)%rhoijp(cplex*lmn2_size,nspden))

! Optional pointers allocations
  if (present(ngrhoij)) then
   if (ngrhoij>0) then
    pawrhoij(irhoij)%ngrhoij=ngrhoij
    allocate(pawrhoij(irhoij)%grhoij(cplex*ngrhoij,lmn2_size,nspden))
   end if
  end if
  if (present(nlmnmix)) then
   if (nlmnmix>0) then
    pawrhoij(irhoij)%lmnmix_sz=nlmnmix
    allocate(pawrhoij(irhoij)%kpawmix(nlmnmix))
   end if
  end if
  if (present(use_rhoij_)) then
   if (use_rhoij_>0) then
    pawrhoij(irhoij)%use_rhoij_=use_rhoij_
    allocate(pawrhoij(irhoij)%rhoij_(cplex*lmn2_size,nspden))
   end if
  end if
  if (present(use_rhoijres)) then
   if (use_rhoijres>0) then
    pawrhoij(irhoij)%use_rhoijres=use_rhoijres
    allocate(pawrhoij(irhoij)%rhoijres(cplex*lmn2_size,nspden))
   end if
  end if

! Intializations to zero (to avoid overflow)
  pawrhoij(irhoij)%rhoijselect(:)=0
  pawrhoij(irhoij)%rhoijp(:,:)=zero

 end do

end subroutine rhoij_alloc
!!***

!!****f* ABINIT/rhoij_free
!! NAME
!! rhoij_free
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2007-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine rhoij_free(pawrhoij)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!arrays
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer :: irhoij,nrhoij

! *************************************************************************

 nrhoij=size(pawrhoij)
 do irhoij=1,nrhoij
  if (associated(pawrhoij(irhoij)%rhoijp)) deallocate(pawrhoij(irhoij)%rhoijp)
  if (associated(pawrhoij(irhoij)%rhoijselect)) deallocate(pawrhoij(irhoij)%rhoijselect)
  if (pawrhoij(irhoij)%ngrhoij>0) deallocate(pawrhoij(irhoij)%grhoij)
  if (pawrhoij(irhoij)%lmnmix_sz>0) deallocate(pawrhoij(irhoij)%kpawmix)
  if (pawrhoij(irhoij)%use_rhoij_>0) deallocate(pawrhoij(irhoij)%rhoij_)
  if (pawrhoij(irhoij)%use_rhoijres>0) deallocate(pawrhoij(irhoij)%rhoijres)
 end do

end subroutine rhoij_free
!!***

!!****f* ABINIT/rhoij_copy
!! NAME
!! rhoij_copy
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2007-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine rhoij_copy(pawrhoij_in,pawrhoij_out)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!arrays
 type(pawrhoij_type),intent(in) :: pawrhoij_in(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij_out(:)

!Local variables-------------------------------
!scalars
 integer :: cplex,ilmn,irhoij,ispden,lmn2_size,lmnmix,ngrhoij,nrhoij_in
 integer :: nrhoij_out,nselect,nspden,use_rhoij_,use_rhoijres

! *************************************************************************

 nrhoij_in=size(pawrhoij_in);nrhoij_out=size(pawrhoij_out)
 if (nrhoij_in/=nrhoij_out) stop "Error in rhoij_copy: wrong sizes ! "

 do irhoij=1,nrhoij_in

  cplex=pawrhoij_in(irhoij)%cplex
  lmn2_size=pawrhoij_in(irhoij)%lmn2_size
  nspden=pawrhoij_in(irhoij)%nspden

! Scalars
  pawrhoij_out(irhoij)%lmn_size=pawrhoij_in(irhoij)%lmn_size
  pawrhoij_out(irhoij)%nspden=nspden
  pawrhoij_out(irhoij)%nsppol=pawrhoij_in(irhoij)%nsppol
  pawrhoij_out(irhoij)%nrhoijsel=pawrhoij_in(irhoij)%nrhoijsel

! Mandatory pointers
  if (pawrhoij_out(irhoij)%lmn2_size/=lmn2_size.or. &
&  pawrhoij_out(irhoij)%cplex/=cplex) then
   pawrhoij_out(irhoij)%lmn2_size=lmn2_size
   deallocate(pawrhoij_out(irhoij)%rhoijselect)
   deallocate(pawrhoij_out(irhoij)%rhoijp)
   allocate(pawrhoij_out(irhoij)%rhoijselect(lmn2_size))
   allocate(pawrhoij_out(irhoij)%rhoijp(cplex*lmn2_size,nspden))
  end if
  nselect=pawrhoij_in(irhoij)%nrhoijsel+0
  pawrhoij_out(irhoij)%rhoijselect(1:nselect)=pawrhoij_in(irhoij)%rhoijselect(1:nselect)+0
  do ispden=1,nspden
   pawrhoij_out(irhoij)%rhoijp(1:cplex*nselect,ispden)=pawrhoij_in(irhoij)%rhoijp(1:cplex*nselect,ispden)+zero
  end do
  if (nselect<lmn2_size) then
   pawrhoij_out(irhoij)%rhoijselect(nselect+1:lmn2_size)=0
   do ispden=1,nspden
    pawrhoij_out(irhoij)%rhoijp(cplex*nselect+1:cplex*lmn2_size,ispden)=zero
   end do
  end if

! Optional pointers
  ngrhoij=pawrhoij_in(irhoij)%ngrhoij
  if (pawrhoij_out(irhoij)%ngrhoij/=ngrhoij.or. &
&  pawrhoij_out(irhoij)%cplex/=cplex) then
   if (pawrhoij_out(irhoij)%ngrhoij>0) deallocate(pawrhoij_out(irhoij)%grhoij)
   if (ngrhoij>0) allocate(pawrhoij_out(irhoij)%grhoij(ngrhoij,lmn2_size,nspden))
   pawrhoij_out(irhoij)%ngrhoij=ngrhoij
  end if
  if (ngrhoij>0) then
   do ispden=1,nspden
    do ilmn=1,lmn2_size
     pawrhoij_out(irhoij)%grhoij(1:cplex*ngrhoij,ilmn,ispden)=pawrhoij_in(irhoij)%grhoij(1:cplex*ngrhoij,ilmn,ispden)
    end do
   end do
  end if
! ---
  lmnmix=pawrhoij_in(irhoij)%lmnmix_sz
  if (pawrhoij_out(irhoij)%lmnmix_sz/=lmnmix) then
   if (pawrhoij_out(irhoij)%lmnmix_sz>0) deallocate(pawrhoij_out(irhoij)%kpawmix)
   if (lmnmix>0) allocate(pawrhoij_out(irhoij)%kpawmix(lmnmix))
   pawrhoij_out(irhoij)%lmnmix_sz=lmnmix
  end if
  if (lmnmix>0) then
   pawrhoij_out(irhoij)%kpawmix(1:lmnmix)=pawrhoij_in(irhoij)%kpawmix(lmnmix)
  end if
! ---
  use_rhoij_=pawrhoij_in(irhoij)%use_rhoij_
  if (pawrhoij_out(irhoij)%use_rhoij_/=use_rhoij_) then
   if (pawrhoij_out(irhoij)%use_rhoij_>0) deallocate(pawrhoij_out(irhoij)%rhoij_)
   if (use_rhoij_>0) allocate(pawrhoij_out(irhoij)%rhoij_(lmn2_size,nspden))
   pawrhoij_out(irhoij)%use_rhoij_=use_rhoij_
  end if
  if (use_rhoij_>0) then
   do ispden=1,nspden
    pawrhoij_out(irhoij)%rhoij_(1:cplex*lmn2_size,ispden)=pawrhoij_in(irhoij)%rhoij_(1:cplex*lmn2_size,ispden)
   end do
  end if
! ---
  use_rhoijres=pawrhoij_in(irhoij)%use_rhoijres
  if (pawrhoij_out(irhoij)%use_rhoijres/=use_rhoijres) then
   if (pawrhoij_out(irhoij)%use_rhoijres>0) deallocate(pawrhoij_out(irhoij)%rhoijres)
   if (use_rhoijres>0) allocate(pawrhoij_out(irhoij)%rhoijres(lmn2_size,nspden))
   pawrhoij_out(irhoij)%use_rhoijres=use_rhoijres
  end if
  if (use_rhoijres>0) then
   do ispden=1,nspden
    pawrhoij_out(irhoij)%rhoijres(1:cplex*lmn2_size,ispden)=pawrhoij_in(irhoij)%rhoijres(1:cplex*lmn2_size,ispden)
   end do
  end if

  pawrhoij_out(irhoij)%cplex=cplex

 end do ! irhoij

end subroutine rhoij_copy
!!***
