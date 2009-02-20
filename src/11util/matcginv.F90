!{\src2tex{textfont=tt}}
!!****f* ABINIT/matcginv
!! NAME
!! matcginv
!!
!! FUNCTION
!! Invert a general matrix of complex elements.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! lda=leading dimension of complex matrix a
!! n=size of complex matrix a
!! a=matrix of complex elements
!! OUTPUT
!! a=inverse of a input matrix
!! SIDE EFFECTS
!! a(lda,n)= array of complex elements, input, inverted at output
!!
!!
!! PARENTS
!!      calc_ffm,calc_rpa_functional,eps1_tc,make_epsm1_driver
!!
!! CHILDREN
!!      cbgmdi,cbgmlu,cgeicd,cgetrf,cgetri,wrtout,zgetrf,zgetri
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine matcginv(a,lda,n)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lda,n
!arrays
 complex(gwpc),intent(inout) :: a(lda,n)

!Local variables-------------------------------
!scalars
 integer :: ierr,istat,nwork
 character(len=500) :: message
!arrays
 integer :: ipvt(n)
 complex(gwpc),allocatable :: work(:)
!no_abirules
#if defined HAVE_IBM_ESSL_OLD
 complex :: rcond
 complex :: det(2)
#else
 real(dp) :: det
 complex(gwpc) :: cdet
#endif

! *************************************************************************
#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'CGETRI' :: cgetri
!DEC$ ATTRIBUTES ALIAS:'CGETRF' :: cgetrf
#endif

#if defined HAVE_IBM_ESSL_OLD
 nwork=200*n
#else
 nwork=n
#endif

 allocate(work(nwork))

#if defined HAVE_IBM_ESSL_OLD
 call cgeicd(a,lda,n,0,rcond,det,work,nwork)
 if(abs(rcond)==zero) then
  write(message,'(10a)')ch10,&
&  ' matcginv : BUG -',ch10,&
&  '  The matrix that has been passed in argument of this subroutine',ch10,&
&  '  is probably either singular or nearly singular.',ch10,&
&  '  The ESSL routine cgeicd failed.',ch10,&
&  '  Action : Contact ABINIT group '
  call wrtout(6,message,'COLL') ; call leave_new('COLL')
 end if
#elif defined HAVE_NEC_ASL
 call cbgmlu(a,lda,n,ipvt,ierr)
 if(ierr/=0) then
  write(message,'(10a)')ch10,&
&  ' matcginv : BUG -',ch10,&
&  '  The matrix that has been passed in argument of this subroutine',ch10,&
&  '  is probably either singular or nearly singular.',ch10,&
&  '  The ASL routine cbgmlu failed.',ch10,&
&  '  Action : Contact ABINIT group '
  call wrtout(6,message,'COLL') ; call leave_new('COLL')
 end if
 call cbgmdi(a,lda,n,ipvt,cdet,det,-1,work,ierr)
 if(ierr/=0) then
  write(message,'(10a)')ch10,&
&  ' matcginv : BUG -',ch10,&
&  '  The matrix that has been passed in argument of this subroutine',ch10,&
&  '  is probably either singular or nearly singular.',ch10,&
&  '  The ASL routine dbgmdi failed.',ch10,&
&  '  Action : Contact ABINIT group '
  call wrtout(6,message,'COLL') ; call leave_new('COLL')
 end if
#else
#if defined HAVE_GW_DPC
 call zgetrf(n,n,a,lda,ipvt,ierr)
#else
 call cgetrf(n,n,a,lda,ipvt,ierr)
#endif
 if(ierr/=0) then
  write(message,'(10a)')ch10,&
&  ' matcginv : BUG -',ch10,&
&  '  The matrix that has been passed in argument of this subroutine',ch10,&
&  '  is probably either singular or nearly singular.',ch10,&
&  '  The LAPACK routine cgetrf failed.',ch10,&
&  '  Action : Contact ABINIT group '
  call wrtout(6,message,'COLL') ; call leave_new('COLL')
 end if
#if defined HAVE_GW_DPC
 call zgetri(n,a,n,ipvt,work,n,ierr)
#else
 call cgetri(n,a,n,ipvt,work,n,ierr)
#endif
 if(ierr/=0) then
  write(message,'(10a)')ch10,&
&  ' matcginv : BUG -',ch10,&
&  '  The matrix that has been passed in argument of this subroutine',ch10,&
&  '  is probably either singular or nearly singular.',ch10,&
&  '  The LAPACK routine cgetri failed.',ch10,&
&  '  Action : Contact ABINIT group '
  call wrtout(6,message,'COLL') ; call leave_new('COLL')
 end if
#endif

 deallocate(work)

end subroutine matcginv
!!***
