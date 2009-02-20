!{\src2tex{textfont=tt}}
!!****f* ABINIT/setshells
!! NAME
!! setshells
!!
!! FUNCTION
!! Set consistently the number of shells, the number of plane-waves,
!! and the energy cut-off
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  gprimd(3,3)=dimensional primitive vectors in reciprocal space
!!  gmet(3,3)=metric tensor in reciprocal space
!!  nsym=number of symmetry operations
!!  symrel(3,3,nsym)=symmetry operations in real space
!!  tag=suffix to account for the different possibilities for these variables (npw, ecut or nsh ..)
!!  ucvol=unit cell volume
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  ecut,npw,nsh=one of them is an input, the two others are output
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!      kpgsph,leave_new,sort_dp,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine setshells(ecut,npw,nsh,nsym,gmet,gprimd,symrel,tag,ucvol)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13recipspace
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(inout) :: npw,nsh
 real(dp),intent(in) :: ucvol
 real(dp),intent(inout) :: ecut
 character(len=*),intent(in) :: tag
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: exchn2n3d,ifound,ig,ii,ish,isym,istat,nproc_fft,npw_found,npwave,npwwrk
 integer :: nsh_found,pad=50
 real(dp) :: ecut_found,ecut_trial,eps,scale=1.3_dp
 logical :: found
 character(len=500) :: msg
 type(MPI_type) :: mpi_enreg
!arrays
 integer :: geq(3)
 integer,allocatable :: gvec(:,:),gvec_sh(:,:),insort(:),npw_sh(:)
 real(dp) :: gctr(3)
 real(dp),allocatable :: gnorm(:),gnorm_sh(:)

!******************************************************************
!BEGIN EXECUTABLE SECTION

#ifdef DEBUG_MODE
 write(msg,'(a)')' setshells : enter '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif
 !
 ! === Check coherence of input variables ecut, npw, and nsh ===
 !
 ! 1-> one at least should be non-null
 if (npw==0.and.nsh==0.and.ecut==0) then
  write(msg,'(12a)')ch10,&
&  ' setshells: ERROR -',ch10,&
&  '  One of the three variables ecut',TRIM(tag),', npw',TRIM(tag),', or nsh',TRIM(tag),&
&  '  must be non-null.',ch10,&
&  '  Action : modify the value of one of these in input file.'
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if
 !
 ! 2-> one and only one should be non-null
 if (npw/=0.and.nsh/=0) then
  write(msg,'(10a)')ch10,&
&  ' setshells: ERROR -',ch10,&
&  '  Only one of the two variables npw',TRIM(tag),' and nsh',TRIM(tag),&
&  '  can be non-null.',ch10,&
&  '  Action : modify the value of one of these in input file.'
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if
 if (ecut>tol6.and.npw/=0) then
  write(msg,'(10a)')ch10,&
&  ' setshells: ERROR -',ch10,&
&  '  Only one of the two variables ecut',TRIM(tag),' and npw',TRIM(tag),&
&  '  can be non-null.',ch10,&
&  '  Action : modify the value of one of these in input file.'
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if
 if (ecut>tol6.and.nsh/=0) then
  write(msg,'(10a)')ch10,&
&  ' setshells: ERROR -',ch10,&
&  '  Only one of the two variables ecut',TRIM(tag),' and nsh',TRIM(tag),&
&  '  can be non-null.',ch10,&
&  '  Action : modify the value of one of these in input file.'
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if
 !
 ! === Calculates an upper bound for npw ===
 ! gctr is center of the g-vector sphere
 gctr(:)=(/zero,zero,zero/)
 if (ecut>tol6) then
  ! The average number of plane-waves in the cutoff sphere is given by: 
  ! npwave = (2*ecut)**(3/2)*ucvol/(6*pi**2)
  ! The upper bound is calculated as npwwrk=int(scale * npwave) + pad
  npwave=NINT(ucvol*(two*ecut)**1.5_dp/(six*pi**2))
  npwwrk=NINT(DBLE(npwave)*scale)+pad
  ecut_trial=ecut
 else if (npw/=0) then
  ! npw is given in the input
  npwwrk=NINT(DBLE(npw)*scale)+pad
  ecut_trial=(six*pi**2*npw/ucvol)**two_thirds/two
 else
  ! If nsh is given in the input
  npwwrk=nsh*18+2*pad
  ecut_trial=(six*pi**2*nsh*18/ucvol)**two_thirds/two
 end if

 allocate(gvec(3,npwwrk),stat=istat)
 ifound=0
 do while(ifound==0)
  write(msg,'(a,f8.2)')' setshells : ecut_trial = ',ecut_trial
  call wrtout(std_out,msg,'COLL')
  exchn2n3d=0 ! For the time being, no exchange of n2 and n3
  mpi_enreg%me_fft=0   ! Sequential
  mpi_enreg%nproc_fft=1 ! Sequential
  mpi_enreg%paral_compil_fft=0 ! Sequential
  call kpgsph(ecut_trial,exchn2n3d,gmet,0,1,1,gvec,gctr,1,mpi_enreg,npwwrk,npw_found)
  allocate(gnorm(npw_found),insort(npw_found))

  do ig=1,npw_found
   insort(ig)=ig
   gnorm(ig)=zero
   do ii=1,3
    gnorm(ig)=gnorm(ig)+(gvec(1,ig)*gprimd(ii,1)+&
&                        gvec(2,ig)*gprimd(ii,2)+&
&                        gvec(3,ig)*gprimd(ii,3))**2
   end do
  end do
  call sort_dp(npw_found,gnorm,insort,tol14)
  allocate(npw_sh(npw_found),gnorm_sh(npw_found),gvec_sh(3,npw_found))
  npw_sh(:)=0
  gnorm_sh(:)=zero
  gvec_sh(:,:)=0
  ! Count the number of shells:
  ! (search for the G-vectors generating the others by symmetry)
  nsh_found=0
  do ig=1,npw_found
   eps=1.d-8*gnorm(ig)
   found=.FALSE.
   ish=1
   do while ((.not.found).and.(ish<=nsh_found))
    if (ABS(gnorm(ig)-gnorm_sh(ish))<=eps) then
     isym=1
     do while ((.not.found).and.(isym<=nsym))
      geq(:)=(symrel(1,:,isym)*gvec(1,insort(ig))+&
&             symrel(2,:,isym)*gvec(2,insort(ig))+&
&             symrel(3,:,isym)*gvec(3,insort(ig)))

      found=((geq(1)==gvec_sh(1,ish)).and.&
&            (geq(2)==gvec_sh(2,ish)).and.&
&            (geq(3)==gvec_sh(3,ish)))
      isym=isym+1
     end do
    end if
    ish=ish+1
   end do
   if (.not.found) then
    nsh_found=nsh_found+1
    gnorm_sh(nsh_found)=gnorm(ig)
    gvec_sh(:,nsh_found)=gvec(:,insort(ig))
    npw_sh(nsh_found)=1
   else
    ish=ish-1
    npw_sh(ish)=npw_sh(ish)+1
   end if
  end do
  ecut_found=two*pi**2*gnorm(npw_found)
! 
  if(ecut>tol6) then
   ! ecut is given in the input
   if (ecut_found<ecut-tol6) then
    write(msg,'(6a,e14.6,9a,e14.6,3a)')ch10,&
&    ' setshells: WARNING -',ch10,&
&    '  The value ecut',TRIM(tag),'=',ecut,' given in the input file leads to',ch10,&
&    '  the same values for nsh',TRIM(tag),' and npw',TRIM(tag),' as ecut',TRIM(tag),'=',ecut_found,ch10,&
&    '  This value will be adopted for the calculation.',ch10
    call wrtout(std_out,msg,'COLL')
   end if
   ifound=1
  elseif(npw/=0) then
   ! If npw is given in the input
   if (npw_found==npw) then
    ecut_found=two*pi**2*gnorm(npw_found)
    ifound=1
   else if(npw_found>npw) then
    npw_found=0
    nsh_found=0
    do while(npw_found<npw)
     nsh_found=nsh_found+1
     npw_found=npw_found+npw_sh(nsh_found)
    end do
    ! check that the shell is closed
    if(npw_found>npw) then
     ! shell not closed
     npw_found=npw_found-npw_sh(nsh_found)
     nsh_found=nsh_found-1
     do while (ABS(gnorm_sh(nsh_found)-gnorm_sh(nsh_found+1))<0.000001)
      npw_found=npw_found-npw_sh(nsh_found)
      nsh_found=nsh_found-1
     end do
     write(msg,'(6a,i6,6a,i6,3a)')ch10,&
&     ' setshells: WARNING -',ch10,&
&     '  The value npw',TRIM(tag),'=',npw,' given in the input file does',&
&     '  not close the shell',ch10,&
&     '  The lower closed-shell is obtained for a value npw',TRIM(tag),'=',npw_found,ch10,&
&     '  This value will be adopted for the calculation.',ch10
     call wrtout(std_out,msg,'COLL')
    end if
    ecut_found=two*pi**2*gnorm(npw_found)
    ifound=1
   end if
  else if (nsh/=0) then
   ! If nsh is given in the input
   if (nsh_found==nsh) then
    ecut_found=two*pi**2*gnorm(npw_found)
    ifound=1
   else if (nsh_found>nsh) then
    npw_found=0
    nsh_found=0
    do ish=1,nsh
     npw_found=npw_found+npw_sh(ish)
     nsh_found=nsh_found+1
    end do
    if (ABS(gnorm_sh(nsh_found)-gnorm_sh(nsh_found+1))<0.000001) then
     do while (ABS(gnorm_sh(nsh_found)-gnorm_sh(nsh_found+1))<0.000001)
      nsh_found=nsh_found+1
      npw_found=npw_found+npw_sh(nsh_found)
     end do
     write(msg,'(6a,i6,6a,i6,3a)')ch10,&
&     ' setshells: WARNING -',ch10,&
&     '  The value nsh',TRIM(tag),'=',nsh,' given in the input file',&
&     ' corresponds to the same',ch10,&
&     '  cut-off energy as for closed-shell upto nsh',TRIM(tag),'=',nsh_found,ch10,&
&     '  This value will be adopted for the calculation.',ch10
     call wrtout(std_out,msg,'COLL')
    end if
    ecut_found=two*pi**2*gnorm(npw_found)
    ifound=1
   end if
  end if

  if (ifound==0) then
   ecut_trial=1.1*ecut_trial
   deallocate(gnorm,gnorm_sh,gvec_sh,insort,npw_sh)
  else
   ecut=ecut_found
   npw=npw_found
   nsh=nsh_found
  end if
 end do !while(ifound==0)

 deallocate(gnorm,gnorm_sh,gvec,gvec_sh,insort,npw_sh)

#ifdef DEBUG_MODE
 write(msg,'(a)')' setshells : exit '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif
 
end subroutine setshells
!!***
