!{\src2tex{textfont=tt}}
!!****f* ABINIT/blok8
!!
!! NAME
!! blok8
!!
!! FUNCTION
!! This routine reads and writes blocks of data in the DDBs.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! choice=(1 => read), (2 => write), (3 => write minimal info )
!! mpert =maximum number of ipert
!! msize=maximum size of the arrays flags and values
!! natom=number of atoms
!! nunit=unit number for the data block file
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! blktyp=type of the block:
!!   0 => total energy
!!   1 => second-order energy derivatives, non-stationary block
!!   2 => second-order energy derivatives, stationary block
!!   3 => third-order energy derivatives
!!   4 => first-order energy derivatives: forces, stresses and polarization
!!   5 => second-order eigenvalue derivatives
!! blkflg(msize)=flag for every matrix element (0=> the element is
!!  not in the data block), (1=> the element is in the data blok)
!! blkqpt(9)=wavevector of the perturbation(s). The elements from
!!  1 to 3 are used if we are dealing with the 2nd derivative of
!!  total energy (only one wavevector), while all elements are
!!  used in case of a third order derivative of total energy
!!  (three wavevector could be present)
!! blknrm(3)=normalization factors for the three allowed wavevectors.
!! blkval(2,msize)=real(dp), complex, value of the
!!  matrix elements that are present in the data block
!!
!! NOTES
!! only executed by one processor.
!!
!! PARENTS
!!      d3output,gstate,mrgddb,rdddb9
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine blok8(blkflg,blknrm,blkqpt,&
&     blktyp,blkval,choice,mband,mpert,msize,natom,nkpt,nunit,blkval2,kpt)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: choice,mband,mpert,msize,natom,nkpt,nunit
 integer,intent(inout) :: blktyp
!arrays
 integer,intent(inout) :: blkflg(msize)
 real(dp),intent(inout) :: blknrm(3),blkqpt(9),blkval(2,msize)
 real(dp),intent(inout),optional :: blkval2(2,msize,mband,nkpt),kpt(3,nkpt)

!Local variables -------------------------
!scalars
 integer :: band,iband,idir1,idir2,idir3,ii,ikpt,index,ipert1,ipert2,ipert3
 integer :: nelmts
 real(dp) :: ai,ar
 character(len=32) :: name
 character(len=500) :: message

! *********************************************************************
 
 if(choice==1)then
! Zero every flag
  blkflg(1:msize)=0
  if(present(blkval2))blkval2(:,:,:,:)=zero
  if(present(kpt))kpt(:,:)=zero

! Read the block type and number of elements
  read(nunit,*)
  read(nunit, '(a32,12x,i8)' )name,nelmts
  if(name==' 2rd derivatives (non-stat.)  - ')then
   blktyp=1
  else if(name==' 2rd derivatives (stationary) - ')then
   blktyp=2
  else if(name==' 3rd derivatives              - ')then
   blktyp=3
  else if(name==' Total energy                 - ')then
   blktyp=0
  else if(name==' 1st derivatives              - ')then
   blktyp=4
  else if(name==' 2rd eigenvalue derivatives   - ')then
   blktyp=5
  else
   write(message, '(a,a,a,a,a,a,a,a)' )&
&   ' blok8 : ERROR -',ch10,&
&   '  The following string appears in the DDB in place of',&
&   ' the block type description :',ch10,name,ch10,&
&   '  Action : check your DDB.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! Read the 2nd derivative block
  if(blktyp==1.or.blktyp==2)then

!  First check if there is enough space to read it
   if(msize<(3*mpert*3*mpert))then
    write(message, '(a,a,a,a,a,i10,a,i10,a,a,a)' )&
&    ' blok8 : ERROR -',ch10,&
&    '  There is not enough space to read a second-derivative block.',&
&    ch10,'  The size provided is only ',msize,' although ',&
&    3*mpert*3*mpert,' is needed.',ch10,&
&    '  Action : increase msize and recompile.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

!  Read the phonon wavevector
   read(nunit, '(4x,3es16.8,f6.1)' )(blkqpt(ii),ii=1,3),blknrm(1)

!  write(6,*)' Blok8 : For graphic purposes, format changed'

!  Read every element
   do ii=1,nelmts
    read(nunit,*)idir1,ipert1,idir2,ipert2,ar,ai
    index=idir1+                                     &
&    3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))
    blkflg(index)=1
    blkval(1,index)=ar
    blkval(2,index)=ai
   end do

!  Read the 3nd derivative block
  else if(blktyp==3)then

!  First check if there is enough space to read it
   if(msize<(3*mpert*3*mpert*3*mpert))then
    write(message, '(a,a,a,a,a,i10,a,i10,a,a,a)' )&
&    ' blok8 : ERROR -',ch10,&
&    '  There is not enough space to read a third-derivative block.',&
&    ch10,'  The size provided is only ',msize,' although ',&
&    3*mpert*3*mpert*3*mpert,' is needed.',ch10,&
&    '  Action : increase msize and recompile.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

!  Read the perturbation wavevectors
   read(nunit,'(4x,3es16.8,f6.1)')(blkqpt(ii),ii=1,3),blknrm(1)
   read(nunit,'(4x,3es16.8,f6.1)')(blkqpt(ii),ii=4,6),blknrm(2)
   read(nunit,'(4x,3es16.8,f6.1)')(blkqpt(ii),ii=7,9),blknrm(3)

!  Read every element
   do ii=1,nelmts
    read(nunit,'(6i4,2d22.14)')&
&    idir1,ipert1,idir2,ipert2,idir3,ipert3,ar,ai
    index=idir1+                                              &
&    3*((ipert1-1)+mpert*((idir2-1)+                 &
&    3*((ipert2-1)+mpert*((idir3-1)+3*(ipert3-1)))))
    blkflg(index)=1
    blkval(1,index)=ar
    blkval(2,index)=ai
   end do

!  Read the total energy
  else if(blktyp==0)then

!  First check if there is enough space to read it
   if(msize<1)then
    write(message, '(a,a,a,a,a,i10,a,a,a,a)' )&
&    ' blok8 : ERROR -',ch10,&
&    '  There is not enough space to read a total energy block.',&
&    ch10,'  The size provided is only ',msize,' although 1',&
&    ' is needed.',ch10,&
&    '  Action : increase msize and recompile.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

!  Read the total energy
   read(nunit,'(2d22.14)')ar,ai
   blkflg(1)=1
   blkval(1,1)=ar
   blkval(2,1)=ai


!  Read the 1st derivative block
  else if (blktyp == 4) then

!  First check if there is enough space to read it
   if (msize < (3*mpert)) then
    write(message, '(a,a,a,a,a,i10,a,i10,a,a,a)' )&
&    ' blok8 : ERROR -',ch10,&
&    '  There is not enough space to read a first-derivative block.',&
&    ch10,'  The size provided is only ',msize,' although ',&
&    3*mpert,' is needed.',ch10,&
&    '  Action : increase msize and recompile.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

!  Read every element
   do ii=1,nelmts
    read(nunit,'(2i4,2d22.14)')&
&    idir1,ipert1,ar,ai
    index=idir1 + 3*(ipert1 - 1)
    blkflg(index)=1
    blkval(1,index)=ar
    blkval(2,index)=ai
   end do

!  Read the 2nd eigenvalue derivative block
  else if(blktyp==5)then

!  First check if there is enough space to read it
   if(msize<(3*mpert*3*mpert))then
    write(message, '(a,a,a,a,a,i10,a,i10,a,a,a)' )&
&    ' blok8 : ERROR -',ch10,&
&    '  There is not enough space to read a second-derivative block.',&
&    ch10,'  The size provided is only ',msize,' although ',&
&    3*mpert*3*mpert*mband*nkpt,' is needed.',ch10,&
&    '  Action : increase msize and recompile.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

!  Read the phonon wavevector
!  read(nunit, '(5x,54a)' ) message
!  write(*,*)'blok9: ',message(49:54)
   read(nunit, '(4x,3es16.8,f6.1)' )(blkqpt(ii),ii=1,3),blknrm(1)

!  Read the K point and band
   if(present(blkval2).and.present(kpt))then
    do ikpt=1,nkpt
     read(nunit, '(9x,3es16.8)')(kpt(ii,ikpt),ii=1,3)
     do iband=1,mband
      read(nunit, '(6x,i3)') band
!     Read every element
      do ii=1,nelmts
       read(nunit,*)idir1,ipert1,idir2,ipert2,ar,ai
       index=idir1+                                     &
&       3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))
       blkflg(index)=1
       blkval2(1,index,iband,ikpt)=ar
       blkval2(2,index,iband,ikpt)=ai
      end do !nelmts
     end do  !band
    end do   !kpt
   end if
  end if

! The other option is to write the block
 else if(choice==2 .or. choice==3)then

! Count the number of elements
  nelmts=0
  do ii=1,msize
   if(blkflg(ii)==1)nelmts=nelmts+1
  end do

! Write the block type and number of elements
  write(nunit,*)' '
  if (blktyp == 0) then
   write(nunit, '(a,i8)' )&
&   ' Total energy                 - # elements :',nelmts
  else if (blktyp==1) then
   write(nunit, '(a,i8)' )&
&   ' 2rd derivatives (non-stat.)  - # elements :',nelmts
  else if(blktyp==2) then
   write(nunit, '(a,i8)' )&
&   ' 2rd derivatives (stationary) - # elements :',nelmts
  else if(blktyp==3) then
   write(nunit, '(a,i8)' )&
&   ' 3rd derivatives              - # elements :',nelmts
  else if (blktyp == 4) then
   write(nunit, '(a,i8)' )&
&   ' 1st derivatives              - # elements :',nelmts
  else if (blktyp == 5) then
   write(nunit, '(a,i8)' )&
&   ' 2rd eigenvalue derivatives   - # elements :',nelmts
  end if

! Write the 2nd derivative block
  if(blktyp==1.or.blktyp==2)then

!  Write the phonon wavevector
   write(nunit, '(a,3es16.8,f6.1)' )' qpt',(blkqpt(ii),ii=1,3),blknrm(1)

!  Write the matrix elements
   if(choice==2)then
    ii=0
    do ipert2=1,mpert
     do idir2=1,3
      do ipert1=1,mpert
       do idir1=1,3
        ii=ii+1
        if(blkflg(ii)==1)then
         write(nunit,'(4i4,2d22.14)')idir1,ipert1,idir2,ipert2,&
&         blkval(1,ii),blkval(2,ii)
        end if
       end do
      end do
     end do
    end do
   end if

!  Write the 3nd derivative block
  else if(blktyp==3)then

!  Write the phonon wavevectors
   write(nunit, '(a,3es16.8,f6.1)' )&
&   ' qpt',(blkqpt(ii),ii=1,3),blknrm(1)
   write(nunit, '(a,3es16.8,f6.1)' )&
&   '    ',(blkqpt(ii),ii=4,6),blknrm(2)
   write(nunit, '(a,3es16.8,f6.1)' )&
&   '    ',(blkqpt(ii),ii=7,9),blknrm(3)

!  Write the matrix elements
   if(choice==2)then
    ii=0
    do ipert3=1,mpert
     do idir3=1,3
      do ipert2=1,mpert
       do idir2=1,3
        do ipert1=1,mpert
         do idir1=1,3
          ii=ii+1
          if(blkflg(ii)==1)then
           write(nunit, '(6i4,2d22.14)' )&
&           idir1,ipert1,idir2,ipert2,idir3,ipert3,&
&           blkval(1,ii),blkval(2,ii)
          end if
         end do
        end do
       end do
      end do
     end do
    end do
   end if

!  Write total energy
  else if (blktyp == 0) then
   if (choice == 2) then
    write(nunit,'(2d22.14)')blkval(1,1),blkval(2,1)
   end if

!  Write the 1st derivative blok
  else if (blktyp == 4) then
   if (choice == 2) then
    ii = 0
    do ipert1 = 1, mpert
     do idir1 = 1, 3
      ii = ii + 1
      if (blkflg(ii) == 1) then
       write(nunit,'(2i4,2d22.14)')idir1,ipert1,&
&       blkval(1,ii),blkval(2,ii)
      end if
     end do
    end do
   end if

  else if (blktyp==5) then
!  Write the phonon wavevector
   write(nunit, '(a,3es16.8,f6.1)' )' qpt',(blkqpt(ii),ii=1,3),blknrm(1)
!  Write the matrix elements
   if(choice==2)then
    if(present(blkval2).and.present(kpt))then
     do ikpt=1,nkpt
      write(nunit,'(a,3es16.8)')' K-point:',(kpt(ii,ikpt),ii=1,3)
      do iband=1,mband
       write(nunit,'(a,i3)')' Band:',iband
       ii=0
       do ipert2=1,mpert
        do idir2=1,3
         do ipert1=1,mpert
          do idir1=1,3
           ii=ii+1
           if(blkflg(ii)==1)then
            write(nunit,'(4i4,2d22.14)')idir1,ipert1,idir2,ipert2,&
&            blkval2(1,ii,iband,ikpt),blkval2(2,ii,iband,ikpt)
           end if
          end do !idir1
         end do  !ipert1
        end do   !idir2
       end do    !ipert2
      end do     !iband
     end do      !ikpt
    end if !blkval2
   end if !choice
  end if !blktyp
 end if !choice again
end subroutine blok8
!!***
