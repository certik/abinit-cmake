!{\src2tex{textfont=tt}}
!!****f* ABINIT/getkgrid
!!
!! NAME
!! getkgrid
!!
!! FUNCTION
!! Compute the grid of k points
!!
!! Note that nkpt can be computed by called this routine with
!! input value of nkpt=0, provided kptopt/=0.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dsifkpt(3)=defines a finer k-point sampling along
!!         the three primitive axis of the reciprocal space
!! iout=unit number for echoed output
!! iscf= ( <= 0 =>non-SCF), >0 => SCF)
!! kptopt=option for the generation of k points
!! msym=default maximal number of symmetries
!! nkpt=number of k points (might be zero, see output description)
!! nsym=number of symetries
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)=symmetry operations in real space in terms
!!  of primitive translations
!! tnons(3,nsym)=nonsymmorphic translations for symmetry operations
!! vacuum(3)=for each direction, 0 if no vacuum, 1 if vacuum
!!
!! OUTPUT
!! kptrlen=length of the smallest real space supercell vector
!!  associated with the lattice of k points.
!! nkpt_computed=number of k points, computed in the present routine
!! If nkpt/=0  the following are also output :
!! kpt(3,nkpt)=reduced coordinates of k points.
!! wtk(nkpt)=weight assigned to each k point.
!!
!!
!! SIDE EFFECTS
!! Input/Output
!! kptrlatt(3,3)=k-point lattice specification
!! nshiftk=actual number of k-point shifts in shiftk
!! shiftk(3,8)=shift vectors for k point generation
!!
!! PARENTS
!!      elphon,getshell,initberry,inkpts,nonlinear,testkgrid
!!
!! CHILDREN
!!      leave_new,mati3inv,matr3inv,metric,smallprim,smpbz,symkpt,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine getkgrid(dsifkpt,iout,iscf,kpt,kptopt,kptrlatt,kptrlen,&
& msym,nkpt,nkpt_computed,nshiftk,nsym,rprimd,shiftk,symafm,&
& symrel,tnons,vacuum,wtk)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13recipspace, except_this_one => getkgrid
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,iscf,kptopt,msym,nkpt,nsym
 integer,intent(inout) :: nshiftk
 integer,intent(out) :: nkpt_computed
 real(dp),intent(out) :: kptrlen
!arrays
 integer,intent(in) :: dsifkpt(3),symafm(msym),symrel(3,3,msym),vacuum(3)
 integer,intent(inout) :: kptrlatt(3,3)
 real(dp),intent(in) :: rprimd(3,3),tnons(3,msym)
 real(dp),intent(inout) :: shiftk(3,8)
 real(dp),intent(out) :: kpt(3,nkpt),wtk(nkpt)

!Local variables-------------------------------
!scalars
 integer :: brav,decreased,idir,ii,ikpt,ishiftk,isym,jshiftk,kshiftk,mkpt
 integer :: nkpt_fullbz,nkptlatt,nshiftk2,nshiftk4,nsym_used,ok,ok_all,option
 integer :: timrev
 real(dp) :: length2,ucvol,ucvol_super
 character(len=500) :: message
!arrays
 integer :: kptrlatt2(3,3),shiftk_pair(8)
 integer,allocatable :: indkpt(:),symrec(:,:,:)
 real(dp) :: cart(3,3),deltak(3,8),dijk(3),dmult(3),fact_vacuum(3),gmet(3,3)
 real(dp) :: gmet_super(3,3),gprimd(3,3),gprimd_super(3,3),klatt2(3,3)
 real(dp) :: klatt3(3,3),kptrlattr(3,3),ktransf(3,3),ktransf_invt(3,3)
 real(dp) :: metmin(3,3),minim(3,3),rmet(3,3),rmet_super(3,3),rprimd_super(3,3)
 real(dp) :: shiftk2(3,8),shiftk3(3,8)
 real(dp),allocatable :: kpt_fullbz(:,:),shiftk4(:,:),spkpt(:,:),wtk_folded(:)
 real(dp),allocatable :: wtk_fullbz(:)

! *************************************************************************

!DEBUG
!write(6,*)' getkgrid : enter '
!stop
!write(6,*)' kptrlen=',kptrlen
!write(6, '(a,9i4)' )' kptrlatt=',kptrlatt(:,:)
!write(6, '(a,i4)' )' nshiftk=',nshiftk
!do ii=1,nshiftk
!write(6,*)' ii, shiftk=',ii,shiftk(:,ii)
!end do
!write(6, '(3es16.6)' )rprimd(:,1)
!write(6, '(3es16.6)' )rprimd(:,2)
!write(6, '(3es16.6)' )rprimd(:,3)
!ENDDEBUG
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 if(kptopt==1)then
  nsym_used=0
! Cannot use antiferromagnetic symmetry operations
! to decrease the number of k points
  do isym=1,nsym
   if(symafm(isym)==1)nsym_used=nsym_used+1
  end do
  allocate(symrec(3,3,nsym_used))
! Get the symmetry matrices in terms of reciprocal basis
  nsym_used=0
  do isym=1,nsym
   if(symafm(isym)==1)then
    nsym_used=nsym_used+1
    call mati3inv(symrel(:,:,isym),symrec(:,:,nsym_used))
   end if
  end do
 else if(kptopt==2)then
! Use only the time-reversal
  nsym_used=1
  allocate(symrec(3,3,1))
  symrec(1:3,1:3,1)=0
  do ii=1,3
   symrec(ii,ii,1)=1
  end do
 end if

 kptrlatt2(:,:)=kptrlatt(:,:)
 nshiftk2=nshiftk
 shiftk2(:,:)=shiftk(:,:)

!Find a primitive k point lattice, if possible, by decreasing the
!number of shifts.
!WARNING : only able to take into account half-lattice vectors,
!but not third of lattice vectors (would be required to reduce
!hexagonal and rhombohedral lattices). This works only for nshiftk=2, 4, or 8.
 if(nshiftk2/=1)then

! Loop to be repeated if there has been a successful reduction of nshiftk2
  do

   decreased=0
   deltak(1,2:nshiftk2)=shiftk2(1,2:nshiftk2)-shiftk2(1,1)
   deltak(2,2:nshiftk2)=shiftk2(2,2:nshiftk2)-shiftk2(2,1)
   deltak(3,2:nshiftk2)=shiftk2(3,2:nshiftk2)-shiftk2(3,1)

!  Try different shifts
   do ishiftk=2,nshiftk2
    dmult(:)=2*deltak(:,ishiftk)
    if(sum(abs( dmult(:)-nint(dmult(:)) ))<tol8)then

     shiftk_pair(:)=0 ; shiftk_pair(1)=1 ; shiftk_pair(ishiftk)=-1

     ok_all=1
     if(nshiftk2>3)then

!     If there is more than two shifts, must be able to group
!     the other shifts by pairs
      do jshiftk=2,nshiftk2-1
       if(jshiftk==ishiftk .or. shiftk_pair(jshiftk)/=0 )cycle
       ok=0
       do kshiftk=jshiftk+1,nshiftk2
        if(kshiftk==ishiftk .or. shiftk_pair(kshiftk)/=0 )cycle
        dijk(:)=deltak(:,kshiftk)-deltak(:,jshiftk)-deltak(:,ishiftk)
        dijk(:)=dijk(:)-nint(dijk(:))
        if(sum(abs(dijk(:)))<tol8)then
         shiftk_pair(kshiftk)=-1
         shiftk_pair(jshiftk)=1
         ok=1
        end if
       end do
       if(ok==0)ok_all=0
      end do

     else if(nshiftk2==3)then
      ok_all=0
      cycle
     end if

!    DEBUG
!    write(6,*)' getkgrid : nshiftk2=',nshiftk2
!    write(6,*)' shiftk_pair(1:nshiftk2)',shiftk_pair(1:nshiftk2)
!    ENDDEBUG

     if(ok_all==1)then

      ktransf(:,:)=0.0_dp
      ktransf(1,1)=1.0_dp
      ktransf(2,2)=1.0_dp
      ktransf(3,3)=1.0_dp
!     Replace one of the unit vectors by the shift vector deltak(:,shiftk).
!     However, must pay attention not to make linear combinations.
!     Also, choose positive sign for first-non-zero value.
      if(abs(deltak(1,ishiftk)-nint(deltak(1,ishiftk)))>tol8)then
       if(deltak(1,ishiftk)>0)ktransf(:,1)= deltak(:,ishiftk)
       if(deltak(1,ishiftk)<0)ktransf(:,1)=-deltak(:,ishiftk)
      else if(abs(deltak(2,ishiftk)-nint(deltak(2,ishiftk)))>tol8)then
       if(deltak(2,ishiftk)>0)ktransf(:,2)= deltak(:,ishiftk)
       if(deltak(2,ishiftk)<0)ktransf(:,2)=-deltak(:,ishiftk)
      else if(abs(deltak(3,ishiftk)-nint(deltak(3,ishiftk)))>tol8)then
       if(deltak(3,ishiftk)>0)ktransf(:,3)= deltak(:,ishiftk)
       if(deltak(3,ishiftk)<0)ktransf(:,3)=-deltak(:,ishiftk)
      end if
!     Copy the integers to real(dp)
      kptrlattr(:,:)=kptrlatt2(:,:)
!     Go to reciprocal space
      call matr3inv(kptrlattr,klatt2)
!     Make the transformation
      do ii=1,3
       klatt3(:,ii)=ktransf(1,ii)*klatt2(:,1)+&
&       ktransf(2,ii)*klatt2(:,2)+&
&       ktransf(3,ii)*klatt2(:,3)
      end do
!     DEBUG
!     write(6,*)' klatt2 in cartesian coordinates :'
!     do ii=1,3
!     cart(:,ii)=klatt2(1,ii)*gprimd(:,1)+&
!     &                 klatt2(2,ii)*gprimd(:,2)+&
!     &                 klatt2(3,ii)*gprimd(:,3)
!     write(6,'(a,3es16.6)' ) 'cart2(:,ii)=',cart(:,ii)
!     end do
!     write(6,*)' klatt3 in cartesian coordinates :'
!     do ii=1,3
!     cart(:,ii)=klatt3(1,ii)*gprimd(:,1)+&
!     &                 klatt3(2,ii)*gprimd(:,2)+&
!     &                 klatt3(3,ii)*gprimd(:,3)
!     write(6,'(a,3es16.6)' ) 'cart3(:,ii)=',cart(:,ii)
!     end do
!     ENDDEBUG
!     Back to real space
      call matr3inv(klatt3,kptrlattr)
!     real(dp) to integer
      kptrlatt2(:,:)=nint(kptrlattr(:,:))

!     Prepare the transformation of the shifts
      call matr3inv(ktransf,ktransf_invt)
      decreased=1
      kshiftk=0
      do jshiftk=1,nshiftk2
       if(shiftk_pair(jshiftk)==1)then
        kshiftk=kshiftk+1
!       Place the shift with index jshiftk in place of the one in kshiftk,
!       also transform the shift from the old to the new coordinate system
        shiftk3(:,kshiftk)=ktransf_invt(1,:)*shiftk2(1,jshiftk)+&
&        ktransf_invt(2,:)*shiftk2(2,jshiftk)+&
&        ktransf_invt(3,:)*shiftk2(3,jshiftk)
       end if
      end do
      nshiftk2=nshiftk2/2
      shiftk2(:,:)=shiftk3(:,:)
      if(kshiftk/=nshiftk2)then
       write(message, '(a,a,a,a)' ) ch10,&
&       ' getkgrid : BUG -',ch10,&
&       '  The search for a primitive k point lattice contains a bug.'
       call wrtout(6,message,'COLL')
       call leave_new('COLL')
      end if
     end if

!    End of the condition for the trial shift to be half a lattice vector
    end if

!   If this trial shift was succesful, must exit the loop on trial ishiftk,
!   and reinitialize the global loop
    if(decreased==1)exit

!   Next trial ishiftk
   end do

   if(decreased==0 .or. nshiftk2==1)exit

!  Infinite loop
  end do

! End nshiftk being 1 or larger
 end if

!Impose shiftk coordinates to be in [0,1[
 do ishiftk=1,nshiftk2
  do ii=1,3
   if(shiftk2(ii,ishiftk)>one-tol8)&
&   shiftk2(ii,ishiftk)=shiftk2(ii,ishiftk)-1.0_dp
   if(shiftk2(ii,ishiftk)<-tol8)&
&   shiftk2(ii,ishiftk)=shiftk2(ii,ishiftk)+1.0_dp
  end do
 end do

!DEBUG
!write(6,*)' getkgrid : afer nkshift reduction '
!write(6, '(a,9i4)' )' kptrlatt2=',kptrlatt2(:,:)
!write(6, '(a,i4)' )' nshiftk2=',nshiftk2
!do ii=1,nshiftk2
!write(6,*)' ii, shiftk2=',ii,shiftk2(:,ii)
!end do
!stop
!ENDDEBUG

!Compute the number of k points in the G-space unit cell
 nkptlatt=kptrlatt2(1,1)*kptrlatt2(2,2)*kptrlatt2(3,3) &
& +kptrlatt2(1,2)*kptrlatt2(2,3)*kptrlatt2(3,1) &
& +kptrlatt2(1,3)*kptrlatt2(2,1)*kptrlatt2(3,2) &
& -kptrlatt2(1,2)*kptrlatt2(2,1)*kptrlatt2(3,3) &
& -kptrlatt2(1,3)*kptrlatt2(2,2)*kptrlatt2(3,1) &
& -kptrlatt2(1,1)*kptrlatt2(2,3)*kptrlatt2(3,2)

!Check whether the number of k points is positive,
!otherwise, change the handedness of kptrlatt2
 if(nkptlatt<=0)then
! DEBUG
! write(6,*)' getkgrid : nkptlatt is negative !'
! ENDDEBUG
  kptrlatt2(:,3)=-kptrlatt2(:,3)
  nkptlatt=-nkptlatt
  do ishiftk=1,nshiftk2
   shiftk2(3,ishiftk)=-shiftk2(3,ishiftk)
  end do
 end if

!Determine the smallest supercell R-vector whose contribution
!is not taken correctly into account in the k point integration.
!Increase enormously the size of the cell when vacuum is present.
 fact_vacuum(:)=1
 if(vacuum(1)==1)fact_vacuum(1)=1000.0_dp
 if(vacuum(2)==1)fact_vacuum(2)=1000.0_dp
 if(vacuum(3)==1)fact_vacuum(3)=1000.0_dp
 do ii=1,3
  rprimd_super(:,ii)=fact_vacuum(1)*rprimd(:,1)*kptrlatt2(1,ii)+&
&  fact_vacuum(2)*rprimd(:,2)*kptrlatt2(2,ii)+&
&  fact_vacuum(3)*rprimd(:,3)*kptrlatt2(3,ii)
 end do
!Warning : here compute quantities connect to the supercell lattice
!DEBUG
!write(6,*)' getkgrid : will enter metric with rprimd_super'
!write(6,*)' getkgrid : nkptlatt=',nkptlatt
!ENDDEBUG
 call metric(gmet_super,gprimd_super,-1,rmet_super,rprimd_super,ucvol_super)
 call smallprim(metmin,minim,rmet_super,rprimd_super)
 length2=min(metmin(1,1),metmin(2,2),metmin(3,3))
 kptrlen=sqrt(length2)
 write(message,'(a,es16.6)' )&
& ' getkgrid : length of smallest supercell vector (bohr)=',kptrlen
 call wrtout(6,message,'COLL')

!If the number of shifts has been decreased,
!determine the set of kptrlatt2 vectors with minimal length (without using fact_vacuum)
!It is worth to determine the minimal set of vectors
!so that the kptrlatt that is output does not seem screwy,
!although correct but surprising.
 if(nshiftk/=nshiftk2)then
  do ii=1,3
   rprimd_super(:,ii)=rprimd(:,1)*kptrlatt2(1,ii)+&
&   rprimd(:,2)*kptrlatt2(2,ii)+&
&   rprimd(:,3)*kptrlatt2(3,ii)
  end do
  call metric(gmet_super,gprimd_super,-1,rmet_super,rprimd_super,ucvol_super)
! Shift vectors in cartesian coordinates (reciprocal space)
  do ishiftk=1,nshiftk2
   shiftk3(:,ishiftk)=gprimd_super(:,1)*shiftk2(1,ishiftk)+&
&   gprimd_super(:,2)*shiftk2(2,ishiftk)+&
&   gprimd_super(:,3)*shiftk2(3,ishiftk)
  end do
  call smallprim(metmin,minim,rmet_super,rprimd_super)
  call metric(gmet_super,gprimd_super,-1,rmet_super,minim,ucvol_super)
! This is the new kptrlatt2
  do ii=1,3
   kptrlatt2(:,ii)=nint(gprimd(1,:)*minim(1,ii)+&
&   gprimd(2,:)*minim(2,ii)+&
&   gprimd(3,:)*minim(3,ii))
  end do
! Shifts in the new set of kptrlatt vectors
  do ishiftk=1,nshiftk2
   shiftk2(:,ishiftk)=minim(1,:)*shiftk3(1,ishiftk)+&
&   minim(2,:)*shiftk3(2,ishiftk)+&
&   minim(3,:)*shiftk3(3,ishiftk)
  end do
 end if

!DEBUG
!write(6,*)' getkgrid : before smpbz '
!stop
!ENDDEBUG

!brav=1 is able to treat all bravais lattices.
 brav=1
 mkpt=nkptlatt*nshiftk2*(sum(dsifkpt(:)) - 2)

!Define the "densified" k-point grid
!*************************************
 nshiftk4 = nshiftk2 + nshiftk2*(sum(dsifkpt(:))-3)
 allocate(shiftk4(3,nshiftk4))
 shiftk4(:,1:nshiftk2) = shiftk2(:,1:nshiftk2)
 ishiftk = nshiftk2
 do idir = 1,3
  if (dsifkpt(idir)>1) then
   do jshiftk = 1,nshiftk2
    dijk(:) = 0.0_dp
    dijk(idir) = 1.0_dp/(dsifkpt(idir)*1.0_dp)
    do ii = 1, dsifkpt(idir)-1
     ishiftk = ishiftk + 1
     shiftk4(:,ishiftk) = dijk(:)*(ii*1.0_dp) + shiftk2(:,jshiftk)
    end do
   end do
  end if
 end do

!DEBUG
!write(100,*)
!write(100,*)'nshiftk4,mkpt: ',nshiftk4,mkpt
!do ishiftk = 1,nshiftk4
!write(100,*)(shiftk4(:,ishiftk))
!end do
!write(100,*)
!stop
!ENDDEBUG

 allocate(spkpt(3,mkpt))
 option=0
 call smpbz(brav,iout,kptrlatt2,mkpt,nkpt_fullbz,nshiftk4,option,shiftk4,spkpt)

!DEBUG
!write(6,*)' getkgrid : after smpbz, write spkpt, kptopt= ',kptopt
!do ikpt=1,nkpt_fullbz
!write(6, '(i4,3es15.5)' )ikpt,spkpt(:,ikpt)
!end do
!ENDDEBUG

 if(kptopt==1 .or. kptopt==2)then

  allocate(indkpt(nkpt_fullbz),kpt_fullbz(3,nkpt_fullbz))
  allocate(wtk_fullbz(nkpt_fullbz),wtk_folded(nkpt_fullbz))

  kpt_fullbz(:,:)=spkpt(:,1:nkpt_fullbz)
  wtk_fullbz(1:nkpt_fullbz)=1.0_dp/dble(nkpt_fullbz)

  timrev=1 ; option=0
  call symkpt(gmet,indkpt,kpt_fullbz,nkpt_fullbz,&
&  nkpt_computed,nsym_used,option,symrec,timrev,wtk_fullbz,wtk_folded)
  deallocate(symrec,wtk_fullbz)

 else if(kptopt==3)then
  nkpt_computed=nkpt_fullbz
 end if

!The number of k points has been computed from kptopt, kptrlatt,
!nshiftk, shiftk,
!and the eventual symmetries, it is presently called nkpt_computed.

!Check that the argument nkpt is coherent with nkpt_computed,
!if nkpt/=0.
 if(nkpt/=nkpt_computed .and. nkpt/=0)then
  write(message,  '(a,a,a,i6,a,a,a,a,a,i6,a,a,a,a,a,a,a)' ) &
&  ' getkgrid : BUG or ERROR -',ch10,&
&  '  The argument nkpt=',nkpt,', does not match',ch10,&
&  '  the number of k points generated by kptopt, kptrlatt, shiftk,',ch10,&
&  '  and the eventual symmetries, that is, nkpt=',nkpt_computed,'.',&
&  ch10,&
&  '  However, note that it might be due to the user,',ch10,&
&  '  if nkpt is explicitely defined in the input file.',ch10,&
&  '  In this case, please check your input file.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if


 if(kptopt==1 .or. kptopt==2)then

  if(nkpt/=0)then
   do ikpt=1,nkpt
    kpt(:,ikpt)=kpt_fullbz(:,indkpt(ikpt))
!   if (dsifkpt(1)==1.and.dsifkpt(2)==1.and.dsifkpt(3)==1) then
    if(iscf>0 .or. iscf==-3 .or. iscf==-1.or.iscf==-2)wtk(ikpt)=wtk_folded(indkpt(ikpt))
!   end if
   end do
  end if
  deallocate(indkpt,kpt_fullbz,spkpt,wtk_folded)

 else if(kptopt==3)then

  if(nkpt/=0)then
   kpt(:,1:nkpt)=spkpt(:,1:nkpt)
!  if (dsifkpt(1)==1.and.dsifkpt(2)==1.and.dsifkpt(3)==1) then
   if(iscf>1 .or. iscf==-3 .or. iscf==-1.or.iscf==-2)wtk(1:nkpt)=1.0_dp/dble(nkpt)
!  end if
  end if
  deallocate(spkpt)

 end if

 kptrlatt(:,:)=kptrlatt2(:,:)
 nshiftk=nshiftk2
 shiftk(:,:)=shiftk2(:,:)

 deallocate(shiftk4)

!DEBUG
!write(6,*)' getkgrid : exit, nkpt_computed= ',nkpt_computed
!ENDDEBUG

end subroutine getkgrid
!!***
