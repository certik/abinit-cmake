!{\src2tex{textfont=tt}}
!!****f* ABINIT/inkpts
!!
!! NAME
!! inkpts
!!
!! FUNCTION
!! Initialize k points (list of k points, weights, storage)
!! for one particular dataset, characterized by jdtset.
!! Note that nkpt can be computed by called this routine with
!! input value of nkpt=0, provided kptopt/=0.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! bravais(11): bravais(1)=iholohedry
!!              bravais(2)=center
!!              bravais(3:11)=coordinates of rprim in the axes
!!               of the conventional bravais lattice (*2 if center/=0)
!! iout=unit number for echoed output
!! iscf= ( <= 0 =>non-SCF), >0 => SCF)
!! jdtset=number of the dataset looked for
!! lenstr=actual length of the string
!! kptopt=option for the generation of k points
!! msym=default maximal number of symmetries
!! nqpt=number of q points (0 or 1)
!! nsym=number of symetries
!! occopt=option for occupation numbers
!! response=0 if GS case, =1 if RF case.
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! string*(*)=character string containing all the input data.
!!  Initialized previously in instrng.
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)=symmetry operations in real space in terms
!!  of primitive translations
!! tnons(3,nsym)=nonsymmorphic translations for symmetry operations
!! vacuum(3)=for each direction, 0 if no vacuum, 1 if vacuum
!!
!! OUTPUT
!! dsifkpt(3)=defines a finer (densified) k-point sampling along
!!         the three primitive axis of the reciprocal space
!! kptnrm=normalisation of k points
!! kptrlatt(3,3)=k-point lattice specification (if kptopt/=0)
!! kptrlen=length of the smallest real space supercell vector
!! nshiftk=actual number of k-point shifts in shiftk (if kptopt/=0)
!! qpt(3)=reduced coordinates of eventual q point.
!! qptnrm=eventual normalisation of q point
!! shiftk(3,8)=shift vectors for k point generation (if kptopt/=0)
!! If nkpt/=0  the following arrays are also output :
!!  istwfk(nkpt)=option parameters that describes the storage of wfs
!!  kpt(3,nkpt)=reduced coordinates of k points.
!!  wtk(nkpt)=weight assigned to each k point.
!!
!! SIDE EFFECTS
!! Input/output:
!! nkpt=number of k points
!!  if non-zero at input, is only an input variable
!!  if zero at input, its actual value will be computed
!!
!! NOTES
!! Warning : this routine can be called with nkpt=0 (in which
!! case it returns the true value of nkpt), which can lead
!! to strange bugs in the debugging proceudre, if
!! one tries to print wtk or istwfk, in this case !
!!
!! PARENTS
!!      invars1,invars2
!!
!! CHILDREN
!!      getkgrid,intagm,leave_new,testkgrid,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine inkpts(bravais,dsifkpt,iout,iscf,istwfk,jdtset,&
& kpt,kptopt,kptnrm,kptrlatt,&
& kptrlen,lenstr,msym,nkpt,nqpt,nshiftk,nsym,&
& occopt,qpt,qptnrm,response,rprimd,&
& shiftk,string,symafm,symrel,tnons,userid,vacuum,wtk)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_12parser
 use interfaces_13recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,iscf,jdtset,kptopt,lenstr,msym,nqpt,nsym,occopt
 integer,intent(in) :: response,userid
 integer,intent(inout) :: nkpt
 integer,intent(out) :: nshiftk
 real(dp),intent(out) :: kptnrm,kptrlen,qptnrm
 character(len=*),intent(in) :: string
!arrays
 integer,intent(in) :: bravais(11),symafm(msym),symrel(3,3,msym),vacuum(3)
 integer,intent(out) :: dsifkpt(3),istwfk(nkpt),kptrlatt(3,3)
 real(dp),intent(in) :: rprimd(3,3),tnons(3,msym)
 real(dp),intent(out) :: kpt(3,nkpt),qpt(3),shiftk(3,8),wtk(nkpt)

!Local variables-------------------------------
!scalars
 integer :: bit0,dkpt,ii,ikpt,jj,jkpt,kk,marr,mkpt,ndiv_small,nkpt_computed
 integer :: nsegment,prtkpt,tread,tread_kptrlatt,tread_ngkpt
 real(dp) :: fraction,norm,ucvol,wtksum
 character(len=30) :: token
 character(len=500) :: message
!arrays
 integer :: bit(3),ngkpt(3)
 integer,allocatable :: ndivk(:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),rmet(3,3),tsec(2)
 real(dp),allocatable :: kptbounds(:,:)
!no_abirules
!Dummy arguments for subroutine 'intagm' to parse input file
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)

! *************************************************************************

!DEBUG
!write(6,*)' inkpts : enter'
!Warning : do not try to print istwfk when the input nkpt=0
!write(6,*)' inkpts : istwfk=',istwfk
!stop
!ENDDEBUG

 call timab(192,1,tsec)

!Compute the maximum size of arrays intarr and dprarr
 marr=max(3*nkpt,3*8)
 allocate(intarr(marr),dprarr(marr))

!Read qpt and qptnrm
 token = 'qpt'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'DPR')
 if(tread==1) qpt(1:3)=dprarr(1:3)
 qptnrm=1.0_dp
 token = 'qptnrm'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) qptnrm=dprarr(1)
 if(qptnrm<tol10)then
  write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&  ' inkpts: ERROR -',ch10,&
&  '  The input variable qptnrm is lower than 1.0d-10,',ch10,&
&  '  while it must be a positive, non-zero number.   ',ch10,&
&  '  Action : correct the qptnrm in the input file.'
  call wrtout(iout,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!Initialize kptrlen
 kptrlen=20.0_dp
 token = 'kptrlen'
 call intagm(dprarr,intarr,jdtset,marr,&
& 1,string(1:lenstr),token,tread,'DPR')
 if(tread==1)kptrlen=dprarr(1)

!Initialize kpt, kptnrm and wtk according to kptopt.
!For kptopt==0, one must have nkpt defined.
 if(kptopt==0)then

  kpt(:,:)=0.0_dp
  token = 'kpt'
  call intagm(dprarr,intarr,jdtset,marr,&
&  3*nkpt,string(1:lenstr),token,tread,'DPR')
  if(tread==1) kpt(:,:)=reshape( dprarr(1:3*nkpt), (/3,nkpt/) )
  kptnrm=1.0_dp
  token = 'kptnrm'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
  if(tread==1) kptnrm=dprarr(1)

! Only read wtk when iscf >0 or iscf=-1 or iscf=-3
! (this last option is for Zach Levine)
! Normalize the k-point weights when occopt/=2
! Check that k point weights add to 1 when occopt==2
  if  (iscf>0.or.iscf==-1.or.iscf==-3)  then

   token = 'wtk'
   call intagm(dprarr,intarr,jdtset,marr,nkpt,string(1:lenstr),token,tread,'DPR')
   if(tread==1) wtk(1:nkpt)=dprarr(1:nkpt)

   wtksum=sum(wtk(:))
   write(message, '(a,i4,a,f12.6)' ) &
&   ' inkpts: Sum of ',nkpt,' k point weights is',wtksum
   call wrtout(06,message,'COLL')
   if (wtksum<1.d-6) then
    write(message, '(6a)' ) ch10,&
&    ' inkpts: ERROR -',ch10,&
&    '  This sum is too close to zero. ',ch10,&
&    '  Action : correct the array wtk in the input file.'
    call wrtout(iout,message,'COLL')
    call wrtout(06,  message,'COLL')
    call leave_new('COLL')
   end if
   if (abs(wtksum-1.0_dp)>1.d-06) then
    if(occopt==2)then
     write(message, '(a,a,a,a,1p,e18.8,a,a,a)' ) ch10,&
&     ' inkpts: ERROR -',ch10,&
&     '  wtksum=',wtksum,' /=1.0 means wts do not add to 1 , while occopt=2.',&
&     ch10,' Action : correct the array wtk in input file.'
     call wrtout(iout,message,'COLL')
     call wrtout(06,  message,'COLL')
     call leave_new('COLL')
    else
     write(message,  '(a,i4,a)' ) &
&     '   With present occopt=',occopt,' , renormalize it to one'
     call wrtout(06,message,'COLL')
     norm=1.0_dp/wtksum
     wtk(1:nkpt)=wtk(1:nkpt)*norm
!    call DSCAL(nkpt,norm,wtk,1)
    end if
   end if
  end if

 else if (kptopt<0) then

! Band structure calculation
  nsegment=abs(kptopt)

  if(iscf/=-2)then
   write(message,  '(a,a,a,a,a,i4,a,a,a)' ) &
&   ' inkpts : ERROR -',ch10,&
&   '  For a negative kptopt, iscf must be -2,',ch10,&
&   '  while it is found to be',iscf,'.',ch10,&
&   '  Action : change the value of iscf in your input file, or change kptopt.'
   call wrtout(06,message,'COLL') 
   call leave_new('COLL')
  end if

  if(marr<3*nsegment+3)then
   marr=3*nsegment+3
   deallocate(dprarr,intarr)
   allocate(dprarr(marr),intarr(marr))
  end if

  allocate(kptbounds(3,nsegment+1),ndivk(nsegment))

  token = 'kptbounds'
  call intagm(dprarr,intarr,jdtset,marr,&
&  3*nsegment+3,string(1:lenstr),token,tread,'DPR')
  if(tread==1)then
   kptbounds(:,:)=reshape( dprarr(1:3*nsegment+3), (/3,nsegment+1/) )
  else
   write(message,  '(a,a,a,a,a,a,a)' ) &
&   ' inkpts : ERROR -',ch10,&
&   '  When kptopt is negative, kptbounds must be initialized ',ch10,&
&   '  in the input file, which is not the case.',ch10,&
&   '  Action : initialize kptbounds in your input file, or change kptopt.'
   call wrtout(06,message,'COLL') 
   call leave_new('COLL')
  end if

  token = 'ndivk'
  call intagm(dprarr,intarr,jdtset,marr,nsegment,string(1:lenstr),&
&  token,tread,'INT')
  if(tread==1)then
   ndivk(1:nsegment)=intarr(1:nsegment)
!  The 1 stand for the first point
   nkpt_computed=1+sum(ndivk(1:nsegment))
  else
!  Calculate ndivk such as the path is normalized
!  Note that if both ndivk and ndivsm are defined in in input file, only ndivk is used !
   token = 'ndivsm'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1)then
    ndiv_small=intarr(1)
    call metric(gmet,gprimd,std_out,rmet,rprimd,ucvol)
    call mknormpath(nsegment+1,kptbounds,gmet,ndiv_small,ndivk,nkpt_computed)
   else
    write(message,'(7a)') &
&    ' inkpts : ERROR -',ch10,&
&    '  When kptopt is negative, ndivsm or ndivk must be initialized ',ch10,&
&    '  in the input file, which is not the case.',ch10,&
&    '  Action : initialize ndivsm or ndivk in your input file, or change kptopt.'
    call wrtout(06,message,'COLL') 
    call leave_new('COLL')
   end if
  end if

! Check that the argument nkpt is coherent with nkpt_computed,
  if(nkpt/=0 .and. nkpt/=nkpt_computed)then
   write(message,  '(a,a,a,i6,a,a,a,a,a,i6,a,a,a,a,a,a,a)' ) &
&   ' inkpts : BUG -',ch10,&
&   '  The argument nkpt=',nkpt,', does not match',ch10,&
&   '  the number of k points generated by kptopt, ndivk, kptbounds,',ch10,&
&   '  and the eventual symmetries, that is, nkpt=',nkpt_computed,'.',&
&   ch10,&
&   '  However, note that it might due to the user,',ch10,&
&   '  if nkpt is explicitely defined in the input file.',ch10,&
&   '  In this case, please check your input file.'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  if (nkpt/=0) then
!  the array kpt has the right dimension and we can generate the k-path
   token = 'kptbounds'
   call intagm(dprarr,intarr,jdtset,marr,&
&   3*nsegment+3,string(1:lenstr),token,tread,'DPR')
   if(tread==1)then
    kptbounds(:,:)=reshape( dprarr(1:3*nsegment+3), (/3,nsegment+1/) )
   else
    write(message,  '(a,a,a,a,a,a,a)' ) &
&    ' inkpts : ERROR -',ch10,&
&    '  When kptopt is negative, kptbounds must be initialized ',ch10,&
&    '  in the input file, which is not the case.',ch10,&
&    '  Action : initialize kptbounds in your input file, or change kptopt.'
    call wrtout(06,message,'COLL')
    call leave_new('COLL')
   end if

!  First k point
   jkpt=1
   kpt(:,1)=kptbounds(:,1)
   do ii=1,nsegment
    dkpt=ndivk(ii)
    do ikpt=1,dkpt
     fraction=dble(ikpt)/dble(dkpt)
     kpt(:,ikpt+jkpt)=fraction *kptbounds(:,ii+1)+&
&     (1.0_dp-fraction)*kptbounds(:,ii)
    end do
    jkpt=jkpt+dkpt
   end do

  end if

  kptnrm=1.0_dp
  token = 'kptnrm'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
  if(tread==1) kptnrm=dprarr(1)

  deallocate(kptbounds,ndivk)

 else if (kptopt>=1 .and. kptopt<=4) then

  ngkpt(:)=0
  token = 'ngkpt'
  call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),&
&  token,tread_ngkpt,'INT')
  if(tread_ngkpt==1)then
   ngkpt(1:3)=intarr(1:3)
   do ii=1,3
    if(ngkpt(ii)<1)then
     write(message,  '(a,a,a,i1,a,a,a,i4,a,a,a)' ) &
&     ' inkpts : ERROR -',ch10,&
&     '  The input variable ngkpt(',ii,') must be strictly positive,',ch10,&
&     '  while it is found to be',ngkpt(ii),'.',ch10,&
&     '  Action : change it in your input file, or change kptopt.'
     call wrtout(06,message,'COLL')
     call leave_new('COLL')
    end if
   end do
  end if

  token = 'kptrlatt'
  call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),&
&  token,tread_kptrlatt,'INT')
  if(tread_kptrlatt==1)&
&  kptrlatt(:,:)=reshape(intarr(1:9), (/3,3/) )

  if(tread_ngkpt==1 .and. tread_kptrlatt==1)then
   write(message,  '(a,a,a,a,a,a,a)' ) &
&   ' inkpts : ERROR -',ch10,&
&   '  The input variables ngkpt and kptrlatt cannot both ',ch10,&
&   '  be defined in the input file.',ch10,&
&   '  Action : change one of ngkpt or kptrlatt in your input file.'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  else if(tread_ngkpt==1)then
   kptrlatt(:,:)=0
   kptrlatt(1,1)=ngkpt(1)
   kptrlatt(2,2)=ngkpt(2)
   kptrlatt(3,3)=ngkpt(3)
  end if

! DEBUG
! write(6,*)' inkpts : after init, kptrlatt='
! write(6, '(9i3)' )kptrlatt(:,:)
! ENDDEBUG

  token = 'nshiftk'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),&
&  token,tread,'INT')
  if(tread==1)nshiftk=intarr(1)

  if(nshiftk<1 .or. nshiftk>8 )then
   write(message,  '(6a,i4,a,a,a)' )ch10, &
&   ' inkpts : ERROR -',ch10,&
&   '  The only allowed values of nshiftk are between 1 and 8,',ch10,&
&   '  while it is found to be',nshiftk,'.',ch10,&
&   '  Action : change the value of nshiftk in your input file, or change kptopt.'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  token = 'shiftk'
  call intagm(dprarr,intarr,jdtset,marr,&
&  3*nshiftk,string(1:lenstr),token,tread,'DPR')
  if(tread==1)then
   shiftk(:,1:nshiftk)=reshape( dprarr(1:3*nshiftk), (/3,nshiftk/) )
  else
   if(nshiftk/=1)then
    write(message,  '(6a,i4,a,a)' )ch10, &
&    ' inkpts : ERROR -',ch10,&
&    '  When nshiftk is not equal to 1, shiftk must be defined in the input file.',ch10,&
&    '  However, shiftk is not defined, while nshiftk=',nshiftk,ch10,&
&    '  Action : change the value of nshiftk in your input file, or define shiftk.'
    call wrtout(06,message,'COLL')
    call leave_new('COLL')
   end if
  end if

  dsifkpt(:)=1
  token = 'dsifkpt'
  call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),&
&  token,tread,'INT')
  if(tread==1)dsifkpt(1:3)=intarr(1:3)

  do ii=1,3
   if(dsifkpt(ii)<1)then
    write(message,  '(4a,i1,a,a,a,i4,a,a,a)' )ch10, &
&    ' inkpts : ERROR -',ch10,&
&    '  The input variable dsifkpt(',ii,') must be strictly positive,',ch10,&
&    '  while it is found to be',dsifkpt(ii),'.',ch10,&
&    '  Action : change it in your input file.'
    call wrtout(06,message,'COLL')
    call leave_new('COLL')
   end if
  end do

  do ii = 1, 3
   if (dsifkpt(ii) > 1) then
    if (nshiftk > 1) then
     write(message,  '(6a,i1,a,i4,a,a,a,a,a)' ) ch10,&
&     ' inkpts : ERROR -',ch10,&
&     '  the variable nshiftk has been found to be > 1',ch10,&
&     '  while dsifkpt(',ii,') = ',dsifkpt(ii),'.',ch10,&
&     '  This is not allowed.',ch10,&
&     '  Action : change nshiftk or dsifkpt in your input file.'
     call wrtout(06,message,'COLL')
     call leave_new('COLL')
    end if
    do jj = 1, 3
     do kk = 1, 3
      if ((kptrlatt(jj,kk) /= 0).and.(jj /= kk)) then
       write(message,  '(6a,i1,a,i4,a,a,a,a,a)' )ch10, &
&       ' inkpts : ERROR -',ch10,&
&       '  the variable kptrlatt has been found to be not diagonal',ch10,&
&       '  while dsifkpt(',ii,') = ',dsifkpt(ii),'.',ch10,&
&       '  This is not allowed.',ch10,&
&       '  Action : change kptrlatt or dsifkpt in your input file.'
       call wrtout(06,message,'COLL')
       call leave_new('COLL')
      end if
     end do
    end do
   end if
  end do

  prtkpt=0
  token = 'prtkpt'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),&
&  token,tread,'INT')
  if(tread==1)prtkpt=intarr(1)

! DEBUG
! write(6,*)' inkpts : before testkgrid, kptrlatt='
! write(6, '(9i3)' )kptrlatt(:,:)
! ENDDEBUG

  if(sum(abs(kptrlatt(:,:)))==0)then

!  The parameters of the k lattice are not known, compute
!  kptrlatt, nshiftk, shiftk.
   call testkgrid(bravais,iout,kptrlatt,kptrlen,&
&   msym,nshiftk,nsym,prtkpt,rprimd,shiftk,symafm,symrel,tnons,vacuum)

  end if

  call getkgrid(dsifkpt,iout,iscf,kpt,kptopt,kptrlatt,kptrlen,&
&  msym,nkpt,nkpt_computed,nshiftk,nsym,rprimd,&
&  shiftk,symafm,symrel,tnons,vacuum,wtk)

! DEBUG
! write(6,*)' inkpts : after getkgrid, nkpt=',nkpt
! ENDDEBUG

  kptnrm=1.0_dp

 else

  write(message,  '(a,a,a,a,a,i4,a,a,a)' ) &
&  ' inkpts : ERROR -',ch10,&
&  '  The only values of kptopt allowed are smaller than 4.',ch10,&
&  '  The input value of kptopt is',kptopt,'.',ch10,&
&  '  Action : change kptopt in your input file.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')

 end if

 if(kptnrm<tol10)then
  write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&  ' inkpts: ERROR -',ch10,&
&  '  The input variable kptnrm is lower than 1.0d-10,',ch10,&
&  '  while it must be a positive, non-zero number.   ',ch10,&
&  '  Action : correct the kptnrm in the input file.'
  call wrtout(iout,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!The k point number has been computed, and, if nkpt/=0, also
!the list of k points.

!Now, determine istwfk, and eventually shift the k points by the
!value of qpt.

!DEBUG
!write(6,*)' inkpts : istwfk before intagm ',istwfk(1:nkpt)
!ENDDEBUG

 if(nkpt/=0)then

  istwfk(1:nkpt)=0
  token = 'istwfk'
  call intagm(dprarr,intarr,jdtset,marr,&
&  nkpt,string(1:lenstr),token,tread,'INT')
  if(tread==1) istwfk(1:nkpt)=intarr(1:nkpt)
! Impose istwfk=1 for RF calculations
  if(response==1)istwfk(1:nkpt)=1
! Impose the same for TFvW
! if(userid==1.or.userid==2)istwfk(1:nkpt)=1
! Also impose istwfk=1 for spinor calculations
  token = 'nspinor'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),&
&  token,tread,'INT')
  if(tread/=0 .and. intarr(1)/=1)istwfk(1:nkpt)=1

  do ikpt=1,nkpt
   if(istwfk(ikpt)==0)then
    bit0=1
!   DEBUG
!   write(6,*)' inkpts : ikpt,kpt(:,ikpt) ',kpt(:,ikpt)
!   write(6,*)' inkpts : kptnrm ',kptnrm
!   write(6,*)' nqpt,response=',nqpt,response
!   write(6,*)' inkpts : qpt(:)',qpt(:)
!   write(6,*)' inkpts : qptnrm',qptnrm
!   ENDDEBUG

    do ii=1,3
     kpoint(:)=kpt(:,ikpt)/kptnrm
     if(nqpt/=0 .and. response==0) kpoint(:)=kpoint(:)+qpt(:)/qptnrm
     if(abs(kpoint(ii))<1.d-10)then
      bit(ii)=0
     else if(abs(kpoint(ii)-0.5_dp)<1.e-10_dp )then
      bit(ii)=1
     else
      bit0=0
     end if
    end do
!   DEBUG
!   write(6,*)' bit0, bit(:)=',bit0,bit(1:3)
!   ENDDEBUG
    if(bit0==0)then
     istwfk(ikpt)=1
    else
!    Note the inversion between bit(2) and bit(3)
     istwfk(ikpt)=2+bit(1)+4*bit(2)+2*bit(3)
    end if
   end if
  end do
  write(message, '(a,a,6i2)' )ch10,&
&  ' inkpts : istwfk preprocessed, gives following first values (max. 6):',&
&  istwfk(1:(min(6,nkpt)))
  call wrtout(06,message,'COLL')
 end if

!If nkpt was to be computed, transfer it from nkpt_computed
 if(nkpt==0)nkpt=nkpt_computed

 deallocate(intarr,dprarr)

 call timab(192,2,tsec)

!DEBUG
!write(6,*)' inkpts : end of subroutine, kptnrm= ',kptnrm
!write(6,*)'    nkpt=',nkpt
!write(6,*)'  istwfk=',istwfk
!stop
!ENDDEBUG

end subroutine inkpts
!!***
