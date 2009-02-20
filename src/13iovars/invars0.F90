!{\src2tex{textfont=tt}}
!!****f* ABINIT/invars0
!! NAME
!! invars0
!!
!! FUNCTION
!! Initialisation phase : prepare the main input subroutine call by
!! reading most of the NO MULTI variables, as well as natom, and ntypat,
!! needed for allocating some input arrays in abinit, and also useri
!! and userr. The variable usewvl is also read here for later reading
!! of input path for the atomic orbital file (if required).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  iout=unit number of output file
!!  lenstr=actual length of string
!!  ndtset= number of datasets to be read; if 0, no multi-dataset mode
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!               one data set.
!!  string*(*)=string of characters containing all input variables and data
!!
!! OUTPUT
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are initialized here :
!!   cpus,jdtset,natom,npsp,ntypat,useri*,userr*
!!  istatr=repetition rate for status file
!!  istatshft=shift of the repetition rate for status file
!!  mxnatom=maximal value of input natom for all the datasets
!!  mxntypat=maximal value of input ntypat for all the datasets
!!  npsp=number of pseudopotentials
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      intagm,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine invars0(dtsets,iout,istatr,istatshft,lenstr,&
& mxnatom,mxntypat,ndtset,ndtset_alloc,npsp,string)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12parser
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,lenstr,ndtset,ndtset_alloc
 integer,intent(out) :: istatr,istatshft,mxnatom,mxntypat,npsp
 character(len=*),intent(in) :: string
!arrays
 type(dataset_type),intent(out) :: dtsets(0:ndtset_alloc)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,idtset,ii,jdtset,marr,natom,tjdtset,tread,treadh,treadm
 integer :: treads
 real(dp) :: cpus
 character(len=30) :: token
 character(len=500) :: message
!arrays
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)

!******************************************************************

 marr=max(ndtset_alloc,2)
 allocate(dprarr(marr),intarr(marr))

!Set up jdtset
 if(ndtset/=0)then

! Default values
  dtsets(0)%jdtset = -1 ! unused value
  dtsets(1:ndtset_alloc)%jdtset=(/ (ii,ii=1,ndtset_alloc) /)

! Read explicitely the jdtset array
  token = 'jdtset'
  call intagm(dprarr,intarr,0,marr,ndtset,string(1:lenstr),&
&  token,tjdtset,'INT')
  if(tjdtset==1) dtsets(1:ndtset)%jdtset=intarr(1:ndtset)

! Read the udtset array
  token = 'udtset'
  call intagm(dprarr,intarr,0,marr,2,string(1:lenstr),&
&  token,tread,'INT')

! jdtset and udtset cannot be defined together
  if(tjdtset==1 .and. tread==1)then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' invars0: ERROR -',ch10,&
&   '  jdtset and udtset cannot be defined both in the input file.',ch10,&
&   '  Action : remove one of them from your input file.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! Check values of udtset
  if(tread==1)then
   if(intarr(1)<1 .or. intarr(1)>9)then
    write(message, '(a,a,a,a,i5,a,a,a)' )ch10,&
&    ' invars0: ERROR -',ch10,&
&    '  udtset(1) must be between 1 and 9, but it is',intarr(1),'.',ch10,&
&    '  Action : change the value of udtset(1) in your input file.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if(intarr(2)<1 .or. intarr(2)>9)then
    write(message, '(a,a,a,a,i5,a,a,a)' )ch10,&
&    ' invars0: ERROR -',ch10,&
&    '  udtset(2) must be between 1 and 9, but it is',intarr(2),'.',ch10,&
&    '  Action : change the value of udtset(2) in your input file.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if(intarr(1)*intarr(2) /= ndtset)then
    write(message, '(a,a,a,a,a,a,i2,a,a,a,i2,a,i3,a,a,a,i3,a,a,a)' )ch10,&
&    ' invars0: ERROR -',ch10,&
&    '  udtset(1)*udtset(2) must be equal to ndtset,',ch10,&
&    '  but it is observed that udtset(1)=',intarr(1),',',ch10,&
&    '  and udtset(2)=',intarr(2),' so that their product is',&
&    intarr(1)*intarr(2),',',ch10,&
&    '  while ndtset is',ndtset,'.',ch10,&
&    '  Action : change udtset or ndtset in your input file.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   idtset=0
   do i1=1,intarr(1)
    do i2=1,intarr(2)
     idtset=idtset+1
     dtsets(idtset)%jdtset=i1*10+i2
    end do
   end do
  end if

! Final check on the jdtset values
  do idtset=1,ndtset
   if(dtsets(idtset)%jdtset<1 .or. dtsets(idtset)%jdtset>99)then
    write(message, '(a,a,a,a,a,a,i3,a,i3,a,a)' )ch10, &
&    ' invars0: ERROR -',ch10,&
&    '  The components of jdtset must be between 1 and 99.',ch10,&
&    '  However, the input value of the component',idtset,&
&    ' of jdtset is',dtsets(idtset)%jdtset,ch10,&
&    '  Action : correct jdtset in your input file.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
  end do

 else
  dtsets(1)%jdtset=0
 end if

 istatr=49
#if defined T3E
!The T3E has slow IO, so the rate for status files must be small by default
 istatr=149
#endif
 token = 'istatr'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) istatr=intarr(1)

 istatshft=1
 token = 'istatshft'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) istatshft=intarr(1)

 cpus=0.0_dp
 token = 'cpus '
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,treads,'DPR')
 if(treads==1) cpus=dprarr(1)
 token = 'cpum '
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,treadm,'DPR')
 if(treadm==1) cpus=dprarr(1)*60.0_dp
 token = 'cpuh '
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,treadh,'DPR')
 if(treadh==1) cpus=dprarr(1)*3600.0_dp
 if(treads+treadm+treadh>1)then
  write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&  ' invars0: ERROR -',ch10,&
&  '  More than one input variable is used to defined the CPU time limit.',&
&  ch10,'  This is not allowed.',ch10,&
&  '  Action : in the input file, suppress either cpus, cpum or cpuh.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 dtsets(:)%cpus=cpus

!Default for natom, ntypat, useri and userr
 dtsets(:)%natom=1
 dtsets(:)%ntypat=1 ; dtsets(0)%ntypat=0    ! Will always echo ntypat
 dtsets(:)%useria=0
 dtsets(:)%userib=0
 dtsets(:)%useric=0
 dtsets(:)%userid=0
 dtsets(:)%userie=0
 dtsets(:)%userra=zero
 dtsets(:)%userrb=zero
 dtsets(:)%userrc=zero
 dtsets(:)%userrd=zero
 dtsets(:)%userre=zero
 dtsets(:)%usewvl = 0

!Loop on datasets, to find natom and mxnatom, as well as useri and userr
 do idtset=1,ndtset_alloc
  jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0

! Read natom from string
  token = 'natom'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
! Might initialize natom from CML file
  if(tread==0)then
   token = '_natom'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
  end if

  if(tread==1)then
   dtsets(idtset)%natom=intarr(1)
  else
   write(message, '(a,a,a,a,i6,a,a)' ) ch10,&
&   ' invars0: ERROR -',ch10,&
&   '  Input natom must be defined, but was absent for dataset',jdtset,ch10,&
&   '  Action : check the input file.'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if
! Check that natom is greater than 0
  if (dtsets(idtset)%natom<=0) then
   write(message, '(a,a,a,a,i12,a,a,i6,a,a,a)' ) ch10,&
&   ' invars0: ERROR -',ch10,&
&   '  Input natom must be > 0, but was ',dtsets(idtset)%natom,ch10,&
&   '  for dataset',jdtset,'. This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(06,  message,'COLL')
   call leave_new('COLL')
  end if

! Read ntypat from string
  token = 'ntypat'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
  if(tread==1)dtsets(idtset)%ntypat=intarr(1)
! Check that ntypat is greater than 0
  if (dtsets(idtset)%ntypat<=0) then
   write(message, '(a,a,a,a,i12,a,a,i6,a,a,a)' ) ch10,&
&   ' invars0: ERROR -',ch10,&
&   '  Input ntypat must be > 0, but was ',dtsets(idtset)%ntypat,ch10,&
&   '  for dataset',jdtset,'. This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(06,  message,'COLL')
   call leave_new('COLL')
  end if

  token = 'ntype'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
  if(tread==1) then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' invars0: ERROR -',ch10,&
&   '  The use of the "ntype" input variable is forbidden since version 4.1 .',ch10,&
&   '  Action : replace "ntype" by "ntypat".'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  token = 'useria'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
  if(tread==1) dtsets(idtset)%useria=intarr(1)
  token = 'userib'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
  if(tread==1) dtsets(idtset)%userib=intarr(1)
  token = 'useric'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
  if(tread==1) dtsets(idtset)%useric=intarr(1)
  token = 'userid'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
  if(tread==1) dtsets(idtset)%userid=intarr(1)
  token = 'userie'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
  if(tread==1) dtsets(idtset)%userie=intarr(1)

  token = 'userra'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
  if(tread==1) dtsets(idtset)%userra=dprarr(1)
  token = 'userrb'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
  if(tread==1) dtsets(idtset)%userrb=dprarr(1)
  token = 'userrc'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
  if(tread==1) dtsets(idtset)%userrc=dprarr(1)
  token = 'userrd'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
  if(tread==1) dtsets(idtset)%userrd=dprarr(1)
  token = 'userre'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
  if(tread==1) dtsets(idtset)%userre=dprarr(1)

  token = 'usewvl'
  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
  if(tread==1) dtsets(idtset)%usewvl=intarr(1)

 end do

!mxnatom =maxval(dtsets(1:ndtset_alloc)%natom)
!mxntypat =maxval(dtsets(1:ndtset_alloc)%ntypat)
!There is a bug in the HP compiler, the following should execute properly
 mxnatom=dtsets(1)%natom ; mxntypat=dtsets(1)%ntypat
 if(ndtset_alloc>1)then
  do ii=2,ndtset_alloc
   mxnatom =max(dtsets(ii)%natom,mxnatom)
   mxntypat =max(dtsets(ii)%ntypat,mxntypat)
  end do
 end if

!Set up npsp
 npsp=mxntypat   ! Default value
 token = 'npsp'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1)then
  npsp=intarr(1)
 else
  if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
    if(dtsets(idtset)%ntypat/=mxntypat)then
     write(message, '(8a,i6,a,i6,a,a,i6,2a)' ) ch10,&
&     ' invars0: ERROR -',ch10,&
&     '  When npsp is not defined, the input variable ntypat must be',ch10,&
&     '  the same for all datasets. However, it has been found that for',ch10,&
&     '  jdtset:',dtsets(idtset)%jdtset,', ntypat=',dtsets(idtset)%ntypat,ch10,&
&     '  differs from the maximum value of ntypat=',mxntypat,ch10,&
&     '  Action : check the input variables npsp and ntypat.'
     call wrtout(06,  message,'COLL')
     call leave_new('COLL')
    end if
   end do
  end if
 end if
 dtsets(0)%npsp = mxntypat   ! Default value
 dtsets(1:ndtset_alloc)%npsp = npsp

!BEGIN TF_CHANGES
 dtsets(0)%paral_rf = 0    ! Default value
 dtsets(0)%ngroup_rf = 1   ! Default value
 do idtset=1,ndtset_alloc
  token = 'paral_rf'
  call intagm(dprarr,intarr,dtsets(idtset)%jdtset,marr,1,string(1:lenstr),token,tread,'INT')
  if(tread==1) then
   dtsets(idtset)%paral_rf=intarr(1)
  else
   dtsets(idtset)%paral_rf=0
  end if

  token = 'ngroup_rf'
  call intagm(dprarr,intarr,dtsets(idtset)%jdtset,marr,1,string(1:lenstr),token,tread,'INT')
  if(tread==1) then
   dtsets(idtset)%ngroup_rf=intarr(1)
  else
   dtsets(idtset)%ngroup_rf=1
  end if

  if(dtsets(idtset)%ngroup_rf<1 .and. dtsets(idtset)%paral_rf==1) then
   write(message, * ) ch10,&
&   ' invars0: ERROR -',ch10,&
&   '  When paral_rf is 1, ngroup_rf must be greater then 0',ch10,&
&   '  However, it has been found that for',ch10,&
&   '  paral_rf=',dtsets(idtset)%paral_rf,ch10,&
&   '  ngroup_rf=',dtsets(idtset)%ngroup_rf,ch10,&
&   '  Action : check the input variables paral_rf and ngroup_rf.'
   call wrtout(06,  message,'COLL')
   call leave_new('COLL')
  end if
 end do

!END TF_CHANGES

 deallocate(dprarr,intarr)

!We allocate the internal array, depending on the computed values.
!WARNING : do not forget to deallocate these arrays in the routine dtsetfree
!WARNING : do not forget to allocate and deallocate these arrays in the routine newsp
!(should make a separate subroutine for allocating/deallocating these records)
 do idtset=0,ndtset_alloc
  allocate(dtsets(idtset)%algalch(mxntypat))
  allocate(dtsets(idtset)%amu(mxntypat))
  allocate(dtsets(idtset)%corecs(mxntypat))
  allocate(dtsets(idtset)%densty(mxntypat,4))
  allocate(dtsets(idtset)%kberry(3,20)) !some times (ex: v1/t60) nberry is not set
  allocate(dtsets(idtset)%jpawu(mxntypat))
  allocate(dtsets(idtset)%lexexch(mxntypat))
  allocate(dtsets(idtset)%lpawu(mxntypat))
  allocate(dtsets(idtset)%mixalch(npsp,mxntypat))
  allocate(dtsets(idtset)%ptcharge(mxntypat))
  allocate(dtsets(idtset)%quadmom(mxntypat))
  allocate(dtsets(idtset)%so_psp(npsp))
  allocate(dtsets(idtset)%shiftk(3,8)) !some times nshiftk is not set
  allocate(dtsets(idtset)%upawu(mxntypat))
  allocate(dtsets(idtset)%ziontypat(mxntypat))
  allocate(dtsets(idtset)%znucl(npsp))
  allocate(dtsets(idtset)%iatfix(3,mxnatom))
  allocate(dtsets(idtset)%ratsph(mxntypat))
  allocate(dtsets(idtset)%spinat(3,mxnatom))
  allocate(dtsets(idtset)%typat(mxnatom))
  allocate(dtsets(idtset)%vel_orig(3,mxnatom))
  allocate(dtsets(idtset)%xred_orig(3,mxnatom))
 end do

end subroutine invars0
!!***
