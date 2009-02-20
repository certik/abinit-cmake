!{\src2tex{textfont=tt}}
!!****f* ABINIT/invars9
!!
!! NAME
!! invars9
!!
!! FUNCTION
!! Open input file for the anaddb code, then
!! reads or echoes the input information.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG,JCC,CL,MVeithen,XW)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! lenstr=actual length of string
!! natom=number of atoms, needed for atifc
!! nunit=unit number for error messages
!! qtol=tolerance for wavevector comparison
!! string*(*)=string of characters containing all input variables and data
!!
!! OUTPUT
!! anaddb_dtset= (derived datatype) contains all the input variables
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      intagm,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine invars9 (anaddb_dtset,lenstr,natom,nunit,qtol,string)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12parser
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: lenstr,natom,nunit
 real(dp),intent(in) :: qtol
 character(len=*),intent(in) :: string
 type(anaddb_dataset_type),intent(out) :: anaddb_dtset

!Local variables -------------------------
!Dummy arguments for subroutine 'intagm' to parse input file
!Set routine version number here:
!scalars
 integer,parameter :: vrsddb=010929
 integer :: iatom,ii,iph1,iph2,iqshft,istrain,jdtset,marr,tread,vrsinddb
 character(len=30) :: token
 character(len=500) :: message
!arrays
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)

!*********************************************************************

!DEBUG
!write(6,*)' invars9 : enter '
!ENDDEBUG

 marr=3
 allocate(intarr(marr),dprarr(marr))

 jdtset=1

 vrsinddb=010929
 token = 'vrsinddb'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) vrsinddb=intarr(1)
 if(vrsinddb/=vrsddb)then
  write(message, '(a,a,a,i6,a,a,a,i6,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  The DDB version number from the input file, vrsinddb= ',vrsinddb,',',ch10,&
&  '  does not agree with the code DDB version number ',vrsddb,'.',ch10,&
&  '  Action : make the input file, the code, and the DDB consistent.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%rfmeth=1
 token = 'rfmeth'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%rfmeth=intarr(1)
 if(anaddb_dtset%rfmeth<1.or.anaddb_dtset%rfmeth>2)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  rfmeth is',anaddb_dtset%rfmeth,', but the only allowed values',ch10,&
&  '  are 1 or 2 . ',ch10,&
&  '  Action : correct rfmeth in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%enunit=0
 token = 'enunit'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%enunit=intarr(1)
 if(anaddb_dtset%enunit<0.or.anaddb_dtset%enunit>2)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  enunit is',anaddb_dtset%enunit,', but the only allowed values',ch10,&
&  '  are 0, 1 or 2 . ',ch10,&
&  '  Action : correct enunit in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%eivec=0
 token = 'eivec'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%eivec=intarr(1)
 if(anaddb_dtset%eivec<0.or.anaddb_dtset%eivec>4)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  eivec is',anaddb_dtset%eivec,', but the only allowed values',ch10,&
&  '  are 0, 1, 2, 3 or 4.',ch10,&
&  '  Action : correct eivec in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%symdynmat=1
 token = 'symdynmat'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%symdynmat=intarr(1)
 if(anaddb_dtset%symdynmat/=0.and.anaddb_dtset%symdynmat/=1)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  symdynmat is',anaddb_dtset%symdynmat,', but the only allowed values',ch10,&
&  '  are 0, or 1.',ch10,&
&  '  Action : correct symdynmat in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%alphon=0
 token = 'alphon'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%alphon=intarr(1)

 anaddb_dtset%asr=1
 token = 'asr'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%asr=intarr(1)
 if(anaddb_dtset%asr<-2.or.anaddb_dtset%asr>4)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  asr is',anaddb_dtset%asr,', but the only allowed values',ch10,&
&  '  are 0, 1, 2, 3, 4, -1 or -2 .',ch10,&
&  '  Action : correct asr in your input file.'
! Note : negative values are allowed when the acoustic sum rule
! is to be applied after the analysis of IFCs
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%chneut=0
 token = 'chneut'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%chneut=intarr(1)
 if(anaddb_dtset%chneut<0.or.anaddb_dtset%chneut>2)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  chneut is',anaddb_dtset%chneut,', but the only allowed values',ch10,&
&  '  are 0, 1 or 2 .',ch10,&
&  '  Action : correct chneut in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%ramansr=0
 token = 'ramansr'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%ramansr=intarr(1)
 if(anaddb_dtset%ramansr<0.or.anaddb_dtset%ramansr>2)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  ramansr is',anaddb_dtset%ramansr,', but the only allowed values',ch10,&
&  '  are 0, 1 or 2 .',ch10,&
&  '  Action : correct ramansr in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%selectz=0
 token = 'selectz'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%selectz=intarr(1)
 if(anaddb_dtset%selectz<0.or.anaddb_dtset%selectz>2)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  selectz is',anaddb_dtset%selectz,', but the only allowed values',ch10,&
&  '  are 0, 1 or 2 .',ch10,&
&  '  Action : correct selectz in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%dieflag=0
 token = 'dieflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%dieflag=intarr(1)
 if(anaddb_dtset%dieflag<0.or.anaddb_dtset%dieflag>4)then
  write(message, '(3a,i8,5a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  dieflag is',anaddb_dtset%dieflag,', but the only allowed values',ch10,&
&  '  are 0, 1, 2, 3 or 4.',ch10,&
&  '  Action : correct dieflag in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%elaflag=0
 token = 'elaflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%elaflag=intarr(1)
 if(anaddb_dtset%elaflag<0.or.anaddb_dtset%elaflag>5)then
  write(message,'(3a,i8,5a)' )&
&  ' invars9 : Error -',ch10,&
&  '  elaflag is',anaddb_dtset%elaflag,', but the only allowed values',ch10,&
&  '  are 0,1,2,3,4 or 5 .',ch10,&
&  '  Action : correct elaflag in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%elphflag=0
 token = 'elphflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%elphflag=intarr(1)
 if(anaddb_dtset%elphflag<0.or.anaddb_dtset%elphflag>1)then
  write(message,'(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : Error -',ch10,&
&  '  elphflag is',anaddb_dtset%elphflag,', but the only allowed values',ch10,&
&  '  are 0, or 1.',ch10,&
&  '  Action : correct elphflag in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!typical value for mustar, but can vary sensibly with the metal
 anaddb_dtset%mustar = 0.1_dp
 token = 'mustar'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) then
  anaddb_dtset%mustar=dprarr(1)
 else if (anaddb_dtset%elphflag==1) then
  write (*,*) ' invars9 : warning. Default value of mustar = ', anaddb_dtset%mustar,&
&  'will be used'
 end if

!typical value for gaussian smearing, but can vary sensibly with the metal
 anaddb_dtset%elphsmear = 0.01_dp
 token = 'elphsmear'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENERGY')
 if(tread==1) then
  anaddb_dtset%elphsmear=dprarr(1)
 else if (anaddb_dtset%elphflag==1) then
  write (*,*) ' invars9 : warning. Default value of elphsmear = ', &
&  anaddb_dtset%elphsmear, 'will be used'
 end if

!typical value for gaussian smearing of a2F function
 anaddb_dtset%a2fsmear = 0.00002_dp
 token = 'a2fsmear'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENERGY')
 if(tread==1) then
  anaddb_dtset%a2fsmear=dprarr(1)
 else if (anaddb_dtset%elphflag==1) then
  write (*,*) ' invars9 : warning. Default value of a2fsmear = ', &
&  anaddb_dtset%a2fsmear, 'will be used'
 end if

 anaddb_dtset%nqpath=0
 token = 'nqpath'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
  anaddb_dtset%nqpath=intarr(1)
 else if (tread==0 .and. anaddb_dtset%elphflag==1) then
  write(message,'(a,a,a,a,a)' )&
&  ' invars9 : Error -',ch10,&
&  '  elphflag is 1 but no nqpath has been specified for phonon linewidths',ch10,&
&  '  Action : specify nqpath and qpath(3,nqpath) in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if(anaddb_dtset%nqpath<0)then
  write(message,'(a,a,a,i8,a,a,a)' )&
&  ' invars9 : Error -',ch10,&
&  '  nqpath is',anaddb_dtset%nqpath,', but must be positive',ch10,&
&  '  Action : correct elphflag in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 if (anaddb_dtset%nqpath > 0) then
  allocate (anaddb_dtset%qpath(3,anaddb_dtset%nqpath))
  if(3*anaddb_dtset%nqpath>marr)then
   marr=3*anaddb_dtset%nqpath
   deallocate(intarr,dprarr)
   allocate(intarr(marr),dprarr(marr))
  end if
  anaddb_dtset%qpath(:,:)=zero
  token = 'qpath'
  call intagm(dprarr,intarr,jdtset,marr,3*anaddb_dtset%nqpath,string(1:lenstr),token,tread,'DPR')
  if(tread==1) then
   anaddb_dtset%qpath(1:3,1:anaddb_dtset%nqpath)=&
&   reshape(dprarr(1:3*anaddb_dtset%nqpath),(/3,anaddb_dtset%nqpath/))
  else
   write(message,'(a,a,a,a,a)')&
&   ' invars9 : Error -',ch10,&
&   '  nqpath is non zero but qpath is absent ',ch10,&
&   '  Action : specify qpath in your input file.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
 end if

!Default is use gaussian integration
 anaddb_dtset%telphint = 1
 token = 'telphint'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%telphint = intarr(1)
 if(anaddb_dtset%telphint < 0 .or. anaddb_dtset%telphint > 2) then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  telphint is',anaddb_dtset%telphint,', but the only allowed values',ch10,&
&  '  are 0 (tetrahedron) or 1 (gaussian) or 2 (set of bands occupied ep_b_min,ep_b_max).',ch10,&
&  '  Action : correct telphint in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Default is 0 - not used unless telphint==2
 anaddb_dtset%ep_b_min = 0
 anaddb_dtset%ep_b_max = 0
 token = 'ep_b_min'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%ep_b_min = intarr(1)
 token = 'ep_b_max'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%ep_b_max = intarr(1)
 if(anaddb_dtset%telphint /= 2 .and. (anaddb_dtset%ep_b_min /= 0 .or. anaddb_dtset%ep_b_max /= 0)) then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : WARNING -',ch10,&
&  '  telphint is',anaddb_dtset%telphint,', but ep_b_min or ep_b_max',ch10,&
&  '  are set non zero. They will not be used'
  call wrtout(6,message,'COLL')
 else if(anaddb_dtset%telphint == 2 .and. (anaddb_dtset%ep_b_min == 0 .or. anaddb_dtset%ep_b_max == 0)) then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : WARNING -',ch10,&
&  '  telphint is',anaddb_dtset%telphint,', but ep_b_min or ep_b_max',ch10,&
&  '  are not both set. ',ch10,&
&  '  Action : set ep_b_min and ep_b_max in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Default is no output for PHDOS
 anaddb_dtset%prtdos=0
 token='prtdos'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%prtdos = intarr(1)
 if(anaddb_dtset%prtdos < 0 .or. anaddb_dtset%prtdos >= 2) then
  write(message, '(3a,i8,5a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  prtdos is',anaddb_dtset%prtdos,', but the only allowed values',ch10,&
&  '  are 0 (no output) or 1 (gaussian method)',ch10,  &
&  '  Action : correct prtdos in your input file.'
  call wrtout(6,message,'COLL') 
  call leave_new('COLL')
 end if
 if(anaddb_dtset%prtdos/=0 .and. anaddb_dtset%ifcflag/=1) then
  write(message, '(5a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  ifcflag must be 1 when the calculation of the phonon DOS is required ',ch10,&
&  '  Action : correct ifcflag in your input file.'
  call wrtout(6,message,'COLL') 
  call leave_new('COLL')
 end if
 
!Default is no output of the nesting factor
 anaddb_dtset%prtnest = 0
 token = 'prtnest'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%prtnest = intarr(1)
 if(anaddb_dtset%prtnest < 0 .or. anaddb_dtset%prtnest > 2) then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  prtnest is',anaddb_dtset%prtnest,', but the only allowed values',ch10,&
&  '  are 0 (no nesting), 1 (XY format) or 2 (XY + Xcrysden format)',ch10,  &
&  '  Action : correct prtnest in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Default is no output for the Fermi Surface
 anaddb_dtset%prtfsurf = 0
 token = 'prtfsurf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%prtfsurf = intarr(1)
 if(anaddb_dtset%prtfsurf < 0 .or. anaddb_dtset%prtfsurf > 2) then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  prtfsurf is',anaddb_dtset%prtfsurf,', but the only allowed values',ch10,&
&  '  are 0 (no output) or 1 (Xcrysden bxsf format)',ch10,  &
&  '  Action : correct prtfsurf in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%kptrlatt(:,:)=0
 if (anaddb_dtset%telphint == 0 .or. anaddb_dtset%prtnest == 1 .or. &
& anaddb_dtset%prtnest == 2 .or. anaddb_dtset%prtfsurf== 1) then
  marr = 9
  deallocate(intarr,dprarr)
  allocate(intarr(marr),dprarr(marr))
  token = 'kptrlatt'
  call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),token,tread,'INT')
  if(tread==1) then
   anaddb_dtset%kptrlatt(1:3,1:3)=reshape(intarr(1:9),(/3,3/))
  else
   write (*,*) ' invars9 :  ERROR : if tetrahedron integration is used, ',              &
&   'or the output of the nesting function/Fermi surface is required, ',          &
&   'you must specify the kptrlatt'
   stop
  end if
 end if

 anaddb_dtset%gkk2exist = 0
 token = 'gkk2exist'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%gkk2exist = intarr(1)
 if(anaddb_dtset%gkk2exist < 0 .or. anaddb_dtset%gkk2exist > 1) then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  gkk2exist is',anaddb_dtset%gkk2exist,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct gkk2exist in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%gkk2write = 0
 token = 'gkk2write'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%gkk2write = intarr(1)
 if(anaddb_dtset%gkk2write < 0 .or. anaddb_dtset%gkk2write > 1) then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  gkk2write is',anaddb_dtset%gkk2write,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct gkk2write in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%gkk_rptexist = 0
 token = 'gkk_rptexist'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%gkk_rptexist = intarr(1)
 if(anaddb_dtset%gkk_rptexist < 0 .or. anaddb_dtset%gkk_rptexist > 1) then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  gkk_rptexist is',anaddb_dtset%gkk_rptexist,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct gkk_rptexist in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%gkk_rptwrite = 0
 token = 'gkk_rptwrite'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%gkk_rptwrite = intarr(1)
 if(anaddb_dtset%gkk_rptwrite < 0 .or. anaddb_dtset%gkk_rptwrite > 1) then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  gkk_rptwrite is',anaddb_dtset%gkk_rptwrite,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct gkk_rptwrite in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%gkqexist = 0
 token = 'gkqexist'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%gkqexist = intarr(1)
 if(anaddb_dtset%gkqexist < 0 .or. anaddb_dtset%gkqexist > 1) then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  gkqexist is',anaddb_dtset%gkqexist,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct gkqexist in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%gkqwrite = 0
 token = 'gkqwrite'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%gkqwrite = intarr(1)
 if(anaddb_dtset%gkqwrite < 0 .or. anaddb_dtset%gkqwrite > 1) then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  gkqwrite is',anaddb_dtset%gkqwrite,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct gkqwrite in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 
 anaddb_dtset%ifltransport = 0
 token = 'ifltransport'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%ifltransport = intarr(1)
 if(anaddb_dtset%ifltransport < 0 .or. anaddb_dtset%ifltransport > 1) then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  ifltransport is',anaddb_dtset%ifltransport,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct ifltransport in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 
 anaddb_dtset%phfrqexist = 0
 token = 'phfrqexist'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%phfrqexist = intarr(1)
 if(anaddb_dtset%phfrqexist < 0 .or. anaddb_dtset%phfrqexist > 1) then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  phfrqexist is',anaddb_dtset%phfrqexist,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct phfrqexist in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%phfrqwrite = 0
 token = 'phfrqwrite'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%phfrqwrite = intarr(1)
 if(anaddb_dtset%phfrqwrite < 0 .or. anaddb_dtset%phfrqwrite > 1) then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  phfrqwrite is',anaddb_dtset%phfrqwrite,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct phfrqwrite in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if


 anaddb_dtset%tkeepbands = 0
 token = 'tkeepbands'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%tkeepbands = intarr(1)
 if(anaddb_dtset%tkeepbands < 0 .or. anaddb_dtset%tkeepbands > 1) then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  tkeepbands is',anaddb_dtset%tkeepbands,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct tkeepbands in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%doscalprod = 0
 token = 'doscalprod'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%doscalprod = intarr(1)
 if(anaddb_dtset%doscalprod < 0 .or. anaddb_dtset%doscalprod > 1) then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  doscalprod is',anaddb_dtset%doscalprod,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct doscalprod in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!By default use the real fermie, but allow
 anaddb_dtset%elph_fermie = zero
 token = 'elph_fermie'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENERGY')
 if(tread==1) then
  anaddb_dtset%elph_fermie=dprarr(1)
 end if


 anaddb_dtset%ifcflag=0
 token = 'ifcflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%ifcflag=intarr(1)
 if(anaddb_dtset%ifcflag<0.or.anaddb_dtset%ifcflag>1)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  ifcflag is',anaddb_dtset%ifcflag,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct ifcflag in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%instrflag=0
 token = 'instrflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%instrflag=intarr(1)
 if(anaddb_dtset%instrflag<0.or.anaddb_dtset%instrflag>1)then
  write(message,'(3a,i8,5a)' )&
&  'invars9 : Error -',ch10,&
&  ' instrflag is',anaddb_dtset%instrflag,', but the only allowed values',ch10,&
&  '  are 0, 1  .',ch10,&
&  '  Action : correct instrflag in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%piezoflag=0
 token = 'piezoflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%piezoflag=intarr(1)
 if(anaddb_dtset%piezoflag<0.or.anaddb_dtset%piezoflag>7)then
  write(message,'(3a,i8,5a)' )&
&  'invars9 : Error -',ch10,&
&  ' piezoflag is',anaddb_dtset%piezoflag,', but the only allowed values',ch10,&
&  '  are 0, 1,2,3,4,5,6,7  .',ch10,&
&  '  Action : correct piezoflag in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%nlflag=0
 token = 'nlflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%nlflag=intarr(1)
 if(anaddb_dtset%nlflag<0.or.anaddb_dtset%nlflag>2)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  nlflag is',anaddb_dtset%nlflag,', but the only allowed values',ch10,&
&  '  are 0, 1 or 2.',ch10,&
&  '  Action : correct nlflag in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%thmflag=0
 token = 'thmflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%thmflag=intarr(1)
 if(anaddb_dtset%thmflag<0.or.anaddb_dtset%thmflag>3)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  thmflag is',anaddb_dtset%thmflag,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct thmflag in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%polflag=0
 token = 'polflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%polflag=intarr(1)
 if(anaddb_dtset%polflag<0.or.anaddb_dtset%polflag>1)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  polflag is',anaddb_dtset%polflag,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct polflag in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%relaxat=0
 token = 'relaxat'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%relaxat=intarr(1)
 if(anaddb_dtset%relaxat < 0.or.anaddb_dtset%relaxat > 1)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  relaxat is',anaddb_dtset%relaxat,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct relaxat in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%relaxstr=0
 token = 'relaxstr'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%relaxstr=intarr(1)
 if(anaddb_dtset%relaxstr<0.or.anaddb_dtset%relaxstr>1)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  relaxstr is',anaddb_dtset%relaxstr,&
&  '  but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct relaxstr in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%natfix=0
 token = 'natfix'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%natfix=intarr(1)
 if(anaddb_dtset%natfix > natom)then
  write(message, '(a,a,a,i8,a,a,i4,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  natfix is',anaddb_dtset%natfix,', which is larger than natom',&
&  ' (=',natom,')',ch10,&
&  '  Action : correct natfix in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if(anaddb_dtset%natfix < 0)then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  natfix is',anaddb_dtset%natfix,', which is < 0',ch10,&
&  '  Action : correct natfix in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%nstrfix=0
 token = 'nstrfix'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%nstrfix=intarr(1)
 if(anaddb_dtset%nstrfix > 6)then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  nstrfix is',anaddb_dtset%nstrfix,', which is larger than 6',ch10,&
&  '  Action : correct nstrfix in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if(anaddb_dtset%nstrfix < 0)then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  nstrfix is',anaddb_dtset%nstrfix,', which is < 0',ch10,&
&  '  Action : correct nstrfix in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 allocate(anaddb_dtset%iatfix(natom))
 anaddb_dtset%iatfix(:) = 0
 if ((anaddb_dtset%relaxat == 1).and.(anaddb_dtset%natfix > 0)) then
  if(natom > marr)then
   marr = natom
   deallocate(intarr,dprarr)
   allocate(intarr(marr),dprarr(marr))
  end if
  token = 'iatfix'
  call intagm(dprarr,intarr,jdtset,marr,anaddb_dtset%natfix,&
&  string(1:lenstr),token,tread,'INT')
  if(tread==1) anaddb_dtset%iatfix(1:anaddb_dtset%natfix) = &
&  intarr(1:anaddb_dtset%natfix)
 end if

 if ((anaddb_dtset%relaxstr == 1).and.(anaddb_dtset%nstrfix > 0)) then
  anaddb_dtset%istrfix(:) = 0
  token = 'istrfix'
  call intagm(dprarr,intarr,jdtset,marr,anaddb_dtset%nstrfix,&
&  string(1:lenstr),token,tread,'INT')
  if(tread==1) anaddb_dtset%istrfix(1:anaddb_dtset%nstrfix) = &
&  intarr(1:anaddb_dtset%nstrfix)
 end if

 anaddb_dtset%targetpol(:) = 0._dp
 token = 'targetpol'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,&
& tread,'DPR')
 if(tread==1) anaddb_dtset%targetpol(1:3) = dprarr(1:3)

 anaddb_dtset%nfreq=1
 token = 'nfreq'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%nfreq=intarr(1)
 if(anaddb_dtset%nfreq<0)then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  nfreq is',anaddb_dtset%nfreq,', which is lower than 0 .',ch10,&
&  '  Action : correct nfreq in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%frmin=zero
 token = 'frmin'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) anaddb_dtset%frmin=dprarr(1)

 anaddb_dtset%frmax=ten
 token = 'frmax'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) anaddb_dtset%frmax=dprarr(1)

 anaddb_dtset%dipdip=1
 token = 'dipdip'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%dipdip=intarr(1)
 if(anaddb_dtset%dipdip<0.or.anaddb_dtset%dipdip>1)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  dipdip is',anaddb_dtset%dipdip,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct dipdip in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%nsphere=0
 token = 'nsphere'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%nsphere=intarr(1)
 if(anaddb_dtset%nsphere<0)then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  nsphere is',anaddb_dtset%nsphere,', which is lower than 0 .',ch10,&
&  '  Action : correct nsphere in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%rifcsph=zero
 token = 'rifcsph'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) then
  anaddb_dtset%rifcsph=dprarr(1)
 end if

 anaddb_dtset%ifcout=0
 token = 'ifcout'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%ifcout=intarr(1)
 if(anaddb_dtset%ifcout<0)then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  ifcout is',anaddb_dtset%ifcout,', which is lower than 0 .',ch10,&
&  '  Action : correct ifcout in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%iavfrq=0
 token = 'iavfrq'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%iavfrq=intarr(1)
 if(anaddb_dtset%iavfrq<0.or.anaddb_dtset%iavfrq>1)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  iavfrq is',anaddb_dtset%iavfrq,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct iavfrq in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%ifcana=0
 token = 'ifcana'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%ifcana=intarr(1)
 if(anaddb_dtset%ifcana<0.or.anaddb_dtset%ifcana>1)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  ifcana is',anaddb_dtset%ifcana,', but the only allowed values',ch10,&
&  '  are 0 or 1 .',ch10,&
&  '  Action : correct ifcana in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%natifc=0
 token = 'natifc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%natifc=intarr(1)
 if(anaddb_dtset%natifc<0)then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  natifc is',anaddb_dtset%natifc,', which is lower than 0 .',ch10,&
&  '  Action : correct natifc in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 allocate(anaddb_dtset%atifc(natom))
 anaddb_dtset%natom=natom
 anaddb_dtset%atifc(:)=0
 if(anaddb_dtset%natifc>=1)then
  if(anaddb_dtset%natifc>marr)then
   marr=anaddb_dtset%natifc
   deallocate(intarr,dprarr)
   allocate(intarr(marr),dprarr(marr))
  end if
  token = 'atifc'
  call intagm(dprarr,intarr,jdtset,marr,anaddb_dtset%natifc,string(1:lenstr),token,&
&  tread,'INT')
  if(tread==1) anaddb_dtset%atifc(1:anaddb_dtset%natifc)=intarr(1:anaddb_dtset%natifc)
 end if

 anaddb_dtset%prtmbm=0
 token = 'prtmbm'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%prtmbm=intarr(1)

 anaddb_dtset%nchan=800
 token = 'nchan'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%nchan=intarr(1)
 if(anaddb_dtset%nchan <0)then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  nchan is',anaddb_dtset%nchan,', which is lower than 0 .',ch10,&
&  '  Action : correct nchan in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%nwchan=10
 token = 'nwchan'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%nwchan=intarr(1)
 if(anaddb_dtset%nwchan<0)then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  nwchan is',anaddb_dtset%nwchan,', which is lower than 0 .',ch10,&
&  '  Action : correct nwchan in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%dosdeltae=1.0/Ha_cmm1
 token = 'dosdeltae'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) anaddb_dtset%dosdeltae=dprarr(1)
 if(anaddb_dtset%dosdeltae<=zero)then
  write(message, '(a,a,a,es14.4,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  dosdeltae is',anaddb_dtset%dosdeltae,', which is lower than 0 .',ch10,&
&  '  Action : correct dosdeltae in your input file.'
  call wrtout(6,message,'COLL') 
  call leave_new('COLL')
 end if

 anaddb_dtset%dossmear=5.0/Ha_cmm1
 token = 'dossmear'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) anaddb_dtset%dossmear=dprarr(1)
 if(anaddb_dtset%dossmear<=zero)then
  write(message, '(a,a,a,es14.4,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  dossmear is',anaddb_dtset%dossmear,', which is lower than 0 .',ch10,&
&  '  Action : correct dossmear in your input file.'
  call wrtout(6,message,'COLL') 
  call leave_new('COLL')
 end if

 anaddb_dtset%dostol=0.25_dp
 token = 'dostol'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) anaddb_dtset%dostol=dprarr(1)
 if(anaddb_dtset%dostol<zero)then
  write(message, '(a,a,a,es14.4,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  dostol is',anaddb_dtset%dostol,', which is lower than 0 .',ch10,&
&  '  Action : correct dostol in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%thmtol=0.25_dp
 token = 'thmtol'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) anaddb_dtset%thmtol=dprarr(1)
 if(anaddb_dtset%thmtol<zero)then
  write(message, '(a,a,a,es14.4,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  thmtol is',anaddb_dtset%thmtol,', which is lower than 0 .',ch10,&
&  '  Action : correct thmtol in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%ntemper=10
 token = 'ntemper'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%ntemper=intarr(1)
 if(anaddb_dtset%ntemper <0)then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  ntemper is',anaddb_dtset%ntemper,', which is lower than 0 .',ch10,&
&  '  Action : correct ntemper in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%temperinc=two
 token = 'temperinc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) anaddb_dtset%temperinc=dprarr(1)
 if(anaddb_dtset%temperinc < zero)then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  temperinc is',anaddb_dtset%temperinc,', which is lower than 0 .',ch10,&
&  '  Action : correct temperinc in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%tempermin=one
 token = 'tempermin'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) anaddb_dtset%tempermin=dprarr(1)
 if(anaddb_dtset%tempermin<-tol12)then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  tempermin is',anaddb_dtset%tempermin,', which is lower than 0 .',ch10,&
&  '  Action : correct tempermin in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%brav=1
 token = 'brav'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%brav=intarr(1)
 if(anaddb_dtset%brav<=0.or.anaddb_dtset%brav>=5)then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  brav is',anaddb_dtset%brav,', but the only allowed values',ch10,&
&  '  are 1,2,3 or 4 .',ch10,&
&  '  Action : correct brav in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%ngqpt(:)=0
 token = 'ngqpt'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%ngqpt(1:3)=intarr(1:3)
 do ii=1,3
  if(anaddb_dtset%ngqpt(ii)<0)then
   write(message, '(a,a,a,i1,a,i8,a,a,a,i1,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ngqpt(',ii,') is',anaddb_dtset%ngqpt(ii),', which is lower than 0 .',ch10,&
&   '  Action : correct ngqpt(',ii,') in your input file.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
 end do

!DEBUG
!write(6,*)' invars9 : ngqpt(:)=',anaddb_dtset%ngqpt(1:3)
!ENDDEBUG

 anaddb_dtset%nqshft=1
 token = 'nqshft'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%nqshft=intarr(1)
 if(anaddb_dtset%nqshft<0 .or. anaddb_dtset%nqshft==3 .or. anaddb_dtset%nqshft>=5 )then
  write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  nqshft is',anaddb_dtset%nqshft,', but the only allowed values',ch10,&
&  '  are 1, 2 or 4 .',ch10,&
&  '  Action : correct nqshft in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 if (anaddb_dtset%nqshft/=0)then
  if(3*anaddb_dtset%nqshft>marr)then
   marr=3*anaddb_dtset%nqshft
   deallocate(intarr,dprarr)
   allocate(intarr(marr),dprarr(marr))
  end if
  anaddb_dtset%q1shft(:,:)=zero
  token = 'q1shft'
  call intagm(dprarr,intarr,jdtset,marr,3*anaddb_dtset%nqshft,string(1:lenstr),token,tread,'DPR')
  if(tread==1) anaddb_dtset%q1shft(1:3,1:anaddb_dtset%nqshft)=&
&  reshape(dprarr(1:3*anaddb_dtset%nqshft),(/3,anaddb_dtset%nqshft/))
 end if

 anaddb_dtset%ng2qpt(:)=0
 token = 'ng2qpt'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%ng2qpt(:)=intarr(1:3)
 do ii=1,3
  if(anaddb_dtset%ng2qpt(ii)<0)then
   write(message, '(a,a,a,i1,a,i8,a,a,a,i1,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ng2qpt(',ii,') is',anaddb_dtset%ng2qpt(ii),', which is lower than 0 .',ch10,&
&   '  Action : correct ng2qpt(',ii,') in your input file.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
 end do

 anaddb_dtset%ngrids=4
 token = 'ngrids'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%ngrids=intarr(1)
 if(anaddb_dtset%ngrids<0)then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  ngrids is',anaddb_dtset%ngrids,', which is lower than 0 .',ch10,&
&  '  Action : correct ngrids in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 anaddb_dtset%q2shft(:)=zero
 token = 'q2shft'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'DPR')
 if(tread==1) anaddb_dtset%q2shft(:)=dprarr(1:3)

 anaddb_dtset%nph1l=0
 token = 'nph1l'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%nph1l=intarr(1)
 if(anaddb_dtset%nph1l<0)then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  nph1l is',anaddb_dtset%nph1l,', which is lower than 0 .',ch10,&
&  '  Action : correct nph1l in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!DEBUG
!write(6,*)' invars9 : before allocate qph1l, qnrml1, nph1l=',anaddb_dtset%nph1l
!stop
!ENDDEBUG

 allocate(anaddb_dtset%qph1l(3,anaddb_dtset%nph1l))
 allocate(anaddb_dtset%qnrml1(anaddb_dtset%nph1l))

 if (anaddb_dtset%nph1l/=0)then
  if(4*anaddb_dtset%nph1l>marr)then
   marr=4*anaddb_dtset%nph1l
   deallocate(intarr,dprarr)
   allocate(intarr(marr),dprarr(marr))
  end if
  anaddb_dtset%qph1l(:,:)=zero
  anaddb_dtset%qnrml1(:)=zero
  token = 'qph1l'
  call intagm(dprarr,intarr,jdtset,marr,4*anaddb_dtset%nph1l,string(1:lenstr),token,tread,'DPR')
  if(tread==1)then
   do iph1=1,anaddb_dtset%nph1l
    do ii=1,3
     anaddb_dtset%qph1l(ii,iph1)=dprarr(ii+(iph1-1)*4)
    end do
    anaddb_dtset%qnrml1(iph1)=dprarr(4+(iph1-1)*4)
    if(abs(anaddb_dtset%qnrml1(iph1))<qtol)then
     write(message, '(a,a,a,a,a)' )&
&     ' invars9 : ERROR -',ch10,&
&     '  The first list of wavevectors should not have non-analytical data.',ch10,&
&     '  Action : correct the first list of wavevectors in the input file.'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')
    end if
   end do
  end if
 end if

 anaddb_dtset%nph2l=0
 token = 'nph2l'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) anaddb_dtset%nph2l=intarr(1)
 if(anaddb_dtset%nph2l<0)then
  write(message, '(a,a,a,i8,a,a,a)' )&
&  ' invars9 : ERROR -',ch10,&
&  '  nph2l is',anaddb_dtset%nph2l,', which is lower than 0 .',ch10,&
&  '  Action : correct nph2l in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!DEBUG
!write(6,*)' invars9 : before allocate qph2l, qnrml2, nph2l=',anaddb_dtset%nph2l
!ENDDEBUG

 allocate(anaddb_dtset%qph2l(3,anaddb_dtset%nph2l))
 allocate(anaddb_dtset%qnrml2(anaddb_dtset%nph2l))

 if (anaddb_dtset%nph2l/=0)then
  if(4*anaddb_dtset%nph2l>marr)then
   marr=4*anaddb_dtset%nph2l
   deallocate(intarr,dprarr)
   allocate(intarr(marr),dprarr(marr))
  end if
  anaddb_dtset%qph2l(:,:)=zero
  anaddb_dtset%qnrml2(:)=zero
  token = 'qph2l'
  call intagm(dprarr,intarr,jdtset,marr,4*anaddb_dtset%nph2l,string(1:lenstr),token,tread,'DPR')
  if(tread==1)then
   do iph2=1,anaddb_dtset%nph2l
    do ii=1,3
     anaddb_dtset%qph2l(ii,iph2)=dprarr(ii+(iph2-1)*4)
    end do
    anaddb_dtset%qnrml2(iph2)=dprarr(4+(iph2-1)*4)
    if(abs(anaddb_dtset%qnrml2(iph2))>qtol)then
     write(message, '(a,a,a,a,a)' )&
&     ' invars9 : ERROR -',ch10,&
&     '  The second list of wavevectors should have only non-analytical data.',ch10,&
&     '  Action : correct the second list of wavevectors in the input file.'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')
    end if
   end do
  end if
 end if

 deallocate(dprarr,intarr)

!DEBUG
!write(6,*)' invars9 : exit '
!ENDDEBUG

end subroutine invars9
!!***
