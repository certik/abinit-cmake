!{\src2tex{textfont=tt}}
!!****f* ABINIT/status
!!
!! NAME
!! status
!!
!! FUNCTION
!! Routine for description of the status of the calculation
!! Eventually open the status file, write different information,
!! and close the file. The output rate is governed by istatr
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG,TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  counter=value of the loop counter at that level
!!  file (optional argument)=name of the status file
!!  istatr=gives the rate of output. The status file will be opened
!!     and written only once every "istatr" calls.
!!     This variable is saved at the fifth call (just after the first
!!     call to invars. In preceeding and subsequent
!!     calls, istatr has no meaning.
!!  level=number of the level of the calling subroutine
!! (1=initialisation in abinit ; 2=driver ; 3=gstate ; 4=move ;
!!  5=brdmin ; 6=scfcv ; 7=vtorho ; 8=vtowfk ; 9=cgwf ; 10=getghc ;
!!  13=respfn ; 14=loper3 ; 15=moldyn ; 23=screening ; 24=cchi0/cchi0q0 ; 
!!  25=sigma ; 26=csigme)
!!  routine=string of 14 characters indicating the status inside the level
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! Warning : The string "routine" can have any size but
!! it is truncated to a size of 14.
!! because of the behaviour of some compilers, the string
!! "routine" should always have 14 characters in the calling subroutine
!!
!! PARENTS
!!      abinit,afterscfloop,brdmin,cchi0,cchi0q0,cgwf,cgwf3,csigme,ctocprj
!!      delocint,driver,eig1fixed,getgh1c,getghc,getgsc,gstate,hessinit
!!      loop3dte,loper3,moldyn,move,mv_3dte,nonlinear,pstate,resp3dte,respfn
!!      rhofermi3,scfcv,scfcv3,screening,setup_positron,sigma,subdiago,suscep
!!      vtorho,vtorho3,vtorhorec,vtorhotf,vtowfk,vtowfk3,wannier,wfkfermi3
!!
!! CHILDREN
!!      leave_new,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine status(counter,filstat,istatr,level,routine)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: counter,istatr,level
 character(len=*),intent(in) :: routine
 character(len=fnlen),intent(in) :: filstat

!Local variables-------------------------------
!scalars
 integer,parameter :: mcounter=2,mlevel=26
 integer,save :: first=1,output_rate=1,shift_rate=1,statnu=0
 integer :: ilevel
 character(len=12) :: headwr
 character(len=500) :: message
!arrays
 integer,save :: active(mlevel),actual_counter(mlevel,mcounter)
 integer,save :: ncounter(mlevel)
 real(dp) :: tsec(2)
 character(len=14),save :: nm_levels(mlevel),nm_routine(mlevel)
 character(len=9),save :: nm_counter(mlevel,mcounter)

!***********************************************************************

 call timab(73,1,tsec)

!Note : all processors have their own file, so no special
!attention must be paid to the parallel case.

!Initialisation
 if(first/=0)then
  if(first==1)then
   first=5
   nm_routine(:)='              '
   active(:)=0
   actual_counter(:,:)=0

!  List of names for each level
   nm_levels(1)   ='abinit        '
   ncounter(1)=0
   nm_counter(1,1)='         '

   nm_levels(2)   ='driver        '
   ncounter(2)=1
   nm_counter(2,1)='jdtset  ='

   nm_levels(3)   ='gstate        '
   ncounter(3)=1
   nm_counter(3,1)='itime   ='

   nm_levels(4)   ='move          '
   ncounter(4)=2
   nm_counter(4,1)='icalls  ='
   nm_counter(4,2)='itime   ='

   nm_levels(5)   ='brdmin/moldyn '
   ncounter(5)=1
   nm_counter(5,1)='itime   ='

   nm_levels(6)   ='scfcv         '
   ncounter(6)=1
   nm_counter(6,1)='istep   ='

   nm_levels(7)   ='vtorho        '
   ncounter(7)=2
   nm_counter(7,1)='isppol  ='
   nm_counter(7,2)='ikpt    ='

   nm_levels(8)   ='vtowfk        '
   ncounter(8)=2
   nm_counter(8,1)='inonsc  ='
   nm_counter(8,2)='iband   ='

   nm_levels(9)   ='cgwf          '
   ncounter(9)=1
   nm_counter(9,1)='iline   ='

   nm_levels(10)   ='getghc        '
   ncounter(10)=0
   nm_counter(10,1)='         '

   nm_levels(13)   ='respfn        '
   ncounter(13)=0
   nm_counter(13,1)='         '

   nm_levels(14)   ='loper3        '
   ncounter(14)=1
   nm_counter(14,1)='respcase='

   nm_levels(16)   ='scfcv3        '
   ncounter(16)=1
   nm_counter(16,1)='istep   ='

   nm_levels(17)   ='vtorho3       '
   ncounter(17)=2
   nm_counter(17,1)='isppol  ='
   nm_counter(17,2)='ikpt    ='

   nm_levels(18)   ='vtowfk3       '
   ncounter(18)=2
   nm_counter(18,1)='inonsc  ='
   nm_counter(18,2)='iband   ='

   nm_levels(19)   ='cgwf3         '
   ncounter(19)=1
   nm_counter(19,1)='iline   ='

   nm_levels(20)   ='nonlinear     '
   ncounter(20)=0
   nm_counter(20,1)='         '

   nm_levels(21)   ='loop3dte      '
   ncounter(21)=2
   nm_counter(21,1)='pert1case ='
   nm_counter(21,2)='pert3case ='

   nm_levels(22)   ='mv_/resp3dte  '
   ncounter(22)=2
   nm_counter(22,2)='ikpt ='

   nm_levels(23)   ='screening     '
   ncounter(23)=1
   nm_counter(23,1)='iqpt    ='

   nm_levels(24)   ='cchi0/cchi0q0 '
   ncounter(24)=1
!  nm_counter(24,1)='isppol  ='
   nm_counter(24,1)='ikpt    ='

   nm_levels(25)   ='sigma         '
   ncounter(25)=1
   nm_counter(25,1)='ikpt_gw ='

   nm_levels(26)   ='csigme        '
   ncounter(26)=1
!  nm_counter(26,1)='isppol  ='
   nm_counter(26,1)='ikpt    ='



  else if(first==5)then

!  The true input variable "output_rate" is only available at the fifth
!  call to "status".
   if(statnu+1==5)then
    first=0
    if(istatr<=0)then
     write(message, '(a,a,a,a,i7,a,a,a,a,a)' ) ch10,&
&     ' status : ERROR -',ch10,&
&     '  the value of the input variable istatr is',istatr,' ,',ch10,&
&     '  while it must be a positive, non-zero number.',ch10,&
&     '  Action : change istatr in your input file.'
     call wrtout(6,message,'COLL')
     call leave_new('PERS')
    end if
    output_rate=istatr
   end if

  end if
 end if

!The true input variable "shift_rate" is only available at the sixth call
 if(statnu+1==6)then
  if(istatr<0 .or. istatr>=output_rate)then
   write(message, '(a,a,a,a,i7,a,a,a,a,a)' ) ch10,&
&   ' status : ERROR -',ch10,&
&   '  the value of the input variable istatshft is',istatr,' ,',ch10,&
&   '  while it must be a positive number smaller than istatr.',ch10,&
&   '  Action : change istatshft in your input file.'
   call wrtout(6,message,'COLL')
   call leave_new('PERS')
  end if
  shift_rate=istatr
 end if

!Check the value of level
 if(level>mlevel)then
  write(message, '(a,a,a,a,i5,a,a,a,i5)' ) ch10,&
&  ' status : BUG -',ch10,&
&  '  The value of level in the calling routine is',level,' ,',ch10,&
&  '  while the maximum allowed value is',mlevel
  call wrtout(6,message,'COLL')
  call leave_new('PERS')
 end if

!Assign the info about the actual routine
 write(unit=nm_routine(level),fmt='(a14)') routine
 if(trim(nm_routine(level))=='exit')then
! The value 2 will be changed to 0 at the end of the routine.
  active(level)=2
 else if(trim(nm_routine(level))=='')then
  active(level)=0
 else
  active(level)=1
 end if

!Assign the info about the actual counter
 if(counter>=0)then
  if(ncounter(level)==1)then
   actual_counter(level,1)=counter
  else if(ncounter(level)==2)then
   actual_counter(level,2)=counter/100
!  The counter number 1 does not allow more than 99 passes
   actual_counter(level,1)=counter-(counter/100)*100
  end if
 end if

!============================================================

!After treatment of present information, output of the status
 statnu=statnu+1

!DEBUG
!write(6,*)' status : statnu=',statnu
!ENDDEBUG

 if( mod(statnu,output_rate)==shift_rate .or. statnu==10 )then

  open (tmp_unit,file=trim(filstat),form='formatted',status='unknown')
  rewind tmp_unit

  write(tmp_unit,*)

  headwr='(a,i4,a,i6 )'
  if(statnu>=100000)   headwr='(a,i4,a,i9 )'
  if(output_rate>=1000)headwr='(a,i6,a,i6 )'
  if(statnu>=100000 .and. output_rate>=1000)   headwr='(a,i6,a,i9 )'
  if(statnu>=100000000)headwr='(a,i6,a,i12)'
  write(tmp_unit,headwr)&
&  ' Status file, with repetition rate',output_rate,&
&  ', status number',statnu
  write(tmp_unit,*)

! Treat every level one after the other
  do ilevel=1,mlevel
!  This level must be activated in order to have a corresponding output
   if(active(ilevel)==1 .or. active(ilevel)==2)then

    write(tmp_unit,'(4a)')&
&    '  Level ',nm_levels(ilevel),' : ',nm_routine(ilevel)

!   Is there a counter for this level ?
    if(ncounter(ilevel)>=1)then

     if(actual_counter(ilevel,1)>0)then
      write(tmp_unit,'(a,a,i5)')'  ',nm_counter(ilevel,1),&
&      actual_counter(ilevel,1)
     end if
     if(ncounter(ilevel)==2)then
      if(actual_counter(ilevel,2)>0)then
       write(tmp_unit,'(a,a,i5)')'  ',nm_counter(ilevel,2),&
&       actual_counter(ilevel,2)
      end if
     end if

    end if

!   End of the check on activation of the level
   end if

!  End of the loop on the levels
  end do

  close (tmp_unit)

! End of the repetition rate check
 end if

 if (active(level)==2)then
  active(level)=0
  nm_routine(level)='              '
 end if

 call timab(73,2,tsec)

end subroutine status
!!***
