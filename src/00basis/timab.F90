!{\src2tex{textfont=tt}}
!!****f* ABINIT/timab
!! NAME
!! timab
!!
!! FUNCTION
!! Timing subroutine.  Calls machine-dependent "timein" which
!! returns elapsed cpu and wall clock times in sec.
!!
!! Depending on value of "option" routine will:
!! (0) zero all accumulators
!! (1) start with new incremental time slice for accumulator n
!!   also increase by one the counter for this accumulator
!! (2) stop time slice; add time to accumulator n
!! (3) not used (use now time_accu)
!! (4) report time slice for accumlator n (not full time accumlated)
!! (5) option to suppress timing (nn should be 0) or reenable it (nn /=0)
!!
!! If, on first entry, subroutine is not being initialized, it
!! will automatically initialize as well as rezero accumulator n.
!! However, initialization SHOULD be done explicitly by the user
!! so that it can be done near the top of his/her main routine.
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nn=index of accumulator (distinguish what is being timed);
!!   not used if option=0
!!  option=see comment above
!!
!! OUTPUT
!!  on option=4 :
!!    tottim(2,nn)=accumulated time for accumulator nn; otherwise
!!     tottim is a dummy variable.
!!    option gives the number of times that the
!!     accumulator has been incremented
!!
!! PARENTS
!!      abinit,acfd_dyson,acfd_intexact,atm2fft,bestwfs,cchi0,cgwf,cgwf3
!!      csigme,dielmt,dieltcel,dotprod_g,dotprod_v,dotprod_vn,dotprodm_v
!!      dotprodm_vn,driver,drivergw,dyfnl3,eltfrhar3,eltfrkin3,eltfrloc3
!!      eltfrnl3,eltfrxc3,eneres3,energy,etotfor,filterpot,forces,forstrnps
!!      fourdp,fourwf,fxphas,getghc,gstate,hartre,hartre1,initylmg
!!      inkpts,invars2,inwffil,inwffil3,iofn1,iofn2,kpgio,kpgsph,leave_test
!!      lobpcgIIwf,lobpcgccIIwf,lobpcgccwf,lobpcgwf,loop3dte,loper3
!!      matrixelmt_g,mean_fftr,meanvalue_g,mkcore,mkffnl,mklocl,mkresi,mkrho
!!      mkrho3,mkvxc3,mkvxcstr3,newkpt,newocc,newrho,newvtr,newvtr3,nhatgrid
!!      nonlinear,nonlop,nonlop_pl,nonlop_ylm,nstdy3,nstwf3,opernla_ylm,orthon
!!      outkss,outscfcv,outwf,pareigocc,pawdij,pawinit,pawmkrhoij
!!      pawdenpot,pawxc,pawxcm,precon,prep_getghc,projbd,pspini,redgr,respfn
!!      rhofermi3,rhohxc_coll,rhotov,rwwf,scfcv,scfcv3,screening,setsym,setvtr
!!      sigma,sqnorm_g,sqnorm_v,sqnormm_v,status,stress,strhar,subdiago,suscep
!!      suscep_dyn,suscep_kxc_dyn,suscep_stat,susk,susk_dyn,susk_dyn_pgg
!!      susk_kxc_dyn,suskmm,suskmm_dyn,suskmm_kxc_dyn,symrhg,symsgcube,tddft
!!      timana,vtorho,vtorho3,vtorhotf,vtowfk,vtowfk3,wfkfermi3
!!      wfsinp,wrtout,xcden,xcpot
!!
!! CHILDREN
!!      leave_new,timein,wrtout
!!
!! SOURCE
#define C_INT INTEGER 
#define C_FLOAT REAL
#define C_LONG_LONG INTEGER*8

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine timab(nn,option,tottim)

 use defs_basis
 use defs_time


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis, except_this_one => timab
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

#ifdef HAVE_PAPI
#include "f90papi.h"
#else
#define PAPI_MAX_STR_LEN 0
#endif

!Arguments ------------------------------------
  !scalars
  integer,intent(in) :: nn,option
  !arrays
  real(dp),intent(out) :: tottim(2)
!Local variables-------------------------------

  !scalars
  real(dp) :: cpu,wall
  character(len=500) :: message
 character*(PAPI_MAX_STR_LEN) papi_errstr
  C_INT retval 
  C_FLOAT :: real_time, proc_time, mflops1
  C_LONG_LONG :: flops1
  !Performance mesurement
! *************************************************************************
!DEBUG
!write(6,*)' timab : enter with  nn, option, timopt',nn,option,timopt
!if(entry==5)stop
!ENDDEBUG
 if (option==5) then
! DEBUG
! write(6,*)' timab : option=5, nn=',nn
! ENDDEBUG
  timopt=nn
 end if
!If timopt was set to zero by a call with option=5, suppress
!all action of this routine (might as well return at this point !)
 if(timopt/=0 .and. option/=5)then
! 
! Check that nn lies in sensible bounds
  if (nn<0.or.nn>mtim) then
   write(message, '(a,a,a,a,i6,a,i8,a)' ) ch10,&
&   ' timab: BUG -',ch10,&
&   '  dim mtim=',mtim,' but input nn=',nn,'.'
   call wrtout(6,message,'PERS')
   call leave_new('PERS')
  end if
#ifdef HAVE_PAPI
! ie for all active options for time  
! if papi analysis has been selected or for initializing papi
! even if papi has not been choose 
  if (option/=3.and.(papiopt==1.or.initpapiopt==1)) then 
   initpapiopt=0
   call PAPIf_flops(real_time, proc_time, flops1, mflops1, retval)
   if (retval.NE.PAPI_OK) then
    PRINT *, 'Problem to initialize papi high level inteface'
    call papif_perror(retval,papi_errstr,retval)
    PRINT *, 'Error code', papi_errstr
   end if ! DEBUG
!  PRINT *, 'flops  = ', flops1, 'mflops= ',  mflops1
   if (flops1 >= HUGE(flops1)) then  
    PRINT *, 'Mflops analysis : Number of floating point instruction Overflow '
    flops(:)=-1            
   end if
  end if
#endif
! 
  if (option==0) then
!  Zero out all accumulators of time and init timers
   acctim(:,:)=0.0d0
   tzero(:,:)=0.0d0
   ncount(:)=0
   flops(:)=0
   papi_acctim(:,:)=0. 
   papi_accflops(:)=0. 
   papi_tzero(:,:)=0. 
  else if (option==1) then
!  initialize timab for nn
!  write(6,*)' timab : enter timein '
   call timein(cpu,wall)
!  write(6,*)' timab : exit timein '
   tzero(1,nn)=cpu
   tzero(2,nn)=wall
!  initialize megaflops for nn
   flops(nn)=flops1
   papi_tzero(1,nn) = proc_time
   papi_tzero(2,nn) = real_time
  else if (option==2) then
!  accumulate time for nn
!  write(6,*)' timab : enter timein '
   call timein(cpu,wall)
!  write(6,*)' timab : exit timein '
   acctim(1,nn)=acctim(1,nn)+cpu -tzero(1,nn)
   acctim(2,nn)=acctim(2,nn)+wall-tzero(2,nn)
!  accumulate time and flops for nn Difference entre 2 calls a Papif_flops 
   papi_acctim(1,nn)=papi_acctim(1,nn)+ proc_time - papi_tzero(1,nn)
   papi_acctim(2,nn)=papi_acctim(2,nn)+ real_time - papi_tzero(2,nn)
   papi_accflops(nn)=papi_accflops(nn)+ flops1- flops(nn) 
   ncount(nn)=ncount(nn)+1
  else if (option==3) then
!  Now use time_accu
   write(message, '(a,a,a,a)' ) ch10,&
&   ' timab: BUG -',ch10,&
&   ' option 3 not valid (use time_accu).'
   call wrtout(6,message,'PERS')
   call leave_new('PERS')
  else if (option==4) then
!  return elapsed time for nn (do not accumulate)
   call timein(cpu,wall)
   tottim(1)=cpu-tzero(1,nn)
   tottim(2)=wall-tzero(2,nn)
!  return ellapsed floating point operationfor nn (do not accumulate)
   papi_tottim(1,nn)= proc_time - papi_tzero(1,nn)
   papi_tottim(2,nn)= real_time - papi_tzero(2,nn)
   papi_totflops(nn)= flops1- flops(nn) 
  else
   write(message, '(a,a,a,a,i10,a)' ) ch10,&
&   ' timab: BUG -',ch10,&
&   '  Input option not valid, =',option,'.'
   call wrtout(6,message,'PERS')
   call leave_new('PERS')
  end if
 end if
!DEBUG
!write(6,*)' timab : exit '
!ENDDEBUG
end subroutine timab
!!***
