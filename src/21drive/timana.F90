!{\src2tex{textfont=tt}}
!!****f* ABINIT/timana
!! NAME
!! timana
!!
!! FUNCTION
!! Analyse the timing, and print in unit ab_out. Some discussion of the
!! number of calls to different routines is also provided, as comments,
!! at the end of the routine, as well as, in the single dataset mode (ndtset<2),
!! a detailed analysis of the time-consuming routines.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell.
!!  nband(nkpt*nsppol)=number of bands at each k point, for each polarization
!!  ndtset=number of datasets
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nkpt=number of k points
!!  npwtot(nkpt)=number of planewaves in basis at this k point
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  timopt= if >0, write short analysis, if <0, write full analysis
!!          if abs(timopt)>=2, do not time the timer
!!
!!  papiopt=1 write an analysis of time and speed execution done with the papi library
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! One can suppress the cpu timer call in timein.f, if line 315
!! of the present routine is uncommented.
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      timab,time_accu,wrtout,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine timana(mpi_enreg,natom,nband,ndtset,nfft,nkpt,npwtot,nsppol,timopt, papiopt)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ndtset,nfft,nkpt,nsppol,papiopt,timopt
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: nband(nkpt*nsppol),npwtot(nkpt)

!Local variables-------------------------------
!scalars
 integer,parameter :: mtim=599
 integer :: allsub,ierr,ii,ikpt,ilink,ilist,isppol,itim,itimab,ltimab,maxii,me
 integer :: nlink,nlist,nothers,nproc,npwmean,npwnbdmean,return_ncount
 integer :: spaceworld,temp_list,totcount,utimab
 real(dp) :: cpunm,lflops,other_cpu,other_wal,percent_limit,subcpu,subwal
 real(dp) :: timab_cpu,timab_wall,wallnm
 character(len=500) :: message
!arrays
 integer :: ncount(mtim),ncount_s(mtim),ndata(mtim)
 integer,allocatable :: list(:)
 real(dp) :: ftimes(2,mtim),ftsec(2),mflops(mtim),nflops(mtim),tcount(2)
 real(dp) :: times(2,mtim),times_s(2,mtim),tsec(2),tsec_s(2)
 character(len=16) :: names(mtim)
!no_abirules
 character(len=*), parameter :: &
       &  format01040 ="('- ',a16,f12.3,f6.1,f11.3,f6.1,i15)",&
       &  format01041 ="('- ',a16,f12.3,f6.1,f11.3,g12.3,i15)",&
       &  format01045 ="('-',i3,a11,f15.3,f6.1,f11.3,f6.1)",&
       &  format01200 ="('- subtotal     ',f15.3,f6.1,f11.3,f6.1)"
01200 format('- subtotal     ',f15.3,f6.1,f11.3,f6.1)

! *************************************************************************

!DEBUG
!write(6,*)' timana : enter, timopt= ',timopt
!if(.true.)stop
!ENDDEBUG

!The means are computed as integers, for later compatibility
 npwmean=0
 npwnbdmean=0
 do isppol=1,nsppol
  do ikpt=1,nkpt
   npwmean=npwmean+npwtot(ikpt)
   npwnbdmean=npwnbdmean+npwtot(ikpt)*nband(ikpt+(isppol-1)*nkpt)
!  DEBUG
!  write(6,*)'ikpt,isppol,nband',ikpt,isppol,nband(ikpt+(isppol-1)*nkpt)
!  ENDDEBUG
  end do
 end do

 npwmean=dble(npwmean)/dble(nkpt*nsppol)
 npwnbdmean=dble(npwnbdmean)/dble(nkpt*nsppol)

!List of timed subroutines, eventual initialisation of the number of data
!Channels 1 to 299 are for optdriver=0 (GS), 1 (RF) and 2 (Suscep), at random
!Channels 300 to 399 are for optdriver=3 (Chieps)
!Channels 400 to 499 are for optdriver=4 (Sigma)
!Channels 500 to 529 are for optdriver=5 (Nonlinear)
!Channels 530 to 549 are for various counters
!Channels 550 to 599 are for PAW
 names(1:mtim)='***             '
 names(1)='abinit          '
 names(2)='fourwf(pot)     '           ;    ndata(2)=2*nfft
 names(3)='fourwf(den)     '           ;    ndata(3)=nfft
 names(4)='projbd          '           ;    ndata(4)=npwnbdmean
 names(5)='ewald           '
 names(6)='setsym          '
 names(7)='vtowfk(ssdiag)  '
 names(8)='vtowfk(contrib) '
 names(9)='fourdp          '           ;    ndata(9)=nfft
 names(10)='hartre          '
 names(11)='xc:pot/=fourdp  '          ;    ndata(11)=nfft*nsppol
 names(12)='mkcore          '
 names(13)='mkresi          '
 names(14)='rwwf            '
 names(15)='pspini          '
 names(16)='mkffnl          '
 names(17)='symrhg(no FFT)  '
 names(18)='vtorho  (1)     '
 names(19)='inwffil         '
 names(20)='scfcv           '
 names(21)='vtorho          '
 names(22)='cgwf            '
 names(23)='kpgsph          '
 names(24)='vtorho  (1)(2)  '
 names(25)='vtorho  (2)     '
 names(26)='vtorho-kpt loop '
 names(27)='vtorho  (4)     '
 names(28)='vtowfk          '
 names(29)='vtorho:MPI      '
 names(30)='vtowfk  (3)     '
 names(31)='vtowfk  (1)     '
 names(32)='gstate          '
 names(33)='gstate->kpgsph  '
 names(34)='gstate  (2)     '
 names(35)='test_ionmov     '
 names(36)='gstate  (3)     '
 names(37)='stress          '
 names(38)='ewald2          '
 names(39)='vtowfk (loop)   '
 names(40)='cgwf-O(npw)     '
 names(41)='abinit(1)       '
 names(42)='abinit(2)       '
 names(43)='invars2m        '
 names(44)='abinit(4)       '
 names(45)='abinit(5)       '
 names(46)='abinit(6)       '
 names(47)='ingeo/symgroup  '
 names(48)='communic.MPI    '
 names(50)='timing timab    '
 names(51)='total timab     '
 names(52)='scfcv-scprqt    '
 names(53)='forces-mkcore   '
 names(54)='scfcv   (1)     '
 names(55)='stress-mkcore   '
 names(56)='scfcv-read      '
 names(57)='rhotov          '
 names(58)='newvtr/newrho   '
 names(59)='energy          '
 names(60)='scfcv   (6)     '
 names(61)='scfcv :synchro  '
 names(62)='kpgio :synchro  '
 names(63)='mkrho :synchro  '
 names(64)='strkin:synchro  '
 names(65)='forstrnps:synchr'
 names(66)='vtorho:synchro  '
 names(67)='inwffil:synchro '
 names(68)='fourwf(total )  '
 names(69)='forces          '
 names(70)='vtorho(symrhg)  '
 names(71)='mkrho :MPIrhor  '
 names(72)='mklocl(2)       '
 names(73)='status          '
 names(74)='newocc          '
 names(75)='nonlop(apply)   '           ;    ndata(75)=npwmean*natom
 names(76)='nonlop(forces)  '           ;    ndata(76)=npwmean*natom
 names(77)='nonlop(forstr)  '           ;    ndata(77)=npwmean*natom
 names(78)='nonlop(dyfrnl)  '
 names(79)='nonlop(ddk)     '
 names(80)='etotfor/=forces '
 names(81)='xc:pot          '    ! rhohxc_coll, except the call to hartre.f
 names(82)='xc:fourdp       '
 names(83)='newvtr/rho(3):io'
 names(84)='suscep          '
 names(85)='suscep:MPI      '
 names(86)='suscep:synchro  '
 names(87)='susk(m:loop(1)  '
 names(88)='susk(m:loop(2)  '
 names(89)='suscep:other    '
 names(90)='dielmt          '
 names(91)='setvtr          '
 names(92)='setvtr:mkcore   '
 names(93)='fourwf(G->r)    '
 names(94)='fourwf(r->G)    '
 names(95)='tddft           '
 names(96)='dieltcel        '
 names(97)='nonlop(total)   '
 names(98)='getghc-other    '
 names(99)='vtorho(4)-mkrho-paw'

 names(100)='driver          '

 names(101)='nstdy3          '
 names(102)='nstwf3          '
 names(108)='vtowfk3(contrib)'
 names(118)='vtorho3 (1)     '
 names(120)='scfcv3          '
 names(121)='vtorho3         '
 names(122)='cgwf3           '
 names(124)='vtorho3 (1)(2)  '
 names(125)='vtorho3 (2)     '
 names(126)='vtorho3-kpt loop'
 names(127)='vtorho3 (4)     '
 names(128)='vtowfk3         '
 names(129)='vtorho3:MPI     '
 names(130)='vtowfk3 (3)     '
 names(131)='vtowfk3 (1)     '
 names(132)='respfn          '
 names(133)='respfn(kpgio)   '
 names(134)='respfn(pspini)  '
 names(135)='respfn(inwffil) '
 names(136)='respfn(frozen)  '
 names(137)='respfn(1st)     '
 names(138)='respfn(init 1st)'
 names(139)='vtowfk3 (loop)  '
 names(140)='cgwf3-O(npw)    '
 names(141)='loper3          '
 names(142)='loper3(kpgio)   '
 names(143)='loper3(getmpw)  '
 names(144)='loper3(inwffil) '
 names(145)='loper3(other)   '
 names(146)='loper3(outwf)   '
 names(147)='nstdy3          '
 names(152)='scfcv3-scprqt   '
 names(154)='scfcv3  (1)     '
 names(157)='rhotov3         '
 names(158)='newvtr3         '
 names(159)='dyfnl3          '
 names(160)='scfcv3  (6)     '
 names(161)='nstdy3:synchro  '
 names(166)='vtorho3:synchro '
 names(181)='mkvxc3          '

 names(191)='invars2         '
 names(192)='inkpts          '

 names(200)='getghc          '
 names(201)='getghc%cgwf     '
 names(202)='getghc%cgwf3    '
 names(203)='getghc%mkresi   '
 names(204)='getghc%lobpcgwf '
 names(205)='getghc%lobpcgccw'
 names(206)='getghc%prep_getg'
 names(207)='getgh1c%cgwf3   '

 names(210)='projbd          '         ;    ndata(210)=npwnbdmean
 names(211)='projbd%cgwf     '
 names(212)='projbd%cgwf3    '
 names(213)='projbd%vtowfk3  '

 names(220)='nonlop%(other)  '
 names(221)='nonlop%getghc   '
 names(222)='nonlop%vtowfk   '
 names(223)='nonlop%energy   '
 names(224)='nonlop%forstrnps'
 names(225)='nonlop%nstwf3   '
 names(226)='nonlop%dyfnl3   '
 names(227)='nonlop%cgwf3 !2  '
 names(228)='nonlop%cgwf3 !5  '
 names(229)='nonlop%outkss   '
 names(230)='nonlop%rhoij    '

 names(231)='suscep_dyn      '
 names(232)='suscep_stat     '

 names(233)='outkss         '
 names(234)='outkss(1)      '
 names(235)='outkss(k-loop1)'
 names(236)='outkss(diago)  '
 names(237)='outkss(k-loop2)'

 names(240)='fourwf%(other)  '
 names(241)='fourwf%getghc   '
 names(242)='fourwf%vtowfk   '
 names(243)='fourwf%mkrho    '
 names(244)='fourwf%cgwf3    '
 names(245)='fourwf%vtowfk3  '
 names(246)='fourwf%mkrho2   '
 names(247)='fourwf%mkrho3   '
 names(248)='fourwf%susk !0   '
 names(249)='fourwf%susk !3   '
 names(250)='fourwf%susk_dyn0'
 names(251)='fourwf%susk_dyn3'
 names(252)='fourwf%suskmm !0 '
 names(253)='fourwf%suskmm !3 '
 names(254)='fourwf%tddft    '
 names(255)='fourwf%outkss  '
 names(256)='fourwf%prep_four'
 names(257)='fourwf%susk _PAW '

 names(260)='fourdp%(other)  '
 names(261)='fourdp%rhotwg%ch'
 names(262)='fourdp%rhotwg%si'
 names(263)='fourdp%ckxcldag '
 names(264)='fourdp%fftwfn%ch'
 names(265)='fourdp%fftwfn%si'

 names(270)='rwwf%(other)    '
 names(271)='rwwf%vtorho     '
 names(272)='rwwf%initwf     '
 names(273)='rwwf%energy     '
 names(274)='rwwf%inwffil(GS)'
 names(275)='rwwf%mkrho      '
 names(276)='rwwf%outwf      '
 names(277)='rwwf%strnps     '
 names(278)='rwwf%tddft      '
 names(279)='rwwf%suscep     '
 names(280)='rwwf%suscep_dyn '
 names(281)='rwwf%inwffil(RF)'
 names(282)='rwwf%mkrho2     '
 names(283)='rwwf%outwf2     '
 names(284)='rwwf%dyfnl3     '
 names(285)='rwwf%mkrho3     '
 names(286)='rwwf%nstwf3     '
 names(287)='rwwf%vtorho3    '
 names(288)='rwwf%vtowfk3    '
 names(289)='rwwf%nstdy3     '

 names(291)='mkrho%gstate    '
 names(292)='mkrho%vtorho    '
 names(293)='mkrho%energy    '
 names(294)='mkrho%respfn    '
 names(298)='mkrho/=         '
 names(299)='mkrho/=+fourwf  '

 names(301)='screening tot   '
 names(302)='screening(1)    '
 names(303)='screening(fftwfn'
 names(304)='KS => QP; [wfrg]'
 names(305)='screening(2)    '
 names(306)='overall q-loop  '
 names(307)='cchi0q0         '
 names(308)='cchi0           '
 names(309)='chi => eps      '
 names(310)='write scr files '
 names(320)='screenin/=fourdp'

 names(331)='cchi0           '
 names(332)='cchi0(rho_tw_g) '
 names(333)='cchi0(assembly) '
!do not use slot 399

 names(401)='sigma           '
 names(402)='sigma(1)        '
 names(403)='sigma(rdkss)    '
 names(404)='sigma(2)        '
 names(405)='sigma(csigme)   '
 names(410)='sigma/=fourdp   '

 names(421)='csigme (tot)    '
 names(422)='csigme(init0)   '
 names(423)='csigme(initq))  '
 names(424)='csigme(SigX)    '
 names(425)='csigme(SigC)    '
 names(426)='csigme(b-loop)  '

 names(501)='nonlinear       '
 names(502)='loop3dte        '
 names(521)='mv_3dte         '
 names(522)='resp3dte        '

 names(530)='lobpcg          '
 names(531)='nonlop%lobpcg   '
 names(532)='xgemm%lobpcg    '
 names(533)='xsum_mpi%lobpcg '
 names(534)='prep_getghc%lobp'
 names(535)='xorthon-xtrsm   '
 names(536)='xprecon%lobpcg  '
 names(537)='prep_fourwf%vtow'
 names(538)='prep_fourwf%mkrh'
 names(539)='prep_fourwf     '

 names(540)='sg_fourwf%fourwf'
 names(541)='back_wf%sg_fourw'
 names(542)='forw_wf%sg_fourw'
 names(543)='alltoall%back_wf'
 names(544)='alltoall%forw_wf'
 names(545)='alltoall%prep_ge'
 names(546)='allgather%prep_g'
 names(547)='alltoall%prep_fo'
 names(548)='allgather%prep_f'
 names(549)='symrhg%mkrho   '

 names(550)='forces:pawatm2ff'
 names(551)='stress:pawatm2ff'
 names(552)='setvtr:pawatm2ff'
 names(553)='pawinit         '
 names(554)='vtowfk:rhoij    '
 names(555)='vtorho:rhoij    '
 names(556)='vtorho:mknhat   '
 names(557)='vtorho:ntild+nhat'
 names(558)='scfcv:mknhat    '
 names(559)='nhatgrid        '
 names(560)='pawdenpot       '
 names(561)='pawdij          '
 names(562)='respfn:pawatm2ff'
 names(563)='dyfro3:pawatm2ff'
 names(564)='scfcv3:mknhat3  '
 names(565)='vtorho3:symrhoij3'
 names(566)='vtorho3:mknhat3 '
 names(567)='vtorho3:tild+hat'

 names(570)='prep_nonlop     '
 names(571)='prep_nonlop%lobp'
 names(572)='prep_nonlop%vtow'
 names(573)='prep_nonlop%fors'

 names(580)='nonlop%prep_nonl'
 names(581)='alltoall%prep_no'
 names(582)='allgather%prep_n'
 names(583)='orthon          '
 names(584)='xcopy%lobpcg    '
 names(585)='subdiago        '
 names(586)='nonlocalpart    '
 names(587)='zheegv-dsyegv   '

 me=mpi_enreg%me
 nproc=mpi_enreg%nproc
 spaceworld=mpi_enreg%world_comm

 if(abs(timopt)==1)then
! Time the timing routine (precision should be better than 3%)
  ltimab=1
  utimab=1000
  maxii=20
! maxii=1    ! Uncomment this line if no timer is provided in timein.f
  do ii=1,20

   call timab(50,1,tsec)
   do itimab=ltimab,utimab
!   The channel 51 is here used as a dummy channel
    call timab(51,1,tsec)
    call timab(51,2,tsec)
   end do
   call timab(50,2,tsec)
   call time_accu(50,return_ncount,tsec, lflops, ftsec)
!  Exit the timing loop if the CPU time is bigger than 0.10 second
!  of if the number of calls is too large.
!  Since the accuracy of the timing is expected to be better than 0.01 sec,
!  gives about 10% accuracy
   if(tsec(1)>0.10_dp)then
    exit
   else
    ltimab=utimab+1
!   Increase the number of timab calls in a block.
!   This small factor of increase allows to have less than
!   0.15 second for this testing
    utimab=(3*utimab)/2
   end if
  end do
! Get the time per combined call timab(*,1,tsec) + timab(*,2,tsec)
  timab_cpu=tsec(1)/utimab
  timab_wall=tsec(2)/utimab
  if(timopt<0 .and. me==0)then
   write(ab_out,*)
   write(ab_out,*)'Test the timer : '
   write(ab_out,*)' a combined call timab(*,1,tsec) + timab(*,2,tsec) is '
   write(ab_out, '(a,es14.4,a,es14.4,a)' )&
&   '- CPU time =',timab_cpu,' sec,    Wall time =',timab_wall,' sec'
  end if
 else
  timab_cpu=0.0_dp
  timab_wall=0.0_dp
 end if

!DEBUG
!write(6,*)' timana : after time the timer  '
!if(.true.)stop
!ENDDEBUG

!Eventually reenable the timab routine
 call timab(1,5,tsec)

!Get overall elapsed cpu and wall clock time
 call timab(1,2,tsec)
 call time_accu(1,return_ncount,tsec,lflops,ftsec)
 ncount(1)=return_ncount
!Sum over all procs

 tsec_s(:)=tsec(:)
 call xsum_mpi(tsec_s,tsec,2,spaceworld,ierr)
!Only the world master writes
 if(me==0) then
  write(ab_out, '(/,a,f13.1,f12.2,f11.3)' ) &
&  '- Total cpu        time (s,m,h):',&
&  tsec(1),tsec(1)/60._dp,tsec(1)/3600._dp
  write(ab_out, '(a,f13.1,f12.2,f11.3)' ) &
&  '- Total wall clock time (s,m,h):',&
&  tsec(2),tsec(2)/60._dp,tsec(2)/3600._dp
 end if

!Get separate time reports from all timed sections
 totcount=0
 do itim=1,mtim
  call time_accu(itim,return_ncount,times(:,itim),nflops(itim), ftimes(:,itim))
  ncount(itim)=return_ncount
  totcount=totcount+return_ncount
! DEBUG
! write(6,*)itim,times(1:2,itim),ncount(itim),totcount
! ENDDEBUG
 end do

!Estimate additional timings.

!Estimate the values associated with timab, put it in channel 51
 ncount(51)=totcount
 times(1,51)=timab_cpu*totcount
 times(2,51)=timab_wall*totcount

!Gather the different parts of nonlop
 ncount(75)=ncount(221)+ncount(223)+ncount(229)
 times(1:2,75)=times(1:2,221)+times(1:2,223)+times(1:2,229)
 nflops(75) = nflops(221) + nflops(223) + nflops( 229)
 ftimes(1:2,75)=ftimes(1:2,221)+ftimes(1:2,223)+ftimes(1:2,229)

 ncount(76)=ncount(222)+ncount(225)+ncount(227)
 times(1:2,76)=times(1:2,222)+times(1:2,225)+times(1:2,227)
 nflops(76)=nflops(222)+nflops(225)+nflops(227)
 ftimes(1:2,76)=ftimes(1:2,222)+ftimes(1:2,225)+ftimes(1:2,227)

 ncount(77)=ncount(224)
 times(1:2,77)=times(1:2,224)
 ftimes(1:2,77)=ftimes(1:2,224)
 nflops(77)=nflops(224)

 ncount(78)=ncount(226)
 times(1:2,78)=times(1:2,226)
 ftimes(1:2,78)=ftimes(1:2,226)
 nflops(78)=nflops(226)

 ncount(79)=ncount(228)
 times(1:2,79)=times(1:2,228)
 ftimes(1:2,79)=ftimes(1:2,228)
 nflops(79)=nflops(228)

 ncount(97)=ncount(75)+ncount(76)+ncount(77)+ncount(78)+ncount(79)+ncount(220)
 times(1:2,97)=times(1:2,75)+times(1:2,76)+times(1:2,77)&
& +times(1:2,78)+times(1:2,79)+times(1:2,220)
 ftimes(1:2,97)=ftimes(1:2,75)+ftimes(1:2,76)+ftimes(1:2,77)&
& +ftimes(1:2,78)+ftimes(1:2,79)+ftimes(1:2,220)
 nflops(97)=nflops(75)+nflops(76)+nflops(77)&
& +nflops(78)+nflops(79)+nflops(220)


!Gather the different parts of fourwf
!NOTE : Should attribute the channel 240 to one of the 4 modes !!!
 ncount(93)=ncount(245)+ncount(247)+ncount(248)+&
& ncount(250)+ncount(252)+ncount(254)
 times(1:2,93)=times(1:2,245)+times(1:2,247)+times(1:2,248)+&
& times(1:2,250)+times(1:2,252)+times(1:2,254)
 ftimes(1:2,93)=ftimes(1:2,245)+ftimes(1:2,247)+ftimes(1:2,248)+&
& ftimes(1:2,250)+ftimes(1:2,252)+ftimes(1:2,254)
 nflops(93)=nflops(245)+nflops(247)+nflops(248)+&
& nflops(250)+nflops(252)+nflops(254)


 ncount(3)=ncount(242)+ncount(243)+ncount(246)
 times(1:2,3)=times(1:2,242)+times(1:2,243)+times(1:2,246)
 ftimes(1:2,3)=ftimes(1:2,242)+ftimes(1:2,243)+ftimes(1:2,246)
 nflops(3)=nflops(242)+nflops(243)+nflops(246)


 ncount(2)=ncount(241)+ncount(244)
 times(1:2,2)=times(1:2,241)+times(1:2,244)
 ftimes(1:2,2)=ftimes(1:2,241)+ftimes(1:2,244)
 nflops(2)=nflops(241)+nflops(244)

 ncount(94)=ncount(249)+ncount(251)+ncount(253)+ncount(257)
 times(1:2,94)=times(1:2,249)+times(1:2,251)+times(1:2,253)+times(1:2,257)
 ftimes(1:2,94)=ftimes(1:2,249)+ftimes(1:2,251)+ftimes(1:2,253)
 nflops(94)=nflops(249)+nflops(251)+nflops(253)

!In the following, the part coming from the prep_fourwf interface
!is added to the total.
 ncount(68)=ncount(2)+ncount(3)+ncount(93)+ncount(94)+ncount(240)+ncount(256)
 times(1:2,68)=times(1:2,2)+times(1:2,3)+times(1:2,93)+times(1:2,94)+times(1:2,240)+times(1:2,256)
 ftimes(1:2,68)=ftimes(1:2,2)+ftimes(1:2,3)+ftimes(1:2,93)+ftimes(1:2,94)+ftimes(1:2,240)+ftimes(1:2,256)
 nflops(68)=nflops(2)+nflops(3)+nflops(93)+nflops(94)+nflops(240)+nflops(256)

!Gather the different parts of prep_nonlop
 ncount(570)=ncount(571)+ncount(572)
 times(1:2,570)=times(1:2,571)+times(1:2,572)
 ftimes(1:2,570)=ftimes(1:2,571)+ftimes(1:2,572)
 nflops(570)=nflops(571)+nflops(572)

!Gather the different parts of prep_fourwf
 ncount(539)=ncount(537)+ncount(538)
 times(1:2,539)=times(1:2,537)+times(1:2,538)
 ftimes(1:2,539)=ftimes(1:2,537)+ftimes(1:2,538)
 nflops(539)=nflops(537)+nflops(538)


!Gather the different parts of fourdp
 ncount(9)=sum(ncount(260:265))
 times(1,9)=sum(times(1,260:265))
 times(2,9)=sum(times(2,260:265))
 ftimes(1,9)=sum(ftimes(1,260:265))
 ftimes(2,9)=sum(ftimes(2,260:265))
 nflops(9)=sum(nflops(260:265))


!Gather the different parts of getghc
 ncount(200)=sum(ncount(201:206))
 times(1,200)=sum(times(1,201:206))
 times(2,200)=sum(times(2,201:206))
 ftimes(1,200)=sum(ftimes(1,201:206))
 ftimes(2,200)=sum(ftimes(2,201:206))
 nflops(200)=sum(nflops(201:206))

!Gather the different parts of projbd
 ncount(210)=sum(ncount(211:213))
 times(1,210)=sum(times(1,211:213))
 times(2,210)=sum(times(2,211:213))
 ftimes(1,210)=sum(ftimes(1,211:213))
 ftimes(2,210)=sum(ftimes(2,211:213))
 nflops(210)=sum(nflops(211:213))

!Gather the different parts of rwwf (wavefunctions read/write)
 ncount(14)=sum(ncount(270:288))
 times(1,14)=sum(times(1,270:288))
 times(2,14)=sum(times(2,270:288))
 ftimes(1,14)=sum(ftimes(1,270:288))
 ftimes(2,14)=sum(ftimes(2,270:288))
 nflops(14)=sum(nflops(270:288))


!Compute xc part of rhohxc and mkvxc3,
!minus the calls to fourdp inside that part
 ncount(11)=ncount(81)+ncount(181)
 times(1:2,11)=times(1:2,81)+times(1:2,181)-times(1:2,82)
 ftimes(1:2,11)=ftimes(1:2,81)+ftimes(1:2,181)-ftimes(1:2,82)
 nflops(11)=nflops(81)+nflops(181)-nflops(82)


!Estimate the complement of getghc (non fourwf, non nonlop)
 ncount(98)=-1
 times(1:2,98)=times(1:2,200)-times(1:2,241)-times(1:2,221)
 ftimes(1:2,98)=ftimes(1:2,200)-ftimes(1:2,241)-ftimes(1:2,221)
 nflops(98)=nflops(200)-nflops(241)-nflops(221)


!Estimate the complement of cgwf (non getghc,projbd)
 ncount(40)=-1
 times(1:2,40)=times(1:2,22)+times(1:2,530)-times(1:2,201)-times(1:2,211)
 ftimes(1:2,40)=ftimes(1:2,22)+ftimes(1:2,530)-ftimes(1:2,201)-ftimes(1:2,211)
 nflops(40)=nflops(22)+nflops(530)-nflops(201)-nflops(211)

!Estimate the complement of cgwf3 (non getghc,projbd,nonlop,fourwf)
 ncount(140)=-1
 times(1:2,140)=times(1:2,122)-times(1:2,202)-times(1:2,212) &
& -times(1:2,227)-times(1:2,228)-times(1:2,244)
 ftimes(1:2,140)=ftimes(1:2,122)-ftimes(1:2,202)-ftimes(1:2,212) &
& -ftimes(1:2,227)-ftimes(1:2,228)-ftimes(1:2,244)
 nflops(140)=nflops(122)-nflops(202)-nflops(212)-nflops(227)-nflops(228)-nflops(244)

!Estimate different complements in vtowfk
!vtowfk(ssdiag) (= vtowfk(loop)    - cgwf )
 ncount(7)=-1
 times(1:2,7)=times(1:2,39)-times(1:2,22)-times(1:2,530)
 ftimes(1:2,7)=ftimes(1:2,39)-ftimes(1:2,22)-ftimes(1:2,530)
 nflops(7)=nflops(39)-nflops(22)-nflops(530)

!vtowfk(contrib) (= vtowfk (3) - nonlop%vtowfk - fourwf%vtowfk )
 ncount(8)=ncount(30)
 times(1:2,8)=times(1:2,30)-times(1:2,222)-times(1:2,242)
 ftimes(1:2,8)=ftimes(1:2,30)-ftimes(1:2,222)-ftimes(1:2,242)
 nflops(8)=nflops(30)-nflops(222)-nflops(242)
!vtowfk (1) = vtowfk - vtowfk(loop) - vtowfk(3)
 ncount(31)=ncount(28)
 times(1:2,31)=times(1:2,28)-times(1:2,39)-times(1:2,30)
 ftimes(1:2,31)=ftimes(1:2,28)-ftimes(1:2,39)-ftimes(1:2,30)
 nflops(31)=nflops(28)-nflops(39)-nflops(30)

!Estimate different complements in vtowfk3
!vtowfk3(contrib) (= vtowfk3(loop) - cgwf - fourwf%vtowfk3 - rwwf%vtowfk3
!- projbd%vtowfk3 )
 ncount(108)=-1
 times(1:2,108)=times(1:2,139)-times(1:2,122)-times(1:2,245)-times(1:2,288)&
& -times(1:2,213)
 ftimes(1:2,108)=ftimes(1:2,139)-ftimes(1:2,122)-ftimes(1:2,245)-ftimes(1:2,288)&
& -ftimes(1:2,213)
 nflops(108)=nflops(139)-nflops(122)-nflops(245)-nflops(288)&
& -nflops(213)

!vtowfk (1) = vtowfk3 - vtowfk3(loop) - vtowfk3 (3)
 ncount(131)=ncount(128)
 times(1:2,131)=times(1:2,128)-times(1:2,139)-times(1:2,130)
 ftimes(1:2,131)=ftimes(1:2,128)-ftimes(1:2,139)-ftimes(1:2,130)
 nflops(131)=nflops(128)-nflops(139)-nflops(130)


!Estimate different complements in vtorho
!vtorho (1) (= vtorho (1,2) - vtorho(2) - vtorho:synchro )
 ncount(18)=ncount(21)
 times(1:2,18)=times(1:2,24)-times(1:2,25)-times(1:2,66)
 ftimes(1:2,18)=ftimes(1:2,24)-ftimes(1:2,25)-ftimes(1:2,66)
 nflops(18)=nflops(24)-nflops(25)-nflops(66)

!vtorho-kpt loop (= vtowfk (2) - vtowfk - rwwf)
 ncount(26)=ncount(25)
 times(1:2,26)=times(1:2,25)-times(1:2,28)-times(1:2,271)
 ftimes(1:2,26)=ftimes(1:2,25)-ftimes(1:2,28)-ftimes(1:2,271)
 nflops(26)=nflops(25)-nflops(28)-nflops(271)

!vtorho(4)-mkrho (= vtorho(4) - mkrho%vtorho- PAW)
 ncount(99)=ncount(27)
 times(1:2,99)=times(1:2,27)-times(1:2,292)
 ftimes(1:2,99)=ftimes(1:2,27)-ftimes(1:2,292)
 nflops(99)=nflops(27)-nflops(292)

!Estimate different complements in vtorho3
!vtorho3 (1) (= vtorho3 (1,2) - vtorho3(2) - vtorho3:synchro )
 ncount(118)=ncount(121)
 times(1:2,118)=times(1:2,124)-times(1:2,125)-times(1:2,166)
 ftimes(1:2,118)=ftimes(1:2,124)-ftimes(1:2,125)-ftimes(1:2,166)
 nflops(118)=nflops(124)-nflops(125)-nflops(166)
!vtorho3-kpt loop (= vtowfk3 (2) - vtowfk3 - rwwf)
 ncount(126)=ncount(125)
 times(1:2,126)=times(1:2,125)-times(1:2,128)-times(1:2,287)
 ftimes(1:2,126)=ftimes(1:2,125)-ftimes(1:2,128)-ftimes(1:2,287)
 nflops(126)=nflops(125)-nflops(128)-nflops(287)

!Estimate complement in mkrho
!mkrho/= (= mkrho/=+fourwf - fourwf%mkrho )
 ncount(298)=ncount(299)
 times(1:2,298)=times(1:2,299)-times(1:2,243)
 ftimes(1:2,298)=ftimes(1:2,299)-ftimes(1:2,243)
 nflops(298)=nflops(299)-nflops(243)


!Estimate complement in loper3
!loper3(other) (= loper3 - loper3(kpgio) - loper3(getmpw) - loper3(inwffil)
!scfcv3 - loper3(outwf) )
 ncount(145)=ncount(141)
 times(1:2,145)=times(1:2,141)-times(1:2,142)-times(1:2,143)-times(1:2,144)&
& -times(1:2,120)-times(1:2,146)
 ftimes(1:2,145)=ftimes(1:2,141)-ftimes(1:2,142)-ftimes(1:2,143)-ftimes(1:2,144)&
& -ftimes(1:2,120)-ftimes(1:2,146)
 nflops(145)=nflops(141)-nflops(142)-nflops(143)-nflops(144)&
& -nflops(120)-nflops(146)

!Estimate complement in respfn
!respfn(init 1st) (= respfn(1st) - loper3 )
 ncount(138)=ncount(137)
 times(1:2,138)=times(1:2,137)-times(1:2,141)
 ftimes(1:2,138)=ftimes(1:2,137)-ftimes(1:2,141)
 nflops(138)=nflops(137)-nflops(141)


!Estimate complement in screening
!screening/=fourdp = screening - fourdp%rhotwg%ch - fourdp%ckxcldag - fourdp%fftwfn%ch
 ncount(320)=ncount(301)
 times(1:2,320)=times(1:2,301)-times(1:2,261)-times(1:2,263)-times(1:2,264)
 ftimes(1:2,320)=ftimes(1:2,301)-ftimes(1:2,261)-ftimes(1:2,263)-ftimes(1:2,264)
 nflops(320)=nflops(301)-nflops(261)-nflops(263)-nflops(264)


!Estimate complement in sigma
!sigma/=fourdp = sigma - fourdp%rhotwg%si - fourdp%fftwfn%si
 ncount(410)=ncount(401)
 times(1:2,410)=times(1:2,401)-times(1:2,262)-times(1:2,265)
 ftimes(1:2,410)=ftimes(1:2,401)-ftimes(1:2,262)-ftimes(1:2,265)
 nflops(410)=nflops(401)-nflops(262)-nflops(265)

!Calculating Gigaflops for all case
 do itim=1,mtim
  mflops(itim)=-2
  if(ftimes(1,itim)/= 0) then
   mflops(itim)=nflops(itim)*1.e-9/ftimes(1,itim)
  else
   mflops(itim)=-1
  end if
 end do

!Warning if the time is negative
 do itim=1,mtim
  if(times(1,itim)<-tol6 .or. times(2,itim)<-tol6 .or. &
&  ncount(itim)<-1 )then
   write(message, '(a,a,a,a,i4,a,a,a,a,es16.6,a,es16.6,a,i6,a,es16.6)' ) ch10,&
&   ' timana : WARNING -',ch10,&
&   '  Timing section #',itim,', name :  ',names(itim),ch10,&
&   '  CPU =',times(1,itim),', Wall=',times(2,itim),', &
&   ncount=',ncount(itim), &
&   ' flops=',nflops(itim)
   call wrtout(6,message,'PERS')
  end if
 end do


!DEBUG
!write(6,*)' timana : after get time reports  '
!if(.true.)stop
!ENDDEBUG

!Major independent code sections
 nlist=70
 allocate(list(nlist))
 list=(/ 2  ,3  ,5  ,6  ,7  ,  8  ,9  ,11 ,12 ,14,&
& 15 ,17 ,18 ,19 ,23 ,  26 ,29 ,31 ,37 ,38,&
& 42 ,50 ,52 ,66 ,69 ,  73 ,75 ,76 ,77 ,83,&
& 85 ,86 ,90 ,93 ,94 ,  96 ,98 ,99 ,108,118,&
& 126,129,130,131,140,  154,161,166,191,210,&
& 236,237,298,320,410,  521,522,553,555,556,&
& 557,559,560,561,562,  563,564,565,566,567  /)

 percent_limit=0.5_dp
 if(timopt<0)percent_limit=-0.1_dp

!In case there is parallelism
 if(me==0 .and. nproc>1 ) then

! First, report times for node 0

! Find normalization to report timing as % total time
  cpunm=100._dp/tsec(1)
  wallnm=100._dp/tsec(2)

! (0) Take care of major independent code sections
! for this account of node 0 timing

  if ( papiopt == 0 ) then
   write(ab_out,  '(a,a,a,a,/,a,a,a)' ) '-',ch10,&
&   '- For major independent code sections,', &
&   ' cpu and wall times (sec),',&
&   '-  as well as % of the time and number of calls for node 0',&
&   '-'
   write(ab_out, '(a,t26,a,t34,a,t42,a,t51,a,t57,a)' )&
&   '- routine','cpu','%','wall','%',' number of calls '


  else if  ( papiopt == 1 ) then
   write(ab_out,  '(a,a,a,a/,a/,a/a)' ) '-',ch10,&
&   '- For major independent code sections,', &
&   ' cpu ( given by papy) and wall times (sec),',&
&   '-  as well as % of the time and number of calls for node 0',&
&   '- by decreasing CPU', &
&   '-'

   write(ab_out, '(a,t26,a,t35,a,t43,a,t51,a,t63,a)' )&
&   '- routine','wall','%', 'cpu', 'Gigaflops', ' number of calls '
  end if



! Sort the list by decreasing CPU time
  do ii=1,nlist
   do ilist=1,nlist-1
    if(times(1,list(ilist))<times(1,list(ilist+1)))then
     temp_list=list(ilist)
     list(ilist)=list(ilist+1)
     list(ilist+1)=temp_list
    end if
   end do
  end do

  subcpu=0.0_dp
  subwal=0.0_dp
  other_cpu=0.0_dp
  other_wal=0.0_dp
  nothers=0

  do ilist=1,nlist
   if( (times(1,list(ilist))*cpunm > percent_limit .or.       &
&   times(2,list(ilist))*wallnm> percent_limit     ).and. &
&   ncount(list(ilist))/=0                                 )then
!   Analyse des temps

    if (  papiopt == 0 ) then

     write(ab_out,format01040)names(list(ilist)),&
&     times(1,list(ilist)),times(1,list(ilist))*cpunm,&
&     times(2,list(ilist)),times(2,list(ilist))*wallnm,&
&     ncount(list(ilist))
    else if (  papiopt == 1 ) then
!    if ( mflops(list(ilist)) > 0.1 ) then
     write(ab_out,format01041)names(list(ilist)),&
&     times(2,list(ilist)),times(2,list(ilist))*wallnm,&
&     ftimes(1,list(ilist)), &
&     mflops(list(ilist)), &
&     ncount(list(ilist))
!    endif
    end if
   else
    nothers=nothers+1
    other_cpu=other_cpu+times(1,list(ilist))
    other_wal=other_wal+times(2,list(ilist))
   end if
   subcpu=subcpu+times(1,list(ilist))
   subwal=subwal+times(2,list(ilist))
  end do

  if  (  papiopt == 0 ) then
   write(ab_out,format01045)nothers,' others  ',&
&   other_cpu,other_cpu*cpunm,&
&   other_wal,other_wal*wallnm
   write(ab_out,'(a)' ) '-'
   write(ab_out,01200) subcpu,subcpu*cpunm,subwal,subwal*wallnm
  else if  (  papiopt == 1 ) then
   write(ab_out,format01045)nothers,' others  ',&
&   other_wal,other_wal*wallnm
   write(ab_out,'(a)' ) '-'
   write(ab_out,01200) subwal,subwal*wallnm
  end if


! writing the result by decreasing execution speed
! first sort the list by decreasing execution speed
  if (papiopt == 1) then
   do ii=1,nlist
    do ilist=1,nlist-1
     if(mflops(list(ilist))<mflops(list(ilist+1)))then
      temp_list=list(ilist)
      list(ilist)=list(ilist+1)
      list(ilist+1)=temp_list
     end if
    end do
   end do
   write(ab_out,  '(a,a,a,a,/,a/,a,a,a)' ) '-',ch10,&
&   '- For major independent code sections,', &
&   ' cpu ( given by papy) and wall times (sec),',&
&   '-  as well as % of the time and number of calls for node 0', &
&   '- by significant decreasing speed ( > 100 MegaFlops )'


   write(ab_out, '(a,t26,a,t35,a,t43,a,t51,a,t63,a)' )&
&   '- routine','wall','%', 'cpu', 'Gigaflops', ' number of calls '
   do ilist=1,nlist
    if ( mflops(list(ilist)) > 0.1 ) then
     write(ab_out,format01041)names(list(ilist)),&
&     times(2,list(ilist)),times(2,list(ilist))*wallnm,&
&     ftimes(1,list(ilist)), &
&     mflops(list(ilist)), &
&     ncount(list(ilist))
    end if
   end do
  end if
 end if

!Now, gather all information
 call xsum_mpi(times,spaceworld,ierr)
 call xsum_mpi(ncount,spaceworld,ierr)
 call xsum_mpi(ftimes,spaceworld,ierr)
 call xsum_mpi(nflops,spaceworld,ierr)

!Only the world master writes
 if(me==0) then

! Find normalization to report timing as % total time
  cpunm=100._dp/tsec(1)
  wallnm=100._dp/tsec(2)


! Calculating Gigaflops for all process
  do itim=1,mtim
   mflops(itim)=-2
   if(ftimes(1,itim)/= 0) then
    mflops(itim)=nflops(itim)*1.e-9/ftimes(1,itim)
   else
    mflops(itim)=-1
   end if
  end do

! _______________________________________

! Write timing output for cpu times

! (1) Take care of major independent code sections
  if ( papiopt == 0 ) then
   write(ab_out,'(/,a,/,a,/)' )&
&   '- For major independent code sections, cpu and wall times (sec),',&
&   '- as well as % of the total time and number of calls '
   write(ab_out,'(a,t27,a,t35,a,t43,a,t52,a,t58,a)')&
&   '- routine','cpu','%','wall','%', ' number of calls '
   write(ab_out,'(a,t27,a,t35,a,t43,a,t52,a,t58,a,t70,a)')&
&   '-        ','   ',' ','    ',' ','  (-1=no count)'
  else if (  papiopt == 1 ) then
   write(ab_out,'(/,a,/,a,/,a/)' )&
&   '- For major independent code sections, cpu given by papi and wall times (sec),',&
&   '- as well as % of the total time and number of calls ', &
&   '- by decreasing CPU '

   write(ab_out, '(a,t26,a,t35,a,t43,a,t50,a,t63,a)' )&
&   '- routine','wall','%', 'cpu', 'Gigaflops', ' number of calls '

   write(ab_out,'(a,t26,a,t35,a,t43,a,t50,a,t63,a)')&
&   '-        ','   ',' ','    ','(-1=no count) ','  (-1=no count)'

  end if

! Sort the list by decreasing CPU time
  do ii=1,nlist
   do ilist=1,nlist-1
    if(times(1,list(ilist))<times(1,list(ilist+1)))then
     temp_list=list(ilist)
     list(ilist)=list(ilist+1)
     list(ilist+1)=temp_list
    end if
   end do
  end do

  subcpu=0.0_dp
  subwal=0.0_dp
  other_cpu=0.0_dp
  other_wal=0.0_dp
  nothers=0

  do ilist=1,nlist
   if( (times(1,list(ilist))*cpunm > percent_limit .or.  &
&   times(2,list(ilist))*wallnm> percent_limit     ) .and. &
&   ncount(list(ilist))/=0                            )then
    if ( papiopt == 0 ) then

     write(ab_out,format01040)names(list(ilist)),&
&     times(1,list(ilist)),times(1,list(ilist))*cpunm,&
&     times(2,list(ilist)),times(2,list(ilist))*wallnm,&
&     ncount(list(ilist))
    else if ( papiopt == 1 ) then
!    Analyse des temps + Gigaflops
!    if ( mflops(list(ilist)) > 0.1 ) then
     write(ab_out,format01041)names(list(ilist)),&
&     times(2,list(ilist)),times(2,list(ilist))*wallnm,&
&     ftimes(1,list(ilist)),  &
&     mflops(list(ilist)),  &
&     ncount(list(ilist))
!    endif
    end if
   else
    nothers=nothers+1
    other_cpu=other_cpu+times(1,list(ilist))
    other_wal=other_wal+times(2,list(ilist))
   end if
   subcpu=subcpu+times(1,list(ilist))
   subwal=subwal+times(2,list(ilist))
  end do

  if ( papiopt == 0 ) then

   write(ab_out,format01045)nothers,' others  ',&
&   other_cpu,other_cpu*cpunm,other_wal,other_wal*wallnm
   write(ab_out, '(/,a,f15.3,f6.1,f11.3,f6.1)' )&
&   '- subtotal     ', subcpu,subcpu*cpunm,subwal,subwal*wallnm
  else if ( papiopt == 1 ) then

   write(ab_out,format01045)nothers,' others  ',&
&   other_wal,other_wal*wallnm
   write(ab_out, '(/,a,f15.3,f6.1,f11.3,f6.1)' )&
&   '- subtotal     ', subwal,subwal*wallnm
  end if

! writing the result by decreasing execution speed
! first sort the list by decreasing execution speed
  if (papiopt == 1) then
   do ii=1,nlist
    do ilist=1,nlist-1
     if(mflops(list(ilist))<mflops(list(ilist+1)))then
      temp_list=list(ilist)
      list(ilist)=list(ilist+1)
      list(ilist+1)=temp_list
     end if
    end do
   end do
   write(ab_out,  '(a,a,a,a,/,a/)' ) '-',ch10,&
&   '- For major independent code sections,', &
&   ' cpu ( given by papy) and wall times (sec),',&
&   '-  as well as % of the time and number of calls ' ,&
&   '- by significant decreasing speed ( > 100 Megaflops)'


   write(ab_out, '(a,t26,a,t35,a,t43,a,t51,a,t63,a)' )&
&   '- routine','wall','%', 'cpu', 'Gigaflops', ' number of calls '

   do ilist=1,nlist
    if (  mflops(list(ilist)) > 0.1 ) then
     write(ab_out,format01041)names(list(ilist)),&
&     times(2,list(ilist)),times(2,list(ilist))*wallnm,&
&     ftimes(1,list(ilist)), &
&     mflops(list(ilist)), &
&     ncount(list(ilist))
    end if
   end do
  end if


  deallocate(list)


! (2) Linked chains

  if(timopt<0)then

   nlink=27
   do ilink=1,nlink

!   TO BE SUPPRESSED FOR v4.5 DELIVERY
!   if(ilink==17)cycle
!   END REINSTALL


    if(ilink==1)then

     nlist=9;allsub=0
     allocate(list(nlist))
     list=(/1,41,42,43,44,45,100,46,50/)
     if(ncount(list(1))/=0)write(ab_out,'(/,a)') ' Partitioning of abinit '

    else if(ilink==2)then

     nlist=7;allsub=0
     allocate(list(nlist))
     list=(/100,32,132,84,301,401,501/)
     if(ncount(list(1))/=0)write(ab_out,'(/,a)') ' Partitioning of driver '

    else if(ilink==3)then

     nlist=7;allsub=0
     allocate(list(nlist))
     list=(/32,33,15,34,35,36,562/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of gstate '

    else if(ilink==4)then

     nlist=19;allsub=0
     allocate(list(nlist))
     list=(/20,54,91,92,55,56,21,80,52,58,59,60,61,233,37,558,559,560,563/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of scfcv '

    else if(ilink==5)then

     nlist=4
     allocate(list(nlist))
     list=(/80,57,53,69/);allsub=sum(ncount(list(:)))
     if(allsub>0)write(ab_out, '(/,a)' ) ' Partitioning of etotfor/rhotov '

    else if(ilink==6)then

!    Original coding
     nlist=14;allsub=0
     allocate(list(nlist))
     list=(/21,18,28,26,66,29,27,232,70,95,271,555,556,557/)
!    DEBUG : Examination of suscep ; to be updated
!    nlist=12
!    allocate(list(nlist))
!    list=(/21,24,25,28,26,66,29,27,87,88,89,90/)
!    ENDDEBUG
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of vtorho '

    else if(ilink==7)then

     nlist=10;allsub=0
     allocate(list(nlist))
     list=(/28,530,22,585,583,222,572,242,537,586/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of vtowfk '

    else if(ilink==8)then

     if(abs(timopt)==3)then
      nlist=12;allsub=0
      allocate(list(nlist))
      list=(/530,204,205,531,571,532,533,534,535,536,584,587/)
      if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of lobpcg '
     else 
      nlist=3;allsub=0
      allocate(list(nlist))
      list=(/530,204,205/)
      if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of lobpcg (light analysis: for a deeper one, use timopt=3 or -3)'
     endif

    else if(ilink==9)then

     nlist=4;allsub=0
     allocate(list(nlist))
     list=(/22,201,40,211/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of cgwf '

    else if(ilink==10)then

     nlist=7;allsub=0
     allocate(list(nlist))
     list=(/132,133,134,135,136,138,141/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of respfn '

    else if(ilink==11)then

     nlist=7;allsub=0
     allocate(list(nlist))
     list=(/141,142,143,144,120,145,146/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of loper3 '

    else if(ilink==12)then

     nlist=9;allsub=0
     allocate(list(nlist))
     list=(/120,154,121,157,152,158,160,147,564/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of scfcv3 '

    else if(ilink==13)then

     nlist=11;allsub=0
     allocate(list(nlist))
     list=(/121,118,128,126,287,166,129,127,565,566,567/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of vtorho3 '

    else if(ilink==14)then

     nlist=8;allsub=0
     allocate(list(nlist))
     list=(/128,131,122,245,288,213,108,130/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of vtowfk3 '

    else if(ilink==15)then

     nlist=8;allsub=0
     allocate(list(nlist))
     list=(/122,140,202,207,212,227,228,244/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of cgwf3 '

    else if(ilink==16)then

     nlist=4;allsub=0
     allocate(list(nlist))
     list=(/200,241,221,98/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of getghc '

    else if(ilink==17)then

     nlist=19;allsub=0
     allocate(list(nlist))
     list=(/68,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of fourwf '

    else if(ilink==18)then

     nlist=7;allsub=0
     allocate(list(nlist))
     list=(/233,234,235,255,229,236,237/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of outkss '

    else if(ilink==19)then

     nlist=10;allsub=0
     allocate(list(nlist))
     list=(/301,302,303,304,305,306,307,308,309,310/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of screening '

    else if(ilink==20)then

     nlist=5;allsub=0
     allocate(list(nlist))
     list=(/401,402,403,404,405/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of sigma  '

    else if(ilink==21)then

     nlist=6;allsub=0
     allocate(list(nlist))
     list=(/421,422,423,424,425,426/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of csigme '

    else if(ilink==22)then

     nlist=18
     allocate(list(nlist))
     list=(/550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567/)
     allsub=sum(ncount(list(:)))
     if(allsub>0)write(ab_out, '(/,a)' ) ' Partitioning of PAW '

    else if(ilink==23)then

     nlist=4;allsub=0
     allocate(list(nlist))
     list=(/534,206,545,546/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of prep_getghc '

    else if(ilink==24)then

     nlist=4;allsub=0
     allocate(list(nlist))
     list=(/539,256,547,548/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of prep_fourwf '

    else if(ilink==25)then

     nlist=4;allsub=0
     allocate(list(nlist))
     list=(/570,580,581,582/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of prep_nonlop '

    else if(ilink==26)then

!    nlist=3;allsub=0
!    allocate(list(nlist))
!    list=(/540,541,542/)
!    if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of sg_fourwf '
     nlist=0;allsub=0
     write(ab_out, '(/,a,/,/)' ) ' Partitioning of sg_fourwf (disabled) '

    else if(ilink==27)then

     nlist=3;allsub=0
     allocate(list(nlist))
     list=(/292,549,538/)
     if(ncount(list(1))/=0)write(ab_out, '(/,a)' ) ' Partitioning of mkrho '

    end if

    if(nlist>0)then
     if(ncount(list(1))/=0.or.allsub>0)then
      subcpu=0.0_dp
      subwal=0.0_dp
      do ilist=1,nlist
       if(ncount(list(ilist))/=0)then
        if (papiopt == 0) then
         write(ab_out,format01040)names(list(ilist)),&
&         times(1,list(ilist)),times(1,list(ilist))*cpunm,&
&         times(2,list(ilist)),times(2,list(ilist))*wallnm,&
&         ncount(list(ilist))
        else if  (papiopt == 1) then
         if ( mflops(list(ilist)) > 0.1) then
          write(ab_out,format01041)names(list(ilist)),&
&          times(2,list(ilist)),times(2,list(ilist))*wallnm,&
&          ftimes(1,list(ilist)), &
&          mflops(list(ilist)),ncount(list(ilist))
         end if
        end if
        if(ilist/=1.or.allsub>0)then
         subcpu=subcpu+times(1,list(ilist))
         subwal=subwal+times(2,list(ilist))
        else
         write(ab_out, '(a)' ) ' '
        end if
       end if
      end do

      if ( papiopt == 0) then
       write(ab_out, '(/,a,f15.3,f6.1,f11.3,f6.1)' )&
&       '- subtotal     ', subcpu,subcpu*cpunm,subwal,subwal*wallnm
      else if ( papiopt == 1) then
       write(ab_out, '(/,a,f15.3,f6.1,f11.3,f6.1)' )&
&       '- subtotal     ', subwal,subwal*wallnm
      end if
 
     end if
     deallocate(list)
    endif

!   End of loop on linked chains
   end do

!  For parallel case
   if(mpi_enreg%paral_compil==1)then
    write(ab_out, '(a,/,a)' )'-',&
&    '-Synchronisation (=leave_test) and MPI calls '
    nlist=14
    allocate(list(nlist))
    list=(/48,29,61,62,63,64,65,66,67,71,85,86,543,544/)
    subcpu=0.0_dp
    subwal=0.0_dp
    if(ncount(list(1))/=0)then
     do ilist=1,nlist
      if(ncount(list(ilist))/=0)then

       if (papiopt == 0) then
        write(ab_out,format01040)names(list(ilist)),&
&        times(1,list(ilist)),times(1,list(ilist))*cpunm,&
&        times(2,list(ilist)),times(2,list(ilist))*wallnm,&
&        ncount(list(ilist))
       else if  (papiopt == 1) then
        if (mflops(list(ilist)) > 0.1 ) then
         write(ab_out,format01041)names(list(ilist)),&
&         times(2,list(ilist)),times(2,list(ilist))*wallnm,&
&         ftimes(1,list(ilist)), &
&         mflops(list(ilist)),ncount(list(ilist))
        end if
       end if !papiopt

       if(ilist/=1)then
        subcpu=subcpu+times(1,list(ilist))
        subwal=subwal+times(2,list(ilist))
       else
        write(ab_out, '(a)' ) '-'
       end if
      end if !ncount
     end do !ilist

     if (papiopt == 0) then
      write(ab_out, '(a,/,a,f15.3,f6.1,f11.3,f6.1)' )'-',&
&      '- subtotal     ', subcpu,subcpu*cpunm,subwal,subwal*wallnm
     else if (papiopt == 1) then
      write(ab_out, '(a,/,a,f15.3,f6.1,f11.3,f6.1)' )'-',&
&      '- subtotal     ', subwal,subwal*wallnm
     end if !papiopt

    end if !ncount
    deallocate(list)
   end if !mpi_enreg%paral_compil

   write(ab_out, '(/,a)' ) ' Additional information '
   nlist=19
   allocate(list(nlist))
   list=(/47,49,51,68,72,73,74,77,78,79,97,82,87,88,93,94,331,332,333/)
   do ilist=1,nlist
    if(ncount(list(ilist))/=0)then
     write(ab_out,format01040)names(list(ilist)),&
&     times(1,list(ilist)),times(1,list(ilist))*cpunm,&
&     times(2,list(ilist)),times(2,list(ilist))*wallnm,&
&     ncount(list(ilist))
    end if
   end do
   deallocate(list)

!  The detailed analysis cannot be done in the multidataset mode
   if(ndtset<2)then
    write(ab_out, '(/,/,a,/,a,/,a)' ) &
&    ' Detailed analysis of some time consuming routines ',&
&    '                          tcpu    ncalls  tcpu/ncalls    ndata tcpu/ncalls/ndata',&
&    '                         (sec)                (msec)              (microsec)'
    nlist=8
    allocate(list(nlist))
    list=(/2,3,9,75,76,77,210,11/)
    do ilist=1,nlist
     if(ncount(list(ilist))/=0)then
      write(ab_out, '(a,a16,f12.3,i10,f12.3,i10,f12.3)' )'- ',names(list(ilist)),&
&      times(1,list(ilist)),ncount(list(ilist)),&
&      1000.0_dp*times(1,list(ilist))/dble(ncount(list(ilist))),&
&      ndata(list(ilist)),&
&      1000000.0_dp*times(1,list(ilist))/dble(ncount(list(ilist))&
&      *dble(ndata(list(ilist))))
     else
      write(ab_out, '(a,a16,f12.3,i10)' )'- ',names(list(ilist)),&
&      times(1,list(ilist)),ncount(list(ilist))
     end if
    end do !ilist
    deallocate(list)
   else
    write(ab_out, '(/,a)' ) &
&    ' timana : in multi dataset mode, the more detailed analysis is not done.'
    write(6, '(/,a)' ) &
&    ' timana : in multi dataset mode, the more detailed analysis is not done.'
   end if !ndtset

!  End the condition of timopt<0
  end if

 end if ! me==0

end subroutine timana

!The number of fourwf and nonlop calls can be computed as follows, in the
!groud-state case, with no reading of wavefunctions (irdwfk==0 and the like),
!and iscf>0 :
!
!1) For fourwf.f
!
!In each cgwf call, there will be
!1 call (isign=+1 and -1) for the first gradient calculation,
!and iline calls for the line minimizations,
!minus the number of ffts skipped because some wfs are sufficiently converged
!(there is a counter for that, see the log file)
!
!There are nband*nkpt*(nstep+2) calls to cgwf presently, where the
!(nstep+2) comes from the number of the presence of 2 nonscf loops
!in the first 2 steps.
!Thus, the number of fourwf calls in cgwf is
!nband*nkpt*(nstep+2)*(1+iline) - nskip_fourwf_in_cgwf
!
!To compute the density (either in vtowfk or in vtorho - by a mkrho call - )
!at each step, there will be nband*nkpt one-way calls,
!minus the number of bands skipped because the occupation number
!is too small (smaller than 1.0d-14). There is another counter for that.
!Thus, the number of fourwf calls for the density is
!nband*nkpt*nstep - nskip_fourwf_for_density
!
!For example, for Si with nline=3, nkpt=2, nband=4, nstep=10, and supposing
!no fourwf calls are skipped, there will be
!at most 4*2*12=96 calls to cgwf, with 4 two-way fft,
!that is 384 two-way ffts,
!and 4*2*10=80 one-way ffts to make the density.
!Altogether 464-nskip one-way ffts at most.
!
!2) For nonlop.f
!
!Presently, there are three different types of call to nonlop :
!for energy and gradient wrt wavefunctions (choice=1), for forces (choice=2),
!and for stresses (choice=3).
!
!In each cgwf call, there will be one nonlop call for two fourwf calls
!(independently of the number of skipped fourwf calls, since
!nonlop is also skipped then). These are the only calls with choice=1.
!Thus the number will be
!nband*nkpt*(nstep+2)*(1+iline) - nskip_fourwf_in_cgwf
!
!The number of choice=2 nonlop calls is equal to the number of fourwf calls
!to make the density, that is
!nband*nkpt*nstep - nskip_fourwf_for_density
!
!The number of choice=8 calls is equal to the number of occupied bands
!at the end of the calculation :
!nband(occupied)*nkpt
!The number of bands skipped then is not counted.
!
!NOTE : the number of fourwf calls is equal to
!the # of nonlop (choice=1) calls + the # of nonlop (choice=2) calls
!!***
