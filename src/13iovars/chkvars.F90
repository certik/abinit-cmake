!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkvars
!! NAME
!! chkvars
!!
!! FUNCTION
!! Examine the input string, to check whether all names are allowed
!!
!! COPYRIGHT
!! Copyright (C) 2007-2008 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  string*(*)=string of character
!!   the string (with upper case) from the input file, to which the CML data are appended
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine chkvars (string)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: string

!Local variables-------------------------------
 character :: blank=' '
!scalars
 integer :: ii,index_blank,index_current,index_endword,problem
 character(len=100) :: list_logicals,list_strings
 character(len=10000) :: list_vars
 character(len=20) :: word
 character(len=500) :: message

!************************************************************************

!DEBUG
!write(6,*)' chkvars : enter '
!write(6,*)trim(string)
!stop
!ENDDEBUG

!Here, list all admitted variable names (max 10 per line, to fix the ideas)
!A
 list_vars=                 ' accesswff acell algalch alpha amu angdeg atvshift awtr' 
!B
 list_vars=trim(list_vars)//' bandpp bdberry bdgw berryopt bmass boxcenter boxcutmin brvltt bxctmindg'
!C
 list_vars=trim(list_vars)//' ceksph charge chkexit chkprim cpus cpum cpuh'
!D
 list_vars=trim(list_vars)//' dedlnn delayperm densty diecut diegap dielam dielng diemac'
 list_vars=trim(list_vars)//' diemix dilatmx dmatpuopt dmatudiag dmatpawu dosdeltae dsifkpt dtion'
!E
 list_vars=trim(list_vars)//' ecut ecuteps ecutsigx ecutsm ecutwfn effmass efield'
 list_vars=trim(list_vars)//' enunit eshift etsfgroups etsfmain exchmix exchn2n3d'
!F
 list_vars=trim(list_vars)//' fband fftalg fftcache fftgw fft_opt_lob fixmom'
 list_vars=trim(list_vars)//' freqremax freqspmax freqsusin freqsuslo friction frzfermi'
!G
 list_vars=trim(list_vars)//' genafm getacfd getcell getsuscep getddk getden getkss getocc getqps getscr'
 list_vars=trim(list_vars)//' getvel getwfk getwfq getxcart getxred get1den get1wf'
 list_vars=trim(list_vars)//' gwcalctyp gwcomp gwencomp gwgamma gwmem gwpara'
   
!I
 list_vars=trim(list_vars)//' iatcon iatfix iatfixx iatfixy iatfixz iatsph'
 list_vars=trim(list_vars)//' iboxcut icoulomb icutcoul idyson ieig2rf ikhxc'
 list_vars=trim(list_vars)//' inclvkb intexact intexact intxc ionmov iprcch'
 list_vars=trim(list_vars)//' iprcel iprctfvw iprcfc irdsuscep irdddk irdqps'
 list_vars=trim(list_vars)//' irdkss irdscr irdwfk irdwfq ird1wf iscf'
 list_vars=trim(list_vars)//' isecur istatr istatshft istwfk ixc ixcpositron'
!J
 list_vars=trim(list_vars)//' jdtset jellslab jpawu'
!K
 list_vars=trim(list_vars)//' kberry kpara kpt kptbounds kptgw'
 list_vars=trim(list_vars)//' kptnrm kptopt kptrlatt kptrlen kssform'
!L
 list_vars=trim(list_vars)//' ldgapp lexexch localrdwf lofwrite lpawu ltypeorb'
!M
 list_vars=trim(list_vars)//' mdftemp mditemp mdwall mffmem mgriddg mixalch'
 list_vars=trim(list_vars)//' mkmem mkqmem mk1mem mqgrid mqgriddg'
!N
 list_vars=trim(list_vars)//' natcon natfix natfixx natfixy natfixz'
 list_vars=trim(list_vars)//' natom natpawu natrd natsph natvshift nband nbandkss'
 list_vars=trim(list_vars)//' nbandsus nbdblock nbdbuf nberry ncenter nconeq'
 list_vars=trim(list_vars)//' nctime ndivk ndivsm ndtset ndyson'
 list_vars=trim(list_vars)//' nfreqim nfreqre nfreqsp nfreqsus ngfft ngfftdg'
 list_vars=trim(list_vars)//' ngkpt ngroup_rf nkpt nkptgw nline nloalg nnos nnsclo'
 list_vars=trim(list_vars)//' nobj nomegasf nomegasrd norb noseinert npack npara npband'
 list_vars=trim(list_vars)//' npfft npkpt npsp npulayit npweps npwkss npwsigx'
 list_vars=trim(list_vars)//' npwwfn nqpt nqptdm nscforder sheps nshiftk'
 list_vars=trim(list_vars)//' nshsigx nshwfn nspden nspinor nsppol nstep nsym'
 list_vars=trim(list_vars)//' ntime ntypalch ntypat numorb nwfshist'
!O
 list_vars=trim(list_vars)//' objaat objbat objaax objbax objan objbn objarf'
 list_vars=trim(list_vars)//' objbrf objaro objbro objatr objbtr occ'
 list_vars=trim(list_vars)//' occopt omegasrdmax optcell optdriver optforces'
 list_vars=trim(list_vars)//' optfreqsus optnlxccc optstress ortalg outputxml'
!P
 list_vars=trim(list_vars)//' paral_kgb paral_rf parareel pawecutdg pawlcutd pawlmix pawmixdg'
 list_vars=trim(list_vars)//' pawnhatxc pawnphi pawntheta pawnzlm pawovlp pawoptmix'
 list_vars=trim(list_vars)//' pawprtdos pawprtvol pawsphmix pawspnorb pawstgylm pawusecp pawxcdev'
 list_vars=trim(list_vars)//' positron ppmfrq ppmodel prepanl prepgkk'
 list_vars=trim(list_vars)//' prtacfd prtbbb prtcml prtden papiopt'
 list_vars=trim(list_vars)//' prtdensph prtdos prtdosm prtefg prteig'
 list_vars=trim(list_vars)//' prtfc prtfsurf prtgeo prtgkk prtkpt'
 list_vars=trim(list_vars)//' prtnabla prtpot prtspcur prtstm prtvha prtvhxc'
 list_vars=trim(list_vars)//' prtvol prtvxc prtwant prtwf prt1dm pspso ptcharge'
!Q
 list_vars=trim(list_vars)//' qmass qprtrb qpt qptdm qptnrm quadmom'
!R
 list_vars=trim(list_vars)//' ratsph rcoord rcut rdmnb'
 list_vars=trim(list_vars)//' recefermi recnpath recnrec recptrott recrcut rectesteg rectolden'
 list_vars=trim(list_vars)//' restartxf rfasr rfatpol rfdir rfelfd rfmgfd rfmeth rfphon'
 list_vars=trim(list_vars)//' rfstrs rfthrd rfuser rf1atpol rf1dir rf1elfd'
 list_vars=trim(list_vars)//' rf1phon rf2atpol rf2dir rf2elfd rf2phon'
 list_vars=trim(list_vars)//' rf3atpol rf3dir rf3elfd rf3phon rhoqpmix rprim rtheta'
!S
 list_vars=trim(list_vars)//' sciss shiftk signperm slabwsrad slabzbeg slabzend soenergy so_psp '
 list_vars=trim(list_vars)//' spbroad spectral spgaxor spgorig spgroup spgroupma spinat splitsigc'
 list_vars=trim(list_vars)//' stmbias strfact strprecon strtarget supercell'
 list_vars=trim(list_vars)//' suskxcrs symafm symchi symrel symmorphi symsigma'
!T
 list_vars=trim(list_vars)//' td_maxene td_mexcit tfkinfunc tl_nprccg tl_radius'
 list_vars=trim(list_vars)//' tfnewton timopt '
 list_vars=trim(list_vars)//' tnons toldfe toldff tolmxf tolrff tolvrs'
 list_vars=trim(list_vars)//' tolwfr tphysel tsmear typat'
!U
 list_vars=trim(list_vars)//' udtset upawu usedmatpu useexexch usepawu usewvl'
 list_vars=trim(list_vars)//' useria userib useric userid userie'
 list_vars=trim(list_vars)//' userra userrb userrc userrd userre useylm'
!V
 list_vars=trim(list_vars)//' vaclst vacnum vacuum vacwidth vcutgeo vel vis vmass vprtrb'
!W
 list_vars=trim(list_vars)//' wfoptalg wtatcon wtk' 
 list_vars=trim(list_vars)//' wvl_cpmult wvl_crmult wvl_fpmult wvl_frmult wvl_hgrid wvl_nprccg'
 list_vars=trim(list_vars)//' w90iniprj w90lplot w90nplot w90prtunk'
!X
 list_vars=trim(list_vars)//' xangst xcart xred'
!Y
!Z
 list_vars=trim(list_vars)//' zcut znucl'

!Extra token, also admitted :
 list_vars=trim(list_vars)//' au Angstr Angstrom Angstroms Bohr Bohrs eV Ha Hartree Hartrees K Ry Rydberg Rydbergs sqrt'
!SIESTA strings
 list_vars=trim(list_vars)//' LatticeConstant '

!Logical input variables
 list_logicals=' SpinPolarized '

!String input variables
 list_strings=' cmlfile XCname '

!Transform to upper case
 call inupper(list_vars)
 call inupper(list_logicals)
 call inupper(list_strings)

!DEBUG
!write(6,*)' len_trim(string)=',len_trim(string)
!ENDDEBUG

 index_current=1
 do ! Infinite do-loop, to identify the presence of each potential variable names
! do ii=1,100 ! Infinite do-loop, to identify the presence of each potential variable names

  if(len_trim(string)<=index_current)exit
  index_blank=index(string(index_current:),blank)+index_current-1

  if(index('ABCDEFGHIJKLMNOPQRSTUVWXYZ',string(index_current:index_current))/=0)then

!  Skip characters like : + or the digits at the end of the word
   do index_endword=index_blank-1,index_current,-1
    if(index('ABCDEFGHIJKLMNOPQRSTUVWXYZ',string(index_endword:index_endword))/=0)exit
   end do

! DEBUG
! write(6,*)' index_current, index_blank, string',index_current,index_blank,string(index_current:index_endword)
! ENDDEBUG

   if(index(list_vars,string(index_current:index_endword)//blank)==0)then

! DEBUG
! write(6,*)' list_logicals,index=:',list_logicals,':',index(list_strings,string(index_current:index_endword)//blank)
! write(6,*)' list_strings,index=:',list_strings,':',index(list_strings,string(index_current:index_endword)//blank)
! ENDDEBUG


!   Treat possible logical input variables
    if(index(list_logicals,string(index_current:index_endword)//blank)/=0)then
     index_blank=index(string(index_current:),blank)+index_current-1
     if(index(' F T ',string(index_blank:index_blank+2))==0)then
      write(message, '(11a)' ) ch10,&
&      ' chkvars : ERROR - ',ch10,&
&      '  Found the token ',string(index_current:index_endword),' in the input file.',ch10,&
&      '  This variable should be given a logical value (T or F), but the following string was found :',&
&      string(index_blank:index_blank+2),ch10,&
&      '  Action : check your input file. You likely misused the input variable.'
      call wrtout(06,message,'COLL')
      call leave_new('COLL')
     else
      index_blank=index_blank+2
     end if
!    Treat possible string input variables
    else if(index(list_strings,string(index_current:index_endword)//blank)/=0)then
!    Every following string is accepted
     index_current=index(string(index_current:),blank)+index_current
     index_blank=index(string(index_current:),blank)+index_current-1

!    If still not admitted, then there is a problem
    else
     write(message, '(10a)' ) ch10,&
&     ' chkvars : ERROR - ',ch10,&
&     '  Found the token ',string(index_current:index_endword),' in the input file.',ch10,&
&     '  This name is not one of the registered input variable names (see the Web list of input variables).',ch10,&
&     '  Action : check your input file. You likely mistyped the input variable.'
     call wrtout(06,message,'COLL')
     call leave_new('COLL')
    end if
   end if
  end if
  index_current=index_blank+1
 end do

!DEBUG
!write(6,*)' chkvars : exit '
!stop
!ENDDEBUG

end subroutine chkvars
!!***
