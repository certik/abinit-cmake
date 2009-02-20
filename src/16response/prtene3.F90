!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtene3
!!
!! NAME
!! prtene3
!!
!! FUNCTION
!! Print components of second derivative of total energy in nice format
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! eberry=energy associated with Berry phase
!! edocc=correction to 2nd-order total energy coming from changes of occupation
!! eeig0=0th-order eigenenergies part of 2nd-order total energy
!! eew=Ewald part of 2nd-order total energy
!! efrhar=hartree frozen-wavefunction part of 2nd-order tot. en.
!! efrkin=kinetic frozen-wavefunction part of 2nd-order tot. en.
!! efrloc=local psp. frozen-wavefunction part of 2nd-order tot. en.
!! efrnl=nonlocal psp. frozen-wavefunction part of 2nd-order tot. en
!! efrx1=xc core corr.(1) frozen-wavefunction part of 2nd-order tot. en
!! efrx2=xc core corr.(2) frozen-wavefunction part of 2nd-order tot. en
!! ehart01=inhomogeneous 1st-order Hartree part of 2nd-order total energy
!!   for strain perturbation only (zero otherwise, and not used)
!! ehart1=1st-order Hartree part of 2nd-order total energy
!! eii=pseudopotential core part of 2nd-order total energy
!! ek0=0th-order kinetic energy part of 2nd-order total energy.
!! ek1=1st-order kinetic energy part of 2nd-order total energy.
!! eloc0=0th-order local (psp+vxc+Hart) part of 2nd-order total energy
!! elpsp1=1st-order local pseudopot. part of 2nd-order total energy.
!! enl0=0th-order nonlocal pseudopot. part of 2nd-order total energy.
!! enl1=1st-order nonlocal pseudopot. part of 2nd-order total energy.
!! exc1=1st-order exchange-correlation part of 2nd-order total energy
!! iout=unit number to which output is written
!! ipert=type of the perturbation
!! natom=number of atoms in unit cell
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! all energies in Hartree
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine prtene3(berryopt,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,& 
&  ehart01,ehart1,eii,ek0,ek1,eloc0,elpsp1,enl0,enl1,exc1,iout,ipert,natom)  

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: berryopt,iout,ipert,natom
 real(dp),intent(in) :: eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1
 real(dp),intent(in) :: efrx2,ehart01,ehart1,eii,ek0,ek1,eloc0,elpsp1,enl0,exc1
 real(dp),intent(inout) :: enl1

!Local variables -------------------------
!scalars
 real(dp) :: enl1_effective,erelax,etotal
 character(len=500) :: message

! *********************************************************************

 if (ipert>=1.and.ipert<=natom) then
  write(message, '(a,a,a)' )&
&  ch10,' Thirteen components of 2nd-order ',&
&  'total energy (hartree) are '
 else if (ipert==natom+1)then
  write(message, '(a,a,a)' )&
&  ch10,' Eight components of 2nd-order ',&
&  'total energy (hartree) are '
 else if (ipert==natom+2 .and. berryopt /=4 )then   
  write(message, '(a,a,a)' )&
&  ch10,' Seven components of 2nd-order ',&
&  'total energy (hartree) are '
 else if (ipert==natom+2 .and. berryopt ==4 )then   
  write(message, '(a,a,a)' )&
&  ch10,' Seven components of 2nd-order ',&
&  'total energy (hartree) are '
 else if (ipert==natom+3 .or. ipert==natom+4)then
  write(message, '(a,a,a)' )&
&  ch10,' Seventeen components of 2nd-order ',&
&  'total energy (hartree) are '
 else if (ipert==natom+5 )then   
  write(message, '(a,a,a)' )&
&  ch10,' Seven components of 2nd-order ',&
&  'total energy (hartree) are '
 end if
 call wrtout(iout,message,'COLL')
 call wrtout(6,message,'COLL')
 write(message, '(a)' )&
& ' 1,2,3: 0th-order hamiltonian combined with 1st-order wavefunctions'
 call wrtout(iout,message,'COLL')
 call wrtout(6,message,'COLL')
 write(message, '(a,es17.8,a,es17.8,a,es17.8)' )&
& '     kin0=',ek0,   ' eigvalue=',eeig0,'  local=',eloc0
 call wrtout(iout,message,'COLL')
 call wrtout(6,message,'COLL')

 if(ipert==natom+3 .or. ipert==natom+4) then
  write(message, '(a)' )&
&  ' 4,5,6,7: 1st-order hamiltonian combined with 1st and 0th-order wfs'
 else
  write(message, '(a)' )&
&  ' 4,5,6: 1st-order hamiltonian combined with 1st and 0th-order wfs'
 end if
 call wrtout(iout,message,'COLL')
 call wrtout(6,message,'COLL')

 if(ipert/=natom+1.and.ipert/=natom+2.and.ipert/=natom+5)then
  write(message, '(a,es17.8,a,es17.8,a,es17.8,a,a)' ) &
&  ' loc psp =',elpsp1,'  Hartree=',ehart1,'     xc=',exc1,ch10,&
&  ' note that "loc psp" includes a xc core correction that could be resolved'
 else if(ipert==natom+1) then
  write(message, '(a,es17.8,a,es17.8,a,es17.8)' ) &
&  '     kin1=',ek1,   '  Hartree=',ehart1,'     xc=',exc1
 else if(ipert==natom+2 .or. ipert==natom+5 ) then
  write(message, '(a,es17.8,a,es17.8,a,es17.8)' ) &
&  '    dotwf=',enl1,  '  Hartree=',ehart1,'     xc=',exc1
 end if
 if(ipert==natom+3 .or. ipert==natom+4) then
  write(message, '(a,es17.8,a,es17.8,a,es17.8,a,a,es17.8)' ) &
&  ' loc psp =',elpsp1,'  Hartree=',ehart1,'     xc=',exc1,ch10,&
&  '     kin1=',ek1
 end if
 call wrtout(iout,message,'COLL')
 call wrtout(6,message,'COLL')

 if(ipert>=1.and.ipert<=natom)then
  erelax=ek0+edocc+eeig0+eloc0+elpsp1+ehart1+exc1+enl0+enl1
 else if(ipert==natom+1.or.ipert==natom+2.or.ipert==natom+5)then
  erelax=ek0+edocc+eeig0+eloc0+ek1+ehart1+exc1+enl0+enl1
  if(1.0_dp+enl1/10.0_dp==1.0_dp)enl1=0.0_dp
 else if(ipert==natom+3.or.ipert==natom+4)then
  erelax=ek0+edocc+eeig0+eloc0+ek1+elpsp1+ehart1+exc1+enl0+enl1
 end if

 enl1_effective=enl1
 if(ipert==natom+2)enl1_effective=0.0_dp
 if(ipert==natom+5)enl1_effective=0.0_dp
 if(ipert==natom+3 .or. ipert==natom+4) then
  write(message, '(a,a,a,es17.8,a,es17.8,a,es17.8)' )&
&  ' 8,9,10: eventually, occupation + non-local contributions',ch10,&
&  '    edocc=',edocc,'     enl0=',enl0,'   enl1=',enl1_effective
 else
  write(message, '(a,a,a,es17.8,a,es17.8,a,es17.8)' )&
&  ' 7,8,9: eventually, occupation + non-local contributions',ch10,&
&  '    edocc=',edocc,'     enl0=',enl0,'   enl1=',enl1_effective
 end if
 call wrtout(iout,message,'COLL')
 call wrtout(6,message,'COLL')

 if(ipert==natom+3 .or. ipert==natom+4) then
  write(message, '(a,a,a,es17.8)' )&
&  ' 1-10 gives the relaxation energy (to be shifted if some occ is /=2.0)',&
&  ch10,'   erelax=',erelax
 else
  write(message, '(a,a,a,es17.8)' )&
&  ' 1-9 gives the relaxation energy (to be shifted if some occ is /=2.0)',&
&  ch10,'   erelax=',erelax
 end if
 call wrtout(iout,message,'COLL')
 call wrtout(6,message,'COLL')

 if(ipert>=1.and.ipert<=natom)then
  write(message, '(a,a)' )&
&  ' 10,11,12 Non-relaxation  contributions : ',&
&  'frozen-wavefunctions and Ewald'
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  write(message, '(a,es17.8,a,es17.8,a,es17.8)' ) &
&  ' fr.local=',efrloc,' fr.nonlo=',efrnl,'  Ewald=',eew
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')

  write(message, '(a,es16.6)' )' prtene3 : non-relax=',efrloc+efrnl+eew
  call wrtout(6,message,'COLL')

  write(message, '(a)' )&
&  ' 13,14 Frozen wf xc core corrections (1) and (2)'
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  write(message, '(a,es17.8,a,es17.8)' ) &
&  ' frxc 1  =',efrx1,'  frxc 2 =',efrx2
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  write(message, '(a)' )' Resulting in : '
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  etotal=erelax+eew+efrloc+efrnl+efrx1+efrx2
  write(message, '(a,e20.10,a,e22.12,a)' ) &
&  ' 2DEtotal=',etotal,' Ha. Also 2DEtotal=',etotal*Ha_eV,' eV'
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  write(message, '(a,es20.10,a,es20.10,a)' ) &
&  '    (2DErelax=',erelax,' Ha. 2DEnonrelax=',etotal-erelax,' Ha)'
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  write(message, '(a,es20.10,a,a)' ) &
&  '    (  non-var. 2DEtotal :',&
&  0.5_dp*(elpsp1+enl1)+eew+efrloc+efrnl+efrx1+efrx2,' Ha)',ch10
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
 else if(ipert==natom+1.or.ipert==natom+2.or.ipert==natom+5)then
  write(message,*)' No Ewald or frozen-wf contrib.:',&
&  ' the relaxation energy is the total one'
  if( berryopt == 4 )then
   write(message,'(a,es20.10)')'Berry phase energy :',eberry
  end if
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  etotal=erelax
  write(message, '(a,e20.10,a,e22.12,a)' ) &
&  ' 2DEtotal=',etotal,' Ha. Also 2DEtotal=',etotal*Ha_eV,' eV'
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  write(message, '(a,es20.10,a)' ) &
&  '    (  non-var. 2DEtotal :',0.5_dp*(ek1+enl1),' Ha)'
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
 else if(ipert==natom+3 .or. ipert==natom+4) then
  write(message, '(a,a)' )&
&  ' 11,12,13 Non-relaxation  contributions : ',&
&  'frozen-wavefunctions and Ewald'
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  write(message, '(a,es17.8,a,es17.8,a,es17.8)' ) &
&  '  fr.hart=',efrhar,'   fr.kin=',efrkin,' fr.loc=',efrloc
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')

  write(message, '(a,a)' )&
&  ' 14,15,16 Non-relaxation  contributions : ',&
&  'frozen-wavefunctions and Ewald'
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  write(message, '(a,es17.8,a,es17.8,a,es17.8)' ) &
&  '  fr.nonl=',efrnl,'    fr.xc=',efrx1,'  Ewald=',eew
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  write(message, '(a,a)' )&
&  ' 17 Non-relaxation  contributions : ',&
&  'pseudopotential core energy'
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  write(message, '(a,es17.8)' ) &
&  '  pspcore=',eii
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  write(message, '(a,es16.6)' )' prtene3 : non-relax=',&
&  efrhar+efrkin+efrloc+efrnl+efrx1+eew
  call wrtout(6,message,'COLL')

  write(message, '(a)' )' Resulting in : '
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  etotal=erelax+efrhar+efrkin+efrloc+efrnl+efrx1+eew+eii
  write(message, '(a,e20.10,a,e22.12,a)' ) &
&  ' 2DEtotal=',etotal,' Ha. Also 2DEtotal=',etotal*Ha_eV,' eV'
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  write(message, '(a,es20.10,a,es20.10,a)' ) &
&  '    (2DErelax=',erelax,' Ha. 2DEnonrelax=',etotal-erelax,' Ha)'
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
  write(message, '(a,es20.10,a,a)' ) &
&  '    (  non-var. 2DEtotal :',&
&  0.5_dp*(elpsp1+enl1+ek1+ehart01)+&
&  efrhar+efrkin+efrloc+efrnl+efrx1+eew+eii,' Ha)',ch10
  call wrtout(iout,message,'COLL')
  call wrtout(6,message,'COLL')
 end if

end subroutine prtene3
!!***
