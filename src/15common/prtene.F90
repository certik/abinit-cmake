!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtene
!!
!! NAME
!! prtene
!!
!! FUNCTION
!! Print components of total energy in nice format
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, LBoeri, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryphase
!!   | kptopt
!!   | occopt
!!   | tphysel="physical" electronic temperature with FD occupations
!!   | tsmear=smearing energy or temperature (if metal)
!!  energies <type(energies_type)>=values of parts of total energy
!!  iout=unit number to which output is written
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine prtene(dtset,energies,iout,usepaw)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,usepaw
 type(dataset_type),intent(in) :: dtset
 type(energies_type),intent(in) :: energies

!Local variables-------------------------------
! Do not modify the length of this string
!scalars
 integer :: mu,optdc
 logical :: directE_avail
 real(dp) :: eent,eintern,einterndc,enevalue,etotal,etotaldc
 character(len=22) :: eneName
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(iout, '(a)') ' prtene : enter '
!ENDDEBUG

 optdc=-1
!if (abs(energies%e_localpsp)>1.e-15_dp.and.abs(energies%e_xcdc)<1.e-15_dp) optdc=0
 if (abs(energies%e_xcdc)<1.e-15_dp) optdc=0
 if (abs(energies%e_localpsp)<1.e-15_dp.and.abs(energies%e_xcdc)>1.e-15_dp) optdc=1
 if (abs(energies%e_localpsp)>1.e-15_dp.and.abs(energies%e_xcdc)>1.e-15_dp) optdc=2
 directE_avail=(usepaw==0.or.dtset%pawspnorb==0.or.dtset%kptopt==1.or.dtset%kptopt==2)

!============= Evaluate some parts of the energy ===========

!Here, treat the case of metals
!In re-smeared case the free energy is defined with tphysel
 if(dtset%occopt>=3 .and. dtset%occopt<=7)then
  if (abs(dtset%tphysel) < tol10) then
   eent=-dtset%tsmear * energies%entropy
  else
   eent=-dtset%tphysel * energies%entropy
  end if
 else
  eent=zero
 end if

 if (optdc==0.or.optdc==2) then
  etotal = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
&  energies%e_localpsp + energies%e_corepsp + eent
  if (dtset%usewvl == 1) etotal = etotal - energies%e_vxc
  if (usepaw==0) etotal = etotal + energies%e_nonlocalpsp
  if (usepaw==1) etotal = etotal + energies%e_paw
  if (dtset%berryopt==4) etotal=etotal+energies%e_elecfield
  etotal = etotal + energies%e_ewald
 end if
 if (optdc>=1) then
  etotaldc = energies%e_eigenvalues - energies%e_hartree + energies%e_xc - &
&  energies%e_xcdc + energies%e_corepsp + eent
  if (usepaw==1) etotaldc = etotaldc + energies%e_pawdc
  if (dtset%berryopt==4) etotaldc = etotaldc + energies%e_elecfield
  etotaldc = etotaldc + energies%e_ewald
 end if

 write(message,'(a,80a)') ch10,('-',mu=1,80)
 call wrtout(iout,message,'COLL')

!============= Printing of Etotal by direct scheme ===========

 if (dtset%icoulomb == 1) then
  write(eneName, "(A)") "    Ion-ion energy  = "
 else
  write(eneName, "(A)") "    Ewald energy    = "
 end if
 enevalue = energies%e_ewald

 if (optdc==0.or.optdc==2) then

  if (directE_avail) then
   write(message, '(2a)' ) &
&   ' Components of total free energy (in Hartree) :',ch10
   call wrtout(iout,message,'COLL')
   write(message, '(5(a,es21.14,a),a,es21.14)' ) &
&   '    Kinetic energy  = ',energies%e_kinetic,ch10,&
&   '    Hartree energy  = ',energies%e_hartree,ch10,&
&   '    XC energy       = ',energies%e_xc,ch10,&
&   eneName,                 enevalue,ch10,&
&   '    PspCore energy  = ',energies%e_corepsp,ch10,&
&   '    Loc. psp. energy= ',energies%e_localpsp
   call wrtout(iout,message,'COLL')
   if (dtset%berryopt == 4) then
    write(message, '(a,es21.14)' ) &
&    '    Electric field  = ',energies%e_elecfield
    call wrtout(iout,message,'COLL')
   end if
   if (usepaw==0) then
    write(message, '(a,es21.14)' ) &
&    '    NL   psp  energy= ',energies%e_nonlocalpsp
   else
    write(message, '(a,es21.14)' ) &
&    '    Spherical terms = ',energies%e_paw
   end if
   call wrtout(iout,message,'COLL')
   if (dtset%usewvl == 1) then
    write(message, '(a,es21.14)' ) &
&    '    XC pot. energy  = ', energies%e_vxc
    call wrtout(iout,message,'COLL')
   end if
   if(dtset%occopt>=3 .and. dtset%occopt<=7) then
    write(message, '(a,es21.14,a,a,a,es21.14)' ) &
&    '    >>>>> Internal E= ',etotal-eent,ch10,ch10,&
&    '    -kT*entropy     = ',eent
    call wrtout(iout,message,'COLL')
   end if
   write(message, '(a,es21.14)' ) &
&   '    >>>>>>>>> Etotal= ',etotal
   call wrtout(iout,message,'COLL')
   if (abs(energies%e_pulay) >=1.0d-12) then
    write(message, '(a,a,es21.14,a,a,es21.14)' ) ch10,&
&    '    Pulay correction= ',energies%e_pulay,ch10,&
&    '    >Etot with Pulay= ',etotal+energies%e_pulay
    call wrtout(iout,message,'COLL')
   end if

  else
   write(message, '(9a)' ) &
&   ' COMMENT: ',ch10,&
&   '  "Direct" decomposition of total free energy cannot be printed out !!!',ch10,&
&   '  PAW contribution due to spin-orbit coupling cannot be evaluated',ch10,&
&   '  without the knowledge of imaginary part of Rhoij atomic occupancies',ch10,&
&   '  (computed only when pawcpxocc=2).'
   call wrtout(iout,message,'COLL')
  end if

 end if
!============= Printing of Etotal by double-counting scheme ===========

 if (optdc>=1) then

  write(message, '(a,a,a,3(a,es21.14,a),a,es21.14)' ) ch10,&
&  ' "Double-counting" decomposition of free energy:',ch10,&
&  '    Band energy     = ',energies%e_eigenvalues,ch10,&
&  eneName,                 enevalue,ch10,&
&  '    PspCore energy  = ',energies%e_corepsp,ch10,&
&  '    Dble-C XC-energy= ',-energies%e_hartree+energies%e_xc-energies%e_xcdc
  call wrtout(iout,message,'COLL')
  if (dtset%berryopt == 4) then
   write(message, '(a,es21.14)' ) &
&   '    Electric field  = ',energies%e_elecfield
   call wrtout(iout,message,'COLL')
  end if
  if (usepaw==1) then
   write(message, '(a,es21.14)' ) &
&   '    Spherical terms = ',energies%e_pawdc
   call wrtout(iout,message,'COLL')
  end if
  if(dtset%occopt>=3 .and. dtset%occopt<=7) then
   write(message, '(a,es21.14,a,a,a,es21.14)' ) &
&   '    >>>>> Internal E= ',etotaldc-eent,ch10,ch10,&
&   '    -kT*entropy     = ',eent
   call wrtout(iout,message,'COLL')
  end if
  write(message, '(a,es21.14)' ) &
&  '    >>>> Etotal (DC)= ',etotaldc
  call wrtout(iout,message,'COLL')
  if (abs(energies%e_pulay) >=1.0d-12) then
   write(message, '(a,a,es21.14,a,a,es21.14)' ) ch10,&
&   '    Pulay correction= ',energies%e_pulay,ch10,&
&   '    >EDC with Pulay= ',etotaldc+energies%e_pulay
   call wrtout(iout,message,'COLL')
  end if
 end if

!======= Additional printing for compatibility  ==========

 if (usepaw==0.and.optdc==0) then
  write(message, '(a,a,a,a,es21.14,a,es18.10)' ) ch10,&
&  ' Other information on the energy :',ch10,&
&  '    Total energy(eV)= ',etotal*Ha_eV,' ; Band energy (Ha)= ',energies%e_eigenvalues
  call wrtout(iout,message,'COLL')
 end if

 if ((optdc==0.or.optdc==2).and.(.not.directE_avail)) then
  write(message, '(a,a,es18.10)' ) ch10,&
&  ' Band energy (Ha)= ',energies%e_eigenvalues
  call wrtout(iout,message,'COLL')
 end if

 if (usepaw==1) then
  if ((optdc==0.or.optdc==2).and.(directE_avail)) then
   write(message, '(a,a,es21.14)' ) ch10,&
&   '  >Total energy in eV           = ',etotal*Ha_eV
   call wrtout(iout,message,'COLL')
  end if
  if (optdc>=1) then
   if (optdc==1) write(message, '(a,a,es21.14)' ) ch10,&
&   '  >Total DC energy in eV        = ',etotaldc*Ha_eV
   if (optdc==2) write(message, '(a,es21.14)' ) &
&   '  >Total DC energy in eV        = ',etotaldc*Ha_eV
   call wrtout(iout,message,'COLL')
  end if
 end if

!=============
 write(message,'(a,80a)')('-',mu=1,80)
 call wrtout(iout,message,'COLL')

end subroutine prtene
!!***
