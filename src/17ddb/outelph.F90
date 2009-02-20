!{\src2tex{textfont=tt}}
!!****f* ABINIT/outelph
!! NAME
!! outelph
!! structured variables.
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  elph_ds  the elph_type structured variable
!!  enunit   from the anaddb dataset 0 ==> Hartree and cm-1;
!!                                   1 ==> meV and Thz;
!!  nfile      unit number for writing
!!
!! OUTPUT
!!  only write
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine outelph(elph_ds,enunit,fname)

 use defs_basis
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_17ddb, except_this_one => outelph
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: enunit
 character(len=fnlen),intent(in) :: fname
 type(elph_type),intent(in) :: elph_ds

!Local variables-------------------------------
!scalars
 integer :: ibranch,ii,iost,iqfull,iqirr,isppol,jj,nfile,qmax,qnest_max
 integer :: qnest_min
 real(dp) :: lambda_q_max,lambda_qbranch_max,lambda_tot,nest_max,nest_min
 real(dp) :: omegalog_q,omegalog_qgrid,tc_macmill
 character(len=500) :: message
!arrays
 integer :: qbranch_max(2)
 integer,allocatable :: q_lambda_max(:),qbranch_lambda_max(:,:)
 real(dp),allocatable :: lambda_q(:,:),nestfactor(:),qirred(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' outelph : enter '
!ENDDEBUG

 if (enunit /= 0 .and. enunit /=1 .and. enunit /=2)  then
  write (message,'(4a,i9)')ch10, ' outelph : BUG-',ch10,&
&  ' enunit should be 0 or 1 or 2 while it is',enunit
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!==========================================================
!write header
!==========================================================
 nfile = 66
 open (unit=nfile,file=fname,form='formatted',status='unknown',iostat=iost)
 if (iost /= 0) then
  write (message,'(2a)')' outelph : ERROR- opening file ',trim(fname)
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 write(message,'(2a,80a,4a,80a)')ch10,' ',('=',ii=1,80),ch10,&
& ' Values of the parameters that define the electron-phonon calculation',ch10,&
& ' ',('=',ii=1,80)
 call wrtout(nfile,message,'COLL')

 write(message,'(a,i10,a,i10,a,i10)')&
& ' nFSkpt    = ',elph_ds%nFSkpt,   ' nFSkptirred = ',elph_ds%nFSkptirred,&
& ' nqpt      = ',elph_ds%nqpt
 call wrtout(nfile,message,'COLL')

 if (elph_ds%nsppol==1) then
  write(message,'(2a,f10.7,a,f10.6,a,f10.7)')ch10,&
&  ' Fermi DOS = ',elph_ds%n0(1),       ' Fermi level = ',elph_ds%fermie,&
&  ' mustar    = ',elph_ds%mustar
  call wrtout(nfile,message,'COLL')
 else if (elph_ds%nsppol==2) then
  write(message,'(2a,f10.7,f10.7,a,f10.6,a,f10.7)')ch10,&
&  ' Fermi DOS (up/dn) = ',elph_ds%n0(1),elph_ds%n0(2),       ' Fermi level = ',elph_ds%fermie,&
&  ' mustar    = ',elph_ds%mustar
  call wrtout(nfile,message,'COLL')
 else 
  stop 'outelph: Error : bad value for nsppol'
 end if

 write(message,'(2a,i10,a,i10,a,i10)')ch10,&
& ' minFSband = ',elph_ds%minFSband,' maxFSband   = ',elph_ds%maxFSband,&
& ' ngkkband  = ',elph_ds%ngkkband
 call wrtout(nfile,message,'COLL')

 write(message,'(80a,a)')('=',ii=1,80),ch10
 call wrtout(nfile,message,'COLL')

!==========================================================
!evaluate lambda and omega_log as a weighted sum over the q grid
!NOTE: in this part of the code atomic units are used
!==========================================================

 allocate (lambda_q(elph_ds%nqptirred,elph_ds%nsppol))
 lambda_q(:,:)=zero
 lambda_tot=zero  ; lambda_q_max=zero
 qmax=0           ; lambda_qbranch_max=zero
 qbranch_max(:)=1 ; omegalog_qgrid=zero

 do iqirr=1,elph_ds%nqptirred
  omegalog_q=zero

  do isppol=1,elph_ds%nsppol
   do ibranch=1,elph_ds%nbranch
!   find Max lambda(q,n)
    if (elph_ds%qgrid_data(iqirr,ibranch,isppol,3) > lambda_qbranch_max) then
     lambda_qbranch_max=elph_ds%qgrid_data(iqirr,ibranch,isppol,3)
     qbranch_max(1)=iqirr
     qbranch_max(2)=ibranch
    end if
    lambda_q(iqirr,isppol)=lambda_q(iqirr,isppol)+elph_ds%qgrid_data(iqirr,ibranch,isppol,3)
    if (abs(elph_ds%qgrid_data(iqirr,ibranch,isppol,1)) <= tol10) cycle
    omegalog_q=omegalog_q + elph_ds%qgrid_data(iqirr,ibranch,isppol,3)*log(abs(elph_ds%qgrid_data(iqirr,ibranch,isppol,1)))
   end do

   lambda_tot=lambda_tot+elph_ds%wtq(elph_ds%qirredtofull(iqirr))*lambda_q(iqirr,isppol)
   omegalog_qgrid=omegalog_qgrid+elph_ds%wtq(elph_ds%qirredtofull(iqirr))*omegalog_q


!  find Max lambda(q)
   if (lambda_q(iqirr,isppol) > lambda_q_max) then
    lambda_q_max=lambda_q(iqirr,isppol)
    qmax=iqirr
   end if
  end do

 end do !iqirr

 omegalog_qgrid=exp(omegalog_qgrid/lambda_tot)

 write (message,'(3a,2(a,es16.8))')                                                                              &
& ' Values of Lambda, Omega_log and Tc obtained using the weighted sum over the input Q-grid',ch10,ch10,&
& ' Isotropic Lambda = ',lambda_tot,'  Input mustar     = ',elph_ds%mustar
 call wrtout(nfile,message,'COLL')

 if (enunit==0) then !use hartree and cm-1
  write (message,'(2a,es16.8,a,es16.8,a)')ch10,                                           &
&  ' Omega_log        = ',omegalog_qgrid,' (Ha) ',omegalog_qgrid*Ha_cmm1,' (cm-1)'
  call wrtout(nfile,message,'COLL')
 else if (enunit==1) then !mev Thz
  write (message,'(2a,es16.8,a,es16.8,a)')ch10,                                                     &
&  ' Omega_log        = ',omegalog_qgrid*Ha_eV/1000._dp,' (meV) ',omegalog_qgrid*Ha_THz,' (THz)'
  call wrtout(nfile,message,'COLL')
 else !hartree,cm-1,mev,Thz,kelvin
  write (message,'(2a,es16.8,a,es16.8,3a,es16.8,a,es16.8,3a,es16.8,a)')ch10,                        &
&  ' Omega_log        = ',omegalog_qgrid,' (Ha)  ',omegalog_qgrid*Ha_cmm1,' (cm-1)',ch10,            &
&  '                  = ',omegalog_qgrid*Ha_eV/1000._dp,' (meV) ',omegalog_qgrid*Ha_THz,' (THz)',ch10,&
&  '                  = ',omegalog_qgrid*Ha_K,' (K) '
  call wrtout(nfile,message,'COLL')
 end if

 tc_macmill = omegalog_qgrid/1.2_dp&
& *exp((-1.04_dp*(one+lambda_tot)) / (lambda_tot-elph_ds%mustar*(one+0.62_dp*lambda_tot)))

 if (enunit==0) then !use hartree and cm-1
  write (message,'(2a,es16.8,a,es16.8,2a)')ch10,&
&  ' MacMillan Tc     = ',tc_macmill,' (Ha) ',tc_macmill*Ha_cmm1,' (cm-1) ',ch10
  call wrtout(nfile,message,'COLL')
 else if (enunit==1) then !use mev and Thz
  write (message,'(2a,es16.8,a,es16.8,2a)')ch10,&
&  ' MacMillan Tc     = ',tc_macmill*Ha_eV/1000._dp,' (meV) ',tc_macmill*Ha_THz,' (THz) ',ch10
  call wrtout(nfile,message,'COLL')
 else !use hartree,cm-1,mev,Thz,kelvin
  write (message,'(2a,es16.8,a,es16.8,3a,es16.8,a,es16.8,3a,es16.8,2a)')ch10,                 &
&  ' MacMillan Tc     = ',tc_macmill,' (Ha)  ',tc_macmill*Ha_cmm1,' (cm-1) ',ch10,            &
&  '                  = ',tc_macmill*Ha_eV/1000._dp,' (meV) ',tc_macmill*Ha_THz,' (THz) ',ch10,&
&  '                  = ',tc_macmill*Ha_K,' (K) ',ch10
  call wrtout(nfile,message,'COLL')
 end if

!==========================================================
!output lambda(q) values for each q point in the irred grid
!==========================================================

 write(message,'(2a)')' Irreducible q-points and corresponding Lambda(q)',ch10
 call wrtout(nfile,message,'COLL')

 do isppol=1,elph_ds%nsppol
  write(message,'(a,i3,2a)')'  === isppol ', isppol,' === ',ch10
  call wrtout(nfile,message,'COLL')
  do iqirr=1,elph_ds%nqptirred
   iqfull=elph_ds%qirredtofull(iqirr)
   write(message,'(i5,a,3(es16.8,1x),a,es16.8,a)')&
     & iqfull,') ',elph_ds%spqpt(:,iqfull),'(',lambda_q(iqirr,isppol),'  )'
   call wrtout(nfile,message,'COLL')
  end do
 end do

!use same indexing as that used for the full q-grid
 qmax=elph_ds%qirredtofull(qmax)
 qbranch_max(1)=elph_ds%qirredtofull(qbranch_max(1))

 write (message,'(2a,es16.8,a,i6,3a,es16.8,a,i6,a,i4)')ch10,            &
& ' Max lambda(q)      = ',lambda_q_max,      ' at qpt ',qmax,')',ch10, &
& ' Max lambda(q,n)    = ',lambda_qbranch_max,' at qpt ',qbranch_max(1),&
& ') and Mode number ',qbranch_max(2)
 call wrtout(nfile,message,'COLL')

!==========================================================
!evaluation of the nesting-factor over the irreducible q grid.
!==========================================================

!fill irreducile q-grid
 allocate (qirred(3,elph_ds%nqptirred))
 qirred(:,:)=zero

 do iqirr=1,elph_ds%nqptirred
  qirred(:,iqirr)=elph_ds%spqpt(:,elph_ds%qirredtofull(iqirr))
 end do

 allocate (nestfactor(elph_ds%nqptirred))
 nestfactor(:)=zero

!NOTE: weights are not normalised, the normalisation factor in reintroduced in bfactor
 call bfactor(elph_ds%nFSkpt,elph_ds%FSkpt,elph_ds%nqptirred,qirred,elph_ds%FSintweight,&
& elph_ds%nFSband,nestfactor)

 deallocate(qirred)

!find Max and min of the nesting factor
!NOTE maxloc and minloc are arrays so they cannot be used in the formatted output
!anyway the size of nestfactor is not so huge!!!
 nest_max=maxval(nestfactor); nest_min=minval(nestfactor)

 qnest_max=0
 do iqirr=1,elph_ds%nqptirred
  if (nestfactor(iqirr)==nest_max) then
   qnest_max=iqirr
   exit
  end if
 end do

 qnest_min=0
 do iqirr=1,elph_ds%nqptirred
  if (nestfactor(iqirr)==nest_min) then
   qnest_min=iqirr
   exit
  end if
 end do

 write (*,*) maxloc(nestfactor),minloc(nestfactor)
 write(message,'(a,(a,es16.8,a,i6,a),a,(a,es16.8,a,i6,a))')ch10,  &
& ' Max nesting factor = ',nest_max,' at qpt ',qnest_max,') ',ch10,&
& ' min nesting factor = ',nest_min,' at qpt ',qnest_min,') '
 call wrtout(nfile,message,'COLL')

!==========================================================
!Write ph-linewidths and lambda(q,n) obtained before the
!Fourier interpolation
!==========================================================

 write (message,'(2a)')ch10,&
& ' Phonon frequencies, linewidths and e-ph coefficients for each irreducible q point '
 call wrtout(nfile,message,'COLL')

 do isppol=1,elph_ds%nsppol
  write (message,'(a,i3,a)') '========= quantities for isppol = ', isppol, ' ================='
  call wrtout(nfile,message,'COLL')
  do iqirr=1,elph_ds%nqptirred
!  same numbering as that used for irred q points
   iqfull=elph_ds%qirredtofull(iqirr)
   write(*,*) 'iqfull = ', iqfull
   write(message,'(64a,i6,a,3(es16.8),3a,es16.8,a,es16.8,2a,es16.8,a,f8.3,65a)')ch10,&
&   ' ',('=',jj=1,60),ch10,&
&   ' qpt ',iqfull,') ',elph_ds%spqpt(:,iqfull),ch10,ch10,&
&   ' Weight    = ',elph_ds%wtq(iqfull),'    Lambda(q,isppol) = ',lambda_q(iqirr,isppol),ch10,&
&   ' Nest fact = ',nestfactor(iqirr),'    (',100*nestfactor(iqirr)/nest_max,' % of max_value )',ch10,&
&   ' ',('=',jj=1,60),ch10,' Mode number    Frequency       Linewidth        Lambda(q,n)'
   call wrtout(nfile,message,'COLL')

!  use units according to enunit
   if (enunit==0 .or. enunit==2) then !hartree and cm-1
    write(message,'(63a)')' ',('-',jj=1,60),ch10,&
    '                  (Ha)             (Ha)'
    call wrtout(nfile,message,'COLL')
    do ibranch=1,elph_ds%nbranch
!    branch index, frequency, linewidth, lamda(q,n) (hartree units)
     write(message,'(i6,5x,3(es16.8,1x))' )ibranch,(elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,3)
     call wrtout(nfile,message,'COLL')
    end do
    write(message,'(63a)')' ',('-',jj=1,60),ch10,&
&    '                 (cm-1)           (cm-1)'
    call wrtout(nfile,message,'COLL')
    do ibranch=1,elph_ds%nbranch
!    branch index, frequency, linewidth (in cm-1)
     write(message,'(i6,5x,2(es16.8,1x))' )ibranch,(Ha_cmm1*elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,2)
     call wrtout(nfile,message,'COLL')
    end do
   end if !hartree and cm-1

   if (enunit==2 .or. enunit==1) then !write also meV Thz and Kelvin
    write(message,'(63a)')' ',('-',jj=1,60),ch10,&
&    '                 (meV)             (meV)'
    call wrtout(nfile,message,'COLL')
    if (enunit == 1 ) then !write also lambda values
     do ibranch=1,elph_ds%nbranch
!     branch index, frequency, linewidth, lamda(q,n) (mev units)
      write(message,'(i6,5x,3(es16.8,1x))' )ibranch,((Ha_eV/1000._dp)*elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,2),&
&      elph_ds%qgrid_data(iqirr,ibranch,isppol,3)
      call wrtout(nfile,message,'COLL')
     end do
    else !do not write lambda values
     do ibranch=1,elph_ds%nbranch
!     branch index, frequency, linewidth (in meV)
      write(message,'(i6,5x,2(es16.8,1x))' )ibranch,((Ha_eV/1000._dp)*elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,2)
      call wrtout(nfile,message,'COLL')
     end do
    end if

    write(message,'(63a)')' ',('-',jj=1,60),ch10,&
&    '                 (Thz)             (Thz)'
    call wrtout(nfile,message,'COLL')
    do ibranch=1,elph_ds%nbranch
!    branch index, frequency, linewidth (in Thz)
     write(message,'(i6,5x,2(es16.8,1x))' )ibranch,(Ha_THz*elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,2)
     call wrtout(nfile,message,'COLL')
    end do

    if (enunit == 2 ) then !kelvin
     write(message,'(63a)')' ',('-',jj=1,60),ch10,&
&     '                  (K)               (K)'
     call wrtout(nfile,message,'COLL')
     do ibranch=1,elph_ds%nbranch
!     branch index, frequency, linewidth (in Kelvin)
      write(message,'(i6,5x,2(es16.8,1x))' )ibranch,(Ha_K*elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,2)
      call wrtout(nfile,message,'COLL')
     end do
    end if !kelvin

   end if  !end write also meV Thz and Kelvin

   write(message,'(62a)')' ',('=',jj=1,60),ch10
   call wrtout(nfile,message,'COLL')

  end do !nqptirred
 end do !nsppol

 deallocate(nestfactor)
 deallocate (lambda_q)

 close (nfile)
!DEBUG
!write(6,*)' out_elph : exit'
!ENDDEBUG

end subroutine outelph
!!***
