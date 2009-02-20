!{\src2tex{textfont=tt}}
!!****f* ABINIT/clnup1
!! NAME
!! clnup1
!!
!!
!! FUNCTION
!! Perform "cleanup" at end of execution of gstate routine.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  acell(3)=length scales of primitive translations (bohr)
!!  dosdeltae=DOS delta of Energy
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  eigen(mband*nkpt*nsppol)=eigenvalues (hartree) for all bands at each k point
!!  enunit=choice for units of output eigenvalues: 0=>hartree,
!!   1=> eV, 2=> hartree and eV
!!  fermie=fermi energy (Hartree)
!!  filnam= character string giving the root to form the name of the DOS file
!!   and to form the name of the output WFK or WFQ file.
!!  fred(3,natom)=d(E)/d(xred) (hartree)
!!  iatfix(3,natom)=0 if not fixed along specified direction, 1 if fixed
!!  iscf=parameter controlling scf or non-scf choice
!!  kptopt=option for the generation of k points
!!  kptns(3,nkpt)=k points in terms of recip primitive translations
!!  mband=maximum number of bands
!!  mkmem=maximum number of k-points in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum number of plane waves
!!  natom=number of atoms in unit cell
!!  nband(nkpt*nsppol)=number of bands
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nstep=desired number of electron iteration steps
!!  occ(maxval(nband(:))*nkpt*nsppol)=occupancies for each band and k point
!!  occopt=option for occupancies
!!  prtdos= if == 1, will print the density of states
!!  prtfor= if >0, will print the forces
!!  prtstm= input variable prtstm
!!  prtvol=control print volume and debugging
!!  resid(mband*nkpt*nsppol)=squared residuals for each band and k point
!!   where resid(n,k)=|<C(n,k)|(H-e(n,k))|C(n,k)>|^2
!!  rhor(nfft,nspden)=electron density (electrons/bohr^3)
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  tphysel="physical" electronic temperature with FD occupations
!!  tsmear=smearing energy or temperature (if metal)
!!  vxcavg=average of vxc potential
!!  wtk(nkpt)=real(dp) array of k-point weights
!!  xred(3,natom)=reduced atomic coordinates
!!
!! OUTPUT
!!  (only print and write to disk)
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      getnel,prteigrs,prtrhomxmn,prtxf,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine clnup1(acell,dosdeltae,dtset,eigen,enunit,&
&  fermie,filnam,fred,iatfix,iscf,kptns,kptopt,mband,mkmem,mpi_enreg,mpw,&
&  natom,nband,nfft,ngfft,nkpt,nspden,nspinor,nsppol,nstep,occ,occopt,prtdos,&
&  prteig,prtfor,prtstm,prtvol,resid,rhor,rprimd,tphysel,tsmear,vxcavg,wtk,xred)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_14occeig
 use interfaces_15common, except_this_one => clnup1
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: enunit,iscf,kptopt,mband,mkmem,mpw,natom,nfft,nkpt
 integer,intent(in) :: nspden,nspinor,nsppol,nstep,occopt,prtdos,prteig,prtfor
 integer,intent(in) :: prtstm,prtvol
 real(dp),intent(in) :: dosdeltae,fermie,tphysel,tsmear,vxcavg
 character(len=fnlen),intent(in) :: filnam
 type(dataset_type),intent(inout) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: iatfix(3,natom),nband(nkpt*nsppol),ngfft(18)
 real(dp),intent(in) :: acell(3),eigen(mband*nkpt*nsppol),fred(3,natom)
 real(dp),intent(in) :: kptns(3,nkpt),resid(mband*nkpt*nsppol)
 real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3),wtk(nkpt),xred(3,natom)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iscf_dum,iwfrc,nnonsc,option,unitdos
 real(dp) :: entropy,grmax,grsum,maxocc,nelect,tolwf
 logical :: tread
 character(len=500) :: message
 character(len=fnlen) :: fildos
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: doccde(:)

! *************************************************************************

!DEBUG
!write(6,*)' clnup1 : enter,  rhor(1,1),rhor(1,2)',rhor(1,1),rhor(1,2)
!stop
!ENDDEBUG

 if(prtstm==0)then
! Write reduced coordinates xred
  write(message, '(a,i5,a)' ) &
&  ' reduced coordinates (array xred) for',natom,' atoms'
  call wrtout(ab_out,message,'COLL')
  do iatom=1,natom
   write(message, '(1x,3f20.12)' ) xred(:,iatom)
   call wrtout(ab_out,message,'COLL')
  end do
 end if

!Write reduced gradients if iscf > 0 and nstep>0 and prtstm==0
 if (iscf>0.and.nstep>0.and.prtstm==0) then

! Compute absolute maximum and root mean square value of gradients
  grmax=0.0_dp
  grsum=0.0_dp
  do iatom=1,natom
   do ii=1,3
!   To be activated in v5.5
!   grmax=max(grmax,abs(fred(ii,iatom)))
    grmax=max(grmax,fred(ii,iatom))
    grsum=grsum+fred(ii,iatom)**2
   end do
  end do
  grsum=sqrt(grsum/dble(3*natom))

  write(message, '(1x,a,1p,e12.4,a,e12.4,a)' ) &
&  'rms dE/dt=',grsum,'; max dE/dt=',grmax,&
&  '; dE/dt below (all hartree)'
  call wrtout(ab_out,message,'COLL')
  do iatom=1,natom
   write(message, '(i5,1x,3f20.12)' ) iatom,fred(1:3,iatom)
   call wrtout(ab_out,message,'COLL')
  end do

 end if

 if(prtstm==0)then

! Compute and write out dimensional cartesian coords and forces:
  write(message,*)' '
  call wrtout(ab_out,message,'COLL')

! (only write forces if iscf > 0 and nstep>0)
  if (iscf<=0.or.nstep<=0.or.prtfor==0) then
   iwfrc=0
  else
   iwfrc=1
  end if

  call prtxf(fred,iatfix,ab_out,iwfrc,natom,rprimd,xred)

! Write length scales
  write(message, '(1x,a,3f16.12,a)' ) 'length scales=',acell,' bohr'
  call wrtout(ab_out,message,'COLL')
  write(message, '(14x,a,3f16.12,a)' ) '=',Bohr_Ang*acell(1:3),' angstroms'
  call wrtout(ab_out,message,'COLL')

 end if

 option=1 ; nnonsc=0 ; tolwf=0.0_dp
 if(iscf<=0 .and. iscf/=-3)option=3
 iscf_dum=iscf
 if(nstep==0)iscf_dum=0
 if(dtset%tfkinfunc==0)then
  call prteigrs(eigen,enunit,fermie,filnam,ab_out,iscf_dum,kptns,kptopt,&
&  mband,nband,nkpt,nnonsc,nsppol,occ,occopt,option,prteig,prtvol,resid,tolwf,vxcavg,wtk)
 endif

!Compute and print location of maximal and minimal density
!DEBUG
!write(6,*) ' clnup1 : before prtrhomxmn '
!write(6,*) ' nspden,  rhor(1,1), rhor(1,2) ',nspden,rhor(1,1), rhor(1,2)
!ENDDEBUG
 call prtrhomxmn(ab_out,mpi_enreg,nfft,ngfft,nspden,2,rhor)

!If needed, print DOS
 if(prtdos==1)then
! Initialize the file
  fildos=trim(filnam)//'_DOS'
  unitdos=tmp_unit
  option=2
  open (unit=unitdos,file=fildos,status='unknown',form='formatted')
  rewind(unitdos)
  maxocc=two/(nspinor*nsppol)  ! Will not work in the fixed moment case
  allocate(doccde(mband*nkpt*nsppol))
  call getnel(doccde,dosdeltae,eigen,entropy,fermie,maxocc,mband,nband,&
&  nelect,nkpt,nsppol,occ,occopt,option,tphysel,tsmear,unitdos,wtk)
  deallocate(doccde)
 end if

!DEBUG
!write(6,*)' clnup1 : exit,  fred=',fred(1,1)
!stop
!ENDDEBUG

end subroutine clnup1
!!***
