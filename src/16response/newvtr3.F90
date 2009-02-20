!{\src2tex{textfont=tt}}
!!****f* ABINIT/newvtr3
!! NAME
!! newvtr3
!!
!! FUNCTION
!! Call prcref to compute preconditioned residual potential and forces,
!! Then, call one of the self-consistency drivers,
!! then update vtrial.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam.
!!  etotal=the total energy obtained from the input vtrial
!!  character(len=fnlen) :: filfft=name of _FFT file
!!  initialized= if 0 the initialization of the RF run is not yet finished
!!  iscf=( <= 0 =>non-SCF), >0 => SCF)
!!   iscf =2 => SCF cycle, simple mixing
!!   iscf =3 => SCF cycle, anderson mixing
!!  isecur=level of security of the computation
!!  istep= number of the step in the SCF cycle
!!  i_rhor(n_index)=indices of the density in the array f_fftgr
!!  i_vresid(n_index)=indices of the potential residuals in the array f_fftgr
!!  i_vrespc(n_index)=indices of the preconditioned potential residuals in the array f_fftgr
!!  i_vtrial(n_index)=indices of the potential in the array f_fftgr
!!  mffmem=governs the number of FFT arrays which are fit in core memory
!!           it is either 1, in which case the array f_fftgr is used,
!!           or 0, in which case the array f_fftgr_disk is used
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  n_fftgr=third dimension of the array f_fftgr
!!  n_index=dimension for indices of potential/density (see i_vresid, ivrespc, i_vtrial...)
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!                                         Use here rhoij residuals (and gradients)
!!  rhor(cplex*nfft,nspden)=array for 1st-order electron density
!!    in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  vresid(cplex*nfft,nspden)=array for the residual of the potential
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  dbl_nnsclo=1 if nnsclo has to be doubled to secure the convergence.
!!
!! SIDE EFFECTS
!!  f_fftgr(cplex*nfft,nspden,n_fftgr*mffmem)=different functions defined
!!    on the fft grid
!!   (see prcref, scfeig, scfopt, and scfcge for a detailed explanation).
!!   If mffmem=0, these data are kept on disk.
!!  vtrial(cplex*nfft,nspden)= at input, it is the trial potential
!!   that gave vresid .
!!       at output, it is an updated trial potential
!!
!! NOTES
!!
!! PARENTS
!!      scfcv3
!!
!! CHILDREN
!!      leave_new,metric,moddiel,scfcge,scfopt,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine newvtr3(cplex,dbl_nnsclo,dielar,etotal,filfft,f_fftgr,&
&  initialized,iscf,isecur,istep,i_rhor,&
&  i_vresid,i_vrespc,i_vtrial,mffmem,mpi_enreg,natom,nfft,ngfft,nspden,&
&  n_fftgr,n_index,paral_kgb,pawrhoij,qphon,rhor,rprimd,vresid,vtrial,xred)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_15common
!End of the abilint section

 implicit none

!Arguments-------------------------------
!scalars
 integer,intent(in) :: cplex,initialized,iscf,isecur,istep,mffmem,n_fftgr
 integer,intent(in) :: n_index,natom,nfft,nspden,paral_kgb
 integer,intent(out) :: dbl_nnsclo
 real(dp),intent(in) :: etotal
 character(len=fnlen),intent(in) :: filfft
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(inout) :: i_rhor(n_index),i_vresid(n_index),i_vrespc(n_index)
 integer,intent(inout) :: i_vtrial(n_index)
 real(dp),intent(in) :: dielar(7),qphon(3),rhor(cplex*nfft,nspden),rprimd(3,3)
 real(dp),intent(in) :: vresid(cplex*nfft,nspden)
 real(dp),intent(inout) :: vtrial(cplex*nfft,nspden),xred(3,natom)
 real(dp),intent(out) :: f_fftgr(cplex*nfft,nspden,n_fftgr*mffmem)
 type(pawrhoij_type),intent(inout) :: pawrhoij(natom)

!Local variables-------------------------------
!scalars
 integer :: i_vresid1,i_vrespc1,ifft,ispden,moved_atm_inside,nfftot,response
 real(dp) :: diemix,ucvol
 character(len=500) :: message
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),tsec(2)
 real(dp),allocatable :: dtn_pc(:,:),f_atm(:,:,:),f_fftgr_disk(:,:,:)
 real(dp),allocatable :: f_paw(:,:),f_paw_disk(:,:),vpaw(:),vrespc(:,:)

! *************************************************************************

!DEBUG
!write(6,*)' newvtr3 : enter '
!stop
!ENDDEBUG

 if(nspden==4)then
  write(6,*)' newvtr3 : does not work yet for nspden=4'
  stop
 end if

 call timab(158,1,tsec)

 nfftot=ngfft(1)*ngfft(2)*ngfft(3)

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 moved_atm_inside=0

 allocate(vrespc(cplex*nfft,nspden),dtn_pc(3,natom))
 allocate(f_atm(3,natom,n_fftgr))

!Model dielectric function preconditioning, or simple mixing
!Get vrespc from vresid
 call moddiel(cplex,dielar,mpi_enreg,nfft,ngfft,nspden,1,paral_kgb,qphon,&
& rprimd,vresid,vrespc)

!------Compute new vtrial and eventual new atomic positions

 i_vresid1=i_vresid(1)
 i_vrespc1=i_vrespc(1)

!Either use the array f_fftgr or the array f_fftgr_disk
 if(mffmem==1)then
  f_fftgr(:,:,i_vresid1)=vresid(:,:)
  f_fftgr(:,:,i_vrespc1)=vrespc(:,:)
 else
! In this case, must first allocate f_fftgr_disk, then take data
! from disk and existing arrays.
  allocate(f_fftgr_disk(cplex*nfft,nspden,n_fftgr))
  if(istep/=1)then
   call timab(83,1,tsec)
   open(unit=tmp_unit,file=filfft,form='unformatted',status='old')
   rewind(tmp_unit)
   read(tmp_unit)f_fftgr_disk
   close(unit=tmp_unit)
   call timab(83,2,tsec)
  end if
  f_fftgr_disk(:,:,i_vresid1)=vresid(:,:)
  f_fftgr_disk(:,:,i_vrespc1)=vrespc(:,:)

 end if

 deallocate(vrespc)

!Unlike for GS, no need to modify the mean of vtrial

 if(iscf==2 .or. iscf==3 .or. iscf==7)then

! Optimize next vtrial using different algorithms, as
! determined by the variable iscf
  if(mffmem==1)then
   call scfopt(cplex,dtn_pc,f_fftgr,f_paw,iscf,istep,i_vrespc,i_vtrial,&
&   moved_atm_inside,mpi_enreg,natom,nfft,0,nspden,n_fftgr,n_index,0,0,vpaw,vtrial,xred)
  else
   call scfopt(cplex,dtn_pc,f_fftgr_disk,f_paw_disk,iscf,istep,i_vrespc,i_vtrial,&
&   moved_atm_inside,mpi_enreg,natom,nfft,0,nspden,n_fftgr,n_index,0,0,vpaw,vtrial,xred)
  end if

 else if(iscf==5 .or. iscf==6)then

! Optimize next vtrial using an algorithm based
! on the conjugate gradient minimization of etotal
  response=1
  if(mffmem==1)then
   call scfcge(cplex,dbl_nnsclo,dtn_pc,etotal,f_atm,&
&   f_fftgr,initialized,iscf,isecur,istep,&
&   i_rhor,i_vresid,i_vrespc,moved_atm_inside,mpi_enreg,&
&   natom,nfft,nfftot,nspden,n_fftgr,n_index,response,rhor,ucvol,vtrial,xred)
  else
   call scfcge(cplex,dbl_nnsclo,dtn_pc,etotal,f_atm,&
&   f_fftgr_disk,initialized,iscf,isecur,istep,&
&   i_rhor,i_vresid,i_vrespc,moved_atm_inside,mpi_enreg,&
&   natom,nfft,nfftot,nspden,n_fftgr,n_index,response,rhor,ucvol,vtrial,xred)
  end if

 else

  write(message, '(a,a,a,a,i5,a)' ) ch10,&
&  ' newvtr3 : BUG -',ch10,&
&  '  Invalid option: iscf =',iscf,'.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')

 end if

!Eventually write the data on disk and deallocate f_fftgr_disk
 if(mffmem==0)then
  call timab(83,1,tsec)
  open(unit=tmp_unit,file=filfft,form='unformatted',status='unknown')
  rewind(tmp_unit)
  write(tmp_unit)f_fftgr_disk
  close(unit=tmp_unit)
  deallocate(f_fftgr_disk)
  call timab(83,2,tsec)
 end if

 deallocate(dtn_pc,f_atm)

 call timab(158,2,tsec)

!DEBUG
!write(6,*)' newvtr3 : exit '
!write(6,*)' newvtr3 : vtrial(1,1)=',vtrial(1,1)
!stop
!ENDDEBUG

end subroutine newvtr3
!!***
