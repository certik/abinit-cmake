!{\src2tex{textfont=tt}}
!!****f* ABINIT/transgrid
!! NAME
!! transgrid
!!
!! FUNCTION
!! Convert a given density (or potential) from the coarse
!! to the fine rectangular grid  -  and vice versa
!! Used in PAW calculations
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cplex=1 if rhor[f] is real, 2 if rhor[f] is complex
!!  mpi_enreg=informations about MPI parallelization
!!  nspden=number of spin-density components
!!  optgrid=+1 to go from the coarse grid towards the fine grid
!!          -1 to go from the fine grid towards the coarse grid
!!  optin= 0: input density/potential is taken from rhor(:,nspden)
!!         1: input density/potential is taken from rhog(:)     (ispden=1)
!!                                              and rhor(:,2:4) (ispden=2,3,4)
!!  optout= 0: output density/potential is given in r space in rhor(:,nspden)
!!          1: output density/potential is given in r space in rhor(:,nspden)
!!                                           and in g space in rhog(:)
!!  pawfgr <type(paw_fgr_type)>=fine rectangular grid parameters
!!    %nfftc=number of points in the coarse FFT box
!!    %nfft =number of points in the fine FFT box
!!    %ngfftc(18)=all needed information about 3D FFT, for the coarse grid
!!    %ngfft(18) =all needed information about 3D FFT, for the fine grid
!!    %coatofin(nfftc)=Index of the points of the coarse grid on the fine grid
!!    %fintocoa(nfft) =Index of the points of the fine grid on the coarse grid
!!    %usefinegrid= 1 if a fine FFT grid is used (0 otherwise)
!!  if optgrid=+1 and optin=1:
!!    rhog(2,nfftc)=Fourier transform of input density/potential on the coarse grid
!!  if optgrid=-1 and optin=1:
!!    rhogf(2,nfftf)=Fourier transform of input density/potential on the fine grid
!!  if optgrid=+1
!!    rhor(cplex*nfftc,nspden)=input density/potential in r space on the coarse grid
!!  if optgrid=-1:
!!    rhorf(cplex*nfftf,nspden)=input density/potential in r space on the fine grid
!!
!! OUTPUT
!!  if optgrid=-1 and optout=1:
!!    rhog(2,nfftc)=Fourier transform of output density/potential on the coarse grid
!!  if optgrid=+1 and optout=1:
!!    rhogf(2,nfftf)=Fourier transform of output density/potential on the fine grid
!!  if optgrid=-1
!!    rhor(cplex*nfftc,nspden)=output density/potential in r space on the coarse grid
!!  if optgrid=+1:
!!    rhorf(cplex*nfftf,nspden)=output density/potential in r space on the fine grid
!!
!! PARENTS
!!      energy,gstate,vtorho
!!
!! CHILDREN
!!      fourdp,indirect_parallel_fourier,leave_new,wrtout,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine transgrid(cplex,mpi_enreg,nspden,optgrid,optin,optout,paral_kgb,pawfgr,rhog,rhogf,rhor,rhorf)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_lib01fftnew
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,nspden,optgrid,optin,optout,paral_kgb
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawfgr_type),intent(in) :: pawfgr
!arrays
 real(dp),intent(inout) :: rhog(2,pawfgr%nfftc),rhogf(2,pawfgr%nfft)
 real(dp),intent(inout) :: rhor(cplex*pawfgr%nfftc,nspden),rhorf(cplex*pawfgr%nfft,nspden)

!Local variables ---------------------------------------
!scalars
 integer :: i1,ispden,nfftc,nfftctot,nfftf,nfftftot
 character(len=500) :: message
!arrays
 integer :: ngfftc(18),ngfftf(18)
 real(dp),allocatable :: vectg(:,:),work(:,:),workfft(:)

! *************************************************************************

!Tests
 if(pawfgr%nfft<pawfgr%nfftc) then
  write(message, '(4a)' )ch10,&
&  ' transgrid : ERROR -',ch10,&
&  '  nfft (fine grid) must be >= nfft (coarse grid).'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Store FFT dimensions
 nfftc=pawfgr%nfftc;ngfftc(:)=pawfgr%ngfftc(:);nfftctot=ngfftc(1)*ngfftc(2)*ngfftc(3)
 nfftf=pawfgr%nfft ;ngfftf(:)=pawfgr%ngfft (:);nfftftot=ngfftf(1)*ngfftf(2)*ngfftf(3)

!If no fine FFT grid is used, this is only a simple transfer
 if (pawfgr%usefinegrid==0) then
  if (optgrid==1) then
   rhorf=rhor
   if (optout==1.and.optin==1) rhogf=rhog
   if (optout==1.and.optin/=1) then
    allocate(workfft(cplex*nfftc));workfft(:)=rhor(:,1)
    call fourdp(cplex,rhogf,workfft,-1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
    deallocate(workfft)
   end if
  end if
  if (optgrid==-1) then
   rhor=rhorf
   if (optout==1.and.optin==1) rhog=rhogf
   if (optout==1.and.optin/=1) then
    allocate(workfft(cplex*nfftc));workfft(:)=rhorf(:,1)
    call fourdp(cplex,rhog,workfft,-1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
    deallocate(workfft)
   end if
  end if
  return
 end if

!====== FROM THE COARSE GRID TOWARDS THE FINE GRID =============
!===============================================================
!Calculate the FT of rhor to have it in the g space on the coarse grid
!Transfer the FT of rhor on the coarse grid towards the fine grid
!Then calculate the FT back to get rhorf on the fine grid
 if (optgrid==1) then

  allocate(work(2,nfftc))

! First spin component
! --------------------------------------------------------------
  if (optout==0) then
!  if optout=0, rhog on the fine grid is temporary (in vectg)
   allocate(vectg(2,nfftf));vectg(:,:)=zero
   if (optin==1) then
    call zerosym(rhog,2,mpi_enreg,ngfftc(1),ngfftc(2),ngfftc(3))
    if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_compil_fft==1) then
     call indirect_parallel_Fourier&
&     (pawfgr%coatofin,vectg,mpi_enreg,ngfftf,ngfftc,nfftf,nfftc,paral_kgb,rhog,nfftctot)
    else
     do i1=1,nfftc
      vectg(:,pawfgr%coatofin(i1))=rhog(:,i1)
     end do
    end if
   else
    allocate(workfft(cplex*nfftc));workfft(:)=rhor(:,1)
    call fourdp(cplex,work,workfft,-1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
    deallocate(workfft)
    call zerosym(work,2,mpi_enreg,ngfftc(1),ngfftc(2),ngfftc(3))
    if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_compil_fft==1) then
     call indirect_parallel_Fourier&
&     (pawfgr%coatofin,vectg,mpi_enreg,ngfftf,ngfftc,nfftf,nfftc,paral_kgb,work,nfftctot)
    else
     do i1=1,nfftc
      vectg(:,pawfgr%coatofin(i1))=work(:,i1)
     end do
    end if
   end if
!  call zerosym(vectg,2,mpi_enreg,ngfftf(1),ngfftf(2),ngfftf(3))
   allocate(workfft(cplex*nfftf))
   call fourdp(cplex,vectg,workfft,1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)
   rhorf(:,1)=workfft(:);deallocate(workfft)
   deallocate(vectg)
  else
!  if optout=1, rhog on the fine grid is saved
   call zerosym(rhog,2,mpi_enreg,ngfftc(1),ngfftc(2),ngfftc(3))
   rhogf(:,:)=zero
   if (optin==1) then
    if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_compil_fft==1) then
     call indirect_parallel_Fourier&
&     (pawfgr%coatofin,rhogf,mpi_enreg,ngfftf,ngfftc,nfftf,nfftc,paral_kgb,rhog,nfftctot)
    else
     do i1=1,nfftc
      rhogf(:,pawfgr%coatofin(i1))=rhog(:,i1)
     end do
    end if
   else
    allocate(workfft(cplex*nfftc));workfft(:)=rhor(:,1)
    call fourdp(cplex,work,workfft,-1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
    deallocate(workfft)
    call zerosym(work,2,mpi_enreg,ngfftc(1),ngfftc(2),ngfftc(3))
    if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_compil_fft==1) then
     call indirect_parallel_Fourier&
&     (pawfgr%coatofin,rhogf,mpi_enreg,ngfftf,ngfftc,nfftf,nfftc,paral_kgb,work,nfftctot)
    else
     do i1=1,nfftc
      rhogf(:,pawfgr%coatofin(i1))=work(:,i1)
     end do
    end if
   end if
!  call zerosym(rhogf,2,mpi_enreg,ngfftf(1),ngfftf(2),ngfftf(3))
   allocate(workfft(cplex*nfftf))
   call fourdp(cplex,rhogf,workfft,1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)
   rhorf(:,1)=workfft(:);deallocate(workfft)
  end if

! Additional spin components
! ----------------------------------------------------
  if (nspden>=2) then
   allocate(vectg(2,nfftf))
   do ispden=2,nspden
    vectg(:,:)=zero
    allocate(workfft(cplex*nfftc));workfft(:)=rhor(:,ispden)
    call fourdp(cplex,work,workfft,-1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
    deallocate(workfft)
    call zerosym(work,2,mpi_enreg,ngfftc(1),ngfftc(2),ngfftc(3))
    if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_compil_fft==1) then
     call indirect_parallel_Fourier&
&     (pawfgr%coatofin,vectg,mpi_enreg,ngfftf,ngfftc,nfftf,nfftc,paral_kgb,work,nfftctot)
    else
     do i1=1,nfftc
      vectg(:,pawfgr%coatofin(i1))=work(:,i1)
     end do
    end if
!   call zerosym(vectg,2,mpi_enreg,ngfftf(1),ngfftf(2),ngfftf(3))
    allocate(workfft(cplex*nfftf))
    call fourdp(cplex,vectg,workfft,1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)
    rhorf(:,ispden)=workfft(:);deallocate(workfft)
   end do
   deallocate(vectg)
  end if

  deallocate(work)


! ====== FROM THE FINE GRID TOWARDS THE COARSE GRID =============
! ==============================================================
! Calculate the FT of rhorf to have it in the g space on the fine grid
! Transfer the FT of rhorf on the fine grid towards the coarse grid
! Then calculate the FT back to get rhor on the coarse grid
 else if (optgrid==-1) then

  allocate(work(2,nfftf))

! First spin component
! --------------------------------------------------------------
  if (optout==0) then
!  if optout=0, rhog on the fine grid is temporary (in vectg)
   allocate(vectg(2,nfftc));vectg(:,:)=zero
   if (optin==1) then
    do i1=1,nfftf
     if (pawfgr%fintocoa(i1)/=0) vectg(:,pawfgr%fintocoa(i1))=rhogf(:,i1)
    end do
   else
    allocate(workfft(cplex*nfftf));workfft(:)=rhorf(:,1)
    call fourdp(cplex,work,workfft,-1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)
    deallocate(workfft)
    if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_compil_fft==1) then
     call indirect_parallel_Fourier&
&     (pawfgr%fintocoa,vectg,mpi_enreg,ngfftc,ngfftf,nfftc,nfftf,paral_kgb,work,nfftftot)
    else
     do i1=1,nfftf
      if (pawfgr%fintocoa(i1)/=0) vectg(:,pawfgr%fintocoa(i1))=work(:,i1)
     end do
    end if
   end if
   call zerosym(vectg,2,mpi_enreg,ngfftc(1),ngfftc(2),ngfftc(3))
   allocate(workfft(cplex*nfftc))
   call fourdp(cplex,vectg,workfft,1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
   rhor(:,1)=workfft(:);deallocate(workfft)
   deallocate(vectg)
  else
!  if optout=1, rhog on the fine grid is saved
   rhog(:,:)=zero
   if (optin==1) then
    do i1=1,nfftf
     if (pawfgr%fintocoa(i1)/=0) rhog(:,pawfgr%fintocoa(i1))=rhogf(:,i1)
    end do
   else
    allocate(workfft(cplex*nfftf));workfft(:)=rhorf(:,1)
    call fourdp(cplex,work,workfft,-1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)
    deallocate(workfft)
    if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_compil_fft==1) then
     call indirect_parallel_Fourier&
&     (pawfgr%fintocoa,rhog,mpi_enreg,ngfftc,ngfftf,nfftc,nfftf,paral_kgb,work,nfftftot)
    else
     do i1=1,nfftf
      if (pawfgr%fintocoa(i1)/=0) rhog(:,pawfgr%fintocoa(i1))=work(:,i1)
     end do
    end if
   end if
   call zerosym(rhog,2,mpi_enreg,ngfftc(1),ngfftc(2),ngfftc(3))
   allocate(workfft(cplex*nfftc))
   call fourdp(cplex,rhog,workfft,1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
   rhor(:,1)=workfft(:);deallocate(workfft)
  end if

! Additional spin components
! ----------------------------------------------------
  if (nspden>=2) then
   allocate(vectg(2,nfftc))
   do ispden=2,nspden
    vectg(:,:)=zero
    allocate(workfft(cplex*nfftf));workfft(:)=rhorf(:,ispden)
    call fourdp(cplex,work,workfft,-1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)
    deallocate(workfft)
    if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_compil_fft==1) then
     call indirect_parallel_Fourier&
&     (pawfgr%fintocoa,vectg,mpi_enreg,ngfftc,ngfftf,nfftc,nfftf,paral_kgb,work,nfftftot)
    else
     do i1=1,nfftf
      if (pawfgr%fintocoa(i1)/=0) vectg(:,pawfgr%fintocoa(i1))=work(:,i1)
     end do
    end if
    call zerosym(vectg,2,mpi_enreg,ngfftc(1),ngfftc(2),ngfftc(3))
    allocate(workfft(cplex*nfftc))
    call fourdp(cplex,vectg,workfft,1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
    rhor(:,ispden)=workfft(:);deallocate(workfft)
   end do
   deallocate(vectg)
  end if

  deallocate(work)

 end if

end subroutine transgrid
!!***
