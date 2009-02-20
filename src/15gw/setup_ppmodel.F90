!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_ppmodel
!! NAME
!! setup_ppmodel
!!
!! FUNCTION
!!  Initialize some values of several arrays of the Er% datastructure 
!!  that are used in case of plasmonpole calculations
!!  Just a wrapper around different plasmonpole routines.
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  paral_kgb=variable related to band parallelism
!!  Qmesh<bz_mesh_type>=the q-mesh used for the inverse dielectric matrix
!!    %nibz=number of irreducible q-points
!!    %ibz(3,%nibz)=the irred q-point
!!  Sp<sigma_parameters>=parameters defining the self-energy calculation TODO to be removed
!!    %npwc=number of G vectors for the correlation part
!!    %ppmodel=the type of  plasmonpole model 
!!  Er<epsilonm1_results>=the inverse dielectric matrix 
!!    %nomega=number of frequencies in $\epsilon^{-1}$
!!    %epsm1=the inverse dielctric matrix 
!!    %omega=frequencies in epsm1
!!  MPI_enreg<MPI_type>=informations about MPI parallelization
!!  ngfftf(18)=contain all needed information about the 3D fine FFT mesh, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  nfftf=the number of points in the FFT mesh (for this processor)
!!  rhor_tot(nfftf)=the total charge in real space
!!
!! OUTPUT
!!  
!!
!! SIDE EFFECTS
!!  PPm<PPmodel_type>: 
!!  == if ppmodel 1 or 2 ==
!!   %omegatw and %bigomegatwsq=PPmodel parameters 
!!  == if ppmodel 3 ==
!!   %omegatw, %bigomegatwsq and %eigpot=PPmodel parameters
!!  == if ppmodel 4 ==
!!   %omegatw and %bigomegatwsq=PPmodel parameters 
!!
!! NOTES
!! TODO: rhor_tot should be replaced by rhog_tot
!! FFT parallelism won"t work 
!! Solve Issue with MPI_enreg
!!
!! PARENTS
!!  
!!
!! CHILDREN
!!  
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine setup_ppmodel(PPm,paral_kgb,Qmesh,Sp,Er,MPI_enreg,&
& nfftf,gvec,ngfftf,gmet,gprimd,rhor_tot,&
& epsm1q,iqiA) !Optional

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => setup_ppmodel
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,paral_kgb
 integer,intent(in),optional :: iqiA
 type(BZ_mesh_type),intent(in) :: Qmesh
 type(Epsilonm1_results),intent(inout) :: Er
 type(MPI_type),intent(inout) :: MPI_enreg
 type(PPmodel_type),intent(inout) :: PPm
 type(Sigma_parameters),intent(in) :: Sp
!arrays
 integer,intent(in) :: gvec(3,Sp%npwc),ngfftf(18)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3)
 real(dp),intent(inout) :: rhor_tot(nfftf)
 complex(gwpc),intent(in),optional :: epsm1q(Sp%npwc,Sp%npwc,Er%nomega)

!Local variables-------------------------------
 character(len=50),parameter :: sub_name='setup_ppmodel.F90'
!scalars
 integer :: istat,master,npwc2,npwc3,nqiA,ppmsize,rank,spaceComm
 logical :: ltest,single_q
 character(len=500) :: msg

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' setup_ppmodel : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 call xcomm_init  (MPI_enreg,spaceComm)  
 call xme_init    (MPI_enreg,rank     )          
 call xmaster_init(MPI_enreg,master   )  
 !
 ! === if iqiA is present, then consider only one qpoint to save memory ===
 ! * This means the object has been already initialized
 nqiA=Qmesh%nibz ; single_q=.FALSE.
 if (PRESENT(epsm1q)) then 
  nqiA=1 ; single_q=.TRUE.
  ltest=PRESENT(iqiA)
  call assert(ltest,'For single q-point mode, also iqiA must be present')
 end if
 !
 ! Allocate plasmonpole parameters 
 ! TODO ppmodel==1 by default, should be set to 0 if AC and CD
 SELECT CASE (Sp%ppmodel)

 CASE (0)
  write(msg,'(a)')' Skipping Plasmompole model calculation' 
  call wrtout(std_out,msg,'COLL') ; RETURN

 CASE (1) 
  ! === Godby-Needs, q-dependency enters only through epsilon^-1 ===
  ppmsize=Sp%npwc**2*nqiA+Sp%npwc**2*nqiA
  if (.not.single_q) then
   call cppm1par(Sp%npwc,nqiA,Er%nomega,Er%epsm1,Er%omega,PPm%bigomegatwsq,PPm%omegatw,PPm%drude_plsmf)
  else
   call cppm1par(Sp%npwc,nqiA,Er%nomega,epsm1q,  Er%omega,PPm%bigomegatwsq,PPm%omegatw,PPm%drude_plsmf)
  end if

 CASE (2)
  ! === Hybertsen-Louie ===
  if (.not.single_q) then
   call cppm2par(paral_kgb,Sp%npwc,nqiA,Er%nomega,Er%epsm1,PPm%bigomegatwsq,PPm%omegatw,&
&   ngfftf,gvec,gprimd,rhor_tot,nfftf,Qmesh,gmet)
  else
   call cppm2par(paral_kgb,Sp%npwc,nqiA,Er%nomega,epsm1q,PPm%bigomegatwsq,PPm%omegatw,&
&   ngfftf,gvec,gprimd,rhor_tot,nfftf,Qmesh,gmet,iqiA)
  end if

 CASE (3)
  ! === von Linden-Horsh model ===
  ! TODO Check better double precision, this routine is in a messy state
  if (.not.single_q) then
   call cppm3par(paral_kgb,Sp%npwc,nqiA,Er%nomega,Er%epsm1,PPm%bigomegatwsq,&
&   PPm%omegatw,ngfftf,gvec,gprimd,rhor_tot,nfftf,PPm%eigpot,Qmesh)
  else
   call cppm3par(paral_kgb,Sp%npwc,nqiA,Er%nomega,epsm1q,PPm%bigomegatwsq,&
&   PPm%omegatw,ngfftf,gvec,gprimd,rhor_tot,nfftf,PPm%eigpot,Qmesh,iqiA)
  end if

 CASE (4)
  ! === Engel Farid ===
  ! TODO Check better double precision, this routine is in a messy state
  if (.not.single_q) then
   call cppm4par(paral_kgb,Sp%npwc,nqiA,Er%epsm1,Er%nomega,PPm%bigomegatwsq,&
&   PPm%omegatw,ngfftf,gvec,gprimd,rhor_tot,nfftf,Qmesh)
  else
   call cppm4par(paral_kgb,Sp%npwc,nqiA,epsm1q,Er%nomega,PPm%bigomegatwsq,&
&   PPm%omegatw,ngfftf,gvec,gprimd,rhor_tot,nfftf,Qmesh,iqiA)
  end if

 CASE DEFAULT
  write(msg,'(6a,i6)')ch10,&
&  ' setup_ppmodel: BUG -',ch10,&
&  '  The argument ppmodel should be 1 or 2,',ch10,&
&  '  however, ppmodel=',Sp%ppmodel
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 END SELECT

#if defined DEBUG_MODE
 write(msg,'(a)')' setup_ppmodel : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine setup_ppmodel
!!***

