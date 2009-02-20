!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_epsm1_driver
!! NAME
!! make_epsm1_driver
!!
!! FUNCTION
!!  Driver routine to calculate the symmetrical inverse dielectric matrix starting
!!  from the irreducible polarizability. This routines consider a single q-point, and 
!!  performs the following tasks:
!!
!!  1) Calculate epsilon^-1 using different approximations
!!      * RPA
!!      * ALDA within TDDFT
!!
!!  2) Use a special treatment of non-Analytic behavior of heads and wings in reciprocal space
!!     calculating these quantities for different small q-directions specified by the user
!!     (Not yet operative)
!!
!!  3) Output Electron energy loss function and Macroscopic dielectric function with and 
!!     without nonlocal field effect (only if the frequency dependency along the real axis has 
!!     been calculated
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  iq=index of the q-point in the array qibz where epsilon^-1 has to be calculated
!!  npwe=Number of G-vectors in chi0.
!!  nqibz=Number of q-points.
!!  qibz(3,nqibz)=q-points in the IBZ.
!!  nI,nJ=Number of rows/columns in chi0_ij (1,1 in collinear case)
!!  nomega=Number of frequencies.
!!  nomegaer=Number of real frequencies.
!!  approx_type=Integer flag defining the type of approximation
!!   == 0 for RPA   ==
!!   == 1 for TDDFT ==
!!  option_test=Only for TDDFT:
!!   == 0 for TESTPARTICLE ==
!!   == 1 for TESTELECTRON ==
!!  Vcp<Coulombian_type>=Structure gathering data on the Coulombian interaction
!!  Dtfil<Datafiles_type)>=variables related to files
!!  FIXME treatment of kxc has to be rewritten.
!!  kxc(npwe*approx_type,npwe*approx_type)=TDDFT kernel, required only if approx_type==1
!!  MPI_enreg=MPI-parallelisation information
!!
!! OUTPUT
!!  Different files are written according to the type of calculation
!!  See also side effects
!!
!! SIDE EFFECTS
!!  chi0(npwe*nI,npwe*nJ,nomega): in input the irreducible polarizability, in output 
!!   the symmetrized inverse dielectric matrix.
!!
!! NOTES
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

subroutine make_epsm1_driver(iq,npwe,nI,nJ,nomega,nomegaer,omega,&
& approx_type,option_test,nqibz,qibz,Vcp,Dtfil,gmet,kxc,MPI_enreg,chi0)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : GW_TOLQ0
 use m_errors, only : assert_eq, assert
 use m_numeric_tools, only : is_zero,print_arr,hermitianize
 use m_io_tools, only : flush_unit, get_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq,nI,nJ,npwe,nomega,nomegaer,nqibz,approx_type,option_test
 type(Coulombian_type),intent(in) :: Vcp
 type(Datafiles_type),intent(in) :: Dtfil
 type(MPI_type),intent(in) :: MPI_enreg
!arrays
 real(dp),intent(in) :: gmet(3,3),qibz(3,nqibz)
 complex(gwpc),intent(in) :: kxc(npwe*approx_type,npwe*approx_type) !FIXME this has to be rewritten LDA is (npwe,1)
 complex(gwpc),intent(in) :: omega(nomega)
 complex(gwpc),intent(inout) :: chi0(npwe*nI,npwe*nJ,nomega)

!Local variables-------------------------------
 character(len=50),parameter :: FILE__='make_epsm1_driver.F90'
!scalars
 integer :: ig1,ig2,io,istat,master,rank,skxc,unt
 real(dp) :: epsilon0,epsilon0_nlf
 logical :: qeq0
 character(len=500) :: msg
 character(len=fnlen) :: fnam
!arrays
 real(dp) :: tsec(2)
 complex(gwpc),allocatable :: chitmp(:,:),epsm_lf(:),epsm_nlf(:)
 complex(gwpc),pointer :: vc_sqrt(:)

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' make_epsm1_driver : enter '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 call timab(309,1,tsec) ! chi/eps
 call xmaster_init(MPI_enreg,master) 
 call xme_init    (MPI_enreg,rank  )          

if (rank==master) then ! presently only master has chi0 in screening

 allocate(epsm_lf(nomega),epsm_nlf(nomega)) 

 ! === vc_sqrt contains vc^{1/2}(q,G) ====
 ! * complex-valued to allow for a possible cutoff
 vc_sqrt => Vcp%vc_sqrt(:,iq)  
 qeq0=(normv(qibz(:,iq),gmet,'G')<GW_TOLQ0)
 if (qeq0) then 
  !write(*,*)' Analyzing long wavelength limit '
  !vc_sqrt => Vcp%vcqs_sqrt(:,1) !for the moment only first point
 end if

 if (nI/=1.or.nJ/=1) then 
  call assert(.FALSE.,'nI or nJ=/1 not yet implemented')
 end if

 select case (approx_type)
  case (0)
   ! === RPA: \tilde\epsilon=(1-Vc^{1/2}*chi0Vc^{1/2}) ===
   allocate(chitmp(npwe,npwe),STAT=istat) 
   if (istat/=0) call memerr(FILE__,'chitmp',npwe**2,'gwpc')
   do io=1,nomega
    do ig2=1,npwe
     chitmp(:,ig2)=-vc_sqrt(:)*chi0(:,ig2,io)*vc_sqrt(ig2)
     chitmp(ig2,ig2)=one+chitmp(ig2,ig2)
    end do
    ! === chi0(io), now contains epsilon(io) ===
    chi0(:,:,io)=chitmp(:,:)
    epsm_nlf(io)=chitmp(1,1)
    write(msg,'(a,i4,a,2f9.4,a)')&
 &   ' Symmetrical epsilon(G,G'') at the ',io,' th omega',omega(io)*Ha_eV,' [eV]'
    call wrtout(std_out,msg,'COLL')
    call print_arr(chi0(:,:,io))
   end do
   deallocate(chitmp)
   !
   ! === Invert epsilon tilde and calculate macroscopic dielectric constant ===
   ! * epsm_lf(w)=1/epsm1(G=0,Gp=0,w). 
   ! * Since G=Gp=0 there is no difference btw symmetrical and not symmetrical 
   do io=1,nomega
    call matcginv(chi0(:,:,io),npwe,npwe)
    epsm_lf(io)=one/chi0(1,1,io)
    write(msg,'(a,i4,a,2f9.4,a)')&
 &   ' Symmetrical epsilon^-1(G,G'') at the ',io,' th omega',omega(io)*Ha_eV,' [eV]'
    call wrtout(std_out,msg,'COLL')
    call print_arr(chi0(:,:,io))
   end do

  case (1)
   ! === Vertex correction from Adiabatic TDDFT ===
   skxc=assert_eq(SIZE(kxc,1),SIZE(kxc,2),'kxc not square',__FILE__,__LINE__) 
   if (skxc/=npwe) STOP 'wrong size in kxc'
   allocate(chitmp(npwe,npwe),STAT=istat) ; if (istat/=0) call memerr(FILE__,'chitmp',npwe**2,'spc')
   do io=1,nomega
    ! === Calculate chi0*fxc ===
    chitmp(:,:)=MATMUL(chi0(:,:,io),kxc(:,:))
    ! === Calculate (1-chi0*Vc-chi0*Kxc) and put it in chitmp ===
    do ig1=1,npwe
     do ig2=1,npwe
      chitmp(ig1,ig2)=-chitmp(ig1,ig2)-vc_sqrt(ig2)**2*chi0(ig1,ig2,io) 
     end do
     chitmp(ig1,ig1)=chitmp(ig1,ig1)+one
    end do
    ! === Invert (1-chi0*Vc-chi0*Kxc) and Multiply by chi0 ===
    call matcginv(chitmp,npwe,npwe)
    chitmp=MATMUL(chitmp,chi0(:,:,io))
    ! === Save result, now chi0 contains chi ===
    chi0(:,:,io)=chitmp(:,:)
    write(std_out,'(a,i2,a,i1,a)')' chi(q= ',iq,',omega= ',io,',G,G")'
    call print_arr(chi0(:,:,io),mode_paral='PERS')
   end do

   select case (option_test)
    case (0) 
     ! === Calculate symmetrized TESTPARTICLE epsilon^-1 ===
     write(msg,'(a)')' calculating TESTPARTICLE epsilon^-1(G,G") = 1 + Vc*chi'
     call wrtout(std_out,msg,'COLL')
     do io=1,nomega
      do ig1=1,npwe
       chi0(ig1,:,io)=(vc_sqrt(ig1)*vc_sqrt(:))*chi0(ig1,:,io)
       chi0(ig1,ig1,io)=one+chi0(ig1,ig1,io)
      end do 
     end do 
    case (1)
     ! === Calculate symmetrized TESTELECTRON epsilon^-1 ===
     write(msg,'(a)')' calculating TESTELECTRON epsilon^-1(G,G") = 1 + (Vc + fxc)*chi'
     call wrtout(std_out,msg,'COLL')
     do io=1,nomega
      chitmp=MATMUL(kxc(:,:),chi0(:,:,io))
      ! === Perform hermitianization (why ?) ===
      call hermitianize(chitmp)
      do ig1=1,npwe
       chi0(ig1,:,io)=(vc_sqrt(ig1)*vc_sqrt(:))*chi0(ig1,:,io)+chitmp(ig1,:)
       chi0(ig1,ig1,io)=one+chi0(ig1,ig1,io)
      end do 
     end do
    case default
     write(msg,'(4a,i3)')ch10,&
      ' make_epsm1_driver : BUG - ',ch10,&
      ' wrong value for option_test = ',option_test
     call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
   end select
   deallocate(chitmp)
   !
   ! === chi0 now contains symmetrical epsm1 ===
   ! * Calculate macroscopic dielectric constant epsm_lf(w)=1/epsm1(G=0,Gp=0,w) ===
   epsm_lf(:)=one/chi0(1,1,:)
   do io=1,nomega
    !write (msg,'(a,i2,a,i1,a)')' Symmetrical epsilon^-1(q=',iq,',omega=',io,',G,G")'
    !call wrtout(std_out,msg,'COLL')
    call print_arr(chi0(:,:,io),mode_paral='PERS')
   end do

 case default
  write(msg,'(4a,i3)')ch10,&
   ' make_epsm1_driver : BUG - ',ch10,&
   ' wrong value for approx_type = ',approx_type 
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end select

end if !master
 call timab(309,2,tsec) ! chi/eps

 ! =============================================
 ! === Now chi0 contains \tilde\epsilon^{-1} ===
 ! =============================================
 ! * Master node writes the following data:
 !   1) ELF 
 !   2) Macroscopic dielectric constant with and without local fields, 
 if (rank==master) then 
  if (qeq0.and.nomega>2) then
   unt=get_unit() ; fnam=TRIM(Dtfil%filnam_ds(4))//'_ELF' ; open(unit=unt,file=fnam)
   write(unt,'(a)'      )'# Electron Energy Loss Function'
   write(unt,'(a,3f9.6)')'# q = ',qibz(:,iq)
   write(unt,'(a)'      )'# Omega [eV]    ELF -Im(1/epsilon_M)'
   do io=1,nomegaer
    write(unt,'(1x,f7.3,7x,e12.4)')REAL(omega(io))*Ha_eV,-AIMAG(chi0(1,1,io))
   end do
   close(unt)
   fnam=TRIM(Dtfil%filnam_ds(4))//'_EM1_LF' ; open(unit=unt,file=fnam)
   write(unt,'(a)'      )'# Macroscopic Dielectric Function with local fields included'
   write(unt,'(a,3f9.6)')'# q = ',qibz(:,iq)
   write(unt,'(a)'      )'# Omega [eV]    Re epsilon_M       Im eps_M '
   do io=1,nomegaer
    write(unt,'(1x,f7.3,7x,2(e12.4,7x))')REAL(omega(io))*Ha_eV,epsm_lf(io)
   end do
   close(unt)
   if (approx_type==0) then 
    fnam=TRIM(dtfil%filnam_ds(4))//'_EM1_NLF' ; open(unit=unt,file=fnam)
    write(unt,'(a)'      )'# Macroscopic Dielectric Function without local fields'
    write(unt,'(a,3f9.6)')'# q = ',qibz(:,iq)
    write(unt,'(a)'      )'# Omega [eV]    Re epsilon_M       IM eps_M '
    do io=1,nomegaer
     write(unt,'(1x,f7.3,7x,2(e12.4,7x))')REAL(omega(io))*Ha_eV,epsm_nlf(io)
    end do
    close(unt)
   end if 
  end if
  if (qeq0) then
   epsilon0=REAL(epsm_lf(1))
   do io=1,nomega
    if (ABS(REAL(omega(io)))<1.e-3.and.ABS(AIMAG(omega(io)))<1.e-3) then
     write(msg,'(1x,a,f8.4)')' dielectric constant = ',epsilon0
     call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
     if (approx_type==0) then 
      epsilon0_nlf=REAL(epsm_nlf(1))
      write(msg,'(1x,a,f8.4,a)')' dielectric constant without local fields = ',epsilon0_nlf,ch10
      call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
     end if 
    end if
   end do 
  end if
 deallocate(epsm_lf,epsm_nlf)
end if !master

#if defined DEBUG_MODE
 write(msg,'(a)')' make_epsm1_driver : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine make_epsm1_driver
!!***
