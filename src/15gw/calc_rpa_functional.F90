!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_rpa_functional
!! NAME
!! calc_rpa_functional
!!
!! FUNCTION
!!  Routine used to calculate the RPA approximation to the correlation energy
!!  from the irreducible polarizability. 
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (FB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  iq=index of the q-point in the array Qmesh%ibz where epsilon^-1 has to be calculated
!!  Ep<Epsilonm1_parameters>=Structure with parameters and dimensions related to the inverse dielectric matrix.
!!  Pvc<Coulombian_type>=Structure gathering data on the Coulombian interaction
!!  Qmesh<BZ_mesh_type>=Data type with information on the q-sampling
!!  Dtfil<Datafiles_type)>=variables related to files
!!  kxc(Ep%npwe,Ep%npwe)=TDDFT kernel, only if Ep%tddft is .TRUE.
!!  MPI_enreg=MPI-parallelisation information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
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

subroutine calc_rpa_functional(iq,Ep,Pvc,Qmesh,Dtfil,gmet,kxc,MPI_enreg,chi0)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : GW_TOLQ0
 use m_errors, only : assert_eq
 use m_numeric_tools, only : is_zero,print_arr,hermitianize
 use m_IO_tools, only : flush_unit, get_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_lib00numeric
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq
 type(BZ_mesh_type),intent(in) :: Qmesh
 type(Coulombian_type),intent(in) :: Pvc
 type(Datafiles_type),intent(in) :: Dtfil
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 real(dp),intent(in) :: gmet(3,3)
 complex(gwpc),intent(in) :: kxc(:,:)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)

!Local variables-------------------------------
 character(len=50),parameter :: FILE__='calc_rpa_functional.F90'
!NEW
!scalars
 integer,parameter :: nlambda=8
 integer :: ig1,ig2,ilambda,io,istat,kinde,master,rank,skxc,unt
 real(dp),save :: ecorr=0.0_dp,ecorr_diag=0.0_dp,ecorr_l(nlambda)=0.0_dp
 real(dp) :: epsilon0,epsilon0_nlf,lambda
 logical :: qeq0
 character(len=500) :: msg
 character(len=fnlen) :: fnam
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: eig(:),rwork(:),z(:),zl(:),zlw(:),zw(:)
 complex(dpc),allocatable :: mtmp(:,:),work(:)
 complex(gwpc),allocatable :: chi0_diag(:),chitmp(:,:),epsm_lf(:),epsm_nlf(:)
 complex(gwpc),pointer :: vc_sqrt(:)

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' calc_rpa_functional : enter '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif
 write(*,'(a)')' calculate_RPA_functional : enter '

 call timab(309,1,tsec) ! chi/eps
 call xmaster_init(MPI_enreg,master) 
 call xme_init    (MPI_enreg,rank  )          

if (rank==master) then ! presently only master has chi0 in screening

 allocate(epsm_lf(Ep%nomega),epsm_nlf(Ep%nomega)) 
 ! vc_sqrt contains vc^{1/2}(q,G), complex-valued to allow for a possible cutoff
 vc_sqrt => Pvc%vc_sqrt(:,iq)  ; qeq0=(normv(Qmesh%ibz(:,iq),gmet,'G')<GW_TOLQ0)

 ! Calculate Gauss-Legendre quadrature knots and weights for the omega integration
 allocate(zw(Ep%nomegaei),z(Ep%nomegaei))
 call coeffs_gausslegint(zero,one,z,zw,Ep%nomegaei)

 ! Calculate Gauss-Legendre quadrature knots and weights for the lambda integration
 allocate(zlw(nlambda),zl(nlambda))
 call coeffs_gausslegint(zero,one,zl,zlw,nlambda)


 allocate(chi0_diag(Ep%npwe))
 allocate(chitmp(Ep%npwe,Ep%npwe),STAT=istat) ; if (istat/=0) call memerr(FILE__,'chitmp',Ep%npwe**2,'gwpc')


!!static
! do ig1=1,Ep%npwe
!  write(10,'(i4,2x,f12.6,2x,f12.6)') ig1, real(sqrt(4.*pi)/vc_sqrt(ig1)), real(1.0 - vc_sqrt(ig1)**2 * chi0(ig1,ig1,1))
! enddo !ig1
! do ig2=1,Ep%npwe
!  do ig1=1,Ep%npwe
!   chitmp(ig1,ig2) = - vc_sqrt(ig1) * vc_sqrt(ig1) * chi0(ig1,ig2,1)
!  enddo !ig1
!  chitmp(ig2,ig2) = chitmp(ig2,ig2) + 1.0_dp
! enddo !ig2
! call matcginv(chitmp(:,:),Ep%npwe,Ep%npwe)
! chitmp(:,:) = matmul( chi0(:,:,1) , chitmp(:,:) )
! do ig1=1,Ep%npwe
!  chitmp(ig1,ig1) = 1.0_dp + vc_sqrt(ig1) * vc_sqrt(ig1) * chitmp(ig1,ig1)
! enddo
! do ig1=1,Ep%npwe
!  write(20,'(i4,2x,f12.6,2x,f12.6)') ig1, real(sqrt(4.*pi)/vc_sqrt(ig1)), real(chitmp(ig1,ig1))
! enddo !ig1


 do io=2,Ep%nomega ! 1,Ep%nomega
  do ig1=1,Ep%npwe
   chi0_diag(ig1) = vc_sqrt(ig1)**2 * chi0(ig1,ig1,io)
  end do

  do ilambda=1,nlambda
   lambda=zl(ilambda)

   do ig2=1,Ep%npwe
    do ig1=1,Ep%npwe
     chitmp(ig1,ig2) = - lambda * vc_sqrt(ig1) * vc_sqrt(ig1) * chi0(ig1,ig2,io)
    end do !ig1
    chitmp(ig2,ig2) = chitmp(ig2,ig2) + 1.0_dp
   end do !ig2
   call matcginv(chitmp(:,:),Ep%npwe,Ep%npwe)
   chitmp(:,:) = matmul( chi0(:,:,io) , chitmp(:,:) )
 
   do ig1=1,Ep%npwe
    chitmp(ig1,ig1) = vc_sqrt(ig1) * vc_sqrt(ig1) * chitmp(ig1,ig1)
   end do
 
   do ig1=1,Ep%npwe
    ecorr_l(ilambda) = ecorr_l(ilambda) &
&      - zw(io-1) / ( z(io-1) * z(io-1) ) * Qmesh%wt(iq) * real( chitmp(ig1,ig1) - chi0_diag(ig1) ) / (2.0_dp * pi )
   end do

  end do ! ilambda

 end do ! io
 ecorr = sum( zlw(:)*ecorr_l(:) ) 


 
 if(iq==Qmesh%nibz) then 
  unt=get_unit() ; fnam=TRIM(Dtfil%filnam_ds(4))//'_RPA' ; open(unit=unt,file=fnam)
  write(unt,'(a,(2x,f14.8))') '#RPA',ecorr
  write(msg,'(2a,(2x,f14.8))') ch10,' RPA energy [Ha] :',ecorr
  call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 end if
 do ilambda=1,nlambda
  write(msg,'(i6,2x,f10.6,2x,e12.6)') ilambda,zl(ilambda),ecorr_l(ilambda)
  call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
  if(iq==Qmesh%nibz) write(unt,'(i6,2x,f10.6,2x,e12.6)') ilambda,zl(ilambda),ecorr_l(ilambda)
 end do




 deallocate(chi0_diag,chitmp)
 deallocate(zl,zlw,z,zw)
end if !master


#if defined DEBUG_MODE
 write(msg,'(a)')' calc_rpa_functional : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine calc_rpa_functional
!!***
