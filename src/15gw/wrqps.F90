!{\src2tex{textfont=tt}}
!!****f* ABINIT/wrqps
!! NAME
!! wrqps
!!
!! FUNCTION
!!  Write the _QPS file containing information on the quasi-particles energies and wavefunctions 
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2008 ABINIT group (FBruneval, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Dtfil= datafile type gathering all the variables related to files names
!!  Sp<Sigma_parameters>=Parameters characterizing the self-energy calculation.  
!!     %nsppol=1 for unpolarized, 2 for spin-polarized
!!     %nbnds=number of bands used for sigma
!!  Sr<Sigma_results>=Structure containing the results of the sigma run.
!!     %en_qp_diago(nbnds,nibz,nsppol)= NEW quasi-particle energies
!!     %eigvec_qp(nbnds,nbnds,nibz,nsppol)= NEW QP amplitudes in the KS basis set 
!!      obtained by diagonalizing H0 + Herm(Sigma).
!!  m_lda_to_qp(nbnds,nbnds,nibz,nsppol)= expansion of the OLD QP amplitudes in terms of KS wavefunctions
!!  Kmesh<Bz_mesh_type>=information on the k-point sampling.
!!     %nibz=number of irreducible k-points
!!     %ibz(3,kibz)=reduced coordinates of the irreducible k-points
!!  nscf=number of self consistent cycles performed
!!  nspden=number of spin-density components
!!
!! OUTPUT
!!  Only writing
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wrqps(Dtfil,Sp,Kmesh,nspden,nscf,nfftot,Sr,m_lda_to_qp,rho_qp)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot,nscf,nspden
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Datafiles_type),intent(in) :: Dtfil
 type(Sigma_parameters),intent(in) :: Sp
 type(Sigma_results),intent(in) :: Sr
!arrays
 real(dp),intent(in) :: rho_qp(nfftot,nspden)
 complex(dpc),intent(in) :: m_lda_to_qp(Sp%nbnds,Sp%nbnds,Kmesh%nibz,Sp%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ib,ik,is,unqps
 character(len=500) :: msg
 character(len=fnlen) :: filnam
!arrays
 complex(dpc),allocatable :: mtmp(:,:,:)

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' wrqps : enter ' 
 call wrtout(std_out,msg,'COLL')
 call flush_unit(std_out)
#endif

 !TODO add the abinit header, it may happens indeed that 
 !the calculation crashes if the density in the QPS file
 !is on a FFT mesh different from that used in the calling routine
 
 unqps=Dtfil%unqps
 filnam=TRIM(Dtfil%filnam_ds(4))//'_QPS'
 open(unit=unqps,file=TRIM(filnam),form='formatted',status='unknown')

 write(msg,'(3a)')ch10,&
& ' writing QP data on file : ',TRIM(Dtfil%filnam_ds(4))//'_QPS'
 call wrtout(std_out,msg,'COLL') 
 call wrtout(ab_out,msg,'COLL')

 write(unqps,*)nscf+1
 write(unqps,*)Kmesh%nibz
 write(unqps,*)Sp%nbnds
 write(unqps,*)Sp%nsppol
 !
 ! === Calculate the new m_lda_to_qp ===
 allocate(mtmp(Sp%nbnds,Sp%nbnds,Sp%nsppol))
 do is=1,Sp%nsppol
  do ik=1,Kmesh%nibz
   mtmp(:,:,is)=MATMUL(m_lda_to_qp(:,:,ik,is),Sr%eigvec_qp(:,:,ik,is))
   write(unqps,*) Kmesh%ibz(:,ik)
   do ib=1,Sp%nbnds
    write(unqps,*)Sr%en_qp_diago(ib,ik,is)
    write(unqps,*)mtmp(:,ib,is)
   end do 
  end do 
 end do 
 deallocate(mtmp)
 !
 ! === Write QP density ===
 write(unqps,*)rho_qp(:,:)
 close(unqps)

#if defined DEBUG_MODE
 write(msg,'(a)')' sigma : SCF data written '
 call wrtout(std_out,msg,'COLL')
 call flush_unit(std_out)
#endif

end subroutine wrqps
!!***
