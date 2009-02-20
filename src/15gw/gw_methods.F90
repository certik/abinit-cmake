!{\src2tex{textfont=tt}}
!!****f* ABINIT/gw_methods
!! NAME
!! gw_methods
!!
!! FUNCTION
!!  This file contains methods acting on the data types used in the 
!!  GW part of abinit. This methods can be used to nullify the pointers 
!!  defined in the data types before starting the calculation and to 
!!  deallocate the memory occupied before exiting.
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  
!!
!! OUTPUT
!!  
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  All the pointer defined in the data types could be nullified 
!!  using the null() functions. Unfortunately null() has been 
!!  introduced in the F95 specifications and this might lead 
!!  to portability problems.
!!
!! TODO 
!!  write other methods to write the content of the data type.
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

subroutine nullify_epsilonm1_parameters(Ep)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(epsilonm1_parameters),intent(inout) :: Ep

!Local variables-------------------------------
 !character(len=500) :: msg                   

! *************************************************************************
 
 nullify(Ep%qlwl   )
 nullify(Ep%omegasf)
 nullify(Ep%omega  ) 

end subroutine nullify_epsilonm1_parameters
!!***


!!****f* ABINIT/destroy_epsilonm1_parameters
!! NAME
!! destroy_epsilonm1_parameters
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine destroy_epsilonm1_parameters(Ep)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(Epsilonm1_parameters),intent(inout) :: Ep

!Local variables-------------------------------
 !character(len=500) :: msg                   

! *************************************************************************
 
 if (associated(Ep%qlwl   ))  deallocate(Ep%qlwl)
 if (associated(Ep%omegasf))  deallocate(Ep%omegasf)
 if (associated(Ep%omega  ))  deallocate(Ep%omega) 

end subroutine destroy_epsilonm1_parameters
!!***

!!****f* ABINIT/nullify_sigma_parameters
!! NAME
!! nullify_sigma_parameters
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine nullify_sigma_parameters(Sp)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(Sigma_parameters),intent(inout) :: Sp

! *************************************************************************

 nullify(Sp%kcalc )
 nullify(Sp%minbnd)
 nullify(Sp%maxbnd)
 nullify(Sp%xkcalc) 
 nullify(Sp%omegasi)
 nullify(Sp%omegasf)

end subroutine nullify_sigma_parameters
!!***


!!****f* ABINIT/destroy_sigma_parameters
!! NAME
!! destroy_sigma_parameters
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine destroy_sigma_parameters(Sp)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(Sigma_parameters),intent(inout) :: Sp

!Local variables-------------------------------
 !character(len=500) :: msg                   

! *************************************************************************

 if (associated(Sp%kcalc ))  deallocate(Sp%kcalc )
 if (associated(Sp%minbnd))  deallocate(Sp%minbnd)
 if (associated(Sp%maxbnd))  deallocate(Sp%maxbnd)

 if (associated(Sp%xkcalc))  deallocate(Sp%xkcalc) 

 if (associated(Sp%omegasi)) deallocate(Sp%omegasi) 
 if (associated(Sp%omegasf)) deallocate(Sp%omegasf) 

end subroutine destroy_sigma_parameters
!!***


!!****f* ABINIT/nullify_sigma_results
!! NAME
!! nullify_sigma_results
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine nullify_sigma_results(Sr)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(Sigma_results),intent(inout) :: Sr

!Local variables-------------------------------
 !character(len=500) :: msg                   

! *************************************************************************

 nullify(Sr%ak         )                 
 nullify(Sr%ame        )
 nullify(Sr%degwgap    )
 nullify(Sr%egwgap     )
 nullify(Sr%en_qp_diago)
 nullify(Sr%e0         )
 nullify(Sr%e0gap      )
 nullify(Sr%sigxme     )
 nullify(Sr%vxcme      )
 nullify(Sr%vUme       )

 nullify(Sr%degw       )
 nullify(Sr%dsigmee0   )
 nullify(Sr%egw        )
 nullify(Sr%eigvec_qp  )
 nullify(Sr%hhartree   )
 nullify(Sr%sigcme     )
 nullify(Sr%sigmee     )
 nullify(Sr%sigcmee0   )
 nullify(Sr%sigcmesi   )
 nullify(Sr%sigcmesrd  )
 nullify(Sr%sigxcme    )
 nullify(Sr%sigxcmesi  )
 nullify(Sr%sigxcmesrd )
 nullify(Sr%ze0        )

 nullify(Sr%omegasrd   )

end subroutine nullify_sigma_results
!!***


!!****f* ABINIT/init_sigma_results
!! NAME
!! init_sigma_results
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! TODO
!!  Write documentation.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine init_sigma_results(Sp,Kmesh,Sr)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(Sigma_parameters),intent(in) :: Sp
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Sigma_results),intent(inout) :: Sr

!Local variables-------------------------------
 !character(len=500) :: msg                   
!scalars
 integer :: b1gw,b2gw,nkibz,mod10

! *************************************************************************

 mod10=MOD(Sp%gwcalctyp,10)
 nkibz=Kmesh%nibz

 ! === min and Max GW band index over k-pts and spin ===
 ! * Used to dimension arrays 
 b1gw=Sp%minbdgw   
 b2gw=Sp%maxbdgw

 Sr%nomega=Sp%nomegasr  !FIXME remove this

 !======================================================
 ! === Allocate arrays in the sigma_results datatype ===
 !======================================================

 ! hhartree(b1,b2,k,s)= <b1,k,s|T+v_{loc}+v_{nl}+v_{H}|b2,k,s>
 allocate(Sr%hhartree(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab)) 
 Sr%hhartree=czero

 ! === QP amplitudes and energies ===
 allocate(Sr%en_qp_diago(Sp%nbnds,Kmesh%nibz,Sp%nsppol))        ; Sr%en_qp_diago(:,:,:)=zero
 allocate(Sr%eigvec_qp(Sp%nbnds,Sp%nbnds,Kmesh%nibz,Sp%nsppol)) ; Sr%eigvec_qp(:,:,:,:)=czero

 ! Dont know if it is better to do this here or in the sigma
 ! * Initialize with KS wavefunctions and energies
 !do ib=1,Sp%nbnds
 ! Sr%en_qp_diago(ib,:,:)=en(:,ib,:)
 ! Sr%eigvec_qp(ib,ib,:,:)=cone
 !end do 

 allocate(Sr%vxcme(b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab))
 allocate(Sr%vUme (b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab))

 allocate(Sr%sigxme (b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab))
 allocate(Sr%sigcme (b1gw:b2gw,Kmesh%nibz,Sr%nomega,Sp%nsppol*Sp%nsig_ab))
 allocate(Sr%sigxcme(b1gw:b2gw,Kmesh%nibz,Sr%nomega,Sp%nsppol*Sp%nsig_ab))
 allocate(Sr%sigcmee0(b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab))
 allocate(Sr%ze0(b1gw:b2gw,Kmesh%nibz,Sp%nsppol))
 allocate(Sr%dsigmee0(b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab))
 allocate(Sr%sigmee(b1gw:b2gw,Kmesh%nibz,Sp%nsppol*Sp%nsig_ab))
 allocate(Sr%degw(b1gw:b2gw,Kmesh%nibz,Sp%nsppol))
 allocate(Sr%e0(Sp%nbnds,kmesh%nibz,Sp%nsppol),Sr%egw(Sp%nbnds,Kmesh%nibz,Sp%nsppol))
 allocate(Sr%e0gap(Kmesh%nibz,Sp%nsppol),Sr%degwgap(Kmesh%nibz,Sp%nsppol),Sr%egwgap(Kmesh%nibz,Sp%nsppol))
 !allocate(Sr%ame(Sp%nbnds,Kmesh%nibz,Sr%nomega),Sr%ak(Kmesh%nibz,Sr%nomega))
 !
 ! These quantities are used to evaluate $\Sigma(E)$ around the KS\QP eigenvalue
 Sr%nomegasrd=Sp%nomegasrd
 allocate(Sr%omegasrd  (b1gw:b2gw,Kmesh%nibz,Sp%nomegasrd,Sp%nsppol))
 allocate(Sr%sigcmesrd (b1gw:b2gw,Kmesh%nibz,Sp%nomegasrd,Sp%nsppol*Sp%nsig_ab))
 allocate(Sr%sigxcmesrd(b1gw:b2gw,Kmesh%nibz,Sp%nomegasrd,Sp%nsppol*Sp%nsig_ab))

 Sr%e0(:,:,:)=zero
 Sr%egw(:,:,:)=czero
 Sr%e0gap(:,:)=zero
 Sr%sigcme(:,:,:,:)=czero
 Sr%sigxme(:,:,:)=czero
 Sr%sigxcme(:,:,:,:)=czero
 Sr%sigcmee0(:,:,:)=czero
 Sr%ze0(:,:,:)=czero
 Sr%dsigmee0(:,:,:)=czero
 Sr%sigmee(:,:,:)=czero
 Sr%omegasrd(:,:,:,:)=czero
 Sr%sigcmesrd(:,:,:,:)=czero
 Sr%sigxcmesrd(:,:,:,:)=czero
 Sr%degw(:,:,:)=czero

 ! === Analytical Continuation ===
 if (mod10==1) then 
  allocate(Sr%sigcmesi (b1gw:b2gw,Kmesh%nibz,Sp%nomegasi,Sp%nsppol*Sp%nsig_ab))
  allocate(Sr%sigxcmesi(b1gw:b2gw,Kmesh%nibz,Sp%nomegasi,Sp%nsppol*Sp%nsig_ab))
  Sr%sigcmesi (:,:,:,:)=czero
  Sr%sigxcmesi(:,:,:,:)=czero
 end if

end subroutine init_sigma_results
!!***


!!****f* ABINIT/destroy_sigma_results
!! NAME
!! destroy_sigma_results
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine destroy_sigma_results(Sr)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(Sigma_results),intent(inout) :: Sr

!Local variables-------------------------------
 !character(len=500) :: msg                   

! *************************************************************************

 if (associated(Sr%ak         ))    deallocate(Sr%ak)                 
 if (associated(Sr%ame        ))    deallocate(Sr%ame)
 if (associated(Sr%degwgap    ))    deallocate(Sr%degwgap)
 if (associated(Sr%egwgap     ))    deallocate(Sr%egwgap)
 if (associated(Sr%en_qp_diago))    deallocate(Sr%en_qp_diago)
 if (associated(Sr%e0         ))    deallocate(Sr%e0)
 if (associated(Sr%e0gap      ))    deallocate(Sr%e0gap)
 if (associated(Sr%sigxme     ))    deallocate(Sr%sigxme)
 if (associated(Sr%vxcme      ))    deallocate(Sr%vxcme)
 if (associated(Sr%vUme       ))    deallocate(Sr%vUme)
 
 if (associated(Sr%degw       ))    deallocate(Sr%degw)
 if (associated(Sr%dsigmee0   ))    deallocate(Sr%dsigmee0)
 if (associated(Sr%egw        ))    deallocate(Sr%egw)
 if (associated(Sr%eigvec_qp  ))    deallocate(Sr%eigvec_qp)
 if (associated(Sr%hhartree   ))    deallocate(Sr%hhartree)
 if (associated(Sr%sigcme     ))    deallocate(Sr%sigcme)
 if (associated(Sr%sigmee     ))    deallocate(Sr%sigmee)
 if (associated(Sr%sigcmee0   ))    deallocate(Sr%sigcmee0)
 if (associated(Sr%sigcmesi   ))    deallocate(Sr%sigcmesi)
 if (associated(Sr%sigcmesrd  ))    deallocate(Sr%sigcmesrd)
 if (associated(Sr%sigxcme    ))    deallocate(Sr%sigxcme)
 if (associated(Sr%sigxcmesi  ))    deallocate(Sr%sigxcmesi)
 if (associated(Sr%sigxcmesrd ))    deallocate(Sr%sigxcmesrd)
 if (associated(Sr%ze0        ))    deallocate(Sr%ze0)
 
 if (associated(Sr%omegasrd   ))    deallocate(Sr%omegasrd)

end subroutine destroy_sigma_results
!!***
