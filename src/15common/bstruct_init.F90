!{\src2tex{textfont=tt}}
!!****f* ABINIT/bstruct_init
!! NAME
!! bstruct_init
!!
!! FUNCTION
!! This subroutine initializes the bandstructure structured datatype
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! bantot=total number of bands (=sum(nband(:))
!! doccde(bantot)=derivative of the occupation numbers with respect to the energy (Ha)
!! eig(bantot)=eigenvalues (hartree)
!! istwfk(nkpt)=parameter that describes the storage of wfs.
!! kptns(3,nkpt)=k points in terms of recip primitive translations
!! nband(nkpt*nsppol)=number of bands
!! nkpt=number of k points
!! npwarr(nkpt)=number of planewaves at each k point
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! occ(bantot)=occupation numbers
!! wtk(nkpt)=weight assigned to each k point
!!
!! OUTPUT
!! bstruct <type(bstruct_type)>=the bandstructure datatype
!!
!! PARENTS
!!      gstate,loper3,newsp,nonlinear,respfn
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine bstruct_init(bantot,bstruct,doccde,eig,istwfk,kptns,&
& nband,nkpt,npwarr,nsppol,occ,wtk)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bantot,nkpt,nsppol
 type(bandstructure_type),intent(out) :: bstruct
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt)
 real(dp),intent(in) :: doccde(bantot),eig(bantot),kptns(3,nkpt),occ(bantot)
 real(dp),intent(in) :: wtk(nkpt)

!Local variables-------------------------------

! *************************************************************************

!Copy the scalars
 bstruct%nkpt=nkpt
 bstruct%nsppol=nsppol
 bstruct%bantot=bantot

!Allocate the components
 allocate(bstruct%nband(nkpt*nsppol), bstruct%istwfk(nkpt))
 allocate(bstruct%npwarr(nkpt),  bstruct%kptns(3,nkpt))
 allocate(bstruct%eig(bantot),   bstruct%occ(bantot))
 allocate(bstruct%doccde(bantot),bstruct%wtk(nkpt))

!Copy the arrays
 bstruct%nband(1:nkpt*nsppol)=nband(1:nkpt*nsppol)
 bstruct%istwfk(1:nkpt)      =istwfk(1:nkpt)
 bstruct%npwarr(1:nkpt)      =npwarr(1:nkpt)
 bstruct%kptns(1:3,1:nkpt)   =kptns(1:3,1:nkpt)
 bstruct%eig(1:bantot)       =eig(1:bantot)
 bstruct%occ(1:bantot)       =occ(1:bantot)
 bstruct%doccde(1:bantot)    =doccde(1:bantot)
 bstruct%wtk(1:nkpt)         =wtk(1:nkpt)

end subroutine bstruct_init
!!***
