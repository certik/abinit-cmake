!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_Qmesh
!! NAME
!! setup_Qmesh
!!
!! FUNCTION
!! Initialize and construct a bz_mesh_type datatype 
!! gathering information on the q-mesh used for GW calculations
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! nqibz=number of irreducible q-points
!! nsym=number of symmetry operations
!! prtvol=verbosity level
!! timrev=1 if time-reversal cannot be used, 2 otherwise
!! qibz(3,nqibz)=irreducible q-points
!! symrec(3,3,nsym)=symmetry operations in reciprocal space
!!
!! OUTPUT
!! Qmesh<bz_mesh_type>=datatype gathering information on the q point sampling. see defs_datatypes.F90
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

subroutine setup_Qmesh(nqibz,Cryst,prtvol,qibz,Qmesh)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_15gw, except_this_one => setup_Qmesh
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqibz,prtvol
 type(BZ_mesh_type),intent(inout) :: Qmesh
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 real(dp),intent(in) :: qibz(3,nqibz)

!Local variables-------------------------------
!scalars
 integer :: iq_bz,iq_ibz,isym,itim,nqbz,nqbzX,nsym,timrev
 real(dp) :: ucvol
 logical :: ltest
 character(len=500) :: msg
!arrays
 integer,allocatable :: qtab(:),qtabi(:),qtabo(:)
 integer,pointer :: symrec(:,:,:)
 real(dp) :: Sq(3),gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: qbz(:,:),wtq(:)

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' setup_qmesh : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 ltest=(Cryst%timrev==1.or.Cryst%timrev==2)
 call assert(ltest,"timrev should be 1 or 2",__FILE__,__LINE__)

 nsym   =  Cryst%nsym
 timrev =  Cryst%timrev
 symrec => Cryst%symrec

 call metric(gmet,gprimd,-1,rmet,Cryst%rprimd,ucvol)
 Qmesh%gmet   = gmet
 Qmesh%gprimd = gprimd

 nqbzX=nqibz*nsym*timrev
 allocate(qbz(3,nqbzX),qtab(nqbzX),qtabo(nqbzX),qtabi(nqbzX),wtq(nqibz))
 qbz(:,:)=0 ; qtab(:)=0 ; qtabo(:)=0 ; qtabi(:)=0

 call identq(qibz,nqibz,nqbzX,REAL(symrec,dp),nsym,timrev,wtq,qbz,qtab,qtabi,qtabo,nqbz)

 do iq_bz=1,nqbz
  isym=qtabo(iq_bz) ; iq_ibz=qtab(iq_bz) ; itim=(3-qtabi(iq_bz))/2
  call dosym(REAL(symrec(:,:,isym),dp),itim,qibz(:,iq_ibz),Sq(:))
  if (ANY(ABS(qbz(:,iq_bz)-Sq(:) )>1.0d-4)) then
   write(msg,'(4a,3f6.3,a,3f6.3,2a,9i3,2a)')ch10,&
&   ' setup_qmesh : ERROR - ',ch10,&
&   ' qpoint ',qbz(:,iq_bz),' is the symmetric of ',qibz(:,iq_ibz),ch10,&
&   ' through operation ',symrec(:,:,isym),ch10,&
&   ' however a non zero umklapp G_o vector is required and this is not yet allowed'
   call wrtout(std_out,msg,'COLL') ; write(*,*) Sq,qbz(:,iq_bz) 
   call leave_new('COLL')
  end if
 end do 
 !
 ! ==== Create data structure to store information on q-points ====
 ! * Dimensions
 Qmesh%nbz    = nqbz
 Qmesh%nibz   = nqibz      
 Qmesh%timrev = timrev
 Qmesh%nsym   = nsym
 !
 ! * Arrays
 allocate(Qmesh%ibz(3,nqibz)) ; Qmesh%ibz(:,:) = qibz(:,1:nqibz) 
 allocate(Qmesh%wt(nqibz))    ; Qmesh%wt(:)    = wtq(1:nqibz)

 allocate(Qmesh%bz(3,nqbz))   ; Qmesh%bz(:,:)  = qbz(:,1:nqbz)
 allocate(Qmesh%tab(nqbz))    ; Qmesh%tab(:)   = qtab (1:nqbz)
 allocate(Qmesh%tabi(nqbz))   ; Qmesh%tabi(:)  = qtabi(1:nqbz)
 allocate(Qmesh%tabo(nqbz))   ; Qmesh%tabo(:)  = qtabo(1:nqbz)

 !
 ! TODO For the time being these arrays are not used however they should be defined just
 ! to be consistent
 nullify(Qmesh%small)
 nullify(Qmesh%tabp)       
 nullify(Qmesh%rottb)
 nullify(Qmesh%rottbm1)

 !call print_bz_mesh(Qmesh,prtvol=prtvol)
 deallocate(qbz,qtab,qtabi,qtabo,wtq)

#if defined DEBUG_MODE
 write(msg,'(a)')' setup_qmesh : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine setup_Qmesh
!!***


!!****f* ABINIT/find_Qmesh
!! NAME
!! find_Qmesh
!!
!! FUNCTION
!!  Find the q-mesh defined as all the possible differences between k-points
!!  Find the irreducible q-points using a special treatment for the Gamma point.
!!  Then call setup_kmesh to initialize the Qmesh datatype
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Cryst<crystal_structure>=datatype gathering info on the unit cell and symmetries
!!    %nsym=number of symmetry operations
!!    %symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  prtvol=verbosity level
!!  timrev=1 if time-reversal cannot be used, 2 otherwise
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  Kmesh<bz_mesh_type>=datatype gathering information on the k-mesh
!!  nqsm=Number of small qs to be used for the long wavelength limit
!!  qsmall(3,nqsm)= the small q-points (reduced coordinates
!!
!! OUTPUT
!!  Qmesh<bz_mesh_type>=datatype gathering information on the q point sampling. see defs_datatypes.F90
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

subroutine find_Qmesh(Cryst,gprimd,Kmesh,nqsm,qsmall,Qmesh,prtvol)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15gw, except_this_one => find_Qmesh
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqsm,prtvol
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(BZ_mesh_type),intent(out) :: Qmesh
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 real(dp),intent(in) :: gprimd(3,3),qsmall(3,nqsm)

!Local variables-------------------------------
!scalars
 integer :: nqibz,nsym,timrev
 logical :: avoid_zero
 character(len=500) :: msg
!arrays
 integer,pointer :: symrec(:,:,:)
 real(dp),allocatable :: qibz(:,:)

! *************************************************************************

 timrev = Cryst%timrev
 nsym = Cryst%nsym
 symrec => Cryst%symrec
 !
 ! === Find number of q-points q=k-kp ===  
 call findnq(Kmesh%nbz,Kmesh%bz,nsym,symrec,nqibz,timrev) 
 !
 ! === Find q-points === 
 allocate(qibz(3,nqibz)) 
 avoid_zero=.TRUE.

 call findq(Kmesh%nbz,Kmesh%bz,nsym,symrec,gprimd,nqibz,qibz,timrev,avoid_zero)
 !
 ! === Create Qmesh datatype starting from the IBZ ===
 ! Here I should call setup_Kmesh, just to keep it simple but I have to 
 ! solve some problems with the small q
 call setup_Qmesh(nqibz,Cryst,prtvol,qibz,Qmesh)
 deallocate(qibz)

 ! === Define the set of small non-zero qs for long wavelength limit ===
 ! TODO check if small q-points are along high symmetry directions
 ! this part might be encapsulated in setup_Qmesh
 Qmesh%nsmall=nqsm
 allocate(Qmesh%small(3,nqsm)) ; Qmesh%small(:,:)=qsmall(:,:)

end subroutine find_Qmesh
!!***
