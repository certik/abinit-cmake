!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_Kmesh
!! NAME
!! setup_Kmesh
!!
!! FUNCTION
!!  Initialize and construct a bz_mesh_type datatype 
!!  gathering information on the k-mesh used for GW calculations
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nkibz=number of irreducible k-points
!!  prtvol=verbosity level
!!  timrev=1 if time-reversal cannot be used, 2 otherwise
!!  kibz(3,nkibz)=irreducible k-points
!!  Cryst<Crystal_structure> = Info on unit cell and its symmetries
!!     %nsym=number of symmetry operations
!!     %symrec(3,3,nsym)=symmetry operations in reciprocal space
!!     %tnons(3,nsym)=fractional translations
!!
!! OUTPUT
!!  Kmesh<bz_mesh_type>=datatype gathering information on the k point sampling. see defs_datatypes.F90
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!  
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine setup_Kmesh(nkibz,kibz,Cryst,Kmesh,prtvol)

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_15gw, except_this_one => setup_Kmesh
!End of the abilint section

 implicit none

!Arguments ------------------------------------
! scalars
! arrays
!scalars
 integer,intent(in) :: nkibz,prtvol
 type(Bz_mesh_type),intent(inout) :: Kmesh
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 real(dp),intent(in) :: kibz(3,nkibz)

!Local variables-------------------------------
!scalars
 integer :: ierr,ik_bz,ik_ibz,isym,nkbz,nkbzX,nsym,timrev
 real(dp) :: shift1,shift2,shift3,ucvol
 logical :: ltest,use_antiferro
 character(len=500) :: msg
!arrays
 integer,allocatable :: ktab(:),ktabi(:),ktabo(:),ktabr(:,:)
 integer,pointer :: symafm(:),symrec(:,:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rm1t(3),rmet(3,3)
 real(dp),allocatable :: kbz(:,:),kbz_wrap(:,:),wtk(:)
 real(dp),pointer :: tnons(:,:)

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' setup_kmesh : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 ! === Initial tests on input ===
 ltest=(Cryst%timrev==1.or.Cryst%timrev==2)
 write(msg,'(a,i4)')'Wrong value for timrev= ',Cryst%timrev 
 call assert(ltest,msg,__FILE__,__LINE__)

 ! === Get symmetries related stuff ==
 nsym   = Cryst%nsym
 timrev = Cryst%timrev
 use_antiferro = Cryst%use_antiferro
 symrec => Cryst%symrec
 symafm => Cryst%symafm
 tnons  => Cryst%tnons

 call metric(gmet,gprimd,-1,rmet,Cryst%rprimd,ucvol)
 Kmesh%gmet   = gmet
 Kmesh%gprimd = gprimd
 !
 ! === Find BZ from IBZ and fill tables === 
 nkbzX=nkibz*nsym*timrev ! Maximum possible number 
 allocate(kbz(3,nkbzX),wtk(nkibz),ktab(nkbzX),ktabi(nkbzX),ktabo(nkbzX))

 call identk(kibz,nkibz,nkbzX,nsym,timrev,symrec,symafm,use_antiferro,kbz,ktab,ktabi,ktabo,nkbz,wtk)
 !
 ! Wrap the BZ points in the interval ]-1/2,1/2]
 !allocate(kbz_wrap(3,nkbz))
 !do ik_bz=1,nkbz
 ! call canon9(kbz(1,ik_bz),kbz_wrap(1,ik_bz),shift1)
 ! call canon9(kbz(2,ik_bz),kbz_wrap(2,ik_bz),shift2)
 ! call canon9(kbz(3,ik_bz),kbz_wrap(3,ik_bz),shift3)
 !end do
 !
 ! ==== Create data structure to store information on k-points ====
 !
 ! === Dimensions ===
 Kmesh%nbz=nkbz         ! Number of points in the full BZ
 Kmesh%nibz=nkibz       ! Number of points in the IBZ
 Kmesh%nsym=nsym        ! Number of operations
 Kmesh%nsmall=0         ! No need to deal with coulombian singularity
 Kmesh%timrev=timrev    ! 2 if time-reversal is used, 1 otherwise 
 ! 
 ! === Arrays ===
 allocate(Kmesh%bz(3,nkbz))   ;  kmesh%bz(:,:)=  kbz(:,1:nkbz)  ! Red. coordinates of points in full BZ.
 allocate(Kmesh%ibz(3,nkibz)) ; kmesh%ibz(:,:)=kibz(:,1:nkibz)  ! Red. coordinates of points in IBZ.
 !allocate(Kmesh%bz(3,nkbz))   ; kmesh%bz(:,:)=kbz_wrap(:,1:nkbz)

 nullify (Kmesh%small)
 !allocate(Kmesh%small(3,Kmesh%nsmall))                         ! Small k-points (not used)

 allocate(Kmesh%tab(nkbz))    ; Kmesh%tab(:)=ktab (1:nkbz)      ! Index of the irred. point in the array IBZ.
 allocate(Kmesh%tabi(nkbz))   ; Kmesh%tabi(:)=ktabi(1:nkbz)     !-1 if time reversal must be used to obtain this point, 
                                                                ! 1 otherwise
 allocate(Kmesh%tabo(nkbz))   ; Kmesh%tabo(:)=ktabo(1:nkbz)     ! Symm. operation that rotates k_IBZ onto \pm k_BZ 
                                                                ! (depending on tabi)
 allocate(Kmesh%wt(nkibz))    ; Kmesh%wt(:)=wtk(1:nkibz)        ! Weight for each k_IBZ

 allocate(Kmesh%rottbm1(nkbz,timrev,nsym))                      ! Index of rotated point (IS)^{-1} kbz 
 allocate(Kmesh%rottb  (nkbz,timrev,nsym))                      ! Index of rotated point IS kbz where I is either the
                                                                ! identity or time-reversal

 call setup_k_rotation(nsym,symrec,timrev,nkbz,Kmesh%bz,Kmesh%rottb,Kmesh%rottbm1)

 allocate(Kmesh%tabp(nkbz)) ! Phase factors for non-symmorphic operations $e{-i2\pik_BZ.\tau}$
 do ik_bz=1,nkbz
  isym=Kmesh%tabo(ik_bz) ; ik_ibz=Kmesh%tab(ik_bz)
  rm1t=MATMUL(TRANSPOSE(symrec(:,:,isym)),tnons(:,isym))
  Kmesh%tabp(ik_bz)=EXP(-(0.,1.)*two_pi*DOT_PRODUCT(kibz(:,ik_ibz),rm1t))
 end do  

 deallocate(kbz,wtk,ktab,ktabi,ktabo)
 !deallocate(kbz_wrap)

 call print_bz_mesh(Kmesh,prtvol=prtvol)
 !call print_bz_mesh(Kmesh,prtvol=10)

#if defined DEBUG_MODE
 write(msg,'(a)')' setup_Kmesh : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine setup_Kmesh
!!***


!!****f* ABINIT/destroy_bz_mesh_type
!! NAME
!! destroy_bz_mesh_type
!!
!! FUNCTION
!! Deallocate a previously allocated data bz_mesh_type structure.
!!
!! INPUTS
!! Kmesh<bz_mesh_type>: the datatype to be destroyed
!!
!! PARENTS
!!
!! CHILDREN
!!
!! OUTPUT
!! See SIDE EFFECTS
!!
!! SIDE EFFECTS
!! Deallocate all the (associated) pointers defined in the bz_mesh_type datatype
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine destroy_bz_mesh_type(Kmesh)

 use defs_basis
 use defs_datatypes 

 implicit none

!Arguments ------------------------------------
!scalars
 type(BZ_mesh_type),intent(inout) :: Kmesh

!Local variables-------------------------------
!scalars
 integer :: istat

! *********************************************************************

 if (associated(Kmesh%ibz    )) deallocate(Kmesh%ibz    )
 if (associated(Kmesh%bz     )) deallocate(Kmesh%bz     )
 if (associated(Kmesh%small  )) deallocate(Kmesh%small,stat=istat)
 if (associated(Kmesh%rottb  )) deallocate(Kmesh%rottb  )
 if (associated(Kmesh%rottbm1)) deallocate(Kmesh%rottbm1)
 if (associated(Kmesh%tab    )) deallocate(Kmesh%tab    )
 if (associated(Kmesh%tabi   )) deallocate(Kmesh%tabi   )
 if (associated(Kmesh%tabo   )) deallocate(Kmesh%tabo   )
 if (associated(Kmesh%tabp   )) deallocate(Kmesh%tabp   )
 if (associated(Kmesh%wt     )) deallocate(Kmesh%wt     )

end subroutine destroy_bz_mesh_type
!!***


!!****f* ABINIT/print_bz_mesh_type
!! NAME
!! print_bz_mesh_type
!!
!! FUNCTION
!! Print the content of a bz_mesh_type datatype
!!
!! INPUTS
!! Kmesh<bz_mesh_type>=the datatype to be printed
!! unit(optional)=the unit number for output 
!! prtvol(optional)=verbosity level
!! mode_paral=either "COLL" or "PERS"
!!
!! PARENTS
!!
!! CHILDREN
!!
!! OUTPUT
!! Only printing 
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine print_BZ_mesh(Kmesh,unit,prtvol,mode_paral)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 type(BZ_mesh_type),intent(in) :: Kmesh

!Local variables-------------------------------
!scalars
 integer,parameter :: nmaxk=50
 integer :: ii,ik,unt,verbose
 character(len=100) :: fmt
 character(len=4) :: mode
 character(len=500) :: msg

! *************************************************************************

 unt=std_out ; if (PRESENT(unit)) unt=unit
 verbose=0   ; if (PRESENT(prtvol)) verbose=prtvol
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 write(msg,'(2a,i5,3a)')ch10,&
& ' Number of points in the IBZ : ',Kmesh%nibz,ch10,&
& ' Reduced Coordinates and Weights : ',ch10
 call wrtout(unt,msg,mode)

 write(fmt,*)'(1x,i5,a,2x,3es16.8,3x,f11.5)'
 do ik=1,Kmesh%nibz
  write(msg,fmt) ik,') ',(Kmesh%ibz(ii,ik),ii=1,3),Kmesh%wt(ik)
  call wrtout(unt,msg,mode)
 end do

 if (Kmesh%timrev==2) then 
  write(msg,'(2a,i2,3a,i5,a)')ch10,&
&  ' together with ',Kmesh%nsym,' symmetry operations and time-reversal ',ch10,&
&  ' have yielded ',Kmesh%nbz,' k-points in the Brillouin Zone (BZ) :'
 else if (Kmesh%timrev==1) then
  write(msg,'(2a,i2,3a,i5,a)')ch10,&
&  ' together with ',Kmesh%nsym,' symmetry operations (time-reversal not used) ',ch10,&
&  ' have yielded ',Kmesh%nbz,' k-points in the Brillouin Zone (BZ) :'
 else 
  write(msg,'(4a,i3)')ch10,&
&  ' print_bz_mesh : BUG -',ch10,&
&  ' wrong value for timrev = ',Kmesh%timrev
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if
 call wrtout(unt,msg,mode)

 write(fmt,*)'(1x,i5,a,2x,3es16.8)'
 do ik=1,Kmesh%nbz
  if (verbose==0 .and. ik>nmaxk) then
   write(msg,'(a)')' prtvol=0, do not print more k-points.'
   call wrtout(unt,msg,mode) ; EXIT
  end if
  write(msg,fmt)ik,') ',(Kmesh%bz(ii,ik),ii=1,3)
  call wrtout(unt,msg,mode)
 end do
 !
 ! === If used, write small-points around Gamma point ===
 if (Kmesh%nsmall/=0) then 
  write(msg,'(a)')' Small points around Gamma '
  call wrtout(unt,msg,mode)
  write(fmt,*)'(1x,i5,a,2x,3es16.8)'
  do ik=1,Kmesh%nsmall
   write(msg,fmt)ik,') ',(Kmesh%small(ii,ik),ii=1,3)
   call wrtout(unt,msg,mode)
  end do
 end if
 !
 ! === Additional printing ===
 if (verbose>=10) then 
  write(msg,'(2a)')ch10,&
&  '                  Irred point -->              Full-point           through:  Symrec  Time-Rev (1=No,-1=Yes) '
  call wrtout(unt,msg,mode)
  write(fmt,*)'(2x,i5,2x,2(3es16.8,2x),i3,2x,i2)'
  do ik=1,Kmesh%nbz
   write(msg,fmt)ik,Kmesh%bz(:,ik),Kmesh%ibz(:,Kmesh%tab(ik)),Kmesh%tabo(ik),Kmesh%tabi(ik)
   call wrtout(unt,msg,mode)
  end do
 end if

end subroutine print_BZ_mesh
!!***


!!****f* ABINIT/setup_k_rotation
!! NAME
!! setup_k_rotation
!!
!! FUNCTION
!! Set up tables indicating rotation of k-points
!!
!! INPUTS
!! kbz(3,nbz)=reduced coordinates of k-points
!! timrev=if 2, take into account time-reversal, 1 otherwise
!! nsym=number of symmetry operations
!! nbz=number of k-points
!! symrec(3,3,nsym)=symmetry operations in reciprocal space (reduced coordinates)
!!
!! OUTPUT
!! krottb(k,I,S)=Index of (IS) k in the array bz
!! krottbm1(k,I,S)=Index of IS^{-1} k
!!
!! NOTES: 
!! Always use intent(inout) because datatype elements are passed to this routine
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

subroutine setup_k_rotation(nsym,symrec,timrev,nbz,kbz,krottb,krottbm1)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => setup_k_rotation
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbz,nsym,timrev
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 integer,intent(inout) :: krottb(nbz,timrev,nsym),krottbm1(nbz,timrev,nsym)
 real(dp),intent(in) :: kbz(3,nbz)

!Local variables ------------------------------
!scalars
 integer :: ik,ikp,isym,itim
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: G0(3)
 real(dp) :: kbase(3),krot(3)

!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' setup_k_rotation : setting up k-rotation tables'
 call wrtout(std_out,msg,'PERS')
#endif
 !
 ! === Set up k-rotation tables ===
 do ik=1,nbz
  kbase(:)=kbz(:,ik)
  do itim=1,timrev
   do isym=1,nsym
    krot(:)=(3-2*itim)*MATMUL(symrec(:,:,isym),kbase)
    found=.FALSE.
    do ikp=1,nbz
     if (is_samek(krot,kbz(:,ikp),G0)) then
      found=.TRUE.
      krottb  (ik ,itim,isym)=ikp
      krottbm1(ikp,itim,isym)=ik
     end if
    end do 
    if (.not.found) then
     write(msg,'(6a,i4,a,i4,2x,2(3f12.6,2a),i3,a,i2)')ch10,&
&     ' setup_k_rotation : ERROR-',ch10,&
&     ' k-mesh not closed',ch10,&
&     ' Initial k-point ',ik,'/',nbz,kbase(:),ch10,&
&     ' Rotated k-point ',krot(:),ch10,&
      ' Through sym ',isym,' and itim ',itim
     call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
    end if
   end do 
  end do 
 end do 

#if defined DEBUG_MODE
 write(msg,'(2a)')' k-rotation tables set up',ch10
 call wrtout(std_out,msg,'PERS')
#endif

end subroutine setup_k_rotation
!!***


!!****f* ABINIT/get_BZ_item
!! NAME
!! get_BZ_item
!!
!! FUNCTION
!! Given the index of a point in the full BZ, report all related information. 
!!
!! INPUTS
!! ikbz=The index of the required point in the BZ
!! Kmesh<bz_mesh_type>=datatype gathering information on the k point sampling. see defs_datatypes.F90
!!
!! OUTPUT
!! kbz(3)=the k-point in reduced coordinated 
!! isym=index of the symmetry required to rotate ik_ibz onto ik_bz
!! itim=2 is time-reversal has to be used, 1 otherwise 
!! ik_ibz=the index of corresponding point in the IBZ
!! phnons=the phase factor for non-symmorphic operations 
!!  i.e e^{-i 2 \pi k_IBZ \cdot R{^-1}t}=e{-i 2\pi k_BZ cdot t}
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

subroutine get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym,itim,ph_mkbzt)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz
 integer,intent(out) :: ik_ibz,isym,itim
 complex(gwpc),intent(out),optional :: ph_mkbzt
 type(BZ_mesh_type),intent(in) :: Kmesh
!arrays
 real(dp),intent(out) :: kbz(3)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *********************************************************************

 if (ik_bz>Kmesh%nbz.or.ik_bz<=0) then 
  write(msg,'(4a,i3)')ch10,&
&  ' get_BZ_item : BUG ',ch10,&
&  ' wrong value for ik_bz: ',ik_bz
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 kbz(:)=Kmesh%bz(:,ik_bz)

 ik_ibz=Kmesh%tab(ik_bz)
 isym=Kmesh%tabo(ik_bz)
 itim=(3-Kmesh%tabi(ik_bz))/2
 if (PRESENT(ph_mkbzt)) ph_mkbzt=Kmesh%tabp(ik_bz)

end subroutine get_BZ_item
!!***


!!****f* ABINIT/get_IBZ_item
!! NAME
!! get_IBZ_item
!!
!! FUNCTION
!! Given the idex of a point in the IBZ, report information on this point
!!
!! INPUTS
!! ik_bz=The index of the required point in the IBZ
!! Kmesh<bz_mesh_type>=datatype gathering information on the k point sampling. see defs_datatypes.F90
!!
!! OUTPUT
!! kibz(3)=the k-point in reduced coordinated 
!! wtk=the weight 
!!
!! TODO 
!!  Add mapping ibz2bz, ibz2star
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

subroutine get_IBZ_item(Kmesh,ik_ibz,kibz,wtk)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz
 real(dp),intent(out) :: wtk
 type(bz_mesh_type),intent(in) :: Kmesh
!arrays
 real(dp),intent(out) :: kibz(3)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *********************************************************************

 if (ik_ibz>Kmesh%nibz.or.ik_ibz<=0) then 
  write(msg,'(4a,i3)')ch10,&
&  ' get_BZ_item : BUG ',ch10,&
&  ' wrong value for ik_ibz: ',ik_ibz
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 kibz(:)=Kmesh%ibz(:,ik_ibz)
 wtk=Kmesh%wt(ik_ibz)

end subroutine get_IBZ_item
!!***


!!****f* ABINIT/get_BZ_diff
!! NAME
!! get_BZ_diff
!!
!! FUNCTION
!! Given two points k1 and k2 where k1 belongs to the BZ, check if the difference
!! k1-k2 still belongs to the BZ reporting useful quantities
!!
!! INPUTS
!!  Kmesh<bz_mesh_type>=datatype gathering information on the k-mesh
!!  k1(3)=the first k-points (supposed to be in the BZ)
!!  k2(3)=the second point
!!
!! OUTPUT
!!  idiff_bz=the idex of k1-k2 in the BZ 
!!  G0(3)=the umklapp G0 vector required to bring k1-k2 back to the BZ 
!!  nfound= the number of points in the BZ that are equal to k1-k2 (should be 1 if everything is OK)
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

subroutine get_BZ_diff(Kmesh,k1,k2,idiff_bz,G0,nfound)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => get_BZ_diff
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: idiff_bz,nfound
 type(bz_mesh_type),intent(in) :: Kmesh
!arrays
 integer,intent(out) :: G0(3)
 real(dp),intent(in) :: k1(3),k2(3)

!Local variables-------------------------------
!scalars
 integer :: ikp
 character(len=500) :: msg
!arrays
 integer :: umklp(3)
 real(dp) :: kdiff(3),ktrial(3)

! *********************************************************************

 if (.not.has_BZ_item(Kmesh,k1,ikp,umklp)) then 
  write(msg,'(4a,3f12.6)')ch10,&
&  ' get_BZ_diff : ERROR ',ch10,&
&  ' first point must be in BZ ',k1(:)
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 kdiff(:)=k1(:)-k2(:) 
 nfound=0 ; idiff_bz=0
 !
 ! === Find p such k1-k2=p+G0 where p in the BZ === 
 do ikp=1,Kmesh%nbz
  ktrial=Kmesh%bz(:,ikp)
  if (is_samek(kdiff,ktrial,umklp)) then 
   idiff_bz=ikp
   G0=umklp
   nfound=nfound+1
  end if
 end do
 !
 ! === Check if p has not found of found more than once ===
 ! For extremely dense meshes, tol1q in defs_basis might be too large!
 if (nfound/=1) then 
  if (nfound==0) then
   write(msg,'(4a)')ch10,&
&   ' get_BZ_diff : WARNING - ',ch10,&
&   ' k1-k2-G0 not found in BZ '
  else 
   write(msg,'(4a,i3)')ch10,&
&   ' get_BZ_diff : COMMENT - ',ch10,&
&   ' Multiple k1-k2-G0 found in BZ ',nfound
  end if
  call wrtout(std_out,msg,'COLL') 
  write(msg,'(5a,3(a,3f12.6,a))')ch10,&
&  ' k1    = ',k1(:),ch10,&
&  ' k2    = ',k2(:),ch10,&
&  ' k1-k2 = ',kdiff(:),ch10
  call wrtout(std_out,msg,'COLL') !; call leave_new('COLL') 
 end if

end subroutine get_BZ_diff
!!***

!!****f* ABINIT/is_samek
!! NAME
!! is_samek
!!
!! FUNCTION
!! Compare two k-points
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!
!! INPUTS
!!  k1(3),k2(3)=two k points to be compared
!!
!! OUTPUT
!! Return .TRUE. if they are the same within a RL vector,
!!        .FALSE. if they are different.
!! G0(3)=if .TRUE. G0(3) is the reciprocal lattice vector such as k1=k2+G0
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

logical function is_samek(k1,k2,G0)

 use defs_basis
 use m_gwdefs, only : GW_TOLQ

 implicit none

!Arguments ------------------------------------
!arrays
 integer,intent(out) :: G0(3)
 real(dp),intent(in) :: k1(3),k2(3)

!Local variables-------------------------------
!scalars
 real(dp) :: f,x

! *************************************************************************

 ! === Statement function definition: f is zero only if x is integer ===
 f(x)=ABS(x-NINT(x))

 is_samek=.FALSE. ; G0(:)=0
 if (f(k1(1)-k2(1))<GW_TOLQ) then
  if (f(k1(2)-k2(2))<GW_TOLQ) then
   if (f(k1(3)-k2(3))<GW_TOLQ) then
    is_samek=.TRUE.
    G0(:)=NINT(k1(:)-k2(:))
   end if
  end if
 end if

end function is_samek
!!***


!!****f* ABINIT/has_BZ_item
!! NAME
!! has_BZ_item
!!
!! FUNCTION
!!  check if item belongs to the BZ  within a reciprocal lattice vector
!!  and return the index number and the reciprocal vector g0.
!!
!! INPUTS
!!  Kmesh<bz_mesh_type>=datatype gathering information on the k-mesh
!!  item(3)=the k-point to be checked
!!
!! OUTPUT
!!  .TRUE. if item is the BZ within a RL vector
!!  ikbz=Index of the k-point in the Kmesh%bz array
!!  G0(3)=Umklapp vector.
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

logical function has_BZ_item(Kmesh,item,ikbz,G0)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => has_BZ_item
!End of the abilint section

 implicit none

!Arguments ------------------------------------
! scalars
! arrays
!scalars
 integer,intent(out) :: ikbz
 type(BZ_mesh_type),intent(in) :: Kmesh
!arrays
 integer,intent(out) :: G0(3)
 real(dp),intent(in) :: item(3)

!Local variables-------------------------------
!scalars
 integer :: ik_bz,yetfound
 character(len=500) :: msg

! *************************************************************************

 has_BZ_item=.FALSE. ; ikbz=0 ; G0(:)=0 ; yetfound=0
 do ik_bz=1,Kmesh%nbz
  if (is_samek(item,Kmesh%bz(:,ik_bz),G0)) then 
   has_BZ_item=.TRUE. 
   ikbz=ik_bz
   yetfound=yetfound+1
   !EXIT
  end if
 end do

 if (yetfound/=0.and.yetfound/=1) then
  write(msg,'(4a)')ch10,' has_BZ_item : BUG ',ch10,' multiple k-points found'   
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

end function has_BZ_item
!!***


!!****f* ABINIT/has_IBZ_item
!! NAME
!! has_IBZ_item
!!
!! FUNCTION
!!  Check if item belongs to the IBZ within a reciprocal lattice vector
!!
!! INPUTS
!!  Kmesh<bz_mesh_type>=datatype gathering information on the k-mesh
!!  item(3)=the k-point to be checked
!!
!! OUTPUT
!!  .TRUE. if item is the IBZ within a RL vector
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

logical function has_IBZ_item(Kmesh,item,ikibz,G0)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => has_IBZ_item
!End of the abilint section

 implicit none

!Arguments ------------------------------------
! scalars
! arrays
!scalars
 integer,intent(out) :: ikibz
 type(BZ_mesh_type),intent(in) :: Kmesh
!arrays
 integer,intent(out) :: G0(3)
 real(dp),intent(in) :: item(3)

!Local variables-------------------------------
!scalars
 integer :: ik_ibz,yetfound
 character(len=500) :: msg

! *************************************************************************

 has_IBZ_item=.FALSE. ; ikibz=0 ; G0(:)=0 ; yetfound=0
 do ik_ibz=1,Kmesh%nibz
  if (is_samek(item,Kmesh%ibz(:,ik_ibz),G0)) then 
   has_IBZ_item=.TRUE. 
   ikibz=ik_ibz
   yetfound=yetfound+1
   !EXIT
  end if
 end do

 if (yetfound/=0.and.yetfound/=1) then
  write(msg,'(4a)')ch10,' has_IBZ_item : BUG ',ch10,' multiple k-points found'   
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

end function has_IBZ_item
!!***
