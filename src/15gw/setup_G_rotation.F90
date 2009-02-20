!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_G_rotation
!! NAME
!! setup_G_rotation
!!
!! FUNCTION
!! Set up tables indicating rotation of G-vectors
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MT, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gvec(3,npw)=coordinates of plane waves, supposed to be ordered in increasing modulus
!! timrev=if 2, take into account time-reversal, 1 otherwise
!! nsym=number of symmetry operations
!! npw=number of planewaves used
!! symrec(3,3,nsym)=symmetry operations in reciprocal space
!!
!! OUTPUT
!!  grottb  (npw,2,nsym)= grottb(G,I,S) is the index of (SI) G in the array gvec; I is the identity or the inversion 
!!  grottbm1(npw,2,nsym)= index no of IS^{-1} G 
!!
!! NOTES: 
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

subroutine setup_G_rotation(only_one_kpt,nsym,symrec,timrev,npw,gvec,g2sh,nsh,shlim,grottb,grottbm1)

 use defs_basis
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nsh,nsym,timrev
 logical,intent(in) :: only_one_kpt
!arrays
 integer,intent(in) :: g2sh(npw),gvec(3,npw),shlim(nsh+1),symrec(3,3,nsym)
 integer,intent(inout) :: grottb(npw,timrev,nsym),grottbm1(npw,timrev,nsym)

!Local variables ------------------------------
!scalars
 integer :: ee,ig1,ig2,ish1,isym,itim,ss
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: gbase(3),grot(3)

!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(2a)')ch10,&
& ' setup_G_rotation : setting up G-rotation tables'
 call wrtout(std_out,msg,'COLL')
 call flush_unit(std_out)
#endif
 !
 ! === Set up G-rotation table ===
 ! * This loop might be CPU consuming in isolated systems
 !   Therefore we skip it in case of one single k-point
 if (only_one_kpt) then
  ! As only_one_kpt is true, the only symmetry needed is identity
  do ig1=1,npw
   grottb  (ig1,1,1)=ig1
   grottbm1(ig1,1,1)=ig1
   !TODO check if also inversion might enter somewhere!!!
  end do
 else 
  ! === Several k-points ===
  do ig1=1,npw
   ish1=g2sh(ig1) ; ss=shlim(ish1) ; ee=shlim(ish1+1)-1
   gbase(:)=gvec(:,ig1)
   do itim=1,timrev
    do isym=1,nsym
     grot=(3-2*itim)*MATMUL(symrec(:,:,isym),gbase)
     found=.FALSE.
     ! * Looping on the shell of ig1 to have better scalling
     do ig2=ss,ee 
      if (ALL(ABS(grot(:)-gvec(:,ig2))==0)) then
       found=.TRUE.
       grottb  (ig1,itim,isym)=ig2
       grottbm1(ig2,itim,isym)=ig1
      end if
     end do 
     if (.not.found) then
      write(msg,'(6a,i5,a,i5,1x,2(3i5,a),a,i3,a,i3)')ch10,&
&      ' setup_G_rotation : ERROR-',ch10,&
&      ' G-shell not closed',ch10,&
&      ' Initial G vector ',ig1,'/',npw,gbase(:),' Rotated G vector ',grot(:),ch10,&
       ' Through sym ',isym,' and itim ',itim
      call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
     end if
    end do 
   end do 
  end do 
 end if !only_one_kpt

#if defined DEBUG_MODE
 write(msg,'(a)')' G-rotation tables set up'
 call wrtout(std_out,msg,'COLL')
 call flush_unit(std_out)
#endif

end subroutine setup_G_rotation
!!***

!!****f* ABINIT/init_Gvectors_type 
!! NAME
!! init_Gvectors_type
!!
!! FUNCTION
!!  Initialize a Gvectors data type
!!
!! NOTES
!!  gvec are supposed to be ordered with increasing norm.

!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine init_Gvectors_type(only_one_kpt,Gsph,Cryst,ng,gvec,gmet,gprimd)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : j_dpc
 use m_errors, only : assert
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => init_Gvectors_type
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng
 logical,intent(in) :: only_one_kpt
 type(Crystal_structure),intent(in) :: Cryst
 type(Gvectors_type),intent(out) :: Gsph
!arrays
 integer,intent(in) :: gvec(3,ng)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ig,istat,isym,nsh,nsym,timrev
 real(dp) :: eps,norm,norm_old
 logical :: ltest
 character(len=500) :: msg
!arrays
 integer :: sg(3)
 integer,allocatable :: shlim(:)
 integer,pointer :: symrec(:,:,:)
 real(dp),allocatable :: shlen(:)
 real(dp),pointer :: tnons(:,:)

!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' init_Gvectors_type : enter '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 ! === Get info on symmetries ===
 nsym   =  Cryst%nsym
 timrev =  Cryst%timrev
 symrec => Cryst%symrec
 tnons  => Cryst%tnons
 !
 ! === Initialize the structure ===
 Gsph%ng     = ng
 Gsph%nsym   = nsym
 Gsph%timrev = timrev

 Gsph%gmet   = gmet
 Gsph%gprimd = gprimd

 allocate(Gsph%gvec(3,ng)) 
 Gsph%gvec(:,:)=gvec(:,:)
 !
 ! === Calculate phase exp{-i2\pi G.\tau} ===
 allocate(Gsph%phmGt(ng,nsym))
 do ig=1,ng
  do isym=1,nsym
   Gsph%phmGt(ig,isym)=EXP(-j_dpc*two_pi*DOT_PRODUCT(gvec(:,ig),tnons(:,isym)))
  end do 
 end do
 !
 ! === Calculate phase phsgt= exp{-i2\pi SG\cdot t} ===
 ! Here we can store only one of this arrays but I have to rewrite screeening!
 allocate(Gsph%phmSGt(ng,nsym))
 do ig=1,ng
  do isym=1,nsym
   sg(:)=MATMUL(symrec(:,:,isym),gvec(:,ig))
   Gsph%phmSGt(ig,isym)=EXP(-j_dpc*two_pi*DOT_PRODUCT(sg,tnons(:,isym)))
  end do
 end do
 !
 ! === Calculate number of shells and relative starting index ===
 ! * Shells are useful to speed up search algorithms see e.g setup_G_rotation.
 ! * The last shell ends at ng+1, thus gvec is supposed to be closed
 ! TODO wrote a checking routine for debugging purpose
 !
 ltest=ALL(gvec(:,1)==0)
 call assert(ltest,'First G should be 0',__FILE__,__LINE__)

 allocate(Gsph%g2sh(ng)) ; Gsph%g2sh(1)=1 ! This table is useful if we dont loop over shell

 ! For each shell, give the index of the initial G-vector 
 allocate(shlim(ng+1))  
 shlim(1)=1 

 ! For each shell, give the radius of the shell
 allocate(shlen(ng))  
 shlen(1)=zero 

 nsh=1 ; norm_old=zero
 do ig=2,ng
  norm=two_pi*SQRT(DOT_PRODUCT(gvec(:,ig),MATMUL(gmet,gvec(:,ig))))
  eps=norm*tol8
  if (ABS(norm-norm_old)>eps) then 
   norm_old=norm
   nsh=nsh+1
   shlim(nsh)=ig
   shlen(nsh)=norm 
  end if
  Gsph%g2sh(ig)=nsh
 end do
 shlim(nsh+1)=ng+1

 ! === Save info on the shells ===
 Gsph%nsh=nsh
 allocate(Gsph%shlim(nsh+1)) ; Gsph%shlim=shlim(1:nsh+1)
 allocate(Gsph%shlen(nsh  )) ; Gsph%shlen=shlen(1:nsh)
 deallocate(shlim,shlen)
 !
 ! === Calculate table for rotated G"s ===
 allocate(Gsph%rottb  (ng,timrev,nsym),STAT=istat)   
 allocate(Gsph%rottbm1(ng,timrev,nsym),STAT=istat) 

 call setup_G_rotation(only_one_kpt,nsym,symrec,timrev,Gsph%ng,Gsph%gvec,&
& Gsph%g2sh,Gsph%nsh,Gsph%shlim,Gsph%rottb,Gsph%rottbm1)

 call print_Gvectors(Gsph,unit=std_out,prtvol=1)

#if defined DEBUG_MODE
 write(msg,'(a)')' init_Gvectors_type : exit '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine init_Gvectors_type
!!***

!!****f* ABINIT/print_Gvectors
!! NAME
!! print_Gvectors
!!
!! FUNCTION
!!  Print the content of a gvectors data type
!!

!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine print_Gvectors(Gsph,unit,prtvol,mode_paral)

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 type(Gvectors_type),intent(in) :: Gsph

!Local variables-------------------------------
!scalars
 integer :: ii,ish,nsc,unt,verbose
 real(dp) :: fact,kin
 character(len=100) :: fmt
 character(len=4) :: mode
 character(len=500) :: msg

! *************************************************************************

 unt=std_out ; if (PRESENT(unit)) unt=unit
 verbose=0   ; if (PRESENT(prtvol)) verbose=prtvol
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 write(msg,'(3a,2(a,i8,a))')ch10,&
& ' ==== Info on the G-sphere ==== ',ch10,&
& '  Number of G vectors ... ',Gsph%ng,ch10,&
& '  Number of shells ...... ',Gsph%nsh,ch10
 call wrtout(unt,msg,mode)

 SELECT CASE (Gsph%timrev) 
 CASE (1)
  write(msg,'(a)')' Time-reversal symmetry is used'
 CASE (2)
  write(msg,'(a)')' Time-reversal symmetry is used'
 CASE DEFAULT
  call assert(.FALSE.,'Wrong value for timrev',__FILE__,__LINE__)
 END SELECT
 call wrtout(unt,msg,mode)

 if (verbose/=0) then 
  fact=half*two_pi**2
  write(msg,'(a)')' Shell   Tot no. of Gs   Cutoff [Ha]'
  call wrtout(unt,msg,mode)
  do ish=1,Gsph%nsh
   nsc=Gsph%shlim(ish+1)-1 
   kin=half*Gsph%shlen(ish)**2
   write(msg,'(2x,i4,10x,i6,5x,f8.3)')ish,nsc,kin
   call wrtout(unt,msg,'COLL')
  end do
  write(msg,'(a)')ch10 ; call wrtout(unt,msg,mode)
 end if

end subroutine print_Gvectors 
!!***

!!****f* ABINIT/destroy_Gvectors
!! NAME
!! destroy_Gvectors
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

subroutine destroy_Gvectors(Gsph)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(Gvectors_type),intent(inout) :: Gsph

! *************************************************************************
 if (ASSOCIATED(Gsph%g2sh   ))  DEALLOCATE(Gsph%g2sh   )
 if (ASSOCIATED(Gsph%gvec   ))  DEALLOCATE(Gsph%gvec   )
 if (ASSOCIATED(Gsph%rottb  ))  DEALLOCATE(Gsph%rottb )
 if (ASSOCIATED(Gsph%rottbm1))  DEALLOCATE(Gsph%rottbm1)
 if (ASSOCIATED(Gsph%shlim  ))  DEALLOCATE(Gsph%shlim  )
 if (ASSOCIATED(Gsph%shlen  ))  DEALLOCATE(Gsph%shlen  )
 if (ASSOCIATED(Gsph%phmGt  ))  DEALLOCATE(Gsph%phmGt  )
 if (ASSOCIATED(Gsph%phmSGt ))  DEALLOCATE(Gsph%phmSGt )

end subroutine destroy_Gvectors
!!***

!!****f* ABINIT/nullify_Gvectors
!! NAME
!! nullify_Gvectors
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

subroutine nullify_Gvectors(Gsph)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(Gvectors_type),intent(inout) :: Gsph

! *************************************************************************
 nullify(Gsph%g2sh   )
 nullify(Gsph%gvec   )
 nullify(Gsph%rottb )
 nullify(Gsph%rottbm1)
 nullify(Gsph%shlim  )
 nullify(Gsph%shlen  )
 nullify(Gsph%phmGt  )
 nullify(Gsph%phmSGt )

end subroutine nullify_Gvectors
!!***
