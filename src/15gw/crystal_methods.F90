!{\src2tex{textfont=tt}}
!!****f* ABINIT/Crystal_methods
!! NAME
!! Crystal_methods
!!
!! FUNCTION
!! This module (?) contains methods operating on the crystal_structure data type.
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
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

!!***

!!****f* ABINIT/init_Crystal_structure
!! NAME
!!  init_Crystal_structure 
!!
!! FUNCTION
!!  Initialize a Crystal_structure data type 
!!
!! INPUTS
!!  natom=number of atom
!!  ntypat=number of type of atoms
!!  nsym=number of symmetry operations
!!  rprimd(3,3)=dimensional lattive vector (real space)
!!  typat(natom)=type of each atom
!!  xred(3,natom)=reduced coordinates of each atom
!!  symrel(3,3) [optional]=symmetry operations in real space
!!  tnons(3,nsym) [optional]=fractional Translations
!!  symafm(nsym) [optional]=  ferromagnetic symmetries
!!  remove_inv [optional]= if .TRUE. the inversion is removed from the set of symmetries
!!  timrev ==2 => take advantage of time-reversal symmetry
!!         ==1 ==> do not use time-reversal symmetry 
!!
!! OUTPUT
!!  Cryst <Crystal_structure>= the data type
!!
!! SOURCE

subroutine init_Crystal_structure(Cryst,natom,ntypat,nsym,rprimd,typat,xred,timrev,use_antiferro,remove_inv,&
& symrel,tnons,symafm) ! optional
    
 use defs_basis
 use defs_datatypes
 use m_numeric_tools, only : set2unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_15gw, except_this_one => init_Crystal_structure
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat,nsym,timrev
 type(Crystal_structure),intent(inout) :: Cryst
 logical,intent(in) :: remove_inv,use_antiferro
!arrays
 integer,intent(in) :: typat(natom)
 integer,optional,intent(in) :: symrel(3,3,nsym),symafm(nsym)
 real(dp),intent(in) :: xred(3,natom),rprimd(3,3)
 real(dp),optional,intent(in) :: tnons(3,nsym)

!Local variables-------------------------------
!scalars
 integer :: iat,indx,itypat,pinv,isym,nsym_noI
 real(dp) :: ucvol
 character(len=500) :: msg      
!arrays
 integer :: symrec(3,3),inversion(3,3)
 real(dp) :: gprimd(3,3),gmet(3,3),rmet(3,3)
 real(dp) :: spinrot(4)
 integer,pointer :: symrel_noI(:,:,:)
 integer,allocatable :: indsym(:,:,:)
 real(dp),pointer :: tnons_noI(:,:)
! *************************************************************************

  Cryst%natom  = natom 
  Cryst%ntypat = ntypat

  call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
  
  Cryst%ucvol  = ucvol
  Cryst%rprimd = rprimd 
  Cryst%rmet   = rmet

  Cryst%angdeg(1)=ACOS(Cryst%rmet(2,3)/SQRT(Cryst%rmet(2,2)*Cryst%rmet(3,3)))/two_pi*360.0d0
  Cryst%angdeg(2)=ACOS(Cryst%rmet(1,3)/SQRT(Cryst%rmet(1,1)*Cryst%rmet(3,3)))/two_pi*360.0d0
  Cryst%angdeg(3)=ACOS(Cryst%rmet(1,2)/SQRT(Cryst%rmet(1,1)*Cryst%rmet(2,2)))/two_pi*360.0d0

  allocate(Cryst%typat(natom),Cryst%xred(3,natom),Cryst%xcart(3,natom)) 
  Cryst%typat=typat 
  Cryst%xred=xred 
  call xredxcart(natom,1,rprimd,Cryst%xcart,Cryst%xred)
  !
  ! === Generate an index table of atoms, in order for them to be used type after type ===
  allocate(Cryst%atindx(natom),Cryst%atindx1(natom),Cryst%nattyp(ntypat))
  indx=1
  do itypat=1,ntypat
   Cryst%nattyp(itypat)=0
   do iat=1,natom
    if (Cryst%typat(iat)==itypat) then
     Cryst%atindx (iat )=indx 
     Cryst%atindx1(indx)=iat
     indx=indx+1
     Cryst%nattyp(itypat)=Cryst%nattyp(itypat)+1
    end if
   end do
  end do

  Cryst%timrev = timrev
  call set2unit(inversion) ; inversion=-inversion

  if (PRESENT(symrel).and.PRESENT(tnons).and.PRESENT(symafm)) then 
   if (.not.remove_inv) then
    ! * Just a copy 
    Cryst%nsym= nsym
    allocate(Cryst%symrel(3,3,nsym),Cryst%symrec(3,3,nsym))
    allocate(Cryst%tnons(3,nsym),Cryst%symafm(nsym))
    Cryst%symrel=symrel 
    Cryst%tnons=tnons
    Cryst%symafm=symafm
    Cryst%use_antiferro = use_antiferro
    Cryst%has_inversion=.FALSE.
    do isym=1,nsym
     call mati3inv(symrel(:,:,isym),symrec)
     Cryst%symrec(:,:,isym)=symrec
     if (ALL(symrel(:,:,isym)==inversion)) Cryst%has_inversion=.TRUE. 
    end do
   else 
    ! * Remove inversion, just to be compatible with old GW implementation
    call remove_inversion(nsym,symrel,tnons,nsym_noI,symrel_noI,tnons_noI,pinv)
    Cryst%nsym=nsym_noI
    allocate(Cryst%symrel(3,3,nsym_noI),Cryst%symrec(3,3,nsym_noI))
    allocate(Cryst%tnons(3,nsym_noI),Cryst%symafm(nsym_noI))
    Cryst%symrel=symrel_noI
    Cryst%tnons=tnons_noI
    Cryst%has_inversion=.FALSE.
    if (ANY(symafm==-1)) then 
     write(msg,'(a)')' First of all solve the problem with inversion before adding ferromagnetic symmetries '
     call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
    end if
    Cryst%symafm=1
    Cryst%use_antiferro=use_antiferro 
    do isym=1,nsym_noI
     call mati3inv(symrel_noI(:,:,isym),symrec)
     Cryst%symrec(:,:,isym)=symrec
    end do
    deallocate(symrel_noI,tnons_noI)
   end if

  else
   ! * Find symmetries symrec,symrel,tnons,symafm
   ! TODO This should be a wrapper around the abinit library whose usage is not so straightforward
  end if

  ! === Obtain a list of rotated atom labels ===
  allocate(indsym(4,Cryst%nsym,natom))
  call symatm(indsym,natom,Cryst%nsym,Cryst%symrec,Cryst%tnons,Cryst%typat,Cryst%xred)

  allocate(Cryst%indsym(4,Cryst%nsym,natom))
  Cryst%indsym=indsym  
  deallocate(indsym)

  ! === Compute rotation in spinor space ===
  allocate(Cryst%spinrot(4,Cryst%nsym))
  do isym=1,Cryst%nsym
   call getspinrot(Cryst%rprimd,spinrot,Cryst%symrel(:,:,isym))
   Cryst%spinrot(:,isym)=spinrot(:)
  end do

end subroutine init_Crystal_structure
!!***

!!****f* ABINIT/init_Crystal_from_Hdr 
!! NAME
!!  init_Crystal_from_Hdr
!!
!! FUNCTION
!!  initialize a Crystal_structure data type starting from the abinit header
!!
!! INPUTS
!!  Hdr<Hdr_type>=the abinit header
!!  timrev ==2 => take advantage of time-reversal symmetry
!!         ==1 ==> do not use time-reversal symmetry 
!!  remove_inv [optional]= if .TRUE. the inversion symmetry is removed from the set of operations
!!  even though it is present in the header
!!
!! OUTPUT
!!  Cryst <Crystal_structure>= the data type filled with data reported in the abinit header 
!!
!! TODO
!!  Add information on the use of time-reversal in the Abinit header.
!!
!! SOURCE

subroutine init_Crystal_from_Hdr(Cryst,Hdr,timrev,remove_inv)

 use defs_basis 
 use defs_datatypes
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15gw, except_this_one => init_Crystal_from_Hdr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(Hdr_type),intent(in) :: Hdr
 type(Crystal_structure),intent(out) :: Cryst 
 integer,intent(in) :: timrev
 logical,optional,intent(in) :: remove_inv

!Local variables-------------------------------
 logical :: rinv,ltest,use_antiferro
! *********************************************************************

 rinv=.FALSE. ; if (PRESENT(remove_inv)) rinv=remove_inv
 use_antiferro=(Hdr%nspden==2.and.Hdr%nsppol==1)

 ! === consistency check ===
 ltest = (timrev==1.or.timrev==2)
 call assert(ltest,"Wrong value for timrev (1|2)",__FILE__,__LINE__)
 if (use_antiferro) then
  ltest = (ANY(Hdr%symafm==-1))
  call assert(ltest,"Wrong combination of nspden, nsppol, symafm.",__FILE__,__LINE__)
 end if

 call init_Crystal_structure(Cryst,Hdr%natom,Hdr%ntypat,Hdr%nsym,Hdr%rprimd,Hdr%typat,Hdr%xred,timrev,use_antiferro,rinv,&
& Hdr%symrel,Hdr%tnons,Hdr%symafm) ! Optional

end subroutine init_Crystal_from_Hdr
!!***

!!****f* ABINIT/nullify_Crystal_structure 
!! NAME
!! nullify_Crystal_structure 
!!
!! FUNCTION
!!  Nullify the dinamic pointers in a Crystal_structure data type
!!
!! SOURCE

subroutine nullify_Crystal_structure(Cryst)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
 type(Crystal_structure),intent(inout) :: Cryst
! *********************************************************************

! integer
 NULLIFY(Cryst%indsym )
 NULLIFY(Cryst%symafm )
 NULLIFY(Cryst%symrec )
 NULLIFY(Cryst%symrel )
 NULLIFY(Cryst%atindx )   
 NULLIFY(Cryst%atindx1)   
 NULLIFY(Cryst%typat  )   
 NULLIFY(Cryst%nattyp )   

 NULLIFY(Cryst%tnons  )
 NULLIFY(Cryst%xcart  )
 NULLIFY(Cryst%xred   )
 NULLIFY(Cryst%spinrot)

end subroutine nullify_Crystal_structure
!!***

!!****f* ABINIT/destroy_Crystal_structure
!! NAME
!!  destroy_Crystal_structure 
!!
!! FUNCTION
!!  Destroy the dinamic pointers in a Crystal_structure data type
!!
!! SOURCE

subroutine destroy_Crystal_structure(Cryst)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(Crystal_structure),intent(inout) :: Cryst

!Local variables-------------------------------
 character(len=500) :: msg
! *********************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')'destroy_Crystal_structure : enter ' 
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

!integer
 if (ASSOCIATED(Cryst%indsym )) DEALLOCATE(Cryst%indsym )   
 if (ASSOCIATED(Cryst%symafm )) DEALLOCATE(Cryst%symafm )   
 if (ASSOCIATED(Cryst%symrec )) DEALLOCATE(Cryst%symrec )   
 if (ASSOCIATED(Cryst%symrel )) DEALLOCATE(Cryst%symrel )   
 if (ASSOCIATED(Cryst%atindx )) DEALLOCATE(Cryst%atindx )   
 if (ASSOCIATED(Cryst%atindx1)) DEALLOCATE(Cryst%atindx1)   
 if (ASSOCIATED(Cryst%typat  )) DEALLOCATE(Cryst%typat  )   
 if (ASSOCIATED(Cryst%nattyp )) DEALLOCATE(Cryst%nattyp )   

!real
 if (ASSOCIATED(Cryst%tnons  )) DEALLOCATE(Cryst%tnons  )   
 if (ASSOCIATED(Cryst%xcart  )) DEALLOCATE(Cryst%xcart  )   
 if (ASSOCIATED(Cryst%xred   )) DEALLOCATE(Cryst%xred   )   
 if (ASSOCIATED(Cryst%spinrot)) DEALLOCATE(Cryst%spinrot)   

#if defined DEBUG_MODE
 write(msg,'(a)')'destroy_Crystal_structure : exit ' 
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine destroy_Crystal_structure
!!***

!!****f* ABINIT/print_Crystal_structure
!! NAME
!!  print_Crystal_structure 
!!
!! FUNCTION
!!  Print the content of Crystal_structure data type
!!
!! SOURCE

subroutine print_Crystal_structure(Cryst,unit,mode_paral,prtvol) 

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => print_Crystal_structure
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral 

!Local variables-------------------------------
 integer :: unt,verb,nu,isym
 character(len=4) :: mode
 character(len=500) :: msg      
! ********************************************************************* 

 unt=std_out ; if (PRESENT(unit))       unt=unit
 verb=0      ; if (PRESENT(prtvol))     verb=prtvol 
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 write(msg,'(2a)')' Real(R) ',&
&  ' space primitive vectors, cartesian coordinates (Bohr):'
 call wrtout(unt,msg,mode)
 do nu=1,3
  write(msg,'(1x,a,i1,a,3f11.7)')&
&  'R(',nu,')=',Cryst%rprimd(:,nu)+tol10  !tol10 is used to be consistent with metric.F90
  call wrtout(unt,msg,mode)
 end do

 write(msg,'(a,1p,e15.7,a)')&
& ' Unit cell volume ucvol=',Cryst%ucvol+tol10,' bohr^3'
 call wrtout(unt,msg,mode)

 write(msg,'(a,3es16.8,a)')&
& ' Angles (23,13,12)=',Cryst%angdeg(1:3),' degrees'
 call wrtout(unt,msg,mode)

 if (Cryst%timrev==1) then 
  write(msg,'(a)')' Time-reversal symmetry is not present '
 else if (Cryst%timrev==2) then 
  write(msg,'(a)')' Time-reversal symmetry is present '
 else 
  call assert(.FALSE.,'Wrong value for timrev (1|2)',__FILE__,__LINE__)
 end if
 call wrtout(unt,msg,mode)

 !if (Cryst%use_antiferro) then 
 ! write(msg,'(a)')' System has magnetic symmetries '
 ! call wrtout(unt,msg,mode)
 !end if

 call print_symmetries(Cryst%nsym,Cryst%symrel,Cryst%tnons,Cryst%symafm,unit=unt,mode_paral=mode)

end subroutine print_Crystal_structure
!!***

subroutine print_symmetries(nsym,symrel,tnons,symafm,unit,mode_paral)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,optional,intent(in) :: unit
 character(len=4),optional,intent(in) :: mode_paral
!arrays
 integer,intent(in) :: symrel(3,3,nsym),symafm(nsym)
 real(dp),intent(in) :: tnons(3,nsym)

!Local variables-------------------------------
 integer :: unt,isym
 character(len=500) :: msg      
 character(len=4) :: mode
! *********************************************************************

 unt=std_out ; if (PRESENT(unit)) unt=unit
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 write(msg,'(2a)')ch10,' Rotations                           Translations     Symafm '
 do isym=1,nsym
  write(msg,'(1x,3(3i3,1x),4x,3(f11.7,1x),6x,i2)')&
&  symrel(:,:,isym),tnons(:,isym),symafm(isym)
  call wrtout(unt,msg,mode)
 end do 

end subroutine print_symmetries 
!!***

!!****f* ABINIT/is_equal_Crystal_structure
!! NAME
!!  is_equal_Crystal_structure 
!!
!! FUNCTION
!!  Compare two Crystal_structure data types 
!!
!! SOURCE

function is_equal_Crystal_structure(Cryst1,Cryst2) result(equal)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
 type(Crystal_structure),intent(in) :: Cryst1,Cryst2 
 integer :: equal

!!Local variables-------------------------------
 real(dp) :: TOL=tol6
 ! *********************************************************************

!TODO check dimensions
 equal=0
 if (Cryst1%natom/=Cryst2%natom .or. Cryst1%ntypat/=Cryst2%ntypat)  equal=1
 !
 ! Check integer pointers
 if ( &
&    ANY (Cryst1%symafm/=Cryst2%symafm) .or.& 
&    ANY (Cryst1%symrec/=Cryst2%symrec) .or.& 
&    ANY (Cryst1%symrel/=Cryst2%symrel)     & 
    ) equal=1

 if (Cryst1%timrev /= Cryst2%timrev) equal=2

 ! Check real arrays
 if ( &
&    ANY (ABS(Cryst1%rprimd-Cryst2%rprimd)>TOL) .or.& 
&    ANY (ABS(Cryst1%rmet-Cryst2%rmet)>TOL) .or.& 
&    ANY (ABS(Cryst1%angdeg-Cryst2%angdeg)>TOL) .or. &
&    ANY (ABS(Cryst1%tnons-Cryst2%tnons)>TOL) &
    ) equal=3

end function is_equal_Crystal_structure 
!!***

!!****f* ABINIT/remove_inversion
!! NAME
!!  remove_inversion
!!
!! FUNCTION
!!  Remove the inversion symmetry from a symmetry set as well 
!!  all the improper rotations (if present)
!!
!! INPUTS
!!  nsym=initial number of symmetries
!!  symrel(3,3,nsym)=Initial set of symmetry operarations in real space
!!  tnons(3,nsym)=Initial fractional translations
!!
!! OUTPUT
!!  nsym_out=Number of symmetries in the set without improper rotation
!!  symrel_out(:,:) [pointer] = output symmetries without improper rotations
!!  tnons_out(:) [pointer] = fraction translations associated to symrel_out
!!  pinv=-1 if the inversion has been removed, 1 otherwise
!!
!! NOTES
!!  Note the use of pointers, memory is allocated inside the procedure and passed back 
!!  to the caller. Thus memory deallocation is relegated to the caller. To be on the safe side
!!  the pointers should be nulified before entering.
!!
!! SOURCE

subroutine remove_inversion(nsym,symrel,tnons,nsym_out,symrel_out,tnons_out,pinv)

 use defs_basis
 use m_numeric_tools, only : is_integer, set2unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(out) :: nsym_out,pinv
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 integer,pointer :: symrel_out(:,:,:)
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),pointer :: tnons_out(:,:)

!Local variables-------------------------------
!scalars
 integer :: i,j,is,is2,is_inv,nsym2
 integer :: is_retained,is_discarded
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: symrel2(3,3,nsym),inversion(3,3),determinant(nsym)
 real(dp) :: tnons2(3,nsym),dtnons(3)
! *********************************************************************

  write(msg,'(4a)')ch10,&
&  ' remove_inversion: WARNING - ',ch10,&
&  ' Removing inversion related symmetrie from initial set '
 call wrtout(std_out,msg,'COLL')
 !
 ! === Find the occurence of the inversion symmetry ===
 call set2unit(inversion) ; inversion=-inversion

 is_inv=0 ; found=.FALSE.
 do while (is_inv<nsym .and. .not.found)
  is_inv=is_inv+1 ; found=ALL(symrel(:,:,is_inv)==inversion)
 end do
 if (found) then
  write(msg,'(a,i3)')' The inversion is symmetry operation no. ',is_inv
 else
  write(msg,'(a)')' The inversion was not found in the symmetries list.'
 end if
 call wrtout(std_out,msg,'COLL')
 !
 ! === Find the symmetries that are related through the inversion symmetry ===
 call symdet(determinant,nsym,symrel)
 nsym2=0
 do is=1,nsym-1
  do is2=is+1,nsym

   dtnons(:)=tnons(:,is2)-tnons(:,is)-tnons(:,is_inv)
   found=ALL(symrel(:,:,is)==-symrel(:,:,is2)).and.is_integer(dtnons,tol8)

   if (found) then
    nsym2=nsym2+1
    ! * Retain symmetries with positive determinant
    if (ALL(tnons(:,is2)<tol8).and.ALL(tnons(:,is)<tol8)) then
     is_retained=is2 ; is_discarded=is
     if (determinant(is)==1) then  
      is_retained=is  ; is_discarded=is2
     end if
    else if (ALL(tnons(:,is2)<tol8)) then
     is_retained=is2 ; is_discarded=is
    else
     is_retained=is ;  is_discarded=is2
    end if

    symrel2(:,:,nsym2)=symrel(:,:,is_retained)
    tnons2   (:,nsym2)=tnons   (:,is_retained)
    write(msg,'(a,i3,a,i3,3a,i3,a)')&
&    ' Symmetry operations no. ',is,' and no. ',is2,&
&    ' are related through the inversion.',ch10,&
&    ' Symmetry operation no. ',is_discarded,' will be suppressed.'
    call wrtout(std_out,msg,'COLL')
   end if ! found

  end do !is2
 end do !is

 if (nsym2/=(nsym/2).or.nsym==1) then
  write(msg,'(a)')' Program uses the original set of symmetries '
  call wrtout(std_out,msg,'COLL')
  nsym_out=nsym
  allocate(symrel_out(3,3,nsym),tnons_out(3,nsym))
  symrel_out(:,:,:)=symrel(:,:,1:nsym)
  tnons_out(:,:)=tnons(:,1:nsym)
  pinv=1
 else
  write(msg,'(a)')' Inversion related operations have been suppressed from symmetries list.'
  call wrtout(std_out,msg,'COLL')
  nsym_out=nsym2
  allocate(symrel_out(3,3,nsym2),tnons_out(3,nsym2))
  symrel_out(:,:,:)=symrel2(:,:,1:nsym2)
  tnons_out(:,:)=tnons(:,1:nsym2)
  pinv=-1
 end if

end subroutine remove_inversion
!!***
