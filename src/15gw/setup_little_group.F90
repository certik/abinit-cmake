!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_little_group
!! NAME
!! setup_little_group
!!
!! FUNCTION
!!  Finds symmetry operations belonging to the little group associated to an external 
!!  point ext_pt and fills symmetry tables.
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! gmet(3,3)=reciprocal space metric (bohr**-2).
!! gvec(3,npwvec) coordinates of G vectors
!! npwvec=number of G vectors
!! Kmesh<BZ_mesh_type> 
!!   %nbz=number of points in the full BZ
!!   %kbz(3,nbz)=points in the full Brillouin Zone
!! Cryst<Crystal_structure>= Info on symmetries and unit cell
!!   %symrec(3,3,nsym)=symmetry operations in reciprocal space
!!   %nsym=number of symmetry operations in the space group
!!   %timrev=if 2 time-reversal can be used; 1 otherwise
!! prtvol=input variable defining the verbosity
!! ext_pt(3)= External point in the Brillouin zone in reduce coordinated 
!! use_umklp=flag to include umklapp G0 vectors in the definition of the little group (0:n0,1:yes)
!!
!! OUTPUT
!! Ltg% <little_group_datatype>.
!!  %ibzq(nbz)= 1 if the kpoint belongs to the IBZ defined by ext_pt, 0 otherwise
!!  %bz2ibz(nbz)= sequential index of the point in the IBZ defined by ext_pt
!!  %ibz2bz(nibz_Ltg) For each nibz_Ltg the correspondind index in the BZ array
!!  %igmG0(npwepstimrev,nsym)= index of the uklapp vector G_o in the FFT array
!!  %flag_umklp(timrev,nsym)= flag for umklapp processes 
!!    1 if operation (IS) requires a G_o to preserve ext_pt, 0 otherwise 
!!  %tab(nbz)=table giving, for each k-point in the BZ (kBZ), the corresponding 
!!   irreducible point (kIBZ) in the irreducible zone defined by the little group of ext_pt,
!!   i.e kBZ= (IS) kIBZ where I is either the inversion or the identity and S is an
!!   operation in the little group defined by ext_pt
!!  %tabo(nbz)=the symmetry operation S in the little group that takes kIBZ to each kBZ
!!  %tabi(nbz)= defines whether inversion has to be considered in the 
!!   relation kBZ=(IS) kIBZ (1 => only S; -1 => -S)  
!!  %preserve(2,nsym)= 1 if ISq=q, 0 otherwise, the first index is for the identity or the time reversal symmetry, 
!!  %wtksym(2,nsym,nbz)= for each kpoint is equal to 1 if the symmetry operation (with or without time reversal)  
!!   must be considered in the calculation of \chi_o, 0 otherwise  
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
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

subroutine setup_little_group(ext_pt,Kmesh,gmet,Cryst,npwvec,gvec,npwe,use_umklp,prtvol,Ltg)

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_13recipspace
 use interfaces_15gw, except_this_one => setup_little_group
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwe,npwvec,prtvol,use_umklp
 type(Crystal_structure),intent(in) :: Cryst
 type(bz_mesh_type),intent(in) :: Kmesh
 type(little_group),intent(inout) :: Ltg
!arrays
 integer,intent(in) :: gvec(3,npwvec)
 real(dp),intent(in) :: ext_pt(3),gmet(3,3)

!Local variables-------------------------------
!scalars
 integer :: dummy_timrev,enough,idx,ig,ige,igpw,ik,ikp,ind,iold,isym,itest,itim
 integer :: nbz,nkibzq,nsym,nsym_Ltg,ntest,option,timrev
 real(dp) :: G0len,kin,mG0len,max_kin
 logical :: found,found_identity,ltest,use_antiferro
 character(len=500) :: msg
!arrays
 integer :: g0(3),gg(3),gmG0(3),identity(3,3),nop(Cryst%timrev),nopg0(2)
 integer :: symxpt(4,2,Cryst%nsym)
 integer,allocatable :: indkpt1(:),symafm_ltg(:),symrec_Ltg(:,:,:)
 integer,pointer :: symafm(:),symrec(:,:,:)
 real(dp) :: knew(3)
 real(dp),allocatable :: ktest(:,:),wtk(:),wtk_folded(:)
 real(dp),pointer :: kbz(:,:)

!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' setup_little_group : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif 
 !
 ! === Initial check ====
 ltest=(Cryst%timrev==1.or.Cryst%timrev==2)
 write(msg,'(a,i3)')'Wrong value for timrev :',Cryst%timrev 
 call assert(ltest,msg,__FILE__,__LINE__)
 !
 ! === Get useful data ===
 nsym   = Cryst%nsym
 timrev = Cryst%timrev
 symrec => Cryst%symrec
 symafm => Cryst%symafm
 use_antiferro = Cryst%use_antiferro

 nbz = Kmesh%nbz 
 kbz => Kmesh%bz(1:3,1:nbz)
 !
 ! === Destroy structure if it already exists ===
 call destroy_little_group(Ltg)
 !
 ! === Store dimensions and useful info ===
 Ltg%nsym_sg  =nsym
 Ltg%timrev   =timrev
 Ltg%nbz      =nbz
 !Ltg%use_umklp=use_umklp ! 1 if umklapp processes are used
 Ltg%ext_pt(:)=ext_pt(:)

 allocate(Ltg%G0(3,timrev,nsym)) 
 allocate(Ltg%ibzq(nbz),Ltg%bz2ibz(nbz))   
 allocate(Ltg%preserve(timrev,nsym),Ltg%wtksym(timrev,nsym,nbz))
 allocate(Ltg%tab(nbz),Ltg%tabi(nbz),Ltg%tabo(nbz))
 allocate(Ltg%flag_umklp(timrev,nsym))
 !
 ! In the old GW implementation we were removing symmetries related by time-reversal and 
 ! sometimes it happened that only the inversion was reported in the KSS file (see outkss.F90).
 identity(:,:)=RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/)) ; found_identity=.FALSE. 
 do isym=1,nsym
  if (ALL(symrec(:,:,isym)==identity)) then 
   found_identity=.TRUE. ; EXIT
  end if
 end do 
 if (.not.found_identity) then 
  write(msg,'(8a)')ch10,&
&  ' setup_little_group : ERROR -',ch10,&
&  ' Only the inversion was found in the set of symmetries read from the KSS file ',ch10,&
&  ' Likely you are using a KSS file generated with an old version of Abinit, ',ch10,&
&  ' To run a GW calculation with an old KSS file, use version < 5.5 '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if 
 !
 ! === Find operations in the little group as well as umklapp vectors G0 === 
 call symq3(nsym,ext_pt,symxpt,symrec,dummy_timrev) 

 Ltg%preserve(:,:)=0 ; Ltg%g0(:,:,:)=0 ; Ltg%flag_umklp(:,:)=0 ; mG0len=zero
 do itim=1,timrev
  do isym=1,nsym

   if (use_antiferro.and.symafm(isym)==-1) CYCLE

   if (symxpt(4,itim,isym)==1) then  !\pm Sq = q+g0
    if (ANY(symxpt(1:3,itim,isym)/=0).and.use_umklp==0) CYCLE ! Exclude non zero G0 vectors
    Ltg%preserve(itim,isym)=1
    g0(:)=symxpt(1:3,itim,isym) ; Ltg%g0(:,itim,isym)=g0(:)
    if (ANY(Ltg%g0(:,itim,isym)/=0)) Ltg%flag_umklp(itim,isym)=1
    ! Max radius to be considered to include all G0s
    G0len=normv(DBLE(g0),gmet,'g') 
    mG0len=MAX(mG0len,G0len) 
   end if

  end do
 end do 

 nop(:)=0 ; nopg0(:)=0
 do itim=1,timrev
  nop  (itim)=SUM(Ltg%preserve  (itim,:))
  nopg0(itim)=SUM(Ltg%flag_umklp(itim,:))
 end do
 nsym_Ltg=SUM(nop(:))
 !
 ! === Store little group operations, include time-reversal if present ===
 Ltg%nsym_Ltg=nsym_Ltg
 allocate(symrec_Ltg(3,3,Ltg%nsym_Ltg))  

 ind=1
 do itim=1,timrev
  do isym=1,nsym
   if (Ltg%preserve(itim,isym)==1) then  
    if (itim==1) symrec_Ltg(:,:,ind)= symrec(:,:,isym)
    if (itim==2) symrec_Ltg(:,:,ind)=-symrec(:,:,isym)
    ind=ind+1
   end if
  end do 
 end do
 !
 ! === Check the closure of the (ferromagnetic) little group ===
 allocate(symafm_ltg(Ltg%nsym_Ltg)) ; symafm_ltg(:)=1
 call chkgrp(Ltg%nsym_Ltg,symafm_ltg,symrec_Ltg)
 deallocate(symafm_ltg)
 !
 ! === Find the irreducible zone associated to ext_pt ===
 ! * Do not use time-reversal since it has been manually introduced previously
 allocate(indkpt1(nbz),wtk_folded(nbz),wtk(nbz))
 wtk(:)=one ; option=0 ; dummy_timrev=0
 call symkpt(gmet,indkpt1,kbz,nbz,nkibzq,Ltg%nsym_Ltg,option,symrec_Ltg,dummy_timrev,wtk,wtk_folded)
 deallocate(indkpt1,wtk) 

 Ltg%nibz_Ltg=nkibzq
 !
 ! === Set up table in the BZ === 
 ! * 0 if the point does not belong to IBZ_xpt, 1 otherwise
 allocate(Ltg%ibz2bz(nkibzq)) 
 Ltg%ibzq(:)=0 ; Ltg%bz2ibz(:)=0 ; Ltg%ibz2bz(:)=0

 ind=0 ; enough=0
 do ik=1,nbz
  if (wtk_folded(ik)>tol8) then
   ind=ind+1
   Ltg%ibzq(ik)=1
   Ltg%bz2ibz(ik) =ind
   Ltg%ibz2bz(ind)=ik
  end if
 end do
 if (ind/=Ltg%nibz_Ltg) STOP " BUG ind/=Ltg%nibz_Ltg"
 !
 ! === Reconstruct full BZ starting from IBZ_q  ===
 ! Calculate appropriate weight for each item (point,symmetry operation,time-reversal)
 Ltg%tab=0 ; Ltg%tabo=0 ; Ltg%tabi=0 ; Ltg%wtksym(:,:,:)=0

 ! === Zero no. of k-points found ===
 ntest=0 ; allocate(ktest(3,nbz)) ; ktest=zero  

 do ik=1,nbz
  if (Ltg%ibzq(ik)/=1) CYCLE
  ! * Loop over symmetry operations S and time-reversal.
  do isym=1,nsym
   do itim=1,timrev

    ! Form IS k only for (IS) pairs in the (ferromagnetic) little group.
    if (use_antiferro.and.symafm(isym)==-1) CYCLE
    if (Ltg%preserve(itim,isym)==0) CYCLE
    knew(:)=(3-2*itim)*MATMUL(symrec(:,:,isym),kbz(:,ik))
    !
    ! === Check whether it has already been found (to within a RL vector) ===
    iold=0
    do itest=1,ntest
     if (is_samek(knew(:),ktest(:,itest),gg)) iold=iold+1
    end do

    if (iold==0) then 
     ! == Found new BZ point ===
     ! For this point the operation (isym,itim) must be considered to reconstruct the full BZ
     Ltg%wtksym(itim,isym,ik)=1
     ntest=ntest+1
     ktest(:,ntest)=knew(:)
     ! 
     ! === Now find knew in the BZ array ===
     found=.FALSE.
     do idx=1,nbz
      if (is_samek(knew(:),kbz(:,idx),gg)) then 
       ! They are the same within a RL vector
       Ltg%tab (idx)=ik
       Ltg%tabo(idx)=isym
       Ltg%tabi(idx)=3-2*itim
       found=.TRUE. ; EXIT
      end if 
     end do
     if (.not.found) then 
      write(msg,'(4a,3f12.6,a)')ch10,&
&      ' setup_little_group : BUG - ',ch10,&
&      ' not able to find the ',knew(:),' in the array BZ '
      call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
     end if
    end if

   end do
  end do
 end do
 deallocate(ktest)

 if (ntest/=nbz) then 
  write(msg,'(4a,i5)')ch10,&
&  ' setup_little_group : BUG - ',ch10,&
&  ' ntest-nbz = ',ntest-nbz
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if 
 if (SUM(Ltg%wtksym)/=nbz) then 
  write(msg,'(4a,i5)')ch10,&
&  ' setup_little_group : BUG - ',ch10,&
&  ' sum(Ltg%wtksym)-nbz = ',SUM(Ltg%wtksym)-nbz
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if 

 nullify(Ltg%igmG0) ; Ltg%max_kin_gmG0=zero
 if (npwe>0) then 
  ! This correspond to the case in which we need to know the index of G-Go in the gvec array
  ! where G is one of the npwe vectors. This is required in screenig but not in sigma.
  ! The drawback is that the effective G sphere used to calculate the oscillators must be smaller 
  ! that gvec if we want to avoid possible aliasing effects. Lifting this constraint would require 
  ! a lot of boring coding. (no need to do this if ext_pt=zero, but oh well)
  allocate(Ltg%igmG0(npwe,timrev,nsym)) ; Ltg%igmG0(:,:,:)=0
  max_kin=zero
  !
  ! === Loop over symmetry operations S and time-reversal ===
  do itim=1,timrev
   do isym=1,nsym

    ! * Form IS k only for (IS) pairs in the little group
    if (use_antiferro.and.symafm(isym)==-1) CYCLE
    if (Ltg%preserve(itim,isym)/=0) then 
     g0(:)=Ltg%g0(:,itim,isym)
     do ige=1,npwe
      gmG0(:)=gvec(:,ige)-g0(:)
      kin=half*normv(DBLE(gmG0),gmet,'g')**2 
      max_kin=MAX(max_kin,kin)

      found=.FALSE.
      do igpw=1,npwvec
       if (ALL(gvec(:,igpw)-gmG0(:)==0)) then 
        Ltg%igmG0(ige,itim,isym)=igpw
        found=.TRUE. ; EXIT
       end if
      end do
      if (.not.found) then 
       write(msg,'(8a,f8.3,2a,3i5,2a)')ch10,&
&       ' setup_little_group : ERROR ',ch10,&
&       ' Not able to found G-G0 in the largest G-spere ',ch10,&
&       ' Decrease the size of epsilon or, if possible, increase ecutwfn (>ecuteps) ',ch10,&
&       ' Minimum required cutoff energy for G-G0 sphere = ',kin,ch10,&
&       ' G0 = ',g0(:)
       call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
      end if 

     end do 
    end if
   end do 
  end do 
  Ltg%max_kin_gmG0=max_kin
 end if
 deallocate(symrec_Ltg)

#if defined DEBUG_MODE
 do ik=1,nbz
  if (ABS(SUM(Ltg%wtksym(1,:,ik)+Ltg%wtksym(2,:,ik))-wtk_folded(ik))>tol6) then 
   write(*,*)' sum(Ltg%wtksym,ik)-wtk_folded(ik) = ',sum(Ltg%wtksym(1,:,ik)+Ltg%wtksym(2,:,ik))-wtk_folded(ik)
   write(*,*)Ltg%wtksym(1,:,ik),Ltg%wtksym(2,:,ik),wtk_folded(ik)
   write(*,*)ik,kbz(:,ik) 
   call leave_new('COLL')
  end if 
 end do
 do ik=1,nbz 
  call dosym(REAL(symrec(:,:,Ltg%tabo(ik)),dp),(3-Ltg%tabi(ik))/2,kbz(:,Ltg%tab(ik)),knew(:))
  if (.not.is_samek(knew,kbz(:,ik),gg)) then 
   write(*,*)knew,kbz(:,ik) 
   write(*,*)Ltg%tabo(ik),Ltg%tabi(ik),Ltg%tab(ik)
   write(*,*)' BUG in setup_little_group'
   call leave_new('COLL')
  end if 
 end do
 write(msg,'(a)')' setup_little_group : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 deallocate(wtk_folded)

end subroutine setup_little_group
!!***


!!****if* ABINIT/nullify_little_group
!! NAME
!! nullify_little_group
!!
!! FUNCTION
!!
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

subroutine nullify_little_group(Ltg)

 use defs_basis
 use defs_datatypes 

 implicit none

!Arguments ------------------------------------
!scalars
 type(Little_group),intent(inout) :: Ltg

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *********************************************************************
 
 nullify(Ltg%g0)
 nullify(Ltg%ibzq)
 nullify(Ltg%bz2ibz)
 nullify(Ltg%ibz2bz)
 nullify(Ltg%igmG0)
 nullify(Ltg%flag_umklp)
 nullify(Ltg%preserve)          
 nullify(Ltg%tab)
 nullify(Ltg%tabo)
 nullify(Ltg%tabi)
 nullify(Ltg%wtksym)

end subroutine nullify_little_group
!!***


!!****if* ABINIT/destroy_little_group
!! NAME
!! destroy_little_group
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

subroutine destroy_little_group(Ltg)

 use defs_basis
 use defs_datatypes 

 implicit none

!Arguments ------------------------------------
!scalars
 type(Little_group),intent(inout) :: Ltg

! *********************************************************************
 if (associated(Ltg%g0      ))    deallocate(Ltg%g0        )
 if (associated(Ltg%ibzq    ))    deallocate(Ltg%ibzq      )
 if (associated(Ltg%bz2ibz  ))    deallocate(Ltg%bz2ibz    )
 if (associated(Ltg%ibz2bz  ))    deallocate(Ltg%ibz2bz    )
 if (associated(Ltg%igmG0    ))   deallocate(Ltg%igmG0     )
 if (associated(Ltg%flag_umklp))  deallocate(Ltg%flag_umklp)
 if (associated(Ltg%preserve))    deallocate(Ltg%preserve  )   
 if (associated(Ltg%tab     ))    deallocate(Ltg%tab       )
 if (associated(Ltg%tabo    ))    deallocate(Ltg%tabo      )
 if (associated(Ltg%tabi    ))    deallocate(Ltg%tabi      )
 if (associated(Ltg%wtksym  ))    deallocate(Ltg%wtksym    )
end subroutine destroy_little_group
!!***


!!****if* ABINIT/print_little_group
!! NAME
!! print_little_group
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

subroutine print_little_group(Ltg,unit,prtvol,mode_paral)

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
 type(Little_group),intent(in) :: Ltg

!Local variables-------------------------------
!scalars
 integer :: ik,itim,unt,verbose
 character(len=4) :: mode
 character(len=500) :: msg
!arrays
 integer :: nop(Ltg%timrev),nopg0(Ltg%timrev)

! *********************************************************************

 unt=std_out ; if (PRESENT(unit)) unt=unit
 verbose=0   ; if (PRESENT(prtvol)) verbose=prtvol
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 write(msg,'(3a,3es16.8,2a,i5,a,i5,2a,i3,a,i3)')&
& ' ==== Little group Info ==== ',ch10,&
& '  External point ',Ltg%ext_pt(:),ch10,&
& '  Number of points in the IBZ defined by little group  ',Ltg%nibz_Ltg,'/',Ltg%nbz,ch10,&
& '  Number of operations in the little group : ',Ltg%nsym_Ltg,'/',Ltg%nsym_sg
 call wrtout(unt,msg,mode)

 nop(:)=0 ; nopg0(:)=0
 do itim=1,Ltg%timrev
  nop  (itim)=SUM(Ltg%preserve  (itim,:))
  nopg0(itim)=SUM(Ltg%flag_umklp(itim,:))
 end do
 do itim=1,Ltg%timrev
  if (itim==1) then 
   write(msg,'(a,2(a,i2,a))')ch10,&
&   '  No time-reversal symmetry +     zero Umklapp vector : ',nop(1)-nopg0(1),ch10,&
&   '  No time-reversal symmetry + non-zero Umklapp vector : ',nopg0(1),ch10
   call wrtout(unt,msg,mode)
  else if (itim==2) then 
   write(msg,'(a,2(a,i2,a))')ch10,&
&   '  time-reversal symmetry    +     zero Umklapp vector : ',nop(2)-nopg0(2),ch10,&
&   '  time-reversal symmetry    + non-zero Umklapp vector : ',nopg0(2),ch10
   call wrtout(unt,msg,mode)
  end if
 end do

end subroutine print_little_group 
!!***
