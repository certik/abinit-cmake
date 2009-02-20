!{\src2tex{textfont=tt}}
!!****f* ABINIT/bands_classification
!! NAME
!! bands_classification
!!
!! FUNCTION
!!  Initialize a Bands_symmetries datatypes containing information
!!  needed to analyze the irreducible representations at a particular k-point
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nspinor=number of spinorial components
!!  nsppol=number of independent polarizations
!!  nbnds=number of bands (supposed to be the same for spin up and down)
!!  kpt(3)=the k-point where the classification of bands is required
!!  nsym=number of space group symmetries
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space (reduced coordinates)
!!  tnons(3,nsym)=fractional translations of the space group
!!  ene_k(nbnds,nsppol)=energies for this k-point
!!  EDIFF_TOL=tolerance below which two states are considered to belong to the same irreducible representation 
!!
!! OUTPUT
!!  ierr=if different from 0, the bands classification cannot be performed, 
!!   The present implementation does NOT work at zone border if the space group of
!!   the crystal is non-symmorphic (non-zero fractionary translations)
!!  Bsym<Bands_Symmetries>= Initialized data type gathering information of the small group
!!   of the k-point as well as the irreducible representations 
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

subroutine init_Bands_Symmetries(BSym,only_trace,nspinor,nsppol,nbnds,kpt,nsym,symrec,tnons,EDIFF_TOL,ene_k,ierr)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_13recipspace
 use interfaces_15gw, except_this_one => init_Bands_Symmetries
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbnds,nspinor,nsppol,nsym
 integer,intent(out) :: ierr
 real(dp),intent(in) :: EDIFF_TOL
 logical,intent(in) :: only_trace
 type(Bands_Symmetries),intent(inout) :: BSym
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 real(dp),intent(in) :: ene_k(nbnds,nsppol),kpt(3),tnons(3,nsym)

!Local variables-------------------------------
!scalars
 integer :: dim_cplx,dummy,ib,ic,iclass,idx,ie,ii,is,isym,itim,kk,nclass,ncplx
 integer :: ncplx_MAX,nelmt,nsym_sgk,timrev_
 logical :: has_inversion,is_symmorphic
 character(len=500) :: msg
!arrays
 integer :: inversion(3,3),ksym(4,2,nsym)
 integer,allocatable :: G0(:,:),bbounds(:,:,:),elements_idx(:,:),nelements(:)
 integer,allocatable :: sgk2symrec(:),smallgk(:,:,:)
 real(dp),allocatable :: tnons_k(:,:)

! *************************************************************************

!#if defined DEBUG_MODE
! write(msg,'(a)')' init_Bands_Symmetries : enter '
! call wrtout(std_out,msg,'COLL') 
! call flush_unit(std_out)
!#endif

 ierr=0 ; timrev_=1 !this has to be tested, for the moment do not include time-reversal
 inversion=RESHAPE((/-1,0,0,0,-1,0,0,0,-1/),(/3,3/)) ; has_inversion=.FALSE.
 !
 ! ==== Copy basic info ====
 BSym%kpt(:)=kpt(:)
 BSym%nspinor=nspinor
 BSym%nsppol=nsppol
 BSym%nbnds=nbnds
 BSym%timrev=timrev_
 BSym%only_trace=only_trace
 Bsym%is_symmorphic=.FALSE.
 BSym%tol_deg=EDIFF_TOL
 !
 ! === Find small group of kpt ===
 call symq3(nsym,kpt,ksym,symrec,dummy)

 nsym_sgk=0 
 do itim=1,timrev_ 
  nsym_sgk=nsym_sgk+SUM(ksym(4,itim,:))
 end do
 allocate(smallgk(3,3,nsym_sgk),sgk2symrec(nsym_sgk),G0(3,nsym_sgk))
 
 idx=0 ; is_symmorphic=.FALSE.
 do itim=1,timrev_
  do isym=1,nsym
   if (ksym(4,itim,isym)==1) then 
    idx=idx+1
    smallgk(:,:,idx)=symrec(:,:,isym)*(3-2*itim) 
    G0(:,idx)=ksym(1:3,itim,isym)
    sgk2symrec(idx)=isym
    if (ALL(symrec(:,:,isym)==inversion)) has_inversion=.TRUE.
    ! If this case symmetry operations with fractional translations might not form a group
    if (ANY(ksym(1:3,itim,isym)/=0).and.(ANY(ABS(tnons(:,isym))>tol6))) is_symmorphic=.TRUE.
   end if
  end do
 end do
 if (is_symmorphic) then 
  write(*,*)' Non-symmorphic small group and zone border '
  write(*,*)' Character analysis not available '
  BSym%is_symmorphic=is_symmorphic 
  ierr=1
  !RETURN
 end if
 BSym%has_inversion=has_inversion
 !
 ! === Find classes and elements ===
 allocate(nelements(nsym_sgk),elements_idx(nsym_sgk,nsym_sgk))
 call get_class(nsym_sgk,smallgk,nclass,nelements,elements_idx)
 !
 ! ==== Save info on small group ====
 BSym%nsym_sgk=nsym_sgk ; BSym%nclass=nclass
 allocate(BSym%nelements(nclass)) ; BSym%nelements(1:nclass)=nelements(1:nclass)
 allocate(BSym%sgk2symrec(nsym_sgk),BSym%G0(3,nsym_sgk))
 kk=0
 do ic=1,nclass 
  do ie=1,nelements(ic)
   kk=kk+1
   BSym%sgk2symrec(kk)=sgk2symrec(elements_idx(ie,ic))
   BSym%G0(:,kk)=G0(:,elements_idx(ie,ic))
  end do
 end do
 deallocate(smallgk,sgk2symrec,nelements,elements_idx,G0)
 !
 ! ==== Initialize degenerate_bands structure ====
 allocate(BSym%ncplx(nsppol),bbounds(2,nbnds,nsppol))
 do is=1,nsppol
  ! ==== Find degenerate_bands ====
  ncplx=1 ; bbounds(:,:,is)=0 ; bbounds(1,1,is)=1
  do ib=2,nbnds
   if (ABS(ene_k(ib,is)-ene_k(ib-1,is))>EDIFF_TOL) then
    bbounds(2,ncplx,is)=ib-1
    ncplx=ncplx+1
    bbounds(1,ncplx,is)=ib
   end if
  end do
  bbounds(2,ncplx,is)=nbnds
  BSym%ncplx(is)=ncplx 
 end do
 !
 ! ==== Initialize Degenerate_Bands structure ====
 ncplx_MAX=MAXVAL(Bsym%ncplx) 
 allocate(BSym%Cplx(ncplx_MAX,nsppol)) ; call nullify_Degenerate_Bands(BSym%Cplx)
 do is=1,nsppol
  do ii=1,Bsym%ncplx(is)
   dim_cplx=bbounds(2,ii,is)-bbounds(1,ii,is)+1
   BSym%Cplx(ii,is)%dim_cplx=dim_cplx
   BSym%Cplx(ii,is)%ib_start=bbounds(1,ii,is)
   BSym%Cplx(ii,is)%ib_end  =bbounds(2,ii,is)
   allocate(BSym%Cplx(ii,is)%trace(nclass))
   if (.not.only_trace) allocate(BSym%Cplx(ii,is)%Rirr(dim_cplx,dim_cplx,nsym_sgk))
   allocate(BSym%Cplx(ii,is)%ene(dim_cplx))
   do ib=0,dim_cplx-1 
    kk=ib+bbounds(1,ii,is)
    BSym%Cplx(ii,is)%ene(ib+1)=ene_k(kk,is)
   end do
  end do
 end do !isppol
 deallocate(bbounds)

!#if defined DEBUG_MODE
! write(msg,'(a)')' init_Bands_Symmetries : exit '
! call wrtout(std_out,msg,'COLL') 
! call flush_unit(std_out)
!#endif
 
end subroutine init_Bands_Symmetries
!!***


!!****if* ABINIT/print_Bands_Symmetries
!! NAME
!! print_Bands_Symmetries
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

subroutine print_Bands_Symmetries(Bsym,unitno,prtvol,mode_paral)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalaras
!scalars
 integer,intent(in),optional :: prtvol,unitno
 character(len=4),intent(in),optional :: mode_paral
 type(Bands_Symmetries),intent(in) :: BSym

!Local variables-------------------------------
!scalars
 integer :: icl,isp,isym,ix,ix1,ix2,nclass,ncplx,nsp,nsym_sgk,unt,verbose
 complex(dpc) :: test
 character(len=4) :: mode
 character(len=500) :: fmt,msg
 type(Degenerate_Bands),pointer :: Cplx
!arrays
 complex(dpc),pointer :: Rirr(:,:),ntrace2(:),trace(:),trace1(:)

! *********************************************************************

 unt=std_out ; if (PRESENT(unitno)) unt=unitno
 verbose=0   ; if (PRESENT(prtvol)) verbose=prtvol
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 nsp=Bsym%nsppol ; nclass=Bsym%nclass ; nsym_sgk=BSym%nsym_sgk 
 write(fmt,*)'(2a,3es16.8,3a,i4,2a,',nsp,'i3,2a,i2,2a,i2,a,,',nclass,'i2,a)'
 write(msg,fmt)ch10,&
& ' ===== Character of bands at k-point: ',BSym%kpt(:),' ===== ',ch10,&
& '  Total number of bands ',BSym%nbnds,ch10,&
& '  Number of set of degenerate states detected ',Bsym%ncplx(1:nsp),ch10,&
& '  Number of operations in the little group ',nsym_sgk,ch10,&
& '  Number of classes ',nclass,' (',(BSym%nelements(icl),icl=1,nclass),' )' 
 call wrtout(unt,msg,mode)

 if (Bsym%is_symmorphic) then 
  write(msg,'(4a)')ch10,&
&  ' Non-symmorphic small group and zone border ',ch10,&
&  ' Character analysis not available '
  call wrtout(unt,msg,mode)
  !RETURN
 end if

 write(fmt,*)'(i3,a,i3,2x,',nclass,'(a,2f4.1,1x),a)'
 do isp=1,Bsym%nsppol
  ncplx=Bsym%ncplx(isp)
  do ix=1,ncplx
   Cplx => BSym%Cplx(ix,isp)
   write(msg,fmt)&
&   Cplx%ib_start,'-',Cplx%ib_end,('|',Cplx%trace(icl),icl=1,nclass),'|'
   call wrtout(unt,msg,mode)
  end do
 end do
 !
 ! === Test basic properties of irreducible representations ===
 !
 ! 1) \sum_R chi^*_a(R)\chi_b(R)= N_R \delta_{ab} 
 allocate(ntrace2(nclass))
 do isp=1,Bsym%nsppol
  ncplx=Bsym%ncplx(isp)
  do ix2=1,ncplx
   ntrace2(:)=BSym%nelements(:)*BSym%Cplx(ix2,isp)%trace(:)
   do ix1=1,ix2
    trace1 => BSym%Cplx(ix1,isp)%trace(:)
    test=DOT_PRODUCT(trace1,ntrace2) ; if (ix1==ix2) test=test-cone
    if (ABS(test)>tol6) write(unt,*)'cx1 cx2',ix1,ix2,test/Bsym%nsym_sgk
   end do
  end do
 end do !isp
 deallocate(ntrace2)

 !test if they are unitary
 !do isp=1,Bsym%nsppol
 ! ncplx=Bsym%ncplx(isp)
 !  do ix1=1,ncplx
 !   do isym=1,nsym_sgk
 !   Rirr => BSym%Cplx(ix1,isp)%Rirr(:,:,isym)
 !   write(unt,*)"===== ",ix1,isym," ====="
 !   write(unt,*)MATMUL(Rirr,TRANSPOSE(CONJG(Rirr)))
 !  end do
 ! end do
 !end do !isp

end subroutine print_Bands_Symmetries


subroutine destroy_Bands_Symmetries(Bsym)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15gw, except_this_one => destroy_Bands_Symmetries
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Bands_Symmetries),intent(inout) :: BSym
!Local variables-------------------------------
!scalars
 integer :: ii,isp
! *************************************************************************

 if (associated(BSym%G0)        ) deallocate(BSym%G0)
 if (associated(Bsym%nelements) ) deallocate(Bsym%nelements)  
 if (associated(Bsym%sgk2symrec)) deallocate(Bsym%sgk2symrec) 
 if (associated(Bsym%ncplx)     ) deallocate(Bsym%ncplx)  

 call destroy_Degenerate_Bands(Bsym%Cplx)
 if (associated(BSym%Cplx )) deallocate(BSym%Cplx)

end subroutine destroy_Bands_Symmetries


subroutine destroy_Degenerate_Bands(Cplx)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(Degenerate_Bands),intent(inout) :: Cplx(:,:)
!Local variables-------------------------------
!scalars
 integer :: ix,isp,nsppol,ncplx_MAX
! *************************************************************************

 nsppol=SIZE(Cplx,DIM=2) ; ncplx_MAX=SIZE(Cplx,DIM=1)
 do isp=1,nsppol
  do ix=1,ncplx_MAX
   if (associated(Cplx(ix,isp)%ene)  ) deallocate(Cplx(ix,isp)%ene)
   if (associated(Cplx(ix,isp)%Rirr) ) deallocate(Cplx(ix,isp)%Rirr)
   if (associated(Cplx(ix,isp)%trace)) deallocate(Cplx(ix,isp)%trace)
  end do
 end do

end subroutine destroy_Degenerate_Bands


subroutine nullify_Bands_Symmetries(BSym)

 use defs_basis
 use defs_datatypes 

 implicit none

!Arguments ------------------------------------
!scalars
 type(Bands_symmetries),intent(inout) :: Bsym
! *************************************************************************

 nullify(BSym%G0)
 nullify(BSym%ncplx)
 nullify(BSym%nelements)
 nullify(BSym%sgk2symrec)
 nullify(BSym%Cplx)

end subroutine nullify_Bands_Symmetries 


subroutine nullify_Degenerate_Bands(Cplx)

 use defs_basis
 use defs_datatypes 

 implicit none

!Arguments ------------------------------------
!scalars
 type(Degenerate_Bands),intent(inout) :: Cplx(:,:)
!Local variables-------------------------------
!scalars
 integer :: ix,isp,nsppol,ncplx_MAX
! *************************************************************************

 nsppol=SIZE(Cplx,DIM=2) ; ncplx_MAX=SIZE(Cplx,DIM=1)
 do isp=1,nsppol
  do ix=1,ncplx_MAX
   nullify(Cplx(ix,isp)%ene)
   nullify(Cplx(ix,isp)%Rirr)
   nullify(Cplx(ix,isp)%trace)
  end do
 end do

end subroutine nullify_Degenerate_Bands


!!***f* ABINIT/get_class
!! NAME
!! get_class
!!
!! FUNCTION
!!  Given a set of nsym 3x3 operations in reciprocal/real space,
!!  which are supposed to form a group, divide the group in classes
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! nsym=number of symmetry operation
!! symrec(3,3,nsym)=the operations
!!
!! OUTPUT
!! nclass=The number of classes
!! nelements(1:nclass)=For each class, the number of elements 
!! elements_idx(ii,1:nclass)=For each class, it reports the elements (ii=1,..,nelements(jclass))
!!
!! SIDE EFFECTS
!!
!! NOTES
!! symafm not yet treated
!! Does not work in case of non-collinear magnetism
!! No test is done to check if the input set forms a group 
!!
!! A class is defined as the set of distinct elements obtained by 
!! considering for each element, S, of the group all its conjugate
!! elements X^-1 S X where X range over all the elements of the group.
!!
!! PARENTS
!!  
!!
!! CHILDREN
!! 
!!
!! SOURCE

subroutine get_class(nsym,symrec,nclass,nelements,elements_idx)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(out) :: nclass
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 integer,intent(out) :: nelements(nsym),elements_idx(nsym,nsym)
!Local variables-------------------------------
!scalars
 integer :: isym,jsym,ksym,ie
 character(len=500) :: msg
!arrays
 logical :: found(nsym),found_identity
 integer :: cjg(3,3),ss(3,3),xx(3,3),xxm1(3,3),test(3,3)
 integer :: identity(3,3)

!************************************************************************

 identity(:,:)=RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/)) ; found_identity=.FALSE. 
 do isym=1,nsym
  if (ALL(symrec(:,:,isym)==identity)) then 
   found_identity=.TRUE. ; ie=isym ; EXIT
  end if
 end do 
 if (.not.found_identity.or.ie/=1) then 
  write(msg,'(6a)')ch10,&
&  ' get_class : ERROR -',ch10,&
&  ' either identity is not present or it is not the first operation ',ch10,&
&  ' check set of symmetry operations '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if 

 nclass=0 ; nelements(:)=0 ; elements_idx(:,:)=0 ; found(:)=.FALSE.
 do isym=1,nsym
  if (.not.found(isym)) then 
   nclass=nclass+1
   ss(:,:)=symrec(:,:,isym)
   ! === Form conjugate ===
   do jsym=1,nsym
    xx(:,:)=symrec(:,:,jsym)
    call mati3inv(xx,xxm1) ; xxm1=TRANSPOSE(xxm1)
    cjg(:,:)=MATMUL(xxm1,MATMUL(ss,xx))
    ! === Is it already found? ===
    do ksym=1,nsym
     test(:,:)=symrec(:,:,ksym)
     if (.not.found(ksym).and.(ALL((test-cjg)==0))) then 
      found(ksym)=.TRUE.
      nelements(nclass)=nelements(nclass)+1
      elements_idx(nelements(nclass),nclass)=ksym
     end if
    end do
   end do
  end if
 end do

end subroutine get_class
!!***
