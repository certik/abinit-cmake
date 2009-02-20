!{\src2tex{textfont=tt}}
!!****f* ABINIT/check_zarot
!! NAME
!! check_zarot
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (the_author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
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

subroutine check_zarot(nsym,symrec,timrev,npwvec,rprimd,gprimd,gvec,psps,pawang,grottb,grottbm1)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwvec,nsym,timrev
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: grottb(npwvec,timrev,nsym),grottbm1(npwvec,timrev,nsym)
 integer,intent(in) :: gvec(3,npwvec),symrec(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: aa,ig,ig_sym,iginv_sym,ii,ilpa,ilpm,isym,itim,jj,ll,lmax,mm,mqmem_
 integer :: nqpt_,optder,sign_l,two_lmaxp1
 real(dp) :: err,max_diff,test,tmp,ylm_sym
 logical :: found
 character(len=500) :: message
 type(MPI_type) :: mpi_enreg_seq
!arrays
 integer :: toinv(nsym),trial(3,3)
 integer,allocatable :: nband(:),npwarr(:)
 real(dp),allocatable :: DS_mmpl(:,:,:),DSinv_mmpl(:,:,:),qptns(:,:),ylm_q(:,:)
 real(dp),allocatable :: ylmgr_q(:,:,:)

! *************************************************************************

 write(message,'(a)')' check_zarot  : enter '
 call wrtout(std_out,message,'COLL') !;call flush_unit(std_out)

 do jj=1,nsym 
  found=.FALSE.
  do ii=1,nsym
   call mati3inv(symrec(:,:,ii),trial) 
   trial=transpose(trial)
   if (ALL(trial==symrec(:,:,jj))) then
    toinv(jj)=ii
    found=.TRUE.
    exit
   end if
  end do
  if (.not. found) stop "inverse not found! "
 end do


 mqmem_=1 ; nqpt_=1 ; optder=0
 allocate(npwarr(mqmem_),qptns(3,mqmem_)) 
 npwarr(:)=npwvec ; qptns(:,:)=zero

 lmax=psps%mpsang-1
 write(*,*)'lmax= ',lmax
 allocate(ylm_q(npwvec*mqmem_,(lmax+1)**2))
 allocate(ylmgr_q(npwvec*mqmem_,3+6*(optder/2),(lmax+1)**2))
 call initmpi_seq(mpi_enreg_seq) ; mpi_enreg_seq%nproc_fft=1 ; mpi_enreg_seq%me_fft=0
 allocate(nband(1)) ; nband=0

 ! Note: dtset%nband and dtset%nsppol are not used in sequential mode
 call initylmg(gprimd,gvec,qptns,mqmem_,mpi_enreg_seq,psps%mpsang,npwvec,nband,nqpt_,npwarr,0,optder,rprimd,1,1,ylm_q,ylmgr_q)

 allocate(DS_mmpl(2*lmax+1,2*lmax+1,lmax+1))
 allocate(DSinv_mmpl(2*lmax+1,2*lmax+1,lmax+1))
 max_diff=zero ; test=zero

 do ig=1,npwvec
  if (ig==1) cycle

  do isym=1,nsym
   do itim=1,timrev

    ig_sym=grottb(ig,itim,isym) !index of IS G
    DS_mmpl(:,:,:)=pawang%zarot(:,:,:,isym)

    iginv_sym=grottbm1(ig,itim,isym) !index of (IS)^-1 G
    DSinv_mmpl(:,:,:)=pawang%zarot(:,:,:,toinv(isym))

    do ll=0,lmax
     do mm=1,2*ll+1
      ilpm=1+ll**2+ll+(mm-1-ll)
      ylm_sym=ylm_q(ig_sym,ilpm) !Ylm(ISG)
      !
      ! here we calculate the symmetric
      tmp=zero
      do aa=1,2*ll+1
       test=MAX(test,ABS(DS_mmpl(aa,mm,ll+1)-DSinv_mmpl(mm,aa,ll+1)))
       ilpa=1+ll**2+ll+(aa-1-ll)
       tmp= tmp+ ylm_q(ig,ilpa)*DS_mmpl(aa,mm,ll+1)
      end do
      if (itim==2) tmp=tmp*(-1)**ll
      err=ABS(tmp-ylm_sym)

      if (err > tol6) then
       write(77,'(6(a,i3),a)')' -- ig: ',ig,' igsym: ',ig_sym,' isym ',isym,' itim:',itim,' ll: ',ll,' mm: ',(mm-1-ll)," --" 
       write(77,*)tmp,ylm_sym,ABS(tmp-ylm_sym)
      end if
      max_diff=max(max_diff,err)

     end do
    end do !itim

   end do  !isym
  end do !sym
 end do !ig

 write(*,*)"MAX DIFF ",max_diff
 write(*,*)"MAX TEST ",test

 deallocate(DS_mmpl,DSinv_mmpl,nband)
 deallocate(npwarr,qptns) 
 deallocate(ylm_q,ylmgr_q)

end subroutine check_zarot
!!***


!!****if* ABINIT/get_TDij0
!! NAME
!! get_TDij0
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

subroutine get_TDij0(npw,natom,ntypat,typat,kpt,gmet,gvec,pawtab,cprj1,cprj2,cg1,cg2,res)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_13recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,npw,ntypat
!arrays
 integer,intent(in) :: gvec(3,npw),typat(natom)
 real(dp),intent(in) :: gmet(3,3),kpt(3)
 real(dp),intent(out) :: res(2)
 complex(gwpc),intent(in) :: cg1(npw),cg2(npw)
 type(cprj_type),intent(in) :: cprj1(natom),cprj2(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iat,ij_size,ilmn,ipw,itypat,jlmn,k0lmn,klmn,l_size,lm_size
 integer :: lmn2_size,lmn_size,mesh_size
 real(dp) :: dij0,fij,im_p,re_p
 complex(gwpc) :: ctmp
!arrays
 integer,allocatable :: indklmn_(:,:)
 real(dp) :: kinpw(npw)
 complex(gwpc) :: vec(npw)

! *************************************************************************

 call mkkin(10000.0_dp,zero,one,gmet,gvec,kinpw,kpt,npw)

 vec(:)=cg2(:)*kinpw(:)
 ctmp=DOT_PRODUCT(cg1,vec)

 res(1)=REAL(ctmp)
 res(2)=AIMAG(ctmp)

 do iat=1,natom
  itypat=typat(iat)
  l_size=pawtab(itypat)%l_size
  lmn_size=pawtab(itypat)%lmn_size
  lmn2_size=pawtab(itypat)%lmn2_size
  ij_size=pawtab(itypat)%ij_size
  allocate(indklmn_(4,lmn2_size)) 
  indklmn_(:,:)=pawtab(itypat)%indklmn(:,:)

  do jlmn=1,lmn_size
   k0lmn=jlmn*(jlmn-1)/2 
   do ilmn=1,jlmn

     re_p=   cprj1(iat)%cp(1,ilmn)*cprj2(iat)%cp(1,jlmn) &
&           +cprj1(iat)%cp(2,ilmn)*cprj2(iat)%cp(2,jlmn) &
&           +cprj1(iat)%cp(1,jlmn)*cprj2(iat)%cp(1,ilmn) &
&           +cprj1(iat)%cp(2,jlmn)*cprj2(iat)%cp(2,ilmn) 

     im_p=   cprj1(iat)%cp(1,ilmn)*cprj2(iat)%cp(2,jlmn) &
&           -cprj1(iat)%cp(2,ilmn)*cprj2(iat)%cp(1,jlmn) &
&           +cprj1(iat)%cp(1,jlmn)*cprj2(iat)%cp(2,ilmn) &
&           -cprj1(iat)%cp(2,jlmn)*cprj2(iat)%cp(1,ilmn)


     klmn=k0lmn+ilmn ; fij=one ; if (jlmn==ilmn) fij=half
     dij0=pawtab(itypat)%dij0(klmn)
     res(1)=res(1)+ fij*re_p*dij0 
     res(2)=res(2)+ fij*im_p*dij0 

   end do !ilmn
  end do !jlmn 
  deallocate(indklmn_)
 end do !iat

end subroutine get_TDij0


subroutine test_PAWH(natom,ntypat,typat,nspinor,Kmesh,gmet,gvec,nfftf,vpsps,usepaw,pawtab,cprj_ibz,&
& b1gw,b2gw,sp,en,vxc,vhartr,Wfs,ngfftf,igfftf,mpi_enreg)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15gw, except_this_one => test_PAWH
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(sigma_parameters),intent(in) :: sp
 type(wavefunctions_information),intent(in) :: Wfs
 type(BZ_mesh_type),intent(in) :: Kmesh
 integer,intent(in) :: igfftf(Wfs%npwwfn),ngfftf(18)
 integer,intent(in) :: b1gw,b2gw,natom,ntypat,nspinor,usepaw,nfftf
 integer,intent(in) :: typat(natom)
 integer,intent(in) :: gvec(3,Wfs%npwwfn)
 real(dp),intent(in) :: gmet(3,3)
 type(cprj_type),intent(in) :: cprj_ibz(natom,nspinor*Kmesh%nibz*sp%nbnds*sp%nsppol)
 type(pawtab_type),intent(in) :: pawtab(ntypat)
 type(MPI_type),intent(inout) :: mpi_enreg
 real(dp),intent(in) :: en(Kmesh%nibz,sp%nbnds,sp%nsppol)
 real(dp),intent(in) :: vpsps(nfftf)
 complex(gwpc),intent(in) :: vxc(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,sp%nsppol)
 complex(gwpC),intent(in) :: vhartr(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,sp%nsppol)

!Local variables-------------------------------
 integer :: is,ik,ib1,ib2,iat
 type(cprj_type),allocatable :: cprj1(:),cprj2(:)
 integer :: dimlmn(natom)
 real(dp) :: test(2),dij0(2)
 complex(gwpc),pointer :: cg1(:),cg2(:)
 complex(gwpc) :: ctmp
 complex(gwpc),allocatable :: wfr1(:),wfr2(:)
! *************************************************************************

 do iat=1,natom
  dimlmn(iat)=pawtab(typat(iat))%lmn_size
 end do
 allocate(cprj1(natom),cprj2(natom))
 ! TODO use new cprj methods!
 !call cprj_alloc(cprj1,0,dimlmn)
 !call cprj_alloc(cprj2,0,dimlmn)

 allocate(wfr1(nfftf),wfr2(nfftf))

 do ik=1,Kmesh%nibz
  write(*,*)' ik   band1 band2 ispin PAW Hamiltonian test '
  do is=1,sp%nsppol  
   do ib1=b1gw,b2gw
    do ib2=ib1,b2gw

     !call cprj1_get(cprj_ibz,ik,ib1,is,1,sp%nbnds,Kmesh%nibz,sp%nsppol,cprj1)
     !call cprj1_get(cprj_ibz,ik,ib2,is,1,sp%nbnds,Kmesh%nibz,sp%nsppol,cprj2)

     cg1 => Wfs%wfg(:,ib1,ik,is)
     cg2 => Wfs%wfg(:,ib2,ik,is)
     call get_TDij0(Wfs%npwwfn,natom,ntypat,typat,Kmesh%ibz(:,ik),gmet,gvec,pawtab,cprj1,cprj2,cg1,cg2,dij0)

     call fft_onewfn(Wfs%paral_kgb,Wfs%nspinor,Wfs%npwwfn,nfftf,cg1,wfr1,igfftf,ngfftf,0,mpi_enreg)
     call fft_onewfn(Wfs%paral_kgb,Wfs%nspinor,Wfs%npwwfn,nfftf,cg2,wfr2,igfftf,ngfftf,0,mpi_enreg)

     ctmp=SUM(CONJG(wfr1(:))*vpsps(:)*wfr2(:))/nfftf
     test(:)=zero
     if (ib1==ib2) then 
      test(1)=en(ik,ib1,is)-(dij0(1)+vhartr(ib1,ib1,ik,is)+vxc(ib1,ib1,ik,is))
      test(1)=test(1)-REAL(ctmp)
     else 
      test(1)=(dij0(1)+REAL(vhartr(ib1,ib2,ik,is)+vxc(ib1,ib2,ik,is)))
      test(2)=(dij0(2)+AIMAG(vhartr(ib1,ib2,ik,is)+vxc(ib1,ib2,ik,is)))
      test(1)=test(1)+REAL(ctmp)
      test(2)=test(2)+AIMAG(ctmp)
     end if


     write(*,'(3i4,2(2x,f10.6))')ik,ib1,is,test(1)*Ha_eV,test(2)*Ha_eV

    end do
   end do
  end do
 end do

 deallocate(wfr1,wfr2)
 !call cprj_free(cprj1) ; deallocate(cprj1)
 !call cprj_free(cprj2) ; deallocate(cprj2)

end subroutine test_PAWH


subroutine hack_symmetries(str,nsym,timrev)

 use defs_basis
 use m_IO_tools, only : get_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(inout) :: nsym,timrev
 character(len=*),intent(in) :: str

!Local variables-------------------------------
 integer :: unt
 logical :: ltemp
 character(len=500) :: msg
! *********************************************************************

 inquire(file='__sym.in__',exist=ltemp)
 if (ltemp) then
  write(msg,'(4a)')ch10,TRIM(str)//': COMMENT - ',ch10,&
&  ' reading symmetry operations from file __sym.in__'
  call wrtout(std_out,msg,'COLL') 
  unt=get_unit() ; open(unit=unt,file='__sym.in__')
  read(unt,*)nsym,timrev ; close(unt)
  write(msg,'(2a,i4,a,i2)')ch10,&
&  ' Code will use first ',nsym,' symmetry operations and time_reversal ',timrev
  call wrtout(std_out,msg,'COLL') 
 end if

end subroutine hack_symmetries 
!!***
