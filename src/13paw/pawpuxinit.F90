!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawpuxinit
!! NAME
!! pawpuxinit
!!
!! FUNCTION
!! Initialize some starting values of several arrays used in
!! PAW+U or local exact-exchange calculations
!!
!!
!! A-define useful indices for LDA+U/local exact-exchange
!! B-Compute overlap between atomic wavefunction
!! C-Compute matrix elements of coulomb interaction (see PRB vol.52 5467)
!!    (angular part computed from Gaunt coefficients)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (BA,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors.
!!
!! INPUTS
!!  dmatpuopt= select expression for the density matrix
!!  exchmix= mixing factor for local exact-exchange
!!  jpawu(ntypat)= value of J
!!  llexexch(ntypat)= value of l on which local exact-exchange applies
!!  llpawu(ntypat)= value of l on which LDA+U applies
!!  indlmn(6,i,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn (for a given atom type)
!!  lmnmax=max. number of (l,m,n) components over all type of psps
!!  ntypat=number of types of atoms in unit cell.
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!     %lmax=Maximum value of angular momentum l+1
!!     %gntselect((2*l_max-1)**2,l_max**2,l_max**2)=
!!                     selection rules for Gaunt coefficients
!!  pawprtvol=output printing level for PAW
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data:
!!  upawu(ntypat)= value of U
!!  useexexch= 0 if no local exact-exchange; 1 if local exact-exchange
!!  usepawu= 0 if no LDA+U; 1 if LDA+U
!!
!! OUTPUT
!!  pawtab <type(pawtab_type)>=paw tabulated data read at start:
!!     %ij_proj=nproj*(nproju+1)/2
!!     %klmntomn(4,lmn2_size) = Array giving im, jm ,in, and jn for each klmn=(ilmn,jlmn)
!!     %lnproju(nproj)= value of ln for projectors on which paw+u/local exact-exchange acts.
!!     %nproju=number of projectors for orbitals on which paw+u/local exact-exchange acts.
!!     %phiphjint(pawtabitypat%ij_proj)=Integral of Phi(:,i)*Phi(:,j) for correlated orbitals.
!!     %usepawu=0 if no LDA+U; 1 if LDA+U
!!     %useexexch=0 if no local exact-exchange; 1 if local exact-exchange
!!     === if usepawu>0
!!     %jpawu= value of J
!!     %upawu= value of U
!!     %vee(2*lpawu+1,:,:,:)=matrix of the screened interaction for correlated orbitals
!!     === if useexexch>0
!!     %fk
!!     %vex(2*lpawu+1,:,:,:)=matrix of the screened interaction for correlated orbitals
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      leave_new,simp_gen,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine pawpuxinit(dmatpuopt,exchmix,jpawu,llexexch,llpawu,indlmn,lmnmax,ntypat,pawang,&
&                      pawprtvol,pawrad,pawtab,upawu,useexexch,usepawu)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: dmatpuopt,lmnmax,ntypat,pawprtvol,useexexch,usepawu
 real(dp),intent(in) :: exchmix
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),llexexch(ntypat),llpawu(ntypat)
 real(dp),intent(in) :: jpawu(ntypat),upawu(ntypat)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(inout) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: icount,il,ilmn,isela,iselb,itemp,itypat,iu,j0lmn,jl,jlmn,ju,klm0u
 integer :: klm0x,klma,klmb,klmn,kln,kln1,kln2,kyc,lcur,lexexch,lkyc,ll,ll1
 integer :: lmexexch,lmkyc,lmn_size,lmpawu,lpawu,m1,m11,m2,m21,m3,m31,m4,m41
 integer :: mesh_size,int_meshsz,mkyc
 real(dp) :: ak,f4of2,f6of2,int1,int2,intg
 character(len=500) :: message
!arrays
 integer :: indlmn_(6,lmnmax)
 real(dp),allocatable :: ff(:),fk(:),gg(:)

! *************************************************************************

 if(useexexch==0.and.usepawu==0) return

!PAW+U and local exact-exchange restriction
 if(useexexch>0.and.usepawu>0)then
  do itypat=1,ntypat
   if (llpawu(itypat)/=llexexch(itypat).and.&
&   llpawu(itypat)/=-1.and.llexexch(itypat)/=-1) then
    write(message, '(8a,i2,3a)' ) ch10,&
&    ' pawpuxinit: ERROR - ',ch10,&
&    '  When PAW+U (usepawu>0) and local exact-exchange (exexch>0)',ch10,&
&    '  are selected together, they must apply on the same',ch10,&
&    '  angular momentum (lpawu/=lexexch forbidden, here for typat=',&
&    itypat,') !',ch10,'  Action: correct your input file.'
    call wrtout(ab_out,message,'COLL')
    call wrtout(06,  message,'COLL')
    call leave_new('COLL')
   end if
  end do
 end if

!Print title
 write(message, '(3a)' ) ch10,ch10,"******************************************"
 if(usepawu==1) then
  write(message, '(3a)' ) trim(message),ch10," LDA+U Method used: FLL"
 else if(usepawu==2) then
  write(message, '(3a)' ) trim(message),ch10," LDA+U Method used: AMF"
 end if
 if(useexexch>0) write(message, '(3a)' ) trim(message),ch10," PAW Local Exact exchange: PBE0"
 write(message, '(3a)' ) trim(message),ch10,"******************************************"
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')

!Loop on atom types
 do itypat=1,ntypat
  indlmn_(:,:)=indlmn(:,:,itypat)
  lmn_size=pawtab(itypat)%lmn_size
  mesh_size=pawrad(itypat)%mesh_size
  int_meshsz=pawrad(itypat)%int_meshsz

! PAW+U data
  if (usepawu>0) then
   pawtab(itypat)%lpawu=llpawu(itypat)
   lcur=llpawu(itypat)
   if(pawtab(itypat)%lpawu==-1) pawtab(itypat)%usepawu=0
   if(pawtab(itypat)%lpawu/=-1) pawtab(itypat)%usepawu=usepawu
   if(lcur/=-1) then
    pawtab(itypat)%upawu=upawu(itypat)
    pawtab(itypat)%jpawu=jpawu(itypat)
   end if
  end if

! Local exact-echange data
  if (useexexch>0) then
   pawtab(itypat)%exchmix=exchmix
   pawtab(itypat)%lexexch=llexexch(itypat)
   lcur=llexexch(itypat)
   if(pawtab(itypat)%lexexch==-1) pawtab(itypat)%useexexch=0
   if(pawtab(itypat)%lexexch/=-1) pawtab(itypat)%useexexch=useexexch
  end if

! Select only atoms with +U
  if(lcur/=-1) then

!  Compute number of projectors for LDA+U/local exact-exchange
   icount=count(indlmn(1,1:lmn_size,itypat)==lcur)
   pawtab(itypat)%nproju=icount/(2*lcur+1)
   if(useexexch>0.and.pawtab(itypat)%nproju>2)  then
    write(message, '(a,a,a,a,a,a)' ) ch10,&
&    ' pawpuxinit : ERROR -',ch10,&
&    '  Error on the number of projectors ',ch10,&
&    '  more than 2 projectors is not allowed for local exact-exchange'
    call wrtout(ab_out,message,'COLL')
    call wrtout(06,  message,'COLL')
    call leave_new('COLL')
   end if
   if(pawtab(itypat)%nproju*(2*lcur+1)/=icount)  then
    write(message, '(4a)' ) ch10,&
&    ' pawpuxinit : BUG -',ch10,&
&    '  Error on the number of projectors '
    call wrtout(06,message,'COLL')
    call leave_new('COLL')
   end if
   write(message, '(a,a,i4,a,a,i4)' ) ch10,&
&   ' pawpuxinit : for species ',itypat,ch10,&
&   '   number of projectors is',pawtab(itypat)%nproju
   call wrtout(6,message,'COLL')

   pawtab(itypat)%ij_proj=pawtab(itypat)%nproju*(pawtab(itypat)%nproju+1)/2

   allocate(pawtab(itypat)%lnproju(pawtab(itypat)%nproju))
   allocate(pawtab(itypat)%phiphjint(pawtab(itypat)%ij_proj))

!  ==================================================
!  A-define useful indexes
!  --------------------------------------------------
   icount=0
   do ilmn=1,lmn_size
    if(indlmn_(1,ilmn)==lcur) then
     icount=icount+1
     itemp=(icount-1)/(2*lcur+1)
     if (itemp*(2*lcur+1)==icount-1) then
      pawtab(itypat)%lnproju(itemp+1)=indlmn_(5,ilmn)
     end if
    end if
   end do

   do jlmn=1,lmn_size
    jl= indlmn_(1,jlmn)
    j0lmn=jlmn*(jlmn-1)/2
    do ilmn=1,jlmn
     il= indlmn_(1,ilmn)
     klmn=j0lmn+ilmn
     pawtab(itypat)%klmntomn(1,klmn)=indlmn_(2,ilmn)+il+1
     pawtab(itypat)%klmntomn(2,klmn)=indlmn_(2,jlmn)+jl+1
     pawtab(itypat)%klmntomn(3,klmn)=indlmn_(3,ilmn)
     pawtab(itypat)%klmntomn(4,klmn)=indlmn_(3,jlmn)
    end do
   end do

!  ==================================================
!  B-PAW+U: overlap between atomic wavefunctions
!  --------------------------------------------------
!  if (usepawu>0) then

   if(dmatpuopt==1) then
    write(message, '(4a)' ) ch10,&
&    ' pawpuxinit : dmatpuopt=1 ',ch10,&
&    '   PAW+U: dens. mat. constructed by projection on atomic wfn inside PAW augm. region(s)'
    call wrtout(6,message,'COLL')
    write(message, '(8a)' ) ch10,&
&    ' pawpuxinit: WARNING: Check that the first partial wave for lpawu:', ch10, &
&    '                      - Is an atomic eigenfunction  ',ch10, &
&    '                      - Is normalized ',ch10, &
&    '                      In other cases, choose dmatpuopt=2'
    call wrtout(6,message,'COLL')
   else if(dmatpuopt==2) then
    write(message, '(6a)' ) ch10,&
&    ' pawpuxinit : dmatpuopt=2 ',ch10,&
&    '   PAW+U: dens. mat. constructed by selecting contribution',ch10,&
&    '          for each angular momentum to the density (inside PAW augm. region(s))'
    call wrtout(6,message,'COLL')
   else if(dmatpuopt==3) then
    write(message, '(a,a,a,a,a,a)' ) ch10,&
&    ' pawpuxinit : dmatpuopt=3 ',ch10,&
&    '    PAW+U: dens. mat. constructed by projection on atomic wfn inside PAW augm. region(s)',ch10,&
&    '           and normalized inside PAW augm. region(s)'
    call wrtout(6,message,'COLL')
    write(message, '(6a)' ) ch10,&
&    ' pawpuxinit: WARNING: Check that the first partial wave for lpawu:', ch10, &
&    '                     is an atomic eigenfunction',ch10, &
&    '                     In the other case, choose dmatpuopt=2'
    call wrtout(6,message,'COLL')
   end if

   allocate(ff(mesh_size))
   icount=0
   do ju=1,pawtab(itypat)%nproju
    do iu=1,ju
     icount=icount+1
     if ((dmatpuopt==1).and.(useexexch==0)) then
      ff(1:mesh_size)=pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(1))&
&      *pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(iu))
      call simp_gen(int1,ff,pawrad(itypat))
      ff(1:mesh_size)=pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(1))&
&      *pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(ju))
      call simp_gen(int2,ff,pawrad(itypat))
      pawtab(itypat)%phiphjint(icount)=int1*int2
     else if((dmatpuopt==2).or.(useexexch>0)) then
      ff(1:mesh_size)=pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(iu))&
&      *pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(ju))
      call simp_gen(int1,ff,pawrad(itypat))
      pawtab(itypat)%phiphjint(icount)=int1
     else if((dmatpuopt==3).and.(useexexch==0)) then
      ff(1:mesh_size)=pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(1))&
&      *pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(1))
      call simp_gen(intg,ff,pawrad(itypat))
      ff(1:mesh_size)=pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(1))&
&      *pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(iu))
      call simp_gen(int1,ff,pawrad(itypat))
      ff(1:mesh_size)=pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(1))&
&      *pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(ju))
      call simp_gen(int2,ff,pawrad(itypat))
      pawtab(itypat)%phiphjint(icount)=int1*int2/intg
     else
      write(message, '(6a)' ) ch10,&
&      ' pawpuxinit : ERROR -',ch10,&
&      '  PAW+U: dmatpuopt has a wrong value !',ch10,&
&      '  Action : change value in input file'
      call wrtout(ab_out,message,'COLL')
      call wrtout(06,  message,'COLL')
      call leave_new('COLL')
     end if
    end do
   end do
   if(pawtab(itypat)%ij_proj/=icount)  then
    write(message, '(4a)' ) ch10,&
&    ' pawpuxinit : BUG -',ch10,&
&    '  Error in the loop for calculating phiphjint '
    call wrtout(06,message,'COLL')
    call leave_new('COLL')
   end if
   deallocate(ff)
   if(abs(pawprtvol)>=2) then
    do icount=1,pawtab(itypat)%ij_proj
     write(message, '(a,a,i2,f9.5,a)' ) ch10,&
&     '  PAW+U: icount, phiphjint(icount)=',icount,pawtab(itypat)%phiphjint(icount)
     call wrtout(06,message,'COLL')
    end do
   end if
!  end if

!  ======================================================================
!  C-PAW+U: Matrix elements of coulomb interaction (see PRB vol.52 5467)
!  1. angular part computed from Gaunt coefficients
!  --------------------------------------------------------------------
   if (usepawu>0) then
    lpawu=lcur

!   a. compute F(k)
!   ---------------------------------------------
    allocate(fk(lpawu+1))
    fk(1)=pawtab(itypat)%upawu
    if(lpawu==2) then
     f4of2=0.625_dp
     fk(2)=pawtab(itypat)%jpawu*14._dp/(One+f4of2)
     fk(3)=fk(2)*f4of2
    else if(lpawu==3) then
     f4of2=0.6681_dp
     f6of2=0.4943_dp
     fk(2)=pawtab(itypat)%jpawu*6435._dp/(286._dp+195._dp*f4of2+250._dp*f6of2)
     fk(3)=fk(2)*f4of2
     fk(4)=fk(2)*f6of2
    else
     write(message, '(a,a,i3,a,a)' ) ch10,&
&     ' pawpuxinit :  ERROR - lpawu=',lpawu,ch10,&
&     '   lpawu not equal to 2 or 3 is not allowed'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')
    end if

!   b. Compute ak and vee.
!   ---------------------------------------------
    pawtab(itypat)%vee=zero
    lmpawu=(lpawu-1)**2+2*(lpawu-1)+1  ! number of m value below correlated orbitals
    klm0u=lmpawu*(lmpawu+1)/2          ! value of klmn just below correlated orbitals
!   --------- 4 loops for interaction matrix
    do m1=-lpawu,lpawu
     m11=m1+lpawu+1
     do m2=-lpawu,m1
      m21=m2+lpawu+1
!     klma= number of pair before correlated orbitals +
!     number of pair for m1 lower than correlated orbitals
!     (m1+lpawu+1)*(lpawu-1) + number of pairs for correlated orbitals
!     before (m1,m2) + number of pair for m2 lower than current value
      klma=klm0u+m11*lmpawu+(m11-1)*m11/2+m21
      do m3=-lpawu,lpawu
       m31=m3+lpawu+1
       do m4=-lpawu,m3
        m41=m4+lpawu+1
        klmb=klm0u+m31*lmpawu+(m31-1)*m31/2+m41
!       --------- loop on k=1,2,3 (4 if f orbitals)
        do kyc=1,2*lpawu+1,2
         lkyc=kyc-1
         lmkyc=(lkyc+1)*(lkyc)+1
         ak=zero
         do mkyc=-lkyc,lkyc,1
          isela=pawang%gntselect(lmkyc+mkyc,klma)
          iselb=pawang%gntselect(lmkyc+mkyc,klmb)
          if (isela>0.and.iselb>0) ak=ak +pawang%realgnt(isela)*pawang%realgnt(iselb)
         end do
!        ----- end loop on k=1,2,3 (4 if f orbitals)
         ak=ak/(two*dfloat(lkyc)+one)
         pawtab(itypat)%vee(m11,m31,m21,m41)=ak*fk(lkyc/2+1)+pawtab(itypat)%vee(m11,m31,m21,m41)
        end do  !kyc
        pawtab(itypat)%vee(m11,m31,m21,m41)=pawtab(itypat)%vee(m11,m31,m21,m41)*four_pi
        pawtab(itypat)%vee(m21,m31,m11,m41)=pawtab(itypat)%vee(m11,m31,m21,m41)
        pawtab(itypat)%vee(m11,m41,m21,m31)=pawtab(itypat)%vee(m11,m31,m21,m41)
        pawtab(itypat)%vee(m21,m41,m11,m31)=pawtab(itypat)%vee(m11,m31,m21,m41)
       end do
      end do
     end do
    end do
    deallocate(fk)

   end if ! usepawu

!  ======================================================================
!  D-Local ex-exchange: Matrix elements of coulomb interaction and Fk
!  ----------------------------------------------------------------------
   if (useexexch>0) then
    lexexch=lcur

!   a. compute F(k)
!   ---------------------------------------------
    pawtab(itypat)%fk=zero
    allocate(ff(mesh_size),gg(mesh_size))
    kln=(pawtab(itypat)%lnproju(1)*( pawtab(itypat)%lnproju(1)+1)/2)
    do ll=1,lexexch+1
     ll1=2*ll-2
     if (int_meshsz<mesh_size) ff(int_meshsz+1:mesh_size)=zero
     ff(1:int_meshsz)=pawtab(itypat)%phiphj(1:int_meshsz,kln)
     call poisson(ff,ll1,intg,pawrad(itypat),gg)
     ff(1)=zero
     ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln)*gg(2:mesh_size))&
&     /pawrad(itypat)%rad(2:mesh_size)
     call simp_gen(intg,ff,pawrad(itypat))
     pawtab(itypat)%fk(1,ll)=intg*(two*ll1+one)
    end do
    if (pawtab(itypat)%nproju==2) then
     kln1=kln+pawtab(itypat)%lnproju(1)
     kln2=kln1+1
     do ll=1,lexexch+1
      ll1=2*ll-2
      if (int_meshsz<mesh_size) ff(int_meshsz+1:mesh_size)=zero
      ff(1:int_meshsz)=pawtab(itypat)%phiphj(1:int_meshsz,kln1)
      call poisson(ff,ll1,intg,pawrad(itypat),gg)
      ff(1)=zero
      ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln1)*gg(2:mesh_size))&
&      /pawrad(itypat)%rad(2:mesh_size)
      call simp_gen(intg,ff,pawrad(itypat))
      pawtab(itypat)%fk(2,ll)=intg*(two*ll1+one)
     end do
     do ll=1,lexexch+1
      ll1=2*ll-2
      if (int_meshsz<mesh_size) ff(int_meshsz+1:mesh_size)=zero
      ff(1:int_meshsz)=pawtab(itypat)%phiphj(1:int_meshsz,kln2)
      call poisson(ff,ll1,intg,pawrad(itypat),gg)
      ff(1)=zero
      ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln2)*gg(2:mesh_size))&
&      /pawrad(itypat)%rad(2:mesh_size)
      call simp_gen(intg,ff,pawrad(itypat))
      pawtab(itypat)%fk(3,ll)=intg*(two*ll1+one)
     end do
     do ll=1,lexexch+1
      ll1=2*ll-2
      if (int_meshsz<mesh_size) ff(int_meshsz+1:mesh_size)=zero
      ff(1:int_meshsz)=pawtab(itypat)%phiphj(1:int_meshsz,kln)
      call poisson(ff,ll1,intg,pawrad(itypat),gg)
      ff(1)=zero
      ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln1)*gg(2:mesh_size))&
&      /pawrad(itypat)%rad(2:mesh_size)
      call simp_gen(intg,ff,pawrad(itypat))
      pawtab(itypat)%fk(4,ll)=intg*(two*ll1+one)
     end do
     do ll=1,lexexch+1
      ll1=2*ll-2
      if (int_meshsz<mesh_size) ff(int_meshsz+1:mesh_size)=zero
      ff(1:int_meshsz)=pawtab(itypat)%phiphj(1:int_meshsz,kln)
      call poisson(ff,ll1,intg,pawrad(itypat),gg)
      ff(1)=zero
      ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln2)*gg(2:mesh_size))&
&      /pawrad(itypat)%rad(2:mesh_size)
      call simp_gen(intg,ff,pawrad(itypat))
      pawtab(itypat)%fk(5,ll)=intg*(two*ll1+one)
     end do
     do ll=1,lexexch+1
      ll1=2*ll-2
      if (int_meshsz<mesh_size) ff(int_meshsz+1:mesh_size)=zero
      ff(1:int_meshsz)=pawtab(itypat)%phiphj(1:int_meshsz,kln1)
      call poisson(ff,ll1,intg,pawrad(itypat),gg)
      ff(1)=zero
      ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln2)*gg(2:mesh_size))&
&      /pawrad(itypat)%rad(2:mesh_size)
      call simp_gen(intg,ff,pawrad(itypat))
      pawtab(itypat)%fk(6,ll)=intg*(two*ll1+one)
     end do
     f4of2=0.6681_dp
     f6of2=0.4943_dp
    end if
    deallocate(ff,gg)

!   b. Compute vex.
!   ---------------------------------------------
    pawtab(itypat)%vex=zero
    lmexexch=(lexexch-1)**2+2*(lexexch-1)+1  ! number of m value below correlated orbitals
    klm0x=lmexexch*(lmexexch+1)/2            ! value of klmn just below correlated orbitals
!   --------- 4 loops for interaction matrix
    do m1=-lexexch,lexexch
     m11=m1+lexexch+1
     do m2=-lexexch,m1
      m21=m2+lexexch+1
!     klma= number of pair before correlated orbitals +
!     number of pair for m1 lower than correlated orbitals
!     (m1+lexexch+1)*(lexexch-1) + number of pairs for correlated orbitals
!     before (m1,m2) + number of pair for m2 lower than current value
      klma=klm0x+m11*lmexexch+(m11-1)*m11/2+m21
      do m3=-lexexch,lexexch
       m31=m3+lexexch+1
       do m4=-lexexch,m3
        m41=m4+lexexch+1
        klmb=klm0x+m31*lmexexch+(m31-1)*m31/2+m41
!       --------- loop on k=1,2,3 (4 if f orbitals)
        do kyc=1,2*lexexch+1,2
         lkyc=kyc-1
         ll=(kyc+1)/2
         lmkyc=(lkyc+1)*(lkyc)+1
         ak=zero
         do mkyc=-lkyc,lkyc,1
          isela=pawang%gntselect(lmkyc+mkyc,klma)
          iselb=pawang%gntselect(lmkyc+mkyc,klmb)
          if (isela>0.and.iselb>0) ak=ak +pawang%realgnt(isela)*pawang%realgnt(iselb)
         end do
!        ----- end loop on k=1,2,3 (4 if f orbitals)
         pawtab(itypat)%vex(m11,m31,m21,m41,ll)=ak/(two*dfloat(lkyc)+one)
        end do  !kyc
        do ll=1,4
         pawtab(itypat)%vex(m11,m31,m21,m41,ll)=pawtab(itypat)%vex(m11,m31,m21,m41,ll)*four_pi
         pawtab(itypat)%vex(m21,m31,m11,m41,ll)=pawtab(itypat)%vex(m11,m31,m21,m41,ll)
         pawtab(itypat)%vex(m11,m41,m21,m31,ll)=pawtab(itypat)%vex(m11,m31,m21,m41,ll)
         pawtab(itypat)%vex(m21,m41,m11,m31,ll)=pawtab(itypat)%vex(m11,m31,m21,m41,ll)
        end do
       end do
      end do
     end do
    end do

   end if !useexexch>0

  end if !lcur/=-1
 end do !end loop on typat

 end subroutine pawpuxinit
!!***
