!{\src2tex{textfont=tt}}
!!****f* ABINIT/smatrix_pawinit
!! NAME
!! smatrix_pawinit
!!
!! FUNCTION
!! Routine which computes paw part of the overlap used to compute LMWF wannier
!!  functions and berryphase
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (BAmadon,FJollet,PHermet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cprj(natom,nspinor*mband*mkmem*nsppol)= <p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                          and each |p_lmn> non-local projector
!!  g1(3)= reciprocal vector to put k1+b inside the BZ. bb=k2-k1=b-G1
!!  ("b" is the true b, so we have to correct bb with G1).
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  mband=maximum number of bands
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nkpt=number of k points.
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  usepaw= 1: use paw framework. 0:do not use paw.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  cm2: Inside sphere part of the overlap needed for constructing wannier function
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      berryphase_new,mlwfovlp
!!
!! CHILDREN
!!      initylmr,leave_new,sbf8,simp_gen,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine smatrix_pawinit(cm2,cprj,ikpt1,ikpt2,&
& g1,gprimd,kpt,mband,mbandw,mkmem,&
& natom,nattyp,nband,&
& nkpt,nspinor,nsppol,ntypat,pawang,pawrad,pawtab,rprimd,&
& typat,usepaw,xred)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13paw, except_this_one => smatrix_pawinit
 use interfaces_13psp
!End of the abilint section

 implicit none

!Arguments---------------------------
! scalars
! arrays
! datatypes
!scalars
 integer,intent(in) :: ikpt1,ikpt2,mband,mbandw,mkmem,natom,nkpt,nspinor,nsppol
 integer,intent(in) :: ntypat,usepaw
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: nattyp(ntypat),nband(nsppol*nkpt),typat(natom)
 real(dp),intent(in) :: gprimd(3,3),kpt(3,nkpt),rprimd(3,3),xred(3,natom)
 real(dp),intent(inout) :: cm2(2,mbandw,mbandw)
 type(cprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables---------------------------
! scalars
! arrays
!scalars
 integer :: iatom,iband,iband1,iband2,icat,icg1,icg2,ii,ikpt,ilmn,ind,intot,ir
 integer :: isel,ispinor,isppol,itypat,j0lmn,jj,jlmn,klm,klmn,kln,ll,lm0,lmax
 integer :: lmin,lmn_size,mesh_size,mm,nband_k
 real(dp) :: a1,a2,arg,b1,b2,bnorm,c1,c2,delta,intg,ppi,ppr,qijbtemp,qijtot,x1
 real(dp) :: x2,xsum,xtemp,xx,yy,zz
 character(len=500) :: message
!arrays
 integer :: g1(3)
 real(dp),parameter :: ili(7)=(/0.d0,-1.d0,0.d0,1.d0,0.d0,-1.d0,0.d0/)
 real(dp),parameter :: ilr(7)=(/1.d0,0.d0,-1.d0,0.d0,1.d0,0.d0,-1.d0/)
 real(dp) :: bb(3),bb1(3),bbn(3),qijb(2),xcart(3,natom)
 real(dp),allocatable :: ff(:),j_bessel(:,:),ylmb(:),ylmrgr_dum(:,:,:)
!no_abirules
 complex(dpc) :: e1,e2,e3

! *************************************************************************

!Test
 if (nspinor==2) then
  write(message,'(a,a,a,a,i3,a,a)')ch10,&
&  ' smatrix_pawinit: ERROR - ',ch10,&
&  '   Not compatible with nspinor=2 !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!
!write(6,*) "compute PAW overlap for k-points",ikpt1,ikpt2
 do iatom=1,natom
  xcart(:,iatom)=rprimd(:,1)*xred(1,iatom)+&
&  rprimd(:,2)*xred(2,iatom)+&
&  rprimd(:,3)*xred(3,iatom)
 end do
!write(6,*) "xcart", xcart
!write(6,*) "xred",xred
!write(6,*) "rprimd",rprimd


!!!!!!!!!!!!!!!!!
!--- Compute intermediate quantities: "b" vector=k2-k1 and its
!normalized value: bbn (and its norm: bnorm)
!compute also Ylm(b).
 allocate(ylmb(pawang%l_size_max*pawang%l_size_max))
 allocate(ylmrgr_dum(1,1,0))
 bb(:)=kpt(:,ikpt2)-kpt(:,ikpt1)+g1(:)
 bb1=bb
!write(6,*) nint(0.5), "nint"
!if((abs(bb(1))-half).le.tol10) then
!bb(1)=bb(1)
!else
!bb(1)=bb(1)-nint(bb(1))
!endif
!if((abs(bb(2))-half).le.tol10) then
!bb(2)=bb(2)
!else
!bb(2)=bb(2)-nint(bb(2))
!endif
!if((abs(bb(3))-half).le.tol10) then
!bb(3)=bb(3)
!else
!bb(3)=bb(3)-nint(bb(3))
!endif
!bb(:)=0.d0
!write(6,*) "bb apres mod",bb(:)
 xx=gprimd(1,1)*bb(1)+gprimd(1,2)*bb(2)+gprimd(1,3)*bb(3)
 yy=gprimd(2,1)*bb(1)+gprimd(2,2)*bb(2)+gprimd(2,3)*bb(3)
 zz=gprimd(3,1)*bb(1)+gprimd(3,2)*bb(2)+gprimd(3,3)*bb(3)
!write(6,*) xx,yy,zz
 bnorm=two_pi*dsqrt(xx**2+yy**2+zz**2)
 if(bnorm<0.00000001) then
! write(6,*) "WARNING: bnorm=",bnorm
  bbn(:)=0.d0
 else
  xx=xx*two_pi
  yy=yy*two_pi
  zz=zz*two_pi
  bb(1)=xx
  bb(2)=yy
  bb(3)=zz
  bbn(1)=xx/bnorm
  bbn(2)=yy/bnorm
  bbn(3)=zz/bnorm
 end if
!debug  bbn=0
!debug  bnorm=0
!bbn has to ne normalized
 call initylmr(pawang%l_size_max,0,1,(/one/),1,bbn(:),ylmb(:),ylmrgr_dum)
!write(6,*) "ylmb(:)",ylmb(:)
!write(6,*) pawang%l_size_max
!write(6,*) "bbn",bbn(:)
!write(6,*) "xx,yy,zz",xx,yy,zz
!write(6,*) "bnorm",bnorm
 deallocate(ylmrgr_dum)

!------- First Compute Qij(b)-
 cm2=zero
 isppol=1
 ispinor=1
 icg1=0
 do iatom=1,natom
! write(6,*) " iatom", iatom
  itypat=typat(iatom)
  lmn_size=pawtab(itypat)%lmn_size
! ---  en coordonnnes reelles cartesiennes (espace reel)
! ---  first radial part(see pawinit)
  mesh_size=pawrad(itypat)%mesh_size
  allocate(j_bessel(mesh_size,pawang%l_size_max))
  icg1=0
  do ii=1,nsppol
   do jj=1,ikpt1-1
    icg1=icg1+nband(jj+(ii-1)*nkpt)
   end do
  end do
  icg2=0
  do ii=1,nsppol
   do jj=1,ikpt2-1
    icg2=icg2+nband(jj+(ii-1)*nkpt)
   end do
  end do
! ---  compute bessel function for (br) for all angular momenta necessary
! ---  and for all value of r.
! ---  they are needed for radial part
! ---  of the integration => j_bessel(ir,:)
  do ir=1,mesh_size
   arg=bnorm*pawrad(itypat)%rad(ir)
   call sbf8(pawang%l_size_max,arg,j_bessel(ir,:))
  end do
! do jlmn=1,pawang%l_size_max
! write(665,*) "j_bessel",j_bessel(1:mesh_size,jlmn)
! enddo
! write(6,*) "bessel function computed"
! ---  Compute \Sum b.R=xsum for future use
  xtemp=0.d0
  do mm=1,3
   xtemp=xtemp+xred(mm,iatom)*bb1(mm)
  end do
  xtemp=xtemp*two_pi
  xsum=0.d0
  do mm=1,3
   xsum=xsum+xcart(mm,iatom)*bbn(mm)*bnorm
  end do
! write(6,*)'xsum',xsum,xtemp,lmn_size

! ---  Loop on jlmn and ilmn
  qijtot=zero
  do jlmn=1,lmn_size
   j0lmn=jlmn*(jlmn-1)/2
   do ilmn=1,jlmn
    klmn=j0lmn+ilmn
    klm=pawtab(itypat)%indklmn(1,klmn);kln=pawtab(itypat)%indklmn(2,klmn)
    lmin=pawtab(itypat)%indklmn(3,klmn);lmax=pawtab(itypat)%indklmn(4,klmn)
!   ---  Sum over angular momenta
!   ---  compute radial part integration for each angular momentum => intg
!   ---  (3j) symbols follows the rule: l belongs to abs(li-lj), li+lj.
    qijb=zero
    do ll=lmin,lmax,2
     lm0=ll*ll+ll+1
     allocate(ff(mesh_size))
     ff(1:mesh_size)=(pawtab(itypat)%phiphj(1:mesh_size,kln)&
&     -pawtab(itypat)%tphitphj(1:mesh_size,kln))&
&     *j_bessel(1:mesh_size,ll+1)
     call simp_gen(intg,ff,pawrad(itypat))
     deallocate(ff)
     qijbtemp=0.d0
     do mm=-ll,ll
      isel=pawang%gntselect(lm0+mm,klm)
      if (isel>0) qijbtemp=qijbtemp&
&      +pawang%realgnt(isel)*ylmb(lm0+mm)
     end do ! mm
!    ---     compute angular part with a summation
!    ---     qijb =\sum_{lm} intg(lm)*qijbtemp
     qijb(1)=qijb(1) +intg*qijbtemp*ilr(ll+1)
     qijb(2)=qijb(2) +intg*qijbtemp*ili(ll+1)
!    if(ilmn==jlmn) write(6,*) "intg, qij",intg,qijbtemp
!    write(666,*) "intg, qij",intg,qijbtemp,ll
    end do ! ll
!   ---  Add exp(-i.b*R) for each atom.
    if(ilmn==jlmn) qijtot=qijtot+qijb(1)
!   if(ilmn==jlmn) write(6,*) "qijtot",qijtot
!   write(6666,*) "qijb avant xsum",qijb(1),qijb(2)
    x1=qijb(1)*dcos(-xsum)-qijb(2)*dsin(-xsum)
    x2=qijb(1)*dsin(-xsum)+qijb(2)*dcos(-xsum)
!   x1 x2 necessary to avoid changing qijb(1) before
!   computing qijb(2)
    qijb(1)=x1
    qijb(2)=x2
!   write(6666,*) "qijb apres xsum",qijb(1),qijb(2)
!   write(666,*) "ilmn,jlmn",ilmn,jlmn,klmn
!   write(666,*) "lmin,lmax",lmin,lmax
!   write(666,*) "qijb",qijb(1),qijb(2)
!   
!   
!   if(ilmn==jlmn) write(6,*) "qij",jlmn,ilmn,qijb(1),qijb(2)
    do iband1=1,mbandw ! limite inferieure a preciser
     do iband2=1,mbandw
!     write(6,*) "iband2",iband2
!     product of (a1+ia2)*(b1-ib2) (minus sign because conjugated)
      ppr=&
!     real part a_1*b_1+a_2*b_2
&      cprj(iatom,iband1+icg1)%cp(1,ilmn)*cprj(iatom,iband2+icg2)%cp(1,jlmn)+&
&      cprj(iatom,iband1+icg1)%cp(2,ilmn)*cprj(iatom,iband2+icg2)%cp(2,jlmn)+&
!     add term on the other triangle  of the matrix
!     qij is the same for this part because phi are real.
&      cprj(iatom,iband1+icg1)%cp(1,jlmn)*cprj(iatom,iband2+icg2)%cp(1,ilmn)+&
&      cprj(iatom,iband1+icg1)%cp(2,jlmn)*cprj(iatom,iband2+icg2)%cp(2,ilmn)
      ppi=&
!     imaginary part a_1*b_2-a_2*b_1
&      cprj(iatom,iband1+icg1)%cp(1,ilmn)*cprj(iatom,iband2+icg2)%cp(2,jlmn)-&
&      cprj(iatom,iband1+icg1)%cp(2,ilmn)*cprj(iatom,iband2+icg2)%cp(1,jlmn)+&
!     add term on the other triangle  of the matrix
&      cprj(iatom,iband1+icg1)%cp(1,jlmn)*cprj(iatom,iband2+icg2)%cp(2,ilmn)-&
&      cprj(iatom,iband1+icg1)%cp(2,jlmn)*cprj(iatom,iband2+icg2)%cp(1,ilmn)
!     delta: diagonal terms are counted twice ! so
!     we need a 0.5 factor for diagonal elements.
      delta=1.d0
!     write(6,*) "ppr  and ppi computed",ikpt1,ikpt2,iband1,iband2
      if(ilmn==jlmn) delta=0.5d0
      cm2(1,iband1,iband2)= cm2(1,iband1,iband2)+ &
&      (qijb(1)*ppr-qijb(2)*ppi)*delta
      cm2(2,iband1,iband2)= cm2(2,iband1,iband2)+ &
&      (qijb(2)*ppr+qijb(1)*ppi)*delta
     end do ! iband2
    end do ! iband1
   end do ! ilmn
  end do ! jlmn
! write(6,*) "final qijtot",qijtot
  deallocate(j_bessel)
 end do ! iatom
 deallocate(ylmb)

!write(6,*) "smatrix_pawinit end"

 end subroutine    smatrix_pawinit
!!***

