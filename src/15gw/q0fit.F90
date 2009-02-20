!{\src2tex{textfont=tt}}
!!****f* ABINIT/q0fit
!! NAME
!! q0fit
!!
!! FUNCTION
!! Calculate chi0 for qq->0 by a fit on the finite qq
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (MC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nq number of q-points
!!  q  vector of q-points
!!  gvec vector of G-points
!!  nomega number of frequencies
!!  omega vector of frequencies
!!  npwvec number of G points
!!  chi0 matrix of chi^0
!!  qcut max modulis of q for the fit
!!  metal logical true if simulating a metal
!!  nop number of symmetry operations
!!  op matrix of symmetry operations
!!  ninv number of inversion
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  
!! OUTPUT
!!  chi0 matrix of chi^0
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine q0fit(nq,q,gvec,nomega,omega,npwvec,chi0,qcut,metal,&
&          nop,op,ninv,gprimd)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_15gw, except_this_one => q0fit
!End of the abilint section

 implicit none
 
!Arguments ------------------------------------
 integer,intent(in) :: nomega,nq
 integer,intent(in) :: npwvec,ninv,nop
 integer,intent(in) :: gvec(3,npwvec)
 real(dp),intent(in) :: qcut
 real(dp),intent(in) :: q(3,nq),op(3,3,nop)
 real(dp),intent(in) :: gprimd(3,3)
 complex(gwpc),intent(inout) :: chi0(nq,npwvec,npwvec,nomega) 
 complex(gwpc),intent(in) :: omega(nomega)
 logical,intent(in) :: metal
 
!Local variables-------------------------------
 integer :: istat,i1,i2,iinv,iop
 integer :: cg,cg1,cq,cfit,cm,cont,tr,c0
 integer :: nfit, iomega,d
 integer :: mq(nq),totq,gp(3)
 integer, allocatable :: sym(:,:,:)
 real(dp),allocatable :: qpg(:,:,:,:)
 real(dp) :: q0(3),sqm,det,qp(3),qpp(3),ga(3),gb(3) 
 real(dp) :: b1(3),b2(3),b3(3)
 real(dp),allocatable :: qq(:,:),sismat(:,:),sisnotre(:),sisresre(:)
 real(dp),allocatable :: sisnotim(:),sisresim(:)
 real(dp) :: rm1(3,3)
 complex(gwpc) :: newchi, chiw(nomega)
 complex(gwpc),allocatable :: fitchi(:,:)
 
! *************************************************************************

 b1=two_pi*gprimd(:,1)
 b2=two_pi*gprimd(:,2)
 b3=two_pi*gprimd(:,3)


! part I : preliminar computations

 write (*,*) 'ATTENTION: calculation fitting the q->0 term'
 write (7,*) 'ATTENTION: calculation fitting the q->0 term'

 allocate(qpg(3,nq,npwvec,nop*ninv))
 allocate(sym(2,nq,nop*ninv))
 sym (:,:,:)=0
 qpg(:,:,:,:)=1000
 totq=0

      allocate(sismat(7,7),stat=istat)
      if(istat/=0) stop 'out of memory'
      allocate(sisnotre(7),stat=istat)
      if(istat/=0) stop 'out of memory'
      allocate(sisresre(7),stat=istat)
      if(istat/=0) stop 'out of memory'
      allocate(sisnotim(7),stat=istat)
      if(istat/=0) stop 'out of memory'
      allocate(sisresim(7),stat=istat)
      if(istat/=0) stop 'out of memory'


 do cq =1,nq
  if(all(abs(q(:,cq))<1.e-3)) then
   c0=cq
   mq(cq)=1 
   qpg(1,cq,1,1)=q(1,cq)*b1(1)+q(2,cq)*b2(1)+q(3,cq)*b3(1)
   qpg(2,cq,1,1)=q(1,cq)*b1(2)+q(2,cq)*b2(2)+q(3,cq)*b3(2)
   qpg(3,cq,1,1)=q(1,cq)*b1(3)+q(2,cq)*b2(3)+q(3,cq)*b3(3)
   do cg=2,npwvec
    qpg(1,cq,cg,1)=gvec(1,cg)*b1(1)+gvec(2,cg)*b2(1)+gvec(3,cg)*b3(1)
    qpg(2,cq,cg,1)=gvec(1,cg)*b1(2)+gvec(2,cg)*b2(2)+gvec(3,cg)*b3(2)
    qpg(3,cq,cg,1)=gvec(1,cg)*b1(3)+gvec(2,cg)*b2(3)+gvec(3,cg)*b3(3)
   enddo
  else
   cm=1
    do  iop=1,nop 
      do  iinv=1,1 !ninv TODO TO be checked
          call dosym(op(1,1,iop),iinv,q(1,cq),qp(1))
          tr=0 ! counter
          qpp(:)=qp(1)*b1(:)+qp(2)*b2(:)+qp(3)*b3(:)
          do cont =cm,1,-1
!            if (all(qpp(:).eq.qpg(:,cq,1,cont))) tr=1
            if (all(abs(qpp(:)-qpg(:,cq,1,cont))<0.001)) tr=1
          enddo
          if (tr.eq.0) then
            sym(1,cq,cm)=iop
            sym(2,cq,cm)=iinv
            do cg=1,npwvec 

             qpg(1,cq,cg,cm)=(gvec(1,cg)+qp(1))*b1(1)+(gvec(2,cg)+qp(2))*b2(1)+&
&	     (gvec(3,cg)+qp(3))*b3(1)
             qpg(2,cq,cg,cm)=(gvec(1,cg)+qp(1))*b1(2)+(gvec(2,cg)+qp(2))*b2(2)+&
&	     (gvec(3,cg)+qp(3))*b3(2)
             qpg(3,cq,cg,cm)=(gvec(1,cg)+qp(1))*b1(3)+(gvec(2,cg)+qp(2))*b2(3)+&
&	     (gvec(3,cg)+qp(3))*b3(3)
          
           enddo
          cm =cm+1
          endif
      enddo
    enddo
   mq(cq)=cm-1

  endif


!debug
! do cg=1,npwvec
!  do cm=1,mq(cq)
!  write (*,*) 'qpg=',qpg(1,cq,1,cm),qpg(2,cq,1,cm),qpg(3,cq,1,cm)
!   write (*,*) 'q=',q(1,cq),q(2,cq),q(3,cq)
!   write (*,*) 'G=',gvec(1,cg),gvec(2,cg),gvec(3,cg)
!   write (*,*) 'mq=',mq(cq)
!  enddo
! enddo
!end debug
 enddo


 do cq=1,nq
   totq=totq+mq(cq)
 enddo
 write (*,*) 'q0fit: number of q-points=',totq 

! part II : fit the head 
 nfit=0  
  
 do cg =1,npwvec
    do cq =1,nq
     do cm=1,mq(cq)
       if (sqrt(qpg(1,cq,cg,cm)**2+qpg(2,cq,cg,cm)**2+qpg(3,cq,cg,cm)**2).le.qcut) then
       if (cq.ne.c0) then
          if (cg.ne.1) then
             write (*,*) 'debug =',qpg (:,cq,cg,cm),cg,cm
             write (*,*) 'debug =',qpg(:,cq,1,1),cq
             write (*,*) 'bug: too large cut-off'
             stop	  
	  endif
          nfit=nfit+1
       endif
       endif
       if (cq==c0) then
          q0(1)=q(1,cq)*b1(1)+q(2,cq)*b2(1)+q(3,cq)*b3(1)
          q0(2)=q(1,cq)*b1(2)+q(2,cq)*b2(2)+q(3,cq)*b3(2)
          q0(3)=q(1,cq)*b1(3)+q(2,cq)*b2(3)+q(3,cq)*b3(3)   
       endif
     enddo
    enddo
 enddo        

 write(*,*) 'q0fit: I am starting the fit'
! debug
 write (*,*) 'nfit=',nfit,qcut
 write (*,*) 'q0=',q0
!end debug
 
 if (nfit.le.10) then
   write (*,*)'ATTENTION: it is impossible to continue the fit'
   write (*,*)'no enought points nothing will be done'
   write (7,*)'ATTENTION: it is impossible to continue the fit'
   write (7,*)'no enought points nothing will be done'
 
 else      
   allocate(qq(3,nfit),stat=istat)
   if(istat/=0) stop 'out of memory'
   allocate(fitchi(nfit,nomega),stat=istat)
   if(istat/=0) stop 'out of memory'  

!open (31,file="w0")
!open (32,file="wip")         
   cfit=1
   cg=1

   do cq =1,nq
     do cm =1,mq(cq)
         if (sqrt(qpg(1,cq,cg,cm)**2+qpg(2,cq,cg,cm)**2+qpg(3,cq,cg,cm)**2).le.qcut) then           
         if (cfit.le.nfit.and. cq.ne.c0 ) then
            qq(:,cfit)=qpg(:,cq,cg,cm)
	    fitchi(cfit,:)=chi0(cq,cg,cg,:)
	    
!           write (*,*)'caio', qq(1,cfit),qq(2,cfit),qq(3,cfit),real(fitchi(cfit,1))
          cfit=cfit+1
         endif
         endif
      enddo
   enddo        

!close (31)
!close (32)
!debug
!  write (*,*)'q0fit: beginning of frequency cycle'
!   write (*,*)'qq=',qq
!   write (*,*)'chi=',fitchi
!end debug
do iomega=1,nomega   

   if (metal.and.(omega(iomega)).eq.cmplx(0,0)) then
     d=7
    else
      d=6  
    endif

      sismat(:,:)=0.0
      sisnotre(:)=0.0  
      sisnotim(:)=0.0  
      sisresre(:)=0.0  
      sisresim(:)=0.0  
!debugging
!    write (*,*) 'I am initializing the arrays'
!end debugging      
    if (metal.and.(omega(iomega)).eq.cmplx(0,0)) then
     do cq=1,nfit 
      sismat(1,1)=sismat(1,1)+1.
      sismat(1,2)=sismat(1,2)+qq(1,cq)*qq(1,cq)
      sismat(1,3)=sismat(1,3)+qq(1,cq)*qq(2,cq)
      sismat(1,4)=sismat(1,4)+qq(1,cq)*qq(3,cq)
      sismat(1,5)=sismat(1,5)+qq(2,cq)*qq(2,cq)
      sismat(1,6)=sismat(1,6)+qq(2,cq)*qq(3,cq)
      sismat(1,7)=sismat(1,7)+qq(3,cq)*qq(3,cq)
      
      sismat(2,1)=sismat(2,1)+qq(1,cq)*qq(1,cq)
      sismat(2,2)=sismat(2,2)+qq(1,cq)*qq(1,cq)*qq(1,cq)*qq(1,cq)
      sismat(2,3)=sismat(2,3)+qq(1,cq)*qq(2,cq)*qq(1,cq)*qq(1,cq)
      sismat(2,4)=sismat(2,4)+qq(1,cq)*qq(3,cq)*qq(1,cq)*qq(1,cq)
      sismat(2,5)=sismat(2,5)+qq(2,cq)*qq(2,cq)*qq(1,cq)*qq(1,cq)
      sismat(2,6)=sismat(2,6)+qq(2,cq)*qq(3,cq)*qq(1,cq)*qq(1,cq)
      sismat(2,7)=sismat(2,7)+qq(3,cq)*qq(3,cq)*qq(1,cq)*qq(1,cq)
      
      sismat(3,1)=sismat(3,1)+qq(1,cq)*qq(2,cq)
      sismat(3,2)=sismat(3,2)+qq(1,cq)*qq(1,cq)*qq(1,cq)*qq(2,cq)
      sismat(3,3)=sismat(3,3)+qq(1,cq)*qq(2,cq)*qq(1,cq)*qq(2,cq)
      sismat(3,4)=sismat(3,4)+qq(1,cq)*qq(3,cq)*qq(1,cq)*qq(2,cq)
      sismat(3,5)=sismat(3,5)+qq(2,cq)*qq(2,cq)*qq(1,cq)*qq(2,cq)
      sismat(3,6)=sismat(3,6)+qq(2,cq)*qq(3,cq)*qq(1,cq)*qq(2,cq)
      sismat(3,7)=sismat(3,7)+qq(3,cq)*qq(3,cq)*qq(1,cq)*qq(2,cq)
      
      sismat(4,1)=sismat(4,1)+qq(1,cq)*qq(3,cq)
      sismat(4,2)=sismat(4,2)+qq(1,cq)*qq(1,cq)*qq(1,cq)*qq(3,cq)
      sismat(4,3)=sismat(4,3)+qq(1,cq)*qq(2,cq)*qq(1,cq)*qq(3,cq)
      sismat(4,4)=sismat(4,4)+qq(1,cq)*qq(3,cq)*qq(1,cq)*qq(3,cq)
      sismat(4,5)=sismat(4,5)+qq(2,cq)*qq(2,cq)*qq(1,cq)*qq(3,cq)
      sismat(4,6)=sismat(4,6)+qq(2,cq)*qq(3,cq)*qq(1,cq)*qq(3,cq)
      sismat(4,7)=sismat(4,7)+qq(3,cq)*qq(3,cq)*qq(1,cq)*qq(3,cq)
      
      sismat(5,1)=sismat(5,1)+qq(2,cq)*qq(2,cq)
      sismat(5,2)=sismat(5,2)+qq(1,cq)*qq(1,cq)*qq(2,cq)*qq(2,cq)
      sismat(5,3)=sismat(5,3)+qq(1,cq)*qq(2,cq)*qq(2,cq)*qq(2,cq)
      sismat(5,4)=sismat(5,4)+qq(1,cq)*qq(3,cq)*qq(2,cq)*qq(2,cq)
      sismat(5,5)=sismat(5,5)+qq(2,cq)*qq(2,cq)*qq(2,cq)*qq(2,cq)
      sismat(5,6)=sismat(5,6)+qq(2,cq)*qq(3,cq)*qq(2,cq)*qq(2,cq)
      sismat(5,7)=sismat(5,7)+qq(3,cq)*qq(3,cq)*qq(2,cq)*qq(2,cq)

      sismat(6,1)=sismat(6,1)+qq(2,cq)*qq(3,cq)
      sismat(6,2)=sismat(6,2)+qq(1,cq)*qq(1,cq)*qq(2,cq)*qq(3,cq)
      sismat(6,3)=sismat(6,3)+qq(1,cq)*qq(2,cq)*qq(2,cq)*qq(3,cq)
      sismat(6,4)=sismat(6,4)+qq(1,cq)*qq(3,cq)*qq(2,cq)*qq(3,cq)
      sismat(6,5)=sismat(6,5)+qq(2,cq)*qq(2,cq)*qq(2,cq)*qq(3,cq)
      sismat(6,6)=sismat(6,6)+qq(2,cq)*qq(3,cq)*qq(2,cq)*qq(3,cq)
      sismat(6,7)=sismat(6,7)+qq(3,cq)*qq(3,cq)*qq(2,cq)*qq(3,cq)
      
      sismat(7,1)=sismat(7,1)+qq(3,cq)*qq(3,cq)
      sismat(7,2)=sismat(7,2)+qq(1,cq)*qq(1,cq)*qq(3,cq)*qq(3,cq)
      sismat(7,3)=sismat(7,3)+qq(1,cq)*qq(2,cq)*qq(3,cq)*qq(3,cq)
      sismat(7,4)=sismat(7,4)+qq(1,cq)*qq(3,cq)*qq(3,cq)*qq(3,cq)
      sismat(7,5)=sismat(7,5)+qq(2,cq)*qq(2,cq)*qq(3,cq)*qq(3,cq)
      sismat(7,6)=sismat(7,6)+qq(2,cq)*qq(3,cq)*qq(3,cq)*qq(3,cq)
      sismat(7,7)=sismat(7,7)+qq(3,cq)*qq(3,cq)*qq(3,cq)*qq(3,cq)
     enddo
      
       do cq=1,nfit
         sisnotre(1)=sisnotre(1)+real(fitchi(cq,iomega))
	 sisnotre(2)=sisnotre(2)+real(fitchi(cq,iomega))*qq(1,cq)*qq(1,cq)
	 sisnotre(3)=sisnotre(3)+real(fitchi(cq,iomega))*qq(1,cq)*qq(2,cq)
	 sisnotre(4)=sisnotre(4)+real(fitchi(cq,iomega))*qq(1,cq)*qq(3,cq)
	 sisnotre(5)=sisnotre(5)+real(fitchi(cq,iomega))*qq(2,cq)*qq(2,cq)
	 sisnotre(6)=sisnotre(6)+real(fitchi(cq,iomega))*qq(2,cq)*qq(3,cq)
	 sisnotre(7)=sisnotre(7)+real(fitchi(cq,iomega))*qq(3,cq)*qq(3,cq)
       
         sisnotim(1)=sisnotim(1)+aimag(fitchi(cq,iomega))
	 sisnotim(2)=sisnotim(2)+aimag(fitchi(cq,iomega))*qq(1,cq)*qq(1,cq)
	 sisnotim(3)=sisnotim(3)+aimag(fitchi(cq,iomega))*qq(1,cq)*qq(2,cq)
	 sisnotim(4)=sisnotim(4)+aimag(fitchi(cq,iomega))*qq(1,cq)*qq(3,cq)
	 sisnotim(5)=sisnotim(5)+aimag(fitchi(cq,iomega))*qq(2,cq)*qq(2,cq)
	 sisnotim(6)=sisnotim(6)+aimag(fitchi(cq,iomega))*qq(2,cq)*qq(3,cq)
	 sisnotim(7)=sisnotim(7)+aimag(fitchi(cq,iomega))*qq(3,cq)*qq(3,cq)
       
       enddo

!debug
!write (*,*) 'mat=',sismat
!end debug


!      call dete(sismat,d,det,sisnotre,sisresre)
  
      call matrginv(sismat,7,d)
      
!debug
!write(*,*) 'mat=',sismat
!write(*,*) 'not=',sisnotre
!write(*,*) 'not=',sisnotim
!write(*,*) 'resa=',sisresre
!sisresre(:)=0.
!end debug 
     
     
       do i1=1,7
         do i2=1,7
           sisresre(i1)=sisresre(i1)+sismat(i1,i2)*sisnotre(i2)
           sisresim(i1)=sisresim(i1)+sismat(i1,i2)*sisnotim(i2)         
	 enddo
       enddo
      
!debug
  write (*,*) 'result=', sisresre
  write (*,*) 'result=', sisresim
!end debug      
      
        newchi=cmplx(sisresre(1)+sisresre(2)*q0(1)*q0(1)&
&                    +sisresre(5)*q0(2)*q0(2)+sisresre(7)*q0(3)*q0(3)&
&                    +2*sisresre(3)*q0(1)*q0(2)+2*sisresre(4)*q0(1)*q0(3)&
&                    +2*sisresre(6)*q0(2)*q0(3),&
&                    sisresim(1)+sisresim(2)*q0(1)*q0(1)&
&                    +sisresim(5)*q0(2)*q0(2)+sisresim(7)*q0(3)*q0(3)&
&                    +2*sisresim(3)*q0(1)*q0(2)+2*sisresim(4)*q0(1)*q0(3)&
&                    +2*sisresim(6)*q0(2)*q0(3))
        do cq=1,nq
	 do cg=1,npwvec
	  if (all(abs(qpg(:,cq,cg,1))<1.0e-3)) chi0(cq,cg,cg,iomega)=newchi 
	 enddo
	enddo
  
  sqm=0
  do cq=1,nfit
   sqm=sqm+abs(fitchi(cq,iomega)-cmplx(sisresre(1)+sisresre(2)*qq(1,cq)*qq(1,cq)&
&                    +sisresre(5)*qq(2,cq)*qq(2,cq)+sisresre(7)*qq(3,cq)*qq(3,cq)&
&                    +2*sisresre(3)*qq(1,cq)*qq(2,cq)+2*sisresre(4)*qq(1,cq)*qq(3,cq)&
&                    +2*sisresre(6)*qq(2,cq)*qq(3,cq),&
&                    sisresim(1)+sisresim(2)*qq(1,cq)*qq(1,cq)&
&                    +sisresim(5)*qq(2,cq)*qq(2,cq)+sisresim(7)*qq(3,cq)*qq(3,cq)&
&                    +2*sisresim(3)*qq(1,cq)*qq(2,cq)+2*sisresim(4)*qq(1,cq)*qq(3,cq)&
&                    +2*sisresim(6)*qq(2,cq)*qq(3,cq)))**2
  enddo


  chiw(iomega)=newchi 
  write (7,*) 'q0fit: fitted omega=',omega(iomega)
  write (7,*) 'q0fit: fitted chi0=',newchi
  write (7,*) 'q0fit: mean square difference=',sqm/nfit
!  write (7,*) 'q0fit: determinant=',det
  write (*,*) 'q0fit: fitted omega=',omega(iomega)
  write (*,*) 'q0fit: fitted chi0=',newchi      
  write (*,*) 'q0fit: mean square difference=',sqm/nfit
!  write (*,*) 'q0fit: determinant=',det
      
   else

  do cq=1,nfit     
      sismat(1,1)=sismat(1,1)+qq(1,cq)*qq(1,cq)*qq(1,cq)*qq(1,cq)
      sismat(1,2)=sismat(1,2)+qq(1,cq)*qq(2,cq)*qq(1,cq)*qq(1,cq)
      sismat(1,3)=sismat(1,3)+qq(1,cq)*qq(3,cq)*qq(1,cq)*qq(1,cq)
      sismat(1,4)=sismat(1,4)+qq(2,cq)*qq(2,cq)*qq(1,cq)*qq(1,cq)
      sismat(1,5)=sismat(1,5)+qq(2,cq)*qq(3,cq)*qq(1,cq)*qq(1,cq)
      sismat(1,6)=sismat(1,6)+qq(3,cq)*qq(3,cq)*qq(1,cq)*qq(1,cq)
      
      sismat(2,1)=sismat(2,1)+qq(1,cq)*qq(1,cq)*qq(1,cq)*qq(2,cq)
      sismat(2,2)=sismat(2,2)+qq(1,cq)*qq(2,cq)*qq(1,cq)*qq(2,cq)
      sismat(2,3)=sismat(2,3)+qq(1,cq)*qq(3,cq)*qq(1,cq)*qq(2,cq)
      sismat(2,4)=sismat(2,4)+qq(2,cq)*qq(2,cq)*qq(1,cq)*qq(2,cq)
      sismat(2,5)=sismat(2,5)+qq(2,cq)*qq(3,cq)*qq(1,cq)*qq(2,cq)
      sismat(2,6)=sismat(2,6)+qq(3,cq)*qq(3,cq)*qq(1,cq)*qq(2,cq)
      
      sismat(3,1)=sismat(3,1)+qq(1,cq)*qq(1,cq)*qq(1,cq)*qq(3,cq)
      sismat(3,2)=sismat(3,2)+qq(1,cq)*qq(2,cq)*qq(1,cq)*qq(3,cq)
      sismat(3,3)=sismat(3,3)+qq(1,cq)*qq(3,cq)*qq(1,cq)*qq(3,cq)
      sismat(3,4)=sismat(3,4)+qq(2,cq)*qq(2,cq)*qq(1,cq)*qq(3,cq)
      sismat(3,5)=sismat(3,5)+qq(2,cq)*qq(3,cq)*qq(1,cq)*qq(3,cq)
      sismat(3,6)=sismat(3,6)+qq(3,cq)*qq(3,cq)*qq(1,cq)*qq(3,cq)
            
      sismat(4,1)=sismat(4,1)+qq(1,cq)*qq(1,cq)*qq(2,cq)*qq(2,cq)
      sismat(4,2)=sismat(4,2)+qq(1,cq)*qq(2,cq)*qq(2,cq)*qq(2,cq)
      sismat(4,3)=sismat(4,3)+qq(1,cq)*qq(3,cq)*qq(2,cq)*qq(2,cq)
      sismat(4,4)=sismat(4,4)+qq(2,cq)*qq(2,cq)*qq(2,cq)*qq(2,cq)
      sismat(4,5)=sismat(4,5)+qq(2,cq)*qq(3,cq)*qq(2,cq)*qq(2,cq)
      sismat(4,6)=sismat(4,6)+qq(3,cq)*qq(3,cq)*qq(2,cq)*qq(2,cq)

      sismat(5,1)=sismat(5,1)+qq(1,cq)*qq(1,cq)*qq(2,cq)*qq(3,cq)
      sismat(5,2)=sismat(5,2)+qq(1,cq)*qq(2,cq)*qq(2,cq)*qq(3,cq)
      sismat(5,3)=sismat(5,3)+qq(1,cq)*qq(3,cq)*qq(2,cq)*qq(3,cq)
      sismat(5,4)=sismat(5,4)+qq(2,cq)*qq(2,cq)*qq(2,cq)*qq(3,cq)
      sismat(5,5)=sismat(5,5)+qq(2,cq)*qq(3,cq)*qq(2,cq)*qq(3,cq)
      sismat(5,6)=sismat(5,6)+qq(3,cq)*qq(3,cq)*qq(2,cq)*qq(3,cq)
      
      sismat(6,1)=sismat(6,1)+qq(1,cq)*qq(1,cq)*qq(3,cq)*qq(3,cq)
      sismat(6,2)=sismat(6,2)+qq(1,cq)*qq(2,cq)*qq(3,cq)*qq(3,cq)
      sismat(6,3)=sismat(6,3)+qq(1,cq)*qq(3,cq)*qq(3,cq)*qq(3,cq)
      sismat(6,4)=sismat(6,4)+qq(2,cq)*qq(2,cq)*qq(3,cq)*qq(3,cq)
      sismat(6,5)=sismat(6,5)+qq(2,cq)*qq(3,cq)*qq(3,cq)*qq(3,cq)
      sismat(6,6)=sismat(6,6)+qq(3,cq)*qq(3,cq)*qq(3,cq)*qq(3,cq)
   enddo

       do cq=1,nfit

	 sisnotre(1)=sisnotre(1)+real(fitchi(cq,iomega))*qq(1,cq)*qq(1,cq)
	 sisnotre(2)=sisnotre(2)+real(fitchi(cq,iomega))*qq(1,cq)*qq(2,cq)
	 sisnotre(3)=sisnotre(3)+real(fitchi(cq,iomega))*qq(1,cq)*qq(3,cq)
	 sisnotre(4)=sisnotre(4)+real(fitchi(cq,iomega))*qq(2,cq)*qq(2,cq)
	 sisnotre(5)=sisnotre(5)+real(fitchi(cq,iomega))*qq(2,cq)*qq(3,cq)
	 sisnotre(6)=sisnotre(6)+real(fitchi(cq,iomega))*qq(3,cq)*qq(3,cq)
       
	 sisnotim(1)=sisnotim(1)+aimag(fitchi(cq,iomega))*qq(1,cq)*qq(1,cq)
	 sisnotim(2)=sisnotim(2)+aimag(fitchi(cq,iomega))*qq(1,cq)*qq(2,cq)
	 sisnotim(3)=sisnotim(3)+aimag(fitchi(cq,iomega))*qq(1,cq)*qq(3,cq)
	 sisnotim(4)=sisnotim(4)+aimag(fitchi(cq,iomega))*qq(2,cq)*qq(2,cq)
	 sisnotim(5)=sisnotim(5)+aimag(fitchi(cq,iomega))*qq(2,cq)*qq(3,cq)
	 sisnotim(6)=sisnotim(6)+aimag(fitchi(cq,iomega))*qq(3,cq)*qq(3,cq)
       
       enddo
!debug
!write (*,*) 'mat=',sismat
!end debug

!      call dete(sismat,d,det,sisnotre,sisresre)
     
      call matrginv(sismat,7,d)
      
!debug
!write (*,*) 'mat=',sismat
!write(*,*) 'not=',sisnotre
!write(*,*) 'not=',sisnotim
!write(*,*) 'resa=',sisresre
!sisresre(:)=0.
!end debug

       do i1=1,6
         do i2=1,6
           sisresre(i1)=sisresre(i1)+sismat(i1,i2)*sisnotre(i2)
           sisresim(i1)=sisresim(i1)+sismat(i1,i2)*sisnotim(i2)         
	 enddo
       enddo
      
!debug
  write (*,*) 'result=',sisresre
  write (*,*) 'result=',sisresim
!end debug


        newchi=cmplx(sisresre(1)*q0(1)*q0(1)+sisresre(4)*q0(2)*q0(2)&
&                    +sisresre(6)*q0(3)*q0(3)&
&                    +2*sisresre(2)*q0(1)*q0(2)+2*sisresre(3)*q0(1)*q0(3)&
&                    +2*sisresre(5)*q0(2)*q0(3),&
&                    sisresim(1)*q0(1)*q0(1)+sisresim(4)*q0(2)*q0(2)&
&                    +sisresim(6)*q0(3)*q0(3)&
&                    +2*sisresim(2)*q0(1)*q0(2)+2*sisresim(3)*q0(1)*q0(3)&
&                    +2*sisresim(5)*q0(2)*q0(3))
        do cq=1,nq
	 do cg=1,npwvec
	  if (all(abs(qpg(:,cq,cg,1))<1.0e-3)) chi0(cq,cg,cg,iomega)=newchi 
	 enddo
	enddo  

  sqm=0
  do cq=1,nfit
   sqm=sqm+abs(fitchi(cq,iomega)-cmplx(sisresre(1)*qq(1,cq)*qq(1,cq)&
&                    +sisresre(4)*qq(2,cq)*qq(2,cq)+sisresre(6)*qq(3,cq)*qq(3,cq)&
&                    +2*sisresre(2)*qq(1,cq)*qq(2,cq)+&
&                    2*sisresre(3)*qq(1,cq)*qq(3,cq)&
&                    +2*sisresre(5)*qq(2,cq)*qq(3,cq),&
&                    sisresim(1)*qq(1,cq)*qq(1,cq)&
&                    +sisresim(4)*qq(2,cq)*qq(2,cq)+sisresim(6)*qq(3,cq)*qq(3,cq)&
&                    +2*sisresim(2)*qq(1,cq)*qq(2,cq)+&
&                    2*sisresim(3)*qq(1,cq)*qq(3,cq)&
&                    +2*sisresim(5)*qq(2,cq)*qq(3,cq)))**2
  enddo

  chiw(iomega)=newchi
  write (7,*) 'q0fit: fitted omega=',omega(iomega)
  write (7,*) 'q0fit: fitted chi0=',newchi
  write (7,*) 'q0fit: mean square difference=',sqm/nfit
!  write (7,*) 'q0fit: determinant=',det
  write (*,*) 'q0fit: fitted omega=',omega(iomega)
  write (*,*) 'q0fit: fitted chi0=',newchi
  write (*,*) 'q0fit: mean square difference=',sqm/nfit
!  write (*,*) 'q0fit: determinant=',det
            
   endif

   


  enddo
   deallocate (fitchi)
   deallocate (qq)
   write (7,'(A, f7.4,A, i5,A, e13.6,A, e13.6)') &
&   'q0fit: summary ',qcut/sqrt(b1(1)**2+b1(2)**2+b1(3)**2),' '&
&   ,nfit,' ',real(chiw(1)),' ',real(chiw(2))
   write (*,'(A, f7.4,A, i5,A, e13.6,A, e13.6)') &
&   'q0fit: summary ',qcut/sqrt(b1(1)**2+b1(2)**2+b1(3)**2),' '&
&   ,nfit,' ',real(chiw(1)),' ',real(chiw(2))


 endif    


! part III : fit the wings
  if (nfit.gt.10) then
  if(.true.) then
! fit the wings
!  c0
   allocate(qq(3,nfit),stat=istat)
   if(istat/=0) stop 'out of memory'
   allocate(fitchi(nfit,nomega),stat=istat)
   if(istat/=0) stop 'out of memory'
! part III.1 : wing 0G
   do cg=2,npwvec
    cfit=1 
 
     do cq=1,nq
     do cm =1,mq(cq)
         if (sqrt(qpg(1,cq,1,cm)**2+qpg(2,cq,1,cm)**2+qpg(3,cq,1,cm)**2).le.qcut) then
         if (cfit.le.nfit.and. cq.ne.c0 ) then
! look for right G'
          rm1(:,:)=op(:,:,sym(1,cq,cm))
          call matrginv(rm1,3,3)
          ga(:)=real(gvec(:,cg))
          call dosym(rm1(1,1),sym(2,cq,cm),ga(1),gb(1))
          gp(:)=int(gb(:))          
          do cg1=1,npwvec
             if (all(gp(:)==gvec(:,cg1))) exit   
          enddo
            qq(:,cfit)=qpg(:,cq,1,cm)
            fitchi(cfit,:)=chi0(cq,1,cg1,:)
!            write (*,*)'pluto',sym(1,cq,cm),sym(2,cq,cm)
!            write (*,*)'pippo',gvec(:,cg),gp,gvec(:,cg1)
!            write (*,*)'ciao', qq(1,cfit),qq(2,cfit),qq(3,cfit),real(fitchi(cfit,1))
            cfit=cfit+1
         endif
         endif
      enddo
      enddo

! faccio il fit
do iomega=1,nomega
 write (*,*)' ci provo'
   if (metal.and.(omega(iomega)).eq.cmplx(0,0)) then
     d=4
   else
     d=3
   endif 
      sismat(:,:)=0.0
      sisnotre(:)=0.0
      sisnotim(:)=0.0
      sisresre(:)=0.0
      sisresim(:)=0.0
!debugging
!    write (*,*) 'I am initializing the arrays'
!end debugging
    if (metal.and.(omega(iomega)).eq.cmplx(0,0)) then
     do cq=1,nfit
      sismat(1,1)=sismat(1,1)+1.
      sismat(1,2)=sismat(1,2)+qq(1,cq)
      sismat(1,3)=sismat(1,3)+qq(2,cq)
      sismat(1,4)=sismat(1,4)+qq(3,cq)

      sismat(2,1)=sismat(2,1)+qq(1,cq)
      sismat(2,2)=sismat(2,2)+qq(1,cq)*qq(1,cq)
      sismat(2,3)=sismat(2,3)+qq(1,cq)*qq(2,cq)
      sismat(2,4)=sismat(2,4)+qq(1,cq)*qq(3,cq)

      sismat(3,1)=sismat(3,1)+qq(2,cq)
      sismat(3,2)=sismat(3,2)+qq(2,cq)*qq(1,cq)
      sismat(3,3)=sismat(3,3)+qq(2,cq)*qq(2,cq)
      sismat(3,4)=sismat(3,4)+qq(2,cq)*qq(3,cq)

      sismat(4,1)=sismat(4,1)+qq(3,cq)
      sismat(4,2)=sismat(4,2)+qq(3,cq)*qq(1,cq)
      sismat(4,3)=sismat(4,3)+qq(3,cq)*qq(2,cq)
      sismat(4,4)=sismat(4,4)+qq(3,cq)*qq(3,cq)
    enddo
       do cq=1,nfit
         sisnotre(1)=sisnotre(1)+real(fitchi(cq,iomega))
         sisnotre(2)=sisnotre(2)+real(fitchi(cq,iomega))*qq(1,cq)
         sisnotre(3)=sisnotre(3)+real(fitchi(cq,iomega))*qq(2,cq)
         sisnotre(4)=sisnotre(4)+real(fitchi(cq,iomega))*qq(3,cq)

         sisnotim(1)=sisnotim(1)+aimag(fitchi(cq,iomega))
         sisnotim(2)=sisnotim(2)+aimag(fitchi(cq,iomega))*qq(1,cq)
         sisnotim(3)=sisnotim(3)+aimag(fitchi(cq,iomega))*qq(2,cq)
         sisnotim(4)=sisnotim(4)+aimag(fitchi(cq,iomega))*qq(3,cq)

       enddo
      write (*,*)'1'
       call matrginv(sismat,7,d)
      write (*,*)'2'
      do i1=1,4
         do i2=1,4
           sisresre(i1)=sisresre(i1)+sismat(i1,i2)*sisnotre(i2)
           sisresim(i1)=sisresim(i1)+sismat(i1,i2)*sisnotim(i2)
         enddo
       enddo

!debug
!  write (*,*) 'result=',sisresre
!  write (*,*) 'result=',sisresim
!end debug


        newchi=cmplx(sisresre(1)+sisresre(2)*q0(1)+sisresre(3)*q0(2)&
&                    +sisresre(4)*q0(3),sisresim(1)+sisresim(2)*q0(1)&
&                    +sisresim(3)*q0(2)+sisresim(4)*q0(3))
        do cq=1,nq
          if (all(abs(qpg(:,cq,1,1))<1.0e-3)) chi0(cq,1,cg,iomega)=newchi
        enddo


 sqm=0
  do cq=1,nfit
   sqm=sqm+abs(fitchi(cq,iomega)-cmplx(sisresre(1)+sisresre(2)*qq(1,cq)&
&                    +sisresre(3)*qq(2,cq)+sisresre(4)*qq(3,cq),&
&                    sisresim(1)+sisresim(2)*qq(1,cq)&
&                    +sisresim(3)*qq(2,cq)+sisresim(4)*qq(3,cq)))**2
  enddo
  write (*,*)'fatto'
  chiw(iomega)=newchi
  write (7,*) 'q0fit: fitted omega=',omega(iomega),'G=0G=',gvec(:,cg)
  write (7,*) 'q0fit: fitted chi0=',newchi
  write (7,*) 'q0fit: mean square difference=',sqm/nfit
!  write (7,*) 'q0fit: determinant=',det
  write (*,*) 'q0fit: fitted omega=',omega(iomega),'G=0G=',gvec(:,cg)
  write (*,*) 'q0fit: fitted chi0=',newchi
  write (*,*) 'q0fit: mean square difference=',sqm/nfit


    else
     do cq=1,nfit

      sismat(1,1)=sismat(1,1)+qq(1,cq)*qq(1,cq)
      sismat(1,2)=sismat(1,2)+qq(1,cq)*qq(2,cq)
      sismat(1,3)=sismat(1,3)+qq(1,cq)*qq(3,cq)

      sismat(2,1)=sismat(2,1)+qq(2,cq)*qq(1,cq)
      sismat(2,2)=sismat(2,2)+qq(2,cq)*qq(2,cq)
      sismat(2,3)=sismat(2,3)+qq(2,cq)*qq(3,cq)

      sismat(3,1)=sismat(3,1)+qq(3,cq)*qq(1,cq)
      sismat(3,2)=sismat(3,2)+qq(3,cq)*qq(2,cq)
      sismat(3,3)=sismat(3,3)+qq(3,cq)*qq(3,cq)
    enddo
       do cq=1,nfit
         sisnotre(1)=sisnotre(1)+real(fitchi(cq,iomega))*qq(1,cq)
         sisnotre(2)=sisnotre(2)+real(fitchi(cq,iomega))*qq(2,cq)
         sisnotre(3)=sisnotre(3)+real(fitchi(cq,iomega))*qq(3,cq)

         sisnotim(1)=sisnotim(1)+aimag(fitchi(cq,iomega))*qq(1,cq)
         sisnotim(2)=sisnotim(2)+aimag(fitchi(cq,iomega))*qq(2,cq)
         sisnotim(3)=sisnotim(3)+aimag(fitchi(cq,iomega))*qq(3,cq)

       enddo

       call matrginv(sismat,7,d)
      do i1=1,3
         do i2=1,3
           sisresre(i1)=sisresre(i1)+sismat(i1,i2)*sisnotre(i2)
           sisresim(i1)=sisresim(i1)+sismat(i1,i2)*sisnotim(i2)
         enddo
       enddo

!debug
!  write (*,*) 'result=',sisresre
!  write (*,*) 'result=',sisresim
!end debug


        newchi=cmplx(sisresre(1)*q0(1)+sisresre(2)*q0(2)&
&                    +sisresre(3)*q0(3),sisresim(1)*q0(1)&
&                    +sisresim(2)*q0(2)+sisresim(3)*q0(3))
       do cq=1,nq
          if (all(abs(qpg(:,cq,1,1))<1.0e-3)) chi0(cq,1,cg,iomega)=newchi
        enddo

 sqm=0
  do cq=1,nfit
   sqm=sqm+abs(fitchi(cq,iomega)-cmplx(sisresre(1)*qq(1,cq)&
&                    +sisresre(2)*qq(2,cq)+sisresre(3)*qq(3,cq),&
&                    sisresim(1)*qq(1,cq)&
&                    +sisresim(2)*qq(2,cq)+sisresim(3)*qq(3,cq)))**2
  enddo

  chiw(iomega)=newchi
  write (7,*) 'q0fit: fitted omega=',omega(iomega),'G=0G=',gvec(:,cg)
  write (7,*) 'q0fit: fitted chi0=',newchi
  write (7,*) 'q0fit: mean square difference=',sqm/nfit
!  write (7,*) 'q0fit: determinant=',det
  write (*,*) 'q0fit: fitted omega=',omega(iomega),'G=0G=',gvec(:,cg)
  write (*,*) 'q0fit: fitted chi0=',newchi
  write (*,*) 'q0fit: mean square difference=',sqm/nfit


    endif


 enddo


  enddo

! part III.2 : wing G0
   do cg=2,npwvec
    cfit=1

     do cq=1,nq
     do cm =1,mq(cq)
         if (sqrt(qpg(1,cq,1,cm)**2+qpg(2,cq,1,cm)**2+qpg(3,cq,1,cm)**2).le.qcut) then
         if (cfit.le.nfit.and. cq.ne.c0 ) then
! look for right G'
          rm1(:,:)=op(:,:,sym(1,cq,cm))
          call matrginv(rm1,3,3)
          ga(:)=real(gvec(:,cg)) 
          call dosym(rm1(1,1),sym(2,cq,cm),ga(1),gb(1))
          gp(:)=int(gb(:))
          do cg1=1,npwvec
             if (all(gp(:)==gvec(:,cg1))) exit
          enddo
            qq(:,cfit)=qpg(:,cq,1,cm)
            fitchi(cfit,:)=chi0(cq,cg1,1,:)
!            write (*,*)'ciao', qq(1,cfit),qq(2,cfit),qq(3,cfit),real(fitchi(cfit,1))
            cfit=cfit+1
         endif
         endif
      enddo
      enddo

! faccio il fit
do iomega=1,nomega

   if (metal.and.(omega(iomega)).eq.cmplx(0,0)) then
     d=4
    else
      d=3
    endif
      sismat(:,:)=0.0
      sisnotre(:)=0.0
      sisnotim(:)=0.0
      sisresre(:)=0.0
      sisresim(:)=0.0
!debugging
!    write (*,*) 'I am initializing the arrays'
!end debugging
    if (metal.and.(omega(iomega)).eq.cmplx(0,0)) then
     do cq=1,nfit
      sismat(1,1)=sismat(1,1)+1.
      sismat(1,2)=sismat(1,2)+qq(1,cq)
      sismat(1,3)=sismat(1,3)+qq(2,cq)
      sismat(1,4)=sismat(1,4)+qq(3,cq)

      sismat(2,1)=sismat(2,1)+qq(1,cq)
      sismat(2,2)=sismat(2,2)+qq(1,cq)*qq(1,cq)
      sismat(2,3)=sismat(2,3)+qq(1,cq)*qq(2,cq)
      sismat(2,4)=sismat(2,4)+qq(1,cq)*qq(3,cq)

      sismat(3,1)=sismat(3,1)+qq(2,cq)
      sismat(3,2)=sismat(3,2)+qq(2,cq)*qq(1,cq)
      sismat(3,3)=sismat(3,3)+qq(2,cq)*qq(2,cq)
      sismat(3,4)=sismat(3,4)+qq(2,cq)*qq(3,cq)

      sismat(4,1)=sismat(4,1)+qq(3,cq)
      sismat(4,2)=sismat(4,2)+qq(3,cq)*qq(1,cq)
      sismat(4,3)=sismat(4,3)+qq(3,cq)*qq(2,cq)
      sismat(4,4)=sismat(4,4)+qq(3,cq)*qq(3,cq)
    enddo
       do cq=1,nfit
         sisnotre(1)=sisnotre(1)+real(fitchi(cq,iomega))
         sisnotre(2)=sisnotre(2)+real(fitchi(cq,iomega))*qq(1,cq)
         sisnotre(3)=sisnotre(3)+real(fitchi(cq,iomega))*qq(2,cq)
         sisnotre(4)=sisnotre(4)+real(fitchi(cq,iomega))*qq(3,cq)

         sisnotim(1)=sisnotim(1)+aimag(fitchi(cq,iomega))
         sisnotim(2)=sisnotim(2)+aimag(fitchi(cq,iomega))*qq(1,cq)
         sisnotim(3)=sisnotim(3)+aimag(fitchi(cq,iomega))*qq(2,cq)
         sisnotim(4)=sisnotim(4)+aimag(fitchi(cq,iomega))*qq(3,cq)

       enddo

       call matrginv(sismat,7,d)

      do i1=1,4
         do i2=1,4
           sisresre(i1)=sisresre(i1)+sismat(i1,i2)*sisnotre(i2)
           sisresim(i1)=sisresim(i1)+sismat(i1,i2)*sisnotim(i2)
         enddo
       enddo

!debug
!  write (*,*) 'result=',sisresre
!  write (*,*) 'result=',sisresim
!end debug


        newchi=cmplx(sisresre(1)+sisresre(2)*q0(1)+sisresre(3)*q0(2)&
&                    +sisresre(4)*q0(3),sisresim(1)+sisresim(2)*q0(1)&
&                    +sisresim(3)*q0(2)+sisresim(4)*q0(3))
        do cq=1,nq
          if (all(abs(qpg(:,cq,1,1))<1.0e-3)) chi0(cq,cg,1,iomega)=newchi
        enddo


 sqm=0
  do cq=1,nfit
   sqm=sqm+abs(fitchi(cq,iomega)-cmplx(sisresre(1)+sisresre(2)*qq(1,cq)&
&                    +sisresre(3)*qq(2,cq)+sisresre(4)*qq(3,cq),&
&                    sisresim(1)+sisresim(2)*qq(1,cq)&
&                    +sisresim(3)*qq(2,cq)+sisresim(4)*qq(3,cq)))**2
  enddo

  chiw(iomega)=newchi
  write (7,*) 'q0fit: fitted omega=',omega(iomega),'G=0G=',gvec(:,cg)
  write (7,*) 'q0fit: fitted chi0=',newchi
  write (7,*) 'q0fit: mean square difference=',sqm/nfit
!  write (7,*) 'q0fit: determinant=',det
  write (*,*) 'q0fit: fitted omega=',omega(iomega),'G=0G=',gvec(:,cg)
  write (*,*) 'q0fit: fitted chi0=',newchi
  write (*,*) 'q0fit: mean square difference=',sqm/nfit


    else
     do cq=1,nfit

      sismat(1,1)=sismat(1,1)+qq(1,cq)*qq(1,cq)
      sismat(1,2)=sismat(1,2)+qq(1,cq)*qq(2,cq)
      sismat(1,3)=sismat(1,3)+qq(1,cq)*qq(3,cq)

      sismat(2,1)=sismat(2,1)+qq(2,cq)*qq(1,cq)
      sismat(2,2)=sismat(2,2)+qq(2,cq)*qq(2,cq)
      sismat(2,3)=sismat(2,3)+qq(2,cq)*qq(3,cq)

      sismat(3,1)=sismat(3,1)+qq(3,cq)*qq(1,cq)
      sismat(3,2)=sismat(3,2)+qq(3,cq)*qq(2,cq)
      sismat(3,3)=sismat(3,3)+qq(3,cq)*qq(3,cq)
    enddo
       do cq=1,nfit
         sisnotre(1)=sisnotre(1)+real(fitchi(cq,iomega))*qq(1,cq)
         sisnotre(2)=sisnotre(2)+real(fitchi(cq,iomega))*qq(2,cq)
         sisnotre(3)=sisnotre(3)+real(fitchi(cq,iomega))*qq(3,cq)

         sisnotim(1)=sisnotim(1)+aimag(fitchi(cq,iomega))*qq(1,cq)
         sisnotim(2)=sisnotim(2)+aimag(fitchi(cq,iomega))*qq(2,cq)
         sisnotim(3)=sisnotim(3)+aimag(fitchi(cq,iomega))*qq(3,cq)

       enddo
      call matrginv(sismat,7,d)
      do i1=1,3
         do i2=1,3
           sisresre(i1)=sisresre(i1)+sismat(i1,i2)*sisnotre(i2)
           sisresim(i1)=sisresim(i1)+sismat(i1,i2)*sisnotim(i2)
         enddo
       enddo

!debug
!  write (*,*) 'result=',sisresre
!  write (*,*) 'result=',sisresim
!end debug


        newchi=cmplx(sisresre(1)*q0(1)+sisresre(2)*q0(2)&
&                    +sisresre(3)*q0(3),sisresim(1)*q0(1)&
&                    +sisresim(2)*q0(2)+sisresim(3)*q0(3))
       do cq=1,nq
          if (all(abs(qpg(:,cq,1,1))<1.0e-3)) chi0(cq,cg,1,iomega)=newchi
        enddo

 sqm=0
  do cq=1,nfit
   sqm=sqm+abs(fitchi(cq,iomega)-cmplx(sisresre(1)*qq(1,cq)&
&                    +sisresre(2)*qq(2,cq)+sisresre(3)*qq(3,cq),&
&                    sisresim(1)*qq(1,cq)&
&                    +sisresim(2)*qq(2,cq)+sisresim(3)*qq(3,cq)))**2
  enddo

  chiw(iomega)=newchi
  write (7,*) 'q0fit: fitted omega=',omega(iomega),'G=0G=',gvec(:,cg)
  write (7,*) 'q0fit: fitted chi0=',newchi
  write (7,*) 'q0fit: mean square difference=',sqm/nfit
!  write (7,*) 'q0fit: determinant=',det
  write (*,*) 'q0fit: fitted omega=',omega(iomega),'G=0G=',gvec(:,cg)
  write (*,*) 'q0fit: fitted chi0=',newchi
  write (*,*) 'q0fit: mean square difference=',sqm/nfit


    endif


 enddo


  enddo




   deallocate (fitchi)
   deallocate (qq)











  else
! wing=0
  write (ab_out,*)'ATTENTION: wing at the moment placed equal to 0'
  write (*,*)'ATTENTION: wing at the moment placed equal to 0'
  
  do cq=1,nq       
      do cg1=1,npwvec
	 do cg=1,npwvec	    
            if (all(abs(q(:,cq))<1.0e-3).and.all(gvec(:,cg).ne.gvec(:,cg1))) then
	      chi0(cq,cg1,cg,:)=cmplx(0.0,0.0)
	    endif
	 enddo  
      enddo
  enddo



 endif	    
 endif
   deallocate (sismat)
   deallocate (sisnotre)
   deallocate (sisresre)
   deallocate (sisnotim)
   deallocate (sisresim)

 deallocate (qpg,sym)

 end subroutine
!!***
