!{\src2tex{textfont=tt}}
!!****f* ABINIT/secinit
!! NAME
!! secinit
!!
!! FUNCTION
!! Subroutine that computes the first projections of the W initial states
!! over the free energy states
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  eigval(nqpt,3*natom)=array containing the eigenvalues, in cm-1
!!  eigvect(nqpt,3*natom,natom,3,2)=array containing the eigenvectors:
!!  f_,g_,z_ subsp = indexes of the bands in the free,global,frozen window
!!  gsize = number of global bands
!!  indwnz(nqpt,nwnn)= index of the W states in the frozen window
!!  ingss(nwnn,natom,3,2)=initial guess of the W states, user-input
!!  iqpt = inex of the current q point
!!  maxqsize = maximum number of bands in the global window
!!  natom = number of atoms in the unit cell
!!  nqpt = number of q points
!!  nwnn = number of LWF
!!  nwnz(nqpt)= number of Wannier states in the frozen window
!!  qpoint(nqpt,3)= list of the q points in fractional coordinates
!!  xred = atomic fractional coordinates
!!  znucl = atomic mass numbers
!!  zsize = (=qsize(:,3)) number of frozen bands
!!
!! OUTPUT
!!  lambda(nqpt,nwnn,maxqsize,2)= overlaps (real) over all states of the initial LWF
!!
!! NOTE
!!  To be re-written to do all q points in the same time
!!  alfa = projection coefficients: alpha(mm,jj)= < u_jj | gin_mm >
!!
!! PARENTS
!!      lwf
!!
!! CHILDREN
!!      zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine secinit(eigval,eigvect,f_subsp,g_subsp,z_subsp,&
&gsize,indwnz,ingss,iqpt,lambda,maxqsize,natom,nqpt,nwnn,nwnz,qpoint,&
&xred,znucl,zsize)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: gsize,iqpt,maxqsize,natom,nqpt,nwnn,zsize
!arrays
 integer,intent(in) :: f_subsp(nqpt,3*natom),g_subsp(nqpt,3*natom),indwnz(nwnn)
 integer,intent(in) :: nwnz(nqpt),z_subsp(nqpt,3*natom)
 real(dp),intent(in) :: eigval(nqpt,3*natom),eigvect(nqpt,3*natom,natom,3,2)
 real(dp),intent(in) :: ingss(nwnn,natom,3,2),qpoint(nqpt,3),xred(natom,3)
 real(dp),intent(in) :: znucl(natom)
 real(dp),intent(out) :: lambda(nqpt,nwnn,maxqsize,2)

!Local variables-------------------------------
!scalars
 integer :: aa,i1,i2,iatom,iband,ier,ii,jj,kk,mm,ok_z_flag,okalmode,tao,tmpcont
 real(dp) :: normvect,normwann,phi,rcov,testnorm,tmpvar
 character(len=2) :: symbol
 character(len=500) :: message
!arrays
 integer :: ilist(gsize)
 real(dp) :: PQ(gsize,gsize,2),P_G(gsize,gsize,2),QPQ(gsize,gsize,2)
 real(dp) :: Q_inner(gsize,gsize,2),alfa(nwnn,maxqsize,2),eigvalqpq(gsize)
 real(dp) :: eigvecqpq(2,gsize,gsize),eigwan(3),gin(nwnn,natom,3,2)
 real(dp) :: tmplist(gsize)
 real(dp),allocatable :: matrx(:,:),zhpev1(:,:),zhpev2(:)

!******************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zhpev
#endif

!DEBUG
!write(*,*) ' secinit : enter'
!write(*,*)
!write(*,*) ' current q point',iqpt,qpoint(iqpt,:)
!write(*,*) ' nwnn=',nwnn
!write(*,*) ' nwnz(iqpt)=',nwnz(iqpt)
!write(*,*) ' gsize=',gsize
!write(*,*) ' zsize=',zsize
!ENDDEBUG

!imposition of the phase factor, like in the case of the triatomic chain
 phi=2*pi*(qpoint(iqpt,1)+qpoint(iqpt,2)+qpoint(iqpt,3))
 do mm=1,nwnn
  do tao=1,natom
   do aa=1,3
!   gin(mm,tao,aa,1)=ingss(mm,tao,aa,1)*cos(2*pi*qpoint(iqpt,aa))+ingss(mm,tao,aa,2)*sin(2*pi*qpoint(iqpt,aa))
!   gin(mm,tao,aa,2)=ingss(mm,tao,aa,2)*cos(2*pi*qpoint(iqpt,aa))-ingss(mm,tao,aa,1)*sin(2*pi*qpoint(iqpt,aa))
    gin(mm,tao,aa,1)=ingss(mm,tao,aa,1)*cos(phi)+ingss(mm,tao,aa,2)*sin(phi)
    gin(mm,tao,aa,2)=ingss(mm,tao,aa,2)*cos(phi)-ingss(mm,tao,aa,1)*sin(phi)
!   write(*,'(a,2f16.11)') 'gin:',gin(mm,tao,aa,1),gin(mm,tao,aa,2)
   end do
  end do
 end do



!determination of lambdas, as:
!alfa_m = < u_j | gin >
!for each of the gin components (1..m..N_m)
!<u_j| are the states in the global window

 do mm=1,nwnn               !loop over LWF states
  tmpvar=zero
  normwann=zero
  do tao=1,natom           !loop over atoms, tao
   do aa=1,3               !loop over directions, alpha,beta,gamma
    normwann=normwann+gin(mm,tao,aa,1)*gin(mm,tao,aa,1)+gin(mm,tao,aa,2)*gin(mm,tao,aa,2)
   end do                   !directions
  end do                    !atoms
! DEBUG
! write(*,'(a,i4,a,f16.11)') 'LWF state=',mm,'normwann',normwann
! ENDDEBUG

  do jj=1,gsize             !loop over global bands
   alfa(mm,jj,1)=zero
   alfa(mm,jj,2)=zero
   normvect=zero
   do tao=1,natom           !loop over atoms, tao
    do aa=1,3               !loop over directions, alpha,beta,gamma
     alfa(mm,jj,1)=alfa(mm,jj,1)+eigvect(iqpt,g_subsp(iqpt,jj),tao,aa,1)*gin(mm,tao,aa,1)+&
&     eigvect(iqpt,g_subsp(iqpt,jj),tao,aa,2)*gin(mm,tao,aa,2)
     alfa(mm,jj,2)=alfa(mm,jj,2)+eigvect(iqpt,g_subsp(iqpt,jj),tao,aa,1)*gin(mm,tao,aa,2)-&
&     eigvect(iqpt,g_subsp(iqpt,jj),tao,aa,2)*gin(mm,tao,aa,1)
     normvect=normvect+eigvect(iqpt,g_subsp(iqpt,jj),tao,aa,1)*eigvect(iqpt,g_subsp(iqpt,jj),tao,aa,1)+&
&     eigvect(iqpt,g_subsp(iqpt,jj),tao,aa,2)*eigvect(iqpt,g_subsp(iqpt,jj),tao,aa,2)
    end do                   !directions
   end do                    !atoms
!  write(*,'(a,2f16.11,a,f16.11)') 'alfas and norms',alfa(mm,jj,1),alfa(mm,jj,2),' norm of eigvect',normvect
   alfa(mm,jj,1)=alfa(mm,jj,1)/dble(sqrt(normvect*normwann))
   alfa(mm,jj,2)=alfa(mm,jj,2)/dble(sqrt(normvect*normwann))
!  write(*,'(a,3i4,3f16.11)') ' The projections: ',iqpt,g_subsp(iqpt,jj),mm,alfa(mm,jj,1),alfa(mm,jj,2),&
!  alfa(mm,jj,1)*alfa(mm,jj,1)+alfa(mm,jj,2)*alfa(mm,jj,2)
   tmpvar=tmpvar+alfa(mm,jj,1)*alfa(mm,jj,1)+alfa(mm,jj,2)*alfa(mm,jj,2)
  end do                     !global bands
! write(*,*) ' Sum of the projections: ',tmpvar
 end do                      !LWF states


!start the computation of the Q_inner*P_G*Q_inner

!DEBUG
!write(*,*) ' we are in secinit,start the computation of the Q_inner*P_G*Q_inner'
!ENDDEBUG

!initialize matrices
 Q_inner(:,:,:)=zero
 P_G(:,:,:)=zero
 PQ(:,:,:)=zero
 QPQ(:,:,:)=zero

!compute Q_inner
 do ii=1,gsize
  if (f_subsp(iqpt,ii)>0) then
   Q_inner(ii,ii,1)=one
  end if
 end do

!DEBUG
!write(*,*) ' we are in secinit, after compute Q_inner'
!ENDDEBUG

!compute P_G= Sum_wann < alfa_ii,mm | alfa_jj,mm>
 do ii=1,gsize           ! loop over global bands
  do jj=1,gsize          ! loop over global bands
   do mm=1,nwnn          ! loop over Wannier states
    P_G(ii,jj,1)=P_G(ii,jj,1)+alfa(mm,ii,1)*alfa(mm,jj,1)+alfa(mm,ii,2)*alfa(mm,jj,2)
    P_G(ii,jj,2)=P_G(ii,jj,2)+alfa(mm,ii,2)*alfa(mm,jj,1)-alfa(mm,ii,1)*alfa(mm,jj,2)
   end do                 !wannier states
!  write(*,'(a,2i4,2f16.11)') 'P_G ',ii,jj,P_G(ii,jj,1),P_G(ii,jj,2)
  end do                  !global bands
 end do                   !global bands

!DEBUG
!write(*,*) ' we are in secinit, after compute P_G'
!ENDDEBUG
!check hermicity of P_G
 do ii=1,gsize           ! loop over global bands
  do jj=1,gsize          ! loop over global bands
   if (abs(P_G(ii,jj,1)-P_G(jj,ii,1))>0.000000001 .and. &
&   abs(P_G(ii,jj,2)+P_G(jj,ii,2))>0.000000001) then
    write (*,*) ' P_G is not hermitian :'
    write (*,'(2i4,2f16.12,a,2f16.12)') ii,jj,P_G(ii,jj,1),P_G(ii,jj,2),'and',P_G(jj,ii,1),P_G(jj,ii,2)
   end if
  end do
 end do

!compute P_G*Q_inner
 do ii=1,gsize           ! loop over lines
  do jj=1,gsize          ! loop over columns
   do kk=1,gsize         ! loop over elements
    PQ(jj,ii,1)=PQ(jj,ii,1)+P_G(jj,kk,1)*Q_inner(kk,ii,1)-P_G(jj,kk,2)*Q_inner(kk,ii,2)
    PQ(jj,ii,2)=PQ(jj,ii,2)+P_G(jj,kk,2)*Q_inner(kk,ii,1)+P_G(jj,kk,1)*Q_inner(kk,ii,2)
   end do
  end do
 end do

!DEBUG
!write(*,*) ' we are in secinit, after compute  P_G*Q_inner'
!ENDDEBUG

!compute Q_inner*PQ_inner
 do ii=1,gsize           ! loop over lines
  do jj=1,gsize          ! loop over columns
   do kk=1,gsize         ! loop over elements
    QPQ(jj,ii,1)=QPQ(jj,ii,1)+Q_inner(jj,kk,1)*PQ(kk,ii,1)-Q_inner(jj,kk,2)*PQ(kk,ii,2)
    QPQ(jj,ii,2)=QPQ(jj,ii,2)+Q_inner(jj,kk,2)*PQ(kk,ii,1)+Q_inner(jj,kk,1)*PQ(kk,ii,2)
   end do
  end do
 end do

!check hermicity of Q_inner*PQ_inner
 do ii=1,gsize           ! loop over global bands
  do jj=1,gsize          ! loop over global bands
   if (abs(QPQ(ii,jj,1)-QPQ(jj,ii,1))>0.000000001 .and. &
&   abs(QPQ(ii,jj,2)+QPQ(jj,ii,2))>0.000000001) then
    write (*,*) ' QPQ is not hermitian :'
    write (*,'(2i4,2f16.12,a,2f16.12)') ii,jj,QPQ(ii,jj,1),QPQ(ii,jj,2),'and',QPQ(jj,ii,1),QPQ(jj,ii,2)
   end if
  end do
 end do

!DEBUG
!write(*,*) ' we are in secinit, after compute Q_inner*PQ_inner '
!write(*,*) ' QPQ is :'
!do ii=1,gsize
!do jj=1,gsize
!write(*,*) QPQ(ii,jj,1),QPQ(ii,jj,2)
!end do
!write(*,*)
!end do
!write(*,*) ' and we try to do some diagonalization ... '
!ENDDEBUG

 ier=0
 ii=1
 allocate(matrx(2,gsize*(gsize+1)/2))
 do i2=1,gsize
  do i1=1,i2
   matrx(1,ii)=QPQ(i2,i1,1)
   matrx(2,ii)=QPQ(i2,i1,2)
   ii=ii+1
  end do
 end do
 allocate(zhpev1(2,2*gsize-1),zhpev2(3*gsize-2))
 call zhpev('V','U',gsize,matrx,eigvalqpq,eigvecqpq,gsize,zhpev1,zhpev2,ier)
 deallocate(matrx,zhpev1,zhpev2)

!DEBUG
!write(*,*) ' secinit : eigvalues of the QPQ matrix'
!do ii=1,gsize
!write(*,*) 'for the eigenvalue',eigvalqpq(ii)
!do jj=1,gsize
!write(*,'(2f12.7)') eigvecqpq(1,jj,ii),eigvecqpq(2,jj,ii)
!end do
!end do
!ENDDEBUG

!write(*,'(a,i5,a,i4,a)') 'at q ',iqpt,' we choose',nwnn-nwnz(iqpt),' Wannier states outside the inner window'

!DEBUG
!if (nwnz(iqpt)<nwnn) then
!!write(*,*) ' testing the norm of the eigenvector of QPQ matrix'
!do ii=1,gsize
!testnorm=0
!! write(*,*) 'corresponding eigvalue',eigvalqpq(ii)
!! write(*,*) 'eigenvector:'
!do jj=1,gsize
!!  write(*,'(a,2f12.8)') 'eigvecqpq, real, imag:',eigvecqpq(1,jj,ii),eigvecqpq(2,jj,ii)
!testnorm=testnorm+eigvecqpq(1,jj,ii)*eigvecqpq(1,jj,ii)+eigvecqpq(2,jj,ii)*eigvecqpq(2,jj,ii)
!end do
!! write(*,*) 'testnorm of eigvectors of QPQ in secinit =',testnorm
!if (abs(testnorm-1.0_dp)>0.1d-8) write(*,*) 'attention to testnorm, ERROR >1',testnorm
!end do
!end if
!ENDDEBUG

!get the first nwnz(iqpt) values corresponding to the eigenvectors

 lambda(iqpt,:,:,:)=zero
 aa=1
 do ii=1,nwnz(iqpt)
! write(*,*) 'aa=',aa
  do kk=aa,gsize
!  write(*,*) 'z_subsp(iqpt,kk)=',z_subsp(iqpt,kk),' kk=',kk
   if (z_subsp(iqpt,kk)>0) then
    lambda(iqpt,ii,kk,1)=one
    aa=kk+1
    exit
   end if
  end do
 end do

 do ii=nwnz(iqpt)+1,nwnn         ! loop over N-Mq states (free Wannier states)
! write(*,*) 'corresponding eigval of QPQ is',eigvalqpq(gsize-ii+1+nwnz(iqpt))
  do kk=1,gsize                  ! loop over global bands
   if (f_subsp(iqpt,kk)>0) then
    lambda(iqpt,ii,kk,1)=eigvecqpq(1,kk,gsize-ii+1+nwnz(iqpt))
    lambda(iqpt,ii,kk,2)=eigvecqpq(2,kk,gsize-ii+1+nwnz(iqpt))
   else
    lambda(iqpt,ii,kk,1)=zero
    lambda(iqpt,ii,kk,2)=zero
   end if
  end do                          !bands
 end do                           !free Wannier states

 eigwan=zero

 do ii=1,nwnn
! write(*,*) 'Wannier state',ii
  do kk=1,gsize
!  write(*,'(a,4f11.6)') 'lambda,eigval',lambda(iqpt,ii,kk,:),lambda(iqpt,ii,kk,1)*lambda(iqpt,ii,kk,1)+&
!  lambda(iqpt,ii,kk,2)*lambda(iqpt,ii,kk,2),eigval(iqpt,g_subsp(iqpt,kk))
   eigwan(ii)=eigwan(ii)+(lambda(iqpt,ii,kk,1)*lambda(iqpt,ii,kk,1)+&
&   lambda(iqpt,ii,kk,2)*lambda(iqpt,ii,kk,2))*eigval(iqpt,g_subsp(iqpt,kk))
  end do
 end do
!DEBUG
!write(*,'(a,3f10.5,3f15.7)') 'secinit eigwan ',qpoint(iqpt,:),eigwan(:)
!ENDDEBUG

 do ii=nwnz(iqpt)+1,nwnn                         ! loop over Wannier
  testnorm=zero
  do jj=1,gsize                       ! loop over global bands
   testnorm=testnorm+lambda(iqpt,ii,jj,1)*lambda(iqpt,ii,jj,1)+lambda(iqpt,ii,jj,2)*lambda(iqpt,ii,jj,2)
  end do
  if (abs(testnorm-1.0_dp)>1.0d-08) write(*,*) 'ERROR with lambda norm at W state',ii,testnorm
 end do                                !Wannier states


!DEBUG
!write(*,*) ' secinit : exit'
!ENDDEBUG

end subroutine secinit
!!***
