!{\src2tex{textfont=tt}}
!!****f* ABINIT/recursion
!! NAME
!! recursion
!! 
!! FUNCTION
!! This routine computes the recursion coefficients and the corresponding 
!! continued fraction to get the density at a point from a fixed potential. 
!! 
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  exppot=exponential of -1/tsmear*vtrial (computed only once in vtorhorec)
!!  coordx, coordy, coordz=coordonnees of the computed point
!!  nrec=order of recursion
!!  fermie=fermi energy (Hartree)
!!  tsmear=temperature (Hartree)
!!  trotter=trotter parameter
!!  ZT_p=fourier transform of the Green krenel (computed only once in vtorhorec)
!!  tol=tolerance criteria for stopping recursion
!!  get_rec_coef=indicate which calcul should be done : 
!!     0 : a and b2 are input, only the continued fraction is calculated, the density is computed in rho_out
!!     1 : a and b2 are computed, then the density is computed in rho_out
!!     2 : a and b2 are computed, then the kinetic energy is computed in rho_out
!!  prtvol=priting volume
!!  mpi_enreg=information about MPI paralelisation
!!  nfft=number of points in FFT grid
!!  ngfft=information about FFT
!!  inf_ucvol=infinitesimal unit cell volume
!!  tim_fourdp=time counter for fourdp
!! 
!! OUTPUT
!!  rho_out=result of the continued fraction multiplied by a multiplicator
!! 
!! SIDE EFFECTS
!!  a, b2 : coefficient given by recursion. Input if get_rec_coef=0, output else
!! 
!! PARENTS
!!      vtorhorec, musolve
!! 
!! CHILDREN
!!      fft3d,wrtout,fourdp,timab
!! 
!! NOTES
!!  at this time :
!!       - exppot should be replaced by ?
!!       - coord should be replaced by ?
!!       - need a rectangular box (rmet diagonal matrix)
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine recursion(exppot,coordx,coordy,coordz,an,bn2, &
&                             rho_out, &
&                             nrec,fermie,tsmear,trotter, &
&                             ZT_p, tol, &
&                             get_rec_coef,prtvol,&
&                             mpi_enreg,nfft,ngfft,rmet,inf_ucvol,tim_fourdp)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
!End of the abilint section

 implicit none

!Arguments -------------------------------
! real(dp), intent(in) :: exppot(0:ngfft(1)-1,0:ngfft(2)-1,0:ngfft(3)-1)
! real(dp), intent(in) :: ZT_p(1:2, 0:ngfft(1)-1,0:ngfft(2)-1,0:ngfft(3)-1)
!scalars
 integer,intent(in) :: coordx,coordy,coordz,get_rec_coef,nfft,nrec,prtvol
 integer,intent(in) :: tim_fourdp,trotter
 real(dp),intent(in) :: fermie,inf_ucvol,tol,tsmear
 real(dp),intent(out) :: rho_out
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: ZT_p(1:2,0:ngfft(1)-1,0:ngfft(2)-1,0:ngfft(3)-1)
 real(dp),intent(in) :: exppot(0:ngfft(1)-1,0:ngfft(2)-1,0:ngfft(3)-1)
 real(dp),intent(in) :: rmet(3,3)
 real(dp),intent(inout) :: an(0:nrec),bn2(0:nrec)

!Local variables-------------------------------
!not used, debugging purpose only
!for debugging purpose, detailled printing only once for density and ekin
!scalars
 integer,parameter :: level=7,minrec=3
 integer,save :: first(1:2)=0
 integer :: cplex,i_trotter,ii,irec,isign,jj,kk,modi,modj,modk
 real(dp) :: arg,bb,beta,mult,prod_b2,rtrotter,switch
 complex(dpc) :: zj
 character(len=500) :: message
!arrays
 real(dp) :: erreur(0:nrec),tsec(2)
 real(dp),allocatable :: Imun(:,:,:),Imunold(:,:,:),Imvn(:,:,:),Zu(:,:,:,:)
 real(dp),allocatable :: Zvtempo(:,:),switchimu(:,:,:),switchu(:,:,:),un(:,:,:)
 real(dp),allocatable :: unold(:,:,:),vn(:,:,:),vtempo(:)
 complex(dpc) :: acc_rho(0:nrec)
 complex(dpc),allocatable :: D(:),Dnew(:),Dold(:),N(:),Nnew(:),Nold(:)

! *************************************************************************

 call timab(186+get_rec_coef,1,tsec)
 
!DEBUG
!!$first = 0 
!!$write(6,*)' '
!!$write(6,*)'enter recursion : echo input variables'
!!$write(6,*)'exppot', exppot(1,1,:)
!!$write(6,*)'coord',coordx,coordy,coordz
!!$write(6,*)'nrec,fermie,tsmear,trotter',nrec,fermie,tsmear,trotter
!!$write(6,*)'ZT_p', Zt_p(1,1,1,:)
!!$write(6,*)'get_rec_coef,prtvol',get_rec_coef,prtvol
!!$write(6,*)'nfft,ngfft',nfft,ngfft
!!$write(6,*)'rmet(:,1)',rmet(:,1)
!!$write(6,*)'inf_ucvol,tim_fourdp',inf_ucvol,tim_fourdp
!ENDDEBUG

!structured debugging if prtvol=-level : print detailled result the first time we enter entropyrec
 if(prtvol==-level)then
  if(first(get_rec_coef)==0)then
   write(message,'(a)')' ' 
   call wrtout(06,message,'COLL')
   write(message,'(a,i6)')' recursion : enter with get_rec_coef = ', get_rec_coef
   call wrtout(06,message,'COLL')
  end if
 end if
 
!initialisation
 if(get_rec_coef/=0)then 
  cplex = get_rec_coef
  allocate(un(0:ngfft(1)-1,0:ngfft(2)-1,0:ngfft(3)-1),vn(0:ngfft(1)-1,0:ngfft(2)-1,0:ngfft(3)-1))
  allocate(unold(0:ngfft(1)-1,0:ngfft(2)-1,0:ngfft(3)-1))
  allocate(switchu(0:ngfft(1)-1,0:ngfft(2)-1,0:ngfft(3)-1),switchimu(0:ngfft(1)-1,0:ngfft(2)-1,0:ngfft(3)-1))
  allocate(Zu(1:2, 0:ngfft(1)-1,0:ngfft(2)-1,0:ngfft(3)-1))
  allocate(vtempo(0:cplex*nfft-1),Zvtempo(1:2,0:nfft-1))
  if(get_rec_coef==2)then  !kinetic energy computation, un and vn must be complex
   allocate(Imun(0:ngfft(1)-1,0:ngfft(2)-1,0:ngfft(3)-1), Imvn(0:ngfft(1)-1,0:ngfft(2)-1,0:ngfft(3)-1))
   allocate(Imunold(0:ngfft(1)-1,0:ngfft(2)-1,0:ngfft(3)-1))
  end if
 end if
 
 beta = 1/tsmear
 if (trotter == 0) then
  rtrotter = 0.5d0
  allocate(D(0:0), N(0:0), Dold(0:0))
  allocate(Nold(0:0), Dnew(0:0), Nnew(0:0))
 else
  rtrotter = real(trotter,dp)
  allocate(D(0:2*trotter-1), N(0:2*trotter-1), Dold(0:2*trotter-1))
  allocate(Nold(0:2*trotter-1), Dnew(0:2*trotter-1), Nnew(0:2*trotter-1))
 end if
 
 if(get_rec_coef/=0)then 
  an = 0.d0
  bn2 = 0.d0
  bn2(0) = 1.d0
  bb = 0.d0
  vn = 0.d0
  unold = 0.d0
  if(get_rec_coef == 1)then !u0 is a Dirac function
   un = 0.d0
   un(coordx,coordy,coordz) = 1/inf_ucvol**0.5d0
  elseif(get_rec_coef==2)then !u0 is a planevawe
   do ii=0,ngfft(1)-1
    do jj=0,ngfft(2)-1
     do kk=0,ngfft(3)-1
      un(ii,jj,kk) = 1/(inf_ucvol*nfft)**0.5d0 * cos(2*pi*( &   !there is no truncation for ekin computation 
&     real(coordx*ii,dp)/real(ngfft(1),dp)+ &
&      real(coordy*jj,dp)/real(ngfft(2),dp)+ &
&      real(coordz*kk,dp)/real(ngfft(3),dp)))
      Imun(ii,jj,kk) = 1/(inf_ucvol*nfft)**0.5d0 * sin(2*pi*( &    !there is no truncation for ekin computation 
&     real(coordx*ii,dp)/real(ngfft(1),dp)+ &
&      real(coordy*jj,dp)/real(ngfft(2),dp)+ &
&      real(coordz*kk,dp)/real(ngfft(3),dp)))
!     false for a non-rectangular cell
     end do
    end do
   end do
   Imvn = 0.d0
   Imunold = 0.d0
  end if
 end if
 acc_rho=dcmplx(0.d0,0.d0)
 N = dcmplx(0.d0,0.d0)
 D = dcmplx(1.d0,0.d0)
!###############################
!calcul d un estimateur d'erreur
 prod_b2 = 1.d0
 erreur = 0.d0
!###############################
 
!##############################################################
!main loop
 maindo : do irec = 0, nrec
  if(get_rec_coef/=0)then
!  #################################################
!  get an and bn2 coef by the lanczos method
   
!  computation of exp(-beta*V/2*p)*un
   vn(:,:,:) = exppot(:,:,:) * un(:,:,:)   !.* in matlab
   if(get_rec_coef==2)then
    Imvn(:,:,:) = exppot(:,:,:) * Imun(:,:,:)
   end if
   
!  convolution with the Green kernel
   do ii = 0, ngfft(1)-1
    do jj = 0, ngfft(2)-1
     do kk = 0, ngfft(3)-1
      vtempo(cplex*(ii+ ngfft(1)*jj+( ngfft(1)* ngfft(2))*kk)) =  vn(ii,jj,kk)
      if(cplex == 2)then
       vtempo(cplex*(ii+ ngfft(1)*jj+( ngfft(1)* ngfft(2))*kk)+1) =  Imvn(ii,jj,kk)
      end if
     end do
    end do
   end do
   isign = -1
   call fourdp(cplex,Zvtempo,vtempo,isign,mpi_enreg,nfft,ngfft,1,tim_fourdp)
   do ii = 0, ngfft(1)-1
    do jj = 0, ngfft(2)-1
     do kk = 0, ngfft(3)-1
      Zu(:,ii,jj,kk) =  Zvtempo(:,ii+ ngfft(1)*jj+( ngfft(1)* ngfft(2))*kk)
     end do
    end do
   end do
   
   switchu(:,:,:) = Zu(1,:,:,:)
   switchimu(:,:,:) = Zu(2,:,:,:)
   Zu(1,:,:,:) = switchu(:,:,:)*ZT_p(1,:,:,:) - switchimu(:,:,:)*ZT_p(2,:,:,:)
   Zu(2,:,:,:) = switchu(:,:,:)*ZT_p(2,:,:,:) + switchimu(:,:,:)*ZT_p(1,:,:,:)
   
   do ii = 0, ngfft(1)-1
    do jj = 0, ngfft(2)-1
     do kk = 0, ngfft(3)-1
      Zvtempo(:,ii+ ngfft(1)*jj+( ngfft(1)* ngfft(2))*kk) = Zu(:,ii,jj,kk)
     end do
    end do
   end do
   isign = 1
   call fourdp(cplex,Zvtempo,vtempo,isign,mpi_enreg,nfft,ngfft,1,tim_fourdp)
   do ii = 0, ngfft(1)-1
    do jj = 0, ngfft(2)-1
     do kk = 0, ngfft(3)-1
      vn(ii,jj,kk)=inf_ucvol*vtempo(cplex*(ii+ ngfft(1)*jj+( ngfft(1)* ngfft(2))*kk))
      if(cplex == 2)then
       Imvn(ii,jj,kk)=inf_ucvol*vtempo(cplex*(ii+ ngfft(1)*jj+( ngfft(1)* ngfft(2))*kk)+1)
      end if
     end do
    end do
   end do
   
!  computation of exp(-beta*V/2*p)*vn
   vn(:,:,:) = exppot(:,:,:) * vn(:,:,:)
   if(get_rec_coef==2)then
    Imvn(:,:,:) = exppot(:,:,:) * Imvn(:,:,:)
   end if
!  multiplication of a and b2 coef by exp(beta*fermie/(2.d0*rtrotter)) must be done in the continued fraction computation
   
!  computation of a and b2
   an(irec) = inf_ucvol*sum(vn(:,:,:)*un(:,:,:))
   if(get_rec_coef==2)then
!   an must be positive real 
    an(irec) = an(irec) + inf_ucvol*sum(Imvn(:,:,:)*Imun(:,:,:))
   end if
   
   if(irec<nrec)then     !we must compute bn2 and prepare for the next iteration
    switchu(:,:,:)=un(:,:,:)
    un(:,:,:)=vn(:,:,:)-an(irec)*un(:,:,:)-bb*unold(:,:,:)
    unold(:,:,:)=switchu(:,:,:)
    bn2(irec+1)=inf_ucvol*sum(un(:,:,:)*un(:,:,:))
    if(get_rec_coef==2)then
     switchu(:,:,:)=Imun(:,:,:)
     Imun(:,:,:)=Imvn(:,:,:)-an(irec)*Imun(:,:,:)-bb*Imunold(:,:,:)
     Imunold(:,:,:)=switchu(:,:,:)
     bn2(irec+1)=bn2(irec+1)+inf_ucvol*sum(Imun(:,:,:)*Imun(:,:,:))
    end if
    bb = sqrt(bn2(irec+1))
    un(:,:,:) = 1/bb*un(:,:,:)
    if(get_rec_coef==2)then
     Imun(:,:,:) = 1/bb*Imun(:,:,:) 
    end if
   end if
   
  end if
  
! ######################################################
! density computation
! density computation is done inside the main looping, juste after the calculus of a and b2, in order to make 
! it possible to stop the recursion at the needed accuracy, without doing more recursion loop than needed - 
! further developpement
  if(get_rec_coef==2)then !kinetic energy computation
   mult = two_pi**two*( &
&   real(coordx**2,dp)/rmet(1,1)+& !   !Gx**2/Lx**2
&  real(coordy**2,dp)/rmet(2,2)+& !   !Gy**2/Ly**2
&  real(coordz**2,dp)/rmet(3,3)) !   !Gy**2/Ly**2
!  false for a non-rectangular cell
  else !density computation
   mult = 1/inf_ucvol * 2.d0   !non-spined system
  end if
  
  if(trotter==0)then
   
   Nnew(0) = (dcmplx(-1.d0 - exp(beta*fermie/(2.d0*rtrotter))*an(irec),0.d0))*N(0) - &
&   dcmplx(exp(beta*fermie/(rtrotter))*bn2(irec),0.d0)*Nold(0)
   Dnew(0) = (dcmplx(-1.d0 - exp(beta*fermie/(2.d0*rtrotter))*an(irec),0.d0))*D(0) - &
&   dcmplx(exp(beta*fermie/(rtrotter))*bn2(irec),0.d0)*Dold(0)
   if(irec==0)then
    Nnew(0) = dcmplx(1.d0,0.d0)
    Dnew(0) = dcmplx( -1.d0 - exp(beta*fermie/(2.d0*rtrotter))*an(0) ,0.d0)
   end if
   Nold(0) = N(0)
   Dold(0) = D(0)
   N(0) = Nnew(0)
   D(0) = Dnew(0)
   acc_rho(irec) = dcmplx(1.d0,0.d0) + N(0)/D(0)
   
!  ###############################
!  calcul d un estimateur d'erreur
   if (irec/=0)then
    prod_b2 = prod_b2 * exp(beta*fermie/(rtrotter)) * bn2(irec)
   end if
   erreur(irec) = abs(prod_b2/(D(0)*Dold(0)))
!  ###############################
   
  else
   
!  ###############################
!  calcul d un estimateur d'erreur
   if (irec/=0)then
    prod_b2 = prod_b2 * exp(beta*fermie/(rtrotter)) * bn2(irec)
   end if
!  ###############################
   
   do i_trotter=0,2*trotter-1
    arg = pi*( 2.d0 * real(i_trotter,dp) + 1.d0 )/( 2.d0 * real(trotter,dp))
    zj = dcmplx( cos(arg) , sin(arg) )
    if(irec==0)then
     Nnew(i_trotter) = dcmplx(1.d0,0.d0)
     Dnew(i_trotter) = zj - dcmplx( exp(beta*fermie/(2.d0*rtrotter))*an(0) ,0.d0)
    else
     Nnew(i_trotter) = (zj - dcmplx(exp(beta*fermie/(2.d0*rtrotter))*an(irec),0.d0))*N(i_trotter) - &
&     dcmplx(exp(beta*fermie/(rtrotter))*bn2(irec),0.d0)*Nold(i_trotter)
     Dnew(i_trotter) = (zj - dcmplx(exp(beta*fermie/(2.d0*rtrotter))*an(irec),0.d0))*D(i_trotter) - &
&     dcmplx(exp(beta*fermie/(rtrotter))*bn2(irec),0.d0)*Dold(i_trotter)
    end if
    Nold(i_trotter) = N(i_trotter)
    Dold(i_trotter) = D(i_trotter)
    N(i_trotter) = Nnew(i_trotter)
    D(i_trotter) = Dnew(i_trotter)
    acc_rho(irec) = acc_rho(irec) + zj*N(i_trotter)/D(i_trotter)
!   ###############################
!   calcul d un estimateur d'erreur
    erreur(irec) = erreur(irec) + abs(prod_b2/(D(i_trotter)*Dold(i_trotter))*2.d0*real(trotter,dp))
!   ###############################
    
   end do
   acc_rho(irec) = dcmplx(1.d0,0.d0) - acc_rho(irec)/dcmplx(2.d0*real(trotter,dp),0.d0)
   
  end if
  
  rho_out = mult *real(acc_rho(irec),dp)
  
  if(irec/=nrec.and.irec>=minrec)then
   if((bn2(irec+1)<tol14).or.(mult*erreur(irec)<tol.and.mult*erreur(irec-1)<tol))then !stop the recursion
    if(get_rec_coef/=0)then 
     bn2(irec+1)=zero
    end if
    exit
   end if
  end if
  
 end do maindo
 
 if(get_rec_coef/=0)then 
  deallocate(un,vn,unold)
  deallocate(switchu,switchimu)
  deallocate(Zu)
  deallocate(vtempo,Zvtempo)
  if(get_rec_coef==2)then
   deallocate(Imun,Imvn,Imunold)
  end if
 end if
 deallocate(D, N, Dold, Nold, Dnew, Nnew)
 
!!$ rho_out = mult *real(acc_rho(nrec),dp)
 
!structured debugging if prtvol=-level : print detailled result the first time we calculate density/kinetic energy
!if(prtvol==-level)then
 if(get_rec_coef/=0)then
  if(first(get_rec_coef)==0)then
   write(message,'(a,3i3)') ' coordonnees',coordx,coordy,coordz
   call wrtout(06,message,'PERS')
   write(message,'(a)')'  irec, densite, erreur_theorique'
   call wrtout(06,message,'PERS')
   do irec=0,nrec
    if(bn2(irec)>tol14)then
     write(message,'(i8,2d10.3)') irec,  0.5d0*mult*real(acc_rho(irec),dp), &
&     mult*erreur(irec)
     call wrtout(06,message,'PERS')
    else
     exit
    end if
   end do
   first(get_rec_coef) = 1
!  DEBUG
!  write(6,*)'a', an
!  write(6,*)'b2', bn2
!  ENDDEBUG
  end if
 end if
!endif

 call timab(186+get_rec_coef,2,tsec)
 
end subroutine recursion
!!***
