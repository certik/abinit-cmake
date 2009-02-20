!{\src2tex{textfont=tt}}
!!****f* ABINIT/getnel
!!
!! NAME
!! getnel
!!
!! FUNCTION
!! Option=1 :
!! Get the total number of electrons nelect, given a trial fermienergy fermie.
!! For this, compute new occupation numbers at each k point,
!! from eigenenergies eigen, according to the
!! smearing scheme defined by occopt (and smearing width tsmear or tphysel).
!!
!! Option=2 :
!! Compute and output the smeared density of states, and the integrated density
!! of states, then write these data
!!
!! Warning : this routine assumes checks have been done in the calling
!! routine, and that the values of the arguments are sensible
!!
!! NOTE
!! in order to speed the calculation, it would be easy to
!! compute the entropy only when the fermi energy is well converged
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG, AF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dosdeltae= DOS delta of Energy (needed if Option=2)
!! eigen(mband*nkpt*nsppol)=eigenvalues (input or init to large number), hartree
!! fermie= fermi energy (Hartree)
!! maxocc=asymptotic maximum occupation number per band
!! mband=maximum number of bands
!! nband(nkpt*nsppol)=number of bands at each k point
!! nkpt=number of k points
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! occopt=option for occupancies, or re-smearing scheme if dblsmr /= 0
!! option=see above
!! tphysel="physical" electronic temperature with FD occupations
!! tsmear=smearing width (or temperature)
!! unitdos=unit number of output of the DOS. Not needed if option==1
!! wtk(nkpt)=k point weights
!!
!! OUTPUT
!! doccde(maxval(nband(:))*nkpt*nsppol)=derivative of occupancies wrt
!!          the energy for each band and k point
!! entropy= entropy associated with the smearing (adimensional)
!! nelect=number of electrons per unit cell
!! occ(maxval(nband(:))*nkpt*nsppol)=occupancies for each band and k point!!
!!
!! NOTES
!! Modified beginning 23/11/2000 by MV
!! Add an additional smearing on top of a FD type, in order to improve k-point
!! convergence: tsmear = 0 and tphysel ~= 2.e-3 corresponds to a small (300K)
!! temperature on the electrons insufficient for convergence purposes.
!! Feed re-smeared "Dirac delta" to the rest of ABINIT with only one parameter,
!! tphysel, which is the physical temperature.
!! encorr = correction to energy for terms of order tsmear^2:
!!       $  E_{phys} = E_{free} - encorr*(E_{int}-E_{free}) + O(tseamr^3)  $
!!
!! PARENTS
!!      clnup1,conducti,conducti_paw,loper3,newocc
!!
!! CHILDREN
!!      dos_hdr_write,leave_new,splfit,spline,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine getnel(doccde,dosdeltae,eigen,entropy,fermie,maxocc,mband,nband,&
&  nelect,nkpt,nsppol,occ,occopt,option,tphysel,tsmear,unitdos,wtk)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_14occeig, except_this_one => getnel
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: mband,nkpt,nsppol,occopt,option,unitdos
 real(dp) :: dosdeltae,fermie,maxocc,tphysel,tsmear
 real(dp),intent(out) :: entropy,nelect
!arrays
 integer :: nband(nkpt*nsppol)
 real(dp) :: eigen(mband*nkpt*nsppol),wtk(nkpt)
 real(dp),intent(out) :: doccde(mband*nkpt*nsppol),occ(mband*nkpt*nsppol)

!Local variables-------------------------------
! dblsmr=flag for re-smearing of delta by type defined by occopt
! nptsdiv2 is the number of integration points, divided by 2.
! nconvd2 = 1/2 number of integration points for convolution
! convlim = limits of convolution integral (in units)
!           corresponds to nconvd2*incconv
! incconv = increment between point on the convolution integral grid
! secmom = second moment of total smearing delta given by occopt.
! smom1  = second moment of first (FD) delta.
! smom2  = second moment of second re-smearing delta.
! thdmom = third moment of total smearing delta given by occopt.
! tmom1  = third moment of first (FD) delta.
! tmom2  = third moment of second re-smearing delta.
! tratio  = ratio tsmear/tphysel for convoluted smearing function
! save values so we can impose recalculation of smdfun when
! the smearing or electronic temperature change between
! datasets
! corresponds roughly to delta_FD (maxFDarg) = 1.0d-100
! smd1    = temp array for smearing function 1
! smd2    = temp array for smearing function 2
!
! return fermi-dirac smearing function analytically
!
! real(dp) :: smdFD
! smdFD (tt) = 1.0_dp / (exp(-tt/2.0_dp) + exp(tt/2.0_dp))**2
!scalars
 integer,parameter :: nptsdiv2=6000
 integer,save :: dblsmr,occopt_prev=-9999
 integer :: algo,bantot,iband,iene,ii,ikpt,index,index_start,isppol,jj,nconvd2
 integer :: nene,nmaxFD,nminFD,prtdos
 real(dp),parameter :: maxFDarg=500.0_dp
 real(dp),save :: convlim,incconv,limit,tphysel_prev=-9999,tsmear_prev=-9999
 real(dp) :: aa,buffer,deltaene,dosdbletot,doshalftot,dostot,dsqrpi,encorr
 real(dp) :: enemax,enemin,enex,expinc,expx22,expxo2,factor,gauss,increm
 real(dp) :: intdostot,resFD,resFD1,resFD2,resFD3,resFD4,resmom,resmom1,resmom2
 real(dp) :: resmom3,resmom4,secmom,smom1,smom2,thdmom,tmom1,tmom2,tmpexpsum
 real(dp) :: tmpsmdfun,tratio,tsmearinv,tt,xx,yp1,ypn
 character(len=500) :: message
!arrays
 real(dp),save :: entfun(-nptsdiv2:nptsdiv2,2),occfun(-nptsdiv2:nptsdiv2,2)
 real(dp),save :: smdfun(-nptsdiv2:nptsdiv2,2),xgrid(-nptsdiv2:nptsdiv2)
 real(dp),allocatable :: arg(:),derfun(:),dos(:),dosdble(:),doshalf(:),ent(:)
 real(dp),allocatable :: entder(:),intdos(:),occder(:),smd(:),smd1(:),smd2(:)
 real(dp),allocatable :: smdder(:),tgrid(:),work(:),workfun(:)

! *************************************************************************

!DEBUG
!if(option==2)then
!write(6,*) ' getnel : enter '
!stop
!end if
!ENDDEBUG

 if(option/=1 .and. option/=2)then
  write(message, '(a,a,a,a,i6,a)' ) ch10,&
&  ' getnel : BUG -',ch10,&
&  '  Option must be either 1 or 2. It is ',option,'.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Initialize the occupation function and generalized entropy function,
!at the beginning, or if occopt changed
 if(occopt_prev/=occopt           .or. &
& abs(tsmear_prev-tsmear)  >tol12 .or. &
& abs(tphysel_prev-tphysel)>tol12       ) then

  occopt_prev=occopt
  tsmear_prev=tsmear
  tphysel_prev=tphysel

! 
! Check whether input values of tphysel tsmear and occopt are consistent
  dblsmr = 0
  if (abs(tphysel)>tol12) then
!  Use re-smearing scheme
   if (abs(tsmear)>tol12) then
    dblsmr = 1
!   Use FD occupations (one smearing) only with "physical" temperature tphysel
   else if (occopt /= 3) then
    write(message, '(a,a,a,a,i6,a)' ) ch10,&
&    ' getnel : ERROR -',ch10,&
&    '  tphysel /= 0, tsmear == 0, but occopt is not = 3, but ',occopt,'.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
  end if
! DEBUG
! write (6,*) 'getnel : input read.'
! write (6,*) '  dblsmr = ', dblsmr
! write (6,*) '  tphysel, tsmear = ', tphysel, tsmear
! ENDDEBUG

  allocate(entder(-nptsdiv2:nptsdiv2),occder(-nptsdiv2:nptsdiv2))
  allocate(smdder(-nptsdiv2:nptsdiv2) )
  allocate(workfun(-nptsdiv2:nptsdiv2),work(-nptsdiv2:nptsdiv2)  )

! Prepare the points on the grid
! limit is the value of the argument that will give 0.0 or 1.0 , with
! less than about 1.0d-15 error for 4<=occopt<=7, and less than about 1.0d-12
! error for occopt==3. It is not worth to compute the function beyond
! that point. Even with a less severe requirement, it is significantly
! larger for occopt==3, with an exponential
! tail, than for the other occupation functions, with a Gaussian tail.
! Note that these values are useful in newocc.f also.
  limit=6.0_dp
  if(occopt==3)limit=30.0_dp
  if(dblsmr /= 0) then
   tratio = tsmear / tphysel
   limit=30.0_dp + 6.0_dp*tratio
  end if

! With nptsdiv2=6000 (thus increm=0.001 for 4<=occopt<=7,
! and increm=0.005 for occopt==3, the O(1/N4) algorithm gives 1.0d-12
! accuracy on the stored values occfun and entfun. These, together
! with smdfun and xgrid, need permanently about 0.67 MB, which is affordable.
  increm=limit/nptsdiv2
  do ii=-nptsdiv2,nptsdiv2
   xgrid(ii)=ii*increm
  end do

! DEBUG
! write(6,*)' getnel : debug, occopt = ',occopt
! write(6,*)' limit,increm,nptsdiv2',limit,increm,nptsdiv2
! do ii=-nptsdiv2,nptsdiv2,nptsdiv2/30
! write(6, '(es10.2,3es22.14)' )xgrid(ii)
! end do
! ENDDEBUG

! ---------------------------------------------------------
! Ordinary (unique) smearing function
! ---------------------------------------------------------
  if (dblsmr == 0) then

!  Compute the unnormalized smeared delta function between -limit and +limit
!  (well, they are actually normalized ...)
   if(occopt==3)then

!   Fermi-Dirac
    do ii=0,nptsdiv2
     xx=xgrid(ii)
     smdfun( ii,1)=0.25_dp/(cosh(xx/2.0_dp)**2)
     smdfun(-ii,1)=smdfun(ii,1)
    end do

   else if(occopt==4 .or. occopt==5)then

!   Cold smearing of Marzari, two values of the "a" parameter being possible
!   first value gives minimization of the bump
    if(occopt==4)aa=-.5634
!   second value gives monotonic occupation function
    if(occopt==5)aa=-.8165
!   DEBUG
!   if(occopt==4)aa=-.2
!   if(occopt==5)aa=-.1
!   ENDDEBUG
    dsqrpi=1.0_dp/sqrt(pi)
    do ii=0,nptsdiv2
     xx=xgrid(ii)
     gauss=dsqrpi*exp(-xx**2)
     smdfun( ii,1)=gauss*(1.5_dp+xx*(-aa*1.5_dp+xx*(-1.0_dp+aa*xx)))
     smdfun(-ii,1)=gauss*(1.5_dp+xx*( aa*1.5_dp+xx*(-1.0_dp-aa*xx)))
    end do

   else if(occopt==6)then

!   First order Hermite-Gaussian of Paxton and Methfessel
    dsqrpi=1.0_dp/sqrt(pi)
    do ii=0,nptsdiv2
     xx=xgrid(ii)
     smdfun( ii,1)=dsqrpi*(1.5_dp-xx**2)*exp(-xx**2)
     smdfun(-ii,1)=smdfun(ii,1)
    end do

   else if(occopt==7)then

!   Gaussian smearing
    dsqrpi=1.0_dp/sqrt(pi)
    do ii=0,nptsdiv2
     xx=xgrid(ii)
     smdfun( ii,1)=dsqrpi*exp(-xx**2)
     smdfun(-ii,1)=smdfun(ii,1)
    end do

   else

    write(message, '(a,a,a,a,i4,a)' ) ch10,&
&    ' getnel : BUG -',ch10,&
&    '  Occopt=',occopt,' is not allowed in getnel. '
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

!  ---------------------------------------------------------
!  smear FD delta with occopt delta calculated in smdfun
!  ---------------------------------------------------------
  else if (dblsmr /= 0) then

   nconvd2 = 6000
   convlim = 10.0_dp
   incconv = convlim / nconvd2

!  DEBUG
!  write(6,*)' getnel : debug, convlim, = ', convlim
!  write(6,*)' getnel : debug, tratio = ', tratio
!  write(6,*)' estimated (temporary) memory usage: ', &
!  &    int(6.0*nconvd2*8.0/1024.0), ' kb '
!  ENDDEBUG

!  store smearing functions in smd1 and smd2
   allocate(smd1(-nconvd2:nconvd2) )
   allocate(smd2(-nconvd2:nconvd2) )
   allocate(tgrid(-nconvd2:nconvd2) )

!  FD function in smd1( ii) and second smearing delta in smd2( ii)
!  
!  smd1(:) contains delta_FD ( x )
   do ii=0,nconvd2
    tgrid(ii)=ii*incconv
    tgrid(-ii)=-tgrid(ii)
    tt=tgrid(ii)
    smd1( ii)=0.25_dp/(cosh(tt/2.0_dp)**2)
    smd1(-ii)=smd1(ii)
   end do

!  check input values of occopt and fill smd2(:) with appropriate data:
!  smd2(:) contains delta_resmear ( x )
   if(occopt == 3) then
    write(message, '(a,a,a,a,a)' ) ch10,&
&    ' getnel : ERROR -',ch10,&
&    '  Occopt=3 is not allowed as a re-smearing.', &
&    ' Use a single FD, or re-smear with a different delta type (faster cutoff). '
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   else if(occopt==4 .or. occopt==5)then
!   Cold smearing of Marzari, two values of the "a" parameter being possible
!   first value gives minimization of the bump
    if(occopt==4)aa=-.5634
!   second value gives monotonic occupation function
    if(occopt==5)aa=-.8165
!   DEBUG
!   if(occopt==4)aa=-.2
!   if(occopt==5)aa=-.1
!   ENDDEBUG
    dsqrpi=1.0_dp/sqrt(pi)
    do ii=0,nconvd2
     tt=tgrid(ii)
     gauss=dsqrpi*exp(-tt**2)
     smd2( ii)=gauss*(1.5_dp+tt*(-aa*1.5_dp+tt*(-1.0_dp+aa*tt)))
     smd2(-ii)=gauss*(1.5_dp+tt*( aa*1.5_dp+tt*(-1.0_dp-aa*tt)))
    end do
   else if(occopt==6)then
    dsqrpi=1.0_dp/sqrt(pi)
    do ii=0,nconvd2
     tt=tgrid(ii)
     smd2( ii)=dsqrpi*(1.5_dp-tt**2)*exp(-tt**2)
     smd2(-ii)=smd2(ii)
    end do
   else if(occopt==7)then
    dsqrpi=1.0_dp/sqrt(pi)
    do ii=0,nconvd2
     tt=tgrid(ii)
     smd2( ii)=dsqrpi*exp(-tt**2)
     smd2(-ii)=smd2(ii)
    end do
   else
    write(message, '(a,a,a,a,i4,a)' ) ch10,&
&    ' getnel : BUG -',ch10,&
&    '  Occopt=',occopt,' is not allowed in getnel. '
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if


!  Use O(1/N4) algorithm from Num Rec (see below)
!  
!  The grid for the convoluted delta is taken (conservatively)
!  to be that for the FD delta ie 6000 pts in [-limit;limit]
!  Smearing functions are given on [-dbllim;dbllim] and the grid must
!  superpose the normal grid on [-limit:limit]
!  The maximal interval for integration of the convolution is
!  [-dbllim+limit+lim(delta2);dbllim-limit-lim(delta2)] =
!  [-dbllim+36;dbllim-36]

!  test the smdFD function for extreme values:
!  do jj=-nptsdiv2,-nptsdiv2
!  do ii=-nconvd2+4,nconvd2
!  call smdFD(xgrid(jj) - tgrid(ii)*tratio, resFD)
!  write (6,*) 'ii jj = ', ii,jj, ' smdFD (', &
!  &    xgrid(jj) - tgrid(ii)*tratio, ') ', resFD
!  end do
!  end do

   expinc = exp(half*incconv*tratio)

!  
!  jj = position of point at which we are calculating smdfun
!  
   do jj=-nptsdiv2,nptsdiv2
!   Do not care about the 8 boundary points,
!   where the values should be extremely small anyway
    smdfun(jj,1)=0.0_dp
!   only add contribution with delta_FD > 1.0d-100
    nmaxFD = floor  (( maxFDarg+xgrid(jj)) / tratio / incconv )
    nmaxFD = min (nmaxFD, nconvd2)
    nminFD = ceiling((-maxFDarg+xgrid(jj)) / tratio / incconv )
    nminFD = max (nminFD, -nconvd2)

!   DEBUG
!   write (6,*) ' getnel : nmaxFD, nminFD = ', nmaxFD, nminFD
!   ENDDEBUG

!   Calculate the Fermi-Dirac distrib at point xgrid(jj)-tgrid(ii)*tratio
    expxo2 = exp (-half*(xgrid(jj) - (nminFD)*incconv*tratio))
    expx22 = expxo2*expxo2
    tmpexpsum = expxo2 / (expx22 + 1.0_dp)
    resFD4 = tmpexpsum * tmpexpsum
    expxo2 = expxo2*expinc
    expx22 = expxo2*expxo2
    tmpexpsum = expxo2 / (expx22 + 1.0_dp)
    resFD3 = tmpexpsum * tmpexpsum
    expxo2 = expxo2*expinc
    expx22 = expxo2*expxo2
    tmpexpsum = expxo2 / (expx22 + 1.0_dp)
    resFD2 = tmpexpsum * tmpexpsum
    expxo2 = expxo2*expinc
    expx22 = expxo2*expxo2
    tmpexpsum = expxo2 / (expx22 + 1.0_dp)
    resFD1 = tmpexpsum * tmpexpsum

!   DEBUG
!   write (6,*) ' getnel : expinc, expxo2 = ', expinc, expxo2
!   ENDDEBUG

!   
!   core contribution to the integral with constant weight (48)
!   
    tmpsmdfun = 0.0_dp
    do ii=nminFD+4,nmaxFD-4
     expxo2 = expxo2*expinc
!    tmpexpsum = 1.0_dp / (expxo2 + 1.0_dp / expxo2 )
     expx22 = expxo2*expxo2
     tmpexpsum = expxo2 / (expx22 + 1.0_dp)
     tmpsmdfun = tmpsmdfun + smd2(ii) * tmpexpsum * tmpexpsum
    end do
!   
!   Add on end contributions for show (both functions smd and smdFD are
!   very small
    smdfun(jj,1)=smdfun(jj,1)       +48.0_dp*tmpsmdfun             &
&    + 31.0_dp*smd2(nminFD+3)*resFD1 -11.0_dp*smd2(nminFD+2)*resFD2 &
&    +  5.0_dp*smd2(nminFD+1)*resFD3 -       smd2(nminFD)*resFD4

    expxo2 = expxo2*expinc
    expx22 = expxo2*expxo2
    tmpexpsum = expxo2 / (expx22 + 1.0_dp)
    resFD1 = tmpexpsum * tmpexpsum
    expxo2 = expxo2*expinc
    expx22 = expxo2*expxo2
    tmpexpsum = expxo2 / (expx22 + 1.0_dp)
    resFD2 = tmpexpsum * tmpexpsum
    expxo2 = expxo2*expinc
    expx22 = expxo2*expxo2
    tmpexpsum = expxo2 / (expx22 + 1.0_dp)
    resFD3 = tmpexpsum * tmpexpsum
    expxo2 = expxo2*expinc
    expx22 = expxo2*expxo2
    tmpexpsum = expxo2 / (expx22 + 1.0_dp)
    resFD4 = tmpexpsum * tmpexpsum

!   Contribution above
    smdfun(jj,1)=smdfun(jj,1)                                      &
&    + 31.0_dp*smd2(nmaxFD-3)*resFD1  -11.0_dp*smd2(nmaxFD-2)*resFD2 &
&    +  5.0_dp*smd2(nmaxFD-1)*resFD3  -       smd2(nmaxFD)*resFD4
    smdfun(jj,1)=incconv*smdfun(jj,1)/48.0_dp

!   DEBUG
!   write (6,*) ' getnel : smdfun(',jj,')= ', smdfun(jj,1)
!   ENDDEBUG
   end do

!  write (6,*) 'getnel : convolution done ! '
!  secmom = second moment of delta in smdfun(:,1)
!  = $\int_{-\infty}^{+\infty} t^2 delta_2 (t) dt$
!  
   secmom = 0.0_dp
   thdmom = 0.0_dp
   resmom4 = xgrid(-nptsdiv2  )*xgrid(-nptsdiv2  )*smdfun(-nptsdiv2  ,  1)
   resmom3 = xgrid(-nptsdiv2+1)*xgrid(-nptsdiv2+1)*smdfun(-nptsdiv2+1,  1)
   resmom2 = xgrid(-nptsdiv2+2)*xgrid(-nptsdiv2+2)*smdfun(-nptsdiv2+2,  1)
   resmom1 = xgrid(-nptsdiv2+3)*xgrid(-nptsdiv2+3)*smdfun(-nptsdiv2+3,  1)
   resmom  = xgrid(-nptsdiv2+4)*xgrid(-nptsdiv2+4)*smdfun(-nptsdiv2+4,  1)
   do ii=-nptsdiv2,nptsdiv2
    secmom = secmom +                                   &
&    ( 17.0_dp*xgrid(ii)  *xgrid(ii)  *smdfun(ii,  1)   &
&    +42.0_dp*xgrid(ii-1)*xgrid(ii-1)*smdfun(ii-1,1)   &
&    -16.0_dp*xgrid(ii-2)*xgrid(ii-2)*smdfun(ii-2,1)   &
&    + 6.0_dp*xgrid(ii-3)*xgrid(ii-3)*smdfun(ii-3,1)   &
&    -       xgrid(ii-4)*xgrid(ii-4)*smdfun(ii-4,1)  )
    resmom4 = resmom3
    resmom3 = resmom2
    resmom2 = resmom1
    resmom1 = resmom
    resmom  = xgrid(ii+1)  *xgrid(ii+1)  *smdfun(ii+1,  1)
   end do
   secmom=increm * secmom / 48.0_dp
!  thdmom=increm * thdmom / 48.0_dp
!  
!  smom1  = second moment of delta in smd1(:)
!  smom2  = second moment of delta in smd2(:)
!  
   smom1  = 0.0_dp
   smom2  = 0.0_dp
   tmom1  = 0.0_dp
   tmom2  = 0.0_dp
   do ii=-nconvd2+4,nconvd2
    smom1 = smom1+                                       &
&    ( 17.0_dp*tgrid(ii)  *tgrid(ii)  *smd1(ii)         &
&    +42.0_dp*tgrid(ii-1)*tgrid(ii-1)*smd1(ii-1)       &
&    -16.0_dp*tgrid(ii-2)*tgrid(ii-2)*smd1(ii-2)       &
&    + 6.0_dp*tgrid(ii-3)*tgrid(ii-3)*smd1(ii-3)       &
&    -       tgrid(ii-4)*tgrid(ii-4)*smd1(ii-4)  )
    smom2 = smom2+                                       &
&    ( 17.0_dp*tgrid(ii)  *tgrid(ii)  *smd2(ii  )     &
&    +42.0_dp*tgrid(ii-1)*tgrid(ii-1)*smd2(ii-1)     &
&    -16.0_dp*tgrid(ii-2)*tgrid(ii-2)*smd2(ii-2)     &
&    + 6.0_dp*tgrid(ii-3)*tgrid(ii-3)*smd2(ii-3)     &
&    -       tgrid(ii-4)*tgrid(ii-4)*smd2(ii-4)  )
   end do
   smom1 =incconv * smom1  / 48.0_dp
   smom2 =incconv * smom2  / 48.0_dp
!  tmom1 =incconv * tmom1  / 48.0_dp
!  tmom2 =incconv * tmom2  / 48.0_dp

   encorr =  smom2*tratio*tratio/secmom

!  DEBUG
!  write (6,*) ' getnel : debug, secmoms = ', secmom, smom1, smom2
!  write (6,*) ' getnel : debug, thdmoms = ', thdmom, tmom1, tmom2
!  write (6,*) ' getnel : encorr = ', encorr
!  ENDDEBUG

   deallocate (tgrid, smd1, smd2)

  end if

! --------------------------------------------------------
! end of smearing function initialisation, dblsmr case
! --------------------------------------------------------

! DEBUG
! write(6,*)' getnel : debug, occopt = ',occopt
! do ii=-nptsdiv2,nptsdiv2,nptsdiv2/30
! write(6, '(es10.2,3es22.14)' )xgrid(ii)
! end do
! ENDDEBUG

! Now that the smeared delta function has been initialized, compute the
! occupation function
  occfun(-nptsdiv2,1)=zero
  entfun(-nptsdiv2,1)=zero

! Different algorithms are possible, corresponding to the formulas
! (4.1.11), (4.1.12) and (4.1.14) in Numerical recipes (pp 107 and 108),
! with respective O(1/N2), O(1/N3), O(1/N4) convergence, where N is the
! number of points in the interval.
  algo=4

  if(algo==2)then

!  Extended trapezoidal rule (4.1.11), taken in a cumulative way
   do ii=-nptsdiv2+1,nptsdiv2
    occfun(ii,1)=occfun(ii-1,1)+increm*(smdfun(ii,1)+smdfun(ii-1,1))/2.0_dp
    entfun(ii,1)=entfun(ii-1,1)+increm*&
&    ( -xgrid(ii)*smdfun(ii,1) -xgrid(ii-1)*smdfun(ii-1,1) )/2.0_dp
   end do

  else if(algo==3)then

!  Derived from (4.1.12). Converges as O(1/N3).
!  Do not care about the following points,
!  where the values are extremely small anyway
   occfun(-nptsdiv2+1,1)=0.0_dp ;   entfun(-nptsdiv2+1,1)=0.0_dp
   do ii=-nptsdiv2+2,nptsdiv2
    occfun(ii,1)=occfun(ii-1,1)+increm*&
&    ( 5.0_dp*smdfun(ii,1) + 8.0_dp*smdfun(ii-1,1) - smdfun(ii-2,1) )/12.0_dp
    entfun(ii,1)=entfun(ii-1,1)+increm*&
&    ( 5.0_dp*(-xgrid(ii)  )*smdfun(ii,1)  &
&    +8.0_dp*(-xgrid(ii-1))*smdfun(ii-1,1)&
&    -      (-xgrid(ii-2))*smdfun(ii-2,1) )/12.0_dp
   end do

  else if(algo==4)then

!  Derived from (4.1.14)- alternative extended Simpsons rule. Converges as O(1/N4).
!  Do not care about the following points,
!  where the values are extremely small anyway
   occfun(-nptsdiv2+1,1)=0.0_dp ;   entfun(-nptsdiv2+1,1)=0.0_dp
   occfun(-nptsdiv2+2,1)=0.0_dp ;   entfun(-nptsdiv2+2,1)=0.0_dp
   occfun(-nptsdiv2+3,1)=0.0_dp ;   entfun(-nptsdiv2+3,1)=0.0_dp
   do ii=-nptsdiv2+4,nptsdiv2
    occfun(ii,1)=occfun(ii-1,1)+increm*&
&    ( 17.0_dp*smdfun(ii,1)  &
&    +42.0_dp*smdfun(ii-1,1)&
&    -16.0_dp*smdfun(ii-2,1)&
&    + 6.0_dp*smdfun(ii-3,1)&
&    -       smdfun(ii-4,1) )/48.0_dp
    entfun(ii,1)=entfun(ii-1,1)+increm*&
&    ( 17.0_dp*(-xgrid(ii)  )*smdfun(ii,1)  &
&    +42.0_dp*(-xgrid(ii-1))*smdfun(ii-1,1)&
&    -16.0_dp*(-xgrid(ii-2))*smdfun(ii-2,1)&
&    + 6.0_dp*(-xgrid(ii-3))*smdfun(ii-3,1)&
&    -       (-xgrid(ii-4))*smdfun(ii-4,1) )/48.0_dp
   end do

!  End of choice between different algorithms for integration
  end if


! Normalize the functions (actually not needed for occopt=3..7)
! DEBUG
! write(message, '(a,es16.8,a,a,a)' ) &
! &  ' getnel : integral of the function is ',occfun(nptsdiv2,1),'.',ch10,&
! &  '    performs now normalization '
! call wrtout(6,message,'COLL')
! ENDDEBUG
  factor=1.0_dp/occfun(nptsdiv2,1)
  smdfun(:,1)=smdfun(:,1)*factor
  occfun(:,1)=occfun(:,1)*factor
  entfun(:,1)=entfun(:,1)*factor

! Compute the cubic spline fitting of the smeared delta function
  yp1=0.0_dp ; ypn=0.0_dp
  workfun(:)=smdfun(:,1)
  call spline(xgrid, workfun, (2*nptsdiv2+1), yp1, ypn, smdder, work)
  smdfun(:,2)=smdder(:)

! Compute the cubic spline fitting of the occupation function
  yp1=0.0_dp ; ypn=0.0_dp
  workfun(:)=occfun(:,1)
  call spline(xgrid, workfun, (2*nptsdiv2+1), yp1, ypn, occder, work)
  occfun(:,2)=occder(:)

! Compute the cubic spline fitting of the entropy function
  yp1=0.0_dp ; ypn=0.0_dp
  workfun(:)=entfun(:,1)
  call spline(xgrid, workfun, (2*nptsdiv2+1), yp1, ypn, entder, work)
  entfun(:,2)=entder(:)

! DEBUG
! write(6,*)' getnel : debug '
! open(unit=102,file="smrfcts.dat",form='formatted',status='unknown')
! rewind(unit=102)
! do ii=-nptsdiv2,nptsdiv2,nptsdiv2/30
! write(102, '(es10.2,3es22.14)' )xgrid(ii),smdfun(ii,1),occfun(ii,1),entfun(ii,1)
! end do
! close (102)
! write(6,*)' getnel : factor = ', factor
! stop
! ENDDEBUG

  deallocate(entder,occder,smdder,work,workfun)

! The initialisation of occfun and entfun is done
 end if

!---------------------------------------------------------------------

!DEBUG
!write(6,*)' getnel : debug  tphysel, tsmear = ', tphysel, tsmear
!ENDDEBUG
 bantot=sum(nband(:))
 if (abs(tphysel)<tol12) then
  tsmearinv=one/tsmear
 else
  tsmearinv=one/tphysel
 end if

 allocate(arg(bantot),derfun(bantot),ent(bantot),smd(bantot))

!
!normal evaluation of occupations and entropy
!
 if(option==1)then

! Compute the arguments of the occupation and entropy functions
  arg(:)=(fermie-eigen(1:bantot))*tsmearinv

! Compute the values of the occupation function, and the entropy function
! Note : splfit also takes care of the points outside of the interval,
! and assign to them the value of the closest extremal point,
! which is what is needed here.

  call splfit(xgrid,doccde,occfun,1,arg,occ,(2*nptsdiv2+1),bantot)
  call splfit(xgrid,derfun,entfun,0,arg,ent,(2*nptsdiv2+1),bantot)

! Normalize occ and ent, and sum number of electrons and entropy
  nelect=zero
  entropy=zero
  index=0
  do isppol=1,nsppol
   do ikpt=1,nkpt
    do iband=1,nband(ikpt+nkpt*(isppol-1))
     index=index+1
     ent(index)=ent(index)*maxocc
     occ(index)=occ(index)*maxocc
     doccde(index)=-doccde(index)*maxocc*tsmearinv
     entropy=entropy+wtk(ikpt)*ent(index)
     nelect=nelect+wtk(ikpt)*occ(index)
    end do
   end do
  end do

! DEBUG
! write(6,*) ' getnel : debug   wtk, occ, eigen = ', wtk, occ, eigen
! END DEBUG

! DEBUG
! write(6,*)xgrid(-nptsdiv2),xgrid(nptsdiv2)
! write(6,*)'fermie',fermie
! do ii=1,bantot
! write(6,*)ii,arg(ii),doccde(ii)
! end do
! write(6,*)'eigen',eigen(:)
! write(6,*)'arg',arg(:)
! write(6,*)'occ',occ(:)
! write(6,*)'nelect',nelect
! stop
! ENDDEBUG

! 
! evaluate DOS for smearing, half smearing, and double.
! 
 else if(option==2)then

  buffer=limit/tsmearinv*.5_dp
  prtdos=1

! Write the header of the DOS file, and also decides the energy range and
! increment
  call dos_hdr_write(buffer,deltaene,dosdeltae,eigen,enemax,enemin,fermie,mband,nband,nene,&
&  nkpt,nsppol,occopt,prtdos,tphysel,tsmear,unitdos)

  allocate(dos(bantot),dosdble(bantot),doshalf(bantot),intdos(bantot) )

  do isppol=1,nsppol

   if (nsppol==2) then
    if(isppol==1) write(message,'(a,16x,a)')  '#','Spin-up DOS'
    if(isppol==2) write(message,'(a,16x,a)')  '#','Spin-dn DOS '
    call wrtout(unitdos,message,'COLL')
   end if
   index_start=0
   if(isppol==2)then
    do ikpt=1,nkpt
     index_start=index_start+nband(ikpt)
    end do
   end if

   enex=enemin
   do iene=1,nene

!   Compute the arguments of the dos and occupation function
    arg(:)=(enex-eigen(1:bantot))*tsmearinv

    call splfit(xgrid,derfun,smdfun,0,arg,dos,(2*nptsdiv2+1),bantot)
    call splfit(xgrid,derfun,occfun,0,arg,intdos,(2*nptsdiv2+1),bantot)
!   Also compute the dos with tsmear halved and doubled
    arg(:)=arg(:)*2.0_dp
    call splfit(xgrid,derfun,smdfun,0,arg,doshalf,(2*nptsdiv2+1),bantot)
!   Since arg was already doubled, must divide by four
    arg(:)=arg(:)*0.25_dp
    call splfit(xgrid,derfun,smdfun,0,arg,dosdble,(2*nptsdiv2+1),bantot)

!   Now, accumulate the contribution from each eigenenergy
    dostot=0.0_dp
    intdostot=0.0_dp
    doshalftot=0.0_dp
    dosdbletot=0.0_dp
    index=index_start

!   DEBUG
!   write(6,*)' eigen, arg, dos, intdos, doshalf, dosdble'
!   ENDDEBUG
    do ikpt=1,nkpt
     do iband=1,nband(ikpt+nkpt*(isppol-1))
      index=index+1
      dostot=dostot+wtk(ikpt)*maxocc*dos(index)*tsmearinv
      intdostot=intdostot+wtk(ikpt)*maxocc*intdos(index)
      doshalftot=doshalftot+wtk(ikpt)*maxocc*doshalf(index)*tsmearinv*2.0_dp
      dosdbletot=dosdbletot+wtk(ikpt)*maxocc*dosdble(index)*tsmearinv*0.5_dp
     end do
    end do

!   Print the data for this energy
    write(message, '(i5,f8.3,2f14.6,2f14.3)' )&
&    iene-1,enex,dostot,intdostot,doshalftot,dosdbletot
    call wrtout(unitdos,message,'COLL')

    enex=enex+deltaene

!   End the loop over the energies
   end do

!  End the loop over isppol
  end do

  deallocate(dos,dosdble,doshalf,intdos)

! Close the DOS file
  close(unitdos)

 end if

 deallocate(arg,derfun,ent,smd)

!DEBUG
!if(option==2)then
!write(6,*) ' getnel : exit '
!stop
!end if
!ENDDEBUG

end subroutine getnel
!!***
