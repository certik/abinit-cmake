!{\src2tex{textfont=tt}}
!!****f* ABINIT/rhohxc
!! NAME
!! rhohxc
!!
!! FUNCTION
!! Start from the density or spin-density, and
!! compute Hartree (if option>=1) and xc correlation potential and energies.
!! Eventually compute xc kernel (if option=-2, 2, 3, 10 or 12).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, MF, GZ, DRH, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | intxc=0 for old quadrature; 1 for new improved quadrature
!!   | ixc= choice of exchange-correlation scheme (see above, and below)
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!! (gsqcut=(boxcut**2)*ecut/(2.d0*(Pi**2))
!!  izero=if 1, unbalanced components of Vhartree(g) have to be set to zero
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat(nfftf,nspden*nhatdim)= -PAW only- compensation density
!!  nhatdim= -PAW only- 0 if nhat array is not used ; 1 otherwise
!!  nhatgr(nfftf,nspden,3*nhatgrdim)= -PAW only- cartesian gradients of compensation density
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  nkxc=second dimension of the kxc array. If /=0,
!!   the exchange-correlation kernel must be computed.
!!  nspden=number of spin-density components
!!  n3xccc=dimension of the xccc3d array (0 or nfft or cplx*nfft).
!!  option=0 for xc only (exc, vxc, strsxc),
!!         1 for Hxc (idem + vhartr) ,
!!         2 for Hxc and kxc (no paramagnetic part if nspden=1)
!!        10 for xc  and kxc with only LDA part (d2Exc/drho^2)
!!        12 for Hxc and kxc with only LDA part (d2Exc/drho^2)
!!         3 for Hxc, kxc and k3xc
!!        -1 return (no xc terms for positron)
!!        -2 for Hxc and kxc (with paramagnetic part if nspden=1)
!!  rhog(2,nfft)=electron density in G space
!!  rhor(nfft,nspden)=electron density in real space in electrons/bohr**3
!!   (total in first half and spin-up in second half if nspden=2)
!!   (total in first comp. and magnetization in comp. 2 to 4 if nspden=4)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!
!! OUTPUT
!!  enxc=returned exchange and correlation energy (hartree).
!!
!!  === Only if abs(option)=3 ===
!!  k3xc(nfft)=third derivative of the XC energy functional of the density,
!!    at each point of the real space grid (only in the LDA,
!!    non-spin-polarized)
!!
!!  === Only if abs(option)=2, -2, 3, 10, 12 ===
!!  kxc(nfft,nkxc)=exchange and correlation kernel
!!                 (returned only if nkxc/=0)
!!   allowed if LDAs (dtset%xclevel=1 or option=(10 or 12)) :
!!    if nspden==1: return kxc(:,1)= d2Exc/drho2
!!       that is 1/2 ( d2Exc/drho_up drho_up + d2Exc/drho_up drho_dn )
!!    if nspden==1: also return kxc(:,2)= d2Exc/drho_up drho_dn
!!    if nspden>=2, return  kxc(:,1)=d2Exc/drho_up drho_up
!!                          kxc(:,2)=d2Exc/drho_up drho_dn
!!                          kxc(:,3)=d2Exc/drho_dn drho_dn
!!   allowed also if GGAs (dtset%xclevel=2 and option=2 or 3)
!!    for the time being, treat all cases as spin-polarized, with nkxc=23
!!    kxc(:,1)= d2Ex/drho_up drho_up
!!    kxc(:,2)= d2Ex/drho_dn drho_dn
!!    kxc(:,3)= dEx/d(abs(grad(rho_up))) / abs(grad(rho_up))
!!    kxc(:,4)= dEx/d(abs(grad(rho_dn))) / abs(grad(rho_dn))
!!    kxc(:,5)= d2Ex/d(abs(grad(rho_up))) drho_up / abs(grad(rho_up))
!!    kxc(:,6)= d2Ex/d(abs(grad(rho_dn))) drho_dn / abs(grad(rho_dn))
!!    kxc(:,7)= 1/abs(grad(rho_up)) * d/drho_up (dEx/d(abs(grad(rho_up))) /abs(grad(rho_up)))
!!    kxc(:,8)= 1/abs(grad(rho_dn)) * d/drho_dn (dEx/d(abs(grad(rho_dn))) /abs(grad(rho_dn)))
!!    kxc(:,9)= d2Ec/drho_up drho_up
!!    kxc(:,10)=d2Ec/drho_up drho_dn
!!    kxc(:,11)=d2Ec/drho_dn drho_dn
!!    kxc(:,12)=dEc/d(abs(grad(rho))) / abs(grad(rho))
!!    kxc(:,13)=d2Ec/d(abs(grad(rho))) drho_up / abs(grad(rho))
!!    kxc(:,14)=d2Ec/d(abs(grad(rho))) drho_dn / abs(grad(rho))
!!    kxc(:,15)=1/abs(grad(rho)) * d/drho (dEc/d(abs(grad(rho))) /abs(grad(rho)))
!!    kxc(:,16)=rho_up
!!    kxc(:,17)=rho_dn
!!    kxc(:,18)=gradx(rho_up)
!!    kxc(:,19)=gradx(rho_dn)
!!    kxc(:,20)=grady(rho_up)
!!    kxc(:,21)=grady(rho_dn)
!!    kxc(:,22)=gradz(rho_up)
!!    kxc(:,23)=gradz(rho_dn)
!!
!!  strsxc(6)= contribution of xc to stress tensor (hartree/bohr^3),
!!   given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).
!!   Explicitely : strsxc(mu,nu) = (1/N) Sum(i=1,N)
!!    ( delta(mu,nu) * [  exc(i)rhotot(i)
!!               - drhoexc_drho(up,i)*rhor(up,i)-drhoexc_drho(dn,i)*rhor(dn,i)]
!!     - gradrho(up,mu)*gradrho(up,nu) * drhoexc_dgradrho(up,i) / gradrho(up,i)
!!     - gradrho(dn,mu)*gradrho(dn,nu) * drhoexc_dgradrho(dn,i) / gradrho(dn,i) )
!!  vhartr(nfft)=Hartree potential (returned if option/=0 and option/=10)
!!  vxc(nfft,nspden)=xc potential
!!    (spin up in first half and spin down in second half if nspden=2)
!!    (v^11, v^22, Re[V^12], Im[V^12] if nspden=4)
!!  vxcavg=<Vxc>=unit cell average of Vxc = (1/ucvol) Int [Vxc(r) d^3 r].
!!
!! NOTES
!! Start from the density, and
!! compute Hartree (if option>=1) and xc correlation potential and energies.
!! Eventually compute xc kernel (if option=-2, 2, 3, 10 or 12).
!! Allows a variety of exchange-correlation functionals
!! according to ixc. Here is a list of allowed values.
!!                                                    subroutine name
!!    0 means no xc applied (usually for testing)
!! *LDA,LSD
!!    1 means new Teter (4/93) with spin-pol option        xcspol
!!    2 means Perdew-Zunger-Ceperley-Alder                 xcpzca
!!    3 means old Teter (4/91) fit to Ceperley-Alder data  xctetr
!!    4 means Wigner                                       xcwign
!!    5 means Hedin-Lundqvist                              xchelu
!!    6 means "X-alpha" xc                                 xcxalp
!!    7 mean Perdew-Wang 92 LSD fit to Ceperley-Alder data xcpbe
!!    8 mean Perdew-Wang 92 LSD , exchange-only            xcpbe
!!    9 mean Perdew-Wang 92 Ex+Ec_RPA  energy              xcpbe
!!   10 means RPA LSD energy (only the energy !!)          xcpbe
!! *GGA
!!   11 means Perdew-Burke-Ernzerhof GGA functional        xcpbe
!!   12 means x-only Perdew-Burke-Ernzerhof GGA functional xcpbe
!!   13 means LDA (ixc==7), except that the xc potential
!!      is given within the van Leeuwen-Baerends GGA       xclb
!!   14 means revPBE GGA functional                        xcpbe
!!   15 means RPBE GGA functional                          xcpbe
!!   16 means HCTH GGA functional                          xchcth
!!   23 means WC GGA functional                            xcpbe
!! *Fermi-Amaldi
!!   20 means Fermi-Amaldi correction
!!   21 means Fermi-Amaldi correction with LDA(ixc=1) kernel
!!   22 means Fermi-Amaldi correction with hybrid BPG kernel
!! *Nanoquanta libxc routines
!!   30 lda_x.c
!!   31 lda_x.c + lda_c_vwn.c
!!   32 lda_x.c + lda_c_pz.c
!!   33 lda_x.c + lda_c_pw.c
!!   34 lda_x.c + lda_c_amgb.c
!!
!! Allow for improved xc quadrature (intxc=1) by using the usual FFT grid
!! as well as another, shifted, grid, and combining both results.
!! Spin-polarization is allowed only with ixc=0, 1, and GGAs until now.
!! Note : this routine has been optimized already. See the end of this routine.
!!
!! PARENTS
!!      rhohxc
!!
!! CHILDREN
!!      dotprod_vn,drivexc,hartre,leave_new,mean_fftr,metric,mkdenpos,timab
!!      wrtout,xcden,xcmult,xcomm_init,xcpot,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine rhohxc(dtset,enxc,gsqcut,izero,kxc,mpi_enreg,nfft,ngfft,&
& nhat,nhatdim,nhatgr,nhatgrdim,nkxc,nspden,n3xccc,option,rhog,rhor,rprimd, &
& strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d,k3xc)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_12spacepar
 use interfaces_13xc, except_this_one => rhohxc
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: izero,n3xccc,nfft,nhatdim,nhatgrdim,nkxc,nspden,option
 integer,intent(in) :: usexcnhat
 real(dp),intent(in) :: gsqcut
 real(dp),intent(out) :: enxc,vxcavg
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: nhat(nfft,nspden*nhatdim)
 real(dp),intent(in) :: nhatgr(nfft,nspden,3*nhatgrdim),rhog(2,nfft)
 real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3),xccc3d(n3xccc)
 real(dp),intent(out) :: kxc(nfft,nkxc),strsxc(6),vhartr(nfft),vxc(nfft,nspden)
 real(dp),intent(out),optional :: k3xc(1:nfft)

!Local variables-------------------------------
!scalars
 integer :: cplex,i1,i2,i3,ierr,ifft,ii,index,inkxc,ipts,ir,ishift,ispden,iwarn
 integer :: ixc,jj,mpts,n1,n2,n3,ndvxc,nfftot,ngr2,ngrad,npts,nspden_eff
 integer :: nspden_updn,nspgrad,nvxcdgr,old_paral_level,optpbe,order,spaceComm
 real(dp),parameter :: mot=-one/3.0_dp
 real(dp) :: coeff,divshft,doti,dstrsxc,dvdn,dvdz,factor,nelect,rhotmp,s1,s2,s3
 real(dp) :: strdiag,strsxc1_tot,strsxc2_tot,strsxc3_tot,strsxc4_tot
 real(dp) :: strsxc5_tot,strsxc6_tot,ucvol
 logical :: test_nhat
 character(len=500) :: message
!arrays
 real(dp) :: dmnorm(3),drho_updn(3),gmet(3,3),gprimd(3,3),qphon(3),rmet(3,3)
 real(dp) :: tsec(2),vxcmean(4)
 real(dp),allocatable :: d2vxcar(:),dnexcdn(:,:),dvxcdgr(:,:),dvxci(:,:)
 real(dp),allocatable :: exci(:),grho2_updn(:,:),m_norm(:),nhat_up(:)
 real(dp),allocatable :: rho_updn(:,:),rhoarr(:),rhocorval(:,:),rhonow(:,:,:)
 real(dp),allocatable :: rspts(:),vxci(:,:),zeta(:)
!no_abirules
#if defined OPENMP
           integer,external :: OMP_GET_NUM_THREADS
#endif

! *************************************************************************

!DEBUG
!write(6,*)' rhohxc : enter with option, nspden ',option,nspden
!stop
!ENDDEBUG

!Check options
 if(nspden/=1 .and. nspden/=2 .and. nspden/=4)then
  write(message, '(a,a,a,a,a,a,i5)' ) ch10,&
&  ' rhohxc :  BUG -',ch10,&
&  '  The only allowed values of nspden are 1, 2, or 4,',ch10,&
&  '  while the argument nspden=',nspden
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if ((option==3).and.(dtset%ixc/=7).and.(dtset%ixc/=3)) then
  write(message, '(a,a,a,a,a,a,i5)' ) ch10,&
&  ' rhohxc :  ERROR -',ch10,&
&  '  Third-order xc kernel can only be computed for ixc = 3 or ixc =7,',ch10,&
&  '  while it is found to be',dtset%ixc
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if(dtset%icoulomb /= 0)then
  write(message, '(a,a,a,a,a,a,i5)' ) ch10,&
&  ' rhohxc : BUG -',ch10,&
&  '  To use non-periodic computation (icoulomb /= 0), ',ch10,&
&  '  use PSolver_rhohxc() instead, while the argument icoulomb=',dtset%icoulomb
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if(nspden==4.and.dtset%xclevel==2.and.(abs(option)==2))then
  write(message, '(4a)' ) ch10,&
&  ' rhohxc :  BUG -',ch10,&
&  '  When nspden==4 and GGA, the absolute value of option cannot be 2 !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Nothing to do if option==-1
 if (option==-1) then
  enxc=zero
  vxc(:,:)=zero
  vxcavg=zero
  strsxc(:)=zero
  return
 end if

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!In this routine, hartre, xcden and xcpot are called for real
!densities and potentials, corresponding to zero wavevector
 cplex=1
 qphon(:)=zero
 iwarn=0
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)

 if(option/=0.and.option/=10)then
  call hartre(cplex,gmet,gsqcut,izero,mpi_enreg,nfft,ngfft,dtset%paral_kgb,qphon,rhog,vhartr)
 end if

!Note : hartre is excluded from the timing
 call timab(81,1,tsec)

 enxc=zero
 vxc(:,:)=zero
 strsxc(:)=zero
 strsxc1_tot=zero
 strsxc2_tot=zero
 strsxc3_tot=zero
 strsxc4_tot=zero
 strsxc5_tot=zero
 strsxc6_tot=zero

 if (dtset%ixc==0) then

  vxcavg=zero
  if (nkxc/=0) kxc(:,:)=zero
! No xc at all is applied (usually for testing)
  write(message, '(a,a,a,a)' ) ch10,&
&  ' rhohxc : WARNING -',ch10,&
&  '  Note that no xc is applied (ixc=0).'
  call wrtout(06,message,'COLL')

 else if (dtset%ixc/=20) then

! ngrad=1 is for LDAs or LSDs, ngrad=2 is for GGAs
  ngrad=1;if(dtset%xclevel==2)ngrad=2
! Test: has a compensation density to be added/substracted (PAW) ?
  test_nhat=((nhatdim==1).and.(usexcnhat==0.or.(ngrad==2.and.nhatgrdim==1)))
! nspden_updn: 1 for non-polarized, 2 for polarized
  nspden_updn=min(nspden,2)
! nspden_eff: effective value of nspden used to compute gradients of density:
! 1 for non-polarized system,
! 2 for collinear polarized system or LDA (can be reduced to a collinear system)
! 4 for non-collinear polarized system and GGA
  nspden_eff=nspden_updn;if (nspden==4.and.ngrad==2) nspden_eff=4

! The different components of dnexcdn will be
! for nspden=1,   dnexcdn(:,1)=d(n.exc)/d(n)
! and if ngrad=2, dnexcdn(:,2)=1/2*1/|grad n_up|*d(n.exc)/d(|grad n_up|)
! +   1/|grad n|*d(n.exc)/d(|grad n|)
! (do not forget : |grad n| /= |grad n_up| + |grad n_down|
! for nspden>=2,          dnexcdn(:,1)=d(n.exc)/d(n_up)
! dnexcdn(:,2)=d(n.exc)/d(n_down)
! and if ngrad=2, dnexcdn(:,3)=1/|grad n_up|*d(n.exc)/d(|grad n_up|)
! dnexcdn(:,4)=1/|grad n_down|*d(n.exc)/d(|grad n_down|)
! dnexcdn(:,5)=1/|grad n|*d(n.exc)/d(|grad n|)
! Note: if nspden=4, n_up=(n+|m|)/2, n_down=(n-|m|)/2
  nspgrad=nspden_updn*ngrad;if(nspden_updn==2.and.ngrad==2)nspgrad=5
  allocate(dnexcdn(nfft,nspgrad))

! Non-collinear magnetism: store norm of magnetization
  if (nspden==4) then
   allocate(m_norm(nfft));m_norm=zero
   if ((usexcnhat==1).or.(nhatdim==0)) then
    m_norm(:)=sqrt(rhor(:,2)**2+rhor(:,3)**2+rhor(:,4)**2)
   else
    m_norm(:)=sqrt((rhor(:,2)-nhat(:,2))**2+(rhor(:,3)-nhat(:,3))**2 &
&    +(rhor(:,4)-nhat(:,4))**2)
   end if
  end if

! rhocorval will contain effective density used to compute gradients:
! - with core density (if NLCC)
! - without compensation density (if PAW under certain conditions)
! - in (up+dn,up) or (n,mx,my,mz) format according to collinearity
! of polarization and use of gradients (GGA)
  if (n3xccc>0.or.test_nhat.or.nspden_eff/=nspden) then
   allocate(rhocorval(nfft,nspden_eff))
   if (nspden==nspden_eff) then
    rhocorval(:,1:nspden)=rhor(:,1:nspden)
   else if (nspden==4) then
    rhocorval(:,1)=rhor(:,1)
    rhocorval(:,2)=half*(rhor(:,1)+m_norm(:))
   else
    rhocorval=zero
   end if
  end if

! Add core electron density to effective density
  if (n3xccc>0) then
   rhocorval(:,1)=rhocorval(:,1)+xccc3d(:)
   if(nspden_eff==2) then
    rhocorval(:,2)=rhocorval(:,2)+half*xccc3d(:)
   end if
  end if

! If PAW, substract compensation density from effective density:
! - if GGA, because nhat gradients are computed separately
! - if nhat does not have to be included in XC
  if (test_nhat) then
   if (nspden==nspden_eff) then
    rhocorval(:,1:nspden)=rhocorval(:,1:nspden)-nhat(:,1:nspden)
   else if (nspden==4) then
    if (usexcnhat==0) then
     rhocorval(:,1)=rhocorval(:,1)-nhat(:,1)
     rhocorval(:,2)=rhocorval(:,2)-half*nhat(:,2)
    else
     allocate(nhat_up(nfft))
     nhat_up(:)=half*(nhat(:,1)+(rhor(:,2)*nhat(:,2) &
&     +rhor(:,3)*nhat(:,3) &
&     +rhor(:,4)*nhat(:,4))/m_norm(ifft))
     rhocorval(:,1)=rhocorval(:,1)-nhat(:,1)
     rhocorval(:,2)=rhocorval(:,2)-nhat_up(:)
    end if
   end if
  end if

! rhonow will contain effective density (and gradients if GGA)
  allocate(rhonow(nfft,nspden_eff,ngrad*ngrad))

! ====================================================================
! Loop on unshifted or shifted grids
  do ishift=0,dtset%intxc

!  Set up density on unshifted or shifted grid (will be in rhonow(:,:,1)),
!  as well as the gradient of the density, also on the unshifted
!  or shifted grid (will be in rhonow(:,:,2:4)), if needed.
   if ((n3xccc==0).and.(.not.test_nhat).and.(nspden_eff==nspden)) then
    call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,dtset%paral_kgb,qphon,rhor,rhonow)
   else if ((ishift>0).and.(test_nhat)) then
    call xcden(cplex,gprimd,0,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,dtset%paral_kgb,qphon,rhocorval,rhonow)
   else
    call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,dtset%paral_kgb,qphon,rhocorval,rhonow)
   end if

!  -PAW+GGA: add "exact" gradients of compensation density
   if (test_nhat.and.usexcnhat==1) then
    if (ishift==0) then
     if (nspden==nspden_eff) then
      rhonow(:,1:nspden,1)=rhocorval(:,1:nspden)+nhat(:,1:nspden)
     else if (nspden==4) then
      rhonow(:,1,1)=rhocorval(:,1)+nhat(:,1)
      rhonow(:,2,1)=rhocorval(:,2)+nhat_up(:)
     end if
    else
     if (nspden==nspden_eff) then
      rhocorval(:,1:nspden)=rhocorval(:,1:nspden)+nhat(:,1:nspden)
     else if (nspden==4) then
      rhocorval(:,1)=rhocorval(:,1)+nhat(:,1)
      rhocorval(:,2)=rhocorval(:,2)+nhat_up(:)
     end if
     call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,1,nspden_eff,dtset%paral_kgb,qphon,rhocorval,rhonow)
    end if
    if (ngrad==2.and.nhatgrdim==1.and.nspden==nspden_eff) then
     do ii=1,3
      jj=ii+1
      do ispden=1,nspden
       do ifft=1,nfft
        rhonow(ifft,ispden,jj)=rhonow(ifft,ispden,jj)+nhatgr(ifft,ispden,ii)
       end do
      end do
     end do
    end if
   end if

!  Deallocate temporary arrays
   if (ishift==dtset%intxc) then
    if (n3xccc>0.or.test_nhat.or.nspden_eff/=nspden) deallocate(rhocorval)
    if (test_nhat.and.nspden/=nspden_eff.and.usexcnhat==1) deallocate(nhat_up)
   end if

!  In case of non-collinear magnetism, extract up and down density and gradients (if GGA)
   if (nspden==4.and.nspden_eff==nspden) then
    if (ngrad==2) then
     do ifft=1,nfft
      dmnorm(1:3)=zero
      if(m_norm(ifft)>rhonow(ifft,1,1)*tol10+tol14) then
       do jj=1,3  ! Compute here nabla(|m|)=(m.nabla(m))/|m|
        do ii=2,4
         dmnorm(jj)=dmnorm(jj)+rhonow(ifft,ii,1+jj)*rhonow(ifft,ii,1)
        end do
       end do
       dmnorm(1:3)=dmnorm(1:3)/m_norm(ifft)
      end if
      rhonow(ifft,2,2)=half*(rhonow(ifft,1,2)+dmnorm(1))
      rhonow(ifft,2,3)=half*(rhonow(ifft,1,3)+dmnorm(2))
      rhonow(ifft,2,4)=half*(rhonow(ifft,1,4)+dmnorm(3))
     end do
    end if
    rhonow(:,2,1)=half*(rhonow(:,1,1)+m_norm(:))
   end if

!  Make the density positive everywhere (but do not care about gradients)
   call mkdenpos(iwarn,nfft,nspden_updn,1,rhonow(:,1:nspden_updn,1))

!  write(6,*) 'rhonow',rhonow

!  Uses a block formulation, in order to save simultaneously
!  CPU time and memory : xc routines
!  are called only once over mpts times, while the amount of allocated
!  space is kept at a low value, even if a lot of different
!  arrays are allocated, for use in different xc functionals.
!  !!!  !$OMP PARALLEL LASTPRIVATE(mpts)
!  $OMP PARALLEL
!  $OMP SINGLE
#if defined OPENMP
   mpts=8000 / OMP_GET_NUM_THREADS()
#else
   mpts=4000
#endif
!  $OMP END SINGLE
!  $OMP END PARALLEL

!  The variable order indicates to which derivative of the energy
!  the computation must be done. Computing exc and vxc needs order=1 .
!  Meaningful values are 1, 2, 3. Lower than 1 is the same as 1, and larger
!  than 3 is the same as 3.
!  order=1 or 2 supported for all LSD and GGA ixc
!  order=3 supported only for ixc=3 and ixc=7
   order=1
   if(option==2.or.option==10.or.option==12)order=2
   if(option==-2)order=-2
   if(option==3)order=3

!  $OMP PARALLEL DO PRIVATE(coeff,drho_updn,dstrsxc,dvxcdgr,dvxci) &
!  $OMP&PRIVATE(d2vxcar,exci,grho2_updn,ifft,index,ipts,ispden) &
!  $OMP&PRIVATE(npts,rhoarr,rho_updn,rspts,strdiag) &
!  $OMP&PRIVATE(s1,s2,s3,vxci,zeta) &
!  $OMP&REDUCTION(+:enxc,strsxc1_tot,strsxc2_tot,strsxc3_tot) &
!  $OMP&REDUCTION(+:strsxc4_tot,strsxc5_tot,strsxc6_tot) &
!  $OMP&SHARED(dnexcdn,dtset%ixc,kxc,mpts,nfft,ngrad,nspden_updn,order) &
!  $OMP&SHARED(rhonow,vxc)
   do ifft=1,nfft,mpts
!   npts=mpts
!   npts is the number of points to be treated in this bunch
    npts=min(nfft-ifft+1,mpts)

!   Allocation of mandatory arguments of drivexc
    allocate(exci(npts))
    allocate(rhoarr(npts),rho_updn(npts,nspden_updn))
    allocate(rspts(npts),vxci(npts,nspden_updn),zeta(npts))
!   Allocation of optional arguments

    call size_dvxc(dtset%ixc,ndvxc,ngr2,nspden_updn,nvxcdgr,order)
    if (ndvxc/=0) allocate(dvxci(npts,ndvxc))
    if (nvxcdgr/=0) allocate(dvxcdgr(npts,nvxcdgr))
    ixc=dtset%ixc
    if ((ixc==3 .or. (ixc>=7 .and. ixc<=15) .or. ixc==23) .and. order==3) allocate(d2vxcar(npts))
    if (ngrad == 2) allocate(grho2_updn(npts,ngr2))

    do ipts=ifft,ifft+npts-1
!    index=ipts-ifft+1 varies from 1 to npts
     index=ipts-ifft+1
     rhoarr(index)=rhonow(ipts,1,1)
     if(nspden_updn==1)then
      rho_updn(index,1)=rhonow(ipts,1,1)*half
      if(dtset%xclevel==2)then
       grho2_updn(index,1)=quarter*(rhonow(ipts,1,2)**2+rhonow(ipts,1,3)**2+rhonow(ipts,1,4)**2)
      end if
     else
      rho_updn(index,1)=rhonow(ipts,2,1)
      rho_updn(index,2)=rhonow(ipts,1,1)-rhonow(ipts,2,1)
      if(dtset%xclevel==2)then
       grho2_updn(index,1)=rhonow(ipts,2,2)**2+   &
&       rhonow(ipts,2,3)**2+   &
&       rhonow(ipts,2,4)**2
       grho2_updn(index,2)=(rhonow(ipts,1,2)-rhonow(ipts,2,2))**2 +   &
&       (rhonow(ipts,1,3)-rhonow(ipts,2,3))**2 +   &
&       (rhonow(ipts,1,4)-rhonow(ipts,2,4))**2
       grho2_updn(index,3)=rhonow(ipts,1,2)**2+   &
&       rhonow(ipts,1,3)**2+   &
&       rhonow(ipts,1,4)**2
      end if
     end if


    end do

!   case with gradient
    if (dtset%xclevel==2)then
     if (order**2 <= 1 .or. ixc == 16) then
      if (ixc /= 13) then
       call drivexc(exci,ixc,npts,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,     &
&       grho2_updn=grho2_updn,vxcgr=dvxcdgr)
      else
       call drivexc(exci,ixc,npts,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,     &
&       grho2_updn=grho2_updn)
      end if
     else if (order /= 3) then
      if (ixc /= 13) then
       call drivexc(exci,ixc,npts,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,     &
&       dvxc=dvxci,grho2_updn=grho2_updn,vxcgr=dvxcdgr)
      else
       call drivexc(exci,ixc,npts,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,     &
&       dvxc=dvxci,grho2_updn=grho2_updn)
      end if
     else if (order == 3) then
      if (ixc /= 13) then
       call drivexc(exci,ixc,npts,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,     &
&       dvxc=dvxci,d2vxc=d2vxcar,grho2_updn=grho2_updn,vxcgr=dvxcdgr)
      else
       call drivexc(exci,ixc,npts,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,     &
&       dvxc=dvxci,d2vxc=d2vxcar,grho2_updn=grho2_updn)
      end if
     end if
!    cases without gradient
    else
     if (order**2 <=1 .or. ixc >= 31 .and. ixc<=34) then
      call drivexc(exci,ixc,npts,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr)
     else if (order==3 .and. (ixc==3 .or. ixc>=7 .and. ixc<=10)) then
      call drivexc(exci,ixc,npts,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,      &
&      dvxc=dvxci,d2vxc=d2vxcar)
     else
      call drivexc(exci,ixc,npts,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,      &
&      dvxc=dvxci)
     end if
    end if

!   Accumulate enxc, strsxc and store vxc (and eventually kxc)
    dstrsxc=zero
    do ipts=ifft,ifft+npts-1
     index=ipts-ifft+1
     enxc=enxc+rhoarr(index)*exci(index)
     dnexcdn(ipts,1)=vxci(index,1)
     if(nspden_updn==1)then
      strdiag=rhoarr(index)*(exci(index)-vxci(index,1))
     else if(nspden_updn==2)then
      dnexcdn(ipts,2)=vxci(index,2)
!     Note : this is not the complete Vxc in the GGA case
      strdiag=rhoarr(index)*exci(index) &
&      -rho_updn(index,1)*vxci(index,1)&
&      -(rhoarr(index)-rho_updn(index,1))*vxci(index,2)
     end if
     dstrsxc=dstrsxc+strdiag

!    For GGAs, additional terms appear
!    (the LB functional does not lead to additional terms)
     if(ngrad==2 .and. dtset%ixc/=13)then

!     Treat explicitely spin up, spin down and total spin for spin-polarized
!     Will exit when ispden=1 is finished if non-spin-polarized
      do ispden=1,3

       if(nspden_updn==1 .and. ispden>=2)exit

!      If the norm of the gradient vanishes, then the different terms vanishes,
!      but the inverse of the gradient diverges, so skip the update.
       if(grho2_updn(index,ispden) < 1.0d-24) then
        dnexcdn(ipts,ispden+nspden_updn)=zero
        cycle
       end if

!      Store the gradient of up, down or total density, depending on ispden and nspden
       if(nspden_updn==1)then
        drho_updn(1:3)=rhonow(ipts,1,2:4)
       else if(ispden==1 .and. nspden_updn==2)then
        drho_updn(1:3)=rhonow(ipts,2,2:4)
       else if(ispden==2 .and. nspden_updn==2)then
        drho_updn(1:3)=rhonow(ipts,1,2:4)-rhonow(ipts,2,2:4)
       else if(ispden==3 .and. nspden_updn==2)then
        drho_updn(1:3)=rhonow(ipts,1,2:4)
       end if

!      Compute the derivative of n.e_xc wrt the
!      spin up, spin down, or total density. In the non-spin-polarized
!      case take the coefficient that will be multiplied by the
!      gradient of the total density
       if(nspden_updn==1)then
!       Definition of dvxcdgr changed in v3.3
        if (nvxcdgr == 3) then
         coeff=half*dvxcdgr(index,1) + dvxcdgr(index,3)
        else
         coeff=half*dvxcdgr(index,1)
        end if
       else if(nspden_updn==2)then
        if (nvxcdgr == 3) then
         coeff=dvxcdgr(index,ispden)
        else if (ispden /= 3) then
         coeff=dvxcdgr(index,ispden)
        else if (ispden == 3) then
         coeff=zero
        end if
       end if
       dnexcdn(ipts,ispden+nspden_updn)=coeff

!      Compute the contribution to the stress tensor
       s1=-drho_updn(1)*drho_updn(1)*coeff
       s2=-drho_updn(2)*drho_updn(2)*coeff
       s3=-drho_updn(3)*drho_updn(3)*coeff
!      The contribution of the next line comes from the part of Vxc
!      obtained from the derivative wrt the gradient
       dstrsxc=dstrsxc+s1+s2+s3
       strsxc1_tot=strsxc1_tot+s1
       strsxc2_tot=strsxc2_tot+s2
       strsxc3_tot=strsxc3_tot+s3
       strsxc4_tot=strsxc4_tot-drho_updn(3)*drho_updn(2)*coeff
       strsxc5_tot=strsxc5_tot-drho_updn(3)*drho_updn(1)*coeff
       strsxc6_tot=strsxc6_tot-drho_updn(2)*drho_updn(1)*coeff

      end do
     end if
    end do

!   Transfer the xc kernel
    if(nkxc/=0 .and. nkxc/=23)then
     if (ndvxc==15.and.(option==10.or.option==12)) then
      if (nkxc>=3) then
       kxc(ifft:ifft+npts-1,1)=dvxci(1:npts,1)+dvxci(1:npts,9)
       kxc(ifft:ifft+npts-1,2)=dvxci(1:npts,10)
       kxc(ifft:ifft+npts-1,3)=dvxci(1:npts,2)+dvxci(1:npts,11)
       if (nkxc>3) kxc(ifft:ifft+npts-1,4:nkxc)=zero
      else
       kxc(ifft:ifft+npts-1,1)=half*(dvxci(1:npts,1)+dvxci(1:npts,9)+dvxci(1:npts,10))
       if (nkxc>1) kxc(ifft:ifft+npts-1,2:nkxc)=zero
      end if
     else if (nkxc<=ndvxc) then
      kxc(ifft:ifft+npts-1,1:nkxc)=dvxci(1:npts,1:nkxc)
     else
      if (allocated(dvxci)) kxc(ifft:ifft+npts-1,1:ndvxc)=dvxci(1:npts,1:ndvxc) !if ndvxc=0 dvxci is not allocated (PMA)
      kxc(ifft:ifft+npts-1,ndvxc+1:nkxc)=zero
     end if
     if (order==3) k3xc(ifft:ifft+npts-1)=d2vxcar(1:npts)
    else if(nkxc==23)then
     if (ndvxc==15) then
      kxc(ifft:ifft+npts-1,1:15)=dvxci(1:npts,1:15)
     else
      kxc(ifft:ifft+npts-1,1:ndvxc)=dvxci(1:npts,1:ndvxc)
      kxc(ifft:ifft+npts-1,ndvxc+1:15)=zero
     end if
     do ispden=1,nspden_updn
      do ii=1,4
       kxc(ifft:ifft+npts-1,13+ispden+2*ii)=rhonow(ifft:ifft+npts-1,ispden,ii)
      end do
     end do
    end if

!   Add the diagonal part to the xc stress
    strsxc1_tot=strsxc1_tot+dstrsxc
    strsxc2_tot=strsxc2_tot+dstrsxc
    strsxc3_tot=strsxc3_tot+dstrsxc

    deallocate(exci,rhoarr,rho_updn,rspts,vxci,zeta)
    if (allocated(dvxci)) deallocate(dvxci)
    if (allocated(dvxcdgr)) deallocate(dvxcdgr)
    if (allocated(d2vxcar)) deallocate(d2vxcar)
    if (allocated(grho2_updn)) deallocate(grho2_updn)


!   End of the loop on blocks of data
   end do
!  $OMP END PARALLEL DO

   strsxc(1)=strsxc1_tot
   strsxc(2)=strsxc2_tot
   strsxc(3)=strsxc3_tot
   strsxc(4)=strsxc4_tot
   strsxc(5)=strsxc5_tot
   strsxc(6)=strsxc6_tot

!  If GGA, multiply the gradient of the density by the proper
!  local partial derivatives of the XC functional
   if(ngrad==2 .and. dtset%ixc/=13)then
    call xcmult (dnexcdn,nfft,ngrad,nspden_eff,nspgrad,rhonow)
   end if

!  Compute contribution from this grid to vxc, and ADD to existing vxc
   if (nspden/=4) then
    call xcpot(cplex,dnexcdn,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,nspgrad,&
&    dtset%paral_kgb,qphon,rhonow,vxc)
   else

!   If non-collinear magnetism, restore potential in proper axis before adding it
    allocate(vxci(nfft,4));vxci=zero
    call xcpot(cplex,dnexcdn,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,nspgrad,&
&    dtset%paral_kgb,qphon,rhonow,vxci)
    if (usexcnhat==1.or.nhatdim==0) then
     do ifft=1,nfft
      dvdn=half*(vxci(ifft,1)+vxci(ifft,2))
      if(m_norm(ifft)>rhor(ifft,1)*tol10+tol14) then
       dvdz=half*(vxci(ifft,1)-vxci(ifft,2))/m_norm(ifft)
       vxc(ifft,1)=vxc(ifft,1)+dvdn+rhor(ifft,4)*dvdz
       vxc(ifft,2)=vxc(ifft,2)+dvdn-rhor(ifft,4)*dvdz
       vxc(ifft,3)=vxc(ifft,3)+rhor(ifft,2)*dvdz
       vxc(ifft,4)=vxc(ifft,4)-rhor(ifft,3)*dvdz
      else
       vxc(ifft,1:2)=vxc(ifft,1:2)+dvdn
      end if
     end do
    else
     do ifft=1,nfft
      dvdn=half*(vxci(ifft,1)+vxci(ifft,2))
      if(m_norm(ifft)>rhor(ifft,1)*tol10+tol14) then
       dvdz=half*(vxci(ifft,1)-vxci(ifft,2))/m_norm(ifft)
       vxc(ifft,1)=vxc(ifft,1)+dvdn+(rhor(ifft,4)-nhat(ifft,4))*dvdz
       vxc(ifft,2)=vxc(ifft,2)+dvdn-(rhor(ifft,4)-nhat(ifft,4))*dvdz
       vxc(ifft,3)=vxc(ifft,3)+(rhor(ifft,2)-nhat(ifft,2))*dvdz
       vxc(ifft,4)=vxc(ifft,4)-(rhor(ifft,3)-nhat(ifft,3))*dvdz
      else
       vxc(ifft,1:2)=vxc(ifft,1:2)+dvdn
      end if
     end do
    end if
    deallocate(vxci)
   end if

!  End loop on unshifted or shifted grids
  end do

! Normalize enxc, strsxc and vxc
  divshft=one/dble(dtset%intxc+1)
  strsxc(:)=strsxc(:)/dble(nfftot)*divshft
  if (dtset%usewvl == 0) then
   enxc=enxc*ucvol/dble(nfftot)*divshft
  else
   enxc = enxc * (dtset%wvl_hgrid / real(2, dp)) ** 3 * divshft
  end if
  do ispden=1,nspden
   do ifft=1,nfft
    vxc(ifft,ispden)=vxc(ifft,ispden)*divshft
   end do
  end do

! XG030514 : MPIWF Should reduce strsxc and enxc
! in the group of WF processors
! Init mpi_comm
  if(mpi_enreg%paral_compil_fft==1)then
   old_paral_level=mpi_enreg%paral_level
   mpi_enreg%paral_level=3
!  call xcomm_init(mpi_enreg,spaceComm)
   if(mpi_enreg%mode_para=='b')then
    spaceComm=mpi_enreg%comm_fft
    call timab(48,1,tsec)
    call xsum_mpi(strsxc,spaceComm ,ierr)
    call xsum_mpi(enxc,spaceComm ,ierr)
    call timab(48,2,tsec)
   end if
   mpi_enreg%paral_level=old_paral_level
  end if

! Compute vxcavg
  call mean_fftr(vxc,vxcmean,mpi_enreg,nfft,nfftot,min(nspden,2))

! write(6,*) 'vxcmean',vxcmean

  if(nspden==1)then
   vxcavg=vxcmean(1)
  else
   vxcavg=half*(vxcmean(1)+vxcmean(2))
  end if
  deallocate(dnexcdn,rhonow)
  if (nspden==4) deallocate(m_norm)

 else if (dtset%ixc/=20) then

! Not an allowed choice for ixc and nspden
  write(message, '(a,a,a,a,2i8,a)' ) ch10,&
&  ' rhohxc: BUG -',ch10,&
&  '  ixc,nspden=',dtset%ixc,nspden,' not supported.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')

 end if

!Treat separately the Fermi-Amaldi correction.
 if (dtset%ixc==20 .or. dtset%ixc==21 .or. dtset%ixc==22) then

! Fermi-Amaldi correction : minus Hartree divided by the
! number of electrons per unit cell. This is not size consistent, but
! interesting for isolated systems with a few electrons.
  nelect=ucvol*rhog(1,1)
  factor=-one/nelect
  vxc(:,1)=factor*vhartr(:)
  if(nspden>=2) then
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(factor,nfft,vhartr,vxc)
   do ifft=1,nfft
    vxc(ifft,2)=factor*vhartr(ifft)
   end do
!  $OMP END PARALLEL DO
  end if

! Compute corresponding xc energy and stress as well as vxcavg
  call dotprod_vn(1,rhor,enxc,doti,mpi_enreg,nfft,nfftot,1,1,vxc,ucvol)
  enxc=half*enxc
  strsxc(1:3)=-enxc/ucvol

! Compute average of vxc (one component only).
  call mean_fftr(vxc,vxcmean,mpi_enreg,nfft,nfftot,1)
  vxcavg = vxcmean(1)
! For ixc=20, the local exchange-correlation kernel is zero, but the Hartree
! kernel will be modified in tddft. No other use of kxc should be made with ixc==20
  if(nkxc/=0 .and. dtset%ixc==20) then
   do inkxc=1,nkxc
!   $OMP PARALLEL DO PRIVATE(ifft) &
!   $OMP&SHARED(inkxc,kxc,nfft,nkxc)
    do ifft=1,nfft
     kxc(ifft,inkxc)=zero
    end do
!   $OMP END PARALLEL DO
   end do
  end if
! For ixc=21 or 22, the LDA (ixc=1) kernel has been computed previously.

 end if

 call timab(81,2,tsec)

!DEBUG
!write(6,*)' rhohxc : enxc=',enxc
!stop
!write(6,*)' rhohxc : strsxc(1:3)=',strsxc(1:3)
!ENDDEBUG

!DEBUG
!write(6,'(a)') ' rhohxc :  '
!write(6,'(a)')&
!& '   ir              rhor(ir)      vxc(ir)       kxc(ir)     xccc3d(ir) '
!n1=ngfft(1)
!n2=ngfft(2)
!n3=ngfft(3)
!do ir=1,nfft
!! if(ir<=11 .or. mod(ir,301)==0 )then
!i3=(ir-1)/n1/n2
!i2=(ir-1-i3*n1*n2)/n1
!i1=ir-1-i3*n1*n2-i2*n1
!if(n3xccc>0)write(message,'(i5,3i3,a,4es14.6)')ir,i1,i2,i3,' ',&
!&   rhor(ir,1),vxc(ir,1),kxc(ir,1),xccc3d(ir)
!if(n3xccc==0)write(message,'(i5,3i3,a,4es14.6)')ir,i1,i2,i3,' ',&
!&   rhor(ir,1),vxc(ir,1),kxc(ir,1)
!call wrtout(06,message,'COLL')
!if(nspden_updn==2)then
!write(message,'(a,2es14.6)')'               ',rhor(ir,2),vxc(ir,2)
!call wrtout(06,message,'COLL')
!end if
!! end if
!end do
!stop
!ENDDEBUG

!DEBUG
!write(6,*)' rhohxc : debug, stop at exit '
!enxc=zero
!kxc(:,:)=zero
!strsxc(:)=zero
!vhartr(:)=zero
!vxc(:,:)=zero
!vxcavg=zero
!stop
!ENDDEBUG

!call leave_new("COLL")

end subroutine rhohxc

!An optimization of the xc part of this routine
!and the subroutines called by this routine
!has been performed already.
!The number of exponentiations, logarithms, divisions etc has been minimized.
!The results of tests 11-18, 21-24 are as follows, on a P6 at 200 MHz.
!All tests are with intxc=0, except the last one.
!A further optimization should focus on intxc=1 .
!Time is in microsec/data/call
!Test    time     Note
!11     2.643    non spin-pol Teter xcspol
!12     3.103    id           Perdew Zunger
!13     2.836    id           Teter xctetr
!14     2.670    id           Wigner
!15     3.815    id           Hedin-Lundqvist
!16     2.365    id           X-alpha
!17     4.476    id           Perdew Wang
!18    12.894    id           GGA PBE
!21     4.001    spin-pol     Teter xcspol
!22     6.012    id           Perdew Wang
!23    12.797    id           GGA PBE
!24    25.953    id           GGA PBE with intxc=1
!
!Tests where also lead to determine the time for some basic operations :
!(these are very approximative numbers)
!div     0.22
!sqrt    0.3
!log     0.5
!exp     1.3
!**dblep 3.0
!The routine invcb, that compute x**(-one/3.0d0), scores 1.4  .

!!***
