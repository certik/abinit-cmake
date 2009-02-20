!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_coulombian
!! NAME
!! setup_coulombian
!!
!! FUNCTION
!! Perform general check and initialize the data type containing information on the cutoff technique
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Qmesh<Bz_mesh_type>= Info on q-point sampling
!!  Kmesh<Bz_mesh_type>= Info on k-point sampling
!!  Gsph<Gvectors_type>= info of G sphere
!!   %gmet(3,3)= metric in reciprocal space
!!   %gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!   %gvec=G vectors
!!  ng=number of G-vectors to be used to describe the coulombian interaction
!!  ngfft(18)=information on the (fine) FFT grid used for the density.
!!
!! OUTPUT
!!  Vcp <type Coulombian_type> datatype gathering information on the coulombian cutoff technique
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

subroutine setup_coulombian(Dtset,Gsph,Qmesh,Kmesh,ng,rprimd,ngfft,MPI_enreg,Vcp)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : GW_TOLQ0
 use m_errors, only : assert
 use m_io_tools, only : flush_unit
#if defined FC_PGI6
 ! Buggy PGI6 doesnt like OPERATOR
 use m_numeric_tools 
#else
 use m_numeric_tools, only : arth, geop, imin_loc, llsfit_svd, l2norm, OPERATOR(.x.)
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_15gw, except_this_one => setup_coulombian
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng
 type(Bz_mesh_type),intent(in) :: Kmesh,Qmesh
 type(Coulombian_type),intent(out) :: Vcp
 type(Dataset_type),intent(in) :: Dtset
 type(Gvectors_type),intent(in) :: Gsph
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!number of sampling points along each axis for the numerical integration
!scalars
 integer,parameter :: nqx=100
 integer :: ig,ii,info,ipara,iq_bz,iq_ibz,iqg,iqx,iqy,iqz,istat,nfft,npar,npt
 integer :: opt_cylinder,opt_surface,test
 real(dp),parameter :: tolq0=1.d-3
 real(dp) :: bz_geometry_factor,bz_plane,check,chisq,dx,integ,q0_vol,q0_volsph
 real(dp) :: qbz_norm,step,ucvol
 logical,parameter :: prior_to_56=.TRUE.
 logical :: ltest
 character(len=500) :: msg
!arrays
 integer :: gam(3,1)
 integer,pointer :: gvec(:,:)
 real(dp) :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),gmet(3,3),gprimd(3,3)
 real(dp) :: qbz_cart(3),qq(3),rmet(3,3)
 real(dp),allocatable :: cov(:,:),par(:),qfit(:,:),qpg(:),sigma(:),var(:)
 real(dp),allocatable :: vcfit(:,:),vcoul(:,:),xx(:),yy(:)

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' setup_coulombian: enter ' 
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif
 !
 ! === Test if the first q-point is zero === this wont work if nqptdm/=0
 !if (normv(Qmesh%ibz(:,1),gmet,'G')<GW_TOLQ0)) STOP 'setup_coulombian, non zero first point '

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 gvec => Gsph%gvec(:,:) 
 !
 ! === Save dimension and other useful quantities in Vcp% ===
 Vcp%nfft      = ngfft(1)*ngfft(2)*ngfft(3) ! Number of points in FFT mesh
 Vcp%ng        = ng                         ! Number of G-vectors in the coulombian matrix elements
 Vcp%nqibz     = Qmesh%nibz                 ! Number of irred q-point 
 Vcp%nqsmall   = Qmesh%nsmall               ! Number of small q-directions to deal with singularity and non Analytic behavior 
 Vcp%rcut      = Dtset%rcut                 ! Cutoff radius for cylinder
 Vcp%hcyl      = zero                       ! Length of finite cylinder (Rozzi"s method, default is Beigi)
 Vcp%ucvol     = ucvol                      ! Unit cell volume

 Vcp%rprimd    = rprimd(:,:)                ! Dimensional direct lattice
 Vcp%boxcenter = zero                       ! boxcenter at the moment is supposed to be at the origin
 !Vcp%boxcenter=Dtset%boxcenter     
 Vcp%vcutgeo   = Dtset%vcutgeo(:)           ! Info on the orientation and extension of the cutoff region
 Vcp%ngfft(:)  = ngfft(:)                   ! Info on FFT mesh
 !
 ! === Define geometry and cutoff radius (if used) ===
 Vcp%mode='None'
 if (Dtset%icutcoul==0) Vcp%mode='SPHERE'
 if (Dtset%icutcoul==1) Vcp%mode='CYLINDER'
 if (Dtset%icutcoul==2) Vcp%mode='SURFACE'
 if (Dtset%icutcoul==3) Vcp%mode='CRYSTAL'
 !
 ! === Calculate Fourier coefficient of coulombian interaction ===
 ! * For the limit q-->0 we consider ng vectors due to a possible anisotropy in case of a cutoff interaction
 allocate(Vcp%vc_sqrt(ng,Vcp%nqibz))     ; Vcp%vc_sqrt  =czero
 !print*,'nqsmal',Qmesh%nsmall ; call flush_unit(std_out)
 allocate(Vcp%vcqs_sqrt(ng,Vcp%nqsmall),STAT=istat) ; Vcp%vcqs_sqrt=czero  !STAT required by buggy g95 if called by sigma
 allocate(vcoul(ng,Vcp%nqibz)) 

 a1=rprimd(:,1) ; b1=two_pi*gprimd(:,1)
 a2=rprimd(:,2) ; b2=two_pi*gprimd(:,2)
 a3=rprimd(:,3) ; b3=two_pi*gprimd(:,3)

 SELECT CASE (TRIM(Vcp%mode)) 

 CASE ('SPHERE')
  ! TODO check that L-d > R_c > d
  ltest=(Vcp%rcut>tol12) 
  write(msg,'(a,f8.4)')' Negative cutoff radius = ',Vcp%rcut
  call assert(ltest,msg,__FILE__,__LINE__)
  Vcp%vcutgeo=zero
  call cutoff_sphere(Qmesh%nibz,Qmesh%ibz,ng,gvec,gmet,Vcp%rcut,vcoul)
  !Vcp%vc_sqrt=vcoul ; Vcp%vc_sqrt=SQRT(Vcp%vc_sqrt) 
  !deallocate(vcoul) ; allocate(vcoul(ng,Qmesh%nsmall))
  !call cutoff_sphere(Qmesh%nsmall,Qmesh%small,ng,gvec,gmet,Vcp%rcut,vcoul)
  !Vcp%vcqs_sqrt=vcoul ; Vcp%vcqs_sqrt=SQRT(Vcp%vcqs_sqrt) 
  !
  ! === Treat the limit q--> 0 ===
  ! * The small cube is approximated by a sphere, while vc(q=0)=2piR**2
  ! * if a single q-point is used, the expression for the volume is exact
  Vcp%i_sz=two_pi*Vcp%rcut**2
  call print_coulombian(Vcp) ; call print_coulombian(Vcp,unit=ab_out) 

 CASE ('CYLINDER')
  test=COUNT(Vcp%vcutgeo(:)/=zero) 
  call assert((test==1),'Wrong cutgeo for cylinder',__FILE__,__LINE__)
  ! === Beigi method is default i.e infinite cylinder of radius rcut ===
  opt_cylinder=1 ; Vcp%hcyl=zero ; Vcp%pdir(:)=0
  do ii=1,3
   check=Vcp%vcutgeo(ii)
   if (ABS(check)>zero) then 
    Vcp%pdir(ii)=1
    if (check<zero) then 
     ! * In case it is possible to use Rozzi method with finite cylinder
     opt_cylinder=2 ; Vcp%hcyl=ABS(check)*SQRT(SUM(rprimd(:,ii)**2))
    end if
   end if
  end do
  test=COUNT(Vcp%pdir==1) 
  call assert((test==1),'Wrong pdir for cylinder',__FILE__,__LINE__)
  if (Vcp%pdir(3)/=1) STOP "not implemented yet"
  call cutoff_cylinder(Qmesh%nibz,Qmesh%ibz,ng,gvec,Vcp%rcut,Vcp%hcyl,Vcp%pdir,&
&  Vcp%boxcenter,rprimd,vcoul,opt_cylinder,MPI_enreg)
  !Vcp%vc_sqrt=vcoul ; Vcp%vc_sqrt=SQRT(Vcp%vc_sqrt) 
  !deallocate(vcoul) ; allocate(vcoul(ng,Qmesh%nsmall))
  !call cutoff_cylinder(Qmesh%nsmall,Qmesh%small,ng,gvec,Vcp%rcut,Vcp%hcyl,Vcp%pdir,&
  !&Vcp%boxcenter,rprimd,vcoul,opt_cylinder,MPI_enreg)
  !Vcp%vcqs_sqrt=vcoul ; Vcp%vcqs_sqrt=SQRT(Vcp%vcqs_sqrt) 
  !
  ! === Treat the limit q--> 0 ===
  if (opt_cylinder==1) then 
   npar=8 ; npt=100 ; gam=RESHAPE((/0,0,0/),(/3,1/))
   allocate(qfit(3,npt),vcfit(1,npt)) ! TODO for the moment only z-axis
   if (Qmesh%nibz==1) STOP
   qfit(:,:)=zero 
   step=half/(npt*(Qmesh%nibz-1))              ; qfit(3,:)=arth(tol6,step,npt)
   !step=(half/(Qmesh%nibz-1)/tol6)**(one/npt)  ; qfit(3,:)=geop(tol6,step,npt)
   call cutoff_cylinder(npt,qfit,1,gam,Vcp%rcut,Vcp%hcyl,Vcp%pdir,Vcp%boxcenter,rprimd,vcfit,opt_cylinder,MPI_enreg)
   allocate(xx(npt),yy(npt),sigma(npt),par(npar),var(npar),cov(npar,npar))
   do ii=1,npt ; xx(ii)=normv(qfit(:,ii),gmet,'G') ; end do
   sigma=one ; yy(:)=vcfit(1,:)
   !call llsfit_svd(xx,yy,sigma,npar,K0fit,chisq,par,var,cov,info)
   !do ii=1,npt
   ! write(99,*)xx(ii),yy(ii),DOT_PRODUCT(par,K0fit(xx(ii),npar))
   !end do
   bz_plane=l2norm(b1.x.b2)
   !integ=K0fit_int(xx(npt),par,npar)
   !write(*,*)' SVD fit : chi-square',chisq 
   !write(*,*)' fit-parameters : ',par
   !write(*,*)' variance ',var
   !write(*,*)' bz_plane ',bz_plane
   !write(*,*)' SCD integ ',integ
   ! Here Im assuming homogeneous mesh
   dx=(xx(2)-xx(1))
   integ=yy(2)*dx*3.0/2.0
   do ii=3,npt-2
    integ=integ+yy(ii)*dx
   end do
   integ=integ+yy(npt-1)*dx*3.0/2.0
   write(*,*)' simple integral',integ
   q0_volsph=(two_pi)**3/(Kmesh%nbz*ucvol) 
   q0_vol=bz_plane*two*xx(npt)
   write(*,*)' q0 sphere : ',q0_volsph,' q0_vol cyl ',q0_vol
   Vcp%i_sz=bz_plane*two*integ/q0_vol
   write(*,*)' spherical approximation ',four_pi*7.44*q0_volsph**(-two_thirds) 
   write(*,*)' Cylindrical cutoff value ',Vcp%i_sz
   !Vcp%i_sz=four_pi*7.44*q0_vol**(-two_thirds) 
   deallocate(xx,yy,sigma,par,var,cov)
  end if
  call print_coulombian(Vcp) ; call print_coulombian(Vcp,unit=ab_out)
 CASE ('SURFACE') 
  ! TODO should check that R=L_Z/2 and fill Pcv%rcut the surface must be along x-y
  test=COUNT(Vcp%vcutgeo/=zero) ; if (test/=2) STOP ' check cutgeo '
  !
  ! === Default is Beigi"s method ===
  opt_surface=1  ; Vcp%alpha(:)=zero 
  if (ANY(Vcp%vcutgeo<zero)) opt_surface=2
  Vcp%pdir(:)=zero
  do ii=1,3
   check=Vcp%vcutgeo(ii)
   if (ABS(test)>zero) then 
    ! * In case it is possible to use Rozzi method with a finite surface along x-y
    Vcp%pdir(ii)=1
    if (test<zero) then 
     Vcp%alpha(ii)=normv(check*rprimd(:,ii),rmet,'R')
    end if
   end if
  end do
  call cutoff_surface(Qmesh%nibz,Qmesh%ibz,ng,gvec,gprimd,gmet,Vcp%rcut,&
&  Vcp%boxcenter,Vcp%pdir,Vcp%alpha,vcoul,opt_surface)
  !Vcp%vc_sqrt=vcoul ; Vcp%vc_sqrt=SQRT(Vcp%vc_sqrt) 
  !deallocate(vcoul) ; allocate(vcoul(ng,Qmesh%nsmall))
  !call cutoff_surface(Qmesh%nsmall,Qmesh%small,ng,gvec,gprimd,gmet,Vcp%rcut,vcoul,opt_surface)
  !Vcp%vcqs_sqrt=vcoul ; Vcp%vcqs_sqrt=SQRT(Vcp%vcqs_sqrt) 
  !
  ! === Treat the limit q--> 0 TODO 
  call print_coulombian(Vcp) ; call print_coulombian(Vcp,unit=ab_out)

 CASE ('CRYSTAL')
  do iq_ibz=1,Qmesh%nibz
   call cvc(Qmesh%nibz,iq_ibz,Qmesh%ibz,ng,gvec,gprimd,vcoul(:,iq_ibz))
  end do
  vcoul=four_pi/(vcoul**2)
  !call clcqpg(ng,gvec,gprimd,Qmesh%ibz,Qmesh%nibz,vcoul)
  !WHERE (vcoul>tol12) vcoul=four_pi/(vcoul**2)
  !Vcp%vc_sqrt=vcoul ; Vcp%vc_sqrt=SQRT(Vcp%vc_sqrt) 
  !deallocate(vcoul) ; allocate(vcoul(ng,Qmesh%nsmall))
  !call clcqpg(ng,gvec,gprimd,Qmesh%nsmall,Qmesh%small,vcoul)
  !Vcp%vcqs_sqrt=vcoul ; Vcp%vcqs_sqrt=SQRT(Vcp%vcqs_sqrt) 
  !
  ! === Treat 1/q^2 singularity === 
  ! * We use the auxiliary function from PRB 75, 205126 (2007)
  q0_vol=(two_pi)**3/(Kmesh%nbz*ucvol) ; bz_geometry_factor=zero
  do iq_bz=1,Qmesh%nbz
   qbz_cart(:)=Qmesh%bz(1,iq_bz)*b1(:)+Qmesh%bz(2,iq_bz)*b2(:)+Qmesh%bz(3,iq_bz)*b3(:)
   qbz_norm=SQRT(SUM(qbz_cart(:)**2))
   if (qbz_norm>tolq0) bz_geometry_factor=bz_geometry_factor-faux(Qmesh%bz(:,iq_bz))
  end do
  do iqx=1,nqx 
   do iqy=1,nqx 
    do iqz=1,nqx 
     qq(1)=(iqx-0.5)/REAL(nqx)-half
     qq(2)=(iqy-0.5)/REAL(nqx)-half
     qq(3)=(iqz-0.5)/REAL(nqx)-half
     bz_geometry_factor=bz_geometry_factor+faux(qq(:))*Qmesh%nbz/nqx**3
    end do
   end do
  end do
  write(msg,'(2a,2x,f8.4)')ch10,&
&  ' integrate q->0 : numerical BZ geometry factor = ',bz_geometry_factor*q0_vol**(2./3.)
  call wrtout(std_out,msg,'COLL') !; call wrtout(ab_out,msg,'COLL')
  Vcp%i_sz=four_pi*bz_geometry_factor  ! Final result stored here
  !
  ! ££ MG: this is to restore the previous implementation, it will facilitate the merge 
  ! Analytic integration of 4pi/q^2 over the volume element:
  ! $4pi/V \int_V d^3q 1/q^2 =4pi bz_geometric_factor V^(-2/3)$
  ! i_sz=4*pi*bz_geometry_factor*q0_vol**(-two_thirds) where q0_vol= V_BZ/N_k
  ! bz_geometry_factor: sphere=7.79, fcc=7.44, sc=6.188, bcc=6.946, wz=5.255  (see gwa.pdf, appendix A.4) 
if (prior_to_56) Vcp%i_sz=four_pi*7.44*q0_vol**(-two_thirds) 
  !call print_coulombian(Vcp) ; call print_coulombian(Vcp,unit=ab_out)

 CASE DEFAULT
  write(msg,'(5a)')ch10,&
&  ' setup_coulombian BUG - ',ch10,&
&  ' undefined cutoff mode ',TRIM(Vcp%mode)
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 END SELECT
 !
 ! === Store final results in complex array ===
 ! * Rozzi"s cutoff can give real negative values 
 Vcp%vc_sqrt=CMPLX(vcoul,zero) ; Vcp%vc_sqrt=SQRT(Vcp%vc_sqrt) 
 call plot_Vc(Vcp,Qmesh,Gsph,ng,vcoul,MPI_enreg)
 deallocate(vcoul)
 !
 ! === Create a table for FFT points inside the cutoff region ===
 if (TRIM(Vcp%mode)/='CRYSTAL') then
  allocate(Vcp%ctab(Vcp%nfft))
  call cutoff_table(Vcp)
 else 
  allocate(Vcp%ctab(1))
 end if

#if defined DEBUG_MODE
write(msg,'(a)')' setup_coulombian: exit ' 
call wrtout(std_out,msg,'COLL') 
call flush_unit(std_out)
#endif

contains !===============================================================

 function faux(qq)


  real(dp),intent(in) :: qq(3)
  real(dp) :: faux

  faux= four*( dot_product(b1,b1) * SIN(two_pi*qq(1)/two)**2 &
&             +dot_product(b2,b2) * SIN(two_pi*qq(2)/two)**2 &
&             +dot_product(b3,b3) * SIN(two_pi*qq(3)/two)**2 &
&            )                                              &
&       +two*( dot_product(b1,b2) * SIN(two_pi*qq(1))*SIN(two_pi*qq(2)) &
&             +dot_product(b2,b3) * SIN(two_pi*qq(2))*SIN(two_pi*qq(3)) &
&             +dot_product(b3,b1) * SIN(two_pi*qq(3))*SIN(two_pi*qq(1)) &
&            )
  faux=(two_pi)**2/faux
 end function faux

 function K0fit(mq,nn) result(vals)



!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15gw
!End of the abilint section

  integer,intent(in) :: nn
  real(dp),intent(in) :: mq
  real(dp) :: vals(nn)
 !Local variables-------------------------------
 !scalars
  integer :: ii
  real(dp) :: mqh
 !arrays
  real(dp),parameter :: cc(7)=(/-0.57721566,0.42278420,0.23069756,&
&                                0.03488590,0.00262698,0.00010750,0.00000740/) 
 ! *************************************************************************

  if (nn>8.or.nn<1) STOP 'not implemented'
  ! === Eq 9.8.5 in Abramovitz ===
  vals(1)=-LOG(mq*half)*I0(mq) 
  mqh=mq*half
  do ii=2,nn
   vals(ii)=cc(ii-1)*mqh**(2*(ii-2))
  end do
 end function K0fit

 function K0fit_int(mq,par,nn) result(integ)



!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15gw
!End of the abilint section

  integer,intent(in) :: nn
  real(dp),intent(in) :: mq
  real(dp) :: integ
  real(dp),intent(in) :: par(nn)
 !Local variables-------------------------------
 !scalars
  integer :: ii,aa
  real(dp) :: mqh
 !arrays
  real(dp),parameter :: cc(7)=(/-0.57721566,0.42278420,0.23069756,&
&                                0.03488590,0.00262698,0.00010750,0.00000740/) 
 ! *************************************************************************

  if (nn>8.or.nn<1) STOP 'not implemented'

  mqh=mq*half
  integ=-par(1)*int_I0ln(mqh)
  ! primitive of polynomial \sum_0^{N/2} cc_{2i} (x/2)^{2*i}
  do ii=2,nn
   aa=(2*(ii-1)+1)
   integ=integ+par(ii)*two*cc(ii-1)*(mqh**aa)/aa
  end do

 end function K0fit_int

 function I0(xx)


  real(dp),intent(in) :: xx
  real(dp) :: I0
  !Local variables-------------------------------
  real(dp) :: tt
  ! *************************************************************************

  ! Eq 9.8.1 of Abramovitz, entering the expansion of K0 -->0
  ! Expansion holds for |x|<3.75, Error<1.6*10D-07 
  tt=xx/3.75
  I0=one+3.5156229*tt**2+3.0899424*tt**4 +1.2067492*tt**6 &
        +0.2659732*tt**8+0.0360768*tt**10+0.0045813*tt**12
 end function I0
 
! Primitive of x^m Ln(x) for m/=-1
  function int_xmln(xx,mm)  result(res)

!Arguments ------------------------------------

   integer,intent(in) :: mm
   real(dp),intent(in) :: xx
   real(dp) :: res
! *********************************************************************

  if (mm==-1) STOP ' invalid value for mm '
  if (xx<=zero) STOP ' invalid value '

  res= (xx**(mm+1))/(mm+1) * (LOG(xx) - one/(mm+1))

  end function int_xmln

! Primitive function of ln(x/2)*I0(x) = sum_0^{N/2} 2^{2s+1} c_{2s} T(x/2,2s) 
! where T(x,s)=\int x^s ln(x)dx 
  function int_I0ln(xx) result(res)

!Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15gw
!End of the abilint section

  real(dp),intent(in) :: xx
  real(dp) :: res
!Local variables-------------------------------
  real(dp) :: yy
! *********************************************************************
 
  yy=xx*half
  res =  (       one*2    *int_xmln(yy,0)  &
&         +3.5156229*2**3 *int_xmln(yy,2)  &  
&         +3.0899424*2**5 *int_xmln(yy,4)  & 
&         +1.2067492*2**7 *int_xmln(yy,6)  &
&         +0.2659732*2**9 *int_xmln(yy,8)  & 
&         +0.0360768*2**11*int_xmln(yy,10) &
&         +0.0045813*2**13*int_xmln(yy,12) &
&        )

  end function int_I0ln

end subroutine setup_coulombian
!!***


!!****f* ABINIT/plot_Vc
!! NAME
!! plot_Vc
!!
!! FUNCTION
!! Plot vccut(q,G) as a function of |q+G|, calculate vc in real space

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

subroutine plot_Vc(Vcp,Qmesh,Gsph,ng,vc,MPI_enreg)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : get_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_15gw, except_this_one => plot_Vc
 use interfaces_lib00numeric
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng
 type(Bz_mesh_type),intent(in) :: Qmesh
 type(Coulombian_type),intent(in) :: Vcp
 type(Gvectors_type),intent(in) :: Gsph
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 real(dp),intent(in) :: vc(ng,Qmesh%nibz)

!Local variables-------------------------------
!scalars
 integer :: dimr,icount,idx_Sm1G,ierr,ig,igs,ii,iq_bz,iq_ibz,iqg,ir,isym,itim
 integer :: master,my_start,my_stop,nprocs,nqbz,nqibz,nr,ntasks,rank,spaceComm
 integer :: unt
 real(dp) :: arg,fact,l1,l2,l3,lmax,step,tmp,vcft,zmax
 character(len=500) :: msg
 character(len=fnlen) :: filnam
!arrays
 integer,allocatable :: insort(:)
 integer,pointer :: gvec(:,:)
 real(dp) :: b1(3),b2(3),b3(3),gmet(3,3),gprimd(3,3),qbz(3),qpgc(3)
 real(dp),allocatable :: kin(:),rr(:,:,:),vcr(:,:),vcr_cut(:,:)

!************************************************************************

 if (TRIM(Vcp%mode)=='CRYSTAL')  RETURN
 call xcomm_init  (MPI_enreg,spaceComm)
 call xme_init    (MPI_enreg,rank     )         
 call xmaster_init(MPI_enreg,master   ) 

 nqibz=Qmesh%nibz ; nqbz=Qmesh%nbz
 gmet(:,:)=Gsph%gmet(:,:) ; gprimd=Gsph%gprimd(:,:) ; gvec => Gsph%gvec(:,:)
 b1(:)=two_pi*gprimd(:,1)
 b2(:)=two_pi*gprimd(:,2)
 b3(:)=two_pi*gprimd(:,3)

 if (rank==master) then 
  allocate(insort(nqibz*ng),kin(nqibz*ng))
  iqg=1
  do iq_ibz=1,nqibz
   do ig=1,ng
    kin(iqg)=normv(Qmesh%ibz(:,iq_ibz)+gvec(:,ig),gmet,'g')
    insort(iqg)=iqg ; iqg=iqg+1
   end do
  end do
  call sort_dp(nqibz*ng,kin,insort,tol14)
  filnam='_VCoulFT_' ; call isfile(filnam,'new')
  unt=get_unit()
  open(unt,file=filnam,status='new',form='formatted')
  write(unt,'(a,i3,a,i6,a)')&
&  '#  |q+G|^2/2    q-points (',nqibz,')        gvec (',ng,' )  vc      vc_cutoff'
  do iqg=1,nqibz*ng
   iq_ibz=(insort(iqg)-1)/ng +1
   ig=(insort(iqg))-(iq_ibz-1)*ng
   write(unt,'(f12.6,2x,3f8.4,2x,3i6,2x,2es14.6)')&
&   kin(iqg),Qmesh%ibz(:,iq_ibz),gvec(:,ig),two_pi/kin(iqg)**2,vc(ig,iq_ibz)
  end do 
  close(unt)
  deallocate(insort,kin)
 end if

 nr=50 ; fact=one/(Vcp%ucvol*nqbz)
 ntasks=nqbz*ng
 call split_work(ntasks,my_start,my_stop)

 l1=SQRT(SUM(Vcp%rprimd(:,1)**2))
 l2=SQRT(SUM(Vcp%rprimd(:,2)**2))
 l3=SQRT(SUM(Vcp%rprimd(:,3)**2))
 lmax=MAX(l1,l2,l3) ; step=lmax/(nr-1)
 allocate(rr(3,nr,3),vcr(nr,3),vcr_cut(nr,3)) 
 !
 ! numb coding just to check cutoff implementation
 rr=zero
 do ii=1,3
  do ir=1,nr
   rr(ii,ir,ii)=(ir-1)*step
  end do
 end do

 vcr=zero ; vcr_cut=zero 
 do iq_bz=1,nqbz
  call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym,itim)
  if (ABS(qbz(1))<0.01) qbz(1)=zero
  if (ABS(qbz(2))<0.01) qbz(2)=zero
  if (ABS(qbz(3))<0.01) qbz(3)=zero
  igs=1 ; if (ALL(qbz(:)==zero)) igs=2
  do ig=igs,ng
   icount=ig+(iq_bz-1)*ng ; if (icount<my_start.or.icount>my_stop) CYCLE
   idx_Sm1G=Gsph%rottbm1(ig,itim,isym) ! IS{^-1}G
   vcft=vc(idx_Sm1G,iq_ibz)
   qpgc(:)=qbz(:)+gvec(:,ig) ; qpgc(:)=b1(:)*qpgc(1)+b2(:)*qpgc(2)+b3(:)*qpgc(3)
   tmp=SQRT(DOT_PRODUCT(qpgc,qpgc)) ; tmp=tmp**2
   do ii=1,3
    do ir=1,nr
     arg=DOT_PRODUCT(rr(:,ir,ii),qpgc)
     vcr_cut(ir,ii)=vcr_cut(ir,ii)+vcft*COS(arg)
     vcr(ir,ii)=vcr(ir,ii)+four_pi/tmp*COS(arg)
    end do
   end do
  end do !ig
 end do !iq_ibz
 call leave_test(MPI_enreg)
 call xsum_master(vcr_cut,master,spaceComm,ierr)
 call xsum_master(vcr    ,master,spaceComm,ierr)
 if (rank==master) then 
  filnam='_VCoulR_' ; call isfile(filnam,'new')
  unt=get_unit()
  open(unt,file=filnam,status='new',form='formatted')
  do ir=1,nr
   write(unt,'(7es18.6)')(ir-1)*step,(fact*vcr(ir,ii),fact*vcr_cut(ir,ii),ii=1,3)
  end do
  close(unt)
 end if
 deallocate(rr,vcr)

end subroutine plot_Vc
!!***


!!****f* ABINIT/print_coulombian
!! NAME
!! print_coulombian
!!
!! FUNCTION
!!  Print the content of a coulombian datatype
!!
!! INPUTS
!!  Vcp<Coulombian_type>=the datatype to be printed
!!  unit(optional)=the unit number for output 
!!  prtvol(optional)=verbosity level
!!  mode_paral=either "COLL" or "PERS"
!!
!! PARENTS
!!
!! CHILDREN
!!
!! OUTPUT
!!  Only printing 
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine print_coulombian(Vcp,unit,prtvol,mode_paral)

 use defs_basis
 use defs_datatypes
 use m_numeric_tools, only : imin_loc


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 type(Coulombian_type),intent(in) :: Vcp

!Local variables-------------------------------
!scalars
 integer :: ii,unt,verbose
 character(len=100) :: fmt
 character(len=4) :: mode
 character(len=500) :: msg

! *************************************************************************

 unt=std_out ; if (PRESENT(unit))       unt=unit
 verbose=0   ; if (PRESENT(prtvol))     verbose=prtvol
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 SELECT CASE (TRIM(Vcp%mode)) 
 CASE ('SPHERE')
  write(msg,'(5a,f8.4,3a,f8.2,3a,3f8.5,2a)')ch10,&
&  ' === Spherical cutoff === ',ch10,ch10,&
&  '  Cutoff radius ......... ',Vcp%rcut,' [Bohr] ',ch10,&
&  '  Volume of the sphere .. ',four_pi/three*Vcp%rcut**3,' [Bohr^3] ',ch10,&
&  '  Sphere centered at .... ',Vcp%boxcenter,' (r.l.u) ',ch10
  call wrtout(unt,msg,mode) 
 CASE ('CYLINDER')
  ii=imin_loc(ABS(Vcp%pdir-1))
  write(msg,'(5a,f8.4,3a,i2,2a,3f8.2,a)')ch10,&
&  ' === Cylindrical cutoff === ',ch10,ch10,&
&  '  Cutoff radius .......... ',Vcp%rcut,' [Bohr] ',ch10,&
&  '  Axis parallel to dir.... ',ii,ch10,&
&  '  Passing through point .. ',Vcp%boxcenter,' (r.l.u) '
  call wrtout(unt,msg,mode) 
  write(msg,'(2a)')'  Infinite length  ....... ',ch10
  if (Vcp%hcyl/=zero) write(msg,'(a,f8.5,2a)')'  Finite length of ....... ',Vcp%hcyl,' [Bohr] ',ch10
  call wrtout(unt,msg,mode) 
 CASE ('SURFACE') 
  write(msg,'(5a,f8.4,3a,3f8.2,2a)')ch10,&
&  ' === Surface cutoff === ',ch10,ch10,&
&  '  Cutoff radius .................... ',Vcp%rcut,' [Bohr] ',ch10,&
&  '  Central plane passing through .... ',Vcp%boxcenter,' (r.l.u) ',ch10
  call wrtout(unt,msg,mode) 
  !write(msg,'(a)')'  Infinite length  .......'
  !if (Vcp%hcyl/=zero) write(msg,'(a,f8.5,a)')'  Finite length of .......',Vcp%hcyl,' [Bohr] '
  !call wrtout(unt,msg,mode) 
 CASE ('CRYSTAL')
  write(msg,'(3a)')ch10,&
&  ' setup_coulombian : cutoff-mode = ',TRIM(Vcp%mode)
  call wrtout(unt,msg,mode) 
 CASE DEFAULT
  write(msg,'(5a)')ch10,&
&  ' setup_coulombian BUG - ',ch10,&
&  ' undefined cutoff mode ',TRIM(Vcp%mode)
  call wrtout(unt,msg,'COLL') ; call leave_new('COLL')
 END SELECT

 !TODO add additional information

end subroutine print_coulombian 
!!***


!!****f* ABINIT/cutoff_table
!! NAME
!! cutoff_table
!!
!! FUNCTION
!!  Create a table for each point in the FFT grid: 0 if the point is outside the cutoff region, 1 otherwise 
!!  This table can be used to redefine the density in the isolated system (see cutoff_density)
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).
!!  Vcp<type Coulombian_type> 
!!   The following quantities are used :
!!   %ngfft(1:3)= FFT dimensions
!!   %nr= number of point in the FFT mesh
!!   %boxcenter=reduced coordinates of the center of the box  (input variable
!!  (FIXME kept from tddft.F90, should write something in doc) 
!!   %cutoffmode= string of characters defining the cutoff mode 
!!   %rcut= cutoff radius 
!!
!! OUTPUT
!!  Vcp%ctab(Vcp%nr)
!!  For each point of the FFT grid gives 
!!   1 if the point in inside the cutoff region
!!   0 otherwise
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

subroutine cutoff_table(Vcp)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Coulombian_type),intent(inout) :: Vcp

!Local variables-------------------------------
!scalars
 integer :: ir,ix,iy,iz,ngfft1,ngfft2,ngfft3
 real(dp) :: ucvol
 character(len=500) :: msg
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rr(3),rxy(3),rz(3)

!************************************************************************

 if (TRIM(Vcp%mode)/='CRYSTAL') RETURN

 call metric(gmet,gprimd,-1,rmet,Vcp%rprimd,ucvol)
 ngfft1=Vcp%ngfft(1) 
 ngfft2=Vcp%ngfft(2) 
 ngfft3=Vcp%ngfft(3)

 Vcp%ctab(:)=1
 do ix=0,ngfft1-1
  rr(1)=DBLE(ix)/ngfft1
  do iy=0,ngfft2-1
   rr(2)=DBLE(iy)/ngfft2
   do iz=0,ngfft3-1
    rr(3)=DBLE(iz)/ngfft3
    ir=1+ix+iy*ngfft1+iz*ngfft1*ngfft2
    rr(:)=rr(:)-Vcp%boxcenter(:)

    SELECT CASE (TRIM(Vcp%mode))
     CASE ('SPHERE')
      if (normv(rr,rmet,'r')>Vcp%rcut) Vcp%ctab(ir)=0
     CASE ('CYLINDER') ! Infinite cylinder
      rxy(1)=rr(1)
      rxy(2)=rr(2)
      rxy(3)=zero
      if (normv(rxy,rmet,'r')>Vcp%rcut) Vcp%ctab(ir)=0
     CASE ('SURFACE')
      rz(:)=zero
      rz(3)=rr(3)
      if (normv(rz,rmet,'r')>Vcp%rcut) Vcp%ctab(ir)=0
     CASE DEFAULT   
      write(msg,'(6a)')ch10,&
&      ' cutoff_table : BUG -',ch10,&
&      ' cutoffmode ',TRIM(Vcp%mode),' not implemented'
      call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
    END SELECT

   end do 
  end do 
 end do

end subroutine cutoff_table
!!***


!!****f* ABINIT/cutoff_density
!! NAME
!! cutoff_density
!!
!! FUNCTION
!!  Modify density in case of calculations with Coulomb cutoff
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ngfft(18)=information on the FFT grid 
!!  (at the moment only the first three elements are used and are checked
!!  against the values contained in the Coulombian_type datatype)
!!  nspden=input variable 
!!  Vcp <type Coulombian_type>, only ctab(ngfft(1)*ngfft(2)*ngfft(3)) is used: 
!!   For each point of the FFT grid gives 
!!   1 if the point in inside the cutoff region, 0 otherwise
!!  ucvol=unit cell volume 
!!
!! OUTPUT
!!  Only printing
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

subroutine cutoff_density(ngfft,nspden,nsppol,Vcp,rhor,MPI_enreg)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nspden,nsppol
 type(Coulombian_type),intent(in) :: Vcp
 type(MPI_type),intent(in) :: MPI_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(ngfft(1)*ngfft(2)*ngfft(3),nspden)

!Local variables-------------------------------
!scalars
 integer :: ir,is,ix,iy,iz,master,ngfft1,ngfft2,ngfft3,rank
 real(dp),parameter :: EPS=1.d-2
 real(dp) :: fact
 character(len=500) :: msg
!arrays
 real(dp) :: nel(nspden),overflow(nspden),rr(3),rxy(3),rz(3)

! *************************************************************************
 
 ! === Do some check ===
 if (ANY(ngfft(1:3)/=Vcp%ngfft(1:3))) stop 'BUG in cutoff_density'
 if (Vcp%mode=='CRYSTAL') RETURN

 call xmaster_init(MPI_enreg,master   ) 
 call xme_init    (MPI_enreg,rank     )          

 ngfft1=ngfft(1) 
 ngfft2=ngfft(2) 
 ngfft3=ngfft(3)
 overflow(:)=zero ; nel(:)=zero

 do ix=0,ngfft1-1
  rr(1)=DBLE(ix)/ngfft1
  do iy=0,ngfft2-1
   rr(2)=DBLE(iy)/ngfft2
   do iz=0,ngfft3-1
    rr(3)=DBLE(iz)/ngfft3
    ir=1+ix+iy*ngfft1+iz*ngfft1*ngfft2

    nel(:)=nel(:)+rhor(ir,:)
    if (Vcp%ctab(ir)==0) then 
     overflow(:)=overflow+rhor(ir,:)
     !rhor(ir,:)=zero
    end if 

   end do 
  end do 
 end do

 nel(:)=nel(:)*Vcp%ucvol/(ngfft1*ngfft2*ngfft3)
 overflow(:)=overflow(:)*Vcp%ucvol/(ngfft1*ngfft2*ngfft3)

 if (rank==master) then
  write(*,*)' Number of electrons inside unit cell  = ',nel(1:nsppol)
  write(*,*)' Charge density outside cutoff region  =',overflow(1:nsppol)
  !
  ! If overflow is larger that few percent of the total charge warn or stop
  ! since one should increase the size of the supercell or the cutoff radius
  ! if (ANY(overflow(:)>EPS*nel(:)))  then 
  !  write(msg,'(4a,f8.5,3a)')ch10,&
  !&  ' cutoff_density : WARNING - ',ch10,&
  !&  '  More than ',eps,' % of the charge density is outside the cutoff region',ch10,&
  !&  '  You should try to increase the cutoff radius '
  !   call wrtout(std_out,msg,'COLL') 
  ! end if 
 end if

 !write(*,*)' restoring removed charge '
 ! Restore charge neutrality, spreading the overflow inside the cutoff region
 ! Using this approach if the charge is zero close to 
 ! the boundary then there should be no spurious effect
 !do is=1,nspden
 ! fact=nel(is)/(nel(is)-overflow(is))
 ! rhor(:,is)=rhor(:,is)*fact
 !end do

end subroutine cutoff_density
!!***


!!****f* ABINIT/destroy_coulombian
!! NAME
!! destroy_coulombian
!!
!! FUNCTION
!!  Destroy a Coulombian_type type 
!!
!! SIDE EFFECTS
!!  Vcp<Coulombian_type>=the datatype to be destroyed
!!
!! PARENTS
!!
!! CHILDREN
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine destroy_Coulombian(Vcp) 

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Coulombian_type),intent(inout) :: Vcp

!Local variables ------------------------------
!scalars
 integer :: istat
 character(len=500) :: msg

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' destroy_coulombian : enter '
 call wrtout(std_out,msg,'COLL')
#endif

 if (associated(Vcp%ctab     ))  deallocate(Vcp%ctab,     STAT=istat)
 if (associated(Vcp%vc_sqrt  ))  deallocate(Vcp%vc_sqrt,  STAT=istat)
 if (associated(Vcp%vcqs_sqrt))  deallocate(Vcp%vcqs_sqrt,STAT=istat)
 
end subroutine destroy_Coulombian
!!***
