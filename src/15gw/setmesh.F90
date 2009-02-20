!{\src2tex{textfont=tt}}
!!****f* ABINIT/setmesh
!!
!! NAME
!! setmesh
!!
!! FUNCTION
!! Calculate the size of the FFT grid for the GW calculation.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, YMN, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mG0= number of shell that must be added to take into account umklapp processes
!!  enforce_sym = flag to enforce a FFT compatible with all the symmetry operations, both rotations and 
!!   fractional translations (should be default)
!!  gmet(3,3) = reciprocal space metric.
!!  gvec(3,npwvec) = G-vectors array.
!!  method = integer flag for FFT grid (see below)
!!  npwvec= number of G vectors in the array gvec max(npwwfn,npwsigx)
!!  npwsigx = size of the dielectric or self-energy matrix.
!!  npwwfn = number of G-vectors in the wavefunctions.
!!  Cryst<Crystal_structure>=data type gathering information on unit cell and symmetries
!!    %nsym = number of symmetry operations 
!!    %symrel(3,3,nsym) = symmetry operations in real space
!!    %tnons(3,nsym) = fractional translations
!!
!! OUTPUT
!! nfftot= ngfft(1)*ngfft(2)*ngfft(3)=number of points on the  FFT grid.
!! ngfft(18)=contain all needed information about 3D FFT,
!!  see ~abinit/doc/input_variables/vargs.htm#ngfft
!!
!! NOTE
!! Four methods have been coded for the calculation of the mesh (see parameter "method" below):
!!
!!  method=0 FFT mesh defined by the user
!!  method=1 roughly takes the FFT box that encloses the larger of the two spheres of radius
!!            aliasing_factor*rwfn and rsigx, where rwfn and rsigx are the radius of the spheres
!!            with npwwfn and npwsigx planewaves respectively. default aliasing_factor is 1
!!  method=2 calculates the optimal FFT grid that allows aliasing only outside the sphere of the
!!            npwsigx planewaves (finer than method=1 with aliasing_factor=1).
!!  method=3 calculates the FFT grid needed to expand the density.
!!            (even finer than method=2, roughly corresponds to method=1 with aliasing_factor=2).
!!
!! PARENTS
!!      mrgscr,screening,sigma
!!
!! CHILDREN
!!      sizefft,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine setmesh(gmet,gvec,ngfft,npwvec,npwsigx,npwwfn,nfftot,method,mG0,Cryst,enforce_sym)

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_15gw, except_this_one => setmesh
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: enforce_sym,method,npwsigx,npwvec,npwwfn
 integer,intent(out) :: nfftot
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 integer,intent(in) :: gvec(3,npwvec),mG0(3)
 integer,intent(inout) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3)

!Local variables ------------------------------
!scalars
 integer :: aliasing_factor,fftalg,ig,ig1,ig1max,ig2,ig2max,ig3,ig3max,ii,idx
 integer :: is,m1,m2,m3,mm1,mm2,mm3,n1,n2,n3,nsym,nt,unt
 real(dp) :: ecuteff,ecutsigx,ecutwfn,g1,g2,g3,gsq,gsqmax,reff,rsigx,rwfn
 logical :: fft_ok,fftfile,ltest
 character(len=500) :: msg
 character(len=fnlen) :: fnam
!arrays
 integer :: fftnons(3),fftsym(3),mdum(3)
 integer,pointer :: symrel(:,:,:)
 real(dp),pointer :: tnons(:,:)

!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' setmesh enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 ! === Test on input ===
 ltest=ALL(mg0>=0)
 write(msg,'(a,3i4)')' called with wrong value of mG0 = ',mG0(:)
 call assert(ltest,msg,__FILE__,__LINE__)

 nsym = Cryst%nsym
 symrel => Cryst%symrel
 tnons  => Cryst%tnons
 !
 ! Calculate the limits of the sphere of npwwfn G-vectors in each direction.
 m1=MAXVAL(ABS(gvec(1,1:npwwfn)))
 m2=MAXVAL(ABS(gvec(2,1:npwwfn)))
 m3=MAXVAL(ABS(gvec(3,1:npwwfn)))
 !
 ! Calculate the limits of the sphere of npsigx G-vectors in each direction.
 ! Ensure that G+G0 will fit into the FFT grid, where G is any of the npwsigx/npweps vectors 
 ! and G0 is (i,j,k) [-nG0shell<i,j,k<nG0shell]. This is required when npwsigx>npwwfn since 
 ! we have to take into account umklapp G0 vectors to evaluate the oscillator matrix elements 
 ! (see rho_tw_g) or to symmetrize these quantities (see also cigfft).
 mm1=MAXVAL(ABS(gvec(1,1:npwsigx)))
 mm2=MAXVAL(ABS(gvec(2,1:npwsigx)))
 mm3=MAXVAL(ABS(gvec(3,1:npwsigx)))
 mm1=mm1+mG0(1) ; mm2=mm2+mG0(2) ; mm3=mm3+mG0(3)

 write(msg,'(2(2a,i8,a,3i6),2a,3i3)')ch10,&
& ' setmesh: npwwfn           = ',npwwfn, '; Max (m1,m2,m3)    = ',m1,m2,m3,ch10,&
& '          npweps/npwsigx   = ',npwsigx,'; Max (mm1,mm2,mm3) = ',mm1,mm2,mm3,ch10,&
& '          mG0 added        = ',mG0(:)
 call wrtout(std_out,msg,'COLL')
 !
 ! Different FFT grids according to method.
 select CASE (method)

  CASE (0) 
   ! FFT mesh defined by user, useful for testing
   !fnam='__fft.in__'
   !inquire(file=fnam,exist=fftfile)
   !if (fftfile) then 
   ! unt=io_unused()
   ! open(file=fnam,unit=unt,form='formatted')
   ! read(unt,*) n1,n2,n3 
   ! close(unt)
   !else 
   ! write(msg,'(5a)')' setmesh : ERROR : ',ch10,' FFT file ',TRIM(fnam),' not found'
   ! call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
   !end if 
   n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
   write(msg,'(3a,3(a,i3))')ch10,&
&   ' setmesh: COMMENT - ',ch10,&
&   ' mesh size enforced by user = ',n1,'x',n2,'x',n3
   call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
   ngfft(1)=n1 
   ngfft(2)=n2 
   ngfft(3)=n3
   ngfft(4)=2*(ngfft(1)/2)+1
   ngfft(5)=2*(ngfft(2)/2)+1
   ngfft(6)=   ngfft(3)
   nfftot=n1*n2*n3
   RETURN

  CASE (1)
   aliasing_factor=1 
   write(msg,'(2a,i3)')ch10,' using method 1 with aliasing_factor = ',aliasing_factor
   call wrtout(std_out,msg,'COLL')
   m1=m1*aliasing_factor
   m2=m2*aliasing_factor
   m3=m3*aliasing_factor
 
  CASE (2,3)
   ! Calculate the radius of the sphere of npwwfn G-vectors.
   ecutwfn=-one
   do ig=1,npwwfn
    g1=REAL(gvec(1,ig))
    g2=REAL(gvec(2,ig))
    g3=REAL(gvec(3,ig))
    gsq=       gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
&         two*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
    ecutwfn=MAX(ecutwfn,gsq)
   end do
   rwfn=SQRT(ecutwfn) ; ecutwfn=two*ecutwfn*pi**2
   ! Calculate the radius of the sphere of npwsigx/npweps G-vectors.
   ecutsigx=-one
   do ig=1,npwsigx
    g1=REAL(gvec(1,ig))
    g2=REAL(gvec(2,ig))
    g3=REAL(gvec(3,ig))
    gsq=      gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
&        two*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
    ecutsigx=MAX(ecutsigx,gsq)
   end do
   rsigx=SQRT(ecutsigx) ; ecutsigx=two*ecutsigx*pi**2
   write(msg,'(a,f7.3,3a,f7.3,a)')&
&    ' calculated ecutwfn           = ',ecutwfn,' [Ha] ',ch10,&
&    ' calculated ecutsigx/ecuteps = ',ecutsigx,' [Ha]'
   call wrtout(std_out,msg,'COLL')
   ! 
   ! In the calculation of the GW self-energy or of the RPA dielectric matrix,
   ! we make products Rho_12(r)=u_1*(r) u_2(r) of wavefunctions whose Fourier
   ! coefficients lie in the sphere of radius rwfn. Such products will have non
   ! vanishing Fourier coefficients in the whole sphere of radius 2*rwfn since:
   ! Rho(G) = \sum_T u_1*(T) u_2(T+G) 
   ! However, we only need the Fourier coefficients of Rho_12 that lie in the sphere 
   ! of radius rsigx. we can thus allow aliasing outside that sphere, so that the FFT box 
   ! may only enclose a sphere of radius:

   !reff=min(rsigx+rwfn,2.0*rwfn) 
   reff=rsigx+rwfn
   !
   ! Extreme case: this yields back the GS FFT grid if full wavefunctions are considered
   if (method==3) reff=two*rwfn
   ecuteff=two*(pi*reff)**2
   write(msg,'(2a,i2,a,f7.3,a)')ch10,&
&   ' using method = ',method,' with ecuteff = ',ecuteff,' [Ha]'
   call wrtout(std_out,msg,'COLL')
   ! 
   ! Search the limits of that sphere in each direction...
   ! FIXME this might be wrong in case of augmentation in the FFT
   gsqmax=reff**2
   ig1max=2*m1+1
   ig2max=2*m2+1
   ig3max=2*m3+1
   ! this is the correct coding
   !£££ ig1max=MAX(2*m1+1,2*mm1+1,mm1+m1+1)
   !£££ ig1max=MAX(2*m1+1,2*mm1+1,mm1+m1+1)
   !£££ ig1max=MAX(2*m1+1,2*mm1+1,mm1+m1+1)
   m1=-1 ; m2=-1 ; m3=-1
   do ig1=0,ig1max
    do ig2=0,ig2max
     do ig3=0,ig3max
      g1=REAL(ig1)
      g2=REAL(ig2)
      g3=REAL(ig3)
      gsq=     gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
&         two*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
      if (gsq>gsqmax+tol6) CYCLE ! tol6 to improve portability
      m1=MAX(m1,ig1)
      m2=MAX(m2,ig2)
      m3=MAX(m3,ig3)
     end do
    end do
   end do

  CASE DEFAULT
   write(msg,'(4a)')ch10,&
&   ' setmesh : BUG- ',ch10,&
&   ' method > 3 or < 0 is not allowed in setmesh ' 
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end select 
 !print*,m1,m2,m3
 !print*,mm1,mm2,mm3
 !
 ! Warning on low npwwfn
 if (m1<mm1 .or. m2<mm2 .or. m3<mm3) then
  write(msg,'(8a)')ch10,&
&  ' setmesh: COMMENT -',ch10,&
&  '  Note that npwwfn is small with respect to npweps or with respect to npwsigx.',ch10,&
&  '  Such a small npwwfn is a waste:',ch10,&
&  '  You could raise npwwfn without loss in cpu time.'
  call wrtout(std_out,msg,'COLL')
 end if
 !
 ! Keep the largest of the m/mm and and find the FFT grid which is compatible 
 ! with the library and, if required, with the symmetry operations.
 m1=MAX(m1,mm1)
 m2=MAX(m2,mm2)
 m3=MAX(m3,mm3)

 if (enforce_sym==0) then 
  ! 
  ! Determine the best size for the FFT grid *without* considering the symm ops.
  ! Ideally n=2*m+1 but this  could not be allowed by the FFT library 
  call sizefft(m1,n1)
  call sizefft(m2,n2)
  call sizefft(m3,n3)
  ! 
  ! Calculate the number of FFT grid points
  nfftot=n1*n2*n3
  ! Check if the FFT is compatible, write ONLY a warning if it breaks the symmetry
  fftnons(1)=n1
  fftnons(2)=n2
  fftnons(3)=n3
  fft_ok=.TRUE.
  rd: do ii=1,3 
   do is=1,nsym 
    nt=divisor_tns(tnons(ii,is))
    if (((fftnons(ii)/nt)*nt) /= fftnons(ii)) then 
     fft_ok=.FALSE.
     exit rd
    end if 
   end do 
  end do rd
  if (.not.fft_ok) then 
   write(msg,'(5a)')ch10,&
&   ' setmesh: WARNING -',ch10,&
&   ' FFT mesh is not compatible with non-symmorphic translations',ch10
   call wrtout(std_out,msg,'COLL') !; call wrtout(ab_out, msg,'COLL')
  end if 
  ! Check only rotations
  if (.not.(check_rot_fft(nsym,symrel,n1,n2,n3))) then
   write(msg,'(5a)')ch10,&
&   ' setmesh: WARNING -',ch10,&
&   ' FFT mesh is not compatible with rotations',ch10
   call wrtout(std_out, msg,'COLL') !; call wrtout(ab_out,msg,'COLL')
  end if 
 else 
  ! 
  ! Determine the best size for the FFT grid considering symm ops.
  ! Ideally n=2*m+1 but this  could not be allowed by the FFT library 
  write(msg,'(2a)')ch10,&
&  ' finding a FFT mesh compatible with all the symmetries'
  call wrtout(std_out,msg,'COLL')
  fftnons(:)=1
  ! Found a FFT mesh compatible with the non-symmorphic operations 
  do ii=1,3
   fftnons(ii)=1
   do is=1,nsym 
    nt=divisor_tns(tnons(ii,is))
    if (((fftnons(ii)/nt)*nt)/=fftnons(ii)) fftnons(ii)=mcm(fftnons(ii),nt)
   end do 
  end do 
  write(msg,'(a,3i3)')' setmesh: divisor mesh',fftnons(:)
  call wrtout(std_out,msg,'COLL')
  ! 
  ! Check if also rotations preserve the grid.
  ! Initial guess from previous m values.
  call sizefft(m1,fftsym(1))
  call sizefft(m2,fftsym(2))
  call sizefft(m3,fftsym(3))
  mdum(1)=m1 ; mdum(2)=m2 ; mdum(3)=m3

  idx=0
  do ! If a FFT division gets too large the code stops in sizefft
   if ( check_rot_fft(nsym,symrel,fftsym(1),fftsym(2),fftsym(3)) .and.&
&       (MOD(fftsym(1),fftnons(1))==0) .and.&
&       (MOD(fftsym(2),fftnons(2))==0) .and.&
&       (MOD(fftsym(3),fftnons(3))==0)&         
&     ) exit 
   ii=MOD(idx,3)+1
   mdum(ii)=mdum(ii)+1
   call sizefft(mdum(ii),fftsym(ii))
   idx=idx+1
  end do 
  ! 
  ! Got a good FFT grid, Calculate the number of FFT grid points
  n1=fftsym(1) ; n2=fftsym(2) ; n3=fftsym(3)
  nfftot=n1*n2*n3

! !DEBUG dont uncomment this part for the time being
  if ( .not.( check_rot_fft(nsym,symrel,n1,n2,n3)).or. &
&       ( MOD(fftsym(1),fftnons(1))/=0) .and.      &
&       ( MOD(fftsym(2),fftnons(2))/=0) .and.      &
&       ( MOD(fftsym(3),fftnons(3))/=0) &
&    ) then 
   write(msg,'(a)')' setmesh : BUG during the generation of a symmetric FFT'
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end if  
! !ENDDEBUG
 end if ! enforce_sym

 write(msg,'(3(a,i3),2a,i8,a)')&
& ' setmesh: FFT mesh size selected  = ',n1,'x',n2,'x',n3,ch10,&
& '          total number of points  = ',nfftot,ch10
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

 ngfft(1)=n1 
 ngfft(2)=n2 
 ngfft(3)=n3
 ngfft(4)=2*(ngfft(1)/2)+1
 ngfft(5)=2*(ngfft(2)/2)+1
 ngfft(6)=   ngfft(3)
 !
 ! Check the value of fftalg ie ngfft(7), 
 ! presently only SG library is allowed, see sizefft.F90
 fftalg=ngfft(7)/100 
 if (fftalg==2 .or. fftalg==3 .or. fftalg==4) then 
  write(msg,'(8a)')ch10,&
&  ' fftmesh : ERROR - ',ch10,&
&  ' Presently only S. Goedecker routines are allowed in GW calculation',ch10,&
&  ' Action : check the value of fftalg (ngfft(7)) in your input file',ch10,&
&  ' or modify setmesh.F90 to be sure that the FFT mesh is compatible with the FFT library '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if 

#if defined DEBUG_MODE
 write(*,*)' ngfft after setmesh ' 
 call print_ngfft(ngfft)
 write(msg,'(a)')' setmesh : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine setmesh
!!***

!!****f* ABINIT/select_divisor_mesh
!! NAME
!! select_divisor_mesh
!!
!! FUNCTION
!!
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

subroutine select_divisor_mesh(n,nsym,tnons)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => select_divisor_mesh
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!Scalars 
!Arrays
!scalars
 integer,intent(in) :: nsym
!arrays
 integer,intent(inout) :: n(3)
 real(dp),intent(in) :: tnons(3,nsym)

!Local variables ------------------------------
!scalars
 integer :: is,ix,nt
 character(len=500) :: msg

!************************************************************************
!BEGIN EXECUTABLE SECTION

 n(:)=1
 do ix=1,3
  n(ix)= 1
  do is=1,nsym
   nt=divisor_tns(tnons(ix,is))
   !if(((n(ix) / nt) * nt) /= n(ix)) n(ix) = mcm(n(ix),nt)
   if (MOD(n(ix),nt)/=0) n(ix)=mcm(n(ix),nt)
  end do
 end do
 write(msg,'(a,i3)')' divisor mesh = ',n(:)
 call wrtout(std_out,msg,'COLL')

end subroutine select_divisor_mesh
!!***


!!****f* ABINIT/divisor_tns
!! NAME
!! divisor_tns
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

 integer function divisor_tns(dd)

 use defs_basis
!returns the denominator of the rational number d

 implicit none

!Arguments ------------------------------------
!Scalars 
!scalars
 real(dp),intent(in) :: dd

!Local variables ------------------------------
!scalars
 integer :: ii

!************************************************************************

 ii=1
 do
  if (ABS(dd*ii-NINT(dd*ii))<0.0001) then
   divisor_tns=ii
   RETURN
  end if
  ii=ii+1
 end do

end function divisor_tns 

 integer function mcm(ii,jj)
!returns the maximum common multiple of ii and jj

 use defs_basis

 implicit none

 integer,intent(in) :: ii,jj

!************************************************************************
 
 mcm=1
 do
  if ( ((mcm/ii)*ii)==mcm .and. ((mcm/jj)*jj)==mcm ) return
  mcm=mcm+1
 end do

end function mcm

 logical function check_rot_fft(nsym,symrel,nr1,nr2,nr3)
!check that the grid is compatible with the given rotations in real space

 use defs_basis

 implicit none

!Arguments
!Scalar
 integer,intent(in) :: nr1,nr2,nr3,nsym 
!Arrays
 integer,intent(in) :: symrel(3,3,nsym)

!local variables
 integer :: is 

!************************************************************************

 !
 ! The grid is compatible with the symmetries (only rotations) if 
 ! for each symmetry, each n_i and n_j ==> $n_i*R_{ij}/n_j$ is an integer
 !
 check_rot_fft=.true.
 do is=1,nsym
  if ( MOD(symrel(2,1,is)*nr2, nr1) /=0 .or. &
&      MOD(symrel(3,1,is)*nr3, nr1) /=0 .or. &
&      MOD(symrel(1,2,is)*nr1, nr2) /=0 .or. &
&      MOD(symrel(3,2,is)*nr3, nr2) /=0 .or. &
&      MOD(symrel(1,3,is)*nr1, nr3) /=0 .or. &
&      MOD(symrel(2,3,is)*nr2, nr3) /=0      &
&    ) then
   check_rot_fft=.false. 
   exit
  end if 
 end do

end function check_rot_fft
!!***
