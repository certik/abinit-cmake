!{\src2tex{textfont=tt}}
!!****f* ABINIT/move
!! NAME
!! move
!!
!! FUNCTION
!! make molecular dynamics updates
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, SE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amass(natom)=mass of each atom, in unit of electronic mass (=amu*1822...)
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cpus= cpu time limit in seconds)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem =number of k points which can fit in memory; set to 0 if use disk
!!   |  angular momentum for nonlocal pseudopotential
!!   | mpw=maximum dimensioned size of npw.
!!   | natom=number of atoms in unit cell
!!   |  except on first call (hartree/bohr); updated on output
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   |      for the "coarse" grid (see NOTES below)
!!   | nkpt=number of k points.
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  ecore=core psp energy (part of total energy) (hartree)
!!  hh=time step for molecular dynamics in atomic time units
!!   (1 atomic time unit=2.418884e-17 seconds)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mpi_enreg=informations about MPI parallelization
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!       for the "fine" grid (see NOTES below)
!!  npwarr(nkpt)=number of planewaves in basis and boundary at this k point.
!!  nspinor=number of spinorial components of the wavefunctions
!!  nattyp(ntypat)= # atoms of each type.
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  rprimd(3,3)=dimensional primitive translations (bohr)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points.
!!
!! SIDE EFFECTS
!! Rest of i/o is related to lda
!!  acell(3)=length scales of primitive translations (bohr)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=array for planewave coefficients of wavefunctions.
!!  densymop_gs <type(dens_sym_operator_type)>=the density symmetrization
!!   operator (ground-state symmetries)
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations (see initberry.f)
!!  ecore=core psp energy (part of total energy) (hartree)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  initialized= if 0 the initialisation of the gstate run is not yet finished
!!  irrzon(nfft**(1-1/nsym),2,nspden/nsppol)=irreducible zone data
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  occ(mband*nkpt*nsppol=occupation number for each band (usually 2) at each k point.
!!  phnons(2,nfft**(1-1/nsym),nspden/nsppol)=nonsymmorphic translation phases
!!  rhog(2,nfftf)=array for Fourier transform of electron density
!!  rhor(nfftf,nspden)=array for electron density in electrons/bohr**3.
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  wffnew,wffnow= struct info for wf disk files.
!!  vel(3,natom)=old value of velocity; updated on output
!!  xred(3,natom)=reduced dimensionless atomic coordinates; updated on output
!!  xred_old(3,natom)=work space for old xred
!!
!! NOTES
!! This subroutine uses the arguments natom, xred, vel, fcart, amass,
!! vis, and dtion (the last two contained in dtset) to make
!! molecular dynamics updates.  The rest of the lengthy
!! argument list supports the underlying lda computation
!! of forces, returned from subroutine scfcv
!!
!! USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      chkexi,leave_new,prtxvf,scfcv,status,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine move(acell,amass,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,&
&  ecore,eigen,hdr,indsym,initialized,irrzon,&
&  kg,mpi_enreg,&
&  nattyp,nfftf,npwarr,nspinor,occ,&
&  pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&  phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,&
&  scf_history,symrec,wffnew,wffnow,vel,wvl,xred,xred_old,ylm,ylmgr)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_15common
 use interfaces_21drive, except_this_one => move
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: pwind_alloc
 integer,intent(inout) :: initialized,nfftf,nspinor
 real(dp),intent(in) :: cpus,ecore
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(dens_sym_operator_type),intent(inout) :: densymop_gs
 type(efield_type),intent(inout) :: dtefield
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(results_gs_type),intent(out) :: results_gs
 type(scf_history_type),intent(inout) :: scf_history
 type(wffile_type),intent(inout) :: wffnew,wffnow
 type(wvl_data),intent(inout) :: wvl
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
 integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
!no_abirules
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem),nattyp(psps%ntypat)
 integer, intent(in) :: npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
 integer, intent(inout) :: symrec(3,3,dtset%nsym)
 real(dp),intent(inout) :: acell(3)
 real(dp), intent(in) :: amass(dtset%natom)
 real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
 real(dp), intent(out) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
 real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), pointer :: rhog(:,:),rhor(:,:)
 real(dp), intent(inout) :: xred(3,dtset%natom),xred_old(3,dtset%natom)
 real(dp), intent(inout) :: vel(3,dtset%natom),rprimd(3,3)
 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawrhoij_type), intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=4
 integer :: counter,iapp,icalls,iexit,ii,itime,jj,kk,ncalls,openexit,prtvel
 integer :: prtvol
 real(dp) :: aa,alfa,bb,cc,dv,dx,em,hh,vis,x0,xm
 character(len=500) :: message
!arrays
 real(dp),allocatable :: fcart(:,:),fnow(:,:),fprev(:,:),fprev2(:,:),fwork(:,:)
 real(dp),allocatable :: vnow(:,:),vprev(:,:),vprev2(:,:),vwork(:,:),xcart(:,:)
 real(dp),allocatable :: xnow(:,:),xprev(:,:),xprev2(:,:),xwork(:,:)

! *********************************************************************

!DEBUG
!write(6,*)' move : enter '
!ENDDEBUG

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 prtvol=dtset%prtvol
 if(prtvol==-level)then
  write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,&
&  ' move : enter '
  call wrtout(06,message,'COLL')
 end if

 vis=dtset%vis

!prtvel=1 makes prtxvf print the velocities
 prtvel=1

 allocate(xcart(3,dtset%natom),fcart(3,dtset%natom))
 allocate(xprev(3,dtset%natom), fprev(3,dtset%natom), vprev(3,dtset%natom) )
 allocate(xprev2(3,dtset%natom),fprev2(3,dtset%natom),vprev2(3,dtset%natom))
 allocate(xnow(3,dtset%natom),  fnow(3,dtset%natom),  vnow(3,dtset%natom)  )
 allocate(xwork(3,dtset%natom), fwork(3,dtset%natom), vwork(3,dtset%natom) )

 fnow(:,:)=0.0_dp

!Compute xcart from xred, and rprimd
 call xredxcart(dtset%natom,1,rprimd,xcart,xred)

 hh = dtset%dtion
 do itime=1,dtset%ntime

  call status(100*itime,dtfil%filstat,iexit,level,'loop itime    ')

! Will initialize the move with a 4th order Runge-Kutta start
  if(itime==1)then
   write(message,'(a)') ' Entering move for first time; initialize'
   call wrtout(06,message,'COLL')
   write(message,'(a)') '  with four calls to grad'
   call wrtout(06,message,'COLL')
   ncalls=4+1
  else
   ncalls=1
  end if

! ncalls is 5 in the initialisation period, and 1 after.
  do icalls=1,ncalls

   counter=100*itime+icalls
   call status(counter,dtfil%filstat,iexit,level,'loop icalls   ')

   write(message, '(a,a,i3,a)' ) ch10,' TIME STEP NUMBER ',itime,&
&   '  ------------------------------------------------------'
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
   if(itime==1) then
    write(message,'(a,i1,a)') &
&    ' (initialisation period : icalls=',icalls,')'
    call wrtout(ab_out,message,'COLL')
    call wrtout(06,  message,'COLL')
   end if

   if(prtvol==-level)then
    write(message,'(a,i4)')&
&    ' move : loop on calls, icalls= ',icalls
    call wrtout(06,message,'COLL')
   end if

!  Check whether exiting was required by the user.
!  if found then beat a hasty exit from the icalls loop
   openexit=1 ; if(dtset%chkexit==0) openexit=0
   call chkexi(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)
   if (iexit/=0) exit

!  Prepare new positions in the variable called xnow.
!  Also compute velocities at these  new positions

!  When the dynamical run is initialized, xcart contains
!  the original position , while after it, xcart contains
!  the previous positions.
!  
   do kk=1,dtset%natom
    do jj=1,3

!    Take first the atoms that are not allowed to move along this direction
!    Warning : implemented in cartesian coordinates
     if (dtset%iatfix(jj,kk)==1) then
      xnow(jj,kk)=xcart(jj,kk)
      vnow(jj,kk)=0.0_dp

     else

!     The first set of positions was the initializing set
      if(itime==1)then

!      First call to grad
       if(icalls==1)then
        xnow(jj,kk)=xcart(jj,kk)
        vnow(jj,kk)=vel(jj,kk)
       else if(icalls==2)then
        dx=hh*vel(jj,kk)
        dv=hh/amass(kk)*(fcart(jj,kk)-vis*vel(jj,kk))
        xnow(jj,kk)=xcart(jj,kk)+.5_dp*dx
        vnow(jj,kk)=vel(jj,kk)+.5_dp*dv
        xwork(jj,kk)=xcart(jj,kk)+sixth*dx
        vwork(jj,kk)=vel(jj,kk)+sixth*dv
       else if(icalls==3)then
        dx=hh*vprev(jj,kk)
        dv=hh/amass(kk)*(fprev(jj,kk)-vis*vprev(jj,kk))
        xnow(jj,kk)=xcart(jj,kk)+.5_dp*dx
        vnow(jj,kk)=vel(jj,kk)+.5_dp*dv
        xwork(jj,kk)=xwork(jj,kk)+third*dx
        vwork(jj,kk)=vwork(jj,kk)+third*dv
       else if(icalls==4)then
        dx=hh*vprev(jj,kk)
        dv=hh/amass(kk)*(fprev(jj,kk)-vis*vprev(jj,kk))
        xnow(jj,kk)=xcart(jj,kk)+dx
        vnow(jj,kk)=vel(jj,kk)+dv
        xwork(jj,kk)=xwork(jj,kk)+third*dx
        vwork(jj,kk)=vwork(jj,kk)+third*dv
       else if(icalls==5)then
        dx=hh*vprev(jj,kk)
        dv=hh/amass(kk)*(fprev(jj,kk)-vis*vprev(jj,kk))
        xnow(jj,kk)=xwork(jj,kk)+sixth*dx
        vnow(jj,kk)=vwork(jj,kk)+sixth*dv
       end if

      else if(abs(vis)<=1.d-8)then
!      If the viscosity is too small, the equations become ill conditioned
!      due to rounding error so do regular Verlet predictor Numerov corrector.
       xnow(jj,kk)=2._dp*xprev(jj,kk)-xprev2(jj,kk)&
&       +hh**2/amass(kk)*fprev(jj,kk)

      else if(abs(vis)>1.d-8)then
!      These equations come from solving m*d2x/dt2+vis*dx/dt=a+b*t+c*t**2
!      analytically under the boundary conditions that x(0)=x0 and
!      x(-h)=xm, and the following is the expression for x(h).
!      a, b and c are determined from our knowledge of the driving forces.
       em=amass(kk)
       aa=fprev(jj,kk)
       bb=(fprev(jj,kk)-fprev2(jj,kk))/hh
       x0=xprev(jj,kk)
       xm=xprev2(jj,kk)
       alfa=exp(-vis*hh/em)
       xnow(jj,kk)=( (-aa*hh*vis**2 +0.5_dp*bb*hh**2*vis**2&
&       +em*bb*hh*vis +x0*vis**3 -xm*vis**3)*alfa&
&       +aa*hh*vis**2 -em*bb*hh*vis&
&       +0.5_dp*bb*hh**2*vis**2 +x0*vis**3         )/vis**3

!      End of choice between initialisation, damped dynamics
!      and non-damped dynamics
!      
      end if

!     End the choice of moving/non-moving atoms
     end if

!    End loops on direction and atoms
    end do
   end do

   call status(counter,dtfil%filstat,iexit,level,'call scfcv    ')

!  Convert input xcart (cartesian) to reduced coordinates xred
   call xredxcart(dtset%natom,-1,rprimd,xnow,xred)

!  Compute LDA forces (big loop) : stored in fnow
   iapp=itime
   if(itime==1 .and. icalls/=5 )iapp=-icalls-1
   call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
&   eigen,hdr,iapp,indsym,initialized,&
&   irrzon,kg,mpi_enreg,&
&   nattyp,nfftf,npwarr,nspinor,occ,&
&   pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&   phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,&
&   scf_history,symrec,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)

!  Update of x, v and f
   do kk=1,dtset%natom
    do jj=1,3

     fnow(jj,kk)=results_gs%fcart(jj,kk)

!    After first step, "fcart" is computed
     if(itime==1 .and. icalls==1)then
      fcart(jj,kk)=fnow(jj,kk)
     end if

!    At the end of the initialisation phase, the original x, v, f
!    are considered as "previous values"
     if(icalls==5)then
      xprev(jj,kk)=xcart(jj,kk)
      fprev(jj,kk)=fcart(jj,kk)
      vprev(jj,kk)=vel(jj,kk)
     end if

!    Uses a corrector to have better value of xnow, and derive vnow.
!    Only update atoms position anf velocity along its allowed directions
     if(itime/=1 .and. dtset%iatfix(jj,kk)==0)then
!     Case when no viscosity is present
      if(abs(vis)<=1.d-8)then
       em=amass(kk)
       aa=fprev(jj,kk)
       bb=(fnow(jj,kk)-fprev2(jj,kk))/(2._dp*hh)
       cc=(fnow(jj,kk)+fprev2(jj,kk)-2._dp*fprev(jj,kk))/(2._dp*hh*hh)
       x0=xprev(jj,kk)
       xm=xprev2(jj,kk)
       xnow(jj,kk)=2._dp*xprev(jj,kk)-xprev2(jj,kk)+hh**2/amass(kk)/12._dp*&
&       (fprev2(jj,kk)+10._dp*fprev(jj,kk)+fnow(jj,kk))
       vnow(jj,kk)=(bb*hh**2)/(3._dp*em)&
&       +1.5_dp*aa*hh/em+(5._dp/12._dp)*cc*hh**3/em&
&       +x0/hh-xm/hh
      else
       em=amass(kk)
       aa=fprev(jj,kk)
       bb=(fnow(jj,kk)-fprev2(jj,kk))/(2._dp*hh)
       cc=(fnow(jj,kk)+fprev2(jj,kk)-2._dp*fprev(jj,kk))/(2._dp*hh*hh)
       x0=xprev(jj,kk)
       xm=xprev2(jj,kk)
       alfa=exp(-vis*hh/em)
       xnow(jj,kk)=((-aa*hh*vis**2+0.5_dp*bb*hh**2*vis**2&
&       -third*cc*hh**3*vis**2+em*bb*hh*vis&
&       -em*cc*hh**2*vis-2._dp*em**2*cc*hh+x0*vis**3-xm*vis**3)*alfa&
&       +aa*hh*vis**2-em*bb*hh*vis+third*cc*hh**3*vis**2&
&       +2._dp*em**2*cc*hh+0.5D0*bb*hh**2*vis**2-em*cc*hh**2*vis+x0*vis**3)&
&       /vis**3
       vnow(jj,kk)=(em*aa*vis**2*alfa-em*aa*vis**2+bb*hh*vis**2*em*alfa&
&       -bb*hh*vis**2*em+cc*hh**2*vis**2*em*alfa-cc*hh**2*vis**2*em&
&       -em**2*bb*vis*alfa+em**2*bb*vis-2._dp*em**2*cc*hh*vis*alfa+&
&       2._dp*em**2*cc*hh*vis+2._dp*em**3*cc*alfa-2._dp*em**3*cc+&
&       vis**3*alfa**2*aa*hh-0.5_dp*vis**3*alfa**2*bb*hh**2+&
&       third*vis**3*alfa**2*cc*hh**3-vis**2*&
&       alfa**2*em*bb*hh+vis**2*alfa**2*em*cc*hh**2+&
&       2._dp*vis*alfa**2*em**2*cc*hh-vis**4*alfa**2*x0+&
&       vis**4*alfa**2*xm)/vis**3/(alfa-1._dp)/em

      end if

     end if

!    After initialisation phase,
!    save the "now" values of x and v in xcart and vel.
     if(itime/=1 .or. icalls==5 )then
      xcart(jj,kk)=xnow(jj,kk)
      fcart(jj,kk)=fnow(jj,kk)
      vel(jj,kk)=vnow(jj,kk)
     end if

!    In all cases, copy "prev" in "prev2" and "now" in "prev"
     xprev2(jj,kk)=xprev(jj,kk)
     fprev2(jj,kk)=fprev(jj,kk)
     vprev2(jj,kk)=vprev(jj,kk)
     xprev(jj,kk)=xnow(jj,kk)
     fprev(jj,kk)=fnow(jj,kk)
     vprev(jj,kk)=vnow(jj,kk)

!    End update loops on atom number and direction
    end do
   end do

   call status(counter,dtfil%filstat,iexit,level,'call prtxvf   ')

!  Print values of x, v, and f, when adequate
   prtvel=1
   if(icalls/=1)prtvel=0
   call prtxvf(fnow,dtset%iatfix,ab_out,dtset%natom,prtvel,vnow,xnow)
   call prtxvf(fnow,dtset%iatfix, 06 ,dtset%natom,prtvel,vnow,xnow)

!  End loop " do icalls=1,ncalls "
  end do

  if(iexit/=0)exit

! End loop on itime
 end do

!Compute xred from xcart
 call xredxcart(dtset%natom,-1,rprimd,xcart,xred)

!Save fcart
 results_gs%fcart(:,:)=fcart(:,:)

 call status(0,dtfil%filstat,iexit,level,'dealloc       ')

 deallocate(xcart,fcart)
 deallocate(xnow,fnow,vnow)
 deallocate(xprev,fprev,vprev)
 deallocate(xprev2,fprev2,vprev2)
 deallocate(xwork,fwork,vwork)

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
  write(message,'(a1,a,a1,a,i1,a)')ch10,&
&  ' move : exit ',&
&  ch10,'  prtvol=-',level,', debugging mode => stop '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

end subroutine move

!!***
