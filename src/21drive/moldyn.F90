!{\src2tex{textfont=tt}}
!!****f* ABINIT/moldyn
!! NAME
!! moldyn
!!
!! FUNCTION
!! perform dynamics on ions according to ionmov (see help file)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, JCC, JYR, SE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amass(natom)=mass of each atom, in unit of electronic mass (=amu*1822...)
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cpus= cpu time limit in seconds
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem =number of k points which can fit in memory; set to 0 if use disk
!!   | mpw=maximum dimensioned size of npw.
!!   | natom=number of atoms in unit cell
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   |      for the "coarse" grid (see NOTES below)
!!   | nkpt=number ofidum k points.
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  ecore=core psp energy (part of total energy) (hartree)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)= # atoms of each type.
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!       for the "fine" grid (see NOTES below)
!!  npwarr(nkpt)=number of planewaves in basis and boundary at this k point.
!!  nspinor=number of spinorial components of the wavefunctions
!!  mxfh=last dimension of the xfhist array
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data paw
!!                                                 tabulated data read at start
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
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
!!  acell(3)=length scales of primitive translations (bohr)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
!!  densymop_gs <type(dens_sym_operator_type)>=the density symmetrization
!!   operator (ground-state symmetries)
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations (see initberry.f)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  initialized= if 0 the initialization of the gstate run is not yet finished
!!  irrzon(nfft**(1-1/nsym),2,nspden/nsppol)=irreducible zone data
!!  nxfh=actual number of (x,f) history pairs, see xfhist array.
!!  occ(mband*nkpt*nsppol)=occup number for each band (often 2) at each k point.
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  phnons(2,nfft**(1-1/nsym),nspden/nsppol)=nonsymmorphic translation phases
!!  rhog(2,nfftf)=array for Fourier transform of electron density
!!  rhor(nfftf,nspden)=array for electron density in electrons/bohr**3.
!!  rprim(3,3)=dimensionless real space primitive translations
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  vel(3,natom)=cartesian velocities at the initialisation; updated on output
!!  wffnew,wffnow=struct info for wf disk files.
!!  xfhist(3,natom+4,2,mxfh)=(x,f) history array,
!!                                 also includes acell, rprim and stress
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)=work space for old xred
!!
!! NOTES
!! For ionmov=6 :
!! Given a starting point xred that is a vector of length 3*natom
!! (reduced nuclei coordinates), a velocity vector (in cartesian
!! coordinates), and unit cell parameters (acell and rprim - without
!! velocities in the present implementation), the
!! Verlet dynamics is performed, using the gradient
!! of the energy (atomic forces and
!! stress : fred or fcart and stress) as calculated by the routine scfcv.
!! Some atoms can be kept fixed, while the propagation of unit cell
!! parameters is only performed if optcell/=0.
!! No more than dtset%ntime steps are performed.
!! The time step is governed by dtion (contained in dtset)
!! Returned quantities are xred, and eventually acell and rprim (new ones!).
!!
!! For ionmov=7 :
!! Block every atom for which the scalar product of velocity and
!! forces is negative, in order to reach the minimum.
!! The convergence requirement on
!! the atomic forces, dtset%tolmxf,  allows an early exit.
!!
!! For  ionmov=8
!! See ionmov=6, but with a nose-hoover thermostat
!! Velocity verlet algorithm : Swope et al JCP 76 (1982) 637
!!
!! For  ionmov=9
!! Uses a Langevin dynamics algorithm :
!! see J. Chelikowsky, J. Phys. D : Appl Phys. 33(2000)R33
!!
!! For  ionmov=12
!!    Application of Gauss' principle of least constraint according
!!    to Fei Zhang's algorithm (J. Chem. Phys. 106, 1997, p.6102);
!!    see also Minary et al. (J. Chem. Phys. 118, 2003, p.2510)
!!
!! NOTE : there are many similarities between this routine
!! and brdmin.f, so that they might have to be maintained
!! together. Some common subroutines might be extracted from them.
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
!!      chkexi,fconv,initylmg,leave_new,metric,mkrdim,prtxvf,scfcv,status
!!      write_header_moldynnetcdf,write_moldynvaluenetcdf,wrtout,xfpack
!!      xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine moldyn(acell,amass,atindx,atindx1,cg,cpus,densymop_gs,&
& dtefield,dtfil,dtset,ecore,eigen,hdr,indsym,initialized,&
& irrzon,kg,mpi_enreg,mxfh,nattyp,nfftf,npwarr,nspinor,nxfh,occ,&
& pawang,pawfgr,pawrad,pawrhoij,pawtab,&
& phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprim,&
& scf_history,symrec,wffnew,wffnow,vel,wvl,xfhist,xred,xred_old,ylm,ylmgr)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_13ionetcdf
 use interfaces_13recipspace
 use interfaces_15common
 use interfaces_16geomoptim
 use interfaces_21drive, except_this_one => moldyn
 use interfaces_lib00numeric
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
! nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
!scalars
 integer,intent(in) :: mxfh,pwind_alloc
 integer,intent(inout) :: initialized,nfftf,nspinor,nxfh
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
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),nattyp(psps%ntypat)
 integer,intent(in) :: npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
 integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
 integer,intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
 integer,intent(inout) :: symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: amass(dtset%natom),pwnsfac(2,pwind_alloc)
 real(dp),intent(inout) :: acell(3)
 real(dp),intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
 real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
 real(dp),intent(inout) :: rprim(3,3),vel(3,dtset%natom)
 real(dp),intent(inout) :: xfhist(3,dtset%natom+4,2,mxfh),xred(3,dtset%natom)
 real(dp),intent(inout) :: xred_old(3,dtset%natom)
 real(dp),intent(inout) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(out) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),pointer :: rhog(:,:),rhor(:,:)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: level=5,lwork=8
 integer :: counter,iapp,iatom,iatom1,iatom2,idim,idir,idir1,idir2,idum=-5,ierr
 integer :: iexit,ii,ionmov,ios,ipos,istopped,itime,itypat,ixfh,jatom,jdir,jj
 integer :: mcfac,nb1,nbdir,ndim,ndim0,nstopped,ntime,openexit,optcell,option
 integer :: prtvel,prtvol,routine
 real(dp),parameter :: esh2=one/six,esh4=esh2/20._dp,esh6=esh4/42._dp
 real(dp),parameter :: esh8=esh6/72._dp,nosetol=tol10,v2tol=tol8
 real(dp) :: a,as,b,cibinose,delxi,diag,dist,distx,disty,distz,dnose,dtion,ekin
 real(dp) :: ekin_corr,etotal_prev,etotal_temp,favg,fsnose,gnose,hamnose
 real(dp) :: hzeronose,ktemp,massvol,maxp1,maxp2,mttk_aloc,mttk_aloc2,mttk_bloc
 real(dp) :: norm_gauss,noseinert,polysh,ran_num1,ran_num2,rescale_vel,s,s1,s2
 real(dp) :: scdot,scprod,sig_gauss,sigma2,snose,sqb,taylor,ucvol,ucvol0
 real(dp) :: ucvol_next,v2gauss,v2nose,vlogv,vtest,xi_nose,xin_nose,xio
 logical :: ready
 character(len=500) :: message
 type(mttk_type) :: mttk_vars
!arrays
 integer,allocatable :: imax_perm(:),natom_type(:),stopped(:)
 real(dp) :: acell0(3),acell_next(3),angle(3),gmet(3,3),gprimd(3,3),mttk_alc(3)
 real(dp) :: mttk_alc2(3),mttk_blc(3),mttk_psh(3),mttk_tv(3,3),mttk_ubox(3,3)
 real(dp) :: mttk_uu(3),mttk_uv(3),mttk_veig(3),mttk_vt(3,3),rmet(3,3)
 real(dp) :: rprim_next(3,3),rprimd(3,3),rprimd0(3,3),rprimd_next(3,3)
 real(dp) :: strtarget(6),work(lwork)
 real(dp),allocatable :: binose(:,:),cinose(:,:),fac_type(:),fcart_m(:,:)
 real(dp),allocatable :: fcart_mold(:,:),finose(:,:),fred_corrected(:,:)
 real(dp),allocatable :: hessin(:,:),hnose(:,:),lang_force(:,:),max_perm(:)
 real(dp),allocatable :: pot_perm(:),ran_force(:,:),vel_nexthalf(:,:)
 real(dp),allocatable :: vel_prevhalf(:,:),vel_temp(:,:),vin(:),vin_next(:)
 real(dp),allocatable :: vin_prev(:),vonose(:,:),vout(:),xcart(:,:)
 real(dp),allocatable :: xcart_next(:,:),xcart_prev(:,:),xred_next(:,:)
 real(dp),allocatable :: xred_prev(:,:)

!************************************************************************
!Beginning of executable session
!***************************************************************************

!XG020711 : this line, tentatively added when setting up v3.4, is incorrect
!nxfh is initialized in the calling routine
!nxfh=0
 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 prtvol=dtset%prtvol
 if(prtvol==-level)then
  write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' moldyn : enter '
  call wrtout(06,message,'COLL')
 end if

 ipos=0
 ntime=dtset%ntime

 optcell=dtset%optcell
 ionmov=dtset%ionmov

!dtion=time step for molecular dynamics in atomic time units
!(1 atomic time unit=2.418884e-17 seconds)
 dtion=dtset%dtion
 noseinert=dtset%noseinert
 strtarget(1:6)=dtset%strtarget(1:6)

 ndim=3*dtset%natom
 if(optcell==1 .or. optcell==4 .or. optcell==5 .or. optcell==6)ndim=ndim+1
 if(optcell==2 .or. optcell==3)ndim=ndim+6
 if(optcell==7 .or. optcell==8 .or. optcell==9)ndim=ndim+3

 allocate(binose(3,dtset%natom),cinose(3,dtset%natom),hnose(3,dtset%natom))
 allocate(finose(3,dtset%natom),vonose(3,dtset%natom),vel_temp(3,dtset%natom))
 allocate(fcart_m(3,dtset%natom),fcart_mold(3,dtset%natom),lang_force(3,dtset%natom))
 allocate(fac_type(psps%ntypat),natom_type(psps%ntypat),ran_force(3,dtset%natom))
 allocate(fred_corrected(3,dtset%natom),stopped(dtset%natom))
 allocate(pot_perm(dtset%natom),max_perm(psps%ntypat),imax_perm(psps%ntypat))
 allocate(vel_nexthalf(3,dtset%natom),vel_prevhalf(3,dtset%natom))
 allocate(vin(ndim),vin_next(ndim))
 allocate(vin_prev(ndim))
 allocate(vout(ndim))
 allocate(xcart(3,dtset%natom),xcart_next(3,dtset%natom))
 allocate(xcart_prev(3,dtset%natom))
 allocate(xred_next(3,dtset%natom),xred_prev(3,dtset%natom))

!Compute dimensional primitive translations rprimd, then metric tensor gmet
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Save initial values
 acell0(:)=acell(:)
 rprimd0(:,:)=rprimd(:,:)
 ucvol0=ucvol

!Initialize input vectors : first vin, then vout
 option=1
 call xfpack(acell,acell0,results_gs%fred,dtset%natom,ndim,&
& dtset%nsym,optcell,option,rprim,rprimd0,&
& strtarget,results_gs%strten,dtset%symrel,ucvol,ucvol0,vin,vout,xred)
 option=3
 call xfpack(acell,acell0,results_gs%fred,dtset%natom,ndim,&
& dtset%nsym,optcell,option,rprim,rprimd0,&
& strtarget,results_gs%strten,dtset%symrel,ucvol,ucvol0,vin,vout,xred)

!Here, set up the matrix of transformation between forces and
!acceleration. Masses must be included here.
!Beside this feature, one could define
!a preconditioner, in which case it should
!be the inverse hessian, like in Broyden. This explains the
!name chosen for this transformation matrix. This would allow
!to find easily the optimal geometry with ionmov=7.
!The default, now implemented, corresponds to the identity matrix
!in cartesian coordinates, which makes use of metric tensor gmet
!in reduced coordinates.
 allocate(hessin(ndim,ndim))
 hessin(:,:)=0.0_dp
 do iatom=1,dtset%natom
  do idir1=1,3
   do idir2=1,3
!   Warning : implemented in reduced coordinates
    if ( dtset%iatfix(idir1,iatom) ==0 .and. dtset%iatfix(idir2,iatom) ==0 )then
     hessin(idir1+3*(iatom-1),idir2+3*(iatom-1))=&
&     gmet(idir1,idir2)/amass(iatom)
    end if
   end do
  end do
 end do
 if(optcell/=0)then
! These values might lead to too large changes in some cases ...
! No "mass" is included here
  diag=dtset%strprecon*30.0_dp/ucvol
  if(optcell==1)diag=diag/3.0_dp
  do idim=3*dtset%natom+1,ndim
   hessin(idim,idim)=diag
  end do
 end if

!-----------------------------------------------------------------------
!
!Iterative procedure (main loop)
!

 do itime=0,ntime

  if(itime/=0)then
   call status(itime,dtfil%filstat,iexit,level,'loop itime    ')
  end if

! Check whether exiting was required by the user.
! If found then beat a hasty exit from loop on itime
  openexit=1 ; if(dtset%chkexit==0) openexit=0
  call chkexi(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)
  if (iexit/=0) exit

  write(message, '(a,a,i4,a)' ) ch10,' MOLDYN STEP NUMBER ',itime,&
&  '  ------------------------------------------------------'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')

  if (ionmov==8.or.ionmov==9.or.ionmov==13) then
!  The temperature is linear between initial and final values
!  It is here converted from Kelvin to Hartree (kb_HaK)
   ktemp=(dtset%mditemp+((dtset%mdftemp-dtset%mditemp)/dble(ntime))*itime)*kb_HaK
  end if

! If not initialisation time step
  if(itime>0)then

!  Shift the data
   etotal_prev=results_gs%etotal
   acell(:)=acell_next(:)
   rprim(:,:)=rprim_next(:,:)
   rprimd(:,:)=rprimd_next(:,:)
   ucvol=ucvol_next
   vel_prevhalf(:,:)=vel_nexthalf(:,:)
   vin_prev(:)=vin(:)
   vin(:)=vin_next(:)
   xcart_prev(:,:)=xcart(:,:)
   xcart(:,:)=xcart_next(:,:)
   xred(:,:)=xred_next(:,:)
  else

!  Get xred, and eventually acell, rprim and rprimd, from current vin
   option=2
   call xfpack(acell,acell0,results_gs%fred,dtset%natom,ndim,&
&   dtset%nsym,optcell,option,rprim,rprimd0,&
&   strtarget,results_gs%strten,dtset%symrel,ucvol,ucvol0,vin,vout,xred)

!  End condition on itime
  end if

  call status(itime,dtfil%filstat,iexit,level,'call scfcv    ')

  if(optcell/=0)then

   call mkrdim(acell,rprim,rprimd)
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!  Write, but only to log file
   write(message, '(a,a,4(a,3es18.10,a),a,es18.10,a,a)' )&
&   ' Unit cell characteristics (before scfcv) :',ch10,&
&   '  acell=',acell(1:3),ch10,&
&   '  rprim=',rprim(1:3,1),ch10,&
&   '        ',rprim(1:3,2),ch10,&
&   '        ',rprim(1:3,3),ch10,&
&   '  ucvol=',ucvol,' Bohr^3',ch10
   call wrtout(6   ,message,'COLL')

  end if

! If metric has changed since the initialization, update the Ylm's
  if (optcell/=0.and.psps%useylm==1.and.itime>0)then
   call status(0,dtfil%filstat,iexit,level,'call initylmg ')
   option=0;if (dtset%iscf>0) option=1
   call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
&   npwarr,dtset%nsppol,option,rprimd,dtfil%unkg,dtfil%unylm,ylm,ylmgr)
  end if

  if((ionmov/=12.and.ionmov/=13).or.&
&  (ionmov==12.and.itime==0).or.(ionmov==13.and.itime==0)&
&  ) then
!  Compute LDA forces (big loop)
   iapp=-1
   if(itime>0)iapp=itime
   call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&
&   dtset,ecore,eigen,hdr,iapp,indsym,initialized,&
&   irrzon,kg,mpi_enreg,&
&   nattyp,nfftf,npwarr,nspinor,occ,&
&   pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&   phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,&
&   scf_history,symrec,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)
  end if

  call status(itime,dtfil%filstat,iexit,level,'call prtxvf   ')

! Output of acell and/or rprim ( and angles ! - should become a routine later)
  if(optcell/=0)then
   angle(1)=acos(rmet(2,3)/sqrt(rmet(2,2)*rmet(3,3)))/two_pi*360.0
   angle(2)=acos(rmet(1,3)/sqrt(rmet(1,1)*rmet(3,3)))/two_pi*360.0
   angle(3)=acos(rmet(1,2)/sqrt(rmet(1,1)*rmet(2,2)))/two_pi*360.0
   write(message, '(a,a,4(a,3es18.10,a))' )&
&   ' Unit cell characteristics :',ch10,&
&   '  acell=',acell(1:3),ch10,&
&   '  rprim=',rprim(1:3,1),ch10,&
&   '        ',rprim(1:3,2),ch10,&
&   '        ',rprim(1:3,3)
   call wrtout(ab_out,message,'COLL')
   call wrtout(6   ,message,'COLL')
   write(message, '(a,es18.10,a,a,a,3es18.10,a,a,a,3f13.8,a)' )&
&   '  ucvol=',ucvol,' Bohr^3',ch10,&
&   '  lengths=',sqrt(rmet(1,1)),sqrt(rmet(2,2)),sqrt(rmet(3,3)),' Bohr',&
&   ch10,'  angles (23,13,12)=',angle(1:3),' degrees'
   call wrtout(ab_out,message,'COLL')
   call wrtout(6   ,message,'COLL')
  end if

! Get rid off mean force on whole unit cell
  do idir=1,3
   favg=sum(results_gs%fred(idir,:))/dble(dtset%natom)
   fred_corrected(idir,:)=results_gs%fred(idir,:)-favg
   if(dtset%jellslab/=0.and.idir==3) fred_corrected(idir,:)=results_gs%fred(idir,:)
  end do

! Update xfhist
  nxfh=nxfh+1
  xfhist(:,1:dtset%natom,1,nxfh)=xred(:,:)
  xfhist(:,dtset%natom+1,1,nxfh)=acell(:)
  xfhist(:,dtset%natom+2:dtset%natom+4,1,nxfh)=rprim(1:3,1:3)
  xfhist(:,1:dtset%natom,2,nxfh)=fred_corrected(:,:)
  xfhist(:,dtset%natom+2,2,nxfh)=results_gs%strten(1:3)
  xfhist(:,dtset%natom+3,2,nxfh)=results_gs%strten(4:6)

! Store computed gradient in vout
  option=3
  call xfpack(acell,acell0,fred_corrected,&
&  dtset%natom,ndim,dtset%nsym,optcell,option,rprim,rprimd0,&
&  strtarget,results_gs%strten,dtset%symrel,ucvol,ucvol0,vin,vout,xred)

! ####### Test case ionmov=9  Langevin dynamics  ##########
  if (ionmov==9) then
   acell_next(:)=acell(:)
   ucvol_next=ucvol
   rprim_next(:,:)=rprim(:,:)
   rprimd_next(:,:)=rprimd(:,:)

   if(itime==0)then

!   Compute twice the kinetic energy of the system, called v2nose
    v2nose=0.0_dp
    do iatom=1,dtset%natom
     do idim=1,3
      v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
     end do
    end do

!   If there is no kinetic energy, use the
    if (v2nose<=v2tol) then
     v2nose=0.0_dp
     do iatom=1,dtset%natom
      do idim=1,3
!      uniformrandom returns a uniform random deviate between 0.0 and 1.0
!      if it were always 0 or 1, then the following expression
!      would give the requested temperature
       vel(idim,iatom)=(1.0_dp-2.0_dp*uniformrandom(idum))*sqrt(ktemp/amass(iatom))
!      Recompute v2nose
       v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
      end do
     end do
    end if

!   Now, rescale the velocities to give the proper temperature
    rescale_vel=sqrt(3.0_dp*dtset%natom*(dtset%mditemp)*kb_HaK/v2nose)
    vel(:,:)=vel(:,:)*rescale_vel
!   Recompute v2nose with the rescaled velocities
    v2nose=0.0_dp
    do iatom=1,dtset%natom
     do idim=1,3
      v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
     end do
    end do
    write(message, '(a)' )&
&    ' Rescaling or initializing velocities to initial temperature'
    call wrtout(ab_out,message,'COLL')
    call wrtout(6   ,message,'COLL')
    write(message, '(a,D12.5,a,D12.5)' )&
&    ' ---  Scaling factor : ',rescale_vel,' Asked T (K) ',dtset%mditemp
    call wrtout(ab_out,message,'COLL')
    call wrtout(6   ,message,'COLL')
    write(message, '(a,D12.5)' )&
&    ' ---  Effective temperature',v2nose/3.0_dp/(kb_HaK*dtset%natom)
    call wrtout(ab_out,message,'COLL')
    call wrtout(6   ,message,'COLL')
!   end if itime==0
   end if

!  This section is devoted to the optional atom permutation (JYR 001114)
!  Two input variables are needed
!  dtset%delayperm : is the interval (in time steps) at which
!  atoms are tentatively permuted
!  default value could be 0
!  dtset%signperm  : is the type of bias for the permutation
!  +1  to favor alternation of species
!  -1  to favor segregation

!  Force no permutation at initial step
   if (itime/=0 .and. dtset%delayperm/=0 ) then
    if (mod(itime,dtset%delayperm)==0) then
!    Try commutation of atoms.
     write(message, '(a)')' Attempt of commutation '
     call wrtout(ab_out,message,'COLL')
     call wrtout(6   ,message,'COLL')
!    Compute a 'permutation potential'
     do iatom=1,dtset%natom
      pot_perm(iatom)=0.0_dp
      do iatom1=1,dtset%natom
       if (iatom1.ne.iatom) then
        distx=xcart(1,iatom)-xcart(1,iatom1)
        distx=distx-acell(1)*nint(distx/acell(1))
        disty=xcart(2,iatom)-xcart(2,iatom1)
        disty=disty-acell(2)*nint(disty/acell(2))
        distz=xcart(3,iatom)-xcart(3,iatom1)
        distz=distz-acell(3)*nint(distz/acell(3))
!       Here we count each atom below 2 angstr as 1, could be customized
        dist=sqrt(distx*distx+disty*disty+distz*distz)/3.7807
        if (dtset%typat(iatom).ne.dtset%typat(iatom1)) then
         mcfac=-1
        else
         mcfac=1
        end if
        if (dist<1.0_dp)  dist=1.0_dp
        pot_perm(iatom)=pot_perm(iatom)+mcfac*(dtset%signperm)*1.0_dp&
&        /exp(log(dist)*6.0_dp)
       end if
      end do
     end do
     write(message, '(a,10f12.5)' )' Perm_pot ',&
&     (pot_perm(iatom1),iatom1=1,dtset%natom)
     call wrtout(ab_out,message,'COLL')
     call wrtout(6   ,message,'COLL')

!    Find the two atoms, of different types, with the highest perm_pot
     max_perm(:)=-1.0d9
     do iatom=1,dtset%natom
      if (pot_perm(iatom) > max_perm(dtset%typat(iatom))) then
       max_perm(dtset%typat(iatom))=pot_perm(iatom)
       imax_perm(dtset%typat(iatom))=iatom
      end if
     end do
!    DEBUG
!    write(message, '(a,10f12.5)' )' max_Perm ',&
!    &      (max_perm(itypat),itypat=1,psps%ntypat)
!    call wrtout(ab_out,message,'COLL')
!    call wrtout(6   ,message,'COLL')
!    write(message, '(a,10i12)' )' imax_Perm ',&
!    &      (imax_perm(itypat),itypat=1,psps%ntypat)
!    call wrtout(ab_out,message,'COLL')
!    call wrtout(6   ,message,'COLL')
!    ENDDEBUG

!    Loop and keep the 2 largest values
     if (max_perm(1)>max_perm(2)) then
      maxp1=max_perm(1)
      maxp2=max_perm(2)
      iatom1=imax_perm(1)
      iatom2=imax_perm(2)
     else
      maxp1=max_perm(2)
      maxp2=max_perm(1)
      iatom1=imax_perm(2)
      iatom2=imax_perm(1)
     end if

     do itypat=3,psps%ntypat
      if (max_perm(itypat)>maxp1) then
       maxp2=maxp1
       iatom2=iatom1
       maxp1=max_perm(itypat)
       iatom1=imax_perm(itypat)
      else if (max_perm(itypat)>maxp2) then
       maxp2=max_perm(itypat)
       iatom2=imax_perm(itypat)
      end if
     end do
     write(message, '(2(a,i5))' )' Will commute atom...',iatom1,'...of type ',&
&     dtset%typat(iatom1)
     call wrtout(ab_out,message,'COLL')
     call wrtout(6   ,message,'COLL')
     write(message, '(2(a,i5))' )'         with atom...',iatom2,'...of type ',&
&     dtset%typat(iatom2)
     call wrtout(ab_out,message,'COLL')
     call wrtout(6   ,message,'COLL')

!    Commute the atoms positions
     distx=xcart(1,iatom1)
     disty=xcart(2,iatom1)
     distz=xcart(3,iatom1)
     xcart(1,iatom1)=xcart(1,iatom2)
     xcart(2,iatom1)=xcart(2,iatom2)
     xcart(3,iatom1)=xcart(3,iatom2)
     xcart(1,iatom2)=distx
     xcart(2,iatom2)=disty
     xcart(3,iatom2)=distz
!    Convert back to xred (reduced coordinates)
     call xredxcart(dtset%natom,-1,rprimd,xcart,xred)

!    Store current total energy
     etotal_temp=results_gs%etotal

!    Compute LDA forces (big loop)
     iapp=-1
     if(itime>0)iapp=itime
     call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&
&     dtset,ecore,eigen,hdr,iapp,indsym,initialized,&
&     irrzon,kg,mpi_enreg,&
&     nattyp,nfftf,npwarr,nspinor,occ,&
&     pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&     phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,&
&     scf_history,symrec,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)

     if (results_gs%etotal > etotal_temp) then

!     Discard the changes
      distx=xcart(1,iatom1)
      disty=xcart(2,iatom1)
      distz=xcart(3,iatom1)
      xcart(1,iatom1)=xcart(1,iatom2)
      xcart(2,iatom1)=xcart(2,iatom2)
      xcart(3,iatom1)=xcart(3,iatom2)
      xcart(1,iatom2)=distx
      xcart(2,iatom2)=disty
      xcart(3,iatom2)=distz

!     Convert back to xred (reduced coordinates)
      call xredxcart(dtset%natom,-1,rprimd,xcart,xred)
      write(message, '(a)' )' Commutation unsuccessful, recomputing the forces'
      call wrtout(ab_out,message,'COLL')
      call wrtout(6   ,message,'COLL')

!     And recompute the forces
      iapp=-1
      if(itime>0)iapp=itime
      call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&
&      dtset,ecore,eigen,hdr,iapp,indsym,initialized,&
&      irrzon,kg,mpi_enreg,&
&      nattyp,nfftf,npwarr,nspinor,occ,&
&      pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&      phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,&
&      scf_history,symrec,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)

     else

      write(message, '(a)')' Commutation successful ! Going on'
      call wrtout(ab_out,message,'COLL')
      call wrtout(6   ,message,'COLL')

!     Get rid of mean force on whole unit cell, but only if no generalized
!     constraints are in effect
      if(dtset%nconeq==0)then
       do idir=1,3
        favg=sum(results_gs%fred(idir,:))/dble(dtset%natom)
        fred_corrected(idir,:)=results_gs%fred(idir,:)-favg
        if(dtset%jellslab/=0.and.idir==3) fred_corrected(idir,:)=results_gs%fred(idir,:)
       end do
      else
       fred_corrected(:,:)=results_gs%fred(:,:)
      end if

!     Update xfhist
      xfhist(:,1:dtset%natom,1,nxfh)=xred(:,:)
      xfhist(:,dtset%natom+1,1,nxfh)=acell(:)
      xfhist(:,dtset%natom+2:dtset%natom+4,1,nxfh)=rprim(1:3,1:3)
      xfhist(:,1:dtset%natom,2,nxfh)=fred_corrected(:,:)
      xfhist(:,dtset%natom+2,2,nxfh)=results_gs%strten(1:3)
      xfhist(:,dtset%natom+3,2,nxfh)=results_gs%strten(4:6)

!     Store computed gradient in vout
      option=3
      call xfpack(acell,acell0,fred_corrected,&
&      dtset%natom,ndim,dtset%nsym,optcell,option,rprim,rprimd0,&
&      strtarget,results_gs%strten,dtset%symrel,ucvol,ucvol0,vin,vout,xred)

     end if


    end if ! if(mod(itime,dtset%delayperm)==0)
   end if ! if(itime/=0 .or. dtset%delayperm/=0)
!  End of the commutation section

!  Specific to Langevin dynamics
!  Initialize an array of random forces
!  No random force at itime=0
!  if (itime==0) then
   if (itime<0) then

    ran_force(:,:)=0.0_dp

   else

    do iatom=1,dtset%natom
!    sig_gauss is the std deviation of the random distribution
     sig_gauss=sqrt(2.0_dp*(dtset%friction)*amass(iatom)*ktemp)
     do idim=1,3
      delxi=2.0_dp
      do while (delxi >= 1.0_dp)
       ran_num1=2.0_dp*uniformrandom(idum)-1.0_dp
       ran_num2=2.0_dp*uniformrandom(idum)-1.0_dp
       delxi=ran_num1*ran_num1+ran_num2*ran_num2
      end do
      ran_force(idim,iatom)=ran_num1*sqrt(-2.0_dp*log(delxi)/delxi)&
&      *sig_gauss/sqrt(dtion)

     end do
    end do
!   DEBUG
!   The distribution should be gaussian
!   delxi=0.0_dp
!   do iatom=1,dtset%natom
!   do idim=1,3
!   delxi=delxi+(ran_force(idim,iatom)*dtion)**2
!   end do
!   end do
!   delxi=delxi/(3.0_dp*dtset%natom)
!   write(message, '(2(a,es22.14))' )' variance =',delxi,'  asked =',&
!   &    2.0_dp*(dtset%friction)*amass(2)*ktemp*dtion
!   call wrtout(ab_out,message,'COLL')
!   call wrtout(6   ,message,'COLL')
!   ENDDEBUG
!   end if itime\=0

   end if

!  DEBUG
!  write(message, '(a)' )' after initializing ran_force'
!  call wrtout(ab_out,message,'COLL')
!  call wrtout(6   ,message,'COLL')
!  ENDDEBUG

   do iatom=1,dtset%natom
    do idim=1,3
     fcart_m(idim,iatom)=results_gs%fcart(idim,iatom)/amass(iatom)
     ran_force(idim,iatom)=ran_force(idim,iatom)/amass(iatom)
    end do
   end do
   lang_force(:,:)=ran_force(:,:)-(dtset%friction)*vel(:,:)+fcart_m(:,:)

!  DEBUG
!  write(message, '(a)' )'before verlet'
!  call wrtout(ab_out,message,'COLL')
!  call wrtout(6   ,message,'COLL')
!  ENDDEBUG

!  Compute next atomic coordinates using Verlet algorithm

!  Impose no change of acell, ucvol, rprim, and rprimd
   acell_next(:)=acell(:)
   ucvol_next=ucvol
   rprim_next(:,:)=rprim(:,:)
   rprimd_next(:,:)=rprimd(:,:)

!  Convert input xred (reduced coordinates) to xcart (cartesian)
   call xredxcart(dtset%natom,1,rprimd,xcart,xred)
!  Uses the velocity
!  
!  If an atom wants to cross the walls, velocity is reversed.
!  
   do iatom=1,dtset%natom
    do idim=1,3
     delxi=xcart(idim,iatom)+dtion*vel(idim,iatom)+ &
&     0.5_dp*dtion*dtion*lang_force(idim,iatom)
     if ( (delxi > (rprimd(idim,idim)+(dtset%mdwall)) ) .or. &
&     (delxi < - (dtset%mdwall)                   )       ) then
      vel(idim,iatom)=-vel(idim,iatom)
      delxi=xcart(idim,iatom)+dtion*vel(idim,iatom)+ &
&      0.5_dp*dtion*dtion*lang_force(idim,iatom)
     end if
     xcart_next(idim,iatom)=delxi
    end do
   end do
   xcart(:,:)=xcart_next(:,:)

!  Convert back to xred_next (reduced coordinates)
   call xredxcart(dtset%natom,-1,rprimd,xcart_next,xred_next)
   call xredxcart(dtset%natom,-1,rprimd,xcart,xred)

   if (itime==0) then
!   no old forces are available at first step
!   Simple update of the velocity
!   first compute vel_nexthalf for next steps
    vel(:,:)=vel(:,:)+dtion*lang_force(:,:)
   else
!   case itime /= 0 normal verlet integration
    vel(:,:)=vel(:,:)+0.5_dp*dtion*(fcart_mold(:,:)+lang_force(:,:))
   end if

!  Store 'current force' as 'old force'
   fcart_mold(:,:)=lang_force(:,:)

!  End of case ionmov =9

!  #####   case ionmov==8  Nose dynamics ###########
  else if (ionmov==8) then

   acell_next(:)=acell(:)
   ucvol_next=ucvol
   rprim_next(:,:)=rprim(:,:)
   rprimd_next(:,:)=rprimd(:,:)

   if(itime==0)then
    snose=0.0_dp
    xi_nose=0.0_dp
!   Compute twice the kinetic energy of the system, called v2nose
    v2nose=0.0_dp
    do iatom=1,dtset%natom
     do idim=1,3
      v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
     end do
    end do

!   If there is no kinetic energy, use a random initial velocity
    if (v2nose<=v2tol) then
     v2nose=0.0_dp
     do iatom=1,dtset%natom
      do idim=1,3
!      uniformrandom returns a uniform random deviate between 0.0 and 1.0
!      if it were always 0 or 1, then the following expression
!      would give the requested temperature
       vel(idim,iatom)=(1.0_dp-2.0_dp*uniformrandom(idum))*&
&       sqrt( (dtset%mditemp) * kb_HaK / amass(iatom) )
!      Recompute v2nose
       v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
      end do
     end do
    end if

!   Now, rescale the velocities to give the proper temperature
    rescale_vel=sqrt(3.0_dp*dtset%natom*(dtset%mditemp)*kb_HaK/v2nose)
    vel(:,:)=vel(:,:)*rescale_vel
!   Recompute v2nose with the rescaled velocities
    v2nose=0.0_dp
    do iatom=1,dtset%natom
     do idim=1,3
      v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
     end do
    end do
    write(message, '(a)' )&
&    ' Rescaling or initializing velocities to initial temperature'
    call wrtout(ab_out,message,'COLL')
    call wrtout(6   ,message,'COLL')
    write(message, '(2(a,es22.14))' )&
&    ' ---  Scaling factor : ',rescale_vel,' Asked T (K) ',dtset%mditemp
    call wrtout(ab_out,message,'COLL')
    call wrtout(6   ,message,'COLL')
    write(message, '(a,es22.14)' )&
&    ' ---  Effective temperature',v2nose/(3.0_dp*dtset%natom*kb_HaK)
    call wrtout(ab_out,message,'COLL')
    call wrtout(6   ,message,'COLL')
   end if


   do iatom=1,dtset%natom
    do idim=1,3
     fcart_m(idim,iatom)=results_gs%fcart(idim,iatom)/amass(iatom)
    end do
   end do

!  First step of velocity verlet algorithm
   gnose=3*dtset%natom

!  Convert input xred (reduced coordinates) to xcart (cartesian)
   call xredxcart(dtset%natom,1,rprimd,xcart,xred)

!  Calculate nose-hoover force on atoms
!  If first iteration, no old force are available, so use present forces
   if (itime==0) fcart_mold(:,:)=fcart_m(:,:)

   finose(:,:)=fcart_mold(:,:)-xi_nose*vel(:,:)
   xcart(:,:)=xcart(:,:)+dtion*(vel(:,:)+dtion*finose(:,:)/2.0_dp)

!  Convert back to xred (reduced coordinates)
   call xredxcart(dtset%natom,-1,rprimd,xcart,xred)

!  Calculate v2nose
   v2nose=0.0_dp
   do iatom=1,dtset%natom
    do idim=1,3
     v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
    end do
   end do
   vel(:,:)=vel(:,:)+dtion*finose(:,:)/2.0_dp

!  Update thermostat
   fsnose=(v2nose-gnose*ktemp)/noseinert
   snose=snose+dtion*(xi_nose+dtion*fsnose/2.0_dp)
   xi_nose=xi_nose+dtion*fsnose/2.0_dp

!  Second step of the velocity Verlet algorithm, uses the 'new forces'
!  Calculate v2nose
   v2nose=0.0_dp
   do iatom=1,dtset%natom
    do idim=1,3
     v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
    end do
   end do
   vel_temp(:,:)=vel(:,:)

   xin_nose=xi_nose

!  Start Newton-Raphson loop
   ready=.false.
   do while (.not.ready)
    xio=xin_nose
    delxi=0.0D0
    vonose(:,:)=vel_temp(:,:)
    hnose(:,:)=-dtion/2.0_dp*(fcart_m(:,:)-xio*vonose(:,:))-(vel(:,:)-vonose(:,:))
    do iatom=1,dtset%natom
     do idim=1,3
      binose(idim,iatom)=vonose(idim,iatom)*dtion/noseinert*amass(iatom) ! a verifier
      delxi=delxi+hnose(idim,iatom)*binose(idim,iatom)
     end do
    end do
    dnose=-(xio*dtion/2.0D0+1.0D0)
    delxi=delxi-dnose*((-v2nose+gnose*ktemp)*dtion/2.0_dp/ &
&    noseinert-(xi_nose-xio))
    delxi=delxi/(-dtion*dtion/2.0_dp*v2nose/noseinert+dnose)

!   hzeronose=-(xio-xi_nose-(v2nose-gnose*ktemp)*dtion/(2.0_dp*noseinert) )
!   cibinose=-v2nose*dtion*dtion/(2.0_dp*noseinert)
!   delxi=(delxi+hzeronose*dnose)/(dnose+cibinose)
!   DEBUG
!   write(message, '(a,es22.14)' )' after delxi',delxi
!   call wrtout(ab_out,message,'COLL')
!   call wrtout(6   ,message,'COLL')
!   ENDDEBUG
    v2nose=0.0_dp

    vel_temp(:,:)=vel_temp(:,:)+(hnose+dtion/2.0_dp*vonose(:,:)*delxi)/dnose
    do iatom=1,dtset%natom
     do idim=1,3
      v2nose=v2nose+vel_temp(idim,iatom)*vel_temp(idim,iatom)*amass(iatom)
     end do
    end do
!   New guess for xi
    xin_nose=xio+delxi


!   DEBUG
!   write(message, '(a,es22.14)' )' v2nose=',v2nose
!   call wrtout(ab_out,message,'COLL')
!   call wrtout(6   ,message,'COLL')
!   ENDDEBUG

    ready=.true.
!   Test for convergence
    iatom=0
    idim=1
    do while((iatom<=dtset%natom).and.(idim<=3).and.ready)
     iatom=iatom+1
     if (iatom>dtset%natom) then
      iatom=1
      idim=idim+1
     end if
     if ((iatom<=dtset%natom) .and.(idim<=3)) then
      if (abs(vel_temp(idim,iatom))<1.0d-50) vel_temp(idim,iatom)=1.0d-50
      if (abs((vel_temp(idim,iatom)-vonose(idim,iatom))/vel_temp(idim,iatom))&
&      >nosetol) ready=.false.
     else
      if (xin_nose<1.0d-50) xin_nose=1.0d-50
      if (abs((xin_nose-xio)/xin_nose)>nosetol) ready=.false.
     end if
    end do   ! end of while

!   Enddo ready
   end do

!  Update velocities to converged value
   vel(:,:)=vel_temp(:,:)
   write(message, '(a,es13.7)' )' converged velocities for T=',ktemp
   call wrtout(ab_out,message,'COLL')
   call wrtout(6   ,message,'COLL')

!  Update thermostat
   xi_nose=xin_nose
   xcart_next(:,:)=xcart(:,:)
!  Convert back to xred_next (reduced coordinates)
   call xredxcart(dtset%natom,-1,rprimd,xcart_next,xred_next)
!  Store 'new force' as 'old force'
   fcart_mold(:,:)=fcart_m(:,:)

!  End of case ionmov =8

!  #####   case ionmov==12 Isokinetic Ensemble ###########
  else if (ionmov==12) then

!  Application of Gauss' principle of least constraint according to Fei Zhang's algorithm (J. Chem. Phys. 106, 1997, p.6102)

   acell_next(:)=acell(:)
   ucvol_next=ucvol
   rprim_next(:,:)=rprim(:,:)
   rprimd_next(:,:)=rprimd(:,:)

!  v2gauss is twice the kinetic energy
   v2gauss=0.0_dp
   do iatom=1,dtset%natom
    do idim=1,3
     v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
    end do
   end do

!  If there is no kinetic energy
   if (v2gauss<=v2tol.and.itime==0) then
!   Maxwell-Boltzman distribution
    v2gauss=zero
    vtest=zero
    do iatom=1,dtset%natom
     do idim=1,3
      vel(idim,iatom)=sqrt(kb_HaK*dtset%mditemp/amass(iatom))*cos(two_pi*uniformrandom(idum))
      vel(idim,iatom)=vel(idim,iatom)*sqrt(-2._dp*log(uniformrandom(idum)))
     end do
    end do

!   Get rid of center-of-mass velocity
    s1=sum(amass(:))
    do idim=1,3
     s2=sum(amass(:)*vel(idim,:))
     vel(idim,:)=vel(idim,:)-s2/s1
    end do

!   Recompute v2gauss
    do iatom=1,dtset%natom
     do idim=1,3
      v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
      vtest=vtest+vel(idim,iatom)/(3._dp*dtset%natom)
     end do
    end do

!   Now rescale the velocities to give the exact temperature
    rescale_vel=sqrt(3._dp*dtset%natom*kb_HaK*dtset%mditemp/v2gauss)
    vel(:,:)=vel(:,:)*rescale_vel

!   Recompute v2gauss with the rescaled velocities
    v2gauss=zero
    do iatom=1,dtset%natom
     do idim=1,3
      v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
     end do
    end do

!   Compute the variance and print
    sigma2=(v2gauss/(3._dp*dtset%natom)-amass(1)*vtest**2)/kb_HaK

    write(message, '(a)' )&
&    ' Rescaling or initializing velocities to initial temperature'
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
    write(message, '(a,d12.5,a,D12.5)' )&
&    ' --- Scaling factor :',rescale_vel,' Asked T (K) ',dtset%mditemp
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
    write(message, '(a,d12.5,a,D12.5)' )&
&    ' --- Effective temperature',v2gauss/(3*dtset%natom*kb_HaK),' From variance', sigma2
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
   end if

   do iatom=1,dtset%natom
    do idim=1,3
     fcart_m(idim,iatom)=results_gs%fcart(idim,iatom)/amass(iatom)
    end do
   end do

!  Convert input xred (reduced coordinates) to xcart (cartesian)
   call xredxcart(dtset%natom,1,rprimd,xcart,xred)

   if(itime==0) then
    vel_nexthalf(:,:)=vel(:,:)
    xcart_next(:,:)=xcart(:,:)
    call xredxcart(dtset%natom,-1,rprimd,xcart_next,xred_next)
   else
!   Computation of vel_nexthalf (4.16 de Ref.1)
!   Computation of a and b (4.13 de Ref.1)
    a=0.0_dp
    b=0.0_dp
    do iatom=1,dtset%natom
     do idim=1,3
      a=a+fcart_m(idim,iatom)*vel(idim,iatom)*amass(iatom)
      b=b+fcart_m(idim,iatom)*fcart_m(idim,iatom)*amass(iatom)
     end do
    end do
    a=a/v2gauss
    b=b/v2gauss
!   Computation of s and scdot
    sqb=sqrt(b)
    as=sqb*dtion/2.
    s1=cosh(as)
    s2=sinh(as)
    s=a*(s1-1.)/b+s2/sqb
    scdot=a*s2/sqb+s1
    vel_nexthalf(:,:)=(vel(:,:)+fcart_m(:,:)*s)/scdot

!   Computation of the next positions
    xcart_next(:,:)=xcart(:,:)+vel_nexthalf(:,:)*dtion

!   Convert back to xred (reduced coordinates)
    call xredxcart(dtset%natom,-1,rprimd,xcart_next,xred_next)

!   Computation of the forces for the new positions
!   Compute LDA forces (big loop)
    iapp=-1
    if(itime>0)iapp=itime
    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
&    eigen,hdr,iapp,indsym,initialized,&
&    irrzon,kg,mpi_enreg,&
&    nattyp,nfftf,npwarr,nspinor,occ,&
&    pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&    phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,&
&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,dtset%natom
     do idim=1,3
      fcart_m(idim,iatom)=results_gs%fcart(idim,iatom)/amass(iatom)
     end do
    end do

!   Computation of vel(:,:) at the next positions
!   Computation of v2gauss
    v2gauss=0.0_dp
    do iatom=1,dtset%natom
     do idim=1,3
      v2gauss=v2gauss+vel_nexthalf(idim,iatom)*vel_nexthalf(idim,iatom)*amass(iatom)
     end do
    end do
!   Calcul de a et b (4.13 de Ref.1)
    a=0.0_dp
    b=0.0_dp
    do iatom=1,dtset%natom
     do idim=1,3
      a=a+fcart_m(idim,iatom)*vel_nexthalf(idim,iatom)*amass(iatom)
      b=b+fcart_m(idim,iatom)*fcart_m(idim,iatom)*amass(iatom)
     end do
    end do
    a=a/v2gauss
    b=b/v2gauss
!   Calcul de s et scdot
    sqb=sqrt(b)
    as=sqb*dtion/2.
    s1=cosh(as)
    s2=sinh(as)
    s=a*(s1-1.)/b+s2/sqb
    scdot=a*s2/sqb+s1
    vel(:,:)=(vel_nexthalf(:,:)+fcart_m(:,:)*s)/scdot
   end if

!  End of case ionmov = 12

!  #####   case ionmov==13 Reversible integrator of Martyna at al.###########
!  There are three sub cases according to the value of optcell
!  optcell=0 means isothermal, optcell==1:homogeneous cell fluctuations
!  optcell=2: full cell fluctuation in addition to temperature control.
  else if (ionmov==13) then
   write(6,*)'moldyn', dtset%nnos, dtset%qmass(:), dtset%bmass, dtset%vmass
   if(itime==0) then
    allocate(mttk_vars%glogs(dtset%nnos),mttk_vars%vlogs(dtset%nnos),&
&    mttk_vars%xlogs(dtset%nnos))
    mttk_vars%glogs(:)=zero; mttk_vars%vlogs(:)=zero; mttk_vars%xlogs(:)=zero
    mttk_vars%vboxg(:,:)=zero
    vlogv=zero
!   v2gauss is twice the kinetic energy
    v2gauss=0.0_dp
    do iatom=1,dtset%natom
     do idim=1,3
      v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
     end do
    end do

!   If there is no kinetic energy
    if (v2gauss<=v2tol.and.itime==0) then
!    Maxwell-Boltzman distribution
     v2gauss=zero
     vtest=zero
     do iatom=1,dtset%natom
      do idim=1,3
       vel(idim,iatom)=sqrt(kb_HaK*dtset%mditemp/amass(iatom))*cos(two_pi*uniformrandom(idum))
       vel(idim,iatom)=vel(idim,iatom)*sqrt(-2._dp*log(uniformrandom(idum)))
      end do
     end do

!    Get rid of center-of-mass velocity
     s1=sum(amass(:))
     do idim=1,3
      s2=sum(amass(:)*vel(idim,:))
      vel(idim,:)=vel(idim,:)-s2/s1
     end do

!    Recompute v2gauss
     do iatom=1,dtset%natom
      do idim=1,3
       v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
       vtest=vtest+vel(idim,iatom)/(3._dp*dtset%natom)
      end do
     end do

!    Now rescale the velocities to give the exact temperature
     rescale_vel=sqrt(3._dp*dtset%natom*kb_HaK*dtset%mditemp/v2gauss)
     vel(:,:)=vel(:,:)*rescale_vel

!    Recompute v2gauss with the rescaled velocities
     v2gauss=zero
     do iatom=1,dtset%natom
      do idim=1,3
       v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
      end do
     end do

!    Compute the variance and print
     sigma2=(v2gauss/(3._dp*dtset%natom)-amass(1)*vtest**2)/kb_HaK

     write(message, '(a)' )&
&     ' Rescaling or initializing velocities to initial temperature'
     call wrtout(ab_out,message,'COLL')
     call wrtout(6,message,'COLL')
     write(message, '(a,d12.5,a,D12.5)' )&
&     ' --- Scaling factor :',rescale_vel,' Asked T (K) ',dtset%mditemp
     call wrtout(ab_out,message,'COLL')
     call wrtout(6,message,'COLL')
     write(message, '(a,d12.5,a,D12.5)' )&
&     ' --- Effective temperature',v2gauss/(3*dtset%natom*kb_HaK),' From variance', sigma2
     call wrtout(ab_out,message,'COLL')
     call wrtout(6,message,'COLL')
    end if
   end if
!  XG070613 : Do not take away the following line , seems needed for the pathscale compiler
   write(6,*)'moldyn',mttk_vars%vboxg(:,:)
!  ##### sub  case optcell==0 Isothermal Ensemble ###########
   if(optcell==0) then
!   There is no evolution of cell
    acell_next(:)=acell(:)
    ucvol_next=ucvol
    rprim_next(:,:)=rprim(:,:)
    rprimd_next(:,:)=rprimd(:,:)
!   Update Thermostat variables and scale velocitie
    call isotemp(amass,dtion,dtset,ekin,ktemp,mttk_vars,vel)
!   Half velocity step
    do idim=1,3
     fcart_m(idim,:)=results_gs%fcart(idim,:)/amass(:)
    end do
    vel_nexthalf(:,:)=vel(:,:)+dtion/two*fcart_m(:,:)
!   New positions
!   Convert input xred (reduced coordinates) to xcart (cartesian)
    call xredxcart(dtset%natom,1,rprimd,xcart,xred)
    xcart_next(:,:)=xcart(:,:)+vel_nexthalf(:,:)*dtion
!   Convert back to xred (reduced coordinates)
    call xredxcart(dtset%natom,-1,rprimd,xcart_next,xred_next)
!   DEBUG
!   write(6,*)' xcart:'
!   do iatom=1,dtset%natom
!   write(6,*)' atom , position=',iatom,xcart_next(:,iatom)
!   enddo
!   ENDDEBUG
!   Computation of the forces for the new positions
!   Compute LDA forces (big loop)
    iapp=-1
    if(itime>0)iapp=itime
    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&
&    dtset,ecore,eigen,hdr,iapp,indsym,initialized,&
&    irrzon,kg,mpi_enreg,&
&    nattyp,nfftf,npwarr,nspinor,occ,&
&    pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&    phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,&
&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
!   Next Half velocity step
    do idim=1,3
     fcart_m(idim,:)=results_gs%fcart(idim,:)/amass(:)
    end do
    vel(:,:)=vel_nexthalf(:,:)+dtion/two*fcart_m(:,:)
!   Update Thermostat variables and velocity
    call isotemp(amass,dtion,dtset,ekin,ktemp,mttk_vars,vel)
    if(itime==0) massvol=ekin+results_gs%etotal
    write(6,*) 'conserved energy:',(ekin+results_gs%etotal)-massvol,ekin,results_gs%etotal
!   End of sub case optcell=0
!   ##### sub  case optcell==1 Isothermal-Isenthalpic Ensemble (homogeneous cell deformation)##########
   else if (optcell==1) then
!   Only homogeneous evolution of cell
!   Evolution of cell we keep rprim constant
    rprim_next(:,:)=rprim(:,:)
!   Update Thermostat variables and velocity
    call isopress(amass,dtion,dtset,ekin,ktemp,results_gs%strten,strtarget,ucvol,mttk_vars,vel,vlogv)
!   Half velocity step
    do idim=1,3
     fcart_m(idim,:)=results_gs%fcart(idim,:)/amass(:)
    end do
    vel_nexthalf(:,:)=vel(:,:)+dtion/two*fcart_m(:,:)
!   New positions
    mttk_aloc=exp(dtion/two*vlogv)
    mttk_aloc2=(vlogv*dtion/two)**2
    polysh=(((esh8*mttk_aloc2+esh6)*mttk_aloc2+esh4)*mttk_aloc2+esh2)*mttk_aloc2+one
    mttk_bloc=mttk_aloc*polysh*dtion
!   Convert input xred (reduced coordinates) to xcart (cartesian)
    call xredxcart(dtset%natom,1,rprimd,xcart,xred)
    xcart_next(:,:)=xcart(:,:)*mttk_aloc**2+vel_nexthalf(:,:)*mttk_bloc
!   Update the volume and related quantities
    acell_next(:)=acell(:)*exp(dtion*vlogv)
!   ucvol=ucvol*exp(dtion*vlogv)
    call mkrdim(acell_next,rprim,rprimd_next)
    call metric(gmet,gprimd,-1,rmet,rprimd_next,ucvol_next)
!   Convert back to xred (reduced coordinates)
    call xredxcart(dtset%natom,-1,rprimd_next,xcart_next,xred_next)
!   Computation of the forces for the new positions
!   Compute LDA forces (big loop)
!   If metric has changed since the initialization, update the Ylm's
    if (optcell/=0.and.psps%useylm==1.and.itime>0)then
     call status(0,dtfil%filstat,iexit,level,'call initylmg ')
     option=0;if (dtset%iscf>0) option=1
     call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
&     npwarr,dtset%nsppol,option,rprimd_next,dtfil%unkg,dtfil%unylm,ylm,ylmgr)
    end if
!   DEBUG
!   write(6,*)' xcart:   (2)'
!   do iatom=1,dtset%natom
!   write(6,*)' atom , position=',iatom,xcart_next(:,iatom)
!   enddo
!   ENDDEBUG
    iapp=-1
    if(itime>0)iapp=itime
    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&
&    dtset,ecore,eigen,hdr,iapp,indsym,initialized,&
&    irrzon,kg,mpi_enreg,&
&    nattyp,nfftf,npwarr,nspinor,occ,&
&    pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&    phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd_next,&
&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
!   Next Half velocity step
    do idim=1,3
     fcart_m(idim,:)=results_gs%fcart(idim,:)/amass(:)
    end do
    vel(:,:)=vel_nexthalf(:,:)+dtion/two*fcart_m(:,:)
!   Update Thermostat variables and velocity
    call isopress(amass,dtion,dtset,ekin,ktemp,results_gs%strten,strtarget,ucvol_next,mttk_vars,vel,vlogv)
    if(itime==0) massvol=ekin+results_gs%etotal
    write(6,*) 'conserved energy:',(ekin+results_gs%etotal)-massvol,ekin,results_gs%etotal
!   End of sub case optcell = 1
!   ##### sub  case optcell==2 Isothermal-Isenthalpic Ensemble (full cell deformation)##########
   else if (optcell==2) then
!   Fisrt half step for extended variables
    call isostress(amass,dtion,dtset,ekin,ktemp,results_gs%strten,strtarget,ucvol,vel,mttk_vars)
!   Half velocity step
    do idim=1,3
     fcart_m(idim,:)=results_gs%fcart(idim,:)/amass(:)
    end do
    vel_nexthalf(:,:)=vel(:,:)+dtion/two*fcart_m(:,:)
!   Convert input xred (reduced coordinates) to xcart (cartesian)
    call xredxcart(dtset%natom,1,rprimd,xcart,xred)
!   New positions
    mttk_vt(:,:)=mttk_vars%vboxg(:,:)
    call dsyev('V','U',3,mttk_vt,3,mttk_veig,work,lwork,ierr)
    mttk_tv(:,:)=transpose(mttk_vt)
    mttk_alc(:)=exp(dtion/two*mttk_veig(:))
    mttk_alc2(:)=(mttk_veig(:)*dtion/two)**2
    mttk_psh(:)=(((esh8*mttk_alc2(:)+esh6)*mttk_alc2(:)+esh4)*mttk_alc2(:)+esh2)*mttk_alc2(:)+one
    mttk_blc(:)=mttk_alc(:)*mttk_psh(:)*dtion
!   Update the positions
    do iatom=1,dtset%natom
     mttk_uu(:)=matmul(mttk_tv,xcart(:,iatom))
     mttk_uv(:)=matmul(mttk_tv,vel_nexthalf(:,iatom))
     mttk_uu(:)=mttk_uu(:)*mttk_alc(:)**2+mttk_uv(:)*mttk_blc(:)
     xcart_next(:,iatom)=matmul(mttk_vt,mttk_uu)
    end do
!   Update the box (rprimd and rprim)
    mttk_ubox(:,:)=matmul(mttk_tv,rprimd)
    do idim=1,3
     mttk_ubox(:,idim)=mttk_ubox(:,idim)*mttk_alc(:)**2
    end do
    rprimd_next(:,:)=matmul(mttk_vt,mttk_ubox)
    do idim=1,3
     rprim_next(idim,:)=rprimd_next(idim,:)/acell(:)
    end do
!   Update the volume
    call metric(gmet,gprimd,-1,rmet,rprimd_next,ucvol)
!   Convert back to xred (reduced coordinates)
    call xredxcart(dtset%natom,-1,rprimd_next,xcart_next,xred_next)
!   Computation of the forces for the new positions
!   If metric has changed since the initialization, update the Ylm's
    if (optcell/=0.and.psps%useylm==1.and.itime>0)then
     call status(0,dtfil%filstat,iexit,level,'call initylmg ')
     option=0;if (dtset%iscf>0) option=1
     call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
&     npwarr,dtset%nsppol,option,rprimd_next,dtfil%unkg,dtfil%unylm,ylm,ylmgr)
    end if

!   DEBUG
!   write(6,*)' xcart: (3)'
!   do iatom=1,dtset%natom
!   write(6,*)' atom , position=',iatom,xcart_next(:,iatom)
!   enddo
!   ENDDEBUG

!   Compute LDA forces (big loop)
    iapp=-1
    if(itime>0)iapp=itime
    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&
&    dtset,ecore,eigen,hdr,iapp,indsym,initialized,&
&    irrzon,kg,mpi_enreg,&
&    nattyp,nfftf,npwarr,nspinor,occ,&
&    pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&    phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd_next,&
&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
!   Next Half velocity step
    do idim=1,3
     fcart_m(idim,:)=results_gs%fcart(idim,:)/amass(:)
    end do
    vel(:,:)=vel_nexthalf(:,:)+dtion/two*fcart_m(:,:)
!   Next half step for extended variables
    call isostress(amass,dtion,dtset,ekin,ktemp,results_gs%strten,strtarget,ucvol,vel,mttk_vars)
    if(itime==0) massvol=ekin+results_gs%etotal
    write(6,*) 'conserved energy:',itime,(ekin+results_gs%etotal)-massvol,ekin,results_gs%etotal
!   Evolution of cell and volumr
    acell_next(:)=acell(:)
    ucvol_next=ucvol
!   End of sub case optcell=2
   else
    write(message, '(a,a,a,a,i12,a,a)' ) ch10,&
&    ' moldyn : BUG -',ch10,&
&    '  Disallowed value for optcell=',optcell,ch10,&
&    '  Allowed values with ionmov=13 : 0 to 2.'
    call wrtout(06,message,'COLL')
    call leave_new('COLL')
   end if
  else

!  ####### Begin case ionmov/=8 and ionmov/=9 ##########

!  Compute next atomic coordinates and cell parameters, using Verlet algorithm
!  First propagate the position, without acceleration
   if(itime/=0)then
    vin_next(:)=2*vin(:)-vin_prev(:)
    taylor=one
   else
!   Initialisation : no vin_prev is available, but the ionic velocity
!   is available, in cartesian coordinates
!   Convert input xred (reduced coordinates) to xcart (cartesian)
    call xredxcart(dtset%natom,1,rprimd,xcart,xred)
!   Uses the velocity
    xcart_next(:,:)=xcart(:,:)+dtion*vel(:,:)
!   Convert back to xred_next (reduced coordinates)
    call xredxcart(dtset%natom,-1,rprimd,xcart_next,xred_next)
!   Impose no change of acell, ucvol, rprim, and rprimd
    acell_next(:)=acell(:)
    ucvol_next=ucvol
    rprim_next(:,:)=rprim(:,:)
    rprimd_next(:,:)=rprimd(:,:)
!   Store all these next values in vin_next
    option=1
    call xfpack(acell_next,acell0,fred_corrected,&
&    dtset%natom,ndim,dtset%nsym,optcell,option,rprim_next,rprimd0,&
&    strtarget,results_gs%strten,dtset%symrel,ucvol_next,ucvol0,vin_next,vout,xred_next)
    taylor=half
   end if
!  Now, take into account the acceleration
   do idim=1,ndim
!   Note the minus sign: the forces are minus the gradients, contained in vout.
    vin_next(:)=vin_next(:)-dtion**2*hessin(:,idim)*vout(idim)*taylor
   end do
!  Implement fixing of atoms : put back old values for fixed components
   do iatom=1,dtset%natom
    do idir=1,3
!    Warning : implemented in reduced coordinates
     if (dtset%iatfix(idir,iatom) == 1) then
      vin_next(idir+(iatom-1)*3)=vin(idir+(iatom-1)*3)
     end if
    end do
   end do

!  Now, compute the velocity at the next half-step
!  Get xred_next, and eventually acell_next, ucvol_next, rprim_next and
!  rprimd_next, from vin_next
   option=2
   call xfpack(acell_next,acell0,results_gs%fred,dtset%natom,ndim,&
&   dtset%nsym,optcell,option,rprim_next,rprimd0,strtarget,results_gs%strten,&
&   dtset%symrel,ucvol_next,ucvol0,vin_next,vout,xred_next)
   if(optcell/=0)then
    call mkrdim(acell_next,rprim_next,rprimd_next)
    call metric(gmet,gprimd,-1,rmet,rprimd_next,ucvol_next)
   else
    rprimd_next(:,:)=rprimd(:,:)
   end if
!  Convert input xred_next (reduced coordinates) to xcart_next (cartesian)
   call xredxcart(dtset%natom,1,rprimd_next,xcart_next,xred_next)
!  Compute the velocity at half of the new step
   vel_nexthalf(:,:)=(xcart_next(:,:)-xcart(:,:))/dtion

!  If needed, compute the velocity at present position
   if(itime/=0)then
    vel(:,:)=(vel_nexthalf(:,:)+vel_prevhalf(:,:))*0.5_dp
   end if

!  End of case ionmov /=8 and /=9
  end if

! Compute the ionic kinetic energy (no cell shape kinetic energy yet)
  ekin=0.0_dp
  do iatom=1,dtset%natom
   do idir=1,3
!   Warning : the fixing of atomis is implemented in reduced
!   coordinates, so that this expression is wrong
    if (dtset%iatfix(idir,iatom) == 0) then
     ekin=ekin+0.5_dp*amass(iatom)*vel(idir,iatom)**2
    end if
   end do
  end do

! Output coordinates, forces and velocities
  prtvel=1
  call prtxvf(results_gs%fcart,dtset%iatfix,ab_out,dtset%natom,prtvel,vel,xcart)
  call prtxvf(results_gs%fcart,dtset%iatfix, 06 ,dtset%natom,prtvel,vel,xcart)

! Here, stop the atoms for which the scalar product of velocity
! and force is negative, and recompute the kinetic energy.
  if(ionmov==7)then
   stopped(:)=0
   do iatom=1,dtset%natom
    scprod=results_gs%fcart(1,iatom)*vel(1,iatom)+&
&    results_gs%fcart(2,iatom)*vel(2,iatom)+&
&    results_gs%fcart(3,iatom)*vel(3,iatom)
    if(scprod<0.0_dp .and. itime/=0)then
     stopped(iatom)=1
!    Shift the velocities of the previous half-step and current half-step,
!    so that the acceleration is correct but the present velocity vanishes.
     vel_prevhalf(:,iatom)=vel_prevhalf(:,iatom)-vel(:,iatom)
     vel_nexthalf(:,iatom)=vel_nexthalf(:,iatom)-vel(:,iatom)
     vel(:,iatom)=0.0_dp
     xcart_next(:,iatom)=xcart(:,iatom)+dtion*vel_nexthalf(:,iatom)
    end if
   end do

!  Establish a list of stopped atoms
   nstopped=sum(stopped(:))

   if(nstopped/=0)then
    write(message,'(a)') ' List of stopped atoms (ionmov=7) :'
    call wrtout(ab_out,message,'COLL')
    istopped=1
    do iatom=1,dtset%natom
     if(stopped(iatom)==1)then
      stopped(istopped)=iatom
      istopped=istopped+1
     end if
    end do
    do ii=1,nstopped,16
     write(message, '(16i4)' )stopped(ii:min(ii+15,nstopped))
     call wrtout(ab_out,message,'COLL')
    end do
!   Now, compute the corrected kinetic energy
!   Generate xred_next from xcart_next
    call xredxcart(dtset%natom,-1,rprimd_next,xcart_next,xred_next)
!   Store xred_next, and eventual acell_next and rprim_next in vin
    option=1
    call xfpack(acell_next,acell0,fred_corrected,&
&    dtset%natom,ndim,dtset%nsym,optcell,option,rprim_next,rprimd0,&
&    strtarget,results_gs%strten,dtset%symrel,ucvol_next,ucvol0,vin_next,vout,xred_next)
    ekin_corr=0.0_dp
    do iatom=1,dtset%natom
     do idir=1,3
!     Warning : the fixing of atomis is implemented in reduced
!     coordinates, so that this expression is wrong
      if (dtset%iatfix(idir,iatom) == 0) then
       ekin_corr=ekin_corr+0.5_dp*amass(iatom)*vel(idir,iatom)**2
      end if
     end do
    end do
!   End of test nstopped/=0
   end if

!  End of test ionmov==7
  end if

  if (ionmov==8) then
!  compute conserved quantity
   hamnose=results_gs%etotal+v2nose/2.0_dp+&
&   (xi_nose**2*noseinert)/2.0_dp+gnose*ktemp*snose
   write(message, '(a,i6,a,es22.14,a)' )&
&   ' At the end of Moldyn step',itime,', ham=',hamnose,' Ha.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(6,message,'COLL')
  end if

! Output total energy in a format that can be captured easily
  write(message, '(a,i6,a,es22.14,a)' )&
&  ' At the end of Moldyn step',itime,', POT.En.=',results_gs%etotal,' Ha.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(6,message,'COLL')
  write(message, '(a,es22.14,a)' )&
&  '                              KIN+POT.En.=',&
&  results_gs%etotal+ekin,' Ha.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(6,message,'COLL')
  if(ionmov==7 .and. nstopped/=0)then
   write(message, '(a,es22.14,a)' )&
&   '                    corrected KIN+POT.En.=',&
&   results_gs%etotal+ekin_corr,' Ha.'
   call wrtout(ab_out,message,'COLL')
  end if
  if(ionmov==9)then
   write(message, '(a,es22.14,2x,es22.14)' )&
&   '           TEMP         TKA =',&
&   ktemp/kb_HaK,ekin/(1.5_dp*dtset%natom*ktemp)
   call wrtout(ab_out,message,'COLL')
   call wrtout(6,message,'COLL')
   write(message, '(a,es22.14)' )&
&   '                       EKIN =',ekin
   call wrtout(ab_out,message,'COLL')
   call wrtout(6,message,'COLL')
  end if

! Writing outpout netcdf
! the output is being done every nctime time
  if ( mpi_enreg%me == 0) then
   nb1 = size(results_gs%strten)
   nbdir = size(xcart,1)
   if (dtset%nctime > 0)then
    if (itime == 0)then
     call  write_header_moldynnetcdf(dtfil, dtset, dtset%natom, nbdir, nb1 )
    end if
    if ( mod (itime, dtset%nctime ) == 0)then
!    at the ipos position vs time in file netcdf
     ipos = ipos +1
     call write_moldynvaluenetcdf(amass, ipos, dtfil, dtset, results_gs%etotal, ekin, &
&     dtset%natom, nbdir, nb1, xcart, vel, results_gs%strten, rprimd, ucvol )
     open(dtfil%unpos,file='POSABIN',status='replace',form='formatted')
     do iatom=1,dtset%natom
      if(iatom==1) then
       write(dtfil%unpos,'(a7,3d18.5)') 'xred  ',(xred_next(idim,iatom),idim=1,3)
      else
       write(dtfil%unpos,'(3d18.5)') (xred_next(idim,iatom),idim=1,3)
      end if
     end do
     do iatom=1,dtset%natom
      if(iatom==1) then
       write(dtfil%unpos,'(a7,3d18.5)') 'vel  ',(vel(idim,iatom),idim=1,3)
      else
       write(dtfil%unpos,'(3d18.5)') (vel(idim,iatom),idim=1,3)
      end if
     end do
     close(dtfil%unpos)
    end if
   end if
  end if

! Check whether forces and stresses are below tolerance; if so, exit
! from the itime loop
  if(ionmov==7)then
   iexit=0
   if(itime==ntime)iexit=1
   call fconv(results_gs%fcart,dtset%iatfix,iexit,itime,dtset%natom,ntime,&
&   optcell,dtset%strfact,strtarget,results_gs%strten,dtset%tolmxf)
   if (iexit/=0) exit
  end if

! Back to another iteration in the itime loop.
! Note that there are "exit" instructions inside the loop.
 end do

!-----------------------------------------------------------------------

 if(ionmov==8)then
  write(message, '(a,i6)')'Nb time steps',nxfh
  call wrtout(ab_out,message,'COLL')
  call wrtout(6,message,'COLL')
  do ierr=1,nxfh
   write(message, '(a,i6)')'step ',ierr
   call wrtout(ab_out,message,'COLL')
   call wrtout(6,message,'COLL')
   do iatom = 1, dtset%natom
    write(message, '(a,i6,3es22.14)')'atom ',iatom,xfhist(1:3,iatom,1,ierr)
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
   end do
  end do
 end if

 deallocate(binose,cinose,hnose)
 deallocate(finose,vonose)
 deallocate(ran_force, lang_force)
 deallocate(fac_type, STAT=ierr)
 deallocate(natom_type,STAT=ierr)
 deallocate(vel_temp,STAT=ierr)
 deallocate(fcart_m,fcart_mold,pot_perm, max_perm, imax_perm)
 deallocate(fred_corrected,hessin,stopped)
 deallocate(vel_nexthalf,vel_prevhalf)
 deallocate(vin,vin_next,vin_prev,vout)
 deallocate(xcart,xcart_next,xcart_prev,xred_next,xred_prev)
 if(ionmov==13) then
  deallocate(mttk_vars%glogs,mttk_vars%vlogs,mttk_vars%xlogs)
 end if
!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
  write(message,'(a1,a,a1,a,i1,a)') ch10,' moldyn : exit ',&
&  ch10,'  prtvol=-',level,', debugging mode => stop '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

end subroutine moldyn
!!***
