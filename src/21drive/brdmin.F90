!{\src2tex{textfont=tt}}
!!****f* ABINIT/brdmin
!! NAME
!! brdmin
!!
!! FUNCTION
!! Given a starting point xred that is a vector of length 3*natom
!! (reduced nuclei coordinates),
!! and unit cell parameters (acell and rprim) the
!! Broyden-Fletcher-Goldfarb-Shanno minimization is performed on the
!! total energy function, using its gradient (atomic forces and
!! stress : fred or fcart and stress) as calculated by the routine scfcv.
!! Some atoms can be kept fixed, while the optimization of unit cell
!! parameters is only performed if optcell/=0.
!! The convergence requirement on
!! the atomic forces, dtset%tolmxf,  allows an early exit.
!! Otherwise no more than dtset%ntime steps are performed.
!! Returned quantities are xred, and eventually acell and rprim (new ones!).
!! Could see Numerical Recipes (Fortran), 1986, page 307.
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, JCC, SE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cpus= cpu time limit in seconds
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem =number of k points which can fit in memory; set to 0 if use disk
!!   |  angular momentum for nonlocal pseudopotential
!!   | mpw=maximum dimensioned size of npw.
!!   | natom=number of atoms in unit cell
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   |      for the "coarse" grid (see NOTES below)
!!   | nkpt=number of k points.
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  ecore=core psp energy (part of total energy) (hartree)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mpi_enreg=informations about MPI parallelization
!!  mxfh=last dimension of the xfhist array
!!  nattyp(ntypat)= # atoms of each type.
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!       for the "fine" grid (see NOTES below)
!!  npwarr(nkpt)=number of planewaves in basis and boundary at this k point.
!!  nspinor=number of spinorial components of the wavefunctions
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
!!  vel(3,natom) is actually dummy in this routine
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points.
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!
!! SIDE EFFECTS
!!  Input/Output
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
!!  wffnew,wffnow=struct info for wf disk files.
!!  xfhist(3,natom+4,2,mxfh)=(x,f) history array,
!!                                 also includes acell, rprim and stress
!!  xred(3,natom)=reduced dimensionless atomic coordinates; updated on output
!!  xred_old(3,natom)=work space for old xred
!!
!! NOTES
!! there are many similarities between this routine
!! and moldyn.f, so that they might have to be maintained
!! together. Some common subroutines might be extracted from them.
!!
!! rhor and rhog are allocated as input with the nfftf size. But in the wavelet
!! case, this size can change and arrays will be reallocated inside the routine.
!! That's why they are pointers.
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
!!      are set equal to (nfft,ngfft,mgfft) in that case.!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      brdene,chkexi,fconv,hessupdt,initylmg,leave_new,metric,mkrdim,prtxvf
!!      scfcv,status,wrtout,xfpack,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine brdmin(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&
& dtset,ecore,eigen,hdr,indsym,initialized,irrzon,&
& kg,mpi_enreg,mxfh,&
& nattyp,nfftf,npwarr,nspinor,nxfh,occ,&
& pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,&
& rhog,rhor,rprim,scf_history,symrec,wffnew,wffnow,vel,wvl,&
& xfhist,xred,xred_old,ylm,ylmgr)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_13recipspace
 use interfaces_15common
 use interfaces_16geomoptim
 use interfaces_21drive, except_this_one => brdmin
!End of the abilint section

 implicit none

!Arguments ------------------------------------
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
 integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
!no_abirules
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 integer, intent(inout)  :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
 integer, intent(in)  :: kg(3,dtset%mpw*dtset%mkmem),nattyp(psps%ntypat)
 integer, intent(in)  :: npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
 integer, intent(inout)  :: symrec(3,3,dtset%nsym)
 real(dp), intent(inout)  :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
 real(dp), intent(inout)  :: acell(3),rprim(3,3)
 real(dp), intent(out)  :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 real(dp), intent(inout)  :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout)  :: phnons(2,dtset%nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
 real(dp), intent(in)  :: pwnsfac(2,pwind_alloc)
 real(dp), pointer  :: rhog(:,:),rhor(:,:)
 real(dp), intent(out)  :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout)  :: xfhist(3,dtset%natom+4,2,mxfh)
 real(dp), intent(inout)  :: xred(3,dtset%natom),xred_old(3,dtset%natom)
 real(dp), intent(in)  :: vel(3,dtset%natom)
 real(dp), intent(inout)  :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(inout)  :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(pawrad_type), intent(in)  :: pawrad(psps%ntypat*psps%usepaw)
 type(pawrhoij_type), intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawtab_type), intent(in)  :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=5
 integer :: counter,cycle_main,hess_ok,iapp,iatom,iatom1,iatom2,idim,idir,idir1
 integer :: idir2,iexit,ii,ios,itime,ixfh,jatom,jdir,jj,ndim,ndim0,openexit
 integer :: optcell,option,prtvel,prtvol,routine
 real(dp) :: DDOT,diag,etotal_prev,favg,ucvol,ucvol0
 character(len=32) :: statusOut
 character(len=500) :: message
!arrays
 real(dp) :: acell0(3),angle(3),gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp) :: rprimd0(3,3)
 real(dp),allocatable :: fred_corrected(:,:),hessin(:,:),vin(:),vin_prev(:)
 real(dp),allocatable :: vout(:),vout_prev(:),xcart(:,:)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 prtvol=dtset%prtvol
 if(prtvol==-level)then
  write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' brdmin : enter '
  call wrtout(06,message,'COLL')
 end if

 optcell=dtset%optcell

 ndim=3*dtset%natom
 if(optcell==1 .or. optcell==4 .or. optcell==5 .or. optcell==6)ndim=ndim+1
 if(optcell==2 .or. optcell==3)ndim=ndim+6
 if(optcell==7 .or. optcell==8 .or. optcell==9)ndim=ndim+3

 allocate(fred_corrected(3,dtset%natom),xcart(3,dtset%natom))
 allocate(vin(ndim),vout(ndim),vin_prev(ndim),vout_prev(ndim))
 allocate(hessin(ndim,ndim))

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
& dtset%strtarget,results_gs%strten,dtset%symrel,ucvol,ucvol0,vin,vout,xred)
 option=3
 call xfpack(acell,acell0,results_gs%fred,dtset%natom,ndim,&
& dtset%nsym,optcell,option,rprim,rprimd0,&
& dtset%strtarget,results_gs%strten,dtset%symrel,ucvol,ucvol0,vin,vout,xred)

!----------------------------------------------------------------------

!Initialise the Hessian matrix using gmet
 call hessinit(dtfil, dtset, hessin, gmet, ndim, ucvol)

!-----------------------------------------------------------------------
!
!Iterative procedure (main loop)

 do itime=0,dtset%ntime

  if(itime/=0)then
   call status(itime,dtfil%filstat,iexit,level,'loop itime    ')
  end if

! Check whether exiting was required by the user.
! If found then beat a hasty exit from loop on itime
  openexit=1 ; if (dtset%chkexit==0) openexit=0
  call chkexi(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)
  if (iexit/=0) exit

  write(message, '(a,a,i3,a)' ) ch10,' BROYDEN STEP NUMBER ',itime,&
&  '  ------------------------------------------------------'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')

  cycle_main=0

! If not initialisation time step
  if(itime>0)then

   if(dtset%ionmov==2 .or. itime==1)then

!   Previous cartesian coordinates
    etotal_prev=results_gs%etotal
    vin_prev(:)=vin(:)

!   New atomic cartesian coordinates are obtained from vin, hessin and vout
    do idim=1,ndim
     vin(:)=vin(:)-hessin(:,idim)*vout(idim)
    end do

!   Previous atomic forces
    vout_prev(:)=vout(:)

   else if(dtset%ionmov==3)then

!   Here the BFGS algorithm, modified to take into account the energy
    call brdene(results_gs%etotal,etotal_prev,hessin,&
&    ndim,vin,vin_prev,vout,vout_prev)

   end if

!  Implement fixing of atoms : put back old values for fixed components
   do iatom=1,dtset%natom
    do idir=1,3
!    Warning : implemented in reduced coordinates
     if ( dtset%iatfix(idir,iatom) == 1) then
      vin(idir+(iatom-1)*3)=vin_prev(idir+(iatom-1)*3)
     end if
    end do
   end do

   call status(itime,dtfil%filstat,iexit,level,'call scfcv    ')

!  If itime==0
  else

   if(nxfh==0)then
!   In case there is no previous history
    call status(itime,dtfil%filstat,iexit,level,'call scfcv_ini')
   else
!   In case there is a previous history : update hessian, then cycle.
    write(message, '(a,a)' ) ch10,&
&    ' Initialize atomic coordinates and forces from previous history.'
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')

!   Loop over previous time steps
    do ixfh=1,nxfh

!    For that time step, get new (x,f) from xfhist
     xred(:,:)     =xfhist(:,1:dtset%natom        ,1,ixfh)
     rprim(1:3,1:3)=xfhist(:,dtset%natom+2:dtset%natom+4,1,ixfh)
     acell(:)      =xfhist(:,dtset%natom+1        ,1,ixfh)
     results_gs%fred(:,:)     =xfhist(:,1:dtset%natom        ,2,ixfh)
!    This use of results_gs is unusual
     results_gs%strten(1:3)   =xfhist(:,dtset%natom+2        ,2,ixfh)
     results_gs%strten(4:6)   =xfhist(:,dtset%natom+3        ,2,ixfh)
!    Transfer it in vin, vout
     option=1
     call xfpack(acell,acell0,results_gs%fred,dtset%natom,ndim,&
&     dtset%nsym,optcell,option,rprim,rprimd0,&
&     dtset%strtarget,results_gs%strten,dtset%symrel,ucvol,ucvol0,vin,vout,xred)
     option=3
     call xfpack(acell,acell0,results_gs%fred,dtset%natom,ndim,&
&     dtset%nsym,optcell,option,rprim,rprimd0,&
&     dtset%strtarget,results_gs%strten,dtset%symrel,ucvol,ucvol0,vin,vout,xred)

!    Get old time step, if any, and update inverse hessian
     if(ixfh/=1)then
      xred(:,:)     =xfhist(:,1:dtset%natom        ,1,ixfh-1)
      rprim(1:3,1:3)=xfhist(:,dtset%natom+2:dtset%natom+4,1,ixfh-1)
      acell(:)      =xfhist(:,dtset%natom+1        ,1,ixfh-1)
      results_gs%fred(:,:)     =xfhist(:,1:dtset%natom        ,2,ixfh-1)
!     This use of results_gs is unusual
      results_gs%strten(1:3)   =xfhist(:,dtset%natom+2        ,2,ixfh-1)
      results_gs%strten(4:6)   =xfhist(:,dtset%natom+3        ,2,ixfh-1)
!     Tranfer it in vin_prev, vout_prev
      option=1
      call xfpack(acell,acell0,results_gs%fred,dtset%natom,ndim,&
&      dtset%nsym,optcell,option,rprim,rprimd0,&
&      dtset%strtarget,results_gs%strten,dtset%symrel,ucvol,ucvol0,vin_prev,vout_prev,xred)
      option=3
      call xfpack(acell,acell0,results_gs%fred,dtset%natom,ndim,&
&      dtset%nsym,optcell,option,rprim,rprimd0,&
&      dtset%strtarget,results_gs%strten,dtset%symrel,ucvol,ucvol0,vin_prev,vout_prev,xred)

      call hessupdt(hessin,dtset%iatfix,dtset%natom,ndim,vin,vin_prev,vout,vout_prev)

     end if

!    End loop over previous time steps
    end do

!   The hessian has been generated, as well as the latest vin and vout
!   so will cycle the main loop
    cycle_main=1

   end if

!  End condition on itime
  end if

! Get xred, and eventually acell, rprim and rprimd, from current vin
  option=2
  call xfpack(acell,acell0,results_gs%fred,dtset%natom,ndim,&
&  dtset%nsym,optcell,option,rprim,rprimd0,&
&  dtset%strtarget,results_gs%strten,dtset%symrel,ucvol,ucvol0,vin,vout,xred)

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

   write(message, '(a)' ) '  xred='
   call wrtout(6   ,message,'COLL')
   do iatom=1,dtset%natom
    write(message, '(3(es18.10,1x))' ) xred(1:3,iatom)
    call wrtout(6   ,message,'COLL')
   end do

  end if

! If metric has changed since the initialization, update the Ylm's
  if (optcell/=0.and.psps%useylm==1.and.itime>0)then
   call status(0,dtfil%filstat,iexit,level,'call initylmg ')
   option=0;if (dtset%iscf>0) option=1
   call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
&   npwarr,dtset%nsppol,option,rprimd,dtfil%unkg,dtfil%unylm,ylm,ylmgr)
  end if

! Convert input xred (reduced coordinates) to xcart (cartesian)
  call xredxcart(dtset%natom,1,rprimd,xcart,xred)

  prtvel=0
  call prtxvf(results_gs%fcart,dtset%iatfix, 06 ,dtset%natom,prtvel,vel,xcart)

  if(cycle_main==1)cycle

! DEBUG
! write(6,*)' brdmin : call scfcv '
! xred(1:3,1)=0.0_dp
! xred(1:3,2)=0.25_dp
! if(itime==0)acell(:)=10.10_dp
! if(itime==1)acell(:)=10.20_dp
! if(itime==2)acell(:)=10.30_dp
! call mkrdim(acell,rprim,rprimd)
! stop
! ENDDEBUG

! Compute LDA forces (big loop)
  iapp=-1
  if(itime>0)iapp=itime
  call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
&  eigen,hdr,iapp,indsym,initialized,&
&  irrzon,kg,mpi_enreg,&
&  nattyp,nfftf,npwarr,nspinor,occ,&
&  pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&  phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,&
&  scf_history,symrec,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)

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

! Again output coordinates and forces (not velocities) and total energy
  prtvel=0
  call prtxvf(results_gs%fcart,dtset%iatfix,ab_out,dtset%natom,prtvel,vel,xcart)
  call prtxvf(results_gs%fcart,dtset%iatfix, 06 ,dtset%natom,prtvel,vel,xcart)

! Output total energy in a format that can be captured easily
  write(message, '(a,i3,a,es22.14,a,a)' )&
&  ' At the end of Broyden step',itime,', total energy=',&
&  results_gs%etotal,' Ha.',ch10
  call wrtout(ab_out,message,'COLL')

! Get rid of mean force on whole unit cell, but only if no generalized
! constraints are in effect
  if(dtset%nconeq==0)then
   do idir=1,3
    favg=sum(results_gs%fred(idir,:))/dble(dtset%natom)
    fred_corrected(idir,:)=results_gs%fred(idir,:)-favg
    if(dtset%jellslab/=0.and.idir==3) fred_corrected(idir,:)=results_gs%fred(idir,:)
   end do
  else
   fred_corrected(:,:)=results_gs%fred(:,:)
  end if

! Update xfhist
  nxfh=nxfh+1
  xfhist(:,1:dtset%natom,1,nxfh)=xred(:,:)
  xfhist(:,dtset%natom+1,1,nxfh)=acell(:)
  xfhist(:,dtset%natom+2:dtset%natom+4,1,nxfh)=rprim(1:3,1:3)
  xfhist(:,1:dtset%natom,2,nxfh)=fred_corrected(:,:)
  xfhist(:,dtset%natom+2,2,nxfh)=results_gs%strten(1:3)
  xfhist(:,dtset%natom+3,2,nxfh)=results_gs%strten(4:6)

! Check whether forces and stresses are below tolerance; if so, exit
! from the itime loop
  iexit=0
  if(itime==dtset%ntime) then
   iexit=1
   statusOut = "Failed"
  else
   statusOut = "OK"
  end if
  call fconv(results_gs%fcart,dtset%iatfix,iexit,itime,dtset%natom,dtset%ntime,&
&  optcell,dtset%strfact,dtset%strtarget,results_gs%strten,dtset%tolmxf)
  if (iexit/=0) then
   exit
  end if

! Store xred, and eventual acell and rprim in vin
  option=1
  call xfpack(acell,acell0,fred_corrected,&
&  dtset%natom,ndim,dtset%nsym,optcell,option,rprim,rprimd0,&
&  dtset%strtarget,results_gs%strten,dtset%symrel,ucvol,ucvol0,vin,vout,xred)

! Store computed gradient in vout
  option=3
  call xfpack(acell,acell0,fred_corrected,&
&  dtset%natom,ndim,dtset%nsym,optcell,option,rprim,rprimd0,&
&  dtset%strtarget,results_gs%strten,dtset%symrel,ucvol,ucvol0,vin,vout,xred)

  if(itime>0)then
!  Update the hessian matrix, by taking into account the
!  current pair (x,f) and the previous one.
   call hessupdt(hessin,dtset%iatfix,dtset%natom,ndim,vin,vin_prev,vout,vout_prev)
  end if

! Back to another iteration in the itime loop.
! Note that there are "exit" instructions inside the loop.
! There is also a "cycle" instruction inside the loop.
 end do

!-----------------------------------------------------------------------

 call status(0,dtfil%filstat,iexit,level,'print hessian ')

!Print inverse hessian matrix "hessin" to log file
!for study or possible re-start use
 write(message, '(a,a,a,a)' ) ch10,&
& ' From Broyden minimization, approx. inverse',&
& ' Hessian(ndim,ndim):',ch10
 call wrtout(06,message,'COLL')
 do jj=1,ndim
  write(message, '(i6)' )jj
  call wrtout(06,message,'COLL')
  do ii=1,ndim,3
   write(message, '(3es22.14)' )hessin(ii:min(ii+2,ndim),jj)
   call wrtout(06,message,'COLL')
  end do
 end do

 deallocate(fred_corrected,hessin,vin,vin_prev,vout,vout_prev,xcart)

!XML output of the status
 if (mpi_enreg%me == 0 .and. dtset%outputXML == 1) then
  write(ab_xml_out, "(A)") '    <geometryMinimisation type="bfgs">'
  write(ab_xml_out, "(A,A,A)") '      <status cvState="', trim(statusOut) , &
&  '" stop-criterion="tolmxf" />'
  write(ab_xml_out, "(A)") '    </geometryMinimisation>'
 end if

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
  write(message,'(a1,a,a1,a,i1,a)') ch10,' brdmin : exit ',&
&  ch10,'  prtvol=-',level,', debugging mode => stop '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

end subroutine brdmin
!!***
