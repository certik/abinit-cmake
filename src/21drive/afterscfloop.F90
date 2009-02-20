!{\src2tex{textfont=tt}}
!!****f* ABINIT/afterscfloop
!! NAME
!! afterscfloop
!!
!! FUNCTION
!! Perform all calculations needed after the SCF loop, independent of the
!! call to scfcv (with or without atomic displacements), and exclusive
!! of print or write purposes, or deallocations.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=wavefunctions
!!                  (may be read from disk instead of input)
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecrpj)=<p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                                 and each |p_lmn> non-local projector
!!  cpus= cpu time limit in seconds
!!  deltae=change in energy between the previous and present SCF cycle
!!  dimcprj(natom*usecprj)=array of dimensions of array cprj
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs (see NOTES at beginning of scfcv)
!!   | mkmem=maximum number of k points in core memory
!!   | mpw = maximum number of plane waves
!!   | natom=number of atoms in cell
!!   | nfft=(effective) number of FFT grid points (for this processor) (see NOTES at beginning of scfcv)
!!   | nkpt=number of k points in Brillouin zone
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetries in space group
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  filapp character(len=fnlen)=generic output root name, with appendix
!!  filfft character(len=fnlen) =temporary FFT file, to be deleted
!!  grewtn(3,natom)=d(Ewald)/d(xred) (hartree)
!!  gsqcut=cutoff value on G**2 for (large) sphere inside FFT box.
!!                       gsqcut=(boxcut**2)*dtset%ecut/(2._dp*(Pi**2)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  indsym(4,nsym,natom)=index showing transformation of atom labels
!!                       under symmetry operations (computed in symatm)
!!  istep=number of the SCF iteration
!!  kg(3,mpw*mkmem)=reduced (integer) coordinates of G vecs in basis sphere
!!  kxc(nfftf,nkxc)=XC kernel
!!  mgfftc= - PAW only - maximum size of 1D FFTs for the "coarse" grid (see NOTES at beginning of scfcv)
!!  mgfftf= - PAW only - maximum size of 1D FFTs for the "fine" grid (see NOTES at beginning of scfcv)
!!  moved_atm_inside: if==1, the atoms are allowed to move.
!!  mpi_enreg=informations about MPI parallelization
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  nattyp(dtset%ntypat)=number of atoms of each type
!!  nfftf= - PAW only - number of FFT grid points for the "fine" grid (see NOTES at beginning of scfcv)
!!  ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!  ngfftf(18)= - PAW only - contain all needed information about 3D FFT  for the "fine" grid
!!  nhat(nfftf,nspden*psps%usepaw)= -PAW only- compensation density
!!  nkxc=dimension of kxc
!!  npwarr(nkpt)=number of planewaves in basis and on boundary for each k
!!  nvresid(nfftf,nspden)=array for the residual of the density/potential
!!  occ(mband*nkpt*nsppol)=occupancies of bands at various k points
!!  optres=0: the potential residual has been computed in scfcv
!!         1: the density residual has been computed in scfcv
!!  optxc=option for XC
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pel(3)=reduced coordinates of the electronic polarization (a. u.)
!!  pel_cg(3) = reduced coordinates of the electronic polarization (a. u.)
!!             computed in the SCF loop
!!  ph1df(2,3*(2*mgfftf+1)*natom)= - PAW only - 1-dim structure factor phases for the "fine" grid
!!      (see NOTES at beginning of scfcv)
!!  pion(3)=reduced coordinates of the ionic polarization (a. u.)
!!  prtfor=1 only if forces have to be printed (0 otherwise)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!  res2=density/potential residual (squared)
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points and spins
!!  residm=maximum value from resid array (except for nbdbuf highest bands)
!!  rhocore= (to be described)
!!  rhog(2,nfftf)=Fourier transform of total electron density (including compensation density in PAW)
!!  rhor(nfftf,nspden)=total electron density (including compensation density in PAW)
!!  rhore= (to be described)
!!  rhototp= (to be described)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  stress_needed=1 if stresses are needed, 0 otherwise
!!  strsxc(6)=xc correction to stress
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  tollist(12)=list of tolerances
!!  usecprj=1 if cprj datastructure has been allocated
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vhartr(nfftf)=Hartree potential
!!  vpsp(nfftf)=array for holding local psp
!!  vxc(nfftf,nspden)=exchange-correlation potential (hartree) in real space
!!  vxcavg=vxc average
!!  wffnow=unit number for current wf disk file
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!   (should be made a pure output quantity)
!!  xred_old(3,natom)=reduced dimensionless atomic coordinates, from input xred
!!  ==== if forces are required ====
!!   diffor=maximal absolute value of changes in the components of
!!          force between the input and the output.
!!   favg(3)=mean of the forces before correction for translational symmetry
!!   fcart(3,natom)=forces in cartesian coordinates (Ha/Bohr)
!!     at input, previous value of forces,
!!     at output, new value.
!!     Note : unlike fred, this array has been corrected by enforcing
!!     the translational symmetry, namely that the sum of force
!!     on all atoms is zero.
!!   fred(3,natom)=symmetrized grtn = d(etotal)/d(xred)
!!   gresid(3,natom)=forces due to the residual of the potential
!!   grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!   grxc(9+3*natom)=d(Exc)/d(xred) if core charges are used
!!   maxfor=maximal absolute value of the output array force.
!!   synlgr(3,natom)=symmetrized gradients of energy due to nonlocal contributions
!!  ==== if stress tensor is required ====
!!   strten(6)=components of the stress tensor (hartree/bohr^3) for the
!!    6 unique components of this symmetric 3x3 tensor:
!!    Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).
!!
!! SIDE EFFECTS
!! computed_forces=1 if forces have been computed, 0 otherwise
!! dtefield <type(efield_type)> = variables related to Berry phase
!!       and electric field calculations (see initberry.f).
!!       In case dtset%berryopt = 4, the overlap matrices computed
!!       in this routine are stored in dtefield%smat in order
!!       to be used in the electric field calculation.
!!  energies <type(energies_type)>=all part of total energy.
!!   | entropy(IN)=entropy due to the occupation number smearing (if metal)
!!   | e_localpsp(IN)=local psp energy (hartree)
!!   | e_eigenvalues(IN)=Sum of the eigenvalues - Band energy (Hartree)
!!   | e_ewald(IN)=Ewald energy (hartree)
!!   | e_hartree(IN)=Hartree part of total energy (hartree units)
!!   | e_corepsp(IN)=psp core-core energy
!!   | e_kinetic(IN)=kinetic energy part of total energy.
!!   | e_nonlocalpsp(IN)=nonlocal pseudopotential part of total energy.
!!   | e_xc(IN)=exchange-correlation energy (hartree)
!!   | e_xcdc(IN)=exchange-correlation double-counting energy (hartree)
!!   | e_paw(IN)=PAW spherical part energy
!!   | e_pawdc(IN)=PAW spherical part double-counting energy
!!   | e_elecfield(OUT)=the term of the energy functional that depends explicitely
!!   |                  on the electric field:  enefield = -ucvol*E*P
!! etotal=total energy, might be correct by improved polarization computation
!! forold(3,natom)=old forces
!! ===== if dtset%iprcch==3 .and. moved_atm_inside==1 =====
!!   ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phases (coarse grid)
!!   ph1df(2,3*(2*mgfftf+1)*natom)=1-dim structure factor phases (fine PAW grid)
!!  wvl <type(wvl_data)>=all wavelets data.
!!
!! NOTES
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      calc_lifetime,dtsetcopy,dtsetfree,elpolariz,forces,forstr,getph
!!      hamiltonianapplication,hartre1,hdr_update,last_orthon,leave_new,metric
!!      nhatgrid,pawmknhat,rhohxc,scprqt,spin_current,status,wrtout,wvl_mkrho
!!      wvl_newvtr,wvl_nl_gradient,wvl_tail_corrections,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine afterscfloop(atindx,atindx1,cg,computed_forces,cprj,cpus,&
& dimcprj,deltae,diffor,dtefield,dtfil,dtset,eigen,energies,etotal,&
& favg,fcart,filapp,filfft,forold,fred,gresid,grewtn,grhf,&
& grxc,gsqcut,hdr,indsym,&
& istep,kg,kxc,maxfor,mgfftc,mgfftf,&
& moved_atm_inside,mpi_enreg,&
& n3xccc,nattyp,&
& nfftf,ngfft,ngfftf,nhat,nkxc,npwarr,nvresid,&
& occ,optres,optxc,paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,pel,pel_cg,&
& ph1d,ph1df,pion,prtfor,psps,pwind,pwind_alloc,pwnsfac,res2,resid,residm,results_gs,&
& rhocore,rhog,rhor,rhore,rhototp,&
& rprimd,stress_needed,strsxc,strten,symrec,synlgr,tollist,usecprj,usexcnhat,&
& vhartr,vpsp,vxc,vxcavg,wffnow,wvl,xccc3d,xred,xred_old,ylm,ylmgr)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13paw
 use interfaces_13recipspace
 use interfaces_13xc
 use interfaces_14iowfdenpot
 use interfaces_14wvl_wfs
 use interfaces_15common
 use interfaces_21drive, except_this_one => afterscfloop
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,mgfftc,mgfftf,moved_atm_inside,n3xccc,nfftf,nkxc
 integer,intent(in) :: optres,optxc,prtfor,pwind_alloc,stress_needed,usecprj
 integer,intent(in) :: usexcnhat
 integer,intent(inout) :: computed_forces
 real(dp),intent(in) :: cpus,deltae,gsqcut,res2,residm
 real(dp),intent(inout) :: diffor,etotal,maxfor,vxcavg
 character(len=fnlen),intent(in) :: filapp,filfft
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(efield_type),intent(inout) :: dtefield
 type(energies_type),intent(inout) :: energies
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(results_gs_type),intent(inout) :: results_gs
 type(wffile_type),intent(inout) :: wffnow
 type(wvl_data),intent(inout) :: wvl
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
 integer,intent(in) :: dimcprj(dtset%natom*usecprj)
 integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),nattyp(dtset%ntypat)
 integer,intent(in) :: ngfft(18),ngfftf(18),npwarr(dtset%nkpt)
 integer,intent(in) :: pwind(pwind_alloc,2,3),symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: cg(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
 real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: grewtn(3,dtset%natom)
 real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol),pel_cg(3)
 real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
 real(dp),intent(in) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: rhocore(nfftf),rhore(nfftf,dtset%nspden),rhototp(nfftf)
 real(dp),intent(in) :: rprimd(3,3),tollist(12),vpsp(nfftf)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout) :: forold(3,dtset%natom)
 real(dp),intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw)
 real(dp),intent(inout) :: nvresid(nfftf,dtset%nspden),pel(3)
 real(dp),intent(inout) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
 real(dp),intent(inout) :: ph1df(2,3*(2*dtset%mgfft+1)*dtset%natom),pion(3)
 real(dp),intent(inout) :: rhog(2,nfftf),rhor(nfftf,dtset%nspden),strsxc(6)
 real(dp),intent(inout) :: vhartr(dtset%nfft),vxc(nfftf,dtset%nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc),xred(3,dtset%natom)
 real(dp),intent(out) :: favg(3),fcart(3,dtset%natom),fred(3,dtset%natom)
 real(dp),intent(out) :: gresid(3,dtset%natom),grhf(3,dtset%natom)
 real(dp),intent(out) :: grxc(3,dtset%natom),kxc(nfftf,nkxc),strten(6)
 real(dp),intent(out) :: synlgr(3,dtset%natom),xred_old(3,dtset%natom)
 type(cprj_type),intent(in) :: cprj(dtset%natom,dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol*usecprj)
 type(paw_ij_type),intent(in) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%ntypat*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=6,response=0
 integer :: bantot,choice,ierr,iexit,ii,irhoij,ispden,nele,nhatgrdim,nsize
 integer :: optfor,optgr0,optgr1,optgr2,optrad,optstr,optxc_str,quit,spaceComm
 integer :: vtrial_shift
 real(dp) :: dum,ekin_sum,ekin_sum_p,enxcsr,epot_sum,epot_sum_p,eproj_sum
 real(dp) :: eproj_sum_p,evpart,evsum,lifetime,rcut_coulomb,rcut_coulomb1,ucvol
 logical :: ex,parallel
 character(len=1) :: datacode
 character(len=500) :: message
 type(dataset_type) :: dtLocal
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),pelev(3),rhodum(1),rmet(3,3)
 real(dp),allocatable :: ehart1(:),grnl(:,:),nhatgr(:,:,:),qphon(:),vhart1(:)
 real(dp),allocatable :: vtrial(:),xcart(:,:)

! *************************************************************************

!DEBUG
!write(6,*)' afterscfloop : enter '
!ENDDEBUG

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Perform the positron lifetime calculation
 if (dtset%positron/=0) then
  call calc_lifetime(dtset%ixcpositron,lifetime,nfftf,dtset%nspden,dtset%positron,rhocore,rhore,rhototp,ucvol)
 end if

!Before leaving the present routine, save the current value of xred.
 xred_old(:,:)=xred(:,:)

!Recompute structure factor phases if atomic positions have changed
 if (moved_atm_inside==1) then
  if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
   call getph(atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred)
  else
   ph1d(:,:)=ph1df(:,:)
  end if
 end if

!----------------------------------------------------------------------
!Wavelet case: update HPsi and transform psi to KS orbitals
!----------------------------------------------------------------------
 if (dtset%usewvl == 1) then

! We deallocate the DIIS arrays, if necessary.
  if (dtset%nwfshist > 0) then
   deallocate(wvl%wfs%psidst)
   nullify(wvl%wfs%psidst)
   deallocate(wvl%wfs%hpsidst)
   nullify(wvl%wfs%hpsidst)
   deallocate(wvl%wfs%ads)
   nullify(wvl%wfs%ads)
  end if

  parallel = (mpi_enreg%nproc > 1)
  if (parallel) then
   datacode = 'D'
  else
   datacode = 'G'
  end if

! Apply vtrial and the projectors to the wavefubctions, computing HPsi.
  allocate(vtrial(dtset%nfft * dtset%nspden))
  call wvl_newvtr(dtset, mpi_enreg, nele, vtrial_shift, vhartr, vpsp, vtrial, vxc)
#if defined HAVE_BIGDFT
  call HamiltonianApplication(parallel, datacode, mpi_enreg%me, mpi_enreg%nproc, &
&  dtset%natom, dtset%ntypat, dtset%typat, dtset%wvl_hgrid,&
&  psps%gth_params%psppar, psps%pspcod, wvl%wfs%nstates, wvl%wfs%mbandp, &
&  occ, dtset%wvl_internal%nSize(1), dtset%wvl_internal%nSize(2), &
&  dtset%wvl_internal%nSize(3), dtset%wvl_internal%fineGrid(1, 1), &
&  dtset%wvl_internal%fineGrid(2, 1), dtset%wvl_internal%fineGrid(1, 2), &
&  dtset%wvl_internal%fineGrid(2, 2), dtset%wvl_internal%fineGrid(1, 3), &
&  dtset%wvl_internal%fineGrid(2, 3), wvl%wfs%keys, wvl%wfs%bounds, &
&  wvl%projectors%keys, wvl%projectors%proj, mpi_enreg%ngatherarr, &
&  mpi_enreg%nscatterarr(mpi_enreg%me, 2), vtrial(1 + vtrial_shift), &
&  wvl%wfs%psi, wvl%wfs%hpsi, ekin_sum, epot_sum, eproj_sum, dtset%nsppol, &
&  wvl%wfs%spinar)
  write(message,'(a)') " done"
  call wrtout(6,message,'COLL')

  energies%e_kinetic = ekin_sum
  energies%e_localpsp = epot_sum - two * energies%e_hartree
  energies%e_nonlocalpsp = eproj_sum
  energies%e_corepsp = zero
  etotal = energies%e_kinetic + energies%e_localpsp + energies%e_nonlocalpsp + &
&  energies%e_hartree + energies%e_xc - energies%e_vxc + &
&  energies%e_ewald + energies%e_corepsp

! transform to KS orbitals
  call last_orthon(mpi_enreg%me, mpi_enreg%nproc, parallel, wvl%wfs%nstates_up, &
&  wvl%wfs%nstates_dn, wvl%wfs%nstates, wvl%wfs%mbandp, &
&  wvl%wfs%keys%nvctr_c, wvl%wfs%keys%nvctr_f, wvl%wfs%mvctrp, &
&  dtset%nsppol, wvl%wfs%psi, wvl%wfs%hpsi, wvl%wfs%psit, occ, evsum, &
&  wvl%wfs%eval)
#else
  write(message, '(a,a,a,a)' ) ch10,&
&  ' afterscfloop: BigDFT library is not compiled.', ch10, &
&  '   Action, used the flag --enable-bigdft when configuring.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
#endif

! Change the density according to the KS projection.
  call wvl_mkrho(dtset, mpi_enreg, occ, rhor, wvl%wfs)

! TODO put it at the end of gstate.
! WVL - maybe compute the tail corrections to energy
  if (dtset%tl_radius > real(0, dp)) then
!  Store xcart for each atom
   allocate(xcart(3, dtset%natom))
   call xredxcart(dtset%natom, 1, rprimd, xcart, xred)
!  Use the tails to improve energy precision.
   call wvl_tail_corrections(dtset, energies, etotal, &
&   mpi_enreg, occ, psps, vtrial, wvl, xcart)
   deallocate(xcart)
  end if

 end if

!----------------------------------------------------------------------
!Polarization Calculation
!----------------------------------------------------------------------

 if(dtset%berryopt/=0)then
  call elpolariz(atindx1,cg,cprj,dtefield,dtfil,dtset,etotal,energies%e_elecfield,gprimd,hdr,&
&  kg,dtset%mband,mgfftf,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,nattyp,dtset%nkpt,&
&  npwarr,dtset%nspinor,dtset%nsppol,psps%ntypat,pawang,pawrad,pawtab,pel,pel_cg,pelev,pion,&
&  psps,pwind,pwind_alloc,pwnsfac,rprimd,ucvol,usecprj,wffnow,xred)
 end if

!######################################################################
!Compute forces (if they were not computed during the elec. iterations)
!and stresses (if requested by user)
!----------------------------------------------------------------------

 optfor=0
 if (computed_forces==0.and.dtset%optforces>0.and.dtset%iscf>0) then
  if (dtset%nstep>0.or.dtfil%ireadwf==1) optfor=1
 end if
 if (optfor>0.or.stress_needed>0) then

! MT on 2007-02-20: strsxc should be recomputed for density mixing (to be tested more carefully)
! if (dtset%iscf>=10.and.stress_needed==1) then
! nhatgrdim=0
! if (psps%usepaw==1.and.usexcnhat>0.and.dtset%xclevel==2) then
! nhatgrdim=1;allocate(nhatgr(nfftf,dtset%nspden,3))
! call pawmknhat(dum,1,0,mpi_enreg,dtset%natom,nfftf,ngfftf,nhatgrdim,dtset%nspden,dtset%ntypat,&
! &                  dtset%paral_kgb,pawang,pawfgrtab,nhatgr,nhat,pawrhoij,pawtab,dtset%typat,ucvol)
! end if
! optxc_str=0
! call rhohxc(dtset,enxcsr,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,ngfftf,&
! &   nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,dtset%nspden,n3xccc,optxc_str,rhog,rhor,rprimd,strsxc,&
! &   usexcnhat,vhartr,vxc,vxcavg,xccc3d)
! if (nhatgrdim>0) deallocate(nhatgr)
! end if

! PAW: eventually, compute g_l(r).Y_lm(r) gradients (if not already done)
  if (psps%usepaw==1) then
   if ((pawfgrtab(1)%nfgd==0).or.&
&   (pawfgrtab(1)%gylmgr_allocated==0.and.dtset%pawstgylm==1).or.&
&   (stress_needed==1.and.pawfgrtab(1)%rfgd_allocated==0.and.dtset%pawstgylm==1).or.&
   (pawfgrtab(1)%rfgd_allocated==0.and.dtset%pawstgylm==0)) then
    optgr0=0;optgr1=dtset%pawstgylm;optgr2=0
    optrad=1-dtset%pawstgylm;if (stress_needed==1) optrad=1
    call status(istep,dtfil%filstat,iexit,level,'call nhatgrid ')
    call nhatgrid(atindx1,gmet,mpi_enreg,dtset%natom,nattyp,nfftf,ngfftf,dtset%ntypat,&
&    optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,rprimd,dtset%typat,ucvol,xred)
   end if
  end if

  if (dtset%usewvl == 0) then
   call forstr(atindx1,cg,diffor,dtset,&
&   eigen,energies,favg,fcart,forold,fred,gresid,grewtn,&
&   grhf,grxc,gsqcut,indsym,&
&   kg,kxc,maxfor,mgfftf,mpi_enreg,&
&   n3xccc,nattyp,nfftf,ngfftf,nhat,nkxc,&
&   npwarr,dtset%nspinor,dtset%ntypat,nvresid,occ,optfor,optres,&
&   paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,pel_cg,ph1d,ph1df,pion,&
&   psps,rhog,rhor,rprimd,stress_needed,strsxc,strten,symrec,synlgr,&
&   ucvol,dtfil%unkg,dtfil%unylm,usexcnhat,vhartr,vpsp,vxc,wffnow,&
&   xccc3d,xred,ylm,ylmgr)
  else
   allocate(xcart(3, dtset%natom))
   call xredxcart(dtset%natom, 1, rprimd, xcart, xred)
   allocate(grnl(3, dtset%natom))
   call wvl_nl_gradient(dtset, grnl, mpi_enreg, occ, psps, rprimd, wvl, xcart)
   deallocate(xcart)
   call forces(atindx1, diffor, dtset, favg, fcart, forold, fred, gresid, grewtn,&
&   grhf, grnl, grxc, gsqcut, indsym, kxc, maxfor, mgfftf, mpi_enreg, &
&   0, n3xccc, nattyp, dtset%nfft, dtset%ngfft, nkxc, &
&   dtset%ntypat, pawtab, ph1d, psps, rhog, rhor, rprimd, symrec, &
&   synlgr, nvresid, vxc, xred)
   deallocate(grnl)
  end if
 end if
 if (optfor==1) computed_forces=1
 if (stress_needed==0) strten(:)=9.9999999999D99
 if (computed_forces==0) fcart(:,:)=9.9999999999D99
 if (dtset%prtstm/=0) strten(:)=zero
 if (dtset%positron/=0) strten(:)=zero

!Berry phase: compute stress tensor imposing a constant potential
 if (dtset%berryopt == 4 .and. stress_needed/=0) then
  write(message,'(a,a,a)')ch10,&
&  ' Stress tensor imposing a constant potential drop across',&
&  '  each lattice vector:'
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

! Compute stress tensor
  call dtsetCopy(dtLocal, dtset)
  dtLocal%berryopt = 0
  optfor=0
  call forstr(atindx1,cg,diffor,dtLocal,&
&  eigen,energies,favg,fcart,forold,fred,gresid,grewtn,&
&  grhf,grxc,gsqcut,indsym,&
&  kg,kxc,maxfor,mgfftf,mpi_enreg,&
&  n3xccc,nattyp,nfftf,&
&  ngfftf,nhat,nkxc,npwarr,&
&  dtset%nspinor,dtset%ntypat,nvresid,occ,optfor,optres,paw_ij,pawang,&
&  pawfgr,pawfgrtab,pawrhoij,pawtab,pel_cg,ph1d,ph1df,pion,psps,&
&  rhog,rhor,rprimd,stress_needed,strsxc,&
&  strten,symrec,synlgr,ucvol,dtfil%unkg,&
&  dtfil%unylm,usexcnhat,vhartr,vpsp,vxc,wffnow,xccc3d,xred,ylm,ylmgr)
  call dtsetFree(dtLocal)
 end if

!If SCF convergence was not reached (for dtset%nstep>0),
!print a warning to the output file (non-dummy arguments: dtset%nstep,
!residm, diffor - infos from tollist have been saved inside )
 choice=3
 call scprqt(choice,cpus,deltae,diffor,dtset,&
& eigen,etotal,favg,fcart,energies%e_fermie,filapp,dtfil%filnam_ds(1),&
& 1,dtset%iscf,istep,dtset%kptns,maxfor,&
& moved_atm_inside,mpi_enreg,dtset%nband,dtset%nkpt,&
& dtset%nstep,occ,optres,prtfor,quit,&
& res2,resid,residm,response,tollist,psps%usepaw,vxcavg,dtset%wtk,xred)

 if(dtset%nfreqsus>0)then
! At the end of the SCF cycle compute the short-range correlation energy
! in RPA-LSD (dtset%ixc==10). Note: after call to rhohxc, vxc indeed correponds to dtset%ixc==10.
  call dtsetCopy(dtLocal, dtset)
  dtLocal%ixc = 10
  nhatgrdim=0
  call rhohxc(dtLocal,enxcsr,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,ngfftf,&
&  nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,dtset%nspden,n3xccc,optxc,rhog,rhor,rprimd,strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d)
  write(message,'(4a,1x,es16.10,2a)') ch10,&
&  ' scfcv: beyond-RPA short-range correlation energy in LSDA ---', ch10,&
&  '  LSD-RPA+ correlation energy=', enxcsr, ' hartree',ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(6,message,'COLL')
! Restore the xc potential.
  if (psps%usepaw==1.and.usexcnhat>0 .and. dtset%xclevel==2) then
   nhatgrdim=1;allocate(nhatgr(nfftf,dtset%nspden,3))
   call pawmknhat(dum,1,0,mpi_enreg,dtset%natom,nfftf,ngfftf,nhatgrdim,dtset%nspden,dtset%ntypat,&
&   dtset%paral_kgb,pawang,pawfgrtab,nhatgr,nhat,pawrhoij,pawtab,dtset%typat,ucvol)
  end if
  call rhohxc(dtset,enxcsr,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,ngfftf,&
&  nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,dtset%nspden,n3xccc,optxc,rhog,rhor,rprimd,strsxc,&
&  usexcnhat,vhartr,vxc,vxcavg,xccc3d)
  if (nhatgrdim>0) deallocate(nhatgr)
! Also check the Hartree energy with a real-space cutoff for the Coulomb interaction
  allocate(ehart1(3),qphon(3),vhart1(nfftf))
  vhart1(:)=zero;qphon(:)=zero;rcut_coulomb=1.d20
  do ii=1,3
   rcut_coulomb=min(rcut_coulomb,sum(rprimd(:,ii)*rprimd(:,ii)))
  end do
  rcut_coulomb=0.5_dp*sqrt(rcut_coulomb)
  write(6,*) '%scfcv: Check cutoff Coulomb interaction, show Hartree energy:'
  write(6,'(1x,a,1x,a,1x,a,6x,a)') 'cutoff radius',&
&  ' with cutoff (with G=0)','(without G=0)','no cutoff'
  rcut_coulomb1=max(0.5*rcut_coulomb,rcut_coulomb-2._dp)
  do while (rcut_coulomb1<min(2.0*rcut_coulomb,rcut_coulomb+2._dp))
   call hartre1(1,gmet,gsqcut,nfftf,ngfftf,dtset%paral_kgb,qphon,rhog,vhart1,ehart1,&
&   rcut_coulomb1,ucvol)
   write(6,'(f14.8,1x,3(3x,es14.7))') rcut_coulomb1,ehart1(3),ehart1(1),ehart1(2)
   rcut_coulomb1=rcut_coulomb1+0.5_dp
  end do
  call hartre1(1,gmet,gsqcut,nfftf,ngfftf,dtset%paral_kgb,qphon,rhog,vhart1,ehart1,&
&  rcut_coulomb,ucvol)
  write(message,'(4a,1x,es16.8,3a,1x,es16.8,2a)') ch10,&
&  ' scfcv: Hartree energy with real-space cutoff for Coulomb interaction ---', ch10,&
&  '   Cutoff radius=', rcut_coulomb, ' bohr',ch10,&
&  '  Hartree energy=', ehart1(3)   , ' hartree',ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(6,message,'COLL')
  call dtsetFree(dtLocal)
  if(allocated(ehart1)) deallocate(ehart1)
  if(allocated(qphon))  deallocate(qphon)
  if(allocated(vhart1))  deallocate(vhart1)
 end if

!get current operator on wavefunctions
 if (dtset%prtspcur == 1) then
  call spin_current(atindx,atindx1,cg,dtfil,dtset,eigen,gmet,gprimd,hdr,kg,mpi_enreg,&
&  nattyp,nfftf,ph1d,psps,rhog,rhor,rmet,symrec,ucvol,wffnow,ylm,ylmgr)

 end if

!Delete eventual _FFT file
 if(dtset%mffmem==0)then
  inquire (file=filfft,exist=ex)
  if(ex)then
   open(unit=tmp_unit,file=filfft,form='unformatted',status='old')
   close(unit=tmp_unit,status='DELETE')
  end if
 end if

!Update the content of the header (evolving variables)
 bantot=hdr%bantot
 call hdr_update(bantot,etotal,energies%e_fermie,hdr,dtset%natom,residm,rprimd,occ,pawrhoij,psps%usepaw,xred)

!XG 070612 : Do not remove this line - needed for the pathscale compiler ?!
!write(6,*)' afterscfloop : fred=',fred(:,:)

 results_gs%energies   = energies
 results_gs%etotal     =etotal
 results_gs%fcart(:,:) =fcart(:,:)
 results_gs%fred(:,:)  =fred(:,:)
 results_gs%gresid(:,:)=gresid(:,:)
 results_gs%grewtn(:,:)=grewtn(:,:)
 results_gs%grxc(:,:)  =grxc(:,:)
 results_gs%pel(1:3)   =pel(1:3)
 results_gs%residm     =residm
 results_gs%strten(1:6)=strten(1:6)
 results_gs%synlgr(:,:)=synlgr(:,:)
 results_gs%vxcavg     =vxcavg
 if (dtset%nstep == 0 .and. dtset%occopt>=3.and.dtset%occopt<=7) then
  results_gs%etotal = results_gs%etotal - dtset%tsmear * results_gs%entropy
 end if

!DEBUG
!write(6,*)' afterscfloop : exit'
!stop
!ENDDEBUG

end subroutine afterscfloop
!!***
