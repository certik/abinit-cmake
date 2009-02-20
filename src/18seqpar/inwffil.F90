!{\src2tex{textfont=tt}}
!!****f* ABINIT/inwffil
!! NAME
!! inwffil
!!
!! FUNCTION
!! Do initialization of wavefunction files.
!! Also call other relevant routines for this initialisation
!!    (initialization of wavefunctions from scratch or from file,
!!     translations of wavefunctions, ...)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, AR, MB, MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ask_accurate= if 1, the wavefunctions and eigenvalues must be
!!    accurate, that is, they must come from a k point that is
!!    symmetric of the needed k point, with a very small tolerance,
!!    the disk file contained sufficient bands to initialize all of them,
!!    the spinor and spin-polarisation characteristics must be identical
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ecut=effective kinetic energy planewave cutoff (hartree), beyond
!!    which the coefficients of plane waves are zero
!!  ecut_eff=effective kinetic energy planewave cutoff (hartree), needed
!!    to generate the sphere of plane wave
!!  exchn2n3d=if 1, n2 and n3 are exchanged
!!  formeig=explained above
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  ireadwf=option parameter described above for wf initialization
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!    to be initialized here.
!!  kg(3,mpw*mkmem)=dimensionless coords of G vecs in basis sphere at k point
!!  kptns(3,nkpt)=reduced coords of k points
!!  localrdwf=(for parallel case) if 1, the wffnm  file is local to each machine
!!  mband=maximum number of bands
!!  mkmem=maximum number of k-points in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum number of planewaves as dimensioned in calling routine
!!  nband(nkpt*nsppol)=number of bands at each k point
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  npwarr(nkpt)=array holding npw for each k point.
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetry elements in space group
!!  occ(mband*nkpt*nsppol)=occupations (from disk or left at their initial value)
!!  optorth= 1 if the WFS have to be orthogonalized; 0 otherwise
!!  prtvol=control print volume and debugging
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,nsym)=symmetry operations in real space in terms
!!   of primitive translations
!!  tnons(3,nsym)=nonsymmorphic translations for symmetry operations
!!  unkg=unit number for storage of basis sphere data: stores indirect
!!   indexing array and integer coordinates for all planewaves in basis
!!   sphere for each k point being considered
!!  unwff1,unwfnow=
!!            unit numbers for files wffnm and wft1nm.
!!  wffnm=name (character data) of file for input wavefunctions.
!!  wft1nm= name (character data) of file used for temporary wf storage.
!!
!! OUTPUT
!!  wff1  = structure information for files wffnm .
!!  wffnow= structure information for wf file wft1nm
!!  if ground state format (formeig=0):
!!    eigen(mband*nkpt*nsppol)=eigenvalues (input or init to large number), (Ha)
!!  if respfn format (formeig=1):
!!    eigen(2*mband*mband*nkpt*nsppol)=matrix of eigenvalues
!!                                     (input or init to large number), (Ha)
!! Conditional output (returned if mkmem/=0):
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=complex wf array
!!    be careful : an array of size cg(2,npw*nspinor), as used
!!    in the response function code, is not enough !
!!  wvl <type(wvl_data)>=all wavelets data.
!!
!! NOTES
!! Detailed description :
!!  Initialize unit wff1%unwff for input of wf data if ireadwf=1
!!  Opens file on unit wffnow%unwff
!!   if the storage on disk is needed (mkmem==0)
!!  Initializes wf data on wffnow%unwff, by calling the appropriate routine.
!!
!! formeig option (format of the eigenvalues and occupations) :
!!   0 => ground-state format (initialisation of
!!        eigenvectors with random numbers, vector of eigenvalues,
!!        occupations are present)
!!   1 => respfn format (initialisation of
!!        eigenvectors with 0 s, hermitian matrix of eigenvalues)
!!
!! ireadwf options:
!!   0 => initialize with random numbers or 0 s
!!   1 => read from disk file wff1, initializing higher bands
!!        with random numbers or 0 s if not provided in disk file
!!
!! The wavefunctions after this initialisation are stored in unit wffnow%unwff
!!
!! WARNINGS
!! The symmetry operations are used to translate the data from one
!! k point to another, symmetric, k point.
!! They can be completely different from the symmetry operations
!! contained on the disk file. No check is performed between the two sets.
!!
!! Occupations will not be modified nor output,
!!  in the present status of this routine.
!!
!! If ground state format (formeig=0) occ(mband*nkpt*nsppol) was output.
!! NOT OUTPUT NOW !
!!
!! PARENTS
!!      gstate,loop3dte,loper3,nonlinear,respfn
!!
!! CHILDREN
!!      handle_ncerr,hdr_check,hdr_clean,hdr_io,hdr_io_etsf,hdr_io_netcdf
!!      ini_wf_netcdf,leave_new,listkk,matr3inv,newkpt,rhoij_copy,timab,wffkg,wffopen
!!      wfsinp,wrtout,wvl_wfsinp,xcomm_init,xdefineoff,xmaster_init,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine inwffil(ask_accurate,cg,dtset,ecut,ecut_eff,eigen,exchn2n3d,&
&           formeig,gmet,hdr,ireadwf,istwfk,kg,kptns,localrdwf,mband,&
&           mkmem,mpi_enreg,mpw,nband,ngfft,nkpt,npwarr,nspden,nspinor,&
&           nsppol,nsym,occ,optorth,psps,prtvol,rprimd,symafm,symrel,tnons,unkg,wff1,&
&           wffnow,unwff1,unwfnow,wffnm,wft1nm, wvl)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes

#if defined HAVE_NETCDF
 use netcdf
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13io_mpi
 use interfaces_13ionetcdf
 use interfaces_14iowfdenpot
 use interfaces_14wfs
 use interfaces_15common
 use interfaces_18seqpar, except_this_one => inwffil
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: ask_accurate,exchn2n3d,formeig,ireadwf,localrdwf,mband,mkmem,mpw
 integer, intent(in) :: nkpt,nspden,nsppol,nsym,optorth,prtvol,unkg,unwff1,unwfnow
 integer, intent(inout) :: nspinor
 real(dp), intent(in) :: ecut,ecut_eff
 character(len=fnlen), intent(in) :: wffnm,wft1nm
 type(MPI_type), intent(inout) :: mpi_enreg
 type(dataset_type), intent(in) :: dtset
 type(hdr_type), intent(inout) :: hdr
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type), intent(out) :: wff1,wffnow
 type(wvl_data), intent(inout) :: wvl
 integer, intent(in) :: istwfk(nkpt),kg(3,mpw*mkmem),nband(nkpt*nsppol),ngfft(18)
 integer, intent(in) :: npwarr(nkpt),symafm(nsym),symrel(3,3,nsym)
 real(dp), intent(out) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp), intent(out) :: eigen((2*mband)**formeig*mband*nkpt*nsppol)
 real(dp), intent(in) :: gmet(3,3),kptns(3,nkpt),rprimd(3,3),tnons(3,nsym)
 real(dp), intent(inout) :: occ(mband*nkpt*nsppol)

!Local variables-------------------------------
 integer,save :: count=0
 integer :: accesswff,accurate,ceksp,date0,debug,doorth,fform,fform_dum,fill
 integer :: headform0,iband,icg1,ierr,ii,ikassoc,ikpt,ikpt0,ikptsp
 integer :: ikptsp0,increase_nkassoc,ios,irec,ispden,isppol,isppol0
 integer :: master,mband0,mband0_rd,mcg,mcg_disk,me,mkmem0,mpw0,mu,nband0_k
 integer :: nkassoc,nkpt0,npw,npw0,nspinor0,nsppol0,rdwr,restart,restartpaw
 integer :: spaceComm,sppoldbl,squeeze,tim_rwwf
 real(dp) :: dksqmax,ecut0
 character(len=fnlen) :: filname
 character(len=500) :: message
 type(hdr_type) :: hdr0
 integer :: mdims(3),ngfft0(18)
 integer,allocatable :: indkk(:,:),indkk0(:,:),istwfk0(:),nband0(:)
 integer,allocatable :: nband0_rd(:),npwarr0(:)
 integer :: ncid_hdr,ncerr,ndim,nvar,natt,uid
 real(dp) :: gmet0(3,3),gprim0(3,3),rprim0(3,3),tsec(2)
 real(dp),allocatable :: cg_disk(:,:),kptns0(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' inwffil: enter'
!count=count+1
!write(std_out,*)' count=',count
!if(count==13)stop
!stop
!ENDDEBUG

!Keep track of total time spent in inwffil
 call timab(19,1,tsec)

!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)

!PATCH inwffil // KPT & FFT spaceComm --> comm_kpt
  !print *, "inwffil", mpi_enreg%paral_compil_kpt, mpi_enreg%paral_compil_fft
 if ((mpi_enreg%paral_compil_kpt==1) .and. &
    &(mpi_enreg%paral_compil_fft==1)) then
     spaceComm = mpi_enreg%commcart_3d
 endif

!Init me
 call xme_init(mpi_enreg,me)
!Init master
 call xmaster_init(mpi_enreg,master)

!Check the validity of formeig
 if(formeig/=0.and.formeig/=1)then
  write(message, '(a,a,a,a,i12,a)' ) ch10,&
&   ' inwffil: BUG -',ch10,&
&   '  formeig=',formeig,' , but the only allowed values are 0 or 1.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 ngfft0(:)=ngfft(:)

!Default value for headform0 (will be needed later, to read wf blocks)
 headform0=0

!If the input data are on disk, determine the kind of restart
 if (ireadwf==1)then

  accesswff=dtset%accesswff
  if(localrdwf==0)accesswff=-1    ! This is in case the wff file must be read by only the master proc

!DEBUG
! write(std_out,*)' inwffil : will open file ',trim(wffnm)
! write(std_out,*)' inwffil : on unit', unwff1
! call flush(6)
!ENDDEBUG
  call WffOpen(accesswff,spaceComm,wffnm,ierr,wff1,master,me,unwff1)

!DEBUG
! write(std_out,*)' inwffil : before hdr_io '
! call flush(6)
!ENDDEBUG

#if defined HAVE_NETCDF
!DEBUG
!if (accesswff == 2) then
!  ncid_hdr = wff1%unwff
!  ncerr = nf90_Inquire(ncid=ncid_hdr,nDimensions=ndim,nVariables=nvar,nAttributes=natt,unlimitedDimId=uid)
!  call handle_ncerr(ncerr, " general Inquire ")
!  write (*,*) 'inwffil : found ndim,nvar,natt,uid = ', ndim,nvar,natt,uid
!end if
!ENDDEBUG
#endif

! Initialize hdr0 (sent to all procs), thanks to reading of wff1
  rdwr=1
  if (wff1%accesswff < 2) then
   call hdr_io(fform_dum,hdr0,rdwr,wff1)
#if defined HAVE_NETCDF
  else if (wff1%accesswff == 2) then
   call hdr_io_netcdf(fform_dum,hdr0,rdwr,wff1)
#endif
#if defined HAVE_ETSF_IO
  else if (wff1%accesswff == 3) then
   call hdr_io_etsf(fform_dum, hdr0, rdwr, wff1%unwff)
#endif
  end if

!DEBUG
! write(std_out,*)' inwffil : before hdr_check '
! call flush(6)
!ENDDEBUG

! Check hdr0 versus hdr
! (and from now on ignore header consistency and write new info
! to header for each file)
  if (dtset%usewvl == 0) then
    ! wait for plane waves.
    fform=2
  else
    ! wait for wavelets.
    fform = 200
  end if
! call hdr_check(fform,fform_dum,hdr,hdr0,'COLL',restart)
  call hdr_check(fform,fform_dum,hdr,hdr0,'PERS',restart,restartpaw)

  nkpt0=hdr0%nkpt
  nsppol0=hdr0%nsppol
  headform0=hdr0%headform


  write(message, '(a,a)' )&
&  '-inwffil : will read wavefunctions from disk file ',trim(wffnm)
  call wrtout(ab_out,message,'COLL')

 else

  restart=1 ; restartpaw=0

! Fill some data concerning an hypothetical file to be read
! This is to allow the safe use of same routines than with ireadwf==1.
  nkpt0=nkpt ; nsppol0=nsppol

 end if ! end ireadwf



 sppoldbl=1
 if(minval(symafm(:))==-1)then
  if(nsppol0==1 .and. nsppol==2)sppoldbl=2
 end if

 allocate(indkk(nkpt*sppoldbl,6),istwfk0(nkpt0),kptns0(3,nkpt0))
 allocate(nband0(nkpt0*nsppol0),npwarr0(nkpt0))

 if(restart==2)then ! restart with translations

  ecut0=hdr0%ecut_eff
  istwfk0(1:nkpt0)=hdr0%istwfk(1:nkpt0)
  kptns0(1:3,1:nkpt0)=hdr0%kptns(1:3,1:nkpt0)
  nband0(1:nkpt0*nsppol0)=hdr0%nband(1:nkpt0*nsppol0)
  ngfft0(1:3)=hdr0%ngfft(1:3)
  npwarr0(1:nkpt0)=hdr0%npwarr(1:nkpt0)
  nspinor0=hdr0%nspinor
  rprim0(:,:)=hdr0%rprimd(:,:)
  mpw0=maxval(npwarr0(:))

! Compute reciprocal space metric gmet for unit cell of disk wf
  call matr3inv(rprim0,gprim0)
  do ii=1,3
   gmet0(:,ii)=gprim0(1,:)*gprim0(1,ii)+&
&              gprim0(2,:)*gprim0(2,ii)+&
&              gprim0(3,:)*gprim0(3,ii)
  end do

! At this stage, the header of the file wff1i%unwff is read, and
! the pointer is ready to read the first wavefunction block.

! Compute k points from input file closest to the output file
  call listkk(dksqmax,gmet0,indkk,kptns0,kptns,nkpt0,nkpt,nsym,sppoldbl,&
&  symafm,symrel,1)

 else if (restart==1) then ! direct restart

! Fill variables that must be the same, as determined by hdr_check.f
! This is to allow the safe use of the same routines than with restart==2.
  nspinor0=nspinor
  ecut0=ecut_eff
  gmet0(:,:)=gmet(:,:)
  istwfk0(:)=istwfk(:)
  kptns0(:,:)=kptns(:,:)
  npwarr0(:)=npwarr(:)
  mpw0=mpw

  do ikpt=1,nkpt
   indkk(ikpt,1)=ikpt
   indkk(ikpt,2:6)=0
  end do
  dksqmax=0.0_dp

! The treatment of nband0 asks for some care
  if(ireadwf==0)then
   nband0(:)=0
  else
   nband0(1:nkpt0*nsppol0)=hdr0%nband(1:nkpt0*nsppol0)
  end if

 end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Before hdr_clean:
! If restartpaw==1, store hdr0%pawrhoij in hdr%pawrhoij; else if restartpaw==0,
! hdr%pawrhoij(:)has been initialized in hdr_init.
 if(restartpaw==1) then
  call rhoij_copy(hdr0%pawrhoij,hdr%pawrhoij)
 end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!At this stage, all the relevant information from the header of the disk file,
!has been exploited, and stored in variables, on all processors.
!It is also contained in hdr0
!(on all processors, except if restart=1 and localrdwf=0,
!in which case it is only on the master)
!These information might be changed later, while processing the
!wavefunction data, and converting it. The variable hdr0 might be kept
!for further checking, or reference, or debugging, but at present,
!it is simpler to close it. The other header, hdr, will be used
!for the new file, if any.

 if(ask_accurate==1)then

! Check whether the accuracy requirements might be fulfilled
  if(ireadwf==0)then
   write(message, '(a,a,a,a,a, a,a,a,a,a, a,a)' ) ch10, &
&   ' inwffil: ERROR ',ch10,&
&   '  The file ',trim(wffnm),' cannot be used to start the ',ch10,&
&   '  present calculation. It was asked that the wavefunctions be accurate,',ch10,&
&   '  but they were not even read.',ch10,&
&   '  Action: use a wf file, with ireadwf/=0.'
   call wrtout(06,  message,'PERS')
   call leave_new('PERS')
  end if
  if(dksqmax>tol12)then
   write(message, '(12a,es16.6,4a)' ) ch10, &
&   ' inwffil: ERROR ',ch10,&
&   '  The file ',trim(wffnm),' cannot be used to start the ',ch10,&
&   '  present calculation. It was asked that the wavefunctions be accurate, but',ch10,&
&   '  at least one of the k points could not be generated from a symmetrical one.',ch10,&
&   '  dksqmax=',dksqmax,ch10,&
&   '  Action: check your wf file and k point input variables',ch10,&
&   '    (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
   call wrtout(06,  message,'PERS')
   call leave_new('PERS')
  end if
  if(nspinor/=nspinor0)then
   write(message, '(a,a,a,a,a, a,a,a,a,a, a,a,2i5,a,a)' ) ch10, &
&   ' inwffil: ERROR ',ch10,&
&   '  The file ',trim(wffnm),' cannot be used to start the ',ch10,&
&   '  present calculation. It was asked that the wavefunctions be accurate, but',ch10,&
&   '  nspinor differs in the file from the actual nspinor.',ch10,&
&   '  nspinor,nspinor0=',nspinor,nspinor0,ch10,&
&   '  Action: check your wf file, and nspinor input variables.'
   call wrtout(06,  message,'PERS')
   call leave_new('PERS')
  end if
  if((nsppol>nsppol0 .and. sppoldbl==1) .or. nsppol<nsppol0 ) then
   write(message, '(a,a,a,a,a, a,a,a,a,a, a,a,3i5,a,a)' ) ch10, &
&   ' inwffil: ERROR ',ch10,&
&   '  The file ',trim(wffnm),' cannot be used to start the ',ch10,&
&   '  present calculation. It was asked that the wavefunctions be accurate, but',ch10,&
&   '  the nsppol variables do not match in the file and in the actual calculation',ch10,&
&   '  nsppol,nsppol,sppoldbl=',nspinor,nspinor0,sppoldbl,ch10,&
&   '  Action: check your wf file, and nsppol input variables.'
   call wrtout(06,  message,'PERS')
   call leave_new('PERS')
  end if

! Now, check the number of bands
  accurate=1
  do isppol=1,nsppol
   do ikpt=1,nkpt
    ikpt0=indkk(ikpt+(isppol-1)*(sppoldbl-1)*nkpt,1)
    ikptsp =ikpt +(isppol-1)*nkpt
    ikptsp0=ikpt0+(isppol-1)*(2-sppoldbl)*nkpt0
    if(nband0(ikptsp0)<nband(ikptsp))accurate=0
   end do
  end do
  if(accurate==0)then
   write(message, '(a,a,a,a,a, a,a,a,a,a, a,a)' ) ch10, &
&   ' inwffil: ERROR ',ch10,&
&   '  The file ',trim(wffnm),' cannot be used to start the ',ch10,&
&   '  present calculation. It was asked that the wavefunctions be accurate,',ch10,&
&   '  but the number of bands differ in the file and in the actual calculation.',ch10,&
&   '  Action: use a wf file with the correct characteristics.'
   call wrtout(06,  message,'PERS')
   call leave_new('PERS')
  end if

 end if

!Not all bands might be read, if not needed to fill the wavefunctions
 mband0=maxval(nband0(1:nkpt0*nsppol0))
 mband0_rd=min(mband0,(mband/nspinor)*nspinor0)

!Open the wf file on unit wffnow%unwff if mkmem==0
!Also allocate temporary array, and write header at top of wffnow%unwff

 if (mkmem==0) then

#if defined HAVE_NETCDF
    if(dtset%accesswff==2) then
        accesswff=2
!  Create empty netCDF file
        ncerr = nf90_create(path=trim(wft1nm), cmode=NF90_CLOBBER, ncid=ncid_hdr)
        call handle_ncerr(ncerr," inwffil create netcdf wavefunction file")
        ncerr = nf90_close(ncid_hdr)
        call handle_ncerr(ncerr," inwffil close netcdf wavefunction file")
    else if (dtset%accesswff == 3) then
     write (std_out,*) "FIXME: ETSF I/O support in inwffil"
    end if
#endif

  write(message, '(a,i4,a,a)' ) &
&  ' inwffil about to open unit',unwfnow,' for file=',trim(wft1nm)
  call wrtout(06,  message,'PERS')
  call WffOpen(dtset%accesswff,spaceComm,wft1nm,ierr,wffnow,master,me,unwfnow)

  rdwr=2 ; fform=2
  if (wffnow%accesswff /= 2) then
   call hdr_io(fform,hdr,rdwr,wffnow)
#if defined HAVE_NETCDF
  else if (wffnow%accesswff == 2) then
   call hdr_io_netcdf(fform,hdr,rdwr,wffnow)

   call ini_wf_netcdf(mpw,wffnow%unwff,formeig)
  else if (wffnow%accesswff == 3) then
   write (std_out,*) "FIXME: ETSF I/O support in inwffil"
#endif
  end if


! Note that the cg_disk array is allocated when cg uses no space (mkmem==0)
! This dimensioning allows to treat all cases
  mcg_disk=max(mpw0*nspinor0*mband0_rd,mpw*nspinor*mband)
  allocate(cg_disk(2,mcg_disk))

! Define offsets, in case of MPI I/O
  call WffKg(wffnow,1)   ! option optkg in wfsinp
  call xdefineOff(formeig,wffnow,mpi_enreg,hdr%nband,hdr%npwarr,hdr%nspinor,hdr%nsppol,hdr%nkpt)
 end if


!****************************************************************************
!If needed, transfer the input wf from disk to core memory
!(in the parallel case, it allows to change localrdwf=0 in localrdwf=1)

 mkmem0=0

 if(mpi_enreg%paral_compil_kpt == 1 .or. mpi_enreg%paral_compil_fft == 1) then
  if(localrdwf==0 .and. mkmem==0)then
   write(message,'(a,a,a,a)')ch10,&
&   ' inwffil: BUG -',ch10,&
&   '  localrdwf==0 and mkmem==0 are not allowed together (yet)'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
 end if

!Here, treat reading wavefunctions with mkmem/=0, first step
 if(ireadwf==1 .and. mkmem/=0)then
!if(restart==1 .and. ireadwf==1 .and. mkmem/=0)then

! Compute table of k point associations. Make a trial choice for nkassoc.
  nkassoc=(nkpt/nkpt0+1)*2
  allocate(indkk0(nkpt0,nkassoc))
! Infinite loops are allowed in F90
  do
   indkk0(:,:)=0
   increase_nkassoc=0
   do ikpt=1,nkpt*sppoldbl
    ikpt0=indkk(ikpt,1)
    do ikassoc=1,nkassoc
     if(indkk0(ikpt0,ikassoc)==0)then
      indkk0(ikpt0,ikassoc)=ikpt
      exit
     end if
     if(nkassoc==ikassoc)increase_nkassoc=1
    end do
    if(increase_nkassoc==1)then
     deallocate(indkk0)
     nkassoc=2*nkassoc
     allocate(indkk0(nkpt0,nkassoc))
     exit
    end if
   end do
   if(increase_nkassoc==0)exit
  end do

!DEBUG
! write(std_out,*)' inwffil: indkk0, nkassoc=',nkassoc
! do ikpt0=1,nkpt0
!  write(std_out,*)' ikpt0,indkk0(ikpt0,1)=',ikpt0,indkk0(ikpt0,1)
! end do
!ENDDEBUG

!DEBUG
! write(std_out,*)' inwffil : indkk(:,1)=',indkk(:,1)
! write(std_out,*)' inwffil : sppoldbl=',sppoldbl
!ENDDEBUG

  squeeze=0
  if(mkmem/=0)then
   allocate(nband0_rd(nkpt0*nsppol0))
   nband0_rd(:)=0
   do isppol=1,nsppol
    do ikpt=1,nkpt
     ikpt0=indkk(ikpt+(isppol-1)*(sppoldbl-1)*nkpt,1)
     isppol0=min(isppol,nsppol0)
     ikptsp =ikpt +(isppol -1)*nkpt
     ikptsp0=ikpt0+(isppol0-1)*(2-sppoldbl)*nkpt0
     nband0_k=min( nband0(ikptsp0)                  , &
&                 (nband (ikptsp )/nspinor)*nspinor0 )
     nband0_rd(ikptsp0)=max(nband0_rd(ikptsp0),nband0_k)

     npw0=npwarr0(ikpt0)
     npw =npwarr (ikpt)
     if(npw0*nspinor0*nband0_k > npw*nspinor*nband(ikptsp))squeeze=1
    end do
   end do
   if(squeeze==1)then
    mcg_disk=mpw0*nspinor0*mband0_rd
    allocate(cg_disk(2,mcg_disk))
   else
    if(mpi_enreg%paral_compil_kpt == 1 .or. mpi_enreg%paral_compil_fft == 1)then
     if(localrdwf==0)then
      mcg_disk=mpw0*nspinor0*mband0_rd
      allocate(cg_disk(2,mcg_disk))
     end if
    end if
   end if
  end if




!DEBUG
! write(std_out,*)' inwffil : nband0_rd=',nband0_rd
! write(message,'(a,i4)')' inwffil : will call wfsinp, squeeze=',squeeze
! call wrtout(6,message,'COLL')
!ENDDEBUG
  if (dtset%usewvl == 0) then
    call wfsinp(cg,cg_disk,ecut,ecut0,ecut_eff,eigen,&
  &  exchn2n3d,formeig,gmet,gmet0,headform0,&
  &  indkk,indkk0,istwfk,istwfk0,kptns,kptns0,localrdwf,&
  &  mband,mband0_rd,mcg_disk,mkmem,mpi_enreg,mpw,mpw0,&
  &  nband,nband0_rd,ngfft,nkassoc,nkpt,nkpt0,npwarr,npwarr0,nspinor,nspinor0,&
  &  nsppol,nsppol0,nsym,occ,optorth,prtvol,restart,rprimd,sppoldbl,squeeze,&
  &  symafm,symrel,tnons,wff1,wffnow)
  else
     ! Read wavefunctions from file.
     call wvl_wfsinp_disk(dtset, hdr0, hdr, mpi_enreg, 1, &
          & hdr%rprimd, wff1, wvl%wfs, hdr%xred)
  end if

! Now, update xyz0 variables, for use in newkpt
  nband0(:)=nband0_rd(:)

#if defined FC_FUJITSU
!! DEBUG by MM for VPP
! write(std_out,*) 'size of nband0 =',lbound(nband0),ubound(nband0)
! write(std_out,*) 'size of nband0_rd =',lbound(nband0_rd),ubound(nband0_rd)
! write(std_out,*) 'nband0, nband0_rd =', nband0(:), ',', nband0_rd(:) !necessary
!! END DEBUG by MM
#endif

! If squeeze, the conversion was done in wfsinp, so no conversion left.
  if(squeeze==1)then
   ecut0=ecut_eff
   gmet0(:,:)=gmet(:,:)
   deallocate(kptns0,istwfk0,nband0,npwarr0)
   allocate(kptns0(3,nkpt),istwfk0(nkpt),nband0(nkpt*nsppol),npwarr0(nkpt))
   kptns0(:,:)=kptns(:,:)
   istwfk0(:)=istwfk(:)
   npwarr0(:)=npwarr(:)
   nband0(:)=0
   do isppol=1,nsppol
    do ikpt=1,nkpt
     ikpt0=indkk(ikpt+(isppol-1)*(sppoldbl-1)*nkpt,1)
     isppol0=min(isppol,nsppol0)
     ikptsp =ikpt +(isppol -1)*nkpt
     ikptsp0=ikpt0+(isppol0-1)*(sppoldbl-1)*nkpt0
     nband0(ikptsp)=(nband0_rd(ikptsp0)/nspinor0)*nspinor

#if defined FC_FUJITSU
!! DEBUG by MM for VPP
!     write(std_out,*) 'nband0(',ikptsp,'), nband0_rd(',ikptsp0,')=',&
!&     nband0(ikptsp),nband0_rd(ikptsp0)
!     write(std_out,*) 'ikpt=',ikpt,' ikpt0=',ikpt0
!     write(std_out,*) 'isppol=',isppol,' isppol0=',isppol0
      write(std_out,*) 'ikptsp=',ikptsp,' ikptsp0=',ikptsp0         !this four lines
      write(std_out,*) 'nspinor=',nspinor,' nspinor0=',nspinor0     ! necessary
      write(std_out,*) 'nband0   (',ikptsp,')=',nband0(ikptsp)      !  for VPP
      write(std_out,*) 'nband0_rd(',ikptsp0,')=',nband0_rd(ikptsp0) !don''t comment!
!! END DEBUG by MM
#endif

    end do
   end do
   do ikpt=1,nkpt
    indkk(ikpt,1)=ikpt
    indkk(ikpt,2:6)=0
   end do
!  This transfer must come after the nband0 transfer
   nspinor0=nspinor
   nkpt0=nkpt
   nsppol0=nsppol
  end if ! end squeeze == 1

! The input wavefunctions have been transferred from disk to core memory
  mkmem0=mkmem

  deallocate(indkk0,nband0_rd)
  if(squeeze==1)deallocate(cg_disk)
  if(mpi_enreg%paral_compil_fft == 1 .or. mpi_enreg%paral_compil_kpt == 1 )then
   if(localrdwf==0)deallocate(cg_disk)
  end if
 else !ireadwf == 0
  if (dtset%usewvl == 1) then
     ! Compute wavefunctions from input guess.
     call wvl_wfsinp_scratch(dtset, mpi_enreg, psps, hdr%rprimd, wvl, hdr%xred)
  end if
 end if

!Clean hdr0
 if (ireadwf==1)then
  if( restart==2 .or. localrdwf==1 .or. master==me)then
   call hdr_clean(hdr0)
  end if
 end if

!****************************************************************************
!Now, treat translation of wavefunctions if wavefunctions are planewaves

 ceksp=0 ; debug=0 ; doorth=1 ; fill=1

!DEBUG
!write(6,*)' inwffil : will call newkpt, me= ',me,wff1%offwff,wffnow%offwff
!write(std_out,*)' inwffil : will call newkpt, me= ',me
!write(std_out,*)' inwffil : restart=',restart
!if(restart==2)stop
!ENDDEBUG

  if (dtset%usewvl == 0) then
     if(mkmem/=0)then
      mcg=mpw*nspinor*mband*mkmem*nsppol
      call newkpt(ceksp,cg,debug,doorth,ecut0,ecut,ecut_eff,eigen,exchn2n3d,&
    &  fill,formeig,gmet0,gmet,headform0,indkk,ab_out,ireadwf,istwfk0,istwfk,kg,kptns0,kptns,&
    &  mband,mcg,mkmem0,mkmem,mpi_enreg,&
    &  mpw0,mpw,nband0,nband,ngfft0,nkpt0,nkpt,npwarr0,npwarr,&
    &  nspinor0,nspinor,nsppol0,nsppol,nsym,occ,optorth,prtvol,&
    &  restart,rprimd,sppoldbl,symafm,symrel,tnons,unkg,wff1,wffnow)
     else
      call newkpt(ceksp,cg_disk,debug,doorth,ecut0,ecut,ecut_eff,eigen,exchn2n3d,&
    &  fill,formeig,gmet0,gmet,headform0,indkk,ab_out,ireadwf,istwfk0,istwfk,kg,kptns0,kptns,&
    &  mband,mcg_disk,mkmem0,mkmem,mpi_enreg,&
    &  mpw0,mpw,nband0,nband,ngfft0,nkpt0,nkpt,npwarr0,npwarr,&
    &  nspinor0,nspinor,nsppol0,nsppol,nsym,occ,optorth,prtvol,&
    &  restart,rprimd,sppoldbl,symafm,symrel,tnons,unkg,wff1,wffnow)
     end if
  end if
!****************************************************************************

 deallocate(indkk,istwfk0,kptns0,nband0,npwarr0)

!If disk file was created, tell its name, also deallocate the big array
 if (mkmem==0) then
  write(message, '(a,a)' ) &
&  ' inwffil: copy wf on disk file ',trim(wft1nm)
  call wrtout(06,  message,'PERS')
  deallocate(cg_disk)
 end if

 write(message,*)
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')

 call timab(19,2,tsec)

!DEBUG
!write(std_out,*)' inwffil: exit, write eigen '
!do ikpt=1,nkpt
! do iband=1,mband
!  write(std_out,*)'ikpt,iband,eigen',ikpt,iband,eigen(iband+(ikpt-1)*mband)
! end do
!end do
!if(formeig==1)stop
!ENDDEBUG

end subroutine inwffil
!!***
