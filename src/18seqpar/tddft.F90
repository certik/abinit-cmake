!{\src2tex{textfont=tt}}
!!****f* ABINIT/tddft
!!    
!! NAME
!! tddft
!!
!! FUNCTION
!! Compute the excitation energies within TDLDA
!! from input wavefunctions, eigenenergies, and band occupations.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG, JYR, MB, MBELAND, SHAMEL)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=wf in G space
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  etotal=total energy of the ground-state (Ha)
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kxc(nfft,nkxc)=exchange-correlation kernel
!!  mband=maximum number of bands
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of the dielectric matrix
!!  mkmem=maximum number of k points in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum allowed value for npw
!!  nfft=(effective) number of FFT grid points (for this processor)
!!       WARNING about parallelization: see below
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!                see ~abinit/doc/input_variables/vargs.htm#ngfft

!!  nkpt=number of k points
!!  nkxc=second dimension of the array kxc (see rhohxc for a description)
!!  npwarr(nkpt)=number of planewaves at each k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(mband*nkpt*nsppol)=
!!          occupation numbers for each band (usually 2.0) at each k point
!!  ucvol=unit cell volume (Bohr**3)
!!  wffnew=unit number for current wf disk file
!!
!! OUTPUT
!!  (only writing)
!!
!! WARNING:
!! This routine should not be parallelized on space for the time being,
!!    because the already existing parallelisation is not the usual one, found
!!    in the majority of ABINIT routines.
!! 
!! NOTES
!! * Only accept nspinor=1, nsppol=1, nkpt=1 (Gamma point), and occopt<3
!!   (insulating occupation numbers).
!!   It is expected to make it work for nsppol=2 in the future.
!!
!! * For the oscillator strengths, see the paper
!!   ''Time-Dependent Density Functional Response Theory of Molecular
!!     systems: Theory, Computational Methods, and Functionals'', by M.E. Casida,
!!   in Recent Developments and Applications of Modern Density Functional
!!   Theory, edited by J.M. Seminario (Elsevier, Amsterdam, 1996).
!!
!! TODO
!! The temporary files should be deleted at the end of this routine
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      chpev,fourdp,fourwf,hartre,hdr_skip,int2char4,leave_new,matr3inv
!!      mpi_barrier,mpi_bcast,mpi_gatherv,mpi_reduce,mpi_scatterv,rdnpw
!!      read_wfrspa,rwwf,sort_dp,sphereboundary,timab,wrtout,xcomm_world
!!      xme_whoiam,xproc_max,xsum_master,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine tddft(cg,dtfil,dtset,eigen,etotal,gmet,gprimd,gsqcut,&
&  kg,kxc,mband,mgfftdiel,mkmem,mpi_enreg,mpw,nfft,ngfftdiel,nkpt,nkxc,&
&  npwarr,nspinor,nsppol,occ,ucvol,wffnew)

 use defs_basis
 use defs_datatypes
#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_13io_mpi
 use interfaces_13xc
 use interfaces_14iowfdenpot
 use interfaces_lib00numeric
 use interfaces_lib01hidempi
 use interfaces_linalg
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer, intent(in) :: mband,mgfftdiel,mkmem,mpw,nfft,nkpt,nkxc,nsppol
 integer, intent(inout) :: nspinor
 real(dp), intent(in) :: etotal,gsqcut,ucvol
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg
 type(wffile_type), intent(inout) :: wffnew
 integer, intent(in) :: kg(3,mpw*mkmem),ngfftdiel(18),npwarr(nkpt)
 real(dp), intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),eigen(mband*nkpt*nsppol)
 real(dp), intent(in) :: gmet(3,3),gprimd(3,3),kxc(nfft,nkxc),occ(mband*nkpt*nsppol)
 
!Local variables-------------------------------
 integer,parameter :: nexcitout=20
 integer :: cplex,i,i1,i2,i3,iband,idir,ier,ierr,iexcit,iexcit1,iexcit2,ifft
 integer :: old_iexcit,ii,jj,isppol,jsppol,isppol1,isppol2,isppol_l,isppol_n
 integer :: isppol_n1,isppol_n2,iocc_n1,iocc_n2,iunocc_n1,iunocc_n2
 integer :: ikpt,index,iocc,iocc1,iocc2,iocc_l,iocc_n,ipw
 integer :: istwf_k,iunocc,iunocc1,iunocc2,iunocc_l,iunocc_n,istate,jexcit
 integer :: jexcit1,jexcit_cbase,master,mcg_disk,me_loc,mu,nband1
 integer :: nband_k(nsppol), nband_occ(nsppol), nband_unocc(nsppol)
 integer :: nstate_k, nstate_occ, nstate_unocc, nexcit_pol(nsppol)
 integer :: nstate_win,ndiel,ndiel1,ndiel2,ndiel3,ndiel4
 integer :: ndiel5,ndiel6,nexcit,nexcit_max,nexcit_win,nfftdiel,nlargest,nnext
 integer :: nnext1,nnext2
 integer :: nproc_loc,npw1,npw_k,pole_approx,sing_trip,spaceComm,tim_fourwf
 integer :: tim_rwwf,save_accesswff
 integer :: rec,recl,idummy,jdummy
 real(dp) :: buffer,buffer_inv,cim,cre,diffeig,eigbnd,eigunocc,emax_win
 real(dp) :: f_sin_trip(2/nsppol),factor,ff,flargest,fnext,fr_invsquare,fr_power
 real(dp) :: fnext1,fnext2
 real(dp) :: normint,product,root1,root2,saa,sab,sbb,si1,si2,si3,sum_exc
 real(dp) :: sum_exc_i,sum_exc_r,sumsqr,sumx,sumy,sumz
 real(dp) :: sum_kernel(2/nsppol)
 real(dp) :: theta,thppi,weight,xx,yy
 logical :: am_master,file_exist
 logical, allocatable :: done_excit(:,:),done_sexc(:),done_sexc2(:)
 character(len=fnlen) :: fname_tdexcit,fname_wf
 character(len=4) :: tag
 character(len=18) :: chain1,chain2
 character(len=500) :: message
 integer,allocatable :: flag_state_win(:),gbound(:,:),indarr(:),index_state(:)
 integer,allocatable :: kg_dum(:,:),kg_k(:,:)
 integer,allocatable :: excit_coords(:,:)
 integer :: count_to_do, count, displ, countmax, displmax
 integer :: ijexcit, ijexcit2, sendcount, f_sing_trip(2/nsppol)
 real(dp) :: sendbuf(5-nsppol)
 real(dp) :: cauchy(7),poscart(3),qphon(3),rprimd(3,3),tsec(2),dummy(2,1)
 real(dp),allocatable :: cg_disk(:,:,:),cwavef(:,:),eexcit(:)
 real(dp),allocatable :: eexcit2(:),eig_dum(:)
 real(dp),allocatable :: matr(:),occ_dum(:)
 real(dp),allocatable :: kxc_for_tddft(:,:,:,:,:,:),omega_tddft_casida(:,:,:,:,:,:,:)
 real(dp),allocatable :: osc_str(:,:),pos(:,:),rhoaug(:,:,:),rhog(:,:)
 real(dp),allocatable :: sexc(:,:),sqrtks(:),vec(:,:,:),vhartr(:),wfprod(:,:,:)
 real(dp) :: omega_tddft_casida_dummy(2/nsppol)
 real(dp),allocatable :: wfraug(:,:,:,:),wfrspa(:,:,:,:),work(:),zhpev1(:,:),intel_compiler(:,:,:)
 real(dp),allocatable :: zhpev2(:)
!no_abirules
#          if defined MPI
           integer :: iproc,jband_loc,nband_per_proc,nbuf
           integer,allocatable :: counts(:),displs(:),recvcounts(:)
           real(dp), allocatable :: recvbuf(:,:)
#          endif 
 
! *************************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zhpev
#endif

!DEBUG
!write(6,*)' tddft : enter '
!call flush(6)
!if(mkmem==0)stop
!ENDDEBUG

!Init mpi_comm
!BEGIN TF_CHANGES
 call xcomm_world(mpi_enreg,spaceComm)
!END TF_CHANGES

 am_master=.true.
 master = 0
!Init ntot proc max
 call xproc_max(nproc_loc,ierr)

!Define who i am
 call xme_whoiam(me_loc)


#if defined MPI

           if (me_loc/=0) then
            am_master=.FALSE.
           endif
           write(message, '(a,i3,a)' ) ' TDDFT ',nproc_loc,' CPU synchronized'
           call wrtout(6,message,'COLL')
           write(message, '(a,3D12.5,a,3D12.5,a,3D12.5)' ) ' gmet ',&
&                             gmet(1,1),gmet(1,2),gmet(1,3),ch10,&
&                             gmet(2,1),gmet(2,2),gmet(2,3),ch10,&
&                             gmet(3,1),gmet(3,2),gmet(3,3)                               
           call wrtout(6,message,'COLL')
#endif 


!COMMENT these values should become arguments
!the two first define the energy window

!XG 020209 : Jean-Yves, note that the definition 
!of td_maxene is in Hartree ...
 emax_win = dtset%td_maxene
!nexcit_max is defined later
!nexcit_max = 600   
 
 call timab(95,1,tsec)

 istwf_k=dtset%istwfk(1)

 if(nkpt/=1 .or. &
&    abs(dtset%kptns(1,1))+abs(dtset%kptns(2,1))+abs(dtset%kptns(3,1))>1.0d-6 )then
  write(message, '(a,a,a,a,a,a,a,a,i4,a,3es14.6,a,a,a,a,a)' )ch10,&
&  ' tddft : ERROR -',ch10,&
&  '  The computation of excited states using TDDFT is only allowed',ch10,&
&  '  with nkpt=1, kpt=(0 0 0), but the following values are input:',ch10,&
&  '  nkpt=',nkpt,', kpt=',dtset%kptns(1:3,1),'.',ch10,&
&  '  Action : in the input file, set nkpt to 1 and kpt to 0 0 0 ,',ch10,&
&  '  or change iscf.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 endif

 if(nspinor/=1)then
  write(message, '(a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&  ' tddft :  ERROR -',ch10,&
&  '  The computation of excited states using TDDFT is restricted',ch10,&
&  '  for the time being to nspinor=1, while input nspinor=2.',ch10,&
&  '  Action : if you want to compute excited states within TDDFT,',ch10,&
&  '  set nsppol to 1 in the input file. Otherwise, do not use iscf=-1.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 endif


 if(nsppol==2 .and. (dtset%ixc==22 .or. dtset%ixc==20))then
  write(message, '(a,a,a,a,a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&  ' tddft :  ERROR -',ch10,&
&  '  The computation of excited states using TDDFT in the spin',ch10,&
&  '  polarized case for the time being cannot be used with ixc=20',ch10,&
&  '  or ixc=22',ch10,&
&  '  Action : if you want to compute excited states within TDDFT,',ch10,&
&  '  set ixc different from 20 or 22. Otherwise, do not use iscf=-1',ch10,&
&  '  with nsppol=2.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 endif


 if(dtset%occopt>2)then
  write(message, '(a,a,a,a,a,a,i2,a,a,a,a,a)' ) ch10,&
&  ' tddft :  ERROR -',ch10,&
&  '  The computation of excited states using TDDFT is only allowed',ch10,&
&  '  with occopt=0, 1, or 2, while input occopt=',dtset%occopt,'.',ch10,&
&  '  Action : if you want to compute excited states within TDDFT,',ch10,&
&  '  set occopt=0, 1, or 2 in the input file. Otherwise, do not use iscf=-1.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 endif

!Examine the occupation numbers, and determine the number of 
!occupied and unoccupied states and band.
!States are numerated as usual in Abinit, before all spin up band
!and after all spin down bands.
!Note that if nsppol==1 nstate=nband_k
 do isppol=1,nsppol
  nband_k(isppol)=dtset%nband(isppol)
  nband_occ(isppol)=0
  do iband=1,nband_k(isppol)
   if(abs(occ(iband+(isppol-1)*nband_k(1))-two/nsppol)<tol6)  &
&             nband_occ(isppol)=nband_occ(isppol)+1
  enddo
  nband_unocc(isppol)=nband_k(isppol)-nband_occ(isppol)
 !next line make no sense if spin flip is taken into account
  nexcit_pol(isppol)=nband_occ(isppol)*nband_unocc(isppol)
 enddo
 nstate_k=nband_k(1)+(nsppol-1)*nband_k(nsppol)
 nstate_occ=nband_occ(1)+(nsppol-1)*nband_occ(nsppol)
 nstate_unocc=nstate_k-nstate_occ
 !next line to be changed if spin fli is taken into account
 nexcit=nexcit_pol(1)+(nsppol-1)*nexcit_pol(nsppol)   
 
!number of plane wave (does it work even for nsppol=2 ??)
 npw_k=npwarr(1)

!mux number of excitations that is taken into account 
 if(dtset%td_mexcit==0)then
  nexcit_max=nexcit
 else
  nexcit_max =dtset%td_mexcit
 endif

!DEBUG
! write(6,*) nband_occ(1),nband_unocc(1)
! write(6,*) nband_occ(nsppol),nband_unocc(nsppol)
!END DEBUG


 if(nsppol==1)then
  write(message, '(a,a,a,a,i4,a,i4,a,a,i4,a,a,a,i6,a)' )ch10,&
& ' *** TDDFT : computation of excited states *** ',ch10,&
& ' Splitting of',dtset%nband(1),' states in',nband_occ(1),' occupied states,',&
&     ' and',nband_unocc(1),' unoccupied states,',ch10,&
& ' giving',nexcit,' excitations.'
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')
 else
  write(message, '(a,a,a,a,i4,a,i4,a,a,i4,a,a,a,i6,a,a,a)' )ch10,&
&  ' *** TDDFT : computation of excited states *** ',ch10,&
&  ' Splitting of',nstate_k,' states in',nstate_occ,' occupied states,',&
&      ' and',nstate_unocc,' unoccupied states,',ch10,&
&  ' giving',nexcit,' excitations. Note that spin flip is not possible actually.',ch10,&
&  ' So the number of excitation is the half of the product of the number of state' 
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')
 endif

!Allocate the matrices to be diagonalized. 
!Use a simple storage mode, to be improved in the future.
 ii=max(nband_occ(1),nband_occ(nsppol))
 jj=max(nband_unocc(1),nband_unocc(nsppol))
 allocate(omega_tddft_casida(ii,jj,nsppol,ii,jj,nsppol,2/nsppol))
 allocate(eexcit(nexcit),sqrtks(nexcit)) 
 allocate(flag_state_win(nstate_k))
 omega_tddft_casida(:,:,:,:,:,:,:)=zero

 
!Fill the diagonal elements with square of differences of KS eigenvalues
!(also not very efficient, but OK for the present first coding)
!Also compute the square root of Kohn-Sham eigenvalue differences
 do isppol=1,nsppol
  do iunocc=1,nband_unocc(isppol)
   eigunocc=eigen(iunocc+nband_occ(isppol)+(isppol-1)*nband_k(1))
   do iocc=1,nband_occ(isppol)
    iexcit=iocc+(isppol-1)*nexcit_pol(1)+nband_occ(isppol)*(iunocc-1)
    diffeig=eigunocc-eigen(iocc+(isppol-1)*nband_k(1))
    do sing_trip=1,2/nsppol
     omega_tddft_casida(iocc,iunocc,isppol,iocc,iunocc,isppol,sing_trip)=diffeig**2
    enddo
    eexcit(iexcit)=diffeig
    sqrtks(iexcit)=sqrt(diffeig)
   enddo
  enddo
 enddo



!Sort the excitation energies : note that the array eexcit is reordered
 allocate(indarr(nexcit))
 indarr(:)=(/ (ii,ii=1,nexcit) /)
 call sort_dp(nexcit,eexcit,indarr,tol14)

!Determine an energy window for the excitations
!to take into account. This is necessary for large systems
  
 nexcit_win = 0
 do iexcit = 1, nexcit
  if ((eexcit(iexcit) < emax_win ).and.(nexcit_win < nexcit_max)) then
   nexcit_win = nexcit_win + 1
   
   !DEBUG
   !write(message,'(a,F12.5,a,a,i2,a,a,i2)') 'excitation energy:', eexcit(indarr(iexcit)),ch10, &
   !&                                        'excitation number:', indarr(iexcit),ch10,         &
   !&                                        'nexcit_win:       ', nexcit_win
   !call wrtout(6,message,'COLL')
   !END DEBUG
   
  endif 
 enddo

!identification of the bands contributing to the 
!nexcit_win  excitations within the window


 nstate_win = 0
 flag_state_win(:) = 0
 do iexcit = 1, nexcit_win
  iexcit1 = indarr(iexcit)
  isppol1 = min((iexcit1-1)/nexcit_pol(1) +1,2)
  iunocc1 = (iexcit1-(isppol1-1)*nexcit_pol(1)-1)/nband_occ(isppol1)+1
  iocc1   = iexcit1-(isppol1-1)*nexcit_pol(1)-(iunocc1-1)*nband_occ(isppol1)
  if (flag_state_win(nband_occ(isppol1)+(isppol1-1)*nband_k(1)+iunocc1)==0) &
&     flag_state_win(nband_occ(isppol1)+(isppol1-1)*nband_k(1)+iunocc1) =1
  if (flag_state_win(iocc1+(isppol1-1)*nband_k(1))==0) &
&     flag_state_win(iocc1+(isppol1-1)*nband_k(1)) =1 
  !DEBUG
  ! write(message,'(a,i2,a,a,i2,i2,a,a,i2,a,a,i2)') 'isppol:', isppol1,ch10, &
  !  &                                       'iocc,iunocc:', iocc1,iunocc1,ch10,         &
  !  &                                       'flag_state_win:', flag_state_win(iocc1+(isppol1-1)*nband_k(1)),  &
  !  &                               ch10,   'flag_state_win:', flag_state_win(nband_occ(isppol1)+(isppol1-1)*nband_k(1)+iunocc1)
  ! call wrtout(6,message,'COLL')
  !END DEBUG

 enddo
 
 do isppol=1,nsppol
  do iband=1,nband_k(isppol)
    nstate_win=nstate_win+flag_state_win(iband+(isppol-1)*nband_k(1))
  enddo
 enddo


 write(message,'(a,a,i5)') ch10,                     &
&  'Nr of states to Fourier transform : ',nstate_win
 call wrtout(6,message,'COLL')

 ndiel1=ngfftdiel(1) ; ndiel2=ngfftdiel(2) ; ndiel3=ngfftdiel(3)
!ndiel4,ndiel5,ndiel6 are FFT dimensions, modified to avoid cache trashing
 ndiel4=ngfftdiel(4) ; ndiel5=ngfftdiel(5) ; ndiel6=ngfftdiel(6)

!The evaluation of integrals, later, needs the following factor
 normint=one/(ucvol*dble(ndiel1*ndiel2*ndiel3))

!Setup the positions in real space for later integration
 call matr3inv(gprimd,rprimd)

 allocate(pos(max(ndiel1,ndiel2,ndiel3),3))

!Select the reduced position of the point with respect to the box center,
!in the interval ]-0.5,0.5]. 
 buffer=0.05_dp ; buffer_inv=one/buffer
 do idir=1,3
  if(idir==1)ndiel=ndiel1
  if(idir==2)ndiel=ndiel2
  if(idir==3)ndiel=ndiel3
  do ii=1,ndiel
!  dtset%boxcenter(3)=reduced coordinates of the center of the box,
!  in view of the computation of the oscillator strength
   pos(ii,idir)=(ii-1)/(one*ndiel)-dtset%boxcenter(idir)
   pos(ii,idir)=pos(ii,idir)-nint(pos(ii,idir)-tol12)
!  The linear behaviour is cut-off when one becomes
!  close to the boundaries : the buffer allows to match smoothly
!  one side of the cell to the other. This is important
!  to get rid of small breakings of symmetry, that are
!  confusing in accurate tests
   if(abs(pos(ii,idir))>half-buffer)then
!   xx is always positive, and goes linearly from 1 to 0
!   in the buffer region
    xx=(half-abs(pos(ii,idir)))*buffer_inv
!   The cut-off is applied to pos(:,:)
    pos(ii,idir)=pos(ii,idir)*xx*(two-xx)
!DEBUG
!   if (idir==1)then
!    write(6,'(i2)') ndiel
!    write(6,'(a,i2,a,F12.5,F12.5)')'idiel : ',ii,'   x : ',pos(ii,idir),&
!&    dtset%boxcenter(idir)
!   endif 
!ENDDEBUG
   endif
  enddo ! ii
 enddo ! idir

!  need to run in MPI I/O case   
  if (wffnew%accesswff == 1 ) then
   save_accesswff=wffnew%accesswff
   wffnew%accesswff = 0
  else
!  Do not store value but set to have save_accesswff /= 1
   save_accesswff = 0
  end if 

 if(am_master)then

!-----------------------------------------------------------
!Do i/o as needed, mkmem==0 means wf and kg info on disk file
!The disk access is only done by master...

  allocate(gbound(2*mgfftdiel+8,2),kg_k(3,npw_k))
  
  if (mkmem==0) then
 
 ! Skip wffnew header
   call hdr_skip(wffnew,ierr)
 ! rewind the kpgsph data file on unit unkg
   rewind (dtfil%unkg)
   mcg_disk=mpw*nspinor*mband
   allocate(cg_disk(2,mcg_disk,nsppol))
 
   ikpt=1 ; isppol=1
!   do isppol=1,nsppol
   call rdnpw(ikpt,isppol,nband_k(isppol),npw_k,nspinor,0,dtfil%unkg)
!   enddo
!   Read k+g data
   read (dtfil%unkg) kg_k(:,1:npw_k)
!   do isppol=1,nsppol
   call sphereboundary(gbound,istwf_k,kg_k,mgfftdiel,npw_k)
!  enddo

  else   !so, mkmem/=0

   ikpt=1
! Only one k point
!  do isppol=1,nsppol
    kg_k(:,1:npw_k)=kg(:,1:npw_k)
    call sphereboundary(gbound,istwf_k,kg_k,mgfftdiel,npw_k)
!  enddo

  end if !mkmem==0

  if(mkmem==0)then
!  Read the wavefunction block for ikpt,isppol
   tim_rwwf=8
   allocate(eig_dum(mband),kg_dum(3,0),occ_dum(mband))
   do isppol=1,nsppol
    call rwwf(cg_disk(:,:,isppol),eig_dum,0,0,0,ikpt,isppol,kg_dum,mband,mcg_disk,mpi_enreg,nband_k(isppol),nband_k(isppol),&
&   npw_k,nspinor,occ_dum,-2,0,tim_rwwf,wffnew)
   enddo
   deallocate(eig_dum,kg_dum,occ_dum)

  end if !mkmem==0

 endif ! am_master

!  need to run in MPI I/O case   
 if ( save_accesswff == 1 )  wffnew%accesswff = 1

!    end of read
!-------------------------------------------------------------

!DEBUG
!write(6,*) mkmem
!write(6,'(a)') 'After reading the wavefunction'      
!call flush(6)
!call wrtout(6,message,'COLL')
!ENDDEBUG

!Use a simple implementation for the computation of the kernel elements
 if (am_master) then
  allocate(cwavef(2,mpw))
  allocate(rhoaug(ndiel4,ndiel5,ndiel6),wfraug(2,ndiel4,ndiel5,ndiel6))
 endif
 allocate(index_state(nstate_k)) 

!If mkmem is not 0, all real-space states are kept in memory
!                0, only 2   are in memory at the same time
 if (mkmem/=0) then
  allocate(wfrspa(ndiel4,ndiel5,ndiel6,nstate_win))
 else
  allocate(wfrspa(ndiel4,ndiel5,ndiel6,2))
 endif

!DEBUG
! write(message,'(a)') 'After allocating wfrspa'
! call wrtout(6,message,'COLL')
!ENDDEBUG

 weight=zero

!Generate states in real space, only for states contributing to excitations in window
 istate=0

 do isppol=1,nsppol
  do iband=1,nband_k(isppol)
   
    if(flag_state_win(iband+(isppol-1)*nband_k(1)) == 1) then
     istate=istate+1
     index_state(iband+(isppol-1)*nband_k(1))=istate

     if (am_master) then
    !Obtain Fourier transform in fft box
     if(mkmem/=0)cwavef(:,1:npw_k)=cg(:,1+(iband-1)*npw_k+(isppol-1)* &
     & (npw_k*nband_k(1)) : iband*npw_k+(isppol-1)*(npw_k*nband_k(1)))
     if(mkmem==0)cwavef(:,1:npw_k)=cg_disk(:,1+(iband-1)*npw_k        &
     & : iband*npw_k,isppol)

!DEBUG
!   write(6,*)' iband : ',iband, ' isppol', isppol, '  -> index ', &
!&            istate,index_state(iband+(isppol-1)*nband_k(1))
!ENDDEBUG

     tim_fourwf=14
!   This call should be made by master, and then the results be sent to the other procs

     call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
      & istwf_k,kg_k,kg_k,mgfftdiel,mpi_enreg,1,ngfftdiel,npw_k,1,ndiel4,ndiel5,ndiel6,&
      & 0,dtset%paral_kgb,tim_fourwf,weight)

!DEBUG
!   write(6, '(a,i5)')' After Fourier proc ',me_loc
!ENDDEBUG

!   Fix the phase, and checks that the wavefunction is real
!   (should be merged with routine fxphas)
     saa=zero ; sab=zero ; sbb=zero
     do i3=1,ndiel3
      do i2=1,ndiel2
       do i1=1,ndiel1
        saa=saa+wfraug(1,i1,i2,i3)**2
        sbb=sbb+wfraug(2,i1,i2,i3)**2
        sab=sab+wfraug(1,i1,i2,i3)*wfraug(2,i1,i2,i3)
       enddo
      enddo
     enddo

     if(sbb>5.0d-9)then

      write(message, '(a,a,a,a,a,a,es20.10,a,i4,a,i2,a)' )ch10,&
&      ' tddft : WARNING -',ch10,&
&      '  The imaginary part of wavefunctions should be practically zero.',&
&      ch10,'  This is not the case, since sbb=',sbb,' for iband=',iband,  &
&      'with sppol=',1,'.'
      call wrtout(6,message,'PERS')
      if(sbb>1.0d-7)then
       call leave_new('COLL')
      endif
     endif

!   Possibility of writing to disk

     if (mkmem/=0) then
      wfrspa(:,:,:,istate)=wfraug(1,:,:,:)
     else
      call int2char4(istate,tag) 
      fname_wf=trim(dtfil%filnam_ds(5))//'_TDWF'//tag
      write(6,*)'fileout : ',fname_wf
      open(tmp_unit,file=fname_wf,form='unformatted',status='unknown')
      !DEBUG
      ! write(6, '(a)')' After opening tmp_unit '
      !END DEBUG
      write(tmp_unit)iband,isppol,eigen(iband+(isppol-1)*nband_k(1))
      !DEBUG
      ! write(6, '(a)')' After writing haed '
      !END DEBUG
!withou the following line I got a 'segmentation fault' (core dumped)
!working with the intel compiler... ???
      wfrspa(:,:,:,1)=wfraug(1,:,:,:)
      write(tmp_unit)wfrspa(:,:,:,1)
      !DEBUG
      ! write(6, '(a)')' After writing wavefunction '
      !END DEBUG
      
      close(tmp_unit)
     endif
    endif !am_master 

!  if mkmem=0 the eigenstates (in r space) with eigenvalues inside
!  the energy window have been stored to disk
!  may cause enormous disk usage, slow down computation. Use only
!  if limited memory is available. This option allows a restart.

   endif

!End loop on iband
  enddo
!End loop on nsppol  
 enddo

 if (am_master) then
  deallocate(gbound,kg_k)
  deallocate(cwavef)
  deallocate(rhoaug,wfraug)
  if(mkmem==0)deallocate(cg_disk)
 endif
 deallocate(flag_state_win)

!Case mkmem /= 0 , send wfrspa from master to world

#if defined MPI
           if (mkmem/=0) then
            nbuf=ndiel4*ndiel5*ndiel6*nstate_win
            call MPI_BCAST(wfrspa,nbuf,MPI_DOUBLE_PRECISION,master,spaceComm,ierr)
           endif
#endif


!DEBUG
!#if defined MPI
!          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!          write(message,'(a)')' after iband loop synchronization done...'
!          call  wrtout(6,message,'COLL')
!#endif
!ENDDEBUG

!Compute the xc kernel, in the form needed for the singlet or triplet
!excitation energy. 
!In the case ixc=20, kxc vanishes, but no change is made here, for simplicity.
!(ixc=20 implemented only in the not spin polyrized case)

!DEBUG
!write(6,*)' tddft : xc kernel '
!do ifft=1,nkxc,41
! write(6,*)ifft,kxc(ifft,1),kxc(ifft,2),kxc(ifft,3)
!enddo
!stop
!ENDDEBUG

 allocate(kxc_for_tddft(ndiel1,ndiel2,ndiel3,nsppol,nsppol,2/nsppol))
 if(dtset%ixc/=22)then
  do isppol=1,nsppol
   do jsppol=1,nsppol
    index=1
    do i3=1,ndiel3
     do i2=1,ndiel2
      do i1=1,ndiel1
       do sing_trip=1,2/nsppol
        kxc_for_tddft(i1,i2,i3,isppol,jsppol,sing_trip)=two/nsppol* &
&                    (kxc(index,isppol+jsppol-1)-(sing_trip-1)*kxc(index,2))
       enddo
       index=index+1
      enddo
     enddo
    enddo
   enddo
  enddo
 else
! This is for the Burke-Petersilka-Gross hybrid, with ixc=22
! However, the implementation in case of spin-polarized system should not be expected to be the correct one !
  do isppol=1,nsppol
   do jsppol=1,nsppol
    index=1
    do i3=1,ndiel3
     do i2=1,ndiel2
      do i1=1,ndiel1
       do sing_trip=1,2/nsppol
        kxc_for_tddft(i1,i2,i3,isppol,jsppol,sing_trip)=((-1)**(sing_trip+1))*kxc(index,2)
       enddo
       index=index+1
      enddo
     enddo
    enddo
   enddo
  enddo
 endif

 pole_approx=0

 allocate(excit_coords(nexcit_win**2,2))
 
!DEBUG
!write(6,*)'before first loop'
!ENDDEBUG

!0000000000000000000000000000000000000000000000000000000
!check if matrix file fname_tdexcit exists on disk
!if the file is present,calculation is a continuation 
 if (am_master) then
  fname_tdexcit=trim(dtfil%filnam_ds(5))//'_TDEXCIT'
  inquire(file=fname_tdexcit,exist=file_exist)
  ! for direct access to excitation file
  if(nsppol==1)then
   inquire(iolength=recl) omega_tddft_casida(1,1,1,1,1,1,1),omega_tddft_casida(1,1,1,1,1,1,1), &
&                        iexcit,jexcit
  else
   inquire(iolength=recl) omega_tddft_casida(1,1,1,1,1,1,1), &
&                         iexcit,jexcit
  endif
  
  open(tmp_unit2, file=fname_tdexcit,form='unformatted', recl=recl, &
&      access='DIRECT')

  allocate(done_excit(nexcit_win,nexcit_win))

  if(file_exist)then
   write(6,*)'TDDFT continues from a previous run'
   rec=0
   do iexcit=1,nexcit_win
    iexcit2 = indarr(iexcit)
    isppol2 = min((iexcit2-1)/nexcit_pol(1) +1,2)
    iunocc2 = (iexcit2-(isppol2-1)*nexcit_pol(1)-1)/nband_occ(isppol2)+1
    iocc2   = iexcit2-(isppol2-1)*nexcit_pol(1)-(iunocc2-1)*nband_occ(isppol2)
    do jexcit=1,nexcit_win
     iexcit1 = indarr(jexcit)
     isppol1 = min((iexcit1-1)/nexcit_pol(1) +1,2)
     iunocc1 = (iexcit1-(isppol1-1)*nexcit_pol(1)-1)/nband_occ(isppol1)+1
     iocc1   = iexcit1-(isppol1-1)*nexcit_pol(1)-(iunocc1-1)*nband_occ(isppol1)

     rec=rec+1 ! record of the entry in the excitation file
       if(nsppol==1)then
        read(tmp_unit2,rec=rec) omega_tddft_casida_dummy(1), &
&                               omega_tddft_casida_dummy(2), &
&                               idummy, jdummy
       else
        read(tmp_unit2,rec=rec) omega_tddft_casida_dummy(1), &
&                               idummy, jdummy
       endif
       done_excit(jexcit,iexcit)= ( idummy /= -1 .and. jdummy /= -1 ) ! if true, eigsqr_singlet and eigsqr_triplet are ok
                                                                    ! and a true is marked in the logical array done_excit
       if (done_excit(jexcit,iexcit)) then
        do sing_trip=1,2/nsppol
         omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,sing_trip)= &
&        omega_tddft_casida_dummy(sing_trip)
        enddo
       endif
      enddo
     enddo
     
  else
   write(*,*)'no excitation matrix on disk'
   write(6,*)'TDDFT starts from scratch'
   done_excit = .false. ! initialize the logical array to false
   rec=0
   do iexcit=1,nexcit_win
    do jexcit=1,nexcit_win
     rec=rec+1
     if(nsppol==1)then
      write(tmp_unit2,rec=rec) zero,zero,-1,-1
     else
      write(tmp_unit2,rec=rec) zero,-1,-1
     endif
    enddo
   enddo
  endif

! Need to list the elements to compute, taking the symmetry into account: valid only for Gamma point but this is already the case
  count_to_do=0
  do iexcit=1,nexcit_win
   do jexcit=1,iexcit
    if (.not. done_excit(jexcit,iexcit)) then
     count_to_do=count_to_do+1
     excit_coords(count_to_do,1)=iexcit
     excit_coords(count_to_do,2)=jexcit
    endif
   enddo
  enddo

  deallocate(done_excit)



#           if defined MPI
            allocate(counts(0:nproc_loc-1))
            allocate(displs(0:nproc_loc-1))
            allocate(recvcounts(0:nproc_loc-1))
            allocate(recvbuf(5-nsppol,nproc_loc-1))

!           Compute limits for load balancing
            do iproc=0,nproc_loc-1
             displs(iproc)=(iproc*count_to_do)/nproc_loc
             counts(iproc)=min(((iproc+1)*count_to_do)/nproc_loc,count_to_do)-displs(iproc)
            enddo
#           endif

 endif ! am_master


#          if defined MPI
           call MPI_BCAST(count_to_do,1,MPI_INTEGER,master,spaceComm,ierr)
#          endif

 displ=(me_loc*count_to_do)/nproc_loc
 count=min(((me_loc+1)*count_to_do)/nproc_loc,count_to_do)-displ
 displmax=((nproc_loc-1)*count_to_do)/nproc_loc
 countmax=count_to_do-displmax

 write(message,'(A,I6)') 'Maximum number of matrix elements per processor = ',countmax
 call wrtout(6,message,'COLL')


#          if defined MPI
!          Need to dispatch the elements to compute to the different processes
           call MPI_Scatterv(excit_coords(1,1),counts,displs,MPI_INTEGER,excit_coords(1,1),count, &
&                            MPI_INTEGER,0,spaceComm,ierr)
           call MPI_Scatterv(excit_coords(1,2),counts,displs,MPI_INTEGER,excit_coords(1,2),count, &
&                  MPI_INTEGER,0,spaceComm,ierr)
#          endif

 nfftdiel=ndiel1*ndiel2*ndiel3
 allocate(wfprod(ndiel1,ndiel2,ndiel3),work(nfftdiel))
 allocate(sexc(3,nexcit_win))
 allocate(done_sexc(nexcit_win))
 allocate(done_sexc2(nexcit_win))
 allocate(rhog(nfftdiel,2),vhartr(nfftdiel))

  sexc(:,:)=zero 
  done_sexc(:)=.false.

!----------------------------------------------------------
!Main double loop

 old_iexcit=0
 do ijexcit=1,countmax
  ! we really loop only through count, but we need to go through countmax
  ! to make sure that all processes execute MPI_Gatherv below
  if (ijexcit <= count) then
   iexcit=excit_coords(ijexcit,1)
   jexcit=excit_coords(ijexcit,2)

   iexcit2 = indarr(iexcit)
   isppol2 = min((iexcit2-1)/nexcit_pol(1) +1,2)
   iunocc2 = (iexcit2-(isppol2-1)*nexcit_pol(1)-1)/nband_occ(isppol2)+1
   iocc2   = iexcit2-(isppol2-1)*nexcit_pol(1)-(iunocc2-1)*nband_occ(isppol2)

   iexcit1 = indarr(jexcit)
   isppol1 = min((iexcit1-1)/nexcit_pol(1) +1,2)
   iunocc1 = (iexcit1-(isppol1-1)*nexcit_pol(1)-1)/nband_occ(isppol1)+1
   iocc1   = iexcit1-(isppol1-1)*nexcit_pol(1)-(iunocc1-1)*nband_occ(isppol1) 

   if (old_iexcit /= iexcit) then
!   We start a new column of the matrix
!DEBUG
!   write(message,'(a,i5,a,i3)')'treating  iexcit =  ',iexcit,&
!                              &'  with proc ',me_loc
!   call wrtout(6,message,'PERS') 
!ENDDEBUG



!   If mkmem ==0, the real-space filled and empty states are read from disk  
    if (mkmem==0) then
!    Read occupied band using internal subroutine
     call read_wfrspa(index_state(iocc2+(isppol2-1)*nband_k(1)), &
                      &dtfil,eigbnd,iband,isppol,1,ndiel4,ndiel5,ndiel6,wfrspa(:,:,:,1))
!    Read unoccupied band
     call read_wfrspa(index_state(iunocc2+nband_occ(isppol2)+(isppol2-1)*nband_k(1)),&
                      &dtfil,eigbnd,iband,isppol,2,ndiel4,ndiel5,ndiel6,wfrspa(:,:,:,2))
    endif

!DEBUG
!   write(message,'(a,i3)')'Multiplicating phi1 phi2, on proc ',me_loc
!   call wrtout(6,message,'PERS')
!ENDDEBUG 

    ifft=1
    do i3=1,ndiel3
     do i2=1,ndiel2
      do i1=1,ndiel1
       if (mkmem/=0) then
        wfprod(i1,i2,i3)=wfrspa(i1,i2,i3,index_state(iocc2+(isppol2-1)*nband_k(1))) &
&                       *wfrspa(i1,i2,i3,index_state(iunocc2+nband_occ(isppol2)+(isppol2-1)*nband_k(1)))
       else
        wfprod(i1,i2,i3)=wfrspa(i1,i2,i3,1)*wfrspa(i1,i2,i3,2)
       endif
       work(ifft)=wfprod(i1,i2,i3)
       ifft=ifft+1
      enddo
     enddo
    enddo
    if (jexcit == 1) then
     do i3=1,ndiel3
      do i2=1,ndiel2
       do i1=1,ndiel1
        do idir=1,3
         poscart(idir)=rprimd(idir,1)*pos(i1,1)+&
&                      rprimd(idir,2)*pos(i2,2)+&
&                      rprimd(idir,3)*pos(i3,3)
         sexc(idir,iexcit)=sexc(idir,iexcit)+poscart(idir)*wfprod(i1,i2,i3)
        enddo
       enddo
      enddo
      done_sexc(iexcit)=.true.
     enddo
    endif

!   For the singlet correction, must compute the hartre potential created
!   by the product of wavefunctions
    cplex=1 ; qphon(:)=zero

!DEBUG
!   write(message,'(a,i3)')'Before Fourdp, on proc ',me_loc
!   call wrtout(6,message,'PERS')
!ENDDEBUG 

    call fourdp(cplex,rhog,work,-1,mpi_enreg,nfftdiel,ngfftdiel,dtset%paral_kgb,0)

!DEBUG
!   write(message,'(a,i3)')'Before Hartree, on proc ',me_loc
!   call wrtout(6,message,'PERS')
!   write(6,*)'CPU ',me_loc,ch10,&
!&            '   cplex : ',cplex,ch10,&
!&            '   gmet(3,3)  : ',gmet(3,3),ch10,&
!&            '   gsqcut : ',gsqcut,ch10,&
!&            '   qphon(1) : ',qphon(1),ch10,&
!&            '   rhog(1,1) :,',rhog(1,1),ch10,&
!&            '   vhartr(1) :,',vhartr(1)  
!ENDDEBUG 

    call hartre(cplex,gmet,gsqcut,0,mpi_enreg,nfftdiel,ngfftdiel,dtset%paral_kgb,qphon,rhog,vhartr)

!DEBUG
!   write(message,'(a,i3)')'After Hartree, on proc ',me_loc
!   call wrtout(6,message,'PERS')
!ENDDEBUG
   endif
   old_iexcit=iexcit

!DEBUG
!  write(6,*)'  treating  iexcit =  ',jexcit
!  write(6,*)'   indarr(iexcit) =',iexcit1,iocc1,iunocc1
!  write(6,*)'   index ',index_state(iocc1+(isppol1-1)*nband_k(1)), &
!&                       index_state(iunocc1+nband_occ(isppol1)+(isppol1-1)*nband_k(1))
!ENDDEBUG

!  if mkmem==0 the filled and empty states for excitation iexcit1 are read from disk  
   if (mkmem==0) then
!   Read occupied and unoccupied state from disk using internal
!   subroutine
    call read_wfrspa(index_state(iocc1+(isppol1-1)*nband_k(1)),&
                     &dtfil,eigbnd,iband,isppol,1,ndiel4,ndiel5,ndiel6,wfrspa(:,:,:,1))
    call read_wfrspa(index_state(iunocc1+nband_occ(isppol1)+(isppol1-1)*nband_k(1)),&
                     &dtfil,eigbnd,iband,isppol,2,ndiel4,ndiel5,ndiel6,wfrspa(:,:,:,2))
   endif !mkmem==0

   if(pole_approx==0 .or. (iunocc1==iunocc2 .and. iocc1==iocc2 .and. isppol1==isppol2))then
    sum_kernel(:)=zero
    f_sing_trip(1)=two/nsppol
    if(nsppol==1) f_sing_trip(2)=zero
!   For Fermi-Amaldi kxc, the xc contribution is -1/2 the Hartree contribution
!   to the triplet state. The following factors combines both contributions.
    if(dtset%ixc==20 .or. dtset%ixc==22)then
     if(nsppol==1)then
      f_sing_trip(1)= one
      f_sing_trip(2)=-one
     endif
    endif
    ifft=1
    do i3=1,ndiel3
     do i2=1,ndiel2
      do i1=1,ndiel1
       if (mkmem/=0) then
        product=wfrspa(i1,i2,i3,index_state(iocc1+(isppol1-1)*nband_k(1))) &
&              *wfrspa(i1,i2,i3,index_state(iunocc1+nband_occ(isppol1)+(isppol1-1)*nband_k(1)))
       else
        product=wfrspa(i1,i2,i3,1)*wfrspa(i1,i2,i3,2)
       endif
       do sing_trip=1,2/nsppol
        sum_kernel(sing_trip)=sum_kernel(sing_trip)+&
&       product*(f_sing_trip(sing_trip)*vhartr(ifft)+kxc_for_tddft(i1,i2,i3,isppol1,isppol2,sing_trip) &
               *wfprod(i1,i2,i3))
       enddo
       ifft=ifft+1
      enddo ! i1
     enddo ! i2
    enddo ! i3

!   The factor two is coherent with the formulas of Vasiliev et al
    factor=two*sqrtks(iexcit1)*sqrtks(iexcit2)*normint
    do sing_trip=1,2/nsppol
     omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,sing_trip)=   &
&     omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,sing_trip)+  &
&     factor*sum_kernel(sing_trip)
    enddo
     
!  End condition of being diagonal element if pole approximation
   endif

!  Continue writing excitation matrix 
   if (am_master) then
    ! the master writes its results to disk
    if(nsppol==1)then
     write(tmp_unit2, rec=(iexcit-1)*nexcit_win+jexcit ) &
&     omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1),&
&     omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,2),&
&     iexcit, jexcit
    else
     write(tmp_unit2, rec=(iexcit-1)*nexcit_win+jexcit ) &
&     omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1), &
&     iexcit, jexcit
    endif
   
!DEBUG
!   if(nsppol==1)then
!    write(6,*)'singlet: ',omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1),&
!&             'iexcit: ',iexcit,'jexcit :',jexcit
!    write(6,*)'triplet: ',omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,2),&
!&             'iexcit: ',iexcit,'jexcit :',jexcit
!   else
!    write(6,*)'excitation: ',omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1),&
!&             'iexcit: ',iexcit,'jexcit :',jexcit
!   endif
!ENDDEBUG

    sendcount=0
   else
    sendcount=5-nsppol

    if(nsppol==1)then
     sendbuf=(/ omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1), &
&               omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,2), &
&               real(iexcit,dp), real(jexcit,dp) /)
    else
     sendbuf=(/ omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1), &
&               real(iexcit,dp), real(jexcit,dp) /)
    endif
   endif ! am_master
  else
   ! ijexcit > count

   ! done with local work, so send message of zero length
   sendcount=0

  endif ! ijexcit <= count

#           if defined MPI
            if (am_master) then

!            Compute displacements and counts for the gathering of the results
             displs(0)=0
             recvcounts(0)=0
             do iproc=1,nproc_loc-1
              recvcounts(iproc)=min(((iproc+1)*count_to_do)/nproc_loc,count_to_do)-(iproc*count_to_do)/nproc_loc
              if (recvcounts(iproc) < countmax .and. ijexcit==countmax) then
               recvcounts(iproc)=0
              else
               recvcounts(iproc)=5-nsppol
              endif
              displs(iproc)=displs(iproc-1)+recvcounts(iproc-1)
             enddo
            endif

            !***********************************************
            !***** I have to ask about that ****************
            !***********************************************
            call MPI_Gatherv(sendbuf,sendcount,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs, &
&                            MPI_DOUBLE_PRECISION,0,spaceComm,ierr)

            if (am_master) then

             ! Extract eigsqr_singlet, eigsqr_triplet, iexcit, jexcit from receive buffer and
             ! write to file
             do ijexcit2=1,sum(recvcounts)/(4-nsppol)
              iexcit=int(recvbuf(4-nsppol,ijexcit2))
              jexcit=int(recvbuf(5-nsppol,ijexcit2))

              iexcit2 = indarr(iexcit)
              isppol2 = min((iexcit2-1)/nexcit_pol(1) +1,2)
              iunocc2 = (iexcit2-(isppol2-1)*nexcit_pol(1)-1)/nband_occ(isppol2)+1
              iocc2   = iexcit2-(isppol2-1)*nexcit_pol(1)-(iunocc2-1)*nband_occ(isppol2)

              iexcit1 = indarr(jexcit)
              isppol1 = min((iexcit1-1)/nexcit_pol(1) +1,2)
              iunocc1 = (iexcit1-(isppol1-1)*nexcit_pol(1)-1)/nband_occ(isppol1)+1
              iocc1   = iexcit1-(isppol1-1)*nexcit_pol(1)-(iunocc1-1)*nband_occ(isppol1)

              do sing_trip=1,2/nsppol
               omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,sing_trip)= &
&                                 recvbuf(sing_trip,ijexcit2)
              enddo

               if(nsppol==1)then
                write(tmp_unit2, rec=(iexcit-1)*nexcit_win+jexcit ) &
&                 omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1), &
&                 omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,2), &
&                 iexcit, jexcit
               else
                write(tmp_unit2, rec=(iexcit-1)*nexcit_win+jexcit ) &
                  omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1), &
&                 iexcit, jexcit
               endif
!DEBUG
!             if(nsppol==1)then
!              write(6,*)'singlet: ',omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1), &
!&                       'iexcit: ',iexcit,'jexcit :',jexcit
!              write(6,*)'triplet: ',omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,2),
!&                       'iexcit: ',iexcit,'jexcit :',jexcit
!             else
!              write(6,*)'excitation: ',omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,1), &
!&                        'iexcit: ',iexcit,'jexcit :',jexcit
!             endif
!ENDDEBUG
             
             enddo
       
            endif
#           endif

!End indices loops
 enddo ! ijexcit

!End of the main double loop 
!--------------------------------------------------------------------


#          if defined MPI
!          sexc needs to be summed here since it used only by master
           call MPI_BARRIER(spaceComm,ierr)
!          call xsum_master(sexc,master,spaceComm,ierr) ! Does not work on some machines
           call xsum_mpi(sexc,spaceComm,ierr)
           done_sexc2=done_sexc
           call MPI_Reduce(done_sexc2,done_sexc,nexcit_win,MPI_LOGICAL,MPI_LOR,&
&                          master,spaceComm,ierr)
#          endif


 if (am_master) then
! We compute sexc again if it was not done. Will only be executed if
! there was a restart from values read from logical unit tmp_unit2.

  do iexcit=1,nexcit_win

  !DEBUG  
  !write(6,*)'do on excitation',iexcit
  !END DEBUG
  
   if (.not.done_sexc(iexcit)) then
     iexcit2 = indarr(iexcit)
     isppol2 = min((iexcit2-1)/nexcit_pol(1) +1,2)
     iunocc2 = (iexcit1-(isppol2-1)*nexcit_pol(1)-1)/nband_occ(isppol2)+1
     iocc2   = iexcit2-(isppol2-1)*nexcit_pol(1)-(iunocc2-1)*nband_occ(isppol2)
    if (mkmem==0) then
!    Read occupied band using internal subroutine
     call read_wfrspa(index_state(iocc2+(isppol2-1)*nband_k(1)),&
                      &dtfil,eigbnd,iband,isppol,1,ndiel4,ndiel5,ndiel6,wfrspa(:,:,:,1))
!    Read unoccupied band
     call read_wfrspa(index_state(iunocc2+nband_occ(isppol2)+(isppol2-1)*nband_k(1)),&
                      &dtfil,eigbnd,iband,isppol,2,ndiel4,ndiel5,ndiel6,wfrspa(:,:,:,2))
    endif
    do i3=1,ndiel3
     do i2=1,ndiel2
      do i1=1,ndiel1
       if (mkmem/=0) then
        wfprod(i1,i2,i3)=wfrspa(i1,i2,i3,index_state(iocc2+(isppol2-1)*nband_k(1))) &
&                       *wfrspa(i1,i2,i3,index_state(iunocc2+nband_occ(isppol2)+ &
&                                                    (isppol2-1)*nband_k(1)))
       else
        wfprod(i1,i2,i3)=wfrspa(i1,i2,i3,1)*wfrspa(i1,i2,i3,2)
       endif
       do idir=1,3
        poscart(idir)=rprimd(idir,1)*pos(i1,1)+&
&                     rprimd(idir,2)*pos(i2,2)+&
&                     rprimd(idir,3)*pos(i3,3)
        sexc(idir,iexcit)=sexc(idir,iexcit)+poscart(idir)*wfprod(i1,i2,i3)
       enddo
      enddo
     enddo
    enddo
   endif
  enddo
 endif

 deallocate(work,rhog,pos)
 deallocate(vhartr,kxc_for_tddft,wfprod)
 deallocate(index_state,excit_coords)  
 deallocate(wfrspa)

!Write the first excitation energies
 write(message, '(a,a,es18.8,a,a,a,a,a,a,a,a,a)' )ch10,&
& '  Ground state total energy (Ha) :',etotal,ch10,ch10,&
& '  Kohn-Sham energy differences,',ch10,&
& '  corresponding total energies and oscillator strengths (X,Y,Z and average)-',&
& ch10,&
& '  (oscillator strengths smaller than 1.e-6 are set to zero)',ch10,&
& '  Transition  (Ha)  and   (eV)   Tot. Ene. (Ha)  Aver     XX       YY       ZZ'
 call wrtout(ab_out,message,'COLL')
 call wrtout(6,message,'COLL')

 if (am_master) then


#           if defined MPI
            deallocate(counts,displs)
            deallocate(recvbuf,recvcounts)
#           endif

  allocate(osc_str(7,nexcit))

  do iexcit=1,nexcit_win
   iexcit2 = indarr(iexcit)
   isppol = min((iexcit2-1)/nexcit_pol(1) +1,2)
   iunocc = (iexcit2-(isppol-1)*nexcit_pol(1)-1)/nband_occ(isppol)+1
   iocc   = iexcit2-(isppol-1)*nexcit_pol(1)-(iunocc-1)*nband_occ(isppol)
    
   osc_str(1,iexcit)=zero
   do idir=1,3
!   One of the factor of two comes from the spin degeneracy,
!   the other comes from Eq.(40) of Casida
    osc_str(idir+1,iexcit)=&
&    (sexc(idir,iexcit)*sqrtks(iexcit2)*normint*ucvol)**2*two*two/nsppol
    osc_str(1,iexcit)=osc_str(1,iexcit)&
&    +osc_str(idir+1,iexcit)*third    
   enddo
   do ii=1,4
    if(abs(osc_str(ii,iexcit))<tol6)osc_str(ii,iexcit)=zero
   enddo
!  Changed, the whole spectrum is written
!  The array eexcit has been reordered previously, the others also                 
   if(nsppol==1)then
    write(message, '(i4,a,i3,2es12.5,es13.5,es11.4,3es9.2)' ) &
&    iocc,'->',iunocc+nband_occ(isppol),            &
&    eexcit(iexcit), eexcit(iexcit)*Ha_eV,eexcit(iexcit)+etotal, &
!    XG 020209 : Jean-Yves, I assume that the printout of sexc is for debugging ?!
!&   osc_str(1:4,iexcit),sexc(1:3,iexcit)
&    osc_str(1:4,iexcit)
   else
    write(message, '(i4,a,i3,a,i1,2es12.5,es13.5,es11.4,3es9.2)' ) &
&    iocc,'->',iunocc+nband_occ(isppol),' s:',isppol,            &
&    eexcit(iexcit), eexcit(iexcit)*Ha_eV,eexcit(iexcit)+etotal, &
&    osc_str(1:4,iexcit)
   endif
   call wrtout(ab_out,message,'COLL')
   call wrtout(6,message,'COLL')

  enddo

! Check of the Sum rule for Casida eq.47,
! only exact if complete basis of excitations, as well as local potentials only.
  sumx=zero
  do iexcit=1,nexcit_win
   sumx=sumx+osc_str(1,iexcit)
  enddo
  write(message, '(a,es16.6)' ) &
&  '  Sum of osc. strength : ',sumx
  call wrtout(ab_out,message,'COLL')
  call wrtout(6,message,'COLL')

!-Diagonalize the excitation matrices----------------------------

  allocate(eexcit2(nexcit_win),vec(2,nexcit_win,nexcit_win))

  do sing_trip=1,2/nsppol

   if(pole_approx==0)then
 
    allocate(matr(nexcit_win*(nexcit_win+1)))
    allocate(zhpev1(2,2*nexcit_win-1),zhpev2(3*nexcit_win-2))
    matr(:)=zero
    ier=0
!DEBUG
!   write(6,*)' after allocation matrices     '
!ENDDEBUG
   
   
!   Store the matrix in proper mode before calling zhpev

    index=1
    do iexcit=1,nexcit_win
     iexcit2 = indarr(iexcit)
     isppol2 = min((iexcit2-1)/nexcit_pol(1) +1,2)
     iunocc2 = (iexcit2-(isppol2-1)*nexcit_pol(1)-1)/nband_occ(isppol2)+1
     iocc2   = iexcit2-(isppol2-1)*nexcit_pol(1)-(iunocc2-1)*nband_occ(isppol2) 
     do jexcit=1,iexcit
      iexcit1 = indarr(jexcit)
      isppol1 = min((iexcit1-1)/nexcit_pol(1) +1,2)
      iunocc1 = (iexcit1-(isppol1-1)*nexcit_pol(1)-1)/nband_occ(isppol1)+1
      iocc1   = iexcit1-(isppol1-1)*nexcit_pol(1)-(iunocc1-1)*nband_occ(isppol1)

!     if (mkmem/=0) then
      matr(index)=omega_tddft_casida(iocc1,iunocc1,isppol1,iocc2,iunocc2,isppol2,sing_trip)
!     endif
      matr(index+1)=zero
      index=index+2
     enddo

    enddo

!DEBUG
!   write(6,*)' after filling    matrices     '
!ENDDEBUG
     
#if defined T3E
    call CHPEV ('V','U',nexcit_win,matr,eexcit2,vec,nexcit_win,zhpev1,&
&    zhpev2,ier)
#else
    call ZHPEV ('V','U',nexcit_win,matr,eexcit2,vec,nexcit_win,zhpev1,&
&    zhpev2,ier)
#endif

    deallocate(matr,zhpev1,zhpev2)
!DEBUG
!  write(6,*)' after deallocating matrices     '
!ENDDEBUG


   else 

    vec(:,:,:)=zero
    do isppol=1, nsppol
     do iunocc=1,nband_k(isppol) 
      do iocc=1,nband_k(isppol) 
       index=iocc+nband_k(isppol)*(iunocc-1)+(isppol-1)*nexcit_pol(1)
       eexcit2(index)=omega_tddft_casida(iocc,iunocc,isppol,iocc,iunocc,isppol,sing_trip)
       vec(1,index,index)=one
      enddo
     enddo
    enddo
 
   endif

!  Compute the excitation energies from the square root of eexcit2
!  eexcit(:)=sqrt(eexcit2(:)
  
   deallocate(eexcit)
   allocate(eexcit(nexcit_win))

   eexcit(:)=sqrt(dabs(eexcit2(:))+tol10**2)
!  Write the first excitation energies
   if(sing_trip==1)then
    if(nsppol==1)then
     write(message, '(a,a,a,a,a,a)' )ch10,&
&     '  TDDFT singlet excitation energies (at most 20 of them are printed),',&
&     ch10,'  and corresponding total energies.                ',ch10,&
&     '  Excit#   (Ha)    and    (eV)    total energy (Ha)    major contributions '
    else
     write(message, '(a,a,a,a,a,a)' )ch10,&
&     '  TDDFT   mixed excitation energies (at most 40 of them are printed),',&
&     ch10,'  and corresponding total energies.                ',ch10,&
&     '  Excit#   (Ha)    and    (eV)    total energy (Ha)    major contributions '
    endif
   else
    write(message, '(a,a,a,a,a,a)' )ch10,&
&    '  TDDFT triplet excitation energies (at most 20 of them are printed),',&
&    ch10,'  and corresponding total energies.                ',ch10,&
&    '  Excit#   (Ha)    and    (eV)    total energy (Ha)    major contributions '
   endif
   call wrtout(ab_out,message,'COLL')
   call wrtout(6,message,'COLL')
   write(6,*)' tddft : before iexcit loop'
   do iexcit=1,min(nexcit_win,nexcitout)*nsppol
    write(6,*)' tddft : iexcit=',iexcit
!   Select largest and next contributions
    flargest=zero ; fnext=zero
    nlargest=zero ; nnext=zero
    if(nsppol==2)then
     fnext1=zero  ;  fnext2=zero
     nnext1=zero  ;  nnext2=zero
    endif
    do jexcit=1,nexcit_win
     ff=vec(1,jexcit,iexcit)**2+vec(2,jexcit,iexcit)**2
     if(ff>flargest+tol12)then
      if(nsppol==2)then
       nnext2=nnext1  ; fnext2=fnext1
       nnext1=nnext   ; fnext1=fnext
      endif
      nnext=nlargest ; fnext=flargest
      nlargest=indarr(jexcit) ; flargest=ff
     else if(ff>fnext+tol12)then
      if(nsppol==2)then
       nnext2=nnext1 ; fnext2=fnext1
       nnext1=nnext  ; fnext1=fnext
      endif
      nnext=indarr(jexcit) ; fnext=ff
     else if(ff>fnext1+tol12 .AND. nsppol==2)then
       nnext2=nnext1  ; fnext2=fnext1
       nnext1=indarr(jexcit) ; fnext1=ff
     else if(ff>fnext2+tol12 .AND. nsppol==2)then
       nnext2=indarr(jexcit) ; fnext2=ff
     endif
    enddo

    isppol_l = min((nlargest-1)/nexcit_pol(1) +1,nsppol)
    iunocc_l = (nlargest-(isppol_l-1)*nexcit_pol(1)-1)/nband_occ(isppol_l)+1
    iocc_l   = nlargest-(isppol_l-1)*nexcit_pol(1)-(iunocc_l-1)*nband_occ(isppol_l)
    isppol_n = min((nnext-1)/nexcit_pol(1) +1,nsppol)
    iunocc_n = (nnext-(isppol_n-1)*nexcit_pol(1)-1)/nband_occ(isppol_n)+1
    iocc_n   = nnext-(isppol_n-1)*nexcit_pol(1)-(iunocc_n-1)*nband_occ(isppol_n)
    if(nsppol==2)then
     isppol_n1 = min((nnext1-1)/nexcit_pol(1) +1,nsppol)
     iunocc_n1 = (nnext1-(isppol_n1-1)*nexcit_pol(1)-1)/nband_occ(isppol_n1)+1
     iocc_n1   = nnext1-(isppol_n1-1)*nexcit_pol(1)-(iunocc_n1-1)*nband_occ(isppol_n1)
     isppol_n2 = min((nnext2-1)/nexcit_pol(1) +1,nsppol)
     iunocc_n2 = (nnext2-(isppol_n2-1)*nexcit_pol(1)-1)/nband_occ(isppol_n2)+1
     iocc_n2   = nnext2-(isppol_n2-1)*nexcit_pol(1)-(iunocc_n2-1)*nband_occ(isppol_n2)
    endif

    if(nsppol==1)then
     write(message,'(i4,es15.5,es14.5,es16.6,f8.2,a,i3,a,i3,a,f6.2,a,i3,a,i3,a)') &
&     iexcit,eexcit(iexcit),&
&     eexcit(iexcit)*Ha_eV,eexcit(iexcit)+etotal,&
&     flargest,'(',iocc_l,'->',iunocc_l+nband_occ(1),')',&
&     fnext,   '(',iocc_n,'->',iunocc_n+nband_occ(1),')'
     call wrtout(ab_out,message,'COLL')
     call wrtout(6,message,'COLL')
    else 
     write(chain1,'(f8.2,a,i3,a,i3,a)')flargest,'(',iocc_l,'->',iunocc_l+nband_occ(isppol_l),')'
     write(chain2,'(f8.2,a,i3,a,i3,a)')fnext,'(',iocc_n,'->',iunocc_n+nband_occ(isppol_n),')'
     if(trim(chain1)==trim(chain2))then
      write(message,'(i4,es15.5,es14.5,es16.6,a,a,a,a)') &
&      iexcit,eexcit(iexcit),&
&      eexcit(iexcit)*Ha_eV,eexcit(iexcit)+etotal,trim(chain1),'(1)',trim(chain2),'(2)'
     else
      write(message,'(i4,es15.5,es14.5,es16.6,a,a,i1,a,a,a,i1,a)') &
&      iexcit,eexcit(iexcit),&
&      eexcit(iexcit)*Ha_eV,eexcit(iexcit)+etotal,trim(chain1),'(',isppol_l,')',&
&      trim(chain2),'(',isppol_n,')'
     endif
     call wrtout(ab_out,message,'COLL')
     call wrtout(6,message,'COLL')
     write(chain1,'(f8.2,a,i3,a,i3,a)')fnext1,'(',iocc_n1,'->',iunocc_n1+nband_occ(isppol_n1),')'
     write(chain2,'(f8.2,a,i3,a,i3,a)')fnext2,'(',iocc_n2,'->',iunocc_n2+nband_occ(isppol_n2),')'
     if(trim(chain1)==trim(chain2))then
      write(message,'(a,a,a,a,a)' ) &
&      '                                                 ',&
&      chain1,'(1)',chain2,'(2)'
     else
      write(message,'(a,a,a,i1,a,a,a,i1,a)' ) &
&      '                                                 ',&
&      chain1,'(',isppol_n1,')',&
&      chain2,'(',isppol_n2,')'
     endif
     call wrtout(ab_out,message,'COLL')
     call wrtout(6,message,'COLL')
    endif
   enddo

!  For each iexcit excitation, compute the oscillator strength (Casida, eq 47)
   write(message, '(a,a,a,a)' )ch10,&
&   '  Oscillator strengths :  (elements smaller than 1.e-6 are set to zero)',ch10,&
&   '  Excit#   (Ha)   Average    XX        YY        ZZ         XY        XZ        YZ'
   call wrtout(ab_out,message,'COLL')
   call wrtout(6,message,'COLL') 
  
   do iexcit=1,nexcit_win

!   One of the factor of two comes from the spin degeneracy,
!   the other comes from Eq.(40) of Casida
    factor=(normint*ucvol)**2*two*two/nsppol
 
    osc_str(:,iexcit)=zero
   
!    factor=(normint*ucvol)**2*two*two/nsppol

    do jexcit=1,nexcit_win
     jexcit_cbase=indarr(jexcit)
     do idir=1,3
      osc_str(idir+1,iexcit)=osc_str(idir+1,iexcit)+ &
&     sexc(idir,jexcit)*sqrtks(jexcit_cbase)*sqrt(factor)* &
&     vec(1,jexcit,iexcit)
     enddo ! idir
    enddo ! jexcit

!   The "standard" definition of the oscillator strength is the square
!   of the matrix elements. 
!   So, instead of the coding 
!   do idir=1,3
!    osc_str(1,iexcit)=osc_str(1,iexcit)+osc_str(idir+1,iexcit)**2*third      
!   enddo
!   I think that the following is more "standard"
!   Now, osc_str(2:4,iexcit) are the X, Y and Z matrix elements, not
!   yet the oscillator strengths
    osc_str(5,iexcit)=osc_str(2,iexcit)*osc_str(3,iexcit)   ! off diag XY
    osc_str(6,iexcit)=osc_str(2,iexcit)*osc_str(4,iexcit)   ! off diag XZ
    osc_str(7,iexcit)=osc_str(3,iexcit)*osc_str(4,iexcit)   ! off diag ZZ
    do idir=1,3
!    Here the X,Y, and Z matrix elements are combined to give diagonal osc. strengths
     osc_str(idir+1,iexcit)=osc_str(idir+1,iexcit)**2        
     osc_str(1,iexcit)=osc_str(1,iexcit)+osc_str(idir+1,iexcit)*third ! compute the trace
    enddo
!   At this stage, osc_str(1,iexcit) is exactly the same as from your coding
!***End of section to be checked

    do ii=1,7
     if(abs(osc_str(ii,iexcit))<tol6)osc_str(ii,iexcit)=zero
    enddo
!   XG 020209 : Jean-Yves, the off-diagonal oscillator strengths
!   can become negative. It is important for automatic
!   checking that the numbers are separated by a blank, even
!   if they are negative. So replace the following format, to have at least one blank.
!   write(message, '(i4,es12.5,es10.3,3es10.3,3es10.3)' ) &
    write(message, '(i4,es12.5,es10.3,3es10.3,3es10.2)' ) &
&    iexcit,eexcit(iexcit),osc_str(1:7,iexcit)
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
   enddo

!  Check of the Sum rule for Casida eq.47, 
!  only exact if complete basis of excitations, as well as local potentials only.
   sumx=zero
   do iexcit=1,nexcit_win
    sumx=sumx+osc_str(1,iexcit)
   enddo
   write(message, '(a,es16.6)' ) &
&   '  Sum of osc. strength : ',sumx
   call wrtout(ab_out,message,'COLL')
   call wrtout(6,message,'COLL')

!  If singlet, compute Cauchy coefficients
   if(sing_trip==1.AND.nsppol==1)then 
    cauchy(:)=zero   
    do iexcit=1,nexcit_win
     fr_invsquare=one/(eexcit(iexcit)**2)
     fr_power=one
     do ii=1,7
      fr_power=fr_power*fr_invsquare
      cauchy(ii)=cauchy(ii)+osc_str(1,iexcit)*fr_power
     enddo
    enddo
    write(message, '(a,es11.3,a,es11.3,a,es11.3,a,a,es11.3,a,es11.3,a,es11.3,a,es11.3)' ) &
&    '  Cauchy coeffs (au) : ( -2)->',cauchy(1),&
&      ', ( -4)->',cauchy(2),', ( -6)->',cauchy(3),ch10,&
&    '    (-8)->',cauchy(4),', (-10)->',cauchy(5),', (-12)->',cauchy(6),', (-14)->',cauchy(7)
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
   endif

! End the loop on singlet or triplet
  enddo

  deallocate(eexcit2,vec)
  deallocate(osc_str)

  close(tmp_unit2,status='delete')

  call timab(95,2,tsec)

 endif  ! end of am_master

 deallocate(omega_tddft_casida,eexcit,sqrtks,sexc,done_sexc)
 deallocate(indarr)
 deallocate(done_sexc2)

 return

end subroutine tddft
!!***
