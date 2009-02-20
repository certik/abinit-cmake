!{\src2tex{textfont=tt}}
!!****f* ABINIT/outwf
!! NAME
!! outwf
!!
!! FUNCTION
!! Conduct output of a "wave-functions" file.
!!  - Compute the maximal residual
!!  - Then open a permanent file wff2 for final output of wf data
!!  - Create a new header for the file.
!!  - Write wave-functions (and energies)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, AR, MB, MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=wavefunction array (storage if nkpt>1)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen( (2*mband)**response *mband*nkpt*nsppol)=
!!                  eigenvalues (hartree) for all bands at each k point
!!  filnam= character string giving the root to form the name of the
!!   output WFK or WFQ file if response==0, otherwise it is the filename.
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kptns(3,nkpt)=k points in terms of recip primitive translations
!!  mband=maximum number of bands
!!  mkmem=maximum number of k-points in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum number of plane waves
!!  mxfh=last dimension of the xfhist array
!!  natom=number of atoms in unit cell
!!  nband=number of bands
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  npwarr(nkpt)=number of planewaves in basis and on boundary for each k
!!  nqpt: if 0 if no q-point shift has to be taken into account, and the
!!   output wavefunction file will be appended with _WFK ;
!!    if 1 a q-point shift has been taken into account, and one should use _WFQ,
!!    except when response=1
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nstep=desired number of electron iteration steps
!!  nxfh=actual number of (x,f) history pairs, see xfhist array.
!!  occ(mband*nkpt*nsppol)=occupations for all bands at each k point
!!  resid(mband*nkpt*nsppol)=squared residuals for each band and k point
!!   where resid(n,k)=|<C(n,k)|(H-e(n,k))|C(n,k)>|^2 for the ground state
!!  response: if == 0, GS wavefunctions , if == 1, RF wavefunctions
!!  wffnow=structure information for current wavefunction (if nkpt>1)
!!  xfhist(3,natom+4,2,mxfh)=(x,f) history array,
!!                                 also includes rprim and stress
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! * The name of the file wff2 might be the same as that of the file wff1.
!! * The routine includes closing wffnow.
!!
!! PARENTS
!!      berryphase_new,gstate,loper3
!!
!! CHILDREN
!!      handle_ncerr,hdr_io,hdr_io_etsf,hdr_io_netcdf,hdr_skip,ini_wf_netcdf
!!      leave_new,leave_test,mpi_barrier,rwwf,timab,wffclose,wffdelete,wffkg
!!      wffoffset,wffopen,wrtout,wvl_write,xcomm_init,xdefineoff,xexch_mpi
!!      xmaster_init,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine outwf(cg,dtfil,dtset,eigen,filnam,hdr,kg,kptns,mband,mkmem,&
 &                mpi_enreg,mpw,mxfh,natom,nband,nfft,ngfft,nkpt,npwarr,&
 &                nqpt,nspinor,nsppol,nstep,nxfh,occ,resid,response,&
 &                wffnow,wfs,xfhist)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes
#if defined HAVE_NETCDF
 use netcdf
#endif
#if defined HAVE_BIGDFT
 use BigDFT_API
#endif

#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_13io_mpi
 use interfaces_13ionetcdf
 use interfaces_14wvl_wfs
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif
!Arguments ------------------------------------
 integer, intent(in) :: mband,mkmem,mpw,mxfh,natom,nfft,nkpt,nqpt,nspinor,nsppol
 integer, intent(in) :: nstep,nxfh,response
 character(len=fnlen), intent(in) :: filnam
 type(MPI_type), intent(inout) :: mpi_enreg
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(in) :: dtset
 type(hdr_type), intent(inout) :: hdr
 type(wffile_type), intent(inout) :: wffnow
 type(wvl_wf_type), intent(in) :: wfs
 integer, intent(in) :: kg(3,mpw*mkmem),nband(nkpt*nsppol),ngfft(18),npwarr(nkpt)
 real(dp), intent(inout) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp), intent(in) :: eigen((2*mband)**response*mband*nkpt*nsppol),kptns(3,nkpt)
 real(dp), intent(in) :: occ(mband*nkpt*nsppol),resid(mband*nkpt*nsppol)
 real(dp), intent(in) :: xfhist(3,natom+4,2,mxfh)

!Local variables-------------------------------
 integer,parameter :: nkpt_max=50
 integer :: accesswff,action,band_index,fform,formeig,headform,iband,ibdkpt,icg
  integer :: icg0,ierr,ii,ikg,ikpt,ios,ipw,isppol,ixfh,jj,master,mcg,mcg_disk,me,me0
 integer :: nband_disk,nband_k,nkpt_eff,nmaster,npw_k,option,optkg,rdwr,sender,source
 integer :: spaceComm,sread,sskip,tim_rwwf
 integer :: xfdim2
 real(dp) :: residk,residm,resims
 logical :: ex,mydata,od,swrite,tmaster,tread
 character(len=fnlen) :: wffnm
 character(len=3) :: ipara
 character(len=500) :: message
 type(wffile_type) :: wff2
 integer,allocatable :: kg_disk(:,:)
 integer :: spaceCommMPIO
 real(dp) :: tsec(2)
!There are two cases of use of the cg_disk : if mkmem==0, or parallel treatment.
 real(dp),allocatable :: cg_disk(:,:),eig_dum(:),eig_k(:),occ_dum(:),occ_k(:)
!no_abirules
#if defined MPI
           !Variables introduced for MPI version
           integer :: ipwnbd
#endif

#if defined HAVE_NETCDF
           !netCDF variables
           integer :: ncid_hdr,ncerr
           integer :: dimg3_id,mband_id,nkpt_id,nspinor_id,nsppol_id
           integer :: complex_id,eigendim_id,mpw_id,mkmem_id
           integer :: kg_id,eigen_id,cg_id
           integer :: nxfh_id, mxfh_id, xfdim2_id, dim2inout_id, dimr3_id,xfhist_id
           integer :: nxfh_tmp,mxfh_tmp,xfdim2_tmp,dim2inout_tmp

#endif

! *************************************************************************
!For readability of the source file, define a "me" variable
!also in the sequential case

!DEBUG
!  write(6,*)' outwf : enter'
!  write(6,*)' outwf : trim(filnam)=',trim(filnam)
!ENDDEBUG
 xfdim2 = natom+4

!Init me
 call xme_init(mpi_enreg,me)
  me0=me
  if (mpi_enreg%mode_para=='b') then
     me = mpi_enreg%me_kpt
  end if
!Define master
 call xmaster_init(mpi_enreg,master)
!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)
 tim_rwwf =0
 source = master
 sread = master
 tmaster=(master==me)
 swrite=tmaster
 sender=-1

!Compute mean square and maximum residual over all bands and k points
!and spins
!(disregard k point weights and occupation numbers here)
 band_index=sum(nband(1:nkpt*nsppol))
 resims=sum(resid(1:band_index))/dble(band_index)

!Find largest residual over bands, k points, and spins,
!except for nbdbuf highest bands
!Already AVAILABLE in hdr ?!
 ibdkpt=1
 residm=0.0_dp
 do isppol=1,nsppol
  do ikpt=1,nkpt
   nband_k=nband(ikpt+(isppol-1)*nkpt)
   nband_k=max(1,nband_k-dtset%nbdbuf)
   residm=max(residm,maxval(resid(ibdkpt:ibdkpt+nband_k-1)))
   ibdkpt=ibdkpt+nband_k
  end do
 end do

 write(message, '(a,1p,e12.4,a,e12.4)' ) &
&   ' Mean square residual over all n,k,spin= ',resims,'; max=',residm
 call wrtout(ab_out,message,'COLL')

 band_index=0
 nkpt_eff=nkpt
 if( (dtset%prtvol==0 .or. dtset%prtvol==1) .and. nkpt_eff>nkpt_max ) nkpt_eff=nkpt_max

!Loop over spin again
 do isppol=1,nsppol
! Give (squared) residuals for all bands at each k
  do ikpt=1,nkpt
   nband_k=nband(ikpt+(isppol-1)*nkpt)
!  Will not print all residuals when prtvol=0 or 1
   if(ikpt<=nkpt_eff)then
!   Find largest residual over all bands for given k point
    residk=maxval(resid(1+band_index:nband_k+band_index))
    write(message, '(1x,3f8.4,3x,i2,1p,e13.5,a)' ) &
&     kptns(1:3,ikpt),isppol,residk,' kpt; spin; max resid(k); each band:'
    call wrtout(ab_out,message,'COLL')
    do ii=0,(nband_k-1)/8
     write(message, '(1p,8e9.2)' ) &
&     (resid(iband+band_index),iband=1+ii*8,min(nband_k,8+ii*8))
     call wrtout(ab_out,message,'COLL')
    end do
   else if(ikpt==nkpt_eff+1)then
    write(message, '(a,a)' ) &
&    ' outwf : prtvol=0 or 1, do not print more k-points.',ch10
    call wrtout(ab_out,message,'COLL')
   end if
   band_index=band_index+nband_k
  end do
 end do

!Will write the wavefunction file only when nstep>0
 if (nstep>0 .and. dtset%prtwf/=0) then

  if(response==0 .and. nqpt==0)wffnm=trim(filnam)//'_WFK'
  if(response==0 .and. nqpt==1)wffnm=trim(filnam)//'_WFQ'
! If RF case, the appendix _1WF has already been added.
  if(response==1)wffnm=trim(filnam)
  if (mpi_enreg%parareel == 1)then
        if (mpi_enreg%paral_compil_mpio == 1 .and. dtset%accesswff == 1 ) then
         wffnm=trim(wffnm)
        else
          if (mpi_enreg%ipara < 10) write(ipara,'(i1)')mpi_enreg%ipara
        if (mpi_enreg%ipara >= 10) write(ipara,'(i2)')mpi_enreg%ipara
        if (mpi_enreg%ipara >= 100) write(ipara,'(i3)')mpi_enreg%ipara
        wffnm=trim(wffnm)//'_'//ipara
       end if
  end if
! Only the master write the file, except if MPI I/O, but the
! full wff dataset should be provided to WffOpen in this case
  accesswff=-1
  if(dtset%accesswff==1) then
    accesswff=1
#if defined HAVE_NETCDF
  else if(dtset%accesswff==2) then
    accesswff=2
!  Create empty netCDF file
    ncerr = nf90_create(path=wffnm, cmode=NF90_CLOBBER, ncid=ncid_hdr)
    call handle_ncerr(ncerr," create netcdf wavefunction file")
    ncerr = nf90_close(ncid_hdr)
    call handle_ncerr(ncerr," close netcdf wavefunction file")
#endif
#if defined HAVE_ETSF_IO
  else if (dtset%accesswff == 3) then
    accesswff = 3
#endif
  end if
!DEBUG
! write (6,*) 'outwf : accesswff = ', accesswff
! write (6,*) 'outwf : wffnow%accesswff = ', wffnow%accesswff,wffnow%kgwff
!ENDDEBUG
  if (accesswff==1.and.dtset%paral_kgb==1) then !Band/fft parallelisation
   spaceCommMPIO=mpi_enreg%commcart
  else
   spaceCommMPIO=spaceComm
  endif
  call WffOpen(accesswff,spaceCommMPIO,wffnm,ierr,wff2,master,me0,dtfil%unwff2)

! Conduct wavefunction output to wff2
  write(message, '(a,a)' )&
&  ' outwf  : write wavefunction to file ',trim(wffnm)
  if (mpi_enreg%parareel == 0) then
          call wrtout(6,message,'COLL')
        else 
        if (master==me) then
                  call wrtout(6,message,'PERS')
        end if
  end if

  allocate(kg_disk(3,mpw))

  mcg=mpw*nspinor*mband*mkmem*nsppol
  mcg_disk=mpw*nspinor*mband
  formeig=0 ; if(response==1)formeig=1
  if(mkmem == 0 )then
   allocate(eig_dum( (2*mband)**formeig * mband),occ_dum(mband))
  end if
  allocate(eig_k( (2*mband)**formeig * mband))
  allocate(occ_k(mband))

#if defined MPI
!BEGIN TF_CHANGES
  call leave_test(mpi_enreg)
!END TF_CHANGES
  ! Compute mband and mpw
  if(mkmem/=0)allocate(cg_disk(2,mcg_disk))
#endif

  if (mkmem==0) then

!  Skip wffnow header
   call hdr_skip(wffnow,ierr)

   allocate(cg_disk(2,mcg_disk))
!  Define offsets, in case of MPI I/O
   call WffKg(wffnow,1)
   call xdefineOff(formeig,wffnow,mpi_enreg,nband,npwarr,nspinor,nsppol,nkpt)
  end if  !mkmem = 0


  band_index=0
  icg=0
  if((mpi_enreg%paralbd==0) .or. (mpi_enreg%paralbd>1))tim_rwwf=6
  if(mpi_enreg%paralbd==1)tim_rwwf=12

! Write header info for new wf file
  rdwr=2
  if (dtset%usewvl == 0) then
    fform=2
  else
    ! Use 200 as radical for naming file format
    ! used by wavelets.
    fform = 200
  end if

   if (wff2%accesswff < 2) then
    call hdr_io(fform,hdr,rdwr,wff2)
    call WffKg(wff2,1)
#if defined HAVE_NETCDF
   else if (wff2%accesswff == 2.and.tmaster) then

!DEBUG
!   write (6,*) 'outwf : entering hdr_io_netcdf'
!ENDDEBUG

    call hdr_io_netcdf(fform,hdr,rdwr,wff2)

    call ini_wf_netcdf(mpw,wff2%unwff,response)
#endif
#if defined HAVE_ETSF_IO
   else if (wff2%accesswff == 3 .and. tmaster) then

!DEBUG
!   write (6,*) 'outwf : entering hdr_io_etsf'
!ENDDEBUG

    call hdr_io_etsf(fform, hdr, rdwr, wff2%unwff)
#endif
   end if

  do isppol=1,nsppol

   ikg=0

   do ikpt=1,nkpt
    nband_k=nband(ikpt+(isppol-1)*nkpt)
    npw_k=npwarr(ikpt)

!   Read the wavefunction block, without the eigenvalues
    if(mkmem==0)then
#if defined MPI
           sread=-1
           if (mpi_enreg%parareel == 0) then
            if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))==0) sread=me
           else
            if(mpi_enreg%proc_distrb_para(mpi_enreg%ipara,ikpt) ==me ) sread=me
           end if
#endif
     if(sread==me)then

      headform=0 ; icg0=0 ; option=-2 ; optkg=1
      call rwwf(cg_disk,eig_dum,formeig,headform,&
&      icg0,ikpt,isppol,kg_disk,mband,mcg_disk,mpi_enreg,nband_k,nband_disk,&
&      npw_k,nspinor,occ_dum,option,optkg,tim_rwwf,wffnow)

      if(nband_k/=nband_disk)then
       write(message, '(a,a,a,a,i4,a,i6,a,a,a,i6,a)' ) ch10,&
&       ' outwf : BUG -',ch10,&
&       '  For k pt number',ikpt,' disk file has',nband_disk,' bands',ch10,&
&       '  but input file gave nband=',nband_k,'.'
       call wrtout(06,message,'PERS')
       call leave_new('PERS')
      end if ! nband check

     end if ! sread==me
    end if ! mkmem

#if defined MPI
    if (dtset%usewvl == 0) then
           if (mpi_enreg%parareel == 0) then
!BEGIN TF_CHANGES
            call MPI_BARRIER(mpi_enreg%spaceComm,ierr)
!END TF_CHANGES
           else
            call MPI_BARRIER(mpi_enreg%kpt_comm_para(mpi_enreg%ipara),ierr)
           end if

          !Must transfer the wavefunctions to the master processor
          !Separate sections for paralbd=1 or other values ; might be merged
           if(mpi_enreg%paralbd==0 .or. mpi_enreg%paralbd>1)then
            if(mpi_enreg%parareel == 0) then
             nmaster=0
             source=minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol))
            else
             nmaster=mpi_enreg%master_group_para
             source=mpi_enreg%proc_distrb_para(mpi_enreg%ipara,ikpt)
            end if

            mydata=.false.
            if(source==me)mydata=.true.
            action=0
          ! I am the master node, and I have the data in cg or cg_disk
            if((tmaster).and.(mydata))action=1
          ! I am not the master, and I have the data => send to master
            if((.not.tmaster).and.(mydata))action=2
          ! I am the master, and I receive the data
            if((tmaster).and.(.not.mydata))action=3

          ! I have the data in cg or cg_disk ( MPI_IO case)
            if (accesswff==1  ) then
               action = 0
               sender=-1
               swrite=.false.
               if (mydata)then
                action=1        
                swrite=.true.
                sender=me
               endif
            endif



          ! I am the master node, and I have the data in cg or cg_disk
          ! I have the data in cg or cg_disk ( MPI_IO case)
            if(action==1)then
             if(mkmem/=0)then
          !   Copy from kg to kg_disk
              kg_disk(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
          !   Copy from cg to cg_disk
              do ipwnbd=1,nband_k*npw_k*nspinor
               cg_disk(1,ipwnbd)=cg(1,ipwnbd+icg)
               cg_disk(2,ipwnbd)=cg(2,ipwnbd+icg)
              end do
             end if
            end if


          ! I am not the master, and I have the data => send to master
          ! I am the master, and I receive the data
            if ( action.eq.2.or.action.eq.3) then
             call timab(48,1,tsec)
             if(mkmem/=0 .and. action==2)then
          !    I am not the master, and I have the data => send to master, not in
          !    place (meaning this part of the array might not be defined in the
          !    master node, so caution is in order that the master does no use
          !    these calls)
              call xexch_mpi(kg(:,1+ikg:npw_k+ikg),3*npw_k,source,kg_disk,nmaster &
               &             ,mpi_enreg%spaceComm,ierr)
              call xexch_mpi(cg(:,icg+1:icg+nband_k*npw_k*nspinor),2*nband_k*npw_k*nspinor &
               &            ,source,cg_disk,nmaster &
               &            ,mpi_enreg%spaceComm,ierr)
             else
              call xexch_mpi(kg_disk,3*npw_k,source,kg_disk,nmaster &
               &             ,mpi_enreg%spaceComm,ierr)
              call xexch_mpi(cg_disk,2*nband_k*npw_k*nspinor,source,cg_disk,nmaster &
              &             ,mpi_enreg%spaceComm,ierr)
             end if
             call timab(48,2,tsec)
            end if


           else if(mpi_enreg%paralbd==1)then
            nmaster=0
#if defined MPI_IO
             sender=-1
            if( accesswff ==1 ) then
              nmaster=mpi_enreg%proc_distrb(ikpt,1,isppol)
              sender=nmaster
            end if
#endif
  
          ! Note the loop over bands
            do iband=1,nband_k

          !  The message passing related to kg is counted as one band
             action=0

          !   I am the master node, and I have the data in cg or cg_disk
             if( mpi_enreg%proc_distrb(ikpt,iband,isppol)==nmaster .and. &
&                 me==nmaster) then
                action=1
          !   I am not the master, and I have the data => send to master
             elseif( mpi_enreg%proc_distrb(ikpt,iband,isppol)==me &
&                 .and. me/=nmaster ) then
                action = 2
        !   I am the master, and I receive the data
             elseif( mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me &
&                 .and. me==nmaster ) then
                action=3
             end if

             if(action==1) then
          !   I am the master node, and I have the data in cg or cg_disk
              if(mkmem/=0)then
          !    Copy from kg to kg_disk
               if(iband==1)kg_disk(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
          !    Copy from cg to cg_disk
               do ipwnbd=1,npw_k*nspinor
                cg_disk(1,(iband-1)*npw_k*nspinor+ipwnbd)= &
&                       cg(1,(iband-1)*npw_k*nspinor+ipwnbd+icg)
                cg_disk(2,(iband-1)*npw_k*nspinor+ipwnbd)= &
&                       cg(2,(iband-1)*npw_k*nspinor+ipwnbd+icg)
               end do
              end if
             end if  ! action=1

             if ( action.eq.2.or.action.eq.3) then
          ! action=2 :  I am not the master, and I have the data => send to master
          ! action=3 :  I am the master, and I receive the data
              call timab(48,1,tsec)
              if ( iband == 1 ) then
               if ( mkmem/=0 .and. action==2 ) then
                call xexch_mpi(kg(:,1+ikg:npw_k+ikg),3*npw_k,mpi_enreg%proc_distrb(ikpt,iband,isppol) &
                &             ,kg_disk,nmaster &
                &             ,mpi_enreg%spaceComm,ierr)
               else
                call xexch_mpi(kg_disk,3*npw_k,mpi_enreg%proc_distrb(ikpt,iband,isppol)  &
                &             ,kg_disk,nmaster &
                &             ,mpi_enreg%spaceComm,ierr)
               end if
              end if       ! iband =1
              ipwnbd=(iband-1)*npw_k*nspinor
              if(mkmem/=0 .and. action==2 )then
               call xexch_mpi( cg(:,ipwnbd+icg+1:ipwnbd+icg+npw_k*nspinor),2*npw_k*nspinor & 
               &              ,mpi_enreg%proc_distrb(ikpt,iband,isppol)                    &
               &              ,cg_disk(:,ipwnbd+1:ipwnbd+npw_k*nspinor),nmaster            &
               &              ,mpi_enreg%spaceComm,ierr)
              else
               call xexch_mpi( cg_disk(:,ipwnbd+1:ipwnbd+npw_k*nspinor),2*npw_k*nspinor    &
               &              ,mpi_enreg%proc_distrb(ikpt,iband,isppol)                    &
               &              ,cg_disk(:,ipwnbd+1:ipwnbd+npw_k*nspinor),nmaster            &
               &              ,mpi_enreg%spaceComm,ierr)
              end if
 
              call timab(48,2,tsec)
             end if        ! action=2 or action=3


             
             if(accesswff ==1 ) then
              !   I have the data in cg or cg_disk
              swrite=.false.
              if (nmaster == me) then
               swrite=.true. 
              end if
             end if

          ! End of loop over bands
            end do

          !End of paralbd=1
           end if
        end if
#endif

!   Only the master will write to disk the final output wf file.
!  in MPI_IO case only swrite will write to disk the final output wf file.
    if(swrite) then
!DEBUG
!    write(6,*) 'outwf : I am master and will write wf file'
!ENDDEBUG
     if(formeig==0)then
      eig_k(1:nband_k)=eigen(1+band_index:nband_k+band_index)
      occ_k(1:nband_k)=occ(1+band_index:nband_k+band_index)
     else
      eig_k(1:2*nband_k*nband_k)=eigen(1+band_index:2*nband_k*nband_k+band_index)
     end if
     option=2
     !if (dtset%prtwf == 2 .and. mkmem/=0) option=4

     if (dtset%usewvl == 0) then
#if defined MPI 
        call rwwf(cg_disk,eig_k,formeig,0,0,ikpt,isppol,kg_disk,mband,mcg_disk,mpi_enreg, &
             & nband_k, nband_k,npw_k,nspinor,occ_k,option,1,tim_rwwf,wff2)
#elif !defined MPI
        if(mkmem==0)then
           call rwwf(cg_disk,eig_k,formeig,0,0,ikpt,isppol,kg_disk,mband,mcg_disk,mpi_enreg, &
                & nband_k,nband_k, npw_k,nspinor,occ_k,2,1,tim_rwwf,wff2)
        else if(mkmem/=0)then
           kg_disk(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
           call rwwf(cg,eig_k,formeig,0,icg,ikpt,isppol,kg_disk,mband,mcg,mpi_enreg,nband_k, &
                & nband_k, npw_k,nspinor,occ_k,option,1,tim_rwwf,wff2)
        end if
#endif
     else
        call wvl_write(dtset, eigen, mpi_enreg, option, hdr%rprimd, &
             & wff2, wfs, hdr%xred)
     end if
    end if

!   The wavefunctions for the present k point and spin are written
    if(response==0)band_index=band_index+nband_k
    if(response==1)band_index=band_index+2*nband_k*nband_k

    if (mkmem/=0) then

     sskip=1
#if defined MPI
     if (dtset%usewvl == 0) then
        sskip=0
        if (mpi_enreg%parareel == 0) then
           if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))==0)sskip=1
        else
           if(mpi_enreg%proc_distrb_para(mpi_enreg%ipara,ikpt)==me)sskip=1
        end if
     end if
#endif
     if(sskip==1)then
      icg=icg+npw_k*nspinor*nband_k
      ikg=ikg+npw_k
     end if

    end if !mkem/=0


#if defined MPI_IO
         call WffOffset(wff2,sender,spaceComm,ierr)
#endif

   end do ! ikpt
  end do ! isppol
  deallocate(kg_disk)
  if(mkmem==0)deallocate(cg_disk)
#if defined MPI
            if(mkmem/=0)deallocate(cg_disk)
#endif

  if(mkmem==0) deallocate(eig_dum,occ_dum)
  deallocate(eig_k,occ_k)

! Write the (x,f) history

     if(me0==0 .and. nxfh>0 .and. response==0)then
   if (wff2%accesswff /= 2) then
#if defined MPI_IO
      if(wff2%accesswff == 1 ) then
        close(unit=wff2%unwff)
!       the file is to be positioned at the terminal point
        open(unit=wff2%unwff,POSITION="APPEND")
        endif
#endif
    write(unit=wff2%unwff)nxfh
    do ixfh=1,nxfh
     write(unit=wff2%unwff)xfhist(:,:,:,ixfh)
    end do
#if defined HAVE_NETCDF
   else if (wff2%accesswff == 2) then
     ncid_hdr=wff2%unwff

! check if nxfh and xfhist are defined
     ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)

     if (ncerr /= NF90_NOERR) then
! need to define everything
       ncerr = nf90_redef (ncid=ncid_hdr)
       call handle_ncerr(ncerr," outwf : going to define mode ")

       ncerr = nf90_def_dim(ncid=ncid_hdr,name="dim2inout",len=2,dimid=dim2inout_id)
       call handle_ncerr(ncerr," outwf : define dim2inout")
       ncerr = nf90_def_dim(ncid=ncid_hdr,name="mxfh",len=mxfh,dimid=mxfh_id)
       call handle_ncerr(ncerr," outwf : define mxfh")
       ncerr = nf90_def_dim(ncid=ncid_hdr,name="nxfh",len=nxfh,dimid=nxfh_id)
       call handle_ncerr(ncerr," outwf : define nxfh")
       ncerr = nf90_def_dim(ncid=ncid_hdr,name="xfdim2",len=xfdim2,dimid=xfdim2_id)
       call handle_ncerr(ncerr," outwf : define xfdim2")

       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="dimr3",dimid=dimr3_id)
       call handle_ncerr(ncerr," outwf : inquire dimr3")

! xfhist(3,natom+4,2,mxfh)
       ncerr = nf90_def_var(ncid=ncid_hdr,name="xfhist",xtype=NF90_DOUBLE,&
          & dimids=(/dimr3_id,xfdim2_id,dim2inout_id,mxfh_id/),varid=xfhist_id)
       call handle_ncerr(ncerr," outwf : define xfhist")

! End define mode and go to data mode
       ncerr = nf90_enddef(ncid=ncid_hdr)
       call handle_ncerr(ncerr," outwf : enddef call ")
     else
! check that the dimensions are correct
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)
       call handle_ncerr(ncerr," outwf : inquire nxfh")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nxfh_id,&
            &   len=nxfh_tmp)
       call handle_ncerr(ncerr,"  outwf : get nxfh")
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="xfdim2",dimid=xfdim2_id)
       call handle_ncerr(ncerr," outwf : inquire xfdim2")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=xfdim2_id,&
            &   len=xfdim2_tmp)
       call handle_ncerr(ncerr,"  outwf : get xfdim2")
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="mxfh",dimid=mxfh_id)
       call handle_ncerr(ncerr," outwf : inquire mxfh")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=mxfh_id,&
            &   len=mxfh_tmp)
       call handle_ncerr(ncerr,"  outwf : get mxfh")
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="dim2inout",dimid=dim2inout_id)
       call handle_ncerr(ncerr," outwf : inquire dim2inout")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=dim2inout_id,&
            &   len=dim2inout_tmp)
       call handle_ncerr(ncerr,"  outwf : get dim2inout")

       ncerr = nf90_inq_varid(ncid=ncid_hdr,name="xfhist",varid=xfhist_id)
       call handle_ncerr(ncerr," outwf : inquire xfhist")

       if (mxfh_tmp /= mxfh .or. dim2inout_tmp /= 2 .or. xfdim2_tmp /= xfdim2) then
        stop
       end if

     end if

! Now fill the data
     ncerr = nf90_put_var(ncid=ncid_hdr,varid=xfhist_id,values=xfhist,&
        & start=(/1,1,1,1/),count=(/3,xfdim2,2,nxfh/))
     call handle_ncerr(ncerr," outwf : fill xfhist")

! end NETCDF definition ifdef
#endif
   end if
! end accesswff if

  end if
! end tmaster if

! Close the wavefunction file (and do NOT delete it !)
  if (wff2%accesswff /= 2) then
    call WffClose(wff2,ierr)
#          if defined HAVE_NETCDF
  else if (wff2%accesswff == 2 .and. tmaster) then
    ncerr = nf90_close(wff2%unwff)
    call handle_ncerr(ncerr," close netcdf wavefunction file")
#          endif
  end if
!End condition of nstep>0
 end if

!Close the temporary data file, if any
 if (mkmem==0) then
  call WffDelete(wffnow,ierr)
 end if

!DEBUG
!write(6,*)' outwf : exit'
!ENDDEBUG

end subroutine outwf
!!***
