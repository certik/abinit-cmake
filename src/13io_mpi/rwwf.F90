!{\src2tex{textfont=tt}}
!!****f* ABINIT/rwwf
!! NAME
!! rwwf
!!
!! FUNCTION
!!  This subroutine reads (different options) or write (option=2) the block of records
!!  related to one k point, and one spin-polarization, that
!!  contains the wavefunctions (as well as the eigenvalues and occupations).
!!  If called with option -1, the records will be skipped.
!!  If called with option -2, only the wavefunctions are read.
!!  The disk file unitwf should have been prepared
!!  outside of this routine, in order to read or write the correct records.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, MVer,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  formeig=format of the eigenvalues
!!   0 => vector of eigenvalues
!!   1 => hermitian matrix of eigenvalues
!!  headform=format of the header of the wf file, also governing the k block format
!!   in case headform=0, use the default (current) format and headform
!!  icg=shift to be given to the location of the cg array
!!  ikpt=index of current k point (only needed for error message)
!!  isppol=spin polarization currently treated (only needed for error message)
!!  mband=maximum number of bands (dimension of cg, eigen and occ)
!!  mcg=dimention of cg
!!  nband=number of bands actually in cg, eigen and occ
!!   (if writing mode : must be larger or equal to nband_disk, only nband_disk
!!     bands are written ;
!!    if reading mode : can be equal, larger or smaller than nband_disk, but
!!     cg, eigen and occ will not be completely filled if nband>nband_disk)
!!  nband_disk=number of bands on the disk file
!!  npw=number of plane waves
!!  nspinor=number of spinorial components of the wavefunctions
!!  option= 2 for write,
!!          1 for reading cg, eigen and occ,
!!         -1 for reading/skipping,
!!         -2 for reading cg only
!!          3 for reading the eigenvalues only
!!          4 for writing a file containing only eigenvalues and occupations
!!         -4 for reading a file written with 4
!!                 (different from 3 which reads a normal option 2 file)
!!  optkg= if 1 , read or write kg_k ; if 0, do not care about kg_k
!!  tim_rwwf=timing code of the calling routine (set to 0 if not attributed)
!!  wff=struct info for wavefunction
!!   | unitwf=unit number for wavefunction
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  cg(2,npw*nspinor*mband)=planewave coefficients of wavefunctions,
!!   input if option=2; output if option=1 or -2
!!  eigen((2*mband)**formeig *mband)=array for holding eigenvalues (hartree)
!!   input if option=2 or 4; output if option=1
!!  kg_k(3,optkg*npw)=k+g data  (only if optkg==1)
!!   input if option=2; output if option=1 or -2
!!  nband_disk=number of bands on disk
!!   input if option=2 or 4; output in the other cases
!!  occ(mband)=array for holding eigenvalues (hartree)
!!   input if option=2 or 4; output if option=1
!!   no meaning if frmeig/=0
!!
!! NOTES
!!  WARNING : occ is not read in the present status of this routine
!!  WARNING : skipping k-bloks is also done in the randac subroutine
!!  WARNING : reading the two first records is also done in the rdnpw routine
!!  WARNING : writing the two first records is also done in the vtowfk3 routine
!!
!! TODO
!!  Some arguments are contained in the wff datastructure, and should be eliminated.
!!  option 3 should be called -3 (reading -> negative option) and others (-1,1) re-shuffled.
!!
!! PARENTS
!!      WffReadEigK,WffReadSkipK,berryphase,compare_interpol,dyfnl3,eltfrkin3
!!      eltfrnl3,energy,forstrnps,getgsc,initwf,ladielmt,lavnl,mkrho,mkrho3
!!      mrggkk,newkpt,outkss,outwant,outwf,overlap_wf,partial_dos_fractions
!!      prctfvw1,prctfvw2,rhofermi3,suscep_dyn,suscep_kxc_dyn,suscep_stat,tddft
!!      uderiv,vtorho,vtorho3,wffile,wfread,wfsinp
!!
!! CHILDREN
!!      etsf_io_electrons_put,etsf_io_low_error_to_str,etsf_io_main_put
!!      etsf_io_basisdata_put,handle_ncerr,leave_new,wffwritenpwrec,wrtout
!!      xderivewrecend,xderivewrecinit,xderivewrite
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine rwwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,nband,&
     & nband_disk,npw,nspinor,occ,option,optkg,tim_rwwf,wff)

 use defs_basis
 use defs_datatypes
#if defined HAVE_NETCDF
 use netcdf
#endif
#if defined HAVE_ETSF_IO
 use etsf_io
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_13io_mpi, except_this_one => rwwf
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: formeig,headform,icg,ikpt,isppol,mband,mcg,nband,npw
 integer,intent(inout) :: nband_disk
 integer,intent(in) :: nspinor,option,optkg,tim_rwwf
 integer,intent(inout), target :: kg_k(3,optkg*npw)
 real(dp),intent(inout), target :: cg(2,mcg),eigen((2*mband)**formeig*mband),occ(mband)
 type(wffile_type),intent(inout) :: wff
  type(MPI_type), intent(inout) :: mpi_enreg

!Local variables-------------------------------
  integer :: iband,ii,indxx,ios,ipw,nband1,npw1,npwso,npwso1,nspinor1,unitwf,nband2,npwtot
 integer :: use_f90,nrec
 character(len=500) :: message
 real(dp) :: tsec(2)
 integer :: ncid_hdr

!no_abirules
#if defined MPI
           integer :: ndim,nvar,natt,uid
           integer :: generic_id
           character(len=500) :: generic_name
           character(len=10) :: tmpstr
#endif
#if defined HAVE_NETCDF
 integer :: ncerr, kg_id, eigen_id, occ_id, cg_id
#endif
#if defined HAVE_ETSF_IO
 character(len=fnlen) :: file_etsf
 character(len=etsf_io_low_error_len) :: errmess
 type(etsf_main) :: main_folder
 type(etsf_basisdata),target :: wave_folder
 type(etsf_electrons),target :: electrons_folder
 logical :: lstat
 type(etsf_io_low_error) :: error
#endif

! *************************************************************************

!DEBUG
!write(6,*)' rwwf : enter, option=',option
!write(6,*)' rwwf : nband,nband_disk=',nband,nband_disk
!write(6,*)' rwwf : wff%accesswff= ', wff%accesswff
!write(6,*)' rwwf : wff%offwff= ', wff%offwff,trim(wff%fname),'par ',wff%me
!if(headform/=40 .and. headform/=0)then
! write(6,*)' rwwf : headform is',headform
! stop
!end if
!ios=0
!stop
!ENDDEBUG

 call timab(270+tim_rwwf,1,tsec)

!Might check that icg+npw*nband*nspinor is smaller than mcg

!Check that nband is smaller than mband, if one will not skip the records.
 if(nband>mband .and. option/=-1)then
  write(message,'(a,a,a,a,a,a,i5,a,i5,a)')ch10,&
&  ' rwwf : BUG -',ch10,&
&  '  One should have nband<=mband',ch10,&
&  '  However, nband=',nband,', and mband=',mband,'.'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if

!Check that formeig is 0 or 1.
 if(formeig/=0 .and. formeig/=1)then
  write(message,'(a,a,a,a,a,a,i5,a)')ch10,&
&  ' rwwf : BUG -',ch10,&
&  '  The argument formeig should be 0 or 1.',ch10,&
&  '  However, formeig=',formeig,'.'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if

!Check the options
 if(option/=1  .and. option/=2 .and. option/=3 .and. option/=4 .and.&
&   option/=-1 .and. option/=-2                .and. option/=-4)then
  write(message,'(a,a,a,a,a,a,i5,a)')ch10,&
&  ' rwwf : BUG -',ch10,&
&  '  The argument option should be 1, 2, 3, -1 or -2.',ch10,&
&  '  However, option=',option,'.'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if

 unitwf = wff%unwff

 ncid_hdr =unitwf

#if defined HAVE_NETCDF
 !DEBUG
 !if (accesswff == 2) then
 ! ncerr = nf90_Inquire(ncid=ncid_hdr,nDimensions=ndim,nVariables=nvar,nAttributes=natt,unlimitedDimId=uid)
 ! call handle_ncerr(ncerr, " general Inquire ")
 ! write (*,*) 'rwwf : found ndim,nvar,natt,uid = ', ndim,nvar,natt,uid
 !
 !! print all names of dims and vars : equivalent to a ncdump -h
 ! do generic_id=1,ndim
 !  ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=generic_id,name=generic_name)
 !  call handle_ncerr(ncerr," inquire dimension ")
 !  write (*,*) 'Found dimension id ', generic_id, generic_name
 ! end do
 ! do generic_id=1,nvar
 !  ncerr = nf90_Inquire_Variable(ncid=ncid_hdr,varid=generic_id,name=generic_name)
 !  call handle_ncerr(ncerr," inquire variable ")
 !  write (*,*) 'Found variable id ', generic_id, generic_name
 ! end do
 !end if
 !ENDDEBUG
#endif

!#if defined HAVE_ETSF_IO
! if (wff%accesswff == 3) then
!   !Initialize filename in case of ETSF file.
!   file_etsf = trim(wff%fname) // '-etsf.nc'
!   write(message, '(a,a)' ) ' rwwf: created file name for ETSF access ', trim(file_etsf)
!   call wrtout(std_out, message, 'COLL')
! end if
!#endif

 use_f90=0
 if( wff%accesswff == 0   .or.                     &
&   (wff%accesswff ==-1 .and. wff%master==wff%me) )use_f90=1

!------------------------------------------------------------------------
! Read

 if (option/=2 .and. option/=4) then

! Proceed with input

! Read the first record, giving npw1,nspinor1,nband_disk,
! also compute npwso1

!DEBUG     
!write(6,*) "headform =", headform
!ENDDEBUG
     
! headform==0 refers to the current headform
  if(headform>=40 .or. headform==0)then
!DEBUG
!write (6,*) ikpt, "Read npw rec : ", npw1, nspinor1, nband_disk
!ENDDEBUG
   call WffReadNpwRec(ios,ikpt,isppol,nband_disk,npw1,nspinor1,wff)
   npwso1=npw1*nspinor1
   if(ios/=0)then
    write(message,'(a,a,a,a,a,a,i4,a,a,i4,a,a,a)')ch10,&
&    ' rwwf: ERROR -',ch10,&
&    '  Reading option of rwwf. Trying to read',ch10,&
&    '  the (npw,nspinor,nband) record of a wf file, unit=',unitwf,ch10,&
&    '  gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&    '  Action: check your input wf file.'
    call wrtout(6,message,'PERS')
    call leave_new('PERS')
   end if ! ios
  else
!  Old format
   if(use_f90==1)then
    read (unitwf,iostat=ios) npwso1,nband_disk
#if defined MPI_IO
           else if(wff%accesswff==1)then
            call xderiveRRecInit(wff,ios)
            call xderiveReadVal(wff,npwso1)
            call xderiveReadVal(wff,nband_disk)
            call xderiveRRecEnd(wff,ios)
#endif
   end if
   if(ios/=0)then
    write(message,'(a,a,a,a,a,a,i4,a,a,i4,a,a,a)')ch10,&
&    ' rwwf: ERROR -',ch10,&
&    '  Reading option of rwwf. Trying to read',ch10,&
&    '  the (npw,nband) record of a wf file, unit=',unitwf,ch10,&
&    '  gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&    '  Action: check your input wf file.'
    call wrtout(6,message,'PERS')
    call leave_new('PERS')
   end if ! ios
  end if ! headform

!DEBUG
!  write(6,*)' rwwf : option=',option
!  stop
!ENDDEBUG

  if(option==1 .or. option==-2)then
!  Will read the wavefunction and/or kg data, so check npw and nspinor
!  headform==0 refers to the current headform
   if(headform>=40 .or. headform==0)then
!   New format
#ifdef MPI_IO
           npwtot = npw
           call xsum_mpi(npwtot,mpi_enreg%commcart,ios)
           if (npwtot/=npw1) then
              write(message,'(a,a,a,a,a,a,i5,a,i5,a)')ch10,&
                   &     ' rwwf : BUG -',ch10,&
                   &     '  Reading option of rwwf. One should have npwtot=npw1',ch10,&
                   &     '  However, npwtot=',npwtot,', and npw1=',npw1,'.'
              call wrtout(6,message,'PERS')
              call leave_new('PERS')
           end if
           npwso1 = npw*nspinor
#else
    if(npw/=npw1)then
     write(message,'(a,a,a,a,a,a,i5,a,i5,a)')ch10,&
&     ' rwwf : BUG -',ch10,&
&     '  Reading option of rwwf. One should have npw=npw1',ch10,&
&     '  However, npw=',npw,', and npw1=',npw1,'.'
     call wrtout(6,message,'PERS')
     call leave_new('PERS')
    end if
#endif
    if(nspinor/=nspinor1)then
     write(message,'(a,a,a,a,a,a,i5,a,i5,a)')ch10,&
&     ' rwwf : BUG -',ch10,&
&     '  Reading option of rwwf. One should have nspinor=nspinor1',ch10,&
&     '  However, nspinor=',nspinor,', and nspinor1=',nspinor1,'.'
     call wrtout(6,message,'PERS')
     call leave_new('PERS')
    end if
   else
!   Old format
!   This is the reading option : check that npw*nspinor and npwso1 are equal.
    if(npw*nspinor/=npwso1)then
     write(message,'(a,a,a,a,a,a,i5,a,i5,a)')ch10,&
&     ' rwwf : BUG -',ch10,&
&     '  Reading option of rwwf. One should have npwso=npwso1',ch10,&
&     '  However, npwso=',npw*nspinor,', and npwso1=',npwso1,'.'
     call wrtout(6,message,'PERS')
     call leave_new('PERS')
    end if
   end if ! headform
  end if ! option==1 .or. option==2


! Read k+g data
! headform==0 refers to the current headform
  if(headform>=40 .or. headform==0)then

   if( (option==1 .or. option==-2 .or. option==3) .and. &
&       optkg/=0 )then
    if(use_f90==1)then
     read(unitwf,iostat=ios) kg_k(1:3,1:npw)
#if defined MPI_IO
           else if(wff%accesswff==1)then
            call xderiveRRecInit(wff,ios)
            call xderiveRead(wff ,kg_k(1:3,1:npw),3,npw,ios,mpi_enreg%commcart)
            call xderiveRRecEnd(wff,ios)
#endif
#if defined HAVE_NETCDF
    else if (wff%accesswff==2) then
     ncid_hdr = unitwf
!    write (*,*) ' rwwf: reading... ncid_hdr = ', ncid_hdr
!  in netcdf file saved as: kg(3,mpw,nkpt)
     ncerr = nf90_inq_varid(ncid=ncid_hdr,name="kg",varid=kg_id)
     call handle_ncerr(ncerr," inquire kg ")
     !  only get part of kg for this kpoint ikpt and sppol isppol
     ncerr = nf90_get_var(ncid=ncid_hdr,varid=kg_id,values=kg_k,start=(/1,1,ikpt/),count=(/3,npw,1/))
     call handle_ncerr(ncerr," get kg_k ")
#endif
#if defined HAVE_ETSF_IO
    else if (wff%accesswff == 3) then
      ! We read reduced_coordinates_of_plane_waves for one k point.
      wave_folder%reduced_coordinates_of_plane_waves%data2D => kg_k
      wave_folder%red_coord_pw__kpoint_access = ikpt
      call etsf_io_basisdata_get(wff%unwff, wave_folder, lstat, error)
      if (.not. lstat) then
        call etsf_io_low_error_to_str(errmess, error)
        write(message, "(A,A,A,A)") ch10, " rwwf: ERROR -", ch10, errmess(1:min(475, len(errmess)))
        call wrtout(06, message, 'COLL')
      end if
#endif
    end if
   else
    call WffReadSkipRec(ios,1,wff)   ! Skip the record
   end if

   if(ios/=0)then
    write(message,'(a,a,a,a,a,a,i4,a,a,i4,a,a,a)')ch10,&
&    ' rwwf: ERROR -',ch10,&
&    '  Reading option of rwwf. Trying to read',ch10,&
&    '  the k+g record of a wf file, unit=',unitwf,ch10,&
&    '  gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&    '  Action: check your input wf file.'
    call wrtout(6,message,'PERS')
    call leave_new('PERS')
   end if ! ios

  end if ! headform


! Read ground-state eigenvalues
  nband1 = min(nband,nband_disk)
  if(formeig==0)then

   if(option==1 .or. option==3 .or. option==-4)then
    if(use_f90==1)then
     read (unitwf,iostat=ios) eigen(1:nband1)
#if defined MPI_IO
           else if(wff%accesswff==1)then

            call xderiveRRecInit(wff,ios)
            call xderiveRead(wff,eigen,nband1,ios)
            call xderiveRRecEnd(wff,ios)
#endif
    end if

!   The reading of occ should be enabled, BUT taking into account
!   of headform of the disk file : occ was NOT present
!   is the disk files with headform=22
!    read (unitwf) eigen(1:nband1),occ(1:nband1)
   else
    call WffReadSkipRec(ios,1,wff)
   end if ! option

   if(ios/=0)then
    write(message,'(a,a,a,a,a,a,i4,a,a,i4,a,a,a)')ch10,&
&    ' rwwf: ERROR -',ch10,&
&    '  Reading option of rwwf. Trying to read',ch10,&
&    '  an eigenvalue record of a wf file, unit=',unitwf,ch10,&
&    '  gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&    '  Action: check your input wf file.'
    call wrtout(6,message,'PERS')
    call leave_new('PERS')
   end if ! ios

  end if ! formeig

! NetCDF command is universal
!  in netcdf file saved as: eigen((2*mband)**formeig*mband,nkpt,nsppol)
  if(wff%accesswff==2 .and. (option == 1 .or. option == 3 .or. option==-4))then
#if defined HAVE_NETCDF
   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="eigen",varid=eigen_id)
   call handle_ncerr(ncerr," inquire eigen ")
   !  only get part of eigen for this kpoint ikpt and sppol isppol
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=eigen_id,values=eigen,start=(/1,ikpt,isppol/),&
     & count=(/(2*nband1)**formeig*nband1,1,1/))
   call handle_ncerr(ncerr," get eigen ")

!  in netcdf file saved as: occ(mband,nkpt,nsppol)
   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="occ",varid=occ_id)
   call handle_ncerr(ncerr," inquire occ ")
   !  only get part of occ for this kpoint ikpt and sppol isppol
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=occ_id,values=occ,start=(/1,ikpt,isppol/),&
     & count=(/nband1,1,1/))
   call handle_ncerr(ncerr," get occ ")
#endif
#if defined HAVE_ETSF_IO
  else if (wff%accesswff == 3) then
    ! We get eigenvalues and occupations
    electrons_folder%eigenvalues%data1D => eigen
    electrons_folder%eigenvalues__kpoint_access = ikpt
    electrons_folder%eigenvalues__spin_access = isppol
    electrons_folder%occupations%data1D => occ
    electrons_folder%occupations__kpoint_access = ikpt
    electrons_folder%occupations__spin_access = isppol
    call etsf_io_electrons_get(wff%unwff, electrons_folder, lstat, error)
    if (.not. lstat) then
      call etsf_io_low_error_to_str(errmess, error)
      write(message, "(A,A,A,A)") ch10, " rwwf: ERROR -", ch10, errmess(1:min(475, len(errmess)))
      call wrtout(06, message, 'COLL')
    end if
#endif
  end if

! Band loop
  indxx=0
  nband1=min(nband,nband_disk)
  if(nband1>0 .and. option/=-1)then
   do iband=1,nband1

!   Read matrix of eigenvalues
    if(formeig==1)then

     if(option==1 .or. option==3 .or. option==-4)then
      if(use_f90==1)then
       read (unitwf,iostat=ios) eigen(1+indxx:2*nband1+indxx)
#if defined MPI_IO
           else if(wff%accesswff==1)then
            call xderiveRRecInit(wff,ios)
            call xderiveRead(wff,eigen(1+indxx:2*nband1+indxx),2*nband1,ios)
            call xderiveRRecEnd(wff,ios)
#endif
      end if
      indxx=indxx+2*nband1
     else
      call WffReadSkipRec(ios,1,wff) ! Skip the record
     end if

     if(ios/=0)then
      write(message,'(a,a,a,a,a,a,i4,a,a,i4,a,a,a)')ch10,&
&      ' rwwf: ERROR -',ch10,&
&      '  Reading option of rwwf. Trying to read',ch10,&
&      '  a RF eigenvalue record of a wf file, unit=',unitwf,ch10,&
&      '  gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&      '  Action: check your input wf file.'
      call wrtout(6,message,'PERS')
      call leave_new('PERS')
     end if

    end if ! formeig==1

!   Read wavefunctions
    if(option==1 .or. option==-2)then

     ipw=(iband-1)*npwso1+icg
     if(use_f90==1)then
      read(unitwf,iostat=ios) cg(1:2,ipw+1:ipw+npwso1)
#if defined MPI_IO
           else if(wff%accesswff==1)then
           call xderiveRRecInit(wff,ios)
           call xderiveRead(wff,cg(1:2,ipw+1:ipw+npwso1),2,npwso1,ios,mpi_enreg%commcart)
           call xderiveRRecEnd(wff,ios)
#endif
     end if

    else if (option/=-4) then
     call WffReadSkipRec (ios,1,wff) ! Skip the record
    end if ! option

    if(ios/=0)then
     write(message,'(a,a,a,a,a,a,i4,a,a,i4,a,a,a)')ch10,&
&     ' rwwf: ERROR -',ch10,&
&     '  Reading option of rwwf. Trying to read',ch10,&
&     '  a RF wf record of a wf file, unit=',unitwf,ch10,&
&     '  gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&     '  Action: check your input wf file.'
     call wrtout(6,message,'PERS')
     call leave_new('PERS')
    end if

   end do ! nband

#if defined HAVE_NETCDF
!  in the preceding do loop nothing is done in the netCDF case,
!    so here read all wf_k as a block
!    for the moment use icg offset and save wfk as one big array, like in abinit.
!    later could separate dimensions and make a sliced extraction of the wfk
  if(wff%accesswff==2 .and. (option == 1 .or. option == -2))then
!  in netcdf file saved as: cg(2,mpw,nspinor,mband,nkpt,nsppol)
   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="cg",varid=cg_id)
   call handle_ncerr(ncerr," inquire cg ")
   !  only get part of cg for this kpoint ikpt and sppol isppol
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=cg_id,values=cg(:,icg+1:icg+npw1*nspinor1*nband1),start=(/1,1,1,1,ikpt,isppol/),&
     &  count=(/2,npw1,nspinor1,nband1,1,1/))
   call handle_ncerr(ncerr," get cg ")
#endif
#if defined HAVE_ETSF_IO
  else if (wff%accesswff == 3) then
     ! We get the coefficients_of_wavefunctions
!!$     main_folder%coefficients_of_wavefunctions%data2D => &
!!$          & cg(:, icg + 1:icg + npw * nspinor * nband)
     ! With g95, the association done above sometime leads to segfaults.
     ! So we allocate a temporary array to store the wfs of our kpt.
    allocate(main_folder%coefficients_of_wavefunctions%data2D(2, &
         & npw * nspinor * nband))
    main_folder%wfs_coeff__kpoint_access = ikpt
    main_folder%wfs_coeff__spin_access = isppol
    main_folder%wfs_coeff__number_of_states = nband
    main_folder%wfs_coeff__number_of_coefficients = npw * nspinor
    call etsf_io_main_get(wff%unwff, main_folder, lstat, error)
    ! Now we copy our values and deallocate the temporary array.
    cg(:, icg + 1:icg + npw * nspinor * nband) = &
         & main_folder%coefficients_of_wavefunctions%data2D
    deallocate(main_folder%coefficients_of_wavefunctions%data2D)
    if (.not. lstat) then
      call etsf_io_low_error_to_str(errmess, error)
      write(message, "(A,A,A,A)") ch10, " rwwf: ERROR -", ch10, errmess(1:min(475, len(errmess)))
      call wrtout(06, message, 'COLL')
    end if
#endif
#if defined HAVE_NETCDF
  end if
#endif


  end if ! nband>1

! If fewer than all bands were read
! wind disk file forward to end of bands for this k point.
! Will have to fill the non-filled bands outside of this routine ...
  if (nband<nband_disk .or. option==-1) then
   nrec=(formeig+1)*(nband_disk-nband)
   if(option==-1)nrec=(formeig+1)*nband_disk
   call WffReadSkipRec(ios,nrec,wff)
  end if

!------------------------------------------------------------------------
! Write

 else if (option==2 .or. option==4) then
     call writewf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,nband,nband_disk,&
& npw,nspinor,occ,option,optkg,wff)

 end if ! option

 call timab(270+tim_rwwf,2,tsec)

!DEBUG
! write(6,*)' rwwf : exit   ',wff%offwff
! stop
!ENDDEBUG

end subroutine rwwf
!!***




!{\src2tex{textfont=tt}}
!!****f* ABINIT/writewf
!! NAME
!! writewf
!!
!! FUNCTION
!!  This subroutine write the block of records
!!  related to one k point, and one spin-polarization, that
!!  contains the wavefunctions (as well as the eigenvalues and occupations).
!!  The disk file unitwf should have been prepared
!!  outside of this routine, in order to read or write the correct records.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, MVer,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cg(2,npw*nspinor*mband)=planewave coefficients of wavefunctions,
!!  eigen((2*mband)**formeig *mband)=array for holding eigenvalues (hartree)
!!  formeig=format of the eigenvalues
!!   0 => vector of eigenvalues
!!   1 => hermitian matrix of eigenvalues
!!  headform=format of the header of the wf file, also governing the k block format
!!   in case headform=0, use the default (current) format and headform
!!  icg=shift to be given to the location of the cg array
!!  ikpt=index of current k point (only needed for error message)
!!  isppol=spin polarization currently treated (only needed for error message)
!!  kg_k(3,optkg*npw)=k+g data  (only if optkg==1)
!!  mband=maximum number of bands (dimension of cg, eigen and occ)
!!  mcg=dimention of cg
!!  nband=number of bands actually in cg, eigen and occ
!!   (must be larger or equal to nband_disk, only nband_disk bands are written)
!!  nband_disk=number of bands on the disk file
!!  npw=number of plane waves
!!  nspinor=number of spinorial components of the wavefunctions
!!  occ(mband)=array for holding eigenvalues (hartree)
!!  option= 2 for write,
!!          4 for writing a file containing only eigenvalues and occupations
!!  optkg= if 1 , read or write kg_k ; if 0, do not care about kg_k
!!  wff=struct info for wavefunction
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  WARNING : occ is not read in the present status of this routine
!!  WARNING : skipping k-bloks is also done in the randac subroutine
!!  WARNING : reading the two first records is also done in the rdnpw routine
!!  WARNING : writing the two first records is also done in the vtowfk3 routine
!!
!! PARENTS
!!      outkss,rwwf
!!
!! CHILDREN
!!      etsf_io_electrons_put,etsf_io_low_error_to_str,etsf_io_main_put
!!      etsf_io_basisdata_put,handle_ncerr,leave_new,wffwritenpwrec,wrtout
!!      xderivewrecend,xderivewrecinit,xderivewrite
!!
!! SOURCE
subroutine writewf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,nband,nband_disk,&
& npw,nspinor,occ,option,optkg,wff)

 use defs_basis
 use defs_datatypes
#if defined HAVE_NETCDF
 use netcdf
#endif
#if defined HAVE_ETSF_IO
 use etsf_io
#endif

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13io_mpi, except_this_one => writewf
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: formeig,headform,icg,ikpt,isppol,mband,mcg,nband,npw
 integer,intent(in) :: nband_disk
 integer,intent(in) :: nspinor,option,optkg
 integer,intent(in), target :: kg_k(3,optkg*npw)
 real(dp),intent(in), target :: cg(2,mcg),eigen((2*mband)**formeig*mband),occ(mband)
  type(wffile_type),intent(inout) :: wff
  type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
 integer :: iband,ii,indxx,ios,ipw,npwso,npwso1,nspinor1 = 0 ,unitwf,nband2
 integer :: use_f90,nrec
 character(len=500) :: message
 real(dp) :: tsec(2)
 integer :: ncid_hdr

!no_abirules
#if defined MPI
           integer :: ndim,nvar,natt,uid
           integer :: generic_id
           character(len=500) :: generic_name
           character(len=10) :: tmpstr
#endif
#if defined HAVE_NETCDF
 integer :: ncerr, kg_id, eigen_id, occ_id, cg_id
#endif
#if defined HAVE_ETSF_IO
 character(len=fnlen) :: file_etsf
 character(len=etsf_io_low_error_len) :: errmess
 type(etsf_main) :: main_folder
 type(etsf_basisdata),target :: wave_folder
 type(etsf_electrons),target :: electrons_folder
 logical :: lstat
 type(etsf_io_low_error) :: error
#endif
#if defined MPI_IO
  integer :: numproc,npwtot,ierr
  integer, allocatable :: arrkg(:), arrcg(:)
#endif


!Check the options
 if(option/=2 .and. option/=4)then
  write(message,'(a,a,a,a,a,a,i5,a)')ch10,&
&  ' writewf : BUG -',ch10,&
&  '  The argument option should be 2 or 4.',ch10,&
&  '  However, option=',option,'.'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if

 unitwf = wff%unwff

 ncid_hdr =unitwf

 use_f90=0
 if( wff%accesswff == 0   .or.                     &
&   (wff%accesswff ==-1 .and. wff%master==wff%me) )use_f90=1

! This is the writing option : check that nband_disk is not larger than nband.
  if(nband<nband_disk)then
   write(message,'(a,a,a,a,a,a,i5,a,i5,a)')ch10,&
&   ' rwwf : BUG -',ch10,&
&   '  Writing option of rwwf. One should have nband<=nband_disk',ch10,&
&   '  However, nband=',nband,', and nband_disk=',nband_disk,'.'
   call wrtout(6,message,'PERS')
   call leave_new('PERS')
  end if

!DEBUG
!write(6,*)' rwwf : write npw,nspinor,nband_disk=',npw,nspinor,nband_disk
!ENDDEBUG


!  not modified for netCDF: no need to add writing of nband_disk,npw,nspinor
#if defined MPI_IO
  call MPI_COMM_SIZE(mpi_enreg%commcart,numproc,ierr)
  npwtot = npw
  call xsum_mpi(npwtot,mpi_enreg%commcart,ierr)
  call WffWriteNpwRec_cs(ios,mpi_enreg,nband_disk,npwtot,nspinor,wff)
#else
  call WffWriteNpwRec(ios,nband_disk,npw,nspinor,wff)
#endif

  if(optkg/=0 .and. option/=4)then
   if(use_f90==1)then
    write(unitwf) kg_k(1:3,1:optkg*npw)
#if defined MPI_IO
           else if(wff%accesswff==1)then
        call xderiveWRecInit_cs(wff,ios,mpi_enreg%me_cart_2d)
        call xderiveWrite(wff,kg_k,3,optkg*npw,ios,mpi_enreg%commcart)
        call xderiveWRecEnd_cs(wff,ios,mpi_enreg%me_cart_2d)
#endif
#if defined HAVE_NETCDF
   else if (wff%accesswff==2)then

!  in netcdf file saved as: kg(3,mpw,nkpt)
! Write data: all dimensions at once.
!   dimensions should have been created in hdr_io_netcdf
!DEBUG
!write (6,*) ' rwwf : writing... about to inquire for kg'
!write (6,*) ' rwwf: ncid_hdr = ', ncid_hdr
!ENDDEBUG
    ncerr = nf90_inq_varid(ncid=ncid_hdr,name="kg",varid=kg_id)
    call handle_ncerr(ncerr," inquire kg ")
    ncerr = nf90_put_var(ncid=ncid_hdr,varid=kg_id,values=kg_k,start=(/1,1,ikpt/),count=(/3,npw,1/))
    call handle_ncerr(ncerr," fill kg")
#endif
#if defined HAVE_ETSF_IO
   else if (wff%accesswff == 3) then
    ! We write reduced_coordinates_of_plane_waves
    wave_folder%reduced_coordinates_of_plane_waves%data2D => kg_k
    wave_folder%red_coord_pw__kpoint_access = ikpt
    wave_folder%red_coord_pw__number_of_coefficients = npw
    call etsf_io_basisdata_put(wff%unwff, wave_folder, lstat, error)
    if (.not. lstat) then
      call etsf_io_low_error_to_str(errmess, error)
      write(message, "(A,A,A,A)") ch10, " rwwf: ERROR -", ch10, errmess(1:min(475, len(errmess)))
      call wrtout(06, message, 'COLL')
    end if
#endif
   end if
   ! end wff%accesswff if
  else
   if(use_f90==1)then
    write(unitwf)                   ! Still skip the record
#if defined MPI_IO
           else if(wff%accesswff==1)then
        call xderiveWRecInit_cs(wff,ios,mpi_enreg%me_cart_2d)
        call xderiveWRecEnd_cs(wff,ios,mpi_enreg%me_cart_2d)
#endif
   end if
  end if

#if defined HAVE_NETCDF
  if(wff%accesswff==2)then
!
! NetCDF commands are universal
!  in netcdf file saved as: eigen((2*mband)**formeig*mband,nkpt,nsppol)
   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="eigen",varid=eigen_id)
   call handle_ncerr(ncerr," inquire eigen ")
   !  only get part of eigen for this kpoint ikpt and sppol isppol
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=eigen_id,values=eigen,start=(/1,ikpt,isppol/),&
     & count=(/(2*nband_disk)**formeig*nband_disk,1,1/))
   call handle_ncerr(ncerr," put eigen ")

!  in netcdf file saved as: occ(mband,nkpt,nsppol)
   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="occ",varid=occ_id)
   call handle_ncerr(ncerr," inquire occ ")
   !  only get part of occ for this kpoint ikpt and sppol isppol
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=occ_id,values=occ,start=(/1,ikpt,isppol/),&
     & count=(/nband_disk,1,1/))
   call handle_ncerr(ncerr," put occ ")

   if (option /= 4) then
!  in netcdf file saved as: cg(2,mpw,nspinor,mband,nkpt,nsppol)
    ncerr = nf90_inq_varid(ncid=ncid_hdr,name="cg",varid=cg_id)
    call handle_ncerr(ncerr," inquire cg ")
    !  only get part of cg for this kpoint ikpt and sppol isppol
    ncerr = nf90_put_var(ncid=ncid_hdr,varid=cg_id,values=cg(:,icg+1:icg+npw*nspinor1*mband),start=(/1,1,1,1,ikpt,isppol/),&
      &  count=(/2,npw,nspinor,nband_disk,1,1/))
    call handle_ncerr(ncerr," put cg ")
   end if
#if defined HAVE_ETSF_IO
  else if (wff%accesswff == 3) then
    ! FIXME: the case formeig /= 0 is not supported.
    ! We write eigenvalues, occupations and coefficients_of_wavefunctions
    electrons_folder%eigenvalues%data1D => eigen
    electrons_folder%eigenvalues__kpoint_access = ikpt
    electrons_folder%eigenvalues__spin_access = isppol
    electrons_folder%eigenvalues__number_of_states = mband
    electrons_folder%occupations%data1D => occ
    electrons_folder%occupations__kpoint_access = ikpt
    electrons_folder%occupations__spin_access = isppol
    electrons_folder%occupations__number_of_states = mband
    call etsf_io_electrons_put(wff%unwff, electrons_folder, lstat, error)
    if (.not. lstat) then
      call etsf_io_low_error_to_str(errmess, error)
      write(message, "(A,A,A,A)") ch10, " rwwf: ERROR -", ch10, errmess(1:min(475, len(errmess)))
      call wrtout(06, message, 'COLL')
    end if
    if (option /= 4) then
!!$       main_folder%coefficients_of_wavefunctions%data2D => &
!!$            & cg(:, icg + 1:icg + npw * nspinor * nband)
       ! See the read access, this direct association sometime leads to
       ! segfaults with g95, so we use a temporary array.
       allocate(main_folder%coefficients_of_wavefunctions%data2D(2, &
            & npw * nspinor * nband))
       main_folder%coefficients_of_wavefunctions%data2D = &
            & cg(:, icg + 1:icg + npw * nspinor * nband)
       main_folder%wfs_coeff__kpoint_access = ikpt
       main_folder%wfs_coeff__spin_access = isppol
       main_folder%wfs_coeff__number_of_states = nband
       main_folder%wfs_coeff__number_of_coefficients = npw * nspinor
       call etsf_io_main_put(wff%unwff, main_folder, lstat, error)
       deallocate(main_folder%coefficients_of_wavefunctions%data2D)
       if (.not. lstat) then
          call etsf_io_low_error_to_str(errmess, error)
          write(message, "(A,A,A,A)") ch10, " rwwf: ERROR -", ch10, errmess(1:min(475, len(errmess)))
          call wrtout(06, message, 'COLL')
       end if
    end if
#endif
  end if
#endif

  if(formeig==0)then
   if(use_f90==1)then
!   Write eigenvalues
    write(unitwf) (eigen(iband),iband=1,nband_disk),&
&                 (occ(iband),iband=1,nband_disk)
#if defined MPI_IO
           else if(wff%accesswff==1)then
        if (mpi_enreg%me_cart_2d == 0) then
            call xderiveWRecInit(wff,ios)
            call xderiveWrite(wff,eigen,nband_disk,ios)
            call xderiveWrite(wff,occ,nband_disk,ios)
            call xderiveWRecEnd(wff,ios)
        end if
        !call xcast_mpi(wff%offwff,0,mpi_enreg%commcart,ios)
        call MPI_BCAST(wff%offwff,1,MPI_INTEGER8,0,mpi_enreg%commcart,ios)
#endif
   end if

   if (option /= 4) then
    npwso=npw*nspinor
    do iband=1,nband_disk
     ipw=(iband-1)*npwso+icg
     if(use_f90==1)then
      write(unitwf) cg(1:2,ipw+1:ipw+npwso)
#if defined MPI_IO
            else if(wff%accesswff==1)then
              call xderiveWRecInit_cs(wff,ios,mpi_enreg%me_cart_2d)
              call xderiveWrite(wff,cg(1:2,ipw+1:ipw+npwso),2,npwso,ios,mpi_enreg%commcart)
              call xderiveWRecEnd_cs(wff,ios,mpi_enreg%me_cart_2d)
#endif
     end if
    end do
   end if

  else if(formeig==1)then
   npwso=npw*nspinor
   nband2=2*nband_disk
   do iband=1,nband_disk
    ipw=(iband-1)*npwso+icg
    ii=(iband-1)*nband2
    if(use_f90==1)then
     write(unitwf) eigen(1+ii:nband2+ii)
     if (option /= 4) then
      write(unitwf) cg(1:2,1+ipw:npwso+ipw)
     end if
#if defined MPI_IO
           else if(wff%accesswff==1)then
            call xderiveWRecInit(wff,ios)
            call xderiveWrite(wff,eigen(ii:ii+nband2),nband2,ios)
            call xderiveWRecEnd(wff,ios)
            if (option /= 4) then
             call xderiveWRecInit(wff,ios)
              call xderiveWrite(wff,cg(1:2,ipw+1:ipw+npwso),2,npwso,ios,arrcg)
             call xderiveWRecEnd(wff,ios)
            end if

#endif
    end if

   end do
  end if ! formeig==0 or 1
#ifdef MPI_IO
  !deallocate(arrkg,arrcg)
#endif
end subroutine writewf
!!***
