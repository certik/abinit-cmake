!{\src2tex{textfont=tt}}
!!****p* ABINIT/mrgscr
!! NAME
!! mrgscr
!!
!! PROGRAM
!! This codes reads partial _SCR files for different q points, then merges them into a single file
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (R. Shaltaf, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  Main program
!!
!! OUTPUT
!!  Only checking and writing
!! 
!! NOTES
!!
!! If the number of SCR files to be merged is equal to 1, the program only checks 
!! the integrity of the file and report the list of q-points that are missing.
!! Note that for GW calculations the list of q-points required depend on the k-mesh 
!! used during the calculation of the KSS file. The list of missing q-points 
!! is calculated assuming that the same mesh is used in sigma.
!!
!! PARENTS
!!
!! CHILDREN
!!      findnq,findq,hdr_check,hdr_clean,hdr_io,hdr_io_netcdf,herald,identk
!!      leave_new,mati3inv,memerr,metric,rdscr,testscr,wrscr,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

program mrgscr

 use defs_basis
 use m_gwdefs, only : GW_TOLQ
 use defs_datatypes
 use defs_infos


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_14iowfdenpot
 use interfaces_15gw
!End of the abilint section

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
 character(len=50),parameter :: FILE__='mrgscr.F90'
!scalars
 integer,parameter :: localrdwf=0,unitem1=41
 integer :: fform,fform1,ifile,ig,ii,io,ios,iq,iq0,iq1,iqW,is,istat,isym,jj
 integer :: nbnds0,nbnds1,nfiles,nkbz,nkbzX,nkibz,nomega0,nomega1,npweps0
 integer :: npweps1,npwwfn0,npwwfn1,nq_merge,nq_miss0,nqibz,nqibz0,nqibz1,nqlwl1,nsym
 integer :: optfil,rdwr,readnetcdf,restart,restartpaw,timrev
 integer :: unt_out,unt1
 real(dp) :: ucvol
 logical :: avoid_zero,found,is_old,use_antiferro
 character(len=20) :: temp_file
 character(len=24) :: codename
 character(len=500) :: msg
 type(Datafiles_type) :: Dtfil0,Dtfil1,Dtfil_out
 type(Dataset_type) :: Dtset
 type(Hdr_type) :: Hdr0,Hdr1,Hdr_scr
 type(MPI_type) :: MPI_enreg
!arrays
 integer,pointer :: gvec_p(:,:)
 integer,allocatable :: gvec0(:,:),gvec1(:,:),ktab(:),ktabi(:),ktabo(:)
 integer,allocatable :: merge_table(:,:),symafm(:),symrec(:,:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),qdiff(3),rmet(3,3),rprimd0(3,3)
 real(dp),pointer :: qibz_p(:,:),qlwl_p(:,:)
 real(dp),allocatable :: kbz(:,:),kibz(:,:),q_merge(:,:),qibz(:,:),qibz0(:,:)
 real(dp),allocatable :: qibz1(:,:),qlwl(:,:),qmiss0(:,:),tnons2(:,:),wtk(:),qsmall_dum(:,:)
 complex(dpc),pointer :: omega_p(:)
 complex(gwpc),allocatable :: lwing_dum(:,:,:),uwing_dum(:,:,:)
 complex(dpc),allocatable :: omega0d(:),omega1d(:)
 complex(gwpc),allocatable :: epsm1(:,:,:,:),omega1(:)
 logical,allocatable :: foundq(:),repetition(:,:)
 character(len=80) :: titem0(2),titem1(2)
 character(len=fnlen),allocatable :: filenames(:)
 character(len=fnlen) :: fname_out,fname1

! *************************************************************************

 codename='MRGSCR'//repeat(' ',18)
!
!Write greetings, and read number of files.
 call herald(codename,abinit_version,std_out)

 write(msg,'(a)')' Enter the number of files to merge :'
 call wrtout(std_out,msg,'COLL')
 read(*,*)nfiles
 write(msg,'(a,i4,a)')' Will read ',nfiles,' partial screening files'
 call wrtout(std_out,msg,'COLL')
 allocate(filenames(nfiles))

 if (nfiles>1) then
! 
! Read files to be merged and check for existence
  write(msg,'(a)')' Enter name for the final output file'
  call wrtout(std_out,msg,'COLL')
  read(*,*) temp_file
  Dtfil_out%filnam_ds(4)=TRIM(temp_file) 

! Here I have a problem with the name
  write(msg,'(2a)')' Partial screening files will be merged into ',TRIM(Dtfil_out%filnam_ds(4))//'_SCR'
  call wrtout(std_out,msg,'COLL')
  inquire(file=TRIM(temp_file)//'_SCR',exist=is_old)
  if (is_old) then 
   write(msg,'(9a)')ch10,' mrgscr : ERROR - ',ch10,&
&   ' found that file ',TRIM(temp_file)//'_SCR',' already exists, cannot overwrite ',&
&   ch10,' use a different name or (if you prefer) remove ',TRIM(temp_file) 
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
  end if 

  do ifile=1,nfiles
   write(msg,'(a,i4)')' Enter name for partial screening file number: ',ifile
   call wrtout(std_out,msg,'COLL') 
   read(*,*) temp_file 
   filenames(ifile)=TRIM(temp_file)
   write(msg,'(2a)')' Will read file ',TRIM(filenames(ifile))
   call wrtout(std_out,msg,'COLL') 
   inquire(file=TRIM(filenames(ifile)),exist=is_old)
   if (.not.is_old) then
    write(msg,'(6a)')ch10,' ERROR -',ch10,&
&    ' file ',TRIM(filenames(ifile)),' does not exist'
    call wrtout(std_out,msg,'COLL') 
    call leave_new('COLL')
   end if
  end do

 else if (nfiles==1) then
! 
! Single mode
  read(*,*) temp_file
  filenames(1)=TRIM(temp_file)
  write(msg,'(7a)')ch10,&
  ' Running single-file mode:',ch10,&
&  ' mrgscr will check the integrity of file: ',TRIM(filenames(1)),ch10,&
&  ' reporting the list of q-points that are missing '
  call wrtout(std_out,msg,'COLL')
 else
  write(msg,'(4a,i4)')ch10,' ERROR -',ch10,&
&  ' number of files should be >0, while it is : ',nfiles
  call wrtout(std_out,msg,'COLL') 
  call leave_new('COLL')
 end if


 do ifile=1,nfiles
! 
! For each file read dimensions, list of q-points and do some check
  open(unit=unitem1,file=TRIM(filenames(ifile)),status='old',form='unformatted',iostat=ios)
  if (ios/=0) then 
   write(msg,'(6a)')ch10,' mrgscr : ERROR -',ch10,&
&   ' opening file ',TRIM(filenames(ifile)),' as old'
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
  end if  
  readnetcdf=0  ; Dtset%accesswff=0 ; rdwr=1

  if (ifile==1) then
!  
!  Special treatment for the first file: Read the header, Hdr0, and calculate the complete 
!  list of q-points. Hdr0 will be kept in memory for the sake of comparison, 
   write(msg,'(7a)') ch10,&
   ' reading and analyzing the header of file ',TRIM(filenames(ifile)),ch10,&
&   ' Calculating the complete list of q-points',ch10,ch10
   call wrtout(std_out,msg,'COLL')

   if      (readnetcdf == 0) then 
    call hdr_io       (fform,Hdr0,rdwr,unitem1)
   else if (readnetcdf == 1) then 
    call hdr_io_netcdf(fform,Hdr0,rdwr,unitem1)
   end if 

   read(unitem1) titem0
   read(unitem1) npweps0,npwwfn0,nbnds0,nqibz0,nomega0
   allocate(gvec0(3,npweps0),qibz0(3,nqibz0),omega0d(nomega0))
   read(unitem1) gvec0(:,1:npweps0)
   read(unitem1) qibz0(:,1:nqibz0)

   write(msg,'(a)')' List of q points : '
   call wrtout(std_out,msg,'COLL')
   do iq0=1,nqibz0
    write(msg,'(i8,3f12.6)') iq0,qibz0(:,iq0)
    call wrtout(std_out,msg,'COLL')
   end do

   read(unitem1) omega0d(1:nomega0)
   close(unitem1)

   rprimd0(:,:)=Hdr0%rprimd(:,:)
   call metric(gmet,gprimd,-1,rmet,rprimd0,ucvol)

   nsym= Hdr0%nsym
   nkibz=Hdr0%nkpt
   
   allocate(symrec(3,3,nsym),symafm(nsym))
   do isym=1,Hdr0%nsym
    call mati3inv(Hdr0%symrel(:,:,isym),symrec(:,:,isym))
    symafm(isym)=Hdr0%symafm(isym)
   end do

   timrev=2 !this should be read from kptopt
   nkbzX=nkibz*nsym*2

   allocate(kibz(3,nkibz),kbz(3,nkbzX))
   allocate(wtk(nkibz),ktab(nkbzX),ktabi(nkbzX),ktabo(nkbzX))
   kibz(:,:)=Hdr0%kptns(:,:)
   use_antiferro = (Hdr0%nspden==2 .and. Hdr0%nsppol==1)

!  FIXME here there is a portability problem with g95
   call identk(kibz,nkibz,nkbzX,nsym,timrev,symrec,symafm,use_antiferro,kbz,ktab,ktabi,ktabo,nkbz,wtk)

   call findnq(nkbz,kbz,nsym,symrec,nqibz,timrev)

!  Here I have a problem if the small point is set equal to zero
   allocate(qibz(3,nqibz)) ; avoid_zero=.TRUE.
   call findq(nkbz,kbz,nsym,symrec,gprimd,nqibz,qibz,timrev,avoid_zero)

   deallocate(kibz,kbz,wtk,ktab,ktabi,ktabo)

   allocate(foundq(nqibz),repetition(nfiles,nqibz),q_merge(3,nqibz))
   foundq(:)=.false. ; repetition(:,:)=.false.
!  
!  merge_table gives, for each q-point to be merged, the file where it is located
!  as well as its sequential index. It might be useful to do the merge point-by-point 
!  thus avoiding the allocation of the entire epsm1 arrays
   allocate(merge_table(nqibz,2))
   merge_table(:,:)=0

!  Find the needed q points
!  initialize the number of q-points missing/to be merged
   nq_miss0=0  ; nq_merge=0
   allocate(qmiss0(3,nqibz-nqibz0))
   do iq=1,nqibz
    do iq0=1,nqibz0
     qdiff(:)=qibz(:,iq)-qibz0(:,iq0)
     if (all(abs(qdiff(:))<GW_TOLQ)) then
      foundq(iq)=.true.
      nq_merge=nq_merge+1
      q_merge(:,nq_merge)=qibz0(:,iq0)
      merge_table(nq_merge,1)=ifile
      merge_table(nq_merge,2)=iq0

      write(msg,'(2a,3f12.6)')ch10,' ... found q-point :',q_merge(:,nq_merge)
      call wrtout(std_out,msg,'COLL')

      exit 
     end if
    end do
    if (.not.foundq(iq)) then
     nq_miss0=nq_miss0+1
     qmiss0(:,nq_miss0)=qibz(:,iq)
    end if
   end do

   if (nfiles==1) then 
!   Only checking, not merging
    write(msg,'(4a,i4,4a)')ch10,' file ',TRIM(filenames(1)),&
&    ' contains',nqibz0,' q-points in the irreducible wedge (IBZ)',&
&    ch10,' q-points [reciprocal lattice units]:',ch10
    call wrtout(std_out,msg,'COLL')
    do iq0=1,nqibz0
     write(msg,'(i5,3f12.6)')iq0,(qibz0(ii,iq0),ii=1,3)
     call wrtout(std_out,msg,'COLL')
    end do
    if (nqibz0<nqibz) then
     write(msg,'(2a,i4,4a)')ch10,' There are ',nq_miss0,&
&     ' points missing in the file',ch10,' List of missing points: ',ch10
     call wrtout(std_out,msg,'COLL')
     do jj=1,nq_miss0
      write(msg,'(i5,3f12.6)')jj,(qmiss0(ii,jj),ii=1,3)
      call wrtout(std_out,msg,'COLL')
     end do
    else if (nqibz0==nqibz) then
     write(msg,'(3a)')ch10,' File is complete, nothing to do !',ch10
     call wrtout(std_out,msg,'COLL') 
     exit
    else if (nqibz0>nqibz) then
     write(msg,'(6a)')ch10,' BUG -',ch10,&
&     ' nqibz0>nqibz in ',TRIM(filenames(ifile)),ch10
     call wrtout(std_out,msg,'COLL') 
     call leave_new('COLL')
    end if
    deallocate(qmiss0)
   end if ! nfiles==1

   write(msg,'(2a)')ch10,' checking and analysing the first file finished'
   call wrtout(std_out,msg,'COLL')

  else if (ifile>1) then  
!  
!  Read and compare the headers of the other files 
   if      (readnetcdf==0) then 
    call hdr_io       (fform1,Hdr1,rdwr,unitem1) 
   else if (readnetcdf==1) then 
    call hdr_io_netcdf(fform1,Hdr1,rdwr,unitem1)
   end if

   read(unitem1) titem1 
   read(unitem1) npweps1,npwwfn1,nbnds1,nqibz1,nomega1
   allocate(gvec1(3,npweps1),qibz1(3,nqibz1),omega1d(nomega1))
   read(unitem1) gvec1(1:3,1:npweps1)
   read(unitem1) qibz1(1:3,1:nqibz1)
   read(unitem1) omega1d(1:nomega1)
   close(unitem1)

   write(msg,'(5a)' )ch10,&
&   ' checking the header of file ',TRIM(filenames(ifile)),&
&   ', comparing to the header of file ',TRIM(filenames(1))
   call wrtout(std_out,msg,'COLL')
   write(msg,&
&   '(a1,80a,2a1,10x,a,3a1,13x,a,25x,a,a1,8x,19a,25x,12a,a1)' )&
&   ch10,('=',ii=1,80),ch10,ch10,' checking file headers for consistency -',&
&   (ch10,ii=1,3),TRIM(filenames(1)),TRIM(filenames(ifile)),ch10,('-',ii=1,19),('-',ii=1,12),ch10
   call wrtout(std_out,msg,'COLL')

   call hdr_check(fform,fform1,Hdr0,Hdr1,'COLL',restart,restartpaw)
   if (restart==0) then
!   FFT grid might be q-point dependent so we stop only when restart==0
    write(msg,'(4a)')ch10,' ERROR -',ch10,' headers are not consistent'
    call wrtout(std_out,msg,'COLL') 
    call leave_new('COLL')
   end if
   if (npweps1/=npweps0) then
    write(msg,'(4a)')ch10,' ERROR -',ch10,& 
&    ' size of epsilon matrix is not equal in the two files'
    call wrtout(std_out,msg,'COLL') 
    call leave_new('COLL')
   end if
   if (npwwfn1/=npwwfn0) then
    write(msg,'(5a)')ch10,' COMMENT -',ch10,&
&    ' number of planewaves used for wavefunctions is not equal',ch10
    call wrtout(std_out,msg,'COLL')
   end if
   if (nbnds1/=nbnds0) then
    write(msg,'(5a)')ch10,' COMMENT -',ch10,&
&    ' number of planewaves used for wavefunctions is not equal',ch10
    call wrtout(std_out,msg,'COLL')
   end if
   if (nomega1/=nomega0) then
    write(msg,'(4a)')ch10,' ERROR -',ch10,&
&    ' number of omega is not equal in the two files '
    call wrtout(std_out,msg,'COLL') 
    call leave_new('COLL')
   end if
   if ( any(abs (real (omega1d-omega0d)) >tol6) .or. &
   any(abs (aimag(omega1d-omega0d)) >tol6)      &
   ) then
    write(msg,'(5a)')ch10,' ERROR -',ch10,&
&    ' frequencies in the first file differ from actual frequencies : ',ch10
    call wrtout(std_out,msg,'COLL') 
    do io=1,nomega0
     write(*,*) omega0d(io),omega1d(io)
    end do 
    call leave_new('COLL')
   end if
   if (any(abs (gvec0(:,:)-gvec1(:,:)) > 0)) then 
    write(msg,'(5a)')ch10,' ERROR -',ch10,&
&    ' incompatible G vector list found',ch10
    call wrtout(std_out,msg,'COLL') 
    do ig=1,npweps0 
     write(*,*)gvec0(:,ig),gvec1(:,ig)
    end do
    call leave_new('COLL')
   end if
!  
!  Now check if the q points reported in this file are correct
   write(msg,'(4a)')ch10,' analyzing existing q points in file ',TRIM(filenames(ifile)),ch10
   call wrtout(std_out,msg,'COLL')

   write(msg,'(a)')' List of q points : '
   call wrtout(std_out,msg,'COLL')
   do iq1=1,nqibz1
    write(msg,'(i8,3f12.6)') iq1,qibz1(:,iq1)
    call wrtout(std_out,msg,'COLL')
   end do

   do iq=1,nqibz
    do iq1=1,nqibz1
     qdiff(:)=qibz(:,iq)-qibz1(:,iq1)
     if (all(abs(qdiff(:))<GW_TOLQ)) then
      if (.not.foundq(iq)) then
       foundq(iq)=.true.
       nq_merge=nq_merge+1
       q_merge(:,nq_merge)=qibz1(:,iq1)
       merge_table(nq_merge,1)=ifile
       merge_table(nq_merge,2)=iq1

       write(msg,'(2a,3f12.6)')ch10,' ... found q-point :',q_merge(:,nq_merge)
       call wrtout(std_out,msg,'COLL')

      else 
!      A q-point is present at least twice
       repetition(ifile,iq1)=.true.
       write(msg, '(4a,3f12.6,4a)')ch10,' WARNING -',ch10,&
&       ' q point',qibz1(:,iq1),' already found in another file',&
&       ch10,' This point will be skipped ',ch10
       call wrtout(std_out,msg,'COLL')
      end if
      if (foundq(iq)) exit
     end if
    end do 
   end do 
   call hdr_clean(Hdr1)
   deallocate(gvec1,qibz1,omega1d)
  end if ! single or multi mode

 end do ! End loop over ifile
!
 if (nfiles>1) then
! Write list of q-points to be merged
  write(msg,'(3a)') ch10,' q-points will be saved in the final file in the following order:',ch10
  call wrtout(std_out,msg,'COLL')
  do jj=1,nq_merge
   write(msg,'(i5,3f12.6)')jj,(q_merge(ii,jj),ii=1,3)
   call wrtout(std_out,msg,'COLL')
  end do
  if (nq_merge<nqibz) then
   write(msg,'(6a)')ch10,' ERROR -',ch10,&
&   ' the partial screening files do not include enough q points',ch10,&
&   ' check the content of the files'
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
  end if
  if (nq_merge>nqibz) then 
   write(msg,'(3a)')' mrgscr BUG -',ch10,'nq_merge /= nqibz'
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
  end if
  deallocate(qibz0)

  write(msg,'(3a)' )ch10,' start merging ...',ch10
  call wrtout(std_out,msg,'COLL')

  iqW=0 ; Dtset%nqptdm=0
  do ifile=1,nfiles
   write(msg,'(4a)') ch10,' merging file ',TRIM(filenames(ifile)),ch10
   call wrtout(std_out,msg,'COLL')
   Dtfil1%filscr=TRIM(filenames(ifile)) ;  optfil=0 !To read \tilde\epsilon^-1
   unt1=Dtfil1%unscr ; fname1=Dtfil1%filscr
   if (optfil==1) then 
    !WARNING filchi0 and unchi0 are not initialized, this part has to be rewritten with methods.
    unt1=Dtfil1%unchi0 ; fname1=Dtfil1%filchi0
   end if
   call testscr(optfil,unt1,fname1,nqibz1,nqlwl1,nomega1,npweps1,npwwfn1,nbnds1,titem1,&
&   fform1,MPI_enreg,localrdwf,qibz_p,qlwl_p,omega_p,gvec_p,Hdr_scr)
   if (associated(qlwl_p)) deallocate(qlwl_p)
   deallocate(qibz_p,omega_p,gvec_p)
   call hdr_clean(hdr_scr)

   allocate(omega1(nomega1),qibz1(3,nqibz1))
   allocate(epsm1(npweps1,npweps1,nomega1,nqibz1),stat=istat)
   if (istat/=0) then
    call memerr(FILE__,'epsm1',npweps1**2*nomega1*nqibz1,'spc')
   end if

   optfil=0 !To read \tilde\epsilon^-1
   allocate(lwing_dum(npweps1,nomega1,1*optfil))
   allocate(uwing_dum(npweps1,nomega1,1*optfil))
   allocate(qlwl(3,1*optfil))
   unt1=Dtfil1%unscr ; fname1=Dtfil1%filscr
   if (optfil==1) then 
    unt1=Dtfil1%unchi0 ; fname1=Dtfil1%filchi0
   end if
   call rdscr(optfil,unt1,fname1,npweps1,nqibz1,nqibz1,nomega1,qibz1,omega1,gmet,&
&   epsm1,MPI_enreg,localrdwf,0,.TRUE.,qlwl,uwing_dum,lwing_dum)
   deallocate(lwing_dum,uwing_dum,qlwl)

   do iq=1,nqibz1
    if (repetition(ifile,iq)) cycle
    iqW=iqW+1

    Dtfil_out%unscr=777
    unt_out=Dtfil_out%unscr
    fname_out=TRIM(Dtfil_out%filnam_ds(4))//'_SCR'
!   shoul fix the problem with the name,

    optfil=0 !To write \tilde\epsilon^-1
    allocate(lwing_dum(npweps0,nomega0,optfil),uwing_dum(npweps0,nomega0,optfil),qsmall_dum(3,optfil))
    call wrscr(iqW,optfil,unt_out,fname_out,Hdr0,Dtset,npweps0,npwwfn0,nbnds0,&
&    nq_merge,0,nomega0,q_merge,omega1,gvec0,gmet,epsm1(:,:,:,iq),titem0,&
&    qsmall_dum,lwing_dum,uwing_dum)
    deallocate(lwing_dum,uwing_dum,qsmall_dum)
   end do

   deallocate(epsm1,omega1,qibz1)
  end do

  write(msg,'(3a)')ch10,' Merging files finished successfully',ch10
  call wrtout(std_out,msg,'COLL')
 else
  write(msg,'(3a)')ch10,' Checking file finished successfully',ch10
  call wrtout(std_out,msg,'COLL')
 end if

 call hdr_clean(Hdr0)

 deallocate(filenames)
 deallocate(foundq,repetition,q_merge,merge_table)
 deallocate(gvec0,omega0d,symrec,symafm)

 end program mrgscr
!!***
