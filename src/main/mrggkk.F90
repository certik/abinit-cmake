!{\src2tex{textfont=tt}}
!!****p* ABINIT/mrggkk
!! NAME
!! mrggkk
!!
!! FUNCTION
!! This program merges a GS file and several 1WF or GKK files for
!! different qvectors and perturbations.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer, MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! NOTES
!! GKK file structure is composed of header records and eigenvalue arrays,
!! in binary or ascii:
!!   GS header = hdr
!!   GS eigenvalues = eigen
!!   number of perturbations = ntot
!!   for each perturbation
!!      1WF header = hdr1
!!      1st order eigenvalues = eigen1
!!
!! PARENTS
!!
!! CHILDREN
!!      hdr_clean,hdr_io,herald,leave_new,rwwf,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

program mrggkk

 use defs_basis
 use defs_datatypes
 use defs_infos


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13io_mpi
 use interfaces_14iowfdenpot
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!Local variables-------------------------------
 ! Only here for call to rwwf
!scalars
 integer :: binascii,fform,formeig,headform,i1wf,icg,igkk,ikpt,ios,isppol,mband
 integer :: mpw,n1wf,nband_disk,ngkk,nspinor,ntot,ntotgkk,option,optkg=0,rdwr
 integer :: rdwrout,tim_rwwf=0,unit1wf=22,unitgkk=24,unitgs=21,unitout=23
 character(len=24) :: codename
 character(len=500) :: message
 character(len=fnlen) :: file1wf,filegkk,filegs,outfile
 type(MPI_type) :: mpi_enreg
 type(hdr_type) :: hdr,hdr1
 type(wffile_type) :: wff_dum
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp),allocatable :: cg(:,:),eigen(:),occ(:)

! *************************************************************************

!dummy wff datatype to feed to rwwf
 wff_dum%unwff = unitgs
 wff_dum%formwff = 0  !scalar eigenvalues
 wff_dum%accesswff = 0
 wff_dum%kgwff = 0
 wff_dum%fname = ""

 codename='MRGGKK'//repeat(' ',18)
 
!write greating,read the file names, etc.
 call herald(codename,abinit_version,std_out)

 write(message,'(17a)')&
& ' Files file format: ',ch10,ch10,&
& '  Name of the output file',ch10,&
& '  Integer flag: 0 --> binary output,   1 --> ascii formatted output',ch10,&
& '  Name of the groud state wavefunction file WF',ch10,&
& '  Number of 1WF, of GKK files, and number of 1WF files in all the GKK files',ch10,&
& '  Names of the 1WF files...',ch10,&
& '  Names of the GKK files...',ch10,ch10,&
& ' Enter name of output file: '
 call wrtout(6,message,'COLL')

!get file with filenames and number of 1wf files
 read(*,*) outfile
 read(*,*) binascii

 read(*,*) filegs
 read(*,*) n1wf,ngkk,ntotgkk

!write(*,*)'echo file names out/in/#1wf/#gkk/#totgkk: ',&
!&           trim(outfile)," ",trim(filegs),n1wf,ngkk,ntotgkk
 write(message,'(7a,i4,2a,i4,2a,i4,a)')&
& ' Output                     = ',trim(outfile),ch10,&
& ' Ground State file          = ',trim(filegs),ch10,&
& ' Number of 1WF files        = ',n1wf,ch10,&
& ' Number of GKK files        = ',ngkk,ch10,&
& ' Total Number of 1WF in GKK = ',ntotgkk,ch10
 call wrtout(6,message,'COLL')


!output without rewinding the file
 if (binascii == 0) then
! open output file
  open(unit=unitout,file=outfile,form='unformatted',iostat=ios)
  rdwrout = 6
 else if (binascii == 1) then
! rdwrout=4 ! use for screen output and change writes of eigen to (*,*)
! MJV 27/5/2008 removed 'new' constraint on gkk files: presume competent user!
  open(unit=unitout,file=outfile,form='formatted',iostat=ios)
  rdwrout = 4
 else
  write(message,'(4a)')ch10,&
&  ' mrggkk : ERROR- ,',ch10,&
&  ' binascii must be 0 or 1'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 if (ios/=0) then
  write(message,'(5a)')&
&  ' mrggkk: ERROR- ',ch10,&
&  ' opening file: ',trim(outfile),' as new'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 rewind (unitout)

!-------------------------------------------------------
!now read and write information for GS file
!-------------------------------------------------------
 
!open GS wf file
 write(message,'(a)')' normal input for GS file'
 call wrtout(6,message,'COLL')

 open(unit=unitgs,file=filegs,form='unformatted',status='old',iostat=ios)
 if (ios/=0) then
  write(message,'(5a)')&
&  ' mrggkk: ERROR- ',ch10,&
&  ' opening file: ',trim(filegs),' as old'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 rewind (unitgs)
!read header of GS wf file
 rdwr = 5
 call hdr_io(fform,hdr,rdwr,unitgs)
 if (fform == 0) then
  write(message,'(5a)')ch10,&
&  ' mrggkk: ERROR- ',ch10,&
&  ' reading header in ',trim(filegs)
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!copy header of GS file to output
 call hdr_io(fform,hdr,rdwrout,unitout)

 write(message,'(a)')' header echoed to output file'
 call wrtout(6,message,'COLL')

!retrieve GS eigenvalues from GS wf file and echo to output
 mband = maxval(hdr%nband)
 mpw = maxval(hdr%npwarr)
 allocate(cg(2,mpw*hdr%nspinor*mband))
 allocate(eigen(mband))
 allocate(kg_k(3,0))
 allocate(occ(mband))
 option = 1
 formeig = 0
 icg = 0
 headform=hdr%headform
 wff_dum%unwff = unitgs
 wff_dum%formwff = formeig !scalar eigenvalues
 do isppol=1,hdr%nsppol
  do ikpt=1,hdr%nkpt
   call rwwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,&
&   mband,mpw*hdr%nspinor*mband,mpi_enreg,hdr%nband(ikpt),nband_disk,&
&   hdr%npwarr(ikpt),hdr%nspinor,occ,option,optkg,tim_rwwf,wff_dum)
   if (binascii==0) then
    write(unitout) eigen(1:hdr%nband(ikpt))
   else
    write(unitout,*) eigen(1:hdr%nband(ikpt))
   end if
  end do
 end do

 deallocate(cg,eigen,kg_k,occ)

!close GS wf file
 close (unitgs)
 call hdr_clean(hdr)

 ntot = n1wf + ntotgkk
 if (binascii==0) then
  write (unitout) ntot
 else
  write (unitout,*) ntot
 end if

!-------------------------------------------------------
!now read and write information for 1WF files
!-------------------------------------------------------

 formeig = 1
 wff_dum%unwff = unit1wf
 wff_dum%formwff = formeig !scalar eigenvalues

 do i1wf=1,n1wf
! for each 1wf file, get name...
  read(*,*) file1wf

! open 1wf file
  write(message,'(a)')' normal input for 1WF file '
  call wrtout(6,message,'COLL')

  open(unit=unit1wf,file=file1wf,form='unformatted',status='old',iostat=ios)
  if (ios/=0) then
   write(message,'(5a)')&
&   ' mrggkk: ERROR- ',ch10,&
&   ' opening file: ',trim(file1wf),' as new'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

  rewind (unit1wf)

! read in header of _WF1 file
  rdwr = 5
  call hdr_io(fform,hdr1,rdwr,unit1wf)
  if (fform == 0) then
   write(message,'(4a,i4,2a)')ch10,&
&   ' mrggkk : ERROR- ',ch10,&
&   ' 1WF header number ',i1wf,ch10,&
&   ' was mis-read. fform == 0'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! copy header of 1WF file to output
! WARNING: cant use normal hdr_io because it rewinds the output file and
! destroys previous content.
  call hdr_io(fform,hdr1,rdwrout,unitout)

! retrieve 1WF <psi_k+q | H | psi_k> from 1wf file and echo to output
  mband = maxval(hdr1%nband)
  mpw = maxval(hdr1%npwarr)
  allocate(cg(2,mpw*hdr1%nspinor*mband))
  allocate(eigen(2*mband*mband))
  allocate(kg_k(3,0))
  allocate(occ(mband))
  option = 1
  headform=hdr1%headform
  do isppol=1,hdr1%nsppol
   do ikpt=1,hdr1%nkpt
!   write (*,*) 'isppol,ikpt = ', isppol,ikpt
    call rwwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,&
&    mband,mpw*hdr1%nspinor*mband,mpi_enreg,hdr1%nband(ikpt),nband_disk,&
&    hdr1%npwarr(ikpt),hdr1%nspinor,occ,option,optkg,tim_rwwf,wff_dum)
    if (binascii==0) then
     write(unitout) eigen(1:2*hdr1%nband(ikpt)**2)
    else
     write(unitout,*) eigen(1:2*hdr1%nband(ikpt)**2)
    end if
   end do
  end do
  
  deallocate (cg,eigen,kg_k,occ)

  close (unit1wf)
! clean header to deallocate everything
  call hdr_clean(hdr1)

 end do

!-------------------------------------------------------
!now read and write information for small GKK files
!-------------------------------------------------------
 formeig = 1
 do igkk=1,ngkk
! for each gkk file, get name...
  read(*,*) filegkk

! open gkk file
  write(message,'(a)')' normal input for GKK file'
  call wrtout(6,message,'COLL')

  open(unit=unitgkk,file=filegkk,form='unformatted',status='old',iostat=ios)
  if (ios/=0) then
   write(message,'(5a)')&
&   ' mrggkk: ERROR- ',ch10,&
&   ' opening file: ',trim(filegkk),' as old'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

  rewind (unitgkk)

! read in header of GS file and eigenvalues
  call hdr_io(fform,hdr,5,unitgkk)
! 
! could force a comparison of header with global header above for consistency
! 

  allocate(eigen(mband))
  write(message,'(a)')'mrggkk : try to reread GS eigenvalues'
  call wrtout(6,message,'COLL')

  do isppol=1,hdr%nsppol
   do ikpt=1,hdr%nkpt
    read (unitgkk) eigen(1:hdr%nband(ikpt))
   end do
  end do
  read(unitgkk) n1wf
  deallocate (eigen)

  allocate(eigen(2*mband*mband))
  do i1wf=1,n1wf
!  read in header of 1WF file
   rdwr = 5
   call hdr_io(fform,hdr1,rdwr,unitgkk)
   if (fform == 0) then
    write(message,'(4a,i4,a)')ch10,&
&    ' mrggkk : ERROR- ',ch10,&
&    ' 1WF header number ',i1wf,' was mis-read. fform == 0'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

!  copy header of 1WF file to output
!  WARNING: cant use normal hdr_io because it rewinds the output file and
!  destroys previous content.
   call hdr_io(fform,hdr1,rdwrout,unitout)

!  retrieve 1WF <psi_k+q | H | psi_k> from gkk file and echo to output
   mband = maxval(hdr1%nband)
   option = 1
   do isppol=1,hdr1%nsppol
    do ikpt=1,hdr1%nkpt
!    write (*,*) 'isppol,ikpt = ', isppol,ikpt
     read (unitgkk) eigen(1:2*hdr1%nband(ikpt)**2)

     if (binascii==0) then
      write (unitout) eigen(1:2*hdr1%nband(ikpt)**2)
     else
      write (unitout,*) eigen(1:2*hdr1%nband(ikpt)**2)
     end if
    end do
   end do
   call hdr_clean(hdr1)
  end do
! end loop over 1wf segments in small gkk file
  deallocate (eigen)

  close (unitgkk)
  call hdr_clean(hdr)

 end do
!end loop over small gkk files

 close (unitout)

 write(message,'(2a)')ch10,' Done'
 call wrtout(6,message,'COLL')

!DEBUG
!allocate (eigen(2*mband*mband))
!write(*,*) 'mrggkk : try to reread the info we just wrote to file'
!open (unit=unitout,file=outfile,form='unformatted',status='old')
!write(*,*) 'mrggkk : try to reread GS header'
!call hdr_io(fform,hdr,5,unitout)
!write(*,*) 'mrggkk : try to reread GS eigenvalues'
!do isppol=1,hdr%nsppol
!do ikpt=1,hdr%nkpt
!read(unitout) eigen(1:hdr%nband(ikpt))
!end do
!end do
!write(*,*) 'mrggkk : try to reread number of compiled 1WF files'
!read(unitout) n1wf
!write(*,*) 'mrggkk : n1wf = ', n1wf
!do i1wf=1,n1wf
!write(*,*) 'mrggkk : try to reread 1st order header '
!call hdr_io(fform,hdr1,5,unitout)
!do isppol=1,hdr1%nsppol
!do ikpt=1,hdr1%nkpt
!write(*,*) 'mrggkk : try to reread 1st order eigenvalues '
!read(unitout) eigen(1:2*hdr1%nband(1)**2)
!end do
!end do
!write(*,*)' eigen(1:16) = ', eigen(1:16)
!end do
!close(unitout)
!write(*,*)'It worked ! Re-reading of gkk file successful'
!deallocate(eigen)
!ENDDEBUG
 end program mrggkk
!!***
