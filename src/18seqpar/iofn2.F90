!{\src2tex{textfont=tt}}
!!****f* ABINIT/iofn2
!! NAME
!! iofn2
!!
!! FUNCTION
!! First, read and echo pseudopotential filenames from unit 05.
!! Then read the pseudopotential header of each psp file,
!!  in order to initialize pspheads(1:npsp).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, FrD, AF, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  npsp=number of pseudopotentials
!!
!! OUTPUT
!!  pspheads(npsp)=<type pspheader_type>=all the important information from the
!!   pseudopotential file headers, as well as the psp file names
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      close_xmlfile,leave_new,leave_test,mpi_bcast,mpi_comm_rank,open_xmlfile
!!      psxml2ab,timab,wrtout,xml_parse
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!BEGIN TF_CHANGES
subroutine iofn2(npsp,pspheads,mpi_enreg)
!END TF_CHANGES

 use defs_basis
 use defs_datatypes
 use m_pseudo_types
 use m_pseudo
#if defined HAVE_XMLF90
 use flib_sax
#endif
#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_13psp
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: npsp
 type(pspheader_type),intent(out) :: pspheads(npsp)
 type(MPI_type), intent(inout) :: mpi_enreg

!Local variables-------------------------------
!In case a xc core correction is to be taken into account,
!the n1xccc value will be given by n1xccc_default. Otherwise it is set to 0.
 integer,parameter :: n1xccc_default=2501
 integer :: idum,ii,ilmax,ios,ipsang,ipsp,lang,lmax,mmax,mpsang,n1xccc,nmesh,pspcod
 integer :: pspso,pspxc,test_paw
 integer :: usexml
 real(dp) :: al,e990,e999,fchrg,hdum(3),qchrg,r1,rchrg,rp,rr,rs
 character(len=fnlen) :: filpsp
 character(len=3) :: testxc
 character(len=500) :: message
 character(len=80) :: pspline
 integer :: nproj(0:3),nprojso(1:3)
 integer,allocatable :: orb(:)
 character(len=70)         :: testxml

!no_abirules
#if defined HAVE_XMLF90
 integer                   :: iostat
 type(xml_t)               :: fxml
 type(pseudo_t), pointer   :: psxml
#endif

#if defined MPI 
          !Variables introduced for MPI version
           integer :: ierr,me,nproc
           integer,allocatable :: list_int(:)
           real(dp) :: tsec(2)
           real(dp),allocatable :: list_dpr(:)
           character(len=fnlen),allocatable :: list_char(:)
           character(len=6) :: tag
#endif

!*************************************************************************

 test_paw=0

#if defined MPI
          !Determine who I am
           call MPI_COMM_RANK(MPI_COMM_WORLD,me,ierr)

           if(me==0) then
#endif

 do ipsp=1,npsp

! Read the name of the psp file
  write(6, '(/,a)' ) &
&   ' iofn2 : Please give name of formatted atomic psp file'
  read (5, '(a)' , iostat=ios ) filpsp
  pspheads(ipsp)%filpsp=trim(filpsp)
! It might be that a file name is missing
  if(ios/=0)then
   write(message, '(a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&   ' iofn2 : ERROR -',ch10,&
&   '  There are not enough names of pseudopotentials',ch10,&
&   '  provided in the files file.',ch10,&
&   '  Action : check first the variable ntypat (and/or npsp) in the input file;',ch10,&
&   '  if they are correct, complete your files file.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

  write(6, '(a,i4,a,a)' ) &
&   ' iofn2 : for atom type',ipsp,' , psp file is ',trim(filpsp)

! Check if the file is written in XML
  usexml = 0
  open (unit=tmp_unit,file=filpsp,form='formatted',status='old')
  rewind (unit=tmp_unit)
    read(tmp_unit,*) testxml
    if(testxml(1:5)=='<?xml')then
      usexml = 1
    else
      usexml = 0
    end if
  close (unit=tmp_unit)

! Read the header of the pseudopotential file
  if( usexml /= 1) then
! Open the psp file
    open (unit=tmp_unit,file=filpsp,form='formatted',status='old')
    rewind (unit=tmp_unit)

! Read the three first lines
    read (tmp_unit, '(a)' )pspheads(ipsp)%title
    read (tmp_unit,*)pspheads(ipsp)%znuclpsp,pspheads(ipsp)%zionpsp,pspheads(ipsp)%pspdat
    read (tmp_unit,*)pspheads(ipsp)%pspcod,pspheads(ipsp)%pspxc,&
&                    pspheads(ipsp)%lmax,idum,mmax
    pspcod=pspheads(ipsp)%pspcod
    lmax=pspheads(ipsp)%lmax
    write(6, '(a,f5.1,a,i4,a,i4)' ) '  read the values zionpsp=',&
&    pspheads(ipsp)%zionpsp,' , pspcod=',pspcod,' , lmax=',lmax

    nproj(0:3)=0 ; nprojso(1:3)=0

    pspheads(ipsp)%xccc=0
    pspheads(ipsp)%pspso=0

  else if( usexml == 1) then

#if defined HAVE_XMLF90
   if(usexml==1)then
    write(message,'(a,a)')  &
&     ' iofn2 : Reading pseudopotential header in XML form from ', trim(filpsp)
    call wrtout(ab_out,message,'COLL')
    call wrtout(06,  message,'COLL')

    call open_xmlfile(filpsp,fxml,iostat)
    if (iostat /=0) stop "Cannot open file"

    call xml_parse(fxml,begin_element,end_element,pcdata_chunk,verbose=.false.)

    psxml => pseudo

    call psxml2ab( psxml,                       &
&                  pspheads(ipsp)%znuclpsp,     &
&                  pspheads(ipsp)%zionpsp,      &
&                  pspheads(ipsp)%pspcod,       &
&                  pspheads(ipsp)%pspxc,        &
&                  pspheads(ipsp)%lmax, 0 )

    pspcod = pspheads(ipsp)%pspcod
    lmax   = pspheads(ipsp)%lmax

    nproj(0:3) = 0 ; nprojso(1:3) = 0
    if( psxml%header%core_corrections .eq. "yes") then
      pspheads(ipsp)%xccc  = n1xccc_default
    else
      pspheads(ipsp)%xccc  = 0
    endif
    pspheads(ipsp)%pspso = 0

    call close_xmlfile(fxml)
   else
#endif
    write(6, '(4a)' ) ch10,&
&     ' iofn2 : In order to use XML pseudopotentials, you need to compile ABINIT',ch10,&
&     '  with the -DXMLF90 preprocessing option, and also to compile the XMLf90 library. Stop.'
    stop
#if defined HAVE_XMLF90
   end if
#endif

  end if

! DEBUG
!  write(6,*) pspheads(ipsp)%znuclpsp
!  write(6,*) pspheads(ipsp)%zionpsp
!  write(6,*) pspheads(ipsp)%pspcod
!  write(6,*) pspheads(ipsp)%pspxc
!  write(6,*) pspheads(ipsp)%lmax
!  stop
! ENDDEBUG

! Initialize nproj, nprojso, pspso, as well as xccc, for each type of psp
  pspheads(ipsp)%GTHradii = zero

  if(pspcod==1 .or. pspcod==4)then

!  Teter format
   do ilmax=0,lmax
    read (tmp_unit,*) lang,e990,e999,nproj(ilmax)
    read (tmp_unit,*)
   end do
   read (tmp_unit,*) rchrg,fchrg,qchrg
   if (fchrg>1.d-15) pspheads(ipsp)%xccc=n1xccc_default

  else if(pspcod==2)then

!  GTH pseudopotentials
   read (tmp_unit,*) pspheads(ipsp)%GTHradii(0) !rloc
   read (tmp_unit,*) pspheads(ipsp)%GTHradii(1),hdum(1),hdum(2)
   if(abs(hdum(1))>1.d-9) nproj(0)=1
   if(abs(hdum(2))>1.d-9) nproj(0)=2
   read (tmp_unit,*) pspheads(ipsp)%GTHradii(2),hdum(3)
   if(abs(hdum(3))>1.d-9) nproj(1)=1

  else if(pspcod==3)then

!  HGH pseudopotentials
   read (tmp_unit,*) pspheads(ipsp)%GTHradii(0) !rloc
   do ilmax=0,lmax
    read (tmp_unit,*) pspheads(ipsp)%GTHradii(ilmax + 1),hdum(1),hdum(2),hdum(3)
    if (abs(hdum(1))>1.d-9)nproj(ilmax)=1
    if (abs(hdum(2))>1.d-9)nproj(ilmax)=2
    if (abs(hdum(3))>1.d-9)nproj(ilmax)=3
    if (ilmax>0.and.ilmax<3) then
     read (tmp_unit,*) hdum(1),hdum(2),hdum(3)
     if (abs(hdum(1))>1.d-9)nprojso(ilmax)=1
     if (abs(hdum(2))>1.d-9)nprojso(ilmax)=2
     if (abs(hdum(3))>1.d-9)nprojso(ilmax)=3
     if(nprojso(ilmax)>0)pspheads(ipsp)%pspso=2
    end if
    if (ilmax==3) then
     read (tmp_unit,*) hdum(1)
     if (abs(hdum(1))>1.d-9)nprojso(3)=1
     if(nprojso(3)>0)pspheads(ipsp)%pspso=2
    end if
   end do

  else if(pspcod==5)then

!  PHONEY pseudopotentials
!  read parameter for Hamman grid
   pspso=1
   read (tmp_unit,fmt=*,err=10,end=10) r1,al,pspso
10 continue
   do ilmax=0,lmax
    read (tmp_unit,*) lang,e990,e999,nproj(ilmax)
    read (tmp_unit,*)
    if (ilmax>0.and.pspso/=1) then
     read (tmp_unit,*) lang,e990,e999,nprojso(ilmax)
     read (tmp_unit,*)
     pspheads(ipsp)%pspso=pspso
!    Meaning of pspso internally to ABINIT has been changed in v5.4
!    So : file must contain pspso 1 , but ABINIT will have pspso 0 .
     if(pspso==1)pspheads(ipsp)%pspso=0
    end if
   end do
   read (tmp_unit,*) rchrg,fchrg,qchrg
   if (fchrg>1.d-15) pspheads(ipsp)%xccc=n1xccc_default

  else if(pspcod==6)then

!  FHI pseudopotentials
   read (tmp_unit, '(a3)') testxc
!  Note : prior to version 2.2, this 4th line started with  4--  ,
!  and no core-correction was available.
   if(testxc/='4--')then
    backspace(tmp_unit)
    read (tmp_unit,*) rchrg,fchrg,qchrg
   else
    fchrg=0.0_dp
   end if
   if (fchrg>1.d-15) pspheads(ipsp)%xccc=n1xccc_default
!  XG020728 : Should take lloc into account ??
   do ilmax=0,lmax
    nproj(ilmax)=1
   end do

  else if(pspcod==7)then

!  PAW pseudopotentials
   test_paw=1;pspheads(ipsp)%pawheader%pawver=1
   read (tmp_unit,'(a80)') pspline;pspline=adjustl(pspline)
   if (pspline(1:3)=="paw".or.pspline(1:3)=="PAW") &
&   read(unit=pspline(4:80),fmt=*) pspheads(ipsp)%pawheader%pawver
   if (pspheads(ipsp)%pawheader%pawver==1) then   ! Compatibility with Abinit v4.2.x
    read (unit=pspline,fmt=*) pspheads(ipsp)%pawheader%basis_size,&
&                         pspheads(ipsp)%pawheader%lmn_size
    allocate(orb(pspheads(ipsp)%pawheader%basis_size));orb(:)=0
    read (tmp_unit,*) (orb(ii), ii=1,pspheads(ipsp)%pawheader%basis_size)
    read (tmp_unit,*);read (tmp_unit,*) pspheads(ipsp)%pawheader%rpaw
    pspheads(ipsp)%pawheader%rshp=pspheads(ipsp)%pawheader%rpaw
    read (tmp_unit,*) pspheads(ipsp)%pawheader%mesh_size
    read (tmp_unit,*) pspheads(ipsp)%pawheader%shape_type
    if (pspheads(ipsp)%pawheader%shape_type==3) pspheads(ipsp)%pawheader%shape_type=-1
   else
    read (tmp_unit,*) pspheads(ipsp)%pawheader%basis_size,&
&                     pspheads(ipsp)%pawheader%lmn_size
    allocate(orb(pspheads(ipsp)%pawheader%basis_size));orb(:)=0
    read (tmp_unit,*) (orb(ii), ii=1,pspheads(ipsp)%pawheader%basis_size)
    pspheads(ipsp)%pawheader%mesh_size=mmax
    read (tmp_unit,*) nmesh
    do ii=1,nmesh;read(tmp_unit,*);end do
    read (tmp_unit,*) pspheads(ipsp)%pawheader%rpaw
    pspheads(ipsp)%pawheader%rshp=pspheads(ipsp)%pawheader%rpaw
    read (tmp_unit,'(a80)') pspline;pspline=adjustl(pspline);print *,pspline
    read(unit=pspline,fmt=*) pspheads(ipsp)%pawheader%shape_type
    if (pspheads(ipsp)%pawheader%pawver==2.and.&
&       pspheads(ipsp)%pawheader%shape_type==3) pspheads(ipsp)%pawheader%shape_type=-1
    if (pspheads(ipsp)%pawheader%pawver>=3.and.pspheads(ipsp)%pawheader%shape_type==-1) then
     rr=zero;read(unit=pspline,fmt=*,err=20,end=20) ii,rr
20   continue
     if (rr>=tol8) pspheads(ipsp)%pawheader%rshp=rr
    end if
   end if
   do ilmax=0,lmax
    do ii=1,pspheads(ipsp)%pawheader%basis_size
     if(orb(ii)==ilmax) nproj(ilmax)=nproj(ilmax)+1
    end do
   end do
   pspheads(ipsp)%pawheader%l_size=2*maxval(orb)+1
   pspheads(ipsp)%xccc=1  ! We suppose apriori that cc is used (but n1xccc is not used in PAW)
   deallocate (orb)

  else if(pspcod==8)then

!  DRH pseudopotentials
   read (tmp_unit,*) rchrg,fchrg,qchrg
   if (fchrg>1.d-15) pspheads(ipsp)%xccc=n1xccc_default
   read(tmp_unit,*) nproj(0:lmax)
   pspso=0
   pspheads(ipsp)%pspso=pspso

  else if(pspcod==9)then

   fchrg = 0.0_dp
   do ilmax = 0, lmax
    nproj(ilmax) = 1
   end do

  else if(pspcod==10)then

!  HGH pseudopotentials, full h/k matrices
   read (tmp_unit,*) pspheads(ipsp)%GTHradii(0) !rloc
   read (tmp_unit,*) idum; if(idum-1/=lmax) stop "in iofn2: nnonloc-1/ = lmax"
   do ilmax=0,lmax
    read (tmp_unit,*) pspheads(ipsp)%GTHradii(ilmax + 1),nproj(ilmax),(hdum(idum),idum=1,nproj(ilmax))
    do idum=2,nproj(ilmax) !skip the rest of h_ij
     read (tmp_unit,*)
    enddo
    if (ilmax==0) cycle
    nprojso(ilmax)=nproj(ilmax)
    if(nprojso(ilmax)>0)then
     pspheads(ipsp)%pspso=2
     do idum=1,nprojso(ilmax) !skip the rest of k_ij
      read (tmp_unit,*)
     enddo
    endif
   end do

  else

   write(message, '(a,a,a,a,i4,a,a,a,a)' ) ch10,&
&   ' iofn2 : ERROR -',ch10,&
&   '  The pseudopotential code (pspcod) read from file is ',pspcod,ch10,&
&   '  This value is not allowed (should be between 1 and 10). ',ch10,&
&   '  Action : use a correct pseudopotential file.'
   call wrtout(6,message,'PERS')
   call leave_new('PERS')

  end if ! pspcod=...

! Store in pspheads
  pspheads(ipsp)%nproj(0:3)=nproj(0:3)
  pspheads(ipsp)%nprojso(1:3)=nprojso(1:3)

  close(tmp_unit)

 end do ! ipsp=1,npsp

!Note that mpsang is the max of 1+lmax, with minimal value 1 (even for local psps, at present)
!mpsang=max(maxval(pspheads(1:npsp)%lmax)+1,1) ! Likely troubles with HP compiler
!n1xccc=maxval(pspheads(1:npsp)%xccc)
 mpsang=1
 n1xccc=pspheads(1)%xccc
 do ii=1,npsp
  mpsang=max(pspheads(ii)%lmax+1,mpsang)
  n1xccc=max(pspheads(ii)%xccc,n1xccc)
 end do

 write(6, '(/,a,i4,a,i4,a)' ) &
& ' iofn2 : deduce mpsang  =',mpsang,', n1xccc  =',n1xccc,'.'

!Test: if one psp is PAW, all must be
 if (test_paw==1) then
  do ipsp=1,npsp
   if (pspheads(ipsp)%pspcod/=7) then
    write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&    ' iofn2 : ERROR -',ch10,&
&    '  One pseudopotential is PAW (pspcod=7) !',ch10,&
&    '  All pseudopotentials must be PAW (this is not the case here) !',ch10,&
&    '  Action : use only PAW pseudopotential files.'
   call wrtout(6,message,'PERS')
   call leave_new('PERS')
   end if
  end do
 end if

#if defined MPI 
!           End the section me==0
            end if
            call timab(48,1,tsec)
            call leave_test(mpi_enreg)

           !Broadcast the characters (file names and titles)
            allocate(list_char(2*npsp))
            list_char(1:npsp)=pspheads(1:npsp)%filpsp
            list_char(npsp+1:2*npsp)=pspheads(1:npsp)%title
            call MPI_BCAST(list_char,2*npsp*fnlen,MPI_CHARACTER,0,&
          &  MPI_COMM_WORLD,ierr)
            pspheads(1:npsp)%filpsp=list_char(1:npsp)
            pspheads(1:npsp)%title=list_char(npsp+1:2*npsp)
            deallocate(list_char)

           !Brodcast the integers
            allocate(list_int(1+13*npsp))
            list_int(1       :  npsp)=pspheads(1:npsp)%nproj(0)
            list_int(1+  npsp:2*npsp)=pspheads(1:npsp)%nproj(1)
            list_int(1+2*npsp:3*npsp)=pspheads(1:npsp)%nproj(2)
            list_int(1+3*npsp:4*npsp)=pspheads(1:npsp)%nproj(3)
            list_int(1+4*npsp:5*npsp)=pspheads(1:npsp)%lmax
            list_int(1+5*npsp:6*npsp)=pspheads(1:npsp)%xccc
            list_int(1+6*npsp:7*npsp)=pspheads(1:npsp)%pspxc
            list_int(1+7*npsp:8*npsp)=pspheads(1:npsp)%pspdat
            list_int(1+8*npsp:9*npsp)=pspheads(1:npsp)%pspcod
            list_int(1+9*npsp:10*npsp)=pspheads(1:npsp)%pspso
            list_int(1+10*npsp:11*npsp)=pspheads(1:npsp)%nprojso(1)
            list_int(1+11*npsp:12*npsp)=pspheads(1:npsp)%nprojso(2)
            list_int(1+12*npsp:13*npsp)=pspheads(1:npsp)%nprojso(3)
            list_int(1+13*npsp)        =test_paw

            call MPI_BCAST(list_int,1+13*npsp,MPI_INTEGER,0,&
          &  MPI_COMM_WORLD,ierr)

            pspheads(1:npsp)%nproj(0) =list_int(1       :  npsp)
            pspheads(1:npsp)%nproj(1) =list_int(1+  npsp:2*npsp)
            pspheads(1:npsp)%nproj(2) =list_int(1+2*npsp:3*npsp)
            pspheads(1:npsp)%nproj(3) =list_int(1+3*npsp:4*npsp)
            pspheads(1:npsp)%lmax     =list_int(1+4*npsp:5*npsp)
            pspheads(1:npsp)%xccc     =list_int(1+5*npsp:6*npsp)
            pspheads(1:npsp)%pspxc    =list_int(1+6*npsp:7*npsp)
            pspheads(1:npsp)%pspdat   =list_int(1+7*npsp:8*npsp)
            pspheads(1:npsp)%pspcod   =list_int(1+8*npsp:9*npsp)
            pspheads(1:npsp)%pspso    =list_int(1+9*npsp:10*npsp)
            pspheads(1:npsp)%nprojso(1)=list_int(1+10*npsp:11*npsp)
            pspheads(1:npsp)%nprojso(2)=list_int(1+11*npsp:12*npsp)
            pspheads(1:npsp)%nprojso(3)=list_int(1+12*npsp:13*npsp)
            test_paw                   =list_int(1+13*npsp)
            deallocate(list_int)

           !Broadcast zionpsp and znuclpsp
            allocate(list_dpr(2*npsp))
            list_dpr(1:npsp)=pspheads(1:npsp)%zionpsp
            list_dpr(1+npsp:2*npsp)=pspheads(1:npsp)%znuclpsp
            call MPI_BCAST(list_dpr,2*npsp,MPI_DOUBLE_PRECISION,0,&
          &  MPI_COMM_WORLD,ierr)
            pspheads(1:npsp)%zionpsp =list_dpr(1:npsp)
            pspheads(1:npsp)%znuclpsp=list_dpr(1+npsp:2*npsp)
            deallocate(list_dpr)

           !Broadcast additional integers for PAW psps (testpaw was spread, previously)
            if (test_paw==1) then
             allocate(list_int(5*npsp))
             list_int(1       :  npsp)=pspheads(1:npsp)%pawheader%basis_size
             list_int(1+  npsp:2*npsp)=pspheads(1:npsp)%pawheader%l_size
             list_int(1+2*npsp:3*npsp)=pspheads(1:npsp)%pawheader%lmn_size
             list_int(1+3*npsp:4*npsp)=pspheads(1:npsp)%pawheader%mesh_size
             list_int(1+4*npsp:5*npsp)=pspheads(1:npsp)%pawheader%pawver
             call MPI_BCAST(list_int,5*npsp,MPI_INTEGER,0,&
           &  MPI_COMM_WORLD,ierr)
             pspheads(1:npsp)%pawheader%basis_size=list_int(1       :  npsp)
             pspheads(1:npsp)%pawheader%l_size    =list_int(1+  npsp:2*npsp)
             pspheads(1:npsp)%pawheader%lmn_size  =list_int(1+2*npsp:3*npsp)
             pspheads(1:npsp)%pawheader%mesh_size =list_int(1+3*npsp:4*npsp)
             pspheads(1:npsp)%pawheader%pawver    =list_int(1+4*npsp:5*npsp)
             deallocate(list_int)
            end if

            call timab(48,2,tsec)
#endif

end subroutine iofn2
!!***
