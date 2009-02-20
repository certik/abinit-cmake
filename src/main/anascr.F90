!{\src2tex{textfont=tt}}
!!****p* ABINIT/anascr
!! NAME
!! anascr
!!
!! FUNCTION
!!  Tool for analysis the inverse dielectric matrix stored in 
!!  the SCR file. It can be used to:
!!  a) Analyse selected columns of $\epsilon^{-1}_{G Gp}(q,\w)$ for a fixed 
!!     q point and $\w$
!!  b) Analyse the frequency dependency of $\epsilon^{-1}_{G Gp}(q,\w)$ for a
!!     fixed q point and (G,Gp) pair 
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Main program 
!!
!! OUTPUT
!!  Write on an external file  
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      hdr_clean,hdr_io,hdr_io_netcdf,herald,int2char4,leave_new,mati3inv
!!      memerr,metric,mkppm,prompt,rdppm,testppm,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

program anascr

 use defs_basis
 use defs_datatypes
 use defs_infos
 use m_IO_tools, only : prompt


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_14iowfdenpot
 use interfaces_15gw
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!Local variables-------------------------------
!scalars
 integer,parameter :: unitem1=111,unitppm=112,unt=211
 integer :: choice,fform,ic,ig,ig1a,ig2a,igp,ii,index,io,ios,ip,iprompt,iq,iqa
 integer :: istat,isym,itask,iwa,nbnds,ncol,ngpi,nomega,npwc,npwc2,npwc3,npwvec
 integer :: npwwfn,nq,nqa,nsym,ppmodel,rdwr,readnetcdf,step
 real(dp) :: domegareal,omegaermax,ucvol
 logical :: filexist,verbose=.true.
 character(len=24) :: codename
 character(len=4) :: tagq,tagw
 character(len=500) :: msg
 character(len=fnlen) :: filem1,filppm,fname,fname_default
 type(hdr_type) :: hdr
!arrays
 integer,allocatable :: gptab(:,:),gptabo(:,:),gvec(:,:),ip2fp(:,:)
 integer,allocatable :: symrec(:,:,:),symrel(:,:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: op(:,:,:),q(:,:),ww(:)
 complex(dpc),allocatable :: epsm1(:,:),epsm1w(:,:,:),omega(:)
 complex,allocatable :: bigomegatwsq(:,:,:),eigvec(:,:,:),em1w(:)
 complex,allocatable :: omegatw(:,:,:)
 character(len=80) :: titem1(2)

! *************************************************************************

 codename='ANASCR'//REPEAT(' ',18)
 call herald(codename,abinit_version,std_out)

 call prompt('Enter your choice (1 for SCR file, 2 for PPM files)',choice)

 if (choice==1) then 
  do   
   call prompt('Enter name of screening file:',filem1)
   filem1=TRIM(filem1)

!  Checking the existence of data file
   inquire(file=filem1,exist=filexist)
   if (.not.filexist) then
    write(msg,'(4a)')ch10,&
&    ' ERROR-  missing data file: ',trim(filem1),ch10
    call wrtout(std_out,msg,'coll') 
    call leave_new('COLL')
   end if

   open(file=filem1,unit=unitem1,status='old',form='unformatted',iostat=ios)
   if (ios/=0) then 
    write(msg,'(4a)')ch10,&
&    ' ERROR- opening file ',trim(filem1),ch10
    call wrtout(std_out,msg,'COLL') 
    call leave_new('COLL')
   end if
   exit
  end do

! Read the header
  rdwr= 1
  readnetcdf= 0 !should make EM1 a netcdf file as well
  if (readnetcdf==0) then
   call hdr_io(fform,hdr,rdwr,unitem1)
  end if
  if (readnetcdf==1) then
   call hdr_io_netcdf(fform,hdr,rdwr,unitem1)
  end if

  if (fform/=1002) then
   write(msg,'(3a)')ch10,&
&   ' ERROR- Unknown file format found ',ch10
   call wrtout(6,msg,'COLL')
   call leave_new('COLL')
  end if 

! Echo part of the header
  rdwr=4 
  call hdr_io(fform,hdr,rdwr,6)

  read(unitem1) titem1
  write(msg,'(2a,1x,5a)') ch10,ch10,titem1(1),ch10,titem1(2),ch10,ch10
  call wrtout(6,msg,'COLL')

  read(unitem1) npwvec,npwwfn,nbnds,nq,nomega

  write(msg,'(a,i8)')' Number of G vectors in epsilon^-1     ',npwvec
  call wrtout(6,msg,'COLL')
  write(msg,'(a,i8)')' Number of G vectors for wavefunctions ',npwwfn
  call wrtout(6,msg,'COLL')
  write(msg,'(a,i8)')' Number of bands                       ',nbnds
  call wrtout(6,msg,'COLL')
  write(msg,'(a,i8)')' Number of q-points                    ',nq
  call wrtout(6,msg,'COLL')
  write(msg,'(a,i8,2a)')' Number of evaluated frequencies       ',nomega,ch10,ch10
  call wrtout(6,msg,'COLL')
  
  allocate(gvec(3,npwvec))
  read(unitem1) gvec(1:3,1:npwvec)

  write(msg,'(3a)')ch10,' G vectors in epsilon^-1 and corresponding indexes:',ch10
  call wrtout(6,msg,'COLL')
  do ig=1,npwvec
   write(*,'(3(1x,i5,a,3x,a,3(i5,2x),a))',advance='no')ig,'*','(',gvec(:,ig),')'
   if (mod(ig,3)==0) write(*,*) 
  end do

  allocate(q(3,nq))
  read(unitem1) q(1:3,1:nq)

  allocate(omega(1:nomega))
  read(unitem1) omega(1:nomega)

  write(*,'(/,a,/,(i3,3f14.6))')&
&  ' q-points [reciprocal lattice units]:',(iq,(q(ii,iq),ii=1,3),iq=1,nq)
! call wrtout(6,msg,'COLL')
  write(*,'(/,a,/,(i3,2f7.2))')&
&  ' frequencies omega [eV]:',(io,omega(io)*Ha_eV,io=1,nomega)
! call wrtout(6,msg,'COLL')

  call hdr_clean(hdr)
! ------------------------------------------------------------------------

  do 

   write(msg,'(2a)')ch10,' What is your choice? Type:'
   call wrtout(6,msg,'COLL')
   write(msg,'(a)')'  0 => exit'
   call wrtout(6,msg,'COLL')
   write(msg,'(a)')'  1 => q-point    (analyze epsilon^-1 in a single q-point)'
   call wrtout(6,msg,'COLL')
   write(msg,'(a)')'  2 => frequency  (analyze frequency dependency of epsilon^-1)'
   call wrtout(6,msg,'COLL')

   read(*,*) itask
   select case(itask)

    case(1) !single q-point analysis

!    call skip_hdr_scr(unitem1,readnetcdf) 

!    Re-Read the header of the file
     rewind(unitem1)
     rdwr=1
!    should make EM1 a netcdf file as well
     call hdr_io(fform,hdr,rdwr,unitem1)
     call hdr_clean(hdr)
     
     read(unitem1) !skip titem1
     read(unitem1) npwvec,npwwfn,nbnds,nq,nomega
     read(unitem1) !skip gvec
     read(unitem1) !skip q
     read(unitem1) !skip omega

     do 
      call prompt('Enter index of the required q-point',iqa)
      if (iqa> nq.or.iqa<=0) then 
       write(msg,'(2a)')ch10,' Wrong value for q-point index'
       call wrtout(6,msg,'COLL')
       CYCLE
      else 
       exit
      end if 
     end do 

     do 
      call prompt(' Enter index of the required frequency',iwa)
      if (iwa > nomega .or. iwa<=0) then
       write(msg,'(2a)')ch10,' Wrong values for frequency index'
       call wrtout(6,msg,'COLL')
       CYCLE
      else 
       exit
      end if 
     end do 

     call int2char4(iqa,tagq)
     call int2char4(iwa,tagw)
     fname_default=trim(filem1)//'_em1_Q'//tagq//'_w'//tagw
     fname=trim(fname_default)

     do 
      inquire(file=fname,exist=filexist)
      if (filexist) then
       write(msg,'(6a)')ch10,&
&       ' file: ',trim(fname),' already exists ! ',ch10,&
&       ' overwrite? (0=default=no, 1=yes)'
       call wrtout(6,msg,'COLL')
       read(*,*)iprompt 
       if (iprompt==1) exit 
       call prompt('Enter the name of a new output file:',fname)
      else 
       exit 
      end if                                                 
     end do 

     open(unit=unt,file=fname,status='unknown',iostat=ios)
     if (ios/=0) then 
      write(msg,'(4a)')ch10,&
&      ' ERROR- opening file ',trim(fname),ch10
      call wrtout(6,msg,'COLL')
      call leave_new('COLL')
     end if

     do  
      write(msg,'(2a,i8,a)')ch10,' Enter step for G  [ Max = ',npwvec,' ]'
      call wrtout(06,msg,'COLL')
      read(*,*)step
      if(step<=0 .or. step > npwvec) then 
       write(msg,'(2a)')ch10,' Wrong values for G index'
       call wrtout(6,msg,'COLL')
       cycle
      end if 
      ncol=int(npwvec/step)+1
      exit
     end do 
     
     write(msg,'(a,i3,a)')' Writing ',ncol,' columns '
     call wrtout(6,msg,'COLL')
!    write(*,'(1x,i4)')(ic,ic=1,npwvec,step)
     
!    write info on the output file
     write(unt,'(2a,/,a,i5)')&
&     '# epsilon^-1 matrix read from file : ',trim(filem1),&
&     '# Number of G vectors : ',npwvec
     if (ncol>1) then 
      write(unt,'(a)')'# Norm of diagonal elements and selected columns '
      write(unt,'(a)',advance='no')'# Row   Diagonal'  
      do ic=1,npwvec,step
       write(unt,'(i8,2x)',advance='no')ic
      end do
      write(unt,'(/)')
     else
      write(unt,'(a)')'# Norm of diagonal elements'
      write(unt,'(a,/)')'# Row   Diagonal'  
     end if

     allocate(epsm1(npwvec,npwvec),stat=istat)
     if(istat/=0) then
      call memerr('anascr','epsm1',npwvec**2,'dpc')
     end if
     
!    read epsilon^-1 tilde
     do iq=1,nq
      do io=1,nomega
       read(unitem1) epsm1(1:npwvec,1:npwvec)
!      find required q-points and frequency 
       if (iq==iqa .and. io==iwa) then
        do ig=1,npwvec
!        modulus of diagonal elements 
         write(unt,'(1x,i4,2x,f8.5)',advance='no')ig,(abs(epsm1(ig,ig)))
         do ic=1,npwvec,step
          write(unt,'(2x,f8.5)',advance='no')(abs(epsm1(ig,ic)))
         end do
         write(unt,*)
        end do
       end if
      end do
     end do 
     
     deallocate(epsm1)
     close(unt)   

    case(2) 

!    call skip_hdr_scr(unitem1,readnetcdf) 

!    Re-Read the header of the file
     rewind(unitem1)
     rdwr=1
!    should make EM1 a netcdf file as well
     call hdr_io(fform,hdr,rdwr,unitem1)
     nsym= hdr%nsym 

     call metric(gmet,gprimd,std_out,rmet,hdr%rprimd,ucvol)

!    calculate sym op in reciprocal space, no symmo 
     allocate(symrec(3,3,nsym))
     do isym=1,nsym
      call mati3inv(hdr%symrel(:,:,isym),symrec(:,:,isym))
     end do

     
     call hdr_clean(hdr)
     
     read(unitem1) !skip titem1
     read(unitem1) npwvec,npwwfn,nbnds,nq,nomega
     read(unitem1) !skip gvec
     read(unitem1) !skip q
     read(unitem1) !skip omega

     allocate(gptab(npwvec,npwvec),gptabo(npwvec,npwvec),stat=istat)
     if(istat/=0) then
      call memerr('anascr','gptab/gptabo',2*npwvec**2,'i4b')
     end if

!    Calculate the independent G-Gp pairs, TODO gptabo could be i1b
!    used to write the independent matrix elements of epsilon^-1 in outeps 
     allocate(ip2fp(2,npwvec**2)) 
     if(istat/=0) then
      call memerr('anascr','ip2fp',2*npwvec**2,'i4b')
     end if
     
     gptab=  0 
     gptabo= 0
     ip2fp=  0

     allocate(epsm1w(npwvec,npwvec,nomega),stat=istat)
     if(istat/=0) then
      call memerr('anascr','epsm1',nomega*npwvec**2,'dpc')
     end if

     if (verbose) write(*,*)' writing only first and last 50 pairs'
     
     do iq=1,nq

      call int2char4(iq,tagq)
      fname_default= trim(filem1)//'_Q'//tagq
      fname= trim(fname_default)
      open(file=fname,unit=999,status='new',form='formatted',iostat=ios)

!     Here we use the small non zero q to deal with the long wavelength limit
      allocate(op(3,3,nsym))
      op(:,:,:)=symrec(:,:,:)
!     FIXME this is BROKEN
!     call findggp(q(:,iq),nsym,op,npwvec,gvec,gmet,ngpi,gptab,gptabo,ip2fp)
      deallocate(op)

      do io=1,nomega
       read(unitem1) epsm1w(1:npwvec,1:npwvec,io)
      end do

      index=0
      do ip=1,ngpi
       if (verbose .and. (ip>50 .and. ip<=ngpi-50)) cycle
       ig = ip2fp(1,ip)
       igp= ip2fp(2,ip)
       index=index+1
       write(999,'(a,i4,a,i8,/,a,3f12.6,/,a,3i6,a,3i6,/,a,/)')&
&       '# index= ',index,' pair number = ',ip,&
&       '# q = ',q(:,iq),&
&       '# G = ',gvec(:,ig),'  G''= ',gvec(:,igp),&
&       '#   w [eV]           Re             Im ' 
       do io=1,nomega
        if (abs(aimag(omega(io)))>tol8) exit
        write(999,'(f8.2,4x,2es16.8)')real(omega(io))*Ha_eV,epsm1w(ig,igp,io)
       end do
       write(999,*)
       write(999,*)
      end do

      close(999)

     end do

     deallocate(epsm1w)

    case(0)
     write(msg,'(2a)')ch10,' Exit requested by user'
     call wrtout(6,msg,'COLL')
     exit

     case default
     cycle

   end select
  end do

  write(msg,'(5a)')ch10,' Analysis completed',ch10,&
&  ' Thank you for using me',ch10
  call wrtout(6,msg,'COLL')                                                                  

  close(unitem1)
  deallocate(omega)
  stop 

 else if (choice==2) then 

  do   

   call prompt('Enter name of PPM file:',filppm)
   filppm=TRIM(filppm)

!  checking the existence of data file
   inquire(file=filppm,exist=filexist)
   if (.not. filexist) then
    write(msg,'(4a)')ch10,&
&    ' error-  missing data file: ',trim(filppm),ch10
    call wrtout(std_out,msg,'coll')
    CYCLE
   end if
   open(file=filppm,unit=unitppm,status='old',form='unformatted',iostat=ios)
   if (ios/=0) then 
    write(msg,'(4a)')ch10,&
&    ' ERROR- opening file ',trim(filppm),ch10
    call wrtout(std_out,msg,'COLL') 
    call leave_new('COLL')
   end if
   EXIT
  end do

  call testppm(unitppm,ppmodel,npwc,npwc2,npwc3,nq,hdr)
  nsym= hdr%nsym 

  write(*,*)' plasmon-pole model : ',ppmodel

  call metric(gmet,gprimd,std_out,rmet,hdr%rprimd,ucvol)

! allocate(symrel(3,3,nsym))
  allocate(gvec(3,npwc),q(3,nq))
  allocate(bigomegatwsq(npwc,npwc2,nq),omegatw(npwc,npwc3,nq))
! this should be allocated only if ppmodel==3
  allocate(eigvec(npwc,npwc,nq))

  call rdppm(unitppm,ppmodel,npwc,npwc2,npwc3,gvec,nq,q,bigomegatwsq,omegatw,eigvec)

! calculate sym op in reciprocal space, no symmo 
  allocate(symrec(3,3,nsym))
  do isym=1,nsym
   call mati3inv(hdr%symrel(:,:,isym),symrec(:,:,isym))
  end do

! DEBUG
! write(*,*)bigomegatwsq,omegatw
! write(*,*)eigvec
! ENDDEBUG

  call prompt('Enter number of frequencies to evaluate eps^-1',nomega)
  if (nomega<=0) nomega= 20
  allocate (em1w(nomega),ww(nomega))

  call prompt('Enter maximum frequency value on the real axis [Ha]',omegaermax)
  domegareal=omegaermax/(nomega-1)

  do io=1,nomega
   ww(io)= (io-1)*domegareal
  end do

! Allocate tables for G vectors
  allocate(gptab(npwc,npwc),gptabo(npwc,npwc),stat=istat)
  if(istat/=0) then
   call memerr('anascr','gptab/gptabo',2*npwc**2,'i4b')
  end if

! Calculate the independent G-Gp pairs, TODO gptabo could be i1b
! used to write the independent matrix elements of epsilon^-1 in outeps 
  allocate(ip2fp(2,npwc**2)) 
  if(istat/=0) then
   call memerr('anascr','ip2fp',2*npwc**2,'i4b')
  end if
  
  gptab=  0 
  gptabo= 0
  ip2fp=  0

  if (verbose) write(*,*)' writing only first and last 50 pairs'

  do iq=1,nq

   call int2char4(iq,tagq)
   fname_default= trim(filppm)//'_Q'//tagq
   fname= trim(fname_default)
   open(file=fname,unit=999,status='new',form='formatted',iostat=ios)

!  Here we use the small non zero q to deal with the long wavelength limit
   allocate(op(3,3,nsym))
   op(:,:,:)=symrec(:,:,:) 
!  FIXME this is broken
!  call findggp(q(:,iq),nsym,op,npwc,gvec,gmet,ngpi,gptab,gptabo,ip2fp)
   deallocate(op)
   index=0

   do ip=1,ngpi
    if (verbose .and. (ip>50 .and. ip<=ngpi-50)) cycle
    ig = ip2fp(1,ip)
    igp= ip2fp(2,ip)
    index=index+1
    call mkppm(ppmodel,ig,igp,iq,nomega,ww,npwc,npwc2,npwc3,gvec,nq,q,bigomegatwsq,omegatw,eigvec,em1w)
    write(999,'(a,i4,a,i8,/,a,3f12.6,/,a,3i6,a,3i6,/,a,/)')&
&    '# index= ',index,' pair number = ',ip,&
&    '# q = ',q(:,iq),&
&    '# G = ',gvec(:,ig),'  G''= ',gvec(:,igp),&
&    '#   w [eV]           Re             Im ' 
    do io=1,nomega
     write(999,'(f8.2,4x,2es16.8)')ww(io)*Ha_eV,em1w(io)
    end do
    write(999,*)
    write(999,*)
   end do 

   close(999)

  end do

  deallocate(em1w,ww)
  deallocate(gvec,q)
  deallocate(bigomegatwsq,omegatw)
  deallocate(gptab,gptabo)
  deallocate(ip2fp) 

! this should be allocated only if ppmodel==3
  deallocate(eigvec)

  stop
 end if 

 end program anascr
!!***
