!{\src2tex{textfont=tt}}
!!****f* ABINIT/wrppm
!! NAME
!! wrppm
!!
!! FUNCTION
!!  input/output plasmon-pole results
!!  This file contains 4 subroutines : 
!! 
!!  wrppm:   Write plasmon-pole parameters PPM on file
!!  testppm: Test a PPM file, return its dimension and rewind the file
!!  rdppm:   Read a PPM file
!!  mkppm:   Calculate the symmetrized inverse dielectric matrix starting 
!!           from the plasmon-pole paramenters
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dtfil=file names (only dtfil%filnam_ds(4) is used)
!!  hdr=header of the file to be written
!!  ppmodel=plasmon-pole type
!!  npwc=number of G vectors in the plasmon-pole model 
!!  npwc2,npwc3=second and third dimension for plasmon-pole parameters 
!!   (different from npwc if ppmodel==3 or ppmodel==4)
!!  gvec(3,npwc)= G vectors in reduced coordinates
!!  nq=number of q-points 
!!  qpoint(3,np)=q-points
!!  bigomegatwsq(npwc,npwc2,nq)= plasmon-parameters (to be better described)
!!  omegatw(npwc,npwc3,nq)= plasmon-pole paramenters (to be better described)
!!  eigvec(npwc,npwc,nq)= eigenvector od the symmetrized inverse dilectric matrix (only if ppmodel==3)
!!
!! OUTPUT
!!  Only write
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  
!!
!! CHILDREN
!! 
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wrppm(dtfil,hdr,ppmodel,npwc,npwc2,npwc3,gvec,nq,qpoint,bigomegatwsq,omegatw,eigvec)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: nq,npwc,npwc2,npwc3,ppmodel
 type(datafiles_type),intent(in) :: dtfil
 type(hdr_type),intent(inout) :: hdr
!arrays 
 integer,intent(in) :: gvec(3,npwc)
 real(dp),intent(in) :: qpoint(3,nq)
 complex(gwpc),intent(in) :: bigomegatwsq(npwc,npwc2,nq)
 complex(gwpc),intent(in) :: omegatw(npwc,npwc3,nq)
 complex(gwpc),intent(in) :: eigvec(npwc,npwc,nq)

!Local variables-------------------------------
 integer :: iq,ii,ig,ig2,ig3,ios,rdwr,readnetcdf
 integer,parameter :: unitpmp=759 
 integer :: fform=1002
 character(len=fnlen) :: filnam
 character(len=500) :: message  

!this quantities should become input 
! integer,intent(in) :: nbndse,npwvece,npwwfne
 integer :: nbndse,npwvece,npwwfne
 character(len=80) :: titem1(2)
 
! *************************************************************************
 
 filnam= trim(dtfil%filnam_ds(4))//'_PMP'

 write(message,'(3a)')' wrppm:  writing plasmon-pole parameters on : ',trim(filnam),ch10
 call wrtout(6,message,'COLL')

!Open file
 open(unitpmp,file=filnam,status='unknown',form='unformatted',iostat=ios)
 if (ios/=0) then 
  write(message,'(4a)')&
&  ' wrppm : ERROR- ',ch10,&
&  '  opening file: ',trim(filnam)
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Write the header
 rdwr= 2
 readnetcdf= 0 ! should make PMP also a netcdf file
 if (readnetcdf==0) then
  call hdr_io(fform,hdr,rdwr,unitpmp)
 else if (readnetcdf==1) then
  call hdr_io_netcdf(fform,hdr,rdwr,unitpmp)
 end if

!General info
 write(unitpmp) ppmodel
 write(unitpmp) npwc,npwc2,npwc3
 write(unitpmp) nq
!here there should be other information related to ppmodel and epms1
!write(unitpmp) titem1
!write(unitpmp) npwvece,npwwfne,nbndse

!Data
 write(unitpmp) gvec(1:3,1:npwc)
 write(unitpmp) (qpoint(:,iq),iq=1,nq)

!Write ppm parameters on file for each q-point
 do iq=1,nq
  write(unitpmp) ((bigomegatwsq(ig,ig2,iq),ig=1,npwc),ig2=1,npwc2)
  write(unitpmp) ((omegatw(ig,ig3,iq),ig=1,npwc),ig3=1,npwc3)
  if (ppmodel==3) then 
   write(unitpmp) ((eigvec(ig,ig2,iq),ig=1,npwc),ig2=1,npwc)
  end if 
 end do 

!Close file and exit
 close(unitpmp)

!DEBUG
!write (std_out,*) ' wrppm : exit'
!stop
!ENDDEBUG

end subroutine wrppm
!!***

!!****f* ABINIT/testppm
!!
!! NAME
!! testppm
!!
!! FUNCTION
!! Test a PPM file, return its dimension and rewind the file
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  unit number of the PPM file 
!!
!! OUTPUT
!!  nq=number of q points
!!  npwc=number of G vectors in the plasmon-pole model 
!!  ppmodel=plasmon-pole type
!!  npwc2,npwc3=second and third dimension for plasmon-pole parameters 
!!
!! PARENTS
!!      
!!
!! CHILDREN
!!      
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine testppm(unitppm,ppmodel,npwc,npwc2,npwc3,nq,hdr)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13io_mpi
!End of the abilint section

 implicit none
    
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unitppm
 integer,intent(out) :: ppmodel,npwc,npwc2,npwc3,nq
 type(hdr_type),intent(out) :: hdr
  
!arrays
!Local variables-------------------------------
!scalars
 integer :: fform,nbndse,nomega,npwvece,npwwfne
 integer :: rdwr,readnetcdf
 character(len=500) :: message
 character(len=80) :: titem1(2)
!arrays

! *************************************************************************

!Read the header
 rdwr= 1
 readnetcdf= 0 ! should make PPM also a netcdf file

 if (readnetcdf==0) then
  call hdr_io(fform,hdr,rdwr,unitppm)
 else if (readnetcdf == 1) then
  call hdr_io(fform,hdr,rdwr,unitppm)
 end if

 if (fform/=1002) then
  write(message,'(5a)')ch10,&
&  ' testppm : ERROR- ',ch10,&
&  ' Unknown file format found ',ch10
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if 

!Echo part of the header
 rdwr= 4 
 call hdr_io(fform,hdr,rdwr,6)

!call hdr_clean(hdr)

!Read dimensions
 read(unitppm,err=100) ppmodel
 read(unitppm,err=100) npwc,npwc2,npwc3
 read(unitppm,err=100) nq

!Read information related to the screening calculation
!read(unitppm,err=100) titem1
!read(unitppm,err=100) npwvece,npwwfne,nbndse,nq,nomega

!Write general info on the screening calculation
 write(message,'(3a)')ch10,' Parameters used to generate the SCR file',ch10
 call wrtout(6,message,'COLL')
 write(message,'(1x,a79,a,1x,a79)')titem1(1)(:79),ch10,titem1(2)(:79)
 call wrtout(6,message,'COLL')
 write(message,'(a,i8)')' dimension of epsilon^-1 matrix   ',npwvece
 call wrtout(6,message,'COLL')
 write(message,'(a,i8)')' number of planewaves used        ',npwwfne
 call wrtout(6,message,'COLL')
 write(message,'(a,i8)')' number of bands used             ',nbndse
 call wrtout(6,message,'COLL')

 rewind(unitppm)
 return

 100 write(message,'(3a)')&
& ' testppm : ERROR - ',ch10,&
& ' file format unknown'
 call wrtout(6,message,'COLL')
 call leave_new('COLL')

 end subroutine testppm
!!***

!!****f* ABINIT/rdppm
!! NAME
!! rdppm
!!
!! FUNCTION
!!  Read a plasmon-pole file 
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  unitppm=unit number of file
!!  ppmodel=plasmon-pole type
!!  npwc=number of G vectors in the plasmon-pole model 
!!  npwc2,npwc3=second and third dimension for plasmon-pole parameters 
!!  nq=number of q-points 
!!
!! OUTPUT
!!  qpoint(3,np)=q-points
!!  gvec(3,npwc)= G vectors
!!  bigomegatwsq(npwc,npwc2,nq)= plasmon-parameters (to be better described)
!!  omegatw(npwc,npwc3,nq)= plasmon-pole paramenters (to be better described)
!!  eigvec(npwc,npwc,nq)= eigenvector od the symmetrized inverse dilectric matrix (only if ppmodel==3)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  
!!
!! CHILDREN
!! 
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine rdppm(unitppm,ppmodel,npwc,npwc2,npwc3,gvec,nq,qpoint,bigomegatwsq,omegatw,eigvec)

    
 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13io_mpi
 use interfaces_14iowfdenpot
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: nq,npwc,npwc2,npwc3,ppmodel,unitppm 
!arrays 
 integer,intent(out) :: gvec(3,npwc)
 real(dp),intent(out) :: qpoint(3,nq)
 complex(gwpc),intent(out) :: bigomegatwsq(npwc,npwc2,nq)
 complex(gwpc),intent(out) :: omegatw(npwc,npwc3,nq)
 complex(gwpc),intent(out) :: eigvec(npwc,npwc,nq)

!Local variables-------------------------------
 integer :: iq,ii,ig,ig2,ig3,ios,rdwr,readnetcdf

 integer :: fform=1002
 type(hdr_type) :: hdr
 character(len=500) :: message  

!this quantities should become input 
! integer,intent(in) :: nbndse,npwvece,npwwfne
! integer :: nbndse,npwvece,npwwfne
! character(len=80) :: titem1(2)
 
! *************************************************************************
 
!DEBUG
 write (std_out,*) ' rdppm : enter'
!ENDDEBUG

!Read the header
 rdwr= 1
 readnetcdf= 0 ! should make PMP also a netcdf file

 if (readnetcdf==0) then
  call hdr_io(fform,hdr,rdwr,unitppm)
 else if (readnetcdf==1) then
  call hdr_io_netcdf(fform,hdr,rdwr,unitppm)
 end if

 call hdr_clean(hdr)

 if(fform/=1002) then
  write(message,'(5a)')ch10,&
&  ' rdppm : ERROR- ',ch10,&
&  ' Unknown file format found ',ch10
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if 

!General info
 read(unitppm) !skip ppmodel
 read(unitppm) !skip npwc,npwc2,npwc3
 read(unitppm) !skip nq
!here there should be other information related to ppmodel and epms1
!read(unitppm) !skip titem1
!read(unitppm) !skip npwvece,npwwfne,nbndse

!Read Data
 read(unitppm) gvec(1:3,1:npwc)
 read(unitppm) (qpoint(:,iq),iq=1,nq)

!Read ppm parameter on file for each q-point
 do iq=1,nq
  read(unitppm) ((bigomegatwsq(ig,ig2,iq),ig=1,npwc),ig2=1,npwc2)
  read(unitppm) ((omegatw(ig,ig3,iq),ig=1,npwc),ig3=1,npwc3)
  if (ppmodel==3) then 
   read(unitppm) ((eigvec(ig,ig2,iq),ig=1,npwc),ig2=1,npwc)
  end if 
 end do 

!Close file and exit
 close(unitppm)

!DEBUG
 write (std_out,*) ' rdppm : exit'
!stop
!ENDDEBUG

end subroutine rdppm
!!***


!!****f* ABINIT/mkppm
!! NAME
!! mkppm
!!
!! FUNCTION
!!  Calculate symmetrized inverse dielectri matrix starting from 
!!  the plasmon-pole paramenters 
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ig,igp= index for the G vectors where the symmetrized inverse dielectric matrix is required
!!  iq= index of the q-point where the symmetrized inverse dielectric matrix is calculated 
!!  ppmodel=plasmon-pole type
!!  npwc=number of G vectors in the plasmon-pole model 
!!  npwc2,npwc3=second and third dimension for plasmon-pole parameters 
!!  nq=number of q-points 
!!  nomega=number of real frequencies to evaluate the symmetrized inverse dielectric matrix using the 
!!   plasmon-pole model 
!!  omega(nomega)= frequencies
!!  qpoint(3,np)=q-points
!!  gvec(3,npwc)= G vectors
!!  bigomegatwsq(npwc,npwc2,nq)= plasmon-parameters (to be better described)
!!  omegatw(npwc,npwc3,nq)= plasmon-pole paramenters (to be better described)
!!  eigvec(npwc,npwc,nq)= eigenvector od the symmetrized inverse dilectric matrix (only if ppmodel==3)
!!
!! OUTPUT
!!  
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Useful to check the accuracy of the plasmon-pole model
!!
!! PARENTS
!! 
!!
!! CHILDREN
!! 
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mkppm(ppmodel,ig,igp,iq,nomega,omega,npwc,npwc2,npwc3,gvec,nq,q,bigomegatwsq,omegatw,eigvec,em1w)
    
 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ppmodel,nq,npwc,npwc2,npwc3,nomega
 integer,intent(in) :: ig,igp,iq 
!arrays
 integer,intent(in) :: gvec(3,npwc)
 real(dp),intent(in) :: q(3,nq)
 real(dp),intent(in) :: omega(nomega) !could be also complex
 complex(gwpc),intent(in) :: bigomegatwsq(npwc,npwc2,nq)
 complex(gwpc),intent(in) :: omegatw(npwc,npwc3,nq)
 complex(gwpc),intent(in) :: eigvec(npwc,npwc,nq)
 complex(gwpc),intent(out) :: em1w(nomega)

!Local variables-------------------------------
 integer :: io,idm
 character(len=500) :: message      
 real(dp) :: repsm1,den
 !this should be read from file, for the moment use default value
 !complex(dpc) :: delta
 real(dp),parameter :: zcut=3.67493260d-03
 complex(gwpc) :: delta
 
! *************************************************************************
 
!DEBUG
!write (std_out,*) ' mkppm : enter'
!ENDDEBUG

 delta=(0.,zcut)
 em1w(:)= (0.,0.)

 if (ppmodel==1 .or. ppmodel==2) then
  do io=1,nomega
   den= omega(io)**2-real(omegatw(ig,igp,iq)**2)
!  if (den**2>zcut**2) then 
   em1w(io)= bigomegatwsq(ig,igp,iq)/den
!  else 
!  FIXME not so sure Im correct here, maybe I need a normalization factor
!  em1w(io)= bigomegatwsq(ig,igp,iq)*&
!  &   ( 1./(omega(io)-real(omegatw(ig,igp,iq))+delta) - 1./omega(io)+real(omegatw(ig,igp,iq))-delta ) 
!  end if 
  end do
  if (ig==igp) em1w= em1w+1.0

 else if (ppmodel==3) then 
  do idm=1,npwc
   do io=1,nomega
    den=omega(io)**2-(omegatw(ig,igp,iq)-delta)**2
    em1w(io)=em1w(io)+eigvec(ig,idm,iq)*conjg(eigvec(igp,idm,iq))*bigomegatwsq(ig,igp,iq)/den
   end do
  end do 

 else if (ppmodel==4) then 
  write(*,*)' ppmodel 4 not yet implemented'
  call leave_new('COLL')
 else 
  write(message,'(5a)')ch10,&
&  ' mkppm : BUG- ',ch10,&
&  ' not allowed value for ppmodel ',ch10
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if 

!DEBUG
!write (std_out,*) ' mkppm : exit'
!ENDDEBUG

end subroutine mkppm
!!***
