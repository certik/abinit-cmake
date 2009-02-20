!{\src2tex{textfont=tt}}
!!****f* ABINIT/testscr
!!
!! NAME
!! testscr
!!
!! FUNCTION
!! Test epsilon-twiddle**-1 file and return its dimension as well as its header
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  unt=Unit number.
!!  fname=File name
!!  localrdwf=input variable (for parallel case). If 1, the SCR file is local to each machine
!!  optfil=0 to read the _SCR  file.
!!         1 to read the _SUSC file.
!!
!! OUTPUT
!!  fform=file format
!!  nbnds_used=number of bands used to calculate the screening (not used, just info)
!!  npweR=number of planewaves used for the screening
!!  npwwfn_used=number of planewaves used during the chi0 run to describe the wavefunctions (only info)
!!  nomegaR=number of frequencies  
!!  nqibzR=number of irred q-points
!!  nqlwlR=Number of q-directions for long-wavelength limit (only if optfil==1)
!!  title=title of the file
!!  Hdr_scr=header of the screening file
!!  qibz_p(3,nqibzR)=q-points reported in the file.
!!  qlwl_p(3,nqlwlR)="small" q-points used for the long-wavelength limit.
!!   warning: In input it should not be neither allocated nor associated 
!!  omega_p(nomegaR)=Calculated frequencies 
!!  gvec_p(3,npweR)=G-vectors in reduced coordinates.
!!
!! TODO
!!  Use file handlers.
!!
!! PARENTS
!!      mrgscr,sigma
!!
!! CHILDREN
!!      hdr_clean,hdr_io
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine testscr(optfil,unt,fname,nqibzR,nqlwlR,nomegaR,npweR,npwwfn_used,nbnds_used,title,&
& fform,MPI_enreg,localrdwf,qibz_p,qlwl_p,omega_p,gvec_p,Hdr_scr)

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert
#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13io_mpi
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: localrdwf,optfil,unt
 integer,intent(out) :: fform,nbnds_used,nomegaR,npweR,npwwfn_used,nqlwlR,nqibzR
 character(len=fnlen),intent(in) :: fname
 type(MPI_type),intent(in) :: MPI_enreg
 type(Hdr_type),intent(out) :: Hdr_scr

!arrays
 integer,pointer :: gvec_p(:,:)
 real(dp),pointer :: qibz_p(:,:),qlwl_p(:,:)
 complex(dpc),pointer :: omega_p(:)
 character(len=80),intent(out) :: title(2)

!Local variables-------------------------------
!scalars
 integer :: ios,rdwr,readnetcdf
 logical :: ltest
 character(len=500) :: msg

!arrays
 integer :: spaceComm,rank,master,ierr
! *************************************************************************

#if defined DEBUG_MODE
!write(msg,'(a)')' testscr : enter'
!call wrtout(std_out,msg,'PERS') 
!call flush_unit(std_out)
#endif
 !
 ! === Test on input variables ===
 ltest=(optfil==0.or.optfil==1) 
 call assert(ltest,'Wrong optfil',__FILE__,__LINE__)
 !
 ! === Initialize MPI related variables ===
 call xcomm_init  (MPI_enreg,spaceComm)   
 call xme_init    (MPI_enreg,rank)            
 call xmaster_init(MPI_enreg,master)   
 !
 ! === Open file ===
 if (rank==master.or.localrdwf==1) then
  write(msg,'(3a)')' testscr : testing file ',TRIM(fname),ch10
  call wrtout(std_out,msg,'COLL')
  open(unit=unt,file=fname,status='old',form='unformatted',iostat=ios)
  if (ios/=0) then 
   write(msg,'(6a)')ch10,&
&   ' testscr : ERROR - ',ch10,&
&   ' opening file ',TRIM(fname),' as old '
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end if  
  !
  ! === Read the header and do consistency check ===
  rdwr=1 ; readnetcdf=0 ! should make EM1 also a netcdf file
  if (readnetcdf==0) then 
   call hdr_io(fform,Hdr_scr,rdwr,unt)
  else if (readnetcdf==1) then 
   STOP "not implemented"
   !call hdr_io(fform,Hdr_scr,rdwr,unt)
  end if
  if (optfil==0) call assert((fform==1002),'Wrong fform for eps^-1',__FILE__,__LINE__)
  if (optfil==1) call assert((fform==1102),'Wrong fform for chi0  ',__FILE__,__LINE__)

  select case (fform)
  case (1002)
   write(msg,'(a)')' epsilon^-1 file in the format em1'
   read(unt,ERR=10)title
   read(unt,ERR=10)npweR,npwwfn_used,nbnds_used,nqibzR,nomegaR
   nqlwlR=0
   allocate(gvec_p(3,npweR))
   read(unt,ERR=10) gvec_p(:,:)
   allocate(qibz_p(3,nqibzR))
   read(unt,ERR=10)qibz_p(:,:)
   nullify(qlwl_p)
   allocate(omega_p(nomegaR))
   read(unt,ERR=10)omega_p(:)

  case (1102)
   write(msg,'(a)')' chi0 file                        '
   read(unt,ERR=10)title
   read(unt,ERR=10)npweR,npwwfn_used,nbnds_used,nqibzR,nomegaR,nqlwlR
   allocate(gvec_p(3,npweR))
   read(unt,ERR=10) gvec_p(:,:)
   allocate(qibz_p(3,nqibzR))
   read(unt,ERR=10)qibz_p(:,:)
   allocate(qlwl_p(3,nqlwlR))
   read(unt,ERR=10)qlwl_p(:,:)
   allocate(omega_p(nomegaR))
   read(unt,ERR=10)omega_p(:)

  case default
   write(msg,'(4a)')ch10,&
&   ' testscr : ERROR - ',ch10,&
&   ' epsilon^-1 file format unknown'
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end select

  call wrtout(std_out,msg,'COLL')
  write(msg,'(1x,a79,a,1x,a79)')title(1)(:79),ch10,title(2)(:79)
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' dimension of matrix       ',npweR
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' number of planewaves used ',npwwfn_used
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' number of bands used      ',nbnds_used
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' number of q-points        ',nqibzR
  call wrtout(std_out,msg,'COLL')
  if (optfil==1) then 
   write(msg,'(a,i8)')' number of q-directions    ',nqlwlR
   call wrtout(std_out,msg,'COLL')
  end if
  write(msg,'(a,i8,a)')' number of omega           ',nomegaR,ch10
  call wrtout(std_out,msg,'COLL')
  close(unt)
 end if !rank or localrdwf 

 if (rank/=master.and.localrdwf==0) then
  !FIXME also title should be casted to avoid problems in parallel if title is used somewhere
  ! requires xcast_mpi_char
  title(1)='not'
  title(2)='important'
 end if
 !
 ! === Cast data to other processors ===
 if (MPI_enreg%nproc>1.and.localrdwf==0) then       
  call xcast_mpi(nqibzR,     master,spaceComm,ierr)
  call xcast_mpi(nqlwlR,     master,spaceComm,ierr)
  call xcast_mpi(nomegaR,    master,spaceComm,ierr)
  call xcast_mpi(npweR,      master,spaceComm,ierr)
  call xcast_mpi(npwwfn_used,master,spaceComm,ierr)
  call xcast_mpi(nbnds_used, master,spaceComm,ierr)
  call xcast_mpi(fform,      master,spaceComm,ierr)
  if (rank/=master.and.localrdwf==0) then 
   ! this proc has not read thus these pointers are not allocated.
   allocate(qibz_p(3,nqibzR),omega_p(nomegaR),gvec_p(3,npweR))
   nullify(qlwl_p) ; if (fform==1102) allocate(qlwl_p(3,nqlwlR))
  end if
  call xcast_mpi(qibz_p, master,spaceComm,ierr)
  call xcast_mpi(omega_p,master,spaceComm,ierr)
  call xcast_mpi(gvec_p, master,spaceComm,ierr)
  if (fform==1102) call xcast_mpi(qlwl_p,master,spaceComm,ierr)
#if defined MPI
  !TODO avoid preprocessor skipping hdr_comm indide the procedure if !MPI
  call hdr_comm(Hdr_scr,master,rank,spaceComm)
#endif
  call leave_test(MPI_enreg) 
 end if

#if defined DEBUG_MODE
 !write(msg,'(a)')' testscr : exit'
 !call wrtout(std_out,msg,'PERS') 
 !call flush_unit(std_out)
#endif
 !
 RETURN
 !
 ! === Something went wrong! ===
10 write(msg,'(4a)')ch10,&
&  ' testscr : ERROR - ',ch10,&
&  ' File seems to be corrupted.'
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')

end subroutine testscr
!!***
