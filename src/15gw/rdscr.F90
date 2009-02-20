!{\src2tex{textfont=tt}}
!!****f* ABINIT/rdscr
!! NAME
!! rdscr
!!
!! FUNCTION
!! Read screening (epsilon-twiddle**-1) file in the SCR format
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, XG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!!  iqiA= (optional) used if only a particular q-point is required. In this case iqiA define the index
!!   of the required q-point in the array qibz(3,nqibz)
!!  nqibzA=number of asked q-points (used to dimension the output arrays) it is equal to nqibz if the full
!!  matrix is required
!!  localrdwf= input variable (for parallel case) if 1, the SCR file is local to each machine
!!  MPI_enreg= datatype gathering information on parallelism
!!  npweA=number of asked planewaves
!!  nomegaA=number of asked frequencies
!!  nqibz=number of q points
!!  verbose= if .true. output more
!!  optfil=0 to write a SCR file containing the symmetrised inverse dielectric matrix
!!         1 to write a CHI0 file containing irreducible polarizability and heads and wings if q==0
!!  nqlwlA=number of small q-directions in head and wings, only if optfil=1
!!  unt=The unit number
!!
!! OUTPUT
!!  qibz(3,nqibz)=coordinates of q points
!!  epsm1(npweA,npweA,nomegaA,nqibz)=epsilon tilde minus 1 matrix, for each frequency and q point
!!  omegaA(nomegaA)=asked frequencies
!!  lwing(npweA,nomegaA,nqlwlA) (optional)=lower and upper wings, first element is the head, only if optfil=1
!!
!! NOTES
!!  If the epsilon matrix read is bigger than npweA x npweA, it will be truncated; 
!!  if it is smaller, an error will occur
!!  
!!  If the number of frequencies asked for is smaller than that reported in the file, the matrix 
!!  will be truncated. If nomegaA > nomegaR an error will occur 
!!
!! PARENTS
!!      csigme,mrgscr,sigma
!!
!! CHILDREN
!!      assert,hdr_clean,hdr_io,hdr_io_netcdf,leave_test,memerr,wrtout
!!      xcast_mpi,xcomm_init,xmaster_init,xme_init
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine rdscr(optfil,unt,fname,npweA,nqibz,nqibzA,nomegaA,qibz,omegaA,gmet,&
& epsm1,MPI_enreg,localrdwf,nqlwlA,verbose,qlwl,uwing,lwing,&
& iqiA) ! Optional

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : GW_TOLQ0
 use m_errors, only : assert,assert_eq
#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: localrdwf,nomegaA,npweA,nqibz,nqibzA,optfil,nqlwlA,unt
 integer,optional,intent(in) :: iqiA 
 logical,intent(in) :: verbose
 character(len=fnlen),intent(in) :: fname
 type(MPI_type),intent(in) :: MPI_enreg
!arrays
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(out) :: qibz(3,nqibz)
 real(dp),intent(inout) ::  qlwl(3,nqlwlA*optfil)
 complex(gwpc),intent(inout) :: epsm1(npweA,npweA,nomegaA,nqibzA),omegaA(nomegaA)
 complex(gwpc),intent(inout) :: lwing(npweA,nomegaA,nqlwlA*optfil) 
 complex(gwpc),intent(inout) :: uwing(npweA,nomegaA,nqlwlA*optfil)
  
!Local variables-------------------------------
!scalars
 integer :: fform,ii,io,ios,iq,istat,iqsm
 integer :: nbndsR,nomegaR,npweR,npwwfnR,nqR,nqsm,nqlwlR
 integer :: rdwr,readnetcdf
 integer :: spaceComm,rank,master,ierr
 real(dp) :: dielconst
 logical :: lqeq0,ltest,i_read,master_cast,read_qslice
 character(len=500) :: msg
 character(len=50),parameter :: FILE__='rdscr.F90'
 type(Hdr_type) :: Hdr_scr
!arrays
 real(dp) :: qd(3,nqibz)
 real(dp),allocatable :: qsmall(:,:)
 character(len=80) :: titem1(2)
 complex(dpc),allocatable :: epsm1d(:,:),omegad(:)
 complex(dpc),allocatable :: wingd(:)
! *************************************************************************

#if defined DEBUG_MODE
!write(msg,'(a)')' rdscr : enter'
!call wrtout(std_out,msg,'PERS') 
!call flush_unit(std_out)
#endif
 !
 ! === Check input ===
 ltest=(optfil==0.or.optfil==1) 
 call assert(ltest,'Wrong optfil',__FILE__,__LINE__)

 call xcomm_init  (MPI_enreg,spaceComm)  
 call xme_init    (MPI_enreg,rank     )          
 call xmaster_init(MPI_enreg,master   )  
 !
 ! === Define IO for parallel execution ===
 i_read     =.FALSE. 
 master_cast=.FALSE.
 if (localrdwf==1.or.rank==master)            i_read=.TRUE.
 if (localrdwf==0.and.MPI_enreg%nproc>1) master_cast=.TRUE.
 !
 ! === Are we reading a _SCR or _CHI0 file? ===
 !unt=Dtfil%unscr ; fname=Dtfil%filscr
 !if (optfil==1) then 
 ! unt=Dtfil%unchi0 ; fname=Dtfil%filchi0
 !end if
 !
 ! === Slice or full array? ===
 read_qslice=.FALSE.
 if (PRESENT(iqiA)) then 
  read_qslice=.TRUE.
  write(msg,'(2a,i5,2a)')ch10,&
&  ' rdscr : reading the slice corresponding to iq = ',iqiA,' from file : ',TRIM(fname)
  call wrtout(std_out,msg,'COLL')
  ltest=(iqiA>0.and.iqiA<=nqibz)
  call assert(ltest,'iqiA out of range',__FILE__,__LINE__)
 end if 

 allocate(omegad(nomegaA))
 
 if (i_read) then
  if (verbose) then 
   write(msg,'(5a)')ch10,&
&   ' rdscr: reading screening (epsilon-twiddle^-1) ',&
&   ' file in the SCR format from : ',TRIM(fname),ch10
   call wrtout(std_out,msg,'COLL')
  end if
  open(unit=unt,file=TRIM(fname),status='old',form='unformatted',iostat=ios)
  if (ios/=0) then 
   write(msg,'(6a)')ch10,&
&   ' rdscr : ERROR -',ch10,&
&   ' opening file ',TRIM(fname),' as old'
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end if  
  ! 
  ! === Read the header of the file ===
  rdwr=1 ; readnetcdf=0 ! should become input parameter and make EM1 a netcdf file as well
  if (readnetcdf==0) then 
   call hdr_io(fform,Hdr_scr,rdwr,unt)
  else if (readnetcdf==1) then 
   call hdr_io_netcdf(fform,Hdr_scr,rdwr,unt)
  end if 

  read(unt,ERR=10)titem1
  if (optfil==0) read(unt,ERR=10)npweR,npwwfnR,nbndsR,nqR,nomegaR
  if (optfil==1) read(unt,ERR=10)npweR,npwwfnR,nbndsR,nqR,nomegaR,nqlwlR

  if (npweR>npweA) then
   write(msg,'(4a,i8,2a,i8)')ch10,&
&   ' rdscr : COMMENT - ',ch10,&
&   ' Total number of G in the matrix       = ',npweR,ch10,&
&   ' Reading a smaller matrix of dimension = ',npweA
   call wrtout(std_out,msg,'COLL')
  end if
  ltest=(npweR>=npweA)
  write(msg,'(4a,i8,2a,i8)')ch10,&
&  ' rdscr : ERROR - ',ch10,&
&  ' dimension of matrix        = ',npweR,ch10,&
&  ' requiring a too big matrix = ',npweA
  call assert(ltest,'Requiring too much frequencies',__FILE__,__LINE__)
  ! === Do some check ===
  if (nqR/=nqibz.or.nqR<nqibzA) then 
   write(msg,'(4a)')ch10,&
&   ' rdscr : ERROR - ',ch10,&
&   ' requiring too much q-points'
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end if 
  ltest=(nomegaR>=nomegaA) 
  call assert(ltest,'Requiring too much frequencies',__FILE__,__LINE__)
  if (optfil==1) then 
   ltest=(nqlwlA>=nqlwlR)
   call assert(ltest,'Requiring too much nqlwl',__FILE__,__LINE__)
  end if

  read(unt,ERR=10) !Skip gvec
  read(unt,ERR=10)qd(1:3,1:nqibz)
  qibz(:,:)=qd(:,:)

  ! === Additional record for chi0 files ===
  if (optfil==1) then 
   allocate(qsmall(3,nqlwlR))
   read(unt,ERR=10)qsmall(1:3,1:nqlwlR)
   qlwl(:,:)=qsmall(:,1:nqlwlA)
   deallocate(qsmall)
   write(*,*)' Small q-directions ',qlwl(:,:)
  end if

  read(unt,ERR=10) omegad(1:nomegaA)
  omegaA(:)=omegad(:)

  if (verbose) then 
   write(msg,'(2a)')ch10,' q-points [reciprocal lattice units]:'
   call wrtout(std_out,msg,'COLL') 
   do iq=1,nqibz 
    write(msg,'(i5,3f12.6)')iq,(qibz(ii,iq),ii=1,3)
    call wrtout(std_out,msg,'COLL') 
   end do
   write(msg,'(2a)')ch10,' frequencies used [eV]:'
   call wrtout(std_out,msg,'COLL') 
   do io=1,nomegaA
    write(msg,'(i3,2f7.2)')io,REAL(omegaA(io))*Ha_eV,AIMAG(omegaA(io))*Ha_eV
   end do
  end if 
 end if ! i_read

 if (master_cast) then
  call xcast_mpi(qibz   ,master,spaceComm,ierr)
  call xcast_mpi(npweR  ,master,spaceComm,ierr)
  call xcast_mpi(omegaA ,master,spaceComm,ierr)
  call xcast_mpi(nomegaR,master,spaceComm,ierr)
  if (optfil==1) call xcast_mpi(qlwl,master,spaceComm,ierr)
 end if
 !
 ! === Now read epsilon^-1 ===
 allocate(epsm1d(npweR,npweR),STAT=istat) 
 if (istat/=0) call memerr(FILE__,'epsm1d',npweR**2,'dpc')
 !
 ! === Two coding for different case just to keep it readable ===
 ! TODO re-merge the two cases.
 SELECT CASE (read_qslice)
 CASE (.TRUE.)
  ! === Read only a slice of the full array (useful if full array is huge) ===
  if (optfil==1) STOP 'not implemented'
  do iq=1,nqibz

   if (iq==iqiA) then 
    lqeq0=(normv(qibz(:,iq),gmet,'G')<GW_TOLQ0)
    do io=1,nomegaA
     if (i_read) read(unt,ERR=10) epsm1d(1:npweR,1:npweR)
     if (master_cast) call xcast_mpi(epsm1d,master,spaceComm,ierr)
     epsm1(1:npweA,1:npweA,io,1)=epsm1d(1:npweA,1:npweA)
    end do
    ! Skip other frequencies
    do io=nomegaA+1,nomegaR ; if (i_read) read(unt,ERR=10) ; end do
   else 
    ! Skip epsm1d(1:npweR,1:npweR)
    do io=1,nomegaR ; if (i_read) read(unt,ERR=10) ; end do
   end if ! iq==iqiA 

  end do !iq
 CASE (.FALSE.)
  ! === Read the entire array ===
  do iq=1,nqibz

   lqeq0=(normv(qibz(:,iq),gmet,'G')<GW_TOLQ0)

   do io=1,nomegaA
    if (i_read) read(unt,ERR=10) epsm1d(1:npweR,1:npweR)
    if (master_cast)  call xcast_mpi(epsm1d,master,spaceComm,ierr)
    epsm1(1:npweA,1:npweA,io,iq)=epsm1d(1:npweA,1:npweA)
   end do
   ! Skip other frequencies
   do io=nomegaA+1,nomegaR ; if (i_read) read(unt,ERR=10) ; end do
   !
   ! === Read heads and wings ===
   if (optfil==1.and.lqeq0) then 
    allocate(wingd(npweR))
    do iqsm=1,nqsm
     do io=1,nomegaA
      if (i_read) read(unt,ERR=10) wingd(1:npweR)
      if (master_cast)  call xcast_mpi(wingd,master,spaceComm,ierr)
      uwing(1:npweA,io,iqsm)=wingd(1:npweA)
      if (i_read) read(unt,ERR=10) wingd(1:npweR)
      if (master_cast)  call xcast_mpi(wingd,master,spaceComm,ierr)
      lwing(1:npweA,io,iqsm)=wingd(1:npweA)
     end do
     do io=nomegaA+1,nomegaR
      if (i_read) read(unt,ERR=10) ! skip uwing
      if (i_read) read(unt,ERR=10) ! skip lwing
     end do ! io
    end do !iqsm
    deallocate(wingd)
   end if

  end do !iqibz
 END SELECT

 if (i_read) close(unt)                    
 if (master_cast) then 
  call leave_test(MPI_enreg) 
 end if

 if (optfil==0.and.verbose) then
  if (ABS(REAL(omegaA(1)))<0.01.and.ABS(AIMAG(omegaA(1)))<0.01) then
   do iq=1,nqibz
    if (ALL(ABS(qibz(:,iq))<0.001)) then
     dielconst=one/REAL(epsm1(1,1,1,iq))
     write(msg,'(a,f6.2)')' dielectric constant with LF = ', dielconst
     call wrtout(std_out,msg,'COLL')
    end if
   end do
  end if
  write(msg,'(a)')' epsilon-twiddle^-1 read '
  call wrtout(std_out,msg,'COLL')
 end if 
 !
 ! === Free memory ===
 if (allocated(epsm1d)) deallocate(epsm1d)
 if (allocated(omegad)) deallocate(omegad)
 if (i_read) call hdr_clean(Hdr_scr)

#if defined DEBUG_MODE
!write(msg,'(a)')' rdscr : exit'
!call wrtout(std_out,msg,'PERS') 
!call flush_unit(std_out)
#endif

 RETURN
 !
 ! === Something went wrong! ===
10 write(msg,'(4a)')ch10,&
&  ' rdscr : ERROR - ',ch10,&
&  ' File seems to be corrupted.'
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')

end subroutine rdscr
!!***
