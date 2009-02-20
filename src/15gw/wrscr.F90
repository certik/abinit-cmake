!{\src2tex{textfont=tt}}
!!****f* ABINIT/wrscr
!! NAME
!! wrscr
!!
!! FUNCTION
!! For a single q-point, write either \tilde epsilon^{-1} on the _SCR file 
!! or chi0 on the _SUSC file according to the value of optfil. 
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR,VO,LR,RWG,R.Shaltaf,MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  unt=the unit number.
!!  epsm1(npwe,npwe,nomega)=screening matrix, for different frequencies, and q-point iq
!!  gvec(3,npwe)=coordinates of G vectors
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!!  Hdr=header of the file to be written
!!  iq=index, in the array qibz, of the q-point where epsm1 is going to be written. 
!!   It is used to decide whether the Hdr must be written or not.
!!  nbnds_used=number of bands used (only for info)
!!  nomega=number of frequencies
!!  npwe=number of plane waves for epsilon (input variable)
!!  npwwfn_used=number of plane waves used to describe the wavefunction (only for info)
!!  nqibz=number of q-points
!!  nqlwl=Number of small q-points used to calculate the long-wavelenght limit (only if optfil==2)
!!  omega(nomega)=set of frequencies
!!  qibz(3,nqibz)=reduced coordinates of irred q-points
!!  qlwl(3,nqlwl)="Small" q-points for the long-wavelenght limit.
!!  title(2)=title (to be described)
!!  optfil=0 to write a _SCR  file containing the symmetrised inverse dielectric matrix
!!         1 to write a _SUSC file containing irreducible polarizability and heads and wings if q==0
!!  lwing(npwe,nomega,Ep%nqlwl)=Lower and upper wings, first element is the head, 
!!
!! OUTPUT
!!  (only writing on file)
!!
!! PARENTS
!!      mrgscr,screening
!!
!! CHILDREN
!!      hdr_io,hdr_io_netcdf
!!
!! TODO 
!! Use file handlers.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wrscr(iq,optfil,unt,fname,Hdr,Dtset,npwe,npwwfn_used,nbnds_used,&
& nqibz,nqlwl,nomega,qibz,omega,gvec,gmet,epsm1,title,&
& qlwl,uwing,lwing)

 use defs_basis
 use m_gwdefs, only : GW_TOLQ, GW_TOLQ0
 use defs_datatypes
 use m_errors, only : assert, assert_eq
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_13io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq,nbnds_used,nomega,npwe,npwwfn_used,nqibz,nqlwl,optfil
 integer,intent(in) :: unt
 type(Dataset_type),intent(in) :: Dtset
 type(Hdr_type),intent(inout) :: Hdr
 character(len=fnlen),intent(in) :: fname
!arrays
 integer,intent(in) :: gvec(3,npwe)
 real(dp),intent(in) :: gmet(3,3),qibz(3,nqibz)
 real(dp),intent(in) :: qlwl(3,nqlwl*optfil)
 complex(gwpc),intent(in) :: epsm1(npwe,npwe,nomega),omega(nomega)
 complex(gwpc),intent(in) :: lwing(npwe,nomega,nqlwl*optfil),uwing(npwe,nomega,nqlwl*optfil)
 character(len=80),intent(in) :: title(2)

!Local variables-------------------------------
 character(len=50),parameter :: FILE__='wrscr.F90'
!scalars
 integer :: fform,has_q0,io,ios,iqdm,iqs,istat,itst
 integer :: nq_loc,nqlwl_loc,rdwr
 logical :: lqeq0,ltest,is_partial,first_qptdm1
 character(len=500) :: msg

!arrays
 real(dp) :: qdiff(3)
 real(dp),allocatable :: qd_temp(:,:)
 complex(dpc),allocatable :: epsm1d(:,:),omegad(:),wingd(:)

! *************************************************************************

#if defined DEBUG_MODE
! write(msg,'(a)')' wrscr : enter '
! call wrtout(std_out,msg,'COLL') 
! call flush_unit(std_out)
#endif
 !
 ! === Define the file format: e^-1 or chi0 ===
 ltest=(optfil==0.or.optfil==1) 
 call assert(ltest,'Wrong value for optfil',__FILE__,__LINE__)
 if (optfil==0) fform=1002
 if (optfil==1) fform=1102
 !
 ! === Check if a partial file has to be written ===
 ! * first_qptdm1 is .TRUE. if we are writing the first point in the qptdm array.
 ! * has_q0 is 1 if Gamma the point is contained in qptdm, if .not.partial always 1
 has_q0=1
 is_partial  =.FALSE. 
 first_qptdm1=.FALSE. 

 if (Dtset%nqptdm>0.and.Dtset%nqptdm<nqibz) then 
  is_partial=.TRUE. 
  qdiff(:)=qibz(:,iq)-Dtset%qptdm(:,1)
  if (ALL(ABS(qdiff(:))<1.0e-3)) first_qptdm1=.TRUE.
  !if (normv(qdiff,gmet,'G')<GW_TOLQ) first_qptdm1=.TRUE.
  has_q0=0
  do iqdm=1,Dtset%nqptdm
   if (normv(Dtset%qptdm(:,iqdm),gmet,'G')<GW_TOLQ0) then 
    has_q0=has_q0+1
   end if
  end do
  ltest=(has_q0==0.or.has_q0==1)
  call assert(ltest,'Too much q0 points',__FILE__,__LINE__)
 end if
 !
 ! === _SCR or _SUSC file? ===
 ! FIXME this is not treated correctly, problem with mrgscr
 !unt=Dtfil%unscr 
 !fname=TRIM(Dtfil%filnam_ds(4))//'_SCR'
 if (optfil==1) then
  !unt=Dtfil%unchi0 
  !fname=TRIM(Dtfil%filnam_ds(4))//'_SUS'
  lqeq0=(normv(qibz(:,iq),gmet,'G')<GW_TOLQ0)
  if (has_q0==1) call assert((nqlwl>=1),'nqlwl must be >=1',__FILE__,__LINE__)
 end if
 !
 ! === If iq=1 or first element in qptdm, open the file and write info ===
 if (iq==1.or.first_qptdm1) then

  if (is_partial) then
   nq_loc=Dtset%nqptdm
   allocate(qd_temp(3,nq_loc))
   qd_temp(:,:)=Dtset%qptdm(:,1:nq_loc) 
  else
   nq_loc=nqibz
   allocate(qd_temp(3,nqibz)) 
   qd_temp(:,:)=qibz(:,:)
  end if

  open(unit=unt,file=TRIM(fname),status='unknown',form='unformatted',iostat=ios)
  if (ios/=0) then 
   write(msg,'(6a)')ch10,&
&   ' wrscr : ERROR -',ch10,&
&   ' opening file ',TRIM(fname),' as unknown-unformatted'
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
  end if  
  ! 
  ! === Write the header of the SCR/SUSC file ===
  rdwr=2
  if (Dtset%accesswff==0) then 
   call hdr_io(fform,Hdr,rdwr,unt)
  else if (Dtset%accesswff==1) then 
   call hdr_io_netcdf(fform,Hdr,rdwr,unt)
  end if
  !if (Dtset%accesswff==3) call hdr_io_etsf(fform,Hdr,rdwr,unt)

  ! === Write dimensions, gvec and q-points ===
  select case (optfil)
   case (0)
    write(unt)title
    write(unt)npwe,npwwfn_used,nbnds_used,nq_loc,nomega
   case (1) 
    nqlwl_loc=0 ; if (has_q0==1) nqlwl_loc=nqlwl
    write(unt)title
    !write(unt)title,nI,nJ,ID,approx_type,tordering,test_type
    write(unt)npwe,npwwfn_used,nbnds_used,nq_loc,nomega,nqlwl_loc 
  end select

  write(unt)gvec(1:3,1:npwe)

  write(unt)qd_temp(1:3,1:nq_loc)
  ! === For chi0 and q-->0 add q-points for heads and wings ===
  ! * WARNING this wont work if qptdm==.TRUE. and the first point is not gamma
  if (optfil==1.and.has_q0==1) write(unt)qlwl(1:3,1:nqlwl)

  allocate(omegad(nomega)) 
  omegad(:)=omega(:)
  write(unt)omegad(1:nomega)
  deallocate(qd_temp,omegad)
 end if ! iq==1 or first point in qptdm
 !
 ! === Write a record for each omega ===
 ! If chi0 and q==0, zero heads and wings as they are written separately.
 allocate(epsm1d(npwe,npwe),STAT=istat) 
 if (istat/=0) call memerr(FILE__,'epsm1d',npwe**2,'gwpc')

 do io=1,nomega
  epsm1d(:,:)=epsm1(:,:,io) !spc ==> dpc
  if (optfil==1.and.lqeq0) then 
   epsm1d(1,:)=czero
   epsm1d(:,1)=czero
  end if
  write(unt)epsm1d(1:npwe,1:npwe)
 end do
 deallocate(epsm1d)
 !
 ! === For chi0 files and q==>0 write heads and wings ===
 if (optfil==1.and.lqeq0) then
  allocate(wingd(npwe))
  do iqs=1,nqlwl
   do io=1,nomega
    wingd(:)=uwing(:,io,iqs)  !spc --> dpc
    write(unt)wingd(1:npwe)
    wingd(:)=lwing(:,io,iqs)  !spc --> dpc
    write(unt)wingd(1:npwe)
   end do
  end do
  deallocate(wingd)
 end if

 ! === If last q-point, close the file ===
 if (iq==nqibz) close(unt)

#if defined DEBUG_MODE
! write(msg,'(a)')' wrscr : exit'
! call wrtout(std_out,msg,'COLL') 
! call flush_unit(std_out)
#endif

end subroutine wrscr
!!***
