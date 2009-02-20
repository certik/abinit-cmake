!{\src2tex{textfont=tt}}
!!****f* ABINIT/kpgio
!! NAME
!! kpgio
!!
!! FUNCTION
!! Do initialization of kg information.
!! Includes opening disk files for kpgsph i/o.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, AR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ecut=kinetic energy planewave cutoff (hartree)
!!  exchn2n3d=if 1, n2 and n3 are exchanged
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kgnam=name of unkg file
!!  kptns(3,nkpt)=reduced coords of k points
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  character(len=4) : mode_paral=either 'COLL' or 'PERS', tells whether
!!   the loop over k points must be done by all processors or not,
!!   in case of parallel execution.
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum number of planewaves as dimensioned in calling routine
!!  nband(nkpt*nsppol)=number of bands at each k point
!!  nkpt=number of k points
!!  nsppol=1 for unpolarized, 2 for polarized
!!  unkg=unit number for storage of basis sphere data: stores indirect
!!   indexing array and integer coordinates for all planewaves in basis
!!   sphere for each k point being considered
!!
!! OUTPUT
!!  npwarr(nkpt)=array holding npw for each k point, taking into account
!!   the effect of istwfk, and the spreading over processors
!!  npwtot(nkpt)=array holding the total number of plane waves for each k point,
!!  kg(3,mpw*mkmem)=dimensionless coords of G vecs in basis sphere at k point
!!
!! NOTES
!! Note that in case of band parallelism, the number of spin-up
!! and spin-down bands must be equal at each k points
!!
!! PARENTS
!!      gstate,initberry3,initmv,loper3,nonlinear,overlap_wf
!!      partial_dos_fractions,respfn,scfcv,suscep,wffile,wfread
!!
!! CHILDREN
!!      kpgsph,leave_test,timab,wrtout,xcomm_init,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine kpgio(ecut,exchn2n3d,gmet,istwfk,kg,kgnam,kptns,mkmem,nband,nkpt,&
& mode_paral,mpi_enreg,mpw,npwarr,npwtot,nsppol,unkg)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_13recipspace, except_this_one => kpgio
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: exchn2n3d,mkmem,mpw,nkpt,nsppol,unkg
 real(dp),intent(in) :: ecut
 character(len=4),intent(in) :: mode_paral
 character(len=fnlen),intent(in) :: kgnam
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol)
 integer,intent(out) :: kg(3,mpw*mkmem),npwarr(nkpt),npwtot(nkpt)
 real(dp),intent(in) :: gmet(3,3),kptns(3,nkpt)

!Local variables-------------------------------
!scalars
 integer :: ierr,ikg,ikpt,istwf_k,me,nband_down,nband_k,npw1,npw2,spaceComm
 character(len=500) :: message
!arrays
 integer,allocatable :: kg_disk(:,:)
 real(dp) :: kpoint(3),tsec(2)

! *************************************************************************

!BEGIN TF_CHANGES
!Define me
 call xme_init(mpi_enreg,me)
!END TF_CHANGES

!PATCH kpgio // KPT & FFT me-->me_kpt
 if ((mpi_enreg%paral_compil_kpt==1) .and. &
& (mpi_enreg%paral_compil_fft==1)) then
  me = mpi_enreg%me_kpt
 end if

 if(mpi_enreg%paral_compil_kpt==1)then
  if((mpi_enreg%paralbd==1) .and. (mode_paral=='PERS')) then
   if(nsppol==2)then
    do ikpt=1,nkpt
     nband_k=nband(ikpt)
     nband_down=nband(ikpt+nkpt)
     if(nband_k/=nband_down)then
      write(message,'(a,a,a,a,a,a,a,a,i4,a,i4,a,a,a,i4,a,a,a)')ch10,&
&      ' kpgio : ERROR -',ch10,&
&      '  Band parallel case, one must have same number',ch10,&
&      '  of spin up and spin down bands, but input is :',ch10,&
&      '  nband(up)=',nband_k,', nband(down)=',nband_down,',',ch10,&
&      '  for ikpt=',ikpt,'.',ch10,&
&      '  Action : correct nband in your input file.'
     end if
    end do
   end if
  end if
 end if
 npwarr(:)=0
 npwtot(:)=0


!In case the information are to be kept on disk, open the disk file
!and allocate the work space
 if (mkmem==0) then
  open (unit=unkg,file=kgnam,form='unformatted',status='unknown')
  rewind (unit=unkg)
  allocate(kg_disk(3,mpw))
 end if

 ikg=0
!Find (k+G) sphere for each k.

 do ikpt=1,nkpt

  nband_k = nband(ikpt)

  if(mpi_enreg%paral_compil_kpt==1)then
   if(mode_paral=='PERS')then
    if (mpi_enreg%parareel == 0) then
!    BEGIN TF_CHANGES
     if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,1:nsppol)-me))/=0) then
!     END TF_CHANGES
      cycle
     end if
    else
     if(mpi_enreg%proc_distrb_para(mpi_enreg%ipara,ikpt)/= mpi_enreg%me) then
      cycle
     end if
    end if
   end if
  end if

  kpoint(:)=kptns(:,ikpt)
  istwf_k=istwfk(ikpt)

  if(mkmem/=0)then
   call kpgsph(ecut,exchn2n3d,gmet,ikg,ikpt,istwf_k,kg,kpoint,&
&   mkmem,mpi_enreg,mpw,npw1)
  else
   call kpgsph(ecut,exchn2n3d,gmet,ikg,ikpt,istwf_k,kg_disk,kpoint,&
&   1,mpi_enreg,mpw,npw1)
  end if

! This section should be simplified !!
  if(mpi_enreg%paral_compil_kpt==0)then
   npwarr(ikpt)=npw1
  else if(mpi_enreg%paral_compil_kpt==1)then
   if (mpi_enreg%parareel == 0) then
!   BEGIN TF_CHANGES
    if((minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,1:nsppol))== &
&    me) .and. (mode_paral=='PERS')) then
!    END TF_CHANGES
     npwarr(ikpt)=npw1
    end if
   else
    if((mpi_enreg%proc_distrb_para(mpi_enreg%ipara,ikpt) &
&    == mpi_enreg%me) .and. (mode_paral=='PERS')) then
     npwarr(ikpt)=npw1
    end if
   end if
  end if

! Make sure npw < nband never happens:
! if (npw1<nband(ikpt)) then
! write(message, '(a,a,a,a,i5,a,3f8.4,a,a,i10,a,i10,a,a,a,a)' )ch10,&
! &   ' kpgio : ERROR -',ch10,&
! &   '  At k point number',ikpt,' k=',(kptns(mu,ikpt),mu=1,3),ch10,&
! &   '  npw=',npw1,' < nband=',nband(ikpt),ch10,&
! &   '  Indicates not enough planewaves for desired number of bands.',ch10,&
! &   '  Action : change either ecut or nband in input file.'
! call wrtout(06,message,mode_paral)
! call leave_new(mode_paral)
! end if

! Find boundary of G sphere for efficient zero padding,
  if (mkmem/=0) then

!  Shift to next section of each array kg
   ikg=ikg+npw1

  else

!  Write G sphere data to disk file when mkmem==1
   write(unkg) npw1
   write(unkg)
   write(unkg) kg_disk(1:3,1:npw1)

  end if

  if(mpi_enreg%paral_compil_kpt==1)then
!  This is subtle : before the MPI_ALLREDUCE, only one
!  processor can contain the correct non-zero npw(:,nkpt)
   if (mpi_enreg%parareel == 0) then
!   BEGIN TF_CHANGES
    if((minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,1:nsppol) &
&    -me))/=0) .and. mode_paral=='PERS') then
!    END TF_CHANGES
     npwarr(ikpt)=0
    end if
   else
    if((mpi_enreg%proc_distrb_para(mpi_enreg%ipara,ikpt) &
&    /= mpi_enreg%me) .and. mode_paral=='PERS') then
     npwarr(ikpt)=0
    end if
   end if
  end if


! End of the loop over k points
 end do

 if(mkmem==0)deallocate(kg_disk)

 if(mpi_enreg%paral_compil_kpt==1)then
  call timab(62,1,tsec)
  if(mode_paral == 'PERS') then
!  BEGIN TF_CHANGES
   call leave_test(mpi_enreg)
!  END TF_CHANGES
!  Recreate full npwarr on all proc.
   call timab(48,1,tsec)

!  PATCH kpgio // KPT & FFT SUM NPWARR"
   if ((mpi_enreg%paral_compil_kpt==1) .and. &
&   (mpi_enreg%paral_compil_fft==1)) then
    call xsum_mpi(npwarr,mpi_enreg%comm_kpt,ierr)
   else
    call xcomm_init(mpi_enreg,spaceComm)
    call xsum_mpi(npwarr,spaceComm,ierr)
   end if

   call timab(48,2,tsec)
  end if
  write(message, '(a)' ) 'kpgio: loop on k-points done in parallel'
  call wrtout(06,message,'COLL')
 end if

!XG030513 MPIWF : now, one should sum npwarr over all processors
!of the WF group, to get npwtot (to be spread on all procs of the
!WF group
 npwtot(:)=npwarr(:)

!Taking into account istwfk
 do ikpt=1,nkpt
  if(istwfk(ikpt)>1)then
   if(istwfk(ikpt)==2)then
    npwtot(ikpt)=2*npwtot(ikpt)-1
   else
    npwtot(ikpt)=2*npwtot(ikpt)
   end if
  end if
 end do

end subroutine kpgio
!!***
