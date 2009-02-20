!{\src2tex{textfont=tt}}
!!****f* ABINIT/outvars
!! NAME
!! outvars
!!
!! FUNCTION
!! Echo variables for the ABINIT code.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  choice= 1 if echo of preprocessed variables, 2 if echo after call driver
!!  dmatpuflag=flag controlling the use of an initial density matrix in PAW+U (max. value over datasets)
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  iout=unit number for echoed output
!!  istatr=repetition rate for status file
!!  istatshft=shift of the repetition rate for status file
!!  mpi_enreg=informations about MPI parallelization
!!  mxlpawu=maximal value of input lpawu for all the datasets
!!  mxmband=maximum number of bands
!!  mxnatom=maximal value of input natom for all the datasets
!!  mxnatpawu=maximal value of number of atoms on which +U is applied for all the datasets
!!  mxnatsph=maximal value of input natsph for all the datasets
!!  mxnatvshift=maximal value of input natvshift for all the datasets
!!  mxncenter=maximal value of input ncenter for all the datasets
!!  mxnconeq=maximal value of input nconeq for all the datasets
!!  mxnkptgw=maximal value of input nkptgw for all the datasets
!!  mxnkpt=maximal value of input nkpt for all the datasets
!!  mxnorb=maximal value of input norb for all the datasets
!!  mxnqptdm=maximal value of input nqptdm for all the datasets
!!  mxnspinor=maximal value of input nspinor for all the datasets
!!  mxnsppol=maximal value of input nsppol for all the datasets
!!  mxnsym=maximum number of symmetries
!!  mxntypat=maximum number of type of atoms
!!  ndtset=number of datasets
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!   one data set. Use for most dimensioned arrays.
!!  npsp=number of atom types
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!  results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including evolving variables
!!  timopt=input variable to modulate the timing
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Note that this routine is called only by the processor me==0 .
!! In consequence, no use of message and wrtout routine.
!! The lines of code needed to output the defaults are preserved
!! (see last section of the routine, but are presently disabled)
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      outvar1,prtocc,prttagm,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine outvars (choice,dmatpuflag,dtsets,iout,istatr,istatshft,mpi_enreg,&
&  mxlpawu,mxmband,mxnatom,mxnatpawu,mxnatsph,mxnatvshift,mxncenter,mxnconeq,mxnkptgw,mxnkpt,&
&  mxnorb,mxnnos,mxnqptdm,mxnspinor,mxnsppol,mxnsym,mxntypat,&
&  ndtset,ndtset_alloc,npsp,pspheads,results_out,timopt)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_13iovars, except_this_one => outvars
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,dmatpuflag,iout,istatr,istatshft,mxlpawu,mxmband
 integer,intent(in) :: mxnatom,mxnatpawu,mxnatsph,mxnatvshift,mxncenter,mxnconeq,mxnkpt
 integer,intent(in) :: mxnkptgw,mxnnos,mxnorb,mxnqptdm,mxnspinor,mxnsppol,mxnsym,mxntypat
 integer,intent(in) :: ndtset,ndtset_alloc,npsp,timopt
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(pspheader_type),intent(in) :: pspheads(npsp)
 type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)

!Local variables-------------------------------
 character(len=*), parameter :: format01110 ="(1x,a9,1x,(t13,8i8) )"
 character(len=*), parameter :: format01150a="(1x,a9,a,1x,(t13,3es16.8))"
 character(len=*), parameter :: format01155 ="(1x,a9,1x,(t13,10i5))"
 character(len=*), parameter :: format01155a="(1x,a9,a,1x,(t13,10i5))"
 character(len=*), parameter :: format01160 ="(1x,a9,1x,(t13,3es18.10)) "
 character(len=*), parameter :: format01160a="(1x,a9,a,1x,(t13,3es18.10)) "
!scalars
 integer,parameter :: nkpt_max=50
 integer :: dmatsize,first,iat,iban,icount,idtset,ii,ikptsp,iscf,jdtset,kptopt,lpawu
 integer :: lpawu1,marr,mu,multi_lpawu,multi_mxnsp,multi_natom,multi_natpawu,multi_ncenter
 integer :: multi_nconeq,multi_nkpt,multi_nnos,multi_norb,multi_nshiftk
 integer :: multi_nspinor,multi_nsppol,multi_nsym,multi_ntypat,multi_occopt,mxnsp,natom,natpawu
 integer :: nban,ncenter,nconeq,ndtset_kptopt,nkpt,nkpt_eff,nnos,norb,nshiftk
 integer :: nspinor,nsppol,nsym,ntypat,occopt,prtvol_glob,response,rfelfd,rfmgfd,rfphon,rfstrs
 integer :: rfuser,timopt_default,tnkpt,usepaw
 character(len=2) :: appen
 character(len=500) :: message
!arrays
 integer,allocatable :: intarr(:,:),jdtset_(:),jdtset_kptopt(:),response_(:)
 real(dp) :: rprimd(3,3)
 real(dp),allocatable :: dprarr(:,:),xangst(:,:),xangst_(:,:,:),xcart(:,:)
 real(dp),allocatable :: xcart_(:,:,:),xred(:,:)

! *************************************************************************

!DEBUG
!write(6,*)' outvars : enter '
!ENDDEBUG

!XG 080819 Warning : do not remove the following (silly) line : 
!this is needed to avoid a compiler bug with g95
 allocate(dprarr(1,1)) ; deallocate(dprarr)

!Set up a 'global' prtvol value
 prtvol_glob=1
 if(sum((dtsets(:)%prtvol)**2)==0)prtvol_glob=0

!Echo all variables which either were input or given default values:

 if(choice==1)then
  write(iout, '(a)' )&
&  ' -outvars: echo values of preprocessed input variables --------'
 else
  write(iout, '(a)' )&
&  ' -outvars: echo values of variables after computation  --------'
 end if

 marr=max(3*mxnatom,3*mxnkptgw,mxnkpt*mxnsppol*mxmband,3*mxnkpt,npsp,mxntypat,&
& 9*mxnsym,3*8,3*mxnatom*mxnconeq,3*mxncenter,3*mxnorb,mxnnos,&
& (2*mxlpawu+1)**2*max(mxnsppol,mxnspinor)*mxnatpawu*dmatpuflag)
 allocate(dprarr(marr,0:ndtset_alloc))
 allocate(intarr(marr,0:ndtset_alloc))

!Set up dimensions : determine whether these are different for different
!datasets.
 multi_natom=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%natom/=dtsets(idtset)%natom)multi_natom=1
  end do
 end if
 if(multi_natom==0)natom=dtsets(1)%natom

 multi_ncenter=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%ncenter/=dtsets(idtset)%ncenter)multi_ncenter=1
  end do
 end if
 if(multi_ncenter==0)ncenter=dtsets(1)%ncenter

 multi_nconeq=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%nconeq/=dtsets(idtset)%nconeq)multi_nconeq=1
  end do
 end if
 if(multi_nconeq==0)nconeq=dtsets(1)%nconeq

 multi_norb=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%norb/=dtsets(idtset)%norb)multi_norb=1
  end do
 end if
 if(multi_norb==0)norb=dtsets(1)%norb

 multi_nnos=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%nnos/=dtsets(idtset)%nnos)multi_nnos=1
  end do
 end if
 if(multi_nnos==0)nnos=dtsets(1)%nnos

 multi_nshiftk=0
 nshiftk=1
 if(sum((dtsets(1:ndtset_alloc)%kptopt)**2)/=0)then
  first=0
  do idtset=1,ndtset_alloc
   kptopt=dtsets(idtset)%kptopt
   if(kptopt>=1)then
    if(first==0)then
     first=1
     nshiftk=dtsets(idtset)%nshiftk
    else
     if(nshiftk/=dtsets(idtset)%nshiftk)multi_nshiftk=1
    end if
   end if
  end do
 end if

 multi_nkpt=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%nkpt/=dtsets(idtset)%nkpt)multi_nkpt=1
  end do
 end if
 if(multi_nkpt==0)nkpt=dtsets(1)%nkpt

 multi_nsppol=0;multi_nspinor=0;multi_mxnsp=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%nsppol/=dtsets(idtset)%nsppol)multi_nsppol=1
   if(dtsets(1)%nspinor/=dtsets(idtset)%nspinor)multi_nspinor=1
   if(dtsets(1)%nsppol*dtsets(1)%nspinor/=dtsets(idtset)%nsppol*dtsets(idtset)%nspinor)multi_mxnsp=1
  end do
 end if
 if(multi_nsppol==0)nsppol=dtsets(1)%nsppol
 if(multi_nspinor==0)nspinor=dtsets(1)%nspinor
 if(multi_mxnsp==0)mxnsp=max(nsppol,nspinor)

 multi_nsym=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%nsym/=dtsets(idtset)%nsym)multi_nsym=1
  end do
 end if
 if(multi_nsym==0)nsym=dtsets(1)%nsym

 multi_ntypat=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%ntypat/=dtsets(idtset)%ntypat)multi_ntypat=1
  end do
 end if
 if(multi_ntypat==0)ntypat=dtsets(1)%ntypat

 multi_occopt=0
 if(ndtset_alloc>1)then
  do idtset=1,ndtset_alloc
   if(dtsets(1)%occopt/=dtsets(idtset)%occopt)multi_occopt=1
  end do
 end if
 if(multi_occopt==0)occopt=dtsets(1)%occopt

 response=0
 allocate(response_(ndtset_alloc))
 response_(:)=0
 do idtset=1,ndtset_alloc
  rfelfd=dtsets(idtset)%rfelfd
  rfmgfd=dtsets(idtset)%rfmgfd
  rfphon=dtsets(idtset)%rfphon
  rfstrs=dtsets(idtset)%rfstrs
  rfuser=dtsets(idtset)%rfuser
  if(rfelfd/=0 .or. rfmgfd/=0 .or. rfphon/=0 .or. rfstrs/=0 .or. rfuser/=0)then
   response_(idtset)=1 ; response=1
  end if
 end do

!Must compute xangst and xcart
 allocate(xangst_(3,mxnatom,0:ndtset_alloc),xcart_(3,mxnatom,0:ndtset_alloc))
 xangst_(:,:,:)=0.0_dp ; xcart_(:,:,:)=0.0_dp
 do idtset=1,ndtset_alloc
  natom=dtsets(idtset)%natom
  allocate(xred(3,natom),xangst(3,natom),xcart(3,natom))
  xred(:,1:natom)=results_out(idtset)%xred(:,1:natom)
  rprimd(:,:)     =results_out(idtset)%rprimd(:,:)
! Compute xcart from xred and rprimd
  call xredxcart(natom,1,rprimd,xcart,xred)
! Compute xangst from xcart
  xangst(:,:)=xcart(:,:)*Bohr_Ang
! Save the data
  xangst_(1:3,1:natom,idtset)=xangst(:,:)
  xcart_(1:3,1:natom,idtset)=xcart(:,:)
  deallocate(xred,xangst,xcart)
 end do

 allocate(jdtset_(0:ndtset_alloc))
 jdtset_(0:ndtset_alloc)=dtsets(0:ndtset_alloc)%jdtset

!Determine whether we are in a PAW run
 usepaw=0;if (maxval(dtsets(0:ndtset_alloc)%usepaw)==1) usepaw=1

!Print variables between acell and natom (by alphabetic order)
 call outvar1(choice,dtsets,iout,istatr,istatshft,&
& jdtset_,mxmband,mxnatom,mxnatsph,mxnatvshift,mxnkptgw,mxnkpt,mxnqptdm,mxnsppol,mxnsym,mxntypat,&
& ndtset,ndtset_alloc,npsp,prtvol_glob,pspheads,response,response_,&
& results_out,usepaw)

!DEBUG
!write(6,*)' outvars : after outvar1 '
!stop
!ENDDEBUG

!Print remaining variables, one at a time

!natsph   need to be printed only if there is some occurence of prtdos==3
 do idtset=1,ndtset_alloc
  if(dtsets(idtset)%prtdos==3)then
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01155a)'natsph',appen,dtsets(idtset)%natsph
  end if
 end do

!nband
 if(multi_nkpt==0.and.multi_nsppol==0.and.multi_occopt==0)then
  if(occopt==2)then
   do idtset=0,ndtset_alloc
    intarr(1:nkpt*nsppol,idtset)=&
&    dtsets(idtset)%nband(1:nkpt*nsppol)
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,nkpt*nsppol,&
&   ndtset_alloc,'nband','INT')
  else
   do idtset=0,ndtset_alloc
    intarr(1,idtset)=dtsets(idtset)%nband(1)
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,&
&   ndtset_alloc,'nband','INT')
  end if
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   occopt=dtsets(idtset)%occopt
!  The quantity of data to be output vary with occopt
   if(occopt==2)then
    write(iout,format01155a)&
&    'nband',appen,dtsets(idtset)%nband(1:dtsets(idtset)%nkpt*dtsets(idtset)%nsppol)
   else
    write(iout,format01155a)'nband',appen,dtsets(idtset)%nband(1)
   end if
  end do
 end if

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%nbandsus
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nbandsus','INT')

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%nbdblock
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nbdblock','INT')

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%nbdbuf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nbdbuf','INT')

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%nberry
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nberry','INT')

 intarr(1,:)=dtsets(:)%ncenter
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ncenter','INT')

 intarr(1,:)=dtsets(:)%nconeq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nconeq','INT')

 intarr(1,:)=dtsets(:)%nctime
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nctime','INT')

 intarr(1,:)=dtsets(:)%ndyson
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ndyson','INT')

 if(ndtset>0) write(iout,format01110) 'ndtset',ndtset

 intarr(1,:)=dtsets(:)%nfreqsus
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nfreqsus','INT')

 intarr(1,:)=dtsets(:)%ngeohist
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ngeohist','INT')

 intarr(1,:)=dtsets(:)%ngfft(1)
 intarr(2,:)=dtsets(:)%ngfft(2)
 intarr(3,:)=dtsets(:)%ngfft(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,ndtset_alloc,'ngfft','INT')

 if (usepaw==1) then
  intarr(1,:)=dtsets(:)%ngfftdg(1)
  intarr(2,:)=dtsets(:)%ngfftdg(2)
  intarr(3,:)=dtsets(:)%ngfftdg(3)
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,ndtset_alloc,'ngfftdg','INT')
 end if

 intarr(1,:)=dtsets(:)%ngroup_rf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ngroup_rf','INT')

 intarr(1,:)=dtsets(:)%nkptgw
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nkptgw','INT')

 intarr(1,:)=dtsets(:)%nkpt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nkpt','INT')

 intarr(1,:)=dtsets(:)%nline
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nline','INT')

 intarr(1,:)=dtsets(:)%nloalg(1)+10*dtsets(:)%nloalg(5)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nloalg','INT')

 intarr(1,:)=dtsets(:)%nnsclo
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nnsclo','INT')

 intarr(1,:)=dtsets(:)%nomegasf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nomegasf','INT')

 intarr(1,:)=dtsets(:)%nomegasi
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nomegasi','INT')

 intarr(1,:)=dtsets(:)%nomegasrd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nomegasrd','INT')

 intarr(1,:)=dtsets(:)%norb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'norb','INT')

 intarr(1,:)=dtsets(:)%nnos
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nnos','INT')

 dprarr(1,:)=dtsets(:)%noseinert
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'noseinert','DPR')

 if(multi_ntypat/=0)then
  intarr(1,:)=dtsets(:)%npsp
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'npsp','INT')
 else if(multi_ntypat==0 .and. ntypat/=npsp)then
  intarr(1,:)=dtsets(:)%npsp
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'npsp','INT')
 end if

!intarr(1,:)=dtsets(:)%npspalch
!call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'npspalch','INT')

 intarr(1,:)=dtsets(:)%npweps
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'npweps','INT')

 intarr(1,:)=dtsets(:)%npwkss
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'npwkss','INT')

 intarr(1,:)=dtsets(:)%npwsigx
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'npwsigx','INT')

 intarr(1,:)=dtsets(:)%npwwfn
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'npwwfn','INT')

 intarr(1,:)=dtsets(:)%nqpt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nqpt','INT')

 intarr(1,:)=dtsets(:)%nqptdm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nqptdm','INT')

 intarr(1,:)=dtsets(:)%nscforder
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nscforder','INT')

 intarr(1,:)=dtsets(:)%nsheps
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nsheps','INT')

!nshiftk
 if(sum((dtsets(1:ndtset_alloc)%kptopt)**2)/=0)then
  ndtset_kptopt=0
  intarr(1:1,0)=dtsets(0)%nshiftk
  allocate(jdtset_kptopt(0:ndtset_alloc))
! Define the set of datasets for which kptopt>0
  do idtset=1,ndtset_alloc
   kptopt=dtsets(idtset)%kptopt
   if(kptopt>0)then
    ndtset_kptopt=ndtset_kptopt+1
    jdtset_kptopt(ndtset_kptopt)=jdtset_(idtset)
    intarr(1:1,ndtset_kptopt)=dtsets(idtset)%nshiftk
   end if
  end do
  if(ndtset_kptopt>0)then
   call prttagm(dprarr,intarr,iout,jdtset_kptopt,2,marr,1,&
&   ndtset_kptopt,'nshiftk','INT')
  end if
  deallocate(jdtset_kptopt)
 end if

 intarr(1,:)=dtsets(:)%nshsigx
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nshsigx','INT')

 intarr(1,:)=dtsets(:)%nshwfn
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nshwfn','INT')

 intarr(1,:)=dtsets(:)%nspden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nspden','INT')

 intarr(1,:)=dtsets(:)%nspinor
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nspinor','INT')

 intarr(1,:)=dtsets(:)%nsppol
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nsppol','INT')

 intarr(1,:)=dtsets(:)%nstep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nstep','INT')

 intarr(1,:)=dtsets(:)%nsym
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nsym','INT')

 intarr(1,:)=dtsets(:)%ntime
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ntime','INT')

 intarr(1,:)=dtsets(:)%ntypalch
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ntypalch','INT')

!intarr(1,:)=dtsets(:)%ntyppure
!call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ntyppure','INT')

 intarr(1,:)=dtsets(:)%ntypat
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ntypat','INT')

!numorb
 if(multi_ncenter==0)then
  do idtset=0,ndtset_alloc
   intarr(1:ncenter,idtset)=dtsets(idtset)%numorb(1:ncenter)
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,ncenter,&
&  ndtset_alloc,'numorb','INT')
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,'(1x,a9,a,1x,(t13,20i3))')&
&   'numorb',appen,dtsets(idtset)%numorb(1:dtsets(idtset)%ncenter)
  end do
 end if

 intarr(1,:)=dtsets(:)%nwfshist
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nwfshist','INT')

!occ
!The use of prttagm for occ if occopt>=2 is not possible because
!the different k-point and spins must be separated on different lines ...
 call prtocc(dtsets,iout,jdtset_,mxmband,mxnkpt,mxnsppol,&
& ndtset_alloc,prtvol_glob,results_out)

 intarr(1,:)=dtsets(:)%occopt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'occopt','INT')

 dprarr(1,:)=dtsets(:)%omegasimax
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'omegasimax','ENE')

 dprarr(1,:)=dtsets(:)%omegasrdmax
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'omegasrdmax','ENE')

 intarr(1,:)=dtsets(:)%optcell
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'optcell','INT')

 intarr(1,:)=dtsets(:)%optdriver
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'optdriver','INT')

 intarr(1,:)=dtsets(:)%optforces
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'optforces','INT')

 intarr(1,:)=dtsets(:)%optstress
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'optstress','INT')

 intarr(1,:)=dtsets(:)%optnlxccc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'optnlxccc','INT')

 intarr(1,:)=dtsets(:)%ortalg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ortalg','INT')

 intarr(1,:)=dtsets(:)%outputxml
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'outputxml','INT')

 intarr(1,:)=dtsets(:)%paral_rf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'paral_rf','INT')

 intarr(1,:)=dtsets(:)%paral_kgb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'paral_kgb','INT')


 if (usepaw==1) then

  dprarr(1,:)=dtsets(:)%pawecutdg
  call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'pawecutdg','ENE')

  intarr(1,:)=dtsets(:)%pawlcutd
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawlcutd','INT')

  intarr(1,:)=dtsets(:)%pawlmix
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawlmix','INT')

  intarr(1,:)=dtsets(:)%pawmixdg
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawmixdg','INT')

  intarr(1,:)=dtsets(:)%pawnhatxc
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawnhatxc','INT')

  intarr(1,:)=dtsets(:)%pawnphi
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawnphi','INT')

  intarr(1,:)=dtsets(:)%pawntheta
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawntheta','INT')

  intarr(1,:)=dtsets(:)%pawnzlm
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawnzlm','INT')

  intarr(1,:)=dtsets(:)%pawoptmix
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawoptmix','INT')

  dprarr(1,:)=dtsets(:)%pawovlp
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawovlp','DPR')

  intarr(1,:)=dtsets(:)%pawprtdos
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawprtdos','INT')

  intarr(1,:)=dtsets(:)%pawprtvol
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawprtvol','INT')

  dprarr(1,:)=dtsets(:)%pawsphmix
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawsphmix','DPR')

  intarr(1,:)=dtsets(:)%pawspnorb
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawspnorb','INT')

  intarr(1,:)=dtsets(:)%pawstgylm
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawstgylm','INT')

  intarr(1,:)=dtsets(:)%pawusecp
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawusecp','INT')
  
  intarr(1,:)=dtsets(:)%pawxcdev
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawxcdev','INT')

  icount=0
  do idtset=1,ndtset_alloc
   if (dtsets(idtset)%usepawu>0) icount=icount+1
  end do
  if(icount>0) then
   intarr(1,:)=dtsets(:)%dmatudiag
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'dmatudiag','INT')
   intarr(1,:)=dtsets(:)%dmatpuopt
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'dmatpuopt','INT')
   intarr(1,:)=dtsets(:)%usedmatpu
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'usedmatpu','INT')
   intarr(1,:)=dtsets(:)%usepawu
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'usepawu','INT')
   if(multi_ntypat==0)then
    do idtset=0,ndtset_alloc
     dprarr(1:ntypat,idtset)=dtsets(idtset)%jpawu(1:ntypat)
    end do
    call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,ndtset_alloc,'jpawu','DPR')
    do idtset=0,ndtset_alloc
     intarr(1:ntypat,idtset)=dtsets(idtset)%lpawu(1:ntypat)
    end do
    call prttagm(dprarr,intarr,iout,jdtset_,2,marr,ntypat,ndtset_alloc,'lpawu','INT')
    do idtset=0,ndtset_alloc
     dprarr(1:ntypat,idtset)=dtsets(idtset)%upawu(1:ntypat)
    end do
    call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,ndtset_alloc,'upawu','DPR')
   else
    do idtset=1,ndtset_alloc
     jdtset=jdtset_(idtset)
     if(jdtset<10)write(appen,'(i1)')jdtset
     if(jdtset>=10)write(appen,'(i2)')jdtset
     write(iout,format01160a)'jpawu',appen,dtsets(idtset)%jpawu(1:dtsets(idtset)%ntypat)
     write(iout,format01155a)'lpawu',appen,dtsets(idtset)%lpawu(1:dtsets(idtset)%ntypat)
     write(iout,format01160a)'upawu',appen,dtsets(idtset)%upawu(1:dtsets(idtset)%ntypat)
    end do
   end if
   if (dmatpuflag==1.and.mxnatpawu>0) then
    multi_lpawu=0;multi_natpawu=0
    lpawu=maxval(dtsets(1)%lpawu(:))
    natpawu=dtsets(1)%natpawu
    if(ndtset_alloc>1)then
     do idtset=1,ndtset_alloc
      lpawu1=maxval(dtsets(idtset)%lpawu(:))
      if(lpawu/=lpawu1) multi_lpawu=1
      if(dtsets(idtset)%natpawu/=natpawu) multi_natpawu=1
     end do
    end if
    if(multi_lpawu==0.and.multi_mxnsp==0.and.multi_natpawu==0)then
     ii=(2*lpawu+1)*(2*lpawu+1)*mxnsp*natpawu
     do idtset=0,ndtset_alloc
      if (lpawu1/=-1) dprarr(1:ii,idtset)=&
&      reshape(dtsets(idtset)%dmatpawu(1:2*lpawu+1,1:2*lpawu+1,1:mxnsp,1:natpawu),(/ii/))
     end do
     call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ii,ndtset_alloc,'dmatpawu','DPR')
    else
     do idtset=1,ndtset_alloc
      lpawu1=maxval(dtsets(idtset)%lpawu)
      if (dtsets(idtset)%usedmatpu/=0.and.lpawu1/=-1) then
       jdtset=jdtset_(idtset)
       if(jdtset<10)write(appen,'(i1)')jdtset
       if(jdtset>=10)write(appen,'(i2)')jdtset
       do iat=1,dtsets(idtset)%natpawu
        do mu=1,dtsets(idtset)%nsppol*dtsets(idtset)%nspinor
         if (iat==1.and.mu==1) then
          write(iout,'(1x,a9,a,es12.4,8(1x,es12.4))') &
&          'dmatpawu',appen,dtsets(idtset)%dmatpawu(1:2*lpawu1+1,1,mu,iat)
         else
          write(iout,'(t13,es12.4,8(1x,es12.4))') dtsets(idtset)%dmatpawu(1:2*lpawu1+1,1,mu,iat)
         end if
         if (lpawu1>0) then
          do ii=2,2*lpawu1+1
           write(iout,'(t13,es12.4,8(1x,es12.4))') dtsets(idtset)%dmatpawu(1:2*lpawu1+1,ii,mu,iat)
          end do
         end if
        end do
       end do
      end if
     end do
    end if
   end if
  end if
  icount=0
  do idtset=1,ndtset_alloc
   if (dtsets(idtset)%useexexch>0) icount=icount+1
  end do
  if(icount>0) then
   intarr(1,:)=dtsets(:)%useexexch
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'useexexch','INT')
   if(multi_ntypat==0)then
    do idtset=0,ndtset_alloc
     intarr(1:ntypat,idtset)=dtsets(idtset)%lexexch(1:ntypat)
    end do
    call prttagm(dprarr,intarr,iout,jdtset_,2,marr,ntypat,ndtset_alloc,'lexexch','INT')
   else
    do idtset=1,ndtset_alloc
     jdtset=jdtset_(idtset)
     if(jdtset<10)write(appen,'(i1)')jdtset
     if(jdtset>=10)write(appen,'(i2)')jdtset
     write(iout,format01155a)'lexexch',appen,dtsets(idtset)%lexexch(1:dtsets(idtset)%ntypat)
    end do
   end if
  end if

! intarr(1,:)=dtsets(:)%ngfftdg(1)
! intarr(2,:)=dtsets(:)%ngfftdg(2)
! intarr(3,:)=dtsets(:)%ngfftdg(3)
! call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,ndtset_alloc,'ngfftdg','INT')

 end if

 intarr(1,:)=dtsets(:)%positron
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'positron','INT')


 dprarr(1,:)=dtsets(:)%ppmfrq
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'ppmfrq','ENE')

 intarr(1,:)=dtsets(:)%ppmodel
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ppmodel','INT')

 intarr(1,:)=dtsets(:)%gwcomp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gwcomp','INT')

 dprarr(1,:)=dtsets(:)%gwencomp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gwencomp','ENE')

 intarr(1,:)=dtsets(:)%gwgamma
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gwgamma','INT')

 intarr(1,:)=dtsets(:)%gwmem
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gwmem','INT')

 intarr(1,:)=dtsets(:)%gwpara
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gwpara','INT')


 intarr(1,:)=dtsets(:)%nfreqmidm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'nfreqmidm','INT')
 
 do idtset=0, ndtset_alloc
  do ii = 1, ntypat
   dprarr(ii,idtset) = dtsets(idtset)%ptcharge(ii)
  end do ! end loop over ntypat
 end do ! end loop over datasets
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,ndtset_alloc,'ptcharge','DPR')


 intarr(1,:)=dtsets(:)%prepanl
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prepanl','INT')

 intarr(1,:)=dtsets(:)%prepgkk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prepgkk','INT')

 intarr(1,:)=dtsets(:)%prtbbb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtbbb','INT')

 intarr(1,:)=dtsets(:)%prtcml
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtcml','INT')

 intarr(1,:)=dtsets(:)%prtcs
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtcs','INT')

 intarr(1,:)=dtsets(:)%prtden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtden','INT')

 intarr(1,:)=dtsets(:)%prtdensph
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtdensph','INT')

 intarr(1,:)=dtsets(:)%prtdos
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtdos','INT')

 intarr(1,:)=dtsets(:)%prtdosm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtdosm','INT')

 intarr(1,:)=dtsets(:)%prtefg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtefg','INT')

 intarr(1,:)=dtsets(:)%prteig
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prteig','INT')
 
 intarr(1,:)=dtsets(:)%prtfc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtfc','INT')

 intarr(1,:)=dtsets(:)%prtfsurf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtfsurf','INT')

 intarr(1,:)=dtsets(:)%prtgeo
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtgeo','INT')

 intarr(1,:)=dtsets(:)%prtgkk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtgkk','INT')

 intarr(1,:)=dtsets(:)%prtnabla
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtnabla','INT')

 intarr(1,:)=dtsets(:)%prtpot
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtpot','INT')

 intarr(1,:)=dtsets(:)%prtpmp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtpmp','INT')

 intarr(1,:)=dtsets(:)%prtspcur
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtspcur','INT')

 intarr(1,:)=dtsets(:)%prtstm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtstm','INT')

 intarr(1,:)=dtsets(:)%prtvha
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtvha','INT')

 intarr(1,:)=dtsets(:)%prtvhxc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtvhxc','INT')

 intarr(1,:)=dtsets(:)%prtvxc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtvxc','INT')

 intarr(1,:)=dtsets(:)%prtvol
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtvol','INT')

 intarr(1,:)=dtsets(:)%prtwant
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtwant','INT')

 intarr(1,:)=dtsets(:)%prtwf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prtwf','INT')

 intarr(1,:)=dtsets(:)%prt1dm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'prt1dm','INT')

 intarr(1,:)=dtsets(:)%ptgroupma
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ptgroupma','INT')

 intarr(1,:)=dtsets(:)%qprtrb(1)
 intarr(2,:)=dtsets(:)%qprtrb(2)
 intarr(3,:)=dtsets(:)%qprtrb(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,ndtset_alloc,'qprtrb','INT')

 dprarr(1,:)=dtsets(:)%qpt(1)
 dprarr(2,:)=dtsets(:)%qpt(2)
 dprarr(3,:)=dtsets(:)%qpt(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,ndtset_alloc,'qpt','DPR')

 do idtset=0, ndtset_alloc
  do ii = 1, ntypat
   dprarr(ii,idtset) = dtsets(idtset)%quadmom(ii)
  end do ! end loop over ntypat
 end do ! end loop over datasets
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,ndtset_alloc,'quadmom','DPR')

 dprarr(1,:)=dtsets(:)%qptnrm
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'qptnrm','DPR')

!Variables used for recursion method
 intarr(1,:)=dtsets(:)%recnpath
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'recnpath','INT')
 intarr(1,:)=dtsets(:)%recnrec
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'recnrec','INT')
 intarr(1,:)=dtsets(:)%recptrott
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'recptrott','INT')
 intarr(1,:)=dtsets(:)%rectesteg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rectesteg','INT')
 dprarr(1,:)=dtsets(:)%recrcut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'recrcut','LEN')
 dprarr(1,:)=dtsets(:)%recefermi
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'recefermi','ENE')
 dprarr(1,:)=dtsets(:)%rectolden
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'rectolden','DPR')

 do idtset=0, ndtset_alloc
  do ii = 1, ntypat
   dprarr(ii,idtset) = dtsets(idtset)%ratsph(ii)
  end do ! end loop over ntypat
 end do ! end loop over datasets
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,ndtset_alloc,'ratsph','LEN')

!rcoord
 if(multi_ncenter==0)then
! Required if for pathscale to avoid failure with -ff_bounds_check
  if (ncenter > 0) then
   do idtset=0,ndtset_alloc
    intarr(1:3*ncenter,idtset)=reshape(dtsets(idtset)%rcoord(1:3,1:ncenter),(/3*ncenter/))
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,3,marr,3*ncenter,ndtset_alloc,'rcoord','DPR')
  end if
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01150a)&
&   'rcoord',appen,dtsets(idtset)%rcoord(1:3,1:dtsets(idtset)%ncenter)
  end do
 end if

 dprarr(1,:)=dtsets(:)%rcut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'rcut','LEN')

 intarr(1,:)=dtsets(:)%rdmnb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rdmnb','INT')
 intarr(1,:)=dtsets(:)%restartxf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'restartxf','INT')

 if(response==1)then

  intarr(1,:)=dtsets(:)%rfasr
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rfasr','INT')

  intarr(1,:)=dtsets(:)%rfatpol(1)
  intarr(2,:)=dtsets(:)%rfatpol(2)
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,ndtset_alloc,'rfatpol','INT')

  intarr(1,:)=dtsets(:)%rfdir(1)
  intarr(2,:)=dtsets(:)%rfdir(2)
  intarr(3,:)=dtsets(:)%rfdir(3)
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,ndtset_alloc,'rfdir','INT')

  intarr(1,:)=dtsets(:)%rfelfd
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rfelfd','INT')

  intarr(1,:)=dtsets(:)%rfmgfd
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rfmgfd','INT')

  intarr(1,:)=dtsets(:)%rfmeth
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rfmeth','INT')

  intarr(1,:)=dtsets(:)%rfphon
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rfphon','INT')

  intarr(1,:)=dtsets(:)%rfstrs
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rfstrs','INT')

  intarr(1,:)=dtsets(:)%rfthrd
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rfthrd','INT')

  intarr(1,:)=dtsets(:)%rfuser
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rfuser','INT')

 end if

 intarr(1,:)=dtsets(:)%rf1atpol(1)
 intarr(2,:)=dtsets(:)%rf1atpol(2)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,ndtset_alloc,'rf1atpol','INT')
 intarr(1,:)=dtsets(:)%rf1dir(1)
 intarr(2,:)=dtsets(:)%rf1dir(2)
 intarr(3,:)=dtsets(:)%rf1dir(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,ndtset_alloc,'rf1dir','INT')
 intarr(1,:)=dtsets(:)%rf1elfd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rf1elfd','INT')
 intarr(1,:)=dtsets(:)%rf1phon
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rf1phon','INT')

 intarr(1,:)=dtsets(:)%rf2atpol(1)
 intarr(2,:)=dtsets(:)%rf2atpol(2)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,ndtset_alloc,'rf2atpol','INT')
 intarr(1,:)=dtsets(:)%rf2dir(1)
 intarr(2,:)=dtsets(:)%rf2dir(2)
 intarr(3,:)=dtsets(:)%rf2dir(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,ndtset_alloc,'rf2dir','INT')
 intarr(1,:)=dtsets(:)%rf2elfd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rf2elfd','INT')
 intarr(1,:)=dtsets(:)%rf2phon
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rf2phon','INT')

 intarr(1,:)=dtsets(:)%rf3atpol(1)
 intarr(2,:)=dtsets(:)%rf3atpol(2)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,ndtset_alloc,'rf3atpol','INT')
 intarr(1,:)=dtsets(:)%rf3dir(1)
 intarr(2,:)=dtsets(:)%rf3dir(2)
 intarr(3,:)=dtsets(:)%rf3dir(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,ndtset_alloc,'rf3dir','INT')
 intarr(1,:)=dtsets(:)%rf3elfd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rf3elfd','INT')
 intarr(1,:)=dtsets(:)%rf3phon
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'rf3phon','INT')

 do idtset=0,ndtset_alloc
  dprarr(1:9,idtset)= reshape(results_out(idtset)%rprim(:,:),(/9/))
  do ii=1,9
   if(abs(dprarr(ii,idtset))<tol12)dprarr(ii,idtset)=zero  ! This is to improve the portability
  end do
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,9,ndtset_alloc,'rprim','DPR')

!rtheta
 if(multi_norb==0)then
! Required if for pathscale to avoid failure with -ff_bounds_check
  if (norb > 0) then
   do idtset=0,ndtset_alloc
    dprarr(1:3*norb,idtset)=reshape(dtsets(idtset)%rtheta(1:3,1:norb),(/3*norb/))
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,3,marr,3*norb,&
&   ndtset_alloc,'rtheta','DPR')
  end if
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01150a)&
&   'rtheta',appen,dtsets(idtset)%rtheta(1:3,1:dtsets(idtset)%norb)
  end do
 end if

!qmass
 if(multi_nnos==0)then
! Required if for pathscale to avoid failure with -ff_bounds_check
  if (nnos > 0) then
   do idtset=0,ndtset_alloc
    dprarr(1:nnos,idtset)=reshape(dtsets(idtset)%qmass(1:nnos),(/nnos/))
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,nnos,&
&   ndtset_alloc,'qmass','DPR')
  end if
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01150a)&
&   'qmass',appen,dtsets(idtset)%qmass(1:dtsets(idtset)%nnos)
  end do
 end if

 if(response==1)then
  dprarr(1,:)=dtsets(:)%sciss
  call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'sciss','ENE')
 end if

!shiftk (printed only when kptopt>0)
 if(sum((dtsets(1:ndtset_alloc)%kptopt)**2)/=0)then
  if(multi_nshiftk==0)then
   allocate(jdtset_kptopt(0:ndtset_alloc))
   ndtset_kptopt=0
   dprarr(1:3*nshiftk,:)=0.0_dp
   do idtset=1,ndtset_alloc
    kptopt=dtsets(idtset)%kptopt
    if(kptopt>0)then
     ndtset_kptopt=ndtset_kptopt+1
     jdtset_kptopt(ndtset_kptopt)=jdtset_(idtset)
     dprarr(1:3*nshiftk,ndtset_kptopt)=&
&     reshape(dtsets(idtset)%shiftk(1:3,1:nshiftk),(/3*nshiftk/) )
    end if
   end do
   if(ndtset_kptopt>0)then
    call prttagm(dprarr,intarr,iout,jdtset_kptopt,1,marr,3*nshiftk,&
&    ndtset_kptopt,'shiftk','DPR')
   end if
   deallocate(jdtset_kptopt)
  else
   do idtset=1,ndtset_alloc
    if(dtsets(idtset)%kptopt>0)then
     jdtset=jdtset_(idtset)
     if(jdtset<10)write(appen,'(i1)')jdtset
     if(jdtset>=10)write(appen,'(i2)')jdtset
     write(iout,format01150a)&
&     'shiftk',appen,dtsets(idtset)%shiftk(1:3,1:dtsets(idtset)%nshiftk)
    end if
   end do
  end if
! End of test to see whether kptopt/=0 for some dataset
 end if

 intarr(1,:)=dtsets(:)%signperm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'signperm','INT')

 dprarr(1,:)=dtsets(:)%slabwsrad
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'slabwsrad','DPR')

 dprarr(1,:)=dtsets(:)%slabzbeg
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'slabzbeg','DPR')

 dprarr(1,:)=dtsets(:)%slabzend
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'slabzend','DPR')

 dprarr(1,:)=dtsets(:)%soenergy
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'soenergy','ENE')

 dprarr(1,:)=dtsets(:)%spbroad
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'spbroad','ENE')

 do idtset=0,ndtset_alloc
  intarr(1:npsp,idtset)=dtsets(idtset)%so_psp(1:npsp)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,npsp,ndtset_alloc,'so_psp','INT')

 intarr(1,:)=dtsets(:)%spgroup
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'spgroup','INT')

!spinat
 if(multi_natom==0)then
  do idtset=0,ndtset_alloc
   dprarr(1:3*natom,idtset)= &
&   reshape( dtsets(idtset)%spinat(1:3,1:natom) , (/3*natom/) )
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3*natom,&
&  ndtset_alloc,'spinat','DPR')
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   if(sum(abs( dtsets(idtset)%spinat(1:3,1:dtsets(idtset)%natom)- &
&   dtsets(0)%spinat(1:3,1:dtsets(idtset)%natom)      )) > tol12 )then
    write(iout,format01150a)'spinat',appen,dtsets(idtset)%spinat(1:3,1:dtsets(idtset)%natom)
   end if
  end do
 end if

 intarr(1,:)=dtsets(:)%splitsigc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'splitsigc',&
& 'INT')

 intarr(1,:)=dtsets(:)%spmeth
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'spmeth','INT')

 dprarr(1,:)=dtsets(:)%stmbias
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'stmbias','DPR')

 dprarr(1,:)=dtsets(:)%strfact
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'strfact','DPR')

 do ii=1,6
  dprarr(ii,:)=dtsets(:)%strtarget(ii)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,6,ndtset_alloc,'strtarget','DPR')

!strten
 if(choice==2)then
  do idtset=1,ndtset_alloc
   iscf=dtsets(idtset)%iscf
   if(iscf>0)then
    if(ndtset>0)then
     jdtset=jdtset_(idtset)
     if(jdtset<10)write(appen,'(i1)')jdtset
     if(jdtset>=10)write(appen,'(i2)')jdtset
     write(iout,format01160a)'strten',appen,results_out(idtset)%strten(:)
    else
     write(iout,format01160)'strten',results_out(idtset)%strten(:)
    end if
   end if
  end do
 end if

 intarr(1,:)=dtsets(:)%supercell(1)
 intarr(2,:)=dtsets(:)%supercell(2)
 intarr(3,:)=dtsets(:)%supercell(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,ndtset_alloc,'supercell','INT')

!symafm
 if(multi_nsym==0)then
  do idtset=0,ndtset_alloc
   intarr(1:nsym,idtset)=dtsets(idtset)%symafm(1:nsym)
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,3,marr,nsym,ndtset_alloc,'symafm','INT')
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01155a)'symafm',appen,dtsets(idtset)%symafm(1:dtsets(idtset)%nsym)
  end do
 end if

 intarr(1,:)=dtsets(:)%symmorphi
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'symmorphi','INT')

!symrel
 if(multi_nsym==0)then
  do idtset=0,ndtset_alloc
   intarr(1:3*3*nsym,idtset)=&
&   reshape(dtsets(idtset)%symrel(1:3,1:3,1:nsym),(/3*3*nsym/) )
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,3,marr,3*3*nsym,&
&  ndtset_alloc,'symrel','INT')
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout, '(1x,a9,a,1x,(t13,3(3i3,1x),4x,3(3i3,1x)))' )&
&   'symrel',appen,dtsets(idtset)%symrel(1:3,1:3,1:dtsets(idtset)%nsym)
  end do
 end if

 intarr(1,:)=dtsets(:)%symsigma
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'symsigma','INT')

 intarr(1,:)=dtsets(:)%symchi
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'symchi','INT')

 dprarr(1,:)=dtsets(:)%td_maxene
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'td_maxene','DPR')

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%td_mexcit
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'td_mexcit','INT')

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%tfkinfunc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'tfkinfunc','INT')

 dprarr(1,:)=dtsets(:)%tfnewton
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'tfnewton','DPR')

 timopt_default=1
!MPI parallel case
 if(mpi_enreg%paral_compil==1)then
  timopt_default=0
 end if
 if(timopt/=timopt_default)write(iout,format01110) 'timopt',timopt

!WVL - tails related variables
 intarr(1,:)=dtsets(:)%tl_nprccg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'tl_nprccg','INT')
 dprarr(1,:)=dtsets(:)%tl_radius
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'tl_radius','DPR')

!tnons
 if(multi_nsym==0)then
  do idtset=0,ndtset_alloc
   dprarr(1:3*nsym,idtset)=&
&   reshape( dtsets(idtset)%tnons(1:3,1:nsym) , (/3*nsym/) )
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,-3,marr,3*nsym,&
&  ndtset_alloc,'tnons','DPR')
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout, '(1x,a9,a,1x,(t13,3f11.7,3x,3f11.7))' )&
&   'tnons',appen,dtsets(idtset)%tnons(1:3,1:dtsets(idtset)%nsym)
  end do
 end if

 dprarr(1,:)=dtsets(:)%toldfe
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'toldfe','ENE')

 dprarr(1,:)=dtsets(:)%toldff
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'toldff','DPR')

 dprarr(1,:)=dtsets(:)%tolrff
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'tolrff','DPR')

 dprarr(1,:)=dtsets(:)%tolmxf
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'tolmxf','DPR')

 dprarr(1,:)=dtsets(:)%tolvrs
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'tolvrs','DPR')

 dprarr(1,:)=dtsets(:)%tolwfr
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'tolwfr','DPR')

 dprarr(1,:)=dtsets(:)%tphysel
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'tphysel','ENE')

 dprarr(1,:)=dtsets(:)%tsmear
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'tsmear','ENE')

!typat
 if(multi_natom==0)then
  do idtset=0,ndtset_alloc
   intarr(1:natom,idtset)=dtsets(idtset)%typat(1:natom)
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,4,marr,natom,&
&  ndtset_alloc,'typat','INT')
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,'(1x,a9,a,1x,(t13,20i3))')&
&   'typat',appen,dtsets(idtset)%typat(1:dtsets(idtset)%natom)
  end do
 end if

 intarr(1,:)=dtsets(:)%useria
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'useria','INT')

 intarr(1,:)=dtsets(:)%userib
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'userib','INT')

 intarr(1,:)=dtsets(:)%useric
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'useric','INT')

 intarr(1,:)=dtsets(:)%userid
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'userid','INT')

 intarr(1,:)=dtsets(:)%userie
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'userie','INT')

 dprarr(1,:)=dtsets(:)%userra
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'userra','DPR')

 dprarr(1,:)=dtsets(:)%userrb
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'userrb','DPR')

 dprarr(1,:)=dtsets(:)%userrc
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'userrc','DPR')

 dprarr(1,:)=dtsets(:)%userrd
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'userrd','DPR')

 dprarr(1,:)=dtsets(:)%userre
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'userre','DPR')

 intarr(1,:)=dtsets(:)%usewvl
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'usewvl','INT')

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%useylm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'useylm','INT')

 dprarr(1,:)=dtsets(:)%vcutgeo(1)
 dprarr(2,:)=dtsets(:)%vcutgeo(2)
 dprarr(3,:)=dtsets(:)%vcutgeo(3)
 call prttagm(dprarr,intarr,iout,jdtset_,3,marr,3,ndtset_alloc,'vcutgeo','DPR')

!vel
 if(multi_natom==0)then
  do idtset=0,ndtset_alloc
   dprarr(1:3*natom,idtset)=reshape(results_out(idtset)%vel(:,1:natom),(/3*natom/) )
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3*natom,&
&  ndtset_alloc,'vel','DPR')
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   if(sum(abs( results_out(idtset)%vel(:,1:dtsets(idtset)%natom)- &
&   results_out(0)     %vel(:,1:dtsets(idtset)%natom)      )) > tol12 )then
    write(iout,format01160a)'vel',appen,results_out(idtset)%vel(:,1:dtsets(idtset)%natom)
   end if
  end do
 end if

 dprarr(1,:)=dtsets(:)%vis
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'vis','DPR')

 dprarr(1,:)=dtsets(:)%vprtrb(1)
 dprarr(2,:)=dtsets(:)%vprtrb(2)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,2,ndtset_alloc,'vprtrb','ENE')

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%wfoptalg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'wfoptalg','INT')

!wtatcon
 if(multi_natom==0 .and. multi_nconeq==0)then
! Required if for pathscale to avoid failure with -ff_bounds_check
  if (nconeq > 0) then
   do idtset=0,ndtset_alloc
    dprarr(1:3*natom*nconeq,idtset)=&
&    reshape(dtsets(idtset)%wtatcon(1:3,1:natom,1:nconeq),(/3*natom*nconeq/) )
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,-2,marr,3*natom*nconeq,&
&   ndtset_alloc,'wtatcon','DPR')
  end if
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   nconeq=dtsets(idtset)%nconeq
   if(nconeq>0)then
    write(iout,format01160a)&
&    'wtatcon',appen,dtsets(idtset)%wtatcon(1:3,1:dtsets(idtset)%natom,1:nconeq)
   end if
  end do
 end if

!wtk
 if(multi_nkpt==0)then
! Might restrict the number of k points to be printed
  tnkpt=0
  nkpt_eff=nkpt
  if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
   nkpt_eff=nkpt_max
   tnkpt=1
  end if
! Modify slightly the value of wtk, thanks to tol12, to improve portability
  do idtset=0,ndtset_alloc
   dprarr(1:nkpt_eff,idtset)=dtsets(idtset)%wtk(1:nkpt_eff)+tol12
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,4,marr,nkpt_eff,&
&  ndtset_alloc,'wtk','DPR')
  if(tnkpt==1)write(iout,'(16x,a)' )&
&  'outvar1 : prtvol=0, do not print more k-points.'
 else
  do idtset=1,ndtset_alloc
   tnkpt=0
   nkpt_eff=dtsets(idtset)%nkpt
   if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
    nkpt_eff=nkpt_max
    tnkpt=1
   end if
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
!  Modify slightly the value of wtk, thanks to tol12, to improve portability
   write(iout,'(1x,a9,a,1x,(t13,6f11.5))')&
&   'wtk',appen,dtsets(idtset)%wtk(1:nkpt_eff)+tol12
   if(tnkpt==1)write(iout,'(16x,a)' )&
&   'outvar1 : prtvol=0, do not print more k-points.'
  end do
 end if

!WVL - wavelets variables
 dprarr(1,:)=dtsets(:)%wvl_crmult
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'wvl_crmult','DPR')
 dprarr(1,:)=dtsets(:)%wvl_frmult
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'wvl_frmult','DPR')
 dprarr(1,:)=dtsets(:)%wvl_cpmult
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'wvl_cpmult','DPR')
 dprarr(1,:)=dtsets(:)%wvl_fpmult
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'wvl_fpmult','DPR')
 dprarr(1,:)=dtsets(:)%wvl_hgrid
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'wvl_hgrid','DPR')
 intarr(1,:)=dtsets(:)%wvl_nprccg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'wvl_nprccg','INT')

!Wannier90 interface related variables
 if(sum(dtsets(1:ndtset_alloc)%prtwant) >1)then
  intarr(1,:)=dtsets(:)%w90iniprj
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'w90iniprj','INT')
  intarr(1,:)=dtsets(:)%w90prtunk
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'w90prtunk','INT')
  intarr(1,:)=dtsets(:)%w90nplot
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'w90nplot','INT')
  do idtset=0,ndtset_alloc
   intarr(1:3,idtset)=reshape(dtsets(idtset)%w90cplot(:),(/3/))
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,ndtset_alloc,'w90cplot','INT')
! w90lplot
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   if(dtsets(idtset)%w90nplot>0) then
    write(iout,'(1x,a9,a,1x,(t13,20i3))')&
&    'w90lplot',appen,dtsets(idtset)%w90lplot(1:dtsets(idtset)%w90nplot)
   end if
  end do !idtset
 end if !prtwant>1

!xangst
 if(multi_natom==0)then
  dprarr(1:3*natom,0:ndtset_alloc)=&
&  reshape(xangst_(1:3,1:natom,0:ndtset_alloc),(/3*natom,ndtset_alloc+1/) )
  call prttagm(dprarr,intarr,iout,jdtset_,-2,marr,3*natom,&
&  ndtset_alloc,'xangst','DPR')
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01160a)'xangst',appen,xangst_(1:3,1:dtsets(idtset)%natom,idtset)
  end do
 end if

!xcart
 if(multi_natom==0)then
  dprarr(1:3*natom,0:ndtset_alloc)=&
&  reshape(xcart_(1:3,1:natom,0:ndtset_alloc),(/3*natom,ndtset_alloc+1/) )
  call prttagm(dprarr,intarr,iout,jdtset_,-2,marr,3*natom,&
&  ndtset_alloc,'xcart','DPR')
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01160a)'xcart',appen,xcart_(1:3,1:dtsets(idtset)%natom,idtset)
  end do
 end if

!xred
 if(multi_natom==0)then
  do idtset=0,ndtset_alloc
   dprarr(1:3*natom,idtset)=&
&   reshape(results_out(idtset)%xred(:,1:natom),(/3*natom/) )
  end do
  call prttagm(dprarr,intarr,iout,jdtset_,-2,marr,3*natom,&
&  ndtset_alloc,'xred','DPR')
 else
  do idtset=1,ndtset_alloc
   jdtset=jdtset_(idtset)
   if(jdtset<10)write(appen,'(i1)')jdtset
   if(jdtset>=10)write(appen,'(i2)')jdtset
   write(iout,format01160a)'xred',appen,results_out(idtset)%xred(:,1:dtsets(idtset)%natom)
  end do
 end if

 dprarr(1,:)=dtsets(:)%zcut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'zcut','ENE')

!ziontypat   ! After all, should always echo this value
 if(sum(dtsets(:)%ntypalch)>0)then   ! After all, should always echo this value ...
  if(multi_ntypat==0)then
   do idtset=0,ndtset_alloc
    dprarr(1:ntypat,idtset)=dtsets(idtset)%ziontypat(1:ntypat)
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,ndtset_alloc,'ziontypat','DPR')
  else
   do idtset=1,ndtset_alloc
    jdtset=jdtset_(idtset)
    if(jdtset<10)write(appen,'(i1)')jdtset
    if(jdtset>=10)write(appen,'(i2)')jdtset
    write(iout,format01160a)'ziontypat',appen,dtsets(idtset)%ziontypat(1:dtsets(idtset)%ntypat)
   end do
  end if
 end if

 do idtset=0,ndtset_alloc
  dprarr(1:npsp,idtset)=dtsets(idtset)%znucl(1:npsp)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,4,marr,npsp,ndtset_alloc,'znucl','DPR')

 deallocate(dprarr,intarr)
 deallocate(jdtset_,response_,xangst_,xcart_)

 write(message,'(a,80a)')ch10,('=',mu=1,80)
 call wrtout(iout,message,'COLL')

!**************************************************************************

!DEBUG
!write(6,*)' outvars : end of subroutine '
!if(.true.)stop
!ENDDEBUG

end subroutine outvars
!!***
