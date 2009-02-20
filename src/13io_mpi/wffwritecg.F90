!{src2tex{texfont=tt}}
!!****f* ABINIT/wffwritecg
!! NAME
!! wffwritecg
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, MVer,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! FUNCTION
!! This procedure write cg in the file .o_WFK using MPI_IO 
!! in case of GS calculation
!!  for one kpt cg are dispatch amoung commcart communicator
!!
!! cg is written like that :
!!   BeginMarker cg (iproc =0, iband =1 ) cg(iproc=1, iband=1) ... cg(iproc=nproc-1, iband=1) EndMarker  
!!  BeginMarker cg (iproc =0, iband =2 ) cg(iproc=1, iband=2) ... cg(iproc=nproc-1, iband=2) EndMarker  
!! ....
!! ...
!!  BeginMarker cg(iproc=0, iband= nband_disk ) ../..  cg(iproc=nproc-1, iband=nband_disk) EndMarker
!! BeginMarker and EndMarker are given the value of the total length 
!! of cg for one band 
!! For writting cg, for improving performance we use a view of the file being written
!!  each for one proc. 
!!
!! INPUTS
!!  wff=struct info for wavefunction
!!  nband_disk =number of bands on disk files to be write
!!  cg(2,npw*nspinor*mband)=planewave coefficients of wavefunctions,
!!  icg=shift to be given to the location of the cg array
!!  mcg=dimention of cg
!! npwso =npw*nspinor 
!! with npw =number of plane waves
!! and nspinor =number of spinotial components of wavefunctions 
!! spaceComm = Communication space where are dispatched the cg ( commcart)
!!
!! OUTPUT
!! ierr1 = Error status
!!
!! PARENTS
!!     writewf
!!
!! CHILDREN
!!
!!SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wffwritecg(wff,cg,mcg, icg,nband_disk,npwso, spaceComm, ierr1)
 use defs_basis
 use defs_datatypes 

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: nband_disk,npwso
 integer,intent(in) :: mcg 
 integer,intent(in) :: icg 
 integer,intent(out) :: ierr1
 real(dp),intent(in):: cg(2,mcg)
 integer, intent(in) :: spaceComm 
!integer :: status1(MPI_STATUS_SIZE)
 

integer :: length(2), depl(2), typ(2)
integer :: length1(3), depl1(3), typ1(3)
integer :: myrank
integer :: nbOct_dp
integer :: ierr 
integer :: filetype, filetype1
integer :: delim_record
integer :: nbproc, iproc
 integer, allocatable:: local_offset(:)
integer :: total_taille, total_offset
integer :: i,j,j1,k
real(dp) ::  buff(nband_disk*2*npwso)
 integer, allocatable :: liste_taille(:)
integer , allocatable :: buffdelim(:)
 integer(abinit_offset) :: offset
 integer(abinit_offset) :: offset_zero = 0
integer :: fh,fh1
integer :: iband
integer :: n3 , np1,kk
integer :: tt, tu, ipw 

 ierr =0
#if defined MPI_IO
 call MPI_COMM_SIZE(spaceComm, nbproc, ierr)
 call MPI_COMM_RANK(spaceComm, myrank, ierr)

 allocate(liste_taille(nbproc))
 n3=2*npwso
call MPI_ALLGATHER(n3,1,MPI_INTEGER,liste_taille,1,MPI_INTEGER,spaceComm,ierr)
 
!First we put cg in a buffer
      kk= 1
         do iband=1,nband_disk
           ipw=(iband-1)*npwso+icg
              do tu=ipw+1, ipw+npwso
            do tt=1, 2
               buff(kk) = cg(tt,tu)
              kk=kk+1
!               PRINT *, 'iband ', iband, 'me ', myrank, 'tt ', tt , &
!      &         'tu ', tu, 'cg ', cg(tt,tu), 'kk ', kk
              enddo
          enddo
         enddo
 

 wff%off_recs = wff%offwff
 allocate(local_offset(nbproc))

local_offset(1)=wff%nbOct_int !begin marker
do iproc=2,nbproc
   local_offset(iproc)=local_offset(iproc-1)+liste_taille(iproc-1)*wff%nbOct_dp
enddo
total_taille=local_offset(nbproc) + liste_taille(nbproc)*wff%nbOct_dp
!end mark
total_taille=total_taille +wff%nbOct_int

!Writing cg in one collective order
! using a view ( better performance)

    depl(1)=0
    length(1) = liste_taille(myrank+1)*wff%nbOct_dp
    typ(1)=MPI_BYTE
    depl(2)=total_taille
    length(2)=0
    typ(2)=MPI_BYTE
  
    call MPI_TYPE_INDEXED(2,length, depl, typ, FILETYPE, ierr)
    call MPI_TYPE_COMMIT(FILETYPE, ierr)

    offset=wff%offwff+local_offset(myrank+1)

!Defining the view corresponding to cg
   call MPI_FILE_SET_VIEW(wff%fhwff, offset, MPI_BYTE, filetype, &
      & "native", MPI_INFO_NULL, ierr)
 
!writing cg
   call MPI_FILE_WRITE_ALL(wff%fhwff, buff, 2*npwso*nband_disk, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)

   call MPI_TYPE_FREE(FILETYPE, ierr)

!writing marks in only one order for performance reason
!first we put all marks in a buffer
    allocate(buffdelim(2*nband_disk))
    j=1
   do iband = 1, nband_disk
     delim_record = total_taille -2*wff%nbOct_int
     buffdelim(j) =  delim_record
     j=j+1
     buffdelim(j) =  delim_record
     j=j+1
   enddo

   depl1(1)=0
   length1(1)= wff%nbOct_int
   typ1(1)= MPI_BYTE
   depl1(2)= total_taille -wff%nbOct_int 
length1(2)= wff%nbOct_int
   typ1(2)= MPI_BYTE
   depl1(3)=total_taille
   length1(3)= 0
   typ1(3)= MPI_BYTE

    call MPI_TYPE_INDEXED(3,length1, depl1, typ1, FILETYPE1, ierr)

    call MPI_TYPE_COMMIT(FILETYPE1, ierr)

!Defining the view corresponding to markers 
   call MPI_FILE_SET_VIEW(wff%fhwff, wff%offwff, MPI_BYTE, filetype1, &
      & "native", MPI_INFO_NULL, ierr)
!Only one process write
   if (myrank == 0 ) then
   call MPI_FILE_WRITE(wff%fhwff, buffdelim, 2*nband_disk, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
   endif
   call MPI_TYPE_FREE(FILETYPE1, ierr)
    deallocate(buffdelim)
    deallocate(local_offset)
    deallocate(liste_taille)

    wff%offwff=wff%offwff+total_taille*nband_disk
 
!Reinit View by the default view
   call MPI_FILE_SET_VIEW(wff%fhwff, offset_zero, MPI_BYTE, MPI_BYTE, &
      & "native", MPI_INFO_NULL, ierr)


#endif
    ierr1=ierr 
   end subroutine wffwritecg
!!***
