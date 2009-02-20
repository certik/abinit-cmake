!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_gkk_qpt_tr
!!
!! NAME
!! get_gkk_qpt_tr
!!
!! FUNCTION
!!  calculate the product of gkk_qpt with the (in) and (out) velocity factors for transport
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (JPC)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  elph_ds
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%nFSkpt = number of kpts included in the FS integration
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%minFSband = index of the lowest FS band
!!    elph_ds%nqpt  = number of Q pts 
!!    elph_ds%gkk_qpt = standrad gkk matrix
!!  to index the GS electronic states :
!!  hdr%nkpt = full number of k points 
!!  nband =full number of bands 
!!  FSfulltoirred = mapping of full FS kpts to irreducible ones
!!   FSfullpqtofull = mapping of k + q to k
!!   FSirredtoGS = mapping of irreducible kpoints to GS set
!! OUTPUT
!! elph_tr_ds%gkk_qpt_trout = out gkk_qpt
!! elph_tr_ds%gkk_qpt_trin = in_gkk_qpt
!! elph_tr_ds%FSelecveloc_sq = avergae FS electronic velocity
!! SIDE EFFECTS
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      read_el_veloc
!!
!! NOTES
!!   
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine get_gkk_qpt_tr(elph_ds,mpi_enreg,nband,hdr,FSfulltoirred, FSirredtoGS,FSfullpqtofull,elph_tr_ds)

  use defs_datatypes
  use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_17ddb, except_this_one => get_gkk_qpt_tr
!End of the abilint section

  implicit none


!Arguments ------------------------------------
  !scalars
  integer,intent(in) :: nband
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  type(hdr_type),intent(in) :: hdr
  type(elph_type),intent(in) :: elph_ds
  type(elph_tr_type):: elph_tr_ds


  integer,intent(in) :: FSfulltoirred(3,elph_ds%nFSkpt), FSirredtoGS(elph_ds%nFSkptirred)
  integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)

!Local variables-------------------------------
  !scalars
  integer:: iqpt,iFSkpt,ikpt,iFSkptq,ikptpq
  integer:: ib1,ib2,fib1,fib2,isppol,ibeff
  real(dp)::eta1,eta2,eta,etain,etaout
  real(dp)::n0(elph_ds%nsppol)
  !arrays
  real(dp):: elvelock(3),elvelockpq(3)
  character(len=fnlen) :: fname
  character(len=500) :: message
  integer::iost,ierr

  real(dp),allocatable :: tmp_gkkin(:,:,:,:,:)
  real(dp),allocatable :: tmp_gkkout(:,:,:,:,:)
  real(dp),allocatable :: tmp_gkk(:,:,:,:,:)

! *********************************************************************

!check inputs
!TODO: should be done at earlier stage of initialization and checking
 write (*,*) ' ngkkband,nFSband,minFSband = ', elph_ds%ngkkband,elph_ds%nFSband,elph_ds%minFSband
 if (elph_ds%ngkkband /= elph_ds%nFSband) then
  write (message,'(a)') 'Error: need to keep electron band dependency in memory for transport calculations'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 elph_tr_ds%unitgkq_trin=40 
 elph_tr_ds%unitgkq_trout=41

!if the gkk_tr already are on disk
 if (elph_ds%gkqexist==1) then

  fname=trim(elph_ds%elph_base_name) // '_GKKQ_trin'
  open (unit=elph_tr_ds%unitgkq_trin,file=fname,access='direct',recl=elph_tr_ds%onegkksize,&
&  form='unformatted',status='old',iostat=iost)
  if (iost /= 0) then
   write (message,'(3a)')' get_gkk_qpt_tr : ERROR- opening file ',trim(fname),' as old'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  fname=trim(elph_ds%elph_base_name) // '_GKKQ_trout'
  open (unit=elph_tr_ds%unitgkq_trout,file=fname,access='direct',recl=elph_tr_ds%onegkksize,&
&  form='unformatted',status='old',iostat=iost)

  if (iost /= 0) then
   write (message,'(3a)')' get_gkk_qpt_tr : ERROR- opening file ',trim(fname),' as old'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  return
 end if

 if (elph_ds%gkqwrite==1) then

  fname=trim(elph_ds%elph_base_name) // '_GKKQ_trin'
  write(6,*)'SIZES ', elph_tr_ds%onegkksize
  open (unit=elph_tr_ds%unitgkq_trin,file=fname,access='direct',recl=elph_tr_ds%onegkksize,&
&  form='unformatted',status='new',iostat=iost)

  if (iost /= 0) then
   write (message,'(3a)')' get_gkk_qpt_tr : ERROR- opening file ',trim(fname),' as new'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if


  fname=trim(elph_ds%elph_base_name) // '_GKKQ_trout'
  open (unit=elph_tr_ds%unitgkq_trout,file=fname,access='direct',recl=elph_tr_ds%onegkksize,&
&  form='unformatted',status='new',iostat=iost)

  if (iost /= 0) then
   write (message,'(3a)')' get_gkk_qpt_tr : ERROR- opening file ',trim(fname),' as new'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if


 end if


 write(6,*)'reading of electronic velocities'
 allocate(elph_tr_ds%el_veloc(hdr%nkpt,nband,3,elph_ds%nsppol))
 call read_el_veloc(mpi_enreg,nband,hdr%nkpt,elph_ds%nsppol,elph_tr_ds)

!write(6,*)elph_tr_ds%el_veloc
 allocate (elph_tr_ds%FSelecveloc_sq(elph_ds%nsppol))
 elph_tr_ds%FSelecveloc_sq=0.
 do isppol=1,elph_ds%nsppol
  do iFSkpt=1,elph_ds%nFSkpt
   ikpt=FSirredtoGS(FSfulltoirred(1,iFSkpt))
   do ib1=1,elph_ds%nFSband
    fib1=ib1+elph_ds%minFSband-1
    elvelock(:)=elph_tr_ds%el_veloc(ikpt,fib1,:,isppol)
    eta2=elvelock(1)*elvelock(1)+elvelock(2)*elvelock(2)+elvelock(3)*elvelock(3)  
!   MJV: removed 2/7/2007 because never used below
!   eta1=eta2*elph_ds%FSintweight(ib1,iFSkpt,isppol)/elph_ds%nFSkpt
    elph_tr_ds%FSelecveloc_sq(isppol)=elph_tr_ds%FSelecveloc_sq(isppol)&
&    +eta2*elph_ds%FSintweight(ib1,iFSkpt,isppol)/elph_ds%nFSkpt
   end do
  end do
  elph_tr_ds%FSelecveloc_sq(isppol) =elph_tr_ds%FSelecveloc_sq(isppol)/elph_ds%n0(isppol)
 end do ! end isppol


!calculates the in and out factor and multiplies

 if (elph_ds%gkqwrite == 0) then
  write (message,'(a)')' get_gkk_qpt_tr : keeping matrices in memory'
  call wrtout(06,message,'COLL')



  do iqpt=1,elph_ds%nqpt

   do isppol=1,elph_ds%nsppol

    do iFSkpt=1,elph_ds%nFSkpt
     ikpt=FSirredtoGS(FSfulltoirred(1,iFSkpt))        !k     
     iFSkptq = FSfullpqtofull(iFSkpt,iqpt)
     ikptpq=FSirredtoGS(FSfulltoirred(1,iFSkptq))     !k'=k+q

     do ib1=1,elph_ds%nFSband
      fib1=ib1+elph_ds%minFSband-1
      elvelock(:)=elph_tr_ds%el_veloc(ikpt,fib1,:,isppol)
!     write(6,*)'elvelock',elvelock
      do ib2=1,elph_ds%nFSband
       ibeff=ib2+(ib1-1)*elph_ds%nFSband

       fib2=ib2+elph_ds%minFSband-1
       elvelockpq(:)= elph_tr_ds%el_veloc(ikptpq,fib2,:,isppol)
!      write(6,*)'elvelockPQ',elvelockpq
       etain=elvelock(1)*elvelockpq(1)+elvelock(2)*elvelockpq(2)+elvelock(3)*elvelockpq(3)
       etaout=elvelock(1)*elvelock(1)+elvelock(2)*elvelock(2)+elvelock(3)*elvelock(3)
       etain=etain/elph_tr_ds%FSelecveloc_sq(isppol)
       etaout=etaout/elph_tr_ds%FSelecveloc_sq(isppol)
       eta=(etaout-etain)
!      write(6,*)'eta',etain,etaout,elph_ds%gkk_qpt(:,ibeff,:,iFSkpt,isppol,iqpt)
       elph_tr_ds%gkk_qpt_trout(:,ibeff,:,iFSkpt,isppol,iqpt)= elph_ds%gkk_qpt(:,ibeff,:,iFSkpt,isppol,iqpt)*etaout
       elph_tr_ds%gkk_qpt_trin(:,ibeff,:,iFSkpt,isppol,iqpt)= elph_ds%gkk_qpt(:,ibeff,:,iFSkpt,isppol,iqpt)*etain
      end do
     end do
    end do ! ik

   end do ! isppol


  end do ! iq
! write(987,*)elph_tr_ds%gkk_qpt_trout

 else if (elph_ds%gkqwrite == 1) then
  allocate (tmp_gkkin (2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,&
&  elph_ds%nFSkpt,elph_ds%nsppol),stat=ierr)
  if (ierr /= 0 ) then
   write (message,'(3a)')' get_gkk_qpt_tr : ERROR- ',ch10,&
&   ' trying to allocate array tmp_gkkin '
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  allocate (tmp_gkkout (2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,&
&  elph_ds%nFSkpt,elph_ds%nsppol),stat=ierr)
  if (ierr /= 0 ) then
   write (message,'(3a)')' get_gkk_qpt_tr : ERROR- ',ch10,&
&   ' trying to allocate array tmp_gkkout '
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  allocate (tmp_gkk (2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,&
&  elph_ds%nFSkpt,elph_ds%nsppol),stat=ierr)
  if (ierr /= 0 ) then
   write (message,'(3a)')' get_gkk_qpt_tr : ERROR- ',ch10,&
&   ' trying to allocate array tmp_gkkout '
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  
  do iqpt=1,elph_ds%nqpt
   write(6,*)'iqpt',iqpt
   read (elph_ds%unitgkq,REC=iqpt) tmp_gkk
   write(6,*)'read iqpt OK'
   do isppol=1,elph_ds%nsppol
    do iFSkpt=1,elph_ds%nFSkpt
     ikpt=FSirredtoGS(FSfulltoirred(1,iFSkpt))        !k     
     iFSkptq = FSfullpqtofull(iFSkpt,iqpt)
     ikptpq=FSirredtoGS(FSfulltoirred(1,iFSkptq))     !k'=k+q

     do ib1=1,elph_ds%nFSband
      fib1=ib1+elph_ds%minFSband-1
      elvelock(:)=elph_tr_ds%el_veloc(ikpt,fib1,:,isppol)

      do ib2=1,elph_ds%nFSband
       ibeff=ib2+(ib1-1)*elph_ds%nFSband

       fib2=ib2+elph_ds%minFSband-1
       elvelockpq(:)= elph_tr_ds%el_veloc(ikptpq,fib2,:,isppol)

       etain=elvelock(1)*elvelockpq(1)+elvelock(2)*elvelockpq(2)+elvelock(3)*elvelockpq(3)
       etaout=elvelock(1)*elvelock(1)+elvelock(2)*elvelock(2)+elvelock(3)*elvelock(3)
       etain=etain/elph_tr_ds%FSelecveloc_sq(isppol)
       etaout=etaout/elph_tr_ds%FSelecveloc_sq(isppol)
       eta=(etaout-etain)

       tmp_gkkout(:,ibeff,:,iFSkpt,isppol)= tmp_gkk(:,ibeff,:,iFSkpt,isppol)*etaout
       tmp_gkkin(:,ibeff,:,iFSkpt,isppol)= tmp_gkk(:,ibeff,:,iFSkpt,isppol)*etain

      end do
     end do
    end do
   end do
   write(6,*)'write 0'
   write (elph_tr_ds%unitgkq_trin,REC=iqpt) tmp_gkkin
   write(6,*)'write in'
   write (elph_tr_ds%unitgkq_trout,REC=iqpt) tmp_gkkout
   write(6,*)'write out'

  end do

  deallocate (tmp_gkk)
  deallocate (tmp_gkkout)
  deallocate (tmp_gkkin)


 end if
 deallocate(elph_tr_ds%el_veloc)
 write(6,*)'out of get_gkk_qpt_tr'
 return
end subroutine get_gkk_qpt_tr
!!***
