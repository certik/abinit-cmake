!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_transitions
!! NAME
!! make_transitions
!!
!! FUNCTION
!!  Calculate transition energies entering the espression for the irreducible polarizability 
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nsspol=1 for spin unpolarized, 2 for spin polarized calculations
!!  nbnds=total number of bands
!!  kmesh<bz_mesh_type>=datatype gathering info on the k-mesh:
!!   | %nbz=number of k-points in the full BZ
!!   | %nibz=number of k-points in the IBZ
!!   | $tab(nkbz)=table giving for each k-point in the BZ, the corresponding irreducible point in the IBZ array
!!   | %bz(3,nkbz)=reduced coordinated of k-points
!!  mG0(3)=integer defining the number of shells in which the umklapp G0 vector has to be found
!!  TOL_DELTA_OCC=tolerance on the difference of the occupation numbers
!!  gwenergy(kmesh%nkibz,nbnds,nsppol)=quasi-particle energies energies 
!!  occ(kmesh%nkibz,nbnds,nsppol)=occupation numbers
!!  chi0alg=integer defining the method used to calculate chi0
!!   0 ==> calculate chi0 using the Adler-Wiser expression
!!   1 ==> use spectral method 
!!  timrev=if 2, time-reversal symmetry is considered; 1 otherwise
!!
!! OUTPUT
!!  
!!
!! SIDE EFFECTS
!!
!! NOTES
!! 
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

subroutine make_transitions(chi0alg,nbnds,nbvw,nsppol,symchi,timrev,TOL_DELTA_OCC,zcut,&
& max_rest,min_rest,my_max_rest,my_min_rest,kmesh,ltg_q,mpi_enreg,mG0,gwenergy,occ,qpoint)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => make_transitions
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chi0alg,nbnds,nbvw,nsppol,symchi,timrev
 real(dp),intent(in) :: TOL_DELTA_OCC,zcut
 real(dp),intent(out) :: max_rest,min_rest
 real(dp),intent(out) :: my_max_rest,my_min_rest
 type(bz_mesh_type),intent(in) :: kmesh
 type(little_group),intent(in) :: ltg_q
 type(MPI_type    ),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: mG0(3)
 real(dp),intent(in) :: gwenergy(kmesh%nibz,nbnds,nsppol)
 real(dp),intent(in) :: occ(kmesh%nibz,nbnds,nsppol),qpoint(3)

!Local variables-------------------------------
!scalars
 integer :: ib1,ib2,ii,ikbz,ik_ibz,ikmq_bz,ikmq_ibz,is,nt,io,ntrans,my_ntrans
 integer :: ier,rank,spaceComm,nprocs,master,iloop
 real(dp) :: delta_ene,delta_occ,spin_fact
 character(len=500) :: msg
!arrays
 integer :: G0(3)
 real(dp) :: kmq(3)

!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' make_transitions : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 if (chi0alg<0 .or. chi0alg>=2) then 
  write(msg,'(4a,i3,a)')ch10,&
&  ' make_transitions : BUG ',ch10,&
&  ' chi0alg = ',chi0alg,' not allowed '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if 
 if (timrev/=1 .and. timrev/=2) then 
  write(msg,'(4a,i3,a)')ch10,&
&  ' make_transitions : BUG ',ch10,&
&  ' timrev = ',timrev,' not allowed'
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if 
 !
 ! Initialize MPI quantities.
 call xcomm_init  (mpi_enreg,spaceComm) 
 call xmaster_init(mpi_enreg,master   )  
 call xme_init    (mpi_enreg,rank     )          
 call xproc_max(nprocs,ier)

 if (nprocs==1 .and. mpi_enreg%gwpara/=0) then 
  write(msg,'(4a,i4)')ch10,&
&  ' make_transitions : BUG : ',ch10,&
&  ' in a sequential run gwpara should be 0 while it is ',mpi_enreg%gwpara 
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if
 !
 ! In the first loop calculate total number of transitions for this q-point
 ! as well min and max transition without taking into account distribution of bands. 
 ! In the second iteration calculate min and Max transition for this processor.
 !
 spin_fact=half ; if (nsppol==2) spin_fact=one
 my_max_rest=smallest_real ; my_min_rest=greatest_real
    max_rest=smallest_real ;    min_rest=greatest_real
 print*,"done"

 do iloop=1,2
  nt=0
 print*,"in loop",iloop
  do ikbz=1,kmesh%nbz

   ik_ibz=kmesh%tab(ikbz)
   kmq(:)=kmesh%bz(:,ikbz)-qpoint(:)

   if (symchi==1) then  
    if (ltg_q%ibzq(ikbz)/=1) cycle 
    ! This point does not belong to the IBZ defined by the little group
   end if 
   !
   ! Find kp=k-q-G0 and also G0 where kp is in the first BZ
   if (.not.has_BZ_item(Kmesh,kmq,ikmq_bz,g0)) then
    ! Stop as the weight 1.0/nkbz is wrong and should be changed in cchi0/cchi0q0
    write(msg,'(4a,2(2a,3f12.6),2a)')ch10,&
&    ' make_transitions : ERROR - ',ch10,&
&    ' kp  = k-q-G0 not found in the BZ mesh',ch10,&
&    ' k   = ',(kmesh%bz(ii,ikbz),ii=1,3),ch10,&
&    ' k-q = ',(kmq(ii),ii=1,3),ch10,&
&    ' weight in cchi0/cchi0q is wrong ' 
    call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
   end if 

   ikmq_ibz=kmesh%tab(ikmq_bz)
   do is=1,nsppol

    do ib1=1,nbnds           
     if (iloop==2 .and. mpi_enreg%gwpara==2) then
      if (mpi_enreg%proc_distrb(ik_ibz,ib1,is)/=rank) cycle
     end if

     do ib2=1,nbnds
      if (timrev==2 .and. ib1<ib2) cycle ! Thanks to time-reversal we gain a factor ~2.
      if (iloop==2 .and. mpi_enreg%gwpara==2) then
       if (mpi_enreg%proc_distrb(ik_ibz,ib2,is)/=rank) cycle
      end if

      ! Take care of "valence-valence" transitions in case of gwpara==2
      if (iloop==2 .and. mpi_enreg%gwpara==2) then
       if (ib1<=nbvw .and. ib2<=nbvw .and. rank/=master) cycle 
      end if

      delta_occ=spin_fact*(occ(ikmq_ibz,ib1,is)-occ(ik_ibz,ib2,is))
      delta_ene=gwenergy(ikmq_ibz,ib1,is)-gwenergy(ik_ibz,ib2,is)

      select CASE (chi0alg)
       CASE (0)  
        ! Adler-Wiser expression.
        ! Skip only if factor due to occupation number is smaller than TOL_DELTA_OCC
        if (abs(delta_occ) < abs(TOL_DELTA_OCC)) cycle
       CASE (1)
        ! Spectral method with time-reversal, only resonant transitions 
        ! This has to changed to include spectral method without time-reversal
        if (delta_ene < -abs(TOL_DELTA_OCC) .or. abs(delta_occ) < abs(TOL_DELTA_OCC)) cycle
      end select 
      !
      ! We have a new transition
      nt=nt+1

      if (iloop==1) then 
       max_rest=MAX(max_rest,zero,delta_ene)
       if (delta_ene>=-tol6) min_rest=MIN(min_rest,delta_ene)
      end if
      if (iloop==2) then 
       my_max_rest=MAX(my_max_rest,zero,delta_ene)
       if (delta_ene>=-tol6) my_min_rest=MIN(my_min_rest,delta_ene)
      end if

     end do 
    end do 
   end do
  end do
  if (iloop==1) ntrans=nt
  if (iloop==2) my_ntrans=nt
 end do !iloop

 !if (mpi_enreg%gwpara==2) then 
 ! call xmin_mpi(my_min_rest,min_rest,spaceComm,ier)
 ! call xmax_mpi(my_max_rest,max_rest,spaceComm,ier)
 !end if

 write(msg,'(2a,i9,2a,f8.3,3a,f8.3,a)')ch10,&
& ' Total number of transitions = ',ntrans,ch10,&
& ' min resonant     = ',min_rest*Ha_eV,' [eV] ',ch10,&
& ' Max resonant     = ',max_rest*Ha_eV,' [eV] '
 call wrtout(std_out,msg,'COLL')
 if (nprocs/=1) then 
  write(msg,'(2a,i9,2a,f8.3,3a,f8.3,a)')ch10,&
&  ' Total number of transitions for this processor= ',my_ntrans,ch10,&
&  ' min resonant     = ',my_min_rest*Ha_eV,' [eV] ',ch10,&
&  ' Max resonant     = ',my_max_rest*Ha_eV,' [eV] '
  call wrtout(std_out,msg,'PERS')
 end if

#if defined DEBUG_MODE
 write(msg,'(a)')' make_transitions : exit'
 call wrtout(std_out,msg,'COLL')
 call flush_unit(std_out)
#endif

end subroutine make_transitions 
!!***

subroutine nullify_transitions(self)

  use defs_basis
  use defs_datatypes

  implicit none

  type(transitions_type),intent(out) :: self

  nullify(self%bands)
  nullify(self%distrb)
  nullify(self%ik_ibz)
  nullify(self%ikmq_ibz)
  nullify(self%ik_bz)
  nullify(self%ikmq_bz)    
  nullify(self%G0) 
  nullify(self%spin)      
  nullify(self%delta_occ)
  nullify(self%qpoint)
  nullify(self%delta_ene) 
  nullify(self%num_w)

end subroutine nullify_transitions 


subroutine init_transitions(self,nkbz,nbnds,nbvw,nsppol,nomega,qpoint,ntrans)

   use defs_basis
   use defs_datatypes

   implicit none

   type(transitions_type),intent(out) :: self
   integer,intent(in) :: nkbz,nbnds,nbvw,nsppol,nomega
   integer,optional,intent(in) :: ntrans
   real(dp),intent(in) :: qpoint(3)
   integer :: nt

#if defined DEBUG_MODE
write(*,*)" init_transitions : enter"
#endif

   self%nbnds=nbnds
   self%nbvw=nbvw
   self%nomega=nomega
   self%nkbz=nkbz
   self%nsppol=nsppol

   self%ntrans=nkbz*nbnds**2*nsppol
   if (present(ntrans)) self%ntrans=ntrans


   nt=self%ntrans
   allocate(self%bands(2,nt)     )  ;  self%bands(:,:)    =0
   allocate(self%distrb(nt)      )  ;  self%bands(:,:)    =-999
   allocate(self%ik_ibz(nt)      )  ;  self%ik_ibz(:)     =0
   allocate(self%ikmq_ibz(nt)    )  ;  self%ikmq_ibz(:)   =0
   allocate(self%ik_bz(nt)       )  ;  self%ik_bz(:)      =0
   allocate(self%ikmq_bz(nt)     )  ;  self%ik_bz(:)      =0
   allocate(self%G0(3,nt)        )  ;  self%G0(:,:)       =0
   allocate(self%spin(2,nt)      )  ;  self%spin(:,:)     =0
   allocate(self%delta_occ(nt)   )  ;  self%delta_occ(:)  =zero
   allocate(self%qpoint(3)       )  ;  self%qpoint(:)     =qpoint(:)
   allocate(self%delta_ene(nt)   )  ;  self%delta_ene(:)  =czero
   allocate(self%num_w(nomega,nt))  ;  self%num_w(:,:)    =czero

#if defined DEBUG_MODE
write(*,*)" init_transitions : exit"
#endif

end subroutine init_transitions 


subroutine destroy_transitions(self)

  use defs_basis
  use defs_datatypes

  implicit none

  type(transitions_type),intent(inout) :: self

   if (associated(self%bands    ))  deallocate(self%bands    )    
   if (associated(self%distrb   ))  deallocate(self%distrb   )    
   if (associated(self%ik_ibz   ))  deallocate(self%ik_ibz   )
   if (associated(self%ikmq_ibz ))  deallocate(self%ikmq_ibz )
   if (associated(self%ik_bz    ))  deallocate(self%ik_bz    )
   if (associated(self%ikmq_bz  ))  deallocate(self%ikmq_bz  ) 
   if (associated(self%G0       ))  deallocate(self%G0       )
   if (associated(self%spin     ))  deallocate(self%spin     ) 
   if (associated(self%delta_occ))  deallocate(self%delta_occ)  
   if (associated(self%qpoint   ))  deallocate(self%qpoint   )  
   if (associated(self%delta_ene))  deallocate(self%delta_ene)   
   if (associated(self%num_w    ))  deallocate(self%num_w    )

end subroutine destroy_transitions 

subroutine copy_transitions(t_in,t_out)

  use defs_basis 
  use defs_datatypes

  implicit none

  type(transitions_type),intent(in) :: t_in
  type(transitions_type),intent(out) :: t_out

  !local
  integer :: nt,nw

  nt=min(t_in%ntrans,t_out%ntrans) 
  nw=min(t_in%nomega,t_out%nomega)

  if (t_in%ntrans/=t_out%ntrans) then 
   write(*,*)" Copying ",nt,"/",t_in%ntrans," transitions from input datatype"
  end if
  if (t_in%ntrans/=t_out%ntrans) then 
   write(*,*)" Copying ",nw,"/",t_in%nomega," frequencies from input datatype"
  end if 

  t_out%nbnds =t_in%nbnds
  t_out%nbvw  =t_in%nbvw
  t_out%nkbz  =t_in%nkbz
  t_out%nsppol=t_in%nsppol

  ! Use whit care!
  t_out%my_min_res=t_in%my_min_res
  t_out%my_max_res=t_in%my_max_res

  ! Copy a smaller number of transitions
  t_out%ntrans=nt
  t_out%nomega=nw

  t_out%bands(:,1:nt)      = t_in%bands(:,1:nt)   
  t_out%distrb(1:nt)       = t_in%distrb(1:nt)
  t_out%ik_ibz(1:nt)       = t_in%ik_ibz(1:nt)
  t_out%ikmq_ibz(1:nt)     = t_in%ikmq_ibz(1:nt)
  t_out%ik_bz(1:nt)        = t_in%ik_bz(1:nt)
  t_out%ikmq_bz(1:nt)      = t_in%ikmq_bz(1:nt) 
  t_out%G0(:,1:nt)         = t_in%G0(:,1:nt)       
  t_out%spin(:,1:nt)       = t_in%spin(:,1:nt)
  t_out%delta_occ(1:nt)    = t_in%delta_occ(1:nt)
  t_out%qpoint(:)          = t_in%qpoint(:)
  t_out%delta_ene(1:nt)    = t_in%delta_ene(1:nt)   
  t_out%num_w(1:nw,1:nt)   = t_in%num_w(1:nw,1:nt)

end subroutine copy_transitions

subroutine print_transitions(self,unit,prtvol)

 use defs_basis
 use defs_datatypes

 implicit none
 
 integer,optional,intent(in) :: unit,prtvol
 type(transitions_type),intent(in) :: self

 integer :: unt,verbose

 unt=std_out
 if (present(unit)) unt=unit
 verbose=0
 if (present(prtvol)) verbose=prtvol

 write(unt,*)" q-point : ",self%qpoint(:)
 write(unt,*)" Total number of transitions : ",self%ntrans 

end subroutine print_transitions


subroutine get_my_extrema(self,my_min_res,my_max_res)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

 type(transitions_type),intent(in) :: self
 real(dp),intent(out) :: my_min_res,my_max_res
 
 integer :: it,my_rank,my_ntrans
 real(dp),allocatable :: ene(:)

 if (all(self%distrb(:)==-999)) then 
  write(*,*)" get_my_extrema : ERROR -",ch10,&
&  " Either BUG while distributing or ",ch10,&
&  " number of processors exceeds number of transitions"
  call leave_new('COLL')
 end if

 call xme_whoiam(my_rank)
 my_ntrans=COUNT(self%distrb==my_rank)
 if (self%my_ntrans/=my_ntrans) stop "BUG in distribution"
 
 allocate(ene(my_ntrans)) ; ene=zero
 !do it=1,self%ntrans
 ! if (self%distrb(it)==my_rank) 
 ! ene(it)=self%delta_ene(it)
 !end do
 ene(:)= PACK(real(self%delta_ene), MASK=self%distrb==my_rank)

 my_min_res = MINVAL(ene, MASK=ene>=zero)
 my_max_res = MAXVAL(ene, MASK=ene>=zero)

end subroutine get_my_extrema 


subroutine split_transitions(self,mpi_enreg)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => split_transitions
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(transitions_type),intent(inout) :: self
 type(MPI_type),intent(in) :: mpi_enreg

!Local variables ------------------------------
!scalars
 integer :: ikibz,ib1,ib2,is,it
 integer :: nproc,ntrans,my_rank
 integer :: my_start,my_stop,my_ntrans,master
 real(dp) :: my_min_res,my_max_res
 character(len=500) :: msg
!arrays

 call xme_whoiam(my_rank)
 call xmaster_init(mpi_enreg,master) 
 nproc =mpi_enreg%nproc
 ntrans=self%ntrans

 self%distrb(:)=-999
 my_ntrans=0

 select case (mpi_enreg%gwpara) 

  case (0)
   !
   ! Sequential run
   !
   write(*,*)" split_transitions : Sequential run, nothing to do"
   self%distrb(:)=my_rank
   my_ntrans=self%ntrans

  case (1)
   !
   ! Parallelization over k-points, each node has the entire set of wavefunctions.
   ! Split the total number of transitions, the data access is optimized in make_transitions.
   !
   call split_work(ntrans,my_start,my_stop)

   self%distrb(my_start:my_stop)=my_rank
   my_ntrans=my_stop-my_start+1

  case (2)
   !
   ! Parallelism over band. Be careful because wavefunctions are separated into 
   ! "valence" and "conduction". "valence" means all the states with occupation smaller
   ! than a TOL_DELTA_OCC value. Remember that "conduction" states are distributed but "valence"
   ! states are on each node. Avoid multiple distribution in case of metallic systems.
   ! Moreover Im assumig that valence index is always on the left, this distribution 
   ! does not work if we dont use time-reversal.
   ! PResently only master takes care of the valence-valence trasitions in metals 

   if (self%nbvw<=0) stop "BUG while distributing"
   !
   do it=1,self%ntrans
    ib1 = self%bands(1,it)
    ib2 = self%bands(2,it)

    if (ib1<=self%nbvw .and. ib2<=self%nbvw) then 
     ! Each processor has these states, for the moment only master than we can optimize
     if (my_rank/=master) cycle
    end if

    if ( any(mpi_enreg%proc_distrb(:,ib1,:)==my_rank) .and. &
&        any(mpi_enreg%proc_distrb(:,ib2,:)==my_rank)       & 
&      ) then 
     self%distrb(it)=my_rank
     my_ntrans=my_ntrans+1
    end if 
   end do

!DEBUG this is just to keep in mind that I cannot parallelize memory over spin
   do ib1=1,size(mpi_enreg%proc_distrb, DIM=2)
    do is=1,size(mpi_enreg%proc_distrb, DIM=3)
     if (any(mpi_enreg%proc_distrb(:,ib1,is)/=mpi_enreg%proc_distrb(1,1,is))) then
      write(*,*)"BUG while distributing wavefunctions"
      call leave_new('COLL')
     end if 
    end do
   end do
!ENDDEBUG

  case DEFAULT
   write(msg,'(4a)')ch10,&
&   " split_transitions : BUG-",ch10,&
&   " called with wrong value of gwpara "
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end select

 self%my_ntrans=my_ntrans

 call get_my_extrema(self,my_min_res,my_max_res)
 self%my_min_res=my_min_res
 self%my_max_res=my_max_res

 write(msg,'(2a,i4,a,i5,a)')ch10,&
& " Processor ",my_rank," will treat ",my_ntrans," transitions "
 call wrtout(std_out,msg,'PERS')
 write(*,*)" min and Max resonant transitions = ",my_min_res,my_max_res

end subroutine split_transitions


subroutine find_my_indeces(self,nomega,omega_mesh,my_w1,my_w2)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(transitions_type),intent(in) :: self
 integer,intent(in) :: nomega
 integer,intent(out) :: my_w1,my_w2
!arrays 
 real(dp) :: omega_mesh(nomega)

!Local variables-------------------------------
!scalars
 integer :: io
 character(len=500) :: msg
!************************************************************************

 ! Note nomega-1, that is because omega_mesh encloses the 
 ! interval made of the possible resonant transitions
 ! This part is very sensitive to changes in setup_mesh
 my_w2=-999
 do io=1,nomega-1
  if (omega_mesh(io) > self%my_max_res) then 
   my_w2=io+1
   exit
  end if 
 end do

 my_w1=-999
 do io=nomega,1,-1
  if (omega_mesh(io)<= self%my_min_res) then ! Check metals
   my_w1=io
   exit
  end if 
 end do 

 if (my_w1==-999 .or. my_w2==-999) then 
  write(msg,"(a,2i4)")" ERROR in find_my_indeces ",my_w1,my_w2
  call wrtout(std_out,msg,'PERS') ; call leave_new('COLL')
 end if 

end subroutine find_my_indeces


subroutine gw2abi(nkibz,nbnds,nsppol,gw_array,abi_array)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkibz,nbnds,nsppol 
!arrays
 real(dp),intent(in) :: gw_array(nkibz,nbnds,nsppol)
 real(dp),intent(out) :: abi_array(nbnds*nkibz*nsppol)

!Local variables-------------------------------
 integer :: is,ik,ib,idx
 character(len=500) :: msg

! *************************************************************************

 if (SIZE(abi_array)/=SIZE(gw_array)) then 
  write(msg,'(4a,i8)')ch10,&
&  ' gw2abi BUG - ',ch10,&
&  ' input and output arrays have different dimension ',SIZE(abi_array),SIZE(gw_array)
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 idx=0
 do is=1,nsppol 
  do ik=1,nkibz
   do ib=1,nbnds
    idx=idx+1
    abi_array(idx)=gw_array(ik,ib,is)
   end do
  end do
 end do

end subroutine gw2abi


subroutine abi2gw(nkibz,nbnds,nsppol,abi_array,gw_array)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkibz,nbnds,nsppol 
!arrays
 real(dp),intent(in) :: abi_array(nbnds*nkibz*nsppol)
 real(dp),intent(out) :: gw_array(nkibz,nbnds,nsppol)

!Local variables-------------------------------
 integer :: is,ik,ib,idx
 character(len=500) :: msg

! *************************************************************************

 if (SIZE(abi_array)/=SIZE(gw_array)) then 
  write(msg,'(4a,i8)')ch10,&
&  ' abi2gw BUG - ',ch10,&
&  ' input and output arrays have different dimension ',SIZE(abi_array),SIZE(gw_array)
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 idx=0
 do is=1,nsppol 
  do ik=1,nkibz
   do ib=1,nbnds
    idx=idx+1
    gw_array(ik,ib,is)=abi_array(idx)
   end do
  end do
 end do

end subroutine abi2gw
