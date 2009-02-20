!{\src2tex{textfont=tt}}
!!****f* ABINIT/init_wf_info_1
!! NAME
!! init_wf_info_1
!!
!! FUNCTION
!!  This module (?) contains methods acting on the wafefunction object
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (FB,MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! igfft0(npwwfn)=index of each plane wave in FFT grid
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nk=number of k points
!! npwwfn=number of plane waves
!! nsppol=number of independent spin polarizations
!! tim_fourdp=4 if called from within screening ; =5 if called from within sigma
!! wfg(npwwfn,my_minb:my_maxb,nk,nsppol)=wavefunctions in reciprocal space treated by this processor.
!! my_minb,my_maxb = min and max band treated by this processor
!! MPI_enreg= datatype containing information on parallelism to be passed to fourdp
!!
!! OUTPUT
!!  wfr(ngfft(1)*ngfft(2)*ngfft(3),my_minb:my_maxb,nk,nsppol)
!!   wavefunctions in real space, for each band, k point and spin
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      fourdp,wrtout,xcomm_init,xmaster_init,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine init_wf_info_1(Wf_info,gwmem,paral_kgb,npwwfn,my_minb,my_maxb,nk,nsppol,nspden,nspinor)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : czero_gw


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: gwmem,paral_kgb,my_maxb,my_minb
 integer,intent(in) :: nk,npwwfn,nsppol,nspden,nspinor
 type(wavefunctions_information),intent(inout) :: Wf_info
!arrays

!Local variables ------------------------------
!scalars
 integer :: istat
 character(len=500) :: msg
 character(len=50),parameter :: sub_name='init_wf_info_1.F90'
!************************************************************************

 Wf_info%gwmem     = gwmem
 Wf_info%my_minb   = my_minb
 Wf_info%my_maxb   = my_maxb
 !TODO FFT should be done here, waiting for better restructuring of sigma and screening
 Wf_info%nk        = nk
 Wf_info%npwwfn    = npwwfn
 Wf_info%nspden    = nspden
 Wf_info%nsppol    = nsppol
 Wf_info%nspinor   = nspinor
 Wf_info%paral_kgb = paral_kgb

 allocate(Wf_info%is_already_stored(my_minb:my_maxb,nk,nsppol))
 Wf_info%is_already_stored(:,:,:)=.FALSE.

 allocate(Wf_info%wfg(npwwfn*nspinor,my_minb:my_maxb,nk,nsppol),stat=istat)
 if (istat/=0) then
  call memerr(sub_name,'wfg',nspinor*npwwfn*(my_maxb-my_minb+1)*nk*nsppol,'gwpc')
 end if
 Wf_info%wfg(:,:,:,:)=czero_gw

 write(msg,'(2a)')ch10,' initialization 1 of the Wf_info object done'
 call wrtout(std_out,msg,'COLL')

end subroutine init_wf_info_1
!!***


!!****f* ABINIT/init_wf_info_2
!! NAME
!!  init_wf_info_2 
!!
!! FUNCTION
!!  Finalize the initialization of the data structure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine init_wf_info_2(Wf_info,igfft0,ngfft)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : czero_gw


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wavefunctions_information),intent(inout) :: Wf_info
!arrays
 integer,intent(in) :: igfft0(Wf_info%npwwfn),ngfft(18)
!Local variables ------------------------------
!scalars
 integer :: istat,my_minb,my_maxb,nkibz,nsppol,nspinor,nfft
 character(len=500) :: msg
 character(len=50),parameter :: sub_name='init_wf_info_2.F90'
!************************************************************************

 Wf_info%nfftot   = ngfft(1)*ngfft(2)*ngfft(3)
 Wf_info%nfft     = Wf_info%nfftot !MG TODO this has to be fixed

 Wf_info%ngfft(:) = ngfft(:)

 allocate(Wf_info%igfft0(Wf_info%npwwfn))
 Wf_info%igfft0(:)=igfft0(:)

 select case (MODULO(Wf_info%gwmem,10))

 case (0)
  ! * Do not allocate the pointer wfr  
  write(msg,'(a)')' init2 : COMMENT- wavefunctions are NOT stored in memory! '
  call wrtout(std_out,msg,'COLL') 

 case (1)
  write(msg,'(a)')' init2 : wavefunctions are stored in memory '
  call wrtout(std_out,msg,'COLL') 
  nsppol  = Wf_info%nsppol
  nspinor = Wf_info%nspinor
  nkibz   = Wf_info%nk
  my_minb = Wf_info%my_minb
  my_maxb = Wf_info%my_maxb
  nfft    = Wf_info%nfft
  allocate(Wf_info%wfr(nfft*nspinor,my_minb:my_maxb,nkibz,nsppol),stat=istat)
  if (istat/=0) then 
   call memerr(sub_name,'wfr',nfft*nspinor*(my_maxb-my_minb+1)*nkibz*nsppol,'gwpc')
  end if
  Wf_info%wfr=czero_gw

 case default
  write(msg,'(4a)')ch10,&
&  ' init_wf_info_2 : BUG -',ch10,&
&  ' Wf_info%gwmem/=x0,x1 not yet implemented '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 END SELECT

end subroutine init_wf_info_2
!!***


!!****f* ABINIT/destroy_wf_info 
!! NAME
!!  destroy_wf_info
!!
!! FUNCTION
!!  Free the dynamic entities in the data type
!!
!! SOURCE

subroutine destroy_wf_info(Wf_info)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(wavefunctions_information),intent(inout) :: Wf_info
!************************************************************************

 if (associated(Wf_info%igfft0))            deallocate(Wf_info%igfft0)
 if (associated(Wf_info%is_already_stored)) deallocate(Wf_info%is_already_stored)
 if (associated(Wf_info%wfg))               deallocate(Wf_info%wfg)
 if (associated(Wf_info%wfr))               deallocate(Wf_info%wfr)

end subroutine destroy_wf_info
!!***


!!****f* ABINIT/reinit_wf_info
!! NAME
!!  reinit_wf_info
!!
!! FUNCTION
!!  reinitialize the storage mode
!!
!! SOURCE

subroutine reinit_wf_info(Wf_info)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(wavefunctions_information),intent(inout) :: Wf_info
!************************************************************************

 Wf_info%is_already_stored(:,:,:)=.FALSE.

end subroutine reinit_wf_info
!!***


!!****f* ABINIT/get_wfr
!! NAME
!!  get_wfr
!!
!! FUNCTION
!!  Get a wave function is real space, either by doing a FFT G-->R
!!  or by just retrieving the data already stored in the data type
!!
!! INPUTS
!!  Wf_info<wavefunctions_information>=the data type
!!  MPI_enreg=Info on the parallelism
!!  ib=bands index
!!  ik=Index of the k-point in the IBZ 
!!  is=spin index
!!
!! OUTPUT
!!  wfr(Wf_info%nfft)=the required wavefunction in real space
!!
!! SOURCE

subroutine get_wfr(Wf_info,MPI_enreg,ib,ik,is,wfr)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => get_wfr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib,ik,is
 type(wavefunctions_information),intent(inout) :: Wf_info
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 complex(gwpc),intent(out) :: wfr(Wf_info%nfft*Wf_info%nspinor)
!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp=5
 integer :: istat,npwwfn,nfft,nspinor
 character(len=500) :: msg
 character(len=50),parameter :: sub_name='get_wfr.F90'
!arrays
!************************************************************************

 if (ib>Wf_info%my_maxb.or.ib<Wf_info%my_minb) then
  write(msg,'(6a,i4,a,2i4)')ch10,&
&  ' get_wfr: BUG- ',ch10,&
&  '  Band index out of range. ',ch10,&
&  '  Requiring ib = ',ib,' while range is ',Wf_info%my_minb,Wf_info%my_maxb 
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 if (ik>Wf_info%nk) then
  write(msg,'(6a,i4,a,i4)')ch10,&
&  ' get_wfr: BUG- ',ch10,&
&  '  k-point index out of range. ',ch10,&
&  '  Requiring ik = ',ik,' while nkpt = ',Wf_info%nk
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 if (is>Wf_info%nsppol) then
  write(msg,'(6a,i4,a,i4)')ch10,&
&  ' get_wfr: BUG- ',ch10,&
&  '  spin index out of range. ',ch10,&
&  '  Requiring is = ',is,' while nsppol = ',Wf_info%nsppol 
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 ! MG Here it would be better using a pointer instead of an array
 ! If wfr is big a lot of time would be wasted to copy %wfr

 if (.not.Wf_info%is_already_stored(ib,ik,is)) then
  npwwfn = Wf_info%npwwfn
  nfft   = Wf_info%nfft
  nspinor= Wf_info%nspinor

  call fft_onewfn(Wf_info%paral_kgb,nspinor,npwwfn,nfft,Wf_info%wfg(:,ib,ik,is),wfr,&
&  Wf_info%igfft0,Wf_info%ngfft,tim_fourdp,MPI_enreg)

  if (Wf_info%gwmem==1) then
   Wf_info%wfr(:,ib,ik,is)=wfr(:)
   Wf_info%is_already_stored(ib,ik,is)=.TRUE.
   !write(*,'(a,3(i4,x),a)')' wavefunction',ib,ik,is,'stored'
  end if

 else 
  ! * wfr is_already_stored, just copy it back
  wfr(:)=Wf_info%wfr(:,ib,ik,is)
 end if 

end subroutine get_wfr
!!***

!!****f* ABINIT/duplicate_wf_info 
!! NAME
!!  duplicate_wf_info
!!
!! FUNCTION
!!  Copy the content of a wavefunctions_information data type 
!!  taking into account the spreading among processors
!!
!! INPUTS
!!  MPI_enreg=information about MPI parallelization
!!  Kmesh<BZ_mesh_type>=Structure reporting information on the k-point sampling 
!!  kcalc
!!
!! OUTPUT
!!  See sides effect
!!
!! SIDES EFFECT
!!  Wf_info=
!!  Wf_info_braket=
!!  
!!
!! SOURCE

subroutine duplicate_wf_info(MPI_enreg,Wf_info,Wf_info_braket,kcalc,Kmesh)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(in) :: MPI_enreg
 type(wavefunctions_information),intent(inout) :: Wf_info,Wf_info_braket
 type(BZ_mesh_type),intent(in) :: Kmesh
!arrays
 integer,intent(in) :: kcalc(Wf_info_braket%nk)
!Local variables ------------------------------
!scalars
 integer :: ib,ik,ikibz,is,my_maxb,my_minb,nsppol,nkcalc
 integer :: me,spaceComm,ikbs_proc_rank,ierr
 character(len=500) :: msg
 character(len=50),parameter :: sub_name='get_wfr.F90'
!arrays
!************************************************************************

 call xcomm_init  (MPI_enreg,spaceComm)
 call xme_init    (MPI_enreg,me       )
 !write(*,*) me,MPI_enreg%proc_distrb(:,:,:)

 ! For band parallelism the processors must communicate
 ! the wavefunctions where the GW corrections are required
 if (MPI_enreg%gwpara==2) then
  write(msg,'(3a,2(a,i3))')ch10,&
&  ' band parallelism : communicating wavefunctions ',ch10,&
&  ' from band = ',Wf_info_braket%my_minb,' up to band = ',Wf_info_braket%my_maxb
  call wrtout(std_out,msg,'PERS')
 end if

 nsppol  = Wf_info_braket%nsppol
 nkcalc  = Wf_info_braket%nk
 my_minb = Wf_info_braket%my_minb
 my_maxb = Wf_info_braket%my_maxb

 do is=1,nsppol
  do ik=1,nkcalc
   ikibz=Kmesh%tab(kcalc(ik))
   !write(*,*)' === me',me,ik,kcalc(ik),ikibz

   do ib=my_minb,my_maxb
    ikbs_proc_rank=MPI_enreg%proc_distrb(ikibz,ib,is)
    !write(*,*)'== me',ib,ikbs_proc_rank
    if (me==ikbs_proc_rank) Wf_info_braket%wfg(:,ib,ik,is)=Wf_info%wfg(:,ib,ikibz,is)
    if (MPI_enreg%gwpara==2) then 
     call xcast_mpi(Wf_info_braket%wfg(:,ib,ik,is),ikbs_proc_rank,spaceComm,ierr)
    end if
   end do

  end do
 end do
 call leave_test(MPI_enreg)

end subroutine duplicate_wf_info
!!***

!!****f* ABINIT/nullify_wf_info
!! NAME
!!  nullify_wf_info
!!
!! FUNCTION
!!  Nullify the dynamic arrays of the data structure
!!
!! SOURCE

subroutine nullify_wf_info(Wf_info)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(wavefunctions_information),intent(inout) :: Wf_info
!************************************************************************

 nullify(Wf_info%igfft0)
 nullify(Wf_info%is_already_stored)
 nullify(Wf_info%wfg)
 nullify(Wf_info%wfr)

end subroutine nullify_wf_info
!!***


!!****f* ABINIT/print_wavefunctions_information
!! NAME
!! print_wavefunctions_information 
!!
!! FUNCTION
!!  Print the content of a wavefunctions_information datatype
!!
!! INPUTS
!!  Wf_info<wavefunctions_information>= the datatype
!!  unitno(optional)=unit number for output
!!  prtvol(optional)=verbosity level
!!  mode_paral(optional): either "COLL" or "PERS"
!!
!! OUTPUT
!!  Only printing 
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

subroutine print_wavefunctions_information(Wf_info,unitno,prtvol,mode_paral)
    
 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wavefunctions_information),intent(in) :: Wf_info
 integer,optional,intent(in) :: unitno,prtvol
 character(len=4),optional,intent(in) :: mode_paral

!local variables-------------------------------
 integer :: verb,unt
 character(len=4) :: mode
 character(len=500) :: msg      
! *************************************************************************

 verb=0      ; if (PRESENT(prtvol))     verb=prtvol
 unt=std_out ; if (PRESENT(unitno))     unt=unitno
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 write(msg,'(2a,3(a,i5,a),a,i5)')&
& ' ==== Content of the Wf_info datatype ==== ',ch10,&
& '   Number of irreducible k-points ........ ',Wf_info%nk,ch10,&
& '   Number of spinorial components ........ ',Wf_info%nspinor,ch10,&
& '   Number of spin-density components ..... ',Wf_info%nspden,ch10,&
& '   Number of spin polarizations .......... ',Wf_info%nsppol
 call wrtout(unt,msg,mode)
 write(msg,'(5(a,i5,a),a,2i5)')&
& '   Number of reciprocal lattice vectors .. ',Wf_info%npwwfn,ch10,&
& '   Memory storage option (gwmem) ......... ',Wf_info%gwmem,ch10,&
& '   Total number of FFT points ....... .... ',Wf_info%nfftot,ch10,&
& '   Number of FFT points treated by me .... ',Wf_info%nfft,ch10,&
& '   Parallelism over k-b-g (paral_kgb) .... ',Wf_info%paral_kgb,ch10,&
& '   min and Max band index stored by me ... ',Wf_info%my_minb,Wf_info%my_maxb
 call wrtout(unt,msg,mode)

 call print_ngfft(Wf_info%ngfft,"FFT mesh for wavefunctions",unt,mode,verb)

end subroutine print_wavefunctions_information 
!!***


!!****f* ABINIT/rotate_wfg
!! NAME
!! rotate_wfg
!!
!! FUNCTION
!! Rotate the Fourier components of a wave function in reciprocal space to obtain
!! the wave function at a symmetric k-point (assuming a nondegenerate state)
!!
!! INPUTS
!! Wf<wavefunctions_information>=datatype gathering information of wavefunctions in the IBZ
!! ik_bz=the BZ k-point asked for
!! iband=the required band index
!! isppol=the spin polarization
!! grottbm1(Wf%npwwfn,Kmesh%timreversal,Kmsh%nsym)=Index of (IS)^{-1}G 
!! phmgt(Wf%npwwfn,Kmesh%nsym)=phase e^{-iG.t}
!!
!! OUTPUT
!! wfg_rot(Wf%npwwfn)=the Fourier coefficients of the wave function at point ik_bz
!! Notes that possible umklapp G0 vectors are not yet supported
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine rotate_wfg(Wfs,Kmesh,iband,ik_bz,isppol,grottbm1,phmGt,wfg_rot)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15gw, except_this_one => rotate_wfg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,iband,isppol
 type(wavefunctions_information),intent(in) :: Wfs
 type(bz_mesh_type),intent(in) :: Kmesh
!arrays
 integer,target,intent(in) :: grottbm1(:,:,:)
 complex(gwpc),intent(in) :: phmGt(:,:)
 complex(gwpc),intent(out) :: wfg_rot(Wfs%npwwfn)
!Local variables ------------------------------
!scalars
 integer :: isym,itim,ik_ibz,ig
 real(dp) :: kbz(3)
 complex(gwpc) :: phnons_k 
!arrays
 integer,pointer :: Sm1G(:)
 complex(gwpc),pointer :: wfg_irr(:)
!************************************************************************

#if defined DEBUG_MODE
 if (Wfs%npwwfn/=SIZE(grottbm1,DIM=1))    STOP ' rotate_wfg : Wf%npwwfn /= grottbm1 DIM 1 '
 if (Kmesh%timrev/=SIZE(grottbm1,DIM=2)) STOP ' rotate_wfg : timrev /= grottbm1 DIM 2 ' 
 if (Kmesh%nsym/=SIZE(grottbm1,DIM=3))   STOP ' rotate_wfg : nsym /= grottbm1 DIM 3 ' 
 if (Wfs%npwwfn/=SIZE(phmGt,DIM=1))      STOP ' rotate_wfg : Wf%npwwfn /= phmGt DIM 1 '
 if (Kmesh%nsym/=SIZE(phmGt,DIM=2))     STOP ' rotate_wfg : nsym /= phmGt DIM 2 ' 
#endif
 
 ! ****WARNING WARNING ****
 !FIXME Umklapp not yet implemented I have to store the umklapp somewhere
 ! but it requires boring modification of the code
 STOP "Still under development"

 call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym,itim,phnons_k)
 write(*,*) " Rotate_wfg WARNING Umklapp not implemented yet "
 !
 ! === Rotate irreducible wave function ===
 ! $u_{Sk}(G) = e^{-i(Sk+G).\tau} u_k(S^{-1}G)$
 ! $u_{-k}(G) = u_k(-G)^*$
 ! $u_{k+G0}(G) = u_k(G+G0)$   
 !
 Sm1G => grottbm1(:,itim,isym)
 wfg_irr => Wfs%wfg(:,iband,ik_ibz,isppol)
 do ig=1,Wfs%npwwfn
  wfg_rot(ig)=wfg_irr(Sm1G(ig))*phmGt(ig,isym)
 end do
 wfg_rot=wfg_rot*phnons_k
 if (itim==2) wfg_rot=CONJG(wfg_rot)

end subroutine rotate_wfg
!!***


!!****f* ABINIT/rotate_wfr
!! NAME
!! rotate_wfr
!! 
!! FUNCTION
!! Rotate the periodic part of a wave function in real space to obtain
!! the lattice-period part of the Bloch function at a symmetric k-point (assuming a nondegenerate state)
!!
!! INPUTS
!! Wfs<wavefunctions_information>=datatype gathering information of wavefunctions in the IBZ
!! ik_bz=the BZ k-point asked for
!! iband=the required band index
!! isppol=the spin polarization
!! irottb(Wf%nfft,Kmesh%timreversal,Kmsh%nsym)=Index of (IS)^{-1}G 
!! phmgt(Wf%npwwfn,Kmesh%nsym)=phase e^{-iG.t}
!!
!! OUTPUT
!! wfr_rot(Wf%npwwfn)=the Fourier coefficients of the wave function at point ik_bz
!! Notes that possible umklapp G0 vectors are not yet supported
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine rotate_wfr(Wfs,Kmesh,iband,ik_bz,isppol,irottb,MPI_enreg,ur_rot)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15gw, except_this_one => rotate_wfr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,iband,isppol
 type(Wavefunctions_information),intent(inout) :: Wfs
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in),target :: irottb(Wfs%nfft,Kmesh%nsym)
 complex(gwpc),intent(out) :: ur_rot(Wfs%nfft)

!Local variables ------------------------------
!scalars
 integer :: isym,itim,ik_ibz,ir
 real(dp) :: kbz(3)
 complex(gwpc) :: ph_mkbzt
!arrays
 integer,pointer :: Rm1rt(:)
 complex(gwpc),allocatable :: ur_irr(:)
!************************************************************************

 STOP "Still under development"
 
 call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym,itim,ph_mkbzt)
 !
 ! === Rotate irreducible wave function ===
 ! $u_{Sk}(r) = e^{-i(Sk.\tau} u_k(R^{-1}(r-\tau))$
 ! $u_{-k}(r) = u_k(r)^*$
 !
 allocate(ur_irr(Wfs%nfft))
 call get_wfr(Wfs,MPI_enreg,iband,ik_ibz,isppol,ur_irr)
 Rm1rt => irottb(:,isym)

 do ir=1,Wfs%nfft
  ur_rot(ir)=ur_irr(Rm1rt(ir)) 
 end do
 ur_rot(:)=ur_rot(:)*ph_mkbzt
 if (itim==2) ur_rot=CONJG(ur_rot)

 deallocate(ur_irr)

end subroutine rotate_wfr
!!***
