!{\src2tex{textfont=tt}}
!!****f* ABINIT/destroy_paw_ij
!! NAME
!! destroy_paw_ij
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  paw_ij(:)<type(paw_ij_type)>=paw arrays given on (i,j) channels
!!
!! OUTPUT
!!  See side effects
!!
!! SIDE EFFECTS
!!  All associated pointers in paw_ij(:) are deallocated
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

subroutine destroy_paw_ij(Paw_ij)
    
 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(Paw_ij_type),intent(inout) :: Paw_ij(:)

!Local variables-------------------------------
 integer :: iat,natom
 character(len=500) :: msg      
 
! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' destroy_paw_ij : enter '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 natom=SIZE(Paw_ij)
 do iat=1,natom
  if (ASSOCIATED(Paw_ij(iat)%dij       )) deallocate(Paw_ij(iat)%dij       )
  if (ASSOCIATED(Paw_ij(iat)%dijhartree)) deallocate(Paw_ij(iat)%dijhartree)
  if (ASSOCIATED(Paw_ij(iat)%dijhat    )) deallocate(Paw_ij(iat)%dijhat    )
  if (ASSOCIATED(Paw_ij(iat)%dijU      )) deallocate(Paw_ij(iat)%dijU      )
  if (ASSOCIATED(Paw_ij(iat)%dijso     )) deallocate(Paw_ij(iat)%dijso     )
  if (ASSOCIATED(Paw_ij(iat)%dijxc     )) deallocate(Paw_ij(iat)%dijxc     )
  if (ASSOCIATED(Paw_ij(iat)%dijxc_val )) deallocate(Paw_ij(iat)%dijxc_val )
  if (ASSOCIATED(Paw_ij(iat)%noccmmp   )) deallocate(Paw_ij(iat)%noccmmp   )
  if (ASSOCIATED(Paw_ij(iat)%nocctot   )) deallocate(Paw_ij(iat)%nocctot   )
  if (ASSOCIATED(Paw_ij(iat)%vpawx     )) deallocate(Paw_ij(iat)%vpawx     )
 end do

#if defined DEBUG_MODE
 write(msg,'(a)')' destroy_paw_ij : exit '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine destroy_Paw_ij
!!***

!!****f* ABINIT/init_paw_ij
!! NAME
!! init_paw_ij
!!
!! FUNCTION
!!  Initialize a Paw_ij data type.
!!
!! INPUTS
!!  cplex=1 if all on-site PAW quantities are real, 2 if they are complex
!!  cplex_dij=1 if dij are real, 2 if they are complex
!!  natom=Number of atoms.
!!  ntypat=Number of types of atoms in cell.
!!  nspinor=Number of spinor components
!!  nsppol=Number of independent spin polarizations.
!!  nspden=Number of spin-density components
!!  pawspnorb=1 if spin-orbit coupling is activated
!!  typat(natom)=Type of each atom
!!  Pawtab(ntypat)<type(pawtab_type)>=PAW tabulated starting data
!!
!! OPTIONAL INPUTS
!!  has_dijhat=1 to allocate Paw_ij%dijhat, 0 otherwise (default)
!!  has_dijxc=1 to allocate Paw_ij%dijxc, 0 otherwise (default)
!!  has_dijxc_val=1 to allocate Paw_ij%dijxc_val, 0 otherwise (default)
!!  has_dijhartree=1 to allocate Paw_ij%dijhartree, 0 otherwise (default)
!!  has_dijso=1 to allocate Paw_ij%dijso, used only if pawspnorb>0. 0 otherwise (default) 
!!  has_dijU=1 to allocate Paw_ij%dijU, used only if Pawtab(itypat)%usepawu>0. 0 otherwise (default).
!!
!! OUTPUT
!!  Paw_ij(natom)<type(paw_ij_type)>=data structure containing PAW arrays given on (i,j) channels.
!!   In output all the basic dimensions are defined and the arrays are allocated  
!!   according to the input variables.
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

subroutine init_paw_ij(Paw_ij,cplex,cplex_dij,nspinor,nsppol,nspden,pawspnorb,natom,ntypat,typat,Pawtab,&
& has_dijhartree,has_dijhat,has_dijxc,has_dijxc_val,has_dijso,has_dijU) ! Optional
    
 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,cplex_dij,nspinor,nspden,nsppol,pawspnorb,natom,ntypat
 integer,optional,intent(in) :: has_dijhat,has_dijxc,has_dijxc_val 
 integer,optional,intent(in) :: has_dijso,has_dijhartree,has_dijU
!arrays
 integer,intent(in) :: typat(natom)
 type(Paw_ij_type),intent(inout) :: Paw_ij(natom)
 type(Pawtab_type),intent(in) :: Pawtab(ntypat)

!Local variables-------------------------------
 integer :: iat,itypat,lmn2_size,ndij
 character(len=500) :: msg      
 
! *************************************************************************

 !allocate(Paw_ij(natom))
 !call nullify_paw_ij(Paw_ij)

 do iat=1,natom
  itypat=typat(iat)
  lmn2_size              =Pawtab(itypat)%lmn2_size
  Paw_ij(iat)%cplex      =cplex
  !Paw_ij(iat)%cplex_dij  =nspinor
  !Paw_ij(iat)%cplex_dij  =MAX(cplex,1+pawspnorb,nspinor)
  Paw_ij(iat)%cplex_dij  =cplex_dij
  Paw_ij(iat)%nspden     =nspden
  Paw_ij(iat)%nsppol     =nsppol
  Paw_ij(iat)%lmn_size   =Pawtab(itypat)%lmn_size
  Paw_ij(iat)%lmn2_size  =lmn2_size
  Paw_ij(iat)%ndij       =MAX(nspinor**2,nspden)
  !Paw_ij(iat)%lmnmix_sz =  do we need this? It seems it is not used anymore and can be removed

  !cplex_dij=Paw_ij(iat)%cplex_dij
  ndij     =Paw_ij(iat)%ndij

  allocate(Paw_ij(iat)%dij(cplex_dij*lmn2_size,ndij))

  ! === Allocation for PAW+U ===
  if (Pawtab(itypat)%usepawu>0) then
   allocate(Paw_ij(iat)%noccmmp(2*Pawtab(itypat)%lpawu+1,2*Pawtab(itypat)%lpawu+1,nspden))
   allocate(Paw_ij(iat)%nocctot(nspden))
  end if

  ! === Allocation for PAW+LEXX ===
  if (Pawtab(itypat)%useexexch>0) then
   ! TODO solve issue with first dimension
   allocate(Paw_ij(iat)%vpawx(1,lmn2_size,nspden))
  end if                                                                                    

  ! ============================
  ! === Optional allocations ===
  ! ============================
  Paw_ij(iat)%has_dijhartree=0
  if (PRESENT(has_dijhartree)) then
   if (has_dijhartree/=0) then
    Paw_ij(iat)%has_dijhartree=1
    allocate(Paw_ij(iat)%dijhartree(cplex*lmn2_size))
   end if
  end if

  Paw_ij(iat)%has_dijhat=0 
  if (PRESENT(has_dijhat)) then
   if (has_dijhat/=0) then
    Paw_ij(iat)%has_dijhat=1
    allocate(Paw_ij(iat)%dijhat(cplex_dij*lmn2_size,ndij))
   end if
  end if

  Paw_ij(iat)%has_dijxc=0
  if (PRESENT(has_dijxc)) then
   if (has_dijxc/=0) then
    Paw_ij(iat)%has_dijxc=1
    allocate(Paw_ij(iat)%dijxc(cplex_dij*lmn2_size,ndij))
   end if
  end if

  Paw_ij(iat)%has_dijxc_val=0
  if (PRESENT(has_dijxc_val)) then
   if (has_dijxc_val/=0) then
    Paw_ij(iat)%has_dijxc_val=1
    allocate(Paw_ij(iat)%dijxc_val(cplex_dij*lmn2_size,ndij))
   end if
  end if

  Paw_ij(iat)%has_dijU=0
  if (PRESENT(has_dijU)) then
   if (has_dijU/=0.and.Pawtab(itypat)%usepawu>0) then
    Paw_ij(iat)%has_dijU=1
    allocate(Paw_ij(iat)%dijU(cplex_dij*lmn2_size,ndij))
   end if
  end if

  Paw_ij(iat)%has_dijso=0
  if (PRESENT(has_dijso)) then
   if (has_dijso/=0.and.pawspnorb>0) then
    Paw_ij(iat)%has_dijso=1
    allocate(Paw_ij(iat)%dijso(cplex_dij*lmn2_size,ndij))
   end if
  end if

 end do !iat 

end subroutine init_paw_ij
!!***


!!****f* ABINIT/destroy_paw_an
!! NAME
!! destroy_paw_an
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Paw_an(:)<type(Paw_an_type)>=various arrays given on ANgular mesh or ANgular moments
!!
!! OUTPUT
!!  See side effects
!!
!! SIDE EFFECTS
!!  All associated pointers in Paw_an(:) are deallocated
!!
!! NOTES
!!  vh1 and vht1 are defined in the data structure but never used. 
!!  Cannot test for association status since these quantities are 
!!  not nulliified before entering the calculation
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

subroutine destroy_Paw_an(Paw_an)
    
 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!arrays
 type(Paw_an_type),intent(inout) :: Paw_an(:)

!Local variables-------------------------------
 integer :: iat,natom,itypat
 character(len=500) :: msg      
 
! *************************************************************************

 natom=SIZE(Paw_an)

! do iat=1,natom
!  itypat=typat(iat)
!  deallocate(Paw_an(iat)%lmselect)
!  !deallocate(Paw_an(iat)%vh1,Paw_an(iat)%vht1)      !TODO nullify these arrays
!  deallocate(paw_an(iat)%vxc1,Paw_an(iat)%vxct1)
!  if (Paw_an(iat)%has_vxcval==1 ) deallocate(Paw_an(iat)%vxc1_val,Paw_an(iat)%vxct1_val)
!  if (Pawtab(itypat)%useexexch>0) deallocate(Paw_an(iat)%vxc_ex)
! end do

 do iat=1,natom
  if (ASSOCIATED(Paw_an(iat)%lmselect )) deallocate(Paw_an(iat)%lmselect )
  if (ASSOCIATED(Paw_an(iat)%vh1      )) deallocate(Paw_an(iat)%vh1      )
  if (ASSOCIATED(Paw_an(iat)%vht1     )) deallocate(Paw_an(iat)%vht1     )
  if (ASSOCIATED(Paw_an(iat)%vxc1     )) deallocate(Paw_an(iat)%vxc1     )
  if (ASSOCIATED(Paw_an(iat)%vxc1_val )) deallocate(Paw_an(iat)%vxc1_val )
  if (ASSOCIATED(Paw_an(iat)%vxct1    )) deallocate(Paw_an(iat)%vxct1    )
  if (ASSOCIATED(Paw_an(iat)%vxct1_val)) deallocate(Paw_an(iat)%vxct1_val)
  if (ASSOCIATED(Paw_an(iat)%vxc_ex   )) deallocate(Paw_an(iat)%vxc_ex   )
 end do !iat 

end subroutine destroy_Paw_an 
!!***

!!****f* ABINIT/init_pawfgr
!! NAME
!! init_pawfgr
!!
!! FUNCTION
!!  Initialize a pawfgr_type datatype, reporting also info on the FFT mesh 
!!  according to the method used (norm-conserving or PAW)
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! k0(3)=input k vector for k+G sphere 
!! Dtset <type(dataset_type)>=all input variables for this dataset
!!   %dilatmx
!!   %usepaw
!!   %natom
!!   %ngfft
!!   %ngfftdg
!!   %nfft
!!   %mgfft
!!   %mgfftdg
!!   %dilatmx
!!   %pawecutdg 
!!   %ecut
!!
!! OUTPUT
!!  ecut_eff=effective energy cutoff (hartree) for coarse planewave basis sphere
!!  ecutdg_eff=effective energy cutoff (hartree) for dense planewave basis sphere
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  mgfft=maximum size of 1D FFTs
!1  ngfftc(18),ngfftf(18)=contain all needed information about 3D FFT, for coarse and dense FFT mesh, respectively.
!!   see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  Pawfgr<pawfgr_type>=For PAW, Fine rectangular GRid parameters and related data
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  is not used and not initialized inside this routine
!!
!! PARENTS
!!  
!!
!! CHILDREN
!!  
!!
!! SOURCE

subroutine init_pawfgr(Dtset,k0,gmet,Pawfgr,mgfftf,nfftf,ecut_eff,ecutdg_eff,&
& gsqcutc_eff,gsqcutf_eff,ngfftc,ngfftf)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13paw
 use interfaces_13recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: nfftf,mgfftf 
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(out) :: ecut_eff,ecutdg_eff
 real(dp),intent(out) :: gsqcutf_eff,gsqcutc_eff
 type(dataset_type),intent(in) :: Dtset
 type(pawfgr_type),intent(out) :: Pawfgr
!arrays
 integer,intent(out) :: ngfftc(18),ngfftf(18)
 real(dp),intent(in) :: k0(3)

!Local variables-------------------------------
 integer :: ii,nfftc_tot,nfftf_tot
 real(dp) :: boxcut,boxcutc
 character(len=500) :: msg

!************************************************************************

 ngfftc(:)=Dtset%ngfft(:)

 SELECT CASE (Dtset%usepaw)
 CASE (0)
  ! === Norm-conserving pseudopotentials ===
  nfftf=Dtset%nfft ; mgfftf=Dtset%mgfft ; ngfftf(:)=Dtset%ngfft(:)
  Pawfgr%usefinegrid=0 ; allocate(Pawfgr%coatofin(0),Pawfgr%fintocoa(0))
  ecut_eff  =Dtset%ecut*Dtset%dilatmx**2
  ecutdg_eff=ecut_eff
 CASE (1)
  ! == PAW calculation ===
  if (Dtset%pawecutdg>=1.0000001_dp*Dtset%ecut) then
   ! === Use fine FFT grid generated according to pawecutdg ===
   nfftf=Dtset%nfftdg ; mgfftf=Dtset%mgfftdg ; ngfftf(:)=Dtset%ngfftdg(:)
   nfftc_tot =ngfftc(1)*ngfftc(2)*ngfftc(3)
   nfftf_tot =ngfftf(1)*ngfftf(2)*ngfftf(3)
   Pawfgr%usefinegrid=1 ; allocate(Pawfgr%coatofin(nfftc_tot),Pawfgr%fintocoa(nfftf_tot))
   call indgrid(Pawfgr%coatofin,Pawfgr%fintocoa,nfftc_tot,nfftf_tot,ngfftc,ngfftf)
  else 
   ! === Do not use fine FFT mesh. Simple transfer that can be done in parallel with only local info ===
   nfftf=Dtset%nfft ; mgfftf=Dtset%mgfft ; ngfftf(:)=Dtset%ngfft(:)
   Pawfgr%usefinegrid=0 ; allocate(Pawfgr%coatofin(Dtset%nfft),Pawfgr%fintocoa(Dtset%nfft))
   do ii=1,Dtset%nfft
    Pawfgr%coatofin(ii)=ii ; Pawfgr%fintocoa(ii)=ii
   end do
  end if
  ! == Store useful dimensions in Pawfgr ===
  Pawfgr%natom=Dtset%natom
  Pawfgr%nfftc=Dtset%nfft ; Pawfgr%mgfftc=Dtset%mgfft ; Pawfgr%ngfftc(:)=Dtset%ngfft(:)
  Pawfgr%nfft=nfftf       ; Pawfgr%mgfft=mgfftf       ; Pawfgr%ngfft (:)=ngfftf(:)
  ecutdg_eff=Dtset%pawecutdg*Dtset%dilatmx**2
  ecut_eff  =Dtset%ecut*Dtset%dilatmx**2  
 CASE DEFAULT
  write(msg,'(4a,i4)')ch10,&
&  ' init_pawfgr : BUG ',ch10,&
&  ' called with wrong value of usepaw: ',Dtset%usepaw
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 END SELECT
 !
 ! === Get boxcut for given gmet, ngfft, and ecut (center at k0) ===
 !     boxcut=ratio of basis sphere diameter to fft box side
 if (Dtset%usepaw==1) then
  write(msg,'(2a)')ch10,' Coarse grid specifications '!(used for wave-functions):'
  call wrtout(std_out,msg,'COLL') !; call wrtout(ab_out,msg,'COLL') 
  call getcut(boxcutc,ecut_eff,gmet,gsqcutc_eff,Dtset%iboxcut,std_out,k0,ngfftc)
  write(msg,'(2a)')ch10,' Fine grid specifications (used for densities):'
  call wrtout(std_out,msg,'COLL') !; call wrtout(ab_out,msg,'COLL') 
  call getcut(boxcut,ecutdg_eff,gmet,gsqcutf_eff,Dtset%iboxcut,std_out,k0,ngfftf)
  !FIXME this is never use, should ask Marc 
  !Pawfgr%gsqcut=gsqcutf_eff 
 else
  call getcut(boxcut,ecut_eff,gmet,gsqcutc_eff,Dtset%iboxcut,std_out,k0,ngfftc)
  gsqcutf_eff=gsqcutc_eff
 end if
 !
 ! === Check that boxcut>=2 if intxc=1; otherwise intxc must be set=0 ===
 if (boxcut<two .and. Dtset%intxc==1) then
  write(msg,'(4a,es12.4,5a)')ch10,&
&  ' init_pawfgr : ERROR -',ch10,&
&  ' boxcut=',boxcut,' is < 2.0  => intxc must be 0;',ch10,&
&  ' Need larger ngfft to use intxc=1.',ch10,&
&  ' Action : you could increase ngfft, or decrease ecut, or put intxc=0.'
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

end subroutine init_pawfgr
!!***

!!****f* ABINIT/nullify_paw_ij
!! NAME
!!  nullify_paw_ij
!!
!! FUNCTION
!!  Nullify the pointers in a paw_ij data type
!!
!! INPUTS
!!  See SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  Paw_ij(:)<type(paw_ij_type)>=PAW arrays given on (i,j) channels. Nullified in output
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

subroutine nullify_paw_ij(Paw_ij)
    
 use defs_basis
 use defs_datatypes 

 implicit none

!Arguments ------------------------------------
!arrays
 type(Paw_ij_type),intent(inout) :: Paw_ij(:)

!Local variables-------------------------------
 integer :: iat,natom 
 
! *************************************************************************

 natom=SIZE(Paw_ij(:))

 do iat=1,natom
  nullify(Paw_ij(iat)%dij)
  nullify(Paw_ij(iat)%dijhartree)
  nullify(Paw_ij(iat)%dijhat)
  nullify(Paw_ij(iat)%dijU)
  nullify(Paw_ij(iat)%dijso)
  nullify(Paw_ij(iat)%dijxc )
  nullify(Paw_ij(iat)%dijxc_val)
  nullify(Paw_ij(iat)%noccmmp)
  nullify(Paw_ij(iat)%nocctot)
  nullify(Paw_ij(iat)%vpawx)
 end do !iat 

end subroutine nullify_paw_ij
!!***

!!****f* ABINIT/nullify_paw_an
!! NAME
!!  nullify_paw_an
!!
!! FUNCTION
!!  Nullify the pointers in a paw_an data type
!!
!! INPUTS
!!  See SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  Paw_an(:)<type(paw_an_type)>=PAW arrays given on ANgular mesh or ANgular moments.
!!  Nullified in output
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

subroutine nullify_paw_an(Paw_an)
    
 use defs_basis
 use defs_datatypes 

 implicit none

!Arguments ------------------------------------
!arrays
 type(Paw_an_type),intent(inout) :: Paw_an(:)

!Local variables-------------------------------
 integer :: iat,natom 
 
! *************************************************************************

 natom=SIZE(Paw_an(:))

 do iat=1,natom
  nullify(Paw_an(iat)%lmselect )
  nullify(Paw_an(iat)%vh1      )
  nullify(Paw_an(iat)%vht1     )
  nullify(Paw_an(iat)%vxc1     )
  nullify(Paw_an(iat)%vxc1_val )
  nullify(Paw_an(iat)%vxct1_val)
  nullify(Paw_an(iat)%vxc_ex   )
 end do !iat 

end subroutine nullify_paw_an
!!***

!!****f* ABINIT/init_paw_an
!! NAME
!!  init_paw_an
!!
!! FUNCTION
!!  Initialize a paw_an data type. 
!!
!! INPUTS
!!  See SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  Paw_an(:)<type(paw_an_type)>=PAW arrays given on ANgular mesh or ANgular moments.
!!  Initialized in output
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

subroutine init_paw_an(natom,ntypat,nspden,cplex,pawxcdev,pawspnorb,typat,Pawang,Pawtab,Paw_an,&
& has_vxcval) ! Optional
    
 use defs_basis
 use defs_datatypes 

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat,cplex,nspden,pawspnorb,pawxcdev
 integer,optional,intent(in) :: has_vxcval
!arrays
 integer,intent(in) :: typat(natom)
 type(Pawang_type),intent(in) :: Pawang
 type(Pawtab_type),intent(in) :: Pawtab(ntypat)
 type(Paw_an_type),intent(inout) :: Paw_an(:)

!Local variables-------------------------------
 integer :: iat,itypat,lm_size,v_size
 
! *************************************************************************

 !allocate(Paw_an(natom))
 !call nullify_paw_an(Paw_an)

 do iat=1,natom
  itypat=typat(iat)
  lm_size                =Pawtab(itypat)%lcut_size**2
  Paw_an(iat)%angl_size  =Pawang%angl_size
  Paw_an(iat)%cplex      =cplex
  Paw_an(iat)%lm_size    =lm_size
  Paw_an(iat)%mesh_size  =Pawtab(itypat)%mesh_size
  Paw_an(iat)%nspden     =nspden

  ! Non-zero LM-moments of "one-center" densities/potentials.
  ! * Filled in pawdenpot.
  allocate(Paw_an(iat)%lmselect(lm_size)) 

  ! xc potential inside the sphere.
  !  * LM-moments of potential if pawxcdev/=0
  !  * (theta,phi) values of potential if pawxcdev=0
  v_size=Paw_an(iat)%lm_size ; if (pawxcdev==0) v_size=Paw_an(iat)%angl_size
  allocate(Paw_an(iat)%vxc1 (cplex*Pawtab(itypat)%mesh_size,v_size,nspden))
  allocate(Paw_an(iat)%vxct1(cplex*Pawtab(itypat)%mesh_size,v_size,nspden))
 
  ! ==========================
  ! === Optional arguments ===
  ! ==========================

  ! xc potential inside PAW spheres generated by valence electrons.
  Paw_an(iat)%has_vxcval=0
  if (PRESENT(has_vxcval)) then 
   if (has_vxcval>0) then
    Paw_an(iat)%has_vxcval=1
    allocate(Paw_an(iat)%vxc1_val (cplex*Pawtab(itypat)%mesh_size,v_size,nspden))
    allocate(Paw_an(iat)%vxct1_val(cplex*Pawtab(itypat)%mesh_size,v_size,nspden))
   end if
  end if

  ! xc potential for local exact exchange inside the sphere.
  if (Pawtab(itypat)%useexexch>0) then 
   allocate(Paw_an(iat)%vxc_ex(cplex*Pawtab(itypat)%mesh_size,v_size,nspden))
  end if

  ! Hartree potential LM-moments inside the sphere
  ! what about vht1?
  if (pawspnorb>0) then 
   allocate(Paw_an(iat)%vh1(cplex*Pawtab(itypat)%mesh_size,1,1))
  end if

 end do !iat

end subroutine init_paw_an
!!***
