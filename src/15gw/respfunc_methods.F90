!{\src2tex{textfont=tt}}
!!****f* ABINIT/respfunc_methods
!! NAME
!! respfunc_methods
!!
!! FUNCTION
!!  This module contains methods acting on the Epsilonm1_parameters data type
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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

!!***
!!****f* ABINIT/nullify_epsilonm1_results
!! NAME
!! nullify_epsilonm1_results
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine nullify_epsilonm1_results(Er)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(Epsilonm1_results),intent(inout) :: Er
! *************************************************************************
 
 nullify(Er%gvec)

 nullify(Er%qibz)
 nullify(Er%qlwl) 

 nullify(Er%epsm1) 
 nullify(Er%lwing)
 nullify(Er%omega)             
 nullify(Er%uwing)       

 !call nullify(Er%Hdr) Hdr is always read but nullification might be useful!

end subroutine nullify_epsilonm1_results
!!***

!!****f* ABINIT/destroy_epsilonm1_results
!! NAME
!! destroy_epsilonm1_results
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine destroy_epsilonm1_results(Er)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_14iowfdenpot
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Epsilonm1_results),intent(inout) :: Er
!Local variables-------------------------------
! *************************************************************************

 if (ASSOCIATED(Er%gvec )) deallocate(Er%gvec)

 if (ASSOCIATED(Er%qibz )) deallocate(Er%qibz) 
 if (ASSOCIATED(Er%qlwl )) deallocate(Er%qlwl) 

 if (ASSOCIATED(Er%epsm1)) deallocate(Er%epsm1) 
 if (ASSOCIATED(Er%lwing)) deallocate(Er%lwing)
 if (ASSOCIATED(Er%omega)) deallocate(Er%omega)             
 if (ASSOCIATED(Er%uwing)) deallocate(Er%uwing)       

 call hdr_clean(Er%Hdr)

end subroutine destroy_epsilonm1_results
!!***

!!****f* ABINIT/print_epsilonm1_results
!! NAME
!!  print_epsilonm1_results
!!
!! FUNCTION
!!  print basic dimensions and most important quantities contained 
!!  in a Epsilonm1_results data type
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine print_epsilonm1_results(Er,unit,prtvol,mode_paral)

 use defs_basis
 use defs_datatypes
 use m_numeric_tools, only : print_arr


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 type(Epsilonm1_results),intent(in) :: Er

!Local variables-------------------------------
 integer :: iomega,iqibz,iqlwl,unt,verbose,rdwr
 character(len=100) :: fmt
 character(len=50) :: rfname,rforder,rfapprox,rftest
 character(len=500) :: msg
 character(len=4) :: mode
! *************************************************************************

 unt=std_out ; if (PRESENT(unit))       unt=unit
 verbose=0   ; if (PRESENT(prtvol))     verbose=prtvol
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 ! === chi0 or \epsilon^{-1} ? ===
 if (Er%ID==0) then
  rfname='Symmetrical Inverse Dielectric Matrix'
 else if (Er%ID==1) then
  rfname='Irreducible Polarizability'
 else  
  write(msg,'(4a,i3)')ch10,&
&  ' print_epsilonm1_results : BUG - ',ch10,&
&  ' Wrong value of Er%ID = ',Er%ID
  call wrtout(unt,msg,'COLL') 
  call leave_new('COLL')
 end if

 ! === If \epsilon^{-1}, define the approximation ===
 rfapprox='None'
 if (Er%ID==0) then
  if (Er%ikxc==0) then 
   rfapprox='RPA'
  else if (Er%ikxc>0) then 
   rfapprox='Static TDDFT'
  else  
   rfapprox='TDDFT'
  end if
 end if

 ! === If TDDFT and \epsilon^{-1}, define the type ===
 rftest='None'
! if (Er%ID==0) then
!  if (Er%test_type==0) then 
!   rftest='TEST-PARTICLE'
!  else if (Er%test_type==1) then 
!   rftest='TEST-ELECTRON'
!  else 
!   write(msg,'(4a,i3)')ch10,&
!&   ' print_epsilonm1_results : BUG - ',ch10,&
!&   ' Wrong value of Er%test_type = ',Er%test_type
!   call wrtout(unt,msg,'COLL') 
!   call leave_new('COLL')
!  end if
! end if

 ! === Define time-ordering ===
 if (Er%Tordering==0) then 
  rforder='Time-Ordered'
 else if (Er%Tordering==1) then 
  rforder='Advanced'
 else if (Er%Tordering==2) then 
  rforder='Retarded'
 else  
  write(msg,'(4a,i3)')ch10,&
&  ' print_epsilonm1_results : BUG - ',ch10,&
&  ' Wrong value of Er%Tordering = ',Er%Tordering
  call wrtout(unt,msg,'COLL') 
  call leave_new('COLL')
 end if

 write(msg,'(3a,4(3a))')ch10,&
& ' ==== Info on the Response Function ==== ',ch10,&
& '  Response Function Type .......... ',TRIM(rfname),ch10,&
& '  Type of Approximation ........... ',TRIM(rfapprox),ch10,&
& '  Type of probing particle ........ ',TRIM(rftest),ch10,&
& '  Time-Ordering ................... ',TRIM(rforder),ch10
 call wrtout(unt,msg,mode) 
 write(msg,'(a,2i4,a,3(a,i4,a),a,3i4,a,a,i4,a)')&
& '  Number of components ............ ',Er%nI,Er%nJ,ch10,&
& '  Number of q-points in the IBZ ... ',Er%nqibz,ch10,&
& '  Number of q-points for q-->0      ',Er%nqlwl,ch10,&
& '  Number of G-vectors ............. ',Er%npwe,ch10,&
& '  Number of frequencies ........... ',Er%nomega,Er%nomega_r,Er%nomega_i,ch10,&
& '  Value of mqmem .................. ',Er%mqmem,ch10
 call wrtout(unt,msg,mode) 

 if (Er%nqlwl/=0) then 
  write(msg,'(a,i3)')' q-points for long wavelength limit: ',Er%nqlwl
  call wrtout(unt,msg,mode) 
  do iqlwl=1,Er%nqlwl
   write(msg,'(1x,i5,a,3es16.8)')iqlwl,') ',Er%qlwl(:,iqlwl)
   call wrtout(unt,msg,mode)
  end do
 end if

 if (verbose>0) then
  ! === If chi0, write out head and wings in the long-wavelenght limit ===
  if (Er%ID==1) then 
   write(msg,'(1x,2a)')' Heads and wings of chi0(G,G'')',ch10
   call wrtout(unt,msg,mode)
   do iqlwl=1,Er%nqlwl
    write(msg,'(1x,a,i2,a)')' chi0(qlwl =',iqlwl,')'
    call wrtout(unt,msg,mode)
    do iomega=1,Er%nomega
     write(msg,'(2x,a,i4,a,2f9.4,a)')&
&     ' Upper and lower wings at the ',iomega,' th omega',Er%omega(iomega)*Ha_eV,' [eV]'
     call wrtout(unt,msg,mode)
     call print_arr(Er%uwing(:,iomega,iqlwl),max_r=9,unit=unt)
     call print_arr(Er%lwing(:,iomega,iqlwl),max_r=9,unit=unt)
    end do
   end do
  end if

  write(msg,'(a,i4)')' Calculated Frequencies: ',Er%nomega
  call wrtout(unt,msg,mode) 
  do iomega=1,Er%nomega 
   write(msg,'(i4,es14.6)')iomega,Er%omega(iomega)*Ha_eV
   call wrtout(unt,msg,mode) 
  end do

  write(msg,'(a,i4)')' Calculated q-points: ',Er%nqibz
  call wrtout(unt,msg,mode) 
  do iqibz=1,Er%nqibz
   write(msg,'(1x,i4,a,3es16.8)')iqibz,') ',Er%qibz(:,iqibz)
   call wrtout(unt,msg,mode)
  end do
 end if ! verbose>0

 rdwr=4
 !call hdr_io_int(Er%fform,Er%Hdr,rdwr,unt)

end subroutine print_epsilonm1_results 
!!***


!!****f* ABINIT/Epsm1_symmetrizer
!! NAME
!!  Epsm1_symmetrizer
!!
!! FUNCTION
!!  Symmetrize the inverse dielectric matrix, namely calculate epsilon^{-1} at a generic 
!!  q-point in the BZ starting from the knowledge of the matrix at a q-point in the IBZ.
!!  The procedure is quite generic and can be used for every two-point function which has 
!!  the same symmetry as the the crystal. 
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nomega=Number of frequencies required. All frequencies from 1 up to nomega are symmetrized.
!!  npw=Number of G vectors in symmetrized matrix, has to be smaller than Er%npwe.
!!  remove_exchange=If .TRUE., return e^{.1}-1 namely remove the exchange part.
!!  Er <Epsilonm1_results>= Data structure containing data on the inverse dielectric matrix.
!!  Gsph<Gvectors_type>=data related to the G-sphere
!!    %grottb
!!    %phmSGt 
!!  Qmesh<BZ_mesh_type>=Info on the q-mesh
!!    %nbz=number if q-points in the BZ
!!    %tab(nbz)=index of the symmeric q-point in the IBZ, for each point in the BZ
!!    %tabo(nbz)=the operation that rotates q_ibz onto \pm q_bz (depending on tabi) 
!!    %tabi(nbz)=-1 if time-reversal has to be considered, 1 otherwise
!!  iq_bz=Index of the q-point in the BZ where epsilon^-1 is required. 
!!
!! OUTPUT
!!  epsm1_qbz(npwc,npwc,nomega)=The inverse dielectric matrix at the q-point defined by iq_bz. 
!!   Exchange part can be subtracted out.
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the 
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere 
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required 
!!  to reconstruct the BZ.
!! 
!!  * Remember the symmetry properties of \tilde\espilon^{-1}
!!    If q_bz=Sq_ibz+G0:
!! 
!!    $\epsilon^{-1}_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau}\epsilon^{-1}_{G1,G2)}(q)
!!
!!    If time-reversal symmetry can be used then :
!!    $\epsilon^{-1}_{G1,G2}(-q_bz) = e^{+i(G1-G2).\tau}\epsilon^{-1}_{-S^{-1}(G1+Go),-S^{-1}(G2+G0)}^*(q)
!!
!! TODO
!!  Symmetrization can be skipped if iq_bz correspond to a point in the IBZ
!!
!! PARENTS
!!      csigme
!!
!! CHILDREN
!!
!! SOURCE

subroutine Epsm1_symmetrizer(iq_bz,nomega,npwc,Er,Gsph,Qmesh,remove_exchange,epsm1_qbz) 

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,nomega,npwc
 logical,intent(in) :: remove_exchange
 type(Epsilonm1_results),intent(in) :: Er
 type(Gvectors_type),intent(in) :: Gsph
 type(BZ_mesh_type),intent(in) :: Qmesh
!arrays
 complex(gwpc),intent(out) :: epsm1_qbz(npwc,npwc,nomega)  

!Local variables-------------------------------
!scalars
 integer :: iomega,ii,jj,iq_ibz,iiq,isymq,iq_loc
 character(len=500) :: msg
!arrays
 integer,pointer :: grottb(:)
 complex(gwpc),pointer :: phmSgt(:)

! *********************************************************************

#if defined DEBUG_MODE
 call assert( (Er%nomega>=nomega),'Too much frequencies required in Epsm1_symmetrizer')
 call assert( (Er%npwe  >=npwc),  'Too much G-vectors required in Epsm1_symmetrizer')
#endif

 ! FIXME here there is a problem with the small q, still cannot use BZ methods
 iq_ibz = Qmesh%tab (iq_bz) 
 isymq  = Qmesh%tabo(iq_bz) 
 iiq    = (3-Qmesh%tabi(iq_bz))/2

 grottb => Gsph%rottb (1:npwc,iiq,isymq)
 phmSgt => Gsph%phmSGt(1:npwc,isymq) 

 ! If out-of-memory, only Er%espm1(:,:,:,1) has been allocated and filled. 
 iq_loc=iq_ibz ; if (Er%mqmem==0) iq_loc=1 

 do iomega=1,nomega
  do ii=1,npwc
   do jj=1,npwc
    epsm1_qbz(grottb(ii),grottb(jj),iomega)=Er%epsm1(ii,jj,iomega,iq_loc)*phmSgt(ii)*CONJG(phmSgt(jj))
   end do
  end do
 end do
 !
 ! === Account for time-reversal ===
 if (iiq==2) epsm1_qbz(:,:,:)=CONJG(epsm1_qbz(:,:,:))

 ! === Subtract the exchange contribution ===
 if (remove_exchange) then
  do iomega=1,nomega
   do ii=1,npwc
    epsm1_qbz(ii,ii,iomega)=epsm1_qbz(ii,ii,iomega)-1.0_gwp
   end do
  end do
 end if

end subroutine Epsm1_symmetrizer
!!***

!!****f* ABINIT/init_Er_from_file
!! NAME
!!  init_Er_from_file
!!
!! FUNCTION
!!  Initialize basic dimensions and some important arrays in an Epsilonm1_results data type
!!  starting either from a _SCR or a _SUSC file.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_Er_from_file(Er,optfil,unt,fname,mqmem,accesswff,localrdwf,MPI_enreg)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => init_Er_from_file
!End of the abilint section

 implicit none

!Arguments ------------------------------------

 integer,intent(in) :: localrdwf,mqmem,optfil,unt,accesswff
 character(len=fnlen),intent(in) :: fname
 type(MPI_type),intent(in) :: MPI_enreg
 type(Epsilonm1_results),intent(inout) :: Er

!Local variables-------------------------------
!scalars
 integer :: iomega,fform,nqibzR,nqlwlR,nomegaR,npweR,npwwfn_used
 integer :: nbnds_used
 character(len=500) :: msg                   
!arrays
 integer,pointer :: gvec_p(:,:)
 real(dp),pointer :: qibz_p(:,:),qlwl_p(:,:)
 complex(dpc),pointer :: omega_p(:)
 character(len=80) :: title(2) 

! *********************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' init_Er_from_file : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 call nullify_epsilonm1_results(Er)

 !this has to be rewritten, for the moment it works although it is quite dirty.
 Er%ID=optfil ! 0 for e^-1 
              ! 1 for chi0

 call testscr(optfil,unt,fname,nqibzR,nqlwlR,nomegaR,npweR,npwwfn_used,nbnds_used,title,&
& fform,MPI_enreg,localrdwf,qibz_p,qlwl_p,omega_p,gvec_p,Er%Hdr)

 ! === Generic Info ===
 Er%fname      =fname
 Er%fform      =fform
 Er%npwe_file  =npweR
 Er%nomega_file=nomegaR

 Er%npwwfn_used=npwwfn_used
 Er%nbnds_used =nbnds_used

!FIXME Hardcoded values but they should be reported by testscr
!BEGIN HARCODED
 Er%nI=1
 Er%nJ=1
 Er%ikxc=0
 Er%Tordering=0     ! TO
 Er%test_type=-1    ! No probing charge, only for TDDFT
 Er%inclvkb=0
!END HARCODED

 Er%npwe=npweR  

 Er%nqibz=nqibzR
 Er%mqmem=mqmem 
 if (mqmem/=0) Er%mqmem=Er%nqibz
 allocate(Er%qibz(3,Er%nqibz))
 Er%qibz(:,:)=qibz_p(:,:)
 deallocate(qibz_p)

 ! this has to be done in a cleaner way.
 Er%nqlwl=nqlwlR
 if (associated(qlwl_p)) then 
  allocate(Er%qlwl(3,nqlwlR))
  Er%qlwl(:,:)=qlwl_p(:,:)
  deallocate(qlwl_p)
 end if

 Er%nomega=nomegaR
 allocate(Er%omega(Er%nomega))
 Er%omega(:)=omega_p(:)
 deallocate(omega_p)

 if (Er%nomega==2) then
  Er%nomega_r=1 
  Er%nomega_i=1
 else
  ! Real frequencies are packed first.
  Er%nomega_r=1
  do iomega=1,Er%nomega
   if (REAL(Er%omega(iomega))>0.001*Ha_eV) Er%nomega_r=iomega
  end do
  Er%nomega_i=Er%nomega-Er%nomega_r
 end if     

 ! === Save G-vectors ===
 allocate(Er%gvec(3,Er%npwe_file))
 Er%gvec(:,:)=gvec_p(:,:)
 deallocate(gvec_p)

#if defined DEBUG_MODE
 write(msg,'(a)')' init_Er_from_file : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine init_Er_from_file
!!***


!!****f* ABINIT/mkdump_Er
!! NAME
!!  mkdump_Er
!!
!! FUNCTION
!!  Dump the content of an Epsilonm1_results data type on file.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkdump_Er(Er,Dtset,Dtfil,unt_dump,fname_dump,accesswff,localrdwf,gmet,MPI_enreg)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => mkdump_Er
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: localrdwf,accesswff,unt_dump
 type(Dataset_type),intent(in) :: Dtset
 type(Datafiles_type),intent(in) :: Dtfil
 type(MPI_type),intent(in) :: MPI_enreg
 type(Epsilonm1_results),intent(inout) :: Er
 character(len=fnlen),intent(in) :: fname_dump
!arrays 
 real(dp),intent(in) :: gmet(3,3)

!Local variables-------------------------------
!scalars
 integer :: iqibz,optfil,mqmem
 character(len=80) :: title(2) 
 character(len=500) :: msg                   
!arrays
 complex(gwpc),pointer :: epsm1(:,:,:)

! *********************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' mkdump_Er : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 ! According to options either perform a direct dumping or 
 ! calculate a new e^-1 and dump the result on a file for 
 ! a subsequent use.

 ! TODO, write function to return title, just for info
 title(1)='SCR file: epsilon^-1'
 title(2)='TESTPARTICLE' 
 optfil=Er%ID
 ! TODO remove Dtset as code-reusing and encapsulation is hard!
 do iqibz=1,Er%nqibz
  epsm1 => Er%epsm1(:,:,:,iqibz)
  call wrscr(iqibz,optfil,unt_dump,fname_dump,Er%Hdr,Dtset,Er%npwe,Er%npwwfn_used,Er%nbnds_used,&
&  Er%nqibz,Er%nqlwl,Er%nomega,Er%qibz,Er%omega,Er%gvec,gmet,epsm1,title,&
&  Er%qlwl,Er%uwing,Er%lwing) 
 end do

 ! Now Er% "belongs" to file fname_dump, thus we destroy and re-initialize the object.
 ! TODO should write a method to copy Er% structures
 ! TODO recheck this part
 !call destroy_Epsilonm1_results(Er)
 !fname=Dtfil%filscr 
 !unt_dump =Dtfil%unscr
 !if (optfil==1) then 
 ! fname=Dtfil%filchi0 
 ! unt_dump =Dtfil%unchi0
 !end if
 !mqmem=Er%mqmem
 !call init_Er_from_file(Er,optfil,unt_dump,fname,mqmem,accesswff,localrdwf,MPI_enreg)

#if defined DEBUG_MODE
 write(msg,'(a)')' mkdump_Er : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine mkdump_Er
!!***


!!****f* ABINIT/get_epsm1
!! NAME
!!  get_epsm1
!!
!! FUNCTION
!!  Working in progress but the main is idea is as follows:
!!
!!  Return the symmetrized inverse dielectric matrix.
!!  This method implements both in-core and the out-of-core solution 
!!  In the later, epsilon^-1 or chi0 are read from file.
!!  It is possible to specify options to retrieve (RPA |TDDDT, [TESTCHARGE|TESTPARTICLE]).
!!  All dimensions are already initialized in the Er% object, this method 
!!  should act as a wrapper around rdscr and make_epsm1_driver. A better 
!!  implementation will be done in the following once the coding of file handlers is completed.
!!
!! INPUTS
!!  Dtfil<Datafiles_type)>=datatype containing filenames
!!  Vcp<Coulombian_type>=Structure gathering data on the Coulombian interaction
!!  iqibzA[optional]=Index of the q-point to be read from file (only for out-of-memory solutions)
!!  accesswff=option definig the file format.
!!  localrdwf=1 if file is local to each machine
!!  option_test
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  kxcg=XC kernel in reciprocal space.
!!  MPI_enreg=informations about MPI parallelization
!!
!! OUTPUT
!!  Er%epsm1
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE


subroutine get_epsm1(Er,Vcp,Dtfil,approx_type,option_test,accesswff,localrdwf,kxcg,gmet,MPI_enreg,iqibzA)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => get_epsm1
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: localrdwf,accesswff,option_test,approx_type
 integer,optional,intent(in) :: iqibzA
 type(Coulombian_type),intent(in) :: Vcp
 type(Datafiles_type),intent(in) :: Dtfil
 type(MPI_type),intent(in) :: MPI_enreg
 type(Epsilonm1_results),intent(inout) :: Er
!arrays
 real(dp),intent(in) :: gmet(3,3)
 !complex(gwpc),intent(in) :: kxcg(npwe*approx_type,npwe*approx_type) 
 !FIXME this has to be rewritten LDA is (npwe,1), should pass first and second_dimension
 complex(gwpc),intent(in) :: kxcg(:,:) 

!Local variables-------------------------------
!scalars
 integer :: optfil,istat,iqibz,unt
 character(len=fnlen) :: fname
 character(len=500) :: msg                   
!arrays
 real(dp) :: qibz_dum(3,Er%nqibz)
 complex(gwpc) :: omega_dum(Er%nomega)
 complex(gwpc),pointer :: epsm1(:,:,:)

! *********************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' get_epsm1 : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 ! === File identifier ===
 ! --> 0 to read epsilon^-1
 ! --> 1 to read chi0
 !optfil=Er%IDfile
 optfil=Er%ID
 unt  =Dtfil%unscr 
 fname=Dtfil%filscr
 if (optfil==1) then 
  unt  =Dtfil%unchi0 
  fname=Dtfil%filchi0
 end if

 select case (Er%mqmem)
  case (0)
   ! === Out-of-core solution ===
   ! TODO clean a bit the treatment of wings and qlwl
   if (associated(Er%qlwl )) deallocate(Er%qlwl )
   if (associated(Er%lwing)) deallocate(Er%lwing)
   if (associated(Er%uwing)) deallocate(Er%uwing)
   if (associated(Er%epsm1)) deallocate(Er%epsm1)
   allocate(Er%qlwl(3,Er%nqlwl*optfil))
   allocate(Er%lwing(Er%npwe,Er%nomega,Er%nqlwl*optfil))
   allocate(Er%uwing(Er%npwe,Er%nomega,Er%nqlwl*optfil))
   allocate(Er%epsm1(Er%npwe,Er%npwe,Er%nomega,1),STAT=istat)
   !if (istat/=0) call memerr(FILE__,'Er%epsm1',Er%npwe**2*Er%nomega*Er%nibz,'gwpc')

   !Here I should read a new file containg a previsously calculated e^-1
   ! but I have to pass the name of the input file!
   !TODO clean rdscr, remove all the inout quantities except wings and epsm1
   ! this an headache due to mrgscr
   call rdscr(optfil,unt,fname,Er%npwe,Er%nqibz,1,Er%nomega,qibz_dum,omega_dum,gmet,&
&   Er%epsm1,MPI_enreg,localrdwf,Er%nqlwl,.FALSE.,Er%qlwl,Er%uwing,Er%lwing,iqiA=iqibzA)

   if (Er%ID==0) then 
     ! === If q-slice of epsilon^-1 has been read then return === 
     call print_epsilonm1_results(Er,unit=std_out)
     RETURN 
   else if (Er%ID==1) then 
    ! === if q-slice if chi0 has been read then make e^-1 ===
    ! FIXME remember the wings and fix issue with master inside make_epsm1_driver!!!!!
    ! FIXME treatment of kxcg is messy!!!!
    STOP "Solve problem with Er%ID"
    epsm1 => Er%epsm1(:,:,:,1)
    call make_epsm1_driver(iqibzA,Er%npwe,Er%nI,Er%nJ,Er%nomega,Er%nomega_r,Er%omega,&
&    approx_type,option_test,Er%nqibz,Er%qibz,Vcp,Dtfil,gmet,kxcg,MPI_enreg,epsm1)
    !Er%ID=0 be careful here, I cannot change ID!!!!!!
   else 
    STOP 'Wrong Er%ID'
   end if

  case default
   ! ========================
   ! === In-core solution ===
   ! ========================
   allocate(Er%qlwl(3,Er%nqlwl*optfil))
   allocate(Er%lwing(Er%npwe,Er%nomega,Er%nqlwl*optfil))
   allocate(Er%uwing(Er%npwe,Er%nomega,Er%nqlwl*optfil))
   allocate(Er%epsm1(Er%npwe,Er%npwe,Er%nomega,Er%nqibz),STAT=istat)
   !if (istat/=0) call memerr(FILE__,'Er%epsm1',Er%npwe**2*Er%nomega*Er%nibz,'gwpc')

   !TODO clean rdscr, remove all the INOUT quantities except wings and epsm1
   call rdscr(optfil,unt,fname,Er%npwe,Er%nqibz,Er%nqibz,Er%nomega,qibz_dum,omega_dum,gmet,&
&   Er%epsm1,MPI_enreg,localrdwf,Er%nqlwl,.TRUE.,Er%qlwl,Er%uwing,Er%lwing)

   if (Er%ID==0) then 
     ! === If epsilon^-1 has been read then return === 
     !call print_epsilonm1_results(Er,unit=std_out)
     RETURN 
   else if (Er%ID==1) then 
    ! === if chi0 has been read then make e^-1 ===
    ! FIXME remember the wings and fix issue with master inside make_epsm1_driver!!!!!
    ! FIXME treatment of kxcg is messy!!!!
    do iqibz=1,Er%nqibz
     epsm1 => Er%epsm1(:,:,:,iqibz)
     call make_epsm1_driver(iqibz,Er%npwe,Er%nI,Er%nJ,Er%nomega,Er%nomega_r,Er%omega,&
&     approx_type,option_test,Er%nqibz,Er%qibz,Vcp,Dtfil,gmet,kxcg,MPI_enreg,epsm1)
    end do
    Er%ID=0 !this is correct, figure how to solve the problem with out-of-core.
   else 
    STOP 'Wrong Er%ID'
   end if

 end select

#if defined DEBUG_MODE
 write(msg,'(a)')' get_epsm1 : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine get_epsm1
!!***


