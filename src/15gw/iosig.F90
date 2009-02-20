!{\src2tex{textfont=tt}}
!!****f* ABINIT/write_sigma_header
!! NAME
!! write_sigma_header
!!
!! FUNCTION
!! input/output sigma results
!! This file contains 4 routines : 
!!  write_sigma_results_header,
!!  write_sigma_results, rdgw,
!!  rdgw
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  Sp=sigma_parameters
!!  Cryst<Crystal_structure>= Info on the Crystal structure
!!  Kmesh<Bz_mesh_type>= Description of the BZ sampling.
!!
!! OUTPUT
!!  (for writing routines, no output) otherwise, should be described
!!
!! NOTES
!!
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      cgemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine write_sigma_results_header(Sp,Er,Cryst,Kmesh,Qmesh)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Bz_mesh_type),intent(in) :: Kmesh,Qmesh
 type(Crystal_structure),intent(in) :: Cryst
 type(Epsilonm1_results),intent(in) :: Er
 type(Sigma_parameters),intent(in) :: Sp

!Local variables-------------------------------
!scalars
 integer :: mod10
 character(len=500) :: msg

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' write_sigma_results_header : enter '
  call wrtout(std_out,msg,'COLL')
  call flush_unit(std_out)
#endif

 write(msg,'(a)')' SIGMA fundamental parameters:'
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

 mod10=MOD(Sp%gwcalctyp,10)
 SELECT CASE (mod10)
 CASE (0)
  write(msg,'(a,i2)')' PLASMON POLE MODEL ',Sp%ppmodel
 CASE (1)
  write(msg,'(a)')' ANALYTIC CONTINUATION'
 CASE (2)
  write(msg,'(a)')' CONTOUR DEFORMATION'
 CASE (5)
  write(msg,'(a)')' Hartree-Fock'
 CASE (6)
  write(msg,'(a)')' Screened Exchange'
 CASE (7)
  write(msg,'(a)')' COHSEX'
 CASE (8)
  write(msg,'(a,i2)')' MODEL GW with PLASMON POLE MODEL ',Sp%ppmodel
 CASE (9)
  write(msg,'(a)')' MODEL GW without PLASMON POLE MODEL'
 CASE DEFAULT
  write(msg,'(4a,i3)')ch10,&
&  ' write_sigma_results_header : BUG ',ch10,&
&  ' wrong value for Sp%gwcalctyp = ',Sp%gwcalctyp 
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 END SELECT
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

 write(msg,'(a,i12)')' number of plane-waves for SigmaX         ',Sp%npwx
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of plane-waves for SigmaC and W   ',Sp%npwc
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of plane-waves for wavefunctions  ',Sp%npwwfn
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of bands                          ',Sp%nbnds
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of independent spin polarizations ',Sp%nsppol
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 !write(msg,'(a,i12)')' number of spinorial components           ',Sp%nspinor
 !call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of k-points in IBZ                ',Kmesh%nibz
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of q-points in IBZ                ',Qmesh%nibz
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of symmetry operations            ',Cryst%nsym
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of k-points in BZ                 ',Kmesh%nbz
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of q-points in BZ                 ',Qmesh%nbz
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of omega for Sigma on real axis   ',Sp%nomegasrd
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,f12.2)')' deltae [eV]                              ',Sp%deltae*Ha_eV
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of omega for Sigma on real axis   ',Sp%nomegasr
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,f12.2)')' max omega for Sigma on real axis  [eV]   ',Sp%omegasrmax*Ha_eV
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

 if (Sp%soenergy>0.1d-4) then 
  write(msg,'(a,f12.2)')' scissor energy [eV]                      ',Sp%soenergy*Ha_eV
  call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 end if

 if (MOD(Sp%gwcalctyp,10)==1) then
  write(msg,'(a,i12)')' number of omega for Sigma on imag axis   ',Sp%nomegasi
  call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,f12.2)')' max omega for Sigma on imag axis  [eV]   ',Sp%omegasimax*Ha_eV
  call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 end if 

 write(msg,'(2a)')ch10,' EPSILON^-1 parameters (SCR file):'
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 !write(std_out,*) titem1(2)(1:79)
 write(msg,'(a,i12)')' dimension of the eps^-1 matrix           ',Er%npwe_file
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 !£write(msg,'(a,i12)')' dimension of the eps^-1 matrix on file   ',Er%npwe_file
 !£call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 !£write(msg,'(a,i12)')' dimension of the eps^-1 matrix used      ',Er%npwe
 !£call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of plane-waves for wavefunctions  ',Er%npwwfn_used
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of bands                          ',Er%nbnds_used
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of q-points in IBZ                ',Qmesh%nibz
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of frequencies                    ',Er%nomega
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of real frequencies               ',Er%nomega_r
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of imag frequencies               ',Er%nomega_i
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

 write(msg,'(3a)')ch10,' matrix elements of self-energy operator (all in [eV])',ch10
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

 if (Sp%gwcalctyp<10) then
  write(msg,'(a)')' Perturbative Calculation'
 else if (Sp%gwcalctyp<20) then
  write(msg,'(a)')' Self-Consistent on Energies only'
 else
  write(msg,'(a)')' Self-Consistent on Energies and Wavefunctions'
 end if
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

#if defined DEBUG_MODE
 write(msg,'(a)')' write_sigma_results_header : exit '
 call wrtout(std_out,msg,'COLL')
 call flush_unit(std_out)
#endif

end subroutine write_sigma_results_header
!!***


!!****f* ABINIT/write_sigma_results
!! NAME
!! write_sigma_results
!!
!! FUNCTION
!! write the final results of the GW calculation
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2008 ABINIT group (FBruneval, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! en_lda(nkibz,sp%nbnds,sp%nsppol)= KS energies
!! ikibz= index of the k-point in the array kibz, where GW corrections are calculated 
!! ikcalc= index of the k-point in the array sp%kcalc 
!! Kmesh<Bz_mesh_type>=Info on the k-point mesh
!! sp=sigma_parameters datatype
!! sr=sigma results datatype
!!
!! OUTPUT
!!  (for writing routines, no output) otherwise, should be described
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      cgemm
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine write_sigma_results(Sp,Sr,ikcalc,ikibz,Kmesh,usepawu,en_lda)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : unt_gw, unt_sig, unt_sgr, unt_sgm


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => write_sigma_results
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikcalc,ikibz,usepawu
 type(Bz_mesh_type),intent(in) :: Kmesh
 type(Sigma_parameters),intent(in) :: Sp
 type(Sigma_results),intent(in) :: Sr
!arrays
 real(dp),intent(in) :: en_lda(Kmesh%nibz,Sp%nbnds,Sp%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ib,io,is
 character(len=500) :: msg
!arrays
 character(len=12) :: tag_spin(2)

! *************************************************************************

 !unt_gw  = 21  ! File with GW corrections.
 !unt_sig = 22  ! Self-energy as a function of frequency.
 !unt_sgr = 23  ! Derivative wrt omega of the Self-energy.
 !unt_sgm = 20  ! Sigma on the Matsubara axis.

 tag_spin=(/'            ','            '/)
 if (Sp%nsppol==2) tag_spin=(/',  SPIN UP  ',',  SPIN DOWN'/)

 do is=1,Sp%nsppol
  write(msg,'(2a,3f8.3,a)')ch10,' k = ',Sp%xkcalc(:,ikcalc),tag_spin(is)
  call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

  write(msg,'(a)')&
&  ' Band     E0 <VxcLDA>   SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E'

  if (usepawu/=0) then
   write(msg,'(a)')&
&   ' Band     E0 <VxcLDA>   <H_U>  SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E'
  end if

  if (Sp%gwcalctyp>=10) then
   write(msg,'(2a)')&
&   ' Band     E_lda   <Vxclda>   E(N-1)  <Hhartree>   SigX  SigC[E(N-1)]',&
&   '    Z     dSigC/dE  Sig[E(N)]  DeltaE  E(N)_pert E(N)_diago'
  end if
  call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

  write(unt_gw,'(3f10.6)')Sp%xkcalc(:,ikcalc)
  write(unt_gw,'(i4)')Sp%maxbnd(ikcalc)-Sp%minbnd(ikcalc)+1

  write(unt_sig,'("# k = ",3f10.6)')Sp%xkcalc(:,ikcalc)
  write(unt_sig,'("# b = ",2i10)')Sp%minbnd(ikcalc),Sp%maxbnd(ikcalc)

  write(unt_sgr,'("# k = ",3f10.6)')Sp%xkcalc(:,ikcalc)
  write(unt_sgr,'("# b = ",2i10)')Sp%minbnd(ikcalc),Sp%maxbnd(ikcalc)

  do ib=Sp%minbnd(ikcalc),Sp%maxbnd(ikcalc)

   if (Sp%gwcalctyp>=10) then
    call print_Sigma_SC(Sr,ikibz,ib,is,Sp,Kmesh,en_lda,usepawu,unit=ab_out)
    call print_Sigma_SC(Sr,ikibz,ib,is,Sp,Kmesh,en_lda,usepawu,unit=std_out,prtvol=1)
   else
    call print_Sigma_perturbative(Sp,Sr,ikibz,ib,is,usepawu,unit=ab_out)
    call print_Sigma_perturbative(Sp,Sr,ikibz,ib,is,usepawu,unit=std_out,prtvol=1)
   end if

   write(unt_gw,'(i6,3f9.4)')         &
&   ib,                               &
&   REAL (Sr%egw(ib,ikibz,is)) *Ha_eV,&
&   REAL (Sr%degw(ib,ikibz,is))*Ha_eV,&
&   AIMAG(Sr%egw(ib,ikibz,is)) *Ha_eV
   
  end do !ib

  if (Sr%e0gap(ikibz,is)**2+Sr%egwgap(ikibz,is)**2+Sr%degwgap(ikibz,is)**2 > tol10) then
   ! Output the direct gap for each spin
   ! If all the gaps are zero, this means that it could not be computed in the calling routine
   write(msg,'(2a,f8.3)')ch10,' E^0_gap       ',Sr%e0gap(ikibz,is)*Ha_eV
   call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,f8.3)')      ' E^GW_gap      ',Sr%egwgap(ikibz,is)*Ha_eV
   call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,f8.3,a)')    ' DeltaE^GW_gap ',Sr%degwgap(ikibz,is)*Ha_eV,ch10
   call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
  end if

  ! === Output spectral function ===
  do io=1,Sr%nomega
   write(unt_sig,'(100(e11.5,2x))')&
&   REAL(Sp%omegasf(io))*Ha_eV,&
&   (REAL(Sr%sigxcme(ib,ikibz,io,is))*Ha_eV,&
&   AIMAG(Sr%sigxcme(ib,ikibz,io,is))*Ha_eV,&
&   one/pi*ABS(AIMAG(Sr%sigcme(ib,ikibz,io,is)))&
&   /( (REAL(Sp%omegasf(io)-Sr%hhartree(ib,ib,ikibz,is)-Sr%sigxcme(ib,ikibz,io,is)))**2&
&     +(AIMAG(Sr%sigcme(ib,ikibz,io,is)))**2) /Ha_eV,&
&   ib=Sp%minbnd(ikcalc),Sp%maxbnd(ikcalc))
  end do

  do ib=Sp%minbnd(ikcalc),Sp%maxbnd(ikcalc)
   write(unt_sgr,'("# ik, ib",2i5)')ikibz,ib
   do io=1,Sr%nomegasrd
    write(unt_sgr,'(100(e11.5,2x))')             &
&    REAL (Sr%omegasrd(ib,ikibz,io,is))   *Ha_eV,&
&    REAL (Sr%sigxcmesrd(ib,ikibz,io,is)) *Ha_eV,&
&    AIMAG(Sr%sigxcmesrd(ib,ikibz,io,is)) *Ha_eV
   end do
  end do

  if (MOD(sp%gwcalctyp,10)==1) then 
   ! For AC write matrix elements of sigma along the imaginary axis
   do ib=Sp%minbnd(ikcalc),Sp%maxbnd(ikcalc)
    write(unt_sgm,'("# ik, ib",2i5)')ikibz,ib
    do io=1,Sp%nomegasi
     write(unt_sgm,'(3(e11.5,2x))')             &
&     AIMAG(Sp%omegasi(io))              *Ha_eV,&
&     REAL (Sr%sigxcmesi(ib,ikibz,io,is))*Ha_eV,&
&     AIMAG(Sr%sigxcmesi(ib,ikibz,io,is))*Ha_eV
    end do
   end do
  end if 

 end do !is

end subroutine write_sigma_results
!!***


!!****f* ABINIT/rdgw
!! NAME
!! rdgw
!!
!! FUNCTION
!!  This subroutine reads the GW corrections from a _GW file
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2008 ABINIT group (FBruneval, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nkibz=number of k-points
!!  nsppol=number of spin
!!  nbnds=number of bands in the present GW calculation
!!  kibz(3,nkibz)= irreducible k-points
!!  nbv(nsppol)= index of the valence band for each spin
!!
!! OUTPUT
!!  gwenergy(nkibz,nbnds,nsppol) : QP energies as read or derived from the data contained
!!   in the external file
!! 
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!      cgemm
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine rdgw(nkibz,nbnds,nbv,nsppol,kibz,gwenergy)

 use defs_basis
 use defs_datatypes
 use m_numeric_tools, only : linfit
 use m_io_tools, only : get_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbnds,nkibz,nsppol
!arrays
 integer,intent(in) :: nbv(nsppol)
 real(dp),intent(in) :: kibz(3,nkibz)
 real(dp),intent(out) :: gwenergy(nkibz,nbnds,nsppol)

!Local variables ------------------------------
!scalars
 integer :: ib,ibr,ik,ikibz,ikr,ios,is,n,nbR,nkR,nsR,unt
 real(dp) :: a,b,degw,egw,smrt
 character(len=500) :: msg
 character(len=fnlen) :: filnam
!arrays
 real(dp) :: gwcorr(nkibz,nbnds,nsppol),k(3)

!************************************************************************

 filnam='in.gw' !TODO should become input
 write(msg,'(2a)')' reading GW corrections from file ',trim(filnam)
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 
 unt=get_unit()
 open(unt,file=TRIM(filnam),status='old',iostat=ios)
 if (ios/=0) then 
  write(msg,'(6a)')ch10,' rdgw: ERROR- ',ch10,  &
&  ' opening file: ',TRIM(filnam),' as old'
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 read(unt,*) nkR,nsR 
 if (nkR/=nkibz) then  
  write(msg,'(4a,i4,a,i4,3a)')ch10,&
&  ' rdgw : WARNING - ,',ch10,&
&  ' Found less k-points than that required ',nkR,'/',nkibz,&
&  ch10,' Some of the k-points will be skipped ',ch10
  call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 end if 
 if (nsR/=nsppol) then  
  write(msg,'(4a)')ch10,&
&  ' rdgw : ERROR - ,',ch10,&
&  ' Found different value of nsspol'
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if 

 gwcorr(:,:,:)=0.0
 do is=1,nsppol
  do ikr=1,nkR
   read(unt,*) k
   read(unt,*) nbR
   ikibz=0
   do ik=1,nkibz
    if (all(abs(k(:)-kibz(:,ik))<0.0001)) ikibz=ik
   end do
   do ib=1,nbR
    read(unt,*) ibr,egw,degw
    if(ibr<=nbnds .and. ikibz/=0) gwcorr(ikibz,ibr,is)=degw/Ha_eV
   end do
  end do 
 end do
 close(unt)
 !
 ! nbv(is) is used to take into account the valence band index for each spin
 do is=1,nsppol
  do ik=1,nkibz
   
   n=nbnds-nbv(is)
   do ib=nbv(is)+1,nbnds
    if (gwcorr(ik,ib,is)==0) then
     n=ib-1-nbv(is)
     if (n>1) then
      write(msg,'(a)')&
&      ' linear extrapolating (conduction) GW corrections beyond the read values'
      call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
      smrt=linfit(n,gwenergy(ik,nbv(is)+1:nbv(is)+n,is),gwcorr(ik,nbv(is)+1:nbv(is)+n,is),a,b)
     else
      write(msg,'(a)')' assuming constant (conduction) GW corrections beyond the read values'
      call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
      a=0 ; b=gwcorr(ik,nbv(is)+n,is)
     end if
     exit !ib
    end if
   end do !ib

   do ib=nbv(is)+n+1,nbnds
    gwcorr(ik,ib,is)= a*gwenergy(ik,ib,is) + b
   end do

   n=nbv(is)
   do ib=nbv(is),1,-1
    if (gwcorr(ik,ib,is)==0) then
     n= nbv(is)-ib
     if (n>1) then
      write(msg,'(a)')' linear extrapolating (valence) GW corrections beyond the read values'
      call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
      smrt=linfit(n,gwenergy(ik,nbv(is)-n+1:nbv(is),is),gwcorr(ik,nbv(is)-n+1:nbv(is),is),a,b)
     else
      write(msg,'(a)')' assuming constant (valence) GW corrections beyond the read values'
      call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
      a=0 ; b=gwcorr(ik,nbv(is),is)
     end if
     exit !ib
    end if
   end do !ib

   do ib=1,nbv(is)-n
    gwcorr(ik,ib,is)=a*gwenergy(ik,ib,is) + b
   end do

  end do !ik
 end do !is

 write(msg,'(a)')' k  s     GW corrections [eV] '
 call wrtout(std_out,msg,'COLL')
 do is=1,nsppol
  do ik=1,nkibz
   write(*,'(i3,1x,i3,10f7.2/50(10x,10f7.2/))') ik,is,(Ha_eV*gwcorr(ik,ib,is),ib=1,nbnds)
  end do
 end do 
 gwenergy(:,:,:)=gwenergy(:,:,:)+gwcorr(:,:,:)

 write(msg,'(a)')' k   s    GW eigenvalues [eV]'
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 do is=1,nsppol
  do ik=1,nkibz
   write(*,     '(i3,7x,10f7.2/50(10x,10f7.2/))') ik,is,(Ha_eV*gwenergy(ik,ib,is),ib=1,nbnds)
   write(ab_out,'(i3,7x,10f7.2/50(10x,10f7.2/))') ik,is,(Ha_eV*gwenergy(ik,ib,is),ib=1,nbnds)
  end do
 end do
 write(std_out,*) ; write(ab_out,*)
end subroutine rdgw
!!***


!!****f* ABINIT/print_Sigma_perturbative 
!! NAME
!! print_Sigma_perturbative
!!
!! FUNCTION
!!  write the results of the GW calculation done with the perturbative approach
!!

!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine print_Sigma_perturbative(Sp,SigR,ik_ibz,iband,isp,usepawu,unit,prtvol,mode_paral)

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband,ik_ibz,isp,usepawu
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 type(Sigma_parameters),intent(in) :: Sp
 type(Sigma_results),intent(in) :: SigR

!Local variables-------------------------------
!scalars
 integer :: unt,verbose
 logical :: ltest
 character(len=4) :: mode
 character(len=500) :: msg

! *********************************************************************

 unt=std_out ; if (PRESENT(unit))       unt=unit
 verbose=0   ; if (PRESENT(prtvol))     verbose=prtvol
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 !write(msg,'(a)')' Band     E0 <VxcLDA>   SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E '
 !call wrtout(unt,msg,mode) 

 if (usepawu==0) then 
  if (Sp%nsig_ab/=1) then
   write(msg,'(i5,9f8.3)')                         & 
&        iband,                                    &
&        SigR%e0          (iband,ik_ibz,1)*Ha_eV,  &
&        SUM(SigR%vxcme   (iband,ik_ibz,:))*Ha_eV, &
&        SUM(SigR%sigxme  (iband,ik_ibz,:))*Ha_eV, &
&   REAL(SUM(SigR%sigcmee0(iband,ik_ibz,:)))*Ha_eV,&
&   REAL(SigR%ze0         (iband,ik_ibz,1)),       &
&   REAL(SUM(SigR%dsigmee0(iband,ik_ibz,:))),      &
&   REAL(SUM(SigR%sigmee  (iband,ik_ibz,:)))*Ha_eV,&
&   REAL(SigR%degw        (iband,ik_ibz,1))*Ha_eV, &
&   REAL(SigR%egw         (iband,ik_ibz,1))*Ha_eV
    call wrtout(unt,msg,mode) 

   if (verbose/=0) then
    write(msg,'(i5,9f8.3)')                          & 
&          iband,                                    &
&          zero,                                     &
&          zero,                                     &
&          zero,                                     &
&    AIMAG(SUM(SigR%sigcmee0(iband,ik_ibz,:)))*Ha_eV,&
&    AIMAG(SigR%ze0         (iband,ik_ibz,1)),       &
&    AIMAG(SUM(SigR%dsigmee0(iband,ik_ibz,:))),      &
&    AIMAG(SUM(SigR%sigmee  (iband,ik_ibz,:)))*Ha_eV,&
&    AIMAG(SigR%degw        (iband,ik_ibz,1))*Ha_eV, &
&    AIMAG(SigR%egw         (iband,ik_ibz,1))*Ha_eV
     call wrtout(unt,msg,mode) 
   end if
  else
   write(msg,'(i5,9f8.3)')                      & 
&        iband,                                 &
&        SigR%e0      (iband,ik_ibz,isp)*Ha_eV, &
&        SigR%vxcme   (iband,ik_ibz,isp)*Ha_eV, &
&        SigR%sigxme  (iband,ik_ibz,isp)*Ha_eV, &
&   REAL(SigR%sigcmee0(iband,ik_ibz,isp))*Ha_eV,&
&   REAL(SigR%ze0     (iband,ik_ibz,isp)),      &
&   REAL(SigR%dsigmee0(iband,ik_ibz,isp)),      &
&   REAL(SigR%sigmee  (iband,ik_ibz,isp))*Ha_eV,&
&   REAL(SigR%degw    (iband,ik_ibz,isp))*Ha_eV,&
&   REAL(SigR%egw     (iband,ik_ibz,isp))*Ha_eV
    call wrtout(unt,msg,mode) 

   if (verbose/=0) then
    write(msg,'(i5,9f8.3)')                       & 
&          iband,                                 &
&          zero,                                  &
&          zero,                                  &
&          zero,                                  &
&    AIMAG(SigR%sigcmee0(iband,ik_ibz,isp))*Ha_eV,&
&    AIMAG(SigR%ze0     (iband,ik_ibz,isp)),      &
&    AIMAG(SigR%dsigmee0(iband,ik_ibz,isp)),      &
&    AIMAG(SigR%sigmee  (iband,ik_ibz,isp))*Ha_eV,&
&    AIMAG(SigR%degw    (iband,ik_ibz,isp))*Ha_eV,&
&    AIMAG(SigR%egw     (iband,ik_ibz,isp))*Ha_eV
     call wrtout(unt,msg,mode) 
   end if
  end if
 else 
  ! === PAW+U+GW calculation ===
  ltest=(Sp%nsig_ab==1)
  call assert(ltest,'LDA+U with spinor not implemented',&
&  __FILE__,__LINE__)
  write(msg,'(i5,10f8.3)')                     & 
&       iband,                                 &
&       SigR%e0      (iband,ik_ibz,isp)*Ha_eV, &
&       SigR%vxcme   (iband,ik_ibz,isp)*Ha_eV, &
&       SigR%vUme    (iband,ik_ibz,isp)*Ha_eV, &
&       SigR%sigxme  (iband,ik_ibz,isp)*Ha_eV, &
&  REAL(SigR%sigcmee0(iband,ik_ibz,isp))*Ha_eV,&
&  REAL(SigR%ze0     (iband,ik_ibz,isp)),      &
&  REAL(SigR%dsigmee0(iband,ik_ibz,isp)),      &
&  REAL(SigR%sigmee  (iband,ik_ibz,isp))*Ha_eV,&
&  REAL(SigR%degw    (iband,ik_ibz,isp))*Ha_eV,&
&  REAL(SigR%egw     (iband,ik_ibz,isp))*Ha_eV
   call wrtout(unt,msg,mode) 

  if (verbose/=0) then
   write(msg,'(i5,10f8.3)')                      & 
&         iband,                                 &
&         zero,                                  &
&         zero,                                  &
&         zero,                                  &
&         zero,                                  &
&   AIMAG(SigR%sigcmee0(iband,ik_ibz,isp))*Ha_eV,&
&   AIMAG(SigR%ze0     (iband,ik_ibz,isp)),      &
&   AIMAG(SigR%dsigmee0(iband,ik_ibz,isp)),      &
&   AIMAG(SigR%sigmee  (iband,ik_ibz,isp))*Ha_eV,&
&   AIMAG(SigR%degw    (iband,ik_ibz,isp))*Ha_eV,&
&   AIMAG(SigR%egw     (iband,ik_ibz,isp))*Ha_eV
    call wrtout(unt,msg,mode) 
  end if
 end if

end subroutine print_Sigma_perturbative 
!!***


!!****f* ABINIT/write_sigma_SC
!! NAME
!! print_Sigma_SC
!!
!! FUNCTION
!! write the results of the GW calculation in case of self-consistency
!!

!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine print_Sigma_SC(SigR,ik_ibz,iband,isp,SigP,Kmesh,en_lda,usepawu,unit,prtvol,mode_paral)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!this should be in the datatype
!scalars
 integer,intent(in) :: iband,ik_ibz,isp,usepawu
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Sigma_parameters),intent(in) :: SigP
 type(Sigma_results),intent(in) :: SigR
!arrays
 real(dp),intent(in) :: en_lda(Kmesh%nibz,SigP%nbnds,SigP%nsppol)

!Local variables-------------------------------
!scalars
 integer :: unt,verbose
 character(len=4) :: mode
 character(len=500) :: msg

! *********************************************************************

 unt=std_out ; if (PRESENT(unit))       unt=unit
 verbose=0   ; if (PRESENT(prtvol))     verbose=prtvol
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

! write(msg,'(a)')&
!& ' Band     E0 <VxcLDA>   SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E'
! call wrtout(unt,msg,mode) 

 if (usepawu==0) then
  if (SigP%nsig_ab/=1) then
   write(msg,'(i5,12(2x,f8.3))')                         & 
&        iband,                                          &
&        en_lda           (ik_ibz,iband,1)*Ha_eV,        &
         !sr%e0           (iband,ikibz,1)*Ha_eV,         &
&        SUM(SigR%vxcme   (iband,ik_ibz,:))*Ha_eV,       &
&        SigR%e0          (iband,ik_ibz,1)*Ha_eV,        &
&   REAL(SUM(SigR%hhartree(iband,iband,ik_ibz,:)))*Ha_eV,&
&        SUM(SigR%sigxme  (iband,ik_ibz,:))*Ha_eV,       &
&   REAL(SUM(SigR%sigcmee0(iband,ik_ibz,:)))*Ha_eV,      &
&   REAL(SigR%ze0         (iband,ik_ibz,1)),             &
&   REAL(SUM(SigR%dsigmee0(iband,ik_ibz,:))),            &
&   REAL(SUM(SigR%sigmee  (iband,ik_ibz,:)))*Ha_eV,      &
&   REAL(SigR%degw        (iband,ik_ibz,1))*Ha_eV,       &
&   REAL(SigR%egw         (iband,ik_ibz,1))*Ha_eV,       &
&        SigR%en_qp_diago (iband,ik_ibz,1)*Ha_eV
   call wrtout(unt,msg,mode) 

   write(msg,'(i5,12(2x,f8.3))')                          & 
&         iband,                                          &
&         zero,                                           &
          !sr%e0          (iband,ikibz,isp)*Ha_eV,        &
&         zero,                                           &
&         zero,                                           &
&   AIMAG(SUM(SigR%hhartree(iband,iband,ik_ibz,:)))*Ha_eV,&
&         zero,                                           &
&   AIMAG(SUM(SigR%sigcmee0(iband,ik_ibz,:)))*Ha_eV,      &
&   AIMAG(SigR%ze0         (iband,ik_ibz,1)),             &
&   AIMAG(SUM(SigR%dsigmee0(iband,ik_ibz,:))),            &
&   AIMAG(SUM(SigR%sigmee  (iband,ik_ibz,:)))*Ha_eV,      &
&   AIMAG(SigR%degw        (iband,ik_ibz,1))*Ha_eV,       &
&   AIMAG(SigR%egw         (iband,ik_ibz,1))*Ha_eV,       &
&         zero
   if (verbose/=0) call wrtout(unt,msg,mode) 
  else
   write(msg,'(i5,12(2x,f8.3))')                         & 
&        iband,                                          &
&        en_lda(ik_ibz,iband,isp)*Ha_eV,                 &
         !sr%e0          (iband,ikibz,isp)*Ha_eV,        &
&        SigR%vxcme      (iband,ik_ibz,isp)*Ha_eV,       &
&        SigR%e0         (iband,ik_ibz,isp)*Ha_eV,       &
&   REAL(SigR%hhartree   (iband,iband,ik_ibz,isp))*Ha_eV,&
&        SigR%sigxme     (iband,ik_ibz,isp)*Ha_eV,       &
&   REAL(SigR%sigcmee0   (iband,ik_ibz,isp))*Ha_eV,      &
&   REAL(SigR%ze0        (iband,ik_ibz,isp)),            &
&   REAL(SigR%dsigmee0   (iband,ik_ibz,isp)),            &
&   REAL(SigR%sigmee     (iband,ik_ibz,isp))*Ha_eV,      &
&   REAL(SigR%degw       (iband,ik_ibz,isp))*Ha_eV,      &
&   REAL(SigR%egw        (iband,ik_ibz,isp))*Ha_eV,      &
&        SigR%en_qp_diago(iband,ik_ibz,isp)*Ha_eV
   call wrtout(unt,msg,mode) 

   write(msg,'(i5,12(2x,f8.3))')                         & 
&         iband,                                         &
&         zero,                                          &
          !sr%e0          (iband,ikibz,isp)*Ha_eV,       &
&         zero,                                          &
&         zero,                                          &
&   AIMAG(SigR%hhartree  (iband,iband,ik_ibz,isp))*Ha_eV,&
&         zero,                                          &
&   AIMAG(SigR%sigcmee0   (iband,ik_ibz,isp))*Ha_eV,     &
&   AIMAG(SigR%ze0        (iband,ik_ibz,isp)),           &
&   AIMAG(SigR%dsigmee0   (iband,ik_ibz,isp)),           &
&   AIMAG(SigR%sigmee     (iband,ik_ibz,isp))*Ha_eV,     &
&   AIMAG(SigR%degw       (iband,ik_ibz,isp))*Ha_eV,     &
&   AIMAG(SigR%egw        (iband,ik_ibz,isp))*Ha_eV,     &
&         zero
   if (verbose/=0) call wrtout(unt,msg,mode) 
  end if
 else 
  ! === PAW+U+GW calculation ===
  STOP "PAW+U+GW not yet implemented"
 end if

end subroutine print_Sigma_SC
!!***


!!****f* ABINIT/print_QP
!! NAME
!! print_QP
!!
!! FUNCTION
!! Print in a nice format (?) the expansion coefficients of the quasiparticle 
!! amplitudes in terms of KS eigenvectors
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nsppol=1 for unpolarized, 2 for spin-polarized.
!!  nbnds=total number of bands considered in the calculation.
!!  nkibz=number of irreducible k-points.
!!  m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)=matrix giving the decomposition of the QP
!!   amplitued in the mainfold generated by the KS wavefunctions 
!!   (i.e $ m_lda_to_qp(ib,jb,k,s) := \langle \psi_{ib,k,s}^{KS}| \psi_{jb,k,s}^{QP}\rangle $
!!  ene_qp(nkibz,nbnds,nsppol)= QP energies for each k-point, band and spin.
!!  ib_start,ib_stop=initial and final band index for QP, only states in this range are printed
!!  prtvol
!!  unit
!!  mode_paral
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine print_QP(nbnds,nkibz,nsppol,m_lda_to_qp,ene_qp,ib_start,ib_stop,unit,prtvol,mode_paral)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalaras
!here nbnds should be first dimension!
!scalars
 integer,intent(in) :: ib_start,ib_stop,nbnds,nkibz,nsppol
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
!arrays
 real(dp),intent(in) :: ene_qp(nkibz,nbnds,nsppol)
 complex(dpc),intent(in) :: m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)

!Local variables-------------------------------
!scalars
 integer :: counter,ib_KS,ib_QP,ii,ikibz,isp,nspace,unt,verbose
 character(len=10) :: bks,bqp,kpt,spin
 character(len=4) :: mode
 character(len=500) :: KS_ket,QP_ket,fspace,msg

! *********************************************************************

 unt=std_out ; if (PRESENT(unit))       unt=unit
 verbose=0   ; if (PRESENT(prtvol))     verbose=prtvol
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 do ikibz=1,nkibz 
  call int2char(ikibz,kpt) 
  write(*,*)' k-point = ',ikibz
  do isp=1,nsppol 
   call int2char(isp,spin) 

   do ib_QP=ib_start,ib_stop
    call int2char(ib_QP,bqp) 
    QP_ket=' |'//TRIM(bqp)//';'//TRIM(spin)//'> ='
    write(*,'(a)',ADVANCE='NO')TRIM(QP_ket)
    nspace=LEN(TRIM(QP_ket)) ; write(fspace,'(a,i5,a)')'(/,',nspace,'x)'
   
    counter=0
    do ib_KS=ib_start,ib_stop
     counter=counter+1
     call int2char(ib_KS,bks) 
     KS_ket=' |'//TRIM(bks)//';'//TRIM(spin)//'>'
     write(KS_ket,fmt=10)m_lda_to_qp(ib_QP,ib_KS,ikibz,isp),TRIM(KS_ket)
     write(*,'(a)',ADVANCE='NO')TRIM(KS_ket)
     if (MOD(counter,5)==0) write(*,fspace,ADVANCE='NO')
    end do
    write(*,*)

   end do
  end do 
 end do

 10 FORMAT (1x,2f7.4,a,1x)

end subroutine print_QP
!!***

