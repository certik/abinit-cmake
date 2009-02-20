!{\src2tex{textfont=tt}}
!!****f* ABINIT/solve_Dyson
!! NAME
!! solve_Dyson
!!
!! FUNCTION
!!  Solve the Dyson equation for the QP energies. Two different methods are coded:
!!  The first one is based on the standard perturbative approach in which the self-energy
!!  is linearly expanded around the previous single-particle energy (KS energy if one-shot)
!!  and the derivative is evaluated by finite differences.
!!  In the second method (AC), the values of the self-energy operator on the real axis are obtained 
!!  by means of an analitic continuation based on the Pade extrapolation. 
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ikcalc=Index of the considered k-point in the Sp%kcalc array.
!!  nomega_sigc=Number of frequencies used to evaluate the correlation part of Sigma.
!!  Sp<Sigma_parameters>=Structure gathering parameters on the calculation of Sigma.
!!     %nbnds=Number of bands in G0.
!!     %nsppol=Number of independent spin polarizations.
!!     %nsig_ab=Numner of components in the self-energy operator.
!!     %minbnd and %maxbnd= min and Max band index for GW correction (for this k-point)
!!     %nomegasr=Number of real frequencies for spectral function.
!!     %nomegasrd=Number of real frequencies used to evalute the derivative of Sigma.
!!     %nomegasi=Number of imaginary frequencies for AC.
!!     %omegasi=Purely imaginary frequencies for AC.
!!     %gwcalctyp=Type of the GW calculation.
!!     %soenergy=Scissor energy
!!  Kmesh<BZ_mesh_type>=Info on the K-mesh for the wavefunctions.
!!     %nibz=Number of points in the IBZ
!!  sigxme_tmp(ib1:ib2,ib1:ib2,nsppol)=Matrix elements of Sigma_x.
!!  sigcme_tmp=(nomega_sigc,ib1:ib2,ib1:ib2,nsppol)=Matrix elements of Sigma_c.
!!  en_qp(nibz,nbnds,nsppol)= KS or QP energies, only used in case of calculation with scissor operator.
!!
!! OUTPUT
!!  Sr<Sigma_results>=Structure containing the matrix elements of the self-energy:
!!     %sigxme(ib1:ib2,jkibz,nsspol)=Diagonal elements of Sigma_x
!!     %sigcmee0(ib1:ib2,jkibz,nsppol)=Matrix elements of Sigma_c at the initial energy E0.
!!     %dsigmee0(jb,ib1:ib2,nsppol)=Derivate of sigma at the energy E0.
!!     %ze0(ib1:ib2,jkibz,is)=Renormalization factor at the energy E0.
!!     %degw(ib1:ib2,jkibz,is)= QP correction  i.e DeltaE_GW=E-E0 
!!     %egw(ib1:ib2,jkibz,is)=QP energy
!!     %sigmee(ib1:ib2,jkibz,is)=Self-energy evaluated at the QP energy.
!!     %sigcme (ib1:ib2,jkibz,io,is)= Sigma_c as a function of frequency.
!!     %sigxcme(ib1:ib2,jkibz,io,is)= Sigma_xc as a function of frequency.
!!     %sigcmesrd (ib1:ib2,jkibz,io,is)= Diagonal matrix elements of \Sigma_c  at frequencies around the KS eigenvalue
!!     %sigxcmesrd(ib1:ib2,jkibz,io,is)= Diagonal matrix elements of \Sigma_xc at frequencies around the KS eigenvalue
!!    where ib1 and ib2 are the band indeces included in the GW calculation for this k-point.
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

subroutine solve_Dyson(ikcalc,nomega_sigc,Sp,Kmesh,sigxme_tmp,sigcme_tmp,en_qp,Sr)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : czero_gw
 use m_numeric_tools, only : linfit, pade, dpade, newrap_step
 use m_errors, only : assert
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => solve_Dyson
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikcalc,nomega_sigc
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Sigma_parameters),intent(in) :: Sp
 type(Sigma_results),intent(inout) :: Sr
!arrays
!no_abirules
 real(dp),intent(in) :: en_qp(Kmesh%nibz,Sp%nbnds,Sp%nsppol)
 complex(gwpc),intent(in) :: &
& sigcme_tmp(nomega_sigc,Sp%minbnd(ikcalc):Sp%maxbnd(ikcalc),Sp%minbnd(ikcalc):Sp%maxbnd(ikcalc),Sp%nsppol*Sp%nsig_ab)
 complex(gwpc),intent(in) :: &
& sigxme_tmp(Sp%minbnd(ikcalc):Sp%maxbnd(ikcalc),Sp%minbnd(ikcalc):Sp%maxbnd(ikcalc),Sp%nsppol*Sp%nsig_ab)

!Local variables-------------------------------
!scalars
 integer,parameter :: NITER_MAX=1000
 integer :: iab,ib,ib1,ib2,ikbz_gw,info,io,ioe0j,is,is_idx,isym,iter,itim,jb
 integer :: jkibz,kb,ld_matrix,lwork,mod10,nsploop
 real(dp) :: alpha,beta,smrt
 complex(dpc) :: ctdpc,dct,dsigc,sigc,zz
 complex(gwpc) :: phase
 logical :: converged,ltest
 character(len=500) :: msg
!arrays
 real(dp) :: kbz_gw(3)
 real(dp),allocatable :: e0pde(:),eig(:),rwork(:),scme(:)
 complex(dpc),allocatable :: hdp(:,:),tmpcdp(:),work(:)
 complex(gwpc),allocatable :: hhartree(:,:,:),htotal(:,:,:)

! *************************************************************************

#ifdef __VMS
!DEC$ ATTRIBUTES ALIAS:'ZHEEV' :: zheev
#endif

#if defined DEBUG_MODE
 write(msg,'(a)')' solve_Dyson : enter '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 mod10=MOD(Sp%gwcalctyp,10)

 ltest=(nomega_sigc==Sp%nomegasr+Sp%nomegasrd)
 if (mod10==1) ltest=(nomega_sigc==Sp%nomegasi)
 call assert(ltest,'Wrong number of frequencies',__FILE__,__LINE__)

 ! Index of the KS or QP energy.
 ioe0j=Sp%nomegasrd/2+1

 ! min and Max band index for GW corrections (for this k-point).
 ib1=Sp%minbnd(ikcalc) 
 ib2=Sp%maxbnd(ikcalc)

 ! Index of this k-point in the IBZ array.
 ikbz_gw=Sp%kcalc(ikcalc) 
 call get_BZ_item(Kmesh,ikbz_gw,kbz_gw,jkibz,isym,itim,phase)

 ! ===========================================================
 ! ==== Solve the Dyson Equation and store results in Sr% ====
 ! ===========================================================

 ! === Save elements or ab components of Sigma_x (hermitian) ===
 ! TODO It should be hermitian also in the spinorial case re-check
 do is=1,Sp%nsppol
  do jb=ib1,ib2

   do iab=1,Sp%nsig_ab
    is_idx=is ; if (Sp%nsig_ab>1) is_idx=iab
    Sr%sigxme(jb,jkibz,is_idx) = REAL(sigxme_tmp(jb,jb,is_idx))
   end do
   if (Sp%nsig_ab>1) then 
    write(std_out,'(i3,4f8.3,a,f8.3)')jb,Sr%sigxme(jb,jkibz,:)*Ha_eV,' Tot ',SUM(Sr%sigxme(jb,jkibz,:))*Ha_eV
   end if

  end do 
 end do 

 if (mod10/=1) then 
  ! =============================
  ! === Perturbative approach ===
  ! =============================
  do is=1,Sp%nsppol
   do jb=ib1,ib2

    ! === Get matrix elements of Sigma_c at energy E0 ===
    ! * SigC(w) is linearly interpolated and the slope alpha is assumed as dSigC/dE
    do iab=1,Sp%nsig_ab
     is_idx=is ; if (Sp%nsig_ab>1) is_idx=iab

     Sr%sigcmee0(jb,jkibz,is_idx) = sigcme_tmp(Sr%nomega+ioe0j,jb,jb,is_idx)

     allocate(scme(Sp%nomegasrd),e0pde(Sp%nomegasrd))
     e0pde(:) = Sr%omegasrd(jb,jkibz,:,is)
     scme(:)  = REAL(sigcme_tmp(Sr%nomega+1:Sr%nomega+Sp%nomegasrd,jb,jb,is_idx))

     if (Sp%nomegasrd==1) then
      smrt=zero ; alpha=zero
     else
      smrt=linfit(Sp%nomegasrd,e0pde(:),scme(:),alpha,beta)
     end if
     if (smrt>0.1/Ha_eV) then
      write(msg,'(6a,i4,a,i3,2a,2(f22.15,2a))')ch10,&
&      ' solve_Dyson : WARNING - ',ch10,&
&      '  Values of Re Sig_c are not linear ',ch10,&
&      '  band index       = ',jb,' spin|component = ',is_idx,ch10,& 
&      '  root mean square = ',smrt,ch10,&
&      '  estimated slope  = ',alpha,ch10,' Omega [eV] SigC [eV]' 
      call wrtout(std_out,msg,'COLL') 
      do io=1,Sp%nomegasrd
       write(std_out,'(2f8.4)')e0pde(io)*Ha_eV,scme(io)*Ha_eV
      end do
     end if
     deallocate(scme,e0pde)
     !   
     ! === Evaluate renormalization factor and QP correction ===
     ! * Z=(1-dSigma/domega(E0))^-1
     ! * DeltaE_GW=E-E0= (Sigma(E0)-V_xc)/(1-dSigma/domega)
     ! * If nspinor==2, this part is done at the end.

     Sr%dsigmee0(jb,jkibz,is_idx)=CMPLX(alpha,zero)

     if (Sp%nsig_ab==1) then

      Sr%ze0(jb,jkibz,is)=one/(one-Sr%dsigmee0(jb,jkibz,is))
      Sr%degw(jb,jkibz,is) = &
&      (Sr%sigxme(jb,jkibz,is)+Sr%sigcmee0(jb,jkibz,is)-Sr%e0(jb,jkibz,is)+Sr%hhartree(jb,jb,jkibz,is))*Sr%ze0(jb,jkibz,is) 

      if (Sp%soenergy>0.1d-4) then
       ! RS: PATCH for GW+scissor: e0 is replaced by en_qp which contains the updated energy eigenvalue
       Sr%degw(jb,jkibz,is)= &
&      (Sr%sigxme(jb,jkibz,is)+Sr%sigcmee0(jb,jkibz,is)-en_qp(jkibz,jb,is)+Sr%hhartree(jb,jb,jkibz,is))*Sr%ze0(jb,jkibz,is) 
      end if

      Sr%egw(jb,jkibz,is)=Sr%e0(jb,jkibz,is)+Sr%degw(jb,jkibz,is)

      if (Sp%soenergy>0.1d-4) then
       ! RS: PATCH for GW+scissor: e0 is replaced by en_qp which contains the updated energy eigenvalue    
       Sr%egw(jb,jkibz,is)=en_qp(jkibz,jb,is)+Sr%degw(jb,jkibz,is)
      end if

      ! === Estimate Sigma at the QP-energy ===
      ! * Sigma(E_qp)=Sigma(E0)+(E_qp-E0)*dSigma/dE
      Sr%sigmee(jb,jkibz,is)= &
&      Sr%sigxme(jb,jkibz,is)+Sr%sigcmee0(jb,jkibz,is)+Sr%degw(jb,jkibz,is)*Sr%dsigmee0(jb,jkibz,is)

      if (Sp%soenergy>0.1d-4) then
       ! RS: here we report the gw corr with respect to e0 in the output file
       Sr%degw(jb,jkibz,is)=Sr%egw(jb,jkibz,is)-Sr%e0(jb,jkibz,is)
      end if

     end if !Sp%nsig_ab==1

     ! Spectrum of Sigma
     do io=1,Sr%nomega
      Sr%sigcme (jb,jkibz,io,is_idx)= sigcme_tmp(io,jb,jb,is_idx)
      Sr%sigxcme(jb,jkibz,io,is_idx)= Sr%sigxme(jb,jkibz,is_idx)+Sr%sigcme(jb,jkibz,io,is_idx)
     end do
     do io=1,Sp%nomegasrd
      Sr%sigcmesrd (jb,jkibz,io,is_idx)= sigcme_tmp(Sr%nomega+io,jb,jb,is_idx)
      Sr%sigxcmesrd(jb,jkibz,io,is_idx)= Sr%sigxme(jb,jkibz,is_idx)+Sr%sigcmesrd(jb,jkibz,io,is_idx)
     end do

    end do !iab

    if (Sp%nsig_ab>1) then
     ltest=(ABS(Sp%soenergy)<0.1d-4)
     call assert(ltest,'Scissor with spinor not coded yet',&
&     __FILE__,__LINE__)
     !TODO this should be allocated with nsppol
     !TODO recheck this part

     ! === Evaluate renormalization factor and QP correction ===
     ! * Z=(1-dSigma/domega(E0))^-1
     ! * DeltaE_GW=E-E0= (Sigma(E0)-V_xc)/(1-dSigma/domega)
     write(std_out,'(a,i2,10f8.3)')' Correlation',jb,Sr%sigcmee0(jb,jkibz,:)*Ha_eV,SUM(Sr%sigcmee0(jb,jkibz,:))*Ha_eV

     Sr%ze0 (jb,jkibz,1) = one/(one-SUM(Sr%dsigmee0(jb,jkibz,:)))
     Sr%degw(jb,jkibz,1) = Sr%ze0(jb,jkibz,1) * &
&    (SUM(Sr%sigxme(jb,jkibz,:)+Sr%sigcmee0(jb,jkibz,:)+Sr%hhartree(jb,jb,jkibz,:))-Sr%e0(jb,jkibz,1))
!&   (SUM(Sr%sigxme(jb,jkibz,:)+Sr%sigcmee0(jb,jkibz,:)+Sr%hhartree(jb,jb,jkibz,:))-two*Sr%e0(jb,jkibz,1))
!&   (SUM(Sr%sigxme(jb,jkibz,:)+Sr%sigcmee0(jb,jkibz,:)-Sr%vxcme(jb,jkibz,:)))*Sr%ze0(jb,jkibz,1) 
       
     !write(77,'(2f8.3)')SUM(Sr%hhartree(jb,jb,jkibz,:))-2*Sr%e0(jb,jkibz,1)+SUM(Sr%vxcme(jb,jkibz,:))

     Sr%egw(jb,jkibz,1)=Sr%e0(jb,jkibz,1)+Sr%degw(jb,jkibz,1)

     ! === Estimate Sigma at the QP-energy ===
     do iab=1,Sp%nsig_ab
      Sr%sigmee(jb,jkibz,iab)= &
&      Sr%sigxme(jb,jkibz,iab)+Sr%sigcmee0(jb,jkibz,iab)+Sr%degw(jb,jkibz,1)*Sr%dsigmee0(jb,jkibz,iab)
     end do

    end if

   end do !jb
  end do !is

 else 
  ! =============================
  ! === Analytic Continuation ===
  ! =============================
  call assert((Sp%nsig_ab==1),'AC with spinor not implemented',&
&  __FILE__,__LINE__)
  do is=1,Sp%nsppol
   do jb=ib1,ib2

    allocate(tmpcdp(Sp%nomegasi))
    ! * Calculate Sigc(E0), dSigc(E0)
    zz=CMPLX(Sr%e0(jb,jkibz,is),zero)
    if (Sp%soenergy>0.1d-4) then
     ! RS: e0 is replaced by en_qp which contains the updated energy eigenvalue
     zz=CMPLX(en_qp(jkibz,jb,is),zero)   
    end if

    ! === Diagonal elements of sigcme_tmp ===
    ! * if zz in 2 or 3 quadrant, avoid poles in the complex plane using Sigma(-iw)=Sigma(iw)*.
    do iab=1,Sp%nsig_ab
     is_idx=is ; if (Sp%nsig_ab>1) is_idx=iab
     if (REAL(zz)>zero) then
      tmpcdp(:)=sigcme_tmp(:,jb,jb,is_idx)
      Sr%sigcmee0(jb,jkibz,is_idx)=  pade(Sp%nomegasi,Sp%omegasi,tmpcdp,zz)
      Sr%dsigmee0(jb,jkibz,is_idx)= dpade(Sp%nomegasi,Sp%omegasi,tmpcdp,zz)
     else
      tmpcdp(:)=CONJG(sigcme_tmp(:,jb,jb,is_idx))
      Sr%sigcmee0(jb,jkibz,is_idx)=  pade(Sp%nomegasi,CONJG(Sp%omegasi),tmpcdp,zz)
      Sr%dsigmee0(jb,jkibz,is_idx)= dpade(Sp%nomegasi,CONJG(Sp%omegasi),tmpcdp,zz)
     end if
    end do !iab

    ! Z=(1-dSigma/domega(E0))^-1
    if (Sp%nsig_ab==1) then
     Sr%ze0(jb,jkibz,is)=one/(one-Sr%dsigmee0(jb,jkibz,is))
    else
     Sr%ze0(jb,jkibz,1)=one/(one-SUM(Sr%dsigmee0(jb,jkibz,:)))
    end if

    ! Find roots of E^0-V_xc-V_U+Sig_x+Sig_c(z)-z, i.e E^qp.
    ! using Newton-Raphson method and starting point E^0
    zz=CMPLX(Sr%e0(jb,jkibz,is),zero)

    if (Sp%soenergy>0.1d-4) then
     ! RS: e0 is replaced by en_qp which contains the updated energy eigenvalue
     zz=CMPLX(en_qp(jkibz,jb,is),0.0)
    end if

    iter=0 ; converged=.FALSE. ; ctdpc=cone
    do while (ABS(ctdpc)>0.0001/Ha_eV.or.iter<NITER_MAX)
     iter=iter+1
     if (REAL(zz)>zero) then
      tmpcdp(:)=sigcme_tmp(:,jb,jb,is)
      sigc =  pade(Sp%nomegasi,Sp%omegasi,tmpcdp,zz)
      dsigc= dpade(Sp%nomegasi,Sp%omegasi,tmpcdp,zz)
     else
      tmpcdp(:)=CONJG(sigcme_tmp(:,jb,jb,is))
      sigc =  pade(Sp%nomegasi,CONJG(Sp%omegasi),tmpcdp,zz)
      dsigc= dpade(Sp%nomegasi,CONJG(Sp%omegasi),tmpcdp,zz)
     end if
     ctdpc=Sr%e0(jb,jkibz,is)-Sr%vxcme(jb,jkibz,is)-Sr%vUme(jb,jkibz,is)+Sr%sigxme(jb,jkibz,is)+sigc-zz
     if (ABS(ctdpc)<0.0001/Ha_eV) then 
      converged=.TRUE.
      EXIT
     end if
     dct=dsigc-one
     zz=newrap_step(zz,ctdpc,dct)
    end do

    if (.not.converged) then 
     write(msg,'(4a,f8.4)')ch10,&
&     ' solve_Dyson : WARNING - ',ch10,&
&     ' problem in converging: ABS(ctdpc)= ',ABS(ctdpc)
     call wrtout(std_out,msg,'COLL') 
    end if
    !   
    ! Store the final result TODO re-shift everything according to efermi
    Sr%egw(jb,jkibz,is)=zz
    Sr%degw(jb,jkibz,is)=Sr%egw(jb,jkibz,is) - Sr%e0(jb,jkibz,is)
    Sr%sigmee(jb,jkibz,is)=Sr%sigxme(jb,jkibz,is) + sigc
    !   
    ! Spectra of Sigma, remember that Sr%nomega does not contains the frequencies used to evaluate the derivative
    ! each frequency is obtained using the pade_expression
    do io=1,Sr%nomega
     zz=Sp%omegasf(io)
     if (REAL(zz)>zero) then
      tmpcdp(:)=sigcme_tmp(:,jb,jb,is)
      Sr%sigcme(jb,jkibz,io,is) = pade(Sp%nomegasi,Sp%omegasi,tmpcdp,zz)
     else
      tmpcdp(:)=CONJG(sigcme_tmp(:,jb,jb,is))
      Sr%sigcme(jb,jkibz,io,is) = pade(Sp%nomegasi,CONJG(Sp%omegasi),tmpcdp,zz)
     end if
     Sr%sigxcme(jb,jkibz,io,is)= Sr%sigxme(jb,jkibz,is)+Sr%sigcme(jb,jkibz,io,is)
    end do
    !   
    ! === Save sigma values along the imaginary axis ===
    do iab=1,Sp%nsig_ab
     is_idx=is ; if (Sp%nsig_ab>1) is_idx=iab
     do io=1,Sp%nomegasi
      Sr%sigcmesi (jb,jkibz,io,is_idx)= sigcme_tmp(io,jb,jb,is_idx)
      Sr%sigxcmesi(jb,jkibz,io,is_idx)= Sr%sigxme(jb,jkibz,is_idx)+Sr%sigcmesi(jb,jkibz,io,is_idx)
     end do
    end do
 
    deallocate(tmpcdp)
    !if (rank==master) call print_Sigma_perturbative(Sr,jkibz,jb,is,Dtset%usepawu,prtvol=1)

   end do !jb
  end do !is
 end if ! Analytical continuation
 !
 ! === Diagonalize the QP Hamiltonian (forced to be Hermitian) ===
 ! * Calculate Sr%en_qp_diago and Sr%eigvec_qp to be written in the QPS file.
 ! TODO in case of AC results are wrong.

 allocate(hhartree(ib1:ib2,ib1:ib2,Sp%nsppol*Sp%nsig_ab))
 hhartree(:,:,:)=Sr%hhartree(ib1:ib2,ib1:ib2,jkibz,:)

 ! If non self-consistent erase all off-diagonal elements
 if (Sp%gwcalctyp<20) then
  do jb=ib1,ib2
   do kb=ib1,ib2
    if (jb==kb) CYCLE
    hhartree(jb,kb,:)=czero_gw
   end do
  end do
 end if

 allocate(htotal(ib1:ib2,ib1:ib2,Sp%nsppol*Sp%nsig_ab))
 htotal(:,:,:) = hhartree(:,:,:)+sigxme_tmp(:,:,:)+sigcme_tmp(Sr%nomega+ioe0j,:,:,:)

 ! === Get the Hermitian part of htotal ===
 ! * In the noncollinear case A_{12}^{ab} = A_{21}^{ba}^* if A is Hermitian.
 nsploop=Sp%nsppol ; if (Sp%nsig_ab/=1) nsploop=2
 do is=1,nsploop
  htotal(:,:,is)= half*(htotal(:,:,is)+TRANSPOSE(CONJG(htotal(:,:,is))))
 end do 

 if (Sp%nsig_ab==4) then
  htotal(:,:,3)= half*(htotal(:,:,3)+TRANSPOSE(CONJG(htotal(:,:,4))))
  htotal(:,:,4)= TRANSPOSE(CONJG(htotal(:,:,3)))
 end if

 ! === Solve Herm(htotal)*U = E*U ===
 ! TODO: call print_QP
 ld_matrix=ib2-ib1+1 ; lwork=2*ld_matrix-1
 allocate(hdp(ld_matrix,ld_matrix),eig(ld_matrix))
 allocate(work(lwork),rwork(3*ld_matrix-2))

 do is=1,Sp%nsppol
  if (Sp%nsig_ab==1) then
   hdp(:,:)=htotal(ib1:ib2,ib1:ib2,is)
  else 
   hdp(:,:)=SUM(htotal(ib1:ib2,ib1:ib2,:),DIM=3)
  end if

  call ZHEEV('V','U',ld_matrix,hdp,ld_matrix,eig,work,lwork,rwork,info)
  write(msg,'(a,i3)')' ZHEEV reported info = ',info
  call assert((info==0),msg,__FILE__,__LINE__)

  Sr%eigvec_qp(ib1:ib2,ib1:ib2,jkibz,is)=hdp(:,:)
  Sr%en_qp_diago(ib1:ib2,jkibz,is)=eig(:)
 end do 

 deallocate(hdp,eig,work,rwork)
 deallocate(htotal,hhartree)

#if defined DEBUG_MODE
 write(msg,'(a)')' solve_Dyson : exit '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine solve_Dyson
!!***
