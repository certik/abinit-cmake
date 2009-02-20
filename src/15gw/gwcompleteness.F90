!{\src2tex{textfont=tt}}
!!****f* ABINIT/completechi0_deltapart
!! NAME
!! completechi0_deltapart
!!
!! FUNCTION
!!  Apply the delta part of the completeness correction to chi0
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (FB, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see
!! ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ik_bz=Index of the k-point in the full BZ whose contribution has to be added and symmetrized.
!!  qzero=.TRUE. is long wave-length limit.
!!  symchi=1 if we are summing over IBZ_q and symmetrization has to be performed.
!!  npwe=Number of G vectors in chi0.
!!  npwvec=MAX number of G.
!!  nomega=Number of frequencies.
!!  nspinor=Number of spinorial components.
!!  nfftot=Total Number of points in the FFT
!!  ngfft(18)=Info on the FFT.
!!  gvec(3,npwvec)=Reduced coordinates of G.
!!  igfft0(npwvec)=Index of each G in the FFT array.
!!  Gsph_wfn=<Gvectors_type>=Info on the G-sphere for wavefunctions.
!!  Ltg_q=<Little_group>= Structure gathering information on the little group of the external q.
!!  green_enhigh_w=Approximated frequency dependent part of the Green function entering equation (TODO put reference)
!!  wfwfg=Fourier components of u_{kb1}.u_{kb2}
!!
!! OUTPUT
!!  See SIDES EFFECTS
!!
!! SIDES EFFECTS 
!!  chi0(npwe,npwe,nomega)= In input chi0 calculated so far, 
!!  In output the "delta part" of the completeness correction is added.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine completechi0_deltapart(ik_bz,qzero,symchi,npwe,npwvec,nomega,nspinor,&
& nfftot,ngfft,gvec,igfft0,Gsph_wfn,Ltg_q,green_enhigh_w,wfwfg,chi0)

 use defs_basis
 use defs_datatypes
 use m_numeric_tools, only : bisect
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,nfftot,nomega,npwe,npwvec,nspinor,symchi
 logical,intent(in) :: qzero
 type(Gvectors_type),intent(in) :: Gsph_wfn
 type(Little_group),intent(in) :: Ltg_q
!arrays
 integer,intent(in) :: gvec(3,npwvec),igfft0(npwvec),ngfft(18)
 complex(dpc),intent(in) :: green_enhigh_w(nomega)
 complex(gwpc),intent(in) :: wfwfg(nfftot*nspinor**2)
 complex(gwpc),intent(inout) :: chi0(npwe,npwe,nomega)

!Local variables ------------------------------
!scalars
 integer,save :: enough=0
 integer :: iSm1_g1mg2,iSm1_g1mg2_fft,ig,ig1mg2,ige,igmgp,igmgpx,igmgpy,igmgpz
 integer :: igp,igs,igstart,ish,ishbsc,isym,itim,nshwfn,outofGwfn
 real(dp) :: difflen
 complex(gwpc) :: phmGt
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: g1mg2(3)
 real(dp) :: gmet(3,3)

!************************************************************************

 if (qzero) then
  igstart=2
 else
  igstart=1
 end if
 outofGwfn=0

 SELECT CASE (symchi)

 CASE (0)
  ! === Do not use symmetries ===
  ! MG: Here there was a problem since one has to make sure 
  ! that each G1-G2 is still in the FFT mesh for each G1 and G2 in chi0 (not always true)
  ! MODULO and igmgp wraps G1-G2 in the FFT box but the Fourier components are not periodic!
  do igp=igstart,npwe
   do ig=igstart,npwe

    g1mg2(:)=gvec(:,ig)-gvec(:,igp)
    if (ANY(g1mg2(:)>ngfft(1:3)/2 .or. g1mg2(:)<-(ngfft(1:3)-1)/2)) then 
     outofGwfn=outofGwfn+1
     CYCLE
    end if

    igmgpx=MODULO(g1mg2(1),ngfft(1))
    igmgpy=MODULO(g1mg2(2),ngfft(2))
    igmgpz=MODULO(g1mg2(3),ngfft(3))
    igmgp= 1 + igmgpx + igmgpy*ngfft(1) + igmgpz*ngfft(1)*ngfft(2)
    !MG this is always false!
    !if (igmgp>nfftot) then
    ! write(*,*) 'completechi0_deltapart: too large G in the extrapolar correction'
    ! cycle
    !end if
    chi0(ig,igp,:) = chi0(ig,igp,:) + wfwfg(igmgp)*green_enhigh_w(:) 
   end do
  end do

 CASE (1)
  !
  ! * <Sk b|e^{-i(G1-G2}.r}|b Sk> = e^{-i(G1-G2).\tau} <k b|e^{-i(S^{-1}(G1-G2).r)|b k>
  ! * green_enhigh_w in invariant under symmetry
  ! * We symmetrize using the operations of the little group of q since this routine
  !   is called inside a sum over IBZ_q, it would be possible to symmetrize 
  !   this term by just summing over the IBZ and rotating the matrix elements.
  ! * Time-reversal does not lead to a complex conjugated since bra and ket are the same.
  !
  nshwfn   = Gsph_wfn%nsh 
  gmet(:,:)= Gsph_wfn%gmet(:,:)

  do igp=igstart,npwe
   do ig=igstart,npwe

    ! === Find g1mg2 in the array gvec, if out of range skip it ===
    ! * Use shells and bisect to find the starting index thus avoiding storing a table (ig1,ig2)
    g1mg2(:)=gvec(:,ig)-gvec(:,igp)
    difflen = two_pi*SQRT(DOT_PRODUCT(g1mg2,MATMUL(gmet,g1mg2)))

    ishbsc = bisect(Gsph_wfn%shlen,difflen)
    if (ishbsc==nshwfn) then 
     outofGwfn=outofGwfn+1
     CYCLE
    end if

    igs=Gsph_wfn%shlim(ishbsc)
    ige=Gsph_wfn%shlim(MIN(ishbsc+2,nshwfn+1))-1

    ig1mg2=igs-1 ; found=.FALSE.
    do while (.not.found .and. ig1mg2<ige)
     ig1mg2=ig1mg2+1
     found=(ALL(gvec(:,ig1mg2)==g1mg2(:)))
    end do
    if (.not.found) then 
     call assert(.FALSE.,'g1mg2 not found in gvec',__FILE__,__LINE__)
    end if

    do itim=1,Ltg_q%timrev
     do isym=1,Ltg_q%nsym_sg

      if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then 
       ! * This operation belongs to the little group and has to be used to reconstruct the BZ.
       ! * Time-reversal in not used to rotate (G1-G2) see comment above.
       phmGt          = Gsph_wfn%phmGt  (ig1mg2,isym) 
       iSm1_g1mg2     = Gsph_wfn%rottbm1(ig1mg2,1,isym)
       iSm1_g1mg2_fft = igfft0(iSm1_g1mg2)

       chi0(ig,igp,:) = chi0(ig,igp,:) + phmGt*wfwfg(iSm1_g1mg2_fft)*green_enhigh_w(:) 
      end if

     end do !isym
    end do !itim

   end do !igp
  end do !ig

 CASE DEFAULT
  call assert(.FALSE.,'Wrong value of symchi',__FILE__,__LINE__)
 END SELECT

 if (outofGwfn/=0) then 
  enough=enough+1
  if (enough<=50) then 
   write(msg,'(3a,i5)')ch10,&
&  ' completechi0_deltapart : WARNING ',&
&  ' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofGwfn
   call wrtout(std_out,msg,'COLL') 
   if (enough==50) then 
    write(msg,'(a)')' ========== Stop writing Warnings =========='
    call wrtout(std_out,msg,'COLL') 
   end if
  end if
 end if

end subroutine completechi0_deltapart
!!***

!!****f* ABINIT/output_chi0sumrule
!! NAME
!! output_chi0sumrule
!!
!! FUNCTION
!!  Calculate and output the value of the sum rule for 
!!  the non-interacting polarizability chi0
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (FB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see
!! ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  (for writing routines, no output)
!!  otherwise, should be described
!!
!! NOTES
!!
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine output_chi0sumrule(qeq0,iq,npwe,omegaplasma,chi0sumrule,epsm1_w0,vc_sqrt)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq,npwe
 real(dp),intent(in) :: omegaplasma
 logical,intent(in) :: qeq0
!arrays
 real(dp),intent(inout) :: chi0sumrule(npwe)
 complex(gwpc),intent(in) :: epsm1_w0(npwe,npwe),vc_sqrt(npwe)

!Local variables ------------------------------
!scalars
 integer :: ig,igstart
 real(dp) :: average,norm
 character(len=500) :: msg

!************************************************************************

 if (qeq0) then 
  igstart=2
 else
  igstart=1
 end if
 !
 ! The sumrule is for
 ! \int d\omega \omega v * Im[ \chi_0(\omega) ] = \pi/2 * w_p^2
 chi0sumrule(igstart:npwe) = chi0sumrule(igstart:npwe) * vc_sqrt(igstart:npwe)**2
 !
 ! Calculate a weighted average of the fulfilment of the sumrule on epsilon
 ! The weight is given according to the significance of each q+G in the
 ! subsequent GW calculation
 ! It is proportional to v * (epsm1 -1 )
 average = zero
 norm    = zero
 do ig=igstart,npwe
  average = average + chi0sumrule(ig) * real( vc_sqrt(ig)**2 * (epsm1_w0(ig,ig) - 1.0_dp ) )
  norm    = norm    +                   real( vc_sqrt(ig)**2 * (epsm1_w0(ig,ig) - 1.0_dp ) )
!  average = average + chi0sumrule(ig) * real(  (epsm1_w0(ig,ig) - 1.0_dp ) )
!  norm    = norm    +                   real(  (epsm1_w0(ig,ig) - 1.0_dp ) )
!  write(203,'(i4,8(2x,e12.6))') ig,1.0_dp/vc_sqrt(ig),chi0sumrule(ig)/ (0.5d0*omegaplasma**2*pi)
 end do
 write(msg,'(1x,a,i4,a,f10.2,2x,a)')&
& ' Average fulfillment of the sum rule on Im[epsilon] for q-point ',&
& iq,' :',average/norm/(0.5_dp*omegaplasma**2*pi)*100.0_dp,'[%]'
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out, msg,'COLL') 

end subroutine output_chi0sumrule
!!***


!!****f* ABINIT/accumulate_chi0sumrule
!! NAME
!! accumulate_chi0sumrule
!!
!! FUNCTION
!!  Accumulate the contribution to the sum rule for Im chi0
!!  arising from a single transition. Eventually symmetrize 
!!  it using the symmetry operations of the little group of q.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (FB, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see
!! ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ik_bz=Index of the k-point in the full BZ whose contribution has to be added and symmetrized.
!!  symchi=1 if we are summing over IBZ_q and symmetrization has to be performed.
!!  npwe=Number of G vectors in chi0.
!!  npwepG0=Number of G vectors in the "enlarged" sphere to treat umklapp.
!!  factor=factor entering the expression.
!!  delta_ene=Transition energy.
!!  Ltg_q=<Little_group>= Structure gathering information on the little group of the external q.
!!  Gsph_epsG0=<Gvectors_type>=Info on the G-sphere for chi0.
!!  rhotwg(npwepG0)=Fouriet transform of u_{b1 k-q} u_{b2 k} in the "enlarged" sphere. 
!!
!! OUTPUT
!!  See SIDES EFFECTS
!!
!! SIDES EFFECTS 
!!  chi0sumrule(npwe)= In input the sum rule calculated so far, 
!!  In output the contribution of this transition is accounted for, and, eventually, symmetrized.
!!  using the symmetry operations of the little group of the external q.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine accumulate_chi0sumrule(ik_bz,symchi,npwe,factor,delta_ene,&
& Ltg_q,Gsph_epsG0,npwepG0,rhotwg,chi0sumrule)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,npwe,npwepG0,symchi
 real(dp),intent(in) :: delta_ene,factor
 type(Gvectors_type),intent(in) :: Gsph_epsG0
 type(Little_group),intent(in) :: Ltg_q
!arrays
 real(dp),intent(inout) :: chi0sumrule(npwe)
 complex(gwpc),intent(in) :: rhotwg(npwepG0)

!Local variables-------------------------------
!scalars
 integer :: isym,itim
 character(len=500) :: msg
!arrays
 integer,allocatable :: Sm1_gmG0(:)
 integer,pointer :: gmG0(:)
 complex(gwpc),allocatable :: rhotwg_sym(:)

!************************************************************************

 ! Accumulating the sum rule on chi0.
 ! Eq.(5.284) in G. D. Mahan Many-Particle Physics 3rd edition

 SELECT CASE (symchi)

 CASE (0)
  chi0sumrule(:)=chi0sumrule(:) + factor*delta_ene*ABS(rhotwg(1:npwe))**2

 CASE (1)

  allocate(rhotwg_sym(npwe))
  allocate(Sm1_gmG0(npwe))

  do itim=1,Ltg_q%timrev
   do isym=1,Ltg_q%nsym_sg

    if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then 
     ! === This operation belongs to the little group and has to be used to reconstruct the BZ ===
     ! * In the following 2 lines mind the slicing (1:npwe)

     gmG0  => Ltg_q%igmG0(1:npwe,itim,isym)  
     Sm1_gmG0(1:npwe)=Gsph_epsG0%rottbm1(gmG0(1:npwe),itim,isym)
     rhotwg_sym(1:npwe)=rhotwg(Sm1_gmG0) 

     chi0sumrule(:)=chi0sumrule(:) + factor*delta_ene*ABS(rhotwg_sym(1:npwe))**2
    end if

   end do !isym
  end do !itim 

  deallocate(rhotwg_sym)
  deallocate(Sm1_gmG0)

 CASE DEFAULT
  write(msg,'(4a)')ch10,&
&  ' accumulate_chi0sumrule : BUG - ',ch10,&
&  ' wrong value for symchi '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 END SELECT

end subroutine accumulate_chi0sumrule
!!***
