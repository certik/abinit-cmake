!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_vHxc_braket
!! NAME
!!  calc_vHxc_braket
!!
!! FUNCTION
!!  Evaluate the matrix elements of $v_H$ and $v_{xc}$ and $v_U$ 
!!  both in case of NC pseudopotentials and PAW (LDA+U, presently, is only available in PAW)
!!  The matrix elements of $v_{xc}$ are calculated with and without the core contribution. 
!!  The later quantity is required in case of GW calculations.
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  b1,b2=min and max band index to be considered in vxc_me and vhartr_me 
!!  gsqcutf_eff=Fourier cutoff on G^2 for "large sphere" of radius double
!!   that of the basis sphere--appropriate for charge density rho(G),Hartree potential, and pseudopotentials
!!  Dtset <type(dataset_type)>=all input variables in this dataset
!!     %nspden= number of spin-density components
!!     %nspinor=number of spinorial components
!!     %nsppol=number of independent spin polarizations
!!  igfftf(npwvec)=indeg of each G vector in the fine FFT mesh
!!  kcalc2ibz(nkcalc)=index in the irred array of each nkcalc k-point to be considered
!!  MPI_enreg=informations about MPI parallelization
!!  nbnds=number of bands treated (note: does not depend on k)
!!  ngfftf(18)contain all needed information about 3D fine FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nfftf=number of points in the fine FFT mesh (for this processor)
!!  nfftf_tot= Total number of points in the fine FFT mesh (for this processor)
!!  nkcalc=number of points to be calculated
!!  nkibz=number of irreducible k-points
!!  npwvec=Max number of planewaves
!!  Pawtab(Dtset%ntypat*Dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  Paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  Pawang <type(pawang_type)>=paw angular mesh and related data
!!  Paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  Pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  Pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  Cprj(Cryst%natom,Dtset%nspinor*nbnds*nkibz*Dtset%nsppol*Dtset%usepaw) <type(Cprj_type)
!!   projected input wave functions <Proj_i|Cnk> with all NL projectors for each k-point in the IBZ.
!!  Cryst<Crystal_structure>=unit cell and symmetries
!!     %natom=number of atoms in the unit cell
!!     %rprimd(3,3)=direct lattice vectors 
!!     %ucvol=unit cell volume
!!     %ntypat= number of type of atoms 
!!     %typat(natom)=type of each atom
!!  vhartr(nfftf)= Hartree potential in real space on the fine FFT mesh
!!  vxc(nfftf,Dtset%nspden)= xc potential in real space on the fine FFT grid
!!  Wf_info <type (wavefunctions_information)>=Structure gathering information on the wavefunctions.
!!  rhor(nfftf,nspden)=density in real space (smooth part if PAW).
!!  rhog(2,nfftf)=density in reciprocal space (smooth part if PAW).
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
!!
!! OUTPUT
!!  vxc_me    (b1gw:b2gw,b1gw:b2gw,nkibz,Dtset%nsppol)=matrix elements of $v_{xc}[nv+nc$.
!!  vxcval_me (b1gw:b2gw,b1gw:b2gw,nkibz,Dtset%nsppol)=matrix elements of $v_{xc}[nv]$.
!!  vhartr_me (b1gw:b2gw,b1gw:b2gw,nkibz,Dtset%nsppol)=matrix elements of $v_H$.
!!  vUpaw_me (b1gw:b2gw,b1gw:b2gw,nkibz,Dtset%nsppol)=matrix elements of  $v_U$.
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  All the quantities ($v_H$, $v_{xc}$ and $\psi$ are evaluated on the "fine" FFT mesh.
!!  In case of calculations with pseudopotials the usual mesh is defined by ecut. 
!!  For PAW calculations the dense FFT grid defined bt ecutdg is used
!!  Besides, in case of PAW, the matrix elements of V_hartree do not contain the onsite 
!!  contributions due to the coulombian potentials generate by ncore and tncore. 
!!  These quantities, as well as the onsite kinetic terms, are stored in Paw_ij%dij0.
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

subroutine calc_vHxc_braket(Dtset,nkibz,nkcalc,kcalc2ibz,b1,b2,nbnds,nsig_ab,npwvec,gsqcutf_eff,&
& nfftf_tot,nfftf,ngfftf,igfftf,Wf_info,vpsp,vhartr,vxc,Cprj,Pawtab,Paw_an,Pawang,&
& Pawfgrtab,Pawrad,Paw_ij,MPI_enreg,Cryst,tim_fourdp,rhor,rhog,usexcnhat,nhat,nhatgr,nhatgrdim,&
& vhartr_me,vxcval_me,vxc_me,vUpaw_me)
    
 use defs_basis
 use defs_datatypes
 use m_errors, only : assert
 use m_io_tools, only : flush_unit, get_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13nonlocal
 use interfaces_13paw
 use interfaces_13xc
 use interfaces_15gw, except_this_one => calc_vHxc_braket
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkcalc,nkibz,nhatgrdim,usexcnhat,nsig_ab
 integer,intent(in) :: b1,b2,nbnds,nfftf_tot,nfftf,tim_fourdp,npwvec
 real(dp),intent(in) :: gsqcutf_eff
 type(Dataset_type),intent(in) :: Dtset
 type(Wavefunctions_information),intent(inout) :: Wf_info
 type(MPI_type),intent(inout) :: MPI_enreg
 type(Pawang_type),intent(in) :: Pawang
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 integer,intent(in) :: ngfftf(18),igfftf(npwvec)
 integer,intent(in) :: kcalc2ibz(nkcalc)
 real(dp),intent(in) :: vhartr(nfftf),vxc(nfftf,Dtset%nspden),vpsp(nfftf)
 real(dp),intent(in) :: rhor(nfftf,Dtset%nspden),rhog(2,nfftf)
 real(dp),intent(in) :: nhat(nfftf,Dtset%nspden*Dtset%usepaw)
 real(dp),intent(in) :: nhatgr(nfftf,Dtset%nspden,3*nhatgrdim)
 complex(dpc),intent(out) :: vhartr_me(b1:b2,b1:b2,nkibz,Dtset%nsppol*nsig_ab)
 complex(dpc),intent(out) :: vxc_me   (b1:b2,b1:b2,nkibz,Dtset%nsppol*nsig_ab)
 complex(dpc),intent(out) :: vxcval_me(b1:b2,b1:b2,nkibz,Dtset%nsppol*nsig_ab)
 complex(dpc),intent(out) :: vUpaw_me (b1:b2,b1:b2,nkibz,Dtset%nsppol*nsig_ab)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Dtset%usepaw)
 type(Cprj_type),intent(in) ::  Cprj(Cryst%natom,Dtset%nspinor*nbnds*nkibz*Dtset%nsppol*Dtset%usepaw)
 type(Paw_an_type),intent(in) :: Paw_an(Cryst%natom)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat)
 type(Paw_ij_type),intent(in) :: Paw_ij(Cryst%natom)
 type(Pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom)

!Local variables-------------------------------
!scalars
 integer :: iat,ikc,ik_ibz,ib,jb,is,unt,shift,ispden,ispinor
 integer :: itypat,lmn_size,j0lmn,jlmn,ilmn,klmn,klmn1,klm,lmn2_size_max
 integer :: lmin,lmax,mm,isel,l_size,lm_size,lmn2_size,izero,im,cplex_dij
 integer :: ils,ilslm,ic,lm0,lpawu,nspinor,nsppol,nspden
 integer :: isp1,isp2,iab,nsploop,nkxc,option,n3xccc_,natom,ibsp1,ibsp2,iab_tr,first=0
 real(dp) :: nfftfm1,fact,DijH
 real(dp) :: enxc_val,vxcval_avg
 real(dp) :: vxc1,vxc1_val
 real(dp) :: re_p,im_p
 logical :: use_fineFFT,ldebug,ltest,use_PAWpU
 character(len=500) :: msg      
 character(len=fnlen) :: fnam
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer,parameter :: spinor_idxs(2,4)=RESHAPE((/1,1,2,2,1,2,2,1/),(/2,4/))
 integer,parameter :: trsp_idx(2:4)=(/2,4,3/)
 integer,allocatable :: dimlmn(:),indklmn_(:,:)
 real(dp) :: tmp_xc(2,nsig_ab),tmp_xcval(2,nsig_ab),tmp_H(2,nsig_ab),tmp_U(2,nsig_ab),dijU(2)
 real(dp) :: strsxc(6)
 real(dp) :: vxc1ab(2),vxc1ab_val(2)
 real(dp) :: rdum(1)
 real(dp),allocatable :: DijhatH(:,:),dijexxc(:,:,:)
 real(dp),allocatable :: prod(:)
 real(dp),allocatable :: kxc_(:,:),vh_(:),xccc3d_(:),vxc_val(:,:)
 complex(dpc) :: tmp(3)
 complex(gwpc),target,allocatable :: wfr1(:),wfr2(:)
 complex(dpc),allocatable :: vxcab(:),vxcab_val(:),u1cjg_u2dpc(:)
 complex(gwpc),pointer :: wfr1up(:),wfr1dwn(:)
 complex(gwpc),pointer :: wfr2up(:),wfr2dwn(:)
 type(Cprj_type),allocatable ::  Cprj_b1ks(:,:),Cprj_b2ks(:,:)
 
! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' calc_vHxc_braket : enter'
 call wrtout(std_out,msg,'COLL')
 call flush_unit(std_out)
#endif

 nspinor= Dtset%nspinor
 nsppol = Dtset%nsppol
 nspden = Dtset%nspden
 if (nspinor==2) WRITE(*,*)" Remember to ADD SO"
 !
 ! * Fake MPI_type for sequential part
 call initmpi_seq(MPI_enreg_seq) 
 MPI_enreg_seq%nproc_fft=1 
 MPI_enreg_seq%me_fft=0
 !
 ! === Evaluate $v_{xc}$ using only the valence charge ====
 write(msg,'(a)')&
& ' calc_vHxc_braket : calling rhohxc to calculate V_xc[n_val] (excluding non-linear core corrections) '
 call wrtout(std_out,msg,'COLL')

 option = 0 ! --> only exc, vxc, strsxc
 nkxc   = 0 ! --> no computation of XC kernel 
 n3xccc_= 0 ! --> no core
 izero  = Dtset%usepaw
 allocate(xccc3d_(n3xccc_),vh_(nfftf),kxc_(nfftf,nkxc),vxc_val(nfftf,nspden))

 call rhohxc(Dtset,enxc_val,gsqcutf_eff,izero,kxc_,MPI_enreg_seq,nfftf,ngfftf,&
& nhat,Dtset%usepaw,nhatgr,nhatgrdim,nkxc,nspden,n3xccc_,option,rhog,rhor,Cryst%rprimd,&
& strsxc,usexcnhat,vh_,vxc_val,vxcval_avg,xccc3d_)

 deallocate(vh_,xccc3d_,kxc_) 

 write(msg,'(a,f8.4,2a,f8.4,2a)')&
& ' E_xc[n_val]  = ',enxc_val,  ' [Ha]',&
& '<V_xc[n_val]> = ',vxcval_avg,' [Ha]',ch10
 call wrtout(std_out,msg,'COLL')

 if (nspinor==2) then 
  ! === Setup of the hermitian operator vxcab ===
  ! * if nsppol==4 vxc contains (v^11, v^22, Re[V^12], Im[V^12].
  ! * Cannot use directly Re and Im since we also need off-diagonal elements.
  allocate(vxcab(nfftf),vxcab_val(nfftf))
  vxcab    (:)=CMPLX(vxc    (:,3),vxc    (:,4))
  vxcab_val(:)=CMPLX(vxc_val(:,3),vxc_val(:,4))
 end if

 use_fineFFT=.FALSE.
 if (ANY(Wf_info%ngfft(1:3)/=ngfftf(1:3))) then  
  use_fineFFT=.TRUE.
  write(msg,'(a,3i4)')' calc_vHxc_braket : using fine FFT grid ',ngfftf(1:3)
  call wrtout(std_out,msg,'COLL')
 end if
 !
 ! ==== Get matrix elements of vxc[n], vxc_val[n_v], vH and vU ===
 ! * First take care of the plane wave part.
 nfftfm1=one/nfftf
 vxc_me   (:,:,:,:)=czero 
 vxcval_me(:,:,:,:)=czero
 vhartr_me(:,:,:,:)=czero
 vUpaw_me (:,:,:,:)=czero
 allocate(wfr1(nfftf*nspinor),wfr2(nfftf*nspinor))
 allocate(u1cjg_u2dpc(nfftf))

 !if (Dtset%gwmem==0.or.Dtset%gwmem==10) then 
 ! If the wavefunctions in real space are not stored,
 ! it is much faster to temporary store them in a small array.
 !
 !MG: Fabien, it was much faster because, in the old implementation,
 ! I was calling get_wfr twice inside the double loop. 
 ! The overhead of get_wfr was small if wfr were in memory 
 ! but, obviously, was quite large is wfr had to be recalculatead
 ! on the fly for each b1,b2 pair. I removed your changes; this new coding 
 ! should fix the problem as it is equivalent to your patch for gwmem,
 ! and requires less memory.
 !end if
 !
 ! === Loop over required k-points ===
 do ikc=1,nkcalc
  ik_ibz=kcalc2ibz(ikc)

  do is=1,nsppol
   do jb=b1,b2
    if (use_fineFFT) then 
     call fft_onewfn(Wf_info%paral_kgb,nspinor,Wf_info%npwwfn,nfftf,Wf_info%wfg(:,jb,ikc,is),wfr2,&
&     igfftf,ngfftf,tim_fourdp,MPI_enreg)
    else 
     call get_wfr(Wf_info,MPI_enreg,jb,ikc,is,wfr2)
    end if

    do ib=b1,jb ! Upper triangle 
     if (Dtset%gwcalctyp<10.and.ib/=jb) CYCLE
     if (use_fineFFT) then 
      call fft_onewfn(Wf_info%paral_kgb,nspinor,Wf_info%npwwfn,nfftf,Wf_info%wfg(:,ib,ikc,is),wfr1,&
&      igfftf,ngfftf,tim_fourdp,MPI_enreg)
     else 
      call get_wfr(Wf_info,MPI_enreg,ib,ikc,is,wfr1)
     end if

     u1cjg_u2dpc(1:nfftf)=CONJG(wfr1(1:nfftf))*wfr2(1:nfftf)

     vxc_me   (ib,jb,ik_ibz,is)=SUM(u1cjg_u2dpc(1:nfftf)*vxc    (1:nfftf,is))*nfftfm1
     vxcval_me(ib,jb,ik_ibz,is)=SUM(u1cjg_u2dpc(1:nfftf)*vxc_val(1:nfftf,is))*nfftfm1
     vhartr_me(ib,jb,ik_ibz,is)=SUM(u1cjg_u2dpc(1:nfftf)*vhartr (1:nfftf))   *nfftfm1

     if (nspinor==2) then 
      !Here I can skip 21 if ib==jb
      ! TODO rewrite everything in double precision complex should modify get_wfr.
      wfr1up  => wfr1(1:nfftf)
      wfr1dwn => wfr1(nfftf+1:2*nfftf)
      wfr2up  => wfr2(1:nfftf)
      wfr2dwn => wfr2(nfftf+1:2*nfftf)

      tmp(1) = SUM(CONJG(wfr1dwn)*      vxc(:,2) *wfr2dwn)*nfftfm1
      tmp(2) = SUM(CONJG(wfr1up )*      vxcab(:) *wfr2dwn)*nfftfm1
      tmp(3) = SUM(CONJG(wfr1dwn)*CONJG(vxcab(:))*wfr2up )*nfftfm1
      vxc_me(ib,jb,ik_ibz,2:4)=tmp(:)

      tmp(1) = SUM(CONJG(wfr1dwn)*      vxc_val(:,2) *wfr2dwn)*nfftfm1
      tmp(2) = SUM(CONJG(wfr1up )*      vxcab_val(:) *wfr2dwn)*nfftfm1
      tmp(3) = SUM(CONJG(wfr1dwn)*CONJG(vxcab_val(:))*wfr2up )*nfftfm1
      vxcval_me(ib,jb,ik_ibz,2:4)=tmp(:) 

      tmp(1) = SUM(CONJG(wfr1dwn)*vhartr(:)*wfr2dwn)*nfftfm1
      vhartr_me(ib,jb,ik_ibz,2  )=tmp(1) 
      vhartr_me(ib,jb,ik_ibz,3:4)=czero 
     end if
     !
     ! === Symmetrize matrix elements to get the lower triangle ===
     ! * In the collinear case, generate the lower triangle by just doing a complex conjugate.
     ! * In the noncollinear case do also a transposition since A_{12}^{ab} = A_{21}^{ba}^*
     !   2-->2, 3-->4, 4-->3
     if (ib/=jb) then 
      vxc_me   (jb,ib,ik_ibz,is)=CONJG(vxc_me   (ib,jb,ik_ibz,is))
      vxcval_me(jb,ib,ik_ibz,is)=CONJG(vxcval_me(ib,jb,ik_ibz,is))
      vhartr_me(jb,ib,ik_ibz,is)=CONJG(vhartr_me(ib,jb,ik_ibz,is))
      if (nspinor==2) then 
       do iab=2,4
        iab_tr=trsp_idx(iab)
        vxc_me   (jb,ib,ik_ibz,iab)=CONJG(vxc_me   (ib,jb,ik_ibz,iab_tr))
        vxcval_me(jb,ib,ik_ibz,iab)=CONJG(vxcval_me(ib,jb,ik_ibz,iab_tr))
        vhartr_me(jb,ib,ik_ibz,iab)=CONJG(vhartr_me(ib,jb,ik_ibz,iab_tr)) 
        ! off-diagonal are zero but oh well!
       end do
      end if
     else 
      ! For ib==jb, just impose Hermiticity 
      vxc_me   (jb,ib,ik_ibz,is)=half*(vxc_me   (ib,jb,ik_ibz,is)+CONJG(vxc_me   (ib,jb,ik_ibz,is)))
      vxcval_me(jb,ib,ik_ibz,is)=half*(vxcval_me(ib,jb,ik_ibz,is)+CONJG(vxcval_me(ib,jb,ik_ibz,is)))
      vhartr_me(jb,ib,ik_ibz,is)=half*(vhartr_me(ib,jb,ik_ibz,is)+CONJG(vhartr_me(ib,jb,ik_ibz,is)))
      if (nspinor==2) then 
       vxc_me   (jb,ib,ik_ibz,2)=half*(vxc_me   (ib,jb,ik_ibz,2)+CONJG(vxc_me   (ib,jb,ik_ibz,2)))
       vxcval_me(jb,ib,ik_ibz,2)=half*(vxcval_me(ib,jb,ik_ibz,2)+CONJG(vxcval_me(ib,jb,ik_ibz,2)))
       vhartr_me(jb,ib,ik_ibz,2)=half*(vhartr_me(ib,jb,ik_ibz,2)+CONJG(vhartr_me(ib,jb,ik_ibz,2)))
      end if
     end if

    end do !ib
   end do !jb

  end do !is
 end do !ikc

 if (nspinor==2) deallocate(vxcab,vxcab_val)
 deallocate(wfr1,wfr2,vxc_val)
 deallocate(u1cjg_u2dpc)
 !
 ! ====================================
 ! ===== Additional terms for PAW =====
 ! ====================================
 if (Dtset%usepaw==1) then 

  ! === Test if needed pointers in paw_ij are associated ===
  ltest=(associated(paw_ij(1)%dijxc).and.associated(paw_ij(1)%dijxc_val))
  call assert(ltest,'dijxc or dijxc_val not associated',__FILE__,__LINE__)

  if (ANY(pawtab(:)%usepawu>0)) then
   ltest=(Dtset%nspinor==1)
   call assert(ltest,'nspinor==2 not compatible with LDA+U',__FILE__,__LINE__)
   ltest=(associated(paw_ij(1)%dijU))
   call assert(ltest,'dijU not associated',__FILE__,__LINE__)
  end if

  if (Dtset%pawspnorb>0) then
   ltest=(associated(paw_ij(1)%dijso))
   call assert(ltest,'dijso not associated',__FILE__,__LINE__)
  end if

  ! TODO terminate the implementation of this routine.
  !call print_paw_ij(Paw_ij,Dtset%pawprtvol)

  natom=Cryst%natom
  nsploop=nspinor**2
  lmn2_size_max=MAXVAL(Pawtab(1:Cryst%ntypat)%lmn2_size) 

  ! === Build Hartree part of Dij_hat ===
  ! TODO this part should be calculate in pawdij
  ! Waiting for Marc"s patch.
  allocate(DijhatH(lmn2_size_max,natom)) ; DijhatH(:,:)=zero

  do iat=1,natom
   itypat=Cryst%typat(iat)
   l_size=Pawtab(itypat)%l_size
   lm_size=Paw_an(iat)%lm_size
   allocate(prod(l_size*l_size)) ; prod(:)=zero
   !
   ! Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
   if (Pawfgrtab(iat)%gylm_allocated==0) then
    if (associated(Pawfgrtab(iat)%gylm)) deallocate(Pawfgrtab(iat)%gylm)
    allocate(Pawfgrtab(iat)%gylm(Pawfgrtab(iat)%nfgd,lm_size))
    Pawfgrtab(iat)%gylm_allocated=2
    call pawgylm(Pawfgrtab(iat)%gylm,rdum,rdum,iat,Pawfgrtab(iat)%ifftsph,&
&    itypat,lm_size,Pawfgrtab(iat)%nfgd,1,0,0,Pawtab(itypat),&
&    Pawfgrtab(iat)%rfgd,Pawfgrtab(iat)%rfgd_allocated)
   end if
   !
   ! === prod_lm = $\int g_l Ylm v_h[tn+nhat] dr$ on the FFT mesh ===
   ! * Note that this quantity does not depend on ij
   ! TODO RECHECK THIS
   do ilslm=1,l_size**2
    do ic=1,Pawfgrtab(iat)%nfgd
     prod(ilslm)=prod(ilslm)+ &
&     (vhartr(Pawfgrtab(iat)%ifftsph(ic))+vpsp(Pawfgrtab(iat)%ifftsph(ic)))*Pawfgrtab(iat)%gylm(ic,ilslm)
     !prod(ilslm)=prod(ilslm)+&
     !& vhartr(Pawfgrtab(iat)%ifftsph(ic))                                  *Pawfgrtab(iat)%gylm(ic,ilslm)
    end do
   end do

   ! === Assembly Dij_hat ====
   lmn2_size=Paw_ij(iat)%lmn2_size
   allocate(indklmn_(6,lmn2_size))
   indklmn_(:,:)=Pawtab(itypat)%indklmn(:,:)
   do klmn=1,lmn2_size
    klm =indklmn_(1,klmn)
    lmin=indklmn_(3,klmn) 
    lmax=indklmn_(4,klmn)
    ! Get $\sum_lm q_ij^l prod_lm $
    do ils=lmin,lmax,2
     lm0=ils**2+ils+1
     do mm=-ils,ils
      ilslm=lm0+mm 
      isel=Pawang%gntselect(lm0+mm,klm)
      if (isel>0) DijhatH(klmn,iat)=DijhatH(klmn,iat)+prod(ilslm)*Pawtab(itypat)%qijl(ilslm,klmn)
     end do
    end do
   end do

   deallocate(prod,indklmn_)
  end do !iat
  ! * Normalization factor due to integration on the FFT mesh
  DijhatH(:,:)=DijhatH(:,:)*Cryst%ucvol/DBLE(nfftf_tot)

  use_PAWpU=ANY(pawtab(:)%usepawu>0) 
  !
  ! === Assemble PAW matrix elements ===
  allocate(dimlmn(natom)) 
  do iat=1,natom
   dimlmn(iat)=Pawtab(Cryst%typat(iat))%lmn_size
  end do
  allocate(Cprj_b1ks(natom,nspinor),Cprj_b2ks(natom,nspinor))
  call cprj_alloc(Cprj_b1ks,0,dimlmn)
  call cprj_alloc(Cprj_b2ks,0,dimlmn)

  ! === Loop over required k-points ===
  do ikc=1,nkcalc
   ik_ibz=kcalc2ibz(ikc)

   do is=1,nsppol
    shift=nspinor*nbnds*nkibz*(is-1)

    do jb=b1,b2
     ! === Load projected wavefunctions for this k-point, spin and band ===
     ! * Cprj are unsorted, full correspondence with xred. See ctocprj.F90!!
     ibsp2=nspinor*nbnds*(ik_ibz-1)+jb+shift-1
     do ispinor=1,nspinor
      ibsp2=ibsp2+1
      do iat=1,natom
       Cprj_b2ks(iat,ispinor)%cp(:,:)=Cprj(iat,ibsp2)%cp(:,:) 
      end do
     end do

     do ib=b1,jb ! Upper triangle 
      if (Dtset%gwcalctyp<10.and.ib/=jb) CYCLE
      ibsp1=nspinor*nbnds*(ik_ibz-1)+ib+shift-1
      do ispinor=1,nspinor
       ibsp1=ibsp1+1
       do iat=1,natom
        Cprj_b1ks(iat,ispinor)%cp(:,:)=Cprj(iat,ibsp1)%cp(:,:) 
       end do
      end do
      !
      ! === Get onsite matrix elements summing over atoms and channels ===
      ! * Spin is external and fixed (1,2) if collinear.
      ! * if noncollinear loop internally over the four components ab.
      tmp_xc   (:,:)=zero 
      tmp_xcval(:,:)=zero 
      tmp_H    (:,:)=zero 
      tmp_U    (:,:)=zero

      do iat=1,natom
       itypat=Cryst%typat(iat)
       lmn_size=Pawtab(itypat)%lmn_size
       cplex_dij=Paw_ij(iat)%cplex_dij
       klmn1=1

       do jlmn=1,lmn_size
        j0lmn=jlmn*(jlmn-1)/2
        do ilmn=1,jlmn
         klmn=j0lmn+ilmn
         ! TODO Be careful, here I assume that the onsite terms ij are symmetric 
         ! should check the spin-orbit case!
         fact=one ; if (ilmn==jlmn) fact=half

         ! === Loop over four components if nspinor==2 ===
         do iab=1,nsploop

          isp1=spinor_idxs(1,iab)
          isp2=spinor_idxs(2,iab)

          re_p=  Cprj_b1ks(iat,isp1)%cp(1,ilmn) * Cprj_b2ks(iat,isp2)%cp(1,jlmn) &
&               +Cprj_b1ks(iat,isp1)%cp(2,ilmn) * Cprj_b2ks(iat,isp2)%cp(2,jlmn) &
&               +Cprj_b1ks(iat,isp1)%cp(1,jlmn) * Cprj_b2ks(iat,isp2)%cp(1,ilmn) &
&               +Cprj_b1ks(iat,isp1)%cp(2,jlmn) * Cprj_b2ks(iat,isp2)%cp(2,ilmn) 

          im_p=  Cprj_b1ks(iat,isp1)%cp(1,ilmn) * Cprj_b2ks(iat,isp2)%cp(2,jlmn) &
&               -Cprj_b1ks(iat,isp1)%cp(2,ilmn) * Cprj_b2ks(iat,isp2)%cp(1,jlmn) &
&               +Cprj_b1ks(iat,isp1)%cp(1,jlmn) * Cprj_b2ks(iat,isp2)%cp(2,ilmn) &
&               -Cprj_b1ks(iat,isp1)%cp(2,jlmn) * Cprj_b2ks(iat,isp2)%cp(1,ilmn) 
         
          ! ==================================================
          ! === Load onsite matrix elements and accumulate ===
          ! ==================================================

          if (nspinor==1) then
           ! * Accumulate vxc[n1+nc] + vxc[n1+tn+nc].
           vxc1 = Paw_ij(iat)%dijxc(klmn,is)
           tmp_xc(1,iab)=tmp_xc(1,iab) + vxc1*re_p*fact
           tmp_xc(2,iab)=tmp_xc(2,iab) + vxc1*im_p*fact

           ! * Accumulate valence-only XC.
           vxc1_val=Paw_ij(iat)%dijxc_val(klmn,is)
           tmp_xcval(1,1)=tmp_xcval(1,1) + vxc1_val*re_p*fact
           tmp_xcval(2,1)=tmp_xcval(2,1) + vxc1_val*im_p*fact

           ! * Accumulate Hartree term of the PAW Hamiltonian ===
           ! FIXME Waiting for Marc"s patch : XG 080818 DONE
           !!DijH=Paw_ij(iat)%veij(klmn)+DijhatH(klmn,iat)
           DijH=Paw_ij(iat)%dijhartree(klmn)
           tmp_H(1,1)=tmp_H(1,1) + DijH*re_p*fact
           tmp_H(2,1)=tmp_H(2,1) + DijH*im_p*fact

           if (use_PAWpU) then 
            ! * Accumulate U term of the PAW Hamiltonian (only onsite AE contribution)
            dijU(1)=Paw_ij(iat)%dijU(klmn,is)
            tmp_U(1,1)=tmp_U(1,1) + dijU(1)*re_p*fact
            tmp_U(2,1)=tmp_U(2,1) + dijU(1)*im_p*fact
           end if

          else 
           ! === Spinorial case ===
           vxc1ab(1) = Paw_ij(iat)%dijxc(klmn1,  iab)
           vxc1ab(2) = Paw_ij(iat)%dijxc(klmn1+1,iab)
           tmp_xc(1,iab) = tmp_xc(1,iab) + (vxc1ab(1)*re_p - vxc1ab(2)*im_p)*fact
           tmp_xc(2,iab) = tmp_xc(2,iab) + (vxc1ab(2)*re_p + vxc1ab(1)*im_p)*fact
 
           vxc1ab_val(1) = Paw_ij(iat)%dijxc_val(klmn1,  iab)
           vxc1ab_val(2) = Paw_ij(iat)%dijxc_val(klmn1+1,iab)
           tmp_xcval(1,iab) = tmp_xcval(1,iab) + (vxc1ab_val(1)*re_p - vxc1ab_val(2)*im_p)*fact
           tmp_xcval(2,iab) = tmp_xcval(2,iab) + (vxc1ab_val(2)*re_p + vxc1ab_val(1)*im_p)*fact

           ! * Accumulate Hartree term of the PAW Hamiltonian ===
           ! * In GW, dijhartree is always real.
           ! FIXME Waiting for Marc"s patch : XG 080818 DONE
           if (iab==1.or.iab==2) then
            DijH = Paw_ij(iat)%dijhartree(klmn) 
            tmp_H(1,iab) = tmp_H(1,iab) + DijH*re_p*fact
            tmp_H(2,iab) = tmp_H(2,iab) + DijH*im_p*fact
           end if

           ! TODO "ADD LDA+U and SO"
           ! check this part

           if (use_PAWpU) then 
            ! * Accumulate U term of the PAW Hamiltonian (only onsite AE contribution)
            dijU(1)=Paw_ij(iat)%dijU(klmn1  ,iab)
            dijU(2)=Paw_ij(iat)%dijU(klmn1+1,iab)
            tmp_U(1,iab) = tmp_U(1,iab) + (dijU(1)*re_p - dijU(2)*im_p)*fact
            tmp_U(2,iab) = tmp_U(2,iab) + (dijU(2)*re_p + dijU(1)*im_p)*fact
           end if

          end if
         end do !iab
         klmn1=klmn1+cplex_dij

        end do !ilmn
       end do !jlmn
      end do !iat

      if (ib==jb) then 
       first=first+1
       if (nspinor==1) then
        if (first==1) write(std_out,'(a)')&
&        ' ik  b1  b2  is  onsite  vxc_tot  vxc_val  vhartr  vU '
        write(std_out,'(3i3,i2,1x,6f8.3)')&
&        ik_ibz,ib,jb,is,&
&        REAL(vxc_me   (ib,jb,ik_ibz,is))*Ha_eV,tmp_xc(1,1)*Ha_eV,tmp_xcval(1,1)*Ha_eV,&
&        REAL(vhartr_me(ib,jb,ik_ibz,is))*Ha_eV,tmp_H (1,1)*Ha_eV,REAL(vUpaw_me(ib,jb,ik_ibz,is))*Ha_eV
       else 
        if (first==1) write(std_out,'(a)')&
&        ' ik  b1  b2  3*vxc_tot(PW) 3*vxc_tot(sph) 3*vcx_val(PW)  3*vxc_val(sph) 2*vhartr(PW) '
        write(*,'(3i3,1x,14f8.3)')&
&        ik_ibz,ib,jb,&
&         REAL(vxc_me   (ib,jb,ik_ibz,1:3))*Ha_eV,tmp_xc   (1,1:3)*Ha_eV,&
&         REAL(vxcval_me(ib,jb,ik_ibz,1:3))*Ha_eV,tmp_xcval(1,1:3)*Ha_eV,&
&         REAL(vhartr_me(ib,jb,ik_ibz,1:2))*Ha_eV
          !,SUM(REAL(vUpaw_me(ib,jb,ik_ibz,:)))*Ha_eV
       end if
      end if
      !
      ! === Add to plane wave contribution ===
      if (nspinor==1) then
       vxc_me   (ib,jb,ik_ibz,is)=vxc_me   (ib,jb,ik_ibz,is)+CMPLX(tmp_xc(1,1),tmp_xc(2,1))
       vxcval_me(ib,jb,ik_ibz,is)=vxcval_me(ib,jb,ik_ibz,is)+CMPLX(tmp_xcval(1,1),tmp_xcval(2,1))
       vhartr_me(ib,jb,ik_ibz,is)=vhartr_me(ib,jb,ik_ibz,is)+CMPLX(tmp_H (1,1),tmp_H (2,1))
       vUpaw_me (ib,jb,ik_ibz,is)=CMPLX(tmp_U(1,1),tmp_U(2,1))
      else
       vxc_me   (ib,jb,ik_ibz,:)=vxc_me   (ib,jb,ik_ibz,:)+CMPLX(tmp_xc(1,:),tmp_xc(2,:))
       vxcval_me(ib,jb,ik_ibz,:)=vxcval_me(ib,jb,ik_ibz,:)+CMPLX(tmp_xcval(1,:),tmp_xcval(2,:))
       !vxcval_me(ib,jb,ik_ibz,:)=vxcval_me(ib,jb,ik_ibz,:)+CMPLX(tmp_xc(1,:),tmp_xc(2,:))
       vhartr_me(ib,jb,ik_ibz,:)=vhartr_me(ib,jb,ik_ibz,:)+CMPLX(tmp_H (1,:),tmp_H (2,:))
       vUpaw_me (ib,jb,ik_ibz,:)=CMPLX(tmp_U(1,:),tmp_U(2,:))
      end if
      !
      ! === Symmetrize matrix elements ===
      ! * In the collinear case, generate the lower triangle by just doing a complex conjugate.
      ! * In the noncollinear case do also a transposition since A_{12}^{ab} = A_{21}^{ba}^*
      if (ib/=jb) then 
       vxc_me   (jb,ib,ik_ibz,is)=CONJG(vxc_me   (ib,jb,ik_ibz,is))
       vxcval_me(jb,ib,ik_ibz,is)=CONJG(vxcval_me(ib,jb,ik_ibz,is)) 
       vhartr_me(jb,ib,ik_ibz,is)=CONJG(vhartr_me(ib,jb,ik_ibz,is))
       vUpaw_me (jb,ib,ik_ibz,is)=CONJG(vUpaw_me (ib,jb,ik_ibz,is))
       if (nspinor==2) then 
        do iab=2,4
         iab_tr=trsp_idx(iab)
         vxc_me   (jb,ib,ik_ibz,iab)=CONJG(vxc_me   (ib,jb,ik_ibz,iab_tr))
         vxcval_me(jb,ib,ik_ibz,iab)=CONJG(vxcval_me(ib,jb,ik_ibz,iab_tr))
         ! off-diagonal are zero but oh well!
         vhartr_me(jb,ib,ik_ibz,iab)=CONJG(vhartr_me(ib,jb,ik_ibz,iab_tr))
         vUpaw_me (jb,ib,ik_ibz,iab)=CONJG(vUpaw_me (ib,jb,ik_ibz,iab_tr))
        end do
       end if
      else 
       ! For ib==jb just impose Hermiticity 
       vxc_me   (jb,ib,ik_ibz,is)=half*(vxc_me   (jb,ib,ik_ibz,is)+CONJG(vxc_me   (jb,ib,ik_ibz,is)))
       vxcval_me(jb,ib,ik_ibz,is)=half*(vxcval_me(jb,ib,ik_ibz,is)+CONJG(vxcval_me(jb,ib,ik_ibz,is)))
       vhartr_me(jb,ib,ik_ibz,is)=half*(vhartr_me(jb,ib,ik_ibz,is)+CONJG(vhartr_me(jb,ib,ik_ibz,is)))
       vUpaw_me (jb,ib,ik_ibz,is)=half*(vUpaw_me (jb,ib,ik_ibz,is)+CONJG(vUpaw_me (jb,ib,ik_ibz,is)))
       if (nspinor==2) then 
        vxc_me   (jb,ib,ik_ibz,2)=half*(vxc_me   (ib,jb,ik_ibz,2)+CONJG(vxc_me   (ib,jb,ik_ibz,2)))
        vxcval_me(jb,ib,ik_ibz,2)=half*(vxcval_me(ib,jb,ik_ibz,2)+CONJG(vxcval_me(ib,jb,ik_ibz,2)))
        vhartr_me(jb,ib,ik_ibz,2)=half*(vhartr_me(ib,jb,ik_ibz,2)+CONJG(vhartr_me(ib,jb,ik_ibz,2)))
        vUpaw_me (jb,ib,ik_ibz,2)=half*(vUpaw_me (jb,ib,ik_ibz,2)+CONJG(vUpaw_me (jb,ib,ik_ibz,2)))
       end if
      end if

     end do !ib
    end do !jb

   end do !is
  end do !ikc

  deallocate(dimlmn,DijhatH)
  call cprj_free(Cprj_b1ks) ; deallocate(Cprj_b1ks)
  call cprj_free(Cprj_b2ks) ; deallocate(Cprj_b2ks)

 end if !PAW

 ldebug=.FALSE.
 if (ldebug) then
  fnam='__vHxcme.dat__' ; call isfile(fnam,'new')
  unt=get_unit() 
  open(unit=unt,file=fnam)

  do ikc=1,nkcalc
   ik_ibz=kcalc2ibz(ikc)
   do is=1,nsppol
    write(msg,'(a,i3,2a)')&
&    ' < nks| Vxc |nks> for k-point ',ikc,ch10,&
&    '  s  n   <Vxc>   <Vxc_val>   <VH>   [Ev] '
    call wrtout(unt,msg,'COLL')
    do ib=b1,b2
     if (nspinor==1) then
      write(msg,'(2i3,3(e14.6,2x))')&
&      is,ib,REAL(vxc_me(ib,ib,ik_ibz,is))*Ha_eV,&
&            REAL(vxcval_me(ib,ib,ik_ibz,is))*Ha_eV,&
&            REAL(vhartr_me(ib,ib,ik_ibz,is))*Ha_eV
     else 
      write(msg,'(2i3,3(e14.6,2x))')&
&      is,ib,SUM(REAL(vxc_me   (ib,ib,ik_ibz,:)))*Ha_eV,&
&            SUM(REAL(vxcval_me(ib,ib,ik_ibz,:)))*Ha_eV,&
&            SUM(REAL(vhartr_me(ib,ib,ik_ibz,:)))*Ha_eV
     end if
     call wrtout(unt,msg,'COLL')
    end do
   end do
  end do 

  close(unt)
 end if

#if defined DEBUG_MODE
 write(msg,'(a)')' calc_vHxc_braket : exit'
 call wrtout(std_out,msg,'COLL')
 call flush_unit(std_out)
#endif 

end subroutine calc_vHxc_braket
!!***
