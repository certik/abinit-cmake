!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_bands_sym_GW
!! NAME
!! get_bands_sym_GW
!!
!! FUNCTION
!!  This routine finds the irreducible representation associated to a set of degenerate bands at a given k-point 
!!  The irreducible representation is obtained by rotating a set of degenerate wave functions using only the symmetry 
!!  operations belonging to the little group of k. Two states are considered to be degenerate if their
!!  energy differs by less than EDIFF_TOL.
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  EDIFF_TOL(optional)= tolerance on the energy difference of two states (if not specified is set to 0.001 eV)
!!  ik_ibz=the index of the irreducible k-point to be analyzed 
!!  nsym=number of operations in space group
!!  nbnds=number of bands at this k-point (same for spin up and down)
!!  nsppol=number of independent spin polarizations
!!  only_trace=if .TRUE. only the trace of a single matrix per class is calculated (standard procedure if
!!   only the symmetry of bands is required). If .FALSE. all the matrices for each irreducible representation
!!   are calculated and stored in BSym
!!  Kmesh<BZ_mesh_type>=datatype gathering information on the k-mesh used 
!!  Wfs(Wavefunctions_information)= structure gathering information on wave functions
!!  mpi_enreg=MPI-parallelisation information (some already initialized,
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space (reduced coordinates)
!!  tnons(3,nsym)=fractional translations
!!  ene_k(nbnds,nsppol)=energies for this k-point
!!  usepaw=1 if PAW
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ntypat=number of type of atoms (onlu for PAW)
!!  irottb(Wfs%nfft,nsym) !this should be included in Wfs   FIXME rewrite everything
!!  dimlmn(natom)=number of lmn channel for each atom
!!  typat(natom)=type of each atom
!!  psps<pseudopotential_type>
!!    %indlmn(6,lmnmax,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn (for each atom type)
!!
!! OUTPUT
!!  BSym<Bands_Symmetries>=structure containing info on the little group of the k-point as well
!!   as the character of the representation associated to each set of degenerate states
!!  ierr=if differ from zero, the symmetry analysis cannot be performed, usually it means that 
!!   k is at zone border and there are non-symmorphic translations (see Notes)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  The method does not work if k is at zone border and the little group contains a non-symmorphic fractional 
!!  translation. The irreducible representation, D(R_t), multiplies wave functions as a row vector:
!!   R_t \psi_a = \sum_b M(R_t)_{ba} \psi_b
!! 
!!  As a net result, if R_t belong to the little group of k (Sk=k+G0), we have:
!!  M_ab(R_t) = e^{-i(k+G0)+\tau} \int e^{iG0.r} u_{ak}^*(r) u_{bk}(R^{-1}(r-\tau))\,dr 
!!
!!  For PAW there are two additional onsite terms involving <phi_i|phi_j(R^{-1}(r)> that can be 
!!  easily evaluated using the rotation matrix of real spherical harmonis, zarot(mp,m,l,R).
!!   Ylm(Rr)= \sum_mp zarot(mp,m,ll,R) Ylmp(r)
!!  Besides zarot(mp,m,l,R)=zarot(m,mp,l,R^{-1}}
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

subroutine get_Bands_Sym_GW(ntypat,natom,typat,nspinor,nsppol,nbnds,nkibz,Kmesh,Wfs,usepaw,pawtab,&
& pawang,psps,dimlmn,cprj_ibz,mpi_enreg,ik_ibz,ene_k,nsym,symrec,tnons,indsym,only_trace,irottb,BSym,ierr,&
& EDIFF_TOL) ! optional

 use defs_basis
 use defs_datatypes
#if defined DEBUG_MODE
 use m_IO_tools, only : flush_unit
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13nonlocal
 use interfaces_15gw, except_this_one => get_Bands_Sym_GW
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!this should be included in Wfs
!also this
!scalars
 integer,intent(in) :: ik_ibz,natom,nbnds,nkibz,nspinor,nsppol,nsym,ntypat
 integer,intent(in) :: usepaw
 integer,intent(out) :: ierr
 real(dp),intent(in),optional :: EDIFF_TOL
 logical,intent(in) :: only_trace
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Bands_Symmetries),intent(out) :: BSym
 type(MPI_type),intent(inout) :: mpi_enreg
 type(Pawang_type),intent(in) :: pawang
 type(Wavefunctions_information),intent(inout) :: Wfs
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: dimlmn(natom),indsym(4,nsym,natom),irottb(Wfs%nfft,nsym)
 integer,intent(in) :: symrec(3,3,nsym),typat(natom)
 real(dp),intent(in) :: ene_k(nbnds,nsppol),tnons(3,nsym)
 type(Cprj_type),intent(in) :: cprj_ibz(natom,nspinor*nbnds*nkibz*nsppol*usepaw)
 type(Pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: dim_cplx,iat,ib,ib1,ib2,ib_end,ib_start,icls,icplx,idx,idx_ks,il
 integer :: ilmn,ilpm,im,ir,isp,isym,isym_cls,itypat,ix,iy,iz,jb1,jb2,jl,jlmn
 integer :: jlpm,jm,k0lmn,klmn,lmax,nbnds_,nclass,ncplx,ncplx_MAX,ngfft1,ngfft2
 integer :: ngfft3,nlmn,nsym_cls,nsym_ltg,shift
 real(dp) :: EDIFF_TOL_,G0r,arg,dmimj,dmjmi,fact,fij,im_p,re_p,sij
 complex(gwpc) :: phmkG0t,trace
 character(len=500) :: msg
!arrays
 integer :: G0(3),R0(3)
 integer,allocatable :: Rm1rt(:)
 real(dp) :: kibz(3),kpG0(3),tmp(2)
 real(dp),allocatable :: DS_mmpl(:,:,:)
 complex(gwpc) :: dum(2)
 complex(gwpc),allocatable :: phase(:,:),ur1(:),ur2(:),ur2_rot(:)
 type(Degenerate_Bands),allocatable :: PAWCplx(:,:)
 type(cprj_type),allocatable :: cprj_b(:,:),cprj_brot(:,:)

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' get_Bands_symmetries : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 if (nspinor/=1) STOP ' get_Bands_symmetres : spinorial case not implemeented '

 EDIFF_TOL_=0.001/Ha_eV ; if (PRESENT(EDIFF_TOL)) EDIFF_TOL_=ABS(EDIFF_TOL)
 !
 ! ==== Initialize Bsym structure ====
 kibz(:)=Kmesh%ibz(:,ik_ibz) ; nbnds_=Wfs%my_maxb-Wfs%my_minb+1 
 if (nbnds_/=nbnds) STOP ' ERROR:  nbnds /= nbnds_' !this works only if bands are not spread
 call init_Bands_Symmetries(BSym,only_trace,nspinor,Wfs%nsppol,nbnds,kibz,nsym,symrec,tnons,EDIFF_TOL_,ene_k,ierr)
 !if (ierr/=0) RETURN !here I have to deal with the symmorphic case
 !
 ! ==== Precalculate phase on FFT mesh ====
 ngfft1=Wfs%ngfft(1) ; ngfft2=Wfs%ngfft(2) ; ngfft3=Wfs%ngfft(3)
 nsym_ltg=Bsym%nsym_sgk 
 allocate(phase(ngfft1*ngfft2*ngfft3,nsym_ltg)) ; phase(:,:)=cone

 do isym=1,nsym_ltg
  G0=BSym%G0(:,isym) ; if (ALL(G0==0)) CYCLE
  do iz=0,ngfft3-1
   do iy=0,ngfft2-1
    do ix=0,ngfft1-1
     G0r=two_pi*(  (G0(1)*ix)/DBLE(ngfft1) &
                  +(G0(2)*iy)/DBLE(ngfft2) &
                  +(G0(3)*iz)/DBLE(ngfft3) )
     ir=1+ix+iy*ngfft1+iz*ngfft1*ngfft2
     phase(ir,isym)=CMPLX(COS(G0r),SIN(G0r))
    end do
   end do
  end do
 end do
 !
 ! ==== Additional computation for PAW ====
 ! The onsite contribution reads:
 ! M_{12} = \sum_{aij} \delta_{l_i,l_j} <p_i|\psi2>^*<p_j|\psi1> S_ij D_^{l_i}{m_i,m_j}(R^{-1})
 ! where D is the rotation matrix of real spherical harmonics (zarot). 

! **** WARNING this part has to be tested ****

 if (usepaw==1) then 
  allocate(cprj_b   (natom,nbnds)) ; call cprj_alloc(cprj_b,   0,dimlmn)
  allocate(cprj_brot(natom,nbnds)) ; call cprj_alloc(cprj_brot,0,dimlmn)
  lmax=pawang%l_max-1 ; allocate(DS_mmpl(2*lmax+1,2*lmax+1,lmax+1))  !l_max is Max l+1
  ncplx_MAX=MAXVAL(Bsym%ncplx) ; allocate(PAWCplx(ncplx_MAX,nsppol)) 
  call nullify_Degenerate_Bands(PAWCplx)
  nclass=BSym%nclass ; nsym_ltg=Bsym%nsym_sgk 

  do isp=1,nsppol
   !
   ! === Retrieve matrix elements of PAW projectors === !TODO to be replaced by cprj_get
   ! TODO here I have to rotate the projector function otherwise results are wrong!!!!
   shift=nspinor*nbnds*nkibz*(isp-1)      
   idx_ks=nspinor*nbnds*(ik_ibz-1)+shift
   do ib=1,nbnds
    do iat=1,natom
     cprj_b(iat,ib)%cp(:,:)=cprj_ibz(iat,idx_ks+ib)%cp(:,:)
    end do
   end do
   !
   ! === Loop over set of degenerate bands for this spin ====
   ncplx=BSym%ncplx(isp)
   do icplx=1,ncplx 
    ib_start=BSym%Cplx(icplx,isp)%ib_start
    ib_end  =BSym%Cplx(icplx,isp)%ib_end
    dim_cplx=BSym%Cplx(icplx,isp)%dim_cplx
    allocate(PAWCplx(icplx,isp)%trace(nclass))

    !FIXME here there is a crash if only_trace is false
    if (.not.only_trace) allocate(PAWCplx(icplx,isp)%Rirr(dim_cplx,dim_cplx,nsym_ltg))
    !
    ! ==== Loop over classes ====
    idx=0
    do icls=1,nclass
     nsym_cls=BSym%nelements(icls)
     ! ==== Loop over elements in each class ====
     do isym_cls=1,nsym_cls
      !
      ! === If only the character is required, do this block only once ===
      idx=idx+1 ; if (only_trace.and.isym_cls/=1) CYCLE
      isym=BSym%sgk2symrec(idx)
      !arg=two_pi*DOT_PRODUCT(kibz(:,isym),R0)
      !phase(1)=COS(arg)
      !phase(2)=SIN(arg)
      !
      ! THIS PART should become a routine
      ! === DS_mmpl is the rotation matrix of real spherical harmonics associated to symrec(:,:,isym) ===
      ! DS_mmp multiply harmonics as row vectors, we need R^{-1} but we read R and invert m,mp in the equation below 
      DS_mmpl(:,:,:)=pawang%zarot(:,:,:,isym)
      call rotate_cprj(ntypat,natom,nbnds,psps,pawang,pawtab,typat,isym,cprj_b,cprj_brot)
      !
      ! === Loop over band indeces ===
      do ib1=ib_start,ib_end
       do ib2=ib_start,ib_end
        if (only_trace.and.ib1/=ib2) CYCLE
        !
        ! === Accumulating atom-centered contributions ===
        tmp(:)=zero
        do iat=1,natom
         nlmn=dimlmn(iat)
         itypat=typat(iat)
         do jlmn=1,nlmn 
          k0lmn=jlmn*(jlmn-1)/2 
          jl=psps%indlmn(1,jlmn,itypat)
          jm=psps%indlmn(2,jlmn,itypat)
          jlpm=1+jl+jm 
          do ilmn=1,jlmn  
           il=psps%indlmn(1,ilmn,itypat)
           if (il/=jl) CYCLE
           im=psps%indlmn(2,ilmn,itypat)
           ilpm=1+il+im 
           dmjmi=DS_mmpl(jlpm,ilpm,jl+1) ; dmimj=DS_mmpl(ilpm,jlpm,jl+1) 

           re_p=  dmjmi* ( cprj_b(iat,ib1)%cp(1,ilmn)*cprj_brot(iat,ib2)%cp(1,jlmn)  &
&                         +cprj_b(iat,ib1)%cp(2,ilmn)*cprj_brot(iat,ib2)%cp(2,jlmn) )&
&                +dmimj* ( cprj_b(iat,ib1)%cp(1,jlmn)*cprj_brot(iat,ib2)%cp(1,ilmn)  &
&                         +cprj_b(iat,ib1)%cp(2,jlmn)*cprj_brot(iat,ib2)%cp(2,ilmn) )

           im_p=  dmjmi* ( cprj_b(iat,ib1)%cp(1,ilmn)*cprj_brot(iat,ib2)%cp(2,jlmn)  &
&                         -cprj_b(iat,ib1)%cp(2,ilmn)*cprj_brot(iat,ib2)%cp(1,jlmn) )&
&                +dmimj* ( cprj_b(iat,ib1)%cp(1,jlmn)*cprj_brot(iat,ib2)%cp(2,ilmn)  &
&                         -cprj_b(iat,ib1)%cp(2,jlmn)*cprj_brot(iat,ib2)%cp(1,ilmn) )

           klmn=k0lmn+ilmn ; sij=pawtab(itypat)%sij(klmn) ; fij=one ; if (jlmn==ilmn) fij=half
           tmp(1)=tmp(1)+fij*sij*re_p ; tmp(2)=tmp(2)+fij*sij*im_p
          end do !ilmn
         end do !jlmn
        end do !iat
        !
        ! === Save values ===
        jb1=ib1-ib_start+1 ; jb2=ib2-ib_start+1 
        PAWCplx(icplx,isp)%Rirr(jb1,jb2,idx)=CMPLX(tmp(1),tmp(2))
       end do !ib2
      end do !ib1
      !
      ! === Calculate trace for each class ===
      if (isym_cls==1) then 
       trace=czero
       do jb1=1,dim_cplx
        trace=trace+PAWCplx(icplx,isp)%Rirr(jb1,jb1,idx)
       end do
       PAWCplx(icplx,isp)%trace(icls)=trace
      end if

     end do !isym_cls
    end do !icl
   end do !icplx
  end do !isp
  call cprj_free(cprj_b   ) ; deallocate(cprj_b   )
  call cprj_free(cprj_brot) ; deallocate(cprj_brot)
  deallocate(DS_mmpl)
 end if
 !
 ! === Evaluate matrix elements on FFT mesh and add PAW term if present ====
 fact=one/Wfs%nfftot ; nclass=BSym%nclass 
 allocate(ur1(Wfs%nfft),ur2(Wfs%nfft),ur2_rot(Wfs%nfft),Rm1rt(Wfs%nfft))
 do isp=1,nsppol
  ncplx=BSym%ncplx(isp)
  do icplx=1,ncplx 
   ! 
   ! ==== Loop over set of degenerate states ====
   ib_start=BSym%Cplx(icplx,isp)%ib_start
   ib_end  =BSym%Cplx(icplx,isp)%ib_end
   dim_cplx=BSym%Cplx(icplx,isp)%dim_cplx
   !
   ! ==== Loop over classes ====
   idx=0
   do icls=1,nclass
    nsym_cls=BSym%nelements(icls)
    ! ==== Loop over elements in each class ====
    do isym_cls=1,nsym_cls
     !
     ! === If only the character is required, do this block once ===
     idx=idx+1 ; if (BSym%only_trace.and.isym_cls/=1) CYCLE
     isym=BSym%sgk2symrec(idx)
     Rm1rt(:)=irottb(:,isym)
     G0(:)=BSym%G0(:,idx)
     kpG0(:)=two_pi*(kibz(:)+G0(:)) 
     arg=-two_pi*DOT_PRODUCT(kpG0,tnons(:,isym))
     phmkG0t=CMPLX(COS(arg),SIN(arg))
     !
     ! === Loop over matrix elements ===
     do ib1=ib_start,ib_end
      call get_wfr(Wfs,mpi_enreg,ib1,ik_ibz,isp,ur1)
      jb1=ib1-ib_start+1 
      do ib2=ib_start,ib_end
       if (BSym%only_trace.and.ib1/=ib2) CYCLE
       call get_wfr(Wfs,mpi_enreg,ib2,ik_ibz,isp,ur2)
       ! ==== Rotate the wave function and apply phase ====
       ! Note that the k-point is the same within a lattice vector
       do ir=1,Wfs%nfft
        ur2_rot(ir)=ur2(Rm1rt(ir))*phase(ir,idx)
       end do
       jb2=ib2-ib_start+1
       BSym%Cplx(icplx,isp)%Rirr(jb1,jb2,idx)=DOT_PRODUCT(ur1,ur2_rot)*fact*phmkG0t
       if (usepaw==1) then 
        BSym%Cplx(icplx,isp)%Rirr(jb1,jb2,idx)=BSym%Cplx(icplx,isp)%Rirr(jb1,jb2,idx)+PAWCplx(icplx,isp)%Rirr(jb1,jb2,idx)
       end if
      end do !ib2
     end do !ib1
     !
     ! === Calculate trace for each class ===
     if (isym_cls==1) then 
      trace=czero
      do jb1=1,dim_cplx
       trace=trace+BSym%Cplx(icplx,isp)%Rirr(jb1,jb1,idx)
      end do
      if (usepaw==1) trace=trace+PAWCplx(icplx,isp)%trace(icls)
      BSym%Cplx(icplx,isp)%trace(icls)=trace
     end if

    end do !isymc_class 
   end do !icls
  end do !icplx
 end do !isp

 deallocate(ur1,ur2,ur2_rot,Rm1rt,phase)
 if (usepaw==1) call destroy_Degenerate_Bands(PAWCplx)

#if defined DEBUG_MODE
 write(msg,'(a)')' get_Bands_symmetres : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine get_Bands_Sym_GW
!!***

!!****f* ABINIT/rotate_cprj
!! NAME
!! rotate_cprj
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
!! ntypat=number of types of atom
!! natom=number of atoms
!! nbnds=number of bands for this k-point ans spin
!! pawang <type(pawang_type)>=paw angular mesh and related data
!! pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!! typat(natom)=type of eahc atom
!! isym=index of the symmetry in the symrec arrays that preserves the given k-point within a reciprocal lattice vector
!! cprj_in(natom,nbnds) <type(cprj_type)>= projected input wave functions <Proj_i|Cnk> with all NL projectors at fixed k-point 
!!
!! OUTPUT
!! cprj_out(natom,nbnds) <type(cprj_type)>= projection of the smooth PAW wave function onto 
!!  projectors centered on equivalent sites of the crystal (non restricted to be in the firs unit cell)
!!  The equivalent site is defined according to the symmetry operation isym. Thus cprj_out contains
!!
!!  cprj_out(a,b)=<p_j^{R^{-1}(L_a-\tau)} | \tpsi_nk> if  R is the isym operation  with fractional translation \tau
!!  L_a is the position of the initial atom inside the first unit cell
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

subroutine rotate_cprj(ntypat,natom,nbnds,psps,pawang,pawtab,typat,isym,cprj_in,cprj_out)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isym,natom,nbnds,ntypat
 type(Pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: typat(natom)
 type(Cprj_type),intent(in) :: cprj_in(natom,nbnds)
 type(Cprj_type),intent(out) :: cprj_out(natom,nbnds)
 type(Pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iat,ib,iband,indexj,itypat,jl,jl0,jlmn,jln,jln0,jlpm,jm,jn,ll,lmax
 integer :: mmp,nlmn
!arrays
 real(dp) :: tmp(2)
 real(dp),allocatable :: DRm1_mmpl(:,:,:)

! *************************************************************************

 lmax=pawang%l_max-1 ; allocate(DRm1_mmpl(2*lmax+1,2*lmax+1,lmax+1))  !l_max is Max l+1
 DRm1_mmpl(:,:,:)=pawang%zarot(:,:,:,isym)

 do iat=1,natom
  itypat=typat(iat)
  nlmn=pawtab(itypat)%lmn_size
  jl0=-1 ; jln0=-1 ; indexj=1
  do jlmn=1,nlmn
   jl  =psps%indlmn(1,jlmn,itypat)
   jm  =psps%indlmn(2,jlmn,itypat)
   jn  =psps%indlmn(3,jlmn,itypat)
   jln =psps%indlmn(5,jlmn,itypat) 
   jlpm=1+jl+jm 
   if (jln/=jln0) indexj=indexj+2*jl0+1

   do iband=1,nbnds
    ! === For each band, calculate contribution due to rotated real spherical harmonics ===
    ! === we have to apply zarot(R^{-1}) so we use the transpose in the m,mp subspace
    !<p_j^{R^{-1}(L_a-\tau)}|\tpsi_nk> = \sum_x D_{x,m_j}^l(R^-1) <p_nlx,\tpsi_nk> 
    tmp(:)=zero
    do mmp=1,2*jl+1
     tmp(1)=tmp(1)+DRm1_mmpl(jlpm,mmp,jl+1)*cprj_in(iat,iband)%cp(1,indexj+mmp)
     tmp(2)=tmp(2)+DRm1_mmpl(jlpm,mmp,jl+1)*cprj_in(iat,iband)%cp(2,indexj+mmp)
    end do
    !
    ! === Save values ===
    cprj_out(iat,iband)%cp(1,jlmn)=tmp(1)
    cprj_out(iat,iband)%cp(2,jlmn)=tmp(2)
   end do !iband
   jl0=jl ; jln0=jln
  end do !jlmn
 end do !iat
 deallocate(DRm1_mmpl)

end subroutine rotate_cprj
!!***
