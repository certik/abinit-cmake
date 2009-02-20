!{\src2tex{textfont=tt}}
!!****f* ABINIT/joint_dos
!! NAME
!! joint_dos
!!
!! FUNCTION
!!  Calculate the joint density of states 
!!   $J(\omega,q)=\sum_{k,b1,b2} \delta(\omega-\epsilon_{k-q,b2}+\epsilon_{k,b1})$
!!  for a given set of external q-points. 
!!  Two different quadrature methods are employed according to the input variable method:
!!   method=1 => simple gaussian broadenig 
!!   method=2 => tetrahedron method
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nqpt=number of q-points for JDOS
!!  Cryst<Crystal_structure>=Info on unit cell and its symmetries
!!    %nsym=number of symmetry operations
!!    %symrec(3,3,nsym)=symmetry operations in reciprocal space (reduced coordinates)
!!    %tnons(3,nsym)=fractional translations
!!    %timrev=2 if time-reversal holds, 1 otherwise
!!  qpt(3,nqpt)=q-points in reduced coordinates
!!  mband=max number of bands over k-points and spin
!!  nkibz=number of irreducible k-points
!!  kibz(3,nkibz)=the irreducible k-points (reduced coordinates)
!!  nsppol=2 for spin-polarized, 1 otherwise
!!  eigen(mband*nkibz*nsppol)=energies
!!  occ(mband*nkibz*nsppol)=occupation numbers
!!  method=1 for gaussian broadening, 2 for tetrahedron 
!!  step=freqency step for JDOS
!!  broad=only for gaussian method, the broadening in Ha 
!!  nband(nkpt*nsppol)=number of bands for each k-point and spin
!!  
!! OUTPUT
!!  Only write
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      assert,destroy_bz_mesh_type,destroy_little_group,flush_unit,get_bz_diff
!!      get_bz_item,leave_new,nullify_little_group,setup_kmesh,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine joint_dos(nqpt,qpt,mband,nkibz,kibz,nsppol,nband,Cryst,eigen,occ,method,step,broad)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit, get_unit
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => joint_dos
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,method,nkibz,nqpt,nsppol
 real(dp),intent(in) :: broad,step
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 integer,intent(in) :: nband(nkibz*nsppol)
 real(dp),intent(in) :: eigen(mband*nkibz*nsppol),kibz(3,nkibz)
 real(dp),intent(in) :: occ(mband*nkibz*nsppol),qpt(3,nqpt)

!Local variables-------------------------------
!scalars
 integer :: ib1,ib2,idx_kb1,idx_kmqb2,iend,ii,ik_bz,ik_ibz,ikmq_bz,ikmq_ibz,io
 integer :: iq,isppol,istart,isym_k,isym_kmq,itim_k,itim_kmq,nband_k,nband_kmq
 integer :: nfound,nkbz,nomega,npwe,npwvec,prtvol,unt,use_umklp
 real(dp),parameter :: TOL_OCC=tol6
 real(dp) :: dossmear,ene1,ene2,gaussfactor,gaussprefactor,gaussval,max_tr,occ1
 real(dp) :: occ2,trans,xx
 logical :: ltest
 character(len=500) :: msg
 character(len=fnlen) :: fnam
 type(BZ_mesh_type) :: Kmesh
 type(Little_group) :: Ltg_q
!arrays
 integer :: G0(3)
 integer,allocatable :: idx(:,:)
 real(dp) :: gmet_dummy(3,3),kbz(3),kmq(3),qq(3)
 real(dp),allocatable :: jdos(:,:,:),omega(:),wtk(:)

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' joint_dos : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif
 !
 ! === Check input ===
 ltest=(Cryst%timrev==1.or.Cryst%timrev==2)
 call assert(ltest,'timrev should be 1 or 2,',__FILE__,__LINE__)
 ltest=(method==1.or.method==2)
 call assert(ltest,'method should be 1 or 2,',__FILE__,__LINE__)
 !
 ! === Starting index for each k-point and spin ===
 allocate(idx(nkibz,nsppol)) ; ii=0
 do isppol=1,nsppol
  do ik_ibz=1,nkibz
   idx(ik_ibz,isppol)=ii
   nband_k=nband(ik_ibz+(isppol-1)*nkibz)
   ii=ii+nband_k
  end do
 end do
 !
 ! === Initialize the k-mesh structure ===
 prtvol=0 
 call setup_kmesh(nkibz,kibz,Cryst,Kmesh,prtvol)
 nkbz=Kmesh%nbz
 allocate(wtk(nkibz)) ; wtk=Kmesh%wt/SUM(Kmesh%wt)

 max_tr=MAXVAL(eigen)-MINVAL(eigen) ; if (method==1) max_tr=max_tr+five*broad+tol6
 nomega=NINT(max_tr/step)+1
 allocate(omega(nomega))
 do io=1,nomega
  omega(io)=step*(io-1)
 end do

 allocate(jdos(nomega,nsppol,nqpt)) ; jdos=zero
 SELECT CASE (method) 
 CASE (1)
  gaussprefactor=one/(dossmear*SQRT(two_pi))
  gaussfactor=one/(sqrt2*dossmear)
  write(msg,'(4a,f8.5,2a,f8.5)')ch10,&
&  ' joint_dos : calculating Joint DOS using gaussian method :',ch10,&
&  ' gaussian smearing [eV] = ',broad*Ha_eV,ch10,&
&  ' energy step       [eV] = ',step*Ha_eV
 CASE (2)
  write(msg,'(2a)')ch10,&
&  ' mkphdos : calculating phonon DOS using tetrahedron method '
  stop "not implemented yet"
 CASE DEFAULT 
  STOP 'Wrong value for method'
 END SELECT
 call wrtout(std_out,msg,'COLL')

 call nullify_little_group(Ltg_q)

 do iq=1,nqpt
  qq(:)=qpt(:,iq) 
  !use_umklp=1 ; npwvec=0 ; npwe=0
  !call setup_little_group(qq,Kmesh,gmet,Cryst,npwvec,gvec,npwe,use_umklp,prtvol,Ltg_q)

  do isppol=1,nsppol
   do ik_bz=1,nkbz
    !
    ! === Get k and k-q in BZ as well as their image in IBZ ===
    call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k)
    call get_BZ_diff(Kmesh,kbz,qq,ikmq_bz,G0,nfound)
    if (nfound/=1) call leave_new('COLL') 
    call get_BZ_item(Kmesh,ikmq_bz,kmq,ikmq_ibz,isym_kmq,itim_kmq)

    nband_k  =nband(ik_ibz  +(isppol-1)*nkibz)
    nband_kmq=nband(ikmq_ibz+(isppol-1)*nkibz)
    idx_kmqb2=idx(ikmq_ibz,isppol)
    idx_kb1  =idx(  ik_ibz,isppol)

    do ib2=1,nband_kmq
     occ2=occ  (idx_kmqb2+ib2)
     ene2=eigen(idx_kmqb2+ib2)
     do ib1=1,nband_k
      occ1=occ  (idx_kb1+ib1)
      if (ABS(occ1*(one-occ2))<TOL_OCC) CYCLE
      ene1=eigen(idx_kb1+ib1)
      trans=ene2-ene1
      ! three*broad should be enough, five should give smooth curves for coarse k-meshes
      iend  =NINT((trans+five*broad)/step) ; if (iend>nomega) iend=nomega
      istart=NINT((trans-five*broad)/step) ; if (istart<=0)  istart=1
      !
      ! === Accumulate ===
      do io=istart,iend
       xx=(omega(io)-trans)*gaussfactor
       gaussval=gaussprefactor*EXP(-xx*xx)
       jdos(io,isppol,iq)=jdos(io,isppol,iq)+gaussval!+wtk(ik_bz)*gaussval
      end do

     end do !ib1
    end do !ib2

   end do !ikbz
  end do !isppol
 end do !iqpt
 !
 ! === Write results ===
 fnam='JDOS' ; call isfile(fnam,'new')
 unt=get_unit()
 open(file=fnam,unit=unt,form='formatted')
 write(unt,'(a)')         '# Joint density of states. All in eV units '
 write(unt,'(a,es16.8,a)')'# Frequency step ',step*Ha_eV,' [eV]'
 write(unt,'(a,es16.8,a)')'# Max Frequency  ',max_tr*Ha_eV,' [eV]'
 write(unt,'(a,i2)')      '# Number of Independendent polarization',nsppol
 write(unt,'(a,i5)')      '# Number of k-points in BZ ',Kmesh%nbz
 write(unt,'(a,i5)')      '# Number of analysed q-points',nqpt
 do iq=1,nqpt
  write(unt,'(a,i4,a,3es16.8)')'# ',iq,') ',qpt(:,iq)
 end do
 write(unt,'(a)')'#'
 write(unt,'(a)')'# omega(io),((jdos(io,isppol,iq),isppol=1,nsppol),iq=1,nqpt)'
 write(unt,'(a)')'#'
 do io=1,nomega
  write(unt,'(es16.8)')omega(io),((jdos(io,isppol,iq),isppol=1,nsppol),iq=1,nqpt)
 end do
 close(unt)
 !
 ! === Free memory ===
 deallocate(idx,omega,jdos,wtk)
 call destroy_BZ_mesh_type(Kmesh)
 call destroy_Little_group(Ltg_q)

#if defined DEBUG_MODE
 write(msg,'(a)')' joint_dos : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine joint_dos
!!***
