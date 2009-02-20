!{\src2tex{textfont=tt}}
!!****f* ABINIT/eig2tot
!! NAME
!! eig2tot
!!
!! FUNCTION
!! This routine calculates the second-order eigenvalues.
!! The output eig2_tot is this quantity for the input k point.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (PB, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  eig_k(mband*nsppol)= 0-order eigenvalues at all K-points: <k,n'|H(0)|k,n'> (hartree)
!!  eig1_k(2*nsppol*mband**2)=matrix of first-order: <k+Q,n'|H(1)|k,n> (hartree) (calculated in cgwf3)
!!  
!!  mband=maximum number of bands
!!  nband_k=number of bands for this k point
!!  nband_kq=number of bands at k+Q point
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!
!! OUTPUT
!!  eig2_diakq(nband_k)=diagonal part of the second-order eigenvalues: E^{(2),diag}_{k,q,j}
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!      dotprod_g
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine eig2tot(clflg,cg1_pert,gh1_pert,eigen0,eigenq,eigen1,eig2nkq,&
&  indsym,istwfk_pert,mband,mk1mem,natom,mpert,nsym,mpi_enreg,mpw1,nkpt_rbz,&
&  nspinor,nsppol,occ_k,qpt,sciss,symq,symrec,symrel,timrev,tkq)

 use defs_basis
 use defs_datatypes
#if defined HAVE_NETCDF
 use netcdf
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12spacepar
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mk1mem,mpert,mpw1,natom,nkpt_rbz,nspinor,nsppol
 integer,intent(in) :: nsym,timrev
 real(dp),intent(in) :: sciss
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: clflg(3,mpert),indsym(4,nsym,natom)
 integer,intent(in) :: istwfk_pert(nkpt_rbz,3,mpert),symq(4,2,nsym)
 integer,intent(in) :: symrec(3,3,nsym),symrel(3,3,nsym),tkq(nkpt_rbz)
 real(dp),intent(in) :: cg1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol,3,mpert)
 real(dp),intent(in) :: eigen0(nkpt_rbz*mband*nsppol)
 real(dp),intent(in) :: eigen1(nkpt_rbz*2*nsppol*mband**2,3,mpert)
 real(dp),intent(in) :: eigenq(nkpt_rbz*mband*nsppol)
 real(dp),intent(in) :: gh1_pert(nkpt_rbz,mband,3,mpert,2,mpw1*nspinor)
 real(dp),intent(in) :: occ_k(mband*nkpt_rbz*nsppol),qpt(3)
 real(dp),intent(out) :: eig2nkq(2,mband*nsppol,nkpt_rbz,3,mpert,3,mpert)

!Local variables-------------------------------
!tolerance for non degenerated levels
!scalars
 integer :: band2tot_index,band_index,bandtot_index,iband,icg1,icg2,idir1,idir2
 integer :: ikpt,ipert1,ipert2,isppol,istwf_k,jband
 real(dp),parameter :: etol=1.0d-3
 real(dp) :: dot2i,dot2r,doti,dotr,eig1_i1,eig1_i2,eig1_r1,eig1_r2,eig2_diai
 real(dp) :: eig2_diar
!arrays
 integer :: blk1flg(3,mpert,3,mpert)
 real(dp) :: cwavef(2,mpw1*nspinor),cwavef2(2,mpw1*nspinor)
 real(dp) :: eig2nkq_tmp(2,3,mpert,3,mpert),gh(2,mpw1*nspinor)
 real(dp) :: gh1(2,mpw1*nspinor)

! *********************************************************************

 band2tot_index =0
 bandtot_index=0
 band_index=0
 icg1=0
 eig2nkq(:,:,:,:,:,:,:) = zero
 blk1flg(:,:,:,:) = 0


 do isppol=1,nsppol
  do ikpt =1,nkpt_rbz
   do ipert1=1,mpert
    do idir1=1,3
     if(clflg(idir1,ipert1)==0)cycle
     istwf_k = istwfk_pert(ikpt,idir1,ipert1)
     do ipert2=1,mpert
      do idir2=1,3
       if(clflg(idir2,ipert2)==0)cycle
       blk1flg(idir1,ipert1,idir2,ipert2)=1
       do iband=1,mband
        eig2_diar = zero
        eig2_diai = zero
        icg2 = (tkq(ikpt) -1)*mpw1*nspinor*mband !does not work with isppol
        cwavef(:,:) = cg1_pert(:,1+(iband-1)*mpw1*nspinor+icg2:iband*mpw1*nspinor+icg2,idir2,ipert2)
        gh1(:,:)    = gh1_pert(ikpt,iband,idir1,ipert1,:,:)
        cwavef2(:,:)= cg1_pert(:,1+(iband-1)*mpw1*nspinor+icg2:iband*mpw1*nspinor+icg2,idir1,ipert1)
        gh(:,:)     = gh1_pert(ikpt,iband,idir2,ipert2,:,:)
!       if(ipert1==ipert2)then
        call dotprod_g(dotr,doti,istwf_k,mpi_enreg,mpw1*nspinor,1,cwavef,gh1)
        call dotprod_g(dot2r,dot2i,istwf_k,mpi_enreg,mpw1*nspinor,1,gh,cwavef2)
        doti = zero
        dot2i = zero
!       else
!       call dotprod_g(dotr,doti,istwf_k,mpi_enreg,mpw1*nspinor,2,cwavef,gh1)
!       call dotprod_g(dot2r,dot2i,istwf_k,mpi_enreg,mpw1*nspinor,2,gh,cwavef2)
!       end if

        do jband=1,mband
         if((abs(eigenq(jband+bandtot_index)-eigen0(iband+bandtot_index))>etol).and.(abs(occ_k(jband+bandtot_index))>tol8)) then  

          eig1_r1 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir1,ipert1)  
          eig1_r2 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir2,ipert2)
          eig1_i1 = eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
          eig1_i2 = - eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir2,ipert2) !the negative sign is from the CC

          if(abs(occ_k(iband+bandtot_index))>tol8) then
           eig2_diar = eig2_diar + (eig1_r1*eig1_r2 - eig1_i1*eig1_i2)/(eigenq(jband+bandtot_index) - eigen0(iband+bandtot_index)) 
           eig2_diai = eig2_diai + (eig1_r1*eig1_i2 + eig1_i1*eig1_r2)/(eigenq(jband+bandtot_index) - eigen0(iband+bandtot_index))
          else
           eig2_diar = eig2_diar + (eig1_r1*eig1_r2 - eig1_i1*eig1_i2)/(eigenq(jband+bandtot_index) -&
&           eigen0(iband+bandtot_index)-sciss)
           eig2_diai = eig2_diai + (eig1_r1*eig1_i2 + eig1_i1*eig1_r2)/(eigenq(jband+bandtot_index) -&
&           eigen0(iband+bandtot_index)-sciss)
          end if !for sciss  NOTE ONE SHOULD USE A BETTER CONDITION
         end if ! on degenerate bands
        end do !jband
        eig2nkq(1,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = half*(dotr + dot2r) - eig2_diar 
        eig2nkq(2,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = - eig2_diai 
       end do !iband
      end do !idir2
     end do !ipert2
    end do  !idir1
   end do   !ipert1
   band2tot_index = band2tot_index + 2*mband**2
   bandtot_index = bandtot_index + mband
   icg1 = icg1 + mpw1*nspinor*mband
  end do    !ikpt
  band_index = band_index + mband
 end do !isppol
 write(*,*) 'eig2tot: sortie',nsppol

!band_index=0
!do isppol=1,nsppol
!do ikpt=1,nkpt_rbz
!do iband=1,mband
!eig2nkq_tmp(:,:,:,:,:) = eig2nkq(:,iband+band_index,ikpt,:,:,:,:)
!call d2sym3(blk1flg,eig2nkq_tmp,indsym,mpert,natom,nsym,qpt,symq,symrec,symrel,timrev)
!eig2nkq(1,iband+band_index,ikpt,:,:,:,:) = eig2nkq_tmp(1,:,:,:,:)
!eig2nkq(2,iband+band_index,ikpt,:,:,:,:) = eig2nkq_tmp(2,:,:,:,:)
!end do
!end do
!band_index = band_index + mband
!end do
end subroutine eig2tot
!!***

