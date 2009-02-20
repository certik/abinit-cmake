!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawgrnl
!!
!! NAME
!! pawgrnl
!!
!! FUNCTION
!! PAW: Add to GRadients of total energy due to non-local term of Hamiltonian
!!      the contribution due to Dij derivatives
!! In particular, compute contribution to forces, stresses, dyn. matrix
!! Remember: Vnl=Sum_ij[|p_i>Dij<p_j|]
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dimnhat=second dimension of array nhat (0 or # of spin components)
!!  dimvtrial=second dimension of array vtrial (1 or # of spin components)
!!  natom=number of atoms in cell
!!  nattyp(ntypat)=array describing how many atoms of each type in cell
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nhat(nfft,dimnhat)=compensation charge density on rectangular grid in real space
!!  nspden=number of spin-density components
!!  nsym=number of symmetries in space group
!!  ntypat=number of types of atoms
!!  optgr= 1 if gradients with respect to atomic position(s) have to be computed
!!  optgr2= 1 if 2nd gradients with respect to atomic position(s) have to be computed
!!  optstr= 1 if gradients with respect to strain(s) have to be computed
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  typat(natom)=type integer for each atom in cell
!!  vtrial(nfft,dimvtrial)= total potential
!!
!! SIDE EFFECTS
!!  At input, this terms contain contribution from non-local projectors derivatives
!!  At output, they are updated with the contribution of Dij derivatives
!!  ==== if optgr=1 ====
!!   grnl(3*natom) =gradients of NL energy wrt atomic coordinates
!!  ==== if optstr=1 ====
!!   nlstr(6) =gradients of NL energy wrt strains
!!  ==== if optgr2=1 ====
!!   dyfrnl(3,3,natom) =2nd gradients of NL energy wrt atomic coordinates
!!
!! PARENTS
!!      dyfnl3,etotfor,forstr
!!
!! CHILDREN
!!      metric,strconv,xcomm_init,xsum_mpi
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pawgrnl(atindx1,dimnhat,dimvtrial,dyfrnl,grnl,mpi_enreg,natom,nattyp,nfft,ngfft,&
&                  nhat,nlstr,nspden,nsym,ntypat,optgr,optgr2,optstr,pawang,pawfgrtab,&
&                  pawrhoij,pawtab,rprimd,symrec,typat,vtrial)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_13paw, except_this_one => pawgrnl
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dimnhat,dimvtrial,natom,nfft,nspden,nsym,ntypat,optgr
 integer,intent(in) :: optgr2,optstr
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: atindx1(natom),nattyp(ntypat),ngfft(18),symrec(3,3,nsym)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: nhat(nfft,dimnhat),rprimd(3,3),vtrial(nfft,dimvtrial)
 real(dp),intent(inout) :: dyfrnl(3,3,natom*optgr2),grnl(3*natom*optgr)
 real(dp),intent(inout) :: nlstr(6*optstr)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iatm,iatom,iatshft,ic,ier,ils,ilslm,irhoij,isel,ishift_gr
 integer :: ishift_gr2,ishift_str,ispden,ispvtr,itypat,jc,klm,klmn,lm_size
 integer :: lm0,lmax,lmin,lmn2_size,mm,mu,mua,mub,mushift,nfftot,nfgd,ngrad,ngradp
 integer :: nsploop,old_paral_level,opt1,opt2,opt3,spaceComm
 real(dp) :: dlt_tmp,dro,dytmp,fact,grhat_x,hatstr_diag,r2,ro,ro_d,ucvol
 character(len=500) :: message
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 real(dp) :: gmet(3,3),gprimd(3,3),hatstr(6),rdum(1),rmet(3,3),work1(3,3)
 real(dp) :: work2(3,3)
 real(dp),allocatable :: grhat_tmp(:),prod(:,:),prodp(:,:),vloc(:)

! *************************************************************************

!DEBUG
!write(6,*)' pawgrnl: enter '
!stop
!ENDDEBUG

!Compatibility test
 if (optgr2==1.and.pawrhoij(1)%ngrhoij==0) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pawgrnl : BUG -',ch10,&
&  '  Inconsistency between variables optgr2 and ngrhoij !'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Initializations and allocations
 ngrad=0;ngradp=0
 ishift_gr=0;ishift_gr2=0;ishift_str=0
 if (optgr==1) then
  ngrad=ngrad+3
  ishift_gr2=ishift_gr2+3
 end if
 if (optgr2==1) then
  ngrad=ngrad+6;ngradp=ngradp+3
 end if
 if (optstr==1) then
  hatstr=zero
  ngrad=ngrad+6
  ishift_gr=ishift_gr+6
  ishift_gr2=ishift_gr2+6
 end if

 allocate(grhat_tmp(ngrad))
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 fact=ucvol/dble(nfftot)
 nsploop=nspden;if (dimvtrial<nspden) nsploop=2

!Loops over types and atoms
 iatshft=0
 do itypat=1,ntypat

  lmn2_size=pawtab(itypat)%lmn2_size
  do iatm=iatshft+1,iatshft+nattyp(itypat)
   iatom=atindx1(iatm)
   lm_size=pawfgrtab(iatom)%l_size**2
   nfgd=pawfgrtab(iatom)%nfgd
   allocate(vloc(nfgd))
   if (ngrad>0) allocate(prod(ngrad,lm_size))
   if (ngradp>0) allocate(prodp(ngradp,lm_size))

   grhat_tmp=zero

!  Eventually compute g_l(r).Y_lm(r) derivatives for the current atom (if not already done)
   if ((optgr==1.or.optstr==1).and.(optgr2/=1)) then
    if (pawfgrtab(iatom)%gylmgr_allocated==0) then
     if (associated(pawfgrtab(iatom)%gylmgr)) deallocate(pawfgrtab(iatom)%gylmgr)
     allocate(pawfgrtab(iatom)%gylmgr(3,pawfgrtab(iatom)%nfgd,lm_size))
     pawfgrtab(iatom)%gylmgr_allocated=2
     call pawgylm(rdum,pawfgrtab(iatom)%gylmgr,rdum,iatom,pawfgrtab(iatom)%ifftsph,itypat,&
&     lm_size,pawfgrtab(iatom)%nfgd,0,1,0,pawtab(itypat),pawfgrtab(iatom)%rfgd,&
&     pawfgrtab(iatom)%rfgd_allocated)
    end if
   end if
   if (optgr2==1) then
    opt1=0;opt2=0;opt3=0
    if (pawfgrtab(iatom)%gylm_allocated==0) then
     if (associated(pawfgrtab(iatom)%gylm)) deallocate(pawfgrtab(iatom)%gylm)
     allocate(pawfgrtab(iatom)%gylm(pawfgrtab(iatom)%nfgd,lm_size))
     pawfgrtab(iatom)%gylm_allocated=2;opt1=1
    end if
    if (pawfgrtab(iatom)%gylmgr_allocated==0) then
     if (associated(pawfgrtab(iatom)%gylmgr)) deallocate(pawfgrtab(iatom)%gylmgr)
     allocate(pawfgrtab(iatom)%gylmgr(3,pawfgrtab(iatom)%nfgd,lm_size))
     pawfgrtab(iatom)%gylmgr_allocated=2;opt2=1
    end if
    if (pawfgrtab(iatom)%gylmgr2_allocated==0) then
     if (associated(pawfgrtab(iatom)%gylmgr2)) deallocate(pawfgrtab(iatom)%gylmgr2)
     allocate(pawfgrtab(iatom)%gylmgr2(6,pawfgrtab(iatom)%nfgd,lm_size))
     pawfgrtab(iatom)%gylmgr2_allocated=2;opt3=1
    end if
    call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,&
&    pawfgrtab(iatom)%gylmgr2,iatom,pawfgrtab(iatom)%ifftsph,&
&    itypat,lm_size,pawfgrtab(iatom)%nfgd,opt1,opt2,opt3,pawtab(itypat),&
&    pawfgrtab(iatom)%rfgd,pawfgrtab(iatom)%rfgd_allocated)
   end if

!  Loop over spin components
   do ispden=1,nsploop

!   ----- Retrieve potential (subtle if nspden=4 ;-)
    if (nspden/=4) then
     ispvtr=min(dimvtrial,ispden)
     do ic=1,nfgd
      vloc(ic)=vtrial(pawfgrtab(iatom)%ifftsph(ic),ispvtr)
     end do
    else
     if (ispden==1) then
      ispvtr=min(dimvtrial,2)
      do ic=1,nfgd
       jc=pawfgrtab(iatom)%ifftsph(ic)
       vloc(ic)=half*(vtrial(jc,1)+vtrial(jc,ispvtr))
      end do
     else if (ispden==4) then
      ispvtr=min(dimvtrial,2)
      do ic=1,nfgd
       jc=pawfgrtab(iatom)%ifftsph(ic)
       vloc(ic)=half*(vtrial(jc,1)-vtrial(jc,ispvtr))
      end do
     else ! ispden=2 or 3
      ispvtr=min(dimvtrial,ispden+1)
      do ic=1,nfgd
       vloc(ic)=vtrial(pawfgrtab(iatom)%ifftsph(ic),ispvtr)
      end do
     end if
    end if

!   ----- Compute temporary projected scalars
    if (ngrad>0) prod=zero
    if (ngradp>0) prodp=zero
!   ==== Contribution to forces ====
    if (optgr==1) then
     do ilslm=1,lm_size
      do ic=1,pawfgrtab(iatom)%nfgd
       do mu=1,3
        prod(mu+ishift_gr,ilslm)=prod(mu+ishift_gr,ilslm)-&
&        vloc(ic)*pawfgrtab(iatom)%gylmgr(mu,ic,ilslm)
       end do
      end do
     end do
    end if
!   ==== Contribution to stresses ====
    if (optstr==1) then
     do ilslm=1,lm_size
      do ic=1,pawfgrtab(iatom)%nfgd
       do mu=1,6
        mua=alpha(mu);mub=beta(mu)
        prod(mu+ishift_str,ilslm)=prod(mu+ishift_str,ilslm) &
&        +half*vloc(ic) &
&        *(pawfgrtab(iatom)%gylmgr(mua,ic,ilslm)*pawfgrtab(iatom)%rfgd(mub,ic)&
&        +pawfgrtab(iatom)%gylmgr(mub,ic,ilslm)*pawfgrtab(iatom)%rfgd(mua,ic))
       end do
      end do
     end do
    end if
!   ==== Contribution to frozen wf part of dyn. matrix ====
    if (optgr2==1) then
     do ilslm=1,lm_size
      do ic=1,pawfgrtab(iatom)%nfgd
       do mu=1,6
        mua=alpha(mu);mub=beta(mu)
        prod(mu+ishift_gr2,ilslm)=prod(mu+ishift_gr2,ilslm) &
&        +vloc(ic)*pawfgrtab(iatom)%gylmgr2(mu,ic,ilslm) &
&        +(pawfgrtab(iatom)%vlocgr(mua,ic)*pawfgrtab(iatom)%gylmgr(mub,ic,ilslm)&
&        +pawfgrtab(iatom)%vlocgr(mub,ic)*pawfgrtab(iatom)%gylmgr(mua,ic,ilslm))
       end do
       do mu=1,3
        prodp(mu,ilslm)=prodp(mu,ilslm) &
&        -vloc(ic)*pawfgrtab(iatom)%gylmgr(mu,ic,ilslm) &
&        -pawfgrtab(iatom)%vlocgr(mu,ic)*pawfgrtab(iatom)%gylm(ic,ilslm)
       end do
      end do
     end do
    end if

!   --- Reduction in case of parallelization ---
    if(mpi_enreg%paral_compil_fft==1)then
     old_paral_level= mpi_enreg%paral_level
     mpi_enreg%paral_level=3
     call xcomm_init(mpi_enreg,spaceComm)
     if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%comm_fft
     if (ngrad>0) call xsum_mpi(prod,spaceComm,ier)
     if (ngradp>0) call xsum_mpi(prodp,spaceComm,ier)
     mpi_enreg%paral_level=old_paral_level
    end if

!   ---- Compute all gradients
    do irhoij=1,pawrhoij(iatom)%nrhoijsel
     klmn=pawrhoij(iatom)%rhoijselect(irhoij)
     klm =pawtab(itypat)%indklmn(1,klmn)
     lmin=pawtab(itypat)%indklmn(3,klmn)
     lmax=pawtab(itypat)%indklmn(4,klmn)
     ro     =pawrhoij(iatom)%rhoijp(irhoij,ispden)
     ro_d   =ro*pawtab(itypat)%dltij(klmn)
     do ils=lmin,lmax,2
      lm0=ils**2+ils+1
      do mm=-ils,ils
       ilslm=lm0+mm;isel=pawang%gntselect(lm0+mm,klm)
       if (isel>0) then
        grhat_x=ro_d*pawtab(itypat)%qijl(ilslm,klmn)
        do mu=1,ngrad
         grhat_tmp(mu)=grhat_tmp(mu)+grhat_x*prod(mu,ilslm)
        end do
       end if
      end do
     end do
    end do ! irhoij

!   ---- Add additional terms for second gradients
    if (optgr2==1) then
     do klmn=1,lmn2_size
      dlt_tmp=pawtab(itypat)%dltij(klmn)
      klm =pawtab(itypat)%indklmn(1,klmn)
      lmin=pawtab(itypat)%indklmn(3,klmn)
      lmax=pawtab(itypat)%indklmn(4,klmn)
      do ils=lmin,lmax,2
       lm0=ils**2+ils+1
       do mm=-ils,ils
        ilslm=lm0+mm;isel=pawang%gntselect(lm0+mm,klm)
        if (isel>0) then
         do mu=1,6
          mua=alpha(mu);mub=beta(mu)
          grhat_tmp(ishift_gr2+mu)=grhat_tmp(ishift_gr2+mu)&
&          +dlt_tmp*pawtab(itypat)%qijl(ilslm,klmn)&
&          *(pawrhoij(iatom)%grhoij(mua,klmn,ispden)*prodp(mub,ilslm)&
&          +pawrhoij(iatom)%grhoij(mub,klmn,ispden)*prodp(mua,ilslm))
         end do
        end if
       end do
      end do
     end do ! klmn
    end if

   end do ! ispden

!  Eventually free temporary space for g_l(r).Y_lm(r) factors
   if (pawfgrtab(iatom)%gylmgr_allocated==2) then
    deallocate(pawfgrtab(iatom)%gylmgr);allocate(pawfgrtab(iatom)%gylmgr(0,0,0))
    pawfgrtab(iatom)%gylmgr_allocated=0
   end if
   if (optgr2==1) then
    if (pawfgrtab(iatom)%gylm_allocated==2) then
     deallocate(pawfgrtab(iatom)%gylm);allocate(pawfgrtab(iatom)%gylm(0,0))
     pawfgrtab(iatom)%gylm_allocated=0
    end if
    if (pawfgrtab(iatom)%gylmgr2_allocated==2) then
     deallocate(pawfgrtab(iatom)%gylmgr2);allocate(pawfgrtab(iatom)%gylmgr2(0,0,0))
     pawfgrtab(iatom)%gylmgr2_allocated=0
    end if
   end if

!  ==== Forces ====
!  Convert from cart. to reduced coordinates
   if (optgr==1) then
    mushift=3*(iatm-1)
    do mu=1,3
     grnl(mu+mushift)=grnl(mu+mushift)&
&     +fact*(rprimd(1,mu)*grhat_tmp(ishift_gr+1)+&
&     rprimd(2,mu)*grhat_tmp(ishift_gr+2)+&
&     rprimd(3,mu)*grhat_tmp(ishift_gr+3))
    end do
   end if
!  ==== Stresses ====
!  Only store contributions (in reduced coordinates)
   if (optstr==1) then
    hatstr(1:6)=hatstr(1:6)+ grhat_tmp(ishift_str+1:ishift_str+6)
   end if
!  ==== Frozen wf part of dyn. matrix ====
!  Convert from cart. to reduced coordinates
   if (optgr2==1) then
    do mu=1,6
     mua=alpha(mu);mub=beta(mu)
     work1(mua,mub)=fact*grhat_tmp(ishift_gr2+mu)
     if (mua/=mub) work1(mub,mua)=work1(mua,mub)
    end do
    do mu=1,3
     work2(:,mu)=rprimd(1,mu)*work1(:,1)+rprimd(2,mu)*work1(:,2)+rprimd(3,mu)*work1(:,3)
    end do
    do mu=1,6
     mua=alpha(mu);mub=beta(mu)
     dytmp=rprimd(1,mua)*work2(1,mub)+rprimd(2,mua)*work2(2,mub)+rprimd(3,mua)*work2(3,mub)
     dyfrnl(mua,mub,iatm)=dyfrnl(mua,mub,iatm)+dytmp
     if (mua/=mub) dyfrnl(mub,mua,iatm)=dyfrnl(mub,mua,iatm)+dytmp
    end do
   end if

!  End loops on types and atoms
   deallocate(vloc)
   if (ngrad>0) deallocate(prod)
   if (ngradp>0) deallocate(prodp)
  end do
  iatshft=iatshft+nattyp(itypat)
 end do

!Deallocate memory
 deallocate(grhat_tmp)

!Convert stresses from reduced to cartesian coordinates
 if (optstr==1) then
! Has to compute int[nhat+*vtrial]
  hatstr_diag=zero
  if (nspden==1.or.dimvtrial==1) then
   do ic=1,nfft
    hatstr_diag=hatstr_diag+vtrial(ic,1)*nhat(ic,1)
   end do
  else if (nspden==2) then
   do ic=1,nfft
    hatstr_diag=hatstr_diag+vtrial(ic,1)*nhat(ic,2)+vtrial(ic,2)*(nhat(ic,1)-nhat(ic,2))
   end do
  else if (nspden==4) then
   do ic=1,nfft
    hatstr_diag=hatstr_diag+half*(vtrial(ic,1)*(nhat(ic,1)+nhat(ic,4)) &
&    +vtrial(ic,2)*(nhat(ic,1)-nhat(ic,4))) &
&    +vtrial(ic,3)*nhat(ic,2)+vtrial(ic,4)*nhat(ic,3)
   end do
  end if
  if(mpi_enreg%paral_compil_fft==1)then
   old_paral_level= mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%comm_fft
   call xsum_mpi(hatstr_diag,spaceComm,ier)
   mpi_enreg%paral_level=old_paral_level
  end if
! Convert hat contribution
  hatstr(1:3)=(hatstr(1:3)+hatstr_diag)/dble(nfftot)
  hatstr(4:6)= hatstr(4:6)/dble(nfftot)
! Add to already computed NL contrib
  nlstr(1:6)=nlstr(1:6)+hatstr(1:6)
! Apply symmetries
  call stresssym(gprimd,nsym,nlstr,symrec)
 end if

!DEBUG
!write(6,*)' pawgrnl: exit '
!stop
!ENDDEBUG

end subroutine pawgrnl
!!***
