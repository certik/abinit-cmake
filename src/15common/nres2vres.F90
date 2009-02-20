!{\src2tex{textfont=tt}}
!!****f* ABINIT/nres2vres
!!
!! NAME
!! nres2vres
!!
!! FUNCTION
!! Convert a density residual into a potential residual
!! using a first order formula:
!!     V^res(r)=dV/dn.n^res(r)
!!             =V_hartree(n^res)(r) + Kxc.n^res(r)
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dtset <type(dataset_type)>=all input variables in this dataset
!!  | icoulomb=0 periodic treatment of Hartree potential, 1 use of Poisson solver
!!  | natom= number of atoms in cell
!!  | nspden=number of spin-density components
!!  | ntypat=number of atom types
!!  | typat(natom)=type (integer) for each atom
!! gsqcut=cutoff value on G**2 for sphere inside fft box
!! izero=if 1, unbalanced components of Vhartree(g) are set to zero
!! kxc(nfft,nkxc)=exchange-correlation kernel, needed only if nkxc>0
!! mpi_enreg=informations about MPI parallelization
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT
!! nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!! nkxc=second dimension of the array kxc, see rhohxc.F90 for a description
!! nresid(nfft,nspden)= the input density residual
!! n3xccc=dimension of the xccc3d array (0 or nfft).
!! optnc=option for non-collinear magnetism (nspden=4):
!!       1: the whole 2x2 Vres matrix is computed
!!       2: only Vres^{11} and Vres^{22} are computed
!! optxc=0 if LDA part of XC kernel has only to be taken into account (even for GGA)
!!       1 if XC kernel has to be fully taken into
!!      -1 if XC kernel does not have to be taken into account
!! pawang <type(pawang_type)>=paw angular mesh and related data
!! pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!! pawrhoij(natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! rhog(2,nfft)=electron density in reciprocal space
!!              (used only if Kxc was not computed before)
!! rhor(nfft,nspden)=electron density in real space
!!                   (used only if Kxc was not computed before)
!! rprimd(3,3)=dimensional primitive translation vectors (bohr)
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!! xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!
!! OUTPUT
!! vresid(nfft,nspden)= the output potential residual
!!
!! PARENTS
!!      etotfor,forstr
!!
!! CHILDREN
!!      atmdata,fourdp,hartre,leave_new,metric,mkvxc3,pawmknhat,rhohxc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine nres2vres(dtset,gsqcut,izero,kxc,mpi_enreg,nfft,ngfft,nhat,&
&                 nkxc,nresid,n3xccc,optnc,optxc,pawang,pawfgrtab,pawrhoij,pawtab,&
&                 rhog,rhor,rprimd,usepaw,usexcnhat,vresid,xccc3d)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_13paw
 use interfaces_13xc
 use interfaces_14poisson
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: izero,n3xccc,nfft,nkxc,optnc,optxc,usepaw,usexcnhat
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: kxc(nfft,nkxc),nresid(nfft,dtset%nspden),rhog(2,nfft)
 real(dp),intent(in) :: rhor(nfft,dtset%nspden),rprimd(3,3),xccc3d(n3xccc)
 real(dp),intent(inout) :: nhat(nfft,dtset%nspden*usepaw)
 real(dp),intent(out) :: vresid(nfft,dtset%nspden)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%ntypat*usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: ifft,ikxc,ispden,nhatgrdim,nkxc_cur,option
 real(dp) :: dum,dvdn,dvdz,energy,fact,m_dot_mres,ucvol,vxcavg
 character(len=500) :: message
!arrays
 real(dp) :: dummy6(6),gmet(3,3),gprimd(3,3),qq(3),rmet(3,3)
 real(dp),allocatable :: dummy(:),kxc_cur(:,:),m_norm(:),nhatgr(:,:,:)
 real(dp),allocatable :: nresg(:,:),nresid_diag(:,:),rhor0(:,:),vhres(:)
 real(dp),allocatable :: vresid_diag(:,:)

! *************************************************************************

!DEBUG
!write(6,*)' nres2vres : enter ',optxc
!ENDDEBUG

!Compatibility tests:
 if(optxc<-1.or.optxc>1)then
  write(message, '(4a)') ch10,&
&  ' nres2vres : BUG -',ch10,&
  '   Wrong value for optxc !'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if
 if((optnc/=1.and.optnc/=2).or.(dtset%nspden/=4.and.optnc/=1))then
  write(message, '(4a)') ch10,&
&  ' nres2vres : BUG -',ch10,&
  '   Wrong value for optnc !'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if
 if(dtset%icoulomb==1.and.optxc/=-1)then
  write(message, '(4a)') ch10,&
&  ' nres2vres : ERROR -',ch10,&
  '   This routine is not compatible with icoulomb==1 and optxc/=-1 !'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if
 if(dtset%nspden==4.and.dtset%xclevel==2.and.optxc==1.and.nkxc/=23)then
  write(message, '(4a)') ch10,&
&  ' nres2vres : ERROR -',ch10,&
  '   Wrong values for optxc and nkxc !'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if

 qq=zero
 nkxc_cur=0
 if (dtset%xclevel==1.or.optxc==0) nkxc_cur=3-2*mod(dtset%nspden,2)
 if (dtset%xclevel==2.and.optxc==1) nkxc_cur=23
 allocate(vhres(nfft))

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Compute density residual in reciprocal space
 if (dtset%icoulomb==0) then
  allocate(nresg(2,nfft),dummy(nfft));dummy(:)=nresid(:,1)
  call fourdp(1,nresg,dummy,-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
  deallocate(dummy)
 end if

!First case: Kxc has already been computed
!-----------------------------------------
 if (nkxc==nkxc_cur.or.optxc==-1) then

! Compute VH(n^res)(r)
  if (dtset%icoulomb == 0) then
   call hartre(1,gmet,gsqcut,izero,mpi_enreg,nfft,ngfft,dtset%paral_kgb,qq,nresg,vhres)
  else
   call PSolver_hartree(dtset,energy,mpi_enreg,nresid(:,1),rprimd,vhres)
  end if

! Compute Kxc(r).n^res(r)
  if (optxc/=-1) then

!  Collinear magnetism or non-polarized
   if (dtset%nspden/=4) then
    call mkvxc3(1,gmet,gsqcut,kxc,mpi_enreg,nfft,ngfft,nkxc,&
&    dtset%nspden,0,2,dtset%paral_kgb,qq,nresid,rprimd,vresid,dummy)

!   Non-collinear magnetism
!   Has to locally "rotate" n^res(r) (according to magnetization),
!   compute V^res(r) and rotate it back
   else
    allocate(nresid_diag(nfft,2),vresid_diag(nfft,2),&
&    rhor0(nfft,dtset%nspden),m_norm(nfft))
!   -- Compute "initial" density
    rhor0(:,:)=rhor(:,:)-nresid(:,:)
!   -- Rotate n^res(r)
    do ifft=1,nfft
     nresid_diag(ifft,1)=nresid(ifft,1)
     m_norm(ifft)=sqrt(rhor0(ifft,2)**2+rhor0(ifft,3)**2+rhor0(ifft,4)**2)
     m_dot_mres=rhor0(ifft,2)*nresid(ifft,2)+rhor0(ifft,3)*nresid(ifft,3) &
&              +rhor0(ifft,4)*nresid(ifft,4)
     if(m_norm(ifft)>tol14)then
      nresid_diag(ifft,2)=half*(nresid_diag(ifft,1)+m_dot_mres/m_norm(ifft))
     else
      nresid_diag(ifft,2)=nresid_diag(ifft,1)
     end if
    end do
!   -- Compute Kxc(r).n^res(r)_rotated
    call mkvxc3(1,gmet,gsqcut,kxc,mpi_enreg,nfft,ngfft,nkxc,&
&    2,0,2,dtset%paral_kgb,qq,nresid_diag,rprimd,vresid_diag,dummy)
    deallocate(nresid_diag)
!   -- Rotate back V^res(r)
    if (optnc==1) then
     do ifft=1,nfft
      dvdn=(vresid_diag(ifft,1)+vresid_diag(ifft,2))*half
      dvdz=(vresid_diag(ifft,1)-vresid_diag(ifft,2))*half
      if(m_norm(ifft)>tol14)then
       fact=dvdz/m_norm(ifft)
       dum=rhor0(ifft,4)*fact
       vresid(ifft,1)=dvdn+dum
       vresid(ifft,2)=dvdn-dum
       vresid(ifft,3)= rhor0(ifft,2)*fact
       vresid(ifft,4)=-rhor0(ifft,3)*fact
      else
       vresid(ifft,1:2)=dvdn
       vresid(ifft,3:4)=zero
      end if
     end do
    else
     do ifft=1,nfft
      dvdn=(vresid_diag(ifft,1)+vresid_diag(ifft,2))*half
      dvdz=(vresid_diag(ifft,1)-vresid_diag(ifft,2))*half
      if(m_norm(ifft)>tol14)then
       dum=dvdz*rhor0(ifft,4)/m_norm(ifft)
       vresid(ifft,1)=dvdn+dum
       vresid(ifft,2)=dvdn-dum
      else
       vresid(ifft,1:2)=dvdn
      end if
     end do
    end if
    deallocate(vresid_diag,rhor0,m_norm)
   end if

  else
   vresid=zero
  end if

 end if

!2nd case: Kxc has to be computed
!--------------------------------
 if (nkxc/=nkxc_cur.and.optxc/=-1) then
! For GGA, has to recompute gradients of nhat
  nhatgrdim=0
  if (usepaw==1.and.usexcnhat>0.and.dtset%xclevel==2.and.dtset%pawnhatxc>0) then
   nhatgrdim=1;allocate(nhatgr(nfft,dtset%nspden,3))
   call pawmknhat(dum,1,0,mpi_enreg,dtset%natom,nfft,ngfft,nhatgrdim,dtset%nspden,dtset%ntypat,&
&   dtset%paral_kgb,pawang,pawfgrtab,nhatgr,nhat,pawrhoij,pawtab,dtset%typat,ucvol)
  end if

! Has to use the "initial" density to compute Kxc
  allocate(rhor0(nfft,dtset%nspden))
  rhor0(:,:)=rhor(:,:)-nresid(:,:)

! Compute VH(n^res) and XC kernel (Kxc) together
  allocate(kxc_cur(nfft,nkxc_cur))
  option=2;if (dtset%xclevel==2.and.optxc==0) option=12
  call rhohxc(dtset,energy,gsqcut,izero,kxc_cur,mpi_enreg,nfft,ngfft,&
&  nhat,usepaw,nhatgr,nhatgrdim,nkxc_cur,dtset%nspden,n3xccc,option,nresg,&
&  rhor0,rprimd,dummy6,usexcnhat,vhres,vresid,vxcavg,xccc3d)  !vresid=work space
  if (dtset%nspden/=4) deallocate(rhor0)
  if (nhatgrdim>0) deallocate(nhatgr)

! Compute Kxc(r).n^res(r)

! Collinear magnetism or non-polarized
  if (dtset%nspden/=4) then
   call mkvxc3(1,gmet,gsqcut,kxc_cur,mpi_enreg,nfft,ngfft,nkxc_cur,&
&   dtset%nspden,0,2,dtset%paral_kgb,qq,nresid,rprimd,vresid,dummy)

!  Non-collinear magnetism
!  Has to locally "rotate" n^res(r) (accroding to magnetization),
!  compute V^res(r) and rotate it back
  else
   allocate(nresid_diag(nfft,2),vresid_diag(nfft,2),m_norm(nfft))
!  -- Rotate n^res(r)
   do ifft=1,nfft
    nresid_diag(ifft,1)=nresid(ifft,1)
    m_norm(ifft)=sqrt(rhor0(ifft,2)**2+rhor0(ifft,3)**2+rhor0(ifft,4)**2)
    m_dot_mres=rhor0(ifft,2)*nresid(ifft,2)+rhor0(ifft,3)*nresid(ifft,3) &
&    +rhor0(ifft,4)*nresid(ifft,4)
    nresid_diag(ifft,2)=half*(nresid_diag(ifft,1)+m_dot_mres/m_norm(ifft))
   end do
!  -- Compute Kxc(r).n^res(r)_rotated
   call mkvxc3(1,gmet,gsqcut,kxc_cur,mpi_enreg,nfft,ngfft,nkxc_cur,&
&   2,0,2,dtset%paral_kgb,qq,nresid_diag,rprimd,vresid_diag,dummy)
   deallocate(nresid_diag)
!  -- Rotate back V^res(r)
   if (optnc==1) then
    do ifft=1,nfft
     dvdn=(vresid_diag(ifft,1)+vresid_diag(ifft,2))*half
     dvdz=(vresid_diag(ifft,1)-vresid_diag(ifft,2))*half
     if(m_norm(ifft)>tol14)then
      fact=dvdz/m_norm(ifft)
      dum=rhor0(ifft,4)*fact
      vresid(ifft,1)=dvdn+dum
      vresid(ifft,2)=dvdn-dum
      vresid(ifft,3)= rhor0(ifft,2)*fact
      vresid(ifft,4)=-rhor0(ifft,3)*fact
     else
      vresid(ifft,1:2)=dvdn
      vresid(ifft,3:4)=zero
     end if
    end do
   else
    do ifft=1,nfft
     dvdn=(vresid_diag(ifft,1)+vresid_diag(ifft,2))*half
     dvdz=(vresid_diag(ifft,1)-vresid_diag(ifft,2))*half
     if(m_norm(ifft)>tol14)then
      dum=dvdz*rhor0(ifft,4)/m_norm(ifft)
      vresid(ifft,1)=dvdn+dum
      vresid(ifft,2)=dvdn-dum
     else
      vresid(ifft,1:2)=dvdn
     end if
    end do
   end if
   deallocate(vresid_diag,m_norm,rhor0)
  end if

  deallocate(kxc_cur)
 end if

!Assemble potential residual: V^res(r)=VH(n^res)(r) + Kxc(r).n^res(r)
!--------------------------------------------------------------------
 do ispden=1,dtset%nspden/optnc
  vresid(:,ispden)=vresid(:,ispden)+vhres(:)
 end do

 if (dtset%icoulomb==0) deallocate(nresg)
 deallocate(vhres)

end subroutine nres2vres

!!***
