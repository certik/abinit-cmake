!{\src2tex{textfont=tt}}
!!****f* ABINIT/nhatgrid
!! NAME
!! nhatgrid
!!
!! FUNCTION
!! Determine parts of the rectangular (fine) grid that are contained
!! inside spheres around atoms (used to compute n_hat density).
!! If corresponding option is selected, compute also g_l(r)*Y_lm(r)
!! (and derivatives) on this grid (g_l=radial shape function).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms in unit cell
!!  optgr0= 1 if g_l(r)*Y_lm(r) are computed
!!  optgr1= 1 if first derivatives of g_l(r)*Y_lm(r) are computed
!!  optgr2= 1 if second derivatives of g_l(r)*Y_lm(r) are computed
!!  optrad= 1 if vectors (r-r_atom) on the fine grid around atoms have to be stored
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  typat(natom)=type (integer) for each atom
!!  ucvol=unit cell volume in bohr**3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  pawfgrtab(natom)%ifftsph(nfgd)=FFT index (fine grid) of a points in paw spheres around each atom
!!  pawfgrtab(natom)%nfgd= number of (fine grid) FFT points in paw spheres around atoms
!!  if (optgr0==1)
!!    pawfgrtab(natom)%gylm(nfgd,l_size**2)= g_l(r)*Y_lm(r) around each atom
!!  if (optgr1==1)
!!    pawfgrtab(natom)%gylmgr(3,nfgd,l_size**2)= derivatives of g_l(r)*Y_lm(r) wrt cart. coordinates
!!  if (optgr2==1)
!!    pawfgrtab(natom)%gylmgr2(6,nfgd,l_size**2)= second derivatives of g_l(r)*Y_lm(r) wrt cart. coordinates
!!  if (optrad==1)
!!    pawfgrtab(natom)%rfgd(3,nfgd)= coordinates of r-r_atom around each atom
!!
!! PARENTS
!!      afterscfloop,respfn,scfcv,screening,sigma
!!
!! CHILDREN
!!      leave_new,pawgylm,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine nhatgrid(atindx1,gmet,mpi_enreg,natom,nattyp,nfft,ngfft,ntypat,&
& optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,rprimd,typat,ucvol,xred)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_13paw, except_this_one => nhatgrid
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,ntypat,optgr0,optgr1,optgr2,optrad
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),nattyp(ntypat),ngfft(18),typat(natom)
 real(dp),intent(in) :: gmet(3,3),rprimd(3,3),xred(3,natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ------------------------------
!scalars
 integer,parameter :: ishift=5
 integer :: i1,i2,i3,iat,iatm,iatom,ifft_local,itypat,ix,iy,iz,izmod,lm_size
 integer :: me_fft,n1,n1a,n1b,n2,n2a,n2b,n3,n3a,n3b,ncmax,nd3,nfftot,nfgd
 real(dp),parameter :: delta=0.99_dp
 real(dp) :: difx,dify,difz,r2,r2shp,rr1,rr2,rr3,rshp,rx,ry,rz
 character(len=500) :: message
!arrays
 integer,allocatable :: ifftsph_tmp(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: rr(:,:)

! *************************************************************************

 call timab(559,1,tsec)

 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3);nd3=n3/mpi_enreg%nproc_fft
 nfftot=n1*n2*n3
 me_fft=mpi_enreg%me_fft

!Loop over types of atom
!-------------------------------------------
 iatm=0
 do itypat=1,ntypat

  rshp=pawtab(itypat)%rshp
  r2shp=1.0000001_dp*rshp**2
  ncmax=1+int(1.2_dp*nfftot*four_pi/(three*ucvol)*rshp**3)
  allocate(ifftsph_tmp(ncmax),rr(3,ncmax))
  rr1=sqrt(r2shp*gmet(1,1))
  rr2=sqrt(r2shp*gmet(2,2))
  rr3=sqrt(r2shp*gmet(3,3))

! Loop over atoms
! -------------------------------------------
  do iat=1,nattyp(itypat)
   iatm=iatm+1;iatom=atindx1(iatm)

!  Define a "box" around each atom
   n1a=int((xred(1,iatom)-rr1+ishift)*n1+delta)-ishift*n1
   n1b=int((xred(1,iatom)+rr1+ishift)*n1      )-ishift*n1
   n2a=int((xred(2,iatom)-rr2+ishift)*n2+delta)-ishift*n2
   n2b=int((xred(2,iatom)+rr2+ishift)*n2      )-ishift*n2
   n3a=int((xred(3,iatom)-rr3+ishift)*n3+delta)-ishift*n3
   n3b=int((xred(3,iatom)+rr3+ishift)*n3      )-ishift*n3

!  ------------------------------------------------------------------
!  A-Identify the fft indexes of the rectangular grid around the atom
!  ------------------------------------------------------------------

   nfgd=0
   do i3=n3a,n3b
    iz=mod(i3+ishift*n3,n3)
    if(iz/nd3==me_fft) then !MPIWF
     izmod=modulo(iz,nd3)
     difz=dble(i3)/dble(n3)-xred(3,iatom)
     do i2=n2a,n2b
      iy=mod(i2+ishift*n2,n2)
      dify=dble(i2)/dble(n2)-xred(2,iatom)
      do i1=n1a,n1b
       ix=mod(i1+ishift*n1,n1)
       difx=dble(i1)/dble(n1)-xred(1,iatom)
       rx=difx*rprimd(1,1)+dify*rprimd(1,2)+difz*rprimd(1,3)
       ry=difx*rprimd(2,1)+dify*rprimd(2,2)+difz*rprimd(2,3)
       rz=difx*rprimd(3,1)+dify*rprimd(3,2)+difz*rprimd(3,3)
       r2=rx**2+ry**2+rz**2

       if (r2 <= r2shp) then
        ifft_local=1+ix+n1*(iy+n2*izmod)
        nfgd=nfgd+1
        if (nfgd>ncmax) then
         write(message, '(a,a,a,a,i4,a,a,i4,a)' ) ch10,&
&         ' nhatgrid : BUG -',ch10,&
&         '   Number of fft points around atom ',iatom,ch10,&
&         '   exceeds max. allowed (',ncmax,') !'
         call wrtout(06,message,'COLL')
         call leave_new('PERS')
        end if

        rr(1,nfgd)=rx
        rr(2,nfgd)=ry
        rr(3,nfgd)=rz
        ifftsph_tmp(nfgd)=ifft_local
       end if

      end do
     end do
    end if
   end do

   lm_size=pawfgrtab(iatom)%l_size**2

!  Allocate arrays defining sphere (and related data) around current atom
   if (associated(pawfgrtab(iatom)%ifftsph))deallocate(pawfgrtab(iatom)%ifftsph)
   allocate(pawfgrtab(iatom)%ifftsph(nfgd))
   pawfgrtab(iatom)%nfgd=nfgd
   pawfgrtab(iatom)%ifftsph(1:nfgd)=ifftsph_tmp(1:nfgd)

   if (optrad==1) then
    if (associated(pawfgrtab(iatom)%rfgd)) deallocate(pawfgrtab(iatom)%rfgd)
    allocate(pawfgrtab(iatom)%rfgd(3,nfgd))
    pawfgrtab(iatom)%rfgd_allocated=1
    pawfgrtab(iatom)%rfgd(1:3,1:nfgd)=rr(1:3,1:nfgd)
   end if

   if (optgr0==1) then
    if (associated(pawfgrtab(iatom)%gylm)) deallocate(pawfgrtab(iatom)%gylm)
    allocate(pawfgrtab(iatom)%gylm(nfgd,lm_size))
    pawfgrtab(iatom)%gylm_allocated=1
   end if

   if (optgr1==1) then
    if (associated(pawfgrtab(iatom)%gylmgr)) deallocate(pawfgrtab(iatom)%gylmgr)
    allocate(pawfgrtab(iatom)%gylmgr(3,nfgd,lm_size))
    pawfgrtab(iatom)%gylmgr_allocated=1
   end if

   if (optgr2==1) then
    if (associated(pawfgrtab(iatom)%gylmgr2)) deallocate(pawfgrtab(iatom)%gylmgr2)
    allocate(pawfgrtab(iatom)%gylmgr2(6,nfgd,lm_size))
    pawfgrtab(iatom)%gylmgr2_allocated=1
   end if

!  ------------------------------------------------------------------
!  B-Calculate g_l(r)*Y_lm(r) for each r around the atom (if requested)
!  ------------------------------------------------------------------

   call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,pawfgrtab(iatom)%gylmgr2,&
&   iatom,ifftsph_tmp(1:nfgd),itypat,lm_size,nfgd,optgr0,optgr1,optgr2,&
&   pawtab(itypat),rr(:,1:nfgd),1)

!  End loops over types/atoms
!  -------------------------------------------
  end do
  deallocate(ifftsph_tmp,rr)
 end do

 call timab(559,2,tsec)

end subroutine nhatgrid

!!***
