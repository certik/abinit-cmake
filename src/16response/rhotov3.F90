!{\src2tex{textfont=tt}}
!!****f* ABINIT/rhotov3
!!
!! NAME
!! rhotov3
!!
!! FUNCTION
!! This routine is called to compute, from a given 1st-order total density
!! the trial (local) 1st-order potential and the residual potential.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG, DRH, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex: if 1, real space 1-order WF on FFT grid are REAL; if 2, COMPLEX
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  idir=direction of atomic displacement (=1,2 or 3 : displacement of
!!    atom ipert along the 1st, 2nd or 3rd axis).
!!  ipert=type of the perturbation
!!  kxc(nfft,nkxc)=exchange-correlation kernel
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat(cplex*nfft,nspden*usepaw)= -PAW only- compensation density
!!  nhat1(cplex*nfft,nspden*usepaw)= -PAW only- 1st-order compensation density
!!  nhat1dim=first dimension of array nhat1 (depends on perturbation type)
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  nspden=number of spin-density components
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used
!!  optene=option for the computation of additional energies
!!  optres=0: the trial potential residual is computed ; the input potential value is kept
!!         1: the new value of the trial potential is computed in place of the input value
!!  paral_kgb=flag controlling (k,g,bands) parallelization
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  qphon(3)=reduced coordinates for the phonon wavelength
!!  rhog(2,nfft)=array for Fourier transform of GS electron density
!!  rhog1(2,nfft)=RF electron density in reciprocal space
!!  rhor(nfft,nspden)=array for GS electron density in electrons/bohr**3.
!!  rhor1(cplex*nfft,nspden)=RF electron density in real space (electrons/bohr**3).
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol=unit cell volume in ($\textrm{bohr}^{3}$)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vpsp1(cplex*nfft)=first-order derivative of the ionic potential
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density, see n3xccc
!!
!! OUTPUT
!!  ==== if optene>=0
!!    ehart01=inhomogeneous 1st-order Hartree part of 2nd-order total energy
!!    ehart1=1st-order Hartree part of 2nd-order total energy
!!    exc1=1st-order exchange-correlation part of 2nd-order total energy
!!  ==== if optene==0.or.2
!!    elpsp1=1st-order local pseudopot. part of 2nd-order total energy.
!!  ==== if optene==1.or.2
!!    To be completed
!!  ==== if optres==0
!!    vresid1(cplex*nfft,nspden)=potential residual
!!    vres2=square of the norm of the residual
!!
!! SIDE EFFECTS
!!  vhartr1(cplex*nfft)=1-order Hartree potential
!!  ==== if optres==1
!!    vtrial1(cplex*nfft,nspden)= new value of 1st-order trial potential
!!
!! PARENTS
!!      scfcv3
!!
!! CHILDREN
!!      dotprod_vn,hartre,hartrestr,leave_new,mkvxc3,mkvxcstr3,sqnorm_v,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine rhotov3(cplex,ehart01,ehart1,elpsp1,exc1,gmet,gprimd,gsqcut,idir,ipert,&
&           kxc,mpi_enreg,natom,nfft,ngfft,nhat,nhat1,nhat1dim,nkxc,nspden,n3xccc,&
&           optene,optres,paral_kgb,pawfgrtab,qphon,rhog,rhog1,rhor,rhor1,&
&           rprimd,ucvol,usepaw,usexcnhat,vhartr1,vpsp1,vresid1,vres2,vtrial1,xccc3d1)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12spacepar
 use interfaces_13xc
 use interfaces_16response, except_this_one => rhotov3
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,n3xccc,natom,nfft,nhat1dim,nkxc,nspden
 integer,intent(in) :: optene,optres,paral_kgb,usepaw,usexcnhat
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(out) :: ehart01,ehart1,elpsp1,exc1,vres2
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),kxc(nfft,nkxc)
 real(dp),intent(in) :: nhat(cplex*nfft,nspden*usepaw)
 real(dp),intent(in) :: nhat1(cplex*nhat1dim,nspden*usepaw),qphon(3)
 real(dp),intent(in) :: rhog(2,nfft),rhog1(2,nfft),rhor(nfft,nspden)
 real(dp),intent(in) :: rhor1(cplex*nfft,nspden),rprimd(3,3),vpsp1(cplex*nfft)
 real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
 real(dp),intent(inout) :: vtrial1(cplex*nfft,nspden)
 real(dp),intent(out) :: vhartr1(cplex*nfft),vresid1(cplex*nfft,nspden)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)

!Local variables-------------------------------
!scalars
 integer :: ifft,ii,ii1,ii2,ispden,jj,jj1,jj2,nfftot,option
 real(dp) :: doti,elpsp10
 character(len=500) :: message
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp) :: tsec(20)
 real(dp),allocatable :: rhor1_wk(:,:),rhor_wk(:,:),vhartr01(:),vxc1(:,:)
 real(dp),allocatable :: vxc10(:,:)

! *********************************************************************

!DEBUG
!write(6,*)' rhotov3 : enter, stop '
!write(6,*)rprimd
!ENDDEBUG

!Tests
 if(nspden==4)then
  write(message, '(a,a,a,a)' )ch10,&
&  ' rhotov3 : ERROR -',ch10,&
&  '  Does not work yet for nspden=4 !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if(optene>0)then
  write(message, '(a,a,a,a)' )ch10,&
&  ' rhotov3 : BUG -',ch10,&
&  '  optene>0 not yet implemented !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 call timab(157,1,tsec)

!Get size of FFT grid
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)

!------ Compute 1st-order Hartree potential (and energy) ----------------------

 call hartre(cplex,gmet,gsqcut,0,mpi_enreg,nfft,ngfft,paral_kgb,qphon,rhog1,vhartr1)
 if (optene>=0) then
  call dotprod_vn(cplex,rhor1,ehart1,doti,mpi_enreg,nfft,nfftot,1,1,vhartr1,ucvol)
 end if

 if (optene>=0) ehart01=zero
 if(ipert==natom+3 .or. ipert==natom+4) then
  allocate(vhartr01(cplex*nfft))
  call hartrestr(gmet,gprimd,gsqcut,idir,ipert,mpi_enreg,natom,nfft,ngfft,paral_kgb,rhog,vhartr01)
  if (optene>=0) then
   call dotprod_vn(cplex,rhor1,ehart01,doti,mpi_enreg,nfft,nfftot,1,1,vhartr01,ucvol)
   ehart01=two*ehart01
   ehart1=ehart1+ehart01
  end if
! Note that there is a factor 2.0_dp difference with the similar GS formula
  vhartr1(:)=vhartr1(:)+vhartr01(:)
  deallocate(vhartr01)
 end if

!------ Compute 1st-order Hartree potential (and energy) ----------------------
!(including the XC core correction)

 option=0;if (optene<0) option=1
 allocate(vxc1(cplex*nfft,nspden),vxc10(cplex*nfft,nspden))

!Eventually substract compensation density from total density
 if (usepaw==1.and.usexcnhat==0) then
  allocate(rhor1_wk(cplex*nfft,nspden))
  if (nhat1dim==nfft) then
   rhor1_wk(:,:)=rhor1(:,:)-nhat1(:,:)
  else
   do ispden=1,nspden
    do ii=1,pawfgrtab(ipert)%nfgd
     jj=pawfgrtab(ipert)%ifftsph(ii)
     ii1=cplex*(ii-1)+1;ii2=cplex*ii
     jj1=cplex*(jj-1)+1;jj2=cplex*jj
     rhor1_wk(jj1:jj2,ispden)=rhor1(jj1:jj2,ispden)-nhat1(ii1:ii2,ispden)
    end do
   end do
  end if
 end if

 if (usepaw==1.and.usexcnhat==0) then
  if(ipert==natom+3 .or. ipert==natom+4) then
   allocate(rhor_wk(cplex*nfft,nspden));rhor_wk(:,:)=rhor(:,:)-nhat(:,:)
   call mkvxcstr3(cplex,gmet,gsqcut,idir,ipert,kxc,mpi_enreg,natom,nfft,ngfft,&
&   nkxc,nspden,n3xccc,option,paral_kgb,qphon,rhor_wk,rhor1_wk,rprimd,vxc10,xccc3d1)
   deallocate(rhor_wk)
  else
   call mkvxc3(cplex,gmet,gsqcut,kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,n3xccc,option,&
&   paral_kgb,qphon,rhor1_wk,rprimd,vxc10,xccc3d1)
  end if
 else
  if(ipert==natom+3 .or. ipert==natom+4) then
   call mkvxcstr3(cplex,gmet,gsqcut,idir,ipert,kxc,mpi_enreg,natom,nfft,ngfft,&
&   nkxc,nspden,n3xccc,option,paral_kgb,qphon,rhor,rhor1,rprimd,vxc10,xccc3d1)
  else
   call mkvxc3(cplex,gmet,gsqcut,kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,n3xccc,option,&
&   paral_kgb,qphon,rhor1,rprimd,vxc10,xccc3d1)
  end if
 end if

 if (optene==0.or.optene==2) then
  call dotprod_vn(cplex,rhor1,elpsp10,doti,mpi_enreg,nfft,nfftot,nspden,1,vxc10,ucvol)
  call dotprod_vn(cplex,rhor1,elpsp1 ,doti,mpi_enreg,nfft,nfftot,1     ,1,vpsp1,ucvol)
! Note that there is a factor 2.0_dp difference with the similar GS formula
  elpsp1=two*(elpsp1+elpsp10)
 end if

!Compute XC contribution exc1 (except the XC core-correction)
 if (optene>=0) then
  option=2
  if (usepaw==1.and.usexcnhat==0) then
   call mkvxc3(cplex,gmet,gsqcut,kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,n3xccc,option,&
&   paral_kgb,qphon,rhor1_wk,rprimd,vxc1,xccc3d1)
  else
   call mkvxc3(cplex,gmet,gsqcut,kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,n3xccc,option,&
&   paral_kgb,qphon,rhor1,rprimd,vxc1,xccc3d1)
  end if
  call dotprod_vn(cplex,rhor1,exc1,doti,mpi_enreg,nfft,nfftot,nspden,1,vxc1,ucvol)
  if(n3xccc/=0 .or. ipert==natom+3 .or. ipert==natom+4)then
   vxc1(:,:)=vxc1(:,:)+vxc10(:,:)
  end if
 else
  vxc1=vxc10
 end if

 deallocate(vxc10)
 if (usepaw==1.and.usexcnhat==0) deallocate(rhor1_wk)


!DEBUG (do not take away)
!Compute NSC energy ensc1 associated with rhor1 in vtrial1, for debugging purposes
!call dotprod_vn(cplex,rhor1,ensc1,doti,nfft,nfftot,nspden,1,vtrial1,ucvol)
!write(6,*)' ek0+eeig0+eloc0=',ek0+eeig0+eloc0
!write(6,*)' ensc1=',ensc1
!Compute NSC energy associated with vtrial1, for debugging purposes
!call dotprod_vn(cplex,rhor1,ensc1,doti,mpi_enreg,nfft,nfftot,nspden,1,vtrial1,ucvol)
!ensc1=ensc1+half*enl1
!write(6,*)' rhotov3 : check NSC energy, diff=',&
!&  ek0+edocc+eeig0+eloc0+enl0+ensc1
!write(6,*)' evarNSC=',ek0+edocc+eeig0+eloc0+enl0
!write(6,*)' ensc1,exc1=',ensc1,exc1
!ENDDEBUG

!Here, vhartr1 contains Hartree potential, vpsp1 contains local psp,
!while vxc1 contain xc potential

!------ Produce residual vector and square of norm of it -------------
!(only if requested ; if optres==0)

 if (optres==0) then
  do ispden=1,nspden
   do ifft=1,cplex*nfft
    vresid1(ifft,ispden)=vhartr1(ifft)+vxc1(ifft,ispden)+vpsp1(ifft)-vtrial1(ifft,ispden)
   end do
  end do

! Compute square norm vres2 of potential residual vresid
  call sqnorm_v(cplex,mpi_enreg,nfft,vres2,nspden,vresid1)

 else

! ------ Produce new value of trial potential-------------
! (only if requested ; if optres==1)

  do ispden=1,nspden
   do ifft=1,cplex*nfft
    vtrial1(ifft,ispden)=vhartr1(ifft)+vxc1(ifft,ispden)+vpsp1(ifft)
   end do
  end do

 end if

 deallocate(vxc1)

 call timab(157,2,tsec)

!DEBUG
!write(6,*)' rhotov3 : exit '
!stop
!ENDDEBUG

end subroutine rhotov3
!!***
