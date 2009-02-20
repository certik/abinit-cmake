!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp7in
!! NAME
!! psp7in
!!
!! FUNCTION
!! Initialize pspcod=7 ("PAW pseudopotentials"):
!! continue to read the corresponding file and compute the form factors
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (GJ, FJ, MT, FB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ixc=exchange-correlation choice from main routine data file
!!  lmax=value of lmax mentioned at the second line of the psp file
!!  lmnmax=max number of (l,m,n) components over all type of psps
!!  lnmax=max number of (l,n)   components over all type of psps
!!  mmax=max number of pts in real space grid (already read in the psp file header)
!!  mqgrid_ff=dimension of q grid for array ffspl (spline fit of nl form factors)
!!  mqgrid_vl=dimension of q grid for array vlspl (spline fit of Vloc)
!!  pawxcdev=choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  pspso= spin orbit signal
!!  qgrid_ff(mqgrid_ff)=values of q on grid (from 0 to qmax) for nl form factors
!!  qgrid_vl(mqgrid_vl)=values of q on grid (from 0 to qmax) for Vloc
!!  xclevel= XC functional level
!!  zion=nominal valence of atom as specified in psp file
!!
!! OUTPUT
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r}) dr]$(hartree)
!!  ffspl(mqgrid_ff,2,lnmax)=form factor f_l(q) and second derivative
!!   from spline fit for each angular momentum and each projector;
!!   if any, spin-orbit components begin at l=mpsang+1
!!  indlmn(6,lmnmax)= array giving l,m,n,lm,ln,s for i=lmn
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  vlspl(mqgrid_vl,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  xcccrc=XC core correction cutoff radius (bohr) from psp file
!!
!! NOTES
!!  Spin-orbit not yet implemented (to be done)
!!
!! PARENTS
!!      pspatm
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine psp7in(epsatm,ffspl,indlmn,ixc,lmax,lmnmax,lnmax,mmax,&
&                  mqgrid_ff,mqgrid_vl,pawrad,pawtab,pawxcdev,&
&                  pspso,qgrid_ff,qgrid_vl,vlspl,xcccrc,xclevel,zion)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13paw
 use interfaces_13psp, except_this_one => psp7in
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,lmax,lmnmax,lnmax,mmax,mqgrid_ff,mqgrid_vl,pawxcdev
 integer,intent(in) :: pspso,xclevel
 real(dp),intent(in) :: zion
 real(dp),intent(out) :: epsatm,xcccrc
 type(pawrad_type),intent(out) :: pawrad
 type(pawtab_type),intent(out) :: pawtab
!arrays
 integer,intent(out) :: indlmn(6,lmnmax)
 real(dp),intent(in) :: qgrid_ff(mqgrid_ff),qgrid_vl(mqgrid_vl)
 real(dp),intent(out) :: ffspl(mqgrid_ff,2,lnmax)
 real(dp),intent(out) :: vlspl(mqgrid_vl,2)

!Local variables ------------------------------
!scalars
 integer,parameter :: reduced_mshsz=2501
 integer :: creatorid,ib,icoremesh,ii,il,ilm,ilmn,iln,imainmesh,imsh,iprojmesh
 integer :: ir,iread1,iread2,ishpfmesh,isnotzero,itest,ivalemesh,ivlocmesh,j0lmn,jj,jlm
 integer :: jlmn,jln,klmn,msz,nmesh,pspversion
 real(dp),parameter :: reduced_rstep=0.00025_dp,rm_vloc=20.0_dp
 real(dp) :: arg,argj1,argj2,intg,intvh,jbes1,jbes2,qcore,qq,rc,rread1,rread2,yp1,ypn
 logical :: reduced_ncor,reduced_nval,reduced_vloc
 character :: blank=' ',numb=' '
 character(len=500) :: message,submessage
 character(len=80) :: pspline
 type(pawang_type) :: pawang_tmp
 type(pawrad_type) :: core_mesh,rcore_mesh,rvale_mesh,rvloc_mesh,shpf_mesh
 type(pawrad_type) :: tproj_mesh,vale_mesh,vloc_mesh
!arrays
 integer,allocatable :: nprj(:),orbitals(:)
 real(dp),allocatable :: ncore(:),nwk(:),r2k(:),rtncor(:),rtnval(:),rvlocr(:)
 real(dp),allocatable :: shpf(:,:),tncore(:),tnvale(:),tproj(:,:),vbare(:),vh(:),vlocr(:)
 real(dp),allocatable :: work1(:),work2(:),work3(:),work4(:)
 logical :: tmp_lmselect(1)
 type(pawrad_type),allocatable :: radmesh(:)

!************************************************************************

!File format of formatted PAW psp input (the 3 first lines
!have already been read in calling -pspatm- routine) :

!(1) title (character) line
!(2) znucl, zion, pspdat
!(3) pspcod, pspxc, lmax, lloc, mmax, r2well
!(4) psp_version, creatorID
!(5) basis_size, lmn_size
!(6) orbitals (for l=1 to basis_size)
!(7) number_of_meshes
!For imsh=1 to number_of_meshes
!(8)  mesh_index, mesh_type ,mesh_size, rad_step[, log_step]
!(9) r_cut(SPH)
!(10) shape_type, r_shape[, shapefunction arguments]
!For iln=1 to basis_size
!(11) comment(character)
!(12) radial mesh index for phi
!(13) phi(r) (for ir=1 to phi_meshsz)
!For iln=1 to basis_size
!(14) comment(character)
!(15) radial mesh index for tphi
!(16) tphi(r) (for ir=1 to phi_mesh_size)
!For iln=1 to basis_size
!(17) comment(character)
!(18) radial mesh index for tproj
!(19) tproj(r) (for ir=1 to proj_mesh_size)
!(20) comment(character)
!(21) radial mesh index for core_density
!(22) core_density (for ir=1 to core_mesh_size)
!(23) comment(character)
!(24) radial mesh index for pseudo_core_density
!(25) tcore_density (for ir=1 to core_mesh_size)
!(26) comment(character)
!(27) Dij0 (for ij=1 to lmn_size*(lmn_size+1)/2)
!(28) comment(character)
!(29) Rhoij0 (for ij=1 to lmn_size*(lmn_size+1)/2)
!(30) comment(character)
!(31) radial mesh index for Vloc, format of Vloc (0=Vbare, 1=VH(tnzc))
!(32) Vloc(r) (for ir=1 to vloc_mesh_size)
!===== Following lines only if shape_type=-1 =====
!For il=1 to 2*max(orbitals)+1
!(33) comment(character)
!(34) radial mesh index for shapefunc
!(35) shapefunc(r)*gnorm(l)*r**l (for ir=1 to shape_mesh_size)
!(36) comment(character)
!(37) radial mesh index for pseudo_valence_density
!(38) tvale(r) (for ir=1 to vale_mesh_size)
!
!Comments:
!* psp_version= ID of PAW_psp version
!4 characters string of the form 'pawn' (with n varying)
!* creatorID= ID of psp generator
!creatorid=1xyz : psp generated from Holzwarth AtomPAW generator version x.yz
!creatorid=2xyz : psp generated from Vanderbilt ultra-soft generator version x.yz
!creatorid=-1: psp for tests (for developpers only)
!* mesh_type= type of radial mesh
!mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!* radial shapefunction type
!shape_type=-1 ; gl(r)=numeric (read from psp file)
!shape_type= 1 ; gl(r)=k(r).r^l; k(r)=exp[-(r/sigma)**lambda]
!shape_type= 2 ; gl(r)=k(r).r^l; k(r)=[sin(pi*r/rshp)/(pi*r/rshp)]**2 if r<=rshp
!shape_type= 3 ; gl(r)=Alpha(1,l)*jl(q(1,l)*r)+Alpha(2,l)*jl(q(2,l)*r) for each l


!==========================================================

!DEBUG
!write(6,*)' psp7in : enter'
!ENDDEBUG

!==========================================================
!Read psp version in line 4 of the header
 pspversion=1
 read (tmp_unit,'(a80)') pspline;pspline=adjustl(pspline)
 if (pspline(1:3)=="paw".or.pspline(1:3)=="PAW") &
& read(unit=pspline(4:80),fmt=*) pspversion
 if (pspversion<1.or.pspversion>4) then
  write(message, '(a,a,a,a,i2,a,a,a)' ) ch10,&
&  ' psp7in: ERROR -',ch10,&
&  '  This version of PAW psp file (',pspversion,') is not compatible with',ch10,&
&  '  current version of Abinit.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!==========================================================
!Read lines 4 to 11 of the header

!Have to maintain compatibility with Abinit v4.2.x
 if (pspversion==1) then
  read (unit=pspline,fmt=*) pawtab%basis_size,pawtab%lmn_size
  allocate(orbitals(pawtab%basis_size))
  read(tmp_unit,*) (orbitals(ib), ib=1,pawtab%basis_size)
  pawtab%l_size=2*maxval(orbitals)+1
  nmesh=3;allocate(radmesh(nmesh))
  read(tmp_unit,'(a80)') pspline
  radmesh(1)%lstep=zero
  read(unit=pspline,fmt=*,err=10,end=10) radmesh(1)%mesh_type,&
&  radmesh(1)%rstep,radmesh(1)%lstep
  10 read(tmp_unit,*) pawtab%rpaw
  read(tmp_unit,*) radmesh(1)%mesh_size,radmesh(2)%mesh_size,&
&  radmesh(3)%mesh_size
  read(tmp_unit,'(a80)') pspline
  pawtab%shape_lambda=-1;pawtab%shape_sigma=1.d99
  read(unit=pspline,fmt=*,err=11,end=11) pawtab%shape_type,&
&  pawtab%shape_lambda,pawtab%shape_sigma
  11 read(tmp_unit,*) creatorid
  if (pawtab%shape_type==3) pawtab%shape_type=-1
  radmesh(2)%mesh_type=radmesh(1)%mesh_type
  radmesh(3)%mesh_type=radmesh(1)%mesh_type
  radmesh(2)%rstep=radmesh(1)%rstep
  radmesh(3)%rstep=radmesh(1)%rstep
  radmesh(2)%lstep=radmesh(1)%lstep
  radmesh(3)%lstep=radmesh(1)%lstep
 else

! Here psp file for Abinit 4.3+
  read (unit=pspline(5:80),fmt=*) creatorid
  read (tmp_unit,*) pawtab%basis_size,pawtab%lmn_size
  allocate(orbitals(pawtab%basis_size))
  read(tmp_unit,*) (orbitals(ib), ib=1,pawtab%basis_size)
  pawtab%l_size=2*maxval(orbitals)+1
  read(tmp_unit,*) nmesh
  allocate(radmesh(nmesh))
  do imsh=1,nmesh
   rread2=zero
   read(tmp_unit,'(a80)') pspline
   read(unit=pspline,fmt=*,err=20,end=20) ii,iread1,iread2,rread1,rread2
   20 continue
   if (ii<=nmesh) then
    radmesh(ii)%mesh_type=iread1
    radmesh(ii)%mesh_size=iread2
    radmesh(ii)%rstep=rread1
    radmesh(ii)%lstep=rread2
   else
    write(message, '(a,a,a,a,a,a)' ) ch10,&
&    ' psp7in: ERROR -',ch10,&
&    '  Index of mesh out of range !',ch10,&
&    '  Action : check your pseudopotential file.'
    call wrtout(ab_out,message,'COLL')
    call wrtout(06,  message,'COLL')
    call leave_new('COLL')
   end if
  end do
  read(tmp_unit,*) pawtab%rpaw
  read(tmp_unit,'(a80)') pspline
  read(unit=pspline,fmt=*) pawtab%shape_type
  pawtab%shape_lambda=-1;pawtab%shape_sigma=1.d99
 end if

!Here reading shapefunction parameters for Abinit 4.3...4.5
 if (pspversion==2) then
  if (pawtab%shape_type==1) read(unit=pspline,fmt=*) ii,pawtab%shape_lambda,pawtab%shape_sigma
  if (pawtab%shape_type==3) pawtab%shape_type=-1
  pawtab%rshp=zero

! Here reading shapefunction parameters for Abinit 4.6+
 else if (pspversion>=3) then
  pawtab%rshp=zero
  if (pawtab%shape_type==-1) read(unit=pspline,fmt=*,err=21,end=21) ii,pawtab%rshp
  if (pawtab%shape_type== 1) read(unit=pspline,fmt=*,err=21,end=21) ii,pawtab%rshp, &
&  pawtab%shape_lambda,pawtab%shape_sigma
  if (pawtab%shape_type== 2) read(unit=pspline,fmt=*,err=21,end=21) ii,pawtab%rshp
  if (pawtab%shape_type== 3) read(unit=pspline,fmt=*,err=21,end=21) ii,pawtab%rshp
 end if
 21 continue

!If shapefunction type is Bessel, deduce here its parameters from rc
 if (pawtab%shape_type== 3) then
  rc=pawtab%rshp;if (rc<1.d-8) rc=pawtab%rpaw
  do il=1,pawtab%l_size
   call shapebes(pawtab%shape_alpha(1:2,il),pawtab%shape_q(1:2,il),il-1,rc)
  end do
 end if

!==========================================================
!Mirror pseudopotential parameters to the output and log files

 write(message,'(a,i1)')' Pseudopotential format is: paw',pspversion
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')
 write(message,'(2(a,i3),a,64i4)') &
& ' basis_size (lnmax)=',pawtab%basis_size,' (lmn_size=',&
& pawtab%lmn_size,'), orbitals=',orbitals(1:pawtab%basis_size)
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')
 write(message,'(a,f11.8)')' Spheres core radius: rc_sph=',pawtab%rpaw
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')
 write(message,'(a,i1,a)')' ',nmesh,' radial meshes are used:'
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')
 do imsh=1,nmesh
  if (radmesh(imsh)%mesh_type==1) &
&  write(message,'(a,i1,a,i4,a,g12.5)') &
&  '  - mesh ',imsh,': r(i)=step*(i-1), size=',radmesh(imsh)%mesh_size,&
&  ' , step=',radmesh(imsh)%rstep
  if (radmesh(imsh)%mesh_type==2) &
&  write(message,'(a,i1,a,i4,2(a,g12.5))') &
&  '  - mesh ',imsh,': r(i)=AA*[exp(BB*(i-1))-1], size=',radmesh(imsh)%mesh_size,&
&  ' , AA=',radmesh(imsh)%rstep,' BB=',radmesh(imsh)%lstep
  if (radmesh(imsh)%mesh_type==3) &
&  write(message,'(a,i1,a,i4,2(a,g12.5))') &
&  '  - mesh ',imsh,': r(i)=AA*exp(BB*(i-2)), size=',radmesh(imsh)%mesh_size,&
&  ' , AA=',radmesh(imsh)%rstep,' BB=',radmesh(imsh)%lstep
  if (radmesh(imsh)%mesh_type==4) &
&  write(message,'(a,i1,a,i4,a,g12.5)') &
&  '  - mesh ',imsh,': r(i)=-AA*ln(1-(i-1)/n), n=size=',radmesh(imsh)%mesh_size,&
&  ' , AA=',radmesh(imsh)%rstep
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
 end do
 if (pawtab%shape_type==-1) then
  write(message,'(a)')&
  ' Shapefunction is NUMERIC type: directly read from atomic data file'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
 end if
 if (pawtab%shape_type==1) then
  write(message,'(2a,a,f6.3,a,i3)')&
&  ' Shapefunction is EXP type: shapef(r)=exp(-(r/sigma)**lambda)',ch10,&
&  '                            with sigma=',pawtab%shape_sigma,' and lambda=',pawtab%shape_lambda
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
 end if
 if (pawtab%shape_type==2) then
  write(message,'(a)')&
  ' Shapefunction is SIN type: shapef(r)=[sin(pi*r/rshp)/(pi*r/rshp)]**2'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
 end if
 if (pawtab%shape_type==3) then
  write(message,'(a)')&
&  ' Shapefunction is BESSEL type: shapef(r,l)=aa(1,l)*jl(q(1,l)*r)+aa(2,l)*jl(q(2,l)*r)'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
 end if
 if (pawtab%rshp<1.d-8) then
  write(message,'(a)') ' Radius for shape functions = sphere core radius'
 else
  write(message,'(a,f11.8)') ' Radius for shape functions = ',pawtab%rshp
 end if
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')

!==========================================================
!Perfom tests

!Spin-orbit not yet implemented
 if (pspso>1) then
  write(message, '(a,a,a,a,a,a)' ) ch10,&
&  ' psp7in: ERROR -',ch10,&
&  '  Spin-orbit not yet implemented for PAW psps !',ch10,&
&  '  Action : check your pseudopotential or input file.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!Are lmax and orbitals compatibles ?
 if (lmax/=maxval(orbitals)) then
  write(message, '(a,a,a,a,a,a)' ) ch10,&
&  ' psp7in: ERROR -',ch10,&
&  '  lmax /= MAX(orbitals) !',ch10,&
&  '  Action: check your pseudopotential file.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!Are lmnmax and pawtab%lmn_size compatibles ?
 if (lmnmax<pawtab%lmn_size) then
  write(message, '(a,a,a,a)' ) ch10,&
&  ' psp7in: BUG -',ch10,&
&  '  lmnmax < lmn size (read from pseudo) !'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!Only mesh_type=1,2, 3 or 4 allowed
 do imsh=1,nmesh
  if (radmesh(imsh)%mesh_type>4) then
   write(message, '(a,a,a,a,a,a)' ) ch10,&
&   ' psp7in: ERROR -',ch10,&
&   '  Only mesh types 1,2, 3 or 4 allowed !',ch10,&
&   '  Action : check your pseudopotential or input file.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
   call leave_new('COLL')
  end if
 end do

!==========================================================
!Initialize radial meshes

 do imsh=1,nmesh
  allocate(radmesh(imsh)%rad (radmesh(imsh)%mesh_size))
  allocate(radmesh(imsh)%radfact(radmesh(imsh)%mesh_size))
  allocate(radmesh(imsh)%simfact(radmesh(imsh)%mesh_size))
  call compmesh(radmesh(imsh),-1._dp)
 end do

!==========================================================
!Read tabulated atomic data

!---------------------------------
!Read wave-functions (phi)
 do ib=1,pawtab%basis_size
  read (tmp_unit,*)
  if (pspversion==1) iread1=1
  if (pspversion>1) read (tmp_unit,*) iread1
  if (ib==1) then
   pawrad%mesh_type=radmesh(iread1)%mesh_type
   pawrad%mesh_size=radmesh(iread1)%mesh_size
   pawrad%rstep=radmesh(iread1)%rstep
   pawrad%lstep=radmesh(iread1)%lstep
   call compmesh(pawrad,pawtab%rpaw)
   imainmesh=iread1
  else if (iread1/=imainmesh) then
   write(message, '(a,a,a,a,a,a)' ) ch10,&
&   ' psp7in: ERROR -',ch10,&
&   '  All Phi and tPhi must be given on the same radial mesh !',ch10,&
&   '  Action: check your pseudopotential file.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
   call leave_new('COLL')
  end if
  read (tmp_unit,*) (pawtab%phi(ir,ib),ir=1,pawrad%mesh_size)
 end do

!---------------------------------
!Read pseudo wave-functions (tphi)
 do ib=1,pawtab%basis_size
  read (tmp_unit,*)
  if (pspversion==1) iread1=1
  if (pspversion>1) read (tmp_unit,*) iread1
  if (iread1/=imainmesh) then
   write(message, '(a,a,a,a,a,a)' ) ch10,&
&   ' psp7in: ERROR -',ch10,&
&   '  All Phi and tPhi must be given on the same radial mesh !',ch10,&
&   '  Action: check your pseudopotential file.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
   call leave_new('COLL')
  end if
  read (tmp_unit,*) (pawtab%tphi(ir,ib),ir=1,pawrad%mesh_size)
 end do
 write(message,'(a,i1)') &
& ' Radial grid used for partial waves is grid ',imainmesh
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')

!---------------------------------
!Read projectors (tproj)
 do ib=1,pawtab%basis_size
  read (tmp_unit,*)
  if (pspversion==1) iread1=2
  if (pspversion>1) read (tmp_unit,*) iread1
  if (ib==1) then
   call copymesh(radmesh(iread1),tproj_mesh)
   iprojmesh=iread1
   allocate(tproj(tproj_mesh%mesh_size,pawtab%basis_size))
  else if (iread1/=iprojmesh) then
   write(message, '(a,a,a,a,a,a)' ) ch10,&
&   ' psp7in: ERROR -',ch10,&
&   '  All tprojectors must be given on the same radial mesh !',ch10,&
&   '  Action: check your pseudopotential file.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
   call leave_new('COLL')
  end if
  read (tmp_unit,*) (tproj(ir,ib),ir=1,tproj_mesh%mesh_size)
 end do
 write(message,'(a,i1)') &
& ' Radial grid used for projectors is grid ',iprojmesh
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')

!---------------------------------
!Read core density (coredens)
 read (tmp_unit,*)
 if (pspversion==1) iread1=1
 if (pspversion>1) read (tmp_unit,*) iread1
 icoremesh=iread1
 call copymesh(radmesh(iread1),core_mesh)
 if ((radmesh(icoremesh)%mesh_type/=pawrad%mesh_type).or.&
& (radmesh(icoremesh)%rstep    /=pawrad%rstep)    .or.&
& (radmesh(icoremesh)%lstep    /=pawrad%lstep)) then
  write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&  ' psp7in: ERROR -',ch10,&
&  '  Ncore must be given on a radial mesh with the same',ch10,&
&  '  type and step(s) than the main radial mesh (mesh for Phi) !',ch10,&
&  '  Action: check your pseudopotential file.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if
 allocate(ncore(core_mesh%mesh_size))
 read (tmp_unit,*) (ncore(ir),ir=1,core_mesh%mesh_size)
 pawtab%coredens(1:pawrad%mesh_size)=ncore(1:pawrad%mesh_size)

!---------------------------------
!Read pseudo core density (tcoredens)
 read (tmp_unit,*)
 if (pspversion==1) iread1=1
 if (pspversion>1) read (tmp_unit,*) iread1
 if (iread1/=icoremesh) then
  write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&  ' psp7in: ERROR -',ch10,&
&  '  Pseudized core density (tNcore) must be given',ch10,&
&  '  on the same radial mesh as core density (Ncore) !',ch10,&
&  '  Action: check your pseudopotential file.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if
 allocate(tncore(core_mesh%mesh_size))
 read (tmp_unit,*) (tncore(ir),ir=1,core_mesh%mesh_size)
 if (maxval(abs(tncore(:)))<tol6) then
  pawtab%usetcore=0
  pawtab%tcoredens(1:pawrad%mesh_size)=zero
 else
  pawtab%usetcore=1
  pawtab%tcoredens(1:pawrad%mesh_size)=tncore(1:pawrad%mesh_size)
 end if
 write(message,'(a,i1)') &
& ' Radial grid used for (t)core density is grid ',icoremesh
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')

!---------------------------------
!Read frozen part of Dij terms (dij0)
 read (tmp_unit,*)
 read (tmp_unit,*) (pawtab%dij0(ib),ib=1,pawtab%lmn2_size)

!---------------------------------
!Read initial guess of rhoij (rhoij0)
 read (tmp_unit,*)
 read (tmp_unit,*) (pawtab%rhoij0(ib),ib=1,pawtab%lmn2_size)

!---------------------------------
!Read local pseudopotential=Vh(tn_zc) or Vbare
 read (tmp_unit,*)
 if (pspversion==1) ivlocmesh=3
 pawtab%vlocopt=1
 if (pspversion==2) then
  read (tmp_unit,*) ivlocmesh
 else if (pspversion>2) then
! read (tmp_unit,fmt=*,err=30,end=30) ivlocmesh,pawtab%vlocopt
  message=blank
  read (tmp_unit,fmt=*) message
  read (message,fmt=*) ivlocmesh
  write(numb,'(i1)')ivlocmesh
  ii=index(message,numb)
  submessage=trim(message(ii+1:))
  if(len_trim(submessage)/=0)then
   do ii=1,len_trim(submessage)
    numb=submessage(ii:ii)
    if(numb==blank)cycle
    jj=index('123456789',numb)
    if(jj<1 .or. jj>9)exit
    pawtab%vlocopt=jj
   end do
  end if
 end if
 call copymesh(radmesh(ivlocmesh),vloc_mesh)
 allocate(vlocr(vloc_mesh%mesh_size))
 read (tmp_unit,*) (vlocr(ir),ir=1,vloc_mesh%mesh_size)
 write(message,'(a,i1)') &
& ' Radial grid used for Vloc is grid ',ivlocmesh
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')

!---------------------------------
!Eventually read "numeric" shapefunctions (if shape_type=-1)
 if (pawtab%shape_type==-1) then
  do il=1,pawtab%l_size
   read (tmp_unit,*)
   if (pspversion==1) iread1=1
   if (pspversion>1) read (tmp_unit,*) iread1
   if (il==1) then
    call copymesh(radmesh(iread1),shpf_mesh)
    ishpfmesh=iread1
    allocate(shpf(shpf_mesh%mesh_size,pawtab%l_size))
   else if (iread1/=ishpfmesh) then
    write(message, '(a,a,a,a,a,a)' ) ch10,&
&    ' psp7in: ERROR -',ch10,&
&    '  All shape functions must be given on the same radial mesh !',ch10,&
&    '  Action: check your pseudopotential file.'
    call wrtout(ab_out,message,'COLL')
    call wrtout(06,  message,'COLL')
    call leave_new('COLL')
   end if
   read (tmp_unit,*) (shpf(ir,il),ir=1,shpf_mesh%mesh_size)
  end do
  write(message,'(a,i1)') &
&  ' Radial grid used for shape functions is grid ',iread1
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')

! Has to spline shape functions if mesh is not the "main" mesh
  if (ishpfmesh/=imainmesh) then
   msz=shpf_mesh%mesh_size
   allocate(work1(msz),work2(msz),work3(msz),work4(pawrad%mesh_size))
   work3(:)=shpf_mesh%rad(:)
   work4(:)=pawrad%rad(:)
   do il=1,pawtab%l_size
    call bound_deriv(shpf(1:msz,il),shpf_mesh,msz,yp1,ypn)
    call spline(work3,shpf(:,il),msz,yp1,ypn,work1,work2)
    call splint(msz,work3,shpf(:,il),work1,pawrad%mesh_size,work4,pawtab%shapefunc(:,il))
   end do
   deallocate(work1,work2,work3,work4)
  else
   pawtab%shapefunc(:,:)=shpf(:,:)
  end if
  deallocate(shpf)
 end if

!---------------------------------
!Read pseudo valence density (if psp version >=4)
 if (pspversion>=4) then
  read (tmp_unit,*)
  read (tmp_unit,*) iread1
  ivalemesh=iread1
  call copymesh(radmesh(iread1),vale_mesh)
  allocate(tnvale(vale_mesh%mesh_size))
  read (tmp_unit,*) (tnvale(ir),ir=1,vale_mesh%mesh_size)
  pawtab%usetvale=1
  write(message,'(a,i1)') &
&  ' Radial grid used for pseudo valence density is grid ',ivalemesh
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
 else
  pawtab%usetvale=0
 end if

!==========================================================
!Perfom tests on meshes

!Are radial meshes for Phi and Vloc compatibles ?
 if (vloc_mesh%rmax<pawrad%rmax) then
  write(message, '(a,a,a,a,a,a)' ) ch10,&
&  ' psp7in: ERROR -',ch10,&
&  '  Rmax for Vloc < Rmax for Phi !',ch10,&
&  '  Action : check your pseudopotential (increase Vloc meshSize).'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!Are mmax and mesh_size for partial waves compatibles ?
 if (mmax/=pawrad%mesh_size) then
  write(message, '(a,a,a,a,a,a)' ) ch10,&
&  ' psp7in: ERROR -',ch10,&
&  '  mmax /= phi_mesh_size in psp file !',ch10,&
&  '  Action: check your pseudopotential file.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!Are radial meshes for (t)Ncore and Phi compatibles ?
 if (core_mesh%mesh_size<pawrad%mesh_size) then
  write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&  ' psp7in: ERROR -',ch10,&
&  '  Mesh size for core density must be equal or larger',ch10,&
&  '  than mesh size for spheres (partial waves) !',ch10,&
&  '  Action : check your pseudopotential (increase Ncore meshSize).'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!Are radial meshes for (t)Nvale and Phi compatibles ?
 if ((pawtab%usetvale==1).and.(vale_mesh%mesh_size<pawrad%mesh_size)) then
  write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&  ' psp7in: ERROR -',ch10,&
&  '  Mesh size for pseudo valence density must be equal or larger',ch10,&
&  '  than mesh size for spheres (partial waves) !',ch10,&
&  '  Action : check your pseudopotential (increase tNvale meshSize).'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!Is PAW radius included inside radial mesh ?
 if (pawtab%rpaw>pawrad%rmax+tol8) then
  write(message, '(a,a,a,a,a,a)' ) ch10,&
&  ' psp7in: ERROR -',ch10,&
&  '  Radius of PAW sphere is outside the radial mesh !',ch10,&
&  '  Action: check your pseudopotential file.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!Max. radius of mesh for Vloc has to be "small" in order to avoid numeric noise ?
 if (vloc_mesh%rmax>rm_vloc) then
  vloc_mesh%mesh_size=ifromr(vloc_mesh,rm_vloc)
  vloc_mesh%rmax=vloc_mesh%rad(vloc_mesh%mesh_size)
  write(message, '(a,a,a,a,f6.2,a,a,a,a,a,i4,a)' ) ch10,&
&  ' psp7in: WARNING -',ch10,&
&  '  Max. radius for Vloc was too large (>',rm_vloc,' a.u.) !',ch10,&
&  '  Numeric noise was possible.',ch10,&
&  '  Mesh size for Vloc has been set to ',vloc_mesh%mesh_size,'.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
 end if

!This test has been disable... MT 2006-25-10
!For Simpson rule, it is better to have odd mesh sizes
!itest=0
!do imsh=1,nmesh
!if (mod(radmesh(imsh)%mesh_size,2)==0.and.radmesh(imsh)%mesh_type==1) itest=1
!end do
!if (itest==1) then
!write(message, '(8a)' ) ch10,&
!&   ' psp7in: WARNING -',ch10,&
!&   '  Regular radial meshes should have odd number of points ',ch10,&
!&   '  for better accuracy of integration sheme (Simpson rule).',ch10,&
!&   '  Althought it''s not compulsory, you should change mesh sizes in psp file.'
!call wrtout(ab_out,message,'COLL')
!call wrtout(06,  message,'COLL')
!end if

!Test the compatibilty between Rpaw and mesh for (t)Phi
 if (pspversion>=3) then
  itest=ifromr(radmesh(imainmesh),pawtab%rpaw)
  if (itest+2>radmesh(imainmesh)%mesh_size) then
   write(message, '(12a)' ) ch10,&
&   ' psp7in: WARNING -',ch10,&
&   '  Atomic data could produce inaccurate results:',ch10,&
&   '    Wavefunctions and pseudo-wavefunctions should',ch10,&
&   '    be given on a radial mesh larger than the PAW',ch10,&
&   '    spheres (at least 2 additional points) !',ch10,&
&   '  Action: check your pseudopotential file.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
  end if
  if (abs(pawtab%rpaw-radmesh(imainmesh)%rad(itest))<tol8) itest=itest-1
  ib=0;isnotzero=0
  do while ((isnotzero==0).and.(ib<pawtab%basis_size))
   ib=ib+1;ir=itest
   do while ((isnotzero==0).and.(ir<radmesh(imainmesh)%mesh_size))
    ir=ir+1;if (abs(pawtab%phi(ir,ib)-pawtab%tphi(ir,ib))>tol8) isnotzero=1
   end do
  end do
  if (isnotzero>0) then
   write(message, '(10a)' ) ch10,&
&   ' psp7in: ERROR -',ch10,&
&   '  Atomic data are inconsistent:',ch10,&
&   '  For r>=r_paw, pseudo wavefunctions are not',ch10,&
&   '  equal to wave functions (Phi(r)/=tPhi(r)) !',ch10,&
&   '  Action: check your pseudopotential file.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
   call leave_new('COLL')
  end if
 else
! For compatibility reasons set PAW radius at the end of mesh (older versions)
  if (pawtab%rpaw/=pawrad%rmax) then
   pawtab%rpaw=pawrad%rmax
   call compmesh(pawrad,-1._dp)
  end if
 end if
!If Vloc is a "Vbare" potential, it has to be localized inside PAW spheres
 if (pawtab%vlocopt==0.and.(vloc_mesh%rmax>pawtab%rpaw+tol10)) then
  write(message, '(8a)' ) ch10,&
&  ' psp7in: ERROR -',ch10,&
&  '  Atomic data are inconsistent:',ch10,&
&  '  Local potential is a "Vbare" potential',ch10,&
&  '  and is not localized inside PAW sphere !'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!==========================================================
!Initialize various parameters in structured data

!Some dims
 pawtab%l_size=2*maxval(orbitals)+1
 pawtab%mesh_size=pawrad%mesh_size
 pawtab%lmn2_size=pawtab%lmn_size*(pawtab%lmn_size+1)/2

!indlmn calculation (indices for (l,m,n) basis)
 allocate(nprj(0:maxval(orbitals)))
 ilmn=0;iln=0;nprj=0
 do ib=1,pawtab%basis_size
  il=orbitals(ib)
  nprj(il)=nprj(il)+1
  iln=iln+1
  do ilm=1,2*il+1
   indlmn(1,ilmn+ilm)=il
   indlmn(2,ilmn+ilm)=ilm-(il+1)
   indlmn(3,ilmn+ilm)=nprj(il)
   indlmn(4,ilmn+ilm)=il*il+ilm
   indlmn(5,ilmn+ilm)=iln
   indlmn(6,ilmn+ilm)=1
  end do
  ilmn=ilmn+2*il+1
 end do
 deallocate(nprj,orbitals)

!Are ilmn (found here) and pawtab%lmn_size compatibles ?
 if (ilmn/=pawtab%lmn_size) then
  write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&  ' psp7in: ERROR -',ch10,&
&  '  Calculated lmn size differs from',ch10,&
&  '  lmn_size read from pseudo !',ch10,&
&  ' Action: check your pseudopotential file.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!==========================================================
!Compute ffspl(q) (and derivatives)

 ffspl=zero
 call psp7nl(ffspl,indlmn,pawtab%lmn_size,lnmax,mqgrid_ff,qgrid_ff,tproj_mesh,tproj)
 deallocate(tproj)

!==========================================================
!Compute eventually compensation charge radius (i.e. radius for shape functions)

 if (pawtab%shape_type>0.and.pawtab%rshp<1.d-8) then
  pawtab%rshp=pawtab%rpaw
 else if (pawtab%shape_type==-1) then
  ir=ifromr(radmesh(imainmesh),pawtab%rpaw)+1;isnotzero=0
  do while ((isnotzero==0).and.(ir>1))
   ir=ir-1;il=0
   do while ((isnotzero==0).and.(il<pawtab%l_size))
    il=il+1;if (pawtab%shapefunc(ir,il)>tol16) isnotzero=1
   end do
  end do
  ir=min(ir+1,radmesh(imainmesh)%mesh_size)
  pawtab%rshp=radmesh(imainmesh)%rad(ir)
  do il=1,pawtab%l_size
   if (pawtab%shapefunc(ir,il)>tol6) then
    write(message, '(a,a,a,a,a,a)' ) ch10,&
&    ' psp7in: ERROR -',ch10,&
&    '  Shape function is not zero at PAW radius !',ch10,&
&    '  Action: check your pseudopotential file.'
    call wrtout(ab_out,message,'COLL')
    call wrtout(06,  message,'COLL')
    call leave_new('COLL')
   end if
  end do
 end if

!==========================================================
!If read Vloc potential is in "Vbare" format,
!translate it into VH(tnzc) format
!VXC difference is missing !!!
 if (pawtab%vlocopt==0) then
  write(message,'(a)') ' Local potential is in "Vbare" format... '
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
! Compute r2.k(r)
  msz=core_mesh%mesh_size;allocate(r2k(msz))
  call pawshpfun(0,core_mesh,intg,pawtab,r2k)
  r2k(1:msz)=r2k(1:msz)*core_mesh%rad(1:msz)**2
! Compute VH[4pi.r2.n(r)=4pi.r2.tncore(r)+(Qcore-Z).r2.k(r)]
  allocate(nwk(core_mesh%mesh_size),vh(core_mesh%mesh_size))
  nwk(:)=tncore(:)*four_pi*core_mesh%rad(:)**2
  call simp_gen(qcore,nwk,core_mesh)
  nwk(1:msz)=nwk(1:msz)-r2k(1:msz)*(qcore+zion)
  call poisson(nwk,0,intg,core_mesh,vh)
  vh(2:msz)=vh(2:msz)/core_mesh%rad(2:msz)
  call deducer0(vh,msz,core_mesh)
  deallocate(nwk)
! Eventually spline Vbare
  allocate(vbare(core_mesh%mesh_size))
  if ((core_mesh%mesh_type/=vloc_mesh%mesh_type).or.&
&  (core_mesh%rstep    /=vloc_mesh%rstep)    .or.&
&  (core_mesh%lstep    /=vloc_mesh%lstep)) then
   msz=core_mesh%mesh_size;if (vloc_mesh%rmax<core_mesh%rmax) msz=ifromr(core_mesh,vloc_mesh%rmax)
   call bound_deriv(vlocr(1:vloc_mesh%mesh_size),vloc_mesh,vloc_mesh%mesh_size,yp1,ypn)
   allocate(work1(vloc_mesh%mesh_size),work2(vloc_mesh%mesh_size))
   call spline(vloc_mesh%rad,vlocr,vloc_mesh%mesh_size,yp1,ypn,work1,work2)
   call splint(vloc_mesh%mesh_size,vloc_mesh%rad,vlocr,work1,msz,core_mesh%rad(1:msz),vbare)
   deallocate(work1,work2)
  else
   msz=min(core_mesh%mesh_size,vloc_mesh%mesh_size)
   vbare(1:msz)=vlocr(1:msz)
  end if
! Build VH(tnzc) from Vbare
  deallocate(vlocr,vloc_mesh%rad,vloc_mesh%radfact,vloc_mesh%simfact)
  call copymesh(core_mesh,vloc_mesh)
  allocate(vlocr(core_mesh%mesh_size))
  vlocr(:)=vbare(:)+vh(:)
  deallocate(vbare,vh)
! Compute <tPhi_i|VH(tnzc)|tPhi_j> and int[VH(tnzc)*Qijhat(r)dr] parts of Dij0
! Note: it is possible as core_mesh and radmesh(imainmesh) have the same steps
  msz=radmesh(imainmesh)%mesh_size;allocate(work1(msz))
  work1(1:msz)=vlocr(1:msz)*r2k(1:msz)
  call simp_gen(intvh,work1,radmesh(imainmesh))
  do jlmn=1,pawtab%lmn_size
   j0lmn=jlmn*(jlmn-1)/2;jlm=indlmn(4,jlmn);jln=indlmn(5,jlmn)
   do ilmn=1,jlmn
    klmn=j0lmn+ilmn;ilm=indlmn(4,ilmn);iln=indlmn(5,ilmn)
    if (jlm==ilm) then
     work1(1:msz)=pawtab%tphi(1:msz,iln)*pawtab%tphi(1:msz,jln)*(vlocr(1:msz)-intvh) &
&     -pawtab%phi (1:msz,iln)*pawtab%phi (1:msz,jln)*intvh
     call simp_gen(intg,work1,radmesh(imainmesh))
     pawtab%dij0(klmn)=pawtab%dij0(klmn)+intg
    end if
   end do
  end do
  deallocate(work1,r2k)
 end if

!==========================================================
!Try to optimize CPU time:
!If Vloc mesh size is big, spline Vloc into a smaller log. mesh

 reduced_vloc=(vloc_mesh%mesh_size>int(reduced_mshsz))
 if (reduced_vloc) then
  msz=vloc_mesh%mesh_size
  rvloc_mesh%mesh_type=3
  rvloc_mesh%mesh_size=reduced_mshsz
  rvloc_mesh%rstep    =reduced_rstep
  rvloc_mesh%lstep    =log(0.9999999_dp*vloc_mesh%rmax/reduced_rstep)/dble(reduced_mshsz-2)
  allocate(rvloc_mesh%rad (reduced_mshsz))
  allocate(rvloc_mesh%radfact(reduced_mshsz))
  allocate(rvloc_mesh%simfact(reduced_mshsz))
  allocate(rvlocr(reduced_mshsz))
  call compmesh(rvloc_mesh,-1._dp)
  call bound_deriv(vlocr(1:msz),vloc_mesh,msz,yp1,ypn)
  allocate(work1(msz),work2(msz),work3(msz))
  work3(:)=vloc_mesh%rad(:)
  call spline(work3,vlocr,msz,yp1,ypn,work1,work2)
  call splint(msz,work3,vlocr,work1,reduced_mshsz,rvloc_mesh%rad,rvlocr)
  deallocate(work1,work2,work3)
 end if

!==========================================================
!Try to optimize CPU time:
!If ncore mesh size is big, spline tncore into a smaller log. mesh

 reduced_ncor=(core_mesh%mesh_size>int(reduced_mshsz)).and.(pawtab%usetcore/=0)
 if (reduced_ncor) then
  msz=core_mesh%mesh_size
  rcore_mesh%mesh_type=3
  rcore_mesh%mesh_size=reduced_mshsz
  rcore_mesh%rstep    =reduced_rstep
  rcore_mesh%lstep    =log(0.9999999_dp*core_mesh%rmax/reduced_rstep)/dble(reduced_mshsz-2)
  allocate(rcore_mesh%rad (reduced_mshsz))
  allocate(rcore_mesh%radfact(reduced_mshsz))
  allocate(rcore_mesh%simfact(reduced_mshsz))
  allocate(rtncor(reduced_mshsz))
  call compmesh(rcore_mesh,-1._dp)
  call bound_deriv(tncore(1:msz),core_mesh,msz,yp1,ypn)
  allocate(work1(msz),work2(msz),work3(msz))
  work3(:)=core_mesh%rad(:)
  call spline(work3,tncore,msz,yp1,ypn,work1,work2)
  call splint(msz,work3,tncore,work1,reduced_mshsz,rcore_mesh%rad,rtncor)
  deallocate(work1,work2,work3)
 end if

!==========================================================
!Try to optimize CPU time:
!If vale mesh size is big, spline tnvale into a smaller log. mesh

 if (pawtab%usetvale==1) then
  reduced_nval=(vale_mesh%mesh_size>int(reduced_mshsz))
  if (reduced_nval) then
   msz=vale_mesh%mesh_size
   rvale_mesh%mesh_type=3
   rvale_mesh%mesh_size=reduced_mshsz
   rvale_mesh%rstep    =reduced_rstep
   rvale_mesh%lstep    =log(0.9999999_dp*vale_mesh%rmax/reduced_rstep)/dble(reduced_mshsz-2)
   allocate(rvale_mesh%rad (reduced_mshsz))
   allocate(rvale_mesh%radfact(reduced_mshsz))
   allocate(rvale_mesh%simfact(reduced_mshsz))
   allocate(rtnval(reduced_mshsz))
   call compmesh(rvale_mesh,-1._dp)
   call bound_deriv(tnvale(1:msz),vale_mesh,msz,yp1,ypn)
   allocate(work1(msz),work2(msz),work3(msz))
   work3(:)=vale_mesh%rad(:)
   call spline(work3,tnvale,msz,yp1,ypn,work1,work2)
   call splint(msz,work3,tnvale,work1,reduced_mshsz,rvale_mesh%rad,rtnval)
   deallocate(work1,work2,work3)
  end if
 else
  reduced_nval=.false.
 end if

!==========================================================
!Compute Vlspl(q) (and second derivative) from Vloc(r)

!Compute Vlspl(q)=q^2.Vloc(q) from vloc(r)
 if (reduced_vloc) then
  call psp7lo(epsatm,mqgrid_vl,qgrid_vl,vlspl(:,1),rvloc_mesh,rvlocr,yp1,ypn,zion)
 else
  call psp7lo(epsatm,mqgrid_vl,qgrid_vl,vlspl(:,1),vloc_mesh,vlocr,yp1,ypn,zion)
 end if
!Compute second derivative of Vlspl(q)
 allocate(work1(mqgrid_vl))
 call spline(qgrid_vl,vlspl(:,1),mqgrid_vl,yp1,ypn,vlspl(:,2),work1)
 deallocate(work1)

!==========================================================
!Compute tcorespl(q) (and second derivative) from tNcore(r)

 pawtab%mqgrid=mqgrid_vl
 xcccrc=core_mesh%rmax

 if (pawtab%usetcore/=0) then
! Compute tcorespl(q)=tNc(q) from tNcore(r)
  if (reduced_ncor) then
   call psp7cg(pawtab%dncdq0,mqgrid_vl,qgrid_vl,pawtab%tcorespl(:,1),rcore_mesh,rtncor,yp1,ypn)
  else
   call psp7cg(pawtab%dncdq0,mqgrid_vl,qgrid_vl,pawtab%tcorespl(:,1),core_mesh,tncore,yp1,ypn)
  end if
! Compute second derivative of tcorespl(q)
  allocate(work1(mqgrid_vl))
  call spline(qgrid_vl,pawtab%tcorespl(:,1),mqgrid_vl,yp1,ypn,pawtab%tcorespl(:,2),work1)
  deallocate (work1)
 else
  pawtab%tcorespl=zero
  pawtab%dncdq0=zero
 end if

!==========================================================
!Compute tvalespl(q) (and second derivative) from tNvale(r)

 if (pawtab%usetvale/=0) then

! Add compensation density to pseudo valence density
  if (reduced_nval) then
   msz=rvale_mesh%mesh_size;allocate(nwk(msz))
   nwk(1:msz)=rtnval(1:msz)*rvale_mesh%rad(1:msz)**2
   call simp_gen(intg,nwk,rvale_mesh);qq=zion/four_pi-intg
   call pawshpfun(0,rvale_mesh,intg,pawtab,nwk)
   rtnval(1:msz)=rtnval(1:msz)+qq*nwk(1:msz)
   deallocate(nwk)
  else
   msz=vale_mesh%mesh_size;allocate(nwk(msz))
   nwk(1:msz)=tnvale(1:msz)*vale_mesh%rad(1:msz)**2
   call simp_gen(intg,nwk,vale_mesh);qq=zion/four_pi-intg
   call pawshpfun(0,vale_mesh,intg,pawtab,nwk)
   tnvale(1:msz)=tnvale(1:msz)+qq*nwk(1:msz)
   deallocate(nwk)
  end if
! Compute tvalespl(q)=tNv(q) from tNvale(r)
  if (reduced_nval) then
   call psp7cg(pawtab%dnvdq0,mqgrid_vl,qgrid_vl,pawtab%tvalespl(:,1),rvale_mesh,rtnval,yp1,ypn)
  else
   call psp7cg(pawtab%dnvdq0,mqgrid_vl,qgrid_vl,pawtab%tvalespl(:,1),vale_mesh,tnvale,yp1,ypn)
  end if
! Compute second derivative of tvalespl(q)
  allocate(work1(mqgrid_vl))
  call spline(qgrid_vl,pawtab%tvalespl(:,1),mqgrid_vl,yp1,ypn,pawtab%tvalespl(:,2),work1)
  deallocate (work1)
 else
  pawtab%dnvdq0=zero
 end if

!==================================================
!Compute Ex-correlation energy for the core density

 allocate(work1(core_mesh%mesh_size),work2(core_mesh%mesh_size))
 allocate(work3(core_mesh%mesh_size))
 work1(:)=zero;tmp_lmselect(1)=.true.
 if (pawxcdev/=0) then
  call pawxcm(ncore,pawtab%exccore,yp1,0,ixc,1,tmp_lmselect,work3,1,4,&
&  pawang_tmp,core_mesh,pawxcdev,work1,1,0,work2,xclevel)
 else
  pawang_tmp%l_size_max=1;pawang_tmp%angl_size=1
  allocate(pawang_tmp%angwgth(1));pawang_tmp%angwgth(1)=1._dp
  allocate(pawang_tmp%ylmr(1,1));pawang_tmp%ylmr(1,1)=1._dp/sqrt(four_pi)
  call pawxc (ncore,pawtab%exccore,yp1,ixc,1,tmp_lmselect,work3,1,4,&
&  pawang_tmp,core_mesh,work1,1,0,work2,xclevel)
  deallocate(pawang_tmp%angwgth,pawang_tmp%ylmr)
 end if
 deallocate(work1,work2,work3)

!==========================================================
!Free temporary allocated space

 do imsh=1,nmesh
  deallocate(radmesh(imsh)%rad)
  deallocate(radmesh(imsh)%radfact)
  deallocate(radmesh(imsh)%simfact)
 end do
 deallocate(radmesh)
 deallocate(vlocr,ncore,tncore)
 deallocate(tproj_mesh%rad,tproj_mesh%radfact,tproj_mesh%simfact)
 deallocate(core_mesh%rad,core_mesh%radfact,core_mesh%simfact)
 deallocate(vloc_mesh%rad,vloc_mesh%radfact,vloc_mesh%simfact)
 if(pawtab%shape_type==-1)deallocate(shpf_mesh%rad, shpf_mesh%radfact,shpf_mesh%simfact)
 if (reduced_vloc) deallocate(rvloc_mesh%rad,rvloc_mesh%radfact,rvloc_mesh%simfact,rvlocr)
 if (reduced_ncor) deallocate(rcore_mesh%rad,rcore_mesh%radfact,rcore_mesh%simfact,rtncor)
 if (reduced_nval) deallocate(rvale_mesh%rad,rvale_mesh%radfact,rvale_mesh%simfact,rtnval)
 if (pspversion>=4) deallocate(tnvale,vale_mesh%rad,vale_mesh%radfact,vale_mesh%simfact)

!DEBUG
!write(6,*)' psp7in : exit'
!ENDDEBUG

end subroutine psp7in
!!***


!!****f* ABINIT/bound_deriv
!! NAME
!! bound_deriv
!!
!! FUNCTION
!! Computes derivatives of a function a boundaries of interval (first and last derivative)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  func(n)= array containing function
!!  mesh <type(pawrad_type)>= radial mesh and related data
!!  nn= size of intervall
!!
!! OUTPUT
!!  yp1,ypn= derivatives of func at r(1) and r(n)
!!
!! PARENTS
!!      psp7in
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine bound_deriv(func,mesh,nn,yp1,ypn)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments----------------------
 integer, intent(in) :: nn
 real(dp), intent(in) :: func(nn)
 real(dp), intent(out) :: yp1,ypn
 type(pawrad_type),intent(in) :: mesh

!*************************************************************************

 if (mesh%radfact(1)>zero) then
  yp1=1._dp/12._dp/mesh%stepint/mesh%radfact(1) &
&  *(-25._dp*func(1)+48._dp*func(2)-36._dp*func(3)+16._dp*func(4)-3._dp*func(5))
 else
  yp1=(func(2)-func(1))/(mesh%rad(2)-mesh%rad(1))
 end if
 ypn=1._dp/12._dp/mesh%stepint &
& *( 3._dp*func(nn-4)-16._dp*func(nn-3)+36._dp*func(nn-2)-48._dp*func(nn-1) &
& +25._dp*func(nn))/mesh%radfact(nn)

end subroutine bound_deriv

!!***
