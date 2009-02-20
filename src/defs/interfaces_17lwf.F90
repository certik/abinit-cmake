!!****m* ABINIT/interfaces_17lwf
!! NAME
!! interfaces_17lwf
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/17lwf
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!! 
!!
!! SOURCE

module interfaces_17lwf

 implicit none

interface
 subroutine bldlwf(grdsize,iout,natom,nqpt,nwnn,qpoint,rcenter,wannvect)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: nwnn
  integer,intent(in) :: grdsize(3)
  real(dp),intent(in) :: qpoint(nqpt,3)
  real(dp),intent(in) :: rcenter(3)
  real(dp) :: wannvect(nqpt,nwnn,natom,3,2)
 end subroutine bldlwf
end interface

interface
 subroutine chkilwf(allerr,alpha,decflg,enwdmax,enwdmin,frozflg,grdsize,irwfl,lenstr,localqmode,&  
  &  mqpt,natom,nqpt,nstom,nwnn,prtvol,rcenter,string,subwdmax,subwdmin,tolomi,znucl)
  use defs_basis
  implicit none
  integer,intent(inout) :: decflg
  integer,intent(inout) :: frozflg
  integer,intent(inout) :: irwfl
  integer,intent(in) :: lenstr
  integer,intent(in) :: mqpt
  integer,intent(in) :: natom
  integer,intent(inout) :: nqpt
  integer,intent(inout) :: nstom
  integer,intent(inout) :: nwnn
  integer,intent(inout) :: prtvol
  real(dp),intent(inout) :: allerr
  real(dp),intent(inout) :: alpha
  real(dp),intent(inout) :: enwdmax
  real(dp),intent(inout) :: enwdmin
  character(len=*) :: string
  real(dp),intent(inout) :: subwdmax
  real(dp),intent(inout) :: subwdmin
  real(dp),intent(inout) :: tolomi
  integer,intent(in) :: grdsize(3)
  real(dp),intent(in) :: localqmode(3)
  real(dp),intent(in) :: rcenter(3)
  real(dp),intent(in) :: znucl(natom)
 end subroutine chkilwf
end interface

interface
 subroutine invars7w(allerr,alpha,decflg,enwdmax,enwdmin,frozflg,grdsize,ingss,irwfl,lenstr,localqmode,&  
  &  mqpt,natom,nstom,nwnn,prtvol,rcenter,string,subwdmax,subwdmin,tolomi,trialq,znucl)
  use defs_basis
  implicit none
  integer,intent(out) :: decflg
  integer,intent(out) :: frozflg
  integer,intent(out) :: irwfl
  integer,intent(in) :: lenstr
  integer,intent(in) :: mqpt
  integer,intent(in) :: natom
  integer,intent(out) :: nstom
  integer,intent(in) :: nwnn
  integer,intent(out) :: prtvol
  integer,intent(out) :: trialq
  real(dp),intent(out) :: allerr
  real(dp),intent(out) :: alpha
  real(dp),intent(out) :: enwdmax
  real(dp),intent(out) :: enwdmin
  character(len=*) :: string
  real(dp),intent(out) :: subwdmax
  real(dp),intent(out) :: subwdmin
  real(dp),intent(out) :: tolomi
  integer,intent(out) :: grdsize(3)
  real(dp),intent(out) :: ingss(nwnn,natom,3,2)
  real(dp),intent(out) :: localqmode(3)
  real(dp),intent(out) :: rcenter(3)
  real(dp),intent(out) :: znucl(natom)
 end subroutine invars7w
end interface

interface
 subroutine overlap_ph(eigvect,g_subsp,maxqsize,mmnkb,natom,nqpt,qneigh,qsize)
  use defs_basis
  implicit none
  integer,intent(in) :: maxqsize
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  real(dp),intent(in) :: eigvect(nqpt,3*natom,natom,3,2)
  integer,intent(in) :: g_subsp(nqpt,3*natom)
  real(dp),intent(out) :: mmnkb(nqpt,6,maxqsize,maxqsize,2)
  integer,intent(in) :: qneigh(nqpt,6)
  integer,intent(in) :: qsize(nqpt,3)
 end subroutine overlap_ph
end interface

interface
 subroutine readeig(acell,eigval,eigvect,ineig,natom,nqpt,qpoint,rprim,typat,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: ineig
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  real(dp),intent(out) :: acell(3)
  real(dp),intent(out) :: eigval(nqpt,3*natom)
  real(dp),intent(out) :: eigvect(nqpt,3*natom,natom,3,2)
  real(dp),intent(out) :: qpoint(nqpt,3)
  real(dp),intent(out) :: rprim(3,3)
  integer,intent(out) :: typat(natom)
  real(dp),intent(out) :: xred(3,natom)
 end subroutine readeig
end interface

interface
 subroutine rwwan(irwfl,iwf,natom,nqpt,nwnn,qpoint,wannvect)
  use defs_basis
  implicit none
  integer,intent(in) :: irwfl
  integer,intent(in) :: iwf
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: nwnn
  real(dp),intent(inout) :: qpoint(nqpt,3)
  real(dp),intent(inout) :: wannvect(nqpt,nwnn,natom,3,2)
 end subroutine rwwan
end interface

interface
 subroutine secinit(eigval,eigvect,f_subsp,g_subsp,z_subsp,&  
  &  gsize,indwnz,ingss,iqpt,lambda,maxqsize,natom,nqpt,nwnn,nwnz,qpoint,&  
  &  xred,znucl,zsize)
  use defs_basis
  implicit none
  integer,intent(in) :: gsize
  integer,intent(in) :: iqpt
  integer,intent(in) :: maxqsize
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: nwnn
  integer,intent(in) :: zsize
  real(dp),intent(in) :: eigval(nqpt,3*natom)
  real(dp),intent(in) :: eigvect(nqpt,3*natom,natom,3,2)
  integer,intent(in) :: f_subsp(nqpt,3*natom)
  integer,intent(in) :: g_subsp(nqpt,3*natom)
  integer,intent(in) :: indwnz(nwnn)
  real(dp),intent(in) :: ingss(nwnn,natom,3,2)
  real(dp),intent(out) :: lambda(nqpt,nwnn,maxqsize,2)
  integer,intent(in) :: nwnz(nqpt)
  real(dp),intent(in) :: qpoint(nqpt,3)
  real(dp),intent(in) :: xred(natom,3)
  integer,intent(in) :: z_subsp(nqpt,3*natom)
  real(dp),intent(in) :: znucl(natom)
 end subroutine secinit
end interface

interface
 subroutine shellin(acell,nqpt,qneigh,qpoint,rprim)
  use defs_basis
  implicit none
  integer,intent(in) :: nqpt
  real(dp),intent(in) :: acell(3)
  integer,intent(out) :: qneigh(nqpt,6)
  real(dp),intent(in) :: qpoint(nqpt,3)
  real(dp),intent(in) :: rprim(3,3)
 end subroutine shellin
end interface

interface
 subroutine wanvec(atmass,eigvect,g_subsp,iout,iwf,lambda,maxqsize,&  
  &  natom,nqpt,nwnn,qpoint,qsize,wannvect)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: iwf
  integer,intent(in) :: maxqsize
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: nwnn
  real(dp),intent(in) :: atmass(natom)
  real(dp),intent(inout) :: eigvect(nqpt,3*natom,natom,3,2)
  integer,intent(in) :: g_subsp(nqpt,3*natom)
  real(dp),intent(in) :: lambda(nqpt,nwnn,maxqsize,2)
  real(dp),intent(inout) :: qpoint(nqpt,3)
  integer,intent(in) :: qsize(nqpt,3)
  real(dp),intent(out) :: wannvect(nqpt,nwnn,natom,3,2)
 end subroutine wanvec
end interface

interface
 subroutine zmnbld(f_subsp,g_subsp,lambda,maxqsize,mmnkb,natom,nqpt,nwnn,nwnz,qneigh,qsize,Zmn)
  use defs_basis
  implicit none
  integer,intent(in) :: maxqsize
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: nwnn
  real(dp),intent(out) :: Zmn(nqpt,maxqsize,maxqsize,2)
  integer,intent(in) :: f_subsp(nqpt,3*natom)
  integer,intent(in) :: g_subsp(nqpt,3*natom)
  real(dp),intent(in) :: lambda(nqpt,nwnn,maxqsize,2)
  real(dp),intent(in) :: mmnkb(nqpt,6,maxqsize,maxqsize,2)
  integer,intent(in) :: nwnz(nqpt)
  integer,intent(in) :: qneigh(nqpt,6)
  integer,intent(in) :: qsize(nqpt,3)
 end subroutine zmnbld
end interface

end module interfaces_17lwf
!!***
