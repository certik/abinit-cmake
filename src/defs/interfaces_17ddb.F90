!!****m* ABINIT/interfaces_17ddb
!! NAME
!! interfaces_17ddb
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/17ddb
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

module interfaces_17ddb

 implicit none

interface
 subroutine alignph(amu,displ,d2cart,mpert,natom,ntypat,phfrq,typat)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(inout) :: displ(2,3*natom,3*natom)
  real(dp),intent(in) :: phfrq(3*natom)
  integer,intent(in) :: typat(natom)
 end subroutine alignph
end interface

interface
 subroutine asrif9(asr,atmfrc,natom,nrpt,rpt,wghatm)
  use defs_basis
  implicit none
  integer,intent(in) :: asr
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine asrif9
end interface

interface
 subroutine axial9(ifccar,vect1,vect2,vect3)
  use defs_basis
  implicit none
  real(dp),intent(in) :: ifccar(3,3)
  real(dp),intent(in) :: vect1(3)
  real(dp),intent(out) :: vect2(3)
  real(dp),intent(out) :: vect3(3)
 end subroutine axial9
end interface

interface
 subroutine bfactor(nkpt,kpt,nqpt,qpt,weight,nband,nestfactor)
  use defs_basis
  implicit none
  integer,intent(in) :: nband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nqpt
  real(dp),intent(in) :: kpt(3,nkpt)
  real(dp),intent(inout) :: nestfactor(nqpt)
  real(dp),intent(in) :: qpt(3,nqpt)
  real(dp),intent(in) :: weight(nband,nkpt)
 end subroutine bfactor
end interface

interface
 subroutine bigbx9(brav,choice,mrpt,ngqpt,nqshft,nrpt,rprim,rpt)
  use defs_basis
  implicit none
  integer,intent(in) :: brav
  integer,intent(in) :: choice
  integer,intent(in) :: mrpt
  integer,intent(in) :: nqshft
  integer,intent(out) :: nrpt
  integer,intent(in) :: ngqpt(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(out) :: rpt(3,mrpt)
 end subroutine bigbx9
end interface

interface
 subroutine canat9(brav,natom,rcan,rprim,trans,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: brav
  integer,intent(in) :: natom
  real(dp),intent(out) :: rcan(3,natom)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(out) :: trans(3,natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine canat9
end interface

interface
 subroutine canct9(acell,gprim,ib,index,irpt,natom,nrpt,&  
  &  rcan,rcart,rprim,rpt)
  use defs_basis
  implicit none
  integer,intent(out) :: ib
  integer,intent(in) :: index
  integer,intent(out) :: irpt
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(out) :: rcart(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
 end subroutine canct9
end interface

interface
 subroutine chki8(inti,intt,name)
  implicit none
  integer,intent(in) :: inti
  integer,intent(in) :: intt
  character(len=6),intent(in) :: name
 end subroutine chki8
end interface

interface
 subroutine chkin9(atifc,natifc,natom)
  implicit none
  integer,intent(in) :: natifc
  integer,intent(in) :: natom
  integer,intent(inout) :: atifc(natom)
 end subroutine chkin9
end interface

interface
 subroutine chkr8(reali,realt,name,tol)
  use defs_basis
  implicit none
  character(len=6),intent(in) :: name
  real(dp),intent(in) :: reali
  real(dp),intent(in) :: realt
  real(dp),intent(in) :: tol
 end subroutine chkr8
end interface

interface
 subroutine chkrp9(brav,rprim)
  use defs_basis
  implicit none
  integer,intent(in) :: brav
  real(dp),intent(in) :: rprim(3,3)
 end subroutine chkrp9
end interface

interface
 subroutine clean_phon_ds(phon_ds)
  use defs_elphon
  implicit none
  type(phon_type),intent(inout) :: phon_ds
 end subroutine clean_phon_ds
end interface

interface
 subroutine cmpar8 (acell,acell8,amu,amu8,dimekb,ecut,ecut8,ekb,ekb8,&  
  &  fullinit,fullinit8,iscf,iscf8,ixc,ixc8,kpt,kpt8,kptnrm,kptnr8,&  
  &  natom,natom8,nband,nband8,ngfft,ngfft8,nkpt,nkpt8,&  
  &  nsppol,nsppo8,nsym,nsym8,ntypat,ntypat8,occ,occ8,&  
  &  occopt,occop8,rprim,rprim8,sciss,sciss8,symrel,symre8,&  
  &  tnons,tnons8,tolwfr,tolwf8,typat,typat8,usepaw,wtk,wtk8,xred,xred8,zion,zion8)
  use defs_basis
  implicit none
  integer,intent(in) :: dimekb
  integer,intent(inout) :: fullinit
  integer,intent(in) :: fullinit8
  integer,intent(inout) :: iscf
  integer,intent(in) :: iscf8
  integer,intent(inout) :: ixc
  integer,intent(in) :: ixc8
  integer,intent(inout) :: natom
  integer,intent(in) :: natom8
  integer,intent(inout) :: nkpt
  integer,intent(in) :: nkpt8
  integer,intent(in) :: nsppo8
  integer,intent(inout) :: nsppol
  integer,intent(inout) :: nsym
  integer,intent(in) :: nsym8
  integer,intent(inout) :: ntypat
  integer,intent(in) :: ntypat8
  integer,intent(in) :: occop8
  integer,intent(inout) :: occopt
  integer,intent(in) :: usepaw
  real(dp),intent(inout) :: ecut
  real(dp),intent(in) :: ecut8
  real(dp),intent(in) :: kptnr8
  real(dp),intent(inout) :: kptnrm
  real(dp),intent(inout) :: sciss
  real(dp),intent(in) :: sciss8
  real(dp),intent(in) :: tolwf8
  real(dp),intent(inout) :: tolwfr
  integer,intent(inout) :: nband(*)
  integer,intent(in) :: nband8(*)
  integer,intent(inout) :: ngfft(18)
  integer,intent(in) :: ngfft8(18)
  integer,intent(in) :: symre8(3,3,*)
  integer,intent(inout) :: symrel(3,3,*)
  integer,intent(inout) :: typat(*)
  integer,intent(in) :: typat8(*)
  real(dp),intent(inout) :: acell(3)
  real(dp),intent(in) :: acell8(3)
  real(dp),intent(inout) :: amu(*)
  real(dp),intent(in) :: amu8(*)
  real(dp),intent(inout) :: ekb(dimekb,*)
  real(dp),intent(in) :: ekb8(dimekb,*)
  real(dp),intent(inout) :: kpt(3,*)
  real(dp),intent(in) :: kpt8(3,*)
  real(dp),intent(inout) :: occ(*)
  real(dp),intent(in) :: occ8(*)
  real(dp),intent(inout) :: rprim(3,3)
  real(dp),intent(in) :: rprim8(3,3)
  real(dp),intent(inout) :: tnons(3,*)
  real(dp),intent(in) :: tnons8(3,*)
  real(dp),intent(inout) :: wtk(*)
  real(dp),intent(in) :: wtk8(*)
  real(dp),intent(inout) :: xred(3,*)
  real(dp),intent(in) :: xred8(3,*)
  real(dp),intent(inout) :: zion(*)
  real(dp),intent(in) :: zion8(*)
 end subroutine cmpar8
end interface

interface
 subroutine complete_gkk(acell,elph_ds,FSfulltofull,FSkpt,gkk_flag,&  
  &  gprimd,indsym,mpert,natom,nqptirred,nsym,qptirred,qpttoqpt,rprimd,&  
  &  spqpt,symrec,symrel,tnons,ucvol,xred)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nqptirred
  integer,intent(in) :: nsym
  type(elph_type),intent(inout) :: elph_ds
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: FSfulltofull(2,nsym,elph_ds%nFSkpt)
  real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt)
  real(dp),intent(in) :: acell(3)
  integer,intent(inout) :: gkk_flag(elph_ds%nbranch,elph_ds%nbranch, &
  &         elph_ds%nFSkpt,elph_ds%nsppol,elph_ds%nqpt)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: qptirred(3,nqptirred)
  integer,intent(in) :: qpttoqpt(2,nsym,elph_ds%nqpt)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: spqpt(3,elph_ds%nqpt)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine complete_gkk
end interface

interface
 subroutine completeperts(elph_ds,FSfulltofull,gkk_flag,h1_mat_el,h1_mat_el_sq,hdr1,&  
  &  indsym,iqptfull,irredpert,natom,nsym,spqpt,symq,symrec,symrel,timrev,tnons)
  use defs_elphon
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iqptfull
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  type(elph_type),intent(in) :: elph_ds
  type(hdr_type),intent(in) :: hdr1
  integer,intent(in) :: FSfulltofull(2,nsym,elph_ds%nFSkpt)
  integer,intent(inout) :: gkk_flag(elph_ds%nbranch,elph_ds%nbranch, &
  &         elph_ds%nFSkpt,elph_ds%nsppol,elph_ds%nqpt)
  real(dp),intent(in) :: h1_mat_el(2,elph_ds%nFSband*elph_ds%nFSband, &
  &         elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
  real(dp),intent(out) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: irredpert(7,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nqpt)
  real(dp),intent(in) :: spqpt(3,elph_ds%nqpt)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
 end subroutine completeperts
end interface

interface
 subroutine diel9(amu,anaddb_dtset,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&  
  &  iout,lst,mpert,natom,nph2l,ntypat,phfrq,qtol,typat,ucvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nph2l
  integer,intent(in) :: ntypat
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  real(dp),intent(in) :: qtol
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(out) :: dielt_rlx(3,3)
  real(dp),intent(inout) :: displ(2,3*natom,3*natom)
  real(dp),intent(out) :: epsinf(3,3)
  real(dp),intent(out) :: fact_oscstr(2,3,3*natom)
  real(dp),intent(in) :: lst(nph2l)
  real(dp),intent(in) :: phfrq(3*natom)
  integer,intent(in) :: typat(natom)
 end subroutine diel9
end interface

interface
 subroutine dist9(acell,dist,gprim,natom,nrpt,rcan,rprim,rpt)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  real(dp),intent(in) :: acell(3)
  real(dp),intent(out) :: dist(natom,natom,nrpt)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
 end subroutine dist9
end interface

interface
 subroutine dtchi(blkval,dchide,dchidt,mpert,natom,ramansr)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: ramansr
  real(dp),intent(in) :: blkval(2,3*mpert*3*mpert*3*mpert)
  real(dp),intent(out) :: dchide(3,3,3)
  real(dp),intent(out) :: dchidt(natom,3,3,3)
 end subroutine dtchi
end interface

interface
 subroutine dtech9(blkval,dielt,iblok,mpert,natom,nblok,zeff)
  use defs_basis
  implicit none
  integer,intent(in) :: iblok
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
  real(dp),intent(out) :: dielt(3,3)
  real(dp),intent(out) :: zeff(3,3,natom)
 end subroutine dtech9
end interface

interface
 subroutine dymfz9(dynmat,natom,nqpt,gprim,option,&  
  &  rcan,spqpt,trans)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: option
  real(dp),intent(inout) :: dynmat(2,3,natom,3,natom,nqpt)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: spqpt(3,nqpt)
  real(dp),intent(in) :: trans(3,natom)
 end subroutine dymfz9
end interface

interface
 subroutine elast9(anaddb_dtset,blkval,elast,iblok,iblok_stress,instrain,iout,mpert,&  
  &  natom,nblok,ucvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iblok
  integer,intent(in) :: iblok_stress
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
  real(dp),intent(out) :: elast(6,6)
  real(dp),intent(in) :: instrain(3*natom,6)
 end subroutine elast9
end interface

interface
 subroutine electrooptic(dchide,dieflag,&  
  &  epsinf,fact_oscstr,natom,phfrq,prtmbm,rsus,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: dieflag
  integer,intent(in) :: natom
  integer,intent(in) :: prtmbm
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: dchide(3,3,3)
  real(dp),intent(in) :: epsinf(3,3)
  real(dp),intent(in) :: fact_oscstr(2,3,3*natom)
  real(dp),intent(in) :: phfrq(3*natom)
  real(dp),intent(in) :: rsus(3*natom,3,3)
 end subroutine electrooptic
end interface

interface
 subroutine eli_app_m_1d (delta_1d,lambda_1d,nmatsu,tc,z_1d)
  use defs_basis
  implicit none
  integer,intent(in) :: nmatsu
  real(dp),intent(in) :: tc
  real(dp),intent(inout) :: delta_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: z_1d(-nmatsu:nmatsu)
 end subroutine eli_app_m_1d
end interface

interface
 subroutine eli_diag_m_1d (delta_1d,lambda_1d,maxeigval,mustar,nmatsu,tc,z_1d)
  use defs_basis
  implicit none
  integer,intent(in) :: nmatsu
  real(dp),intent(out) :: maxeigval
  real(dp),intent(in) :: mustar
  real(dp),intent(in) :: tc
  real(dp),intent(inout) :: delta_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: z_1d(-nmatsu:nmatsu)
 end subroutine eli_diag_m_1d
end interface

interface
 subroutine eli_lambda_1d (a2f_1d,elph_ds,lambda_1d,nmatsu,tc)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: nmatsu
  type(elph_type),intent(in) :: elph_ds
  real(dp),intent(in) :: tc
  real(dp),intent(in) :: a2f_1d(elph_ds%na2f)
  real(dp),intent(out) :: lambda_1d(-nmatsu:nmatsu)
 end subroutine eli_lambda_1d
end interface

interface
 subroutine eli_m_iter_1d (delta_1d,lambda_1d,maxeigval,nmatsu,tc,z_1d)
  use defs_basis
  implicit none
  integer,intent(in) :: nmatsu
  real(dp),intent(out) :: maxeigval
  real(dp),intent(in) :: tc
  real(dp),intent(inout) :: delta_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: z_1d(-nmatsu:nmatsu)
 end subroutine eli_m_iter_1d
end interface

interface
 subroutine eli_z_1d (lambda_1d,nmatsu,tc,z_1d)
  use defs_basis
  implicit none
  integer,intent(in) :: nmatsu
  real(dp),intent(in) :: tc
  real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu)
  real(dp),intent(out) :: z_1d(-nmatsu:nmatsu)
 end subroutine eli_z_1d
end interface

interface
 subroutine eliashberg_1d(a2f_1d,elph_ds,mustar,n0,natom,nsym,phon_ds)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  type(elph_type),intent(in) :: elph_ds
  real(dp),intent(in) :: mustar
  type(phon_type),intent(in) :: phon_ds
  real(dp),intent(in) :: a2f_1d(elph_ds%na2f)
  real(dp),intent(in) :: n0(elph_ds%nsppol)
 end subroutine eliashberg_1d
end interface

interface
 subroutine elphon(a2fsmear,acell_in,amu,atmfrc,blkflg,blkqpt,blknrm,blktyp,blkval,&  
  &  brav,ddkfilename,dielt,dipdip,dyewq0,elphsmear,elph_fermie,elph_base_name,enunit,ep_b_min,ep_b_max,&  
  &  gkk2exist,gkk2write,gkk_rptexist,gkk_rptwrite,gkqexist,gkqwrite,&  
  &  phfrqexist,phfrqwrite,prtfsurf,prtnest,gmet,gprim,&  
  &  ifcflag,ifltransport,indsym,kptrlatt,mpert,mpi_enreg,msym,&  
  &  mustar,natom,nblok,ngqpt,nqpath,nqshft,nrpt,nsym,ntypat,qpath,&  
  &  q1shft,rcan,rmet,rprim_in,rpt,symrec,symrel,telphint,tkeepbands,&  
  &  doscalprod,tnons,trans,typat,ucvol,unitgkk,wghatm,xred,zeff)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: brav
  integer,intent(in) :: dipdip
  integer,intent(in) :: doscalprod
  integer,intent(in) :: enunit
  integer,intent(in) :: ep_b_max
  integer,intent(in) :: ep_b_min
  integer,intent(in) :: gkk2exist
  integer,intent(in) :: gkk2write
  integer,intent(in) :: gkk_rptexist
  integer,intent(in) :: gkk_rptwrite
  integer,intent(in) :: gkqexist
  integer,intent(in) :: gkqwrite
  integer,intent(in) :: ifcflag
  integer,intent(in) :: ifltransport
  integer,intent(in) :: mpert
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  integer,intent(in) :: nqpath
  integer,intent(inout) :: nqshft
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: phfrqexist
  integer,intent(in) :: phfrqwrite
  integer,intent(in) :: prtfsurf
  integer,intent(in) :: prtnest
  integer :: telphint
  integer,intent(in) :: tkeepbands
  integer,intent(in) :: unitgkk
  real(dp),intent(in) :: a2fsmear
  character(len=fnlen),intent(in) :: ddkfilename
  character(len=fnlen) :: elph_base_name
  real(dp),intent(in) :: elph_fermie
  real(dp),intent(in) :: elphsmear
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: mustar
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: kptrlatt(3,3)
  integer,intent(in) :: ngqpt(3)
  real(dp),intent(in) :: acell_in(3)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(in) :: atmfrc(2,3,natom,3,natom,nrpt)
  integer,intent(in) :: blkflg(3*mpert*3*mpert,nblok)
  real(dp),intent(in) :: blknrm(1,nblok)
  real(dp),intent(in) :: blkqpt(9,nblok)
  integer,intent(in) :: blktyp(nblok)
  real(dp),intent(in) :: blkval(2,3*mpert*3*mpert,nblok)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(in) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: q1shft(3,4)
  real(dp),intent(in) :: qpath(3,nqpath)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim_in(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(in) :: trans(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine elphon
end interface

interface
 subroutine ewald9(acell,dielt,dyew,gmet,gprim,natom,&  
  &  qphon,rmet,rprim,sumg0,ucvol,xred,zeff)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: sumg0
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(out) :: dyew(2,3,natom,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine ewald9
end interface

interface
 subroutine ftgam (wghatm,gam_qpt,gam_rpt,gprim,natom,nqpt,nrpt,qtor,rpt,spqpt)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: nrpt
  integer,intent(in) :: qtor
  real(dp),intent(inout) :: gam_qpt(2,3*natom*3*natom,nqpt)
  real(dp),intent(inout) :: gam_rpt(2,3*natom*3*natom,nrpt)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: spqpt(3,nqpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine ftgam
end interface

interface
 subroutine ftgkk (wghatm,gkk_qpt,gkk_rpt,gkqwrite,gkrwrite,gprim,iFSkpt0,&  
  &  natom,nFSkpt,ngkkband,nkpt_used,nqpt,nrpt,nsppol,&  
  &  qtor,rpt,spqpt,unit_gkk_rpt,unitgkq)
  use defs_basis
  implicit none
  integer,intent(in) :: gkqwrite
  integer,intent(in) :: gkrwrite
  integer,intent(in) :: iFSkpt0
  integer,intent(in) :: nFSkpt
  integer,intent(in) :: natom
  integer,intent(in) :: ngkkband
  integer,intent(in) :: nkpt_used
  integer,intent(in) :: nqpt
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: qtor
  integer,intent(in) :: unit_gkk_rpt
  integer,intent(in) :: unitgkq
  real(dp),intent(inout) :: gkk_qpt(2,ngkkband*ngkkband,3*natom*3*natom,nkpt_used,nsppol,nqpt)
  real(dp),intent(inout) :: gkk_rpt(2,ngkkband*ngkkband,3*natom*3*natom,nkpt_used,nsppol,nrpt)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: spqpt(3,nqpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine ftgkk
end interface

interface
 subroutine ftiaf9(atmfrc,dynmat,gprim,natom,nqpt,&  
  &  nrpt,qtor,rpt,spqpt,wghatm)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: nrpt
  integer,intent(in) :: qtor
  real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(inout) :: dynmat(2,3,natom,3,natom,nqpt)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: spqpt(3,nqpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine ftiaf9
end interface

interface
 subroutine fxgkkphase(elph_ds,gkk_flag,h1_mat_el,iqptfull)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptfull
  type(elph_type),intent(in) :: elph_ds
  integer,intent(in) :: gkk_flag(elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nqpt)
  real(dp),intent(inout) :: h1_mat_el(2,elph_ds%nFSband,elph_ds%nFSband, &
  &         elph_ds%nbranch,elph_ds%nFSkpt)
 end subroutine fxgkkphase
end interface

interface
 subroutine gamma9(gamma,qphon,qphnrm,qtol)
  use defs_basis
  implicit none
  integer,intent(out) :: gamma
  real(dp),intent(in) :: qphnrm
  real(dp),intent(in) :: qtol
  real(dp),intent(in) :: qphon(3)
 end subroutine gamma9
end interface

interface
 subroutine get_all_gkk2(acell,amu,atmfrc,dielt,dipdip,dyewq0,&  
  &  elph_ds,FSkptirred,FSkpt,&  
  &  ftwghtgkk,gmet,gprim,indsym,mpert,msym,&  
  &  natom,nrpt,nsym,ntypat,&  
  &  onegkksize,phon_ds,rcan,rmet,rprim,rprimd,&  
  &  rpt,spqpt,symrel,trans,typat,ucvol,&  
  &  wghatm,xred,zeff)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: dipdip
  integer,intent(in) :: mpert
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: onegkksize
  type(elph_type),intent(inout) :: elph_ds
  type(phon_type),intent(inout) :: phon_ds
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt)
  real(dp),intent(in) :: FSkptirred(3,elph_ds%nFSkptirred)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(in) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(in) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: ftwghtgkk(natom,nrpt)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: spqpt(3,elph_ds%nqpt)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: trans(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine get_all_gkk2
end interface

interface
 subroutine get_all_gkq (acell,amu,elph_ds,FSfullpqtofull,FSfulltofull,FSkpt,FSintweight,&  
  &  gkk_flag,gprimd,indsym,mpert,natom,nband,nqptirred,nsym,ntypat,n1wf,onegkksize,phon_ds,&  
  &  qptirred,qpttoqpt,rprimd,spqpt,symrec,symrel,timrev,tnons,typat,ucvol,unitgkk,xred)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: n1wf
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(inout) :: nqptirred
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: onegkksize
  integer,intent(in) :: timrev
  integer,intent(in) :: unitgkk
  type(elph_type),intent(inout) :: elph_ds
  type(phon_type),intent(inout) :: phon_ds
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
  integer,intent(in) :: FSfulltofull(2,nsym,elph_ds%nFSkpt)
  real(dp),intent(in) :: FSintweight(elph_ds%nFSband,elph_ds%nFSkpt)
  real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(ntypat)
  integer,intent(inout) :: gkk_flag(elph_ds%nbranch,elph_ds%nbranch, &
  &         elph_ds%nFSkpt,elph_ds%nsppol,elph_ds%nqpt)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(inout) :: qptirred(3,n1wf)
  integer,intent(in) :: qpttoqpt(2,nsym,elph_ds%nqpt)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: spqpt(3,elph_ds%nqpt)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine get_all_gkq
end interface

interface
 subroutine get_all_gkr (elph_ds,gprim,natom,nrpt,onegkksize,rpt,spqpt,wghatm)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  integer,intent(in) :: onegkksize
  type(elph_type),intent(inout) :: elph_ds
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: spqpt(3,elph_ds%nqpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine get_all_gkr
end interface

interface
 subroutine get_fs_kpts(eigenGS,elph_ds,FSkptflag,gaussig,hdr)
  use defs_elphon
  use defs_basis
  use defs_datatypes
  implicit none
  type(elph_type),intent(inout) :: elph_ds
  real(dp),intent(out) :: gaussig
  type(hdr_type),intent(in) :: hdr
  integer,intent(out) :: FSkptflag(hdr%nkpt)
  real(dp),intent(in) :: eigenGS(hdr%nband(1),hdr%nkpt,hdr%nsppol)
 end subroutine get_fs_kpts
end interface

interface
 subroutine get_gkk_qpt_tr(elph_ds,mpi_enreg,nband,hdr,FSfulltoirred, FSirredtoGS,FSfullpqtofull,elph_tr_ds)
  use defs_elphon
  use defs_datatypes
  implicit none
  integer,intent(in) :: nband
  type(elph_type),intent(in) :: elph_ds
  type(elph_tr_type) :: elph_tr_ds
  type(hdr_type),intent(in) :: hdr
  type(mpi_type), intent(inout) :: mpi_enreg
  integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
  integer,intent(in) :: FSfulltoirred(3,elph_ds%nFSkpt)
  integer,intent(in) :: FSirredtoGS(elph_ds%nFSkptirred)
 end subroutine get_gkk_qpt_tr
end interface

interface
 subroutine gtblk9(blkflg,blknrm,blkqpt,blktyp,iblok,mpert,msize,&  
  &  natom,nblok,qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)
  use defs_basis
  implicit none
  integer,intent(out) :: iblok
  integer,intent(in) :: mpert
  integer,intent(in) :: msize
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  integer,intent(in) :: rftyp
  real(dp),intent(in) :: qtol
  integer,intent(in) :: rfelfd(4)
  integer,intent(in) :: rfphon(4)
  integer,intent(in) :: rfstrs(4)
  integer,intent(in) :: blkflg(msize,nblok)
  real(dp),intent(in) :: blknrm(3,nblok)
  real(dp),intent(in) :: blkqpt(3,3,nblok)
  integer,intent(in) :: blktyp(nblok)
  real(dp),intent(inout) :: qphnrm(3)
  real(dp),intent(inout) :: qphon(3,3)
 end subroutine gtblk9
end interface

interface
 subroutine gtdyn9(acell,atmfrc,dielt,dipdip,&  
  &  dyewq0,d2cart,gmet,gprim,mpert,natom,&  
  &  nrpt,qphnrm,qpt,rcan,rmet,rprim,rpt,&  
  &  trans,ucvol,wghatm,xred,zeff)
  use defs_basis
  implicit none
  integer,intent(in) :: dipdip
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  real(dp),intent(in) :: qphnrm
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(inout) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: qpt(3)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: trans(3,natom)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine gtdyn9
end interface

interface
 subroutine hybrid9(acell,asr,atmfrc,dielt,dipdip,dyew,dyewq0,&  
  &  gmet,gprim,iout,natom,nrpt,rcan,rmet,&  
  &  rprim,rpt,ucvol,wghatm,xred,zeff)
  use defs_basis
  implicit none
  integer,intent(in) :: asr
  integer,intent(in) :: dipdip
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(inout) :: dielt(3,3)
  real(dp),intent(inout) :: dyew(2,3,natom,3,natom)
  real(dp),intent(inout) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(inout) :: zeff(3,3,natom)
 end subroutine hybrid9
end interface

interface
 subroutine ifclo9(ifccar,ifcloc,vect1,vect2,vect3)
  use defs_basis
  implicit none
  real(dp),intent(in) :: ifccar(3,3)
  real(dp),intent(out) :: ifcloc(3,3)
  real(dp),intent(in) :: vect1(3)
  real(dp),intent(in) :: vect2(3)
  real(dp),intent(in) :: vect3(3)
 end subroutine ifclo9
end interface

interface
 subroutine init8(dscrpt,filnam,mddb,nddb)
  use defs_basis
  implicit none
  integer,intent(in) :: mddb
  integer,intent(out) :: nddb
  character(len=fnlen),intent(out) :: dscrpt
  character(len=fnlen),intent(out) :: filnam(mddb+1)
 end subroutine init8
end interface

interface
 subroutine init9(filnam)
  use defs_basis
  implicit none
  character(len=fnlen),intent(out) :: filnam(7)
 end subroutine init9
end interface

interface
 subroutine inpphon(displ,pheigval,pheigvec,phfrq,phon_ds,qpt)
  use defs_elphon
  use defs_basis
  implicit none
  type(phon_type),intent(inout) :: phon_ds
  real(dp),intent(out) :: displ(2,3*phon_ds%natom,3*phon_ds%natom)
  real(dp),intent(out) :: pheigval(3*phon_ds%natom)
  real(dp),intent(out) :: pheigvec(2*3*phon_ds%natom*3*phon_ds%natom)
  real(dp),intent(out) :: phfrq(3*phon_ds%natom)
  real(dp),intent(inout) :: qpt(3)
 end subroutine inpphon
end interface

interface
 subroutine instr9(blkval,iblok,instrain,iout,mpert,natom,nblok)
  use defs_basis
  implicit none
  integer,intent(in) :: iblok
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
  real(dp),intent(out) :: instrain(3*natom,6)
 end subroutine instr9
end interface

interface
 subroutine integrate_gamma(elph_ds,FSfullpqtofull,gprim,n0,natom,nrpt,rpt,spqpt,wghatm)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: n0(elph_ds%nsppol)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: spqpt(3,elph_ds%nqpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine integrate_gamma
end interface

interface
 subroutine integrate_gamma_tr(elph_ds,FSfullpqtofull,gprim,natom,nrpt,rpt,spqpt,wghatm,elph_tr_ds)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  type(elph_type),intent(inout) :: elph_ds
  type(elph_tr_type) :: elph_tr_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: spqpt(3,elph_ds%nqpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine integrate_gamma_tr
end interface

interface
 subroutine interpolate_gkk (acell,amu,atmfrc,dielt,dipdip,&  
  &  dyewq0,elph_ds,FSkptirred,FSkpt,ftwghtgkk,&  
  &  gmet,gprim,indsym,mpert,msym,natom,&  
  &  nrpt,nsym,ntypat,phon_ds,rcan,rmet,rprim,rprimd,rpt,spqpt,&  
  &  symrel,trans,typat,ucvol,wghatm,xred,zeff)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: dipdip
  integer,intent(in) :: mpert
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  type(elph_type),intent(inout) :: elph_ds
  type(phon_type),intent(inout) :: phon_ds
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt)
  real(dp),intent(in) :: FSkptirred(3,elph_ds%nFSkptirred)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(in) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(in) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: ftwghtgkk(natom,nrpt)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: spqpt(3,elph_ds%nqpt)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: trans(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine interpolate_gkk
end interface

interface
 subroutine interpolate_phfrq (acell,amu,atmfrc,dielt,dipdip,&  
  &  dyewq0,FSkptirred,FSkpt,gmet,&  
  &  gprim,indsym,mpert,msym,natom,&  
  &  nbranch,nFSband,nFSkpt,nFSkptirred,&  
  &  nrpt,nsym,ntypat,phfrq,&  
  &  rcan,rmet,rprim,rprimd,rpt,&  
  &  symrel,trans,typat,ucvol,&  
  &  wghatm,xred,zeff)
  use defs_basis
  implicit none
  integer,intent(in) :: dipdip
  integer,intent(in) :: mpert
  integer,intent(in) :: msym
  integer,intent(in) :: nFSband
  integer,intent(in) :: nFSkpt
  integer,intent(in) :: nFSkptirred
  integer,intent(in) :: natom
  integer,intent(in) :: nbranch
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: FSkpt(3,nFSkpt)
  real(dp),intent(in) :: FSkptirred(3,nFSkptirred)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(in) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(in) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(out) :: phfrq(3*natom,nFSkpt,nFSkptirred)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: trans(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine interpolate_phfrq
end interface

interface
 subroutine invars9 (anaddb_dtset,lenstr,natom,nunit,qtol,string)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: lenstr
  integer,intent(in) :: natom
  integer,intent(in) :: nunit
  type(anaddb_dataset_type),intent(out) :: anaddb_dtset
  real(dp),intent(in) :: qtol
  character(len=*),intent(in) :: string
 end subroutine invars9
end interface

interface
 subroutine mk_irredpert(hdr1,iatom,idir,indsym,iqptfull,irredpert,&  
  &  natom,nbranch,nqpt,nsym,qpt,qtimrev,symq,symrel)
  use defs_datatypes
  implicit none
  integer,intent(in) :: iatom
  integer,intent(in) :: idir
  integer,intent(in) :: iqptfull
  integer,intent(in) :: natom
  integer,intent(in) :: nbranch
  integer,intent(in) :: nqpt
  integer,intent(in) :: nsym
  integer,intent(in) :: qtimrev
  type(hdr_type),intent(in) :: hdr1
  integer,intent(in) :: qpt(3)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(out) :: irredpert(7,nbranch,nbranch,nqpt)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine mk_irredpert
end interface

interface
 subroutine mka2f(a2f_1d,dos_phon,elph_ds,FSintweight,FSirredtofull,FSirredwtk,FSkpt,gprim,gprimd,mustar,&  
  &  n0,natom,nrpt,phon_ds,rpt,wghatm)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  type(elph_type),intent(inout) :: elph_ds
  real(dp),intent(in) :: mustar
  type(phon_type),intent(inout) :: phon_ds
  real(dp),intent(in) :: FSintweight(elph_ds%nFSband,elph_ds%nFSkpt)
  integer,intent(in) :: FSirredtofull(elph_ds%nFSkptirred)
  real(dp),intent(in) :: FSirredwtk(elph_ds%nFSkptirred)
  real(dp),intent(inout) :: FSkpt(3,elph_ds%nFSkpt)
  real(dp),intent(out) :: a2f_1d(elph_ds%na2f)
  real(dp),intent(out) :: dos_phon(elph_ds%na2f)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: n0(elph_ds%nsppol)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine mka2f
end interface

interface
 subroutine mka2fQgrid(elph_ds,fname)
  use defs_elphon
  use defs_basis
  implicit none
  type(elph_type),intent(in) :: elph_ds
  character(len=fnlen),intent(in) :: fname
 end subroutine mka2fQgrid
end interface

interface
 subroutine mka2f_tr(elph_ds,FSkpt,gprim,n0,ucvol,natom,nrpt,phon_ds,rpt,wghatm,elph_tr_ds)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  type(elph_type),intent(inout) :: elph_ds
  type(elph_tr_type) :: elph_tr_ds
  type(phon_type),intent(inout) :: phon_ds
  real(dp),intent(in) :: ucvol
  real(dp),intent(inout) :: FSkpt(3,elph_ds%nFSkpt)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: n0(elph_ds%nsppol)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine mka2f_tr
end interface

interface
 subroutine mkFSkgrid (FSirredwtk,FSirredtofull,tmpFSfulltoirred,&  
  &  tmpFSfulltofull,tmpFSkpt,FSkptirred,nFSkptirred,nFSkpt,&  
  &  nsym,symrec,timrev)
  use defs_basis
  implicit none
  integer,intent(out) :: nFSkpt
  integer,intent(in) :: nFSkptirred
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  integer,intent(out) :: FSirredtofull(nFSkptirred)
  real(dp),intent(out) :: FSirredwtk(nFSkptirred)
  real(dp),intent(in) :: FSkptirred(3,nFSkptirred)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(out) :: tmpFSfulltofull(2,nsym,2*nFSkptirred*nsym)
  integer,intent(out) :: tmpFSfulltoirred(3,2*nFSkptirred*nsym)
  real(dp),intent(out) :: tmpFSkpt(3,2*nFSkptirred*nsym)
 end subroutine mkFSkgrid
end interface

interface
 subroutine mkfsqgrid(FSkpt,FStoqpt,nFSkpt,nFSqpt,tmpFSqpt)
  use defs_basis
  implicit none
  integer,intent(in) :: nFSkpt
  integer,intent(out) :: nFSqpt
  real(dp),intent(in) :: FSkpt(3,nFSkpt)
  integer,intent(out) :: FStoqpt(nFSkpt,nFSkpt)
  real(dp),intent(out) :: tmpFSqpt(3,nFSkpt*nFSkpt)
 end subroutine mkfsqgrid
end interface

interface
 subroutine mkifc9(acell,amu,anaddb_dtset,atmfrc,&  
  &  blkflg,blknrm,blkqpt,blktyp,blkval,dielt,displ,dyewq0,&  
  &  d2cart,eigval,eigvec,gmet,gprim,&  
  &  indsym,iout,mpert,msym,natom,nblok,nrpt,&  
  &  nsym,ntypat,phfrq,rcan,rmet,rprim,&  
  &  rpt,symrec,symrel,tcpui,tnons,trans,twalli,typat,&  
  &  ucvol,wghatm,xred,zeff)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  integer,intent(inout) :: nrpt
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  real(dp),intent(in) :: tcpui
  real(dp),intent(in) :: twalli
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: blkflg(*)
  integer,intent(in) :: blktyp(*)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(out) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(in) :: blknrm(3,*)
  real(dp),intent(in) :: blkqpt(9,*)
  real(dp),intent(in) :: blkval(2,*)
  real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(inout) :: dielt(3,3)
  real(dp),intent(out) :: displ(2*3*natom*3*natom)
  real(dp),intent(out) :: dyewq0(3,3,natom)
  real(dp),intent(out) :: eigval(*)
  real(dp),intent(out) :: eigvec(*)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym*natom)
  real(dp),intent(out) :: phfrq(*)
  real(dp),intent(out) :: rcan(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(out) :: rpt(3,nrpt)
  integer,intent(in) :: symrec(3,3,msym)
  integer,intent(in) :: symrel(3,3,msym)
  real(dp),intent(in) :: tnons(3,msym)
  real(dp),intent(out) :: trans(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(out) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(inout) :: zeff(3,3,natom)
 end subroutine mkifc9
end interface

interface
 subroutine mkkptrank (kpt,nkpt,rank,invrank)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpt
  integer,intent(out) :: invrank(16000000)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(out) :: rank(nkpt)
 end subroutine mkkptrank
end interface

interface
 subroutine mknesting(nkpt,kpt,kptrlatt,nband,weight,nsegment_in,npoint_in,&  
  &  qpath_vertices_in,elph_base_name,gprim,gprimd,rprim,brav,prtnest)
  use defs_basis
  implicit none
  integer,intent(in) :: brav
  integer,intent(in) :: nband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsegment_in
  integer,intent(in) :: prtnest
  character(len=fnlen),intent(in) :: elph_base_name
  integer,intent(in) :: kptrlatt(3,3)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: npoint_in(nsegment_in)
  real(dp),intent(in) :: qpath_vertices_in(3,nsegment_in+1)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: weight(nband,nkpt)
 end subroutine mknesting
end interface

interface
 subroutine mkph_linwid(elph_ds,FSintweight,FSfulltoirred,FSirredtofull,&  
  &  gmet,gprim,gprimd,n0,natom,npoint_in,nrpt,nsegment_in,nsym,phon_ds,&  
  &  qpath_vertices_in,qpttoqpt,rpt,spqpt,wghatm)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsegment_in
  integer,intent(in) :: nsym
  type(elph_type),intent(inout) :: elph_ds
  type(phon_type),intent(inout) :: phon_ds
  integer,intent(in) :: FSfulltoirred(3,elph_ds%nFSkpt)
  real(dp),intent(in) :: FSintweight(elph_ds%nFSband,elph_ds%nFSkpt,elph_ds%nsppol)
  integer,intent(in) :: FSirredtofull(elph_ds%nFSkptirred)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: n0(elph_ds%nsppol)
  integer,intent(in) :: npoint_in(nsegment_in)
  real(dp),intent(in) :: qpath_vertices_in(3,nsegment_in+1)
  integer,intent(in) :: qpttoqpt(2,nsym,elph_ds%nqpt)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: spqpt(3,elph_ds%nqpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine mkph_linwid
end interface

interface
 subroutine mkphdos(acell,amu,anaddb_dtset,atmfrc,dielt,dyewq0,filname,gmet,gprim,indsym,&  
  &  mpert,msym,natom,nrpt,nsym,ntypat,rcan,rmet,rprim,rpt,symrec,symrel,tcpui,&  
  &  trans,twalli,typat,ucvol,wghatm,xred,zeff)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  character(len=fnlen),intent(in) :: filname
  real(dp),intent(in) :: tcpui
  real(dp),intent(in) :: twalli
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(inout) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: trans(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine mkphdos
end interface

interface
 subroutine mkqptequiv (FSfullpqtofull,FSkpt,FSkptirred,nFSkpt,nFSkptirred,nqpt,nsym,&  
  &  qpttoqpt,spqpt,symrec)
  use defs_basis
  implicit none
  integer,intent(in) :: nFSkpt
  integer,intent(in) :: nFSkptirred
  integer,intent(in) :: nqpt
  integer,intent(in) :: nsym
  integer,intent(out) :: FSfullpqtofull(nFSkpt,nqpt)
  real(dp),intent(in) :: FSkpt(3,nFSkpt)
  real(dp),intent(in) :: FSkptirred(3,nFSkptirred)
  integer,intent(out) :: qpttoqpt(2,nsym,nqpt)
  real(dp),intent(in) :: spqpt(3,nqpt)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine mkqptequiv
end interface

interface
 subroutine nanal9(dyew,dynmat,iqpt,natom,nqpt,plus)
  use defs_basis
  implicit none
  integer,intent(in) :: iqpt
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: plus
  real(dp),intent(in) :: dyew(2,3,natom,3,natom)
  real(dp),intent(out) :: dynmat(2,3,natom,3,natom,nqpt)
 end subroutine nanal9
end interface

interface
 subroutine nmsq_gam (accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&  
  &  FSintweight,FSkpt,gkk_qpt_tmp,h1_mat_el_sq,iqptfull,phfrq_tmp,spqpt,wf)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptfull
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
  real(dp),intent(in) :: FSintweight(elph_ds%nFSband,elph_ds%nFSkpt,elph_ds%nsppol)
  real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt)
  real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(inout) :: gkk_qpt_tmp(2,elph_ds%ngkkband*elph_ds%ngkkband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
  real(dp),intent(in) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
  real(dp),intent(in) :: phfrq_tmp(elph_ds%nbranch)
  real(dp),intent(in) :: spqpt(3,elph_ds%nqpt)
  real(dp),intent(in) :: wf(elph_ds%nbranch)
 end subroutine nmsq_gam
end interface

interface
 subroutine nmsq_gam_sumFS(accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&  
  &  FSintweight,FSkpt,gkk_qpt_tmp,h1_mat_el_sq,iqptfull,phfrq_tmp,spqpt,wf)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptfull
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
  real(dp),intent(in) :: FSintweight(elph_ds%nFSband,elph_ds%nFSkpt,elph_ds%nsppol)
  real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt)
  real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(inout) :: gkk_qpt_tmp(2,elph_ds%ngkkband*elph_ds%ngkkband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
  real(dp),intent(in) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband*elph_ds%nbranch, &
  &         elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
  real(dp),intent(in) :: phfrq_tmp(elph_ds%nbranch)
  real(dp),intent(in) :: spqpt(3,elph_ds%nqpt)
  real(dp),intent(in) :: wf(elph_ds%nbranch)
 end subroutine nmsq_gam_sumFS
end interface

interface
 subroutine nmsq_pure_gkk(accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,FSintweight,FSkpt,gkk_qpt_tmp,&  
  &  h1_mat_el_sq,iqptfull,phfrq_tmp,spqpt,wf)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptfull
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
  real(dp),intent(in) :: FSintweight(elph_ds%nFSband,elph_ds%nFSkpt,elph_ds%nsppol)
  real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt)
  real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(inout) :: gkk_qpt_tmp(2,elph_ds%ngkkband*elph_ds%ngkkband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
  real(dp),intent(in) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
  real(dp),intent(in) :: phfrq_tmp(elph_ds%nbranch)
  real(dp),intent(in) :: spqpt(3,elph_ds%nqpt)
  real(dp),intent(in) :: wf(elph_ds%nbranch)
 end subroutine nmsq_pure_gkk
end interface

interface
 subroutine nmsq_pure_gkk_sumfs(accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,FSintweight,FSkpt,gkk_qpt_tmp,&  
  &  h1_mat_el_sq,iqptfull,phfrq_tmp,spqpt,wf)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptfull
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
  real(dp),intent(in) :: FSintweight(elph_ds%nFSband,elph_ds%nFSkpt,elph_ds%nsppol)
  real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt)
  real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(inout) :: gkk_qpt_tmp(2,elph_ds%ngkkband*elph_ds%ngkkband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
  real(dp),intent(in) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
  real(dp),intent(in) :: phfrq_tmp(elph_ds%nbranch)
  real(dp),intent(in) :: spqpt(3,elph_ds%nqpt)
  real(dp),intent(in) :: wf(elph_ds%nbranch)
 end subroutine nmsq_pure_gkk_sumfs
end interface

interface
 subroutine normsq_gkq(displ_red,eigvec,elph_ds,FSfullpqtofull,&  
  &  FSintweight,FSkpt,h1_mat_el_sq,iqptfull,phfrq_tmp,spqpt,wf,qdata)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptfull
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
  real(dp),intent(in) :: FSintweight(elph_ds%nFSband,elph_ds%nFSkpt)
  real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt)
  real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(in) :: h1_mat_el_sq(2,elph_ds%nFSband,elph_ds%nFSband, &
  &         elph_ds%nbranch,elph_ds%nbranch,elph_ds%nFSkpt)
  real(dp),intent(in) :: phfrq_tmp(elph_ds%nbranch)
  real(dp),intent(out) :: qdata(elph_ds%nbranch,elph_ds%nsppol,3)
  real(dp),intent(in) :: spqpt(3,elph_ds%nqpt)
  real(dp),intent(in) :: wf(elph_ds%nbranch)
 end subroutine normsq_gkq
end interface

interface
 subroutine order_fs_kpts(FSkptflag,FSkptirrank,FSkptirred,FSirredtoGS,hdr,nFSkptirred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nFSkptirred
  type(hdr_type),intent(in) :: hdr
  integer,intent(out) :: FSirredtoGS(nFSkptirred)
  integer,intent(in) :: FSkptflag(hdr%nkpt)
  integer,intent(out) :: FSkptirrank(nFSkptirred)
  real(dp),intent(inout) :: FSkptirred(3,nFSkptirred)
 end subroutine order_fs_kpts
end interface

interface
 subroutine outelph(elph_ds,enunit,fname)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: enunit
  type(elph_type),intent(in) :: elph_ds
  character(len=fnlen),intent(in) :: fname
 end subroutine outelph
end interface

interface
 subroutine outlwf9 (acell,iodyn,msym,natom,nph1l,nsym,ntypat,rprim,symrel,typat,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: iodyn
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: nph1l
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
  integer,intent(in) :: symrel(3,3,msym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine outlwf9
end interface

interface
 subroutine outvars9 (anaddb_dtset,nunit)
  use defs_datatypes
  implicit none
  integer,intent(in) :: nunit
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
 end subroutine outvars9
end interface

interface
 subroutine piezo9(anaddb_dtset,blkval,dielt_rlx,elast,iblok,instrain,iout,mpert,&  
  &  natom,nblok,piezo,ucvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iblok
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
  real(dp),intent(in) :: dielt_rlx(3,3)
  real(dp),intent(in) :: elast(6,6)
  real(dp),intent(in) :: instrain(3*natom,6)
  real(dp),intent(out) :: piezo(6,3)
 end subroutine piezo9
end interface

interface
 subroutine printxsf(n1,n2,n3,datagrid,basis,origin,nunit,option)
  use defs_basis
  implicit none
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: nunit
  integer,intent(in) :: option
  real(dp),intent(in) :: basis(3,3)
  real(dp),intent(in) :: datagrid(n1*n2*n3)
  real(dp),intent(in) :: origin(3)
 end subroutine printxsf
end interface

interface
 subroutine ramansus(d2cart,dchide,dchidt,displ,mpert,&  
  &  natom,phfrq,qphon,qphnrm,rsus,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  real(dp),intent(in) :: qphnrm
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(in) :: dchide(3,3,3)
  real(dp),intent(in) :: dchidt(natom,3,3,3)
  real(dp),intent(in) :: displ(2,3*natom,3*natom)
  real(dp),intent(in) :: phfrq(3*natom)
  real(dp),intent(inout) :: qphon(3)
  real(dp),intent(out) :: rsus(3*natom,3,3)
 end subroutine ramansus
end interface

interface
 subroutine rchkGSheader (hdr,natom,nband,unitgkk)
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(out) :: nband
  integer,intent(in) :: unitgkk
  type(hdr_type),intent(out) :: hdr
 end subroutine rchkGSheader
end interface

interface
 subroutine rdddb9(acell,atifc,amu,blkflg,blknrm,blkqpt,&  
  &  blktyp,blkval,ddbun,dimekb,filnam,gmet,gprim,indsym,iout,&  
  &  mband,mpert,msize,msym,&  
  &  natifc,natom,nblok,nkpt,nsym,ntypat,&  
  &  occopt,rmet,rprim,symq,symrec,symrel,thmflag,&  
  &  tnons,typat,ucvol,usepaw,xcart,xred,zion,blkval2,kpnt)
  use defs_basis
  implicit none
  integer,intent(in) :: ddbun
  integer,intent(in) :: dimekb
  integer,intent(in) :: iout
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: msize
  integer,intent(in) :: msym
  integer,intent(in) :: natifc
  integer,intent(inout) :: natom
  integer,intent(inout) :: nblok
  integer,intent(inout) :: nkpt
  integer,intent(inout) :: nsym
  integer,intent(inout) :: ntypat
  integer,intent(inout) :: occopt
  integer,intent(in) :: thmflag
  integer,intent(in) :: usepaw
  character(len=fnlen),intent(in) :: filnam
  real(dp),intent(out) :: ucvol
  integer,intent(out) :: symq(4,2,*)
  real(dp),intent(out) :: acell(3)
  real(dp),intent(out) :: amu(ntypat)
  integer,intent(inout) :: atifc(natom)
  integer,intent(out) :: blkflg(msize,nblok)
  real(dp),intent(out) :: blknrm(3,nblok)
  real(dp),intent(out) :: blkqpt(9,nblok)
  integer,intent(out) :: blktyp(nblok)
  real(dp),intent(out) :: blkval(2,msize,nblok)
  real(dp),intent(out),optional :: blkval2(2,msize,mband,nkpt,nblok)
  real(dp),intent(out) :: gmet(3,3)
  real(dp),intent(out) :: gprim(3,3)
  integer,intent(out) :: indsym(4,nsym,natom)
  real(dp),intent(out),optional :: kpnt(3,nkpt,nblok)
  real(dp),intent(out) :: rmet(3,3)
  real(dp),intent(out) :: rprim(3,3)
  integer,intent(out) :: symrec(3,3,msym)
  integer,intent(out) :: symrel(3,3,msym)
  real(dp),intent(out) :: tnons(3,msym)
  integer,intent(out) :: typat(natom)
  real(dp),intent(out) :: xcart(3,natom)
  real(dp),intent(out) :: xred(3,natom)
  real(dp),intent(out) :: zion(ntypat)
 end subroutine rdddb9
end interface

interface
 subroutine read_el_veloc(mpi_enreg,nbandtot,nkpttot,nsppol_in,elph_tr_ds)
  use defs_elphon
  use defs_datatypes
  implicit none
  integer,intent(in) :: nbandtot
  integer,intent(in) :: nkpttot
  integer,intent(in) :: nsppol_in
  type(elph_tr_type) :: elph_tr_ds
  type(mpi_type) :: mpi_enreg
 end subroutine read_el_veloc
end interface

interface
 subroutine read_gkk(amu,elph_ds,FSfullpqtofull,FSfulltofull,FSintweight,FSkpt,&  
  &  gkk_flag,gprimd,indsym,n1wf,natom,nband,nqptirred,nsym,ntypat,&  
  &  phon_ds,qptirred,rprimd,spqpt,symrec,symrel,timrev,tnons,typat,ucvol,unitgkk)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: n1wf
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(out) :: nqptirred
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: timrev
  integer,intent(in) :: unitgkk
  type(elph_type),intent(inout) :: elph_ds
  type(phon_type),intent(inout) :: phon_ds
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
  integer,intent(in) :: FSfulltofull(2,nsym,elph_ds%nFSkpt)
  real(dp),intent(in) :: FSintweight(elph_ds%nFSband,elph_ds%nFSkpt)
  real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt)
  real(dp),intent(in) :: amu(ntypat)
  integer,intent(out) :: gkk_flag(elph_ds%nbranch,elph_ds%nbranch, &
  &         elph_ds%nFSkpt,elph_ds%nsppol,elph_ds%nqpt)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(out) :: qptirred(3,n1wf)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: spqpt(3,elph_ds%nqpt)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(in) :: typat(natom)
 end subroutine read_gkk
end interface

interface
 subroutine relaxpol(blkflg,blkval,etotal,fred,iatfix,indsym,iout,istrfix,&  
  &  mpert,msize,msym,natfix,natom,nstrfix,nsym,ntypat,pel,relaxat,relaxstr,&  
  &  rprimd,strten,symrel,targetpol,typat,ucvol,xcart,xred,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: msize
  integer,intent(in) :: msym
  integer,intent(in) :: natfix
  integer,intent(in) :: natom
  integer,intent(in) :: nstrfix
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: relaxat
  integer,intent(in) :: relaxstr
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: istrfix(6)
  integer,intent(in) :: blkflg(msize)
  real(dp),intent(inout) :: blkval(2,msize)
  real(dp),intent(in) :: fred(3,natom)
  integer,intent(in) :: iatfix(natom)
  integer,intent(in) :: indsym(4,msym,natom)
  real(dp),intent(in) :: pel(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: strten(6)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(inout) :: targetpol(3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xcart(3,natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine relaxpol
end interface

interface
 subroutine rsiaf9(acell,atifc,atmfrc,dielt,dipdip,dyewq0,&  
  &  gmet,gprim,ifcana,ifcout,iout,natom,nrpt,nsphere,rcan,&  
  &  rifcsph,rmet,rprim,rpt,tcpui,twalli,ucvol,wghatm,zeff)
  use defs_basis
  implicit none
  integer,intent(in) :: dipdip
  integer,intent(in) :: ifcana
  integer,intent(in) :: ifcout
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsphere
  real(dp),intent(in) :: rifcsph
  real(dp),intent(in) :: tcpui
  real(dp),intent(in) :: twalli
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  integer,intent(in) :: atifc(natom)
  real(dp),intent(in) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(in) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(inout) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine rsiaf9
end interface

interface
 subroutine setup_phon_ds(phon_ds,dipdip,mpert,nsym,natom,ntypat,nrpt,&  
  &  ucvol,indsym,symrel,typat,acell,amu,atmfrc,dielt,dyewq0,gprim,gmet,&  
  &  xred,zeff,rcan,rmet,rprim,rprimd,rpt,trans,wghatm)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: dipdip
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  type(phon_type) :: phon_ds
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(in) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(in) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: trans(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine setup_phon_ds
end interface

interface
 subroutine sortph(displ,filnam,&  
  &  iout,natom,phfrq,qphnrm,qphon,udispl,ufreq)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  integer,intent(in) :: udispl
  integer,intent(in) :: ufreq
  character(len=fnlen),intent(in) :: filnam
  real(dp),intent(in) :: qphnrm
  real(dp),intent(in) :: displ(2,3*natom,3*natom)
  real(dp),intent(in) :: phfrq(3*natom)
  real(dp),intent(in) :: qphon(3)
 end subroutine sortph
end interface

interface
 subroutine sym_gkk(acell,FSfulltofull,FSkpt,gkk_flag,gkk_qpt,&  
  &  gprim,indsym,mpert,natom,nbranch,nFSband,nFSkpt,nqpt,nqptirred,nsym,&  
  &  qptirred,rprim,spqpt,symrec,symrel,tnons,ucvol,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: nFSband
  integer,intent(in) :: nFSkpt
  integer,intent(in) :: natom
  integer,intent(in) :: nbranch
  integer,intent(in) :: nqpt
  integer,intent(in) :: nqptirred
  integer,intent(in) :: nsym
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: FSfulltofull(2,nsym,nFSkpt)
  real(dp),intent(in) :: FSkpt(3,nFSkpt)
  real(dp),intent(in) :: acell(3)
  integer,intent(inout) :: gkk_flag(nbranch,nFSkpt,nqpt)
  real(dp),intent(inout) :: gkk_qpt(2,nbranch,nFSband,nFSband,nFSkpt,nqpt)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: qptirred(3,nqptirred)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: spqpt(3,nqpt)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine sym_gkk
end interface

interface
 subroutine symdm9(acell,blkflg,blknrm,blkqpt,blktyp,blkval,&  
  &  dynmat,gprim,indsym,mpert,natom,nblok,nqpt,nsym,rfmeth,&  
  &  rprim,spqpt,symrec,symrel,xred,tnons,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  integer,intent(in) :: nqpt
  integer,intent(in) :: nsym
  integer,intent(in) :: rfmeth
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  integer,intent(in) :: blkflg(3,mpert,3,mpert,nblok)
  real(dp),intent(in) :: blknrm(3,nblok)
  real(dp),intent(in) :: blkqpt(9,nblok)
  integer,intent(in) :: blktyp(nblok)
  real(dp),intent(in) :: blkval(2,3*mpert*3*mpert,nblok)
  real(dp),intent(out) :: dynmat(2,3,natom,3,natom,nqpt)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: spqpt(3,nqpt)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine symdm9
end interface

interface
 subroutine symgamma(elph_ds,FSfulltofull,h1_mat_el,&  
  &  indsym,iqptfull,natom,nsym,symq,symrec)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptfull
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  type(elph_type),intent(in) :: elph_ds
  integer,intent(in) :: FSfulltofull(2,nsym,elph_ds%nFSkpt)
  real(dp),intent(inout) :: h1_mat_el(2,elph_ds%nFSband,elph_ds%nFSband, &
  &         elph_ds%nbranch,elph_ds%nFSkpt)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine symgamma
end interface

interface
 subroutine test_ftgkk(elph_ds,gprim,natom,nrpt,rpt,spqpt,wghatm)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  type(elph_type),intent(inout) :: elph_ds
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: spqpt(3,elph_ds%nqpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine test_ftgkk
end interface

interface
 subroutine thm9(acell,amu,anaddb_dtset,atmfrc,dielt,displ,&  
  &  dyewq0,d2cart,eigval,eigvec,gmet,gprim,indsym,iout,mpert,msym,natom,&  
  &  nrpt,nsym,ntypat,phfrq,rcan,rmet,rprim,rpt,symrec,symrel,tcpui,&  
  &  trans,twalli,typat,ucvol,wghatm,xred,zeff,themflag, filname, udispl, ufreq)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in),optional :: themflag
  integer,intent(in),optional :: udispl
  integer,intent(in),optional :: ufreq
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  character(len=fnlen),optional :: filname
  real(dp),intent(in) :: tcpui
  real(dp),intent(in) :: twalli
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(out) :: displ(2*3*natom*3*natom)
  real(dp),intent(inout) :: dyewq0(3,3,natom)
  real(dp),intent(out) :: eigval(3*natom)
  real(dp),intent(out) :: eigvec(2,3,natom,3*natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(out) :: phfrq(3*natom)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: trans(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine thm9
end interface

interface
 subroutine thmeig(acell,amu,blkval2,eigvec,filnam,kpnt,mband,msize,natom,nkpt,nqpt,ntemper,ntypat,&  
  &  phfreq,qphon,rprim,temperinc,tempermin,typat,wtq,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: msize
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nqpt
  integer,intent(in) :: ntemper
  integer,intent(in) :: ntypat
  character(len=fnlen),intent(in) :: filnam
  real(dp),intent(in) :: temperinc
  real(dp),intent(in) :: tempermin
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(in) :: blkval2(2,msize,mband,nkpt,nqpt)
  real(dp),intent(in) :: eigvec(2,3,natom,3*natom,nqpt)
  real(dp),intent(in) :: kpnt(3,nkpt,nqpt)
  real(dp),intent(in) :: phfreq(3*natom,nqpt)
  real(dp),intent(in) :: qphon(9,nqpt)
  real(dp),intent(in) :: rprim(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wtq(3,nqpt)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine thmeig
end interface

interface
 subroutine wght9(brav,gprim,natom,ngqpt,nqpt,nqshft,&  
  &  nrpt,qshft,rcan,rpt,wghatm)
  use defs_basis
  implicit none
  integer,intent(in) :: brav
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: nqshft
  integer,intent(in) :: nrpt
  integer,intent(inout) :: ngqpt(9)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: qshft(3,4)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(out) :: wghatm(natom,natom,nrpt)
 end subroutine wght9
end interface

interface
 subroutine wstoconv(num,red,shift)
  use defs_basis
  implicit none
  real(dp),intent(in) :: num
  real(dp),intent(out) :: red
  real(dp),intent(out) :: shift
 end subroutine wstoconv
end interface

end module interfaces_17ddb
!!***
