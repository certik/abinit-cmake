!{\src2tex{textfont=tt}}
!!****f* ABINIT/rdddb9
!!
!! NAME
!! rdddb9
!!
!! FUNCTION
!! This routine reads the derivative database entirely,
!! for use in ppddb9, and performs some checks and symmetrisation
!! At the end, the whole DDB is in central memory, contained in the
!! array blkval(2,msize,nblok).
!! The information on it is contained in the four arrays
!! blkflg(msize,nblok) : blok flag for each element
!! blkqpt(9,nblok)  : blok wavevector (unnormalized)
!! blknrm(3,nblok)  : blok wavevector normalization
!! blktyp(nblok)    : blok type
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! atifc(natom) = atifc(ia) equals 1 if the analysis of ifc
!!  has to be done for atom ia; otherwise 0
!! ddbun = unit number for DDB io
!! dimekb=dimension of ekb (for the time being, only for norm-
!!                          conserving psps)
!! iout=unit number for output of formatted data
!! filnam(6) names of input or output files
!! mband=maximum number of bands
!! mpert =maximum number of ipert
!! msize=maximum size of data blocks
!! msym =maximum number of symmetry elements in space group
!! natifc = number of atoms for which the analysis of ifc is done
!! natom = number of atoms
!! ntypat=number of atom types
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!! acell(3)=length scales of cell (bohr)
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! blkflg(msize,nblok)= flag of existence for each element of the DDB
!! blknrm(3,nblok)  : blok wavevector normalization
!! blkqpt(9,nblok)  : blok wavevector (unnormalized)
!! blktyp(nblok)    : blok type
!! blkval(2,msize,nblok)= value of each complex element of the DDB
!! gmet(3,3)=reciprocal space metric tensor in bohr**-2
!! gprim(3,3)=dimensionless reciprocal space primitive translations
!! indsym(4,nsym,natom)=indirect indexing array for symmetries
!! natom=number of atoms in cell
!! nblok= number of bloks in the DDB
!! nsym=number of space group symmetries
!! occopt=occupation option
!! rmet(3,3)=metric tensor in real space (bohr^2)
!! rprim(3,3)= primitive translation vectors
!! symq(4,2,nsym)= (integer) three first numbers define the G vector ;
!!   fourth number is zero if the q-vector is not preserved,
!!   second index is about time-reversal symmetry
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!! tnons(3,nsym)=fractional nonsymmorphic translations
!! typat(natom)=type integer for each atom in cell
!! ucvol=unit cell volume in bohr**3
!! xcart(3,natom)=atomic cartesian coordinates
!! xred(3,natom)=fractional dimensionless atomic coordinates
!! zion(ntypat)=charge on each type of atom (real number)
!!
!! NOTES
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      blok8,cart29,chkin9,d2sym3,d3sym,ioddb8,leave_new,mati3inv,matr3inv
!!      metric,mkrdim,nlopt,psddb8,symatm,symq3,timein,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine rdddb9(acell,atifc,amu,blkflg,blknrm,blkqpt,&
& blktyp,blkval,ddbun,dimekb,filnam,gmet,gprim,indsym,iout,&
& mband,mpert,msize,msym,&
& natifc,natom,nblok,nkpt,nsym,ntypat,&
& occopt,rmet,rprim,symq,symrec,symrel,thmflag,&
& tnons,typat,ucvol,usepaw,xcart,xred,zion,blkval2,kpnt)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13recipspace
 use interfaces_16response
 use interfaces_17ddb, except_this_one => rdddb9
!End of the abilint section

 implicit none

!Arguments ------------------------------------
! NOTE: these are used for dimensioning and then re-assigned in ioddb8.
!   This is almost definitely bad practice. In particular
!    it should be indsym(4,msym,natom),
!   and
!    the allocation allocate(kpt(3,nkpt)) is strange
!scalars
 integer,intent(in) :: ddbun,dimekb,iout,mband,mpert,msize,msym,natifc,thmflag
 integer,intent(in) :: usepaw
 integer,intent(inout) :: natom,nblok,nkpt,nsym,ntypat,occopt
 real(dp),intent(out) :: ucvol
 character(len=fnlen),intent(in) :: filnam
!arrays
 integer,intent(inout) :: atifc(natom)
 integer,intent(out) :: blkflg(msize,nblok),blktyp(nblok),indsym(4,nsym,natom)
 integer,intent(out) :: symq(4,2,*),symrec(3,3,msym),symrel(3,3,msym)
 integer,intent(out) :: typat(natom)
 real(dp),intent(out) :: acell(3),amu(ntypat),blknrm(3,nblok),blkqpt(9,nblok)
 real(dp),intent(out) :: blkval(2,msize,nblok),gmet(3,3),gprim(3,3),rmet(3,3)
 real(dp),intent(out) :: rprim(3,3),tnons(3,msym),xcart(3,natom),xred(3,natom)
 real(dp),intent(out) :: zion(ntypat)
 real(dp),intent(out),optional :: blkval2(2,msize,mband,nkpt,nblok)
 real(dp),intent(out),optional :: kpnt(3,nkpt,nblok)

!Local variables -------------------------
!mtyplo=maximum number of type, locally
!scalars
 integer,parameter :: msppol=2,mtyplo=6
 integer :: choice,fullinit,iblok,ii,index,intxc,ipert,iscf,isym,ixc,lmnmax
 integer :: lnmax,mu,nsize,nspden,nspinor,nsppol,nu,nunit,timrev,useylm,vrsddb
 real(dp) :: dilatmx,ecut,ecutsm,kptnrm,sciss,tcpu,tcpui,tolwfr,tphysel,tsmear
 real(dp) :: twall,twalli
 character(len=500) :: message
 character(len=fnlen) :: dummy
!arrays
 integer :: ngfft(18)
 integer,allocatable :: car3flg(:,:,:,:,:,:),carflg(:,:,:,:),indlmn(:,:,:)
 integer,allocatable :: nband(:),pspso(:),symafm(:),tmpflg(:,:,:,:,:,:)
 real(dp) :: gprimd(3,3),qpt(3),rprimd(3,3)
 real(dp),allocatable :: d2cart(:,:,:,:,:),d3cart(:,:,:,:,:,:,:),ekb(:,:)
 real(dp),allocatable :: kpt(:,:),occ(:),spinat(:,:),tmpval(:,:,:,:,:,:,:)
 real(dp),allocatable :: wtk(:),znucl(:)

! *********************************************************************

 call timein(tcpui,twalli)

!Read the DDB information

 vrsddb=010929

!Check the value of usepaw
 if (usepaw==1) then
  write(message, '(a,a,a,a)' ) ch10,&
&  ' psddb8: BUG -',ch10,&
&  '  Paw not yet allowed !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!The checking of pseudopotentials is not done presently
!so that dimensions are fake
 lmnmax=dimekb
 allocate(ekb(dimekb,ntypat),indlmn(6,lmnmax,ntypat),pspso(ntypat))

 allocate(kpt(3,nkpt))
 allocate(nband(nkpt),occ(nkpt*mband*msppol),spinat(3,natom))
 allocate(symafm(msym),wtk(nkpt),znucl(ntypat))

!Open the input derivative database file
!and read the preliminary information
 choice=1
 nunit=ddbun
!Note that in this call, mkpt has been replaced by nkpt,
!mtypat by ntypat, and matom by natom.
!It is suspected that some problem when preparing v3.2 are linked to this
!(the final v3.2 works)
 call ioddb8 (choice,dummy,filnam,natom,mband,&
& nkpt,msym,ntypat,nunit,vrsddb,&
& acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
& natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
& rprim,sciss,spinat,symafm,symrel,tnons,tolwfr,tphysel,tsmear,&
& typat,wtk,xred,zion,znucl)

!Compute different matrices in real and reciprocal space, also
!checks whether ucvol is positive.
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

!Obtain reciprocal space primitive transl g from inverse trans of r
!(Unlike in abinit, gprim is used throughout ifc; should be changed, later)
 call matr3inv(rprim,gprim)

!Generate atom positions in cartesian coordinates
 call xredxcart(natom,1,rprimd,xcart,xred)

!DEBUG
!FORCE ALL SYMMETRIES FOR FCC Al
!nsym = 48
!symrel(:,:,1:nsym) = reshape( (/&
!&  1, 0, 0,    0, 1, 0,   0, 0, 1,     -1, 0, 0  ,  0, -1, 0  , 0, 0, -1    ,  &
!&  0, -1, 1,   0, -1, 0,  1, -1, 0,     0, 1, -1  , 0, 1, 0  , -1, 1, 0    ,  &
!& -1, 0, 0,   -1, 0, 1,  -1, 1, 0,      1, 0, 0  ,  1, 0, -1  , 1, -1, 0    ,  &
!&  0, 1, -1,   1, 0, -1,  0, 0, -1    , 0, -1, 1  ,-1, 0, 1  ,  0, 0, 1    ,  &
!& -1, 0, 0,   -1, 1, 0,  -1, 0, 1    ,  1, 0, 0  ,  1, -1, 0  , 1, 0, -1    ,  &
!&  0, -1, 1,   1, -1, 0,  0, -1, 0    , 0, 1, -1  ,-1, 1, 0  ,  0, 1, 0    ,  &
!&  1, 0, 0,    0, 0, 1,   0, 1, 0    , -1, 0, 0  ,  0, 0, -1  , 0, -1, 0    ,  &
!&  0, 1, -1,   0, 0, -1,  1, 0, -1    , 0, -1, 1  , 0, 0, 1  , -1, 0, 1    ,  &
!& -1, 0, 1,   -1, 1, 0,  -1, 0, 0    ,  1, 0, -1  , 1, -1, 0  , 1, 0, 0    ,  &
!&  0, -1, 0,   1, -1, 0,  0, -1, 1    , 0, 1, 0  , -1, 1, 0  ,  0, 1, -1    ,  &
!&  1, 0, -1,   0, 0, -1,  0, 1, -1    ,-1, 0, 1  ,  0, 0, 1  ,  0, -1, 1    ,  &
!&  0, 1, 0,    0, 0, 1,   1, 0, 0    ,  0, -1, 0  , 0, 0, -1  ,-1, 0, 0    ,  &
!&  1, 0, -1,   0, 1, -1,  0, 0, -1    ,-1, 0, 1  ,  0, -1, 1  , 0, 0, 1    ,  &
!&  0, -1, 0,   0, -1, 1,  1, -1, 0    , 0, 1, 0  ,  0, 1, -1  ,-1, 1, 0    ,  &
!& -1, 0, 1,   -1, 0, 0,  -1, 1, 0    ,  1, 0, -1  , 1, 0, 0  ,  1, -1, 0    ,  &
!&  0, 1, 0,    1, 0, 0,   0, 0, 1    ,  0, -1, 0  ,-1, 0, 0  ,  0, 0, -1    ,  &
!&  0, 0, -1,   0, 1, -1,  1, 0, -1    , 0, 0, 1  ,  0, -1, 1  ,-1, 0, 1    ,  &
!&  1, -1, 0,   0, -1, 1,  0, -1, 0    ,-1, 1, 0  ,  0, 1, -1  , 0, 1, 0    ,  &
!&  0, 0, 1,    1, 0, 0,   0, 1, 0    ,  0, 0, -1  ,-1, 0, 0  ,  0, -1, 0    ,  &
!& -1, 1, 0,   -1, 0, 0,  -1, 0, 1    ,  1, -1, 0  , 1, 0, 0  ,  1, 0, -1    ,  &
!&  0, 0, 1,    0, 1, 0,   1, 0, 0    ,  0, 0, -1  , 0, -1, 0  ,-1, 0, 0    ,  &
!&  1, -1, 0,   0, -1, 0,  0, -1, 1    ,-1, 1, 0  ,  0, 1, 0  ,  0, 1, -1    ,  &
!&  0, 0, -1,   1, 0, -1,  0, 1, -1    , 0, 0, 1  , -1, 0, 1  ,  0, -1, 1    ,  &
!& -1, 1, 0,   -1, 0, 1,  -1, 0, 0    ,  1, -1, 0  , 1, 0, -1  , 1, 0, 0     /), (/3,3,nsym/))
!ENDDEBUG

!Transposed inversion of the symmetry matrices, for use in
!the reciprocal space
 do isym=1,nsym
  call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
 end do

!SYMATM generates for all the atoms and all the symmetries, the atom
!on which the referenced one is sent and also the translation bringing
!back this atom to the referenced unit cell
 call symatm(indsym,natom,nsym,symrec,tnons,typat,xred)

!Read the psp information of the input DDB
 lmnmax=dimekb;lnmax=dimekb;useylm=0
 call psddb8 (choice,dimekb,ekb,fullinit,indlmn,lmnmax,lnmax,&
& nblok,ntypat,nunit,pspso,usepaw,useylm,vrsddb)

!Check the correctness of some input parameters,
!and perform small treatment if needed.
 call chkin9(atifc,natifc,natom)

!Read the blocks from the input database, and close it.
 write(message, '(a,a,a,i5,a)' )ch10,ch10,&
& ' rdddb9 : read ',nblok,' blocks from the input DDB '
 call wrtout(6,message,'COLL')
 choice=1
 nunit=ddbun
 do iblok=1,nblok
! index=1+nsize*(iblok-1)
  if (thmflag==3) then
   call blok8(blkflg(:,iblok),blknrm(:,iblok),blkqpt(:,iblok),&
&   blktyp(iblok),blkval(:,:,iblok),choice,mband,mpert,msize,&
&   natom,nkpt,nunit,blkval2(:,:,:,:,iblok),kpnt(:,:,iblok))
  else
   call blok8(blkflg(:,iblok),blknrm(:,iblok),blkqpt(:,iblok),&
&   blktyp(iblok),blkval(:,:,iblok),choice,mband,mpert,msize,&
&   natom,nkpt,nunit)
  end if

! Here complete the matrix by symmetrisation of the
! existing elements
  if(blktyp(iblok)==1 .or. blktyp(iblok)==2) then
   qpt(1)=blkqpt(1,iblok)/blknrm(1,iblok)
   qpt(2)=blkqpt(2,iblok)/blknrm(1,iblok)
   qpt(3)=blkqpt(3,iblok)/blknrm(1,iblok)

!  Examine the symmetries of the q wavevector
   call symq3(nsym,qpt,symq,symrec,timrev)

   nsize=3*mpert*3*mpert
   allocate(tmpflg(3,mpert,3,mpert,1,1))
   allocate(tmpval(2,3,mpert,3,mpert,1,1))
   tmpflg(:,:,:,:,1,1) = reshape(blkflg(1:nsize,iblok),&
&   shape = (/3,mpert,3,mpert/))
   tmpval(1,:,:,:,:,1,1) = reshape(blkval(1,1:nsize,iblok),&
&   shape = (/3,mpert,3,mpert/))
   tmpval(2,:,:,:,:,1,1) = reshape(blkval(2,1:nsize,iblok),&
&   shape = (/3,mpert,3,mpert/))

!  Then apply symmetry operations
   call d2sym3(tmpflg,tmpval,indsym,mpert,&
&   natom,nsym,qpt,symq,symrec,symrel,timrev)

!  Transform the dynamical matrix in cartesian coordinates
   allocate(carflg(3,mpert,3,mpert),d2cart(2,3,mpert,3,mpert))
   call cart29(tmpflg,tmpval,carflg,d2cart,&
&   gprimd,1,mpert,natom,1,nsize,ntypat,rprimd,typat,ucvol,zion)

   blkflg(1:nsize,iblok) = reshape(carflg,shape = (/3*mpert*3*mpert/))
   blkval(1,1:nsize,iblok) = reshape(d2cart(1,:,:,:,:),&
&   shape = (/3*mpert*3*mpert/))
   blkval(2,1:nsize,iblok) = reshape(d2cart(2,:,:,:,:),&
&   shape = (/3*mpert*3*mpert/))


   deallocate(carflg,d2cart)
   deallocate(tmpflg,tmpval)

  else if (blktyp(iblok) == 3) then

   nsize=3*mpert*3*mpert*3*mpert
   allocate(tmpflg(3,mpert,3,mpert,3,mpert))
   allocate(tmpval(2,3,mpert,3,mpert,3,mpert))

   tmpflg(:,:,:,:,:,:) = reshape(blkflg(1:nsize,iblok),&
&   shape = (/3,mpert,3,mpert,3,mpert/))
   tmpval(1,:,:,:,:,:,:) = reshape(blkval(1,1:nsize,iblok),&
&   shape = (/3,mpert,3,mpert,3,mpert/))
   tmpval(2,:,:,:,:,:,:) = reshape(blkval(2,1:nsize,iblok),&
&   shape = (/3,mpert,3,mpert,3,mpert/))

   call d3sym(tmpflg,tmpval,indsym,mpert,natom,nsym,&
&   symrec,symrel)

   allocate(d3cart(2,3,mpert,3,mpert,3,mpert))
   allocate(car3flg(3,mpert,3,mpert,3,mpert))
   call nlopt(tmpflg,car3flg,tmpval,d3cart,gprimd,mpert,natom,rprimd,ucvol)

   blkflg(1:nsize,iblok) = reshape(car3flg, shape = (/3*mpert*3*mpert*3*mpert/))
   blkval(1,1:nsize,iblok) = reshape(d3cart(1,:,:,:,:,:,:),&
&   shape = (/3*mpert*3*mpert*3*mpert/))
   blkval(2,1:nsize,iblok) = reshape(d3cart(2,:,:,:,:,:,:),&
&   shape = (/3*mpert*3*mpert*3*mpert/))


   deallocate(d3cart,car3flg)
   deallocate(tmpflg,tmpval)
  end if

 end do

 close(ddbun)

 write(message,'(a)' )' Now the whole DDB is in central memory '
 call wrtout(6,message,'COLL')
 call wrtout(iout,message,'COLL')

 deallocate(ekb,indlmn,kpt,nband,occ,pspso,spinat,symafm,wtk,znucl)

!DEBUG
!write(6,*)' rdddb9 : exit'
!ENDDEBUG

end subroutine rdddb9
!!***
