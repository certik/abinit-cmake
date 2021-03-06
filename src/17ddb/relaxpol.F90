!{\src2tex{textfont=tt}}
!!****f* ABINIT/relaxpol
!!
!! NAME
!! relaxpol
!!
!! FUNCTION
!! 1) Compute polarization in cartesian coordinates
!! 2) Structural relaxation at fixed polarization: this routine
!!    solves the linear system of equations Eq.(13)
!!    of Na Sai et al., PRB 66, 104108 (2002).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! blkflg(msize) = flag for every matrix element (0=> the element
!!   is not in the data block), (1=> the element is in the data blok)
!! blkval(2,msize) = matrix that contains the second-order energy
!!   derivatives
!! etotal = Kohn-Sham energy at zero electric field
!! fred(3,natom) = -1 times the forces in reduced coordinates
!! iatfix(natom) = indices of the atoms that are held fixed in the relaxation
!! indsym(4,msym,natom) = indirect indexing array for atom labels
!! iout = unit number for output
!! istrfix(6) = indices of the elements of the strain tensor that
!!   are held fixed in the relaxation
!!      1 = xx
!!      2 = yy
!!      3 = zz
!!      4 = yz & zy
!!      5 = xz & zx
!!      6 = xy & yx
!! mpert = maximum number of ipert
!! msize = dimension of blkflg and blkval
!! msym = maximal number of symmetries
!! natfix = number of atoms that are held fixed in the relaxation
!! natom = number of atoms in the unit cell
!! nstrfix = number of elements of the strain tensor that
!!   are held fixed in the relaxation
!! nsym = number of symmetry elements in space group
!! ntypat = number of atom types
!! pel(3) = electronic polarization not taking into account
!!   the factor 1/ucvol
!! relaxat = 1: relax atomic positions
!!         = 0: do not relax atomic positions
!! relaxstr = 1: relax cell parameters
!!          = 0: do not relax cell parameters
!! rprimd(3,3) = dimensional primitive vectors
!! strten(6) = stress tensor in cartesian coordinates
!! symrel(3,3,nsym) = symmetry operations in real space in terms
!!   of primitive translations
!! targetpol(3) = target value of the polarization
!! typat(natom) = type of each atom
!! ucvol = unit cell volume in bohr**3
!! xcart(3,natom) = cartesian coordinates of the atoms in the unit cell
!! xred(3,natom) = reduced coordinates of the atoms in the unit cell
!! zion(ntypat) = charge on each type of atom
!!
!! OUTPUT
!!
!! NOTES
!! - The elements of the dynamical matrix stored in blkval
!!   are symmetrized before computing the new atomic positions and
!!   cell parameters.
!! - In case relaxat = 0 and relaxstr = 0, the routine only
!!   computes the polarization in cartesian coordinates.
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      leave_new,matr3inv,polcart,symdyma,wrtout,xredxcart,dzgedi,dzgefa
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine relaxpol(blkflg,blkval,etotal,fred,iatfix,indsym,iout,istrfix,&
& mpert,msize,msym,natfix,natom,nstrfix,nsym,ntypat,pel,relaxat,relaxstr,&
& rprimd,strten,symrel,targetpol,typat,ucvol,xcart,xred,zion)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_16response
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iout,mpert,msize,msym,natfix,natom,nstrfix,nsym,ntypat
 integer,intent(in) :: relaxat,relaxstr
 real(dp),intent(in) :: etotal,ucvol
!arrays
 integer,intent(in) :: blkflg(msize),iatfix(natom),indsym(4,msym,natom)
 integer,intent(in) :: istrfix(6),symrel(3,3,nsym),typat(natom)
 real(dp),intent(in) :: fred(3,natom),pel(3),rprimd(3,3),strten(6)
 real(dp),intent(in) :: xcart(3,natom),xred(3,natom),zion(ntypat)
 real(dp),intent(inout) :: blkval(2,msize),targetpol(3)

!Local variables -------------------------
!scalars
 integer :: flag,iatom,idir,ii,index,index1,index_tild,info,ipert,istrain
 integer :: itypat,jatom,jdir,jj,job,jpert,option,polunit,posi,posj,sizef
 real(dp) :: e1,fmax,poltmp,sigmax,tol,value
 character(len=500) :: message
!arrays
 integer :: irelaxstrain(6)
 integer,allocatable :: ipvt(:),irelaxat(:),rfpert(:,:)
 real(dp) :: acell_new(3),delta_eta(6),delta_xcart(3,natom),det(2,2),diffpol(3)
 real(dp) :: diffsig(6),favg(3),gprimd(3,3),lambda(3),pel_cart(3),pelev(3)
 real(dp) :: pion(3),pion_cart(3),pol(3),ptot_cart(3),qphon(3),rprim(3,3)
 real(dp) :: rprim_new(3,3),rprimd_new(3,3),sigelfd(6),strainmat(3,3)
 real(dp) :: xcart_new(3,natom),xred_new(3,natom)
 real(dp),allocatable :: cfac(:,:),delta(:),dymati(:),fcart(:,:),fcmat(:,:,:)
 real(dp),allocatable :: fdiff(:,:),felfd(:,:),ifcmat(:,:,:),vec(:),zgwork(:,:)

! *********************************************************************

!DEBUG
!write(6,*)'relaxpol: enter'
!write(6,*)'targetpol = ',targetpol(:)
!do ii = 1, nsym
!write(6,*)symrel(:,:,ii)
!end do
!stop
!ENDDEBUG

!Check if some degrees of freedom remain fixed during the
!optimization

 allocate(irelaxat(natom))
 irelaxat(:) = 1   ; irelaxstrain(:) = 1
 if (natfix > 0) then
  do ii = 1, natfix
   iatom = iatfix(ii)
   if ((iatom > natom).or.(iatom < 0)) then
    write(message, '(a,a,a,i4,a,i4,a,a,a,a,a)')&
&    ' relaxpol : ERROR -',ch10,&
&    '  The value of iatfix(',ii,') is',iatom,', which is not allowed.',ch10,&
&    '  iatfix must be larger than 0 and smaller than natom.',ch10,&
&    '  Action : correct iatfix in your input file.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   irelaxat(iatom) = 0
  end do
 end if

 if (nstrfix > 0) then
  do ii = 1, nstrfix
   istrain = istrfix(ii)
   if ((istrain > 6).or.(istrain < 0)) then
    write(message, '(a,a,a,i4,a,i4,a,a,a,a,a)')&
&    ' relaxpol : ERROR -',ch10,&
&    '  istrfix(',ii,') is',istrain,', which is not allowed.',ch10,&
&    '  istrfix must be larger than 0 and smaller than 6.',ch10,&
&    '  Action : correct istrfix in your input file.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   irelaxstrain(istrain) = 0
  end do
 end if


 allocate(rfpert(mpert,3),cfac(mpert,mpert))
 call matr3inv(rprimd,gprimd)

!Compute the size of the matrix that contains the second-order derivatives

 sizef = 3
 rfpert(:,:) = 0
 rfpert(natom+2,1:3) = 1
 if (relaxat == 1) then
  do iatom = 1, natom
   if (irelaxat(iatom) == 1) then
    sizef = sizef + 3
    rfpert(iatom,1:3) = 1
   end if
  end do
 end if
 ii = natom + 2
 if (relaxstr == 1) then
  istrain = 0
  do ipert = (natom+3), (natom+4)
   do idir = 1, 3
    istrain = istrain + 1
    if (irelaxstrain(istrain) == 1) then
     sizef = sizef + 1
     rfpert(ipert,idir) = 1
    end if
   end do
  end do
 end if


 allocate(fcmat(2,sizef,sizef),ifcmat(2,sizef,sizef),vec(sizef))
 allocate(delta(sizef),ipvt(sizef),zgwork(2,sizef))
 allocate(fcart(3,natom),felfd(3,natom),fdiff(3,natom))

!Build the vector that stores the forces, sigma and the polarization

 vec(:) = 0._dp
 posi = 0

 if (relaxat == 1) then

! Note conversion to cartesian coordinates (bohr) AND
! negation to make a force out of a gradient
! Also subtract off average force from each force
! component to avoid spurious drifting of atoms across cell.
  favg(:) = zero
  do iatom = 1, natom
   do idir = 1, 3
    fcart(idir,iatom) = -(gprimd(idir,1)*fred(1,iatom) + &
&    gprimd(idir,2)*fred(2,iatom) + &
&    gprimd(idir,3)*fred(3,iatom))
    favg(idir) = favg(idir) + fcart(idir,iatom)
   end do
  end do
  favg(:) = favg(:)/dble(natom)
  do iatom = 1, natom
   fcart(:,iatom) = fcart(:,iatom) - favg(:)
  end do

  do iatom = 1, natom
   if (irelaxat(iatom) == 1) then
    do idir = 1, 3
     posi = posi + 1
     vec(posi) = fcart(idir,iatom)
    end do
   end if
  end do

 end if    ! relaxat == 1

!DEBUG
!write(*,*)'Forces in cartesian coords'
!do iatom = 1, natom
!write(*,'(3(2x,e16.9))')(fcart(idir,iatom),idir = 1, 3)
!end do
!stop
!ENDDEBUG

!Transform target polarization to atomic units
 targetpol(:) = targetpol(:)*((Bohr_Ang*1.0d-10)**2)/e_Cb

!Compute ionic polarization
 pion(:) = zero
 do iatom = 1, natom
  itypat = typat(iatom)
  do idir = 1, 3
   poltmp = zion(itypat)*xred(idir,iatom)
   poltmp = poltmp - two*nint(poltmp/two)   ! fold into [-1,1]
   pion(idir) = pion(idir) + poltmp
  end do
 end do
 do idir = 1, 3
  pion(idir) = pion(idir) - two*nint(pion(idir)/two) ! fold into [-1,1]
 end do

!Transform the polarization to cartesian coordinates
 polunit = 3
 pelev=zero
 call polcart(pel,pel_cart,pelev,pion,pion_cart,polunit,&
 ptot_cart,rprimd,ucvol,iout)

 do idir = 1, 3
  posi = posi + 1
  vec(posi) = ptot_cart(idir) - targetpol(idir)
 end do


 if (relaxstr == 1) then
  do istrain = 1, 6
   if (irelaxstrain(istrain) == 1) then
    posi = posi + 1
    vec(posi) = -1._dp*strten(istrain)*ucvol
   end if
  end do
 end if


!Symmetrize the dynamical matrix

 allocate(dymati(2*3*natom*3*natom))  ! A vector of this shape is required
!by the symdyma routine
 do ipert = 1, natom
  do idir = 1, 3
   do jpert = 1, natom
    do jdir = 1, 3
     index  = jdir +3*((jpert - 1) + mpert*((idir - 1) + 3*(ipert - 1)))
     index1 = jdir +3*((jpert - 1) + natom*((idir - 1) + 3*(ipert - 1)))
     dymati(2*index1 - 1) = blkval(1,index)
     dymati(2*index1    ) = blkval(2,index)
    end do
   end do
  end do
 end do

 qphon(:) = zero
 call symdyma(dymati,indsym,msym,natom,nsym,qphon,rprimd,symrel,xred)

 do ipert = 1, natom
  do idir = 1, 3
   do jpert = 1, natom
    do jdir = 1, 3
     index  = jdir +3*((jpert - 1) + mpert*((idir - 1) + 3*(ipert - 1)))
     index1 = jdir +3*((jpert - 1) + natom*((idir - 1) + 3*(ipert - 1)))
     blkval(1,index) = dymati(2*index1 - 1)
     blkval(2,index) = dymati(2*index1    )
    end do
   end do
  end do
 end do

 deallocate(dymati)

!Define conversion factors for blkval
 cfac(:,:) = 1._dp
 cfac(1:natom,natom+2) = -1._dp/ucvol
 cfac(natom+2,1:natom) = -1._dp/ucvol
 cfac(natom+3:natom+4,natom+2) = -1._dp
 cfac(natom+2,natom+3:natom+4) = -1._dp


!Build the matrix that contains the second-order derivatives
!ipert = natom + 1 corresponds to the ddk perturbation, that
!is not needed; so skip it

 fcmat(:,:,:) = zero

 posi = 0
 flag = 0  ! When fcmat has been build, flag = 0 if all elements
!were available. Otherwise, it will be 1.
!In case one element is missing, check if
!it can be obtained by changing the order of the
!perturbations

 do ipert = 1, mpert
  do idir = 1, 3
   if (rfpert(ipert,idir) == 1) then
    posi = posi + 1
    posj = 0

    do jpert = 1, mpert
     do jdir = 1, 3
      if (rfpert(jpert,jdir) == 1) then
       index = jdir +3*((jpert - 1) + mpert*((idir - 1) + 3*(ipert - 1)))
       index_tild = idir +3*((ipert - 1) + mpert*((jdir - 1) + 3*(jpert - 1)))
       posj = posj + 1
       if ((ipert /= natom + 2).or.(jpert /= natom + 2)) then
        if (blkflg(index) == 1) then
         fcmat(:,posi,posj) = blkval(:,index)*cfac(ipert,jpert)
        else if (blkflg(index_tild) == 1) then
         fcmat(:,posi,posj) = blkval(:,index_tild)*cfac(ipert,jpert)
         blkval(:,index) = blkval(:,index_tild)
        else
         flag = 1
         write(6,'(a,4(2x,i3))')'relaxpol: could not find element:',&
&         idir,ipert,jdir,jpert
        end if
       end if
!      DEBUG
!      write(100,'(4(2x,i3),5x,f16.9)')idir,ipert,jdir,jpert,&
!      & fcmat(1,posi,posj)
!      ENDDEBUG
      end if
     end do
    end do

   end if
  end do
 end do

 if (flag == 1) then
  write(message, '(a,a,a,a,a,a,i2,a,i2,a,a,a,a)' )ch10,&
&  ' relaxpol : ERROR -',ch10,&
&  '  Some of the second order derivatives required to deal with the case',ch10,&
&  '  relaxat = ',relaxat,', relaxstr = ', relaxstr, ch10,&
&  '  are missing in the DDB.',ch10,&
&  '  Action : correct your DDB or change your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if


!Compute the inverse of the force constant matrix

 if ((relaxat /= 0).or.(relaxstr /= 0)) then

  job = 1          ! compute inverse only
  ifcmat(:,:,:) = fcmat(:,:,:)
  call dzgefa(ifcmat,sizef,sizef,ipvt,info)
  call dzgedi(ifcmat,sizef,sizef,ipvt,det,zgwork,job)

! DEBUG
! write(100,*)'relaxat = ',relaxat
! write(100,*)'relaxstr = ',relaxstr
! write(100,*)'irelaxat = '
! write(100,*)irelaxat(:)
! write(100,*)'irelaxstrain = '
! write(100,*)irelaxstrain(:)
! write(100,*)'sizef = ',sizef
! write(100,*)'targetpol ='
! write(100,*)targetpol(:)
! do ipert = 1, sizef
! do jpert = 1, sizef
! write(100,'(2(2x,i3),2x,e16.9)')ipert,jpert,fcmat(1,ipert,jpert)
! end do
! end do
! stop
! ENDDEBUG

! Compute \delta R, \delta \eta and \lambda
  delta(:) = 0._dp
  do ipert = 1, sizef
   do jpert = 1, sizef
    delta(ipert) = delta(ipert) + ifcmat(1,ipert,jpert)*vec(jpert)
   end do
  end do


! Update atomic positions
  posi = 0
  if (relaxat == 1) then

   delta_xcart(:,:) = 0._dp
   xcart_new(:,:) = 0._dp
   do iatom = 1, natom
    if (irelaxat(iatom) == 1) then
     do idir = 1, 3
      posi = posi + 1
      delta_xcart(idir,iatom) = delta(posi)
     end do
    end if
   end do

!  Drop unsignificant digits in order to eleminate numerical noise
   tol = 10000000._dp
   do iatom = 1, natom
    do idir = 1, 3
     value = delta_xcart(idir,iatom)
     ii = log10(abs(value))
     if (ii <= 0) then
      ii = abs(ii) + 1
      value = 1._dp*int(tol*value*10**ii)/(tol*10**ii)
     else
      value = 1._dp*int(tol*value/(10**ii))*(10**ii)/tol
     end if
     delta_xcart(idir,iatom) = value
    end do
   end do

   xcart_new(:,:) = xcart(:,:) + delta_xcart(:,:)
   option = -1
   call xredxcart(natom,option,rprimd,xcart_new,xred_new)

  end if           ! relaxat == 1

! Compute lambda and the value of the energy functional
! $F - \lambda \cdot P$

  e1 = etotal
  do idir = 1, 3
   posi = posi + 1
   lambda(idir) = delta(posi)
   e1 = e1 - lambda(idir)*ptot_cart(idir)
  end do


! Update cell parameters

  if (relaxstr == 1) then

   delta_eta(:) = 0._dp
   do istrain = 1, 6
    if (irelaxstrain(istrain) == 1) then
     posi = posi + 1
     delta_eta(istrain) = delta(posi)
    end if
   end do

   do istrain = 1, 3
    strainmat(istrain,istrain) = delta_eta(istrain)
   end do
   strainmat(2,3) = delta_eta(4)/2._dp ; strainmat(3,2) = delta_eta(4)/2._dp
   strainmat(1,3) = delta_eta(5)/2._dp ; strainmat(3,1) = delta_eta(5)/2._dp
   strainmat(2,1) = delta_eta(6)/2._dp ; strainmat(1,2) = delta_eta(6)/2._dp

   rprimd_new(:,:) = 0._dp
   do idir = 1, 3
    do jdir = 1, 3
     do ii = 1, 3
      rprimd_new(jdir,idir) = rprimd_new(jdir,idir) + &
&      rprimd(ii,idir)*strainmat(ii,jdir)
     end do
    end do
   end do
   rprimd_new(:,:) = rprimd_new(:,:) + rprimd(:,:)

   acell_new(:) = 0._dp
   do idir = 1, 3
    do jdir = 1, 3
     acell_new(idir) = acell_new(idir) + &
&     rprimd_new(jdir,idir)*rprimd_new(jdir,idir)
    end do
    acell_new(idir) = sqrt(acell_new(idir))
    rprim(:,idir) = rprimd_new(:,idir)/acell_new(idir)
   end do

  end if          ! relaxstr == 1

! Write out the results

  write(iout,*)
  write(iout,'(a,80a,a)') ch10,('=',ii=1,80),ch10
  write(iout,*)
  write(iout,*)'Relaxation of the geometry at fixed polarization:'
  write(iout,*)

  write(iout,'(a,3(2x,f16.9))')' Lambda = ',(lambda(idir),idir = 1, 3)
  write(iout,'(a,e16.9)')' Value of the energy functional E_1 = ',e1
  write(iout,*)
  write(iout,*)'Difference between actual value of the Polarization (C/m^2)'
  write(iout,*)'and the target value:'
  diffpol(:) = (ptot_cart(:) - &
&  targetpol(:))*e_Cb/((Bohr_Ang*1.0d-10)**2)
  write(iout,'(3(3x,f16.9))')(diffpol(idir),idir = 1, 3)



  if (relaxat == 1) then

!  Compute the forces induced on the atoms by the electric field
!  The strength of the field is determined by lambda
   felfd(:,:) = zero
   do iatom = 1, natom
    do idir = 1, 3
     do jdir = 1, 3
      index = idir +3*((iatom - 1) + mpert*((jdir - 1) + 3*(natom + 1)))
      felfd(idir,iatom) = felfd(idir,iatom) - &
&      lambda(jdir)*blkval(1,index)/ucvol
     end do
    end do
   end do

!  Compute remaining forces and write them out

   fdiff(:,:) = fcart(:,:) - felfd(:,:)
   write(iout,*)
   write(iout,*)'Difference between the Hellmann-Feynman forces'
   write(iout,*)'and the forces induced by the electric field'
   write(iout,*)'(cartesian coordinates, hartree/bohr)'
   fmax = zero
   do iatom = 1, natom
    write(iout,'(3(3x,es16.9))')(fdiff(idir,iatom),idir = 1, 3)
    do idir = 1, 3
     if (abs(fdiff(idir,iatom)) > fmax) fmax = abs(fdiff(idir,iatom))
    end do
   end do
   write(iout,'(a,3x,es16.9)')' fmax = ',fmax

   write(iout,*)
   write(iout,*)'Change of cartesian coordinates (delta_xcart):'
   do iatom = 1, natom
    write(iout,'(5x,i3,3(2x,f16.9))')iatom,&
&    (delta_xcart(idir,iatom),idir = 1, 3)
   end do
   write(iout,*)
   write(iout,*)'New cartesian coordinates (xcart_new):'
   write(iout,*)'  xcart'
   do iatom = 1, natom
    write(iout,'(3(3x,d22.14))')(xcart_new(idir,iatom),idir = 1, 3)
   end do
   write(iout,*)
   write(iout,*)'New reduced coordinates (xred_new):'
   write(iout,*)'  xred'
   do iatom = 1, natom
    write(iout,'(3(3x,d22.14))')(xred_new(idir,iatom),idir = 1, 3)
   end do

  end if         ! relaxat == 1

  if (relaxstr == 1) then

!  Compute the stresses induced by the electric field
   sigelfd(:) = zero
   istrain = 0
   do ipert = 1, 2
    jpert = natom + 2 + ipert
    do idir = 1, 3
     istrain = istrain + 1
     do jdir = 1, 3
      index = idir +3*((jpert - 1) + mpert*((jdir - 1) + 3*(natom + 1)))
      sigelfd(istrain) = sigelfd(istrain) + &
&      lambda(jdir)*blkval(1,index)
     end do
     sigelfd(istrain) = sigelfd(istrain)/ucvol
    end do
   end do

!  Compute the remaining stresses and write them out
   diffsig(:) = strten(:) - sigelfd(:)
   sigmax = zero
   do istrain = 1, 6
    if (abs(diffsig(istrain)) > sigmax) sigmax = abs(diffsig(istrain))
   end do
   write(iout,*)
   write(iout,*)'Difference between the Hellmann-Feynman stresses'
   write(iout,*)'and the stresses induced by the electric field'
   write(iout,*)'(cartesian coordinates, hartree/bohr^3)'
   write(iout,'(2x,a,f16.9,5x,a,f16.9)')'diffsig(1) = ',diffsig(1),&
&   'diffsig(4) = ',diffsig(4)
   write(iout,'(2x,a,f16.9,5x,a,f16.9)')'diffsig(2) = ',diffsig(2),&
&   'diffsig(5) = ',diffsig(5)
   write(iout,'(2x,a,f16.9,5x,a,f16.9)')'diffsig(3) = ',diffsig(3),&
&   'diffsig(6) = ',diffsig(6)
   write(iout,'(a,3x,es16.9)')' sigmax = ',sigmax

!  DEBUG
!  do istrain = 1, 6
!  write(*,'(2(5x,e16.9))')strten(istrain),sigelfd(istrain)
!  end do
!  stop
!  ENDDEBUG

   write(iout,*)
   write(iout,*)'Induced strain (delta_eta):'
   write(iout,'(2x,a,f16.9,5x,a,f16.9)')'delta_eta(1) = ',delta_eta(1),&
&   'delta_eta(4) = ',delta_eta(4)
   write(iout,'(2x,a,f16.9,5x,a,f16.9)')'delta_eta(2) = ',delta_eta(2),&
&   'delta_eta(5) = ',delta_eta(5)
   write(iout,'(2x,a,f16.9,5x,a,f16.9)')'delta_eta(3) = ',delta_eta(3),&
&   'delta_eta(6) = ',delta_eta(6)

   write(iout,*)
   write(iout,*)'New lattice constants (acell_new):'
   write(iout,*)'  acell'
   write(iout,'(3(2x,d22.14))')(acell_new(idir),idir = 1, 3)

   write(iout,*)
   write(iout,*)'New primitive vectors (rprim_new):'
   write(iout,*)'  rprim'
   write(iout,'(3(2x,d22.14))')(rprim(idir,1),idir = 1, 3)
   write(iout,'(3(2x,d22.14))')(rprim(idir,2),idir = 1, 3)
   write(iout,'(3(2x,d22.14))')(rprim(idir,3),idir = 1, 3)

  end if         ! relaxstr /= 0

 end if    !  (relaxat /= 0).or.(relaxstr /= 0)

 deallocate(cfac,fdiff,felfd,delta,fcart,fcmat,ifcmat,ipvt,rfpert,vec,zgwork)
 deallocate(irelaxat)

end subroutine relaxpol
!!***
