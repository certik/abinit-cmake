!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkphdos
!!
!! NAME
!! mkphdos
!!
!! FUNCTION
!! This routine calculates the phonon density of states as well as 
!! the contributions associated to the different type of atoms in the unit cell.
!! Two methods are implemented: gaussian method and linear interpolation based on 
!! tetrahedrons.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3) =length scales by which rprim is to be multiplied
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! anaddb_dtset= (derived datatype) contains all the input variables
!! atmfrc(2,3,natom,3,natom,nrpt) = Interatomic Forces in real space
!! dielt(3,3)=dielectric tensor
!! dyewq0(3,3,natom)=Ewald part of the dynamical matrix, at q=0.
!! filnam=prefix for the output file containing the PHDOS, presently not used 
!! gmet(3,3)= metric tensor in reciprocal space.
!! gprim(3,3)=normalized coordinates in reciprocal space
!! indsym(4,nsym,natom)=label given by subroutine symatm, indicating atom
!!  label which gets rotated into given atom by given symmetry
!!  (first three elements are related primitive translation--
!!  see symatm where this is computed)
!! mpert =maximum number of ipert
!! msym = maximum number of symmetries
!! natom=number of atoms in the unit cell
!! nrpt=number of R points in the Big Box
!! nsym=number of symmetries
!! ntypat=number of atom types
!! rcan(3,natom)=atomic position in canonical coordinates
!! rmet(3,3)=metric tensor in real space.
!! rprim(3,3)=dimensionless primitive translations in real space
!! rpt(3,nprt)=canonical coordinates of the R points in the unit cell
!!  These coordinates are normalized (=> * acell(3)!!)
!! symrec(3,3,nsym)=symmetry operations
!! symrel(3,3,nsym)=symmetry operations
!! tcpui=initial cpu time
!! trans(3,natom)=atomic translations : xred = rcan + trans
!! twalli=initial wall clock time
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume
!! wghatm(natom,natom,nrpt)=weights associated to a pair of atoms and to a R vector
!! xred(3,natom)= relative coords of atoms in unit cell (dimensionless)
!! zeff(3,3,natom)=effective charge on each atom, versus electric
!!  field and atomic displacement
!!
!! OUTPUT
!!
!! NOTES
!! 
!! On the use of the q-grids : 
!! Two different q-meshes are used in this subroutine. The first one is the coarse 
!! mesh where the interatomic forces have been calculated during the DFPT run. 
!! This q-grid is used to obtain an initial guess for the max and min frequency 
!! value of the phonon spectrum. These values are, indeed, required to dimension 
!! the array containing the PHDOS. The second (dense) grid is used to perform the 
!! PHDOS calculation. If the Fourier interpolation on the secons dense q-grid 
!! generates a phonon frequency outside the initially calculated frequency mesh,
!! the mesh is enlarged and the calculation is restarted.
!!
!! PARENTS
!!      
!!
!! CHILDREN
!!      
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mkphdos(acell,amu,anaddb_dtset,atmfrc,dielt,dyewq0,filname,gmet,gprim,indsym,&
& mpert,msym,natom,nrpt,nsym,ntypat,rcan,rmet,rprim,rpt,symrec,symrel,tcpui,  &
& trans,twalli,typat,ucvol,wghatm,xred,zeff)

 use defs_basis
 use defs_datatypes
 use m_IO_tools, only : get_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13recipspace
 use interfaces_14occeig
 use interfaces_16response
 use interfaces_17ddb, except_this_one => mkphdos
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,msym,natom,nrpt,nsym,ntypat
 real(dp),intent(in) :: tcpui,twalli,ucvol
 character(len=fnlen),intent(in) :: filname
 type(anaddb_dataset_type),intent(in) :: anaddb_dtset
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrec(3,3,nsym),symrel(3,3,nsym)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: acell(3),amu(ntypat),dielt(3,3),gmet(3,3),gprim(3,3)
 real(dp),intent(in) :: rcan(3,natom),rmet(3,3),rprim(3,3),rpt(3,nrpt)
 real(dp),intent(in) :: trans(3,natom),wghatm(natom,natom,nrpt),xred(3,natom)
 real(dp),intent(in) :: zeff(3,3,natom)
 real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt),dyewq0(3,3,natom)

!Local variables -------------------------
 character(len=50),parameter :: sub_name='mkphdos.F90'
!scalars
 integer :: facbrv,iat,idir,ii,imesh,imode,io,iq,istat,itype,mtetra,nkpt_fullbz
 integer :: nmesh,nomega,nqbz,nqibz,nqpt_max,nqshft,ntetra_ibz,option,timrev
 integer :: unt
 real(dp),parameter :: tol5=0.00001_dp
 real(dp) :: Ha_meV,bzvol,cfact,dossmear,dum,gaussfactor,gaussprefactor
 real(dp) :: gaussval,lbound,max_occ,norm,omega_max,omega_min,omega_step,pnorm
 real(dp) :: qphnrm,tcpu,twall,ubound,vtetra,xx
 logical :: out_of_bounds
 character(len=3) :: unitname
 character(len=50) :: frmt
 character(len=500) :: msg
 character(len=fnlen) :: fnam_loc
!arrays
 integer :: igqpt2(3),qptrlatt(3,3)
 integer,allocatable :: bz2ibz(:),ibz2bz(:),ngqpt(:,:),tetra_full(:,:,:)
 integer,allocatable :: tetra_mult(:),tetra_wrap(:,:,:)
 real(dp) :: d2cart(2,3,mpert,3,mpert),displ(2*3*natom*3*natom),eigval(3*natom)
 real(dp) :: eigvec(2,3,natom,3*natom),gprimd(3,3),kpq(3,1),phfrq(3*natom)
 real(dp) :: qlatt(3,3),qphon(3),rlatt(3,3),rprimd(3,3)
 real(dp),allocatable :: dtweightde(:,:),full_eigvec(:,:,:,:,:),full_phfrq(:,:)
 real(dp),allocatable :: kpt_fullbz(:,:),omega(:),phdos(:),phdos_int(:)
 real(dp),allocatable :: pjdos(:,:,:),pjdos_int(:,:,:),pjdos_typ(:,:)
 real(dp),allocatable :: pjdos_typ_int(:,:),pjdos_xyz_typ(:,:,:),qbz(:,:)
 real(dp),allocatable :: qibz(:,:),qshft(:,:),tmp_phfrq(:),tweight(:,:),wtq(:)
 real(dp),allocatable :: wtq_folded(:),wtqibz(:)

! *********************************************************************

#ifdef DEBUG_MODE 
 write(std_out,*)' mkphdos : enter '
!call flush(std_out)
#endif

 if (anaddb_dtset%prtdos/=1.and.anaddb_dtset%prtdos/=2) then 
  write(msg,'(6a,i5)')ch10,&
&  ' mkphdos : BUG -',ch10,&
&  ' The argument anaddb%prtdos should be 1 or 2 ',ch10,&
&  ' however prtdos = ',anaddb_dtset%prtdos 
  call wrtout(std_out,msg,'COLL') 
  call leave_new('COLL')
 end if 
 if (anaddb_dtset%dosdeltae<=zero) then 
  write(msg,'(6a,es16.8)')ch10,&
&  ' mkphdos : BUG -',ch10,&
&  ' The argument anaddb%dosdeltae should be positive ',ch10,&
&  ' however dosdeltae = ',anaddb_dtset%dosdeltae 
  call wrtout(std_out,msg,'COLL') 
  call leave_new('COLL')
 end if 
 if (anaddb_dtset%prtdos==1.and.anaddb_dtset%dossmear<=zero) then 
  write(msg,'(6a,es16.8)')ch10,&
&  ' mkphdos : BUG -',ch10,&
&  ' The argument anaddb%dossmear should be positive ',ch10,&
&  ' however dossmear = ',anaddb_dtset%dosdeltae 
  call wrtout(std_out,msg,'COLL') 
  call leave_new('COLL')
 end if 
!
!If timing is needed at some point ...
 call timein(tcpu,twall)
 write(msg,'(a,2(a,f11.4),a)')ch10,&
& ' mkphdos : begin at tcpu',tcpu-tcpui,' and twall',twall-twalli,' sec'
 call wrtout(std_out,msg,'COLL')

 call mkrdim(acell,rprim,rprimd)
 call matr3inv(rprimd,gprimd)
 bzvol=ABS ( gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
& -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
& +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))
!
!=== Parameters defining the gaussian approximant ===
 Ha_mev=Ha_eV*1000
 omega_step=anaddb_dtset%dosdeltae

 if (anaddb_dtset%prtdos==1) then 
  dossmear=anaddb_dtset%dossmear
  gaussprefactor=one/(dossmear*sqrt(two_pi))
  gaussfactor=one/(sqrt2*dossmear)
  write(msg,'(4a,f8.5,2a,f8.5)')ch10,&
&  ' mkphdos : calculating phonon DOS using gaussian method :',ch10,&
&  '  gaussian smearing [meV] = ',dossmear*Ha_meV,ch10,&
&  '  frequency step    [meV] = ',omega_step*Ha_meV
 else if (anaddb_dtset%prtdos==2) then 
  write(msg,'(2a)')ch10,&
&  ' mkphdos : calculating phonon DOS using tetrahedron method '
 end if 
 call wrtout(std_out,msg,'COLL')
!
!=== Initial lower and upper bound of the phonon spectrum ===
 lbound=greatest_real 
 ubound=smallest_real
!
!Save memory during the generation of the q-mesh in the full BZ  
!Take into account the type of Bravais lattice
 facbrv=1
 if (anaddb_dtset%brav==2) facbrv=2
 if (anaddb_dtset%brav==3) facbrv=4

 nmesh=2 ; allocate(ngqpt(3,nmesh))
 do imesh=1,nmesh

  if (imesh==1) then 
!  
!  === Coarse q-mesh used during RF calculation ===
   ngqpt(:,imesh)=anaddb_dtset%ngqpt(1:3)
   nqshft=anaddb_dtset%nqshft 
   allocate(qshft(3,nqshft))
!  this has to be fixed  there is a small inconsistency in the dimension of q1shft
   qshft(:,1:nqshft)=anaddb_dtset%q1shft(:,1:nqshft)
  else 
!  
!  Dense mesh used fo Fourier interpolation. 
   ngqpt(1:3,imesh)=anaddb_dtset%ng2qpt(1:3)
   nqshft=1 !always 1 
   allocate(qshft(3,nqshft))
   qshft(:,1)=anaddb_dtset%q2shft(:)  !small inconsistency in the dimension of q1shft
  end if 

  nqpt_max=(ngqpt(1,imesh)*ngqpt(2,imesh)*ngqpt(3,imesh)*nqshft)/facbrv
  allocate(qibz(3,nqpt_max),qbz(3,nqpt_max))

  qptrlatt(:,:)=0
  qptrlatt(1,1)=ngqpt(1,imesh)
  qptrlatt(2,2)=ngqpt(2,imesh)
  qptrlatt(3,3)=ngqpt(3,imesh)
  option=1 

! here I noticed a problem in the declaration of q1shft in the anaddb datatype 
! FIXME we write on unit 6 just to avoid problem with automatic tests
  call smpbz(anaddb_dtset%brav,6,qptrlatt,nqpt_max,nqbz,nqshft,option,qshft,qbz)
! 
! Reduce the number of such points by symmetrization.
  allocate(ibz2bz(nqbz),wtq(nqbz),wtq_folded(nqbz))
  wtq(:)=one/nqbz         ! Weights sum up to one
  timrev=1 ; option=1     ! TODO timrev should be input 

  call symkpt(gmet,ibz2bz,qbz,nqbz,nqibz,nsym,option,symrec,timrev,wtq,wtq_folded)

  allocate(wtqibz(nqibz))
  do iq=1,nqibz
   wtqibz(iq)=wtq_folded(ibz2bz(iq))
   qibz(:,iq)=qbz(:,ibz2bz(iq))
  end do
  deallocate(wtq_folded,qshft)

  if (anaddb_dtset%prtdos==2.and.imesh==2) then
!  
!  Second mesh with tetrahedron method
!  convert kptrlatt to double and invert, qlatt here refer to the shortest qpt vectors
   rlatt(:,:)=qptrlatt(:,:)
   call matr3inv(rlatt,qlatt)

   allocate(qshft(3,nqshft))
   qshft(:,1)=anaddb_dtset%q2shft(:)  !small inconsistency in the dimension of q1shft
   nkpt_fullbz=nqbz 
   allocate(bz2ibz(nkpt_fullbz))
   allocate(kpt_fullbz(3,nkpt_fullbz))
!  
!  Make full kpoint grid and get equivalence to irred kpoints.
!  This routines scales badly wrt nkpt_fullbz, should introduce checl on the norm.
   call get_full_kgrid(bz2ibz,qlatt,qibz,kpt_fullbz,qptrlatt,nqibz,nkpt_fullbz,nqshft,nsym,qshft,symrel)
!  
!  Get tetrahedra, ie indexes of the full q-points at their summits
!  tetra_full(:,1,i) contains the irred qpt  number
!  tetra_full(:,2,i) contains the full  qpt number
!  tetra_wrap(:,:,i) contains a flag to wrap q-points outside the IBZ (+-1) to get the irreducible tetrahedra
!  the number of equivalent tetrahedra is counted in tetra_mult and the inequivalent few (ntetra < mtetra) are 
!  packed into the beginning of tetra_full
   mtetra=6*nqbz
   allocate(tetra_full(4,2,mtetra),tetra_wrap(3,4,mtetra),tetra_mult(mtetra))
   
   call get_tetra(bz2ibz,gprimd,qlatt,kpt_fullbz,mtetra,nqbz,ntetra_ibz,tetra_full,tetra_mult,tetra_wrap,vtetra)
   write(*,*)' Number of irreducible tetrahedrons = ',ntetra_ibz
   deallocate(bz2ibz)
!  
!  === Arrays Used to store the entire spectrum, Required to calculate tetra weights ===
   allocate(full_phfrq(3*natom,nqibz),full_eigvec(2,3,natom,3*natom,nqibz),stat=istat)
   if (istat/=0) call memerr(sub_name,'full_eigvec',18*natom**2*nqibz,'dp')
  end if 

  do 
!  
!  This infinite loop is used to be sure that the frequency mesh is large enough to contain 
!  the entire phonon spectrum. The mesh is enlarged if, during the Fourier interpolation,
!  a phonon frequency turns out to be outside the interval [omega_min:omega_max]
!  
   if (allocated(omega)) deallocate(omega)
   if (allocated(phdos)) deallocate(phdos)
   if (allocated(pjdos)) deallocate(pjdos)
   out_of_bounds=.FALSE.
   omega_min=lbound ; if (ABS(omega_min)<tol5) omega_min=-tol5
   omega_max=ubound 
   nomega=NINT((omega_max-omega_min)/omega_step)+1
   if (imesh/=1) then 
    write(*,*)'nomega = ',nomega,' omega_min [cm-1] =',omega_min*Ha_cmm1,&
&    '                    omega_max [cm-1] =',omega_max*Ha_cmm1
   end if 
   allocate(omega(nomega))
   do io=1,nomega
    omega(io)=omega_min+omega_step*(io-1)
   end do
   allocate(phdos(nomega),pjdos(nomega,3,natom))
   phdos(:)=zero ; pjdos(:,:,:)=zero
!  
!  === Sum over irreducible q-points ===
   do iq=1,nqibz
    qphon(:)=qibz(:,iq) ; qphnrm=one
!   
!   Get d2cart using interatomic forces and the long-range coulomb interaction through Ewald summation
    call gtdyn9(acell,atmfrc,dielt,anaddb_dtset%dipdip,dyewq0,d2cart,gmet,gprim,&
&    mpert,natom,nrpt,qphnrm,qphon,rcan,rmet,rprim,rpt,trans,ucvol,wghatm,xred,zeff)
!   
!   Get eigenvectors and eigenvalues of the dynamical matrix, eigvec are normalized to one
    call phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,mpert,msym,natom,nsym,ntypat,&
&    phfrq,qphnrm,qphon,rprimd,anaddb_dtset%symdynmat,symrel,typat,ucvol,xred)
    
    dum=MINVAL(phfrq) ; omega_min=MIN(omega_min,dum)
    dum=MAXVAL(phfrq) ; omega_max=MAX(omega_max,dum)
    if (omega_min<lbound .or. omega_max>ubound) out_of_bounds=.TRUE. 

    if (imesh>1.and..not.out_of_bounds) then
     select case (anaddb_dtset%prtdos)
      case (1)
!      
!      === Accumulate PHDOS and PJDOS using gaussian method ===
       do imode=1,3*natom 
        do io=1,nomega
         xx=(omega(io)-phfrq(imode))*gaussfactor
         gaussval=gaussprefactor*exp(-xx*xx)
         phdos(io)=phdos(io) + wtqibz(iq)*gaussval
         do iat=1,natom
          do idir=1,3
           pnorm=eigvec(1,idir,iat,imode)**2+eigvec(2,idir,iat,imode)**2
           pjdos(io,idir,iat)=pjdos(io,idir,iat)+ pnorm*wtqibz(iq)*gaussval
          end do
         end do
        end do 
       end do 
      case (2) 
!      
!      Use tetrahedrons. In this case save phonon frequencies and eigenvectors
!      since the summations is done after the loops over the two meshes.
       full_phfrq(:,iq)=phfrq(:)
       full_eigvec(:,:,:,:,iq)=eigvec(:,:,:,:)
     end select

    end if !Second mesh and not out of boundaries
   end do !Irred q-points

   if (out_of_bounds) then 
    ubound=omega_max+ABS(omega_max/ten)
    lbound=omega_min-ABS(omega_min/ten)
    write(msg,'(6a)')ch10,&
&    ' mkphdos : COMMENT : ',ch10,&
&    ' at least one phonon frequency falls outside the frequency mesh chosen',ch10,&
&    ' restarting the calculation with a larger frequency mesh ' 
    if (imesh>1) call wrtout(std_out,msg,'COLL')
   else
    EXIT !Exit the infinite loop
   end if 
  end do !Infinite loop

  deallocate(ibz2bz,qibz,qbz,wtq,wtqibz)
 end do !imesh
 deallocate(ngqpt)

 allocate(phdos_int(nomega)) 
 phdos_int(:)=zero 

 if (anaddb_dtset%prtdos==2) then 
! 
! Do the integrations with tetrahedrons. All the data are contained in full_phfrq and full_eigvec. 
! lbound and ubound contain the entire spectrum calculated on the dense mesh. 
  allocate(tmp_phfrq(nqibz)) 
  allocate(tweight(nqibz,nomega),dtweightde(nqibz,nomega))
  allocate(pjdos_int(nomega,3,natom)) 
  phdos(:)=zero ; pjdos(:,:,:)=zero ; pjdos_int(:,:,:)=zero
  max_occ=one 

  do imode=1,3*natom 
   tmp_phfrq(:)=full_phfrq(imode,:)
!  
!  === Calculate general integration weights at each irred kpoint as in Blochl et al PRB 49 16223 ===
   call get_tetra_weight(tmp_phfrq,lbound,ubound,max_occ,mtetra,nomega,nqibz,&
&   ntetra_ibz,bzvol,tetra_full,tetra_mult,tweight,dtweightde,vtetra)

   do io=1,nomega
    do iq=1,nqibz
     phdos(io)=phdos(io)+dtweightde(iq,io)
     phdos_int(io)=phdos_int(io)+tweight(iq,io)
     do iat=1,natom
      do idir=1,3
       pnorm=full_eigvec(1,idir,iat,imode,iq)**2 + full_eigvec(2,idir,iat,imode,iq)**2
       pjdos(io,idir,iat)=pjdos(io,idir,iat) + pnorm*dtweightde(iq,io)
       pjdos_int(io,idir,iat)=pjdos_int(io,idir,iat) + pnorm*tweight(iq,io)         
      end do
     end do
    end do
   end do

  end do 
  deallocate(tmp_phfrq)
  deallocate(tweight,dtweightde)
 end if 
!
!=== Convert everything into meV and calculate IPDOS ===
 cfact=Ha_eV*1000 ; unitname='meV'

 omega(:)=omega(:)*cfact
 omega_step=omega_step*cfact
 phdos(:)=phdos(:)/cfact
 pjdos(:,:,:)=pjdos(:,:,:)/cfact

 allocate(pjdos_xyz_typ(nomega,3,ntypat))
 pjdos_xyz_typ(:,:,:)=zero
 allocate(pjdos_typ(nomega,ntypat),pjdos_typ_int(nomega,ntypat))
 pjdos_typ(:,:)=zero ; pjdos_typ_int(:,:)=zero

 do iat=1,natom 
  itype=typat(iat)
  do io=1,nomega
   pjdos_xyz_typ(io,:,itype)=pjdos_xyz_typ(io,:,itype)+pjdos(io,:,iat)
   pjdos_typ(io,itype)=pjdos_typ(io,itype)+sum(pjdos(io,:,iat))
  end do
  if (anaddb_dtset%prtdos==2) then 
   do io=1,nomega
    pjdos_typ_int(io,itype)=pjdos_typ_int(io,itype)+sum(pjdos(io,:,iat))
   end do
  end if 
 end do

 if (anaddb_dtset%prtdos==1) then 
! Evaluate IDOS using simple simpson integration
! TODO should avoid the simpson rule using derf.F90, just to be consistent
  call simpson_int(nomega,omega_step,phdos,phdos_int)
  do itype=1,ntypat
   call simpson_int(nomega,omega_step,pjdos_typ(:,itype),pjdos_typ_int(:,itype))
  end do
 end if 
!
!Open external file and write results

!TODO Here I have to rationalize how to write all this stuff!!

!write(frmt,*)'(es16.8,',ntypat,'(3es16.8))'
!do io=1,nomega
!write(888,frmt)omega(io),((pjdos_xyz_typ(io,ii,itype),ii=1,3),itype=1,ntypat)
!end do
!
!fnam_loc=trim(filname)//'_PHDOS'
 fnam_loc='PHDOS'
 call isfile(fnam_loc,'new')
 unt=get_unit()
 open(unit=unt,file=trim(fnam_loc),form='formatted',status='new')

 write(msg,'(3a)')'# ',ch10,'# Phonon density of states and projected DOS generated by anaddb'
 call wrtout(unt,msg,'COLL')
 if (anaddb_dtset%prtdos==1) then 
  write(msg,'(a,es16.8,2a,i4)')'# Gaussian method with smearing = ',dossmear*cfact,unitname,', nqibz =',nqibz
 else 
  write(unt,'(a,i5,a,i4)')'# Tetrahedron method, number if irreducible tetrahedrons = ',ntetra_ibz,', nqibz= ',nqibz
 end if
 call wrtout(unt,msg,'COLL')
!write(msg,'(4a,i5,5a,es16.8,a,es16.8,3a,es16.8,3a)')&
!& ' The DOS (in modes/',unitname,'/cell) and integrated DOS (in modes/cell) are computed,'&
!& ' at ',nomega,' energies (in ',unitname,') covering the interval ',ch10,&
!& ' between ',omega_min*cfact,' and ',omega_max*cfact,' ',unitname,' by steps of ',omega_step*cfact,' ',unitname,'.'
!call wrtout(unt,msg,'COLL')
 write(msg,'(5a)')'# ',ch10,'# omega     PHDOS    IPHDOS   PJDOS[1]  IPJDOS[1] ...  ',ch10,'# '
 call wrtout(unt,msg,'COLL')
 write(frmt,*)'(',ntypat,'(2es16.8))'
 do io=1,nomega
  write(unt,'(3es16.8)',advance='NO')omega(io),phdos(io),phdos_int(io) 
  do itype=1,ntypat
   write(unt,frmt,advance='NO')pjdos_typ(io,itype),pjdos_typ_int(io,itype)
  end do 
  write(unt,*)
 end do
 close(unt)

 deallocate(omega)
 deallocate(phdos,phdos_int)
 deallocate(pjdos,pjdos_typ,pjdos_typ_int)
 deallocate(pjdos_xyz_typ)

 if (anaddb_dtset%prtdos==2) then
  deallocate(tetra_full,tetra_wrap,tetra_mult)
  deallocate(full_phfrq,full_eigvec,pjdos_int)
 end if

#ifdef DEBUG_MODE
 write(std_out,*)' phdos : exit '
!call flush(std_out)
#endif

end subroutine mkphdos
!!***

