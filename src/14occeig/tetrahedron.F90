!{\src2tex{textfont=tt}}
!!****f* ABINIT/tetrahedron
!! NAME
!! tetrahedron
!!
!! FUNCTION
!! calculate DOS by tetrahedron method
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MVer,XG,SM,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dos_fractions= (to be completed by MVerstra)
!!  dos_fractions_m= m-decompozed (to be completed by MVerstra)
!!  dos_fractions_paw1= contribution to dos fractions from the PAW partial waves (phi)
!!  dos_fractions_pawt1= contribution to dos fractions from the PAW pseudo partial waves (phi_tild)
!!  dtset     structured datatype, from which one uses :
!!   kpt(3,nkpt)  =irreducible kpoints
!!   kptrlatt(3,3)=lattice vectors for full kpoint grid
!!   nband        =number of electronic bands for each kpoint
!!   nkpt         =number of irreducible kpoints
!!   nshiftk      =number of kpoint grid shifts
!!   nsym         =number of symmetries
!!   pawprtdos    =option to output the individual contributions to the partial DOS (0, 1 or 2)
!!   shiftk(3,nshiftk)=kpoint shifts
!!   symrel(3,3,nsym)=symmetry matrices in real space
!!   typat(ntypat)=types of atoms
!!   usepaw       =option for PAW
!!  eigen(mband*nkpt*nsppol)=eigenvalues at irred kpoints
!!  fermie=Fermi energy
!!  fildata=name of the DOS output file
!!  m_dos_flag=option for the m-contributions to the partial DOS
!!  ndosfraction= (to be completed by MVerstra)
!!  paw_dos_flag= option for partial dos in PAW
!!  rprimd(3,3)  =real space unit cell vectors
!!
!! OUTPUT
!!  (no explicit output)
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      dos_hdr_write,get_dos_1band,get_dos_1band_m,get_full_kgrid,get_tetra
!!      get_tetra_weight,int2char4,leave_new,matr3inv,metric,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine tetrahedron (dos_fractions,dos_fractions_m,dos_fractions_paw1,dos_fractions_pawt1,&
& dtset,fermie,eigen,fildata,mbesslang,m_dos_flag,ndosfraction,paw_dos_flag,rprimd)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13recipspace
 use interfaces_14occeig, except_this_one => tetrahedron
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: m_dos_flag,mbesslang,ndosfraction,paw_dos_flag
 real(dp),intent(in) :: fermie
 character(len=fnlen),intent(in) :: fildata
 type(dataset_type),intent(in) :: dtset
!arrays
 real(dp),intent(in) :: dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction)
 real(dp),intent(in) :: dos_fractions_m(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*mbesslang*m_dos_flag)
 real(dp),intent(in) :: dos_fractions_paw1(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*paw_dos_flag)
 real(dp),intent(in) :: dos_fractions_pawt1(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*paw_dos_flag)
 real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol),rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: iat,iband,iene,ifract,ikpt,ioffset_eig,isppol,mband,mtetra,natsph
 integer :: nband,nene,nfiner,nkpt,nkpt_fullbz,nsppol,ntetra,occopt,prtdos
 integer :: unitdos
 real(dp) :: buffer,deltaene,enemax,enemin,enex,eps,integral_DOS,max_occ,rcvol
 real(dp) :: tolfermi,tphysel,tsmear,ucvol,vv
 logical :: bigDOS
 character(len=4) :: tag
 character(len=500) :: frmt,message
 character(len=fnlen) :: bsxfname,tmpfil
!arrays
 integer,allocatable :: indkpt(:),tetra_full(:,:,:),tetra_mult(:)
 integer,allocatable :: tetra_wrap(:,:,:),unitdos_arr(:)
 real(dp) :: gmet(3,3),gprimd(3,3),k1(3),k2(3),k3(3),klatt(3,3),rlatt(3,3)
 real(dp) :: rmet(3,3)
 real(dp),allocatable :: dtweightde(:,:),integ_dos(:,:,:),integ_dos_m(:,:,:)
 real(dp),allocatable :: kpt_fullbz(:,:),partial_dos(:,:,:)
 real(dp),allocatable :: partial_dos_m(:,:,:),tmp_eigen(:),total_dos(:,:)
 real(dp),allocatable :: total_dos_m(:,:),total_dos_paw1(:,:)
 real(dp),allocatable :: total_dos_pawt1(:,:),total_integ_dos(:,:)
 real(dp),allocatable :: total_integ_dos_m(:,:),tweight(:,:)

! *********************************************************************

!DEBUG
!write (6,*) 'tetrahedron : enter'
!ENDDEBUG
!DEBUG
!write (6,*) ' dos_fractions = ', dos_fractions
!ENDDEBUG

!m-decompozed DOS not compatible with PAW-decompozed DOS 
 if(m_dos_flag==1.and.paw_dos_flag==1) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' tetrahedron : ERROR -',ch10,&
&  '  m-decompozed DOS not compatible with PAW-decompozed DOS !'
  call wrtout(6,message,'COLL')
  call leave_new('PERS')
 end if

!Refuse only 1 kpoint: the algorithms are no longer valid. DOH !
 if (dtset%nkpt == 1) then
  write (6, '(a,a)') 'tetrahedron: you need at least 2 kpoints',&
&  ' to use the method'
  return
 end if

!Do not support nshiftk > 1: lattice must be decomposed into boxes
!and this is not always possible (I think) with bizzare shiftks
!normally at this point we have incorporated everything into
!kptrlatt, and only 1 shift is needed (in particular for MP grids).

 if (dtset%nshiftk > 1) then
  write (6, *) 'tetrahedron: problem with a composite k-point grid.'
  write (6, *) '  Only simple lattices are supported. ',&
&  ' Action : use nshiftk=1.'
  write (6, *) ' tetrahedron: skip subroutine.'
  write (6, *) '  dtset%nshiftk, dtset%shiftk = ', dtset%nshiftk, dtset%shiftk
  write (6, *) '  dtset%kptrlatt= ', dtset%kptrlatt
  return
 end if

!Refuse nband different for different kpoints

 do isppol=1,dtset%nsppol
  do ikpt=1,dtset%nkpt
   if ( dtset%nband(dtset%nkpt*(isppol-1) + ikpt) /= dtset%nband(1) ) then
    write (6, *) ' tetrahedron : Error ',&
&    ' nband must be the same for all kpoints'
    write (6, *) '  nband = ', dtset%nband
    write (6, *) ' tetrahedron: skip subroutine.'
    return
   end if
  end do
 end do

!Should be correct formula to calculate nkpt_fullbz
 nkpt_fullbz= dtset%kptrlatt(1,1)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,3) &
& +dtset%kptrlatt(1,2)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,1) &
& +dtset%kptrlatt(1,3)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,2) &
& -dtset%kptrlatt(1,2)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,3) &
& -dtset%kptrlatt(1,3)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,1) &
& -dtset%kptrlatt(1,1)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,2)
 nkpt_fullbz = nkpt_fullbz*dtset%nshiftk

!DEBUG
!write (6,*) ' tetrahedron : nkpt_fullbz = ', nkpt_fullbz
!ENDEBUG

 if (nkpt_fullbz==0) then
  write(6, *) ' tetrahedron : no homogeneous grid ',&
&  'of k-points is defined ...'
  write(6, *)'  in order to obtain the DOS using ',&
&  'the tetrahedron method,'
  write(6, *)'  you need to re-define ngkpt or kptrlatt.'
  write (6, *) ' tetrahedron: skip subroutine.'
  return
 end if

!DEBUG
!write (6,*) 'tetrahedron : passed input testing'
!ENDDEBUG

!Make klatt
 rlatt(:,:) = dtset%kptrlatt(:,:)
 call matr3inv(rlatt,klatt)

!DEBUG
!write (6,*) 'tetrahedron : klatt '
!write (6,'(3(E16.8,1x))') klatt
!ENDDEBUG

!Get metric tensors
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 rcvol = abs (gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
& -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
& +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))

 allocate (indkpt(nkpt_fullbz))
 allocate (kpt_fullbz(3,nkpt_fullbz))

!Make full kpoint grid and get equivalence to irred kpoints
 call get_full_kgrid(indkpt,klatt,dtset%kpt,kpt_fullbz,dtset%kptrlatt,&
& dtset%nkpt,nkpt_fullbz,dtset%nshiftk,dtset%nsym,dtset%shiftk,&
& dtset%symrel)

!Matthieu : this is for the next version. And one should
!define corresponding input variables. Morover,
!is this the proper place ??
!Should be discussed with you.
!call prtxcfermsurf(eigen,fermie,gprimd,klatt,indkpt,&
!&        kpt_fullbz,dtset%mband,dtset%nkpt,&
!&        nkpt_fullbz,dtset%nsppol,dtset%shiftk(:,1))

!DEBUG
!write (6,*) 'tetrahedron : passed get_full_kgrid'
!ENDDEBUG

!Get tetrahedra, ie indexes of the full kpoints at their summits
 mtetra = 6 * nkpt_fullbz
 ntetra = mtetra
!tetra_full(:,1,i) contains the irred kpt  number
!tetra_full(:,2,i) contains the fullbz_kpt number
!tetra_wrap(:,:,i) contains the a flag to wrap kpoints outside the IBZ (+-1)
!to get the irreducible tetrahedra, the number of equivalent
!tetrahedras is counted in tetra_mult and the inequivalent few
!(ntetra < mtetra) are packed into the beginning of tetra_full
 allocate (tetra_full(4,2,mtetra),tetra_wrap(3,4,mtetra))
 allocate(tetra_mult(mtetra))
 call get_tetra (indkpt,gprimd,klatt,kpt_fullbz,mtetra,nkpt_fullbz,ntetra,&
& tetra_full,tetra_mult,tetra_wrap,vv)

!DEBUG
!write (6,*) 'tetrahedron : passed get_tetra ntetra irred = ', ntetra, &
!&            ' over ', mtetra
!ENDDEBUG

 natsph=dtset%natsph

!Open the DOS file
 unitdos=tmp_unit
 if (dtset%prtdos == 2) then
  open (unit=unitdos,file=fildata,status='unknown',form='formatted')
  rewind(unitdos)
 else if (dtset%prtdos == 3) then
  allocate (unitdos_arr(natsph))
  do iat=1,natsph
   unitdos_arr(iat)=tmp_unit+iat-1
   call int2char4(dtset%iatsph(iat),tag)
   tmpfil = trim(fildata)//'_AT'//tag
   open (unit=unitdos_arr(iat),file=tmpfil,&
&   status='unknown',form='formatted')
   rewind(unitdos_arr(iat))
   write (6,*) 'opened file : ', tmpfil, 'unit', unitdos_arr(iat)
  end do
 end if

!Write the header of the DOS file, and determine the energy range and spacing
 prtdos=dtset%prtdos
 buffer=0.01_dp ! Size of the buffer around the min and max ranges
 if (dtset%prtdos == 2) then
  call dos_hdr_write(buffer,deltaene,dtset%dosdeltae,eigen,enemax,enemin,fermie,dtset%mband,&
&  dtset%nband,nene,dtset%nkpt,dtset%nsppol,dtset%occopt,prtdos,&
&  dtset%tphysel,dtset%tsmear,unitdos)
 else if (dtset%prtdos == 3) then
  do iat=1,natsph
   call dos_hdr_write(buffer,deltaene,dtset%dosdeltae,eigen,enemax,enemin,fermie,dtset%mband,&
&   dtset%nband,nene,dtset%nkpt,dtset%nsppol,dtset%occopt,prtdos,&
&   dtset%tphysel,dtset%tsmear,unitdos_arr(iat))
  end do
 end if

!DEBUG
!write (6,*) 'tetrahedron : after dos_hdr_write'
!ENDDEBUG

!DEBUG
!write (6,*) ' nene,ndosfraction,dtset%mband,dtset%nkpt,nene',&
!&              nene,ndosfraction,dtset%mband,dtset%nkpt,nene
!ENDDEBUG

 allocate(partial_dos(nene,ndosfraction,dtset%mband))
 allocate(integ_dos(nene,ndosfraction,dtset%mband))
 allocate(total_dos(nene,ndosfraction))
 allocate(total_integ_dos(nene,ndosfraction))
 allocate(tweight(dtset%nkpt,nene),dtweightde(dtset%nkpt,nene))
 if (paw_dos_flag==1) then
  allocate(total_dos_paw1(nene,ndosfraction))
  allocate(total_dos_pawt1(nene,ndosfraction))
 end if
 if (m_dos_flag==1) then
  allocate(partial_dos_m(nene,ndosfraction*mbesslang,dtset%mband))
  allocate(integ_dos_m(nene,ndosfraction*mbesslang,dtset%mband))
  allocate(total_dos_m(nene,ndosfraction*mbesslang))
  allocate(total_integ_dos_m(nene,ndosfraction*mbesslang))
 end if
 
!Get maximum occupation value (2 or 1)
 if (dtset%nspinor == 1 .and. dtset%nsppol == 1) then
  max_occ = two
 else
  max_occ = one
 end if


!-------------------------------------------------------------------
!For each spin polarisation and band, interpolate band over kpoints
!calculate integration weights and DOS contib from
!-------------------------------------------------------------------
 do isppol = 1, dtset%nsppol

  total_dos (:,:) = zero
  total_integ_dos (:,:) = zero
  if (paw_dos_flag==1) then
   total_dos_paw1(:,:)=zero;total_dos_pawt1(:,:)=zero
  end if

  if (dtset%nsppol==2) then
   if(isppol==1) write(message,'(a,16x,a)')  '#','Spin-up DOS'
   if(isppol==2) write(message,'(a,16x,a)')  '#','Spin-dn DOS'
   if (dtset%prtdos == 2) then
    call wrtout(unitdos,message,'COLL')
   else if (dtset%prtdos == 3) then
    do iat=1,natsph
     call wrtout(unitdos_arr(iat),message,'COLL')
    end do
   end if
  end if

  allocate(tmp_eigen(dtset%nkpt))
  do iband = 1, dtset%nband(1)

!  DEBUG
!  write (6,*) ' tetrahedron : isppol, iband = ', isppol, iband
!  write (6,*) 'tetra_full(:,:,2) ', tetra_full(:,:,2)
!  ENDDEBUG

   ioffset_eig = dtset%mband*dtset%nkpt*(isppol-1)

   tmp_eigen(:) = zero
!  For each band get its contribution
   do ikpt=1,dtset%nkpt
    tmp_eigen(ikpt) = eigen(ioffset_eig + dtset%mband*(ikpt-1) + iband)
   end do
!  
!  calculate general integration weights at each irred kpoint as in Blochl et al
!  PRB 49 16223
!  
   call get_tetra_weight(tmp_eigen,enemin,enemax,&
&   max_occ,mtetra,nene,dtset%nkpt,ntetra,rcvol,&
&   tetra_full,tetra_mult,tweight,dtweightde,vv)

!  
!  determine fermi surface contribution from this band,
!  interpolated in irred tetrahedra, on a finer subgrid of klatt
!  NOT YET FUNCTIONAL
!  if(dtset%prtfsurf==1)then
!  tolfermi = dtset%userrc
!  nfiner = dtset%userrb
!  call get_fsurf_1band(dtset,tmp_eigen,fermie,klatt,kpt_fullbz,&
!  &  mtetra,nfiner,nkpt_fullbz,ntetra,tetra_full,tetra_mult,tetra_wrap,tolfermi)
!  end if
!  
!  calculate DOS and integrated DOS projected with the input dos_fractions
!  
   if (paw_dos_flag==1) then
    call get_dos_1band (dos_fractions_paw1(:,iband,isppol,:),enemin,enemax,&
&    integ_dos(:,:,iband),nene,dtset%nkpt,ndosfraction,&
&    partial_dos(:,:,iband),tweight,dtweightde)
    do ifract=1,ndosfraction
     do iene=1,nene
      total_dos_paw1(iene,ifract)=total_dos_paw1(iene,ifract)+partial_dos(iene,ifract,iband)
     end do
    end do
    call get_dos_1band (dos_fractions_pawt1(:,iband,isppol,:),enemin,enemax,&
&    integ_dos(:,:,iband),nene,dtset%nkpt,ndosfraction,&
&    partial_dos(:,:,iband),tweight,dtweightde)
    do ifract=1,ndosfraction
     do iene=1,nene
      total_dos_pawt1(iene,ifract)=total_dos_pawt1(iene,ifract)+partial_dos(iene,ifract,iband)
     end do
    end do
   end if

   call get_dos_1band (dos_fractions(:,iband,isppol,:),enemin,enemax,&
&   integ_dos(:,:,iband),nene,dtset%nkpt,ndosfraction,&
&   partial_dos(:,:,iband),tweight,dtweightde)

   if (m_dos_flag==1) then
    call get_dos_1band_m (dos_fractions_m(:,iband,isppol,:),enemin,enemax,&
&    integ_dos_m(:,:,iband),nene,dtset%nkpt,ndosfraction*mbesslang,&
&    partial_dos_m(:,:,iband),tweight,dtweightde)
   end if

!  DEBUG
!  write (6,*) 'tetrahedron : passed get_dos_1band for iband = ', iband
!  ENDDEBUG

!  Add to total dos
   do ifract=1,ndosfraction
    do iene=1,nene
     total_dos(iene,ifract) = total_dos(iene,ifract) + &
&     partial_dos(iene,ifract,iband)
     total_integ_dos(iene,ifract) = total_integ_dos(iene,ifract) + &
&     integ_dos(iene,ifract,iband)
    end do
   end do
   if (m_dos_flag==1) then
    do ifract=1,ndosfraction*mbesslang
     do iene=1,nene
      total_dos_m(iene,ifract) = total_dos_m(iene,ifract) + &
&      partial_dos_m(iene,ifract,iband)
      total_integ_dos_m(iene,ifract) = total_integ_dos_m(iene,ifract) + &
&      integ_dos_m(iene,ifract,iband)
     end do
    end do
   end if
   
  end do ! iband

  deallocate(tmp_eigen)

  write (6,*) 'about to write to the DOS file '
! Write the DOS value in the DOS file
  enex=enemin
  integral_DOS=zero
  if(prtdos==2)then
   write(message, '(a)' )'#index  energy(Ha)     DOS  integrated DOS'
   call wrtout(unitdos,message,'COLL')
  else if(prtdos==3)then

   bigDOS=(maxval(total_dos)>999._dp)
   if (paw_dos_flag/=1) then
    do iat=1,natsph
     write(message, '(3a,i5,a,i5,a,a,es16.6,3a)' ) &
&     '# Local DOS (columns 3-7) and integrated local DOS (columns 8-12),',ch10,&
&     '# for atom number iat=',iat,'  iatom=',dtset%iatsph(iat),ch10,&
&     '# inside sphere of radius ratsph=',dtset%ratsph(dtset%typat(dtset%iatsph(iat))),' Bohr.',ch10,"#"
     call wrtout(unitdos_arr(iat),message,'COLL')
     if (dtset%usepaw==1.and.dtset%pawprtdos==2) then
      write(message, '(3a)' ) &
&      '# PAW: note that only all-electron on-site part has been used to compute DOS !',ch10,"#"
      call wrtout(unitdos_arr(iat),message,'COLL')
     end if
     if (bigDOS) then
      write(message, '(a,a)' ) &
&      '#index  energy(Ha)   l=0       l=1       l=2       l=3       l=4',&
&      '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
     else
      write(message, '(a,a)' ) &
&      '#index  energy(Ha)  l=0      l=1      l=2      l=3      l=4',&
&      '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
     end if
     if (m_dos_flag==1) then
      write(message, '(7a)' ) trim(message),'          ',&
&      '  lm=0 0',&
&      '  lm=1-1  lm=1 0  lm=1 1',&
&      '  lm=2-2  lm=2-1  lm=2 0  lm=2 1  lm=2 2',&
&      '  lm=3-3  lm=3-2  lm=3-1  lm=3 0  lm=3 1  lm=3 2  lm=3 3',&
&      '  lm=4-4  lm=4-3  lm=4-2  lm=4-1  lm=4 0  lm=4 1  lm=4 2  lm=4 3  lm=4 4'
     end if
     call wrtout(unitdos_arr(iat),message,'COLL')
    end do
   else
    do iat=1,natsph
     write(message, '(9a,i5,a,i5,a,a,es16.6,3a)' ) &
&     '# Local DOS (columns 3-7),',ch10,&
&     '#  plane-waves contrib. to DOS (columns 8-12),',ch10,&
&     '#  AE on-site  contrib. to DOS (columns 13-17),',ch10,&
&     '# -PS on-site  contrib. to DOS (columns 18-22),',ch10,&
&     '# for atom number iat=',iat,'  iatom=',dtset%iatsph(iat),ch10,&
&     '# inside sphere of radius ratsph=',dtset%ratsph(dtset%typat(dtset%iatsph(iat))),' Bohr.',ch10,"#"
     call wrtout(unitdos_arr(iat),message,'COLL')
     if (bigDOS) then
      write(message, '(4a)' ) &
&      '#index  energy(Ha)   l=0       l=1       l=2       l=3       l=4',&
&      '       (PW)  l=0       l=1       l=2       l=3       l=4',&
&      '      (Phi)  l=0       l=1       l=2       l=3       l=4',&
&      '     (tPhi)  l=0       l=1       l=2       l=3       l=4'
     else
      write(message, '(4a)' ) &
&      '#index  energy(Ha)  l=0      l=1      l=2      l=3      l=4',&
&      '       (PW) l=0      l=1      l=2      l=3      l=4',&
&      '      (Phi) l=0      l=1      l=2      l=3      l=4',&
&      '     (tPhi) l=0      l=1      l=2      l=3      l=4'
     end if
     call wrtout(unitdos_arr(iat),message,'COLL')
    end do
   end if
  end if

! write(message, *) '#index energy DOS for ',ndosfraction,' channels',&
! &            '  integrated DOS '
! call wrtout(unitdos,message,'COLL')

  if(prtdos==2)then
   do iene=1,nene
!   Print the data for this energy. Note the upper limit, to be
!   consistent with the format. The use of "E" format is not adequate,
!   for portability of the self-testing procedure.
!   write(message, '(i5,f9.4,f14.6)' ) iene-1,enex,total_dos(iene,:)
    write(message, '(i6,f11.5,5f10.4,10x,5f8.2)') iene-1, enex, &
&    min(total_dos(iene,:),9999.9999_dp),&
&    total_integ_dos(iene,:)
    call wrtout(unitdos,message,'COLL')
    enex=enex+deltaene
   end do
  else if(prtdos==3)then
   if (paw_dos_flag/=1) then
    if (bigDOS) then
     frmt='(i6,f11.5,5f10.4,10x,5f8.2,10x,25f8.2)'
    else
     frmt='(i6,f11.5,5f9.4 ,10x,5f8.2,10x,25f8.2)'
    end if
    do iene=1,nene
     do iat=1,natsph
!     DEBUG
!     write(message, '(i6,5f9.4,10x,5f7.2))') iene-1, &
!     &     total_dos(iene,iat*5-2)
!     ENDDEBUG
!     Note the upper limit, to be
!     consistent with the format. The use of "E" format is not adequate,
!     for portability of the self-testing procedure.
      if (m_dos_flag==0) then
       write(message,fmt=frmt) iene-1, enex, &
&       min(total_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang),9999.9999_dp), &
&       total_integ_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang)
      else
       write(message,fmt=frmt) iene-1, enex, &
&       min(total_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang),9999.9999_dp),&
&       total_integ_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang),&
&       min(total_dos_m(iene,(iat-1)*mbesslang**2+1:iat*mbesslang**2),9999.9999_dp)
      end if

      call wrtout(unitdos_arr(iat),message,'COLL')
     end do
     enex=enex+deltaene
    end do
   else
    if (bigDOS) then
     frmt='(i6,f11.5,5f10.4,3(6x,5f10.4))'
    else
     frmt='(i6,f11.5,5f9.4 ,3(6x,5f9.4 ))'
    end if
    do iene=1,nene
     do iat=1,natsph
      write(message,fmt=frmt) iene-1, enex, &
&      min(total_dos(iene,iat*5-4:iat*5),9999.9999_dp),&
&      min(total_dos(iene,iat*5-4:iat*5)-total_dos_paw1(iene,iat*5-4:iat*5)&
&      +total_dos_pawt1(iene,iat*5-4:iat*5),9999.9999_dp),&
&      min(total_dos_paw1(iene,iat*5-4:iat*5),9999.9999_dp),&
&      min(total_dos_pawt1(iene,iat*5-4:iat*5),9999.9999_dp)
      call wrtout(unitdos_arr(iat),message,'COLL')
     end do
     enex=enex+deltaene
    end do
   end if
  end if

! integral_DOS=integral_DOS+deltaene*sum(total_dos(iene,:))
  integral_DOS=sum(total_integ_dos(nene,:))
  write(message, '(a,es16.8)' ) ' tetrahedron : integrate to',integral_DOS
  call wrtout(6,message,'COLL')

 end do ! isppol

 close(unitdos)
 if(prtdos==3)deallocate (unitdos_arr)
 deallocate (tetra_full,tetra_wrap,tetra_mult)
 deallocate (indkpt)
 deallocate (kpt_fullbz)
 deallocate (partial_dos)
 deallocate (integ_dos)
 deallocate (total_dos)
 deallocate (total_integ_dos)
 deallocate (tweight,dtweightde)
 if (m_dos_flag==1) deallocate(partial_dos_m,integ_dos_m,total_dos_m,total_integ_dos_m)
 if (paw_dos_flag==1) deallocate(total_dos_paw1,total_dos_pawt1)

!DEBUG
!write(6,*)' tetrahedron : exit '
!ENDDEBUG

end subroutine tetrahedron
!!***
