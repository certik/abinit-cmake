!{\src2tex{textfont=tt}}
!!****f* ABINIT/nonlinear
!! NAME
!! nonlinear
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations of
!! non linear response functions.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (MVeithen, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  codvsn = code version
!!  dtfil <type(datafiles_type)> = variables related to files
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  etotal = new total energy (no meaning at output)
!!  iexit= exit flag
!!  mband = maximum number of bands
!!  mgfft = maximum single fft dimension
!!  mkmem = maximum number of k points which can fit in core memory
!!  mpi_enreg=informations about MPI pnarallelization
!!  mpw   = maximum number of planewaves in basis sphere (large number)
!!  natom = number of atoms in unit cell
!!  nfft  = (effective) number of FFT grid points (for this processor)
!!  nkpt  = number of k points
!!  nspden = number of spin-density components
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  nsym   = number of symmetry elements in space group
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  xred(3,natom) = reduced atomic coordinates
!!
!! OUTPUT
!!
!!  npwtot(nkpt) = total number of plane waves at each k point
!!
!! SIDE EFFECTS
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!
!! TODO
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      bstruct_clean,bstruct_init,d3output,d3sym,distrb2,fourdp,getcut
!!      getkgrid,getshell,hdr_clean,hdr_init,hdr_update,initmv,inwffil,ioarr
!!      ioddb8,kpgio,leave_new,loop3dte,mkcore,nlopt,psddb8,pspini,rhohxc
!!      setsym,setup1,status,symzat,sytens,timab,wffclose,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine nonlinear(codvsn,dtfil,dtset,etotal,iexit,&
&  mband,mgfft,mkmem,mpi_enreg,mpw,natom,nfft,nkpt,npwtot,nspden,&
&  nspinor,nsppol,nsym,occ,pawrad,pawtab,psps,xred)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_13psp
 use interfaces_13recipspace
 use interfaces_13xc
 use interfaces_14iowfdenpot
 use interfaces_15common
 use interfaces_16response
 use interfaces_18seqpar
 use interfaces_21drive, except_this_one => nonlinear
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iexit,mband,mgfft,mkmem,mpw,nfft
 integer,intent(inout) :: natom,nkpt,nspden,nspinor,nsppol,nsym
 real(dp),intent(inout) :: etotal
 character(len=6),intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 integer,intent(out) :: npwtot(nkpt)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol),xred(3,natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat,psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat,psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: fform=2,level=20,response=1
 integer :: accessfil,ask_accurate,bantot,choice,dum_nshiftk,fformr=52,flag
 integer :: formeig,fullinit,gencond,gscase,i1dir,i1pert,i2dir,i2pert,i3dir
 integer :: i3pert,ierr,ii,ireadwf,jj,kk,mkmem_max,mpert,n1,n2,n3,n3xccc,nblok
 integer :: nkpt3,nkxc,nneigh,option,optorth,rdwr,rdwrpaw,vrsddb
 real(dp) :: boxcut,ecore,ecut_eff,enxc,fermie,gsqcut,gsqcut_eff,gsqcutdg_eff
 real(dp) :: rdum,residm,tolwfr,ucvol,vxcavg
 character(len=500) :: message
 character(len=fnlen) :: ddbnm,dscrpt
 type(bandstructure_type) :: bstruct
 type(dens_sym_operator_type) :: densymop_gs
 type(hdr_type) :: hdr
 type(wffile_type) :: wffgs,wfftgs
 type(wvl_data) :: wvl
!arrays
 integer :: dum_dsifkpt(3),dum_kptrlatt(3,3),dum_vacuum(3),perm(6)
 integer,allocatable :: blkflg(:,:,:,:,:,:),carflg(:,:,:,:,:,:),cgindex(:,:)
 integer,allocatable :: indsym(:,:,:),irrzon(:,:,:),kg(:,:),kneigh(:,:)
 integer,allocatable :: kptindex(:,:),npwarr(:),pwind(:,:,:),rf1pert(:)
 integer,allocatable :: rf2pert(:),rf3pert(:),rfpert(:,:,:,:,:,:),symrec(:,:,:)
 real(dp) :: dum_shiftk(3,dtset%nshiftk),dummy2(6),gmet(3,3),gprimd(3,3),k0(3)
 real(dp) :: rmet(3,3),rprimd(3,3),strsxc(6),tsec(2)
 real(dp),allocatable :: amass(:),cg(:,:),d3cart(:,:,:,:,:,:,:)
 real(dp),allocatable :: d3lo(:,:,:,:,:,:,:),doccde(:),dum_kptns(:,:)
 real(dp),allocatable :: dum_wtk(:),dyfrx2(:,:,:),eigen(:),grxc(:,:),k3xc(:)
 real(dp),allocatable :: kpt3(:,:),kxc(:,:),mvwtk(:,:),phnons(:,:,:),rhog(:,:)
 real(dp),allocatable :: rhor(:,:),vhartr(:),vxc(:,:),work(:),xccc3d(:)
 character(len=fnlen) :: tmpfil(15)
 type(pawrhoij_type),allocatable :: pawrhoij(:)

! ***********************************************************************

!DEBUG
!write(6,*)'nonlinear.f : enter'
!stop
!ENDDEBUG

 call timab(501,1,tsec)
 call status(0,dtfil%filstat,iexit,level,'enter         ')

!
!If dtset%accesswff == 2 set all array outputs to netcdf format
!
 accessfil = 0
 if (dtset%accesswff == 2) then
  accessfil = 1
 end if
 if (dtset%accesswff == 3) then
  accessfil = 3
 end if

 mpi_enreg%paralbd=0
 mpi_enreg%parareel=0
 mpi_enreg%me_fft=0
 mpi_enreg%nproc_fft=1
 mpi_enreg%paral_fft=0
 mpi_enreg%paral_level=2

 if (mpi_enreg%paral_compil_kpt == 1) then
  allocate(mpi_enreg%proc_distrb(nkpt,mband,nsppol))
  call distrb2(mband, dtset%nband, nkpt, nsppol, mpi_enreg)
 end if


!Check if the perturbations asked in the input file
!can be computed

 if (((dtset%rf1phon == 1).and.(dtset%rf2phon == 1)).or. &
& ((dtset%rf1phon == 1).and.(dtset%rf3phon == 1)).or. &
& ((dtset%rf2phon == 1).and.(dtset%rf3phon == 1))) then

  write(message,'(10(a))') ch10,&
&  ' nonlinear : ERROR - ',ch10,&
&  '  You have asked for a third-order derivative with respect to',ch10,&
&  '  2 or more atomic displacements.',ch10,&
&  '  This is not allowed yet.',ch10,&
&  '  Action : change rf1phon, rf2phon or rf3phon in your input file.'
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')

 end if

!Define the set of admitted perturbations taking into account
!the possible permutations

 mpert=natom+6
 allocate(blkflg(3,mpert,3,mpert,3,mpert))
 allocate(carflg(3,mpert,3,mpert,3,mpert))
 allocate(rfpert(3,mpert,3,mpert,3,mpert))
 allocate(rf1pert(mpert),rf2pert(mpert),rf3pert(mpert))
 allocate(d3lo(2,3,mpert,3,mpert,3,mpert))
 allocate(d3cart(2,3,mpert,3,mpert,3,mpert))
 blkflg(:,:,:,:,:,:) = 0
 d3lo(:,:,:,:,:,:,:) = 0_dp
 rfpert(:,:,:,:,:,:) = 0
 rf1pert(:) = 0 ; rf2pert(:) = 0 ; rf3pert(:) = 0

 if (dtset%rf1phon==1) rf1pert(dtset%rf1atpol(1):dtset%rf1atpol(2))=1
 if (dtset%rf2phon==1) rf2pert(dtset%rf2atpol(1):dtset%rf2atpol(2))=1
 if (dtset%rf3phon==1) rf3pert(dtset%rf3atpol(1):dtset%rf3atpol(2))=1
 if (dtset%rf1elfd/=0) rf1pert(natom+2)=1
 if (dtset%rf2elfd/=0) rf2pert(natom+2)=1
 if (dtset%rf3elfd/=0) rf3pert(natom+2)=1

 do i1pert = 1, mpert
  do i1dir = 1, 3
   do i2pert = 1, mpert
    do i2dir = 1, 3
     do i3pert = 1, mpert
      do i3dir = 1, 3
       perm(1) = rf1pert(i1pert)*dtset%rf1dir(i1dir)* &
&       rf2pert(i2pert)*dtset%rf2dir(i2dir)*rf3pert(i3pert)*dtset%rf3dir(i3dir)
       perm(2) = rf1pert(i1pert)*dtset%rf1dir(i1dir)* &
&       rf2pert(i3pert)*dtset%rf2dir(i3dir)*rf3pert(i2pert)*dtset%rf3dir(i2dir)
       perm(3) = rf1pert(i2pert)*dtset%rf1dir(i2dir)* &
&       rf2pert(i1pert)*dtset%rf2dir(i1dir)*rf3pert(i3pert)*dtset%rf3dir(i3dir)
       perm(4) = rf1pert(i2pert)*dtset%rf1dir(i2dir)* &
&       rf2pert(i3pert)*dtset%rf2dir(i3dir)*rf3pert(i1pert)*dtset%rf3dir(i1dir)
       perm(5) = rf1pert(i3pert)*dtset%rf1dir(i3dir)* &
&       rf2pert(i2pert)*dtset%rf2dir(i2dir)*rf3pert(i1pert)*dtset%rf3dir(i1dir)
       perm(6) = rf1pert(i3pert)*dtset%rf1dir(i3dir)* &
&       rf2pert(i1pert)*dtset%rf2dir(i1dir)*rf3pert(i2pert)*dtset%rf3dir(i2dir)
       if (sum(perm(:)) > 0) rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 1
      end do
     end do
    end do
   end do
  end do
 end do

!DEBUG
!do i1pert = 1, natom + 2
!do i1dir = 1, 3
!do i2pert = 1, natom + 2
!do i2dir = 1, 3
!do i3pert = 1, natom + 2
!do i3dir = 1,3
!if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)/=0) then
!write(100,'(6(2x,i3),5x,i3)')i1pert,i1dir,i2pert,i2dir,i3pert,i3dir,&
!&    rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
!end if
!end do
!end do
!end do
!end do
!end do
!end do
!stop
!ENDDEBUG

!Determine the symmetrical perturbations

 allocate(irrzon(dtset%nfft**(1-1/nsym),2,nspden/nsppol))
 allocate(phnons(2,dtset%nfft**(1-1/nsym),nspden/nsppol))
 allocate(indsym(4,nsym,natom),symrec(3,3,nsym))
 call status(0,dtfil%filstat,iexit,level,'call setsym   ')
 call setsym(densymop_gs,indsym,irrzon,dtset%iscf,natom,&
& nfft,dtset%ngfft,nspden,nsppol,nsym,&
& phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)

 call status(0,dtfil%filstat,iexit,level,'call symzat   ')
 call symzat(indsym,natom,nsym,dtset%symrel,dtset%tnons,xred)

 call status(0,dtfil%filstat,iexit,level,'call sytens   ')
 call sytens(indsym,mpert,natom,nsym,rfpert,symrec,dtset%symrel)

 write(message, '(a,a,a,a,a)' ) ch10, &
& ' The list of irreducible elements of the Raman and non-linear',&
& ch10,' optical susceptibility tensors is:',ch10
 call wrtout(ab_out,message,'COLL')

 write(message,'(12x,a)')&
& 'i1pert  i1dir   i2pert  i2dir   i3pert  i3dir'
 call wrtout(ab_out,message,'COLL')

 n1 = 0
 do i1pert = 1, natom + 2
  do i1dir = 1, 3
   do i2pert = 1, natom + 2
    do i2dir = 1, 3
     do i3pert = 1, natom + 2
      do i3dir = 1,3
       if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
        n1 = n1 + 1
        write(message,'(2x,i4,a,6(5x,i3))') n1,')', &
&        i1pert,i1dir,i2pert,i2dir,i3pert,i3dir
        call wrtout(ab_out,message,'COLL')
       else if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==-2) then
        blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 1
       end if
      end do
     end do
    end do
   end do
  end do
 end do
 write(message,'(a,a)') ch10,ch10
 call wrtout(ab_out,message,'COLL')

!Create names for the temporary files based on dtfil%filnam_ds(5)
!by appending adequate string (see respfn.f).
 tmpfil(1) =trim(dtfil%filnam_ds(5))//'_1WF1'
 tmpfil(2) =trim(dtfil%filnam_ds(5))//'_1WF2'
 tmpfil(3) =trim(dtfil%filnam_ds(5))//'_KG'
 tmpfil(4) =trim(dtfil%filnam_ds(5))//'_KGQ'
 tmpfil(5) =trim(dtfil%filnam_ds(5))//'_KG1'
 tmpfil(6) =trim(dtfil%filnam_ds(5))//'_DUM'
 tmpfil(7) =trim(dtfil%filnam_ds(5))//'_WFGS'
 tmpfil(8) =trim(dtfil%filnam_ds(5))//'_WFKQ'
 tmpfil(10)=trim(dtfil%filnam_ds(5))//'_YLM'
 tmpfil(11)=trim(dtfil%filnam_ds(5))//'_YLM1'
 tmpfil(12)=trim(dtfil%filnam_ds(5))//'_PAW'
 tmpfil(13)=trim(dtfil%filnam_ds(5))//'_PAW1'
 tmpfil(14)=trim(dtfil%filnam_ds(5))//'_PAWQ'
 tmpfil(15)=trim(dtfil%filnam_ds(5))//'_HD1'

!Set up for iterations
 ecut_eff= (dtset%ecut) * (dtset%dilatmx) **2
 allocate(amass(natom))
 call status(0,dtfil%filstat,iexit,level,'call setup1   ')
 call setup1(dtset%acell_orig,amass,dtset%amu,bantot,&
& ecut_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcut_eff,dtset%iboxcut,dtset%intxc,dtset%ionmov,&
& natom,dtset%nband,dtset%ngfft,dtset%ngfft,nkpt,dtset%nqpt,nsppol,nsym,psps%ntypat,&
& dtset%qptn,response,rmet,dtset%rprim_orig,rprimd,dtset%typat,ucvol,psps%usepaw)

!Set up the basis sphere of planewaves
 allocate(kg(3,mpw*dtset%mk1mem),npwarr(nkpt))
 call status(0,dtfil%filstat,iexit,level,'call kpgio    ')
 call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg,tmpfil(3),&
& dtset%kptns,mkmem,dtset%nband,nkpt,'PERS',mpi_enreg,&
& mpw,npwarr,npwtot,nsppol,dtfil%unkg)

!Recompute first large sphere cut-off gsqcut,
!without taking into account dilatmx
 k0(:)=0.0_dp
 call status(0,dtfil%filstat,iexit,level,'call getcut   ')
 call getcut(boxcut,dtset%ecut,gmet,gsqcut,dtset%iboxcut,6,k0,dtset%ngfft)

!Open and read pseudopotential files
 ecore = 0_dp
 call status(0,dtfil%filstat,iexit,level,'call pspini   ')
 call pspini(dtset,ecore,gencond,gsqcut_eff,gsqcutdg_eff,level,&
& pawrad,pawtab,psps,rprimd)

!Initialize band structure datatype
 allocate(doccde(bantot),eigen(bantot))
 doccde(:)=zero ; eigen(:)=zero
 call bstruct_init(bantot,bstruct,doccde,eigen,dtset%istwfk,dtset%kptns,&
& dtset%nband,nkpt,npwarr,nsppol,dtset%occ_orig,dtset%wtk)
 deallocate(doccde,eigen)

!Initialize header
 gscase=0
 call hdr_init(bstruct,codvsn,dtset,hdr,pawtab,gscase,psps)

!Update header, with evolving variables, when available
!Here, rprimd, xred and occ are available
 residm=hdr%residm ; fermie=hdr%fermie
 call hdr_update(bantot,etotal,fermie,hdr,natom,&
& residm,rprimd,occ,pawrhoij,psps%usepaw,xred)

!Clean band structure datatype (should use it more in the future !)
 call bstruct_clean(bstruct)

!Read ground-state wavefunctions
 allocate(cg(2,dtset%mpw*dtset%nspinor*mband*dtset%mkmem*dtset%nsppol))
 allocate(eigen(mband*dtset%nkpt*dtset%nsppol))
 optorth=1;if (psps%usepaw==1) optorth=0
 ireadwf=1 ; formeig=0
 eigen(:)=0_dp ; ask_accurate=1
 call status(0,dtfil%filstat,iexit,level,'call inwffil  ')
 call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen,dtset%exchn2n3d,&
& formeig,gmet,hdr,ireadwf,dtset%istwfk,kg,dtset%kptns,&
& dtset%localrdwf,mband,dtset%mkmem,mpi_enreg,mpw,&
& dtset%nband,dtset%ngfft,dtset%nkpt,&
& npwarr,dtset%nspden,dtset%nspinor,dtset%nsppol,dtset%nsym,&
& occ,optorth,psps,dtset%prtvol,rprimd,dtset%symafm,dtset%symrel,dtset%tnons,&
& dtfil%unkg,wffgs,wfftgs,dtfil%unwffgs,dtfil%unwftgs,&
& dtfil%fnamewffk,tmpfil(7),wvl)

 if (ireadwf==1) then
  call WffClose(wffgs,ierr)
 end if

 deallocate(eigen)


 allocate(rhog(2,nfft),rhor(nfft,nspden))
!Get the ground state chagre density

 if (dtset%getden /= 0) then

  rdwr = 1 ; rdwrpaw = 0
! set to 1 for netcdf
  call status(0,dtfil%filstat,iexit,level,'call ioarr    ')
  call ioarr(accessfil,rhor, dtset, etotal,fformr,dtfil%fildensin,hdr, mpi_enreg, &
&  nfft,pawrhoij,rdwr,rdwrpaw,dtset%ngfft)
! Compute up+down rho(G) by fft
  allocate(work(nfft))
  work(:)=rhor(:,1)
  call status(0,dtfil%filstat,iexit,level,'call fourdp   ')
  call fourdp(1,rhog,work,-1,mpi_enreg,nfft,dtset%ngfft,dtset%paral_kgb,0)
  deallocate(work)

 end if



!Compute core electron density xccc3d
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 allocate(grxc(3,natom),vxc(nfft,nspden),vhartr(nfft))
 n3xccc=0
 if (psps%n1xccc/=0) n3xccc=nfft
 allocate(xccc3d(n3xccc))
 if (psps%n1xccc/=0) then
  option=1
  allocate(dyfrx2(3,3,natom))
  call status(0,dtfil%filstat,iexit,level,'call mkcore   ')
  call mkcore(dummy2,dyfrx2,grxc,mpi_enreg,natom,nfft,nspden,psps%ntypat,&
&  n1,psps%n1xccc,n2,n3,option,rprimd,dtset%typat,ucvol,vxc,&
&  psps%xcccrc,psps%xccc1d,xccc3d,xred)
  deallocate(dyfrx2)
 end if


!Comput kxc (second- and third-order exchange-correlation kernel)
 option=3
 nkxc=2*nspden-1
 if(dtset%xclevel==2) nkxc=23
 allocate(kxc(nfft,nkxc),k3xc(nfft))

 call status(0,dtfil%filstat,iexit,level,'call rhohxc   ')
 call rhohxc(dtset,enxc,gsqcut,psps%usepaw,kxc,mpi_enreg,nfft,dtset%ngfft,&
& work,0,work,0,nkxc,nspden,n3xccc,option,rhog,rhor,rprimd,strsxc,1,&
& vhartr,vxc,vxcavg,xccc3d,k3xc)


 deallocate(vhartr,vxc,xccc3d)

!Initializa finite difference calculation of the ddk

 call status(0,dtfil%filstat,iexit,level,'call getshell ')


 nkpt3 = 0

!Prepare first call to getkgrid (obtain number of k points in FBZ)
 dum_dsifkpt(:) = 1
 dum_kptrlatt(:,:) = dtset%kptrlatt(:,:)
 dum_nshiftk = dtset%nshiftk
 dum_shiftk(:,:) = zero
 dum_shiftk(:,1:dtset%nshiftk) = dtset%shiftk(:,1:dtset%nshiftk)
 dum_vacuum(:) = 0

 allocate(dum_kptns(3,0),dum_wtk(0))
 call getkgrid(dum_dsifkpt,ab_out,dtset%iscf,dum_kptns,3,dum_kptrlatt,&
& rdum,dtset%nsym,0,nkpt3,dum_nshiftk,dtset%nsym,&
& rprimd,dum_shiftk,dtset%symafm,dtset%symrel,dtset%tnons,&
& dum_vacuum,dum_wtk)
 deallocate(dum_kptns,dum_wtk)

 write (*,*) 'nonlinear : nkpt, nkpt3 = ',nkpt,nkpt3

 allocate(kneigh(30,nkpt),kptindex(2,nkpt3),mvwtk(30,nkpt))
 allocate(kpt3(3,nkpt3))
 call getshell(gmet,kneigh,kptindex,dtset%kptopt,&
& dtset%kptrlatt,dtset%kptns,kpt3,mkmem,mkmem_max,mpi_enreg,mvwtk,&
& nkpt,nkpt3,nneigh,dtset%nshiftk,rmet,rprimd,dtset%shiftk,dtset%wtk)

 allocate(pwind(mpw,nneigh,mkmem),cgindex(nkpt,nsppol))
 if (mpi_enreg%paral_compil_kpt == 1) then
  allocate(mpi_enreg%kptdstrb(mpi_enreg%nproc,nneigh,mkmem_max))
 end if
 call status(0,dtfil%filstat,iexit,level,'call initmv   ')
 call initmv(cgindex,dtfil,dtset,gmet,kg,kneigh,kptindex,&
& kpt3,mband,mkmem,mkmem_max,mpi_enreg,mpw,dtset%nband,nkpt,&
& nkpt3,nneigh,npwarr,nsppol,occ,pwind)

 call status(0,dtfil%filstat,iexit,level,'call loop3dte ')
 call loop3dte(blkflg,cg,cgindex,dtfil,dtset,d3lo,etotal,gmet,gprimd,gsqcut,&
& gsqcut_eff,hdr,kg,kneigh,kptindex,kpt3,kxc,k3xc,mband,mgfft,&
& mkmem,mkmem_max,dtset%mk1mem,mpert,mpi_enreg,mpw,mvwtk,natom,nfft,&
& nkpt,nkpt3,nkxc,nneigh,nspinor,nsppol,npwarr,occ,psps,pwind,&
& rfpert,rmet,rprimd,tmpfil,ucvol,xred)

 write(message,'(a,a,a)')ch10,&
& ' --- Third order energy calculation completed --- ',ch10
 call wrtout(ab_out,message,'COLL')

!Complete missing elements using symmetry operations

 call status(0,dtfil%filstat,iexit,level,'call d3sym    ')
 call d3sym(blkflg,d3lo,indsym,mpert,natom,nsym,&
& symrec,dtset%symrel)


!Open the formatted derivative database file, and write the
!preliminary information

 if (mpi_enreg%me == 0) then

  call status(0,dtfil%filstat,iexit,level,'call ioddb8   ')
  vrsddb=010929
  dscrpt=' Note : temporary (transfer) database '
  choice=2
  ddbnm=trim(dtfil%filnam_ds(4))//'_DDB'
! tolwfr must be initialized here, but it is a dummy value
  tolwfr=1.0_dp
  call ioddb8 (choice,dscrpt,ddbnm,natom,mband,&
&  nkpt,nsym,psps%ntypat,dtfil%unddb,vrsddb,&
&  dtset%acell_orig,dtset%amu,dtset%dilatmx,dtset%ecut,dtset%ecutsm,&
&  dtset%intxc,dtset%iscf,dtset%ixc,dtset%kpt,dtset%kptnrm,&
&  natom,dtset%nband,dtset%ngfft,nkpt,nspden,nspinor,&
&  nsppol,nsym,psps%ntypat,occ,dtset%occopt,&
&  dtset%rprim_orig,dtset%sciss,dtset%spinat,dtset%symafm,dtset%symrel,&
&  dtset%tnons,tolwfr,dtset%tphysel,dtset%tsmear,&
&  dtset%typat,dtset%wtk,xred,psps%ziontypat,dtset%znucl)

  nblok=1 ; fullinit=1
  call psddb8 (choice,psps%dimekb,psps%ekb,fullinit,psps%indlmn,&
&  psps%lmnmax,psps%lnmax,nblok,&
&  psps%ntypat,dtfil%unddb,psps%pspso,psps%usepaw,psps%useylm,vrsddb)

! Call main output routine
  call d3output(blkflg,d3lo,dtset%mband,mpert,natom,dtset%nkpt,dtfil%unddb)

! Close DDB
  close(dtfil%unddb)

! Compute tensors related to third-order derivatives
  call nlopt(blkflg,carflg,d3lo,d3cart,gprimd,mpert,natom,rprimd,ucvol)

  if ((rf1pert(natom+2)==1).and.(rf2pert(natom+2)==1).and. &
&  (rf3pert(natom+2)==1)) then

   flag = 1
   i1pert = natom+2

   d3cart(:,:,i1pert,:,i1pert,:,i1pert) = &
&   d3cart(:,:,i1pert,:,i1pert,:,i1pert)*16*(pi**2)*(Bohr_Ang**2)*1.0d-8*eps0/e_Cb

   write(ab_out,*)ch10
   write(ab_out,*)' Non-linear optical susceptibility tensor d (pm/V)'
   write(ab_out,*)' in cartesian coordinates'
   write(ab_out,*)'  i1dir  i2dir  i3dir             d'

   do i1dir = 1, 3
    do i2dir = 1, 3
     do i3dir = 1, 3

      write(ab_out,'(3(5x,i2),5x,f16.9)') i1dir,i2dir,i3dir,&
&      d3cart(1,i1dir,i1pert,i2dir,i1pert,i3dir,i1pert)

      if ((blkflg(i1dir,i1pert,i2dir,i1pert,i3dir,i1pert)/=1).or.&
&      (carflg(i1dir,i1pert,i2dir,i1pert,i3dir,i1pert)/=1)) flag = 0

     end do
    end do
   end do

   if (flag == 0) then
    write(message,'(a,a,a,a,a,a)')ch10,&
&    ' d3output: WARNING -',ch10,&
&    '  matrix of third-order energies incomplete,',ch10,&
&    '  non-linear optical coefficients may be wrong.'
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
   end if

  end if  ! rf1pert,rf2pert,rf3pert

  if (((maxval(rf1pert(1:natom))/=0).and.(rf2pert(natom+2)/=0).and. &
&  (rf3pert(natom+2)/=0)).or.&
  ((maxval(rf2pert(1:natom))/=0).and.(rf1pert(natom+2)/=0).and. &
&  (rf3pert(natom+2)/=0)).or.&
  ((maxval(rf3pert(1:natom))/=0).and.(rf2pert(natom+2)/=0).and. &
&  (rf1pert(natom+2)/=0))) then
!  Perform a check if all relevant elements are available

   flag = 1
   do i1pert = 1, natom
    do i1dir = 1, 3
     do i2dir = 1, 3
      do i3dir = 1, 3
       if ((blkflg(i1dir,i1pert,i2dir,natom+2,i3dir,natom+2) /= 1).or.&
       (blkflg(i1dir,natom+2,i2dir,i1pert,i3dir,natom+2) /= 1).or.&
       (blkflg(i1dir,natom+2,i2dir,natom+2,i3dir,i1pert) /= 1)) flag = 0
       if ((carflg(i1dir,i1pert,i2dir,natom+2,i3dir,natom+2) /= 1).or.&
       (carflg(i1dir,natom+2,i2dir,i1pert,i3dir,natom+2) /= 1).or.&
       (carflg(i1dir,natom+2,i2dir,natom+2,i3dir,i1pert) /= 1)) flag = 0
      end do
     end do
    end do
   end do

   write(ab_out,*)ch10
   write(ab_out,*)' First-order change in the electronic dielectric '
   write(ab_out,*)' susceptibility tensor (Bohr^-1)'
   write(ab_out,*)' induced by an atomic displacement'
   write(ab_out,*)'  atom  displacement'

   do i1pert = 1,natom
    do i1dir = 1,3

     write(ab_out,'(1x,i4,9x,i2,3(3x,f16.9))')i1pert,i1dir,&
&     d3cart(1,i1dir,i1pert,1,natom+2,:,natom+2)
     write(ab_out,'(16x,3(3x,f16.9))')&
&     d3cart(1,i1dir,i1pert,2,natom+2,:,natom+2)
     write(ab_out,'(16x,3(3x,f16.9))')&
&     d3cart(1,i1dir,i1pert,3,natom+2,:,natom+2)

    end do

    write(ab_out,*)

   end do

   if (flag == 0) then
    write(message,'(a,a,a,a,a,a)')ch10,&
&    ' d3output: WARNING -',ch10,&
&    '  matrix of third-order energies incomplete,',ch10,&
&    '  changes in the dielectric susceptibility may be wrong.'
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
   end if


  end if  ! rf1pert,rf2pert,rf3pert

 end if   ! mpi_enreg%me

 deallocate(blkflg,carflg,cg,d3lo,d3cart,rhog,rhor)
 deallocate(kneigh,kptindex,kpt3,mvwtk,pwind,cgindex)
 deallocate(rf1pert,rf2pert,rf3pert,rfpert)
 deallocate(amass,grxc,kg,kxc,k3xc,indsym,npwarr,symrec,irrzon,phnons)


 if (mpi_enreg%paral_compil_kpt == 1) then
  deallocate(mpi_enreg%proc_distrb,mpi_enreg%kptdstrb)
 end if

!Clean the header
 call hdr_clean(hdr)

 call status(0,dtfil%filstat,iexit,level,' exit         ')
 call timab(501,2,tsec)

end subroutine nonlinear
!!***
