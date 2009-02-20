!{\src2tex{textfont=tt}}
!!****f* ABINIT/mlwfovlp
!! NAME
!! mlwfovlp
!!
!! FUNCTION
!! Routine which computes overlap M_{mn}(k,b) and projection A_{mn}(k)
!! for Wannier code (www.wannier.org f90 version).
!! Various file are written (wannier90.*) which can be used to run a
!! separate wannier calculation with the wannier90 code.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (BAmadon,FJollet,TRangel)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
!!  cprj(natom,nspinor*mband*mkmem*nsppol)= <p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                          and each |p_lmn> non-local projector
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  fermie= Fermi energy
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mgfftc=maximum size of 1D FFTs (coarse grid)
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw.
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nfft=(effective) number of FFT grid points (for this processor) (see NOTES at beginning of scfcv)
!!  ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  ucvol=unit cell volume (bohr**3)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      smatrix_pawinit,mlwfovlp_setup,mlwfovlp_pw
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine mlwfovlp(cg,cprj,dtset,ecut,eigen,fermie,gprimd,gsqcut,kg,&
& mband,mgfftc,mkmem,mpi_enreg,mpsang,mpw,natom,&
& nattyp,nfft,ngfft,nkpt,npwarr,nspden,nspinor,nsppol,ntypat,&
& pawang,pawrad,pawtab,prtvol,psps,rprimd,ucvol,xred)

 use defs_basis
 use defs_datatypes
 use defs_wannier90


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_13paw
 use interfaces_15common, except_this_one => mlwfovlp
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mgfftc,mkmem,mpsang,mpw,natom,nfft,nkpt,nspden
 integer,intent(in) :: nsppol,ntypat,prtvol
 integer,intent(inout) :: nspinor
 real(dp),intent(in) :: ecut,fermie,gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer :: kg(3,mpw*mkmem),nattyp(ntypat),ngfft(18),npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),gprimd(3,3),rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom)
 type(cprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: band_index,basis_size,cplex,i,iatom,iband,iband1,iband2,icg,ig,ikg
 integer :: ikpt,ikpt1,ikpt2,ilmn,intot,iscf,ispden,ispinor,isppol,itypat
 integer :: iun_plot,j,jband,jband1,jband2,jj,jj1,jj1tmp,jj2,jj2tmp,jj3,jj3tmp
 integer :: jj4,k,l,lmn_size,lplot,lwanniersetup,mbandw,mesh_size,mgfft,n1
 integer :: n1tmp,n2,n2tmp,n3,n3tmp,n4,n5,n6,nband_inc,nband_k,ndosfraction
 integer :: nntot,npw_k,npwin,num_bands,num_nnmax,nwan,partial_dos_flag,spacing
 integer :: tim_fourwf
 real(dp) :: fatchar,intg,uniformrandom,weight,x1,x2,xnorm,xnormb,xtemp
 logical :: gamma_only,lwannierrun,spinors
 character(len=2) :: symbol
 character(len=20) :: wfnname
 character(len=500) :: message
 character(len=9) :: seed_name
!arrays
 integer :: g1temp(3),gtemp(3),ngkpt(3)
 integer,allocatable :: g1(:,:,:),gbound(:,:),iwav(:,:,:),kg_k(:,:),ovikp(:,:)
 integer,allocatable :: proj_l(:),proj_m(:),proj_radial(:)
 real(dp) :: cgpaw(2,mpw*nspinor*mband*mkmem*nsppol),real_lattice(3,3)
 real(dp) :: recip_lattice(3,3),spreadw(3)
 real(dp),allocatable :: cm1(:,:,:,:,:),cm2_paw(:,:,:),cwavef(:,:)
 real(dp),allocatable :: denpot(:,:,:),dos_fractions(:,:,:,:)
 real(dp),allocatable :: eigenvalues_w(:,:),ff(:),fofgout(:,:),fofr(:,:,:,:)
 real(dp),allocatable :: proj_site(:,:),proj_x(:,:),proj_z(:,:),proj_zona(:)
 real(dp),allocatable :: wann_centres(:,:),wann_spreads(:),xcart(:,:)
 complex(dpc),allocatable :: M_matrix(:,:,:,:),U_matrix(:,:,:)
 complex(dpc),allocatable :: U_matrix_opt(:,:,:)
 complex(dpc),pointer :: A_matrix(:,:,:)
 logical,allocatable :: band_in(:),lwindow(:,:)
 character(len=3),allocatable :: atom_symbols(:)

!************************************************************************

!TODO 
!MAGNETISM

 write(message, '(a,a)' ) ch10,&
& '---------------------------------------------------------------'
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')
 write(message, '(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)' ) ch10,&
& '  Calculation of overlap and call to wannier90 library (beta version)',ch10,&
& '  to obtain maximally localized wannier functions ',ch10,&
& '  - wannier90.win is a mandatory secondary input',ch10,&
& '  - wannier90.wout is the output for the library',ch10,&
& '  - wannier90.amn contains projections',ch10,&
& '  - wannier90random.amn contains random projections',ch10,&
& '  - wannier90.mmn contains the overlap',ch10,&
& '  - wannier90.eig contains the eigenvalues'
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')
 write(message, '(a,a)' ) ch10,&
& '---------------------------------------------------------------'
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')

 write(message, '(a,a,a,a)' ) ch10,&
& '   mlwfovlp:  you should give k-point in the full brillouin zone ',ch10,&
& '   with explicit k-points (or kptopt=3) and istwfk 1'
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')
!********************* Keywords for options
 lwanniersetup=1 ! 1 is mandatory ( 0 is for debug)
!to use lwanniersetup=0, one would need to define which bands to exclude.
 lwannierrun=.TRUE.   ! 0 and 1 will be possible 

 gamma_only=.false. !not yet implemented
 spinors=.false. !not yet implemented
 
 if(nsppol==2) then
  write(6,*) "Warning: nsppol==2 is not supported in mlwfovlp"
  stop
 end if
 if(nspinor==2) then
  write(6,*) "Warning: nspinor==2 is not supported in mlwfovlp"
  stop
 end if
 if(mkmem==0) then
  write(6,*) "Warning: mkmem==0 is not supported in mlwfovlp"
  stop
 end if
 
!get lattice parameters in wannier90 format
 do i=1, 3
  real_lattice(:,i)=Bohr_Ang*rprimd(i,:)
  recip_lattice(:,i)=two_pi*gprimd(i,:)/Bohr_Ang
 end do 
 
!********************* Allocations.
 num_nnmax=12 !limit fixed for compact structure in wannier_setup.
 allocate(g1(3,nkpt,num_nnmax),ovikp(nkpt,num_nnmax))
 allocate(atom_symbols(natom))
 allocate(xcart(3,natom))
 allocate(band_in(mband))
 allocate(proj_site(3,mband),proj_l(mband),proj_m(mband),proj_radial(mband))
 allocate(proj_x(3,mband),proj_z(3,mband),proj_zona(mband))



!********************* Wannier setup 
 if(psps%npsp/=psps%ntypat) then
  write(6,*) "prb npsp"
  stop
 end if
 nullify(A_matrix)
 call mlwfovlp_setup(atom_symbols,band_in,dtset,eigen,gamma_only,&
& g1,gprimd,lwanniersetup,mband,mbandw,mkmem,mpw,natom,nattyp,nband_inc,nkpt,&
& nntot,nsppol,nspinor,ntypat,num_bands,num_nnmax,nwan,ovikp,&
& proj_l,proj_m,proj_radial,proj_site,proj_x,proj_z,proj_zona,&
& real_lattice,recip_lattice,rprimd,spinors,xcart,xred)

!write(300,*) 'proj_l ',proj_l
!write(300,*) 'proj_m ',proj_m
!write(300,*) 'proj_radial ',proj_radial
!write(300,*) 'proj_site ',proj_site
!write(300,*) 'proj_x ',proj_x
!write(300,*) 'proj_z ',proj_z
!write(300,*) 'proj_zona ',proj_zona



 write(message, '(a,a,a,a)' ) ch10,&
& '   mlwfovlp :  mlwfovlp_setup done -',ch10,&
& '   see wannier90.wout for details.'
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')

!********************* Allocate Matrices for wannier_run
 if(lwannierrun) then
  allocate(eigenvalues_w(num_bands,nkpt))
  allocate(M_matrix(num_bands,num_bands,nntot,nkpt))
 end if

!********************* Write Eigenvalues 
 open(unit=444,file='wannier90.eig',form='formatted',status='unknown')
 band_index=zero
 do isppol=1,nsppol
  do ikpt=1,nkpt
   nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
   jband=0
   do iband=1,mbandw
    if(band_in(iband)) then
     jband=jband+1
     write(444, '(2i6,4x,f10.5)' ) jband,ikpt,Ha_eV*eigen(iband+band_index)
     eigenvalues_w(jband,ikpt)=Ha_eV*eigen(iband+band_index)
    end if
   end do
   band_index=band_index+nband_k
  end do
 end do
 close(444)
 write(message, '(a,a)' ) ch10,&
& '   mlwfovlp :  eigenvalues written'
 call wrtout(06,  message,'COLL')

!********************* compute PW Contribution
 allocate(cm1(2,nkpt,nntot,mbandw,mbandw))
 allocate(iwav(nsppol,nkpt,mband))

 call mlwfovlp_pw(cg,cm1,dtset,g1,iwav,kg,mband,mband,mkmem,mpsang,mpw,natom,&
& nfft,ngfft,nkpt,nntot,npwarr,nspden,nspinor,nsppol,ntypat,ovikp,prtvol)

 write(message, '(a,a)' ) ch10,&
& '   mlwfovlp : PW part of overlap computed   '
 call wrtout(06,  message,'COLL')
!********************* compute PAW Contribution and add it to PW contribution
 111 continue
 if(psps%usepaw==1) then
  write(message, '(a,a)' ) ch10,&
&  '** smatrix_pawinit : PAW part of overlap  '
  call wrtout(06,  message,'COLL')
  allocate(cm2_paw(2,mbandw,mbandw))
  do ikpt1=1,nkpt
   write(message, '(a,i6)' ) &
&   '   compute PAW part of overlaps for k-point',ikpt1  
   call wrtout(06,  message,'COLL')
   do intot=1,nntot
    ikpt2= ovikp(ikpt1,intot)
    g1temp(:)=g1(:,ikpt1,intot)
    call smatrix_pawinit(cm2_paw,cprj,ikpt1,ikpt2,g1temp,gprimd,&
&    dtset%kpt,mbandw,mbandw,&
&    mkmem,natom,nattyp,dtset%nband,&
&    nkpt,nspinor,nsppol,dtset%ntypat,pawang,pawrad,pawtab,rprimd,&
&    dtset%typat,psps%usepaw,xred)
    cm1(:,ikpt1,intot,:,:)=cm1(:,ikpt1,intot,:,:)+four_pi*cm2_paw(:,:,:)
!   cm1(:,ikpt1,intot,:,:)=four_pi*cm2_paw(:,:,:)
   end do ! intot
  end do ! ikpt1
  deallocate(cm2_paw)
  write(message, '(a,a)' ) ch10,&
&  '   mlwfovlp : PAW part of overlap computed '
  call wrtout(06,  message,'COLL')
 end if ! usepaw

!********************* write overlap for separate calculation of wannier functions

 open(unit=221,file='wannier90.mmn',form='formatted',status='unknown')
 write(221,*) "nnkp version 90"
 write(221,*) num_bands,nkpt,nntot
 do ikpt1=1,nkpt
  do intot=1,nntot
   write(221,'(2i6,3x,3x,3i5)') ikpt1,ovikp(ikpt1,intot),(g1(jj,ikpt1,intot),jj=1,3)
   jband2=0
   do iband2=1,mbandw ! the first index is faster
    if(band_in(iband2)) then
     jband2=jband2+1
     jband1=0
     do iband1=1,mbandw
      if(band_in(iband1)) then
       jband1=jband1+1
       write(221,*) cm1(1,ikpt1,intot,iband1,iband2),cm1(2,ikpt1,intot,iband1,iband2)
       M_matrix(jband1,jband2,intot,ikpt1)=&
&       cmplx(cm1(1,ikpt1,intot,iband1,iband2),cm1(2,ikpt1,intot,iband1,iband2))
!      write(2211,*) ikpt1,intot,iband1,iband2
!      write(2211,*) cm1(1,ikpt1,intot,iband1,iband2),cm1(2,ikpt1,intot,iband1,iband2)
      end if
     end do
    end if
   end do
  end do
 end do
 close(221)
 write(message, '(a)' )  '   wannier90.mmn written'
 call wrtout(06,  message,'COLL')
 deallocate(ovikp,g1,cm1)

!******************** Calculate projections ************

 allocate(A_matrix(num_bands,nwan,nkpt))
 if(dtset%w90iniprj/=0) then
  call mlwfovlp_proj(A_matrix,band_in,cg,cprj,dtset,eigen,gprimd,kg,&
&  dtset%istwfk,iwav,dtset%w90iniprj,mband,mkmem,mpi_enreg,mpw,natom,&
&  nattyp,nkpt,npwarr,nspden,&
&  nspinor,nsppol,ntypat,num_bands,nwan,pawtab,proj_l,proj_m,&
&  proj_radial,proj_site,proj_x,proj_z,proj_zona,prtvol,psps,ucvol)
  write(message, '(a,a,a,a)' ) ch10,&
&  '   mlwfovlp_setup :  mlwfovlp_proj done -',ch10,&
&  '   Projectors computed.'
  call wrtout(06,  message,'COLL')

 else
  write(message,*) " Warning: mlwfovlp: no projection computed:", &
&  "choose at least random projections"
  call wrtout(6,message,'COLL')
 end if !dtset%w90iniprj/=0


 deallocate(proj_site,proj_l,proj_m,proj_radial)
 deallocate(proj_x,proj_z,proj_zona)


!********************* write files for wannier function plot
 
  spacing = 2


 if(psps%usepaw==1.and.dtset%w90prtunk==1) then
  write(message, '( a,a,a,a,a,a,a,a,a)')ch10,&
  &"   Warning: The UNK matrices will not contain the correct wavefunctions ",ch10,&
  &"   since we are just writing the plane wave contribution.",ch10,&
  &"   The contribution from inside the spheres is missing. ",ch10,&
  &"   However, these files can be used for plotting purposes",ch10
  call wrtout(06,  message,'COLL')
 end if

 if( dtset%w90prtunk==1) then
 write(message, '( a,a,a,a,a,a,a,a,i3,a,a)')ch10,&
 &"   UNK files will be written.",ch10,&
 &"   Warning: In order to reduce the size of the files",ch10,&
 &"   we are not writting all the wavefunctions but",ch10, &
 &"   we are writing every", spacing,"records.",ch10
  call wrtout(06,  message,'COLL')

  allocate(kg_k(3,mpw))
  n1=ngfft(1)
  n2=ngfft(2)
  n3=ngfft(3)
  n4=ngfft(4)
  n5=ngfft(5)
  n6=ngfft(6)

  cplex=1
  ikg=0 
  mgfft=mgfftc ! error
  do isppol=1,nsppol
   do ikpt=1,nkpt
    npw_k=npwarr(ikpt)
    kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
    allocate(denpot(cplex*n4,n5,n6),cwavef(2,npw_k),&
&    fofr(2,n4,n5,n6),gbound(2*mgfft+8,2))
    allocate(fofgout(2,npw_k))
    iun_plot=1000+ikpt
    write(wfnname,222) ikpt, isppol
!   open (unit=iun_plot, file=wfnname,form='formatted')
    open(unit=iun_plot, file=wfnname,form='unformatted')
    222   format ('UNK',i5.5,'.',i1)
!   optimizing grid for UNK files
    n1tmp = n1/spacing
    n2tmp = n2/spacing
    n3tmp = n3/spacing
    if( mod(n1,spacing) /= 0) then
     n1tmp = n1tmp + 1
    end if
    if( mod(n2,spacing) /= 0) then
     n2tmp = n2tmp + 1
    end if
    if( mod(n3,spacing) /= 0) then
     n3tmp = n3tmp + 1
    end if
!   write(iun_plot,*) n1tmp,n2tmp,n3tmp,ikpt,nband_inc
    write(iun_plot) n1tmp,n2tmp,n3tmp,ikpt,nband_inc
!   gbound=zero
    call sphereboundary(gbound,dtset%istwfk(ikpt),kg_k,mgfft,npw_k)
    write(6,*) "  writes UNK file for ikpt=",ikpt,dtset%istwfk(ikpt)
    denpot(:,:,:)=zero
    weight = one
    do iband=1,mbandw
     if(band_in(iband)) then
      do ig=1,npw_k*nspinor
       cwavef(1,ig)=cg(1,ig+iwav(isppol,ikpt,iband))
       cwavef(2,ig)=cg(2,ig+iwav(isppol,ikpt,iband))
      end do
      tim_fourwf=0
      call fourwf(cplex,denpot,cwavef,fofgout,fofr,&
&      gbound,gbound,dtset%istwfk(ikpt),kg_k,kg_k,mgfft,&
&      mpi_enreg,1,ngfft,npw_k,npw_k,n4,n5,n6,0,dtset%paral_kgb,&
&      tim_fourwf,weight)
!     do jj3=1,n3,spacing
!     do jj2=1,n2,spacing
!     do jj1=1,n1,spacing
!     write(iun_plot,*) fofr(1,jj1,jj2,jj3),&
!     & fofr(2,jj1,jj2,jj3)
!     end do !jj1
!     end do !jj2
!     end do !jj3
!     unformatted (must be one record)
      write(iun_plot) (((fofr(1,jj1,jj2,jj3),fofr(2,jj1,jj2,jj3),&
&      jj1=1,n1,spacing),jj2=1,n2,spacing),jj3=1,n3,spacing)

     end if !iband
    end do ! iband
    deallocate(cwavef,fofr,gbound,denpot,fofgout)
    ikg=ikg+npw_k
    close(iun_plot)
   end do  ! ikpt
  end do  ! nsppol
  deallocate(kg_k)
 end if !dtset%w90prtunk
 deallocate(iwav)


!********************* compute wannier function if lwannierrun==1

 if(lwannierrun) then
  write(message, '(a,a,a,a,a,a)' ) ch10,&
&  '** mlwfovlp :   call wannier90 library subroutine wannier_run ',ch10,&
&  '   Calculation is running         ',ch10,&
&  '   see wannier90.wout for details.'
  call wrtout(06,  message,'COLL')
  if(lwanniersetup.ne.1) stop
  allocate(U_matrix(nwan,nwan,nkpt))
  allocate(U_matrix_opt(num_bands,nwan,nkpt))
  allocate(lwindow(num_bands,nkpt))
  allocate(wann_centres(3,nwan))
  allocate(wann_spreads(nwan))
  seed_name="wannier90"
! write(6,*) seed_name
! write(6,*) ngkpt
  ngkpt(1)=dtset%kptrlatt(1,1)
  ngkpt(2)=dtset%kptrlatt(2,2) ! ajouter test de verif que kptrlatt est bien diagonal
  ngkpt(3)=dtset%kptrlatt(3,3)
! write(6,*) nkpt
! write(6,*) rprimd*Bohr_Ang
! write(6,*) two_pi*gprimd/Bohr_Ang
! write(6,*) mband
! write(6,*) "nwan",nwan
! write(6,*) nntot
! write(6,*) natom
! write(6,*) atom_symbols
! write(6,*) xcart
! write(6,*) num_bands,num_bands,nntot,nkpt
! write(6,*) wann_spreads
! wann_spreads=2
! do i=1, nkpt
! do j=1, nntot
! write(6,*) i,j
! do k=1, num_bands
! do l=1, num_bands
! write(6,*) "m",M_matrix(l,k,j,i)
! enddo
! enddo
! enddo
! enddo


#if defined HAVE_WANNIER90
  call wannier_run(seed_name,ngkpt,nkpt,&            !input
& real_lattice,recip_lattice,dtset%kpt,num_bands,& !input
& nwan,nntot,natom,atom_symbols,&                  !input
& xcart*Bohr_Ang,gamma_only,M_matrix,A_matrix,eigenvalues_w,& !input
& U_matrix,U_matrix_opt,lwindow_loc=lwindow,wann_centres_loc=wann_centres,&     !output
& wann_spreads_loc=wann_spreads,spread_loc=spreadw)                            !output
#endif
  deallocate(M_matrix,A_matrix,eigenvalues_w,wann_centres,wann_spreads)
 end if !lwannierrun
 write(message, '(a,a,a,a,a)' ) ch10,&
& '   mlwfovlp :  mlwfovlp_run completed -',ch10,&
& '   see wannier90.wout for details.',ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')


!********************* deallocation 

 deallocate(lwindow,band_in)
 deallocate(U_matrix,U_matrix_opt)
 deallocate(atom_symbols)
 deallocate(xcart)

!DEBUG
!write(6,*)' mlwfovlp : exit'
!stop
!ENDDEBUG

 end subroutine mlwfovlp
!!***
