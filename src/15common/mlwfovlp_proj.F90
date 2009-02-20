!{\src2tex{textfont=tt}}
!!****f* ABINIT/mlwfovlp_proj
!! NAME
!! mlwfovlp_proj
!!
!! FUNCTION
!! Routine which computes projection A_{mn}(k)
!! for Wannier code (www.wannier.org f90 version).
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (BAmadon,FJollet,TRangel,drh)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*nkpt*nsppol)=planewave coefficients of wavefunctions
!!  cprj(natom,nspinor*mband*mkmem*nsppol)= <p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                          and each |p_lmn> non-local projector
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  istwfk=option parameter that describes the storage of wfs
!!  iwav(nsppol,nkpt,mbandw): shift for pw components in cg.
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  lproj= flag 0: no projections, 1: random projections, 2: projections on projectors
!!              3: projections on wavefunctions.
!!  mband=maximum number of bands
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum dimensioned size of npw.
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nkpt=number of k points.
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  num_bands=number of bands actually used to construct the wannier function
!!  nwan= number of wannier fonctions (read in wannier90.win).
!!  proj_l(mband)= angular part of the projection function (quantum number l)
!!  proj_m(mband)= angular part of the projection function (quantum number m)
!!  proj_radial(mband)= radial part of the projection.
!!  proj_site(3,mband)= site of the projection.
!!  proj_x(3,mband)= x axis for the projection.
!!  proj_z(3,mband)= z axis for the projection.
!!  proj_zona(mband)= extension of the radial part.
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!! OUTPUT
!!  A_matrix(num_bands,nwan,nkpt)= Matrix of projections needed by wannier_run
!!  ( also wannier90random.amn is written)
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      mlwfovlp_setup
!!
!! CHILDREN
!!      (none)
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine mlwfovlp_proj(A_matrix,band_in,cg,cprj,dtset,eigen,gprimd,kg,&
&istwfk,iwav,lproj,mband,mkmem,mpi_enreg,mpw,natom,nattyp,&
&nkpt,npwarr,nspden,nspinor,&
&nsppol,ntypat,num_bands,nwan,pawtab,proj_l,proj_m,proj_radial,&
&proj_site,proj_x,proj_z,proj_zona,prtvol,psps,ucvol)

 use defs_basis
 use defs_datatypes
 use defs_wannier90


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_15common, except_this_one => mlwfovlp_proj
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lproj,mband,mkmem,mpw,natom,nkpt,nspden,nspinor,nsppol
 integer,intent(in) :: ntypat,num_bands,nwan,prtvol
 complex(dpc),parameter :: c0=(0._dp,0._dp),c1=(1._dp,0._dp),ci=(0._dp,1._dp)
 complex(dpc) :: cstr_fact
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer :: iwav(nsppol,nkpt,mband),nattyp(ntypat)
 integer,intent(in) :: istwfk(nkpt),kg(3,mpw*mkmem),npwarr(nkpt),proj_l(mband)
 integer,intent(in) :: proj_m(mband),proj_radial(mband)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),gprimd(3,3),proj_site(3,mband)
 real(dp),intent(in) :: proj_x(3,mband),proj_z(3,mband),proj_zona(mband)
 complex(dpc),intent(out) :: A_matrix(num_bands,nwan,nkpt)
 logical,intent(in) :: band_in(mband)
 type(cprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: i,iatom,iatprjn,iband,iband1,iband2,ibg,ibnd1,icat,icg,icg_shift
 integer :: idum,ii,ikg,ikpt,il,il2,ilmn,index,inversion_flag,ipw,ispden
 integer :: ispinor,isppol,itypat,iwan,j,jband,jj,jj1,k,libprjn,lm,lmax,lmax2
 integer :: lmn_size,mesh_size,mr,natprjn,nband_k,nbprjn,ndosfraction,npw_k
 integer :: partial_dos_flag,prtdos,sumtmp
 real(dp),parameter :: qtol=2.0d-8
 real(dp) :: alat,anorm,arg,cg1,cg2,doti,dotr,norm,norm_error,norm_error_bar
 real(dp) :: norm_error_tmp,ucvol,x1,x2,xnorm,xnormb,xx,ylm,yy,zz
 complex(dpc) :: ZDOTC,amn_tmp,gf_trial,lphase,sk
 character(len=500) :: message
 character(len=fnlen) :: fildata
 type(pseudopotential_type) :: pspsphi
!arrays
 integer :: kg_k(3,mpw)
 integer,allocatable :: indgk(:),lprjn(:),npprjn(:)
 real(dp) :: a1(3),a2(3),a3(3),gk_rot(3),kpg(3),kpt(3),umat(3,3),utmp(3,3)
 real(dp),allocatable :: amn(:,:,:,:),amn2(:,:,:,:,:,:,:),cgr(:,:)
 real(dp),allocatable :: chiphimphitint(:,:),dos_fractions(:,:,:,:),ff(:)
 real(dp),allocatable :: gsum2(:),kpg2(:),pjPsi(:,:,:,:,:,:,:),radial(:)
 complex(dpc),allocatable :: cgc(:,:),gf(:,:),gft_lm(:),orb_lm(:,:,:)
 complex(dpc),allocatable :: wan_lm(:,:),ylmcp(:)
!no_abirules
!Tables 3.1 & 3.2, User guide
 integer,save :: orb_l_defs(-5:3)=(/2,2,1,1,1,0,1,2,3/) 
 integer,parameter :: mtransfo(0:3,7)=&
&  reshape((/1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,-2,-1,2,1,0,0,0,-1,1,2,-2,-3,3/),(/4,7/))

!************************************************************************


!TODO AND POSSIBLE BUGS:
!MAGNETISM

 if ((lproj/=1).and.(lproj/=4)) then
  write(message, '(a,a,a,a,a,a)' ) ch10,&
&  '  mlwfovlp_proj : ERROR -',ch10,&
&  '  Value of lproj no allowed ',ch10,&
&  '  Action : change lproj.'
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

 write(message, '(a,a)' )ch10,&
& '** mlwfovlp_proj:  compute A_matrix of initial guess for wannier functions'
 call wrtout(6,message,'COLL')

!********************* Write Random projectors
 if(lproj==1) then
  open(unit=220,file='wannier90random.amn',form='formatted',status='unknown')
  idum=123456
! Compute random projections
  allocate(amn(2,mband,nkpt,nwan))
  amn=zero
  write(220,*) "Random Projections from Abinit"
  write(220,*) num_bands,nkpt,nwan
  do ikpt=1,nkpt
   do iband1=1,mband
    xnormb=0.d0
    do iband2=1,nwan
     x1=uniformrandom(idum)
     x2=uniformrandom(idum)
     xnorm=sqrt(x1**2+x2**2)
     xnormb=xnormb+xnorm
     amn(1,iband1,ikpt,iband2)=x1
     amn(2,iband1,ikpt,iband2)=x2
    end do
    do iband2=1,nwan
     amn(1,iband1,ikpt,iband2)=amn(1,iband1,ikpt,iband2)/xnormb
     amn(2,iband1,ikpt,iband2)=amn(1,iband1,ikpt,iband2)/xnormb
    end do
   end do
  end do
  do ikpt=1,nkpt
   do iband2=1,nwan
    jband=0
    do iband1=1,mband
     if(band_in(iband1)) then
      jband=jband+1
      write(220,'(3i6,3x,3x,2f11.7)')iband1,iband2,ikpt,amn(1,iband1,ikpt,iband2),amn(2,iband1,ikpt,iband2)
      if(jband.gt.num_bands) then
       write(message, '(a,a,a,a,a,a)' ) ch10,&
&       '  mlwfovlp_proj : ERROR -',ch10,&
&       '  Value of jband is above num_bands ',ch10,&
&       '  Action : contact Abinit group'
       call wrtout(06,  message,'COLL')
       call leave_new('COLL')
      end if
      A_matrix(jband,iband2,ikpt)=cmplx(amn(1,iband1,ikpt,iband2),amn(2,iband1,ikpt,iband2))
     end if
    end do
   end do
  end do
  close(220)
  deallocate(amn)
 end if


!*************** computes projection  from PROJECTORS ********************
 if(lproj==2) then  !! if LPROJPRJ
! ----- set values for projections --------------------- ! INPUT
! nbprjn:number of  different l-values for projectors
! lprjn: value of l for each projectors par ordre croissant
! npprjn: number of projectors for each lprjn
  natprjn=1  ! atoms with wannier functions are first
  if(natprjn/=1) then ! in this case lprjn should depend on iatprjn
   stop
  end if
  nbprjn=2
  allocate(lprjn(nbprjn))
  lprjn(1)=0
  lprjn(2)=1
  allocate(npprjn(0:lprjn(nbprjn)))
  npprjn(0)=1
  npprjn(1)=1
! --- test coherence of nbprjn and nwan
  sumtmp=0
  do iatprjn=1,natprjn
   do libprjn=0,lprjn(nbprjn)
    sumtmp=sumtmp+(2*libprjn+1)*npprjn(libprjn)
   end do
  end do
  if(sumtmp/=nwan) then
   write(6,*) "Number of Wannier orbitals is not equal to number of projections"
   write(6,*) "Action: check values of lprjn,npprjn % nwan"
   write(6,*) "nwan, sumtmp=",nwan,sumtmp
   stop
  end if
! --- end test of coherence
  allocate(amn2(2,natom,nsppol,nkpt,mband,nspinor,nwan))
  if(psps%usepaw==1) then
   open(unit=219,file='wannier90.amn',form='formatted',status='unknown')
   amn2=zero
   ibg=0
   do isppol=1,nsppol
    do ikpt=1,nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
     do iband=1,nband_k
!     write(6,*)"amn2",iband,ibg,ikpt
      do ispinor=1,nspinor
       icat=1
       do itypat=1,dtset%ntypat
        lmn_size=pawtab(itypat)%lmn_size
        do iatom=icat,icat+nattyp(itypat)-1
         jj1=0
         do ilmn=1,lmn_size
          if(iatom.le.natprjn) then
!          do iwan=1,nwan
           do libprjn=0,lprjn(nbprjn)
!           if (psps%indlmn(1,ilmn,itypat)==proj_l(iwan)) then
!           if (psps%indlmn(2,ilmn,itypat)==mtransfo(proj_l(iwan),proj_m(iwan))) then
            if (psps%indlmn(1,ilmn,itypat)==libprjn) then
             if (psps%indlmn(3,ilmn,itypat)<=npprjn(libprjn)) then
              if(band_in(iband)) then
               jj1=jj1+1
               if(jj1>nwan) then
                write(6,*) "number of wannier orbitals is lower than lmn_size"
                write(6,*) jj1,nwan
                stop
               end if
               amn2(1,iatom,isppol,ikpt,iband,ispinor,jj1)=cprj(iatom,iband+ibg)%cp(1,ilmn)
               amn2(2,iatom,isppol,ikpt,iband,ispinor,jj1)=cprj(iatom,iband+ibg)%cp(2,ilmn)
              end if
             end if
            end if
           end do ! libprjn
!          endif
!          endif
!          enddo ! iwan
          end if ! natprjn
         end do !ilmn
        end do ! iatom
        icat=icat+nattyp(itypat)
       end do ! itypat
      end do ! ispinor
     end do !iband
     ibg=ibg+nband_k*nspinor
!    write(6,*)'amn2b',iband,ibg,ikpt
    end do !ikpt
   end do ! isppol

!  -----------------------  write      projections    --------------------

   write(219,*) "Projections from Abinit : mband,nkpt,nwan. indices: iband1,iwan,ikpt"
   write(219,*) num_bands,nkpt,nwan
   do ikpt=1,nkpt
    do iband2=1,nwan
     jband=0
     do iband1=1,mband
      if(band_in(iband1)) then
       jband=jband+1
!      write(219,*) iband1,ikpt,iband2
       write(219,'(3i6,3x,3x,2f18.14)')jband,iband2,ikpt,amn2(1,1,1,ikpt,iband1,1,iband2),amn2(2,1,1,ikpt,iband1,1,iband2)
!      write(219,*)iband1,iband2,ikpt,amn(1,iband1,ikpt,iband2),amn(2,iband1,ikpt,iband2)
       A_matrix(jband,iband2,ikpt)=&
&       cmplx(amn2(1,1,1,ikpt,iband1,1,iband2),amn2(2,1,1,ikpt,iband1,1,iband2))
      end if
     end do
    end do
   end do

  end if !usepaw
  close(219)
  deallocate(amn2)
  deallocate(npprjn,lprjn)

 end if ! lproj==2

!#########################################################
 if( lproj == 4) then !based on .win file
! #########################################################
! obtain lmax and lmax2
  lmax=0
  do iwan=1,nwan
   lmax=max(lmax,orb_l_defs(proj_l(iwan)))
  end do !iwan
  lmax2=(lmax+1)**2

! Allocate arrays
  allocate(gf(mpw,nwan),gft_lm(lmax2),gsum2(nwan))
  allocate(kpg2(mpw),radial(lmax2))
  allocate(wan_lm(lmax2,nwan),ylmcp(lmax2))

! get ylmfac, factor used for rotations and hybrid orbitals
  call get_wan_lm(wan_lm,lmax,lmax2,mband,nwan,proj_l,proj_m,proj_x,proj_z)

  norm_error=zero
  norm_error_bar=zero
  ikg=0
  icg=0

  do ikpt=1, mkmem
   write(message, '(a,i6)' ) &
&   '   compute projections for k-point=',ikpt
   call wrtout(06,  message,'COLL')

!  Initialize variables
   npw_k=npwarr(ikpt)
   gsum2(:)=0.d0
   gf(:,:) = (0.d0,0.d0)
   kpt(:)=dtset%kpt(:,ikpt)
   kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)

   do ipw=1, npw_k
    kpg(1)= (kpt(1) + real(kg_k(1,ipw),dp))     !k+G
    kpg(2)= (kpt(2) + real(kg_k(2,ipw),dp))
    kpg(3)= (kpt(3) + real(kg_k(3,ipw),dp))

!   Calculate modulus of k+G
    xx=gprimd(1,1)*kpg(1)+gprimd(1,2)*kpg(2)+gprimd(1,3)*kpg(3)
    yy=gprimd(2,1)*kpg(1)+gprimd(2,2)*kpg(2)+gprimd(2,3)*kpg(3)
    zz=gprimd(3,1)*kpg(1)+gprimd(3,2)*kpg(2)+gprimd(3,3)*kpg(3)
    kpg2(ipw)= two_pi*sqrt(xx**2+yy**2+zz**2)

!   Complex Y_lm for k+G
    if(lmax==0) then
     ylmcp(1)=c1/sqrt(four_pi)
    else
     call ylm_cmplx(lmax,ylmcp,xx,yy,zz)
    end if

    do iwan=1,nwan
!    obtain radial part
     call mlwfovlp_radial(proj_zona(iwan),lmax,lmax2,radial,proj_radial(iwan)&
&     ,kpg2(ipw))

!    scale complex representation of projector orbital with radial functions
!    of appropriate l
     gft_lm(:)=radial(:)*wan_lm(:,iwan)

!    complex structure factor for projector orbital position
     arg = ( kpg(1)*proj_site(1,iwan) + kpg(2)*proj_site(2,iwan) + &
&     kpg(3)*proj_site(3,iwan) ) * 2*pi
     cstr_fact = cmplx(cos(arg), -sin(arg) )

!    obtain guiding functions
     gf(ipw,iwan)=cstr_fact*dot_product(ylmcp(:),gft_lm)

     gsum2(iwan)=gsum2(iwan)+real(gf(ipw,iwan))**2+aimag(gf(ipw,iwan))**2
    end do !iwan
   end do !ipw

   do iwan=1,nwan
    gsum2(iwan)=16._dp*pi**2*gsum2(iwan)/ucvol
    gf(:,iwan)=gf(:,iwan)/sqrt(gsum2(iwan))
    norm_error=max(abs(gsum2(iwan)-one),norm_error)
    norm_error_bar=norm_error_bar+(gsum2(iwan)-one)**2
   end do !iwan

!  ! Guiding functions are computed.
!  compute overlaps of gaussian projectors and wave functions
   do iwan=1,nwan
    jband=0
    do iband=1,mband
     if(band_in(iband)) then
      icg_shift=npw_k*nspinor*(iband-1)+icg
      jband=jband+1
      amn_tmp=cmplx(0.d0,0.d0)
      do ipw=1,npw_k
       amn_tmp=amn_tmp+gf(ipw,iwan)*cmplx(cg(1,ipw+icg_shift),-cg(2,ipw+icg_shift))
      end do !ipw
      A_matrix(jband,iwan,ikpt)=amn_tmp
     end if !band_in
    end do !iband
   end do !iwan
   icg=icg+npw_k*nspinor*mband
   ikg=ikg+npw_k
  end do !ikpt

  norm_error_bar=sqrt(norm_error_bar/real(nkpt*nwan,dp))
  if(norm_error>0.05_dp) then
   write(message, '(6a,f6.3,a,f6.3,10a)' )ch10,&
&   ' mlwfovlp_proj : WARNING',ch10,&
&   '  normalization error for wannier projectors',ch10,&
&   '  is',norm_error_bar,' (average) and',norm_error,' (max).',ch10,&
&   '  this may indicate more cell-to-cell overlap of the radial functions',ch10,&
&   '  than you want.',ch10,&
&   '  Action : modify zona (inverse range of radial functions)',ch10,&
   '  under "begin projectors" in wannier90.win file',ch10
   call wrtout(6,message,'COLL')
  end if

! Deallocations
  deallocate(gf,gft_lm,gsum2)
  deallocate(kpg2,radial)
  deallocate(wan_lm,ylmcp)

! -----------------------  write      projections    --------------------
  open(unit=219,file='wannier90.amn',form='formatted',status='unknown')

  write(219,*) 'Projections from Abinit : mband,nkpt,nwan. indices: iband1,iwan,ikpt'
  write(219,*) num_bands,nkpt,nwan
  do ikpt=1,nkpt
   do iwan=1,nwan
    jband=0
    do iband=1,mband
     if(band_in(iband)) then
      jband=jband+1
      write(219,'(3i6,13x,3x,2f18.14)')jband,iwan,ikpt,A_matrix(jband,iwan,ikpt)
     end if !band_in
    end do !iband
   end do !iwan
  end do !ikpt
  close(219)

  write(message, '(a)' ) &
&  '   wannier90.amn written'
  call wrtout(06,  message,'COLL')

 end if !lproj==4

 contains

!! NAME
!! get_wan_lm
!!
!! FUNCTION
!! Routine that produces a factor by which the initial
!! guess of functions will be multiplied for the Wannier90 interface.
!! It is just used if there are rotations, or if the functions required
!! are linear combinations of the ylm real functions.
!! Example,
!! For a function G(r)= 1/2 s + 1/3 px - 1/2 pz
!!   it would produce a matrix of the following form:
!!   [1/2,-1/2,1/3,0,0...0]
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (TRangel,drh)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  lmax2=number of ylm functions
!!  mband=maximum number of bands
!!  nwan = number of wannier functions
!!  proj_l(mband)= angular part of the projection function (quantum number l)
!!  proj_m(mband)= angular part of the projection function (quantum number m)
!!  proj_x(3,mband)= x axis for the projection.
!!  proj_z(3,mband)= z axis for the projection.
!!
!! OUTPUT
!!  ylmfac=matrix containig a factor for ylm hybrid orbitals
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      mlwfovlp_proj
!!
!! CHILDREN
!!      leave_new,rotmat,wrtout,ylm_cmplx,zgesv
!!
!! SOURCE

  subroutine get_wan_lm(wan_lm,lmax,lmax2,mband,nwan,proj_l,proj_m,proj_x,proj_z)



!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_lib00numeric
!End of the abilint section

  implicit none

! Arguments
  integer, intent(in):: lmax,lmax2,nwan,mband
! arrays
  integer,save :: orb_idx(16)=(/1,3,4,2,7,8,6,9,5,13,14,12,15,11,16,10/) !Tab3.1
  integer,intent(in) :: proj_l(mband),proj_m(mband)
  real(dp),intent(in) :: proj_x(3,mband),proj_z(3,mband)
  complex(dp),intent(out)::wan_lm(lmax2,nwan)

! local variables
  integer :: idum,il,info,ir,ll,lm,lmc,mm,mr
  real(dp):: onem,test

! arrays
  integer:: ipiv(lmax2)
  real(dp)::r(3,lmax2),rp(3,lmax2)
  real(dp)::rs2,rs3,rs6,rs12
  real(dp)::mly_trial(lmax2,lmax2),ylm_trial(lmax2,lmax2),ylm_wan(lmax2)
  complex(dp)::crot(lmax2,lmax2),ctor(lmax2,lmax2),orb_lm(lmax2,-5:3,7)
  complex(dp):: ylmcp(lmax2)
  complex(dp):: ylmc_rr(lmax2,lmax2),ylmc_rr_save(lmax2,lmax2)
  complex(dp):: ylmc_rrinv(lmax2,lmax2),ylmc_rp(lmax2,lmax2)
  complex(dp),parameter :: c0=(0._dp,0._dp),c1=(1._dp,0._dp),ci=(0._dp,1._dp)


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! constants for linear combinations of ylm's
  rs2=1._dp/sqrt(2._dp)
  rs3=1._dp/sqrt(3._dp)
  rs6=1._dp/sqrt(6._dp)
  rs12=1._dp/sqrt(12._dp)

! complex lm coefficients for real spherical harmonics in conventional order
! s, py,pz,px, dxy,dyz,dz2,dxz,dx2-y2, fy(3x2-y2),fxyz,fyz2,fz3,fxz2,
! fz(x2-y2),fx(x2-3y2)
  ctor(:,:)=c0
  do ll=0,lmax
   mm=0
   lm= ll**2+ll+mm+1
   ctor(lm,lm)=c1
   if(ll>0) then
    onem=one
    do mm=1,ll
     onem=-onem !(-1^mm)
     lm= ll**2+ll+mm+1
     lmc=ll**2+ll-mm+1
     ctor(lm ,lm )=rs2*c1
     ctor(lmc,lm )=onem*rs2*c1
     ctor(lm ,lmc)=rs2*ci
     ctor(lmc,lmc)=-onem*rs2*ci
    end do
   end if
  end do

  lm=0
  do ll=0,lmax
   do mm=-ll,ll
    lm=lm+1
    ctor(:,lm)=ctor(:,lm)*conjg(ci)**ll
   end do !mm
  end do !ll


! coefficients for basic wannier orbitals in Table 3.1 order
  orb_lm(:,:,:)=c0
  ii=0
  do ll=0,lmax
   do mr=1,2*ll+1
    ii=ii+1
    orb_lm(:,ll,mr)=ctor(:,orb_idx(ii))
   end do
  end do

! coefficients for linear combinations in table 3.2 order
  if(lmax>=1) then
!  s            px
   orb_lm(:,-1,1)=rs2*ctor(:,1)+rs2*ctor(:,4)
   orb_lm(:,-1,2)=rs2*ctor(:,1)-rs2*ctor(:,4)
!  s            px            py
   orb_lm(:,-2,1)=rs3*ctor(:,1)-rs6*ctor(:,4)+rs2*ctor(:,2)
   orb_lm(:,-2,2)=rs3*ctor(:,1)-rs6*ctor(:,4)-rs2*ctor(:,2)
   orb_lm(:,-2,3)=rs3*ctor(:,1)+2._dp*rs6*ctor(:,4)
!  s        px        py        pz
   orb_lm(:,-3,1)=half*(ctor(:,1)+ctor(:,4)+ctor(:,2)+ctor(:,3))
   orb_lm(:,-3,2)=half*(ctor(:,1)+ctor(:,4)-ctor(:,2)-ctor(:,3))
   orb_lm(:,-3,3)=half*(ctor(:,1)-ctor(:,4)+ctor(:,2)-ctor(:,3))
   orb_lm(:,-3,4)=half*(ctor(:,1)-ctor(:,4)-ctor(:,2)+ctor(:,3))
  end if
  if(lmax>=2) then
!  s            px            py
   orb_lm(:,-4,1)=rs3*ctor(:,1)-rs6*ctor(:,4)+rs2*ctor(:,2)
   orb_lm(:,-4,2)=rs3*ctor(:,1)-rs6*ctor(:,4)-rs2*ctor(:,2)
   orb_lm(:,-4,3)=rs3*ctor(:,1)+2._dp*rs6*ctor(:,4)
!  pz           dz2
   orb_lm(:,-4,4)= rs2*ctor(:,3)+rs2*ctor(:,7)
   orb_lm(:,-4,5)=-rs2*ctor(:,3)+rs2*ctor(:,7)
!  s            px            dz2         dx2-y2
   orb_lm(:,-5,1)=rs6*ctor(:,1)-rs2*ctor(:,4)-rs12*ctor(:,7)+half*ctor(:,9)
   orb_lm(:,-5,2)=rs6*ctor(:,1)+rs2*ctor(:,4)-rs12*ctor(:,7)+half*ctor(:,9)
!  s            py            dz2         dx2-y2
   orb_lm(:,-5,3)=rs6*ctor(:,1)-rs2*ctor(:,2)-rs12*ctor(:,7)-half*ctor(:,9)
   orb_lm(:,-5,4)=rs6*ctor(:,1)+rs2*ctor(:,2)-rs12*ctor(:,7)-half*ctor(:,9)
!  s            pz           dz2
   orb_lm(:,-5,5)=rs6*ctor(:,1)-rs2*ctor(:,3)+rs3*ctor(:,7)
   orb_lm(:,-5,6)=rs6*ctor(:,1)+rs2*ctor(:,3)+rs3*ctor(:,7)
  end if

! stuff complex wannier orbital coefficient array
  do iwan=1,nwan
   wan_lm(:,iwan)=orb_lm(:,proj_l(iwan),proj_m(iwan))
  end do

! setup to rotate wan_lm to new axes if called for
! skip if only s projectors are used
  if ( lmax>0 ) then
!  generate a set of nr=lmax2 random vectors
!  idum=123456
   do ir=1,lmax2
    do ii=1,3
     r(ii,ir) = uniformrandom(idum)-0.5d0
    end do !ii
    call ylm_cmplx(lmax,ylmcp,r(1,ir),r(2,ir),r(3,ir))
    ylmc_rr(ir,:)=conjg(ylmcp(:))
    ylmc_rr_save(ir,:)=conjg(ylmcp(:))
   end do !ir

   ylmc_rrinv(:,:)=c0
   do ii=1,lmax2
    ylmc_rrinv(ii,ii)=c1
   end do !ii
!  calculate inverse of ylmc(ir,lm) matrix
   call ZGESV(lmax2,lmax2,ylmc_rr,lmax2,ipiv,ylmc_rrinv,lmax2,info)

!  check that r points are independent (ie., that matrix inversion wasn't
!  too close to singular)
   ylmc_rr=matmul(ylmc_rrinv,ylmc_rr_save)
   test=zero
   do ii=1,lmax2
    ylmc_rr(ii,ii)=ylmc_rr(ii,ii)-c1
    do jj=1,lmax2
     test=max(abs(ylmc_rr(ii,jj)),test)
    end do !ii
   end do !jj
   if(test>tol8) then
    write(message, '(8a)' )ch10,&
&    ' mlwfovlp_get_wan_lm : ERROR',ch10,&
&    '  matrix inversion error for wannier rotations',ch10,&
&    '  random vectors r(j,1:nr) are not all independent !! ',ch10,&
&    '  Action : re-seed uniformrandom or maybe just try again'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if !test>tol8

!  end of the preliminaries, now to the rotations of the wannier orbitals
   do iwan=1,nwan
!   don't bother for s orbitals
    if(proj_l(iwan)==0) cycle
!   check for default axes and cycle if found
    if(proj_z(1,iwan)==zero .and. proj_z(2,iwan)==zero .and.&
&    proj_z(3,iwan)== one .and. proj_x(1,iwan)==one .and.&
&    proj_x(2,iwan)==zero .and. proj_x(3,iwan)==zero) cycle

!   get the u matrix that rotates the reference frame
    call rotmat(proj_x(:,iwan),proj_z(:,iwan),inversion_flag,umat)

!   find rotated r-vectors. Optional inversion
!   operation is an extension of the wannier90 axis-setting options
!   which only allow for proper axis rotations
    if(inversion_flag==1) then
     rp(:,:)= -matmul ( umat(:,:),  r(:,:) )
    else
     rp(:,:) = matmul ( umat(:,:) , r(:,:) )
    end if !inversion_flag

    do ir=1,lmax2
!    get the ylm representation of the rotated vectors
     call ylm_cmplx(lmax,ylmcp,rp(1,ir),rp(2,ir),rp(3,ir))
     ylmc_rp(ir,:)=conjg(ylmcp(:))
    end do !ir
!   the matrix product sum(ir) ylmc_rrinv(lm,ir)*ylmc_rp(ir,lm') gives the
!   the complex lmXlm matrix representation of the coordinate rotation
    crot(:,:)=matmul(ylmc_rrinv(:,:),ylmc_rp(:,:))

!   now rotate the current wannier orbital
    ylmcp(:)=matmul(crot(:,:),wan_lm(:,iwan))
    wan_lm(:,iwan)=ylmcp(:)
   end do !iwan
  end if !lmax>0
 end subroutine get_wan_lm

end subroutine mlwfovlp_proj
!!***
