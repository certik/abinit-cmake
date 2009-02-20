!{\src2tex{textfont=tt}}
!!****f* ABINIT/overlap_wf
!! NAME
!! overlap_wf
!!
!! FUNCTION
!! returns overlap between energy eigen states and auxiliary wave function
!! (this file is a modified version of wffile.F90)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (JB)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! Needs an unformatted wave function from abinit.
!! exchn2n3d=if 1, n2 and n3 are exchanged
!! headform=format of the wf file
!!
!! natom = number of atoms in the unit cell
!! nbands = size of e_kpt
!! nr1,nr2,nr3 = grid size (nr1 x nr2 x nr3 = filrho dimension)
!! ntypat = number of atom type
!! ucvol = unit cell volume (> 0)
!! densfileformat = flag for the format of the density file:
!!       0 = ASCII
!!       1 = binary
!! denval = density value exported by interpol, to be wrote in the output file
!! filrho = name of the density file (ASCII or binary)
!! filtau = name of the atomic position file (Xmol format)
!! rprim = orientation of the unit cell axes
!! cbandpick = bandindex for the wf
!! ckpt = kpoint index for the wf
!! csppol = spin polarization
!!
!! OUTPUT
!!  coverlap = < auxiliary wf | energy eigen states >
!!  e_kpt = all the energy eigen values at the particular kpoint
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      localorb_S
!!
!! CHILDREN
!!      fourwf,getkpgnorm,hdr_skip,kpgio,metric,rwwf,sphereboundary
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine overlap_wf(cwave1,e_kpt,exchn2n3d,csppol,cbandpick,ckpt,&
     & ecut,headform,istwfk,kpt,nband,nbands,nkpt,npwarr,&
     & nr1,nr2,nr3,nspinor,nsppol,paral_kgb,rprimd,coverlap)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_13recipspace
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: cbandpick,ckpt,csppol,exchn2n3d,headform,nbands,nkpt,nr1
 integer,intent(in) :: nr2,nr3,nspinor,nsppol,paral_kgb
 real(dp),intent(in) :: ecut
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt),npwarr(nkpt)
 real(dp),intent(in) :: kpt(3,nkpt),rprimd(3,3)
 real(dp),intent(out) :: e_kpt(nbands)
!no_abirules
 complex(dp),intent(out) :: coverlap
 complex(dp),intent(in) :: cwave1(nr1,nr2,nr3)

!Local variables-------------------------------
       character(*), parameter :: inputfile='cut.in'
       character*(fnlen) :: ylmnam=" "
!scalars
 integer,save :: tim_fourwf=0,tim_rwwf=0
 integer :: cband,cgshift,ckpt1,ckpt2,cplex,cspinor,formeig,i,i1,i10,i11,i12,i2
 integer :: i3,i4,i5,i6,i7,i8,i9,ia,iband,iband1,iband2,icg,ichoice,ierr,ifile
 integer :: ii,ii1,ii2,ii3,ikpt,ilang,init_prefact,insmet,ioffkg,ios,iout
 integer :: iprompt,ipw,ir1,ir2,ir3,ispden,isppol,istop1,ivect,ix,ix1,ix2,ix3
 integer :: ixfh,ixint,iy,iy1,iy2,iy3,iz,iz1,iz2,iz3,j,j1,j2,k,k1,k2,l,m,m1
 integer :: mband,mbess,mcg,mgfft,mkmem,mlang,mpw,n,n4,n5,n6,nband_disk,nfit
 integer :: nkpts,nplwv,npw_k,npwout,nradint,nshift,nsize,nspden,nstart,nxfh
 integer :: oldcband,oldckpt,oldcspinor,oldcsppol,option,prtsphere,select_exit
 integer :: unkg,unylm=0
 real(dp) :: alpha,arg,bessargmax,bessint_delta,efermi,energy,hx,hy,hz,kpgmax
 real(dp) :: normtot,ratsph,re1,re2,re3,re4,re5,re6,re7,re8,rkptx,rkpty,rkptz
 real(dp) :: rlx,rly,rlz,rnorm,rnorm1,rnorm2,rx,ry,rz,tmpi,tmpr,tpi,ucvol
 real(dp) :: weight,x,xnow,y,ynow,z,znow
 character(len=10) :: string
 character(len=4) :: mode_paral
 character(len=500) :: message
 character(len=fnlen) :: kgnam,output,output1
 character :: outputchar
 type(mpi_type) :: mpi_enreg
 type(wffile_type) :: wff
!arrays
 integer :: ngfft(18)
 integer,allocatable :: gbound(:,:),iindex(:),kg(:,:),kg_dum(:,:),kg_k(:,:)
 integer,allocatable :: npwarr1(:),npwtot1(:)
 real(dp) :: acell(3),gmet(3,3),gprimd(3,3),oldkpt(3),rmet(3,3)
 real(dp),allocatable :: bess_fit(:,:,:),bess_spl(:,:),bess_spl_der(:,:)
 real(dp),allocatable :: cg(:,:),cgcband(:,:),cgcband1(:,:),denpot(:,:,:)
 real(dp),allocatable :: eigen1(:),fofgin(:,:),fofgout(:,:),fofr(:,:,:,:)
 real(dp),allocatable :: fofr1(:,:,:,:),kpgnorm(:),occ1(:),rint(:)
 real(dp),allocatable :: spl_bessint(:,:),sum_1atom_1ll(:,:),x_bess(:)
 real(dp),allocatable :: xfhist(:,:,:,:),ylm(:,:),ylm_k(:,:)
 character(len=fnlen),allocatable :: filename(:)
!no_abirules
 complex(dp) :: cmp1,cmp2,cmp3,cmp4,cmpwan1,cmpwan2,cmpwan3,csum,csum1,csum2,csum3

! *************************************************************************

 cspinor = 1

!begin executable section
!ios = 0
 mpi_enreg%paralbd=0

 formeig=0
 oldckpt=0
 oldcband=0
 oldcsppol=0
 oldcspinor=0

 iout=-1
 call metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

 allocate (kg_dum(3,0))

!#############################################################################

 tpi = 8.0*atan(1.0)

 mband=maxval(nband)
 mpw=maxval(npwarr)
 mcg=mpw*nspinor*mband

 allocate(cg(2,mcg),eigen1((2*mband)**formeig*mband),&
 occ1(mband))

!==========================================================================
!necessary procedures for the fft subroutine:
!==========================================================================
 mpi_enreg%paralbd=0
 mpi_enreg%nproc_fft=1
 mpi_enreg%me=0
 mpi_enreg%me_fft=0
 mpi_enreg%paral_fft=0
 mpi_enreg%paral_compil_kpt=0
 mpi_enreg%paral_compil_fft=0
 mpi_enreg%fft_option_lob=1

 ngfft(1)=nr1
 ngfft(2)=nr2
 ngfft(3)=nr3

 if (mod(nr1,2)==0)then
  ngfft(4)=nr1+1
 else
  ngfft(4)=nr1
 end if
 if(mod(nr2,2)==0)then
  ngfft(5)=nr2+1
 else
  ngfft(5)=nr2
 end if
 ngfft(6)=nr3
 ngfft(7)=111
 ngfft(8)=256

 mode_paral='pers'
 mkmem=nkpt
 mgfft=maxval(ngfft(1:3))

 allocate(npwarr1(nkpt),kg(3,mpw*mkmem),npwtot1(nkpt))
!create positions index for pw
 call kpgio(ecut,exchn2n3d,gmet,istwfk,kg,kgnam,kpt,mkmem,nband,nkpt,&
 mode_paral,mpi_enreg,mpw,npwarr1,npwtot1,nsppol,unkg)

!additional allocation:
 n4=       ngfft(4)
 n5=       ngfft(5)
 n6=       ngfft(6)

 e_kpt = 0.0
!-------------------------------------------------------------------
!reading wave function from "_wfk" file :
!-------------------------------------------------------------------
 cband = cbandpick
 cg = 0.0
 eigen1 = 0.0
 occ1 = 0.0
 rewind(19)

 wff%unwff=19
 wff%accesswff=0
 call hdr_skip(wff,ierr)

!call hdr_skip(19)

 do isppol=1,csppol
! write(*,*)'ckpt',ckpt
  do ikpt=1,nkpt

   if(isppol==csppol .and. ikpt==ckpt)then
    option=1
!   write(*,*)'option',option,cband
   else
    option=-1
   end if
   call rwwf(cg,eigen1,formeig,headform,0,ikpt,isppol,kg_dum,&
&   mband,mcg,mpi_enreg,nband(ikpt),nband_disk,&
&   npwarr(ikpt),nspinor,occ1,option,0,tim_rwwf,wff)

   if(isppol==csppol.and.ikpt==ckpt)then
    re1 = 0.0
    re2 = 0.0
    do ii1 = 1,mband
     energy=eigen1(ii1)
     e_kpt(ii1) = energy
    end do
   end if

   if(option==1)exit       ! when the target wf has been read,
!  exit the wf file reading
  end do
  if(option==1)exit
 end do

 ioffkg=0
 do ikpt=1,ckpt-1
  ioffkg=ioffkg+npwarr1(ikpt)
 end do
 npw_k=npwarr(ckpt)

 allocate(gbound(2*mgfft+8,2),kg_k(3,npw_k))
 allocate(kpgnorm(npw_k))

 kg_k(:,1:npw_k)=kg(:,1+ioffkg:npw_k+ioffkg)
 call getkpgnorm(gprimd,kpt(:,ckpt),kg_k,kpgnorm,npw_k)
 call sphereboundary(gbound,istwfk(ckpt),kg_k,mgfft,npw_k)
 n4=ngfft(4)
 n5=ngfft(5)
 n6=ngfft(6)
!cplex=0
 cplex=1
 cgshift=(cband-1)*npw_k*nspinor + (cspinor-1)*npw_k

 allocate(cgcband(2,npw_k))
 allocate(cgcband1(2,npw_k))
 allocate(denpot(cplex*n4,n5,n6))
 allocate(fofgout(2,npw_k))
 allocate(fofr(2,n4,n5,n6))
 allocate(fofr1(2,n4,n5,n6))

 cgcband(:,1:npw_k)=cg(:,cgshift+1:cgshift+npw_k)

 call fourwf(cplex,denpot,cgcband,fofgout,fofr,&
& gbound,gbound,&
& istwfk(ckpt),kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,&
& npw_k,n4,n5,n6,0,paral_kgb,tim_fourwf,weight)

!Swaping for FFT :
 do k=1,nr3
  if(mod(real(nr3),2.0).eq.0.0)then
   if(k < ((nr3/2)+1))n=k+(nr3/2)
   if(k > (nr3/2))n=k-(nr3/2)
  else
   if(k < (nr3-1)/2+1)n=k+(nr3+1)/2
   if(k > (nr3-1)/2)n=k-(nr3-1)/2
  end if
  do j=1,nr2
   if(mod(real(nr2),2.0).eq.0.0)then
    if(j < ((nr2/2)+1))m1=j+(nr2/2)
    if(j > (nr2/2))m1=j-(nr2/2)
   else
    if(j < (nr2-1)/2+1)m1=j+(nr2+1)/2
    if(j > (nr2-1)/2)m1=j-(nr2-1)/2
   end if
   do i=1,nr1
    if(mod(real(nr1),2.0).eq.0.0)then
     if(i < ((nr1/2)+1))l=i+(nr1/2)
     if(i > (nr1/2))l=i-(nr1/2)
    else
     if(i < (nr1-1)/2+1)l=i+(nr1+1)/2
     if(i > (nr1-1)/2)l=i-(nr1-1)/2
    end if

    fofr1(1,l,m1,n) = real(cwave1(i,j,k))
    fofr1(2,l,m1,n) = aimag(cwave1(i,j,k))

   end do
  end do
 end do

 call fourwf(cplex,denpot,cgcband1,fofgout,fofr1,&
& gbound,gbound,&
& istwfk(ckpt),kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,&
& npw_k,n4,n5,n6,3,paral_kgb,tim_fourwf,weight)

 coverlap = 0.0
 do i = 1,npw_k
  re1 = fofgout(1,i)
  re2 = fofgout(2,i)
  re3 = cgcband(1,i)
  re4 = cgcband(2,i)
  coverlap = coverlap + cmplx(re1*re3 + re2*re4, re1*re4 - re3*re2)
 end do

 deallocate(cgcband)
 deallocate(cgcband1)
 deallocate(denpot)
 deallocate(fofgout,fofr)
 deallocate(fofr1)

 deallocate(gbound,kg_k)
 deallocate(kpgnorm)

 return

end subroutine overlap_wf
!!***
