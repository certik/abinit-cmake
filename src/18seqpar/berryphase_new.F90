!{\src2tex{textfont=tt}}
!!****f* ABINIT/berryphase_new
!! NAME
!! berryphase_new
!!
!! FUNCTION
!! This routine computes the Berry Phase polarization
!!  and the finite difference expression of the ddk.
!!  [see for example Na Sai et al., PRB 66, 104108 (2002)]
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT  group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions
!! cprj(natom,nspinor*mband*mkmem*nsppol*usecrpj)=<p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                                and each |p_lmn> non-local projector
!! dtfil <type(datafiles_type)>=variables related to files
!! dtset <type(dataset_type)>=all input variables in this dataset
!! gprimd(3,3)=reciprocal space dimensional primitive translations
!! hdr <type(hdr_type)>=the header of wf, den and pot files
!! kg(3,mpw*mkmem)=reduced planewave coordinates
!! mband=maximum number of bands
!! mkmem=number of k points which can fit in memory; set to 0 if use disk
!! mpi_enreg=informations about MPI parallelization
!! mpw=maximum dimensioned size of npw
!! natom=number of atoms in cell
!! nattyp(ntypat)= # atoms of each type.
!! nkpt=number of k points
!! npwarr(nkpt)=number of planewaves in basis at this k point
!! nspinor=number of spinorial components of the wavefunctions
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! ntypat=number of types of atoms in unit cell
!! nkpt=number of k-points
!! option = 1: compute Berryphase polarization
!!          2: compute finite difference expression of the ddk
!!          3: compute polarization & ddk
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! pelev(3)= expectation value polarization term (PAW only) in cartesian coordinates
!! pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!! pwind_alloc = first dimension of pwind
!! pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! typat(natom)=type integer for each atom in cell
!! ucvol=unit cell volume in bohr**3.
!! unit_out= unit for output of the results (usually the .out file of ABINIT)
!!   The option unit_out = 0 is allowed. In this case, no information is written
!!   to the output file but only to the log file.
!! usecprj=1 if cprj datastructure has been allocated
!! usepaw= 1: use paw framework. 0:do not use paw.
!! wffnow=struct info for wf disk file
!! xred(3,natom)=reduced atomic coordinates
!! zion(ntypat)=valence charge of each type of atom
!!
!! OUTPUT
!! pel(3) = reduced coordinates of the electronic polarization (a. u.)
!! pion(3)= reduced coordinates of the ionic polarization (a. u.)
!!
!! SIDE EFFECTS
!! Input/Output
!! dtefield <type(efield_type)> = variables related to Berry phase
!!       and electric field calculations (see initberry.f).
!!       In case berryopt = 4, the overlap matrices computed
!!       in this routine are stored in dtefield%smat in order
!!       to be used in the electric field calculation.
!!
!! TODO
!!  - Use the analytical relation between the overlap matrices
!!    S(k,k+dk) and S(k+dk,k) to avoid to recompute them
!!    when ifor = 2.
!!
!! NOTES
!! - pel and pion do not take into account the factor 1/ucvol
!! - In case of a ddk calculation, the eigenvalues are not computed.
!! - The ddk computed by this routine should not be used to
!!   compute the electronic dielectric tensor.
!!
!! PARENTS
!!      elpolariz,scfcv
!!
!! CHILDREN
!!      appdig,leave_new,mpi_recv,mpi_send,outwf,polcart,rhophi,smatrix
!!      smatrix_paw,smatrix_pawinit,wrtout,xcomm_world,xme_init,xredxcart
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine berryphase_new(cg,cprj,dtefield,dtfil,dtset,&
&  gprimd,hdr,kg,mband,&
&  mkmem,mpi_enreg,mpw,natom,nattyp,npwarr,nspinor,nsppol,ntypat,&
&  nkpt,option,pawang,pawrad,pawtab,pel,pelev,pion,pwind,pwind_alloc,pwnsfac,&
&  rprimd,typat,ucvol,unit_out,usecprj,usepaw,wffnow,xred,zion)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes
#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13paw
 use interfaces_15common
 use interfaces_16response
 use interfaces_18seqpar, except_this_one => berryphase_new
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer, intent(in) :: mband,mkmem,mpw,natom,nkpt,nspinor,nsppol,ntypat,option
 integer, intent(in) :: pwind_alloc,unit_out,usecprj,usepaw
 real(dp), intent(in) :: ucvol
 type(MPI_type), intent(inout) :: mpi_enreg
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(in) :: dtset
 type(efield_type), intent(inout) :: dtefield
 type(hdr_type), intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang
 type(wffile_type), intent(inout) :: wffnow
!arrays
 integer, intent(in) :: kg(3,mpw*mkmem)
 integer, intent(in) :: nattyp(ntypat),npwarr(nkpt),pwind(pwind_alloc,2,3)
 integer, intent(in) :: typat(natom)
 real(dp), intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),gprimd(3,3)
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
 real(dp), intent(in) :: rprimd(3,3),zion(ntypat)
 real(dp), intent(inout) :: xred(3,natom)
 real(dp), intent(out) :: pel(3),pelev(3),pion(3)
 type(pawrad_type), intent(in) :: pawrad(ntypat*usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)
 type(cprj_type),intent(in) ::  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)

!Local variables -------------------------
 integer :: count,count1,ddkflag,dest,formeig
 integer :: iatom,icg,icg1,idir,idirrot,idum,idum1,ikpt1i_sp,iikpt
 integer :: ierr,ifor,iforrot,ii,ikg,ikpt,ikpti,ikpt1,ikpt2,ikpt_loc
 integer :: inibz,ikpt1i,ikpt2i,idum2
 integer :: isppol,istr,itrs,itypat,jkpt,jkpti,jkstr,job,jsppol,maxbd,mcg,mcg1_k
 integer :: minbd,mxfh,nband_k,nfor,npw_k1,npw_k2,nqpt,nstep,nxfh,orig,pertcase
 integer :: response,shiftbd,source,spaceComm,tag
 integer :: g1(3)
 real(dp) :: det_mod,dkinv,dphase,dtm_real,dtm_imag,fac,gmod,phase0
 real(dp) :: pol,polbtot,polion,politot,poltot,rho
 character(len=fnlen) :: fiwf1o,wff2nm
 character(len=500) :: message
 type(wvl_wf_type) :: wfs
 integer,allocatable :: pwind_k(:),sflag_k(:)
 real(dp) :: det_average(2),dk(3),dtm_k(2),gpard(3),pel_cart(3),pion_cart(3)
 real(dp) :: polb(nsppol),ptot_cart(3),rel_string(2),xcart(3,natom)
 real(dp),allocatable :: buffer(:,:),buffer1(:),buffer2(:)
 real(dp),allocatable :: cg1(:,:),cg1_k(:,:),cgq(:,:)
 real(dp),allocatable :: smat_k_paw(:,:,:)
 real(dp),allocatable :: det_string(:,:),dtm(:,:),eig_dum(:),occ_dum(:)
 real(dp),allocatable :: polberry(:),resid(:),pwnsfac_k(:,:)
 real(dp),allocatable :: smat_inv(:,:,:),smat_k(:,:,:)
 real(dp),allocatable :: xfhist(:,:,:,:)

!BEGIN TF_CHANGES
 integer :: me
!END TF_CHANGES

!no_abirules
#if defined MPI
             integer :: mpi_status(MPI_STATUS_SIZE)
!BEGIN TF_CHANGES
             call xcomm_world(mpi_enreg,spaceComm)
!END TF_CHANGES
#endif

! ***********************************************************************

!DEBUG
!write(*,*)' berryphase_new : enter '
!do ii=1,pwind_alloc
! write(6,*)ii,pwnsfac(:,ii)
!end do
!stop
!ENDDEBUG

!BEGIN TF_CHANGES
!Define me
 call xme_init(mpi_enreg,me)
!END TF_CHANGES

 allocate(pwind_k(mpw),pwnsfac_k(4,mpw),sflag_k(dtefield%nband_occ))
 pwind_k(:) = 0
 pwnsfac_k(1,:) = 1.0_dp ! bra real
 pwnsfac_k(2,:) = 0.0_dp ! bra imag
 pwnsfac_k(3,:) = 1.0_dp ! ket real
 pwnsfac_k(4,:) = 0.0_dp ! ket imag

 if (nspinor == 2) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' berryphase : ERROR -',ch10,&
&  '  This routine does not yet work for nspinor = 2'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 if (maxval(dtset%istwfk(:)) /= 1) then
  write(message, '(a,a,a,a,a,a)' )ch10,&
&  ' berryphase : BUG -',ch10,&
&  '  This routine does not work yet with istwfk /= 1.',ch10,&
&  '  This should have been tested previously ...'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 if (usecprj/=usepaw) then
  write(message, '(4a)' )ch10,&
&  ' berryphase : BUG -',ch10,&
&  '  cprj datastructure has no been allocated !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 mcg = mpw*nspinor*mband*nsppol*mkmem
 mcg1_k = mpw*nspinor*mband
 shiftbd = 1
 if (option > 1) then
  allocate(cg1(2,mpw*nspinor*mband*mkmem*nsppol))
  allocate(eig_dum(2*mband*mband*nkpt*nsppol),occ_dum(mband*nkpt*nsppol))
  eig_dum(:) = zero
  occ_dum(:) = dtefield%sdeg
 end if

 allocate(dtm(2,dtefield%fnkpt*nsppol))
 allocate(cg1_k(2,mcg1_k))

 pel(:) = zero     ; pion(:) = zero

 minbd = 1   ;  maxbd = dtefield%nband_occ

 do idir = 1, 3

  dtm(:,:) = zero

  if (dtset%rfdir(idir) == 1) then

   if (abs(dtefield%efield_dot(idir)) < tol12) dtefield%sflag(:,:,:,idir) = 0

!  Check whether the polarization or the ddk must be computed

!  nfor = 1 : to compute P, I only need the WF at k + dk
!  nfor = 2 : to compute the ddk, I need the WF at k + dk and k - dk
!  dkinv    : +-1/2dk

  if (option > 1) then

   ddkflag = 1
   nfor = 2
   job = 1
   cg1(:,:) = zero
   if (option == 3) job = 11

  else if (option == 1) then

   ddkflag = 0
   nfor = 1
   job = 10

  end if

  dk(:) = dtefield%dkvecs(:,idir)
  gpard(:) = dk(1)*gprimd(:,1) + dk(2)*gprimd(:,2) + dk(3)*gprimd(:,3)
  gmod = sqrt(dot_product(gpard,gpard))
  if (option > 1) dkinv = one/(two*dk(idir))

  write(message,'(a,a,a,3f9.5,a,a,3f9.5,a)')ch10,&
&  ' Computing the polarization (Berry phase) for reciprocal vector:',ch10,&
&  dk(:),' (in reduced coordinates)',ch10,&
&  gpard(1:3),' (in cartesian coordinates - atomic units)'
  call wrtout(6,message,'COLL')
  if (unit_out /= 0) then
   call wrtout(unit_out,message,'COLL')
  end if

  write(message,'(a,i5,a,a,i5)')&
&  ' Number of strings: ',dtefield%nstr(idir),ch10,&
&  ' Number of k points in string:', dtefield%nkstr(idir)
  call wrtout(6,message,'COLL')
  if (unit_out /= 0) then
   call wrtout(unit_out,message,'COLL')
  end if

  if ((option == 2).or.(option == 3)) then

   write(message,'(a,a,a,3f9.5,a,a,3f9.5,a)')ch10,&
&   ' Computing the ddk (Berry phase) for reciprocal vector:',ch10,&
&   dk(:),' (in reduced coordinates)',ch10,&
&   gpard(1:3),' (in cartesian coordinates - atomic units)'
   call wrtout(6,message,'COLL')
   if (unit_out /= 0) then
    call wrtout(unit_out,message,'COLL')
   end if

  end if

  do ifor = 1, nfor

    if (ifor == 2) then
     dk(:) = -1_dp*dk(:)
     job = 1   ! only the inverse of the overlap matrix is required
     dkinv = -1_dp*dkinv
    end if


!   Compute the determinant and/or the inverse of the overlap matrix
!   for each pair of k-points < u_nk | u_nk+dk >

    icg = 0 ; icg1 = 0
    allocate(smat_k(2,dtefield%nband_occ,dtefield%nband_occ))
    allocate(smat_inv(2,dtefield%nband_occ,dtefield%nband_occ))

    ikpt_loc = 0

    do isppol = 1, nsppol

     ikpt1 = 0
     if (mpi_enreg%paral_compil_kpt == 0) ikpt_loc = 0
     do while (ikpt_loc < dtefield%fmkmem_max)

      if (ikpt_loc < dtefield%fmkmem) ikpt1 = ikpt1 + 1
      if ((ikpt1 > dtefield%fnkpt).and.(ikpt_loc < dtefield%fmkmem)) exit
      ikpt1i = dtefield%indkk_f2ibz(ikpt1,1)
      nband_k = dtset%nband(ikpt1i + (isppol-1)*dtset%nkpt)

!DEBUG
!     Please keep that debuggin feature
!     write(6, '(a,4i4)' )' berryphase_new : ikpt1,isppol,idir,ifor=',ikpt1,isppol,idir,ifor
!ENDDEBUG

#if defined MPI
!BEGIN TF_CHANGES
           if ((minval(abs(mpi_enreg%proc_distrb(ikpt1i,1:nband_k,isppol) -&
&               me)) /= 0).and.(ikpt_loc <= dtefield%fmkmem)) then
!END TF_CHANGES
            cycle
           end if
#endif

      ikpt_loc = ikpt_loc + 1

      inibz=0
      if (dtset%kptns(1,ikpt1i) == dtefield%fkptns(1,ikpt1) .and. &
&         dtset%kptns(2,ikpt1i) == dtefield%fkptns(2,ikpt1) .and. &
&         dtset%kptns(3,ikpt1i) == dtefield%fkptns(3,ikpt1)) inibz=1

      ikg = dtefield%fkgindex(ikpt1)
      ikpt2 = dtefield%ikpt_dk(ikpt1,ifor,idir)
      ikpt2i = dtefield%indkk_f2ibz(ikpt2,1)
      itrs = 0
      if (dtefield%indkk_f2ibz(ikpt1,6) == 1 ) itrs = itrs + 1
      if (dtefield%indkk_f2ibz(ikpt2,6) == 1 ) itrs = itrs + 10

      npw_k1 = npwarr(ikpt1i)
      npw_k2 = npwarr(ikpt2i)

!     ji: the loop is over the FBZ, but sflag and smat only apply to the IBZ
      if (dtset%berryopt == 4 .and. inibz == 1) then
       ikpt1i_sp=ikpt1i+(isppol-1)*dtset%nkpt
       smat_k(:,:,:) = dtefield%smat(:,:,:,ikpt1i_sp,ifor,idir)
      else
       smat_k(:,:,:) = zero
      end if

      pwind_k(1:npw_k1) = pwind(ikg+1:ikg+npw_k1,ifor,idir)
      pwnsfac_k(1,1:npw_k1) = pwnsfac(1,ikg+1:ikg+npw_k1)
      pwnsfac_k(2,1:npw_k1) = pwnsfac(2,ikg+1:ikg+npw_k1)

!DEBUG
!     write(6,*)' berryphase_new : dtset%berryopt,inibz,ikpt1i,isppol,dtset%nkpt,ifor,idir',dtset%berryopt,inibz,ikpt1i,isppol,dtset%nkpt,ifor,idir
!     write(6, '(a,4i4)' )' berryphase_new : sflag_k(:)=',sflag_k(:)
!ENDDEBUG

      if (dtset%berryopt == 4 .and. inibz == 1) then
       ikpt1i_sp=ikpt1i+(isppol-1)*dtset%nkpt
       sflag_k(:) = dtefield%sflag(:,ikpt1i_sp,ifor,idir)
      else
       sflag_k(:) = 0
      end if

!DEBUG
!     write(6, '(a,4i4)' )' berryphase_new : sflag_k(:)=',sflag_k(:)
!ENDDEBUG

!DEBUG
!write(*,'(a,7i4)')'me, idir,ifor, ikpt_loc, ikpt1, isppol = ',&
!& mpi_enreg%me,idir,ifor,ikpt_loc,ikpt1,isppol
!write(*,'(a,10i3)')'pwind_k(1:10) = ',pwind_k(1:10)
!ENDDEBUG

#if defined MPI

           count = npw_k2*nspinor*nband_k
           allocate(cgq(2,count))
           source = mpi_enreg%proc_distrb(ikpt2i,1,isppol)

           do jsppol = 1, nsppol
            do jkpt = 1, dtefield%fnkpt

             jkpti = dtefield%indkk_f2ibz(jkpt,1)

!BEGIN TF_CHANGES
             if ((jkpt == ikpt2).and.(source /= me).and.&
&               (ikpt_loc <= dtefield%fmkmem).and.(jsppol == isppol)) then
!END TF_CHANGES

              allocate(buffer(2,npw_k2))
              tag = jkpt + (jsppol - 1)*dtefield%fnkpt

              call MPI_RECV(buffer,2*npw_k2,MPI_DOUBLE_PRECISION,&
&                            source,tag,spaceComm,mpi_status,ierr)
              pwnsfac_k(3,1:npw_k2) = buffer(1,1:npw_k2)
              pwnsfac_k(4,1:npw_k2) = buffer(2,1:npw_k2)
              deallocate(buffer)

             end if

             do dest = 1, mpi_enreg%nproc

              iikpt=ikpt_loc+dtefield%fmkmem_max*nsppol
!BEGIN TF_CHANGES
              if ((minval(abs(mpi_enreg%proc_distrb(jkpti,1:nband_k,jsppol)&
&                 - me)) == 0).and.&
&                 (mpi_enreg%kptdstrb(dest,ifor+2*(idir-1),iikpt)&
&                 == jkpt + (jsppol - 1)*dtefield%fnkpt)) then
               if (((dest-1) /= me)) then
!END TF_CHANGES

                tag = jkpt + (jsppol - 1)*dtefield%fnkpt
                count1 = npwarr(jkpti)
                allocate(buffer(2,count1))
                idum = dtefield%fkgindex(jkpt)
                buffer(1,1:count1)  = pwnsfac(1,idum+1:idum+count1)
                buffer(2,1:count1)  = pwnsfac(2,idum+1:idum+count1)

                call MPI_SEND(buffer,2*count1,MPI_DOUBLE_PRECISION,&
&                              (dest-1),tag,spaceComm,mpi_status,ierr)
                deallocate(buffer)

               else

                idum = dtefield%fkgindex(jkpt)
                pwnsfac_k(3,1:npw_k2) = pwnsfac(1,idum+1:idum+npw_k2)
                pwnsfac_k(4,1:npw_k2) = pwnsfac(2,idum+1:idum+npw_k2)

               end if

              end if

             end do   ! loop over dest

            end do   ! jkpt
           end do    ! jsppol

           do jsppol = 1, nsppol

            do jkpt = 1, nkpt
!BEGIN TF_CHANGES
             if ((jkpt == ikpt2i).and.(source /= me).and.&
&                (ikpt_loc <= dtefield%fmkmem).and.(jsppol == isppol)) then
!END TF_CHANGES

              tag = jkpt + (jsppol - 1)*nkpt
              call MPI_RECV(cgq,2*count,MPI_DOUBLE_PRECISION,&
&                           source,tag,spaceComm,mpi_status,ierr)
             end if

!----------------------------------------------------------------------------
!--------------- Here: send the WF to ALL the cpus that need it -------------
!----------------------------------------------------------------------------

             do dest = 1, mpi_enreg%nproc
!BEGIN TF_CHANGES
              if ((minval(abs(mpi_enreg%proc_distrb(jkpt,1:mband,jsppol)&
&                 - me)) == 0).and.&
&                 (mpi_enreg%kptdstrb(dest,ifor+2*(idir-1),ikpt_loc) == &
&                 jkpt + (jsppol - 1)*nkpt)) then

               icg1 = dtefield%cgindex(jkpt,jsppol)

               if (((dest-1) /= me)) then
!END TF_CHANGES

                tag = jkpt + (jsppol - 1)*nkpt
                count1 = npwarr(jkpt)*nband_k*nspinor
                allocate(buffer(2,count1))
                buffer(:,1:count1)  = cg(:,icg1+1:icg1+count1)

                call MPI_SEND(buffer,2*count1,MPI_DOUBLE_PRECISION,&
&                             (dest-1),tag,spaceComm,mpi_status,ierr)
                deallocate(buffer)

               else

                cgq(:,1:count)  = cg(:,icg1+1:icg1+count)

               end if

              end if

             end do   ! loop over dest

            end do   ! loop over jkpt

           end do   ! loop over jsppol

           if (ikpt_loc > dtefield%fmkmem) then
            deallocate(cgq)
            cycle
           end if

           icg1 = 0
           icg = dtefield%cgindex(ikpt1i,isppol)
           call smatrix(cg,cgq,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,job,maxbd,&
&                       mband,mcg,count,mcg1_k,minbd,&
&                       mpw,dtefield%nband_occ,&
&                       npw_k1,npw_k2,nspinor,nsppol,pwind_k,pwnsfac_k,sflag_k,&
&                       shiftbd,smat_inv,smat_k)

           deallocate(cgq)

#else

           icg = dtefield%cgindex(ikpt1i,isppol)
           icg1 = dtefield%cgindex(ikpt2i,isppol)

           idum = dtefield%fkgindex(ikpt2)
           pwnsfac_k(3,1:npw_k2) = pwnsfac(1,idum+1:idum+npw_k2)
           pwnsfac_k(4,1:npw_k2) = pwnsfac(2,idum+1:idum+npw_k2)

           if(usepaw.ne.1) then
           call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,job,maxbd,&
&                       mband,mcg,mcg,mcg1_k,minbd,&
&                       mpw,dtefield%nband_occ,&
&                       npw_k1,npw_k2,nspinor,nsppol,pwind_k,pwnsfac_k,sflag_k,&
&                       shiftbd,smat_inv,smat_k)
           else
           allocate(smat_k_paw(2,dtefield%nband_occ,dtefield%nband_occ))
!     ikpt1 could not have the correct ordering, choose ikpt1i in this
!     case
           g1=0.d0
           call smatrix_pawinit(smat_k_paw,cprj,ikpt1,ikpt2,g1,gprimd,&
&                                dtset%kpt,dtefield%nband_occ,dtefield%nband_occ,&
&                               mkmem,natom,nattyp,dtset%nband,&
&                               nkpt,nspinor,nsppol,ntypat,pawang,pawrad,pawtab,rprimd,&
&                               typat,usepaw,xred)
           smat_k=four_pi*smat_k_paw

           call smatrix_paw(cg,cg,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,job,maxbd,&
&                       mband,mcg,mcg,mcg1_k,minbd,&
&                       mpw,dtefield%nband_occ,&
&                       npw_k1,npw_k2,nspinor,nsppol,pwind_k,pwnsfac_k,sflag_k,&
&                       shiftbd,smat_inv,smat_k)

           deallocate(smat_k_paw)
           endif

#endif

      if ((job == 10).or.(job == 11)) then

       if (sqrt(dtm_k(1)*dtm_k(1) + dtm_k(2)*dtm_k(2)) < tol12) then
        write(message,'(a,a,a,a,i5,a,a,a)')ch10,&
&        ' berryphase_new : BUG - ',ch10,&
&        '  For k-point #',ikpt1,',',ch10,&
&        '  the determinant of the overlap matrix is found to be 0.'
        call wrtout(06,message,'PERS')
        call leave_new('PERS')
       end if

       dtm(1,ikpt1+(isppol-1)*dtefield%fnkpt) = dtm_k(1)
       dtm(2,ikpt1+(isppol-1)*dtefield%fnkpt) = dtm_k(2)
      end if

      if (dtset%berryopt == 4 .and. inibz == 1) then
       ikpt1i_sp=ikpt1i+(isppol-1)*dtset%nkpt
       dtefield%smat(:,:,:,ikpt1i_sp,ifor,idir) = &
&         smat_k(:,:,:)
       dtefield%sflag(:,ikpt1i_sp,ifor,idir) = &
&         sflag_k(:)
      end if

      if (((job == 1).or.(job == 11)).and.(inibz == 1)) then
       cg1(:,icg + 1: icg + npw_k1*nband_k*nspinor) = &
          cg1(:,icg + 1:icg + npw_k1*nband_k*nspinor) + &
          dkinv*cg1_k(:,1:npw_k1*nband_k*nspinor)
      end if

     end do ! close loop over k-points

    end do  ! close loop over nsppol

    deallocate(smat_inv,smat_k)

   end do   ! close loop over ifor

#if defined MPI
           allocate(buffer1(2*dtefield%fnkpt*nsppol))
           allocate(buffer2(2*dtefield%fnkpt*nsppol))
           count = 2*dtefield%fnkpt*nsppol
           buffer1(:) = reshape(dtm(:,:),(/count/))
!BEGIN TF_CHANGES
           call xsum_mpi(buffer1,buffer2,count,mpi_enreg%spaceComm,ierr)
!END TF_CHANGES
           dtm(:,:) = reshape(buffer2(:),(/2,dtefield%fnkpt*nsppol/))
           deallocate(buffer1,buffer2)
#endif

!===========================================================================

!  Compute the Berry phase polarization

   if ((option == 1).or.(option == 3)) then

!   Compute the electronic Berry phase

    write(message,'(a,a)')ch10,&
&    ' Compute the electronic contribution to polarization'
    call wrtout(6,message,'COLL')

    allocate(det_string(2,dtefield%nstr(idir)))
    allocate(polberry(dtefield%nstr(idir)))
    write(message,'(a,10x,a,10x,a)')ch10,&
&    'istr','polberry(istr)'
    call wrtout(6,message,'COLL')

    polbtot = zero
    polb(:) = zero
    do isppol = 1, nsppol

     det_string(1,:) = one ; det_string(2,:) = zero
     dtm_k(:) = one
     det_average(:) = zero

     do istr = 1, dtefield%nstr(idir)
      do jkstr = 1, dtefield%nkstr(idir)

       ikpt=dtefield%idxkstr(jkstr,istr,idir)

       dtm_real=dtm(1,ikpt+(isppol-1)*dtefield%fnkpt)
       dtm_imag=dtm(2,ikpt+(isppol-1)*dtefield%fnkpt)

       dtm_k(1) = det_string(1,istr)*dtm_real - &
&         det_string(2,istr)*dtm_imag
       dtm_k(2) = det_string(1,istr)*dtm_imag + &
&         det_string(2,istr)*dtm_real
       det_string(1:2,istr) = dtm_k(1:2)

      end do

      det_average(:) = det_average(:) + &
&        det_string(:,istr)/dble(dtefield%nstr(idir))

     end do

!    First berry phase that corresponds to det_average
!     phase0 = atan2(det_average(2),det_average(1))
     call rhophi(det_average,phase0,rho)
     det_mod = det_average(1)**2+det_average(2)**2

!    Then berry phase that corresponds to each string relative to the average
     do istr = 1, dtefield%nstr(idir)

      rel_string(1) = (det_string(1,istr)*det_average(1) + &
          det_string(2,istr)*det_average(2))/det_mod
      rel_string(2) = (det_string(2,istr)*det_average(1) - &
          det_string(1,istr)*det_average(2))/det_mod
!     dphase = atan2(rel_string(2),rel_string(1))
      call rhophi(rel_string,dphase,rho)
      polberry(istr) = dtefield%sdeg*(phase0 + dphase)/two_pi
      polb(isppol) = polb(isppol) + polberry(istr)/dtefield%nstr(idir)

      write(message,'(10x,i4,7x,e15.9)')istr,polberry(istr)
      call wrtout(6,message,'COLL')

     end do

     write(message,'(9x,a,7x,e15.9,1x,a,i4,a,a)')&
&     'total',polb(isppol),'(isppol=',isppol,')',ch10
     call wrtout(6,message,'COLL')

     polbtot = polbtot + polb(isppol)

    end do    ! isppol

!   Fold into interval [-1,1]
    polbtot = polbtot - 2_dp*nint(polbtot/2_dp)
    pel(idir) = polbtot

    deallocate(det_string,polberry)

!==========================================================================

!   Compute the ionic Berry phase

    call xredxcart(natom,1,rprimd,xcart,xred)
    politot = zero
    write(message,'(a)')' Compute the ionic contributions'
    call wrtout(6,message,'COLL')

    write(message,'(a,2x,a,2x,a,15x,a)')ch10,&
&    'itom', 'itypat', 'polion'
    call wrtout(6,message,'COLL')

    do iatom = 1, natom
     itypat = typat(iatom)

!    polion = zion(itypat)*dtefield%nkstr(idir)*&
!&            dot_product(xcart(1:3,iatom),gpard(1:3))

!    The ionic phase can be computed much easier
     polion = zion(itypat)*xred(idir,iatom)

!    Fold into interval (-1,1)
     polion = polion - 2_dp*nint(polion/2_dp)
     politot = politot + polion
     write(message,'(2x,i2,5x,i2,10x,e15.9)') iatom,itypat,polion
     call wrtout(6,message,'COLL')
    end do

!   Fold into interval [-1,1] again
    politot = politot - 2_dp*nint(politot/2_dp)
    pion(idir) = politot

    write(message,'(9x,a,7x,es19.9)') 'total',politot
    call wrtout(6,message,'COLL')


!==========================================================================

!   Compute the total polarization

    poltot = politot + polbtot

    write(message,'(a,a)')ch10,&
&    ' Summary of the results'
    call wrtout(6,message,'COLL')
    if (unit_out /= 0) then
     call wrtout(unit_out,message,'COLL')
    end if

    write(message,'(a,es19.9)')&
&    ' Electronic Berry phase ' ,polbtot
    call wrtout(6,message,'COLL')
    if (unit_out /= 0) then
     call wrtout(unit_out,message,'COLL')
    end if

    write(message,'(a,es19.9)') &
&    '            Ionic phase ', politot
    call wrtout(6,message,'COLL')
    if (unit_out /= 0) then
     call wrtout(unit_out,message,'COLL')
    end if

    write(message,'(a,es19.9)') &
&    '            Total phase ', poltot
    call wrtout(6,message,'COLL')
    if (unit_out /= 0) then
     call wrtout(unit_out,message,'COLL')
    end if

    poltot = poltot - 2.0_dp*nint(poltot/2._dp)
    write(message,'(a,es19.9)') &
&    '    Remapping in [-1,1] ', poltot
    call wrtout(6,message,'COLL')
    if (unit_out /= 0) then
     call wrtout(unit_out,message,'COLL')
    end if

!   Transform the phase into a polarization
    fac = 1._dp/(gmod*dtefield%nkstr(idir))
    fac = fac/ucvol
    pol = fac*poltot

    write(message,'(a,a,es19.9,a,a,a,es19.9,a,a)')ch10,&
&    '           Polarization ', pol,' (a.u. of charge)/bohr^2',ch10,&
&    '           Polarization ', pol*(e_Cb)/(Bohr_Ang*1d-10)**2,&
&         ' C/m^2',ch10
    call wrtout(6,message,'COLL')
    if (unit_out /= 0) then
     call wrtout(unit_out,message,'COLL')
    end if

   end if   ! option == 1 or option == 3


!  Write the ddk WF to a file

   if ((option == 2).or.(option == 3)) then

    pertcase = idir + 3*natom
    response = 1
    wff2nm=trim(dtfil%filnam_ds(4))//'_1WF'
    call appdig(pertcase,wff2nm,fiwf1o)
    mxfh = 0 ; nxfh = 0 ; nqpt = 1 ; nstep = 1
    allocate(xfhist(3,natom+4,2,mxfh),resid(mband*nkpt*nsppol))
    xfhist(:,:,:,:) = zero
    resid(:) = zero

    call outwf(cg1,dtfil,dtset,eig_dum,fiwf1o,hdr,kg,dtset%kptns,&
&    mband,mkmem,mpi_enreg,mpw,mxfh,natom,dtset%nband,dtset%nfft,&
&    dtset%ngfft,nkpt,npwarr,nqpt,nspinor,nsppol,nstep,&
&    nxfh,occ_dum,resid,response,wffnow,wfs,xfhist)
    deallocate(xfhist,resid)

   end if  ! option == 2 or option == 3

  end if   ! rfdir(idir) == 1

 end do    ! Close loop over idir

! Compute polarization in cartesian coordinates
 if ((dtset%rfdir(1) == 1).and.(dtset%rfdir(2) == 1).and.&
&    (dtset%rfdir(3) == 1)) then

 if(usepaw.ne.1) pelev=zero

  call polcart(pel,pel_cart,pelev,pion,pion_cart,3,&
&              ptot_cart,rprimd,ucvol,unit_out)
  call polcart(pel,pel_cart,pelev,pion,pion_cart,3,&
&              ptot_cart,rprimd,ucvol,6)

 end if

 deallocate(dtm,pwind_k,pwnsfac_k,sflag_k)
 deallocate(cg1_k)
 if (option > 1) deallocate(cg1,eig_dum,occ_dum)


end subroutine berryphase_new
!!***
