!{\src2tex{textfont=tt}}
!!****f* ABINIT/setnoccmmp
!! NAME
!! setnoccmmp
!!
!! FUNCTION
!! PAW+U only:
!! Compute density matrix nocc_{m,m}
!! or
!! Impose value of density matrix using dmatpawu input array, then symetrize it.
!!
!! noccmmp^{\sigma}_{m,m'}=\sum_{ni,nj}[\rho^{\sigma}_{ni,nj}*phiphjint_{ni,nj}]
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (BA,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  compute_dmat= flag: if 1, nocc_{m,mp} is computed
!!  dimdmat=first dimension of dmatpawu array
!!  dmatpawu(dimdmat,dimdmat,nsppol*nspinor,natpawu)=input density matrix to be copied into noccmpp
!!  dmatudiag= flag controlling the use of diagonalization:
!!             0: no diagonalization of nocc_{m,mp}
!!             1: diagonalized nocc_{m,mp} matrix is printed
!!             2: dmatpawu matrix is expressed in the basis where nocc_(m,mp} is diagonal
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  natom=number of atoms in cell
!!  natpawu=number of atoms on which PAW+U is applied
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=number of independant spin components
!!  nsym=number of symmetry elements in space group
!!  ntypat=number of atom types
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  spinat(3,matom)=initial spin of each atom, in unit of hbar/2
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  typat(natom)=type for each atom
!!  impose_dmat= flag: if 1, nocc_{m,mp} is replaced by dmatpawu
!!  useexexch=1 if local-exact-exchange is activated
!!  usepawu=1 if PAW+U is activated
!!
!! OUTPUT
!!   paw_ij(natom)%noccmmp(2*pawtab(itypat)%lpawu+1,2*pawtab(itypat)%lpawu+1,nspden)=density matrix
!!
!! NOTES
!! nocc_{m,mp} is stored as: noccmmp(:,:,1)=   n^{up,up}_{m,mp}
!!                           noccmmp(:,:,2)=   n^{dn,dn}_{m,mp}
!!                           noccmmp(:,:,3)=Re[n^{up,dn}_{m,mp}]
!!                           noccmmp(:,:,4)=Im[n^{up,dn}_{m,mp}]
!!
!! Also ready for future diagonalization of the occupation matrix.
!!
!! PARENTS
!!      pawdenpot,pawprt,scfcv
!!
!! CHILDREN
!!      dgemm,dsyev,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine setnoccmmp(compute_dmat,dimdmat,dmatpawu,dmatudiag,impose_dmat,indsym,natom,natpawu,&
&                     nspden,nspinor,nsppol,nsym,ntypat,paw_ij,pawang,pawprtvol,pawrhoij,pawtab,&
&                     spinat,symafm,typat,useexexch,usepawu)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: compute_dmat,dimdmat,dmatudiag,impose_dmat,natom,natpawu
 integer,intent(in) :: nspden,nspinor,nsppol,nsym,ntypat,useexexch,usepawu
 type(pawang_type),intent(in) :: pawang
 integer,intent(in) :: pawprtvol
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symafm(nsym),typat(natom)
 real(dp),intent(in) :: dmatpawu(dimdmat,dimdmat,nspinor*nsppol,natpawu*impose_dmat)
 real(dp),intent(in) :: spinat(3,natom)
 type(paw_ij_type),intent(inout) :: paw_ij(natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: at_indx,dmatudiag_loc,iafm,iatom,iatpawu,icount,ilm,im1,im2,in1,in2,info,irot,ispden
 integer :: irhoij,itypat,jlm,jm,klmn,lcur,ldim,lmax,lmin,lpawu,lwork,nsploop
 real(dp) :: mnorm,mx,my,mz,ntot,nup,ndn,ro,snorm,sx,sy,szp,szm,zarot2
 logical :: antiferro
 character(len=500) :: message
!arrays
 integer :: nsym_used(2)
 real(dp) :: sumocc(2)
 real(dp),allocatable :: hdp(:,:,:),work(:),noccmmptemp(:,:),noccmmp_tmp(:,:,:),eig(:)
 real(dp),allocatable :: rwork(:),hdp2(:,:,:)
 character(len=9),parameter :: dspin(6)=(/"up       ","down     ","up-up    ","down-down","Re[up-dn]","Im[up-dn]"/)
!no_abirules
  type noccmmp_at
   real(dp),pointer :: noccmmp(:,:,:)
  end type
  type(noccmmp_at),allocatable :: tmp(:)

! *********************************************************************

!Tests
 if (usepawu>0.and.useexexch>0) then
  write(message, '(4a)' ) ch10,&
&    ' setnoccmmp: BUG - ',ch10,&
&    '  usepawu>0 and useexexch>0 not allowed !'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
 end if
 if (impose_dmat/=0.and.dimdmat==0) then
  write(message, '(4a)' ) ch10,&
&    ' setnoccmmp: BUG - ',ch10,&
&    '   dmatpawu must be allocated when impose_dmat/=0 !'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
 end if

!Some inits
 if (usepawu==0.and.useexexch==0) return
 antiferro=(nspden==2.and.nsppol==1)
 dmatudiag_loc=dmatudiag
 if (dmatudiag==2.and.(dimdmat==0.or.impose_dmat==0)) dmatudiag_loc=1

!If needed, store dmatpu in suitable format in tmp%noccmmp
 if (usepawu>0.and.impose_dmat/=0) then
  iatpawu=0
  allocate(tmp(natom))
  do iatom=1,natom
   lpawu=pawtab(typat(iatom))%lpawu
   if (lpawu/=-1) then
    iatpawu=iatpawu+1
    if (nspden/=4) then
     allocate(tmp(iatom)%noccmmp(2*lpawu+1,2*lpawu+1,nsppol))
     tmp(iatom)%noccmmp(1:2*lpawu+1,1:2*lpawu+1,1:nsppol)=&
&     dmatpawu(1:2*lpawu+1,1:2*lpawu+1,1:nsppol,iatpawu)
    else
     allocate(tmp(iatom)%noccmmp(2*lpawu+1,2*lpawu+1,nspden))
     snorm=sqrt(spinat(1,natom)**2+spinat(1,iatom)**2+spinat(3,iatom)**2)
     if (snorm>tol12) then
      sx=half*spinat(1,iatom)/snorm
      sy=half*spinat(2,iatom)/snorm
      szp=half*(one+spinat(3,iatom)/snorm)
      szm=half*(one-spinat(3,iatom)/snorm)
     else
      sx=zero;sy=zero
      szp=one;szm=zero
     end if
     do im2=1,2*lpawu+1
      do im1=1,2*lpawu+1
       nup=dmatpawu(im1,im2,1,iatpawu);ndn=dmatpawu(im1,im2,2,iatpawu)
       tmp(iatom)%noccmmp(im1,im2,1)=nup*szp+ndn*szm
       tmp(iatom)%noccmmp(im1,im2,2)=nup*szm+ndn*szp
       tmp(iatom)%noccmmp(im1,im2,3)=(nup-ndn)*sx
       tmp(iatom)%noccmmp(im1,im2,4)=(ndn-nup)*sy
      end do
     end do
    end if
   end if
  end do
 endif  ! impose_dmat/=0

!Print message
 if (usepawu>0.and.impose_dmat/=0) then
  if (dmatudiag_loc/=2) then
   write(message,'(6a)') ch10,'Occupation matrix for correlated orbitals is kept constant',ch10,&
&                             'and equal to dmatpawu from input file !',ch10,&
&                             '----------------------------------------------------------'
  else
   write(message,'(6a)') ch10,'Occupation matrix for correlated orbitals is imposed',ch10,&
&                             'and equal to dmatpawu in the diagonal basis !',ch10,&
&                             '----------------------------------------------------------'
  end if
  call wrtout(6,message,'COLL')
 end if

 if (usepawu>0.and.dmatudiag_loc/=0) then
  write(message,'(4a)') ch10,'Diagonalized occupation matrix "noccmmp" is printed !',ch10,&
&                            '-------------------------------------------------------------'
  call wrtout(6,message,'COLL')
 endif

!Loops over atoms
 do iatom=1,natom
  itypat=typat(iatom)
  if (useexexch>0) then
   lcur=pawtab(itypat)%lexexch
  else if (usepawu>0) then
   lcur=pawtab(itypat)%lpawu
  end if
  if (lcur/=-1) then

!  ########################################################################################
!  # Compute nocc_mmp
!  ########################################################################################
   if ((usepawu>0.and.compute_dmat/=0).or.useexexch>0) then

    paw_ij(iatom)%noccmmp(:,:,:)=zero

!   Loop over spin components
    do ispden=1,nspden
     allocate(noccmmptemp(2*lcur+1,2*lcur+1));noccmmptemp(:,:)=zero
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
      klmn=pawrhoij(iatom)%rhoijselect(irhoij)
      im1=pawtab(itypat)%klmntomn(1,klmn)
      im2=pawtab(itypat)%klmntomn(2,klmn)
      in1=pawtab(itypat)%klmntomn(3,klmn)
      in2=pawtab(itypat)%klmntomn(4,klmn)
      lmin=pawtab(itypat)%indklmn(3,klmn)
      lmax=pawtab(itypat)%indklmn(4,klmn)
      if (nspden==1) then
       ro=half*pawrhoij(iatom)%rhoijp(irhoij,1)
      else if (nspden==2) then
       ro=pawrhoij(iatom)%rhoijp(irhoij,ispden)
      else
!      Non-collinear magnetism: transfer rhoij from (n,m) storage to n^{alpha,beta}
       if (ispden==1) then
        ro=half*(pawrhoij(iatom)%rhoijp(irhoij,1)+pawrhoij(iatom)%rhoijp(irhoij,4))
       else if (ispden==2) then
        ro=half*(pawrhoij(iatom)%rhoijp(irhoij,1)-pawrhoij(iatom)%rhoijp(irhoij,4))
       else if (ispden==3) then
        ro=half*pawrhoij(iatom)%rhoijp(irhoij,2)
       else
        ro=-half*pawrhoij(iatom)%rhoijp(irhoij,3)
       end if
      end if
      if(lmin==0.and.lmax==2*lcur) then
       icount=in1+(in2*(in2-1))/2
       if(pawtab(itypat)%ij_proj<icount)  then
        write(message, '(4a)' ) ch10,&
&        '  setnoccmmp : BUG -',ch10,&
&        '  PAW+U: Problem in the loop for calculating noccmmp !',ch10
        call wrtout(6,message,'COLL')
        call leave_new('COLL')
       end if
       if(in1/=in2) then
        if(im2<=im1) then
         noccmmptemp(im1,im2)=noccmmptemp(im1,im2)+ro*pawtab(itypat)%phiphjint(icount)
        end if
       end if
       if(im2>=im1) then
        paw_ij(iatom)%noccmmp(im1,im2,ispden)=paw_ij(iatom)%noccmmp(im1,im2,ispden) &
&                                            +ro*pawtab(itypat)%phiphjint(icount)
       end if
      end if
     end do ! irhoij
     do im2=1,2*lcur+1
      do im1=1,im2
       paw_ij(iatom)%noccmmp(im1,im2,ispden)=paw_ij(iatom)%noccmmp(im1,im2,ispden) &
&                                           +noccmmptemp(im2,im1)
      end do
     end do
     do im1=1,2*lcur+1
      do im2=1,im1
       paw_ij(iatom)%noccmmp(im1,im2,ispden)=paw_ij(iatom)%noccmmp(im2,im1,ispden)
      end do
     end do
     deallocate(noccmmptemp)
    end do ! ispden

!   Compute total number of electrons per spin
    paw_ij(iatom)%nocctot(:)=zero
    do ispden=1,nspden
     do im1=1,2*lcur+1
      paw_ij(iatom)%nocctot(ispden)=paw_ij(iatom)%nocctot(ispden) &
&                                  +paw_ij(iatom)%noccmmp(im1,im1,ispden)
     end do
    end do

!   Printing of new nocc_mmp
    if (usepawu>0) write(message, '(2a)' )  ch10,'========== LDA+U DATA =================================================== '
    if (useexexch>0) write(message, '(2a)' )ch10,'======= Local ex-exchange (PBE0) DATA =================================== '
    call wrtout(6,message,'COLL')
    write(message,'(2a,i5,a,i4,a)') ch10,"====== For Atom", iatom,&
&    ", occupations for correlated orbitals. l =",lcur,ch10
    call wrtout(6,message,'COLL')
    if(nspden==2) then
     do ispden=1,2
      write(message,'(a,i4,3a,f10.5)') "Atom", iatom,". Occupations for spin ",&
&      trim(dspin(ispden))," =",paw_ij(iatom)%nocctot(ispden)
      call wrtout(6,message,'COLL')
     end do
     write(message,'(a,i4,a,2x,e16.8)') "=> On atom",iatom,", local Mag. is  ",&
&     paw_ij(iatom)%nocctot(2)-paw_ij(iatom)%nocctot(1)
     call wrtout(6,message,'COLL')
    end if
    if(nspden==4) then
     mx= two*paw_ij(iatom)%nocctot(3)
     my=-two*paw_ij(iatom)%nocctot(4)
     mz=paw_ij(iatom)%nocctot(1)-paw_ij(iatom)%nocctot(2)
     ntot=paw_ij(iatom)%nocctot(1)+paw_ij(iatom)%nocctot(2)
     mnorm=sqrt(mx*mx+my*my+mz*mz)
     write(message,'(a,i4,a,2x,e16.8)') "=> On atom",iatom,", local Mag. x is ",mx
     call wrtout(6,message,'COLL')
     write(message,'(14x,a,2x,e16.8)') "  local Mag. y is ",my
     call wrtout(6,message,'COLL')
     write(message,'(14x,a,2x,e16.8)') "  local Mag. z is ",mz
     call wrtout(6,message,'COLL')
     write(message,'(14x,a,2x,e16.8)') "  norm of Mag. is ",mnorm
     call wrtout(6,message,'COLL')
     write(message,'(14x,a,2x,f10.5)') "  occ. for spin up is ",half*(ntot+mnorm)
     call wrtout(6,message,'COLL')
     write(message,'(14x,a,2x,f10.5)') "  occ. for spin dn is ",half*(ntot-mnorm)
     call wrtout(6,message,'COLL')
    end if
    write(message,'(2a)') ch10,"== Calculated occupation matrix for correlated orbitals:"
    call wrtout(6,message,'COLL')
    do ispden=1,nspden
     write(message,'(3a)') ch10,"Calculated occupation matrix for component ",trim(dspin(ispden+2*(nspden/4)))
     call wrtout(6,message,'COLL')
     do im1=1,lcur*2+1
      write(message,'(12(1x,9(1x,f10.5)))') (paw_ij(iatom)%noccmmp(im1,im2,ispden),im2=1,lcur*2+1)
      call wrtout(6,message,'COLL')
     end do
    end do

   end if ! impose_dmat==0

!  ########################################################################################
!  # Diagonalize nocc_mmp
!  ########################################################################################
   if(usepawu>0.and.dmatudiag_loc>0) then
    lpawu=lcur
    ldim=2*lpawu+1
    lwork=3*ldim-1
    allocate(hdp(ldim,ldim,nspden),noccmmp_tmp(ldim,ldim,nspden),eig(ldim),work(lwork))
    allocate(rwork(3*ldim-1))
    allocate(hdp2(ldim,ldim,nspden))
    hdp=zero
    eig=zero

!   Diagonalization of nocc_mmp
    do ispden=1,nspden
     noccmmp_tmp(:,:,ispden)=paw_ij(iatom)%noccmmp(:,:,ispden)
     call DSYEV('v','u',ldim,noccmmp_tmp(:,:,ispden),ldim,eig,work,lwork,info)
     if(info.ne.0) then
      write(message,'(4a)') ch10,'Error info.ne.0 after DSYEV in diagnoccmmp',ch10
      call wrtout(6,message,'COLL')
      call leave_new('COLL')
     endif
     do ilm=1,2*lpawu+1
      hdp(ilm,ilm,ispden)=eig(ilm)
     enddo
    end do ! ispden

!   Printing of diagonalized matrix and eigenvectors
    do ispden=1,nspden
     write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,' == Diagonalized Occupation matrix'
     if (nspden==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
     if (nspden==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
     if (nspden==4) write(message,fmt='(4a)')     trim(message)," for component ", &
&                   trim(dspin(ispden+2*(nspden/4)))," =="
     call wrtout(6,message,'COLL')
     do ilm=1,2*lpawu+1
      write(message,'(12(1x,9(1x,f10.5)))') (hdp(ilm,jlm,ispden),jlm=1,2*lpawu+1)
      call wrtout(6,message,'COLL')
     end do
    end do ! ispden
    if(abs(pawprtvol)>=1) then
     do ispden=1,nspden
      write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,' == Eigenvectors'
      if (nspden==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
      if (nspden==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
      if (nspden==4) write(message,fmt='(4a)')     trim(message)," for component ", &
&                   trim(dspin(ispden+2*(nspden/4)))," =="
      call wrtout(6,message,'COLL')
      do ilm=1,2*lpawu+1
       write(message,'(12(1x,9(1x,f10.5)))') (noccmmp_tmp(ilm,jlm,ispden),jlm=1,2*lpawu+1)
       call wrtout(6,message,'COLL')
      end do
     end do ! ispden
    endif

!   Back rotation of diagonalized matrix and printing
    do ispden=1,nspden
     call dgemm('n','t',ldim,ldim,ldim,one,hdp(:,:,ispden),&
&               ldim,noccmmp_tmp(:,:,ispden),ldim,zero,hdp2(:,:,ispden),ldim)
     call dgemm('n','n',ldim,ldim,ldim,one,noccmmp_tmp(:,:,ispden),&
&               ldim,hdp2(:,:,ispden),ldim,zero,hdp(:,:,ispden),ldim)
    enddo ! ispden
    if(abs(pawprtvol)>=1) then
     do ispden=1,nspden
      write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,&
&      ' == Rotated back diagonalized matrix'
      if (nspden==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
      if (nspden==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
      if (nspden==4) write(message,fmt='(4a)')     trim(message)," for component ", &
&                    trim(dspin(ispden+2*(nspden/4)))," =="
      call wrtout(6,message,'COLL')
      do ilm=1,2*lpawu+1
       write(message,'(12(1x,9(1x,f10.5)))') (hdp(ilm,jlm,ispden),jlm=1,2*lpawu+1)
       call wrtout(6,message,'COLL')
      end do
     end do ! ispden
    end if

   endif ! dmatudiag_loc

!  ########################################################################################
!  # Impose value of nocc_mmp from dmatpu; symetrize it
!  ########################################################################################
   if (usepawu>0.and.impose_dmat/=0) then

    lpawu=lcur
    nsploop=nsppol;if (nspden==4) nsploop=4

!   Loop over spin components
    do ispden=1,nsploop

!    Loops over components of nocc_mmp
     do jlm=1,2*lpawu+1
      do ilm=1,2*lpawu+1

       if(nsym.gt.1) then

        nsym_used(1:2)=0
        sumocc(1:2)=zero

!       Accumulate values of nocc_mmp over symmetries
        do irot=1,nsym
         if ((symafm(irot)/=1).and.(.not.antiferro)) cycle
         iafm=1;if ((antiferro).and.(symafm(irot)==-1)) iafm=2
         nsym_used(iafm)=nsym_used(iafm)+1
         at_indx=indsym(4,irot,iatom)
         do im2=1,2*lpawu+1
          do im1=1,2*lpawu+1
           sumocc(iafm)=sumocc(iafm)+tmp(at_indx)%noccmmp(ilm,jlm,ispden) &
&                     *pawang%zarot(im1,ilm,lpawu+1,irot)&
&                     *pawang%zarot(im2,jlm,lpawu+1,irot)
          end do
         end do
        end do ! End loop over symmetries

!       Store new values of nocc_mmp
        paw_ij(iatom)%noccmmp(ilm,jlm,ispden)=sumocc(1)/nsym_used(1)

!       Antiferromagnetic case: has to fill up "down" component of nocc_mmp
        if (antiferro.and.nsym_used(2)>0) paw_ij(iatom)%noccmmp(ilm,jlm,2)=sumocc(2)/nsym_used(2)

       else  ! nsym=1

!       Case without symetries
        paw_ij(iatom)%noccmmp(ilm,jlm,ispden)= tmp(iatom)%noccmmp(ilm,jlm,ispden)
       end if

      end do !ilm
     end do !jlm
    end do ! ispden

!   Printing of new nocc_mmp
    do ispden=1,nspden
     if(dmatudiag_loc==2) then
      write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,&
 &    ' == Imposed occupation matrix (in the basis of diagonalization!!)'
     else
      write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,&
 &    ' == Imposed occupation matrix'
     endif
     if (nspden==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
     if (nspden==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
     if (nspden==4) write(message,fmt='(4a)')     trim(message)," for component ", &
&                   trim(dspin(ispden+2*(nspden/4)))," =="
     call wrtout(6,message,'COLL')
     do ilm=1,2*lpawu+1
      write(message,'(12(1x,9(1x,f10.5)))') (paw_ij(iatom)%noccmmp(ilm,jlm,ispden),jlm=1,2*lpawu+1)
      call wrtout(6,message,'COLL')
     end do
    end do

   endif ! impose_dmat/=0

!  ########################################################################################
!  # Rotate imposed occupation matrix in the non-diagonal basis
!  ########################################################################################
   if (usepawu>0.and.impose_dmat/=0.and.dmatudiag_loc==2) then

    lpawu=lcur

!   Rotation of imposed nocc_mmp
    hdp(:,:,:)=paw_ij(iatom)%noccmmp(:,:,:)
    do ispden=1,nspden
     call dgemm('n','t',ldim,ldim,ldim,one,hdp(:,:,ispden),&
&               ldim,noccmmp_tmp(:,:,ispden),ldim,zero,hdp2(:,:,ispden),ldim)
     call dgemm('n','n',ldim,ldim,ldim,one,noccmmp_tmp(:,:,ispden),&
&               ldim,hdp2(:,:,ispden),ldim,zero,hdp(:,:,ispden),ldim)
    enddo ! ispden
    paw_ij(iatom)%noccmmp(:,:,:)=hdp(:,:,:)

!   Printing of rotated imposed matrix
    do ispden=1,nspden
     write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,&
&    ' == Imposed density matrix in original basis'
     if (nspden==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
     if (nspden==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
     if (nspden==4) write(message,fmt='(4a)')     trim(message)," for component ", &
&                   trim(dspin(ispden+2*(nspden/4)))," =="
     call wrtout(6,message,'COLL')
     do ilm=1,2*lpawu+1
      write(message,'(12(1x,9(1x,f10.5)))') (hdp(ilm,jlm,ispden),jlm=1,2*lpawu+1)
      call wrtout(6,message,'COLL')
     end do
    end do ! ispden

   endif ! dmatudiag_loc==2

   if (usepawu>0.and.dmatudiag_loc>0) deallocate(hdp,noccmmp_tmp,eig,work,rwork,hdp2)

  endif ! lcur
 enddo ! iatom

!Memory deallocation
 if (usepawu>0.and.impose_dmat/=0) then
  do iatom=1,natom
   lpawu=pawtab(typat(iatom))%lpawu
   if (lpawu/=-1) deallocate(tmp(iatom)%noccmmp)
  end do
  deallocate(tmp)
 endif

end subroutine setnoccmmp
!!***
