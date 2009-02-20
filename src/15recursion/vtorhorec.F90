!{\src2tex{textfont=tt}}
!!****f* ABINIT/vtorhorec
!! NAME
!! vtorhorec
!! 
!! FUNCTION
!! This routine computes the new density from a fixed potential (vtrial)
!! using a recursion method
!! 
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  densymop_gs <type(dens_sym_operator_type)>=the density symmetrization
!!   operator (ground-state symmetries)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  irrzon(nfft**(1-1/nsym),2,nspden/nsppol)=irreducible zone data
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfftf=number of fft grid points
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetry elements in space group
!!  phnons(2,nfft**(1-1/nsym),nspden/nsppol)=nonsymmorphic translation phases
!!  ucvol=unit cell volume in bohr**3.
!!  vtrial(nfft,nspden)=INPUT Vtrial(r).
!!  rmet=define the metric : rprimd*(transpose(rprimd)) 
!! 
!! OUTPUT
!!  ek=kinetic energy part of total energy.
!!  enl=nonlocal pseudopotential part of total energy.
!!  entropy=entropy due to the occupation number smearing (if metal)
!!  e_eigenvalues=Sum of the eigenvalues - Band energy (Hartree)
!!  fermie=fermi energy (Hartree)
!!  grnl(3*natom)=stores grads of nonlocal energy wrt length scales
!!   (3x3 tensor) and grads wrt atomic coordinates (3*natom)
!! 
!! SIDE EFFECTS
!!  rhog(2,nfft)=array for Fourier transform of electron density
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!! 
!! PARENTS
!!      scfcv
!! 
!! CHILDREN
!!      getngrec, green_kernel, recursion, entropyrec, fermisolverec, xcomm_init, xsum_mpi, 
!!      wrtout, symrhg, leave_new, pre_scatter, xmax_mpi, status, timab, time_accu, xallgatherv_mpi
!! 
!! NOTES
!!  at this time :
!!       - must change the choosen entropy (ie entropy = entropy_test). 
!!       - must deactivate computation of ekin, and use instead a computation of free energy. 
!!
!!       - symetrie usage not implemented (densymop_gs, irrzon not used and nsym should be 1)
!!       - natom seems totaly unuseful in the method
!!       - spin-polarized not implemented (nsppol must be 1, nspden ?)
!!       - phnons ?
!!       - need a rectangular box (ngfft(1)=ngfft(2)=ngfft(3))
!!
!!       - enl and grnl are not computed (set to 0)
!! 
!! 
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine vtorhorec(densymop_gs,dtfil,dtset,&
&  ek,enl,entropy,e_eigenvalues,fermie,&
&  grnl,irrzon,mpi_enreg,natom,nfftf,nspden,nsppol,nsym,phnons,&
&  rhog, rhor,ucvol, vtrial, rmet, quit,get_ek,get_entropy)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_15recursion, except_this_one => vtorhorec
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: get_ek,get_entropy,natom,nfftf,nspden,nsppol,nsym
 integer,intent(inout) :: quit
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: e_eigenvalues,ek,enl,entropy,fermie
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(dens_sym_operator_type),intent(in) :: densymop_gs
!arrays
 integer,intent(in) :: irrzon((dtset%nfft)**(1-1/nsym),2,nspden/nsppol)
 real(dp),intent(in) :: phnons(2,(dtset%nfft)**(1-1/nsym),nspden/nsppol)
 real(dp),intent(in) :: rmet(3,3),vtrial(nfftf,nspden)
 real(dp),intent(inout) :: rhog(2,nfftf)
 real(dp),intent(out) :: grnl(3*natom),rhor(nfftf,dtset%nspden)

!Local variables-------------------------------
!pour calcul de la troncature ; a garder ??? 
!for parallel version
!pour calcul de la troncature ; a garder ??? 
!scalars
 integer,parameter :: level=7
 integer,save :: first=1,max_pt,max_pt_group_band,max_pt_group_fft,min_pt
 integer,save :: min_pt_group_band,min_pt_group_fft,nfftrec,ngfftrec(18)=0
 integer,save :: ntranche,ntranche_group_fft,reste
 integer :: affich,calcul_den,calcul_fermie,ierr,iexit,ii,ik,iklocal,index
 integer :: ipoint,ipointlocal,isign,jj,kk,kx,ky,kz,largest_ngfft,max_k
 integer :: maxpow11,maxpow2,maxpow3,maxpow5,maxpow7,me_band,me_fft,me_rec
 integer :: min_k,mk,mmsrch,modi,modj,modk,n1,n2,n3,n_pt_integ_entropy,nfftot
 integer :: nk,nn,nproc_band,nproc_fft,nproc_rec,nrec,nrec_k,ntranche_k
 integer :: old_paral_level,prtvol,reste_k,return_count,return_count2,spaceComm
 integer :: tim_fourdp,trotter,xx,yy,zz
 real(dp),save :: fermie_local
 real(dp) :: beta,drho,drhomax,dummypotmin,ekink,entropy1,entropy2,entropy3
 real(dp) :: entropy4,entropy_local_test,entropy_local_test1
 real(dp) :: entropy_local_test2,entropy_local_test3,entropy_local_test4
 real(dp) :: entropy_test,entropy_test1,entropy_test2,entropy_test3
 real(dp) :: entropy_test4,entropylocal,entropylocal1,entropylocal2
 real(dp) :: entropylocal3,entropylocal4,free_energy,free_energy1,free_energy2
 real(dp) :: free_energy3,free_energy4,free_energy_local,free_energy_local1
 real(dp) :: free_energy_local2,free_energy_local3,free_energy_local4
 real(dp) :: free_energy_local_test,free_energy_test,inf_ucvol,intrhov
 real(dp) :: long_troncat,nelect,potmin,rtroncat,rtrotter,toldrho,tolrec,tsmear
 real(dp) :: xmax
 character(len=4) :: tag
 character(len=500) :: message,timename
 type(MPI_type) :: mpi_enreg_rec
!arrays
 integer :: displs(1:mpi_enreg%nproc_band),ngfft_ek(18)
 integer :: recvcounts(1:mpi_enreg%nproc_band)
 integer,allocatable :: srch(:)
 real(dp) :: inf_rmet(3,3)
 real(dp) :: pot(0:dtset%ngfft(1)-1,0:dtset%ngfft(2)-1,0:dtset%ngfft(3)-1,1)
 real(dp) :: smallpot(0:dtset%ngfft(1)-1,0:dtset%ngfft(2)-1,0:dtset%ngfft(3)/dtset%ngfft(10)-1,1)
 real(dp) :: tsec(2),tsec2(2)
 real(dp),allocatable :: T_p(:,:,:),ZT_p(:,:,:,:),ZT_p_G(:,:,:,:)
 real(dp),allocatable :: ZT_ptempo(:,:),a_k(:,:),alocal(:,:),b2_k(:,:)
 real(dp),allocatable :: b2local(:,:),exppot(:,:,:),exppotloc(:,:,:)
 real(dp),allocatable :: rholocal(:)

! *************************************************************************

!DEBUG echo input variable
!!$write(6,*)' '
!!$write(6,*) 'enter vtorhorec : echo input variable'
!!$write(6,*) 'ek', ek
!!$write(6,*) 'enl', enl
!!$write(6,*) 'entropy', entropy
!!$write(6,*) 'fermi', fermie_local
!!$write(6,*) 'grnl', grnl(1), grnl(2*natom+1)
!!$write(6,*) 'irrzon',irrzon(1,1,1),irrzon(1,2,1),irrzon(dtset%nfft**(1-1/nsym),1,1)
!!$write(6,*) 'natom',natom
!!$write(6,*) 'nfftf',nfftf
!!$write(6,*) 'nspden',nspden
!!$write(6,*) 'nsppol',nsppol
!!$write(6,*) 'nsym',nsym
!!$write(6,*) 'phnons',phnons(1,1,1),phnons(2,1,1),phnons(1,dtset%nfft**(1-1/nsym),1)
!!$write(6,*) 'ucvol',ucvol
!!$write(6,*) 'vtrial',vtrial(1,1),vtrial(3,1)
!!$write(6,*) 'rmet', rmet(:,1)
!!$write(6,*) 'quitrec', quit
!!$write(6,*) 'get_ek', get_ek
!!$write(6,*) 'get_entropy', get_entropy
!!$
!!$call leave_new('COLL')
!ENDDEBUG
 
 call timab(182,1,tsec2)
 call timab(185,1,tsec2)
 
 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if prtvol=-level
 prtvol=dtset%prtvol
 if(prtvol==-level)then
  write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' vtorho : enter '
  call wrtout(06,message,'PERS')
 end if
 
!DEBUG
!call int2char4(mpi_enreg%me,tag)
!timename='time-analysis'//'_P-'//tag
!if(first==1)then
!open(71,file=timename)
!else
!open(71,file=timename, position='append')
!endif
!ENDDEBUG
 
!################################################################################################## 
!parameters for the recursion method 
 trotter = dtset%recptrott
 nrec = dtset%recnrec
!mk = dtset%useric
 nrec_k = nrec 
!dtset%userid !anyway, the recursion is stopped when values are converged
!: taking the same number of recursion shouldn't be a problem. 
!userid is now used in entropyrec.
 if (dtset%recnpath<=0) then
  n_pt_integ_entropy = 25
 else
  n_pt_integ_entropy = dtset%recnpath
 end if
 rtroncat = dtset%recrcut
 toldrho = dtset%rectolden
 tolrec = toldrho*1.d-2
 if(first==1)then
  fermie_local = dtset%recefermi !initial guess for fermie
 end if
 
!for allowing only kinetic energy computation or to freeze fermie (not implemented)
 calcul_den = 1
 calcul_fermie = 1
 
!##################################################################################################  
!initialisation
 enl = 0.d0
 grnl = 0.d0
 
 nelect = dtset%nelect
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 nfftot=dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)
 me_rec=mpi_enreg%me ; nproc_rec=mpi_enreg%nproc !this should use only fft and band- parallelism - no kpt parallelism
 me_fft = mpi_enreg%me_fft ; nproc_fft = mpi_enreg%nproc_fft
 me_band = mpi_enreg%me_band ; nproc_band = mpi_enreg%nproc_band
 
 mpi_enreg_rec=mpi_enreg
 mpi_enreg_rec%nproc_fft=1; mpi_enreg_rec%me_fft=0 ! no fft-parallelism in the recursion

!DEBUG (echo input for parallel version)
!write(6,*)' '
!write(6,*)'me_rec',me_rec
!write(6,*)'me_fft',mpi_enreg%me_fft
!write(6,*)'me_band',mpi_enreg%me_band
!write(6,*)'me',mpi_enreg%me
!write(6,*)'nproc',mpi_enreg%nproc
!write(6,*)'nproc_fft',mpi_enreg%nproc_fft
!write(6,*)'nproc_band',mpi_enreg%nproc_band
!write(6,*)'num_group_fft',mpi_enreg%num_group_fft
!write(6,*)'num_group',mpi_enreg%num_group
!write(6,*)'spaceComm',mpi_enreg%spaceComm 
!!write(6,*)'mpi_enreg%fft_comm(mpi_enreg%num_group_fft)',mpi_enreg%fft_comm(mpi_enreg%num_group_fft)
!write(6,*)'comm_fft',mpi_enreg%comm_fft
!write(6,*)'paral_compil_respfn',mpi_enreg%paral_compil_respfn 
!ENDDEBUG
 
 smallpot = reshape(vtrial,(/ n1,n2,n3/nproc_fft,1 /))
 pot(:,:,0:n3/nproc_rec-1,1) = smallpot(:,:,:,1)
 call pre_scatter(smallpot,pot,n1,n2,n3,mpi_enreg,'gather')
!!$ write(6,*)'smallpot', smallpot(1,1,:)
!!$ write(6,*)' '
!!$ write(6,*)'pot',pot(1,1,:)
!call leave_new('COLL')
 if(dtset%rectesteg==1)then !electron gas
  pot = 0.d0
 end if
 tsmear = dtset%tsmear
 beta = 1/tsmear
 if (trotter == 0) then
  rtrotter = 0.5d0
 else
  rtrotter = real(trotter,dp)
 end if
!infinitesimal metric
 do ii =1,3
  inf_rmet(ii,:) = one/dble(dtset%ngfft(ii))*rmet(ii,:)/dble(dtset%ngfft(1:3))
 end do
 inf_ucvol=ucvol/nfftot
 
!DEBUG : echo the values 
!if(prtvol==-level)then
!write(message,'(a)') '  version 11 '
!call wrtout(06,message,'COLL')
!write(message,'(a,d)') '  temp (Hartree)    ', tsmear
!call wrtout(06,message,'COLL')
!write(message,'(a,d)') '  beta (1/Hartree)  ', beta
!call wrtout(06,message,'COLL')
!write(message,'(a,d)') '  mu                ', fermie_local
!call wrtout(06,message,'COLL')
!write(message,'(a,d)') '  betamu            ',  beta*fermie_local
!call wrtout(06,message,'COLL')
!write(message,'(a,d)') '  nelect            ', nelect
!call wrtout(06,message,'COLL')
!write(message,'(a,i)')'  trotter            ', trotter
!call wrtout(06,message,'COLL')
!write(message,'(a,i)')'  nrec               ', nrec
!call wrtout(06,message,'COLL')
!write(message,'(a,3i5)')'  nb_point             ', n1, n2, n3
!call wrtout(06,message,'COLL')
!write(message,'(a,3d12.3)')'  pas_maillage       ',sqrt(inf_rmet(1,1)),sqrt(inf_rmet(2,2)),sqrt(inf_rmet(3,3))
!call wrtout(06,message,'COLL')
!write(message,'(a,d)')'  rtroncat           ',rtroncat
!call wrtout(06,message,'COLL')
!write(message,'(a,i)')'  Ec max             ', mk
!call wrtout(06,message,'COLL')
!write(message,'(a,i)')'  nrec_k             ',nrec_k
!call wrtout(06,message,'COLL')
!write(message,'(a,i)')'  calcul densite     ', calcul_den 
!call wrtout(06,message,'COLL')
!write(message,'(a,i)')'  calcul entropie    ', get_entropy 
!call wrtout(06,message,'COLL')
!write(message,'(a,i)')'  calcul ekin        ', get_ek
!call wrtout(06,message,'COLL')
!write(message,'(a,i)')'  calcul mu          ', calcul_fermie
!call wrtout(06,message,'COLL')
!write(message,'(a,i)')'  niveau paralellisme', mpi_enreg%paral_level
!call wrtout(06,message,'COLL')
!end if
!ENDDEBUG 
 
 
!################################################################################################
!computation of the fourier transform of the Green kernel (only once)
 
!troncation of the box : determine new dimension
!the method is similar to the one used in getng (except that ecut and xboxcutmin give no constraint, 
!and symmetries are not handled)
 if(first == 1)then
  if ( rtroncat>tol14) then  !default value rtroncat = 0.d0 means no troncation 
   long_troncat = 2*rtroncat + sqrt(beta/rtrotter) !sqrt(beta/trotter) for guess - should be modified
   call getngrec(dtset%ngfft,inf_rmet,ngfftrec,nfftrec,rtroncat+0.25d0*sqrt(beta/rtrotter)) !1/4*sqrt(beta/trotter) for guess - should be modified
  else
   ngfftrec(:) = dtset%ngfft(:) 
!  For now, recursion method doesn't use paralelism on FFT - which would require a great number of processors 
   nfftrec=ngfftrec(1)*ngfftrec(2)*ngfftrec(3)  
   ngfftrec(9)=0              ! paral
   ngfftrec(10)=1             ! nproc_rec
   ngfftrec(11)=0             ! me_rec
   ngfftrec(12)= ngfftrec(2)  ! n2proc
   ngfftrec(13)= ngfftrec(3)  ! n3proc
   nfftrec=ngfftrec(1)*ngfftrec(2)*ngfftrec(3)   
  end if
  
! DEBUG
! if(prtvol==-level)then
! long_troncat = sqrt(inf_rmet(1,1))* ngfftrec(1)
! write(message,'(a,3i5)')'  n_cell               ', ngfftrec(1:3)
! call wrtout(06,message,'COLL')
! write(message,'(a,d)')'  long_troncat       ', long_troncat
! call wrtout(06,message,'COLL')
! end if
! ENDDEBUG
  
 end if
 
!that part could use fft-parallelism
 allocate(ZT_p(1:2,0: ngfftrec(1)-1,0: ngfftrec(2)-1,0: ngfftrec(3)-1))
 allocate(ZT_ptempo(1:2,0: nfftrec-1))
 call green_kernel(ZT_ptempo,inf_rmet,inf_ucvol,rtrotter/beta,mpi_enreg_rec,ngfftrec,nfftrec,prtvol)
 do ii = 0, ngfftrec(1)-1
  do jj = 0, ngfftrec(2)-1
   do kk = 0, ngfftrec(3)-1
    ZT_p(:,ii,jj,kk) =  real(nfftrec,dp)*ZT_ptempo(:,ii+ ngfftrec(1)*jj+( ngfftrec(1)* ngfftrec(2))*kk)
   end do
  end do
 end do
 deallocate(ZT_ptempo)
!end computation of the fourier transform of the Green kernel
 
!###################################################################################
!computation of exp( -beta*pot/(4.d0*rtrotter) ) (only once)
 
 allocate(exppot(0:n1-1,0:n2-1,0:n3-1))
 do ii = 0,n1-1
  do jj = 0,n2-1
   do kk = 0,n3-1
    exppot(ii,jj,kk) = exp( -beta*pot(ii,jj,kk,1)/(4.d0*rtrotter) )
   end do
  end do
 end do
 
!end computation of exp( -beta*pot/(4.d0*rtrotter) )

!################################################################################### 
!determining which point will compute that proc

!paralelism using the fft group
 ntranche_group_fft = nfftf
 min_pt_group_fft = me_fft*ntranche_group_fft
 max_pt_group_fft = (me_fft+1)*ntranche_group_fft-1

!DEBUG
!!$if(dtset%userrd >= 10.d0)then !use symetrie
!!$n1=n1/2
!!$n2=n2/2
!!$n3=n3/2 
!!$ntranche_group_fft =  ntranche_group_fft/8
!!$min_pt_group_fft=me_fft*ntranche_group_fft
!!$max_pt_group_fft = (me_fft+1)*ntranche_group_fft-1
!!$endif
!ENDDEBUG
 
!paralelism using the band communicator (not used in the recursion)
!if(first==1)then
 if(nproc_band/=0)then
  reste = modulo(ntranche_group_fft,nproc_band)
  ntranche = ntranche_group_fft/nproc_band
  displs(1) = 0
  do ii=0,nproc_band-1
   if (ii<reste) then
    recvcounts(ii+1) = ntranche + 1
    if(ii+1/=nproc_band)displs(ii+2) = displs(ii+1) + recvcounts(ii+1)
   else
    recvcounts(ii+1) = ntranche
    if(ii+1/=nproc_band)displs(ii+2) = displs(ii+1) + recvcounts(ii+1)
   end if
  end do
  ntranche = recvcounts(me_band+1)
  min_pt_group_band = displs(me_band+1)
  max_pt_group_band = displs(me_band+1) + ntranche - 1
  min_pt =  min_pt_group_fft + min_pt_group_band
  max_pt =  min_pt_group_fft + max_pt_group_band
 else ! no fft parallelism implies that nproc_band = 0
  ntranche = ntranche_group_fft
  min_pt_group_band = 0
  max_pt_group_band = ntranche - 1
  min_pt =  min_pt_group_fft + min_pt_group_band
  max_pt =  min_pt_group_fft + max_pt_group_band
 end if
!DEBUG
!write(6,*)'  min_pt_group_band',  min_pt_group_band
!write(6,*)'  max_pt_group_band',  max_pt_group_band
!write(6,*)'ntranche',ntranche, ntranche_group_fft
!write(6,*)'minmax',min_pt,max_pt
!ENDDEBUG
!end if
 
 call timab(185,2,tsec2)
 
!###################################################################################
!main loop
 
 noden : if(calcul_den/=0)then

  allocate(rholocal(1:ntranche))
  allocate(alocal(0:nrec,1:ntranche),b2local(0:nrec,1:ntranche))
  
  ipointlocal = 1
  
  graou2 : do kk = 0,n3-1
   do jj = 0,n2-1
    do ii = 0,n1-1
     ipoint = ii+jj*n1+kk*n1*n2
     if ((ipoint>=min_pt).and.(ipoint<=max_pt))then !computation done by that proc
      if ( ngfftrec(1) >= n1 .and. ngfftrec(2) >= n2 .and. ngfftrec(2) >= n2 ) then
       call status(ipointlocal,dtfil%filstat,iexit,level,'call recursion')
!      DEBUG       
!      call timab(182,1,tsec2)
!      ENDDEBUG

!      DEBUG
!      !$write(6,*)' '
!      !$write(6,*)'before recursion : echo input variables'
!      !$write(6,*)'exppot', exppot(1,1,:)
!      !$write(6,*)'coord',ii,jj,kk
!      !$write(6,*)'nrec,fermie,tsmear,trotter',nrec,fermie,tsmear,trotter
!      !$write(6,*)'ZT_p', Zt_p(1,1,1,:)
!      !$write(6,*)'get_rec_coef,prtvol',1,prtvol
!      !$write(6,*)'nfft,ngfft',nfftrec,ngfftrec
!      !$write(6,*)'rmet(:,1)',rmet(:,1)
!      !$write(6,*)'inf_ucvol,tim_fourdp',inf_ucvol,0
!      ENDDEBUG
       
       tim_fourdp=6
       call recursion(exppot,ii,jj,kk, &
&       alocal(:,ipointlocal), &
&       b2local(:,ipointlocal), &
&       rholocal(ipointlocal),&
&       nrec, fermie_local,tsmear,trotter, &
&       ZT_p, &
&       tolrec, &
&       1,prtvol,&
&       mpi_enreg_rec,nfftrec,ngfftrec,rmet,inf_ucvol,tim_fourdp) 
!      DEBUG
!      call timab(182,4,tsec2)
!      call timab(182,2,tsec2)
!      if(modulo(ipointlocal,200)==1)then
!      call  time_accu(187,return_count2,tsec2)
!      write(6,*)'recursion',ipointlocal,return_count2,tsec2(1)
!      end if
!      ENDDEBUG
      else !we use a troncation
       allocate(exppotloc(0:ngfftrec(1)-1,0:ngfftrec(2)-1,0:ngfftrec(3)-1))
       do xx = 0,ngfftrec(1)-1
        do yy = 0,ngfftrec(2)-1
         do zz = 0,ngfftrec(3)-1         
          modi = modulo(ii+xx-ngfftrec(ii)/2,n1)
          modj = modulo(jj+yy-ngfftrec(ii)/2,n2)
          modk = modulo(kk+zz-ngfftrec(ii)/2,n2)
          exppotloc(xx,yy,zz) = exppot(modi,modj,modk)
         end do
        end do
       end do
       call status(ipointlocal,dtfil%filstat,iexit,level,'call recursion')
!      DEBUG
!      call timab(182,1,tsec2)
!      ENDDEBUG
       tim_fourdp=6
       call recursion(exppotloc,ngfftrec(1)/2,ngfftrec(2)/2,ngfftrec(3)/2, &
&       alocal(:,ipointlocal), &
&       b2local(:,ipointlocal), &
&       rholocal(ipointlocal),&
&       nrec, fermie_local,tsmear,trotter, &
&       ZT_p, &
&       tolrec, &
&       1,prtvol,&
&       mpi_enreg_rec,nfftrec,ngfftrec,rmet,inf_ucvol,tim_fourdp)
!      DEBUG
!      call timab(182,4,tsec2)
!      call timab(182,2,tsec2)
!      if(modulo(ipointlocal,200)==1)then
!      call  time_accu(187,return_count2,tsec2)
!      write(6,*)'recursion',ipointlocal,return_count2,tsec2(1)
!      end if
!      ENDDEBUG
       deallocate(exppotloc)
      end if
      
      ipointlocal =   ipointlocal + 1
      
     end if
    end do
   end do
  end do graou2
  
! DEBUG
! call time_accu(182,return_count2,tsec2)
! write(71,*)'fin recursion', return_count2,tsec2(1)
! ENDDEBUG
  
! ###################################
! fermi energy computation

  nomu : if(calcul_fermie/=0)then
!  trouver le mu convenable 
   call status(0,dtfil%filstat,iexit,level,'call fermisolv')
!  DEBUG
!  call timab(183,1,tsec2)
!  ENDDEBUG
   tim_fourdp=0
   call fermisolverec(fermie_local,rholocal,alocal,b2local,nrec, &
&   tsmear,trotter,nelect, &
&   1.d-10,100, &
&   ntranche,mpi_enreg_rec, me_rec,&
&   rmet,inf_ucvol,tim_fourdp)
!  DEBUG
!  call timab(183,2,tsec2)
!  call time_accu(189,return_count2,tsec2)
!  write(6,*)'fermisolve',tsec2
!  ENDDEBUG
  end if nomu

! DEBUG   
! if (prtvol==-level) then
! write(message,'(a)') ' '
! call wrtout(06,message,'COLL')
! write(message,'(a,d)')'  mu                 ',fermie_local
! call wrtout(06,message,'COLL')
! end if
! ENDDEBUG
! end fermi energy computation
  
! ########################
! entropy computation  
  noentropie : if(get_entropy==1)then 
   entropy = zero
   entropy_test = zero
   entropy1 = zero
   entropy2 = zero
   entropy3 = zero
   entropy4 = zero
   entropy_test1 = zero
   entropy_test2 = zero
   entropy_test3 = zero
   entropy_test4 = zero
!  !$   ipointlocal = 1
!  !$   graou2prime : do kk = 0,n3-1
!  !$    do jj = 0,n2-1
!  !$     do ii = 0,n1-1
!  !$      ipoint = ii+jj*n1+kk*n1*n2
!  !$      if ((ipoint>=min_pt).and.(ipoint<=max_pt))then !computation done by that proc
!  !$!      DEBUG
!  !$!       call timab(184,1,tsec2)
!  !$!      ENDDEBUG
!  !$       call entropyrec(exp(beta*fermie/(2.d0*rtrotter))*alocal(:,ipointlocal), &
!  !$&            exp(beta*fermie/rtrotter)*b2local(:,ipointlocal), &
!  !$&            nrec,trotter,entropylocal,1/ucvol,&
!  !$&            prtvol)
!  !$write(6,*)'ipointlocal',ipointlocal,ii,jj,kk
!  !$       entropy = entropy + entropylocal
!  !$!      DEBUG
!  !$!       call timab(184,4,tsec2)
!  !$!       if(modulo(ipoint,500)==0)write(71,*)'entropyrec',ipoint,tsec2
!  !$!       call timab(184,2,tsec2)
!  !$!      ENDDEBUG
!  !$        ipointlocal =   ipointlocal + 1
!  !$      end if
!  !$     end do
!  !$    end do
!  !$   end do graou2prime

!  seek for the min of the path integral
   potmin = minval(pot)
!  DEBUG
!  do ii=1,10
!  entropy = zero
!  potmin=min(-log(real(2.d0*2.d0,dp))/beta+fermie,potmin)
!  write(6,*)'potmin',potmin,exp(-beta/(2.d0*rtrotter)*(potmin-fermie_local))
!  ENDDEBUG


   do ipoint = 1,ntranche
!   DEBUG
!   call timab(184,1,tsec2)
!   ENDDEBUG
!   write(6,*) 'coord', ipoint
    call entropyrec(exp(beta*fermie_local/(2.d0*rtrotter))*alocal(:,ipoint), &
&    exp(beta*fermie_local/rtrotter)*b2local(:,ipoint), &
&    nrec,trotter,entropylocal,2.d0,&
&    prtvol,n_pt_integ_entropy,exp(-beta/(2.d0*rtrotter)*(potmin-fermie_local)),&
&    entropylocal1,entropylocal2,entropylocal3,entropylocal4)
    entropy = entropy + entropylocal
!   !$    entropy_test = entropy_test + entropy_local_test
!   DEBUG
!   ENDDEBUG
    entropy1 = entropy1 + entropylocal1
    entropy2 = entropy2 + entropylocal2
    entropy3 = entropy3 + entropylocal3
    entropy4 = entropy4 + entropylocal4
!   !$    entropy_test1 = entropy_test1 + entropy_local_test1
!   !$    entropy_test2 = entropy_test2 + entropy_local_test2
!   !$    entropy_test3 = entropy_test3 + entropy_local_test3
!   !$    entropy_test4 = entropy_test4 + entropy_local_test4
!   DEBUG
!   ENDDEBUG
!   DEBUG
!   call timab(184,4,tsec2)
!   if(modulo(ipoint,500)==0)write(71,*)'entropyrec',ipoint,tsec2
!   call timab(184,2,tsec2)
!   ENDDEBUG
   end do

   mpi_enreg%paral_level=3
   call xsum_mpi(entropy,mpi_enreg%commcart ,ierr)
   mpi_enreg%paral_level=2

   mpi_enreg%paral_level=3
   call xsum_mpi(entropy_test,mpi_enreg%commcart ,ierr)
   mpi_enreg%paral_level=2

   mpi_enreg%paral_level=3
   call xsum_mpi(entropy1,mpi_enreg%commcart ,ierr)
   call xsum_mpi(entropy2,mpi_enreg%commcart ,ierr)
   call xsum_mpi(entropy3,mpi_enreg%commcart ,ierr)
   call xsum_mpi(entropy4,mpi_enreg%commcart ,ierr)
   call xsum_mpi(entropy_test1,mpi_enreg%commcart ,ierr)
   call xsum_mpi(entropy_test2,mpi_enreg%commcart ,ierr)
   call xsum_mpi(entropy_test3,mpi_enreg%commcart ,ierr)
   call xsum_mpi(entropy_test4,mpi_enreg%commcart ,ierr)
   mpi_enreg%paral_level=2
   
!  DEBUG
!  call time_accu(184,return_count2,tsec2)
!  write(71,*)'fin entropirec', return_count2,tsec2(1)
!  ENDDEBUG
   
!  DEBUG
!  if (prtvol==-level) then
!  write(message,'(a)') ' '
!  call wrtout(06,message,'COLL')
!  write(message,'(a,d)')'  entropy            ',entropy
!  call wrtout(06,message,'COLL')
!  end if
!  write(6,*)'entropy',potmin,entropy,entropy_test
!  write(6,*)'diff entropy', entropy-entropy_test
!  write(6,*)' '
!  write(6,*)'entropy, horizontal path',entropy1
!  write(6,*)'entropy, xmax path',entropy2
!  write(6,*)'entropy, xmin path',entropy3
!  write(6,*)'entropy, zero path',entropy4
!  write(6,*)' '
!  write(6,*)'test, horizontal path',entropy_test1
!  write(6,*)'test, xmax path',entropy_test2
!  write(6,*)'test, xmin path',entropy_test3
!  write(6,*)'test, zero path',entropy_test4
!  ENDDEBUG
!  !$entropy = entropy_test
!  enddo

   free_energy = zero
   free_energy1 = zero
   free_energy2 = zero
   free_energy3 = zero
   free_energy4 = zero
   potmin = minval(pot)
!  DEBUG
!  do ii=1,10
!  free_energy = zero
!  potmin=-log(real(2*10,dp))/beta+fermie +1.d0
!  write(6,*)'potmin',potmin,exp(-beta/(2.d0*rtrotter)*(potmin-fermie_local))
!  ENDDEBUG


   do ipoint = 1,ntranche
!   DEBUG
!   call timab(184,1,tsec2)
!   ENDDEBUG
!   write(6,*) 'coord', ipoint
    call free_energyrec(exp(beta*fermie_local/(2.d0*rtrotter))*alocal(:,ipoint), &
&    exp(beta*fermie_local/rtrotter)*b2local(:,ipoint), &
&    nrec,trotter,free_energy_local,2.d0,& !/ucvol,&
&    prtvol,n_pt_integ_entropy,exp(-beta/(2.d0*rtrotter)*(potmin-fermie_local)),&
&    free_energy_local1,free_energy_local2,free_energy_local3,free_energy_local4)
    free_energy = free_energy + free_energy_local
    free_energy1 = free_energy1 + free_energy_local1
    free_energy2 = free_energy2 + free_energy_local2
    free_energy3 = free_energy3 + free_energy_local3
    free_energy4 = free_energy4 + free_energy_local4
!   DEBUG
!   call timab(184,4,tsec2)
!   if(modulo(ipoint,500)==0)write(71,*)'entropyrec',ipoint,tsec2
!   call timab(184,2,tsec2)
!   ENDDEBUG
   end do

   mpi_enreg%paral_level=3
   call xsum_mpi(free_energy,mpi_enreg%commcart ,ierr)
   mpi_enreg%paral_level=2

   mpi_enreg%paral_level=3
   call xsum_mpi(free_energy1,mpi_enreg%commcart ,ierr)
   call xsum_mpi(free_energy2,mpi_enreg%commcart ,ierr)
   call xsum_mpi(free_energy3,mpi_enreg%commcart ,ierr)
   call xsum_mpi(free_energy4,mpi_enreg%commcart ,ierr)
   mpi_enreg%paral_level=2
   
!  DEBUG
!  call time_accu(184,return_count2,tsec2)
!  write(71,*)'fin entropirec', return_count2,tsec2(1)
!  ENDDEBUG
   
!  DEBUG
!  if (prtvol==-level) then
!  write(message,'(a)') ' '
!  call wrtout(06,message,'COLL')
!  write(message,'(a,d)')'  free_energy            ',free_energy
!  call wrtout(06,message,'COLL')
!  end if
!  write(6,*)' '
!  write(6,*)'free_energy',free_energy
!  write(6,*)' '
!  write(6,*)'test, horizontal path',free_energy1
!  write(6,*)'test, xmax path',free_energy2
!  write(6,*)'test, xmin path',free_energy3
!  write(6,*)'test, zero path',free_energy4
!  ENDDEBUG
!  enddo

   e_eigenvalues=tsmear*(entropy-free_energy) + fermie_local*nelect

  end if noentropie
  
! end entropy computation
  
  deallocate(alocal, b2local)
  
 end if noden

!DEBUG
!!$if(dtset%userrd >= 10.d0)then !use symetrie
!!$n1=2*n1
!!$n2=2*n2
!!$n3=2*n3
!!$endif
!ENDDEBUG

 
!###########################################
!ek computation
 
 noekin : if(get_ek==1)then
  ek = zero

  ngfft_ek(:) = dtset%ngfft(:) 
! For now, recursion method doesn't use paralelism on FFT - which would require a great number of processors 
  ngfft_ek(9)=0              ! paral
  ngfft_ek(10)=1             ! nproc_rec
  ngfft_ek(11)=0             ! me_rec
  ngfft_ek(12)= ngfft_ek(2)  ! n2proc
  ngfft_ek(13)= ngfft_ek(3)  ! n3proc

! !$write(6,*)'green kernel (2)'
  allocate(ZT_p_G(1:2,0:n1-1,0:n2-1,0:n3-1))

  if(.false.)then !(ngfftrec(1)<n1 .or. ngfftrec(2)<n2 .or. ngfftrec(3)<n3)then ! T_p should not be truncated in ekin computation > new computation of ZT_p
   allocate(ZT_ptempo(1:2,0:nfftot-1))
   call green_kernel(ZT_ptempo,inf_rmet,inf_ucvol,rtrotter/beta,mpi_enreg_rec,dtset%ngfft,nfftot,0)
   do ii = 0,n1-1
    do jj = 0,n2-1
     do kk = 0,n3-1
      ZT_p_G(:,ii,jj,kk) =  real(nfftot,dp)*ZT_ptempo(:,ii+n1*jj+n1*n2*kk)
     end do
    end do
   end do
   deallocate(ZT_ptempo)
  else
   ZT_p_G(:,:,:,:) = ZT_p(:,:,:,:)
  end if
  
  iklocal = 1
  if(mk>minval((dtset%ngfft(1:3)-1)/2))then 
   mk = minval((dtset%ngfft(1:3)-1)/2)
  end if
  nk = 2*mk + 1
  
  reste_k = modulo(nk**3,nproc_rec) ! no npband paralellism
  ntranche_k = nk**3/nproc_rec
  if (me_rec<reste_k) then
   ntranche_k = ntranche_k + 1
   min_k = me_rec*ntranche_k
   max_k = (me_rec+1)*ntranche_k-1
  else
   min_k = me_rec*ntranche_k + reste_k
   max_k = (me_rec+1)*ntranche_k-1 + reste_k
  end if
  allocate(a_k(0:nrec_k,1:ntranche_k),b2_k(0:nrec_k,1:ntranche_k))
  
  graouekin : do kx = -mk,mk
   do ky =  -mk,mk
    do kz =  -mk,mk
     ik = (kz+mk)+(ky+mk)*nk+(kx+mk)*nk**2
     if ((ik>=min_k).and.(ik<=max_k))then !computation done by that proc for that k 
      if(kx/=0.or.ky/=0.or.kz/=0)then   !to be modified
!      DEBUG
!      call timab(185,1,tsec2)
!      ENDDEBUG
       tim_fourdp=7
       call recursion(exppot,kx,ky,kz, &
&       a_k(:,iklocal), &
&       b2_k(:,iklocal), &
&       ekink, &
&       nrec_k, fermie_local,tsmear,trotter, &
&       ZT_p_G, &
&       tolrec, &
&       2,prtvol,&
&       mpi_enreg_rec,nfftot,ngfft_ek,rmet,inf_ucvol,tim_fourdp)

       if ((kx == -mk .and. ky == 0 .and. kz == 0).or.&
&       (kx == 0 .and. ky == -mk .and. kz == 0).or.&
&       (kx == 0 .and. ky == 0 .and. kz == -mk).or.&
&       (kx == mk .and. ky == 0 .and. kz == 0).or.&
&       (kx == 0 .and. ky == mk .and. kz == 0).or.&
&       (kx == 0 .and. ky == 0 .and. kz == mk)) then
        write(6,*)'k',kx,ky,kz
        write(6,*)'ekink', ekink
       end if
       
       ek = ek + ekink
!      DEBUG
!      call timab(185,4,tsec2)
!      call timab(185,2,tsec2)
!      if(modulo(ik,500)==0)write(71,*)'recursion for ekin',ik,tsec2
!      ENDDEBUG
      end if
      iklocal = iklocal+1  
     end if
    end do
   end do
  end do graouekin
  
  mpi_enreg%paral_level=3
  call xsum_mpi(ek,mpi_enreg%commcart  ,ierr)
  mpi_enreg%paral_level=2
  
! DEBUG
! call time_accu(185,return_count2,tsec2)
! write(71,*)'fin recursion for ekin', return_count2,tsec2(1)
! DEBUG
! DEBUG
! if (prtvol==-level) then
! write(message,'(a)') ' '
! call wrtout(06,message,'COLL')
! write(message,'(a,d)')'  ekin               ',ek
! call wrtout(06,message,'COLL')
! end if
! ENDDEBUG
  deallocate(a_k,b2_k) 
  deallocate(ZT_p_G) 

 end if noekin
 
!end ek computation
 
 deallocate(ZT_p) 
 deallocate(exppot)
 
 noden2 : if(calcul_den/=0)then
  
! check if the convergence is reached for rho
! DEBUG
! if(mpi_enreg%me==0) then
! write(6,*) 'rho, rholoc(0) for proc 0', rhor(1,1),rholocal(0)
! endif
! ENDDEBUG
  drho = maxval(abs(rhor(min_pt_group_band+1:max_pt_group_band+1,1)-rholocal(:)))
  drhomax = drho
  mpi_enreg%paral_level=3
  call xmax_mpi(drho,drhomax,mpi_enreg%commcart,ierr)
  mpi_enreg%paral_level=2
  if(drhomax<toldrho)then
   quit=quit+1
  else
   quit=0
  end if
! DEBUG
! write(6,*)'convergence',drhomax,quit
! if (drhomax<1.d-15)then
! write(6,*)' '
! write(6,*)'new rho'
! write(6,*)rholocal
! write(6,*)' '
! write(6,*)'old rho'
! write(6,*)rhor(min_pt_group_band+1:max_pt_group_band+1,1)
! end if
! ENDDEBUG
  
  rhor(min_pt_group_band+1:max_pt_group_band+1,1)=rholocal(:)
  if(nproc_band /= 0)then
   call xallgatherv_mpi(rholocal,ntranche,rhor(:,1),recvcounts,displs,mpi_enreg%comm_band,ierr)
  end if
  
  intrhov = inf_ucvol*sum(rholocal*vtrial(min_pt_group_band+1:max_pt_group_band+1,1))
  mpi_enreg%paral_level=3
  call xsum_mpi(intrhov,mpi_enreg%commcart ,ierr)
  mpi_enreg%paral_level=2
! DEBUG
! write(6,*)
! write(6,*)'ek, int(rho*V), ek+int(rho*V)', ek, intrhov,  ek+ intrhov
! write(6,*)'kT*S, kT*sum(ln(...)), diff', tsmear*entropy, tsmear*free_energy, tsmear*(entropy-free_energy) 
! write(6,*)'kT(S-sum(ln(...))) + mu*nelect', tsmear*(entropy-free_energy) + fermie_local*nelect
! write(6,*)'e_eigenvalues', e_eigenvalues
! write(6,*)
! ENDDEBUG
  
  ek = e_eigenvalues-intrhov
  fermie = fermie_local
  
  deallocate( rholocal)
  
 end if noden2
 
!call timab(21,2,tsec)
 
!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
  write(message,'(a)') '  rhor '
  call wrtout(06,message,'PERS')
  write(6,*)rhor   !pas la bonne methode
  write(message,'(a)') ' '
  call wrtout(06,message,'COLL')
! call time_accu(21,return_count,tsec)
  write(message,'(a,2d10.3)')'  temps recursion    ',tsec
  call wrtout(06,message,'COLL')
  write(message,'(a1,a,a1,a,i1,a)') ch10,' vtorho : exit ',&
&  ch10,'  prtvol=-',level,', debugging mode => stop '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if


!call time_accu(21,return_count,tsec)
!write(message,'(a,2d)')'  temps recursion    ',tsec
!call wrtout(06,message,'COLL')
 
 if(first==1)then
  first=0
 end if

!DEBUG 
!close(71)
!ENDDEBUG
 
 call timab(182,2,tsec2)
 
end subroutine vtorhorec
!!***
