!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtrhomxmn
!! NAME
!! prtrhomxmn
!!
!! FUNCTION
!! If option==1, compute the maximum and minimum of the density (and spin-polarization
!! if nspden==2), and print it.
!! If option==2, also compute and print the second maximum or minimum
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  iout=unit for output file
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  option, see above
!!  rhor(nfft,nspden)=electron density (electrons/bohr^3)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  The tolerance tol12 aims at giving a machine-independent ordering.
!!  (this trick is used in bonds.f, listkk.f, prtrhomxmn.f and rsiaf9.f)
!!
!! PARENTS
!!      clnup1,mkrho,vtorho
!!
!! CHILDREN
!!      leave_new,wrtout,xcomm_init,xmaster_init_fft,xmax_mpi,xmin_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine prtrhomxmn(iout,mpi_enreg,nfft,ngfft,nspden,option,rhor)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,nfft,nspden,option
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,idffz,ierr,ifft,ii,iitems,imn1,imn2,imx1,imx2,index1
 integer :: index2,indsign,iproc,isign,ispden,master,n1,n2,n3,ndffz,nitems
 integer :: old_paral_level,resulti,spaceComm
 real(dp) :: resultr,temp,value1,value2
 character(len=500) :: message
!arrays
 integer,allocatable :: index(:,:,:),index_fft(:,:,:,:)
 real(dp) :: rhomn1(4),rhomn2(4),rhomx1(4),rhomx2(4),ri_rhomn1(3,4)
 real(dp) :: ri_rhomn2(3,4),ri_rhomx1(3,4),ri_rhomx2(3,4),ri_zetmn1(3,2)
 real(dp) :: ri_zetmn2(3,2),ri_zetmx1(3,2),ri_zetmx2(3,2),zeta(2),zetmn1(2)
 real(dp) :: zetmn2(2),zetmx1(2),zetmx2(2)
 real(dp),allocatable :: array(:),coord(:,:,:,:),value(:,:,:)
 real(dp),allocatable :: value_fft(:,:,:,:)

! *************************************************************************

!DEBUG
!write(6,*) ' prhomxmn : enter '
!write(6,*) ' nspden, rhor(1,1), rhor(1,2) ',nspden,rhor(1,1), rhor(1,2)
!ENDDEBUG
 if(option/=1 .and. option/=2)then
  write(message, '(a,a,a,a,i8,a)' ) ch10,&
&  ' prtrhomxmn : BUG -',ch10,&
&  '  Option must be 1 or 2, while it is',option,'.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!--------------------------------------------------------------------------
!One has to determine the maximum and minimum (etc...) values
!over all space, and then output it, as well as to identify
!the point at which it occurs ...
!This will require a bit of data exchange, and correct indirect indexing ...

!For the local processor, find different items :
!maximum and minimum total electron density and locations
!and also spin-polarisation and magnetisation
!also keep the second maximal or minimal value
 if(nspden==1)nitems=1   ! Simply the total density
 if(nspden==2)nitems=5   ! Total density, spin up, spin down, magnetization, zeta
 if(nspden==4)nitems=6   ! Total density, x, y, z, magnetization, zeta

 allocate(value(2,2,nitems),index(2,2,nitems),array(nfft))

 do iitems=1,nitems

! Copy the correct values into the array
! First set of items : the density, for each spin component
  if(iitems<=nspden)then
   array(:)=rhor(:,iitems)
  end if
! Case nspden==2, some computation to be done
  if(nspden==2)then
   if(iitems==3)then ! Spin down
    array(:)=rhor(:,1)-rhor(:,2)
   else if(iitems==4)then  ! Magnetization
    array(:)=2*rhor(:,2)-rhor(:,1)
   else if(iitems==5)then  ! zeta = relative magnetization
    array(:)=(2*rhor(:,2)-rhor(:,1))/rhor(:,1)
   end if
!  Case nspden==4, some other computation to be done
  else if(nspden==4)then
   if(iitems==5)then ! Magnetization
    array(:)=sqrt(rhor(:,2)**2+rhor(:,3)**2+rhor(:,4)**2)
   else if(iitems==6)then ! zeta = relative magnetization
    array(:)=(sqrt(rhor(:,2)**2+rhor(:,3)**2+rhor(:,4)**2))/rhor(:,1)
   end if
  end if

! DEBUG
! write(6,*) ' iitems,array(1:2)=',iitems,array(1:2)
! ENDDEBUG

  do indsign=1,2 ! Find alternatively the maximum and the minimum
   isign=3-2*indsign

!  Initialize the two first values
   value1=array(1) ; value2=array(2)
   index1=1 ; index2=2

!  Ordering, if needed
   if( isign*(value2+tol12) > isign*(value1)) then
    temp=value2 ; value2=value1 ; value1=temp  
    index1=2 ; index2=1
   end if

!  DEBUG
!  write(6,*) ' value1,value2,index1,index2=',value1,value2,index1,index2
!  ENDDEBUG

!  Loop over all points
   do ifft=3,nfft    

    temp=array(ifft)
!   Compares it to the second value
    if( isign*(temp+tol12) > isign*value2 ) then
!    Compare it to the first value
     if( isign*(temp+tol12) > isign*value1 ) then
      value2=value1 ; index2=index1
      value1=temp   ; index1=ifft
     else
      value2=temp   ; index2=ifft
     end if     
    end if

   end do ! ifft

   value(1,indsign,iitems)=value1
   value(2,indsign,iitems)=value2
   index(1,indsign,iitems)=index1
   index(2,indsign,iitems)=index2

!  DEBUG
!  write(6,*) ' value1,value2,index1,index2=',value1,value2,index1,index2
!  ENDDEBUG

  end do ! isign

 end do ! iitems

 deallocate(array)

!-------------------------------------------------------------------
!Enter section for FFT parallel case

 if (mpi_enreg%paral_fft == 1 .or. associated(mpi_enreg%nscatterarr)) then

! Communicate all data to all processors with only two global communications
  old_paral_level=mpi_enreg%paral_level
! DC 20071107:
! These ifs should be useless if xcomm_init was used
! correctly without tricks to change paral_level.
! Please correct this.
  if (.not.associated(mpi_enreg%nscatterarr)) then
   mpi_enreg%paral_level=3
  end if
  call xcomm_init(mpi_enreg,spaceComm)
  if (.not.associated(mpi_enreg%nscatterarr)) then
   if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%comm_fft
  end if
  allocate(value_fft(2,2,nitems,mpi_enreg%nproc_fft))
  allocate(index_fft(2,2,nitems,mpi_enreg%nproc_fft))
  value_fft(:,:,:,:)=zero
  index_fft(:,:,:,:)=zero
  value_fft(:,:,:,mpi_enreg%me_fft + 1)=value(:,:,:)
  index_fft(:,:,:,mpi_enreg%me_fft + 1)=index(:,:,:)
  call xsum_mpi(value_fft,spaceComm,ierr)
  call xsum_mpi(index_fft,spaceComm,ierr)

! Determine the global optimum and second optimum for each item
  do iitems=1,nitems
   do indsign=1,2 ! Find alternatively the maximum and the minimum
    isign=3-2*indsign

!   Initialisation
    value1=value_fft(1,indsign,iitems,1)
    value2=value_fft(2,indsign,iitems,1)
    index1=index_fft(1,indsign,iitems,1)
    index2=index_fft(2,indsign,iitems,1)

!   Loop
    do iproc=1,mpi_enreg%nproc_fft
     do ii=1,2
      if(iproc>1 .or. ii==2)then

       temp=value_fft(ii,indsign,iitems,iproc)
!      Compares it to the second value
       if( isign*(temp+tol12) > isign*value2 ) then
!       Compare it to the first value
        if( isign*(temp+tol12) > isign*value1 ) then
         value2=value1 ; index2=index1
         value1=temp   ; index1=index_fft(ii,indsign,iitems,iproc)
        else
         value2=temp   ; index2=index_fft(ii,indsign,iitems,iproc)
        end if
       end if

      end if ! if(iproc>1 .or. ii==2)
     end do ! ii
    end do ! iproc

    value(1,indsign,iitems)=value1
    value(2,indsign,iitems)=value2
    index(1,indsign,iitems)=index1
    index(2,indsign,iitems)=index2

   end do ! isign
  end do ! iitems

  deallocate(value_fft,index_fft)

  mpi_enreg%paral_level=old_paral_level

 end if ! if (mpi_enreg%paral_fft == 1)

 call xmaster_init_fft(mpi_enreg,master)

!-------------------------------------------------------------------

!Determines the reduced coordinates of the min and max for each item 
 allocate(coord(3,2,2,nitems))
 do iitems=1,nitems
  do indsign=1,2
   do ii=1,2
    index1=index(ii,indsign,iitems)
    i3=(index1-1)/n1/n2
    i2=(index1-1-i3*n1*n2)/n1
    i1=index1-1-i3*n1*n2-i2*n1
    coord(1,ii,indsign,iitems)=dble(i1)/dble(n1)+tol12
    coord(2,ii,indsign,iitems)=dble(i2)/dble(n2)+tol12
    coord(3,ii,indsign,iitems)=dble(i3)/dble(n3)+tol12
!   DEBUG
!   write(6,*)' ii,indsign,iitems,coord(1:3)=',ii,indsign,iitems,coord(:,ii,indsign,iitems)
!   ENDDEBUG
   end do
  end do
 end do

!-------------------------------------------------------------------------
!Output
 if (mpi_enreg%me==master) then
  if(.true.)then
   do iitems=1,nitems

    if(iitems==1) write(message,'(a)')' Total charge density [el/Bohr^3]'
    if(nspden==2)then
     if(iitems==2) write(message,'(a)')' Spin up density      [el/Bohr^3]'
     if(iitems==3) write(message,'(a)')' Spin down density    [el/Bohr^3]'
     if(iitems==4) write(message,'(a)')' Magnetization (spin up - spin down) [el/Bohr^3]'
     if(iitems==5) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
    else if(nspden==4)then
     if(iitems==2) write(message,'(a)')' x component of magnetization [el/Bohr^3]'
     if(iitems==3) write(message,'(a)')' y component of magnetization [el/Bohr^3]'
     if(iitems==4) write(message,'(a)')' z component of magnetization [el/Bohr^3]'
     if(iitems==5) write(message,'(a)')' Magnetization (spin up - spin down) [el/Bohr^3]'
     if(iitems==6) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
    end if
    call wrtout(iout,message,'COLL')

    write(message,'(a,es13.4,a,3f10.4)') ',     Maximum= ',&
&    value(1,1,iitems),'  at reduced coord.',coord(:,1,1,iitems)
    call wrtout(iout,message,'COLL')
    if(option==2)then
     write(message,'(a,es13.4,a,3f10.4)')',Next maximum= ',&
&     value(2,1,iitems),'  at reduced coord.',coord(:,2,1,iitems)
     call wrtout(iout,message,'COLL')
    end if
    write(message,'(a,es13.4,a,3f10.4)') ',     Minimum= ',&
&    value(1,2,iitems),'  at reduced coord.',coord(:,1,2,iitems)
    call wrtout(iout,message,'COLL')
    if(option==2)then
     write(message,'(a,es13.4,a,3f10.4)')',Next minimum= ',&
&     value(2,2,iitems),'  at reduced coord.',coord(:,2,2,iitems)
     call wrtout(iout,message,'COLL')
    end if

   end do ! iitems
  end if

  if(.false.)then

   write(message, '(a,a,1p,e12.4,a,0p,3f8.4)' ) ch10,&
&   ',Min el dens=',value(1,2,1),&
&   ' el/bohr^3 at reduced coord.',coord(:,1,2,1)
   call wrtout(iout,message,'COLL')
   if(option==2)then
    write(message, '(a,1p,e12.4,a,0p,3f8.4)' ) &
&    ',   next min=',value(2,2,1),&
&    ' el/bohr^3 at reduced coord.',coord(:,2,2,1)
    call wrtout(iout,message,'COLL')
   end if
   write(message, '(a,1p,e12.4,a,0p,3f8.4)' )&
&   ',Max el dens=',value(1,1,1),&
&   ' el/bohr^3 at reduced coord.',coord(:,1,1,1)
   call wrtout(iout,message,'COLL')
   if(option==2)then
    write(message, '(a,1p,e12.4,a,0p,3f8.4)' )&
&    ',   next max=',value(2,1,1),&
&    ' el/bohr^3 at reduced coord.',coord(:,2,1,1)
    call wrtout(iout,message,'COLL')
   end if

   if(nspden>=2)then
    write(message, '(a,a,1p,e12.4,a,0p,3f8.4)' ) ch10,&
&    ',Min spin pol zeta=',value(1,2,4+nspden/2),&
&    ' at reduced coord.',coord(:,1,2,4+nspden/2)
    call wrtout(iout,message,'COLL')
    if(option==2)then
     write(message, '(a,1p,e12.4,a,0p,3f8.4)' )&
&     ',         next min=',value(2,2,4+nspden/2),&
&     ' at reduced coord.',coord(:,2,2,4+nspden/2)
     call wrtout(iout,message,'COLL')
    end if
    write(message, '(a,1p,e12.4,a,0p,3f8.4)' )&
&    ',Max spin pol zeta=',value(1,1,4+nspden/2),&
&    ' at reduced coord.',coord(:,1,1,4+nspden/2)
    call wrtout(iout,message,'COLL')
    if(option==2)then
     write(message, '(a,1p,e12.4,a,0p,3f8.4)' )&
&     ',         next max=',value(2,1,4+nspden/2),&
&     ' at reduced coord.',coord(:,2,1,4+nspden/2)
     call wrtout(iout,message,'COLL')
    end if
   end if ! nspden

  end if ! second section always true
  
  if(nspden==2 .and. .false.)then
   write(message,'(a)')&
&   '                               Position in reduced coord.       (  x         y         z )'
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Total  el-den) : [el/Bohr^3]',&
&   rhomn1(1),'  at',ri_rhomn1(1,1),ri_rhomn1(2,1),ri_rhomn1(3,1)
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Spin-up   den) : [el/Bohr^3]',&
&   rhomn1(2),'  at',ri_rhomn1(1,2),ri_rhomn1(2,2),ri_rhomn1(3,2)
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Spin-down den) : [el/Bohr^3]',&
&   zetmn1(1),'  at',ri_zetmn1(1,1),ri_zetmn1(2,1),ri_zetmn1(3,1)
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Spin pol zeta) :   [m/|m|]  ',&
&   zetmn1(2),'  at',ri_zetmn1(1,2),ri_zetmn1(2,2),ri_zetmn1(3,2)
   call wrtout(iout,message,'COLL')
   if(option==2)then                                                            
    write(message,'(a,es13.4,a,3f10.4)')' Next minimum (Total  el-den) : [el/Bohr^3]',&
&    rhomn2(1),'  at',ri_rhomn2(1,1),ri_rhomn2(2,1),ri_rhomn2(3,1)
    call wrtout(iout,message,'COLL')
    write(message,'(a,es13.4,a,3f10.4)')' Next minimum (Spin-up   den) : [el/Bohr^3]',&
&    rhomn2(2),'  at',ri_rhomn2(1,2),ri_rhomn2(2,2),ri_rhomn2(3,2)
    call wrtout(iout,message,'COLL')
    write(message,'(a,es13.4,a,3f10.4)')' Next minimum (Spin-down den) : [el/Bohr^3]',&
&    zetmn2(1),'  at',ri_zetmn2(1,1),ri_zetmn2(2,1),ri_zetmn2(3,1)
    call wrtout(iout,message,'COLL')
    write(message,'(a,es13.4,a,3f10.4)')' Next minimum (Spin pol zeta) :   [m/|m|]  ',&
&    zetmn2(2),'  at',ri_zetmn2(1,2),ri_zetmn2(2,2),ri_zetmn2(3,2)
    call wrtout(iout,message,'COLL')
   end if                                                                       
   write(message,*)' '
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Total  el-den) : [el/Bohr^3]',&
&   rhomx1(1),'  at',ri_rhomx1(1,1),ri_rhomx1(2,1),ri_rhomx1(3,1)
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Spin-up   den) : [el/Bohr^3]',&
&   rhomx1(2),'  at',ri_rhomx1(1,2),ri_rhomx1(2,2),ri_rhomx1(3,2)
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Spin-down den) : [el/Bohr^3]',&
&   zetmx1(1),'  at',ri_zetmx1(1,1),ri_zetmx1(2,1),ri_zetmx1(3,1)
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Spin pol zeta) :   [m/|m|]  ',&
&   zetmx1(2),'  at',ri_zetmx1(1,2),ri_zetmx1(2,2),ri_zetmx1(3,2)
   call wrtout(iout,message,'COLL')
   if(option==2)then                                                            
    write(message,'(a,es13.4,a,3f10.4)')' Next maximum (Total  el-den) : [el/Bohr^3]',&
&    rhomx2(1),'  at',ri_rhomx2(1,1),ri_rhomx2(2,1),ri_rhomx2(3,1)
    call wrtout(iout,message,'COLL')
    write(message,'(a,es13.4,a,3f10.4)')' Next maximum (Spin-up   den) : [el/Bohr^3]',&
&    rhomx2(2),'  at',ri_rhomx2(1,2),ri_rhomx2(2,2),ri_rhomx2(3,2)
    call wrtout(iout,message,'COLL')
    write(message,'(a,es13.4,a,3f10.4)')' Next maximum (Spin-down den) : [el/Bohr^3]',&
&    zetmx2(1),'  at',ri_zetmx2(1,1),ri_zetmx2(2,1),ri_zetmx2(3,1)
    call wrtout(iout,message,'COLL')
    write(message,'(a,es13.4,a,3f10.4)')' Next maximum (Spin pol zeta) :   [m/|m|]  ',&
&    zetmx2(2),'  at',ri_zetmx2(1,2),ri_zetmx2(2,2),ri_zetmx2(3,2)
    call wrtout(iout,message,'COLL')
   end if                                                                       
  end if                                                                        
  
  if(nspden==4 .and. .false.)then                                                             
   write(message,'(a)')&
&   '                               Position in reduced coord.       (  x         y         z )'
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Total  el-den) : [el/Bohr^3]',&
&   rhomn1(1),'  at',ri_rhomn1(1,1),ri_rhomn1(2,1),ri_rhomn1(3,1)
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Magnetizat.-x) :   [m/|m|]  ',&
&   rhomn1(2),'  at',ri_rhomn1(1,2),ri_rhomn1(2,2),ri_rhomn1(3,2)
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Magnetizat.-y) :   [m/|m|]  ',&
&   rhomn1(3),'  at',ri_rhomn1(1,3),ri_rhomn1(2,3),ri_rhomn1(3,3)
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Magnetizat.-z) :   [m/|m|]  ',&
&   rhomn1(4),'  at',ri_rhomn1(1,4),ri_rhomn1(2,4),ri_rhomn1(3,4)
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Spin pol zeta) :   [m/|m|]  ',&
&   zetmn1(1),'  at',ri_zetmn1(1,1),ri_zetmn1(2,1),ri_zetmn1(3,1)
   call wrtout(iout,message,'COLL')
   if(option==2)then                                                            
    write(message,'(a,es13.4,a,3f10.4)')' Next-Minimum (Total  el-den) : [el/Bohr^3]',&
&    rhomn2(1),'  at',ri_rhomn2(1,1),ri_rhomn2(2,1),ri_rhomn2(3,1)
    call wrtout(iout,message,'COLL')
    write(message,'(a,es13.4,a,3f10.4)')' Next-Minimum (Magnetizat.-x) :   [m/|m|]  ',&
&    rhomn2(2),'  at',ri_rhomn2(1,2),ri_rhomn2(2,2),ri_rhomn2(3,2)
    call wrtout(iout,message,'COLL')
    write(message,'(a,es13.4,a,3f10.4)')' Next-Minimum (Magnetizat.-y) :   [m/|m|]  ',&
&    rhomn2(3),'  at',ri_rhomn2(1,3),ri_rhomn2(2,3),ri_rhomn2(3,3)
    call wrtout(iout,message,'COLL')
    write(message,'(a,es13.4,a,3f10.4)')' Next-Minimum (Magnetizat.-z) :   [m/|m|]  ',&
&    rhomn2(4),'  at',ri_rhomn2(1,4),ri_rhomn2(2,4),ri_rhomn2(3,4)
    call wrtout(iout,message,'COLL')
    write(message,'(a,es13.4,a,3f10.4)')' Next-Minimum (Spin pol zeta) :   [m/|m|]  ',&
&    zetmn2(1),'  at',ri_zetmn2(1,1),ri_zetmn2(2,1),ri_zetmn2(3,1)
    call wrtout(iout,message,'COLL')
   end if                                                                       
   write(message,*)' '
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Total  el-den) : [el/Bohr^3]',&
&   rhomx1(1),'  at',ri_rhomx1(1,1),ri_rhomx1(2,1),ri_rhomx1(3,1)
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Magnetizat.-x) :   [m/|m|]  ',&
&   rhomx1(2),'  at',ri_rhomx1(1,2),ri_rhomx1(2,2),ri_rhomx1(3,2)
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Magnetizat.-y) :   [m/|m|]  ',&
&   rhomx1(3),'  at',ri_rhomx1(1,3),ri_rhomx1(2,3),ri_rhomx1(3,3)
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Magnetizat.-z) :   [m/|m|]  ',&
&   rhomx1(4),'  at',ri_rhomx1(1,4),ri_rhomx1(2,4),ri_rhomx1(3,4)
   call wrtout(iout,message,'COLL')
   write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Spin pol zeta) :   [m/|m|]  ',&
&   zetmx1(1),'  at',ri_zetmx1(1,1),ri_zetmx1(2,1),ri_zetmx1(3,1)
   call wrtout(iout,message,'COLL')
   if(option==2)then                                                            
    write(message,'(a,es13.4,a,3f10.4)')' Next-Maximum (Total  el-den) : [el/Bohr^3]',&
&    rhomx2(1),'  at',ri_rhomx2(1,1),ri_rhomx2(2,1),ri_rhomx2(3,1)
    call wrtout(iout,message,'COLL')
    write(message,'(a,es13.4,a,3f10.4)')' Next-Maximum (Magnetizat.-x) :   [m/|m|]  ',&
&    rhomx2(2),'  at',ri_rhomx2(1,2),ri_rhomx2(2,2),ri_rhomx2(3,2)
    call wrtout(iout,message,'COLL')
    write(message,'(a,es13.4,a,3f10.4)')' Next-Maximum (Magnetizat.-y) :   [m/|m|]  ',&
&    rhomx2(3),'  at',ri_rhomx2(1,3),ri_rhomx2(2,3),ri_rhomx2(3,3)
    call wrtout(iout,message,'COLL')
    write(message,'(a,es13.4,a,3f10.4)')' Next-Maximum (Magnetizat.-z) :   [m/|m|]  ',&
&    rhomx2(4),'  at',ri_rhomx2(1,4),ri_rhomx2(2,4),ri_rhomx2(3,4)
    call wrtout(iout,message,'COLL')
    write(message,'(a,es13.4,a,3f10.4)')' Next-Maximum (Spin pol zeta) :   [m/|m|]  ',&
&    zetmx2(1),'  at',ri_zetmx2(1,1),ri_zetmx2(2,1),ri_zetmx2(3,1)
    call wrtout(iout,message,'COLL')
   end if
  end if
 end if

 deallocate(coord,value,index)

end subroutine prtrhomxmn
!!***
