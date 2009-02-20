!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pseudo
!! NAME
!! m_pseudo
!!
!! FUNCTION
!! This module reads a pseudopotential file written in XML.
!! A full example of the building up of a data structure using
!! the SAX paradigm.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (JJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_pseudo

#if defined HAVE_XMLF90
use flib_sax
use m_pseudo_types         ! Data types
#endif

implicit none

private

#if defined HAVE_XMLF90
!
! It defines the routines that are called from xml_parser in response
! to particular events.
!
public  :: begin_element, end_element, pcdata_chunk
private :: die

logical, private  :: in_vps = .false. , in_radfunc = .false.
logical, private  :: in_semilocal = .false. , in_header = .false.
logical, private  :: in_coreCharge = .false. , in_data = .false.
logical, private  :: in_valenceCharge = .false.
logical, private  :: in_pseudowavefun = .false. , in_pswf = .false.

integer, private, save  :: ndata

type(pseudo_t), public, target, save :: pseudo
type(grid_t), private, save        :: grid
type(grid_t), private, save        :: global_grid
!
! Pointers to make it easier to manage the data
!
type(header_t), private, pointer   :: hp
type(vps_t), private, pointer      :: pp
type(pswf_t), private, pointer     :: pw
type(radfunc_t), private, pointer  :: rp

CONTAINS  !===========================================================

!----------------------------------------------------------------------
subroutine begin_element(name,attributes)



!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

character(len=*), intent(in)    :: name
type(dictionary_t), intent(in)  :: attributes

character(len=100)  :: value
integer             :: status


select case(name)

      case ("pseudo")
         pseudo%npots  = 0
         pseudo%npswfs = 0
         global_grid%npts = 0
!         call get_value(attributes,"version",value,status)
!         if (value == "0.5") then
!            print *, "Processing a PSEUDO version 0.5 XML file"
!            pseudo%npots  = 0
!            pseudo%npswfs = 0
!            global_grid%npts = 0
!         else
!            print *, "Can only work with PSEUDO version 0.5 XML files"
!            STOP
!         end if

      case ("header")
         in_header = .true.
         hp => pseudo%header

         call get_value(attributes,"symbol",hp%symbol,status)
         if (status /= 0 ) call die("Cannot determine atomic symbol")

         call get_value(attributes,"atomic-number",value,status)
         if (status /= 0 ) call die("Cannot determine atomic-number")
         read(unit=value,fmt=*) hp%atomicnumber

         call get_value(attributes,"zval",value,status)
         if (status /= 0 ) call die("Cannot determine zval")
         read(unit=value,fmt=*) hp%zval

         call get_value(attributes,"xc-functional-type", &
                        hp%xcfunctionaltype,status)
         if (status /= 0 ) &
            call die("Cannot determine xc-functional-type")

         call get_value(attributes,"xc-functional-parametrization", &
                        hp%xcfunctionalparametrization,status)
         if (status /= 0 ) &
            call die("Cannot determine xc-functional-parametrization ")

         call get_value(attributes,"creator",hp%creator,status)
         if (status /= 0 ) hp%creator="unknown"

         call get_value(attributes,"date",hp%date,status)
         if (status /= 0 ) hp%date="unknown"

         call get_value(attributes,"flavor",hp%flavor,status)
         if (status /= 0 ) hp%flavor="unknown"

         call get_value(attributes,"relativistic",value,status)
         if (status /= 0 ) value = "no"
         hp%relativistic = (value == "yes")

         call get_value(attributes,"polarized",value,status)
         if (status /= 0 ) value = "no"
         hp%polarized = (value == "yes")

         call get_value(attributes,"core-corrections", &
                                    hp%core_corrections,status)
         if (status /= 0 ) hp%core_corrections = "nc"

      case ("vps")
         in_vps = .true.

         pseudo%npots = pseudo%npots + 1
         pp => pseudo%pot(pseudo%npots)
         rp => pp%V                       ! Pointer to radial function

         call get_value(attributes,"l",pp%l,status)
         if (status /= 0 ) call die("Cannot determine l for Vps")

         call get_value(attributes,"principal-n",value,status)
         if (status /= 0 ) call die("Cannot determine n for Vps")
         read(unit=value,fmt=*) pp%n

         call get_value(attributes,"cutoff",value,status)
         if (status /= 0 ) call die("Cannot determine cutoff for Vps")
         read(unit=value,fmt=*) pp%cutoff

         call get_value(attributes,"occupation",value,status)
         if (status /= 0 ) call die("Cannot determine occupation for Vps")
         read(unit=value,fmt=*) pp%occupation

         call get_value(attributes,"spin",value,status)
         if (status /= 0 ) call die("Cannot determine spin for Vps")
         read(unit=value,fmt=*) pp%spin

      case ("grid")

         call get_value(attributes,"type",grid%type,status)
         if (status /= 0 ) call die("Cannot determine grid type")

         call get_value(attributes,"npts",value,status)
         if (status /= 0 ) call die("Cannot determine grid npts")
         read(unit=value,fmt=*) grid%npts

         call get_value(attributes,"scale",value,status)
         if (status /= 0 ) call die("Cannot determine grid scale")
         read(unit=value,fmt=*) grid%scale

         call get_value(attributes,"step",value,status)
         if (status /= 0 ) call die("Cannot determine grid step")
         read(unit=value,fmt=*) grid%step

         !
         ! In this way we allow for a private grid for each radfunc,
         ! or for a global grid specification
         !
         if (in_radfunc) then
            rp%grid = grid
         else
            global_grid = grid
         end if

      case ("data")
         in_data = .true.
         if (rp%grid%npts == 0) STOP "Grid not specified correctly"
         allocate(rp%data(rp%grid%npts))
         ndata = 0             ! To start the build up

      case ("radfunc")
         in_radfunc = .true.
         rp%grid = global_grid     ! Might be empty
                                   ! There should then be a local grid element
                                   ! read later

      case ("pseudocore-charge")
         in_coreCharge = .true.
         rp => pseudo%core_charge

      case ("valence-charge")
         in_valenceCharge = .true.
         rp => pseudo%valence_charge

      case ("semilocal")
         in_semilocal = .true.

         call get_value(attributes,"npots-down",value,status)
         if (status /= 0 ) call die("Cannot determine npots-down")
         read(unit=value,fmt=*) pseudo%npots_down

         call get_value(attributes,"npots-up",value,status)
         if (status /= 0 ) call die("Cannot determine npots-up")
         read(unit=value,fmt=*) pseudo%npots_up

      case ("pseudowave-functions")
         in_pseudowavefun = .true.

      case ("pswf")
         in_pswf = .true.

         pseudo%npswfs = pseudo%npswfs + 1

         pw => pseudo%pswf(pseudo%npswfs)
         rp => pw%V                       ! Pointer to radial function

         call get_value(attributes,"l",pw%l,status)
         if (status /= 0 ) call die("Cannot determine l for Vps")

         call get_value(attributes,"principal-n",value,status)
         if (status /= 0 ) call die("Cannot determine n for Vps")
         read(unit=value,fmt=*) pw%n

         call get_value(attributes,"spin",value,status)
         if (status /= 0 ) call die("Cannot determine spin for Vps")
         read(unit=value,fmt=*) pw%spin

end select

end subroutine begin_element
!----------------------------------------------------------------------

subroutine end_element(name)


character(len=*), intent(in)     :: name

select case(name)

      case ("vps")
         in_vps = .false.

      case ("radfunc")
         in_radfunc = .false.

      case ("data")
      !
      ! We are done filling up the radfunc data
      ! Check that we got the advertised number of items
      !
         in_data = .false.
         if (ndata /= size(rp%data)) STOP "npts mismatch"

      case ("pseudocore-charge")
         in_coreCharge = .false.

      case ("valence-charge")
         in_valenceCharge = .false.

      case ("semilocal")
         in_semilocal = .false.

      case ("pseudowave-functions")
         in_pseudowavefun = .false.

      case ("pswf")
         in_pswf = .false.

      case ("pseudo")
!         call dump_pseudo(pseudo)

end select

end subroutine end_element
!----------------------------------------------------------------------

subroutine pcdata_chunk(chunk)


character(len=*), intent(in) :: chunk


if (len_trim(chunk) == 0) RETURN     ! skip empty chunk

if (in_data) then
!
! Note that we know where we need to put it through the pointer rp...
!
      call build_data_array(chunk,rp%data,ndata)

else if (in_header) then
      !
      ! There should not be any pcdata in header in this version...

      print *, "Header data:"
      print *, trim(chunk)

end if

end subroutine pcdata_chunk
!----------------------------------------------------------------------

subroutine die(str)


      character(len=*), intent(in), optional   :: str
      if (present(str)) then
         write(unit=0,fmt="(a)") trim(str)
      end if
      write(unit=0,fmt="(a)") "Stopping Program"
      stop
end subroutine die

#endif

end module m_pseudo
!!***
