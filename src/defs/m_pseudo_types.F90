!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pseudo_types
!! NAME
!! m_pseudo_types
!!
!! FUNCTION
!!  (to be completed)
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

module m_pseudo_types

implicit none
!
! Data structures for a prototype pseudopotential
!
integer, parameter, private    :: MAXN_POTS = 8
integer, parameter, private    :: dp = selected_real_kind(14)
!
public  :: dump_pseudo
!
!-----------------------------------------------------------
type, public :: grid_t
!
!     It should be possible to represent both log and linear
!     grids with a few parameters here.
!
      character(len=20)              :: type
      real(kind=dp)                  :: scale
      real(kind=dp)                  :: step
      integer                        :: npts
end type grid_t
!
type, public :: radfunc_t
      type(grid_t)                            :: grid
      real(kind=dp), dimension(:), pointer    :: data
end type radfunc_t

type, public :: vps_t
      character(len=1)               :: l
      integer                        :: n
      integer                        :: spin
      real(kind=dp)                  :: occupation
      real(kind=dp)                  :: cutoff
      type(radfunc_t)                :: V
end type vps_t

type, public :: pswf_t
      character(len=1)               :: l
      integer                        :: n
      integer                        :: spin
      type(radfunc_t)                :: V
end type pswf_t

type, public :: header_t
        character(len=2)        :: symbol
        real(kind=dp)           :: atomicnumber
        real(kind=dp)           :: zval
        character(len=10)       :: creator
        character(len=10)       :: date
        character(len=40)       :: flavor
        logical                 :: relativistic
        logical                 :: polarized
        character(len=30)       :: xcfunctionaltype
        character(len=30)       :: xcfunctionalparametrization
        character(len=4)        :: core_corrections
end type header_t

type, public :: pseudo_t
      type(header_t)                     :: header
      integer                            :: npots
      integer                            :: npswfs
      integer                            :: npots_down
      integer                            :: npots_up
      type(vps_t), dimension(MAXN_POTS)  :: pot
      type(pswf_t), dimension(MAXN_POTS) :: pswf
      type(radfunc_t)                    :: core_charge
      type(radfunc_t)                    :: valence_charge
end type pseudo_t


CONTAINS !===============================================

subroutine dump_pseudo(pseudo)


type(pseudo_t), intent(in), target   :: pseudo

integer  :: i
type(vps_t), pointer :: pp
type(pswf_t), pointer :: pw
type(radfunc_t), pointer :: rp

print *, "---PSEUDO data:"

print *, pseudo%header%xcfunctionaltype
do i = 1, pseudo%npots
      pp =>  pseudo%pot(i)
      rp =>  pseudo%pot(i)%V
      print *, "VPS ", i, " angular momentum: ", pp%l
      print *, "                 n: ", pp%n
      print *, "                 occupation: ", pp%occupation
      print *, "                 cutoff: ", pp%cutoff
      print *, "                 spin: ", pp%spin
      print *, "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)
end do
do i = 1, pseudo%npswfs
      pw =>  pseudo%pswf(i)
      rp =>  pseudo%pswf(i)%V
      print *, "PSWF ", i, " angular momentum: ", pw%l
      print *, "                 n: ", pw%n
      print *, "                 spin: ", pw%spin
      print *, "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)
end do
rp => pseudo%valence_charge
print *, "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)
rp => pseudo%core_charge
print *, "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)

end subroutine dump_pseudo

end module m_pseudo_types
!!***
