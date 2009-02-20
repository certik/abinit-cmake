subroutine wrtout_myproc(unit,message)

  use defs_basis

  implicit none

  !Arguments ------------------------------------
  integer,intent(in) :: unit
  character(len=500),intent(inout) :: message

  write(unit, *) "[AB] ", trim(message)
end subroutine wrtout_myproc

subroutine leave_myproc()

  implicit none

  stop
end subroutine leave_myproc

subroutine timab(nn,option,tottim)

  use defs_basis
  use defs_time

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: nn,option
  !arrays
  real(dp),intent(out) :: tottim(2)

end subroutine timab
