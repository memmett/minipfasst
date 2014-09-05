!
! Copyright (c) 2014, Matthew Emmett and Michael Minion.  All rights reserved.
!

! Transfer (interpolate, restrict) routines.

module transfer
  use pf_mod_dtype
  implicit none
contains

  subroutine interpolate(qF, qG, fine, crse, t)
    real(pfdp),     intent(in   ) :: t, qG(:)
    real(pfdp),     intent(  out) :: qF(:)
    type(pf_level), intent(in   ) :: fine, crse
    stop
  end subroutine interpolate

  subroutine restrict(qF, qG, fine, crse, t)
    real(pfdp),     intent(in   ) :: t, qF(:)
    real(pfdp),     intent(  out) :: qG(:)
    type(pf_level), intent(in   ) :: fine, crse
    stop
  end subroutine restrict
end module transfer
