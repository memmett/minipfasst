!
! Copyright (C) 2014 Matthew Emmett and Michael Minion.
!

!
! Dummy user routines.
!

module user
  use pf_mod_dtype
  implicit none
contains

  subroutine user_setup(lev)
    type(pf_level), intent(inout) :: lev
  end subroutine user_setup

  subroutine user_destroy(lev)
    type(pf_level), intent(inout) :: lev
  end subroutine user_destroy
end module user
