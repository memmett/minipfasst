!
! Copyright (c) 2014, Matthew Emmett and Michael Minion.  All rights reserved.
!

module hooks
  use pf_mod_dtype
  implicit none
contains

  subroutine echo_error(pf, level)
    use feval, only: shapiro
    type(pf_pfasst), intent(inout) :: pf
    type(pf_level),  intent(inout) :: level

  end subroutine echo_error

end module hooks
