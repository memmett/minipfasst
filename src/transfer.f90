!
! Copyright (c) 2014, Matthew Emmett and Michael Minion.
!

! Dummy transfer (interpolate, restrict) routines.

module transfer
  use pf_mod_dtype
  implicit none
contains
  subroutine interpolate(qF, qG, levelF, levelG, t)
    real(pfdp),     intent(inout) :: qF(:)
    real(pfdp),     intent(in   ) :: qG(:), t
    type(pf_level), intent(in   ) :: levelF, levelG
    qF = qG
  end subroutine interpolate
  subroutine restrict(qF, qG, levelF, levelG, t)
    real(pfdp),     intent(in   ) :: qF(:), t
    real(pfdp),     intent(inout) :: qG(:)
    type(pf_level), intent(in   ) :: levelF, levelG
    qG = qF
  end subroutine restrict
end module transfer
