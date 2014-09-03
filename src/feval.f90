!
! Copyright (c) 2014, Matthew Emmett and Michael Minion.
!

! Dummy RHS routines

module feval
  use pf_mod_dtype
  implicit none
contains

  ! Evaluate the explicit function at y, t.
  subroutine f1eval(y, t, level, f1)
    real(pfdp),     intent(in   ) :: y(:), t
    type(pf_level), intent(in   ) :: level
    real(pfdp),     intent(  out) :: f1(:)
    f1 = 0
  end subroutine f1eval

  ! Evaluate the implicit function at y, t.
  subroutine f2eval(y, t, level, f2)
    real(pfdp),     intent(in   ) :: y(:), t
    type(pf_level), intent(in   ) :: level
    real(pfdp),     intent(  out) :: f2(:)
    f2 = 0
  end subroutine f2eval

  ! Solve for y and return f2 also.
  subroutine f2comp(y, t, dt, rhs, level, f2)
    real(pfdp),     intent(inout) :: y(:), f2(:)
    real(pfdp),     intent(in   ) :: rhs(:)
    real(pfdp),     intent(in   ) :: t, dt
    type(pf_level), intent(in   ) :: level
    y  = rhs
    f2 = 0
  end subroutine f2comp

end module feval
