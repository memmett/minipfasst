!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
!   1. Redistributions of source code must retain the above copyright
!      notice, this list of conditions and the following disclaimer.
!
!   2. Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions and the following disclaimer in
!      the documentation and/or other materials provided with the
!      distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
