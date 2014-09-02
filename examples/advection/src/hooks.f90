!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

module hooks
  use pf_mod_dtype
  implicit none
contains

  subroutine echo_error(pf, level)
    use iso_c_binding
    use feval, only: exact
    type(pf_pfasst), intent(inout) :: pf
    type(pf_level),  intent(inout) :: level

    ! real(c_double) :: yexact(level%nvars)
    ! real(pfdp), pointer :: qend(:)

    ! qend => array1(level%qend)

    ! call exact(state%t0+state%dt, yexact)
    ! print '("error: step: ",i3.3," iter: ",i4.3," error: ",es14.7)', &
    !      state%step+1, state%iter, maxval(abs(qend-yexact))
  end subroutine echo_error

  subroutine echo_residual(pf, level)
    use iso_c_binding
    use pf_mod_utils
    type(pf_pfasst), intent(inout) :: pf
    type(pf_level),  intent(inout) :: level

    ! real(pfdp), pointer :: r(:)

    ! r => array1(level%R(level%nnodes-1))

    ! print '("resid: step: ",i3.3," iter: ",i4.3," level: ",i2.2," resid: ",es14.7)', &
    !      state%step+1, state%iter, level%level, maxval(abs(r))
  end subroutine echo_residual

end module hooks
