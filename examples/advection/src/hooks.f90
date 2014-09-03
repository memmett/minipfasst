!
! Copyright (c) 2014, Matthew Emmett and Michael Minion.  All rights reserved.
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

    real(c_double) :: yexact(level%ndofs)

    call exact(pf%t0+pf%dt, yexact)
    print '("error: step: ",i3.3," iter: ",i4.3," error: ",es14.7)', &
         pf%step+1, pf%iter, maxval(abs(level%qend-yexact))
  end subroutine echo_error

end module hooks
