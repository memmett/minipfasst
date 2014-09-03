!
! Copyright (C) 2014 Matthew Emmett and Michael Minion.
!

!
! User workspace for advection/diffusion example.
!

module user
  use iso_c_binding
  use pf_mod_dtype
  implicit none

  include 'fftw3.f03'

  real(pfdp), parameter :: &
       Lx     = 1.0_pfdp, &        ! domain size
       v      = 1.0_pfdp, &        ! velocity
       nu     = 0.02_pfdp, &       ! viscosity
       t00    = 0.15_pfdp           ! initial time for exact solution

  real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp
  real(pfdp), parameter :: pi = 3.141592653589793_pfdp

contains

  subroutine user_setup(lev)
    type(pf_level), intent(inout) :: lev

    integer :: i
    real(pfdp) :: kx
    type(c_ptr) :: wk

    ! create in-place, complex fft plans
    wk = fftw_alloc_complex(int(lev%ndofs, c_size_t))
    call c_f_pointer(wk, lev%user%wk, [lev%ndofs])

    lev%user%ffft = fftw_plan_dft_1d(lev%ndofs, &
         lev%user%wk, lev%user%wk, FFTW_FORWARD, FFTW_ESTIMATE)
    lev%user%ifft = fftw_plan_dft_1d(lev%ndofs, &
         lev%user%wk, lev%user%wk, FFTW_BACKWARD, FFTW_ESTIMATE)

    ! create operators
    allocate(lev%user%ddx(lev%ndofs))
    allocate(lev%user%lap(lev%ndofs))
    do i = 1, lev%ndofs
       if (i <= lev%ndofs/2+1) then
          kx = two_pi / Lx * dble(i-1)
       else
          kx = two_pi / Lx * dble(-lev%ndofs + i - 1)
       end if

       lev%user%ddx(i) = (0.0_pfdp, 1.0_pfdp) * kx

       if (kx**2 < 1e-13) then
          lev%user%lap(i) = 0.0_pfdp
       else
          lev%user%lap(i) = -kx**2
       end if
    end do
  end subroutine user_setup

  subroutine user_destroy(lev)
    type(pf_level), intent(inout) :: lev
    deallocate(lev%user%wk)
    deallocate(lev%user%ddx)
    deallocate(lev%user%lap)
    call fftw_destroy_plan(lev%user%ffft)
    call fftw_destroy_plan(lev%user%ifft)
  end subroutine user_destroy

end module user
