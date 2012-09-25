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

! RHS routines for advection/diffusion example.

module feval
  use iso_c_binding
  implicit none
  include 'fftw3.f03'

  real(kind=8), parameter :: &
       Lx     = 1.0d0, &        ! domain size
       v      = 1.0d0, &        ! velocity
       nu     = 0.02d0, &       ! viscosity
       t00    = 1.0d0           ! initial time for exact solution

  real(kind=8), parameter :: pi = 3.141592653589793d0
  real(kind=8), parameter :: two_pi = 6.2831853071795862d0

  type :: ad_level_t
     type(c_ptr) :: ffft, ifft
     complex(kind=8), pointer :: wk(:)              ! work space
     complex(kind=8), pointer :: ddx(:), lap(:)     ! operators
  end type ad_level_t

  type(ad_level_t), pointer :: levels(:)

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine feval_init(nlevels, nvars)

    integer, intent(in) :: nlevels, nvars(nlevels)

    integer :: i, l
    type(c_ptr) :: wk
    real(kind=8) :: kx

    allocate(levels(nlevels))
    do l = 1, nlevels

       ! create in-place, complex fft plans
       wk = fftw_alloc_complex(int(nvars(l), c_size_t))
       call c_f_pointer(wk, levels(l)%wk, [nvars(l)])

       levels(l)%ffft = fftw_plan_dft_1d(nvars(l), &
            levels(l)%wk, levels(l)%wk, FFTW_FORWARD, FFTW_ESTIMATE)
       levels(l)%ifft = fftw_plan_dft_1d(nvars(l), &
            levels(l)%wk, levels(l)%wk, FFTW_BACKWARD, FFTW_ESTIMATE)

       ! create operators
       allocate(levels(l)%ddx(nvars(l)))
       allocate(levels(l)%lap(nvars(l)))
       do i = 1, nvars(l)
          if (i <= nvars(l)/2+1) then
             kx = two_pi / Lx * dble(i-1)
          else
             kx = two_pi / Lx * dble(-nvars(l) + i - 1)
          end if

          levels(l)%ddx(i) = (0.0d0, 1.0d0) * kx

          if (kx**2 < 1e-13) then
             levels(l)%lap(i) = 0.0d0
          else
             levels(l)%lap(i) = -kx**2
          end if
       end do
    end do
  end subroutine feval_init

  subroutine feval_finalize()

    integer :: l

    do l = 1, size(levels)
       deallocate(levels(l)%wk)
       deallocate(levels(l)%ddx)
       deallocate(levels(l)%lap)
    end do

    deallocate(levels)

  end subroutine feval_finalize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set initial condition.
  subroutine initial(q0)
    double precision, intent(inout) :: q0(:)

    call exact(0.0d0, q0)
  end subroutine initial

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exact(t, yex)
    double precision, intent(in)  :: t
    double precision, intent(out) :: yex(:)

    integer :: nvars, i, ii
    double precision :: x

    nvars = size(yex)
    yex   = 0.0d0

    do ii = -2, 2

       do i = 1, nvars
          x = Lx*dble(i-nvars/2-1)/dble(nvars) + ii*Lx - t*v
          yex(i) = yex(i) + 1.0/(4.0*pi*nu*(t+t00))**(0.5)*dexp(-x**2/(4.0*nu*(t+t00)))
       end do

    end do

  end subroutine exact

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the explicit function at y, t.
  subroutine eval_f1(y, t, level, f1)

    double precision, intent(in)  :: y(:), t
    integer,          intent(in)  :: level
    double precision, intent(out) :: f1(:)

    complex(kind=8), pointer :: wk(:)

    wk => levels(level)%wk

    wk = y
    call fftw_execute_dft(levels(level)%ffft, wk, wk)
    wk = -v * levels(level)%ddx * wk / size(wk)
    call fftw_execute_dft(levels(level)%ifft, wk, wk)

    f1 = real(wk)
  end subroutine eval_f1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the implicit function at y, t.
  subroutine eval_f2(y, t, level, f2)

    double precision, intent(in)  :: y(:), t
    integer,          intent(in)  :: level
    double precision, intent(out) :: f2(:)

    complex(kind=8), pointer :: wk(:)

    wk => levels(level)%wk

    wk = y
    call fftw_execute_dft(levels(level)%ffft, wk, wk)
    wk = nu * levels(level)%lap * wk / size(wk)
    call fftw_execute_dft(levels(level)%ifft, wk, wk)

    f2 = real(wk)
  end subroutine eval_f2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solve for y and return f2 also.
  subroutine comp_f2(y, t, dt, rhs, level, f2)
    double precision, intent(inout) :: y(:), f2(:)
    double precision, intent(in)    :: rhs(:)
    double precision, intent(in)    :: t, dt
    integer,          intent(in)    :: level

    complex(kind=8), pointer :: wk(:)

    wk => levels(level)%wk

    wk = rhs
    call fftw_execute_dft(levels(level)%ffft, wk, wk)
    wk = wk / (1.0d0 - nu*dt*levels(level)%lap) / size(wk)
    call fftw_execute_dft(levels(level)%ifft, wk, wk)

    y  = real(wk)
    f2 = (y - rhs) / dt
  end subroutine comp_f2

end module feval
