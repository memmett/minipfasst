!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! RHS routines for advection/diffusion example.

module feval
  use pf_mod_dtype
  use iso_c_binding
  implicit none
  include 'fftw3.f03'

  real(pfdp), parameter :: &
       Lx     = 1.0_pfdp, &        ! domain size
       v      = 1.0_pfdp, &        ! velocity
       nu     = 0.02_pfdp, &       ! viscosity
       t00    = 0.15_pfdp           ! initial time for exact solution

  real(pfdp), parameter :: pi = 3.141592653589793_pfdp
  real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp

  type :: workspace
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:)              ! work space
     complex(pfdp), pointer :: ddx(:), lap(:)     ! operators
  end type workspace

  type(workspace), allocatable :: workspaces(:)

contains

  subroutine feval_create_workspaces(ndofs)
    integer, intent(in   ) :: ndofs(:)

    integer    :: nlevs, l, i
    real(pfdp) :: kx

    type(c_ptr) :: wk

    nlevs = size(ndofs)
    allocate(workspaces(nlevs))

    do l = 1, nlevs

       ! create in-place, complex fft plans
       wk = fftw_alloc_complex(int(ndofs(l), c_size_t))
       call c_f_pointer(wk, workspaces(l)%wk, [ndofs(l)])

       workspaces(l)%ffft = fftw_plan_dft_1d(ndofs(l), &
            workspaces(l)%wk, workspaces(l)%wk, FFTW_FORWARD, FFTW_ESTIMATE)
       workspaces(l)%ifft = fftw_plan_dft_1d(ndofs(l), &
            workspaces(l)%wk, workspaces(l)%wk, FFTW_BACKWARD, FFTW_ESTIMATE)

       ! create operators
       allocate(workspaces(l)%ddx(ndofs(l)))
       allocate(workspaces(l)%lap(ndofs(l)))
       do i = 1, ndofs(l)
          if (i <= ndofs(l)/2+1) then
             kx = two_pi / Lx * dble(i-1)
          else
             kx = two_pi / Lx * dble(-ndofs(l) + i - 1)
          end if

          workspaces(l)%ddx(i) = (0.0_pfdp, 1.0_pfdp) * kx

          if (kx**2 < 1e-13) then
             workspaces(l)%lap(i) = 0.0_pfdp
          else
             workspaces(l)%lap(i) = -kx**2
          end if
       end do
    end do

  end subroutine feval_create_workspaces

  subroutine feval_destroy_workspaces()
    integer :: l

    do l = 1, size(workspaces)
       deallocate(workspaces(l)%wk)
       deallocate(workspaces(l)%ddx)
       deallocate(workspaces(l)%lap)
       call fftw_destroy_plan(workspaces(l)%ffft)
       call fftw_destroy_plan(workspaces(l)%ifft)
    end do

    deallocate(workspaces)
  end subroutine feval_destroy_workspaces

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set initial condition.
  subroutine initial(q0)
    real(pfdp), intent(inout) :: q0(:)
    call exact(0.0_pfdp, q0)
  end subroutine initial

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exact(t, yex)
    real(pfdp), intent(in   ) :: t
    real(pfdp), intent(  out) :: yex(:)

    integer    :: nvars, i, ii, nbox
    real(pfdp) :: tol, x

    nvars = size(yex)
    yex   = 0

    if (nu .gt. 0) then
       do i = 1, nvars
          x = Lx*dble(i-nvars/2-1)/dble(nvars) - t*v
          yex(i) = yex(i) + dcos(2.0_pfdp*pi*x)*dexp(-4.0_pfdp*pi*pi*nu*t)
       end do
    else

       ! decide how many images so that contribution is neglible
       tol  = 1e-16
       nbox = 1 + ceiling( sqrt( -(4.0*t00)*log((4.0*pi*(t00))**(0.5)*tol) ))

       do ii = -nbox, nbox
          do i = 1, nvars
             x = Lx*dble(i-nvars/2-1)/dble(nvars) + ii*Lx - t*v
             yex(i) = yex(i) + 1.0/(4.0*pi*t00)**(0.5)*dexp(-x**2/(4.0*t00))
          end do
       end do
    end if


  end subroutine exact

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the explicit function at y, t.
  subroutine f1eval(y, t, level, f1)
    real(pfdp),     intent(in   ) :: y(:), t
    type(pf_level), intent(in   ) :: level
    real(pfdp),     intent(  out) :: f1(:)

    complex(pfdp),   pointer :: wk(:)

    wk => workspaces(level%level)%wk

    wk = y
    call fftw_execute_dft(workspaces(level%level)%ffft, wk, wk)
    wk = -v * workspaces(level%level)%ddx * wk / size(wk)
    call fftw_execute_dft(workspaces(level%level)%ifft, wk, wk)

    f1 = real(wk)

  end subroutine f1eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the implicit function at y, t.
  subroutine f2eval(y, t, level, f2)
    real(pfdp),     intent(in   ) :: y(:), t
    type(pf_level), intent(in   ) :: level
    real(pfdp),     intent(  out) :: f2(:)

    complex(pfdp),   pointer :: wk(:)

    wk => workspaces(level%level)%wk

    wk = y
    call fftw_execute_dft(workspaces(level%level)%ffft, wk, wk)
    wk = nu * workspaces(level%level)%lap * wk / size(wk)
    call fftw_execute_dft(workspaces(level%level)%ifft, wk, wk)

    f2 = real(wk)
  end subroutine f2eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solve for y and return f2 also.
  subroutine f2comp(y, t, dt, rhs, level, f2)
    real(pfdp),     intent(inout) :: y(:), f2(:)
    real(pfdp),     intent(in   ) :: rhs(:)
    real(pfdp),     intent(in   ) :: t, dt
    type(pf_level), intent(in   ) :: level

    complex(pfdp),   pointer :: wk(:)

    wk => workspaces(level%level)%wk

    wk = rhs
    call fftw_execute_dft(workspaces(level%level)%ffft, wk, wk)
    wk = wk / (1.0_pfdp - nu*dt*workspaces(level%level)%lap) / size(wk)
    call fftw_execute_dft(workspaces(level%level)%ifft, wk, wk)

    y  = real(wk)
    f2 = (y - rhs) / dt
  end subroutine f2comp

end module feval
