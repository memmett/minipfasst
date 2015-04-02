!
! Copyright (c) 2014, Matthew Emmett and Michael Minion.  All rights reserved.
!

! RHS routines for advection/diffusion example.

module feval
  use pf_mod_dtype
  use multigrid, only: lx, spatial_order, fill_bc, multigrid_v_cycle
  implicit none

  real(pfdp), parameter :: &
       v      = 1.0_pfdp, &        ! velocity
       nu     = 0.02_pfdp, &       ! viscosity
       t00    = 0.15_pfdp          ! initial time for exact solution

  real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp
  real(pfdp), parameter :: pi = 3.141592653589793_pfdp

contains

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
  subroutine f1eval(y, t, lev, f1)
    real(pfdp),     intent(in   ) :: y(:), t
    type(pf_level), intent(in   ) :: lev
    real(pfdp),     intent(  out) :: f1(:)

    integer :: nx
    real(pfdp) :: dx, ybc(1-spatial_order:size(y)+spatial_order), y_x(size(y))

    nx = size(y)
    dx = lx/dble(nx)

    call fill_bc(y, ybc, nx, spatial_order)
    y_x = ( ybc(2:nx+1) - ybc(0:nx-1) ) / (2*dx)
    f1 = -v * y_x
  end subroutine f1eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the implicit function at y, t.
  subroutine f2eval(y, t, lev, f2)
    real(pfdp),     intent(in   ) :: y(:), t
    type(pf_level), intent(in   ) :: lev
    real(pfdp),     intent(  out) :: f2(:)

    integer :: nx
    real(pfdp) :: dx, ybc(1-spatial_order:size(y)+spatial_order), lap(size(y))

    nx = size(y)
    dx = lx/dble(nx)

    call fill_bc(y, ybc, nx, spatial_order)
    lap = ( ybc(0:nx-1) - 2*ybc(1:nx) + ybc(2:nx+1) ) / (dx*dx)
    f2 = nu*lap
  end subroutine f2eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solve for y and return f2 also.
  subroutine f2comp(y, t, dt, rhs, lev, f2)
    real(pfdp),     intent(inout) :: y(:), f2(:)
    real(pfdp),     intent(in   ) :: rhs(:)
    real(pfdp),     intent(in   ) :: t, dt
    type(pf_level), intent(in   ) :: lev

    real(pfdp) :: maxresidual
    integer :: k

    do k = 1, 2
       call multigrid_v_cycle(y, t, nu*dt, rhs, 1, maxresidual)
       print *, k, maxresidual
    end do

    f2 = (y - rhs) / dt
  end subroutine f2comp

end module feval
