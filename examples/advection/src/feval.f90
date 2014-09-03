!
! Copyright (c) 2014, Matthew Emmett and Michael Minion.  All rights reserved.
!

! RHS routines for advection/diffusion example.

module feval
  use pf_mod_dtype
  use user
  implicit none
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

    complex(pfdp),   pointer :: wk(:)

    wk => lev%user%wk

    wk = y
    call fftw_execute_dft(lev%user%ffft, wk, wk)
    wk = -v * lev%user%ddx * wk / size(wk)
    call fftw_execute_dft(lev%user%ifft, wk, wk)

    f1 = real(wk)

  end subroutine f1eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the implicit function at y, t.
  subroutine f2eval(y, t, lev, f2)
    real(pfdp),     intent(in   ) :: y(:), t
    type(pf_level), intent(in   ) :: lev
    real(pfdp),     intent(  out) :: f2(:)

    complex(pfdp),   pointer :: wk(:)

    wk => lev%user%wk

    wk = y
    call fftw_execute_dft(lev%user%ffft, wk, wk)
    wk = nu * lev%user%lap * wk / size(wk)
    call fftw_execute_dft(lev%user%ifft, wk, wk)

    f2 = real(wk)
  end subroutine f2eval

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solve for y and return f2 also.
  subroutine f2comp(y, t, dt, rhs, lev, f2)
    real(pfdp),     intent(inout) :: y(:), f2(:)
    real(pfdp),     intent(in   ) :: rhs(:)
    real(pfdp),     intent(in   ) :: t, dt
    type(pf_level), intent(in   ) :: lev

    complex(pfdp),   pointer :: wk(:)

    wk => lev%user%wk

    wk = rhs
    call fftw_execute_dft(lev%user%ffft, wk, wk)
    wk = wk / (1.0_pfdp - nu*dt*lev%user%lap) / size(wk)
    call fftw_execute_dft(lev%user%ifft, wk, wk)

    y  = real(wk)
    f2 = (y - rhs) / dt
  end subroutine f2comp

end module feval
