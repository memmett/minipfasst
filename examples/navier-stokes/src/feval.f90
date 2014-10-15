!
! Copyright (c) 2014, Matthew Emmett and Michael Minion.  All rights reserved.
!

! See the README.org file for details.

module feval
  use pf_mod_dtype
  implicit none

contains

  subroutine f1eval(y, t, lev, f)
    real(pfdp),     intent(in   ) :: y(:), t
    type(pf_level), intent(inout) :: lev
    real(pfdp),     intent(  out) :: f(:)

    integer :: i, j, k, nx, ny, nz

    complex(pfdp), dimension(lev%user%nx,lev%user%ny,lev%user%nz) :: &
         u, v, w, uhat, vhat, what, wk1, wk2, wk3, nluhat, nlvhat, nlwhat

    complex(pfdp) :: kx(lev%user%nx), ky(lev%user%ny), kz(lev%user%nz), phat

    real(pfdp) :: k2, nu, scale, diffop

    ! shortcuts
    nx = lev%user%nx; ny = lev%user%ny; nz = lev%user%nz
    kx = lev%user%kx; ky = lev%user%ky; kz = lev%user%kz
    nu = lev%user%nu; scale = lev%user%scale

    ! unpack
    call unpack3(y, uhat, vhat, what)

    call dfftw_execute_dft_(lev%user%bft, uhat, u)
    call dfftw_execute_dft_(lev%user%bft, vhat, v)
    call dfftw_execute_dft_(lev%user%bft, what, w)

    do k=1,nz ; do j=1,ny ; do i=1,nx
       wk1(i,j,k) = uhat(i,j,k) * kx(i)
       wk2(i,j,k) = uhat(i,j,k) * ky(j)
       wk3(i,j,k) = uhat(i,j,k) * kz(k)
    end do; end do ; end do
    call dfftw_execute_dft_(lev%user%bft, wk1, wk1)
    call dfftw_execute_dft_(lev%user%bft, wk2, wk2)
    call dfftw_execute_dft_(lev%user%bft, wk3, wk3)
    nluhat = (u * wk1 + v * wk2 + w * wk3) * scale
    call dfftw_execute_dft_(lev%user%fft, nluhat, nluhat)

    do k=1,nz ; do j=1,ny ; do i=1,nx
       wk1(i,j,k) = vhat(i,j,k) * kx(i)
       wk2(i,j,k) = vhat(i,j,k) * ky(j)
       wk3(i,j,k) = vhat(i,j,k) * kz(k)
    end do; end do ; end do
    call dfftw_execute_dft_(lev%user%bft, wk1, wk1)
    call dfftw_execute_dft_(lev%user%bft, wk2, wk2)
    call dfftw_execute_dft_(lev%user%bft, wk3, wk3)
    nlvhat = (u * wk1 + v * wk2 + w * wk3) * scale
    call dfftw_execute_dft_(lev%user%fft, nlvhat, nlvhat)

    do k=1,nz ; do j=1,ny ; do i=1,nx
       wk1(i,j,k) = what(i,j,k) * kx(i)
       wk2(i,j,k) = what(i,j,k) * ky(j)
       wk3(i,j,k) = what(i,j,k) * kz(k)
    end do; end do ; end do
    call dfftw_execute_dft_(lev%user%bft, wk1, wk1)
    call dfftw_execute_dft_(lev%user%bft, wk2, wk2)
    call dfftw_execute_dft_(lev%user%bft, wk3, wk3)
    nlwhat = (u * wk1 + v * wk2 + w * wk3) * scale
    call dfftw_execute_dft_(lev%user%fft, nlwhat, nlwhat)

    ! evaluate
    do k=1,nz
       do j=1,ny
          do i=1,nx
             k2 = realpart(kx(i)*kx(i) + ky(j)*ky(j) + kz(k)*kz(k))
             diffop = nu * k2

             ! pressure
             phat = -( kx(i)*nluhat(i,j,k) + ky(j)*nlvhat(i,j,k) + kz(k)*nlwhat(i,j,k) ) / (k2 + 0.1d0**13)

             uhat(i,j,k) = - nluhat(i,j,k) - kx(i)*phat
             vhat(i,j,k) = - nlvhat(i,j,k) - ky(j)*phat
             what(i,j,k) = - nlwhat(i,j,k) - kz(k)*phat

          end do
       end do
    end do

    call pack3(f, uhat, vhat, what)
  end subroutine f1eval

  subroutine f2eval(y, t, lev, f)
    real(pfdp),     intent(in   ) :: y(:), t
    type(pf_level), intent(inout) :: lev
    real(pfdp),     intent(  out) :: f(:)

    integer :: i, j, k, nx, ny, nz

    complex(pfdp), dimension(lev%user%nx,lev%user%ny,lev%user%nz) :: &
         uhat, vhat, what

    complex(pfdp) :: kx(lev%user%nx), ky(lev%user%ny), kz(lev%user%nz)

    real(pfdp) :: k2, nu, scale, diffop

    ! shortcuts
    nx = lev%user%nx; ny = lev%user%ny; nz = lev%user%nz
    kx = lev%user%kx; ky = lev%user%ky; kz = lev%user%kz
    nu = lev%user%nu; scale = lev%user%scale

    ! unpack
    call unpack3(y, uhat, vhat, what)

    do k=1,nz
       do j=1,ny
          do i=1,nx
             k2 = realpart(kx(i)*kx(i) + ky(j)*ky(j) + kz(k)*kz(k))
             diffop = nu * k2

             uhat(i,j,k) = diffop * uhat(i,j,k)
             vhat(i,j,k) = diffop * vhat(i,j,k)
             what(i,j,k) = diffop * what(i,j,k)

          end do
       end do
    end do
    call pack3(f, uhat, vhat, what)
  end subroutine f2eval

  subroutine f2comp(y, t, dt, rhs, lev, f)
    real(pfdp),     intent(inout) :: y(:)
    real(pfdp),     intent(in   ) :: rhs(:), t, dt
    type(pf_level), intent(inout) :: lev
    real(pfdp),     intent(  out) :: f(:)

    integer :: i, j, k, nx, ny, nz

    complex(pfdp), dimension(lev%user%nx,lev%user%ny,lev%user%nz) :: &
         uhat, vhat, what, u0hat, v0hat, w0hat

    complex(pfdp) :: kx(lev%user%nx), ky(lev%user%ny), kz(lev%user%nz), invop

    real(pfdp) :: k2, nu, scale

    ! shortcuts
    nx = lev%user%nx; ny = lev%user%ny; nz = lev%user%nz
    kx = lev%user%kx; ky = lev%user%ky; kz = lev%user%kz
    nu = lev%user%nu; scale = lev%user%scale

    ! unpack
    call unpack3(rhs, u0hat, v0hat, w0hat)
    call unpack3(y, uhat, vhat, what)

    ! solve
    do k=1,nz
       do j=1,ny
          do i=1,nx
             k2 = realpart(kx(i)*kx(i) + ky(j)*ky(j) + kz(k)*kz(k))
             invop = 1.d0 / (1 - dt * nu * k2 )

             ! diffusion solve
             uhat(i,j,k) = invop * u0hat(i,j,k)
             vhat(i,j,k) = invop * v0hat(i,j,k)
             what(i,j,k) = invop * w0hat(i,j,k)
          end do
       end do
    end do

    call pack3(y, uhat, vhat, what)
    f = (y-rhs) / dt
  end subroutine f2comp

  ! Evaluate F(U).
  subroutine impl_eval(y, t, lev, f)
    real(pfdp),     intent(in   ) :: y(:), t
    type(pf_level), intent(inout) :: lev
    real(pfdp),     intent(  out) :: f(:)

    integer :: i, j, k, nx, ny, nz

    complex(pfdp), dimension(lev%user%nx,lev%user%ny,lev%user%nz) :: &
         u, v, w, uhat, vhat, what, wk1, wk2, wk3, nluhat, nlvhat, nlwhat

    complex(pfdp) :: kx(lev%user%nx), ky(lev%user%ny), kz(lev%user%nz), phat

    real(pfdp) :: k2, nu, scale, diffop

    ! shortcuts
    nx = lev%user%nx; ny = lev%user%ny; nz = lev%user%nz
    kx = lev%user%kx; ky = lev%user%ky; kz = lev%user%kz
    nu = lev%user%nu; scale = lev%user%scale

    ! unpack
    call unpack3(y, uhat, vhat, what)

    call dfftw_execute_dft_(lev%user%bft, uhat, u)
    call dfftw_execute_dft_(lev%user%bft, vhat, v)
    call dfftw_execute_dft_(lev%user%bft, what, w)

    do k=1,nz ; do j=1,ny ; do i=1,nx
       wk1(i,j,k) = uhat(i,j,k) * kx(i)
       wk2(i,j,k) = uhat(i,j,k) * ky(j)
       wk3(i,j,k) = uhat(i,j,k) * kz(k)
    end do; end do ; end do
    call dfftw_execute_dft_(lev%user%bft, wk1, wk1)
    call dfftw_execute_dft_(lev%user%bft, wk2, wk2)
    call dfftw_execute_dft_(lev%user%bft, wk3, wk3)
    nluhat = (u * wk1 + v * wk2 + w * wk3) * scale
    call dfftw_execute_dft_(lev%user%fft, nluhat, nluhat)

    do k=1,nz ; do j=1,ny ; do i=1,nx
       wk1(i,j,k) = vhat(i,j,k) * kx(i)
       wk2(i,j,k) = vhat(i,j,k) * ky(j)
       wk3(i,j,k) = vhat(i,j,k) * kz(k)
    end do; end do ; end do
    call dfftw_execute_dft_(lev%user%bft, wk1, wk1)
    call dfftw_execute_dft_(lev%user%bft, wk2, wk2)
    call dfftw_execute_dft_(lev%user%bft, wk3, wk3)
    nlvhat = (u * wk1 + v * wk2 + w * wk3) * scale
    call dfftw_execute_dft_(lev%user%fft, nlvhat, nlvhat)

    do k=1,nz ; do j=1,ny ; do i=1,nx
       wk1(i,j,k) = what(i,j,k) * kx(i)
       wk2(i,j,k) = what(i,j,k) * ky(j)
       wk3(i,j,k) = what(i,j,k) * kz(k)
    end do; end do ; end do
    call dfftw_execute_dft_(lev%user%bft, wk1, wk1)
    call dfftw_execute_dft_(lev%user%bft, wk2, wk2)
    call dfftw_execute_dft_(lev%user%bft, wk3, wk3)
    nlwhat = (u * wk1 + v * wk2 + w * wk3) * scale
    call dfftw_execute_dft_(lev%user%fft, nlwhat, nlwhat)

    ! evaluate
    do k=1,nz
       do j=1,ny
          do i=1,nx
             k2 = realpart(kx(i)*kx(i) + ky(j)*ky(j) + kz(k)*kz(k))
             diffop = nu * k2

             ! pressure
             phat = -( kx(i)*nluhat(i,j,k) + ky(j)*nlvhat(i,j,k) + kz(k)*nlwhat(i,j,k) ) / (k2 + 0.1d0**13)

             uhat(i,j,k) = - nluhat(i,j,k) - kx(i)*phat + diffop * uhat(i,j,k)
             vhat(i,j,k) = - nlvhat(i,j,k) - ky(j)*phat + diffop * vhat(i,j,k)
             what(i,j,k) = - nlwhat(i,j,k) - kz(k)*phat + diffop * what(i,j,k)

          end do
       end do
    end do

    call pack3(f, uhat, vhat, what)
  end subroutine impl_eval

  ! Solve U + a F(U) = RHS.
  subroutine impl_solve(y, t, a, rhs, lev, step)
    real(pfdp),     intent(inout) :: y(:)
    real(pfdp),     intent(in   ) :: rhs(:)
    real(pfdp),     intent(in   ) :: t, a
    type(pf_level), intent(inout) :: lev
    integer,        intent(in   ) :: step

    integer :: iter, i, j, k, nx, ny, nz

    complex(pfdp), dimension(lev%user%nx,lev%user%ny,lev%user%nz) :: &
         u, v, w, uhat, vhat, what, u0hat, v0hat, w0hat, wk1, wk2, wk3, &
         nluhat, nlvhat, nlwhat

    complex(pfdp) :: kx(lev%user%nx), ky(lev%user%ny), kz(lev%user%nz), invop, phat

    real(pfdp) :: k2, nu, residual, scale

    ! shortcuts
    nx = lev%user%nx; ny = lev%user%ny; nz = lev%user%nz
    kx = lev%user%kx; ky = lev%user%ky; kz = lev%user%kz
    nu = lev%user%nu; scale = lev%user%scale

    ! unpack
    call unpack3(rhs, u0hat, v0hat, w0hat)
    call unpack3(y, uhat, vhat, what)

    ! fixed point iteration
    iter = 1
    residual = 1
    do while (residual > lev%user%tol)

       call dfftw_execute_dft_(lev%user%bft, uhat, u)
       call dfftw_execute_dft_(lev%user%bft, vhat, v)
       call dfftw_execute_dft_(lev%user%bft, what, w)

       do k=1,nz ; do j=1,ny ; do i=1,nx
          wk1(i,j,k) = uhat(i,j,k) * kx(i)
          wk2(i,j,k) = uhat(i,j,k) * ky(j)
          wk3(i,j,k) = uhat(i,j,k) * kz(k)
       end do; end do ; end do
       call dfftw_execute_dft_(lev%user%bft, wk1, wk1)
       call dfftw_execute_dft_(lev%user%bft, wk2, wk2)
       call dfftw_execute_dft_(lev%user%bft, wk3, wk3)
       nluhat = (u * wk1 + v * wk2 + w * wk3) * scale
       call dfftw_execute_dft_(lev%user%fft, nluhat, nluhat)

       do k=1,nz ; do j=1,ny ; do i=1,nx
          wk1(i,j,k) = vhat(i,j,k) * kx(i)
          wk2(i,j,k) = vhat(i,j,k) * ky(j)
          wk3(i,j,k) = vhat(i,j,k) * kz(k)
       end do; end do ; end do
       call dfftw_execute_dft_(lev%user%bft, wk1, wk1)
       call dfftw_execute_dft_(lev%user%bft, wk2, wk2)
       call dfftw_execute_dft_(lev%user%bft, wk3, wk3)
       nlvhat = (u * wk1 + v * wk2 + w * wk3) * scale
       call dfftw_execute_dft_(lev%user%fft, nlvhat, nlvhat)

       do k=1,nz ; do j=1,ny ; do i=1,nx
          wk1(i,j,k) = what(i,j,k) * kx(i)
          wk2(i,j,k) = what(i,j,k) * ky(j)
          wk3(i,j,k) = what(i,j,k) * kz(k)
       end do; end do ; end do
       call dfftw_execute_dft_(lev%user%bft, wk1, wk1)
       call dfftw_execute_dft_(lev%user%bft, wk2, wk2)
       call dfftw_execute_dft_(lev%user%bft, wk3, wk3)
       nlwhat = (u * wk1 + v * wk2 + w * wk3) * scale
       call dfftw_execute_dft_(lev%user%fft, nlwhat, nlwhat)

       ! solve
       do k=1,nz
          do j=1,ny
             do i=1,nx
                k2 = realpart(kx(i)*kx(i) + ky(j)*ky(j) + kz(k)*kz(k))
                invop = 1.d0 / (1 + a * nu * k2 )

                ! save for residual
                wk1(i,j,k) = uhat(i,j,k)
                wk2(i,j,k) = vhat(i,j,k)
                wk3(i,j,k) = what(i,j,k)

                ! pressure
                phat = -( kx(i)*nluhat(i,j,k) + ky(j)*nlvhat(i,j,k) + kz(k)*nlwhat(i,j,k) ) / (k2 + 0.1d0**13)

                ! diffusion solve
                uhat(i,j,k) = invop * ( u0hat(i,j,k) + a * nluhat(i,j,k) + a * kx(i)*phat )
                vhat(i,j,k) = invop * ( v0hat(i,j,k) + a * nlvhat(i,j,k) + a * ky(j)*phat )
                what(i,j,k) = invop * ( w0hat(i,j,k) + a * nlwhat(i,j,k) + a * kz(k)*phat )
             end do
          end do
       end do

       residual = maxval(abs(wk1 - uhat)) + maxval(abs(wk2 - vhat)) + maxval(abs(wk3 - what))
       print '(a20,i6,i3,i3,es20.12,es8.1)', 'residual', step, lev%level, iter, residual, lev%user%tol
       iter = iter + 1
    end do

    call pack3(y, uhat, vhat, what)
  end subroutine impl_solve

  ! subroutine fft3(plan, wk, scale, u, v, w, uhat, vhat, what)
  !   integer(8),    intent(in   )                   :: plan
  !   real(pfdp),    intent(in   )                   :: scale
  !   complex(pfdp), intent(in   ), dimension(:,:,:) :: u, v, w
  !   complex(pfdp), intent(  out), dimension(:,:,:) :: wk, uhat, vhat, what
  !   wk = u * scale; call dfftw_execute_dft_(plan, wk, uhat)
  !   wk = v * scale; call dfftw_execute_dft_(plan, wk, vhat)
  !   wk = w * scale; call dfftw_execute_dft_(plan, wk, what)
  ! end subroutine fft3

  subroutine pack3(y, u, v, w)
    real(pfdp),    intent(  out)                   :: y(:)
    complex(pfdp), intent(in   ), dimension(:,:,:) :: u, v, w

    integer :: idx, i, j, k, nx, ny, nz

    nx = size(u, dim=1); ny = size(u, dim=2); nz = size(u, dim=3)

    if (size(y) /= 2*3*nx*ny*nz) stop "ERROR: Size mismatch in pack3."

    idx = 1
    do k = 1, nz; do j = 1, ny; do i = 1, nx
       y(idx) = realpart(u(i,j,k)); idx = idx + 1
       y(idx) = imagpart(u(i,j,k)); idx = idx + 1
       y(idx) = realpart(v(i,j,k)); idx = idx + 1
       y(idx) = imagpart(v(i,j,k)); idx = idx + 1
       y(idx) = realpart(w(i,j,k)); idx = idx + 1
       y(idx) = imagpart(w(i,j,k)); idx = idx + 1
    end do; end do; end do
  end subroutine pack3

  subroutine unpack3(y, u, v, w)
    real(pfdp),    intent(in   )                   :: y(:)
    complex(pfdp), intent(  out), dimension(:,:,:) :: u, v, w

    integer :: idx, i, j, k, nx, ny, nz

    nx = size(u, dim=1); ny = size(u, dim=2); nz = size(u, dim=3)

    if (size(y) /= 2*3*nx*ny*nz) stop "ERROR: Size mismatch in unpack3."

    idx = 1
    do k = 1, nz; do j = 1, ny; do i = 1, nx
       u(i,j,k) = complex(y(idx), y(idx+1)); idx = idx + 2
       v(i,j,k) = complex(y(idx), y(idx+1)); idx = idx + 2
       w(i,j,k) = complex(y(idx), y(idx+1)); idx = idx + 2
    end do; end do; end do
  end subroutine unpack3

end module feval
