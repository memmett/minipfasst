!
! Copyright (c) 2014, Matthew Emmett and Michael Minion.  All rights reserved.
!

! See the README.org file for details.

module feval
  use pf_mod_dtype
  implicit none

  real(pfdp), parameter :: pi=3.14159265358979323846264338327950288419716939937510_pfdp

contains

  ! Evaluate F(U).
  subroutine impl_eval(y, t, lev, f)
    real(pfdp),     intent(in   ) :: y(:), t
    type(pf_level), intent(inout) :: lev
    real(pfdp),     intent(  out) :: f(:)

    integer :: i, j, k, nx, ny, nz

    complex(pfdp), dimension(lev%user%nx,lev%user%ny,lev%user%nz) :: &
         u, v, w, uhat, vhat, what, utmp, vtmp, wtmp, &
         ux, uy, uz, vx, vy, vz, wx, wy, wz, nluhat, nlvhat, nlwhat, phat

    complex(pfdp) :: kx(lev%user%nx), ky(lev%user%ny), kz(lev%user%nz), diffop

    real(pfdp) :: nu, scale

    ! shortcuts
    nx = lev%user%nx; ny = lev%user%ny; nz = lev%user%nz
    kx = lev%user%kx; ky = lev%user%ky; kz = lev%user%kz
    nu = lev%user%nu; scale = lev%user%scale

    ! unpack and transform
    call unpack3(y, u, v, w)
    call fft3(lev%user%fft, lev%user%wk, scale, u, v, w, uhat, vhat, what)

    ! gradients
    do k=1,nz ; do j=1,ny ; do i=1,nx
       utmp(i,j,k) = uhat(i,j,k) * kx(i)
       vtmp(i,j,k) = uhat(i,j,k) * ky(j)
       wtmp(i,j,k) = uhat(i,j,k) * kz(k)
    end do; end do ; end do
    call fft3(lev%user%bft, lev%user%wk, 1.d0, utmp, vtmp, wtmp, ux, uy, uz)

    do k=1,nz ; do j=1,ny ; do i=1,nx
       utmp(i,j,k) = vhat(i,j,k) * kx(i)
       vtmp(i,j,k) = vhat(i,j,k) * ky(j)
       wtmp(i,j,k) = vhat(i,j,k) * kz(k)
    end do; end do ; end do
    call fft3(lev%user%bft, lev%user%wk, 1.d0, utmp, vtmp, wtmp, vx, vy, vz)

    do k=1,nz ; do j=1,ny ; do i=1,nx
       utmp(i,j,k) = what(i,j,k) * kx(i)
       vtmp(i,j,k) = what(i,j,k) * ky(j)
       wtmp(i,j,k) = what(i,j,k) * kz(k)
    end do; end do ; end do
    call fft3(lev%user%bft, lev%user%wk, 1.d0, utmp, vtmp, wtmp, wx, wy, wz)

    ! nonlinear terms
    utmp = u * ux + v * uy + w * uz
    vtmp = u * vx + v * vy + w * vz
    wtmp = u * wx + v * wy + w * wz
    call fft3(lev%user%fft, lev%user%wk, scale, utmp, vtmp, wtmp, nluhat, nlvhat, nlwhat)

    ! pressure
    do k=1,nz
       do j=1,ny
          do i=1,nx
             phat(i,j,k) = -1.0d0 * ( &
                  + kx(i)*nluhat(i,j,k) &
                  + ky(j)*nlvhat(i,j,k) &
                  + kz(k)*nlwhat(i,j,k) ) &
                  / (kx(i)*kx(i) + ky(j)*ky(j) + kz(k)*kz(k) + 0.1d0**13)
          end do
       end do
    end do

    ! evaluate
    do k=1,nz
       do j=1,ny
          do i=1,nx
             diffop = nu * ( kx(i)*kx(i) + ky(j)*ky(j) + kz(k)*kz(k) )
             uhat(i,j,k) = - nluhat(i,j,k) - kx(i)*phat(i,j,k) + diffop * uhat(i,j,k)
             vhat(i,j,k) = - nlvhat(i,j,k) - ky(j)*phat(i,j,k) + diffop * vhat(i,j,k)
             what(i,j,k) = - nlwhat(i,j,k) - kz(k)*phat(i,j,k) + diffop * what(i,j,k)
          end do
       end do
    end do

    call fft3(lev%user%bft, lev%user%wk, 1.d0, uhat, vhat, what, u, v, w)

    call pack3(f, u, v, w)
  end subroutine impl_eval

  ! Solve U + a F(U) = RHS.
  subroutine impl_solve(y, t, a, rhs, lev)
    real(pfdp),     intent(inout) :: y(:)
    real(pfdp),     intent(in   ) :: rhs(:)
    real(pfdp),     intent(in   ) :: t, a
    type(pf_level), intent(inout) :: lev

    integer :: i, j, k, nx, ny, nz

    complex(pfdp), dimension(lev%user%nx,lev%user%ny,lev%user%nz) :: &
         u, v, w, uhat, vhat, what, u0hat, v0hat, w0hat, utmp, vtmp, wtmp, &
         ux, uy, uz, vx, vy, vz, wx, wy, wz, nluhat, nlvhat, nlwhat, phat

    complex(pfdp) :: kx(lev%user%nx), ky(lev%user%ny), kz(lev%user%nz), invop

    real(pfdp) :: nu, residual, scale

    ! shortcuts
    nx = lev%user%nx; ny = lev%user%ny; nz = lev%user%nz
    kx = lev%user%kx; ky = lev%user%ky; kz = lev%user%kz
    nu = lev%user%nu; scale = lev%user%scale

    ! unpack and transform rhs
    call unpack3(rhs, utmp, vtmp, wtmp)
    call fft3(lev%user%fft, lev%user%wk, scale, utmp, vtmp, wtmp, u0hat, v0hat, w0hat)

    ! unpack and transform initial guess
    call unpack3(y, u, v, w)
    call fft3(lev%user%fft, lev%user%wk, scale, u, v, w, uhat, vhat, what)

    residual = 1
    do while (residual > lev%user%tol)

       ! gradients
       do k=1,nz ; do j=1,ny ; do i=1,nx
          utmp(i,j,k) = uhat(i,j,k) * kx(i)
          vtmp(i,j,k) = uhat(i,j,k) * ky(j)
          wtmp(i,j,k) = uhat(i,j,k) * kz(k)
       end do; end do ; end do
       call fft3(lev%user%bft, lev%user%wk, 1.d0, utmp, vtmp, wtmp, ux, uy, uz)

       do k=1,nz ; do j=1,ny ; do i=1,nx
          utmp(i,j,k) = vhat(i,j,k) * kx(i)
          vtmp(i,j,k) = vhat(i,j,k) * ky(j)
          wtmp(i,j,k) = vhat(i,j,k) * kz(k)
       end do; end do ; end do
       call fft3(lev%user%bft, lev%user%wk, 1.d0, utmp, vtmp, wtmp, vx, vy, vz)

       do k=1,nz ; do j=1,ny ; do i=1,nx
          utmp(i,j,k) = what(i,j,k) * kx(i)
          vtmp(i,j,k) = what(i,j,k) * ky(j)
          wtmp(i,j,k) = what(i,j,k) * kz(k)
       end do; end do ; end do
       call fft3(lev%user%bft, lev%user%wk, 1.d0, utmp, vtmp, wtmp, wx, wy, wz)

       ! nonlinear terms
       utmp = u * ux + v * uy + w * uz
       vtmp = u * vx + v * vy + w * vz
       wtmp = u * wx + v * wy + w * wz
       call fft3(lev%user%fft, lev%user%wk, scale, utmp, vtmp, wtmp, nluhat, nlvhat, nlwhat)

       ! scaled pressure
       do k=1,nz
          do j=1,ny
             do i=1,nx
                phat(i,j,k) = -1.0d0 * ( &
                     + kx(i)*nluhat(i,j,k) &
                     + ky(j)*nlvhat(i,j,k) &
                     + kz(k)*nlwhat(i,j,k) ) &
                     / (kx(i)*kx(i) + ky(j)*ky(j) + kz(k)*kz(k) + 0.1d0**13)
             end do
          end do
       end do

       ! solve
       do k=1,nz
          do j=1,ny
             do i=1,nx
                invop = 1.d0 / (1 + a * nu * ( kx(i)*kx(i) + ky(j)*ky(j) + kz(k)*kz(k) ) )
                uhat(i,j,k) = invop * ( u0hat(i,j,k) + a * nluhat(i,j,k) + a * kx(i)*phat(i,j,k) )
                vhat(i,j,k) = invop * ( v0hat(i,j,k) + a * nlvhat(i,j,k) + a * ky(j)*phat(i,j,k) )
                what(i,j,k) = invop * ( w0hat(i,j,k) + a * nlwhat(i,j,k) + a * kz(k)*phat(i,j,k) )
             end do
          end do
       end do

       utmp = u; vtmp = v; wtmp = w
       call fft3(lev%user%bft, lev%user%wk, 1.d0, uhat, vhat, what, u, v, w)

       residual = maxval(abs(utmp - u)) + maxval(abs(vtmp - v)) + maxval(abs(wtmp - w))
       print *, '  residual', residual
    end do

    call pack3(y, u, v, w)
  end subroutine impl_solve

  subroutine fft3(plan, wk, scale, u, v, w, uhat, vhat, what)
    integer(8),    intent(in   )                   :: plan
    real(pfdp),    intent(in   )                   :: scale
    complex(pfdp), intent(in   ), dimension(:,:,:) :: u, v, w
    complex(pfdp), intent(  out), dimension(:,:,:) :: wk, uhat, vhat, what
    wk = u * scale; call dfftw_execute_dft_(plan, wk, uhat)
    wk = v * scale; call dfftw_execute_dft_(plan, wk, vhat)
    wk = w * scale; call dfftw_execute_dft_(plan, wk, what)
  end subroutine fft3

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
