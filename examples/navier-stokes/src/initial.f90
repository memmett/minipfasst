!
! Copyright (c) 2014, Matthew Emmett and Michael Minion.  All rights reserved.
!

module initial
  use feval
  implicit none
contains

  subroutine initialcondition(q0, nx, nu)
    real(pfdp), intent(  out) :: q0(:)
    real(pfdp), intent(in   ) :: nu
    integer,    intent(in   ) :: nx
    ! call shapiro(q0, 0.d0, nx, nu)
    call taylor_green(q0, 0.d0, nx, nu)
  end subroutine initialcondition

  subroutine shapiro(yex, t, nx, nu)
    use probin, only: npts
    real(pfdp), intent(in   ) :: t, nu
    real(pfdp), intent(  out) :: yex(:)
    integer,    intent(in   ) :: nx

    integer    :: i, j, k
    real(pfdp) :: sqrt3, expt, cosx, sinx, cosy, siny, cosz, sinz

    complex(pfdp), dimension(nx,nx,nx) :: u, v, w
    real(pfdp), dimension(nx)          :: x

    do i = -nx/2, nx/2-1
       x(nx/2+i+1)= 2*pi*dble(i)/nx
    end do

    sqrt3 = sqrt(3.d0)
    expt  = exp(-3.d0 * t * nu)

    do k = 1, nx
       cosz = cos(x(k))
       sinz = sin(x(k))
       do j = 1, nx
          cosy = cos(x(j))
          siny = sin(x(j))
          do i = 1, nx
             cosx = cos(x(i))
             sinx = sin(x(i))

             u(i,j,k) = -0.5 * ( sqrt3 * cosx * siny * sinz + sinx * cosy * cosz ) * expt
             v(i,j,k) =  0.5 * ( sqrt3 * sinx * cosy * sinz - cosx * siny * cosz ) * expt
             w(i,j,k) = cosx * cosy * sinz * expt
          end do
       end do
    end do

    call pack3(yex, u, v, w)
  end subroutine shapiro

  subroutine taylor_green(yex, t, nx, nu)
    use probin, only: npts
    real(pfdp), intent(in   ) :: t, nu
    real(pfdp), intent(  out) :: yex(:)
    integer,    intent(in   ) :: nx

    integer    :: i, j, k
    real(pfdp) :: ucoeff, vcoeff, wcoeff, theta, cosx, sinx, cosy, siny, cosz, sinz

    complex(pfdp), dimension(nx,nx,nx) :: u, v, w
    real(pfdp), dimension(nx)          :: x

    do i = -nx/2, nx/2-1
       x(nx/2+i+1)= 2*pi*dble(i)/nx
    end do

    theta  = 0.d0
    ucoeff = 2.d0/sqrt(3.d0) * sin(theta + 2*pi/3)
    vcoeff = 2.d0/sqrt(3.d0) * sin(theta - 2*pi/3)
    wcoeff = 2.d0/sqrt(3.d0) * sin(theta)

    do k = 1, nx
       cosz = cos(x(k))
       sinz = sin(x(k))
       do j = 1, nx
          cosy = cos(x(j))
          siny = sin(x(j))
          do i = 1, nx
             cosx = cos(x(i))
             sinx = sin(x(i))

             u(i,j,k) = ucoeff * sinx * cosy * cosz
             v(i,j,k) = vcoeff * cosx * siny * cosz
             w(i,j,k) = wcoeff * cosx * cosy * sinz
          end do
       end do
    end do

    call pack3(yex, u, v, w)
  end subroutine taylor_green

  ! Evaluate F(U).
  subroutine vorticity(y, lev, vort)
    real(pfdp),     intent(in   ) :: y(:)
    type(pf_level), intent(inout) :: lev
    real(pfdp),     intent(  out) :: vort(:,:,:)

    integer :: i, j, k, nx, ny, nz

    complex(pfdp), dimension(lev%user%nx,lev%user%ny,lev%user%nz) :: &
         u, v, w, uhat, vhat, what, utmp, vtmp, wtmp, ux, uy, uz, vx, vy, vz, wx, wy, wz

    complex(pfdp) :: kx(lev%user%nx), ky(lev%user%ny), kz(lev%user%nz)
    real(pfdp)    :: nu, scale

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

    vort = realpart((wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2)
  end subroutine vorticity

end module initial
