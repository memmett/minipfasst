!
! Copyright (c) 2014, Matthew Emmett and Michael Minion.  All rights reserved.
!

module initial
  use feval
  implicit none

  real(pfdp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510_pfdp

contains

  subroutine shapiro(plan, q, t, nx, nu)
    use probin, only: npts
    real(pfdp), intent(in   ) :: t, nu
    real(pfdp), intent(  out) :: q(:)
    integer,    intent(in   ) :: nx
    integer(8), intent(in   ) :: plan

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

    u = u / nx**3; v = v / nx**3; w = w / nx**3

    call dfftw_execute_dft_(plan, u, u)
    call dfftw_execute_dft_(plan, v, v)
    call dfftw_execute_dft_(plan, w, w)

    call pack3(q, u, v, w)
  end subroutine shapiro

  subroutine taylor_green(q, t, nx, nu)
    use probin, only: npts
    real(pfdp), intent(in   ) :: t, nu
    real(pfdp), intent(  out) :: q(:)
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

    call pack3(q, u, v, w)
  end subroutine taylor_green

  ! Evaluate F(U).
  subroutine vorticity(y, lev, vort)
    real(pfdp),     intent(in   ) :: y(:)
    type(pf_level), intent(inout) :: lev
    real(pfdp),     intent(  out) :: vort(:,:,:)

    integer :: i, j, k, nx, ny, nz

    complex(pfdp), dimension(lev%user%nx,lev%user%ny,lev%user%nz) :: &
         uhat, vhat, what, uy, uz, vx, vz, wx, wy

    complex(pfdp) :: kx(lev%user%nx), ky(lev%user%ny), kz(lev%user%nz)

    ! shortcuts
    nx = lev%user%nx; ny = lev%user%ny; nz = lev%user%nz
    kx = lev%user%kx; ky = lev%user%ky; kz = lev%user%kz

    ! unpack
    call unpack3(y, uhat, vhat, what)

    do k=1,nz ; do j=1,ny ; do i=1,nx
       uy(i,j,k) = uhat(i,j,k) * ky(j)
       uz(i,j,k) = uhat(i,j,k) * kz(k)
    end do; end do ; end do
    call dfftw_execute_dft_(lev%user%bft, uy, uy)
    call dfftw_execute_dft_(lev%user%bft, uz, uz)

    do k=1,nz ; do j=1,ny ; do i=1,nx
       vx(i,j,k) = vhat(i,j,k) * kx(i)
       vz(i,j,k) = vhat(i,j,k) * kz(k)
    end do; end do ; end do
    call dfftw_execute_dft_(lev%user%bft, vx, vx)
    call dfftw_execute_dft_(lev%user%bft, vz, vz)

    do k=1,nz ; do j=1,ny ; do i=1,nx
       wx(i,j,k) = what(i,j,k) * kx(i)
       wy(i,j,k) = what(i,j,k) * ky(j)
    end do; end do ; end do
    call dfftw_execute_dft_(lev%user%bft, wx, wx)
    call dfftw_execute_dft_(lev%user%bft, wy, wy)

    vort = realpart((wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2)
  end subroutine vorticity

end module initial
