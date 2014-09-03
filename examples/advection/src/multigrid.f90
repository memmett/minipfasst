module multigrid
  use pf_mod_dtype
  implicit none

  real(pfdp), parameter :: lx = 1
  integer, parameter :: spatial_order = 2

contains

  subroutine relax(y, nx, t, nudt, rhs, nrelax)
    integer,        intent(in   ) :: nx
    real(pfdp),     intent(inout) :: y(nx)
    real(pfdp),     intent(in   ) :: rhs(nx), t, nudt
    integer,        intent(in   ) :: nrelax

    integer    :: i, n
    real(pfdp) :: dx, sig, sigi
    real(pfdp) :: ybc(1-spatial_order:nx+spatial_order)

    dx = lx/dble(nx)

    sig = nudt/(dx*dx)
    sigi = (dx*dx)/(nudt)

    ! Jacobi
    do n = 1,nrelax
       ybc(1:Nx) = y
       call set_bc(ybc, nx, spatial_order)
       y = ( rhs + sig * (ybc(2:Nx+1) + ybc(0:Nx-1)) ) / ( one + two*sig )
    end do

    ! Gauss-Seidel
    ! ybc(1:nx) = y
    ! do n = 1, nrelax
    !    call set_bc(ybc, nx, spatial_order)
    !    do i = 1,nx
    !       ybc(i) = ( rhs(i) + sig * ( ybc(i-1) + ybc(i+1) ) ) / ( one + two*sig )
    !    end do
    ! end do
    ! y = ybc(1:nx)
  end subroutine relax

  recursive subroutine multigrid_v_cycle(y, t, nudt, rhs, nrelax, maxresidual)
    real(pfdp), intent(inout) :: y(:)
    real(pfdp), intent(in   ) :: rhs(:), t, nudt
    integer,    intent(in   ) :: nrelax
    real(pfdp), intent(inout) :: maxresidual

    integer :: i, nx
    real(pfdp), allocatable :: res(:), resc(:), corrf(:), corr(:)

    nx = size(y)

    if (nx .lt. 8)  then
       do i = 1,4
          call relax(y, nx, t, nudt, rhs, nrelax)
       end do
    else
       allocate(res(nx),resc(nx/2),corrf(nx),corr(nx/2))

       call relax(y, nx, t, nudt, rhs, nrelax)
       call residual(y, nx, t, nudt, rhs, res)
       call coarsen(res, nx, resc)

       corr = 0
       call multigrid_v_cycle(corr, t, nudt, resc, nrelax, maxresidual)

       call interp(corrf, nx, corr)
       y = y + corrf
       call relax(y, nx, t, nudt, rhs, nrelax)
       call residual(y, nx, t, nudt, rhs, res)

       maxresidual = maxval(abs(res))
       deallocate(res,resc,corrf,corr)
    end if
  end subroutine multigrid_v_cycle

  subroutine coarsen(y, nx, yc)
    integer,    intent(in   ) :: nx
    real(pfdp), intent(in   ) :: y(nx)
    real(pfdp), intent(  out) :: yc(nx/2)

    real(pfdp) :: ybc(1-spatial_order:nx+spatial_order)

    call fill_bc(y, ybc, nx, spatial_order)
    yc = half * ybc(1:nx-1:2) + quarter * ( ybc(0:nx-2:2) + ybc(2:nx:2) )
  end subroutine coarsen

  subroutine interp(y, nx, yc)
    integer,    intent(in   ) :: nx
    real(pfdp), intent(  out) :: y(nx)
    real(pfdp), intent(in   ) :: yc(nx/2)

    real(pfdp) :: ybc(1-spatial_order:nx/2+spatial_order)

    call fill_bc(yc, ybc, nx/2, spatial_order)
    y(1:nx-1:2) = ybc(1:nx/2)
    y(2:nx:2)   = (-ybc(0:nx/2-1) + 9.0_pfdp * (ybc(1:nx/2) + ybc(2:nx/2+1)) - ybc(3:nx/2+2))/16.0_pfdp
  end subroutine interp

  subroutine residual(y, nx, t, nudt, rhs, res)
    integer,        intent(in   ) :: nx
    real(pfdp),     intent(in   ) :: y(nx), rhs(nx), t, nudt
    real(pfdp),     intent(  out) :: res(nx)

    real(pfdp) :: dx, ybc(1-spatial_order:nx+spatial_order), lap(nx)

    dx = lx/dble(nx)

    call fill_bc(y, ybc, nx, spatial_order)
    lap = ( ybc(0:nx-1) - 2*ybc(1:nx) + ybc(2:nx+1) ) / (dx*dx)
    res = rhs - (y - nudt*lap)
  end subroutine residual

  subroutine fill_bc(y, ybc, nx, nbc)
    integer,    intent(in   ) :: nx, nbc
    real(pfdp), intent(in   ) :: y(1:nx)
    real(pfdp), intent(  out) :: ybc(1-nbc:nx+nbc)

    ybc(1:Nx)=y
    call set_bc(ybc,Nx,Nbc)
  end subroutine fill_bc

  subroutine set_bc(ybc, nx, nbc)
    integer,    intent(in   ) :: nx, nbc
    real(pfdp), intent(inout) :: ybc(1-Nbc:Nx+Nbc)

    ! print *, 'pre'
    ! print *, '0', ybc(1:nbc+1)
    ! print *, '1', ybc(nx-nbc:nx)

    ybc(1-nbc:0) = ybc(nx-nbc+1:nx)
    ybc(nx+1:nx+nbc) = ybc(1:nbc)

    ! print *, 'post'
    ! print *, '0', ybc(:nbc+1)
    ! print *, '1', ybc(nx-nbc:)

    ! do i = 1, Nbc
    !    ybc(1-i)=-ybc(1+i)
    ! end do
    ! ybc(Nx+1) = 0
    ! do i = 1, Nbc-1
    !    ybc(Nx+1+i)=-ybc(Nx+1-i)
    ! end do
  end subroutine set_bc

end module multigrid
