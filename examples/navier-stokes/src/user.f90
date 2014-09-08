!
! Copyright (C) 2014 Matthew Emmett and Michael Minion.
!

module fftw
  use iso_c_binding
  include 'fftw3.f03'
end module fftw

module user
  use pf_mod_dtype
  implicit none
contains

  subroutine user_setup(lev)
    use probin, only: npts, nu
    use fftw, only: FFTW_FORWARD, FFTW_BACKWARD, FFTW_ESTIMATE
    type(pf_level), intent(inout) :: lev

    integer :: i, nx, ny, nz

    complex(pfdp), allocatable :: wk(:,:,:)

    lev%user%nx = npts(lev%level)
    lev%user%ny = npts(lev%level)
    lev%user%nz = npts(lev%level)

    nx = lev%user%nx
    ny = lev%user%ny
    nz = lev%user%nz

    allocate(wk(nx,ny,nz))
    allocate(lev%user%kx(nx), lev%user%ky(ny), lev%user%kz(nz))

    call dfftw_plan_dft_3d_(lev%user%fft, nx, ny, nz, wk, wk, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_3d_(lev%user%bft, nx, ny, nz, wk, wk, FFTW_BACKWARD, FFTW_ESTIMATE)

    do i = 1, nx/2+1
       lev%user%kx(i) = (0.0d0,1.0d0) * dble(i-1)
    end do
    lev%user%kx(1+nx/2) = 0.0d0
    do i = 1, nx/2 -1
       lev%user%kx(i+1+nx/2) = -lev%user%kx(1-i+nx/2)
    end do

    lev%user%ky = lev%user%kx
    lev%user%kz = lev%user%kx

    lev%user%scale = 1.0d0/dble(nx*ny*nz)
    lev%user%tol   = 1.d-8
    lev%user%nu    = nu
  end subroutine user_setup

  subroutine user_destroy(lev)
    type(pf_level), intent(inout) :: lev
    call dfftw_destroy_plan_(lev%user%fft)
    call dfftw_destroy_plan_(lev%user%bft)
  end subroutine user_destroy
end module user
