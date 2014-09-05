!
! Copyright (C) 2014 Matthew Emmett and Michael Minion.
!

module user_dtype
  use pf_mod_config
  implicit none

  type :: pf_user
     integer    :: nx, ny, nz
     real(pfdp) :: nu, scale, tol
     integer(8) :: fft, bft
     complex(pfdp), dimension(:,:,:), allocatable :: wk
     complex(pfdp), dimension(:), allocatable     :: kx, ky, kz
  end type pf_user
end module user_dtype
