!
! Copyright (C) 2014 Matthew Emmett and Michael Minion.
!

!
! User workspace for advection/diffusion example.
!

module user_dtype
  use pf_mod_config
  use iso_c_binding, only: c_ptr
  implicit none

  type :: pf_user
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:)              ! work space
     complex(pfdp), pointer :: ddx(:), lap(:)     ! operators
  end type pf_user
end module user_dtype
