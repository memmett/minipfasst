!
! Copyright (C) 2014 Matthew Emmett and Michael Minion.
!

!
! Sweeper declaration for IMEX sweeper.
!

module sweeper_dtype
  use pf_mod_config
  implicit none

  integer, parameter :: npieces = 1

  type :: pf_sweeper
     real(pfdp), allocatable :: &
          impl_mat(:,:)         ! implicit integration matrix
  end type pf_sweeper
end module sweeper_dtype
