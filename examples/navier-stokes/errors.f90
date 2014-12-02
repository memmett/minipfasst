program main
  use pfasst
  use feval
  use transfer
  use hooks
  use probin

  type(pf_level)          :: fine, crse
  real(pfdp), allocatable :: q_fine(:), q_rstr(:), q_crse(:)

  character(len=128) :: fname
  integer            :: i

  !
  ! build a mock pfasst levels
  !
  npts(1:2) = [ 128, 256 ]

  crse%level = 1
  crse%ndofs = 2 * 3 * npts(1)**3
  call user_setup(crse)

  fine%level = 2
  fine%ndofs = 2 * 3 * npts(2)**3
  call user_setup(fine)

  allocate(q_fine(fine%ndofs),q_crse(crse%ndofs),q_rstr(crse%ndofs))

  !
  ! load fine solution and restrict
  !
  call get_command_argument(1, fname)
  call read_velocity(trim(fname), q_fine, fine)
  call restrict(q_fine, q_rstr, fine, crse, 0.d0)

  !
  ! load coarse solutions and echo errors
  !
  do i = 2, command_argument_count()
     call get_command_argument(i, fname)
     call read_velocity(trim(fname), q_crse, crse)
     print *, fname, maxval(abs(q_rstr-q_crse))
  end do
end program main
