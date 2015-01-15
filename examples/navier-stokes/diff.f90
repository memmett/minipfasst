program main
  use pfasst
  use feval
  use transfer
  use hooks
  use probin

  type(pf_level)          :: fine, crse
  real(pfdp), allocatable :: q_fine(:), q_itrp(:), q_crse(:)

  character(len=128) :: fname
  integer            :: i

  !
  ! build mock pfasst levels
  !
  npts(1:2) = [ 128, 256 ]

  crse%level = 1
  crse%ndofs = 2 * 3 * npts(1)**3
  call user_setup(crse)

  fine%level = 2
  fine%ndofs = 2 * 3 * npts(2)**3
  call user_setup(fine)

  allocate(q_fine(fine%ndofs),q_crse(crse%ndofs),q_itrp(fine%ndofs))

  !
  ! load and compute difference
  !
  call get_command_argument(1, fname)
  call read_velocity(trim(fname), q_fine, fine)

  call get_command_argument(2, fname)
  call read_velocity(trim(fname), q_crse, crse)
  call interpolate(q_itrp, q_crse, fine, crse, 0.d0)

  q_fine = q_fine - q_intrp

  call get_command_argument(3, fname)
  ! call dump_velocity(trim(fname), q_fine, fine)
  call dump_velocity_c('.' // c_null_char, trim(fname) // c_null_char, &
       fine%user%nx, fine%user%ny, fine%user%nz, q_fine)

end program main
