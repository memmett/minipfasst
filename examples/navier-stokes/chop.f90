program chop
  use pfasst
  use feval
  use transfer
  use hooks
  use probin

  type(pf_level) :: fine, crse
  real(pfdp), allocatable :: q_fine(:), q_crse(:)

  npts = [ 128, 256, 256 ]
  output = './'

  ! build a mock pfasst levels
  fine%level = 2
  fine%ndofs = 2 * 3 * 256**3
  call user_setup(fine)

  crse%level = 1
  crse%ndofs = 2 * 3 * 128**3
  call user_setup(crse)

  allocate(q_fine(fine%ndofs),q_crse(crse%ndofs))

  ! read and restrict fine.dat
  call read_velocity("fine.dat", q_fine, fine)
  call restrict(q_fine, q_crse, fine, crse, 0.d0)

  ! interpolate and write fine_crse.dat
  q_fine = 0
  call interpolate(q_fine, q_crse, fine, crse, 0.d0)
  call dump_velocity("fine_crse", q_fine, fine)

end program chop
