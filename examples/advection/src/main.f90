!
! Simple example of using MINIPFASST.
!

program main
  use pfasst
  use pf_mod_mpi, only: MPI_COMM_WORLD
  use feval, only: initial
  use hooks

  implicit none

  integer, parameter :: maxlevs = 3

  type(pf_pfasst) :: pf
  type(pf_comm)   :: comm
  integer         :: err, ndofs(maxlevs), nnodes(maxlevs)
  real(pfdp)      :: dt

  real(pfdp), allocatable :: q0(:)

  !
  ! initialize mpi
  !

  call mpi_init(err)
  if (err .ne. 0) &
       stop "ERROR: Can't initialize MPI."

  !
  ! initialize pfasst
  !

  dt = 0.01_pfdp
  ndofs  = [ 128, 256, 512 ]   ! number of dofs on the time/space levels
  nnodes = [ 2, 3, 5 ]       ! number of sdc nodes on time/space levels

  pf%qtype  = SDC_UNIFORM + SDC_NO_LEFT
  pf%niters = 4

  call pf_mpi_create(comm, MPI_COMM_WORLD)
  call pf_pfasst_create(pf, comm, maxlevs)

  if (pf%nlevels > 1) then
     pf%levels(1)%nsweeps = 2
  end if

  call pf_mpi_setup(comm, pf)
  call pf_pfasst_setup(pf, ndofs, nnodes)

  !
  ! compute initial condition, add hooks, run
  !

  allocate(q0(ndofs(maxlevs)))
  call initial(q0)

  call pf_print_options(pf, show_mats=.true.)
  call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)
  call pf_pfasst_run(pf, q0, dt, 0.0d0, nsteps_in=4)

  !
  ! cleanup
  !

  deallocate(q0)
  call pf_pfasst_destroy(pf)
  call pf_mpi_destroy(comm)
  call mpi_finalize(err)

end program main
