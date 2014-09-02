!
! Simple example of using LIBPFASST.
!

program main
  use pfasst
  use pf_mod_mpi, only: MPI_COMM_WORLD
  use feval
  use hooks
  use transfer

  implicit none

  integer, parameter :: maxlevs = 3

  type(pf_pfasst) :: pf
  type(pf_comm)   :: comm
  integer         :: err, ndofs(maxlevs), nnodes(maxlevs), l
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

  ndofs  = [ 32, 64, 128 ]   ! number of dofs on the time/space levels
  nnodes = [ 2, 3, 5 ]       ! number of sdc nodes on time/space levels
  dt     = 0.01_pfdp

  call pf_mpi_create(comm, MPI_COMM_WORLD)
  call pf_pfasst_create(pf, comm, maxlevs)

  pf%qtype  = SDC_GAUSS_LOBATTO
  pf%niters = 4

  if (pf%nlevels > 1) then
     pf%levels(1)%nsweeps = 2
  end if

  call pf_mpi_setup(comm, pf)
  call pf_pfasst_setup(pf, ndofs, nnodes)

  call feval_create_workspaces(ndofs)

  !
  ! compute initial condition, add hooks, run
  !

  allocate(q0(ndofs(maxlevs)))
  call initial(q0)

  ! call pf_print_options(pf, show_mats=.true.)

  ! call pf_add_hook(pf, pf%nlevels, PF_POST_ITERATION, echo_error)
  call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)
  ! call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_residual)
  call pf_pfasst_run(pf, q0, dt, 0.0d0, nsteps=4)

  !
  ! cleanup
  !

  deallocate(q0)
  call feval_destroy_workspaces()

  call pf_pfasst_destroy(pf)
  call pf_mpi_destroy(comm)
  call mpi_finalize(err)
  call fftw_cleanup()

end program main
