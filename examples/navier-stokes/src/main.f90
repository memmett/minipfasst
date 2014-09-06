!
! Simple example of using MINIPFASST.
!

program main
  use pfasst
  use probin
  use initial
  implicit none

  real(pfdp), allocatable :: q0(:)
  character(len=256)      :: probin_fname

  integer :: err

  if (command_argument_count() >= 1) then
     call get_command_argument(1, value=probin_fname)
  else
     probin_fname = "probin.nml"
  end if
  call probin_init(probin_fname)

  call dfftw_init_threads(err)
  call dfftw_plan_with_nthreads(nthreads)

  call mpi_init(err)
  if (err .ne. 0) &
       stop "ERROR: Can't initialize MPI."

  allocate(q0(ndofs(maxlevs)))
  call initialcondition(q0, npts(maxlevs), nu)
  ! allocate(q0(ndofs(1)))
  ! call initial(q0, npts(1), nu)
  if (midpoint) then
     call solve_midpoint(q0)
  else
     call solve_pfasst(q0)
  end if
  deallocate(q0)

  call mpi_finalize(err)

contains

  subroutine solve_pfasst(q0)
    use pf_mod_mpi, only: MPI_COMM_WORLD
    use hooks
    implicit none
    real(pfdp), intent(in   ) :: q0(:)

    type(pf_pfasst) :: pf
    type(pf_comm)   :: comm

    print *, '==> pfasst'

    call pf_mpi_create(comm, MPI_COMM_WORLD)
    call pf_pfasst_create(pf, comm, maxlevs)

    pf%qtype  = SDC_GAUSS_LOBATTO
    pf%niters = niters

    if (pf%nlevels > 1) then
       pf%levels(1)%nsweeps = 2
    end if

    call pf_mpi_setup(comm, pf)
    call pf_pfasst_setup(pf, ndofs, nnodes)

    print *, 'npts:   ', pf%levels(:)%user%nx
    print *, 'ndofs:  ', pf%levels(:)%ndofs
    print *, 'nnodes: ', pf%levels(:)%nnodes
    print *, 'nsweeps:', pf%levels(:)%nsweeps

    call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)
    call pf_pfasst_run(pf, q0, dt, 0.0d0, nsteps_in=nsteps)

    call pf_pfasst_destroy(pf)
    call pf_mpi_destroy(comm)

  end subroutine solve_pfasst

  subroutine solve_midpoint(q0)
    use feval
    use hooks
    implicit none
    real(pfdp), intent(in   ) :: q0(:)

    integer        :: n
    real(pfdp)     :: t, q1(size(q0)), rhs(size(q0))
    type(pf_level) :: lev
    character(128) :: fname

    lev%level = maxlevs
    lev%ndofs = ndofs(maxlevs)
    call user_setup(lev)

    print *, '==> implicit midpoint method'

    q1 = q0
    t  = 0
    do n = 1, nsteps
       print *, 'step:', n
       rhs = q1
       call impl_solve(q1, t, -dt/2, rhs, lev)
       q1 = 2*q1 - rhs
       t = t + dt

       ! call shapiro(rhs, t, npts(maxlevs), nu)
       ! print *, '  error   ', maxval(abs(q1-rhs))

       write (fname, "('velocity_s',i0.5)") n
       call dump_velocity(fname, q1, lev)

       write (fname, "('vorticity_s',i0.5)") n
       call dump_vorticity(fname, q1, lev)
    end do
  end subroutine solve_midpoint

end program main
