!
! Simple example of using MINIPFASST.
!

program main
  use pfasst
  use probin
  use initial
  implicit none

  character(len=256) :: probin_fname
  integer            :: err

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

  if (midpoint) then
     call solve_midpoint()
  else
     call solve_pfasst()
  end if

  call dfftw_cleanup_threads()
  call dfftw_cleanup()
  call mpi_finalize(err)

contains

  subroutine solve_pfasst()
    use pf_mod_mpi, only: MPI_COMM_WORLD
    use hooks
    implicit none

    type(pf_pfasst) :: pf
    type(pf_comm)   :: comm

    print *, '==> pfasst'

    call pf_mpi_create(comm, MPI_COMM_WORLD)
    call pf_pfasst_create(pf, comm, maxlevs)

    pf%qtype  = SDC_GAUSS_RADAU
    pf%niters = niters

    if (pf%nlevels > 1) then
       pf%levels(1)%nsweeps = 2
    end if

    pf%levels(:)%Finterp = .true.

    call pf_mpi_setup(comm, pf)
    call pf_pfasst_setup(pf, ndofs, nnodes)

    print *, 'npts:   ', pf%levels(:)%user%nx
    print *, 'ndofs:  ', pf%levels(:)%ndofs
    print *, 'nnodes: ', pf%levels(:)%nnodes
    print *, 'nsweeps:', pf%levels(:)%nsweeps

  ! allocate(q0(ndofs(maxlevs)))
  ! call initialcondition(q0, npts(maxlevs), nu)
  ! ! allocate(q0(ndofs(1)))
  ! ! call initial(q0, npts(1), nu)

  ! deallocate(q0)

    call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)
    ! call pf_pfasst_run(pf, q0, dt, 0.0d0, nsteps_in=nsteps)

    call pf_pfasst_destroy(pf)
    call pf_mpi_destroy(comm)

  end subroutine solve_pfasst

  subroutine solve_midpoint()
    use feval
    use hooks
    implicit none

    integer        :: n
    real(pfdp)     :: t
    type(pf_level) :: lev
    character(128) :: fname

    real(pfdp), allocatable :: q(:), rhs(:)

    lev%level = maxlevs
    lev%ndofs = ndofs(maxlevs)
    call user_setup(lev)

    print *, '==> implicit midpoint method'

    allocate(q(lev%ndofs), rhs(lev%ndofs))

    call shapiro(lev%user%fft, q, t, lev%user%nx, nu)

    do n = 1, nsteps
       print *, 'step:', n
       t   = (n-1) * dt
       rhs = q
       call impl_solve(q, t, -dt/2, rhs, lev)
       q = 2*q - rhs

       call shapiro(lev%user%fft, rhs, t, npts(maxlevs), nu)
       print *, '  error   ', maxval(abs(q-rhs))

       write (fname, "('velocity_s',i0.5)") n
       call dump_velocity(fname, q, lev)
       write (fname, "('vorticity_s',i0.5)") n
       call dump_vorticity(fname, q, lev)
    end do

    call user_destroy(lev)
  end subroutine solve_midpoint

end program main
