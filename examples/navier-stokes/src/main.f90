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

    real(pfdp), allocatable :: q0(:)
    type(pf_level), pointer :: fine

    print *, '==> pfasst'

    call pf_mpi_create(comm, MPI_COMM_WORLD)
    call pf_pfasst_create(pf, comm, nlevs)

    pf%qtype  = SDC_GAUSS_RADAU
    pf%niters = niters
    pf%echo_timings = .true.

    if (pf%nlevels > 1) then
       pf%levels(1)%nsweeps = 2
    end if

    call pf_mpi_setup(comm, pf)
    call pf_pfasst_setup(pf, ndofs, nnodes)

    print *, 'npts:   ', pf%levels(:)%user%nx
    print *, 'ndofs:  ', pf%levels(:)%ndofs
    print *, 'nnodes: ', pf%levels(:)%nnodes
    print *, 'nsweeps:', pf%levels(:)%nsweeps

    pf%levels(:)%user%tol = 0.01

    fine => pf%levels(pf%nlevels)
    allocate(q0(fine%ndofs))
    if (exact) then
       call shapiro(fine%user%fft, q0, 0.d0, fine%user%nx, nu)
       call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error_hook)
    else
       call taylor_green(fine%user%fft, q0, 0.d0, fine%user%nx, nu)
    end if

    call pf_add_hook(pf, -1, PF_POST_STEP, echo_enstrophy_hook)
    call pf_pfasst_run(pf, q0, dt, 0.0d0, nsteps_in=nsteps)

    call pf_pfasst_destroy(pf)
    call pf_mpi_destroy(comm)

  end subroutine solve_pfasst

  subroutine solve_midpoint()
    use feval
    use hooks
    implicit none

    integer        :: n
    real(pfdp)     :: t, ens
    type(pf_level) :: lev
    character(128) :: fname

    real(pfdp), allocatable :: q(:), rhs(:)

    lev%level = 1
    lev%ndofs = ndofs(1)
    call user_setup(lev)

    print *, '==> implicit midpoint method'

    allocate(q(lev%ndofs), rhs(lev%ndofs))

    if (exact) then
       call shapiro(lev%user%fft, q, 0.d0, lev%user%nx, nu)
    else
       call taylor_green(lev%user%fft, q, 0.d0, lev%user%nx, nu)
    end if

    do n = 1, nsteps
       print *, 'step:', n
       t   = (n-1) * dt
       rhs = q
       call impl_solve(q, t, -dt/2, rhs, lev, n)
       q = 2*q - rhs

       if (exact) then
          call shapiro(lev%user%fft, rhs, t+dt, npts(1), nu)
          print *, '  error    ', maxval(abs(q-rhs))
       end if

       write (fname, "('velocity_s',i0.5)") n
       call dump_velocity(fname, q, lev)
       write (fname, "('vorticity_s',i0.5)") n
       call dump_vorticity(fname, q, lev)

       call enstrophy(q, lev, ens)
       print *, '  enstrophy', n, ens
    end do

    call user_destroy(lev)
  end subroutine solve_midpoint

end program main
