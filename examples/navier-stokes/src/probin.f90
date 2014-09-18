module probin
  use pf_mod_dtype

  integer, parameter :: maxlevs = 3

  integer, save           :: nlevs, nsteps, niters, nskip, nthreads
  integer, save           :: npts(maxlevs), ndofs(maxlevs), nnodes(maxlevs)
  double precision, save  :: nu, dt
  character(len=64), save :: input, output
  logical, save           :: midpoint, exact

contains

  subroutine probin_init(filename)
    character(len=*), intent(in) :: filename

    integer :: un

    namelist /prbin/ midpoint, exact
    namelist /prbin/ nu, dt, nsteps
    namelist /prbin/ input, output, nskip
    namelist /prbin/ nlevs, npts, nnodes, niters
    namelist /prbin/ nthreads

    !
    ! defaults
    !

    input  = ''
    output = ''

    nlevs    = 2
    nnodes   = [ 2, 3, 5 ]
    npts     = [ 16, 32, 64 ]
    niters   = 8
    nsteps   = -1
    nskip    = 1
    nthreads = 1
    nspace   = 1
    ntime    = -1

    nu       = 7.6d-4
    dt       = 0.0001d0

    midpoint = .false.
    exact    = .false.

    !
    ! read
    !

    un = 66
    open(unit=un, file=filename, status='old', action='read')
    read(unit=un, nml=prbin)
    close(unit=un)

    ndofs  = 2 * 3 * npts**3

    if (nsteps < 0) stop "nsteps must be set"

  end subroutine probin_init

end module probin
