!
! Copyright (C) 2014 Matthew Emmett and Michael Minion.
!

module pf_mod_dtype
  use pf_mod_config
  use sweeper_dtype
  use user_dtype
  implicit none

  integer, parameter :: PF_MAX_HOOKS = 32

  real(pfdp), parameter :: ZERO  = 0.0_pfdp
  real(pfdp), parameter :: ONE   = 1.0_pfdp
  real(pfdp), parameter :: TWO   = 2.0_pfdp
  real(pfdp), parameter :: HALF  = 0.5_pfdp

  integer, parameter :: SDC_GAUSS_LOBATTO   = 1
  integer, parameter :: SDC_GAUSS_RADAU     = 2
  integer, parameter :: SDC_CLENSHAW_CURTIS = 3
  integer, parameter :: SDC_UNIFORM         = 4
  integer, parameter :: SDC_GAUSS_LEGENDRE  = 5
  integer, parameter :: SDC_PROPER_NODES    = 2**8
  integer, parameter :: SDC_COMPOSITE_NODES = 2**9
  integer, parameter :: SDC_NO_LEFT         = 2**10

  type :: pf_hook
     procedure(pf_hook_p), pointer, nopass :: proc
  end type pf_hook

  type :: pf_level
     integer     :: ndofs = -1          ! number of degrees-of-freedom
     integer     :: nnodes = -1         ! number of sdc nodes
     integer     :: nsweeps = 1         ! number of sdc sweeps to perform
     integer     :: level = -1          ! level number (1 is the coarsest)
     logical     :: Finterp = .false.   ! interpolate functions instead of solutions

     real(pfdp)  :: residual

     real(pfdp), allocatable :: &
          Q(:,:), &                     ! unknowns at sdc nodes
          pQ(:,:), &                    ! unknowns at sdc nodes, previous sweep
          F(:,:,:), &                   ! function values at sdc nodes
          pF(:,:,:), &                  ! function values at sdc nodes, previous sweep
          q0(:), &                      ! initial condition
          qend(:), &                    ! end solution
          send(:), &                    ! send buffer
          recv(:), &                    ! recv buffer
          nodes(:), &                   ! sdc nodes
          qmat(:,:), &                  ! integration matrix (0 to node)
          smat(:,:), &                  ! integration matrix (node to node)
          tmat(:,:), &                  ! time interpolation matrix
          rmat(:,:), &                  ! time restriction matrix
          tau(:,:)                      ! fas corrections

     integer, allocatable :: &
          nflags(:)                     ! sdc node flags

     type(pf_sweeper) :: sweeper
     type(pf_user) :: user
  end type pf_level

  type :: pf_comm
     integer :: nproc = -1              ! total number of processors
     integer :: comm = -1               ! communicator
     integer, allocatable :: &
          recvreq(:), &                 ! receive requests (indexed by level)
          sendreq(:)                    ! send requests (indexed by level)
     integer :: statreq                 ! status send request
  end type pf_comm

  type :: pf_pfasst
     integer :: nlevels = -1            ! number of pfasst levels
     integer :: niters  = 5             ! number of iterations
     integer :: qtype   = SDC_GAUSS_LOBATTO
     integer :: rank    = -1            ! rank of current processor

     ! state
     integer    :: step, iter
     real(pfdp) :: t0, dt

     type(pf_hook), allocatable :: hooks(:,:,:)
     integer,       allocatable :: nhooks(:,:)

     ! real(pfdp) :: abs_res_tol = 0.d0
     ! real(pfdp) :: rel_res_tol = 0.d0

     logical :: Pipeline_G =  .false.
     logical :: PFASST_pred = .false.

     integer     :: taui0 = -999999     ! cutoff for tau inclusion

     type(pf_level), pointer :: levels(:)
     type(pf_comm),  pointer :: comm
  end type pf_pfasst

  interface
     subroutine pf_hook_p(pf, level)
       use iso_c_binding
       import pf_pfasst, pf_level
       type(pf_pfasst), intent(inout) :: pf
       type(pf_level),  intent(inout) :: level
     end subroutine pf_hook_p
  end interface

end module pf_mod_dtype
