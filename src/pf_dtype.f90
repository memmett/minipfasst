!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
! 
!   1. Redistributions of source code must retain the above copyright
!      notice, this list of conditions and the following disclaimer.
! 
!   2. Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions and the following disclaimer in
!      the documentation and/or other materials provided with the
!      distribution.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! 

module pf_mod_dtype
  use iso_c_binding
  implicit none


  ! state type
  type :: pf_state_t
     double precision :: t0, dt
     integer          :: cycle, step, iter, nsteps
  end type pf_state_t


  ! level type
  type :: pf_level_t
     integer     :: nvars = -1          ! number of variables (dofs)
     integer     :: nnodes = -1         ! number of sdc nodes
     integer     :: nsweeps = 1         ! number of sdc sweeps to perform
     integer     :: level = -1          ! level number (1=finest)

     ! arrays
     double precision, pointer :: &
          qSDC(:,:), &                  ! unknowns at sdc nodes
          fSDC(:,:,:), &                ! imex functions values at sdc nodes
          tau(:,:), &                   ! fas correction
          q0(:), &                      ! initial condition
          qend(:), &                    ! end value
          send(:), &                    ! send buffer
          recv(:), &                    ! recv buffer
          nodes(:), &                   ! sdc nodes
          qmat(:,:), &                  ! integration matrix
          smat(:,:,:) => null(), &      ! sdc matrices (allocated by the sweeper)
          tmat(:,:)                     ! time interpolation matrix

     logical :: allocated = .false.
  end type pf_level_t


  ! pfasst communicator
  type :: pf_comm_t
     integer :: nproc = -1              ! total number of processors

     ! mpi
     integer :: comm = -1               ! communicator
     integer, pointer :: &
          recvreq(:), &                 ! receive requests (indexed by level)
          sendreq(:)                    ! send requests (indexed by level)
  end type pf_comm_t


  ! pfasst type
  type :: pf_pfasst_t
     integer :: nlevels = -1            ! number of pfasst levels
     integer :: niters  = 5             ! number of iterations
     integer :: qtype   = 1             ! type of quadrature nodes
     integer :: rank    = -1            ! rank of current processor

     ! pf objects
     type(pf_state_t)          :: state
     type(pf_level_t), pointer :: levels(:)
     type(pf_comm_t), pointer  :: comm
  end type pf_pfasst_t

end module pf_mod_dtype
