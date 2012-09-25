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

! This module implements MPI communications.

module pf_mod_mpi
  include "mpif.h"
end module pf_mod_mpi

module pf_mod_comm
  use pf_mod_dtype

  implicit none

  interface create
     module procedure pf_comm_create
  end interface create

  interface setup
     module procedure pf_comm_setup
  end interface setup

  interface destroy
     module procedure pf_comm_destroy
  end interface destroy

contains

  ! Create an MPI based PFASST communicator
  subroutine pf_comm_create(pf_comm, mpi_comm)
    type(pf_comm_t), intent(out) :: pf_comm
    integer,         intent(in)  :: mpi_comm

    integer :: ierror

    pf_comm%comm = mpi_comm
    call mpi_comm_size(mpi_comm, pf_comm%nproc, ierror)
  end subroutine pf_comm_create

  ! Setup
  subroutine pf_comm_setup(pf_comm, pf)
    type(pf_comm_t),   intent(inout) :: pf_comm
    type(pf_pfasst_t), intent(inout) :: pf

    integer :: ierror

    call mpi_comm_rank(pf_comm%comm, pf%rank, ierror)

    allocate(pf_comm%recvreq(pf%nlevels))
    allocate(pf_comm%sendreq(pf%nlevels))
  end subroutine pf_comm_setup

  ! Destroy
  subroutine pf_comm_destroy(pf_comm)
    type(pf_comm_t), intent(inout) :: pf_comm

    deallocate(pf_comm%recvreq)
    deallocate(pf_comm%sendreq)
  end subroutine pf_comm_destroy

  ! Post
  subroutine post(pf, level, tag)
    use pf_mod_mpi, only: MPI_REAL8

    type(pf_pfasst_t), intent(in)    :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag

    integer :: ierror

    if (pf%comm%nproc > 1 .and. pf%rank > 0) then
       call mpi_irecv(level%recv, level%nvars, MPI_REAL8, &
            pf%rank-1, tag, pf%comm%comm, pf%comm%recvreq(level%level), ierror)
    end if
  end subroutine post

  ! Receive
  subroutine recv(pf, level, tag, blocking)
    use pf_mod_mpi, only: MPI_REAL8, MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking

    integer :: ierror, stat(MPI_STATUS_SIZE)

    if (pf%rank > 0) then
       if (blocking) then
          call mpi_recv(level%recv, level%nvars, MPI_REAL8, &
               pf%rank-1, tag, pf%comm%comm, stat, ierror)
          if (ierror .ne. 0) then
             print *, 'WARNING: MPI ERROR DURING RECEIVE', ierror
          end if
       else
          call mpi_wait(pf%comm%recvreq(level%level), stat, ierror)
       end if

       level%q0 = level%recv
    end if
  end subroutine recv

  ! Send
  subroutine send(pf, level, tag, blocking)
    use pf_mod_mpi, only: MPI_REAL8, MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking

    integer :: ierror, stat(MPI_STATUS_SIZE)

    if (pf%rank < pf%comm%nproc-1) then
       level%send = level%qend

       if (blocking) then
          call mpi_send(level%send, level%nvars, MPI_REAL8, &
               pf%rank+1, tag, pf%comm%comm, stat, ierror)
       else
          ! XXX: check previous status
          call mpi_isend(level%send, level%nvars, MPI_REAL8, &
               pf%rank+1, tag, pf%comm%comm, pf%comm%sendreq(level%level), ierror)
       end if
    end if
  end subroutine send

  ! Broadcast
  subroutine broadcast(pf, y, nvar, root)
    use pf_mod_mpi, only: MPI_REAL8

    type(pf_pfasst_t), intent(inout) :: pf
    real(kind=8),      intent(in)    :: y(nvar)
    integer,           intent(in)    :: nvar, root

    integer :: ierror

    call mpi_bcast(y, nvar, MPI_REAL8, root, pf%comm%comm, ierror)
  end subroutine broadcast

end module pf_mod_comm
