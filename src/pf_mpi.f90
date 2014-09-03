!
! Copyright (C) 2014 Matthew Emmett and Michael Minion.
!

module pf_mod_mpi
  include "mpif.h"
end module pf_mod_mpi

module pf_mod_comm_mpi
  use pf_mod_dtype
  implicit none

contains

  !
  ! Create an MPI based PFASST communicator using the MPI
  ! communicator *mpi_comm*.
  !
  subroutine pf_mpi_create(comm, mpi_comm)
    type(pf_comm), intent(  out) :: comm
    integer,       intent(in   ) :: mpi_comm

    integer :: ierror

    comm%comm = mpi_comm
    call mpi_comm_size(mpi_comm, comm%nproc, ierror)
  end subroutine pf_mpi_create

  !
  ! Setup the PFASST communicator.
  !
  ! This should be called soon after adding levels to the PFASST
  ! controller **pf**.
  !
  subroutine pf_mpi_setup(comm, pf)
    use pf_mod_mpi, only: MPI_REQUEST_NULL

    type(pf_comm),   intent(inout) :: comm
    type(pf_pfasst), intent(inout) :: pf

    integer :: ierror

    call mpi_comm_rank(comm%comm, pf%rank, ierror)

    allocate(comm%recvreq(pf%nlevels))
    allocate(comm%sendreq(pf%nlevels))

    comm%sendreq = MPI_REQUEST_NULL
    comm%statreq = -66
  end subroutine pf_mpi_setup

  !
  ! Destroy the PFASST communicator.
  !
  subroutine pf_mpi_destroy(comm)
    type(pf_comm), intent(inout) :: comm

    deallocate(comm%recvreq)
    deallocate(comm%sendreq)
  end subroutine pf_mpi_destroy

  !
  ! Post receive requests.
  !
  subroutine pf_mpi_post(pf, level, tag)
    use pf_mod_mpi, only: MPI_REAL8

    type(pf_pfasst), intent(in   ) :: pf
    type(pf_level),  intent(inout) :: level
    integer,         intent(in   ) :: tag

    integer :: ierror

    call mpi_irecv(level%recv, level%ndofs, MPI_REAL8, &
                   modulo(pf%rank-1, pf%comm%nproc), tag, pf%comm%comm, pf%comm%recvreq(level%level), ierror)
  end subroutine pf_mpi_post

  !
  ! Send/receive solution.
  !
  subroutine pf_mpi_send(pf, level, tag, blocking)
    use pf_mod_mpi, only: MPI_REAL8, MPI_STATUS_SIZE

    type(pf_pfasst), intent(inout) :: pf
    type(pf_level),  intent(inout) :: level
    integer,         intent(in   ) :: tag
    logical,         intent(in   ) :: blocking

    integer :: ierror, stat(MPI_STATUS_SIZE)

    if (blocking) then
       level%send = level%qend
       call mpi_send(level%send, level%ndofs, MPI_REAL8, &
                     modulo(pf%rank+1, pf%comm%nproc), tag, pf%comm%comm, stat, ierror)
    else
       call mpi_wait(pf%comm%sendreq(level%level), stat, ierror)
       level%send = level%qend
       call mpi_isend(level%send, level%ndofs, MPI_REAL8, &
                      modulo(pf%rank+1, pf%comm%nproc), tag, pf%comm%comm, pf%comm%sendreq(level%level), ierror)
    end if
  end subroutine pf_mpi_send

  subroutine pf_mpi_recv(pf, level, tag, blocking)
    use pf_mod_mpi, only: MPI_REAL8, MPI_STATUS_SIZE
    type(pf_pfasst), intent(inout) :: pf
    type(pf_level),  intent(inout) :: level
    integer,         intent(in   ) :: tag
    logical,         intent(in   ) :: blocking

    integer :: ierror, stat(MPI_STATUS_SIZE)

    if (blocking) then
       call mpi_recv(level%recv, level%ndofs, MPI_REAL8, &
                     modulo(pf%rank-1, pf%comm%nproc), tag, pf%comm%comm, stat, ierror)
    else
       call mpi_wait(pf%comm%recvreq(level%level), stat, ierror)
    end if

    if (ierror .ne. 0) &
         print *, pf%rank, 'warning: mpi error during receive', ierror

    level%q0 = level%recv
  end subroutine pf_mpi_recv

  !
  ! Misc
  !
  subroutine pf_mpi_wait(pf, level)
    use pf_mod_mpi, only: MPI_STATUS_SIZE
    type(pf_pfasst), intent(in   ) :: pf
    integer,         intent(in   ) :: level
    integer :: ierror, stat(MPI_STATUS_SIZE)
    call mpi_wait(pf%comm%sendreq(level), stat, ierror)
  end subroutine pf_mpi_wait

  subroutine pf_mpi_broadcast(pf, y, nvar, root)
    use pf_mod_mpi, only: MPI_REAL8
    type(pf_pfasst), intent(inout) :: pf
    real(pfdp),      intent(in   ) :: y(nvar)
    integer,         intent(in   ) :: nvar, root
    integer :: ierror
    call mpi_bcast(y, nvar, MPI_REAL8, root, pf%comm%comm, ierror)
  end subroutine pf_mpi_broadcast

end module pf_mod_comm_mpi
