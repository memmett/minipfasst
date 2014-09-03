!
! Copyright (C) 2012, 2013 Matthew Emmett and Michael Minion.
!
! This file is part of LIBPFASST.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.
!

! Parallel PFASST routines.

module pf_mod_parallel
  use pf_mod_dtype
  use pf_mod_interpolate
  use pf_mod_restrict
  use pf_mod_utils
  use pf_mod_hooks
  use pf_mod_comm_mpi
  use sweeper
  implicit none
contains

  !
  ! Predictor.
  !
  ! Spreads the fine initial condition (F%q0) to all levels and all
  ! nodes.  If we're running with more than one processor, performs
  ! sweeps on the coarsest level.
  !
  ! No time communication is performed during the predictor.
  !
  subroutine pf_predictor(pf, dt)
    type(pf_pfasst), intent(inout) :: pf
    real(pfdp),      intent(in   ) :: dt

    type(pf_level), pointer :: fine, crse
    integer    :: j, k, l
    real(pfdp) :: t0k

    call call_hooks(pf, 1, PF_PRE_PREDICTOR)

    call spreadq0(pf%levels(pf%nlevels), pf%t0)

    if (pf%nlevels > 1) then

       do l = pf%nlevels, 2, -1
          fine => pf%levels(l); crse => pf%levels(l-1)
          call restrict_time_space_fas(pf, pf%t0, dt, fine, crse)
          call save(crse)
          crse%q0 = crse%Q(:,1)
       end do

       crse => pf%levels(1)

       if (pf%comm%nproc > 1) then

          if (pf%Pipeline_G .and. (pf%levels(1)%nsweeps > 1)) then

             !  this is the weird choice.  we burn in without communication, then do extra sweeps
             do k = 1, pf%rank + 1
                pf%iter = -k

                ! Get new initial value (skip on first iteration)
                if (k > 1) then
                   crse%q0 = crse%qend
                   if (.not. pf%PFASST_pred) then
                      call spreadq0(crse, pf%t0)
                   end if
                end if

                call call_hooks(pf, crse%level, PF_PRE_SWEEP)
                call sweep(pf, crse, pf%t0, dt)
                ! call pf_residual(pf, crse, dt)
                call call_hooks(pf, crse%level, PF_POST_SWEEP)
             end do

             ! now we have mimicked the burn in and we must do pipe-lined sweeps
             do k = 1, crse%nsweeps-1
                pf%iter =-(pf%rank + 1) -k

                !  Get new initial conditions
                call pf_recv(pf, crse, crse%level*20000+pf%rank, .true.)
                call call_hooks(pf, crse%level, PF_PRE_SWEEP)
                call sweep(pf, crse, pf%t0, dt)
                call call_hooks(pf, crse%level, PF_POST_SWEEP)
                !  Send forward
                call pf_send(pf, crse,  crse%level*20000+pf%rank+1, .true.)
             end do
             ! call pf_residual(pf, crse, dt)

          else

             ! normal predictor burn in
             do k = 1, pf%rank + 1
                pf%iter = -k
                t0k = pf%t0-(pf%rank)*dt + (k-1)*dt

                ! get new initial value (skip on first iteration)
                if (k > 1) then
                   crse%q0 = crse%qend
                   if (.not. pf%PFASST_pred) then
                      call spreadq0(crse, t0k)
                   end if
                end if

                call call_hooks(pf, crse%level, PF_PRE_SWEEP)
                do j = 1, crse%nsweeps
                   call sweep(pf, crse, t0k, dt)
                end do
                ! call pf_residual(pf, crse, dt)
                call call_hooks(pf, crse%level, PF_POST_SWEEP)
             end do
          end if

       else

          ! Single processor... sweep on coarse and return to fine level.

          do k = 1, pf%rank + 1
             pf%iter = -k
             t0k = pf%t0-(pf%rank)*dt + (k-1)*dt

             call call_hooks(pf, crse%level, PF_PRE_SWEEP)
             do j = 1, crse%nsweeps
                call sweep(pf, crse, t0k, dt)
             end do
             ! call pf_residual(pf, crse, dt)
             call call_hooks(pf, crse%level, PF_POST_SWEEP)
          end do

       end if

       ! return to fine level

       do l = 2, pf%nlevels-1
          fine => pf%levels(l); crse => pf%levels(l-1)
          call interpolate_time_space(pf, pf%t0, dt, fine, crse, crse%Finterp)
          fine%q0 = fine%Q(:,1)

          call call_hooks(pf, l, PF_PRE_SWEEP)
          do j = 1, fine%nsweeps
             call sweep(pf, fine, pf%t0, dt)
          end do
          call call_hooks(pf, l, PF_POST_SWEEP)
       end do

       fine => pf%levels(pf%nlevels); crse => pf%levels(pf%nlevels-1)
       call interpolate_time_space(pf, pf%t0, dt, fine, crse, crse%Finterp)
       fine%q0 = fine%Q(:,1)

    end if

    call call_hooks(pf, -1, PF_POST_PREDICTOR)
  end subroutine pf_predictor

  !
  ! Run in parallel using PFASST.
  !
  subroutine pf_pfasst_run(pf, q0, dt, tend, nsteps_in)
    type(pf_pfasst), intent(inout)           :: pf
    real(pfdp),      intent(in   )           :: q0(:), dt, tend
    integer,         intent(in   ), optional :: nsteps_in

    type(pf_level), pointer :: finest

    integer :: nsteps, nblocks
    integer :: n, k, b, l

    if (present(nsteps_in)) then
       nsteps = nsteps_in
    else
       nsteps = ceiling(1.0*tend/dt)
    end if

    nblocks = nsteps / pf%comm%nproc

    pf%comm%statreq  = -66
    pf%levels(pf%nlevels)%q0 = q0
    pf%dt = dt

    finest => pf%levels(pf%nlevels)

    do b = 1, nblocks
       pf%step = (b-1) * pf%comm%nproc + pf%rank
       pf%t0   = pf%step * dt

       ! XXX: get rid of dt here...
       call pf_predictor(pf, dt)

       do k = 1, pf%niters
          pf%iter = k
          do l = 2, pf%nlevels
             call pf_post(pf, pf%levels(l), l*10000+pf%iter)
          end do
          call pf_v_cycle(pf, pf%nlevels)
       end do

       if (b < nblocks) then
          call pf_mpi_wait(pf, finest%level)
          finest%send = finest%qend
          call pf_broadcast(pf, finest%send, finest%ndofs, pf%comm%nproc-1)
          finest%q0 = finest%send
       end if
    end do
  end subroutine pf_pfasst_run

  !
  ! Execute a V-cycle
  !
  recursive subroutine pf_v_cycle(pf, level)
    type(pf_pfasst), intent(inout) :: pf
    integer,         intent(in   ) :: level

    if (level == 1) then
       call pf_v_cycle_bottom(pf, level)
    else
       call pf_v_cycle_down(pf, level)
       call pf_v_cycle(pf, level-1)
       call pf_v_cycle_up(pf, level)
    end if
  end subroutine pf_v_cycle

  recursive subroutine pf_v_cycle_bottom(pf, level)
    type(pf_pfasst), intent(inout) :: pf
    integer,         intent(in   ) :: level

    type(pf_level), pointer :: crse
    integer :: j

    crse => pf%levels(level)
    call pf_recv(pf, crse, crse%level*10000+pf%iter, .true.)
    call call_hooks(pf, crse%level, PF_PRE_SWEEP)
    do j = 1, crse%nsweeps
       call sweep(pf, crse, pf%t0, pf%dt)
    end do
    call call_hooks(pf, crse%level, PF_POST_SWEEP)
    call pf_send(pf, crse, crse%level*10000+pf%iter, .true.)
  end subroutine pf_v_cycle_bottom

  recursive subroutine pf_v_cycle_down(pf, level)
    type(pf_pfasst), intent(inout) :: pf
    integer,         intent(in   ) :: level

    type(pf_level), pointer :: crse, fine
    integer :: j

    fine => pf%levels(level)
    crse => pf%levels(level-1)

    call call_hooks(pf, fine%level, PF_PRE_SWEEP)
    do j = 1, fine%nsweeps
       call sweep(pf, fine, pf%t0, pf%dt)
    end do
    ! call pf_residual(pf, fine, pf%dt)
    call call_hooks(pf, fine%level, PF_POST_SWEEP)
    call pf_send(pf, fine, fine%level*10000+pf%iter, .false.)
    call restrict_time_space_fas(pf, pf%t0, pf%dt, fine, crse)
    call save(crse)

  end subroutine pf_v_cycle_down

  recursive subroutine pf_v_cycle_up(pf, level)
    type(pf_pfasst), intent(inout) :: pf
    integer,         intent(in   ) :: level

    type(pf_level), pointer :: crse, fine
    integer :: j

    fine => pf%levels(level)
    crse => pf%levels(level-1)

    call interpolate_time_space(pf, pf%t0, pf%dt, fine, crse, crse%Finterp)
    call pf_recv(pf, fine, fine%level*10000+pf%iter, .false.)

    if (pf%rank /= 0) then
       call interpolate_q0(pf, fine, crse)
    end if

    if (fine%level < pf%nlevels) then
       call call_hooks(pf, fine%level, PF_PRE_SWEEP)
       do j = 1, fine%nsweeps
          call sweep(pf, fine, pf%t0, pf%dt)
       end do
       ! call pf_residual(pf, F, dt)
       call call_hooks(pf, fine%level, PF_POST_SWEEP)
    end if

  end subroutine pf_v_cycle_up

  !
  ! Communication helpers
  !
  subroutine pf_post(pf, level, tag)
    type(pf_pfasst), intent(in   ) :: pf
    type(pf_level),  intent(inout) :: level
    integer,         intent(in   ) :: tag
    if (pf%rank /= 0) then
       call pf_mpi_post(pf, level, tag)
    end if
  end subroutine pf_post

  subroutine pf_send(pf, level, tag, blocking)
    type(pf_pfasst), intent(inout) :: pf
    type(pf_level),  intent(inout) :: level
    integer,         intent(in   ) :: tag
    logical,         intent(in   ) :: blocking
    if (pf%rank /= pf%comm%nproc-1) then
       call pf_mpi_send(pf, level, tag, blocking)
    end if
  end subroutine pf_send

  subroutine pf_recv(pf, level, tag, blocking)
    type(pf_pfasst), intent(inout) :: pf
    type(pf_level),  intent(inout) :: level
    integer,         intent(in   ) :: tag
    logical,         intent(in   ) :: blocking
    if (pf%rank /= 0) then
       call pf_mpi_recv(pf, level, tag, blocking)
    end if
  end subroutine pf_recv

  subroutine pf_broadcast(pf, y, nvar, root)
    type(pf_pfasst), intent(inout) :: pf
    real(pfdp)  ,    intent(in   ) :: y(nvar)
    integer,         intent(in   ) :: nvar, root
    call pf_mpi_broadcast(pf, y, nvar, root)
  end subroutine pf_broadcast

end module pf_mod_parallel
