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

! Parallel PFASST routines.

module pf_mod_parallel
  implicit none
contains
  
  ! Run in parallel using PFASST.
  subroutine pfasst_run(pf, q0, dt, tend, nsteps, qend)
    use pf_mod_comm
    use pf_mod_interpolate
    use pf_mod_restrict
    use pf_mod_sweep
    use pf_mod_utils
    use transfer, only: interpolate

    type(pf_pfasst_t), intent(inout) :: pf
    double precision,  intent(in)    :: q0(:)
    double precision,  intent(in)    :: dt, tend
    double precision,  intent(inout), optional :: qend(:)
    integer,           intent(in), optional    :: nsteps

    double precision :: t0
    integer    :: nblock, b, k, j, l
    type(pf_level_t), pointer :: F, G

    pf%state%t0   = 0.0d0
    pf%state%dt   = dt
    pf%state%iter = -1

    !!!! set initial conditions
    F => pf%levels(1)
    F%q0 = q0

    if(present(nsteps)) then
       pf%state%nsteps = nsteps
    else
       pf%state%nsteps = ceiling(1.0*tend/dt)
    end if

    nblock = pf%state%nsteps/pf%comm%nproc

    !!!! time "block" loop
    do b = 1, nblock
       pf%state%step = pf%rank + (b-1)*pf%comm%nproc
       pf%state%t0   = pf%state%step * dt

       t0 = pf%state%t0

       !!!! predictor loop
       F => pf%levels(1)
       call spreadq0(F, t0)

          do l = 1, pf%nlevels-1
             F => pf%levels(l); G => pf%levels(l+1)
             call restrict_time_space_fas(pf, t0, dt, F, G)
             G%q0 = G%qSDC(:, 1)
          end do

       if (pf%comm%nproc > 1) then

          G => pf%levels(pf%nlevels)
          do k = 1, pf%rank + 1
             pf%state%iter = -k

             ! get new initial value (skip on first iteration)
             if (k > 1) &
                  call recv(pf, G, k-1, .true.)

             do j = 1, G%nsweeps
                call sweep(pf, t0, dt, G)
             end do
             call send(pf, G, k, .true.)
          end do

          ! interpolate
          do l = pf%nlevels, 2, -1
             F => pf%levels(l-1)
             G => pf%levels(l)
             call interpolate_time_space(pf, t0, dt, F, G)
             F%q0 = F%qSDC(:, 1)
          end do

       end if


       !!!! pfasst iterations (v-cycle)
       do k = 1, pf%niters
          pf%state%iter  = k
          pf%state%cycle = 1

          ! post receive requests
          do l = 1, pf%nlevels-1
             F => pf%levels(l)
             call post(pf, F, F%level*100+k)
          end do

          !! go down the v-cycle
          do l = 1, pf%nlevels-1
             pf%state%cycle = pf%state%cycle + 1
             F => pf%levels(l)
             G => pf%levels(l+1)

             do j = 1, F%nsweeps
                call sweep(pf, t0, dt, F)
             end do
             call send(pf, F, F%level*100+k, .false.)
             call restrict_time_space_fas(pf, t0, dt, F, G)
          end do

          !! bottom
          pf%state%cycle = pf%state%cycle + 1
          F => pf%levels(pf%nlevels)

          call recv(pf, F, F%level*100+k, .true.)
          do j = 1, F%nsweeps
             call sweep(pf, t0, dt, F)
          end do
          call send(pf, F, F%level*100+k, .true.)

          !! go up the v-cycle
          do l = pf%nlevels, 2, -1
             pf%state%cycle = pf%state%cycle + 1
             F => pf%levels(l-1)
             G => pf%levels(l)

             call interpolate_time_space(pf, t0, dt, F, G)

             call recv(pf, F, F%level*100+k, .false.)

             if (pf%rank > 0) then
                ! XXX: correction...
                ! XXX: using qSDC(1) for now but this is redundant (beginning of next sweep)
                G%qSDC(:, 1) = G%q0
                call interpolate(F%qSDC(:, 1), G%qSDC(:, 1), F%level, G%level)
                F%q0 = F%qSDC(:, 1)
             end if

             if (l > 2) then
                do j = 1, F%nsweeps
                   call sweep(pf, t0, dt, F)
                   call echo_error(pf, F, pf%state)
                end do
             end if
          end do

          F => pf%levels(1)
          call echo_error(pf, F, pf%state)
       end do

       ! broadcast fine qend (non-pipelined time loop)
       F%send = F%qend

       F => pf%levels(1)
       if (nblock > 1) &
            call broadcast(pf, F%send, F%nvars, pf%comm%nproc-1)

       F%q0 = F%send
    end do

    pf%state%iter = -1 

    if (present(qend)) then
       F => pf%levels(1)
       qend = F%qend
    end if
  end subroutine pfasst_run

end module pf_mod_parallel
