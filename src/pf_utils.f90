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

module pf_mod_utils
  use pf_mod_dtype
  implicit none
contains

  ! Build the time interpolation matrix
  subroutine pf_time_interpolation_matrix(nodesF, nnodesF, nodesG, nnodesG, &
       tmat)
    implicit none
    integer,          intent(in)  :: nnodesF, nnodesG
    double precision, intent(in)  :: nodesF(0:nnodesF-1), nodesG(0:nnodesG-1)
    double precision, intent(out) :: tmat(0:nnodesG-1,0:nnodesF-nnodesG-1)

    integer :: i, j, k
    double precision :: xi, num, den

    if ((nnodesF-1)/(nnodesG-1) .ne. 2) then
       stop 'time ratio must be 2 when building time interpolation matrix'
    end if

    do i = 0, nnodesF-nnodesG-1

       xi = nodesF(i*2+1)

       do j = 0, nnodesG-1

          den = 1.0d0
          num = 1.0d0

          do k = 0, nnodesG-1
             if (k == j) cycle
             den = den * (nodesG(j) - nodesG(k))
             num = num * (xi        - nodesG(k))
          end do

          ! store in transpose order
          tmat(j, i) = num/den

       end do

    end do

  end subroutine pf_time_interpolation_matrix

  ! Spread initial condition
  subroutine spreadq0(F, t0)
    !use feval, only: eval_f1, eval_f2
    use pf_mod_sweep
    type(pf_level_t), intent(inout) :: F
    double precision, intent(in)    :: t0

    integer :: m, p

    F%qSDC(:, 1) = F%q0
    call sdceval(t0, 1, F)

    do m = 2, F%nnodes
       f%qSDC(:, m) = F%qSDC(:, 1)
       do p = 1, npieces
          F%fSDC(:, m, p) = F%fSDC(:, 1, p)
       end do
    end do

  end subroutine spreadq0

  subroutine echo_error(pf, level, state)
    use feval, only: exact
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state

    real(c_double) :: yexact(level%nvars)

    call exact(state%t0+state%dt, yexact)
    print '("error: step: ",i3," iter: ",i3," error: ",es14.7)', &
         state%step+1, state%iter, maxval(abs(level%qSDC(:, level%nnodes)-yexact))
  end subroutine echo_error

end module pf_mod_utils
