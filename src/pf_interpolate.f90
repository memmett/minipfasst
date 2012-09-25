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

module pf_mod_interpolate
  implicit none
contains

  ! Interpolate (in time and space) G to F.
  !
  ! Interpolation is done by interpolating increments.  The fine
  ! function values are re-evaluated after interpolation.
  subroutine interpolate_time_space(pf, t0, dt, F, G)
    use pf_mod_dtype
    use pf_mod_restrict
    use pf_mod_sweep
    use transfer

    type(pf_pfasst_t), intent(inout) :: pf
    double precision,  intent(in)    :: t0, dt
    type(pf_level_t),  intent(inout) :: F, G

    integer          :: m, n, trat
    double precision :: tm(F%nnodes)
    double precision :: &
         qSDCFx(G%nvars, G%nnodes), &
         delG(G%nvars, G%nnodes), &
         delGF(F%nvars, G%nnodes)
    
    ! interpolate increments
    call restrict_time_space(F%qSDC, qSDCFx, F%nnodes, G%nnodes, F, G)

    do m = 1, G%nnodes
       delG(:, m) = G%qSDC(:, m) - qSDCFx(:, m)
       delGF(:, m) = 0.0d0
       call interpolate(delGF(:, m), delG(:, m), F%level, G%level)
    end do

    trat = (F%nnodes-1) / (G%nnodes-1)
    if (trat .ne. 2) stop

    do m = 1, G%nnodes
       F%qSDC(:, trat*m-1) = F%qSDC(:, trat*m-1) + delGF(:, m)
    end do

    if (F%nnodes > G%nnodes) then
       do n = 1, F%nnodes - G%nnodes
          do m = 1, G%nnodes
             F%qSDC(:, 2*n) = F%qSDC(:, 2*n) + F%tmat(m, n) * delGF(:, m)
          end do
       end do
    end if

    ! recompute fs
    tm = t0 + dt*F%nodes
    do m = 1, F%nnodes
       call sdceval(tm(m), m, F)
    end do
  end subroutine interpolate_time_space

end module pf_mod_interpolate
