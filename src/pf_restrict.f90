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

! Restriction and FAS routines.

module pf_mod_restrict
  implicit none
contains

  ! Restrict (in time and space) F to G and set G's FAS correction.
  !
  ! The coarse function values are re-evaluated after restriction.
  subroutine restrict_time_space_fas(pf, t0, dt, F, G)
    use pf_mod_dtype
    use pf_mod_utils
    use pf_mod_sweep
    use transfer, only: restrict

    type(pf_pfasst_t), intent(inout) :: pf
    double precision,  intent(in)    :: t0, dt
    type(pf_level_t),  intent(inout) :: F, G

    integer :: m, mc, trat
    double precision :: tm(G%nnodes)
    double precision :: &
         CofG(G%nvars, G%nnodes-1), &    ! coarse integral of coarse function values
         FofF(F%nvars, F%nnodes-1), &    ! fine integral of fine function values
         CofF(G%nvars, G%nnodes-1), &    ! coarse integral of restricted fine function values
         tmp(G%nvars)

    trat = (F%nnodes - 1) / (G%nnodes - 1)

    ! note: even if the number of variables and nodes is the same, we
    ! should still compute the fas correction since the function
    ! evaluations may be different

    !!!! restrict qs and recompute fs
    call restrict_time_space(F%qSDC, G%qSDC, F%nnodes, G%nnodes, F, G)

    tm = t0 + dt*G%nodes
    do m = 1, G%nnodes
       call sdceval(tm(m), m, G)
    end do

    !!!! bring down fas correction from level above
    G%tau = 0.0d0

    if (associated(F%tau)) then
       ! restrict fine fas corrections and sum between coarse nodes
       do m = 1, F%nnodes-1
          mc = int(ceiling(1.0d0*m/trat))
          call restrict(F%tau(:, m), tmp, F%level, G%level)
          G%tau(:, mc) = G%tau(:, mc) + tmp
       end do
    end if

    ! fas correction
    call sdc_integrate(G%qSDC, G%fSDC, dt, G, CofG)
    call sdc_integrate(F%qSDC, F%fSDC, dt, F, FofF)

    CofF = 0.0d0

    ! restrict fine function values and sum between coarse nodes
    do m = 1, F%nnodes-1
       mc = int(ceiling(1.0d0*m/trat))
       call restrict(FofF(:, m), tmp, F%level, G%level)
       CofF(:, mc) = CofF(:, mc) + tmp
    end do

    G%tau = G%tau + CofF - CofG
  end subroutine restrict_time_space_fas


  ! Restrict qSDCF to qSDCG.
  subroutine restrict_time_space(qSDCF, qSDCG, nnodesF, nnodesG, F, G)
    use pf_mod_dtype
    use transfer, only: restrict

    integer,          intent(in)    :: nnodesF, nnodesG
    type(pf_level_t), intent(in)    :: F 
    type(pf_level_t), intent(inout) :: G
    double precision, intent(inout) :: qSDCF(F%nvars, nnodesF)
    double precision, intent(inout) :: qSDCG(G%nvars, nnodesG)

    integer :: m, trat

    if (nnodesG > 1) then
       trat = (nnodesF-1)/(nnodesG-1)
    else
       trat = 1
    end if

    do m = 1, nnodesG
       call restrict(qSDCF(:, trat*(m-1)+1), qSDCG(:, m), F%level, G%level)
    end do
  end subroutine restrict_time_space

end module pf_mod_restrict

