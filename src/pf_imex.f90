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

module pf_mod_sweep
  use pf_mod_dtype
  implicit none
  integer, parameter :: npieces = 2
contains

  ! Perform on SDC sweep on level F and set qend appropriately.
  subroutine sweep(pf, t0, dt, F)
    use feval, only : eval_f1, eval_f2, comp_f2

    type(pf_pfasst_t), intent(inout) :: pf
    real(c_double),    intent(in)    :: dt, t0
    type(pf_level_t),  intent(inout) :: F

    integer :: m, n
    double precision :: t, dtsdc(1:F%nnodes-1)
    double precision, pointer :: S(:,:), rhs(:)

    allocate(S(F%nvars,F%nnodes-1))
    allocate(rhs(F%nvars))

    ! compute integrals and add fas correction
    do m = 1, F%nnodes-1
       S(:, m) = 0.0d0
       do n = 1, F%nnodes
          S(:, m) = S(:, m) + dt * ( F%smat(m,n,1) * F%fSDC(:,n,1) + F%smat(m,n,2) * F%fSDC(:,n,2) )
       end do
       if (associated(F%tau)) then
          S(:, m) = S(:, m) + F%tau(:, m)
       end if
    end do

    ! do the time-stepping
    F%qSDC(:, 1) = F%q0
    call eval_f1(F%qSDC(:, 1), t0, F%level, F%fSDC(:, 1, 1))
    call eval_f2(F%qSDC(:, 1), t0, F%level, F%fSDC(:, 1, 2))

    t = t0
    dtsdc = dt * (F%nodes(2:F%nnodes) - F%nodes(1:F%nnodes-1))
    do m = 1, F%nnodes-1
       t = t + dtsdc(m)

       rhs = F%qSDC(:, m) + dtsdc(m) * F%fSDC(:, m, 1) + S(:, m)

       call comp_f2(F%qSDC(:, m+1), t, dtsdc(m), rhs, F%level, F%fSDC(:, m+1, 2))
       call eval_f1(F%qSDC(:, m+1), t, F%level, F%fSDC(:, m+1, 1))
    end do

    F%qend = F%qSDC(:, F%nnodes)

    ! done
    deallocate(rhs)
    deallocate(S)
  end subroutine sweep

  ! Evaluate function values
  subroutine sdceval(t, m, F)
    use feval, only: eval_f1, eval_f2

    double precision,       intent(in)    :: t
    integer,          intent(in)    :: m
    type(pf_level_t), intent(inout) :: F

    call eval_f1(F%qSDC(:, m), t, F%level, F%fSDC(:, m, 1))
    call eval_f2(F%qSDC(:, m), t, F%level, F%fSDC(:, m, 2))
  end subroutine sdceval

  ! Initialize smats
  subroutine sdcinit(F)
    use pf_mod_dtype
    type(pf_level_t), intent(inout) :: F

    double precision :: dsdc(F%nnodes-1)
    integer          :: m

    allocate(F%smat(F%nnodes-1, F%nnodes, npieces))

    F%smat(:,:,1) = F%qmat
    F%smat(:,:,2) = F%qmat

    dsdc = F%nodes(2:F%nnodes) - F%nodes(1:F%nnodes-1)
    do m = 1, F%nnodes-1
       F%smat(m,m,1)   = F%smat(m,m,1)   - dsdc(m)
       F%smat(m,m+1,2) = F%smat(m,m+1,2) - dsdc(m)
    end do
  end subroutine sdcinit

  ! Compute SDC integral
  subroutine sdc_integrate(qSDC,fSDC, dt, F, fintSDC)
    type(pf_level_t),  intent(in)    :: F
    double precision,  intent(in)    :: qSDC(F%nvars,F%nnodes,npieces) ! Solution
    double precision,  intent(in)    :: fSDC(F%nvars,F%nnodes,npieces) ! Function values
    double precision,  intent(inout) :: fintSDC(F%nvars,F%nnodes-1)    ! Integrals of f
    double precision,  intent(in)    :: dt

    integer :: n, m, p

    fintSDC = 0.0d0

    do n = 1, F%nnodes-1
       do m = 1, F%nnodes
          do p = 1, npieces
             fintSDC(:, n) = fintSDC(:, n) + dt * F%qmat(n, m) * fSDC(:, m, p)
          end do
       end do
    end do
  end subroutine sdc_integrate

end module pf_mod_sweep

