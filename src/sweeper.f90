!
! Copyright (C) 2014 Matthew Emmett.
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

module sweeper
  use pf_mod_dtype
  use feval
  implicit none

contains

  ! Perform on SDC sweep on level lev and set qend appropriately.
  subroutine sweep(pf, lev, t0, dt)
    type(pf_pfasst), intent(inout) :: pf
    real(pfdp),      intent(in   ) :: dt, t0
    type(pf_level),  intent(inout) :: lev

    integer     :: m, n
    real(pfdp)  :: t, dtsdc(1:lev%nnodes-1), rhs(lev%ndofs), S(lev%ndofs,lev%nnodes-1)


    ! compute integrals and add fas correction
    S = 0
    do m = 1, lev%nnodes-1
       S(:,m) = S(:,m) + dt * lev%sweeper%expl_mat(m,n) * lev%F(:,n,1)
       S(:,m) = S(:,m) + dt * lev%sweeper%impl_mat(m,n) * lev%F(:,n,2)
    end do

    if (allocated(lev%tau)) then
       S = S + lev%tau
    end if

    ! do the time-stepping
    lev%Q(:,1) = lev%q0

    call f1eval(lev%Q(:,1), t0, lev, lev%F(:,1,1))
    call f2eval(lev%Q(:,1), t0, lev, lev%F(:,1,2))

    t = t0
    dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))
    do m = 1, lev%nnodes-1
       t = t + dtsdc(m)

       rhs = lev%Q(:,m) + dtsdc(m) * lev%F(:,m,1) + S(:,m)

       call f2comp(lev%Q(:,m+1), t, dtsdc(m), rhs, lev, lev%F(:,m+1,2))
       call f1eval(lev%Q(:,m+1), t, lev, lev%F(:,m+1,1))
    end do

    lev%qend = lev%Q(:,lev%nnodes)
  end subroutine sweep

  ! Evaluate function values
  subroutine evaluate(lev, t, m)
    real(pfdp),       intent(in   ) :: t
    integer,          intent(in   ) :: m
    type(pf_level), intent(inout) :: lev

    call f1eval(lev%Q(:,m), t, lev, lev%F(:,m,1))
    call f2eval(lev%Q(:,m), t, lev, lev%F(:,m,2))
  end subroutine evaluate

  ! Initialize matrices
  subroutine initialize(lev)
    type(pf_level), intent(inout) :: lev

    real(pfdp) :: dsdc(lev%nnodes-1)
    integer    :: m, nnodes

    nnodes = lev%nnodes
    allocate(lev%sweeper%expl_mat(nnodes-1,nnodes))  !  S-FE
    allocate(lev%sweeper%impl_mat(nnodes-1,nnodes))  !  S-BE

    lev%sweeper%expl_mat = lev%smat
    lev%sweeper%impl_mat = lev%smat

    dsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       lev%sweeper%expl_mat(m,m)   = lev%sweeper%expl_mat(m,m)   - dsdc(m)
       lev%sweeper%impl_mat(m,m+1) = lev%sweeper%impl_mat(m,m+1) - dsdc(m)
    end do
  end subroutine initialize

  ! Compute SDC integral
  subroutine integrate(lev, qSDC, fSDC, dt, fintSDC)
    type(pf_level), intent(in   ) :: lev
    real(pfdp),     intent(in   ) :: qSDC(:,:), fSDC(:,:,:)
    real(pfdp),     intent(in   ) :: dt
    real(pfdp),     intent(inout) :: fintSDC(:,:)

    integer :: n, m, p

    fintSDC = 0

    do n = 1, lev%nnodes-1
       do m = 1, lev%nnodes
          do p = 1, npieces
             fintSDC(:,n) = fintSDC(:,n) + dt * lev%smat(n,m) * fSDC(:,m,p)
          end do
       end do
    end do
  end subroutine integrate

  ! subroutine pf_imex_destroy(sweeper)
  !   type(pf_sweeper_t), intent(inout) :: sweeper
  !   type(pf_imex_t), pointer :: imex
  !   call c_f_pointer(sweeper%sweeperctx, imex)
  !   deallocate(lev%sweeper%impl_mat)
  !   deallocate(lev%sweeper%expl_mat)
  ! end subroutine pf_imex_destroy

end module sweeper
