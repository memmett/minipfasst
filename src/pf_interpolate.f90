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

module pf_mod_interpolate
  use pf_mod_dtype
  use pf_mod_restrict
  use pf_mod_hooks
  use pf_mod_utils
  use sweeper
  use transfer
  implicit none
contains

  !
  ! Interpolate (in time and space) levG to levF.
  !
  ! Interpolation is done by interpolating increments.  The fine
  ! function values are re-evaluated after interpolation.
  !
  subroutine interpolate_time_space(pf, t0, dt, levF, levG, Finterp)
    type(pf_pfasst), intent(inout)           :: pf
    real(pfdp),      intent(in   )           :: t0, dt
    type(pf_level),  intent(inout)           :: levF, levG
    logical,         intent(in   ), optional :: Finterp !  if true, then do interp on f not q

    integer    :: m, p
    real(pfdp) :: tF(levF%nnodes)
    real(pfdp) :: tG(levG%nnodes)
    logical    :: Finterp_loc

    real(pfdp) ::  delG(levG%ndofs,levG%nnodes)   ! coarse in time and space
    real(pfdp) ::  delGF(levF%ndofs,levG%nnodes)  ! coarse in time but fine in space

    call call_hooks(pf, levF%level, PF_PRE_INTERP_ALL)

    ! set time at coarse and fine nodes
    tG = t0 + dt * levG%nodes
    tF = t0 + dt * levF%nodes

    delG = levG%Q - levG%pQ

    ! interpolate q
    do m = 1, levG%nnodes
       call interpolate(delGF(:,m), delG(:,m), levF, levG, tG(m))
    end do

    ! interpolate corrections
    call pf_apply_mat(levF%Q, 1.0_pfdp, levF%tmat, delGF, .false.)

    Finterp_loc = .false.
    if(present(Finterp)) then
       if  (Finterp)  then
          Finterp_loc = .true.
       end if
    end if

    if (Finterp_loc) then
       ! interpolate f
       do p = 1, size(levG%F, dim=3)
          delG = levG%F(:,:,p) - levG%pF(:,:,p)

          do m = 1, levG%nnodes
            call interpolate(delGF(:,m), delG(:,m), levF, levG, tG(m))
         end do

         ! interpolate corrections  in time
         call pf_apply_mat(levF%F(:,:,p), 1.0_pfdp, levF%tmat, delGF, .false.)

       end do
    else
       ! recompute fs
       do m = 1, levF%nnodes
          call evaluate(levF, tF(m), m)
       end do
    end if

    ! reset qend so that it is up to date
    levF%qend = levF%Q(:,levF%nnodes)

    call call_hooks(pf, levF%level, PF_POST_INTERP_ALL)
  end subroutine interpolate_time_space

  subroutine interpolate_q0(pf, levF, levG)
    !  Use to update the fine initial condition from increment

    type(pf_pfasst), intent(inout) :: pf
    type(pf_level),  intent(inout) :: levF, levG

    real(pfdp) :: delG(levG%ndofs), delF(levF%ndofs)

    call call_hooks(pf, levF%level, PF_PRE_INTERP_Q0)

    call restrict(levF%q0, delG, levF, levG, pf%t0)
    delG = delG - levG%q0

    call interpolate(delF, delG, levF, levG, pf%t0)
    levF%q0 = levF%q0 + delF

    call call_hooks(pf, levF%level, PF_POST_INTERP_Q0)
  end subroutine interpolate_q0

end module pf_mod_interpolate
