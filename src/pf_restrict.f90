!
! Copyright (C) 2014 Matthew Emmett and Michael Minion.
!

!
! Restriction and FAS routines.
!

module pf_mod_restrict
  use pf_mod_dtype
  use pf_mod_hooks
  use sweeper
  use transfer
  implicit none
contains


  !
  ! Restrict (in time and space) qF to qG.
  !
  subroutine restrict_sdc(fine, crse, qF, qG, integral, tF)
    use pf_mod_utils, only: pf_apply_mat

    type(pf_level), intent(inout) :: fine, crse
    real(pfdp),     intent(inout) :: qF(:,:), qG(:,:)
    logical,        intent(in   ) :: integral
    real(pfdp),     intent(in   ) :: tF(:)

    real(pfdp) :: qFr(crse%ndofs,fine%nnodes)
    integer :: m

    if (integral) then

       do m = 1, fine%nnodes-1
          call restrict(qF(:,m), qFr(:,m), fine, crse, tF(m))
       end do

       ! when restricting '0 to node' integral terms, skip the first
       ! entry since it is zero
       call pf_apply_mat(qG, 1.d0, fine%rmat(2:,2:), qFr)

    else

       do m = 1, fine%nnodes
          call restrict(qF(:,m), qFr(:,m), fine, crse, tF(m))
       end do

       call pf_apply_mat(qG, 1.d0, fine%rmat, qFr)

    end if

  end subroutine restrict_sdc


  !
  ! Restrict (in time and space) F to G and set G's FAS correction.
  !
  ! The coarse function values are re-evaluated after restriction.
  ! Note that even if the number of variables and nodes is the same,
  ! we should still compute the FAS correction since the function
  ! evaluations may be different.
  !
  subroutine restrict_time_space_fas(pf, t0, dt, fine, crse)
    type(pf_pfasst), intent(inout) :: pf
    real(pfdp),      intent(in   ) :: t0, dt
    type(pf_level),  intent(inout) :: fine, crse

    integer    :: m
    real(pfdp) :: tG(crse%nnodes)
    real(pfdp) :: tF(fine%nnodes)
    real(pfdp) :: &
         tmpG(crse%ndofs,crse%nnodes), &    ! coarse integral of coarse function values
         tmpF(fine%ndofs,fine%nnodes), &    ! fine integral of fine function values
         tmpFr(crse%ndofs,crse%nnodes)      ! coarse integral of restricted fine function values

    call call_hooks(pf, fine%level, PF_PRE_RESTRICT_ALL)

    !
    ! restrict q's and recompute f's
    !
    tG = t0 + dt * crse%nodes
    tF = t0 + dt * fine%nodes

    call restrict_sdc(fine, crse, fine%Q, crse%Q, .false., tF)

    do m = 1, crse%nnodes
       call evaluate(crse, tG(m), m)
    end do

    !
    ! fas correction
    !
    crse%tau = 0

    if (pf%iter >= pf%taui0)  then

       ! compute '0 to node' integral on the coarse level
       call integrate(crse, crse%Q, crse%F, dt, tmpG)
       do m = 2, crse%nnodes-1
          tmpG(:,m) = tmpG(:,m) + tmpG(:,m-1)
       end do

       ! compute '0 to node' integral on the fine level
       call integrate(fine, fine%Q, fine%F, dt, tmpF)
       if (allocated(fine%tau)) tmpF = tmpF + fine%tau
       do m = 2, fine%nnodes-1
          tmpF(:,m) = tmpF(:,m) + tmpF(:,m-1)
       end do

       ! restrict '0 to node' integral on the fine level
       call restrict_sdc(fine, crse, tmpF, tmpFr, .true., tF)

       ! compute 'node to node' tau correction
       crse%tau(:,1) = tmpFr(:,1) - tmpG(:,1)
       do m = 2, crse%nnodes-1
          crse%tau(:,m) = tmpFr(:,m) - tmpFr(:,m-1) - ( tmpG(:,m) - tmpG(:,m-1) )
       end do
    end if

    call call_hooks(pf, fine%level, PF_POST_RESTRICT_ALL)

  end subroutine restrict_time_space_fas

end module pf_mod_restrict
