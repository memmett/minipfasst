!
! Copyright (C) 2014 Matthew Emmett and Michael Minion.
!

module pf_mod_utils
  use pf_mod_dtype
  implicit none
contains

  !
  ! Build time interpolation matrix.
  !
  subroutine pf_time_interpolation_matrix(nodesF, nnodesF, nodesG, nnodesG, tmat)
    integer,    intent(in)  :: nnodesF, nnodesG
    real(pfdp), intent(in)  :: nodesF(0:nnodesF-1), nodesG(0:nnodesG-1)
    real(pfdp), intent(out) :: tmat(0:nnodesF-1,0:nnodesG-1)

    integer    :: i, j, k
    real(pfdp) :: xi, num, den

    do i = 0, nnodesF-1
       xi = nodesF(i)

       do j = 0, nnodesG-1
          den = 1.0_pfdp
          num = 1.0_pfdp

          do k = 0, nnodesG-1
             if (k == j) cycle
             den = den * (nodesG(j) - nodesG(k))
             num = num * (xi        - nodesG(k))
          end do

          tmat(i, j) = num/den
       end do
    end do
  end subroutine pf_time_interpolation_matrix


  !
  ! Spread initial condition.
  !
  subroutine spreadq0(lev, t0)
    use sweeper, only: evaluate

    type(pf_level), intent(inout) :: lev
    real(pfdp),     intent(in   ) :: t0

    integer :: m, p

    lev%q(:,1) = lev%q0
    call evaluate(lev, t0, 1)

    do m = 2, lev%nnodes
       lev%q(:,m) = lev%q(:,1)
       do p = 1, size(lev%f, dim=3)
          lev%f(:,m,p) = lev%f(:,1,p)
       end do
    end do
  end subroutine spreadq0


  !
  ! Save current Q and F.
  !
  subroutine save(lev)
    type(pf_level), intent(inout) :: lev

    if (lev%Finterp) then
       if (allocated(lev%pF)) then
          lev%pf = lev%f
       end if
    else
       if (allocated(lev%pQ)) then
          lev%pq = lev%q
       end if
    end if
  end subroutine save


  !
  ! Apply an interpolation matrix (tmat or rmat) to src.
  !
  subroutine pf_apply_mat(dst, a, mat, src, zero)
    real(pfdp), intent(  out)           :: dst(:, :)
    real(pfdp), intent(in   )           :: a, mat(:, :), src(:,:)
    logical,    intent(in   ), optional :: zero

    logical :: lzero
    integer :: n, m, i, j

    lzero = .true.; if (present(zero)) lzero = zero

    n = size(mat, dim=1)
    m = size(mat, dim=2)

    do i = 1, n
       if (lzero) dst(:,i) = 0
       do j = 1, m
          dst(:,i) = dst(:,i) + a * mat(i, j) * src(:,j)
       end do
    end do
  end subroutine pf_apply_mat

end module pf_mod_utils
