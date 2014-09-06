!
! Copyright (C) 2014 Matthew Emmett.
!

!
! Routines to compute various double precision SDC nodes and matrices.
!
! Computations are performed in quad precision.
!

module sdc_mod_quadrature
  use sdc_mod_poly
  use pf_mod_dtype
  implicit none
contains

  logical function not_proper(flags, node)
    integer, intent(in   ) :: flags(:)
    integer, intent(in   ) :: node

    not_proper = .not. btest(flags(node), 0)
  end function not_proper


  !
  ! Compute high precision quadrature nodes.
  !
  subroutine sdc_qnodes(qnodes, flags, qtype, nnodes)
    integer,    intent(in   ) :: nnodes
    integer,    intent(in   ) :: qtype
    real(pfqp), intent(  out) :: qnodes(nnodes)
    integer,    intent(  out) :: flags(nnodes)

    integer :: j, degree
    real(pfqp), allocatable :: roots(:)
    real(pfqp), allocatable :: coeffs(:), coeffs2(:)

    real(pfqp), parameter :: pi = 3.141592653589793115997963468544185161590576171875_pfqp

    flags = 0

    select case(qtype)

    case (SDC_GAUSS_LEGENDRE)

       degree = nnodes - 2
       allocate(roots(degree))
       allocate(coeffs(degree+1))

       call poly_legendre(coeffs, degree)
       call poly_roots(roots, coeffs, degree)

       qnodes(1)          = 0.0_pfqp
       qnodes(2:nnodes-1) = 0.5_pfqp * (1.0_pfqp + roots)
       qnodes(nnodes)     = 1.0_pfqp

       deallocate(coeffs)
       deallocate(roots)

       do j = 2, nnodes-1
          flags(j) = ibset(flags(j), 0)
       end do

    case (SDC_GAUSS_LOBATTO)

       degree = nnodes - 1
       allocate(roots(degree-1))
       allocate(coeffs(degree+1))

       call poly_legendre(coeffs, degree)
       call poly_diff(coeffs, degree)
       call poly_roots(roots, coeffs(:degree), degree-1)

       qnodes(1)          = 0.0_pfqp
       qnodes(2:nnodes-1) = 0.5_pfqp * (1.0_pfqp + roots)
       qnodes(nnodes)     = 1.0_pfqp

       deallocate(coeffs)
       deallocate(roots)

       do j = 1, nnodes
          flags(j) = ibset(flags(j), 0)
       end do

    case (SDC_GAUSS_RADAU)

       degree = nnodes - 1
       allocate(roots(degree))
       allocate(coeffs(degree+1))
       allocate(coeffs2(degree))

       call poly_legendre(coeffs, degree)
       call poly_legendre(coeffs2, degree-1)
       coeffs(:degree) = coeffs(:degree) + coeffs2
       call poly_roots(roots, coeffs, degree)

       qnodes(1)      = 0.0_pfqp
       do j = 2, nnodes-1
          qnodes(j) = 0.5_pfqp * (1.0_pfqp - roots(nnodes+1-j))
       end do
       qnodes(nnodes) = 1.0_pfqp

       deallocate(coeffs2)
       deallocate(coeffs)
       deallocate(roots)

       do j = 2, nnodes
          flags(j) = ibset(flags(j), 0)
       end do

    case (SDC_CLENSHAW_CURTIS)

       do j = 0, nnodes-1
          qnodes(j+1) = 0.5_pfqp * (1.0_pfqp - cos(j * pi / (nnodes-1)))
       end do

       do j = 1, nnodes
          flags(j) = ibset(flags(j), 0)
       end do

    case (SDC_UNIFORM)

       do j = 0, nnodes-1
          qnodes(j+1) = j * (1.0_pfqp / (nnodes-1))
       end do

       do j = 1, nnodes
          flags(j) = ibset(flags(j), 0)
       end do

    case default

       stop "ERROR: Invalid qtype in sdc_quadrature.f90."

    end select

  end subroutine sdc_qnodes


  !
  ! Compute SDC Q and S matrices given source and destination nodes.
  !
  subroutine sdc_qmats(qmat, smat, dst, src, flags, ndst, nsrc)
    integer,    intent(in   ) :: ndst, nsrc
    real(pfqp), intent(in   ) :: dst(ndst), src(nsrc)
    real(pfdp), intent(  out) :: qmat(ndst-1, nsrc), smat(ndst-1, nsrc)
    integer,    intent(in   ) :: flags(nsrc)

    integer  :: i, j, m
    real(pfqp) :: q, s, den, p(0:nsrc)

    qmat = 0.0_pfdp
    smat = 0.0_pfdp

    ! construct qmat and smat
    do i = 1, nsrc

       if (not_proper(flags, i)) cycle

       ! construct interpolating polynomial coefficients
       p    = 0.0_pfqp
       p(0) = 1.0_pfqp
       do m = 1, nsrc
          if (not_proper(flags, m) .or. m == i) cycle
          p = eoshift(p, -1) - src(m) * p
       end do

       den = poly_eval(p, nsrc, src(i))

       call poly_int(p, nsrc)

       ! evaluate integrals
       do j = 2, ndst
          q = poly_eval(p, nsrc, dst(j)) - poly_eval(p, nsrc,   0.0_pfqp)
          s = poly_eval(p, nsrc, dst(j)) - poly_eval(p, nsrc, dst(j-1))

          qmat(j-1,i) = real(q / den, pfdp)
          smat(j-1,i) = real(s / den, pfdp)
       end do
    end do

  end subroutine sdc_qmats

end module sdc_mod_quadrature
