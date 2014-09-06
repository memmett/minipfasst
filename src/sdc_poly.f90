!
! Copyright (C) 2014 Matthew Emmett.
!

!
! Polynomial manipulation routines.
!
! A polynomial p
!
!   p(x) = a_n x^n + ... + a_2 x^2 + a_1 x + a_0
!
! is stored as a Fortran array p(0:n) according to
!
!   p = [ a_0, a_1, ..., a_n ].
!

module sdc_mod_poly
  use pf_mod_config
  implicit none

  private :: qsort_partition

  interface poly_eval
     module procedure poly_eval
     module procedure poly_eval_complex
  end interface

contains

  !
  ! Evaluate polynomial.
  !
  real(pfqp) function poly_eval(p, n, x) result(v)
    integer,    intent(in   ) :: n
    real(pfqp), intent(in   ) :: p(0:n), x

    integer :: j

    v = p(n)
    do j = n-1, 0, -1
       v = x * v + p(j)
    end do
  end function

  complex(pfqp) function poly_eval_complex(p, n, x) result(v)
    integer,       intent(in   ) :: n
    real(pfqp),    intent(in   ) :: p(0:n)
    complex(pfqp), intent(in   ) :: x

    integer :: j

    v = p(n)
    do j = n-1, 0, -1
       v = x * v + p(j)
    end do
  end function


  !
  ! Differentiate polynomial (in place)
  !
  subroutine poly_diff(p, n)
    integer,    intent(in   ) :: n
    real(pfqp), intent(inout) :: p(0:n)

    integer  :: j
    real(pfqp) :: pp(0:n)

    pp = 0.0_pfqp

    do j = 1, n
       pp(j-1) = j * p(j)
    end do

    p = pp
  end subroutine poly_diff


  !
  ! Integrate polynomial (in place)
  !
  subroutine poly_int(p, n)
    integer,    intent(in   ) :: n
    real(pfqp), intent(inout) :: p(0:n)

    integer  :: j
    real(pfqp) :: pp(0:n)

    pp = 0.0_pfqp

    do j = 0, n-1
       pp(j+1) = p(j) / (j+1)
    end do

    p = pp
  end subroutine poly_int


  !
  ! Compute Legendre polynomial coefficients using Bonnet's recursion
  ! formula.
  !
  subroutine poly_legendre(p, n)
    integer,    intent(in   ) :: n
    real(pfqp), intent(  out) :: p(0:n)

    real(pfqp), dimension(0:n) :: p0, p1, p2
    integer :: j, m

    if (n == 0) then
       p = [ 1.0_pfqp ]
       return
    end if

    if (n == 1) then
       p = [ 0.0_pfqp, 1.0_pfqp ]
       return
    end if

    p0 = 0.0_pfqp; p1 = 0.0_pfqp; p2 = 0.0_pfqp

    p0(0) = 1.0_pfqp
    p1(1) = 1.0_pfqp

    ! (n + 1) P_{n+1} = (2n + 1) x P_{n} - n P_{n-1}
    do m = 1, n-1
       do j = 1, n
          p2(j) = ( (2*m + 1) * p1(j-1) - m * p0(j) ) / (m + 1)
       end do
       p2(0) = - m * p0(0) / (m + 1)

       p0 = p1
       p1 = p2
    end do

    p = p2
  end subroutine poly_legendre


  !
  ! Compute polynomial roots using the Durand-Kerner algorithm.
  !
  ! The roots are assumed to be real.
  !
  subroutine poly_roots(roots, p0, n)
    integer,    intent(in   ) :: n
    real(pfqp), intent(  out) :: roots(n)
    real(pfqp), intent(in   ) :: p0(0:n)

    integer     :: i, j, k
    complex(pfqp) :: num, den, z0(n), z1(n)
    real(pfqp)    :: p(0:n)

    real(pfqp), parameter :: eps = 0.1_pfqp**15

    p = p0 / p0(n)

    ! initial guess
    do i = 1, n
       z0(i) = (0.4_pfqp, 0.9_pfqp)**i
    end do

    ! durand-kerner-weierstrass iterations
    z1 = z0
    do k = 1, 100
       do i = 1, n

          ! evaluate poly at z0(i)
          num = poly_eval(p, n, z0(i))

          ! evaluate denominator
          den = 1.0_pfqp
          do j = 1, n
             if (j == i) cycle
             den = den * (z0(i) - z0(j))
          end do

          ! update
          z0(i) = z0(i) - num / den
       end do

       ! converged?
       if (sum(abs(z0 - z1)) < eps) exit

       z1 = z0
    end do

    ! print *, k, 'iterations used to reach', eps

    roots = real(z0)
    where (abs(roots) < eps) roots = 0.0_pfqp
    call qsort(roots)

  end subroutine poly_roots


  !
  ! Sort (inplace) using the quick sort algorithm.
  !
  ! Adapted from http://www.fortran.com/qsort_c.f95.
  !
  recursive subroutine qsort(a)
    real(pfqp), intent(inout) :: a(:)
    integer :: iq

    if (size(a) > 1) then
       call qsort_partition(a, iq)
       call qsort(a(:iq-1))
       call qsort(a(iq:))
    end if
  end subroutine qsort

  subroutine qsort_partition(a, marker)
    real(pfqp), intent(inout) :: a(:)
    integer,    intent(  out) :: marker

    integer  :: i, j
    real(pfqp) :: temp, x

    x = a(1)
    i = 0
    j = size(a) + 1

    do
       j = j-1
       do
          if (a(j) <= x) exit
          j = j-1
       end do

       i = i+1
       do
          if (a(i) >= x) exit
          i = i+1
       end do

       if (i < j) then
          temp = a(i)
          a(i) = a(j)
          a(j) = temp
       else if (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do
  end subroutine qsort_partition

end module sdc_mod_poly
