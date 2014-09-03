!
! Copyright (c) 2014, Matthew Emmett and Michael Minion.  All rights reserved.
!

! Transfer (interpolate, restrict) routines.

module transfer
  use feval
  implicit none
contains

  subroutine interpolate(qF, qG, fine, crse, t)
    real(pfdp),     intent(in   ) :: t, qG(:)
    real(pfdp),     intent(  out) :: qF(:)
    type(pf_level), intent(in   ) :: fine, crse

    complex(kind=8), pointer :: wkF(:), wkG(:)

    integer :: nvarF, nvarG, xrat

    nvarF = size(qF)
    nvarG = size(qG)
    xrat  = nvarF / nvarG

    if (xrat == 1) then
       qF = qG
       return
    endif

    wkF => fine%user%wk
    wkG => crse%user%wk

    wkG = qG
    call fftw_execute_dft(crse%user%ffft, wkG, wkG)
    wkG = wkG / nvarG

    wkF = 0.0d0
    wkF(1:nvarG/2) = wkG(1:nvarG/2)
    wkF(nvarF-nvarG/2+2:nvarF) = wkG(nvarG/2+2:nvarG)

    call fftw_execute_dft(fine%user%ifft, wkF, wkF)

    qF = real(wkF)
  end subroutine interpolate

  subroutine restrict(qF, qG, fine, crse, t)
    real(pfdp),     intent(in   ) :: t, qF(:)
    real(pfdp),     intent(  out) :: qG(:)
    type(pf_level), intent(in   ) :: fine, crse

    integer :: nvarF, nvarG, xrat

    nvarF = size(qF)
    nvarG = size(qG)
    xrat  = nvarF / nvarG

    qG = qF(::xrat)
  end subroutine restrict
end module transfer
