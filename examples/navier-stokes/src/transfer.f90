!
! Copyright (c) 2014, Matthew Emmett and Michael Minion.  All rights reserved.
!

! Transfer (interpolate, restrict) routines.

module transfer
  use pf_mod_dtype
  use feval, only: pack3, unpack3, fft3
  implicit none
contains

  subroutine interpolate(qF, qG, fine, crse, t)
    real(pfdp),     intent(in   ) :: t, qG(:)
    real(pfdp),     intent(  out) :: qF(:)
    type(pf_level), intent(inout) :: fine, crse

    integer :: i, j, k, ii, jj, kk, nG, nF

    complex(pfdp), dimension(fine%user%nx,fine%user%ny,fine%user%nz) :: fu, fv, fw, fuhat, fvhat, fwhat
    complex(pfdp), dimension(crse%user%nx,crse%user%ny,crse%user%nz) :: cu, cv, cw, cuhat, cvhat, cwhat

    call unpack3(qG, cu, cv, cw)
    call fft3(crse%user%fft, crse%user%wk, crse%user%scale, cu, cv, cw, cuhat, cvhat, cwhat)

    nG = crse%user%nx
    nF = fine%user%nx

    fuhat = 0
    fvhat = 0
    fwhat = 0

    do k = 1, nG
       if (k <= nG/2) then
          kk = k
       else if (k > nG/2+1) then
          kk = nF - nG + k
       else
          cycle
       end if

       do j = 1, nG
          if (j <= nG/2) then
             jj = j
          else if (j > nG/2+1) then
             jj = nF - nG + j
          else
             cycle
          end if

          do i = 1, nG
             if (i <= nG/2) then
                ii = i
             else if (i > nG/2+1) then
                ii = nF - nG + i
             else
                cycle
             end if

             fuhat(ii,jj,kk) = cuhat(i,j,k)
             fvhat(ii,jj,kk) = cvhat(i,j,k)
             fwhat(ii,jj,kk) = cwhat(i,j,k)
          end do
       end do
    end do
    call fft3(fine%user%bft, fine%user%wk, 1.d0, fuhat, fvhat, fwhat, fu, fv, fw)
    call pack3(qF, fu, fv, fw)
  end subroutine interpolate

  subroutine restrict(qF, qG, fine, crse, t)
    real(pfdp),     intent(in   ) :: t, qF(:)
    real(pfdp),     intent(  out) :: qG(:)
    type(pf_level), intent(in   ) :: fine, crse

    integer :: i, j, k

    complex(pfdp), dimension(fine%user%nx,fine%user%ny,fine%user%nz) :: fu, fv, fw
    complex(pfdp), dimension(crse%user%nx,crse%user%ny,crse%user%nz) :: cu, cv, cw

    call unpack3(qF, fu, fv, fw)
    do k = 1, crse%user%nz
       do j = 1, crse%user%ny
          do i = 1, crse%user%nx
             cu(i,j,k) = fu(2*(i-1)+1,2*(j-1)+1,2*(k-1)+1)
             cv(i,j,k) = fv(2*(i-1)+1,2*(j-1)+1,2*(k-1)+1)
             cw(i,j,k) = fw(2*(i-1)+1,2*(j-1)+1,2*(k-1)+1)
          end do
       end do
    end do
    call pack3(qG, cu, cv, cw)
  end subroutine restrict
end module transfer
