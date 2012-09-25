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

module pf_mod_pfasst
  use pf_mod_dtype

  implicit none

  interface create
     module procedure pf_pfasst_create
  end interface create

  interface setup
     module procedure pf_pfasst_setup
  end interface setup

  interface destroy
     module procedure pf_pfasst_destroy
  end interface destroy

contains

  ! Create a PFASST object
  !
  ! Passing the 'nvars' and 'nnodes' arrays is optional, but these
  ! should be set appropriately before calling setup.
  subroutine pf_pfasst_create(pf, pf_comm, nlevels, nvars, nnodes)
    type(pf_pfasst_t), intent(inout)         :: pf
    type(pf_comm_t),   intent(inout), target :: pf_comm
    integer, intent(in)                      :: nlevels
    integer, intent(in), optional            :: nvars(nlevels), nnodes(nlevels)

    integer :: l

    pf%comm => pf_comm

    pf%nlevels = nlevels
    allocate(pf%levels(nlevels))
    
    if (present(nvars)) then
       do l = 1, nlevels
          pf%levels(l)%nvars = nvars(l)
       end do
    end if

    if (present(nnodes)) then
       do l = 1, nlevels
          pf%levels(l)%nnodes = nnodes(l)
       end do
    end if
  end subroutine pf_pfasst_create

  ! Setup (allocate) PFASST object
  subroutine pf_pfasst_setup(pf)
    use pf_mod_quadrature
    use pf_mod_sweep
    use pf_mod_utils

    type(pf_pfasst_t), intent(inout)        :: pf

    integer :: l, nvars, nnodes
    type(pf_level_t), pointer :: F, G

    if (pf%rank < 0) then
       stop 'invalid PF rank: did you call setup correctly?'
    end if

    do l = 1, pf%nlevels
       F => pf%levels(l)

       F%level = l
       nvars   = F%nvars
       nnodes  = F%nnodes

       allocate(F%q0(nvars))
       allocate(F%qend(nvars))
       allocate(F%send(nvars))
       allocate(F%recv(nvars))
       allocate(F%qSDC(nvars,nnodes))
       allocate(F%fSDC(nvars,nnodes,npieces))
       allocate(F%nodes(nnodes))
       allocate(F%qmat(nnodes-1,nnodes))

       if (l > 1) then
          allocate(F%tau(F%nvars, nnodes-1))
       else
          nullify(F%tau)
       end if

       call pf_quadrature(pf%qtype, nnodes, pf%levels(1)%nnodes, F%nodes, F%qmat)
       call sdcinit(F)

       F%allocated = .true.
    end do

    do l = 1, pf%nlevels-1
       F => pf%levels(l)
       G => pf%levels(l+1)

       allocate(F%tmat(G%nnodes,F%nnodes-G%nnodes))
       call pf_time_interpolation_matrix(F%nodes, F%nnodes, &
            G%nodes, G%nnodes, F%tmat)
    end do
  end subroutine pf_pfasst_setup

  ! Deallocate PFASST object
  subroutine pf_pfasst_destroy(pf)
    use pf_mod_sweep, only: npieces
    type(pf_pfasst_t), intent(inout) :: pf

    integer :: l
    type(pf_level_t), pointer :: level

    do l = 1, pf%nlevels
       level => pf%levels(l)

       if (level%allocated) then
          deallocate(level%q0)
          deallocate(level%qend)
          deallocate(level%send)
          deallocate(level%recv)
          deallocate(level%nodes)
          deallocate(level%qmat)

          if (associated(level%smat)) then
             deallocate(level%smat)
          end if

          deallocate(level%qSDC)
          deallocate(level%fSDC)

          if (l > 1) then
             deallocate(level%tau)
          end if

          if (l < pf%nlevels) &
               deallocate(level%tmat)

          level%allocated = .false.
       end if
    end do

    deallocate(pf%levels)
  end subroutine pf_pfasst_destroy

end module pf_mod_pfasst
