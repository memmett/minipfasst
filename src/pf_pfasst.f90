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

module pf_mod_pfasst
  use pf_mod_dtype
  use sweeper
  implicit none
contains

  !
  ! Create a PFASST object
  !
  subroutine pf_pfasst_create(pf, comm, nlevels, fname, nocmd)
    use pf_mod_hooks, only: PF_MAX_HOOK

    use pf_mod_options
    type(pf_pfasst),  intent(inout)           :: pf
    type(pf_comm),    intent(inout), target   :: comm
    integer,          intent(in   ), optional :: nlevels
    character(len=*), intent(in   ), optional :: fname
    logical,          intent(in   ), optional :: nocmd

    logical :: read_cmd
    integer :: l

    if (present(nlevels)) pf%nlevels = nlevels

    ! gather some input from a file and command line
    read_cmd = .true.
    if (present(nocmd)) then
         if (nocmd)  read_cmd = .false.
    end if
    if (present(fname))  then
       call pf_read_opts(pf, read_cmd, fname)
    else
       if (read_cmd)   call pf_read_opts(pf, read_cmd)
    end if

    pf%comm => comm

    allocate(pf%levels(pf%nlevels))
    allocate(pf%hooks(pf%nlevels, PF_MAX_HOOK, PF_MAX_HOOKS))
    allocate(pf%nhooks(pf%nlevels, PF_MAX_HOOK))
    pf%nhooks = 0
  end subroutine pf_pfasst_create

  !
  ! Setup PFASST object
  !
  subroutine pf_pfasst_setup(pf, ndofs, nnodes)
    use pf_mod_utils
    type(pf_pfasst), intent(inout) :: pf
    integer,         intent(in   ) :: ndofs(:), nnodes(:)

    integer :: l, nnodesF, nnodesG

    if (pf%rank < 0) then
       stop 'invalid PF rank: did you call setup correctly?'
    end if

    do l = 1, pf%nlevels
       pf%levels(l)%ndofs = ndofs(l)
       pf%levels(l)%nnodes = nnodes(l)
       pf%levels(l)%level = l
    end do

    do l = 1, pf%nlevels
       call pf_level_setup(pf, pf%levels(l))
    end do

    do l = pf%nlevels, 2, -1
       nnodesF = pf%levels(l)%nnodes
       nnodesG = pf%levels(l-1)%nnodes

       allocate(pf%levels(l)%tmat(nnodesF,nnodesG))
       call pf_time_interpolation_matrix(pf%levels(l)%nodes, nnodesF, pf%levels(l-1)%nodes, nnodesG, pf%levels(l)%tmat)

       allocate(pf%levels(l)%rmat(nnodesG,nnodesF))
       call pf_time_interpolation_matrix(pf%levels(l-1)%nodes, nnodesG, pf%levels(l)%nodes, nnodesF, pf%levels(l)%rmat)
    end do

  end subroutine pf_pfasst_setup

  !
  ! Setup (allocate) PFASST level
  !
  ! If the level is already setup, calling this again will allocate
  ! (or deallocate) tau appropriately.
  !
  subroutine pf_level_setup(pf, lev)
    use pf_mod_quadrature
    type(pf_pfasst), intent(in   ) :: pf
    type(pf_level),  intent(inout) :: lev

    integer :: m, p, ndofs, nnodes

    !
    ! do some sanity checks
    !

    if (lev%ndofs <= 0) stop "ERROR: Invalid nvars/dofs (pf_pfasst.f90)."
    if (lev%nnodes <= 0) stop "ERROR: Invalid nnodes (pf_pfasst.f90)."
    if (lev%nsweeps <= 0) stop "ERROR: Invalid nsweeps (pf_pfasst.f90)."
    ! if (.not. associated(lev%encap)) stop "ERROR: Missing encapsulation (pf_pfasst.f90)."
    ! if (.not. associated(lev%interpolate)) stop "ERROR: Missing spatial interpolation (pf_pfasst.f90)."
    ! if (.not. associated(lev%restrict)) stop "ERROR: Missing spatial restriction (pf_pfasst.f90)."

    ndofs   = lev%ndofs
    nnodes  = lev%nnodes

    ! lev%residual = -1.0_pfdp

    !
    ! allocate tau
    !
    if ((lev%level < pf%nlevels)) then
       allocate(lev%tau(ndofs,nnodes-1))
    end if

    !
    ! allocate flat buffers (q0, send, and recv)
    !
    allocate(lev%q0(ndofs))
    allocate(lev%send(ndofs))
    allocate(lev%recv(ndofs))

    !
    ! nodes, flags, and integration matrices
    !
    allocate(lev%nodes(nnodes))
    allocate(lev%nflags(nnodes))
    allocate(lev%smat(nnodes-1,nnodes))
    allocate(lev%qmat(nnodes-1,nnodes))

    if (btest(pf%qtype, 8)) then
       call pf_quadrature(pf%qtype, nnodes, pf%levels(1)%nnodes, &
            lev%nodes, lev%nflags, lev%smat, lev%qmat)
    else
       call pf_quadrature(pf%qtype, nnodes, pf%levels(pf%nlevels)%nnodes, &
            lev%nodes, lev%nflags, lev%smat, lev%qmat)
    end if

    call initialize(lev)

    !
    ! allocate Q and F
    !
    allocate(lev%Q(ndofs,nnodes))
    allocate(lev%F(ndofs,nnodes,npieces))

    ! !
    ! ! allocate S, I, and R
    ! !
    ! allocate(lev%S(ndofs,nnodes-1))
    ! allocate(lev%I(ndofs,nnodes-1))
    ! allocate(lev%R(ndofs,nnodes-1))

    !
    ! allocate pQ and pF
    !
    if (lev%level < pf%nlevels) then

       if (lev%Finterp) then
          ! store F and Q(1) only
          ! Changed by MM Dec. 20, 2013 to allocate all pQ as well
          !
          allocate(lev%pF(ndofs,nnodes,npieces))
          allocate(lev%pQ(ndofs,nnodes))
       else
          ! store Q
          allocate(lev%pQ(ndofs,nnodes))
       end if

    end if

    !
    ! allocate Qend
    !
    allocate(lev%qend(ndofs))

  end subroutine pf_level_setup


  !
  ! Deallocate PFASST object
  !
  subroutine pf_pfasst_destroy(pf)
    type(pf_pfasst), intent(inout) :: pf

    integer :: l

    do l = 1, pf%nlevels
       call pf_level_destroy(pf%levels(l))
    end do

    deallocate(pf%levels)
    deallocate(pf%hooks)
    deallocate(pf%nhooks)
  end subroutine pf_pfasst_destroy


  !
  ! Deallocate PFASST level
  !
  subroutine pf_level_destroy(lev)
    type(pf_level), intent(inout) :: lev

    integer :: m, p

    ! flat buffers
    deallocate(lev%q0)
    deallocate(lev%send)
    deallocate(lev%recv)

    ! nodes, flags, and integration matrices
    deallocate(lev%nodes)
    deallocate(lev%nflags)
    deallocate(lev%qmat)
    deallocate(lev%smat)

    ! Q and F
    deallocate(lev%Q)
    deallocate(lev%F)

    ! S, I, and R
    ! deallocate(lev%S)
    ! deallocate(lev%I)
    ! deallocate(lev%R)

    ! pQ and pF
    if (allocated(lev%pQ)) then
       if (lev%Finterp) then
          deallocate(lev%pF)
          deallocate(lev%pQ)
       else
          deallocate(lev%pQ)
       end if
    end if

    ! qend
    deallocate(lev%qend)

    ! tau
    if (allocated(lev%tau)) then
       deallocate(lev%tau)
    end if

    ! other
    if (allocated(lev%tmat)) then
       deallocate(lev%tmat)
    end if

    if (allocated(lev%rmat)) then
       deallocate(lev%rmat)
    end if

    ! kill the sweeper
    ! call lev%sweeper%destroy(lev%sweeper)
    ! deallocate(lev%sweeper)

  end subroutine pf_level_destroy

end module pf_mod_pfasst
