!
! Copyright (C) 2014 Matthew Emmett and Michael Minion.
!

module pf_mod_pfasst
  use pf_mod_dtype
  use sweeper
  use user
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

    if (present(nlevels)) pf%nlevels = nlevels

    read_cmd = .true.
    if (present(nocmd)) then
         if (nocmd) read_cmd = .false.
    end if
    if (present(fname))  then
       call pf_read_opts(pf, read_cmd, fname)
    else
       call pf_read_opts(pf, read_cmd)
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

    integer :: l
    type(pf_level), pointer :: fine, crse

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
       fine => pf%levels(l)
       crse => pf%levels(l-1)

       allocate(fine%tmat(fine%nnodes,crse%nnodes))
       allocate(fine%rmat(crse%nnodes,fine%nnodes))
       call pf_time_interpolation_matrix(fine%nodes, fine%nnodes, crse%nodes, crse%nnodes, fine%tmat)
       call pf_time_interpolation_matrix(crse%nodes, crse%nnodes, fine%nodes, fine%nnodes, fine%rmat)
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

    integer :: ndofs, nnodes

    if (lev%ndofs <= 0)   stop "ERROR: Invalid nvars/dofs (pf_pfasst.f90)."
    if (lev%nnodes <= 0)  stop "ERROR: Invalid nnodes (pf_pfasst.f90)."
    if (lev%nsweeps <= 0) stop "ERROR: Invalid nsweeps (pf_pfasst.f90)."

    ndofs  = lev%ndofs
    nnodes = lev%nnodes

    lev%residual = -1

    if (lev%level < pf%nlevels) then
       allocate(lev%tau(ndofs,nnodes-1))
    end if

    allocate(lev%q0(ndofs))
    allocate(lev%send(ndofs))
    allocate(lev%recv(ndofs))
    allocate(lev%nodes(nnodes))
    allocate(lev%nflags(nnodes))
    allocate(lev%smat(nnodes-1,nnodes))
    allocate(lev%qmat(nnodes-1,nnodes))
    allocate(lev%qend(ndofs))
    allocate(lev%Q(ndofs,nnodes))
    allocate(lev%F(ndofs,nnodes,npieces))

    if (lev%level < pf%nlevels) then
       if (lev%Finterp) then
          allocate(lev%pF(ndofs,nnodes,npieces))
          allocate(lev%pQ(ndofs,nnodes))
       else
          allocate(lev%pQ(ndofs,nnodes))
       end if
    end if

    if (btest(pf%qtype, 8)) then
       call pf_quadrature(pf%qtype, nnodes, pf%levels(1)%nnodes, &
            lev%nodes, lev%nflags, lev%smat, lev%qmat)
    else
       call pf_quadrature(pf%qtype, nnodes, pf%levels(pf%nlevels)%nnodes, &
            lev%nodes, lev%nflags, lev%smat, lev%qmat)
    end if

    call sweeper_setup(lev)
    call user_setup(lev)
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

    deallocate(lev%q0)
    deallocate(lev%send)
    deallocate(lev%recv)
    deallocate(lev%nodes)
    deallocate(lev%nflags)
    deallocate(lev%qmat)
    deallocate(lev%smat)
    deallocate(lev%Q)
    deallocate(lev%F)
    if (allocated(lev%pQ)) then
       if (lev%Finterp) then
          deallocate(lev%pF)
          deallocate(lev%pQ)
       else
          deallocate(lev%pQ)
       end if
    end if

    deallocate(lev%qend)
    if (allocated(lev%tau)) then
       deallocate(lev%tau)
    end if
    if (allocated(lev%tmat)) then
       deallocate(lev%tmat)
    end if
    if (allocated(lev%rmat)) then
       deallocate(lev%rmat)
    end if

    call sweeper_destroy(lev)
    call user_destroy(lev)
  end subroutine pf_level_destroy

end module pf_mod_pfasst
