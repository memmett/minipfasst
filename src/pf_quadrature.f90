!
! Copyright (C) 2014 Matthew Emmett and Michael Minion.
!

module pf_mod_quadrature
  use pf_mod_dtype
  use sdc_mod_quadrature
  implicit none
contains

  subroutine pf_quadrature(qtype_in, nnodes, nnodes0, nodes, nflags, smat, qmat)
    integer,    intent(in)  :: qtype_in, nnodes, nnodes0
    real(pfdp), intent(out) :: nodes(nnodes)
    real(pfdp), intent(out) :: smat(nnodes-1,nnodes), qmat(nnodes-1,nnodes)
    integer,    intent(out) :: nflags(nnodes)

    real(pfqp) :: qnodes0(nnodes0), qnodes(nnodes), dt
    real(pfdp) :: qmat0(nnodes0-1,nnodes0), smat0(nnodes0-1,nnodes0)
    integer    :: flags0(nnodes0)

    integer :: qtype, i, r, refine
    logical :: composite, proper, no_left

    proper    = btest(qtype_in, 8)
    composite = btest(qtype_in, 9)
    no_left   = btest(qtype_in, 10)

    qmat = 0
    smat = 0
    flags0 = 0
    nflags = 0

    qtype = qtype_in
    if (proper)    qtype = qtype - SDC_PROPER_NODES
    if (composite) qtype = qtype - SDC_COMPOSITE_NODES
    if (no_left)   qtype = qtype - SDC_NO_LEFT


    if (composite) then

       ! nodes are given by repeating the coarsest set of nodes.  note
       ! that in this case nnodes0 corresponds to the coarsest number
       ! of nodes.

       refine = (nnodes - 1) / (nnodes0 - 1)

       call sdc_qnodes(qnodes0, flags0, qtype-SDC_COMPOSITE_NODES, nnodes0)
       call sdc_qmats(qmat0, smat0, qnodes0, qnodes0, flags0, nnodes0, nnodes0)

       dt = 1.0_pfqp / refine
       do i = 1, refine
          r = (i-1)*(nnodes0-1)+1
          qnodes(r:r+nnodes0) = dt * ((i-1) + qnodes0)
          smat(r:r+nnodes0,r:r+nnodes0-1) = smat0 / refine
       end do

       nodes = real(qnodes, pfdp)

    else if (proper) then

       ! nodes are given by proper quadrature rules

       call sdc_qnodes(qnodes, nflags, qtype-SDC_PROPER_NODES, nnodes)
       nodes = real(qnodes, pfdp)

       call sdc_qmats(qmat, smat, qnodes, qnodes, nflags, nnodes, nnodes)

    else

       ! nodes are given by refining the finest set of nodes.  note
       ! that in this case nnodes0 corresponds to the finest number of
       ! nodes.

       refine = (nnodes0 - 1) / (nnodes - 1)

       call sdc_qnodes(qnodes0, flags0, qtype, nnodes0)

       qnodes = qnodes0(::refine)
       nodes  = real(qnodes, pfdp)
       nflags = flags0(::refine)

       if (no_left) nflags(1) = 0

       call sdc_qmats(qmat, smat, qnodes, qnodes, nflags, nnodes, nnodes)

    end if


    if (all(nodes == 0.0d0)) then
       stop 'ERROR: pf_quadrature: invalid SDC nnodes.'
    end if

  end subroutine pf_quadrature

end module pf_mod_quadrature
