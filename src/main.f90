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

program fpfasst
  use pf_mod_dtype
  use pf_mod_comm
  use pf_mod_pfasst
  use pf_mod_parallel
  use pf_mod_mpi, only: MPI_COMM_WORLD
  use feval

  implicit none

  type(pf_pfasst_t) :: pf
  type(pf_comm_t)   :: tcomm
  integer           :: ierror, nlevs, nvars(3), nnodes(3)
  double precision  :: dt
  double precision, pointer :: q0(:)

 
  !!!! initialize mpi
  call mpi_init(ierror)
  if (ierror .ne. 0) &
       stop "ERROR: Can't initialize MPI."


  !!!! initialize pfasst
  nvars  = [ 64, 32, 16 ]
  nnodes = [ 9, 5, 3 ]
  dt     = 0.01d0
  nlevs  = 3

  call create(tcomm, MPI_COMM_WORLD)
  call create(pf, tcomm, nlevs, nvars(1:nlevs), nnodes(1:nlevs))

  pf%niters  = 12
  pf%qtype   = 1

  if (nlevs > 1) then
     pf%levels(nlevs)%nsweeps = 2
  end if

  call setup(tcomm, pf)
  call setup(pf)

  !!!! initialize advection/diffusion
  call feval_init(size(nvars), nvars)
  allocate(q0(nvars(1)))
  call initial(q0)

  !!!! run
  call pfasst_run(pf, q0, dt, 0.0d0, tcomm%nproc)


  !!!! done
  deallocate(q0)
  call destroy(pf)
  call destroy(tcomm)
  call feval_finalize()
  call mpi_finalize(ierror)

end program fpfasst
