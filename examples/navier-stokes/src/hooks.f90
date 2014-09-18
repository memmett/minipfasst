!
! Copyright (c) 2014, Matthew Emmett and Michael Minion.  All rights reserved.
!

module hooks
  use iso_c_binding
  use pf_mod_dtype
  implicit none

  interface
     subroutine dump_velocity_c(dname, fname, nx, ny, nz, array) bind(c)
       use iso_c_binding
       character(c_char), intent(in   )        :: dname, fname
       integer(c_int),    intent(in   ), value :: nx, ny, nz
       real(c_double),    intent(in   )        :: array(2*3*nx*ny*nz)
     end subroutine dump_velocity_c

     subroutine read_velocity_c(fname, nx, ny, nz, array) bind(c)
       use iso_c_binding
       character(c_char), intent(in   )        :: fname
       integer(c_int),    intent(in   ), value :: nx, ny, nz
       real(c_double),    intent(in   )        :: array(2*3*nx*ny*nz)
     end subroutine read_velocity_c

     subroutine dump_vorticity_c(dname, fname, nx, ny, nz, array) bind(c)
       use iso_c_binding
       character(c_char), intent(in   )        :: dname, fname
       integer(c_int),    intent(in   ), value :: nx, ny, nz
       real(c_double),    intent(in   )        :: array(nx*ny*nz)
     end subroutine dump_vorticity_c
  end interface

contains

  subroutine dump_velocity(fname, q, lev)
    use probin, only: output
    use feval
    character(*),   intent(in   ) :: fname
    real(pfdp),     intent(in   ) :: q(:)
    type(pf_level), intent(inout) :: lev

    complex(pfdp), dimension(lev%user%nx,lev%user%ny,lev%user%nz) :: u, v, w
    real(pfdp), dimension(lev%ndofs) :: qdump

    call unpack3(q, u, w, v)
    call dfftw_execute_dft_(lev%user%bft, u, u)
    call dfftw_execute_dft_(lev%user%bft, v, v)
    call dfftw_execute_dft_(lev%user%bft, w, w)
    call pack3(qdump, u, v, w)

    call dump_velocity_c(trim(output) // c_null_char, trim(fname) // c_null_char, &
         lev%user%nx, lev%user%ny, lev%user%nz, qdump)
  end subroutine dump_velocity

  subroutine read_velocity(fname, q, lev)
    use feval
    character(*),   intent(in   ) :: fname
    real(pfdp),     intent(  out) :: q(:)
    type(pf_level), intent(inout) :: lev

    complex(pfdp), dimension(lev%user%nx,lev%user%ny,lev%user%nz) :: u, v, w
    real(pfdp), dimension(lev%ndofs) :: qdump

    call read_velocity_c(trim(fname) // c_null_char, &
         lev%user%nx, lev%user%ny, lev%user%nz, qdump)

    qdump = qdump * lev%user%scale
    call unpack3(qdump, u, w, v)
    call dfftw_execute_dft_(lev%user%fft, u, u)
    call dfftw_execute_dft_(lev%user%fft, v, v)
    call dfftw_execute_dft_(lev%user%fft, w, w)
    call pack3(q, u, v, w)
  end subroutine read_velocity

  subroutine dump_vorticity(fname, q, lev)
    use probin, only: output
    use initial
    character(*),   intent(in   ) :: fname
    type(pf_level), intent(inout) :: lev
    real(pfdp),     intent(in   ) :: q(:)

    real(pfdp) :: vsq(lev%user%nx,lev%user%ny,lev%user%nz)
    call vorticity2(q, lev, vsq)
    call dump_vorticity_c(trim(output) // c_null_char, trim(fname) // c_null_char, &
         lev%user%nx, lev%user%ny, lev%user%nz, vsq)
  end subroutine dump_vorticity

  subroutine echo_error_hook(pf, level)
    use probin, only: npts, nu
    use initial, only: shapiro
    type(pf_pfasst), intent(inout) :: pf
    type(pf_level),  intent(inout) :: level

    real(pfdp) :: qex(size(level%qend))

    call shapiro(level%user%fft, qex, pf%t0+pf%dt, npts(level%level), nu)
    print *, '-> error   ', pf%step, pf%iter, level%level, maxval(abs(level%qend-qex))
  end subroutine echo_error_hook

  subroutine echo_enstrophy_hook(pf, level)
    use probin, only: npts, nu
    use initial, only: enstrophy
    type(pf_pfasst), intent(inout) :: pf
    type(pf_level),  intent(inout) :: level
    real(pfdp) :: ens
    if (level%level /= pf%nlevels) return
    call enstrophy(level%qend, level, ens)
    print *, '  enstrophy', pf%step, ens
  end subroutine echo_enstrophy_hook

end module hooks
