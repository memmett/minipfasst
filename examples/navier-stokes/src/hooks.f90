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
  end interface

  interface
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
    character(*),   intent(in   ) :: fname
    real(pfdp),     intent(in   ) :: q(:)
    type(pf_level), intent(inout) :: lev
    call dump_velocity_c(trim(output) // c_null_char, trim(fname) // c_null_char, &
         lev%user%nx, lev%user%ny, lev%user%nz, q)
  end subroutine dump_velocity

  subroutine dump_vorticity(fname, q, lev)
    use probin, only: output
    use initial
    character(*),   intent(in   ) :: fname
    type(pf_level), intent(inout) :: lev
    real(pfdp),     intent(in   ) :: q(:)

    real(pfdp) :: vort(lev%user%nx,lev%user%ny,lev%user%nz)
    call vorticity(q, lev, vort)
    call dump_vorticity_c(trim(output) // c_null_char, trim(fname) // c_null_char, &
         lev%user%nx, lev%user%ny, lev%user%nz, vort)
  end subroutine dump_vorticity

  subroutine echo_error(pf, level)
    use probin, only: npts, nu
    use initial, only: shapiro
    type(pf_pfasst), intent(inout) :: pf
    type(pf_level),  intent(inout) :: level

    real(pfdp) :: qex(size(level%qend))

    call shapiro(qex, pf%t0+pf%dt, npts(level%level), nu)
    print *, '  error   ', maxval(abs(level%qend-qex))
  end subroutine echo_error

end module hooks