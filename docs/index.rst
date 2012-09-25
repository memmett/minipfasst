MINIPFASST
==========

MINIPFASST is a lightweight implementation of the Parallel Full
Approximation Scheme in Space and Time (PFASST) algorithm.


Software overview
^^^^^^^^^^^^^^^^^


PFASST types
------------

MINIPFASST is controlled by instantiating and manipulating two main
derived types:

#. ``pf_pfasst_t`` - the main PFASST type.  This contains basic
                     parameters, flags, and a list of the PFASST
                     levels (which are instances of ``pf_level_t``).

#. ``pf_comm_t`` - the 'time' communicator.  This contains information
                   about the time (MPI or pthreads) communicator.

These derived types are defined in ``src/pf_dtype.f90``.


Time evolution
--------------

MINIFASST evolves systems of ODEs (perhaps obtained from a MOL
discretisation of a PDE) using an IMEX SDC time-stepper that uses
forward/backward Euler sub-stepping.

To use the IMEX time-stepper, you must implement the explicit
``eval_f1`` subroutine, and the implicit ``eval_f2`` and ``comp_f2``
subroutines.  These must be in your ``feval`` module.  The ``comp_f2``
routine should perform a backward-Euler solve: :math:`y - dt f_2(y) =
rhs`.


Transfer functions
------------------

To transfer solutions between the various time/space discretisations,
you must implement the ``interpolate`` and ``restrict`` subroutines in
your ``transfer`` module.  These routines interpolate or restrict
solutions between spatial discretisations (which you define).
F90PFASST takes care of temporal interpolation/restriction for you
(using polynomials for interpolation and injection for restriction).

The user defined function evaluation and transfer routines are passed
integers (``levelF`` and/or ``levelG``) that specify which PFASST
level they are operating on (with level 1 being the finest level).
