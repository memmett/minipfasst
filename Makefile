#
# Makefile for MINIPFASST
#

F90    ?= mpif90
FFLAGS  = -Jbuild -Ibuild
FFLAGS += -Wall -fbacktrace -g
#FFLAGS += -O2

LDFLAGS = -L/home/memmett/opt/lib
LIBS    = -lfftw3

SRC = $(wildcard src/*.f90)
OBJ = $(addprefix build/, $(notdir $(SRC:.f90=.o)))

pfasst-advdiff: $(OBJ)
	$(F90) $(OBJ) -o pfasst-advdiff $(LDFLAGS) $(LIBS)

clean:
	rm -rf build pfasst-advdiff

build/%.o: src/%.f90
	@mkdir -p build
	$(F90) $(FFLAGS) -c $< $(OUTPUT_OPTION)


# dependencies
build/main.o: $(addprefix build/, feval.o pf_mpi.o pf_pfasst.o pf_parallel.o)
build/pf_imex.o: $(addprefix build/, pf_dtype.o feval.o)
build/pf_interpolate.o: $(addprefix build/, pf_dtype.o pf_restrict.o pf_imex.o transfer.o)
build/pf_mpi.o: $(addprefix build/, pf_dtype.o)
build/pf_parallel.o: $(addprefix build/, pf_mpi.o pf_interpolate.o pf_restrict.o pf_imex.o pf_utils.o transfer.o)
build/pf_pfasst.o: $(addprefix build/, pf_dtype.o pf_quadrature.o pf_imex.o pf_utils.o)
build/pf_quadrature.o: $(addprefix build/, pf_dtype.o)
build/pf_restrict.o: $(addprefix build/, pf_dtype.o pf_utils.o pf_imex.o transfer.o)
build/pf_utils.o: $(addprefix build/, pf_dtype.o pf_imex.o)
build/transfer.o: $(addprefix build/, feval.o)
