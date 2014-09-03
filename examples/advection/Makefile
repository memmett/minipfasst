#
# Makefile for advection example.
#

MINIPFASST ?= ../..

EXE = main.exe
VPATHS = src
include $(MINIPFASST)/Makefile.defaults

FSRC   += src/main.f90 src/feval.f90 src/transfer.f90 src/hooks.f90 src/multigrid.f90

FFLAGS  += -Wall -Wno-unused-dummy-argument -g
LDFLAGS += -g

all: $(EXE)

include $(MINIPFASST)/Makefile.rules

#
# dependencies
#

build/user.o: build/pf_dtype.o
build/feval.o: build/user.o build/multigrid.o
build/hooks.o: build/feval.o
build/sweeper.o: build/feval.o
build/transfer.o: build/feval.o
build/main.o:  build/feval.o build/transfer.o build/hooks.o build/pfasst.o
