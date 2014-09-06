#
# Makefile for advection example.
#

MINIPFASST ?= ../..

EXE = main.exe
VPATHS = src
FFTW   = $(MINIPFASST)/fftw3
include $(MINIPFASST)/Makefile.defaults

FSRC    += src/main.f90 src/feval.f90 src/hooks.f90 src/probin.f90 src/initial.f90
CSRC    += src/dump.c
FFLAGS  += -Wall -Wno-unused-dummy-argument -g -I$(FFTW)/include
LDFLAGS += -g -L$(FFTW)/lib -lfftw3_threads -lfftw3

#FFLAGS  += -O3 -g -I$(FFTW)/include
#LDFLAGS += -L$(FFTW)/lib -lfftw3

all: $(EXE)

include $(MINIPFASST)/Makefile.rules

#
# dependencies
#

build/initial.o: build/feval.o
build/user.o: build/pf_dtype.o build/probin.o
build/feval.o: build/user.o build/probin.o
build/hooks.o: build/feval.o build/initial.o
build/sweeper.o: build/feval.o
build/transfer.o: build/feval.o
build/main.o:  build/feval.o build/transfer.o build/hooks.o build/pfasst.o build/initial.o
