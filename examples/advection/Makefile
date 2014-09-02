#
# Makefile for advection example.
#

MINIPFASST ?= ../..
FFTW3      ?= ../../fftw3

EXE = main.exe
VPATHS = src
include $(MINIPFASST)/Makefile.defaults

FSRC   += src/main.f90 src/feval.f90 src/transfer.f90 src/hooks.f90

FFLAGS  += -g -I$(FFTW3)/include
LDFLAGS += -g -L$(FFTW3)/lib -lfftw3

all: $(EXE)

include $(MINIPFASST)/Makefile.rules

#
# dependencies
#

build/transfer.o: build/feval.o
build/main.o:  build/feval.o build/transfer.o build/hooks.o
