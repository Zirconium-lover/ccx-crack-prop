all: ccx_2.22

# Compilers
CC = gcc
FC = gfortran

# MKL (берётся из окружения; при необходимости можно задать здесь вручную)
MKLROOT ?= /mnt/pool/3/opt/intel/oneapi/mkl/2024.2
MKLINC   = -I$(MKLROOT)/include
MKLLIBD  = $(MKLROOT)/lib/intel64
MKLLIBS  = -L$(MKLLIBD) -Wl,-rpath,$(MKLLIBD) \
           -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread \
           -lgomp -lpthread -lm -ldl

# Flags
CFLAGS = -Wall -O3 -march=native -mtune=native \
         -DARCH="Linux" -DSPOOLES -DPARDISO -DARPACK -DMATRIXSTORAGE -DNETWORKOUT \
         -I../spooles.2.2 $(MKLINC) -fopenmp

FFLAGS = -Wall -O3 -march=native -mtune=native -fopenmp

.c.o:
	$(CC) $(CFLAGS) -c $<
.f.o:
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.22.c
OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)

ccx_2.22.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

# ccx_2.22.o соберётся по правилу .c.o автоматически
ccx_2.22: ccx_2.22.o ccx_2.22.a
	$(FC) $(FFLAGS) -o $@ ccx_2.22.o ccx_2.22.a \
	../spooles.2.2/spooles.a \
	../ARPACK/build/libarpack.a \
	$(MKLLIBS)

clean:
	rm -f *.o *.mod *.a ccx_2.22
