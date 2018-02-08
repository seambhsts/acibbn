# -----------------------------------------------------------------
# Makefile template
# used by configure to create Makefile
# -----------------------------------------------------------------

#SHELL	= /bin/sh

FC	= gfortran
LIB=/usr/local/Cellar/openblas/0.2.20_1/lib/libopenblas.a

#FFLAGS	= -Ofast
#FFLAGS  = -g -Wall -fbacktrace -ffpe-trap=zero,overflow,underflow -fno-automatic -fno-second-underscore -fno-range-check -lpthread #-freal-8-real-1
#FFLAGS  = -g -Wall -fno-automatic -fno-second-underscore -fno-range-check -lpthread
FFLAGS  = -Wall -framework Accelerate -mtune=corei7 -march=native -flto -fno-automatic -fno-second-underscore -fno-range-check -m64 -lpthread -O2

OBJBBN	= machine.o dbesk_slatec.o quadpack_double.o rowmap.o root_rc.o parthenope_110.o main_110.o 

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@
%.o : %.f
	$(FC) -c $(FFLAGS) $< -o $@

BBN2k17:  $(OBJBBN)
	$(FC) $(FFLAGS) -o $@ $(OBJBBN) $(LIB)

default: BBN2k17 
all: clean BBN2k17

libbesl.a: 
	make -C ./amos
	
libquad.a:
	make -C ./quad
	

clean:
	-rm -f *.o *.mod
	
