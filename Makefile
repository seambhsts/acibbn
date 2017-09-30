# -----------------------------------------------------------------
# Makefile template
# used by configure to create Makefile
# -----------------------------------------------------------------

#SHELL	= /bin/sh

FC	= gfortran
#FFLAGS	= -Ofast
FFLAGS  = -mtune=corei7 -march=native -flto -m64 -lpthread -O2

OBJBBN	= cbesk_slatec.o quadpack_double.o dvode_f90_m.o root_rc.o parthenope_110.o main_110.o 

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

BBN2k17:  $(OBJBBN) 
	$(FC) $(FFLAGS) -o $@ $(OBJBBN) 

default: BBN2k17 
all: clean BBN2k17

libbesl.a: 
	make -C ./amos
	
libquad.a:
	make -C ./quad
	

clean:
	-rm -f *.o *.mod
	
