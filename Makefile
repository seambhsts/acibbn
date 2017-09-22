# -----------------------------------------------------------------
# Makefile template
# used by configure to create Makefile
# -----------------------------------------------------------------

SHELL	= /bin/sh

FC	= gfortran
#FFLAGS	= -Ofast
FFLAGS  = -mtune=corei7 -march=native -flto -m64 -lpthread -O2
BESSELOBJDIR = ./amos

OBJBESSEL = $(BESSELOBJDIR)/*
OBJBBN	= dvode_f90_m.o dqagi.o cbesk.o root_rc.o parthenope_110.o main_110.o 
.SUFFIXES: .f90 .f

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

incpara = card_110


default: BBN2k17 

BBN2k17: $(OBJBBN)
	$(FC) $(FFLAGS) -o $@ $(OBJBBN)

$(OBJBBN): $(incpara)

all: clean BBN2k17

clean:
	-rm -f *.o *.mod
