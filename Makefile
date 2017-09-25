# -----------------------------------------------------------------
# Makefile template
# used by configure to create Makefile
# -----------------------------------------------------------------

SHELL	= /bin/sh

FC	= gfortran
#FFLAGS	= -Ofast
FFLAGS  = -mtune=corei7 -march=native -flto -m64 -lpthread -O2

OBJBESSEL = $(BESSELOBJDIR)/*
OBJBBN	= dvode_f90_m.o dqagie.o xgetua.o cunk1.o gamln.o crati.o cuoik.o cwrsk.o cunk2.o cunik.o cunhj.o cmlri.o cairy.o cacai.o cuni2.o cuni1.o cacon.o cuchk.o cs1s2.o ckscl.o cshch.o cseri.o cbunk.o cbinu.o cbuni.o casyi.o cbknu.o xersav.o xerctl.o xerprt.o s88fmt.o xerabt.o fdump.o j4save.o xerrwv.o xerror.o r1mach.o i1mach.o d1mach.o dqelg.o dqpsrt.o dqk15i.o dqagi.o cbesk.o root_rc.o parthenope_110.o main_110.o 

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

default: BBN2k17 

BBN2k17: $(OBJBBN)
	$(FC) $(FFLAGS) -o $@ $(OBJBBN)

all: clean BBN2k17

clean:
	-rm -f *.o *.mod
