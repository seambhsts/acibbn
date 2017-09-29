# -----------------------------------------------------------------
# Makefile template
# used by configure to create Makefile
# -----------------------------------------------------------------

#SHELL	= /bin/sh

FC	= gfortran
#FFLAGS	= -Ofast
FFLAGS  = -mtune=corei7 -march=native -flto -m64 -lpthread -O2

OBJQUAD = dqagie.o xgetua.o xersav.o xerctl.o xerprt.o s88fmt.o xerabt.o fdump.o j4save.o xerrwv.o xerror.o r1mach.o i1mach.o d1mach.o dqelg.o dqpsrt.o dqk15i.o dqagi.o

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

libquad.a: $(OBJQUAD)
	ar rv libquad.a $(OBJQUAD)
	ranlib libquad.a
	

default: libquad.a 
all: clean libquad.a

clean:
	-rm -f *.o *.mod
