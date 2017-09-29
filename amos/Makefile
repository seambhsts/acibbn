# -----------------------------------------------------------------
# Makefile template
# used by configure to create Makefile
# -----------------------------------------------------------------

#SHELL	= /bin/sh

FC	= gfortran
#FFLAGS	= -Ofast
FFLAGS  = -mtune=corei7 -march=native -flto -m64 -lpthread -O2

OBJBESSEL = cunk1.o gamln.o crati.o cuoik.o cwrsk.o cunk2.o cunik.o cunhj.o cmlri.o cairy.o cacai.o cuni2.o cuni1.o cacon.o cuchk.o cs1s2.o ckscl.o cshch.o cseri.o cbunk.o cbinu.o cbuni.o casyi.o cbknu.o cbesk.o i1mach.o

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

libbesl.a: $(OBJBESSEL)
	ar rv libbesl.a $(OBJBESSEL)
	ranlib libbesl.a
	

default: libbesl.a 
all: clean libbesl.a

clean:
	-rm -f *.o *.mod
