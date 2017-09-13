# -----------------------------------------------------------------
# Makefile template
# used by configure to create Makefile
# -----------------------------------------------------------------

SHELL	= /bin/sh

FC	= gfortran
#FFLAGS	= -Ofast
FFLAGS  = -mtune=corei7 -march=native-flto -m64 -lpthread -O2

OBJIBBN	= bbn2k.o dvode.o dqagi.o cbesk.o

incpara = card_110


default: bbn2k 

Ibbn2k: $(OBJIBBN)
	$(FC) $(FFLAGS) -o $@ $(OBJIBBN)

$(OBJIBBN): $(incpara)


clean:
	-rm -f *.o
