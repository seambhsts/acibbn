# -----------------------------------------------------------------
# Makefile template
# used by configure to create Makefile
# -----------------------------------------------------------------

SHELL	= /bin/sh

FC	= gfortran
#FFLAGS	= -Ofast
FFLAGS  = -mtune=corei7 -march=native-flto -m64 -lpthread -O2

OBJIBBN	= main_110.o parthenope_110.o dvode.o dqagi.o cbesk.o zeroin.o

incpara = card_110


default: bbn2k17 

Ibbn2k: $(OBJIBBN)
	$(FC) $(FFLAGS) -o $@ $(OBJIBBN)

$(OBJIBBN): $(incpara)


clean:
	-rm -f *.o
