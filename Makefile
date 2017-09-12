# -----------------------------------------------------------------
# Makefile template
# used by configure to create Makefile
# -----------------------------------------------------------------

SHELL	= /bin/sh

FC	= xlf
FFLAGS	= -O4

OBJIBBN	= ibbn2k2.o nuclrate.o nudcupl.o daux.o

incpara = param.in


default: Ibbn2k 

Ibbn2k: $(OBJIBBN)
	$(FC) $(FFLAGS) -o $@ $(OBJIBBN)

$(OBJIBBN): $(incpara)


clean:
	-rm -f *.o
