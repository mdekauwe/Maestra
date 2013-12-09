##
## Makefile for building a Maestra executable.
## Martin De Kauwe, 20/06/2011
##
PROG =	maestra

SRCS =	getmet.f90 inout.f90 maestra.f90 maestcom.f90 metcom.f90 physiol.f90 \
        radn.f90 unstor.f90 utils.f90

OBJS =	$(SRCS:.f90=.o)

LIBS =	

INCLS = 

# choose compiler
FC=gfortran
#FC=/opt/local/bin/gfortran-mp-4.4
FFLAGS = -g -Wno-unused-dummy-argument -fbounds-check -ffixed-line-length-none -Wall -Wextra -ffpe-trap=zero,overflow,underflow 

all: $(PROG)

$(PROG):	$(OBJS)
			$(FC) $(FFLAGS) -o $@ $(OBJS) $(LIBS) $(INCLS)

clean:
			rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
			$(FC) $(FFLAGS) -c $<

getmet.o:   maestcom.o metcom.o
inout.o:    maestcom.o
maestra.o:  maestcom.o metcom.o
physiol.o:  maestcom.o metcom.o
radn.o:     maestcom.o
unstor.o:   maestcom.o
utils.o:    maestcom.o

