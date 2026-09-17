TFEMPATH = ../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODINC) $<

SRC = *.f90

OBJS =

MAINP = settling2_ma41 stokes11_ma41 stokes11_ma41_metis stokes11_ma57 stokes11_pardiso_sym stokes11_pardiso_unsym stokes11_sk_metis

all: objects programs tags

objects: $(OBJS)

programs: $(MAINP)

tags: $(SRC)
	ctags --fortran-kinds=+i $(SRC) $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

settling2_ma41: $(OBJS) settling2_ma41.o
	$(FC) $(FFLAGS) -o $@ settling2_ma41.o $(OBJS) $(LIBS)

stokes11_ma41: $(OBJS) stokes11_ma41.o
	$(FC) $(FFLAGS) -o $@ stokes11_ma41.o $(OBJS) $(LIBS)

stokes11_ma41_metis: $(OBJS) stokes11_ma41_metis.o
	$(FC) $(FFLAGS) -o $@ stokes11_ma41_metis.o $(OBJS) $(LIBS)

stokes11_ma57: $(OBJS) stokes11_ma57.o
	$(FC) $(FFLAGS) -o $@ stokes11_ma57.o $(OBJS) $(LIBS)

stokes11_pardiso_sym: $(OBJS) stokes11_pardiso_sym.o
	$(FC) $(FFLAGS) -o $@ stokes11_pardiso_sym.o $(OBJS) $(LIBS)

stokes11_pardiso_unsym: $(OBJS) stokes11_pardiso_unsym.o
	$(FC) $(FFLAGS) -o $@ stokes11_pardiso_unsym.o $(OBJS) $(LIBS)

stokes11_sk_metis: $(OBJS) stokes11_sk_metis.o
	$(FC) $(FFLAGS) -o $@ stokes11_sk_metis.o $(OBJS) $(LIBS)

clean:
	$(RM)  settling2_ma41 stokes11_ma41 stokes11_ma41_metis stokes11_ma57 stokes11_pardiso_sym stokes11_pardiso_unsym stokes11_sk_metis *.o *.mod *.bin tags
