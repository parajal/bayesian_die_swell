TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90

OBJS = bcf_elements.o bcf_elements_DG_2D.o bcf_globals.o stochastic_models_2D.o stochastic_models_defs.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

bcf_elements.o: bcf_elements_DG_2D.o

bcf_elements_DG_2D.o: bcf_globals.o stochastic_models_2D.o

bcf_globals.o: stochastic_models_defs.o

stochastic_models_2D.o: stochastic_models_defs.o

clean:
	$(RM) *.o tags
