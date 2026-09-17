TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90

OBJS = convection_diffusion_elements.o convection_diffusion_elements_embedded_boundary.o convection_diffusion_elements_embedded_interface.o convection_diffusion_elements_generic.o convection_diffusion_globals.o convection_diffusion_set_globals.o laplace_elements_ale.o poisson_elements.o poisson_globals.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

convection_diffusion_elements.o: convection_diffusion_elements_embedded_boundary.o convection_diffusion_elements_embedded_interface.o convection_diffusion_elements_generic.o poisson_elements.o

convection_diffusion_elements_embedded_boundary.o: convection_diffusion_globals.o convection_diffusion_set_globals.o

convection_diffusion_elements_embedded_interface.o: convection_diffusion_elements_embedded_boundary.o convection_diffusion_globals.o convection_diffusion_set_globals.o

convection_diffusion_elements_generic.o: convection_diffusion_globals.o poisson_elements.o poisson_globals.o

convection_diffusion_globals.o: poisson_globals.o

convection_diffusion_set_globals.o: convection_diffusion_elements_generic.o poisson_elements.o

laplace_elements_ale.o: convection_diffusion_set_globals.o poisson_globals.o

poisson_elements.o: poisson_globals.o

clean:
	$(RM) *.o tags
