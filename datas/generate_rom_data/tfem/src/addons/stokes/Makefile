TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<


.F90.o:
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90 *.F90

OBJS = devss_elements.o devss_elements_2D.o devss_elements_3D.o devss_elements_generic.o devss_set_globals.o mixed_stokes_elements_2D.o stokes_elements.o stokes_elements_2D_curve.o stokes_elements_3D_surface.o stokes_elements_embedded_boundary.o stokes_elements_embedded_interface.o stokes_elements_embedded_particle.o stokes_elements_integrate.o stokes_elements_preconditioning.o stokes_elements_stabilized.o stokes_globals.o stokes_set_globals.o streamfunction_elements.o stokes_elements_generic.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

devss_elements.o: devss_elements_2D.o devss_elements_3D.o devss_elements_generic.o

devss_elements_2D.o: devss_set_globals.o stokes_globals.o stokes_set_globals.o

devss_elements_3D.o: devss_set_globals.o stokes_globals.o stokes_set_globals.o

devss_elements_generic.o: stokes_globals.o stokes_set_globals.o

devss_set_globals.o: devss_elements_generic.o

mixed_stokes_elements_2D.o: stokes_globals.o stokes_set_globals.o

stokes_elements.o: stokes_elements_2D_curve.o stokes_elements_3D_surface.o stokes_elements_embedded_boundary.o stokes_elements_embedded_interface.o stokes_elements_embedded_particle.o stokes_elements_generic.o stokes_elements_integrate.o stokes_elements_preconditioning.o stokes_elements_stabilized.o

stokes_elements_2D_curve.o: stokes_elements_generic.o stokes_globals.o stokes_set_globals.o

stokes_elements_3D_surface.o: stokes_elements_generic.o stokes_globals.o stokes_set_globals.o

stokes_elements_embedded_boundary.o: stokes_elements_generic.o stokes_globals.o stokes_set_globals.o

stokes_elements_embedded_interface.o: stokes_elements_embedded_boundary.o stokes_elements_generic.o stokes_globals.o stokes_set_globals.o

stokes_elements_embedded_particle.o: stokes_elements_embedded_boundary.o stokes_globals.o stokes_set_globals.o

stokes_elements_integrate.o: stokes_elements_generic.o stokes_globals.o stokes_set_globals.o

stokes_elements_preconditioning.o: stokes_globals.o stokes_set_globals.o

stokes_elements_stabilized.o: stokes_globals.o stokes_set_globals.o

stokes_set_globals.o: stokes_elements_generic.o

streamfunction_elements.o: stokes_globals.o stokes_set_globals.o

stokes_elements_generic.o: stokes_globals.o

clean:
	$(RM) *.o tags
