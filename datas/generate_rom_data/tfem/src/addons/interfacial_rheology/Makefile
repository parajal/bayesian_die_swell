TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<


.F90.o:
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90 *.F90

OBJS = interfacial_rheology_elements_2D_curve.o interfacial_rheology_elements_generic.o interfacial_rheology_globals.o interfacial_rheology_elements_3D_surface.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

interfacial_rheology_elements_2D_curve.o: interfacial_rheology_elements_generic.o interfacial_rheology_globals.o

interfacial_rheology_elements_generic.o: interfacial_rheology_globals.o

interfacial_rheology_elements_3D_surface.o: interfacial_rheology_elements_generic.o interfacial_rheology_globals.o

clean:
	$(RM) *.o tags
