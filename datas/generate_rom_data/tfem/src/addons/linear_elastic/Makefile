TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90

OBJS = linear_elastic_elements.o linear_elastic_elements_2D_curve.o linear_elastic_elements_3D_surface.o linear_elastic_elements_generic.o linear_elastic_elements_integrate.o linear_elastic_elements_stabilized.o linear_elastic_globals.o linear_elastic_material.o linear_elastic_set_globals.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

linear_elastic_elements.o: linear_elastic_elements_2D_curve.o linear_elastic_elements_3D_surface.o linear_elastic_elements_generic.o linear_elastic_elements_integrate.o linear_elastic_elements_stabilized.o

linear_elastic_elements_2D_curve.o: linear_elastic_elements_generic.o linear_elastic_globals.o linear_elastic_material.o linear_elastic_set_globals.o

linear_elastic_elements_3D_surface.o: linear_elastic_elements_generic.o linear_elastic_globals.o linear_elastic_material.o linear_elastic_set_globals.o

linear_elastic_elements_generic.o: linear_elastic_globals.o linear_elastic_material.o

linear_elastic_elements_integrate.o: linear_elastic_elements_generic.o linear_elastic_globals.o linear_elastic_material.o linear_elastic_set_globals.o

linear_elastic_elements_stabilized.o: linear_elastic_globals.o linear_elastic_material.o linear_elastic_set_globals.o

linear_elastic_globals.o: linear_elastic_material.o

linear_elastic_set_globals.o: linear_elastic_elements_generic.o

clean:
	$(RM) *.o tags
