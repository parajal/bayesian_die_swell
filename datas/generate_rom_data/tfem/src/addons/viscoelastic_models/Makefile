TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<


.F90.o:
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90 *.F90

OBJS = viscoelastic_models.o viscoelastic_models_2D.o viscoelastic_models_2D_b.o viscoelastic_models_2D_log.o viscoelastic_models_adapted_lambda.o viscoelastic_models_defs.o viscoelastic_models_plastic_strain.o viscoelastic_models_3D.o viscoelastic_models_3D_b.o viscoelastic_models_3D_log.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

viscoelastic_models.o: viscoelastic_models_2D.o viscoelastic_models_2D_b.o viscoelastic_models_2D_log.o viscoelastic_models_3D.o viscoelastic_models_3D_b.o viscoelastic_models_3D_log.o

viscoelastic_models_2D.o: viscoelastic_models_adapted_lambda.o viscoelastic_models_defs.o viscoelastic_models_plastic_strain.o

viscoelastic_models_2D_b.o: viscoelastic_models_2D.o viscoelastic_models_adapted_lambda.o viscoelastic_models_defs.o viscoelastic_models_plastic_strain.o

viscoelastic_models_2D_log.o: viscoelastic_models_2D.o viscoelastic_models_adapted_lambda.o viscoelastic_models_defs.o

viscoelastic_models_adapted_lambda.o: viscoelastic_models_defs.o

viscoelastic_models_plastic_strain.o: viscoelastic_models_adapted_lambda.o viscoelastic_models_defs.o

viscoelastic_models_3D.o: viscoelastic_models_adapted_lambda.o viscoelastic_models_defs.o viscoelastic_models_plastic_strain.o

viscoelastic_models_3D_b.o: viscoelastic_models_3D.o viscoelastic_models_adapted_lambda.o viscoelastic_models_defs.o viscoelastic_models_plastic_strain.o

viscoelastic_models_3D_log.o: viscoelastic_models_3D.o viscoelastic_models_adapted_lambda.o viscoelastic_models_defs.o viscoelastic_models_plastic_strain.o

clean:
	$(RM) *.o tags
