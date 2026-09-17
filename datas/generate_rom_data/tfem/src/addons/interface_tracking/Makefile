TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90

OBJS = curvature_gd_ma57.o interface_grid_deformation.o interface_tracking_elements.o interface_tracking_globals.o poisson_gd_ma57.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

curvature_gd_ma57.o: interface_grid_deformation.o

interface_grid_deformation.o: interface_tracking_elements.o interface_tracking_globals.o

interface_tracking_elements.o: interface_tracking_globals.o

poisson_gd_ma57.o: interface_grid_deformation.o

clean:
	$(RM) *.o tags
