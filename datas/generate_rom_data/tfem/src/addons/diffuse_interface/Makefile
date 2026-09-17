TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90

OBJS = diffuse_interface_elements.o diffuse_interface_elements_fully_implicit.o diffuse_interface_elements_generic.o diffuse_interface_globals.o diffuse_interface_tp_coupled_elements_generic.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

diffuse_interface_elements.o: diffuse_interface_elements_fully_implicit.o diffuse_interface_elements_generic.o diffuse_interface_tp_coupled_elements_generic.o

diffuse_interface_elements_fully_implicit.o: diffuse_interface_globals.o

diffuse_interface_elements_generic.o: diffuse_interface_globals.o

diffuse_interface_tp_coupled_elements_generic.o: diffuse_interface_globals.o

clean:
	$(RM) *.o tags
