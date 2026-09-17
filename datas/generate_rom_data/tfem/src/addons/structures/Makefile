TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90

OBJS = structures_elements.o structures_elements_beam1.o structures_elements_beam2.o structures_elements_beam3.o structures_elements_truss.o structures_globals.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

structures_elements.o: structures_elements_beam1.o structures_elements_beam2.o structures_elements_beam3.o structures_elements_truss.o

structures_elements_beam1.o: structures_globals.o

structures_elements_beam2.o: structures_globals.o

structures_elements_beam3.o: structures_globals.o

structures_elements_truss.o: structures_globals.o

clean:
	$(RM) *.o tags
