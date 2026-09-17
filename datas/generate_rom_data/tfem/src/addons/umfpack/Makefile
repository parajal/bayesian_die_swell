TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<


.F90.o:
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90 *.F90

OBJS = umfpack_interfaces.o umfpack.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

umfpack.o: umfpack_interfaces.o

clean:
	$(RM) *.o tags
