TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90

OBJS = hsl_mc75.o hsl_solve1_mi20.o hsl_solve2_mi20.o hsl_solve3_mi20.o hsl_solve_diag.o hsl_solve_mi20.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

hsl_solve1_mi20.o: hsl_solve_diag.o

hsl_solve2_mi20.o: hsl_solve_diag.o

hsl_solve3_mi20.o: hsl_solve_diag.o

clean:
	$(RM) *.o tags
