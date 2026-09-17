TFEMPATH = ../..
include $(TFEMPATH)/Mdefs.mk
include $(TFEMPATH)/external/Mdefs_external.mk

.SUFFIXES:

.SUFFIXES: .o .f90

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

.SUFFIXES: .o .f

.f.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90 *.f

OBJS = distdot.o blassm.o ccn.o color.o dsepart.o formats.o ilut.o infofun.o inout.o iters.o matvec.o unary.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(SRC)

clean:
	$(RM) *.o tags $(LIBCREATE)

dirclean:
	$(RM) *.o tags
