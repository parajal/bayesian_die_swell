TFEMPATH = ../..
include $(TFEMPATH)/Mdefs.mk
include $(TFEMPATH)/external/Mdefs_external.mk

.SUFFIXES:

.SUFFIXES: .o .f90

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90

OBJS = nr.o nrtype.o nrutil.o odeint.o rkck.o rkqs.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(SRC)

nr.o: nrtype.o nrtype.o nrtype.o nrtype.o

nrutil.o: nrtype.o

odeint.o: nrtype.o nrtype.o nrtype.o nrutil.o

rkck.o: nrtype.o nrutil.o

rkqs.o: nrtype.o nr.o nrutil.o

clean:
	$(RM) *.o tags $(LIBCREATE)

dirclean:
	$(RM) *.o tags
