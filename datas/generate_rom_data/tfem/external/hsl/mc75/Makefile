TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk
include $(TFEMPATH)/external/Mdefs_external.mk

.SUFFIXES:

.SUFFIXES: .o .f

.f.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f

OBJS = ma48d.o ma50d.o mc13d.o mc21d.o mc29d.o mc59d.o mc71d.o mc75d.o

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
