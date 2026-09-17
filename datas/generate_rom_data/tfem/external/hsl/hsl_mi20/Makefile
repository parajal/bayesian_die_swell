TFEMPATH = ../../..
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

OBJS = hsl_ma48_double.o hsl_mc65_double.o hsl_mi20d.o hsl_zb01_double.o hsl_zb01_integer.o hsl_zd11_double.o fa14.o mc71d.o mi21d.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(SRC)

hsl_ma48_double.o: hsl_zb01_double.o hsl_zb01_integer.o hsl_zd11_double.o

hsl_mc65_double.o: hsl_zd11_double.o

hsl_mi20d.o: hsl_mc65_double.o hsl_ma48_double.o hsl_zd11_double.o

clean:
	$(RM) *.o tags $(LIBCREATE)

dirclean:
	$(RM) *.o tags
