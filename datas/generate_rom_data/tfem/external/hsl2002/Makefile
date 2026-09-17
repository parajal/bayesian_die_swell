TFEMPATH = ../..
include $(TFEMPATH)/Mdefs.mk
include Mdefs_hsl.mk

.SUFFIXES:

.SUFFIXES: .o .f

.f.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f

OBJS = ea23ad.o fd05adh.o id05adh.o ma41ad.o ma57ad.o mc21ad.o mc29ad.o mc34ad.o mc47ad.o mc59ad.o mc71ad.o

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
