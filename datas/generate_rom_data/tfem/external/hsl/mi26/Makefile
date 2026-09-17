TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk
include $(TFEMPATH)/external/Mdefs_external.mk

.SUFFIXES:

.SUFFIXES: .o .f

.f.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f

OBJS = fd15d.o mi26d.o

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
