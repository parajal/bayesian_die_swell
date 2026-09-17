TFEMPATH = ../../..
include $(TFEMPATH)/Mdefs.mk
include Mdefs_io_utils.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90

OBJS = gbt_utils.o gmsh_utils.o ideas_utils.o io_utils.o printinfo.o printtofile.o sep_utils.o streamio.o tec_utils.o tfem_utils.o vtk_utils.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

gmsh_utils.o: streamio.o

io_utils.o: gbt_utils.o gmsh_utils.o ideas_utils.o printinfo.o printtofile.o sep_utils.o tec_utils.o tfem_utils.o vtk_utils.o

vtk_utils.o: streamio.o

clean:
	$(RM) *.o tags
