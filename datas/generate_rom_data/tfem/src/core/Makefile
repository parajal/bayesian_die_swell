TFEMPATH = ../..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODWRITE) $(FMODINC) $<

SRC = *.f90

OBJS = array_defs.o eig2D3D.o element_defs.o element_utils.o eltree.o eltree_basic.o eltree_interface.o eltree_volume.o elvector.o elvector_defs.o execute_command_line.o gauss.o gauss_defs.o gauss_macro.o gauss_mapped.o gauss_spectral.o gauss_standard.o gauss_standard_numeric.o gauss_standard_numeric2.o glob_defs.o intmesh.o kind_defs.o limits.o math_defs.o mesh.o meshgen.o meshgen_basic.o meshgen_basic_common.o meshgen_basic_hexahedron.o meshgen_basic_line.o meshgen_basic_quadrilateral.o meshgen_construct.o meshgen_extra.o meshgen_gmsh_index.o meshgen_objects.o meshgen_parts.o misc.o pos_array.o postprocessing.o problem.o problem_defs.o set_optional.o shapefunc.o shapefunc_highorder.o shapefunc_macro.o shapefunc_spectral.o shapefunc_standard.o shapefunc_standard_global.o sparse.o submesh.o system.o system_connection.o system_constraint.o system_defs.o system_dependency.o system_dg.o system_matrix.o system_vector.o tensor.o tfem.o tfem_elem.o vector.o vector_defs.o

all: $(LIBCREATE) tags

$(LIBCREATE): $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

tags: $(SRC)
	ctags --fortran-kinds=+i $(SRC)

array_defs.o: kind_defs.o

eig2D3D.o: kind_defs.o

element_defs.o: array_defs.o eltree_basic.o elvector_defs.o kind_defs.o limits.o mesh.o problem_defs.o set_optional.o system_defs.o vector_defs.o

element_utils.o: element_defs.o elvector.o kind_defs.o vector.o

eltree.o: eltree_basic.o eltree_interface.o eltree_volume.o intmesh.o submesh.o

eltree_basic.o: glob_defs.o intmesh.o kind_defs.o limits.o set_optional.o submesh.o

eltree_interface.o: eltree_basic.o glob_defs.o kind_defs.o meshgen.o

eltree_volume.o: eltree_basic.o kind_defs.o meshgen.o

elvector.o: element_defs.o elvector_defs.o kind_defs.o mesh.o problem_defs.o set_optional.o system_defs.o vector_defs.o

elvector_defs.o: array_defs.o kind_defs.o

gauss.o: gauss_defs.o gauss_macro.o gauss_mapped.o gauss_spectral.o gauss_standard.o gauss_standard_numeric.o glob_defs.o limits.o set_optional.o

gauss_defs.o: kind_defs.o

gauss_macro.o: gauss_defs.o gauss_standard.o kind_defs.o

gauss_mapped.o: gauss_defs.o gauss_standard_numeric.o kind_defs.o

gauss_spectral.o: gauss_defs.o kind_defs.o

gauss_standard.o: gauss_defs.o kind_defs.o

gauss_standard_numeric.o: gauss_defs.o gauss_standard_numeric2.o kind_defs.o limits.o

gauss_standard_numeric2.o: gauss_defs.o kind_defs.o

glob_defs.o: kind_defs.o

intmesh.o: glob_defs.o limits.o mesh.o meshgen_extra.o set_optional.o

limits.o: kind_defs.o

math_defs.o: kind_defs.o

mesh.o: array_defs.o glob_defs.o limits.o set_optional.o

meshgen.o: kind_defs.o meshgen_basic.o meshgen_construct.o meshgen_extra.o meshgen_objects.o meshgen_parts.o

meshgen_basic.o: kind_defs.o meshgen_basic_common.o meshgen_basic_hexahedron.o meshgen_basic_line.o meshgen_basic_quadrilateral.o

meshgen_basic_common.o: kind_defs.o limits.o mesh.o

meshgen_basic_hexahedron.o: glob_defs.o meshgen_basic_common.o

meshgen_basic_line.o: gauss_spectral.o kind_defs.o meshgen_basic_common.o

meshgen_basic_quadrilateral.o: gauss_spectral.o glob_defs.o kind_defs.o meshgen_basic_common.o

meshgen_construct.o: gauss.o glob_defs.o kind_defs.o limits.o mesh.o misc.o set_optional.o

meshgen_extra.o: array_defs.o glob_defs.o kind_defs.o limits.o mesh.o set_optional.o shapefunc_standard.o

meshgen_objects.o: array_defs.o glob_defs.o kind_defs.o limits.o mesh.o set_optional.o shapefunc.o

meshgen_parts.o: glob_defs.o limits.o mesh.o meshgen_construct.o meshgen_gmsh_index.o meshgen_objects.o misc.o

misc.o: kind_defs.o

pos_array.o: kind_defs.o mesh.o problem_defs.o

postprocessing.o: element_defs.o kind_defs.o mesh.o pos_array.o problem.o problem_defs.o set_optional.o system_defs.o system_matrix.o vector_defs.o

problem.o: glob_defs.o limits.o math_defs.o mesh.o misc.o pos_array.o problem_defs.o set_optional.o

problem_defs.o: array_defs.o kind_defs.o mesh.o set_optional.o sparse.o system_defs.o

set_optional.o: kind_defs.o

shapefunc.o: glob_defs.o limits.o mesh.o set_optional.o shapefunc_highorder.o shapefunc_macro.o shapefunc_spectral.o shapefunc_standard.o shapefunc_standard_global.o

shapefunc_highorder.o: kind_defs.o shapefunc_standard.o

shapefunc_macro.o: kind_defs.o mesh.o

shapefunc_spectral.o: gauss_spectral.o kind_defs.o set_optional.o shapefunc_highorder.o

shapefunc_standard.o: kind_defs.o limits.o

shapefunc_standard_global.o: kind_defs.o

sparse.o: kind_defs.o misc.o

submesh.o: glob_defs.o limits.o mesh.o meshgen_extra.o set_optional.o

system.o: system_connection.o system_constraint.o system_dependency.o system_dg.o system_matrix.o system_vector.o

system_connection.o: element_defs.o kind_defs.o mesh.o misc.o problem_defs.o set_optional.o sparse.o system_defs.o system_matrix.o

system_constraint.o: element_defs.o kind_defs.o mesh.o problem_defs.o set_optional.o sparse.o system_defs.o system_matrix.o

system_defs.o: kind_defs.o sparse.o

system_dependency.o: element_defs.o glob_defs.o mesh.o problem_defs.o set_optional.o sparse.o system_defs.o system_matrix.o

system_dg.o: element_defs.o kind_defs.o mesh.o problem_defs.o sparse.o system_defs.o system_matrix.o

system_matrix.o: element_defs.o glob_defs.o kind_defs.o math_defs.o mesh.o pos_array.o problem_defs.o shapefunc.o sparse.o system_defs.o system_vector.o

system_vector.o: element_defs.o kind_defs.o mesh.o misc.o pos_array.o problem_defs.o set_optional.o sparse.o system_defs.o

tensor.o: kind_defs.o

tfem.o: element_defs.o eltree.o elvector.o kind_defs.o mesh.o meshgen.o postprocessing.o problem.o system.o vector.o

tfem_elem.o: element_defs.o element_utils.o eltree.o elvector.o gauss.o kind_defs.o mesh.o problem.o shapefunc.o system.o vector.o

vector.o: element_defs.o kind_defs.o mesh.o pos_array.o problem_defs.o set_optional.o system_defs.o vector_defs.o

vector_defs.o: kind_defs.o

clean:
	$(RM) *.o $(MODPATHWRITE)/*.mod tags $(LIBCREATE)
