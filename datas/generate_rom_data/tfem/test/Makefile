TFEMPATH = ..
include $(TFEMPATH)/Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODINC) $<

SRC = *.f90

OBJS =

MAINP = bcf_couette compressible_fluid_piston1 compressible_fluid_piston2 conv_diff_supg convection_diffusion_1d_spectral convection_diffusion_elvector convection_diffusion_interface convection_diffusion_tensor core_DG_2d core_findrefcoor_ho_triangle core_findrefcoor_spectral2d core_poisson core_poisson_objects core_poisson_objects_build core_poisson_objects_weak core_testje dfm_poiseuille dfm_poiseuille_flow diffuse_interface_2D elastic_3D elastic_stokes_2D generalized_newtonian_newton_raphson generalized_stokes_brownian generalized_stokes_channel generalized_stokes_meshgen_extra_eltree generalized_stokes_meshgen_extra_submesh generalized_stokes_stress_work inertia_CN inertia_channel_2D3D inertia_channel_2D3Daxi inertia_ns2D interface_tracking_2D interface_tracking_2D_impl linear_elastic_2D3Daxi linear_elastic_2D3Daximixed linear_elastic_3D linear_elastic_3D_thermal meshgen_blended_mesh meshgen_blended_mesh2 meshgen_cylinder_in_box meshgen_extra_edgesfaces meshgen_extra_objects meshgen_extra_objects2 meshgen_extra_objects3 meshgen_extra_poisson_objects2 meshgen_extra_poisson_objects3 meshgen_extra_testje_1d meshgen_extra_testmesh poisson_2D poisson_2d_spectral poisson_3D poisson_highorder_2D_triangle stokes_delete_create stokes_embedded_interface stokes_embedded_particle stokes_embeddedbc stokes_global_constraint stokes_global_constraint_elementset stokes_layers stokes_meshgen_extra_connection stokes_meshgen_extra_connection2 stokes_meshgen_extra_eltree stokes_meshgen_extra_submesh stokes_mixed stokes_physqmask stokes_point_nodeset_constraint stokes_spectral_2d stokes_transformation_nodesub stokes_transformation_normalvec structures_beam_1D structures_beam_2D structures_beam_3D structures_truss_3D surface_advection_stokes viscoelastic_3Dchannel viscoelastic_btensor_impl viscoelastic_channel2D3Daxi_impl viscoelastic_channel_prism viscoelastic_concyl_expl viscoelastic_couette_DG viscoelastic_expl viscoelastic_expl_logc viscoelastic_flow_egp_gammap viscoelastic_flow_rotating_cylinder_egp_gammap viscoelastic_fully_implicit_b viscoelastic_impl viscoelastic_implicit_b viscoelastic_models_startup2D_b viscoelastic_models_startup3D_b viscoelastic_models_steady2D viscoelastic_models_steady2D_b viscoelastic_models_steady2Daxi viscoelastic_models_steady3D viscoelastic_models_steady3D_b viscoelastic_small_lambda viscoelastic_sphere_srm0 viscoelastic_sphere_srm1 viscoelastic_sphere_srm2 viscoelastic_steady_flow_b viscoelastic_steady_flow_b_egp_mm viscoelastic_steady_poiseuille_egp viscoelastic_steady_rotating_cylinder_egp viscoelastic_steady_rotating_cylinder_egp2 viscoelastic_stress_work viscoelastic_xfem_tale

all: objects programs tags

objects: $(OBJS)

programs: $(MAINP)

tags: $(SRC)
	ctags --fortran-kinds=+i $(SRC) $(TFEMPATH)/src/core/*.[fF]90 $(TFEMPATH)/src/addons/*/*.[fF]90

bcf_couette: $(OBJS) bcf_couette.o
	$(FC) $(FFLAGS) -o $@ bcf_couette.o $(OBJS) $(LIBS)

compressible_fluid_piston1: $(OBJS) compressible_fluid_piston1.o
	$(FC) $(FFLAGS) -o $@ compressible_fluid_piston1.o $(OBJS) $(LIBS)

compressible_fluid_piston2: $(OBJS) compressible_fluid_piston2.o
	$(FC) $(FFLAGS) -o $@ compressible_fluid_piston2.o $(OBJS) $(LIBS)

conv_diff_supg: $(OBJS) conv_diff_supg.o
	$(FC) $(FFLAGS) -o $@ conv_diff_supg.o $(OBJS) $(LIBS)

convection_diffusion_1d_spectral: $(OBJS) convection_diffusion_1d_spectral.o
	$(FC) $(FFLAGS) -o $@ convection_diffusion_1d_spectral.o $(OBJS) $(LIBS)

convection_diffusion_elvector: $(OBJS) convection_diffusion_elvector.o
	$(FC) $(FFLAGS) -o $@ convection_diffusion_elvector.o $(OBJS) $(LIBS)

convection_diffusion_interface: $(OBJS) convection_diffusion_interface.o
	$(FC) $(FFLAGS) -o $@ convection_diffusion_interface.o $(OBJS) $(LIBS)

convection_diffusion_tensor: $(OBJS) convection_diffusion_tensor.o
	$(FC) $(FFLAGS) -o $@ convection_diffusion_tensor.o $(OBJS) $(LIBS)

core_DG_2d: $(OBJS) core_DG_2d.o
	$(FC) $(FFLAGS) -o $@ core_DG_2d.o $(OBJS) $(LIBS)

core_findrefcoor_ho_triangle: $(OBJS) core_findrefcoor_ho_triangle.o
	$(FC) $(FFLAGS) -o $@ core_findrefcoor_ho_triangle.o $(OBJS) $(LIBS)

core_findrefcoor_spectral2d: $(OBJS) core_findrefcoor_spectral2d.o
	$(FC) $(FFLAGS) -o $@ core_findrefcoor_spectral2d.o $(OBJS) $(LIBS)

core_poisson: $(OBJS) core_poisson.o
	$(FC) $(FFLAGS) -o $@ core_poisson.o $(OBJS) $(LIBS)

core_poisson_objects: $(OBJS) core_poisson_objects.o
	$(FC) $(FFLAGS) -o $@ core_poisson_objects.o $(OBJS) $(LIBS)

core_poisson_objects_build: $(OBJS) core_poisson_objects_build.o
	$(FC) $(FFLAGS) -o $@ core_poisson_objects_build.o $(OBJS) $(LIBS)

core_poisson_objects_weak: $(OBJS) core_poisson_objects_weak.o
	$(FC) $(FFLAGS) -o $@ core_poisson_objects_weak.o $(OBJS) $(LIBS)

core_testje: $(OBJS) core_testje.o
	$(FC) $(FFLAGS) -o $@ core_testje.o $(OBJS) $(LIBS)

dfm_poiseuille: $(OBJS) dfm_poiseuille.o
	$(FC) $(FFLAGS) -o $@ dfm_poiseuille.o $(OBJS) $(LIBS)

dfm_poiseuille_flow: $(OBJS) dfm_poiseuille_flow.o
	$(FC) $(FFLAGS) -o $@ dfm_poiseuille_flow.o $(OBJS) $(LIBS)

diffuse_interface_2D: $(OBJS) diffuse_interface_2D.o
	$(FC) $(FFLAGS) -o $@ diffuse_interface_2D.o $(OBJS) $(LIBS)

elastic_3D: $(OBJS) elastic_3D.o
	$(FC) $(FFLAGS) -o $@ elastic_3D.o $(OBJS) $(LIBS)

elastic_stokes_2D: $(OBJS) elastic_stokes_2D.o
	$(FC) $(FFLAGS) -o $@ elastic_stokes_2D.o $(OBJS) $(LIBS)

generalized_newtonian_newton_raphson: $(OBJS) generalized_newtonian_newton_raphson.o
	$(FC) $(FFLAGS) -o $@ generalized_newtonian_newton_raphson.o $(OBJS) $(LIBS)

generalized_stokes_brownian: $(OBJS) generalized_stokes_brownian.o
	$(FC) $(FFLAGS) -o $@ generalized_stokes_brownian.o $(OBJS) $(LIBS)

generalized_stokes_channel: $(OBJS) generalized_stokes_channel.o
	$(FC) $(FFLAGS) -o $@ generalized_stokes_channel.o $(OBJS) $(LIBS)

generalized_stokes_meshgen_extra_eltree: $(OBJS) generalized_stokes_meshgen_extra_eltree.o
	$(FC) $(FFLAGS) -o $@ generalized_stokes_meshgen_extra_eltree.o $(OBJS) $(LIBS)

generalized_stokes_meshgen_extra_submesh: $(OBJS) generalized_stokes_meshgen_extra_submesh.o
	$(FC) $(FFLAGS) -o $@ generalized_stokes_meshgen_extra_submesh.o $(OBJS) $(LIBS)

generalized_stokes_stress_work: $(OBJS) generalized_stokes_stress_work.o
	$(FC) $(FFLAGS) -o $@ generalized_stokes_stress_work.o $(OBJS) $(LIBS)

inertia_CN: $(OBJS) inertia_CN.o
	$(FC) $(FFLAGS) -o $@ inertia_CN.o $(OBJS) $(LIBS)

inertia_channel_2D3D: $(OBJS) inertia_channel_2D3D.o
	$(FC) $(FFLAGS) -o $@ inertia_channel_2D3D.o $(OBJS) $(LIBS)

inertia_channel_2D3Daxi: $(OBJS) inertia_channel_2D3Daxi.o
	$(FC) $(FFLAGS) -o $@ inertia_channel_2D3Daxi.o $(OBJS) $(LIBS)

inertia_ns2D: $(OBJS) inertia_ns2D.o
	$(FC) $(FFLAGS) -o $@ inertia_ns2D.o $(OBJS) $(LIBS)

interface_tracking_2D: $(OBJS) interface_tracking_2D.o
	$(FC) $(FFLAGS) -o $@ interface_tracking_2D.o $(OBJS) $(LIBS)

interface_tracking_2D_impl: $(OBJS) interface_tracking_2D_impl.o
	$(FC) $(FFLAGS) -o $@ interface_tracking_2D_impl.o $(OBJS) $(LIBS)

linear_elastic_2D3Daxi: $(OBJS) linear_elastic_2D3Daxi.o
	$(FC) $(FFLAGS) -o $@ linear_elastic_2D3Daxi.o $(OBJS) $(LIBS)

linear_elastic_2D3Daximixed: $(OBJS) linear_elastic_2D3Daximixed.o
	$(FC) $(FFLAGS) -o $@ linear_elastic_2D3Daximixed.o $(OBJS) $(LIBS)

linear_elastic_3D: $(OBJS) linear_elastic_3D.o
	$(FC) $(FFLAGS) -o $@ linear_elastic_3D.o $(OBJS) $(LIBS)

linear_elastic_3D_thermal: $(OBJS) linear_elastic_3D_thermal.o
	$(FC) $(FFLAGS) -o $@ linear_elastic_3D_thermal.o $(OBJS) $(LIBS)

meshgen_blended_mesh: $(OBJS) meshgen_blended_mesh.o
	$(FC) $(FFLAGS) -o $@ meshgen_blended_mesh.o $(OBJS) $(LIBS)

meshgen_blended_mesh2: $(OBJS) meshgen_blended_mesh2.o
	$(FC) $(FFLAGS) -o $@ meshgen_blended_mesh2.o $(OBJS) $(LIBS)

meshgen_cylinder_in_box: $(OBJS) meshgen_cylinder_in_box.o
	$(FC) $(FFLAGS) -o $@ meshgen_cylinder_in_box.o $(OBJS) $(LIBS)

meshgen_extra_edgesfaces: $(OBJS) meshgen_extra_edgesfaces.o
	$(FC) $(FFLAGS) -o $@ meshgen_extra_edgesfaces.o $(OBJS) $(LIBS)

meshgen_extra_objects: $(OBJS) meshgen_extra_objects.o
	$(FC) $(FFLAGS) -o $@ meshgen_extra_objects.o $(OBJS) $(LIBS)

meshgen_extra_objects2: $(OBJS) meshgen_extra_objects2.o
	$(FC) $(FFLAGS) -o $@ meshgen_extra_objects2.o $(OBJS) $(LIBS)

meshgen_extra_objects3: $(OBJS) meshgen_extra_objects3.o
	$(FC) $(FFLAGS) -o $@ meshgen_extra_objects3.o $(OBJS) $(LIBS)

meshgen_extra_poisson_objects2: $(OBJS) meshgen_extra_poisson_objects2.o
	$(FC) $(FFLAGS) -o $@ meshgen_extra_poisson_objects2.o $(OBJS) $(LIBS)

meshgen_extra_poisson_objects3: $(OBJS) meshgen_extra_poisson_objects3.o
	$(FC) $(FFLAGS) -o $@ meshgen_extra_poisson_objects3.o $(OBJS) $(LIBS)

meshgen_extra_testje_1d: $(OBJS) meshgen_extra_testje_1d.o
	$(FC) $(FFLAGS) -o $@ meshgen_extra_testje_1d.o $(OBJS) $(LIBS)

meshgen_extra_testmesh: $(OBJS) meshgen_extra_testmesh.o
	$(FC) $(FFLAGS) -o $@ meshgen_extra_testmesh.o $(OBJS) $(LIBS)

poisson_2D: $(OBJS) poisson_2D.o
	$(FC) $(FFLAGS) -o $@ poisson_2D.o $(OBJS) $(LIBS)

poisson_2d_spectral: $(OBJS) poisson_2d_spectral.o
	$(FC) $(FFLAGS) -o $@ poisson_2d_spectral.o $(OBJS) $(LIBS)

poisson_3D: $(OBJS) poisson_3D.o
	$(FC) $(FFLAGS) -o $@ poisson_3D.o $(OBJS) $(LIBS)

poisson_highorder_2D_triangle: $(OBJS) poisson_highorder_2D_triangle.o
	$(FC) $(FFLAGS) -o $@ poisson_highorder_2D_triangle.o $(OBJS) $(LIBS)

stokes_delete_create: $(OBJS) stokes_delete_create.o
	$(FC) $(FFLAGS) -o $@ stokes_delete_create.o $(OBJS) $(LIBS)

stokes_embedded_interface: $(OBJS) stokes_embedded_interface.o
	$(FC) $(FFLAGS) -o $@ stokes_embedded_interface.o $(OBJS) $(LIBS)

stokes_embedded_particle: $(OBJS) stokes_embedded_particle.o
	$(FC) $(FFLAGS) -o $@ stokes_embedded_particle.o $(OBJS) $(LIBS)

stokes_embeddedbc: $(OBJS) stokes_embeddedbc.o
	$(FC) $(FFLAGS) -o $@ stokes_embeddedbc.o $(OBJS) $(LIBS)

stokes_global_constraint: $(OBJS) stokes_global_constraint.o
	$(FC) $(FFLAGS) -o $@ stokes_global_constraint.o $(OBJS) $(LIBS)

stokes_global_constraint_elementset: $(OBJS) stokes_global_constraint_elementset.o
	$(FC) $(FFLAGS) -o $@ stokes_global_constraint_elementset.o $(OBJS) $(LIBS)

stokes_layers: $(OBJS) stokes_layers.o
	$(FC) $(FFLAGS) -o $@ stokes_layers.o $(OBJS) $(LIBS)

stokes_meshgen_extra_connection: $(OBJS) stokes_meshgen_extra_connection.o
	$(FC) $(FFLAGS) -o $@ stokes_meshgen_extra_connection.o $(OBJS) $(LIBS)

stokes_meshgen_extra_connection2: $(OBJS) stokes_meshgen_extra_connection2.o
	$(FC) $(FFLAGS) -o $@ stokes_meshgen_extra_connection2.o $(OBJS) $(LIBS)

stokes_meshgen_extra_eltree: $(OBJS) stokes_meshgen_extra_eltree.o
	$(FC) $(FFLAGS) -o $@ stokes_meshgen_extra_eltree.o $(OBJS) $(LIBS)

stokes_meshgen_extra_submesh: $(OBJS) stokes_meshgen_extra_submesh.o
	$(FC) $(FFLAGS) -o $@ stokes_meshgen_extra_submesh.o $(OBJS) $(LIBS)

stokes_mixed: $(OBJS) stokes_mixed.o
	$(FC) $(FFLAGS) -o $@ stokes_mixed.o $(OBJS) $(LIBS)

stokes_physqmask: $(OBJS) stokes_physqmask.o
	$(FC) $(FFLAGS) -o $@ stokes_physqmask.o $(OBJS) $(LIBS)

stokes_point_nodeset_constraint: $(OBJS) stokes_point_nodeset_constraint.o
	$(FC) $(FFLAGS) -o $@ stokes_point_nodeset_constraint.o $(OBJS) $(LIBS)

stokes_spectral_2d: $(OBJS) stokes_spectral_2d.o
	$(FC) $(FFLAGS) -o $@ stokes_spectral_2d.o $(OBJS) $(LIBS)

stokes_transformation_nodesub: $(OBJS) stokes_transformation_nodesub.o
	$(FC) $(FFLAGS) -o $@ stokes_transformation_nodesub.o $(OBJS) $(LIBS)

stokes_transformation_normalvec: $(OBJS) stokes_transformation_normalvec.o
	$(FC) $(FFLAGS) -o $@ stokes_transformation_normalvec.o $(OBJS) $(LIBS)

structures_beam_1D: $(OBJS) structures_beam_1D.o
	$(FC) $(FFLAGS) -o $@ structures_beam_1D.o $(OBJS) $(LIBS)

structures_beam_2D: $(OBJS) structures_beam_2D.o
	$(FC) $(FFLAGS) -o $@ structures_beam_2D.o $(OBJS) $(LIBS)

structures_beam_3D: $(OBJS) structures_beam_3D.o
	$(FC) $(FFLAGS) -o $@ structures_beam_3D.o $(OBJS) $(LIBS)

structures_truss_3D: $(OBJS) structures_truss_3D.o
	$(FC) $(FFLAGS) -o $@ structures_truss_3D.o $(OBJS) $(LIBS)

surface_advection_stokes: $(OBJS) surface_advection_stokes.o
	$(FC) $(FFLAGS) -o $@ surface_advection_stokes.o $(OBJS) $(LIBS)

viscoelastic_3Dchannel: $(OBJS) viscoelastic_3Dchannel.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_3Dchannel.o $(OBJS) $(LIBS)

viscoelastic_btensor_impl: $(OBJS) viscoelastic_btensor_impl.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_btensor_impl.o $(OBJS) $(LIBS)

viscoelastic_channel2D3Daxi_impl: $(OBJS) viscoelastic_channel2D3Daxi_impl.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_channel2D3Daxi_impl.o $(OBJS) $(LIBS)

viscoelastic_channel_prism: $(OBJS) viscoelastic_channel_prism.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_channel_prism.o $(OBJS) $(LIBS)

viscoelastic_concyl_expl: $(OBJS) viscoelastic_concyl_expl.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_concyl_expl.o $(OBJS) $(LIBS)

viscoelastic_couette_DG: $(OBJS) viscoelastic_couette_DG.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_couette_DG.o $(OBJS) $(LIBS)

viscoelastic_expl: $(OBJS) viscoelastic_expl.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_expl.o $(OBJS) $(LIBS)

viscoelastic_expl_logc: $(OBJS) viscoelastic_expl_logc.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_expl_logc.o $(OBJS) $(LIBS)

viscoelastic_flow_egp_gammap: $(OBJS) viscoelastic_flow_egp_gammap.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_flow_egp_gammap.o $(OBJS) $(LIBS)

viscoelastic_flow_rotating_cylinder_egp_gammap: $(OBJS) viscoelastic_flow_rotating_cylinder_egp_gammap.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_flow_rotating_cylinder_egp_gammap.o $(OBJS) $(LIBS)

viscoelastic_fully_implicit_b: $(OBJS) viscoelastic_fully_implicit_b.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_fully_implicit_b.o $(OBJS) $(LIBS)

viscoelastic_impl: $(OBJS) viscoelastic_impl.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_impl.o $(OBJS) $(LIBS)

viscoelastic_implicit_b: $(OBJS) viscoelastic_implicit_b.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_implicit_b.o $(OBJS) $(LIBS)

viscoelastic_models_startup2D_b: $(OBJS) viscoelastic_models_startup2D_b.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_models_startup2D_b.o $(OBJS) $(LIBS)

viscoelastic_models_startup3D_b: $(OBJS) viscoelastic_models_startup3D_b.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_models_startup3D_b.o $(OBJS) $(LIBS)

viscoelastic_models_steady2D: $(OBJS) viscoelastic_models_steady2D.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_models_steady2D.o $(OBJS) $(LIBS)

viscoelastic_models_steady2D_b: $(OBJS) viscoelastic_models_steady2D_b.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_models_steady2D_b.o $(OBJS) $(LIBS)

viscoelastic_models_steady2Daxi: $(OBJS) viscoelastic_models_steady2Daxi.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_models_steady2Daxi.o $(OBJS) $(LIBS)

viscoelastic_models_steady3D: $(OBJS) viscoelastic_models_steady3D.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_models_steady3D.o $(OBJS) $(LIBS)

viscoelastic_models_steady3D_b: $(OBJS) viscoelastic_models_steady3D_b.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_models_steady3D_b.o $(OBJS) $(LIBS)

viscoelastic_small_lambda: $(OBJS) viscoelastic_small_lambda.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_small_lambda.o $(OBJS) $(LIBS)

viscoelastic_sphere_srm0: $(OBJS) viscoelastic_sphere_srm0.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_sphere_srm0.o $(OBJS) $(LIBS)

viscoelastic_sphere_srm1: $(OBJS) viscoelastic_sphere_srm1.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_sphere_srm1.o $(OBJS) $(LIBS)

viscoelastic_sphere_srm2: $(OBJS) viscoelastic_sphere_srm2.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_sphere_srm2.o $(OBJS) $(LIBS)

viscoelastic_steady_flow_b: $(OBJS) viscoelastic_steady_flow_b.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_steady_flow_b.o $(OBJS) $(LIBS)

viscoelastic_steady_flow_b_egp_mm: $(OBJS) viscoelastic_steady_flow_b_egp_mm.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_steady_flow_b_egp_mm.o $(OBJS) $(LIBS)

viscoelastic_steady_poiseuille_egp: $(OBJS) viscoelastic_steady_poiseuille_egp.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_steady_poiseuille_egp.o $(OBJS) $(LIBS)

viscoelastic_steady_rotating_cylinder_egp: $(OBJS) viscoelastic_steady_rotating_cylinder_egp.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_steady_rotating_cylinder_egp.o $(OBJS) $(LIBS)

viscoelastic_steady_rotating_cylinder_egp2: $(OBJS) viscoelastic_steady_rotating_cylinder_egp2.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_steady_rotating_cylinder_egp2.o $(OBJS) $(LIBS)

viscoelastic_stress_work: $(OBJS) viscoelastic_stress_work.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_stress_work.o $(OBJS) $(LIBS)

viscoelastic_xfem_tale: $(OBJS) viscoelastic_xfem_tale.o
	$(FC) $(FFLAGS) -o $@ viscoelastic_xfem_tale.o $(OBJS) $(LIBS)

clean:
	$(RM)  bcf_couette compressible_fluid_piston1 compressible_fluid_piston2 conv_diff_supg convection_diffusion_1d_spectral convection_diffusion_elvector convection_diffusion_interface convection_diffusion_tensor core_DG_2d core_findrefcoor_ho_triangle core_findrefcoor_spectral2d core_poisson core_poisson_objects core_poisson_objects_build core_poisson_objects_weak core_testje dfm_poiseuille dfm_poiseuille_flow diffuse_interface_2D elastic_3D elastic_stokes_2D generalized_newtonian_newton_raphson generalized_stokes_brownian generalized_stokes_channel generalized_stokes_meshgen_extra_eltree generalized_stokes_meshgen_extra_submesh generalized_stokes_stress_work inertia_CN inertia_channel_2D3D inertia_channel_2D3Daxi inertia_ns2D interface_tracking_2D interface_tracking_2D_impl linear_elastic_2D3Daxi linear_elastic_2D3Daximixed linear_elastic_3D linear_elastic_3D_thermal meshgen_blended_mesh meshgen_blended_mesh2 meshgen_cylinder_in_box meshgen_extra_edgesfaces meshgen_extra_objects meshgen_extra_objects2 meshgen_extra_objects3 meshgen_extra_poisson_objects2 meshgen_extra_poisson_objects3 meshgen_extra_testje_1d meshgen_extra_testmesh poisson_2D poisson_2d_spectral poisson_3D poisson_highorder_2D_triangle stokes_delete_create stokes_embedded_interface stokes_embedded_particle stokes_embeddedbc stokes_global_constraint stokes_global_constraint_elementset stokes_layers stokes_meshgen_extra_connection stokes_meshgen_extra_connection2 stokes_meshgen_extra_eltree stokes_meshgen_extra_submesh stokes_mixed stokes_physqmask stokes_point_nodeset_constraint stokes_spectral_2d stokes_transformation_nodesub stokes_transformation_normalvec structures_beam_1D structures_beam_2D structures_beam_3D structures_truss_3D surface_advection_stokes viscoelastic_3Dchannel viscoelastic_btensor_impl viscoelastic_channel2D3Daxi_impl viscoelastic_channel_prism viscoelastic_concyl_expl viscoelastic_couette_DG viscoelastic_expl viscoelastic_expl_logc viscoelastic_flow_egp_gammap viscoelastic_flow_rotating_cylinder_egp_gammap viscoelastic_fully_implicit_b viscoelastic_impl viscoelastic_implicit_b viscoelastic_models_startup2D_b viscoelastic_models_startup3D_b viscoelastic_models_steady2D viscoelastic_models_steady2D_b viscoelastic_models_steady2Daxi viscoelastic_models_steady3D viscoelastic_models_steady3D_b viscoelastic_small_lambda viscoelastic_sphere_srm0 viscoelastic_sphere_srm1 viscoelastic_sphere_srm2 viscoelastic_steady_flow_b viscoelastic_steady_flow_b_egp_mm viscoelastic_steady_poiseuille_egp viscoelastic_steady_rotating_cylinder_egp viscoelastic_steady_rotating_cylinder_egp2 viscoelastic_stress_work viscoelastic_xfem_tale *.o *.mod *.bin tags
