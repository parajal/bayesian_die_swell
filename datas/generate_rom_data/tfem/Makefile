#
# Root Makefile of tfem
#

#####################
# Main make targets #
#####################

all: allc doc

allc: common hsl basic_physics $(TFEMADDONS)

common: tfemcore io_utils figplot timer
	ctags --fortran-kinds=+i src/core/*.[fF]90 src/addons/*/*.[fF]90

basic_physics: convection_diffusion stokes

doc: userguide

physics: basic_physics \
         viscoelastic bcf diffuse_interface elastic inertia \
         generalized_stokes viscoelastic_fh surface_advection \
         grid_deformation compressible_fluid projection structures \
         interface_tracking convection_diffusion_supg dfm \
         surface_convection_diffusion linear_elastic interfacial_rheology \
         material_data

clean: tfemcoreclean tfemtestclean io_utilsclean figplotclean \
       timerclean viscoelastic_modelsclean \
       stokesclean viscoelasticclean hslclean sk_solveclean \
       mumpsclean pardisoclean umfpackclean feastclean bcfclean \
       diffuse_interfaceclean shapefunc_extraclean elasticclean \
       convection_diffusionclean inertiaclean \
       metis5clean generalized_stokesclean mkl_gmresclean zigguratclean \
       particle_trackingclean viscoelastic_fhclean \
       hsl_extraclean surface_advectionclean \
       grid_deformationclean supg_utilsclean compressible_fluidclean \
       projectionclean structuresclean interface_trackingclean \
       convection_diffusion_supgclean dfmclean \
       surface_convection_diffusionclean linear_elasticclean \
       interfacial_rheologyclean material_dataclean

allclean: clean libhslclean libhsldummyclean externalclean

tfemcore: Mdefs.mk
	( cd src/core ; make )

tfemcoreclean:
	( cd src/core ; make clean )

tfemtest:
	( cd test ; make )

tfemtestclean:
	( cd test ; make clean )

Mdefs.mk:
	@echo
	@echo Before compiling TFEM you should have an appropriate file
	@echo Mdefs.mk available. Copy the file Mdefs_distr.mk to
	@echo Mdefs.mk and adapt the file to your own computer/compiler.
	@echo
	@exit 1

developer:
	( cd bin; ln -s fmm mm )


#################
# Documentation #
#################

userguide:
	( cd doc/ug ; make )

userguideclean:
	( cd doc/ug ; make clean )


############################################
# Library for default solvers from HSL2002 #
############################################

external/hsl2002/ea23ad.f:
	@echo
	@echo Source files for HSL2002 are not present.
	@echo Unpack the file hsl2002.tar.gz in the TFEM root.
	@echo
	@exit 1

libhsl: Mdefs.mk external/hsl2002/ea23ad.f
	( cd external/hsl2002 ; make )

libhslclean:
	( cd external/hsl2002 ; make clean )

libhsldirclean:
	( cd external/hsl2002 ; make dirclean )

##########
# ADDONS #
##########

SLEEP ?= 1

hsl: tfemcore
	sleep $(SLEEP)
	( cd src/addons/hsl ; make )

hslclean:
	( cd src/addons/hsl ; make clean )

mumps: tfemcore
	sleep $(SLEEP)
	( cd src/addons/mumps ; make )

mumpsclean:
	( cd src/addons/mumps ; make clean )

sk_solve: tfemcore
	sleep $(SLEEP)
	( cd src/addons/sk_solve ; make )

sk_solveclean:
	( cd src/addons/sk_solve ; make clean )

mkl_gmres: tfemcore
	sleep $(SLEEP)
	( cd src/addons/mkl_gmres ; make )

mkl_gmresclean:
	( cd src/addons/mkl_gmres ; make clean )

pardiso: tfemcore
	sleep $(SLEEP)
	( cd src/addons/pardiso ; make )

pardisoclean:
	( cd src/addons/pardiso ; make clean )

umfpack: tfemcore
	sleep $(SLEEP)
	( cd src/addons/umfpack ; make )

umfpackclean:
	( cd src/addons/umfpack ; make clean )

feast: tfemcore
	sleep $(SLEEP)
	( cd src/addons/feast ; make )

feastclean:
	( cd src/addons/feast ; make clean )

io_utils: tfemcore
	sleep $(SLEEP)
	( cd src/addons/io_utils ; make )

io_utilsclean:
	( cd src/addons/io_utils ; make clean )

figplot: tfemcore
	sleep $(SLEEP)
	( cd src/addons/figplot ; make )

figplotclean:
	( cd src/addons/figplot ; make clean )

libhsldummy:
	sleep $(SLEEP)
	( cd src/addons/dummylibs/hsl ; make )

libhsldummyclean:
	( cd src/addons/dummylibs/hsl ; make clean )

timer:
	sleep $(SLEEP)
	( cd src/addons/timer ; make )

timerclean:
	( cd src/addons/timer ; make clean )

stokes: tfemcore
	sleep $(SLEEP)
	( cd src/addons/stokes ; make )

stokesclean:
	( cd src/addons/stokes ; make clean )

viscoelastic_models: tfemcore
	sleep $(SLEEP)
	( cd src/addons/viscoelastic_models ; make )

viscoelastic_modelsclean:
	( cd src/addons/viscoelastic_models ; make clean )

viscoelastic: stokes viscoelastic_models supg_utils
	sleep $(SLEEP)
	( cd src/addons/viscoelastic ; make )

viscoelasticclean:
	( cd src/addons/viscoelastic ; make clean )

convection_diffusion: tfemcore
	sleep $(SLEEP)
	( cd src/addons/convection_diffusion ; make )

convection_diffusionclean:
	( cd src/addons/convection_diffusion ; make clean )

convection_diffusion_supg: stokes supg_utils
	sleep $(SLEEP)
	( cd src/addons/convection_diffusion_supg ; make )

convection_diffusion_supgclean:
	( cd src/addons/convection_diffusion_supg ; make clean )

elastic: tfemcore
	sleep $(SLEEP)
	( cd src/addons/elastic ; make )

elasticclean:
	( cd src/addons/elastic ; make clean )

bcf: viscoelastic shapefunc_extra
	sleep $(SLEEP)
	( cd src/addons/bcf ; make )

bcfclean:
	( cd src/addons/bcf ; make clean )

diffuse_interface: stokes
	sleep $(SLEEP)
	( cd src/addons/diffuse_interface ; make )

diffuse_interfaceclean:
	( cd src/addons/diffuse_interface ; make clean )

shapefunc_extra:
	sleep $(SLEEP)
	( cd src/addons/shapefunc_extra ; make )

shapefunc_extraclean:
	( cd src/addons/shapefunc_extra ; make clean )

inertia: stokes
	sleep $(SLEEP)
	( cd src/addons/inertia ; make )

inertiaclean:
	( cd src/addons/inertia ; make clean )

metis5:
	sleep $(SLEEP)
	( cd src/addons/metis5 ; make )

metis5clean:
	( cd src/addons/metis5 ; make clean )

generalized_stokes: stokes ziggurat
	sleep $(SLEEP)
	( cd src/addons/generalized_stokes ; make )

generalized_stokesclean:
	( cd src/addons/generalized_stokes ; make clean )

ziggurat:
	sleep $(SLEEP)
	( cd src/addons/ziggurat ; make )

zigguratclean:
	( cd src/addons/ziggurat ; make clean )

particle_tracking:
	sleep $(SLEEP)
	( cd src/addons/particle_tracking ; make )

particle_trackingclean:
	( cd src/addons/particle_tracking ; make clean )

viscoelastic_fh: generalized_stokes viscoelastic
	sleep $(SLEEP)
	( cd src/addons/viscoelastic_fh ; make )

viscoelastic_fhclean:
	( cd src/addons/viscoelastic_fh ; make clean )

hsl_extra: tfemcore
	sleep $(SLEEP)
	( cd src/addons/hsl_extra ; make )

hsl_extraclean:
	( cd src/addons/hsl_extra ; make clean )

surface_advection: supg_utils
	sleep $(SLEEP)
	( cd src/addons/surface_advection ; make )

surface_advectionclean:
	( cd src/addons/surface_advection ; make clean )

grid_deformation: convection_diffusion figplot
	sleep $(SLEEP)
	( cd src/addons/grid_deformation ; make )

grid_deformationclean:
	( cd src/addons/grid_deformation ; make clean )

supg_utils:
	sleep $(SLEEP)
	( cd src/addons/supg_utils ; make )

supg_utilsclean:
	( cd src/addons/supg_utils ; make clean )

compressible_fluid: stokes
	sleep $(SLEEP)
	( cd src/addons/compressible_fluid ; make )

compressible_fluidclean:
	( cd src/addons/compressible_fluid ; make clean )

projection:
	sleep $(SLEEP)
	( cd src/addons/projection ; make )

projectionclean:
	( cd src/addons/projection ; make clean )

structures: shapefunc_extra
	sleep $(SLEEP)
	( cd src/addons/structures ; make )

structuresclean:
	( cd src/addons/structures ; make clean )

interface_tracking:
	sleep $(SLEEP)
	( cd src/addons/interface_tracking ; make )

interface_trackingclean:
	( cd src/addons/interface_tracking ; make clean )

surface_convection_diffusion:
	sleep $(SLEEP)
	( cd src/addons/surface_convection_diffusion ; make )

surface_convection_diffusionclean:
	( cd src/addons/surface_convection_diffusion ; make clean )

dfm: stokes shapefunc_extra
	sleep $(SLEEP)
	( cd src/addons/dfm ; make )

dfmclean:
	( cd src/addons/dfm ; make clean )

linear_elastic: tfemcore
	sleep $(SLEEP)
	( cd src/addons/linear_elastic ; make )

linear_elasticclean:
	( cd src/addons/linear_elastic ; make clean )

interfacial_rheology:
	sleep $(SLEEP)
	( cd src/addons/interfacial_rheology ; make )

interfacial_rheologyclean:
	( cd src/addons/interfacial_rheology ; make clean )

material_data: viscoelastic_models
	sleep $(SLEEP)
	( cd src/addons/material_data ; make )

material_dataclean:
	( cd src/addons/material_data ; make clean )


#############################
# External library routines #
#############################

externalclean: libhsl3clean hslfakemetisclean nrclean sparskit2clean

externaldirclean: libhsl3dirclean hslfakemetisdirclean nrdirclean \
  sparskit2dirclean

libhsl2: Mdefs.mk libhslclean
	sleep $(SLEEP)
	( cd external/hsl/ea23 ; make )
	sleep $(SLEEP)
	( cd external/hsl/ma41 ; make )
	sleep $(SLEEP)
	( cd external/hsl/ma57 ; make )

libhsl2clean:
	( cd external/hsl/ea23 ; make clean )
	( cd external/hsl/ma41 ; make clean )
	( cd external/hsl/ma57 ; make clean )

libhsl2dirclean:
	( cd external/hsl/ea23 ; make dirclean )
	( cd external/hsl/ma41 ; make dirclean )
	( cd external/hsl/ma57 ; make dirclean )

hslfakemetis: Mdefs.mk
	sleep $(SLEEP)
	( cd external/hsl/fakemetis ; make )

hslfakemetisclean:
	( cd external/hsl/fakemetis ; make clean )

hslfakemetisdirclean:
	( cd external/hsl/fakemetis ; make dirclean )

libhsl3: Mdefs.mk libhsl2clean libhsl2
	sleep $(SLEEP)
	( cd external/hsl/mc75 ; make )
	sleep $(SLEEP)
	( cd external/hsl/hsl_mi20 ; make )
	sleep $(SLEEP)
	( cd external/hsl/mi24 ; make )
	sleep $(SLEEP)
	( cd external/hsl/mi26 ; make )

libhsl3clean: libhsl2clean
	( cd external/hsl/mc75 ; make clean )
	( cd external/hsl/hsl_mi20 ; make clean )
	( cd external/hsl/mi24 ; make clean )
	( cd external/hsl/mi26 ; make clean )

libhsl3dirclean: libhsl2dirclean
	( cd external/hsl/mc75 ; make dirclean )
	( cd external/hsl/hsl_mi20 ; make dirclean )
	( cd external/hsl/mi24 ; make dirclean )
	( cd external/hsl/mi26 ; make dirclean )

nr: Mdefs.mk
	sleep $(SLEEP)
	( cd external/nr ; make )

nrclean:
	( cd external/nr ; make clean )

nrdirclean:
	( cd external/nr ; make dirclean )

sparskit2: Mdefs.mk
	sleep $(SLEEP)
	( cd external/sparskit2 ; make )

sparskit2clean:
	( cd external/sparskit2 ; make clean )

sparskit2dirclean:
	( cd external/sparskit2 ; make dirclean )

