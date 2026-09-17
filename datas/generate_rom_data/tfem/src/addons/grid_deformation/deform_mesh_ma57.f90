! Performing grid deformation according to:

! Decheng Wan and Stefan Turek, "Fictitious boundary and moving mesh methods
! for the numerical simulation of rigid particulate flows", J. Comp. Phys.,
! 222 (2007) 28-56

! M. Grajewski, M. Kester, and S. Turek., "Mathematical and numerical analysis
! of a robust and efficient grid deformation method in the finite element
! context", SIAM Journal on Scientific Computing, 31 (2009) 1539-1557


! The deform_mesh routine in this modules assumes:
!    1) Optional: periodical boundaries defined by collocated constraints in the
!       user routine constraints_definition. Also the routine move_object_points
!       needs to be defined.
!    2) Zero valued Neumann boundary conditions on other parts of the boundary.
!    3) A single value specified in a point. The mesh should contain at least
!       one point.
!    4) MA57 solver.

! NOTES:
! - Make a special version by modifying the source below to change 4).


module deform_mesh_ma57_m

  use tfem_m
  use hsl_ma57_m
  use convection_diffusion_elements_m
  use meshgen_grid_deformation_m

  implicit none

contains


! routine for performing grid deformation

  subroutine deform_mesh_ma57 ( mesh, object, problem, coefficients, &
    sysmatrix, f_monitor, numgridsteps, solver_options, lu, &
    constraint_definition, move_object_points, project_boundary, buildmatrix, &
    plotfigs )

    use timer_m
    use figplot_m

!   On entry, this contains the initial mesh.
!   On exit, the coordinates in the object (mesh%objects(object)coor) contains
!   the coordinates of deformed mesh.
    type(mesh_t), intent(inout) :: mesh

!   The object that is being used for obtaining the deformed coordinates.
!   If it doesn't exist (like object=0) the object is created.
!   If it does exist, but the actual size is insufficient it is replaced.
    integer, intent(inout) :: object

!   The Poisson problem defined on the domain.
!   Is created if empty on call and kept at output
    type(problem_t), intent(inout) :: problem

!   The coefficients for the Poisson problem. At least the following needs
!   to be filled:
!      call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )
!      coefficients%i = 0
!      coefficients%i(1) = uintpl
!      coefficients%i(10) = gauss
!      coefficients%r = 0
    type(coefficients_t), intent(in) :: coefficients

!   The system matrix of the Poisson problem defined on the domain.
!   Is created if empty on call and kept at output
    type(sysmatrix_t), intent(inout) :: sysmatrix

!   Number of pseudo time steps on the time interval [0,1]
    integer, intent(in) :: numgridsteps

!   The monitor function in the nodes of the mesh. The value of the monitor
!   function gives the relative size (area/volume) of elements with respect to
!   each other. The absolute value is not important.
    real(dp), dimension(:), intent(in) :: f_monitor

!   solver options for the MA57 solver
    type(solver_options_ma57_t), intent(in), optional :: solver_options

!   LU factorization (MA57) of the system matrix
!   When present, lu is used to store the LU factorization. The next time the
!   routine is called with lu present only a backsubstition is done. This
!   happens even when the matrix has changed. If lu is not present or
!   has been deallocated using delete, a full factorization is performed.
    type(lu_ma57_t), intent(inout), optional :: lu

!   Define constraints for the Poisson problem.
!   NOTE: move_object_points also needs to be present.
    optional constraint_definition
    interface
      subroutine constraint_definition ( mesh, input_probdef )
        use mesh_m
        use problem_defs_m
        implicit none
        type(mesh_t), intent(inout) :: mesh
        type(input_probdef_t), intent(inout) :: input_probdef
      end subroutine constraint_definition
    end interface

!   move object points before intersection with the mesh
!   NOTE: this only affects the coordinates before intersecting. After the time
!   step the coordinates are updated correctly. This is for handling periodical
!   domains.
    optional move_object_points
    interface
      subroutine move_object_points ( mesh, object )
        use mesh_m
        implicit none
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: object
      end subroutine move_object_points
    end interface

!   correct for points moved (numerically) out of the domain at the Neumann
!   boundaries after a time step.
    optional project_boundary
    interface
      subroutine project_boundary ( mesh, object )
        use mesh_m
        implicit none
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: object
      end subroutine project_boundary
    end interface

!   if assigned the value .false. the assembling of the matrix is not done.
    logical, intent(in), optional :: buildmatrix

!   if assigned the value .false. some plots with figplot are made
    logical, intent(in), optional :: plotfigs


!   Local definitions

    type(input_probdef_t) :: input_probdef
    type(oldvectors_t) :: oldvectors
    type(sysvector_t), target :: sol
    type(sysvector_t) :: rhsd
    type(plot_options_t) :: plot_options

!   The area and monitor function (both unscaled and scaled)
    type(vector_t), target :: garea, fmon

!   The grid velocity
    type(vector_t) :: grid_velocity

    integer :: elgrp

    real(dp) :: resultsum(2)

!   Output some plots if .true.
    logical :: plot


    plot = set_optional( variable=plotfigs, default=.false. )


!   test constraints

    if ( present(constraint_definition) ) then
      if ( .not. present(move_object_points) ) then
        write(*,'(/2(a/))') 'Error deform_mesh_ma57: ', &
          ' For constraints move_object_points needs to be present.'
        stop
      end if
    end if

!   object for computing the deformed mesh

    if ( object < 1 .or. object > mesh%nobjects ) then

!     add object

      call add_to_mesh ( mesh, object='coordinates', coor=mesh%coor, &
        warn_mesh_parts=.false. )
      object = mesh%nobjects
      call fill_mesh_parts_objects ( mesh, object1=object )

    else if ( mesh%objects(object)%nnodes /= mesh%nnodes .and. &
              mesh%objects(object)%ndim /= mesh%ndim ) then

!     replace with object of the correct size

      call add_to_mesh ( mesh, object='coordinates', coor=mesh%coor, &
        replace=object, warn_mesh_parts=.false. )
      call fill_mesh_parts_objects ( mesh, object1=object )

    else

!     initialize coordinates of the object

      mesh%objects(object)%coor = mesh%coor

    end if

!   problem definition

    if ( .not. problem%created ) then

      call create_input_probdef ( mesh, input_probdef, nvec=2 )

      do elgrp = 1, mesh%nelgrp
        input_probdef%elementdof(elgrp)%a = 1
        input_probdef%vec_elementdof(elgrp)%a(:,1) = 1
        input_probdef%vec_elementdof(elgrp)%a(:,2) = mesh%ndim
      end do

!     Dirichlet in P1
      call define_essential ( mesh, input_probdef, point=1 )

      if ( present(constraint_definition) ) then
        call constraint_definition ( mesh, input_probdef )
      end if

      call problem_definition ( input_probdef, mesh, problem )

      call delete ( input_probdef )

    end if

!   create system vectors (solution and right-hand side)

    call create ( problem, sol, rhsd )

!   fill solution vector with zero essential boundary conditions

    call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )

!   create oldvectors

    call create_oldvectors ( oldvectors, nsysvec=1, nvec=2 )


    call tic

!   area function

    call create_vector ( problem, garea, vec=1 )

!   derive area function

    call derive_vector ( mesh, problem, garea, elemsub=poisson_area_deriv, &
      coefficients=coefficients, oldvectors=oldvectors )

    oldvectors%v(1)%p => garea

    call integrate ( mesh, problem, resultsum, &
      elemsub=integrate_inverse_griddef_elem, &
      coefficients=coefficients, oldvectors=oldvectors )

!   scaled g such that
!         /
!         | 1/g dx = Omega
!         /Omega
!   or average(1/g)=1.

    garea%u = resultsum(1)/resultsum(2)*garea%u

    if ( plot ) then

      call plot_color_contour ( plot_options, mesh, problem, 'garea.fig', &
        vector=garea )

    end if


    call toc ('area function')


!   monitor function

    call create_vector ( problem, fmon, vec=1 )

    fmon%u = f_monitor

    if ( plot ) then

      call plot_color_contour ( plot_options, mesh, problem, 'fmon.fig', &
        vector=fmon )

    end if

    oldvectors%v(1)%p => fmon

    call integrate ( mesh, problem, resultsum, &
      elemsub=integrate_inverse_griddef_elem, &
      coefficients=coefficients, oldvectors=oldvectors )

!   scaled f such that
!         /
!         | 1/f dx = Omega
!         /Omega
!   or average(1/f)=1.
!   Note that
!         /
!         | ( 1/f - 1/g ) dx = 0
!         /Omega
!   for consistency with the zero flux on the boundary.

    fmon%u = resultsum(1)/resultsum(2)*fmon%u

    if ( plot ) then

      call plot_color_contour ( plot_options, mesh, problem, &
        'fmon_scaled.fig', vector=fmon )

    end if

    call toc ('monitor function')


    if ( .not. sysmatrix%allocated_data ) then

!     create system matrix

      if ( present(constraint_definition) ) then

        call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
          symmetric=.true. )
        call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
        call finalize_sysmatrix_structure ( sysmatrix )

      else

        call create_sysmatrix_structure ( sysmatrix, mesh, problem, &
          symmetric=.true. )

      end if

      call create_sysmatrix_data ( sysmatrix )

    end if

    oldvectors%v(1)%p => garea
    oldvectors%v(2)%p => fmon

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=poisson_griddef_elem, coefficients=coefficients, &
      oldvectors=oldvectors, buildmatrix=buildmatrix )

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      elemsub=poisson_node_conn, addmatvec=.true., coefficients=coefficients  )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    oldvectors%v(1)%p => null()
    oldvectors%v(2)%p => null()

    call toc ('build system')


!   MA57 solver

    call solve_system_ma57 ( sysmatrix, rhsd, sol, lu, solver_options )

    if ( plot ) then

      call plot_color_contour ( plot_options, mesh, problem, 'sol.fig', &
        sysvector=sol )

    end if

    call delete ( rhsd )

    call toc ('solve')


!   derive mesh velocity

    call create_vector ( problem, grid_velocity, vec=2 )

    oldvectors%s(1)%p => sol

    call derive_vector ( mesh, problem, grid_velocity, elemsub=poisson_deriv, &
                         coefficients=coefficients, oldvectors=oldvectors )

    oldvectors%s(1)%p => null()

    call delete ( sol )

    if ( plot ) then

      call plot_vector ( plot_options, mesh, problem, 'grid_velocity.fig', &
        vector=grid_velocity )
      call plot_color_contour ( plot_options, mesh, problem, &
        'grid_velocity_color_u.fig', vector=grid_velocity, degfd=1 )
      call plot_color_contour ( plot_options, mesh, problem, &
        'grid_velocity_color_v.fig', vector=grid_velocity, degfd=2 )

    end if

    call toc ('derive gridvelocity')


!   compute deformed mesh

    call find_grid_coordinates ( mesh, object, coefficients, &
      grid_velocity, fmon=fmon, garea=garea, numgridsteps=numgridsteps, &
      project_boundary=project_boundary, move_object_points=move_object_points )

    call delete ( grid_velocity, fmon, garea )

    call toc ('compute deformed mesh')


!   delete some remaining structures

    call delete ( oldvectors )


  end subroutine deform_mesh_ma57


end module deform_mesh_ma57_m
