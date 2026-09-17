! A left built-in torsion shaft problem consisting of three parts. 3D beams.
! Different polar second moments of area on parts of the shaft
! using three groups of elements with difference coefficients.
! External torques moments in all points.
! P1/P3/P1 interpolation for axial displacement, transverse displacement
! and axial rotation, respectively.

program structures17

  use tfem_m
  use hsl_ma57_m
  use structures_elements_m
  use structures_post_m
  use io_utils_m

  implicit none


! constants

  integer, parameter :: &
    ninti = 3,  & ! number of Gauss points
    ne(3) = [ 1, 1, 1 ]  ! number of elements in each group

  real(dp), parameter :: &
    Emod = 1._dp, &  ! E modulus
    Gmod = 1._dp, &  ! G modulus
    Ac = 1._dp,   &  ! cross-sectional area
    Izz = 1._dp,  &  ! second moment of area Izz
    Iyy = 2._dp,  &  ! second moment of area Iyy
!   Following example for exam question 7 of Mechanics 4RA00 08-11-2019.
    Jrr(3) = [ 1._dp, 2._dp, 3._dp ],  &  ! polar moment of area J in each part
    L(3) = [ 1._dp, 1._dp, 1._dp ], &  ! length of the beam parts
    T(3) = [ 6._dp, -7._dp, 3._dp ]  ! external torques at end of the parts

! definitions

  type(mesh_t) :: mesh1, mesh2, mesh3, mesh4, mesh5, mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd, reacf, load
  type(meshgen_options_t) :: mesh_options
  type(coefficients_t) :: coefficients(3)

  integer :: i


! fill coefficients

  call create_coefficients ( coefficients(1), ncoefi=100, ncoefr=50 )

  coefficients(1)%i = 0
  coefficients(1)%i(2) = ninti
  coefficients(1)%i(3) = 0 ! P3 Hermite interpolation
  coefficients(1)%i(4) = 0 ! constant distributed force
  coefficients(1)%i(12) = 2 ! P1 interpolation for u
  coefficients(1)%i(13) = 2 ! P1 interpolation for theta

  coefficients(1)%r = 0
  coefficients(1)%r(1) = Ac
  coefficients(1)%r(2) = Emod
  coefficients(1)%r(3) = Izz
  coefficients(1)%r(8:10) = [ 0._dp, 0._dp, 1._dp ] ! local z-axis direction
  coefficients(1)%r(11) = Iyy
  coefficients(1)%r(12) = 0   ! Izy = 0
  coefficients(1)%r(13) = Jrr(1)
  coefficients(1)%r(14) = Gmod

  coefficients(2) = coefficients(1)
  coefficients(2)%r(13) = Jrr(2)

  coefficients(3) = coefficients(1)
  coefficients(3)%r(13) = Jrr(3)


! create mesh

  mesh_options%elshape = 1 ! two-node line elements

  mesh_options%nx = ne(1)  ! number of elements, equidistant
  mesh_options%lx = L(1)   ! length of interval

  call line1d ( mesh1, mesh_options )

  mesh_options%nx = ne(2)  ! number of elements, equidistant
  mesh_options%lx = L(2)   ! length of interval
  mesh_options%ox = L(1)   ! origin of interval

  call line1d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh3, point1=2, point2=1, &
    nogroupmerge=.true. )

  mesh_options%nx = ne(3)  ! number of elements, equidistant
  mesh_options%lx = L(3)   ! length of interval
  mesh_options%ox = L(1)+L(2) ! origin of interval

  call line1d ( mesh4, mesh_options )

  call mesh_merge ( mesh3, mesh4, mesh5, point1=3, point2=1, &
    nogroupmerge=.true. )

  call mesh_convert ( mesh5, mesh, coordinates=[1,0,0] )

  call fill_mesh_parts ( mesh )

  call delete ( mesh1, mesh2, mesh3, mesh4, mesh5 )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  do i = 1, 3
    input_probdef%elementdof(i)%a = 6 ! (u,v,w,phix,phiy,phiz) in two end nodes
    input_probdef%vec_elementdof(i)%a = 1      ! scalar quantity
    input_probdef%vec_elementdof(i)%a(:,2) = 3 ! vector quantity
  end do

! at x=0

  call define_essential ( mesh, input_probdef, point=1 )

  call problem_definition ( input_probdef, mesh, problem )

  call printinfo ( mesh, printlevel=4 )
  call printinfo ( problem, printlevel=4 )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd, reacf, load )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )

! fill rhs vector with discrete forces/monents at end points of beam sections:

  load%u = 0

  do i = 1, 3
    call fill_sysvector ( mesh, problem, load, point=i+1, degfd=4, value=T(i) )
  end do

  call copy ( load, rhsd )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=beam_elem3, &
    order='ND', mcoefficients=coefficients, addvec=.true. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call reaction_forces ( problem, sysmatrix, sol, rhsd, reacf )

! post processing

  call postprocessing_beam3 ( mesh, problem, mcoefficients=coefficients, &
    mshfilename='structures17.mso', sol=sol, load=load, reacf=reacf )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, reacf, load )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program structures17

