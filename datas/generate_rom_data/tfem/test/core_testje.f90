module element_functions_m

  use kind_defs_m

  implicit none

contains

  function func ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    select case(nr)
      case(1)
        func = 0
      case(2)
        func = (x(1) + 1)**2
      case default
        write(*,'(/a/)') 'Error: wrong function number: '
        stop
    end select

  end function func

end module element_functions_m

module element_m

  use kind_defs_m
  use mesh_m, only: mesh_t
  use problem_m, only: problem_t
  use element_defs_m, only: coefficients_t, oldvectors_t
  use element_functions_m

  implicit none

  logical :: diag = .true.

contains

  subroutine element ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i

    elemmat = 1
    if ( diag ) then
      do i = 1, size(elemmat,1)
        elemmat(i,i) = 2
      end do
    end if
    elemvec = 2

  end subroutine element

  subroutine bounelement ( mesh, problem, curve, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: curve, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    elemvec = 3

  end subroutine bounelement

  subroutine elementc ( mesh, problem, constr, elem, node, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemmat2, elemmatadd, &
    elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    elemmat(1,1)=5
    elemmat = 1
    elemmat2(1,1)=3
    elemmat2 = 1
    elemmatadd = 1
    elemvec = 1
    elemvecadd = 1

  end subroutine elementc

end module element_m

program core_testje

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use problem_m
  use system_m
  use element_m
  use hsl_ma41_m
  use vector_m

  implicit none

  type(meshgen_options_t) :: mesh_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: ipd
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(vector_t) :: physvector1, physvector2
  type(lu_ma41_t) :: lu
  type(solver_options_ma41_t) :: solver_options

  mesh_options%elshape = 6
  mesh_options%nx = 2
  mesh_options%ny = 2

  call quadrilateral2d ( mesh, mesh_options )

! other parts of the mesh

  call add_to_mesh ( mesh, curve=[1,2] )    ! curve 5
  call add_to_mesh ( mesh, curve=[-3] )     ! curve 6
  call add_to_mesh ( mesh, curve=[6,-5] )   ! curve 7
  call add_to_mesh ( mesh, curve=[-4] )     ! curve 8

  call glue_mesh ( mesh, curve1=1, curve2=6 )
  call glue_mesh ( mesh, curve1=2, curve2=8 )

  call fill_mesh_parts ( mesh )

  call create_input_probdef ( mesh, ipd, nvec=2, nphysq=2 )

!  ipd%elementdof(1)%a = (/3,2,3,2,3,2,3,2,1/)
  ipd%vec_elementdof(1)%a = &
    reshape ( [2,2,2,2,2,2,2,2,0, 1,0,1,0,1,0,1,0,1], [9,2] )
  ipd%physq = [1,2]

  call define_essential &
    ( mesh, ipd, degsfd=[2], curve1=1, curve2=2, step=2, exclude=1 )

  call define_essential &
    ( mesh, ipd, physq=2, degfd=[1], curve1=3 )

  call define_constraint &
     ( mesh, ipd, physq=1, curve1=2, curve2=8, nodedof=1, full2=.true., &
      discretization='collocation', exclude=3 )

!  call define_constraint &
!   ( mesh, ipd, physq=1, curve1=1, nglobalc=1 )

  call problem_definition ( ipd, mesh, problem )

  call create_sysvector ( problem, sol )

  sol%u = 0

  call fill_sysvector ( mesh, problem, sol, &
    physq=2, degfd=1, curve1=1, curve2=2, step=2, value=2._dp )

  call fill_sysvector ( mesh, problem, sol, &
    degfd=1, curve1=3, value=3._dp )

  call fill_sysvector ( mesh, problem, sol, &
    degfd=2, curve1=3, func=func, funcnr=2 )

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

  call create_sysvector ( problem, rhsd )

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=element, &
    order='ND' )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    elemsub=elementc, order='ND', addmatvec=.true. )


!    physqcol=(/1,2/), buildmatrix=.true., buildvector=.true. )
!    physqrow=(/1,2/), physqcol=(/1,2/), buildmatrix=.true., buildvector=.true. )
!    physqrow=(/1/), physqcol=(/1/), buildmatrix=.true., buildvector=.true. )
!  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=element, &
!    physqrow=(/2/), physqcol=(/2/), buildmatrix=.true., buildvector=.true., &
!    addmatvec=.true. )
!  diag = .false.
!  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=element, &
!    physqrow=(/1/), physqcol=(/2/), buildmatrix=.true., buildvector=.false., &
!    addmatvec=.true. )
!  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=element, &
!    physqrow=(/2/), physqcol=(/1/), buildmatrix=.true., buildvector=.false., &
!    addmatvec=.true. )

  call check ( sysmatrix )

  call add_boundary_elements ( mesh, problem, rhsd, elemsub=bounelement, &
    curve=2, physq=[1] )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd  )

  !solver_options%printlevel = 1

  solver_options%integer_storage = 1.6

  call solve_system_ma41 ( sysmatrix, rhsd, sol, lu, &
    solver_options=solver_options  )
  !call delete ( lu )
  !call solve_system_ma41 ( sysmatrix, rhsd, sol, lu )

  call create_vector ( problem, physvector1, physq=2 )
  call extract_physvector ( mesh, problem, sol, physvector1 )
  call create_vector ( problem, physvector2, physq=2, elementwise=.true. )
  call extract_physvector ( mesh, problem, sol, physvector2 )

  call printa

  call delete ( problem )
  call delete ( ipd )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( lu )
  call delete ( physvector1, physvector2 )

contains

  subroutine printa

  integer :: i

  print *, 'mesh%ndim'
  print *, mesh%ndim
  print *, 'mesh%nnodes'
  print *, mesh%nnodes
  print *, 'mesh%nelem'
  print *, mesh%nelem
  print *, 'mesh%nelgrp'
  print *, mesh%nelgrp
  print *, 'mesh%npoints'
  print *, mesh%npoints
  print *, 'mesh%ncurves'
  print *, mesh%ncurves
  print *, 'mesh%element(1)%elshape'
  print *, mesh%element(1)%elshape
  print *, 'mesh%element(1)%globalshape'
  print *, mesh%element(1)%globalshape
  print *, 'mesh%element(1)%ndim'
  print *, mesh%element(1)%ndim
  print *, 'mesh%element(1)%numnod'
  print *, mesh%element(1)%numnod
  print *, 'mesh%element(1)%sidnumvert'
  print *, mesh%element(1)%sidnumvert
  print *, 'mesh%element(1)%sidvert'
  print *, mesh%element(1)%sidvert
  print *, 'mesh%elnumnod'
  print *, mesh%elnumnod
  print *, 'mesh%grpnumel'
  print *, mesh%grpnumel
  print *, 'mesh%topology(1)%a'
  print *, mesh%topology(1)%a
  print *, 'mesh%nodnumel'
  print *, mesh%nodnumel
  print *, 'mesh%nodelem(:,1)'
  print *, mesh%nodelem(:,1)
  print *, 'mesh%nodelem(:,2)'
  print *, mesh%nodelem(:,2)
  print *, 'mesh%nodnumnod'
  print *, mesh%nodnumnod
  print *, 'mesh%nodnod'
  print *, mesh%nodnod
  print *, 'mesh%points(:mesh%npoints)'
  print *, mesh%points(:mesh%npoints)
  print *, 'mesh%curves(:mesh%ncurves)%ndim'
  print *, mesh%curves(:mesh%ncurves)%ndim
  print *, 'mesh%curves(:mesh%ncurves)%nnodes'
  print *, mesh%curves(:mesh%ncurves)%nnodes
  print *, 'mesh%curves(:mesh%ncurves)%nelem'
  print *, mesh%curves(:mesh%ncurves)%nelem
  print *, 'mesh%curves(:mesh%ncurves)%elnumnod'
  print *, mesh%curves(:mesh%ncurves)%elnumnod
  print *, 'mesh%curves(:mesh%ncurves)%element%globalshape'
  print *, mesh%curves(:mesh%ncurves)%element%globalshape
  print *, 'mesh%curves(:mesh%ncurves)%element%elshape'
  print *, mesh%curves(:mesh%ncurves)%element%elshape
  print *, 'mesh%curves(:mesh%ncurves)%element%numnod'
  print *, mesh%curves(:mesh%ncurves)%element%numnod
  print *, 'mesh%curves%nodes'
  do i = 1, mesh%ncurves
    print *, mesh%curves(i)%nodes
  end do
  print *, 'mesh%curves%topology(:,:,1)'
  do i = 1, mesh%ncurves
    print *, mesh%curves(i)%topology(:,:,1)
  end do
  print *, 'mesh%curves%topology(:,:,2)'
  do i = 1, mesh%ncurves
    print *, mesh%curves(i)%topology(:,:,2)
  end do
  print *, 'mesh%curves%nc%nodnumel'
  do i = 1, mesh%ncurves
    print *, mesh%curves(i)%nc%nodnumel
  end do
  print *, 'mesh%curves%nc%nodelem'
  do i = 1, mesh%ncurves
    print *, mesh%curves(i)%nc%nodelem
  end do
  print *, 'mesh%curves%nc%nodnumnod'
  do i = 1, mesh%ncurves
    print *, mesh%curves(i)%nc%nodnumnod
  end do
  print *, 'mesh%curves%nc%nodnod'
  do i = 1, mesh%ncurves
    print *, mesh%curves(i)%nc%nodnod
  end do
  print *, 'mesh%sidelem(1)%a(:,:,1)'
  print *, mesh%sidelem(1)%a(:,:,1)
  print *, 'mesh%sidelem(1)%a(:,:,2)'
  print *, mesh%sidelem(1)%a(:,:,2)
  print *, 'mesh%sidelem(1)%a(:,:,3)'
  print *, mesh%sidelem(1)%a(:,:,3)
  print *, 'mesh%sidelem(1)%a(:,:,4)'
  print *, mesh%sidelem(1)%a(:,:,4)
  print *, 'mesh%ngluepoints'
  print *, mesh%ngluepoints
  print *, 'mesh%gluepoints(:,1)'
  print *, mesh%gluepoints(:,1)
  print *, 'mesh%gluepoints(:,2)'
  print *, mesh%gluepoints(:,2)
  print *, 'mesh%coor'
  print *, mesh%coor

  print *, 'problem%numdegfd'
  print *, problem%numdegfd
  print *, 'problem%elnumdegfd(1)%a'
  print *, problem%elnumdegfd(1)%a
  print *, 'ipd%vec_elementdof(1)%a'
  print *, ipd%vec_elementdof(1)%a
  print *, 'ipd%nphysq'
  print *, ipd%nphysq
  print *, 'ipd%physq'
  print *, ipd%physq
  print *, 'problem%vec_elnumdegfd(1)%a'
  print *, problem%vec_elnumdegfd(1)%a
  print *, 'problem%vec_numdegfd'
  print *, problem%vec_numdegfd
  print *, 'problem%nodnumdegfd'
  print *, problem%nodnumdegfd
  print *, 'problem%vec_nodnumdegfd'
  print *, problem%vec_nodnumdegfd
  print *, 'problem%maxnoddegfd'
  print *, problem%maxnoddegfd
  print *, 'problem%numessnodes'
  print *, problem%numessnodes
  print *, 'problem%essnodes(:,1)'
  print *, problem%essnodes(:,1)
  print *, 'problem%essnodes(:,2)'
  print *, problem%essnodes(:,2)
  print *, 'problem%numessdegfd'
  print *, problem%numessdegfd
  print *, 'problem%degfdperm(:,1)'
  print *, problem%degfdperm(:,1)
  print *, 'problem%degfdperm(:,2)'
  print *, problem%degfdperm(:,2)

  print *, 'sol%u'
  print *, sol%u

  print *, 'sysmatrix%Suu%ia'
  print *, sysmatrix%Suu%ia
  print *, 'sysmatrix%Sup%ia'
  print *, sysmatrix%Sup%ia
  print *, 'sysmatrix%Spu%ia'
  print *, sysmatrix%Spu%ia
  print *, 'sysmatrix%Spp%ia'
  print *, sysmatrix%Spp%ia

  print *, 'sysmatrix%Suu%nnz'
  print *, sysmatrix%Suu%nnz
  print *, 'sysmatrix%Sup%nnz'
  print *, sysmatrix%Sup%nnz
  print *, 'sysmatrix%Spu%nnz'
  print *, sysmatrix%Spu%nnz
  print *, 'sysmatrix%Spp%nnz'
  print *, sysmatrix%Spp%nnz

  print *, 'sysmatrix%Suu%a'
  print *, sysmatrix%Suu%a
  print *, 'sysmatrix%Suu%ja'
  print *, sysmatrix%Suu%ja

  print *, 'sysmatrix%Sup%a'
  print *, sysmatrix%Sup%a
  print *, 'sysmatrix%Sup%ja'
  print *, sysmatrix%Sup%ja

  print *, 'sysmatrix%Spu%a'
  print *, sysmatrix%Spu%a
  print *, 'sysmatrix%Spu%ja'
  print *, sysmatrix%Spu%ja

  print *, 'sysmatrix%Spp%a'
  print *, sysmatrix%Spp%a
  print *, 'sysmatrix%Spp%ja'
  print *, sysmatrix%Spp%ja

  print *, 'rhsd%u'
  print *, rhsd%u

  print *, 'physvector1%u'
  print *, physvector1%u

  print *, 'physvector2%u'
  print *, physvector2%u

  end subroutine printa

end program core_testje

