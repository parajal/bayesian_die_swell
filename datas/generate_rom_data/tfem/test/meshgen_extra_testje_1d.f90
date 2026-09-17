module element_m

  use kind_defs_m
  use mesh_m, only: mesh_t
  use problem_m, only: problem_t
  use element_defs_m, only: coefficients_t, oldvectors_t

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

  subroutine bounelement ( mesh, problem, point, matrix, vector, &
    coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: point
    logical, intent(in) :: matrix, vector
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    elemvec = 3

  end subroutine bounelement

end module element_m

program meshgen_extra_testje_1d

  use tfem_m
  use hsl_ma41_m
  use element_m

  implicit none

  type(meshgen_options_t) :: mesh_options
  type(mesh_t) :: mesh1, mesh2, mesh
  type(input_probdef_t) :: ipd
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(vector_t) :: physvector1, physvector2, physvector3
  type(lu_ma41_t) :: lu

  mesh_options%nx = 10
  mesh_options%elshape = 2

  call line1d ( mesh1, mesh_options )

  mesh_options%ox = 1

  call line1d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh, point1=2, point2=1, mergegroup1=1 )
  call delete ( mesh1, mesh2 )

  !call glue_mesh ( mesh, point1=1, point2=3 )

  call fill_mesh_parts ( mesh )

  call create_input_probdef ( mesh, ipd, nvec=2, nphysq=2 )

  ipd%vec_elementdof(1)%a = reshape ( [2,2,2, 1,0,1], [3,2] )
  ipd%physq = [1,2]

  call define_essential ( mesh, ipd, point=1 )

  call problem_definition ( ipd, mesh, problem )

  call create_sysvector ( problem, sol )

  sol%u = 0

  call fill_sysvector ( mesh, problem, sol, &
    physq=1, point=1, value=2._dp )

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

  call create_sysvector ( problem, rhsd )

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=element, &
    order='ND' )

  call check ( sysmatrix )

  call add_boundary_elements_point ( mesh, problem, rhsd, point=3, &
    elemsub=bounelement, physq=[1] )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma41 ( sysmatrix, rhsd, sol, lu )

  call create_vector ( problem, physvector1, physq=2 )
  call extract_physvector ( mesh, problem, sol, physvector1 )
  call create_vector ( problem, physvector2, physq=2, elementwise=.true. )
  call extract_physvector ( mesh, problem, sol, physvector2 )
  call create_vector ( problem, physvector3, physq=1 )
  call extract_physvector ( mesh, problem, sol, physvector3 )

  call printa

  call delete ( problem )
  call delete ( ipd )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( lu )
  call delete ( physvector1, physvector2, physvector3 )

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
  if ( mesh%ngluepoints > 0 ) then
   print *, 'mesh%gluepoints(:,1)'
   print *, mesh%gluepoints(:,1)
   print *, 'mesh%gluepoints(:,2)'
   print *, mesh%gluepoints(:,2)
  end if
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

  print *, 'physvector3%u'
  print *, physvector3%u

  end subroutine printa

end program meshgen_extra_testje_1d

