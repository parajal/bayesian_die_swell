module subs13_m

  use tfem_elem_m

  implicit none

  integer, parameter :: nparticles=1

  real(dp) :: xp(nparticles,3), rp(nparticles)

contains


! subroutine for freely floating particles (3D)

  subroutine elementc_3D ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer :: nsurfaceconstr
    real(dp) :: xr(1,3), r(3)

    nsurfaceconstr = problem%constraints(constr)%geometry1

    xr(1,:) = mesh%coor(mesh%surfaces(nsurfaceconstr)%nodes(node),:)

!   set shape function in the point

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0

    end if

!   connection through collocation
!   three constraints (vectorial)
!      _   _      _     _   _
!      u - up - omega x r = 0
!
!   where r is the vector from the centerpoint of the particle to the boundary
!
!   or in components
!
!      u - up - omega_y x rz + omega_z x ry= 0
!      v - vp - omega_z x rx + omega_x x rz= 0
!      w - wp - omega_x x ry + omega_y x rx= 0
!
!   the matrix A therefore becomes
!
!     A =  [ phi    0     0]
!          [   0  phi     0]
!          [   0    0   phi]
!
!   and the matrix A_add of the additional unknowns (up,vp,wp,omega_x,
!   omega_y, omega_z)
!
!     A_add =  [ -1  0  0   0 -rz   ry ]
!              [  0 -1  0  rz   0  -rx ]
!              [  0  0 -1 -ry  rx    0 ]

    if ( matrix ) then

      elemmat(:,:) = 0._dp
      elemmat(1,1) = 1._dp
      elemmat(2,2) = 1._dp
      elemmat(3,3) = 1._dp

      r = xr(1,:) - xp(1,:)

      elemmatadd(1,:) = [ -1._dp,   0._dp,   0._dp,  0._dp,  -r(3),   r(2) ]
      elemmatadd(2,:) = [  0._dp,  -1._dp,   0._dp,   r(3),  0._dp,  -r(1) ]
      elemmatadd(3,:) = [  0._dp,   0._dp,  -1._dp,  -r(2),   r(1),  0._dp ]

    end if

  end subroutine elementc_3D

end module subs13_m
