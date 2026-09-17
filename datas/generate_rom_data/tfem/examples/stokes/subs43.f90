module subs43_m

use tfem_elem_m

  implicit none

  real(dp) :: force(3), torque(3)

contains

! element for rigid-body motion

  subroutine elementc ( mesh, problem, constr, elem, node, &
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

    integer :: nsurfconstr
    real(dp) :: xr(1,3)

    nsurfconstr = problem%constraints(constr)%geometry1

    xr(1,:) = mesh%coor(mesh%surfaces(nsurfconstr)%nodes(node),:)

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0
      if ( first ) then
        elemvecadd(1:3) = force
        elemvecadd(4:6) = torque
      end if

    end if

!   connection through collocation

    elemmat(1,1) = 1
    elemmat(1,2) = 0
    elemmat(1,3) = 0

    elemmat(2,1) = 0
    elemmat(2,2) = 1
    elemmat(2,3) = 0

    elemmat(3,1) = 0
    elemmat(3,2) = 0
    elemmat(3,3) = 1

    elemmatadd(1,:) = [ -1._dp,  0._dp,  0._dp,    0._dp,  xr(1,3), -xr(1,2) ]
    elemmatadd(2,:) = [  0._dp, -1._dp,  0._dp, -xr(1,3),    0._dp,  xr(1,1) ]
    elemmatadd(3,:) = [  0._dp,  0._dp, -1._dp,  xr(1,2), -xr(1,1),  0._dp   ]

  end subroutine elementc

! element for normal velocity coupling

  subroutine stokes_constr_normal_vel ( mesh, problem, constr, elem, &
    node, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemmat2, elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer, parameter :: ninti = 16, p = 8, nodalp = 6, ndf = 6, ndim = 3, &
                          ndfl = 6

    integer :: i, j, surface
    real(dp) :: phi(ninti,ndf), psi(ninti,ndfl), xig(ninti,2), x(nodalp,ndim)
    real(dp) :: wg(ninti), dphi(ninti,ndf,2), normal(ninti,ndim)
    real(dp) :: dxdxis(ninti,ndim,2), surfl(ninti)
    real(dp) :: tmp_x(ndfl,ndf), tmp_y(ndfl,ndf), tmp_z(ndfl,ndf)

!   connection through elements (weak)

    surface = problem%constraints(constr)%geometry1

!   set Gauss integration and shape function

    call Gauss_Legendre_numeric_triangle ( p, xig, wg )
    call shape_triangle_P2 ( xig, phi, dphi )
    call shape_triangle_P2 ( xig, psi )

!   compute geometry of element deformed

    call get_coordinates_geometry ( mesh, elem, x, surface=surface )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, normal )

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0

    end if

    if ( matrix ) then

      do i = 1, ndfl
        do j = 1, ndf
          tmp_x(i,j) = sum ( psi(:,i) * phi(:,j) * surfl * normal(:,1)* wg )
          tmp_y(i,j) = sum ( psi(:,i) * phi(:,j) * surfl * normal(:,2)* wg )
          tmp_z(i,j) = sum ( psi(:,i) * phi(:,j) * surfl * normal(:,3)* wg )
        end do
      end do

      elemmat = 0
      elemmat(1:ndfl,1:ndf) = tmp_x
      elemmat(1:ndfl,ndf+1:2*ndf) = tmp_y
      elemmat(1:ndfl,2*ndf+1:3*ndf) = tmp_z

      elemmat2 = - elemmat

    end if

  end subroutine stokes_constr_normal_vel


! element for connection: slip on particle surface

  subroutine elementcn ( mesh, problem, conn, elem, node, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat11, elemmat12, &
    elemmat21, elemmat22, elemvec1, elemvec2 )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: conn, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat11, elemmat12, &
      elemmat21, elemmat22
    real(dp), intent(out), dimension(:) :: elemvec1, elemvec2

    integer, parameter :: ninti = 16, p = 8, nodalp = 6, ndf = 6, ndim = 3

    integer ::  i, j, surface

    real(dp) :: phi1(ninti,ndf), phi2(ninti,ndf), xig(ninti,2), x(nodalp,ndim)
    real(dp) :: wg(ninti), dphi(ninti,ndf,2), normal(ninti,ndim)
    real(dp) :: dxdxis(ninti,ndim,2), surfl(ninti), Inn(ninti,ndim,ndim)
    real(dp) :: etaslip

    surface = problem%connections(conn)%geometry1

!   set Gauss integration and shape function

    call Gauss_Legendre_numeric_triangle ( p, xig, wg )
    call shape_triangle_P2 ( xig, phi1, dphi )
    call shape_triangle_P2 ( xig, phi2 )

!   compute geometry of element deformed

    call get_coordinates_geometry ( mesh, elem, x, surface=surface )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, normal )

    if ( vector ) then

      elemvec1 = 0
      elemvec2 = 0

    end if

    if ( matrix ) then

!     diagonal blocks of I-nn

      do i = 1, ndim
        Inn(:,i,i) = 1 - normal(:,i)**2
      end do

!     off-diagonal blocks of I-nn

      do i = 1, ndim-1
        do j = i+1, ndim
          Inn(:,i,j) = - normal(:,i)*normal(:,j)
          Inn(:,j,i) = Inn(:,i,j)
        end do
      end do

      etaslip = coefficients%r(501)

      elemmat11 = 0
      elemmat12 = 0
      elemmat21 = 0
      elemmat22 = 0

      do i = 1, ndf
        do j = 1, ndf

          elemmat11(i,j) = &
              etaslip * sum ( phi1(:,i) * phi1(:,j) * surfl * Inn(:,1,1) * wg )
          elemmat11(i,j+ndf) = &
              etaslip * sum ( phi1(:,i) * phi1(:,j) * surfl * Inn(:,1,2) * wg )
          elemmat11(i,j+2*ndf) = &
              etaslip * sum ( phi1(:,i) * phi1(:,j) * surfl * Inn(:,1,3) * wg )
          elemmat11(i+ndf,j) = &
              etaslip * sum ( phi1(:,i) * phi1(:,j) * surfl * Inn(:,2,1) * wg )
          elemmat11(i+ndf,j+ndf) = &
              etaslip * sum ( phi1(:,i) * phi1(:,j) * surfl * Inn(:,2,2) * wg )
          elemmat11(i+ndf,j+2*ndf) = &
              etaslip * sum ( phi1(:,i) * phi1(:,j) * surfl * Inn(:,2,3) * wg )
          elemmat11(i+2*ndf,j) = &
              etaslip * sum ( phi1(:,i) * phi1(:,j) * surfl * Inn(:,3,1) * wg )
          elemmat11(i+2*ndf,j+ndf) = &
              etaslip * sum ( phi1(:,i) * phi1(:,j) * surfl * Inn(:,3,2) * wg )
          elemmat11(i+2*ndf,j+2*ndf) =  &
             etaslip * sum ( phi1(:,i) * phi1(:,j) * surfl * Inn(:,3,3) * wg )

          elemmat12(i,j) = &
              -etaslip * sum ( phi1(:,i) * phi2(:,j) * surfl * Inn(:,1,1) * wg )
          elemmat12(i,j+ndf) = &
              -etaslip * sum ( phi1(:,i) * phi2(:,j) * surfl * Inn(:,1,2) * wg )
          elemmat12(i,j+2*ndf) = &
              -etaslip * sum ( phi1(:,i) * phi2(:,j) * surfl * Inn(:,1,3) * wg )
          elemmat12(i+ndf,j) = &
              -etaslip * sum ( phi1(:,i) * phi2(:,j) * surfl * Inn(:,2,1) * wg )
          elemmat12(i+ndf,j+ndf) = &
              -etaslip * sum ( phi1(:,i) * phi2(:,j) * surfl * Inn(:,2,2) * wg )
          elemmat12(i+ndf,j+2*ndf) = &
              -etaslip * sum ( phi1(:,i) * phi2(:,j) * surfl * Inn(:,2,3) * wg )
          elemmat12(i+2*ndf,j) = &
              -etaslip * sum ( phi1(:,i) * phi2(:,j) * surfl * Inn(:,3,1) * wg )
          elemmat12(i+2*ndf,j+ndf) = &
              -etaslip * sum ( phi1(:,i) * phi2(:,j) * surfl * Inn(:,3,2) * wg )
          elemmat12(i+2*ndf,j+2*ndf) = &
              -etaslip * sum ( phi1(:,i) * phi2(:,j) * surfl * Inn(:,3,3) * wg )

          elemmat22(i,j) = &
              etaslip * sum ( phi2(:,i) * phi2(:,j) * surfl * Inn(:,1,1) * wg )
          elemmat22(i,j+ndf) = &
              etaslip * sum ( phi2(:,i) * phi2(:,j) * surfl * Inn(:,1,2) * wg )
          elemmat22(i,j+2*ndf) = &
              etaslip * sum ( phi2(:,i) * phi2(:,j) * surfl * Inn(:,1,3) * wg )
          elemmat22(i+ndf,j) = &
              etaslip * sum ( phi2(:,i) * phi2(:,j) * surfl * Inn(:,2,1) * wg )
          elemmat22(i+ndf,j+ndf) = &
              etaslip * sum ( phi2(:,i) * phi2(:,j) * surfl * Inn(:,2,2) * wg )
          elemmat22(i+ndf,j+2*ndf) = &
              etaslip * sum ( phi2(:,i) * phi2(:,j) * surfl * Inn(:,2,3) * wg )
          elemmat22(i+2*ndf,j) = &
              etaslip * sum ( phi2(:,i) * phi2(:,j) * surfl * Inn(:,3,1) * wg )
          elemmat22(i+2*ndf,j+ndf) = &
              etaslip * sum ( phi2(:,i) * phi2(:,j) * surfl * Inn(:,3,2) * wg )
          elemmat22(i+2*ndf,j+2*ndf) = &
              etaslip * sum ( phi2(:,i) * phi2(:,j) * surfl * Inn(:,3,3) * wg )

        end do
      end do

      elemmat21(1:3*ndf,1:3*ndf) = transpose(elemmat12(1:3*ndf,1:3*ndf))

    end if

  end subroutine elementcn

end module subs43_m
