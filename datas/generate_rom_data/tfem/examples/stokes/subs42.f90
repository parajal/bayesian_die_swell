module subs42_m

  use tfem_elem_m
  use math_defs_m

  implicit none

  real(dp) :: force(2), torque(1)

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

    integer :: ncurveconstr
    real(dp) :: xr(1,2)

    ncurveconstr = problem%constraints(constr)%geometry1

    xr(1,:) = mesh%coor(mesh%curves(ncurveconstr)%nodes(node),:)

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0
      if ( first ) then
        elemvecadd(1:2) = force
        elemvecadd(3) = torque(1)
      end if

    end if

!   connection through collocation

    elemmat(1,1) = 1
    elemmat(1,2) = 0

    elemmat(2,1) = 0
    elemmat(2,2) = 1

    elemmatadd(1,:) = [ -1._dp,  0._dp, -xr(1,2) ]
    elemmatadd(2,:) = [  0._dp, -1._dp,  xr(1,1) ]

  end subroutine elementc

! element for connection: slip on the particle boundary

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

    integer, parameter :: ndf = 3, nodalp = 3, ninti = 3, ndim = 2

    integer :: i, j, curve, coorsys
    real(dp) :: phi1(ninti,ndf), phi2(ninti,ndf), xr(ninti,1), x(nodalp,ndim), wg(ninti)
    real(dp) :: dphi(ninti,ndf,1), normal(ninti,ndim), xg(ninti,ndim)
    real(dp) :: dxdxis(ninti,ndim,1), curvel(ninti), Inn(ninti,ndim,ndim)
    real(dp) :: etaslip

    coorsys = coefficients%i(23)

!   first curve for connection

    curve = problem%connections(conn)%geometry1

!   Gauss points

    call Gauss_Legendre_line ( ninti, xr(:,1), wg )

!   shape function of the quadratic curve in the integration points

    call shape_line_P2 ( xr, phi1, dphi )
    call shape_line_P2 ( xr, phi2 )

!   compute geometry of element deformed

    call get_coordinates_geometry ( mesh, elem, x, curve=curve )

    call isoparametric_deformation_curved ( x, dphi, dxdxis, curvel, normal )

    call isoparametric_coordinates ( x, phi1, xg )

    if ( coorsys == 1 ) then
      curvel = 2 * pi * xg(:,2) * curvel
    end if

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

!     first velocity component only

      do i = 1, ndf
        do j = 1, ndf
          elemmat11(i,j) = etaslip * sum ( phi1(:,i) * phi1(:,j) * curvel * Inn(:,1,1) * wg )
          elemmat11(i,j+ndf) = etaslip * sum ( phi1(:,i) * phi1(:,j) * curvel * Inn(:,1,2) * wg )
          elemmat11(i+ndf,j) = etaslip * sum ( phi1(:,i) * phi1(:,j) * curvel * Inn(:,2,1) * wg )
          elemmat11(i+ndf,j+ndf) = etaslip * sum ( phi1(:,i) * phi1(:,j) * curvel * Inn(:,2,2) * wg )

          elemmat12(i,j) = -etaslip * sum ( phi1(:,i) * phi2(:,j) * curvel * Inn(:,1,1) * wg )
          elemmat12(i,j+ndf) = -etaslip * sum ( phi1(:,i) * phi2(:,j) * curvel * Inn(:,1,2) * wg )
          elemmat12(i+ndf,j) = -etaslip * sum ( phi1(:,i) * phi2(:,j) * curvel * Inn(:,2,1) * wg )
          elemmat12(i+ndf,j+ndf) = -etaslip * sum ( phi1(:,i) * phi2(:,j) * curvel * Inn(:,2,2) * wg )

          elemmat22(i,j) = etaslip * sum ( phi2(:,i) * phi2(:,j) * curvel * Inn(:,1,1) * wg )
          elemmat22(i,j+ndf) = etaslip * sum ( phi2(:,i) * phi2(:,j) * curvel * Inn(:,1,2) * wg )
          elemmat22(i+ndf,j) = etaslip * sum ( phi2(:,i) * phi2(:,j) * curvel * Inn(:,2,1) * wg )
          elemmat22(i+ndf,j+ndf) = etaslip * sum ( phi2(:,i) * phi2(:,j) * curvel * Inn(:,2,2) * wg )
        end do
      end do

      elemmat21(1:2*ndf,1:2*ndf) = transpose(elemmat12(1:2*ndf,1:2*ndf))

    end if

  end subroutine elementcn

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

    integer, parameter :: ndf = 3, ndfl = 3, nodalp = 3, ninti = 3, ndim = 2

    integer :: i, j, curve, coorsys
    real(dp) :: phi(ninti,ndf), psi(ninti,ndfl), xr(ninti,1), x(nodalp,ndim), wg(ninti)
    real(dp) :: dphi(ninti,ndf,1), normal(ninti,ndim), xg(ninti,ndim)
    real(dp) :: dxdxis(ninti,ndim,1), curvel(ninti)
    real(dp) :: tmp_x(ndfl,ndf), tmp_y(ndfl,ndf)

    coorsys = coefficients%i(23)

!   first curve for connection

    curve = problem%constraints(constr)%geometry1

!   Gauss points

    call Gauss_Legendre_line ( ninti, xr(:,1), wg )

!   shape function of the quadratic curve in the integration points

    call shape_line_P2 ( xr, phi, dphi )
    call shape_line_P2 ( xr, psi )

!   compute geometry of element deformed

    call get_coordinates_geometry ( mesh, elem, x, curve=curve )

    call isoparametric_deformation_curved ( x, dphi, dxdxis, curvel, normal )

    call isoparametric_coordinates ( x, phi, xg )

    if ( coorsys == 1 ) then
      curvel = 2 * pi * xg(:,2) * curvel
    end if

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0

    end if

    if ( matrix ) then

      do i = 1, ndfl
        do j = 1, ndf
          tmp_x(i,j) = sum ( psi(:,i) * phi(:,j) * curvel * normal(:,1)* wg )
          tmp_y(i,j) = sum ( psi(:,i) * phi(:,j) * curvel * normal(:,2)* wg )
        end do
      end do

      elemmat = 0
      elemmat(1:ndfl,1:ndf) = tmp_x
      elemmat(1:ndfl,ndf+1:2*ndf) = tmp_y

      elemmat2 = - elemmat

    end if

  end subroutine stokes_constr_normal_vel

end module subs42_m
