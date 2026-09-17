
! Copyright (C) 2020-2020 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Compute material parameters from the input.
! Only isotropic materials are currently available.

! Here (lambda,mu) are the Lame parameters:
!
!               nu * E                   E
!  lambda = --------------,  mu = G = -------
!           (1+nu)(1-2*nu)            2(1+nu)
!
! where E is Young's modulus and nu Poisson's ratio. Note, G is the shear
! modulus. Also computed is the bulk modulus K, given by
!
!           E
!   K = ---------
!       3(1-2*nu)
!
! For the plane-stress case (2D, Cartesian) lambda and 3*K need to be
! replaced by lambda_bar and 2*K_bar with
!
!                nu * E                 E
!  lambda_bar = ---------,   K_bar = -------
!               (1-nu**2)            2(1-nu)
!
! For thermal expansion alpha is the linear coefficient of thermal expansion.

module linear_elastic_material_m

  use tfem_elem_m

  implicit none

  type lemodel_t
    logical :: plane_stress = .false.
    logical :: incompressible = .false.
    logical :: thermal_expansion = .false.
    real(dp) :: Emod, nu, lambda, mu, lambda_bar, Kmod, Kmod_bar, alpha
  end type lemodel_t

contains

! Compute linear elastic material parameters

  subroutine set_linear_elastic_material ( coefficients, coorsys, lemodel )

    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: coorsys
!
!   output of material parameters.
    type(lemodel_t), intent(out) :: lemodel

    real(dp) :: Emod, nu, alpha

    Emod = coefficients%r(1)
    nu = coefficients%r(3)
    alpha = coefficients%r(5)

    lemodel%Emod = Emod
    lemodel%nu = nu
    lemodel%mu = Emod / ( 2*(1+nu) )

    if ( coorsys == 0 .and. coefficients%i(67) == 0 .and. &
         coefficients%i(3) == 1 ) then

!     plane stress (lambda and Kmod not set!)

      lemodel%plane_stress = .true.
      lemodel%incompressible = .false.
      lemodel%lambda_bar = nu * Emod / (1-nu**2)
      lemodel%Kmod_bar = Emod / ( 2*(1-nu) )

    else

!     plane strain (2D,Cartesion), axisymmetric of 3D (lambda_bar not set!)

      lemodel%plane_stress = .false.
      if ( coefficients%i(4) == 1 ) then
!       incompressible material (lambda, Kmod not set)
        lemodel%incompressible = .true.
        lemodel%nu = 0.5_dp
        lemodel%mu = Emod / 3
      else
        lemodel%incompressible = .false.
        lemodel%lambda = nu * Emod / ( (1+nu)*(1-2*nu) )
        lemodel%Kmod = Emod / ( 3*(1-2*nu) )
      end if

    end if

    if ( coefficients%i(5) == 1 ) then

      lemodel%thermal_expansion = .true.
      lemodel%alpha = alpha

    else

      lemodel%thermal_expansion = .false.

    end if

  end subroutine set_linear_elastic_material

end module linear_elastic_material_m

