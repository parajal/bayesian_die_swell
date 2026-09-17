
! Copyright (C) 2007-2011 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the convection-diffusion-reaction equation:
!
!          dc
!  gamma ( -- + u.nabla c ) - nabla ( alpha nabla c ) + beta c = f
!          dt
!
! where u is a velocity vector and alpha, beta and gamma are coefficients
! that may depend on position.
! The coefficient alpha can be either a scalar or a tensor.

module convection_diffusion_elements_m

  use tfem_elem_m
  use poisson_elements_m
  use convection_diffusion_elements_generic_m
  use convection_diffusion_elements_embedded_boundary_m
  use convection_diffusion_elements_embedded_interface_m

  implicit none

end module convection_diffusion_elements_m

