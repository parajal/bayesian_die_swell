
! Copyright (C) 2004-2012 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines routines for input/output of mesh, problem, ...

module io_utils_m

  use tfem_utils_m  ! I/O for core tfem structures etc.
  use sep_utils_m   ! I/O utils for sepran
  use gbt_utils_m   ! I/O utils for Gambit (Neutral file format)
  use gmsh_utils_m  ! I/O utils for Gmsh (msh file format)
  use ideas_utils_m ! I/O utils for NX-Ideas (Universal file format)
  use vtk_utils_m   ! I/O utils for VTK (legacy vtk file format)
  use tec_utils_m   ! I/O utils for Tecplot
  use printtofile_m ! Printing of data in vector/sysvectors to a file
  use printinfo_m   ! Printing of info from mesh, problem,... to standard output

  implicit none

end module io_utils_m
