
! module with additional functions for the extrudate swell problem with xfem

module functions_xfem_m

  use tfem_elem_m

  implicit none
  save

! local pointers to data in main program
  type(mesh_t), pointer :: lmesh_sf_adv => null()
  type(problem_t), pointer :: lproblem_sf_adv => null()
  type(sysvector_t), pointer :: lsol_sf_adv_pred => null()


! initial height
  real(dp) :: h0

! mesh used for mapping of reference element
  type(mesh_t), pointer :: lmesh => null()

  integer :: lelem, lelgrp

contains

! levelset function for the interface

  function levelset ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1)) :: levelset

    integer :: i, elem, nod(2)

    real(dp) :: x0, x1, x2, y0, y1, y2, yp, xtmp(2), u(2)
    logical :: find


    do i = 1, size(x,1)

      x0 = x(i,1)
      y0 = x(i,2)

      if ( x0 <= 0._dp ) then

!       inside the die

        levelset(i) = h0 - y0

      else

!       outside the die

        find  = .false.

        do elem = 1, lmesh_sf_adv%nelem

          nod = lmesh_sf_adv%topology(1)%a(1:2,elem)
          xtmp = lmesh_sf_adv%coor(nod,1)

          x1 = xtmp(1); x2 = xtmp(2)

          if ( (x0 >= x1) .and. (x0 <= x2) ) then

!           get surface height for element elem

            u = lsol_sf_adv_pred%u( &
              lproblem_sf_adv%degfdperm(lproblem_sf_adv%nodnumdegfd(nod)+1,2) )
            y1 = u(1); y2 = u(2)

!           height at position x0
            yp = y1 + (y2-y1)*(x0-x1)/(x2-x1)

            levelset(i) = yp - y0

            find = .true.

            exit

          end if

        end do

        if ( .not. find ) then
          print*, 'Error: Finding levelset failed.'
          print*, i, x0, y0
          print*, 'Program Stop!'
          stop
        end if

      end if

    end do


  end function levelset


! function for mapping reference coordinates to the real coordinates
! in building the eltree within an element.

  function mapcoor ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1),size(x,2)) :: mapcoor

    real(dp) :: phi(size(x,1),9), xnod(9,2)

    call shape_quad_Q2 ( x, phi )

    call get_coordinates ( lmesh, lelgrp, lelem, xnod )

    mapcoor = matmul ( phi, xnod ) ! isoparametric

  end function mapcoor


! function for the incremental vertical displacement of the interface
! Used for computing the temporary ALE mesh (md=mesh displacement).

  function func_md ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func_md

    real(dp) :: disp(1), coor(1,size(x))

    if ( x(1) <= 0 ) then
      func_md = 0._dp
    else
      coor(1,:) = x
      disp = levelset( coor )
      func_md = disp(1)
    end if

  end function func_md


end module functions_xfem_m
