!********************************************************************!
!********************************************************************!
!                                                                    !
!   work_module -- A module for doing generic PDE calcs              !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Version history:                                                 !
!                    Program created: 04Oct19                - raw54 !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Known issues:                                                    !
!                    No known issues.                        - raw54 !
!                                                                    !
!********************************************************************!
!********************************************************************!

module work_module

  ! Declare modules

  use precision
  use mesh_module
  use flow_module
  use equations_module

  implicit none

contains

  !*******************************************************************
  !*******************************************************************

  subroutine work_primitives()

    ! Loop over all cells accumulating, gathering primitives

    do k = 1, mesh%nck
    do j = 1, mesh%ncj
    do i = 1, mesh%nci

       flow%q(:,i,j,k) = u2q(flow%u(:,i,j,k))

    end do
    end do
    end do

    ! Return to calling program

    return

  end subroutine work_primitives

  !*******************************************************************
  !*******************************************************************

  subroutine work_gradients()

    ! Declare variables

    real(kind=RP) :: qx(npdes), qy(npdes), qz(npdes)

    ! Zero the gradients here
    
    flow%qp = 0d0

    ! Loop over all cells accumulating the gradient fluxes

    do k = 1, mesh%nck-1
    do j = 1, mesh%ncj-1
    do i = 1, mesh%nci-1

       ! Get the gradient fluxes for this equation

       qx = 0.5d0 * ( flow%q(:,i,j,k) + flow%q(:,i+1,j,k) ) * mesh%nor(1,i,j,k)
       qy = 0.5d0 * ( flow%q(:,i,j,k) + flow%q(:,i,j+1,k) ) * mesh%nor(2,i,j,k)
       qz = 0.5d0 * ( flow%q(:,i,j,k) + flow%q(:,i,j,k+1) ) * mesh%nor(3,i,j,k)

       ! And accumulate them into the residual
       
       flow%qp(:,1,i  ,j,k) = flow%qp(:,1,i  ,j,k) + qx(:)
       flow%qp(:,1,i+1,j,k) = flow%qp(:,1,i+1,j,k) - qx(:)

       flow%qp(:,2,i,j  ,k) = flow%qp(:,2,i,j  ,k) + qy(:)
       flow%qp(:,2,i,j+1,k) = flow%qp(:,2,i,j+1,k) - qy(:)

       flow%qp(:,3,i,j,k  ) = flow%qp(:,3,i,j,k  ) + qz(:)
       flow%qp(:,3,i,j,k+1) = flow%qp(:,3,i,j,k+1) - qz(:)

    end do
    end do
    end do

    ! Correct the qp calculation at the boundaries
    
    do k = 1, mesh%nck
    do j = 1, mesh%ncj
       flow%qp(:,1,       1,j,k) = flow%qp(:,1,         2,j,k)
       flow%qp(:,1,mesh%nci,j,k) = flow%qp(:,1,mesh%nci-1,j,k)
       flow%qp(:,2,       1,j,k) = flow%qp(:,2,         2,j,k)
       flow%qp(:,2,mesh%nci,j,k) = flow%qp(:,2,mesh%nci-1,j,k)
       flow%qp(:,3,       1,j,k) = flow%qp(:,3,         2,j,k)
       flow%qp(:,3,mesh%nci,j,k) = flow%qp(:,3,mesh%nci-1,j,k)
    end do
    end do
    do k = 1, mesh%nck
    do i = 1, mesh%nci
       flow%qp(:,1,i,       1,k) = flow%qp(:,1,i,         2,k)
       flow%qp(:,1,i,mesh%ncj,k) = flow%qp(:,1,i,mesh%ncj-1,k)
       flow%qp(:,2,i,       1,k) = flow%qp(:,2,i,         2,k)
       flow%qp(:,2,i,mesh%ncj,k) = flow%qp(:,2,i,mesh%ncj-1,k)
       flow%qp(:,3,i,       1,k) = flow%qp(:,3,i,         2,k)
       flow%qp(:,3,i,mesh%ncj,k) = flow%qp(:,3,i,mesh%ncj-1,k)
    end do
    end do
    do j = 1, mesh%ncj
    do i = 1, mesh%nci
       flow%qp(:,1,i,j,       1) = flow%qp(:,1,i,j,         2)
       flow%qp(:,1,i,j,mesh%nck) = flow%qp(:,1,i,j,mesh%nck-1)
       flow%qp(:,2,i,j,       1) = flow%qp(:,2,i,j,         2)
       flow%qp(:,2,i,j,mesh%nck) = flow%qp(:,2,i,j,mesh%nck-1)
       flow%qp(:,3,i,j,       1) = flow%qp(:,3,         i,j,2)
       flow%qp(:,3,i,j,mesh%nck) = flow%qp(:,3,i,j,mesh%nck-1)
    end do
    end do

    ! Divide all the qps by the volume

    do k = 1, mesh%nck
    do j = 1, mesh%ncj
    do i = 1, mesh%nci
    
       flow%qp(:,:,i,j,k) = flow%qp(:,:,i,j,k) / mesh%vol(1,i,j,k)

    end do
    end do
    end do

    ! Return to calling program

    return

  end subroutine work_gradients


  !*******************************************************************
  !*******************************************************************

  subroutine work_fluxes()

    ! Declare variables

    real(kind=RP) :: f(npdes), g(npdes), h(npdes)

    ! Zero the residuals here
    
    flow%res = 0d0

    ! Loop over all cells accumulating the fluxes

    do k = 1, mesh%nck-1
    do j = 1, mesh%ncj-1
    do i = 1, mesh%nci-1

       ! Get the fluxes for this equation

       f = xflux(flow%u(:,i  ,j,k), flow%q(:,i  ,j,k), flow%qp(:,:,i  ,j,k), &
                 flow%u(:,i+1,j,k), flow%q(:,i+1,j,k), flow%qp(:,:,i+1,j,k), &
                 mesh%nor(:,i,j,k))
       g = yflux(flow%u(:,i,j  ,k), flow%q(:,i,j  ,k), flow%qp(:,:,i,j  ,k), &
                 flow%u(:,i,j+1,k), flow%q(:,i,j+1,k), flow%qp(:,:,i,j+1,k), &
                 mesh%nor(:,i,j,k))
       h = zflux(flow%u(:,i,j,k  ), flow%q(:,i,j,k  ), flow%qp(:,:,i,j,k  ), &
                 flow%u(:,i,j,k+1), flow%q(:,i,j,k+1), flow%qp(:,:,i,j,k+1), &
                 mesh%nor(:,i,j,k))

       ! And accumulate them into the residual
       
       flow%res(:,i  ,j  ,k  ) = flow%res(:,i  ,j  ,k  ) - f(:)
       flow%res(:,i+1,j  ,k  ) = flow%res(:,i+1,j  ,k  ) + f(:)

       flow%res(:,i  ,j  ,k  ) = flow%res(:,i  ,j  ,k  ) - g(:)
       flow%res(:,i  ,j+1,k  ) = flow%res(:,i  ,j+1,k  ) + g(:)

       flow%res(:,i  ,j  ,k  ) = flow%res(:,i  ,j  ,k  ) - h(:)
       flow%res(:,i  ,j  ,k+1) = flow%res(:,i  ,j  ,k+1) + h(:)

    end do
    end do
    end do

    ! Return to calling program

    return

  end subroutine work_fluxes

  !*******************************************************************
  !*******************************************************************

  subroutine work_timestep()

    ! Find the maximum stable timestep of each cell

    do k = 1, mesh%nck
    do j = 1, mesh%ncj
    do i = 1, mesh%nci

       mesh%dt(:,i,j,k) = 0.25d0*mesh%dx(1,i,j,k) / lambda(flow%u(:,i,j,k))

    end do
    end do
    end do

    ! Return to calling program

    return

  end subroutine work_timestep

  !*******************************************************************
  !*******************************************************************

  subroutine work_update()

    ! Update the conserved variables by the residuals

    do k = 1, mesh%nck
    do j = 1, mesh%ncj
    do i = 1, mesh%nci

!!       flow%u(:,i,j,k) = flow%u(:,i,j,k) + mesh%dt(1,i,j,k) * flow%res(:,i,j,k) / mesh%vol(1,i,j,k)
       flow%u(:,i,j,k) = flow%u(:,i,j,k) + 0.00005d0 * flow%res(:,i,j,k) / mesh%vol(1,i,j,k)

    end do
    end do
    end do

    ! Return to calling program

    return

  end subroutine work_update

  !*******************************************************************
  !*******************************************************************

  subroutine work_bcs()

    ! Inflow and outflow boundary (x_min = const, x_max = const)

    do k = 1, mesh%nck
    do j = 1, mesh%ncj

       flow%u(:,       1,j,k) = bc_type_a(flow%u(:,         2,j,k))
       flow%u(:,mesh%nci,j,k) = bc_type_b(flow%u(:,mesh%nci-1,j,k))

    end do
    end do

    ! Zero gradient boundaries (y_min = const, y_max = const)

    do k = 1, mesh%nck
    do i = 1, mesh%nci

       flow%u(:,i,       1,k) = bc_type_c(flow%u(:,i,         2,k),(/ 0d0, -1d0, 0d0 /))
       flow%u(:,i,mesh%ncj,k) = bc_type_c(flow%u(:,i,mesh%ncj-1,k),(/ 0d0,  1d0, 0d0 /))

    end do
    end do

    ! Zero gradient boundaries (z_min = const, z_max = const)

    do j = 1, mesh%ncj
    do i = 1, mesh%nci

       flow%u(:,i,j,       1) = bc_type_c(flow%u(:,i,j,         2),(/ 0d0, 0d0, -1d0 /))
       flow%u(:,i,j,mesh%nck) = bc_type_c(flow%u(:,i,j,mesh%nck-1),(/ 0d0, 0d0,  1d0 /))

    end do
    end do

    ! Return to calling program

    return

  end subroutine work_bcs

  !*******************************************************************
  !*******************************************************************

  subroutine work_sources()

    implicit none

    ! Loop over all terms gathering the sources

    do k = 1, mesh%nck
    do j = 1, mesh%ncj
    do i = 1, mesh%nci

       flow%res(:,i,j,k) = flow%res(:,i,j,k) + sources(flow%q(:,i,j,k), mesh%vol(:,i,j,k))

    end do
    end do
    end do

    ! Return to calling program

    return

  end subroutine work_sources


end module work_module

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module work_module                                        !
!                                                                    !
!********************************************************************!
!********************************************************************!

