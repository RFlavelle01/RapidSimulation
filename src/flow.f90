!********************************************************************!
!********************************************************************!
!                                                                    !
!   flow_module -- A module for storing flow data                    !
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

module flow_module

  ! Declare modules

  use precision
  use mesh_module
  use equations_module

  implicit none

  ! Set up counters for the run time

  integer(kind=IP) :: it

  ! Set up a type to contain the flow data

  type flow_type

     ! Array to store the conserved and primitive variables

     real   (kind=RP), allocatable :: u(:,:,:,:), q(:,:,:,:)

     ! Array to store the residuals

     real   (kind=RP), allocatable :: res(:,:,:,:)

     ! Array to store the gradients of the primitives

     real   (kind=RP), allocatable :: qp(:,:,:,:,:)

  end type flow_type

  ! And declare a variable flow

  type(flow_type) :: flow

contains

  !*******************************************************************
  !*******************************************************************

  subroutine flow_constructor()

    ! Allocate memory for the various flow variables based on the mesh

    allocate(flow%u(npdes, mesh%nci, mesh%ncj, mesh%nck))
    allocate(flow%q(npdes, mesh%nci, mesh%ncj, mesh%nck))

    allocate(flow%res(npdes,    mesh%nci, mesh%ncj, mesh%nck))
    allocate(flow%qp( npdes, 3, mesh%nci, mesh%ncj, mesh%nck))

    ! Set up the initial solution in u

    call initial_solution()

    ! Return to calling program

    return

  end subroutine flow_constructor

  !*******************************************************************
  !*******************************************************************

  subroutine initial_solution()

    ! Declare variables

    real   (kind=RP) :: t, d, p, vx, vy, vz

    ! Set up the initial solution

    do k = 1, mesh%nck
    do j = 1, mesh%ncj
    do i = 1, mesh%nci

       t = 280d0

       d = 1.266d0
       p = d*rgas*t

       vx = 0d0
       vy = 0d0
       vz = 0d0

       flow%u(1,i,j,k) = d
       flow%u(2,i,j,k) = d * vx
       flow%u(3,i,j,k) = d * vy
       flow%u(4,i,j,k) = d * vz
       flow%u(5,i,j,k) = p*gmi + 0.5d0 * d * (vx*vx + vy*vy + vz*vz)

    end do
    end do
    end do

    ! Return to calling program

    return

  end subroutine initial_solution


end module flow_module

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module flow_module                                        !
!                                                                    !
!********************************************************************!
!********************************************************************!

