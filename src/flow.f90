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

     !Array to store a switch for inside or outside IBM field

     real (kind=RP), allocatable :: ibm(:,:,:,:)

     !Array to store the force error used to calc the force

     real (kind=RP), allocatable :: ferror(:,:,:,:)
     
     !Array to store the force to be applied
     
     real (kind=RP), allocatable :: force(:,:,:,:)


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

    ! ibm and forcing allocation

    allocate(flow%ibm(1, mesh%nci, mesh%ncj, mesh%nck))
    allocate(flow%ferror(5, mesh%nci, mesh%ncj, mesh%nck))
    allocate(flow%force(5, mesh%nci, mesh%ncj, mesh%nck))

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

       ! initial forcing conditions

       flow%force(1,i,j,k) = 0
       flow%force(2,i,j,k) = 0
       flow%force(3,i,j,k) = 0
       flow%force(4,i,j,k) = 0
       flow%force(5,i,j,k) = 0

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

!********************************************************************!
!********************************************************************!

  subroutine ibm_constructor()
    
    ! Declare variables

    real(kind=Rp) :: r2,rc

    ! set the value of the ibm array to 0 outside and 1 inside

    do k = 1, mesh%nck
    do j = 1, mesh%ncj
    do i = 1, mesh%nci

       ! Get the coordinates
       
       r2 = mesh%x(1,i,j,k)**2 + mesh%x(2,i,j,k)**2 + mesh%x(3,i,j,k)**2
       
       ! Test for inside critical radius 1/4 PI

       rc =1d0*atan(1d0)

      ! write(6,*) mesh%x(1,i,j,k), mesh%x(2,i,j,k), mesh%x(3,i,j,k), r2, rc

       if ( r2 <= rc**2) then
          
          flow%ibm(1,i,j,k) = 1.0
       else
          flow%ibm(1,i,j,k) = 0.0
       end if
      
      ! write(6,*) mesh%x(1,i,j,k), mesh%x(2,i,j,k),mesh%x(3,i,j,k), r2,rc**2,flow%ibm(1,i,j,k)
       
     end do
     end do
     end do

end subroutine ibm_constructor

   subroutine force_constructor

   !Declare Variables
     
     Implicit none

   real (Kind=RP) :: alpha,beta,vel

   
   !calculate the force in all axis

   do k= 1, mesh%nck
   do j= 1, mesh%ncj
   do i= 1, mesh%nci

   !Set values for vel, alpha and beta
   !vel is the desired velocity on the boundary. vel=0 for the no slip condition
   !Current values of Alpha and Beta taken from Goldstein for turbulent flow condition (inital estimate)
 
   vel = 0d0
   alpha =-200d0
   beta =-1.5d0   

   !apply continuous forcing approach 

   !Calculate the force integral error

   flow%ferror(2,i,j,k) = flow%ferror(2,i,j,k) + (flow%u(2,i,j,k) / flow%u(1,i,j,k) - vel)
   flow%ferror(3,i,j,k) = flow%ferror(3,i,j,k) + (flow%u(3,i,j,k) / flow%u(1,i,j,k) - vel)
   flow%ferror(4,i,j,k) = flow%ferror(4,i,j,k) + (flow%u(4,i,j,k) / flow%u(1,i,j,k) - vel)

  !Calculate the force to be applied     
  !ensure force is applied only inside the body (multiply by ibm)

   flow%force(2,i,j,k) = (alpha * flow%ferror(2,i,j,k) + beta * (flow%u(2,i,j,k) / flow%u(1,i,j,k) - vel)) * flow%ibm(1,i,j,k)
   flow%force(3,i,j,k) = (alpha * flow%ferror(3,i,j,k) + beta * (flow%u(3,i,j,k) / flow%u(1,i,j,k) - vel)) * flow%ibm(1,i,j,k)
   flow%force(4,i,j,k) = (alpha * flow%ferror(4,i,j,k) + beta * (flow%u(4,i,j,k) / flow%u(1,i,j,k) - vel)) * flow%ibm(1,i,j,k)
       
   end do
   end do
   end do

end subroutine force_constructor

end module flow_module

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module flow_module                                        !
!                                                                    !
!********************************************************************!
!********************************************************************!

