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
     
     !Array to store wall distance
     
     real (kind=RP), allocatable :: d(:,:,:,:)

     ! Array to store interpolation technique

     real (kind=RP), allocatable :: inter(:,:,:,:)
     
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
    allocate(flow%ferror(npdes, mesh%nci, mesh%ncj, mesh%nck))
    allocate(flow%force(npdes, mesh%nci, mesh%ncj, mesh%nck))

    ! wall distance calculation

    allocate(flow%d(2, mesh%nci,mesh%ncj,mesh%nck))
        
    ! Interpolation techniques
    allocate (flow%inter(1, mesh%nci,mesh%ncj,mesh%nck))
    
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

    ! Wall distance variables
    
    real   (kind=RP) :: u,c,q,r

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

       ! Wall distance initial condition

       u = 0d0
       c = 0d0
       q = 0d0
       r = 0d0

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

       ! pde numbers 6 to 9 relate to wall distance calculation

       flow%u(6,i,j,k) = u
       flow%u(7,i,j,k) = c
       flow%u(8,i,j,k) = q
       flow%u(9,i,j,k) = r

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

       rc =1.0d0*atan(1d0)

      ! write(6,*) mesh%x(1,i,j,k), mesh%x(2,i,j,k), mesh%x(3,i,j,k), r2, rc

      if ( r2 <= rc**2) then
          
          flow%ibm(1,i,j,k) = 1.0d0
       
      else
      
          flow%ibm(1,i,j,k) = 0.0d0
       
      end if
      
     ! write(6,*) mesh%x(1,i,j,k), mesh%x(2,i,j,k),mesh%x(3,i,j,k), r2,rc**2,flow%ibm(1,i,j,k)
      
     end do
     end do
     end do

end subroutine ibm_constructor

!*******************************************************************
!*******************************************************************
     
    subroutine interpolation_constructor()
      
      ! Declare variables
      
      real(kind=RP) :: delta, r2, rc

      do k = 1, mesh%nck
      do j = 1, mesh%ncj
      do i = 1, mesh%nci
      
      ! Set a value for delta on to the IBM constructor to attain the outer value of 

      delta = 0.005d0

      ! Get the coordinates

      r2 = mesh%x(1,i,j,k)**2 + mesh%x(2,i,j,k)**2 + mesh%x(3,i,j,k)**2

      ! Test for the critical radius 1/4 PI + Delta

      rc = 0.5d0*atan(1d0) + delta

      if ( r2 <= rc**2) then

      flow%inter(1,i,j,k) = 1.0d0

      else

      flow%inter(1,i,j,k) = 0.0d0
      
      end if
 
      end do
      end do
      end do
      
      end subroutine interpolation_constructor
      
!*******************************************************************
!*******************************************************************


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
   alpha = -200d0
   beta = -1.5d0   

   !apply continuous forcing approach 

   !Calculate the force integral error

   flow%ferror(2,i,j,k) = flow%ferror(2,i,j,k) + (flow%u(2,i,j,k) / flow%u(1,i,j,k) - vel)
   flow%ferror(3,i,j,k) = flow%ferror(3,i,j,k) + (flow%u(3,i,j,k) / flow%u(1,i,j,k) - vel)
   flow%ferror(4,i,j,k) = flow%ferror(4,i,j,k) + (flow%u(4,i,j,k) / flow%u(1,i,j,k) - vel)
   
   ! Wall distance forcing error
   flow%ferror(6,i,j,k) = flow%ferror(6,i,j,k) + (flow%u(6,i,j,k) - 0d0)

   !Calculate the force to be applied     

   flow%force(2,i,j,k) = (alpha * flow%ferror(2,i,j,k) + beta * (flow%u(2,i,j,k) / flow%u(1,i,j,k) - vel)) * flow%ibm(1,i,j,k)
   flow%force(3,i,j,k) = (alpha * flow%ferror(3,i,j,k) + beta * (flow%u(3,i,j,k) / flow%u(1,i,j,k) - vel)) * flow%ibm(1,i,j,k)
   flow%force(4,i,j,k) = (alpha * flow%ferror(4,i,j,k) + beta * (flow%u(4,i,j,k) / flow%u(1,i,j,k) - vel)) * flow%ibm(1,i,j,k)
   
   ! Wall distance boundary force

   ! Set u to = 0 by applying a feedback force

   flow%force(6,i,j,k) = (alpha * flow%ferror(6,i,j,k) + beta * (flow%u(6,i,j,k) / flow%u(1,i,j,k) - 0)) *flow%ibm(1,i,j,k)
     
   end do
   end do
   end do

end subroutine force_constructor

 !*******************************************************************
 !*******************************************************************

  subroutine work_walldist
   
    ! Declare Variables
    
    real(kind=rp) :: Y
    
    do k= 1, mesh%nck
    do j= 1, mesh%ncj
    do i= 1, mesh%nci
  
       Y =  sqrt(flow%u(7,i,j,k)**2 + flow%u(8,i,j,k)**2 + flow%u(9,i,j,k)**2 + 2d0*flow%u(6,i,j,k))
     
       flow%d(1,i,j,k) = -sqrt(flow%u(7,i,j,k)**2 + flow%u(8,i,j,k)**2 + flow%u(9,i,j,k)**2) + Y
    
       flow%d(2,i,j,k) = -sqrt(flow%u(7,i,j,k)**2 + flow%u(8,i,j,k)**2 + flow%u(9,i,j,k)**2) - Y

    end do
    end do
    end do

  end subroutine work_walldist

  end module flow_module


!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module flow_module                                        !
!                                                                    !
!********************************************************************!
!********************************************************************!

