!********************************************************************!
!********************************************************************!
!                                                                    !
!   rectilinear -- A code to solve NS equations for IBM              !
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

program rectilinear

  ! Declare modules

  use precision
  use mesh_module
  use flow_module
  use work_module

  implicit none

  ! Set up the mesh

  call mesh_constructor()

  ! Set up the flow

  call flow_constructor()

  ! ******************************* !
  ! START OF THE TIME STEPPING LOOP !
  ! ******************************* !
  
  do it = 1, 2000

     ! Apply the BCs

     call work_bcs()

     ! Get the primitives from the conserved

     call work_primitives()

     ! Compute the gradients of the primitives

     call work_gradients()

     ! Compute the intercell fluxes

     call work_fluxes()

     ! Compute any source terms

     call work_sources()

     ! Set the time steps

     call work_timestep()

     ! Update the solution

     call work_update()

     ! Print the timestep number

     write(6,*) it

  end do
  
  ! ******************************* !
  !  END OF THE TIME STEPPING LOOP  !
  ! ******************************* !


  ! Plot and have a look

!!$  do k = 2, mesh%nck-1
!!$  do j = 2, mesh%ncj-1
!!$  do i = 2, mesh%nci-1
  do k = 1, mesh%nck
  do j = 1, mesh%ncj
  do i = 1, mesh%nci
     write(19,*) mesh%x(1,i,j,k), mesh%x(2,i,j,k), mesh%x(3,i,j,k), flow%q(5,i,j,k), flow%qp(5,1,i,j,k)
  end do
  end do
  end do



end program rectilinear

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of program rectilinear                                       !
!                                                                    !
!********************************************************************!
!********************************************************************!

