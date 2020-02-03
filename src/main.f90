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
  use io_module

  implicit none

  ! Set up the mesh

  call mesh_constructor()

  ! Set up the flow

  call flow_constructor()

  ! Set up the IBM array

  call ibm_constructor()

  ! ******************************* !
  ! START OF THE TIME STEPPING LOOP !
  ! ******************************* !
  
  do it = 1, 10000

     ! Apply the BCs

     call work_bcs()

     ! Get the primitives from the conserved

     call work_primitives()

     ! Compute the gradients of the primitives

     call work_gradients()

     ! Compute the intercell fluxes

     call work_fluxes()

     ! Call IBM force constructor 
     
     call force_constructor

     ! Compute any source terms

     call work_sources()

     ! Compute the nearest wall distance d
     
     call work_walldist()
      
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

  call write_plot3d()
  
  do k=1, mesh%nck
  do j=1, mesh%ncj
  do i=1, mesh%nci

     !write(19,*) mesh%x(1,i,j,k), mesh%x(2,i,j,k), mesh%x(3,i,j,k), flow%inter(1,i,j,k), flow%ibm(1,i,j,k) 
     write(19, *) flow%inter(1,i,j,k), flow%ibm(1,i,j,k)

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

