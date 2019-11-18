!********************************************************************!
!********************************************************************!
!                                                                    !
!   io_module -- Disk I/O Routines                                   !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Version history:                                                 !
!                    Program created: 08Nov19                - raw54 !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Known issues:                                                    !
!                    No known issues.                        - raw54 !
!                                                                    !
!********************************************************************!
!********************************************************************!

module io_module

  ! Declare modules

  use precision
  use mesh_module
  use flow_module

  implicit none

  ! Declare variables

contains

  !*******************************************************************
  !*******************************************************************

  subroutine write_ascii()

    implicit none

    ! Write out data

    do k = 2, mesh%nck-1
    do j = 2, mesh%ncj-1
    do i = 2, mesh%nci-1
       write(19,*) mesh%x(1,i,j,k), mesh%x(2,i,j,k), mesh%x(3,i,j,k), flow%q(1,i,j,k), flow%qp(2,1,i,j,k)
    end do
    end do
    end do

    ! Return to calling program

    return

  end subroutine write_ascii

  !*******************************************************************
  !*******************************************************************

  subroutine write_plot3d()

    implicit none

    ! Declare local variables

    integer(kind=IP) :: ipde

    ! Write the mesh

    open(unit=13, file="flow.xyz", form="unformatted")

    write(13) 1
    write(13) mesh%nci-2, mesh%ncj-2, mesh%nck-2

    write(13) (((mesh%x(1,i,j,k), i=2, mesh%nci-1), j=2, mesh%ncj-1), k=2, mesh%nck-1), &
              (((mesh%x(2,i,j,k), i=2, mesh%nci-1), j=2, mesh%ncj-1), k=2, mesh%nck-1), &
              (((mesh%x(3,i,j,k), i=2, mesh%nci-1), j=2, mesh%ncj-1), k=2, mesh%nck-1)
    close(13)

    ! Write the flow

    open(unit=14, file="flow.q", form="unformatted")

    write(14) 1
    write(14) mesh%nci-2, mesh%ncj-2, mesh%nck-2
    write(14) 0d0, 0d0, 0d0, 0d0
    write(14) ((((flow%u(ipde,i,j,k), i=2, mesh%nci-1), j=2, mesh%ncj-1), k=2, mesh%nck-1), ipde=1,5)
    
    ! Close flow file
    
    close(14)

    ! Return to calling program

    return

  end subroutine write_plot3d

  !*******************************************************************
  !*******************************************************************

  subroutine read_restart()

    implicit none

    ! Declare local variables

    integer(kind=IP) :: ipde

    ! No need to read the mesh, since that is taken from mesh.dat

    ! Read the flow

    open(unit=14, file="flow.dat", form="unformatted")

    read(14) ((((flow%u(ipde,i,j,k), i=1, mesh%nci), j=1, mesh%ncj), k=1, mesh%nck), ipde=1,5)
    
    ! Close flow file
    
    close(14)

    ! Return to calling program

    return

  end subroutine read_restart


  !*******************************************************************
  !*******************************************************************

  subroutine write_restart()

    implicit none

    ! Declare local variables

    integer(kind=IP) :: ipde

    ! No need to write the mesh, since that is taken from mesh.dat

    ! Write the flow

    open(unit=14, file="flow.dat", form="unformatted")

    write(14) ((((flow%u(ipde,i,j,k), i=1, mesh%nci), j=1, mesh%ncj), k=1, mesh%nck), ipde=1,5)
    
    ! Close flow file
    
    close(14)

    ! Return to calling program

    return

  end subroutine write_restart

end module io_module

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module io_module                                          !
!                                                                    !
!********************************************************************!
!********************************************************************!

