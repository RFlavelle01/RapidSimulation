!********************************************************************!
!********************************************************************!
!                                                                    !
!   mesh_module -- A module which contains mesh code                 !
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

module mesh_module

  ! Declare modules

  use precision

  implicit none

  ! Set the file the mesh will be stored in

  character(len=CP), parameter :: meshfile = 'mesh.dat'

  ! Set the loop variables
  
  integer (kind=IP) :: i, j, k

  ! Set up a type to contain the mesh data

  type mesh_type

     ! Size of the mesh in the three dimensions

     integer (kind=IP) :: nci, ncj, nck

     ! Set the cell centred points

     real    (kind=RP), allocatable :: x(:,:,:,:)

     ! Set the normals for the three relevant faces

     real    (kind=RP), allocatable :: nor(:,:,:,:)
     
     ! And set the volumes

     real    (kind=RP), allocatable :: vol(:,:,:,:)

     ! Set the minimum dx and dt for each cell

     real    (kind=RP), allocatable :: dx(:,:,:,:), dt(:,:,:,:)

  end type mesh_type

  ! Declare a mesh

  type(mesh_type) :: mesh

contains

  !*******************************************************************
  !*******************************************************************

  subroutine mesh_constructor()

    ! Declare variables

    integer (kind=IP) :: ni, nj, nk

    real    (kind=RP), allocatable :: xo(:), yo(:), zo(:)

    ! Open the file containing the mesh data

    open(file=meshfile, unit=11)

    ! Read the mesh size from the top of the file

    read(11,*) ni, nj, nk

    ! Then allocate memory for storing local points

    allocate(xo(ni+2), yo(nj+2), zo(nk+2))

    ! Read data into outer variables

    do i = 2, ni+1
       read(11,*) xo(i)
    end do
    do j = 2, nj+1
       read(11,*) yo(j)
    end do
    do k = 2, nk+1
       read(11,*) zo(k)
    end do

    ! Then set up the ghost cells via reflection

    xo(1) = 2d0*xo(2) - xo(3)
    yo(1) = 2d0*yo(2) - yo(3)
    zo(1) = 2d0*zo(2) - zo(3)

    xo(ni+2) = 2d0*xo(ni+1) - xo(ni)
    yo(nj+2) = 2d0*yo(nj+1) - yo(nj)
    zo(nk+2) = 2d0*zo(nk+1) - zo(nk)

    ! Set variable size for cells, including ghost cells

    mesh%nci = ni + 1
    mesh%ncj = nj + 1
    mesh%nck = nk + 1

    ! Allocate memory for storing main mesh

    allocate(mesh%x(3,mesh%nci,mesh%ncj,mesh%nck))

    allocate(mesh%nor(3,mesh%nci,mesh%ncj,mesh%nck))
    allocate(mesh%vol(1,mesh%nci,mesh%ncj,mesh%nck))

    allocate(mesh%dx(1,mesh%nci,mesh%ncj,mesh%nck))

    allocate(mesh%dt(1,mesh%nci,mesh%ncj,mesh%nck))

    ! Put data into cell-centred format

    do k = 1, mesh%nck
    do j = 1, mesh%ncj
    do i = 1, mesh%nci

       mesh%x(1,i,j,k) = 0.5d0 * ( xo(i)+xo(i+1) )
       mesh%x(2,i,j,k) = 0.5d0 * ( yo(j)+yo(j+1) )
       mesh%x(3,i,j,k) = 0.5d0 * ( zo(k)+zo(k+1) )

    end do
    end do
    end do

    ! Set up the normals 

    do k = 1, mesh%nck
    do j = 1, mesh%ncj
    do i = 1, mesh%nci

       mesh%nor(1,i,j,k) = ( yo(j+1)-yo(j) )*( zo(k+1)-zo(k) )
       mesh%nor(2,i,j,k) = ( xo(i+1)-xo(i) )*( zo(k+1)-zo(k) )
       mesh%nor(3,i,j,k) = ( xo(i+1)-xo(i) )*( yo(j+1)-yo(j) )

    end do
    end do
    end do

    ! Set up the cell volumes

    do k = 1, mesh%nck
    do j = 1, mesh%ncj
    do i = 1, mesh%nci

       mesh%vol(1,i,j,k) = ( xo(i+1)-xo(i) )*( yo(j+1)-yo(j) )*( zo(k+1)-zo(k) )

    end do
    end do
    end do

    ! Set up the minimum cell sizes

    do k = 1, mesh%nck
    do j = 1, mesh%ncj
    do i = 1, mesh%nci

       mesh%dx(1,i,j,k) = min( xo(i+1)-xo(i), yo(j+1)-yo(j), zo(k+1)-zo(k) )

    end do
    end do
    end do

    ! Deallocate local memory

    deallocate(xo, yo, zo)

    ! Close the file containing the mesh data

    close(11)

    ! Return to calling program

    return

  end subroutine mesh_constructor

end module mesh_module

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module mesh_data                                          !
!                                                                    !
!********************************************************************!
!********************************************************************!
