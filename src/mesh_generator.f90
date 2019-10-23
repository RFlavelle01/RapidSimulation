!********************************************************************!
!********************************************************************!
!                                                                    !
!   regular_rectilinear -- Generates a simple test mesh              !
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

program regular_rectilinear

  ! Declare modules

  use precision

  implicit none

  ! Declare variables

  integer(kind=IP) :: nx, ny, nz

  integer(kind=IP) :: ix, iy, iz

  real   (kind=RP) :: lx, ly, lz

  real   (kind=RP), allocatable :: xo(:), yo(:), zo(:)

  ! Set the size of the mesh

  nx = 200
  ny = 80
  nz = 80

  ! Set the domain of the mesh

  lx = 12d0 * 4d0 * atan(1d0)
  ly = 2d0 * 4d0 * atan(1d0)
  lz = 2d0 * 4d0 * atan(1d0)

  ! Allocate memory

  allocate(xo(nx), yo(nx), zo(nx))

  ! Set up the coordinates

  do ix = 1, nx
     xo(ix) = (ix-1)*lx / dble(nx-1) - 0.5d0*lx
  end do
  do iy = 1, ny
     yo(iy) = (iy-1)*ly / dble(ny-1) - 0.5d0*ly
  end do
  do iz = 1, nz
     zo(iz) = (iz-1)*lz / dble(nz-1) - 0.5d0*lz
  end do

  ! Write out the coordinates

  open(unit=12,file='mesh.dat')

  write(12,*) nx, ny, nz

  do ix = 1, nx
     write(12,*) xo(ix)
  end do

  do iy = 1, ny
     write(12,*) yo(iy)
  end do

  do iz = 1, nz
     write(12,*) zo(iz)
  end do

  close(12)

end program regular_rectilinear

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of program regular_rectilinear                               !
!                                                                    !
!********************************************************************!
!********************************************************************!

