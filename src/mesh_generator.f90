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

  !nx = 199
  !ny = 79
  !nz = 79
  
  ! Set the domain of the mesh

  lx = 12d0 * 4d0 * atan(1d0)
  ly = 2d0 * 4d0 * atan(1d0)
  lz = 2d0 * 4d0 * atan(1d0)

  ! Allocate memory

  allocate(xo(nx), yo(nx), zo(nx))

  ! Set up the coordinates

  do ix = 1, nx
    xo(ix) = (ix-1)*lx / dble(nx-1) - 0.5d0*lx
    !xo(ix) = (ix-dble((nx+1)/2))**3 * (1/dble((nx+1)/2 -1)**3 * 0.5d0 *lx) 
  
  !if ( ix <= (nx+1)/4) then
      
   ! xo(ix) = (ix-1)*(lx/3)/((nx+1)*0.25d0-1)-(lx*0.5d0)

  !else if (ix > (nx+1)/4 .AND. ix <= (nx+1)*3/4) then

   !  xo(ix) = ( ix-((nx+1)/4))*(lx/3)/((nx+1)*0.25d0)*0.5d0 - (lx*0.5d0/3)

  !else 

   !  xo(ix) = (ix-((nx+1)*3/4))*(lx/3)/((nx+1)*0.25d0-1)+(lx/6)

  !end if

  end do

  do iy = 1, ny
   yo(iy) = (iy-1)*ly / dble(ny-1) - 0.5d0*ly
   
   !if ( iy <= (ny+1)/4) then
      
    !yo(iy) = (iy-1)*(ly/3)/((ny+1)*0.25d0-1)-(ly*0.5d0)

  !else if (iy > (ny+1)/4 .AND. iy <= (ny+1)*3/4) then

   !  yo(iy) = (iy-((ny+1)/4))*(ly/3)/((ny+1)*0.25d0)*0.5d0 - (ly*0.5d0/3)

  !else 

   !  yo(iy) = (iy-((ny+1)*3/4))*(ly/3)/((ny+1)*0.25d0-1)+(ly/6)

  !end if   
  end do


  do iz = 1, nz
     zo(iz) = (iz-1)*lz / dble(nz-1) - 0.5d0*lz
    
  !if ( iz <= (nz+1)/4) then
      
   !zo(iz) = (iz-1)*(lz/3)/((nz+1)*0.25d0-1)-(lz*0.5d0)

  !else if (iz > (nz+1)/4 .AND. iz <= (nz+1)*3/4) then

   !  zo(iz) = (iz-((nz+1)/4))*(lz/3)/((nz+1)*0.25d0)*0.5d0 - (lz*0.5d0/3)

  !else 

   !  zo(iz) = (iz-((nz+1)*3/4))*(lz/3)/((nz+1)*0.25d0-1)+(lz/6)

  !end if 

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

