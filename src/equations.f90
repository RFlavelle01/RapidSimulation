!********************************************************************!
!********************************************************************!
!                                                                    !
!   equations_module -- Contains the equations to be solved          !
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

module equations_module

  ! Declare modules
  
  use precision
  
  implicit none

  ! Set the number of equations - increased form 5 to 9 for wall distance computation

  integer(kind=IP), parameter :: npdes = 9

  ! Set some major flow parameters

  real   (kind=RP), parameter :: gam = 1.4d0
  real   (kind=RP), parameter :: gm1 = gam - 1d0
  real   (kind=RP), parameter :: gmi = 1d0 / gm1
  real   (kind=RP), parameter :: ggi = gam / gm1
  real   (kind=RP), parameter :: ggo = gm1 / gam

  real   (kind=RP), parameter :: rgas = 287.1d0

  real   (kind=RP), parameter :: nul = 0.0000155d0

  real   (kind=RP), parameter :: pr  = 0.71d0

  ! Set some major geometric parameters

  real   (kind=RP), parameter :: pi = 4d0 * atan(1d0)


contains

  !*******************************************************************
  !*******************************************************************

  function u2q(u) result(q)

    implicit none

    ! Declare variables

    real   (kind=RP) :: u(npdes), q(npdes)

    ! Set the primitive variables from the conserved variables

    q(1) = u(1)
    q(2) = u(2) / q(1)
    q(3) = u(3) / q(1)
    q(4) = u(4) / q(1)
    q(5) = gm1*(u(5) - 0.5d0*q(1)*(q(2)*q(2)+q(3)*q(3)+q(4)*q(4)))

    !Wall distances

    q(6) = u(6)
    q(7) = u(7)
    q(8) = u(8)
    q(9) = u(9)

    ! Return cleanly

    return

  end function u2q

  !*******************************************************************
  !*******************************************************************

  function lambda(u) result(lam)

    implicit none

    ! Declare variables

    real   (kind=RP) :: u(npdes), lam

    real   (kind=RP) :: ri, ui, vi, wi, pi, ai

    ! Move convenient form

    ri = u(1)
    ui = u(2) / ri
    vi = u(3) / ri
    wi = u(4) / ri
    pi = gm1*(u(5) - 0.5d0*ri*(ui*ui + vi*vi + wi*wi))

    ! Compute speed of sound

    ai = sqrt(gam*pi/ri)

    ! Set the maximum wavespeed

    lam = sqrt(ui*ui + vi*vi + wi*wi) + ai

    ! Return cleanly

    return

  end function lambda


  !*******************************************************************
  !*******************************************************************

  function xflux(uL, qL, qpL, uR, qR, qpR, n, dist) result(f)

    implicit none

    ! Declare variables

    real(kind=RP) :: f(npdes), n(3)
    real(kind=RP) :: uL(npdes), qL(npdes), qpL(npdes,3)
    real(kind=RP) :: uR(npdes), qR(npdes), qpR(npdes,3)

    real(kind=RP) :: q(npdes), qp(npdes,3), smij(3,3), tij(3,3), qx, tx, mul

    real(kind=RP) :: HL, HR, RT, umix, vmix, wmix, Hmix, amix, smax, vall, smth

    integer(kind=IP) :: ipde, i, j

    ! Turbulence variables

    real(kind=RP) :: L, mut, rij(3,3), qtx, sumsmij, dist(1)

    ! Set the wave speeds in the relevant direction (Rusanov)

    HL = (uL(5) + qL(5)) / qL(1)
    HR = (uR(5) + qR(5)) / qR(1)

    RT = sqrt(qR(1) / qL(1))

    umix = (qL(2) + RT*qR(2)) / (1d0 + RT)
    vmix = (qL(3) + RT*qR(3)) / (1d0 + RT)
    wmix = (qL(4) + RT*qR(4)) / (1d0 + RT)

    Hmix = (HL    + RT*HR   ) / (1d0 + RT)

    vall = sqrt(umix*umix+vmix*vmix+wmix*wmix)
    amix = sqrt(gm1*(Hmix - 0.5d0*(vall*vall)))

    smax = abs(vall) + amix
    
    smth = 0.50d0

    ! Set the fluxes

    f(1) = 0.5d0*(qL(1)*qL(2)               + (qR(1)*qR(2)              )) - smth*smax*(uR(1)-uL(1))
    f(2) = 0.5d0*(qL(1)*qL(2)*qL(2) + qL(5) + (qR(1)*qR(2)*qR(2) + qR(5))) - smth*smax*(uR(2)-uL(2))
    f(3) = 0.5d0*(qL(1)*qL(2)*qL(3)         + (qR(1)*qR(2)*qR(3)        )) - smth*smax*(uR(3)-uL(3))
    f(4) = 0.5d0*(qL(1)*qL(2)*qL(4)         + (qR(1)*qR(2)*qR(4)        )) - smth*smax*(uR(4)-uL(4))
    f(5) = 0.5d0*(qL(1)*qL(2)*HL            + (qR(1)*qR(2)*HR           )) - smth*smax*(uR(5)-uL(5))


    !Wall distance

    f(6) = -(0.5d0*(qL(7)+qR(7)))
    f(7) = -0.5d0*(qL(6)+qR(6))
    f(8) = 0d0
    f(9) = 0d0

    ! Get the mean gradients

    do ipde = 1, npdes
       q(ipde) = 0.5d0 * (qL(ipde) + qR(ipde))
    end do
    do i = 1, 3
    do ipde = 1, npdes
       qp(ipde,i) = 0.5d0 * (qpL(ipde,i) + qpR(ipde,i))
    end do
    end do
    
    ! Get the mean modifed rate of strain tensor

    do j = 1, 3
    do i = 1, 3
       smij(i,j) = 0.5d0 * (qp(i+1,j) + qp(j+1,i)) - (1d0/3d0)*(qp(1+1,1)+qp(2+1,2)+qp(3+1,3))
    end do
    end do

    ! Compute the viscosity

    mul = q(1)*nul

    ! Compute the viscous stresses

    do j = 1, 3
    do i = 1, 3
       tij(i,j) = 2d0 * mul * smij(i,j)
    end do
    end do

    ! Compute the thermal gradient

    tx = ggi * (qp(5,1)*q(1) - qp(1,1)*q(5)) / (q(1)*q(1))

    ! Find the diffusion of heat

    qx = -(mul / pr) * tx

    ! Turbulence modelling

    ! Set the length scale
     
    L = 0.25 * pi
    !L = min(0.41*dist(1), 0.25 *pi)
    
    ! Compute turbulent eddy viscosity

    sumsmij = (smij(1,1)**2 + smij(2,1)**2 + smij(1,2)**2 + smij(2,2)**2 *smij(3,1)**2 + smij(3,2)**2 + smij(3,3)**2)
    
    mut  = q(1) * sqrt( L**2 * 2 *sumsmij)

    ! Calculate the Reynolds stress tensor
    
    do j = 1, 3
    do i = 1, 3
       rij(i,j) = 2 * mut * smij(i,j)
    end do
    end do
    
    ! Find the turbulent diffusion of heat
    
    qtx = -(mut / pr) * tx

    ! Add on the viscous terms + zero equation turbulence terms

    f(1) = f(1)
    f(2) = f(2) - tij(1,1) !- rij(1,1)
    f(3) = f(3) - tij(1,2) !- rij(1,2)
    f(4) = f(4) - tij(1,3) !- rij(1,3)
    f(5) = f(5) - q(2)*tij(1,1) - q(3)*tij(1,2) - q(4)*tij(1,3) + qx !+ qtx - q(2)*rij(1,1) - q(3)*rij(1,2) - q(4)*rij(1,3)
    f(6) = f(6)
    f(7) = f(7)
    f(8) = f(8)
    f(9) = f(9)

    ! And dot with the normals

    do ipde = 1, npdes
       f(ipde) = f(ipde) * n(1)
    end do

    ! Return cleanly

    return

  end function xflux

  !*******************************************************************
  !*******************************************************************

  function yflux(uL, qL, qpL, uR, qR, qpR, n, dist) result(g)

    implicit none

    ! Declare variables

    real(kind=RP) :: g(npdes), n(3)
    real(kind=RP) :: uL(npdes), qL(npdes), qpL(npdes,3)
    real(kind=RP) :: uR(npdes), qR(npdes), qpR(npdes,3)

    real(kind=RP) :: q(npdes), qp(npdes,3), smij(3,3), tij(3,3), qy, ty, mul

    real(kind=RP) :: HL, HR, RT, umix, vmix, wmix, Hmix, amix, smax, vall, smth

    integer(kind=IP) :: ipde, i, j

    ! Turbulence variables

    real(kind=RP) :: L, mut, rij(3,3), qty, sumsmij, dist(1)

    ! Set the wave speeds in the relevant direction (Rusanov)

    HL = (uL(5) + qL(5)) / qL(1)
    HR = (uR(5) + qR(5)) / qR(1)

    RT = sqrt(qR(1) / qL(1))

    umix = (qL(2) + RT*qR(2)) / (1d0 + RT)
    vmix = (qL(3) + RT*qR(3)) / (1d0 + RT)
    wmix = (qL(4) + RT*qR(4)) / (1d0 + RT)

    Hmix = (HL    + RT*HR   ) / (1d0 + RT)

    vall = sqrt(umix*umix+vmix*vmix+wmix*wmix)
    amix = sqrt(gm1*(Hmix - 0.5d0*(vall*vall)))

    smax = abs(vall) + amix

    smth = 0.50d0

    ! Set the fluxes

    g(1) = 0.5d0*(qL(1)*qL(3)               + (qR(1)*qR(3)              )) - smth*smax*(uR(1)-uL(1))
    g(2) = 0.5d0*(qL(1)*qL(3)*qL(2)         + (qR(1)*qR(3)*qR(2)        )) - smth*smax*(uR(2)-uL(2))
    g(3) = 0.5d0*(qL(1)*qL(3)*qL(3) + qL(5) + (qR(1)*qR(3)*qR(3) + qR(5))) - smth*smax*(uR(3)-uL(3))
    g(4) = 0.5d0*(qL(1)*qL(3)*qL(4)         + (qR(1)*qR(3)*qR(4)        )) - smth*smax*(uR(4)-uL(4))
    g(5) = 0.5d0*(qL(1)*qL(3)*HL            + (qR(1)*qR(3)*HR           )) - smth*smax*(uR(5)-uL(5))
    
    !Wall distances 
   
    g(6) = -(0.5d0*(qL(8)+qR(8)))
    g(7) = 0d0
    g(8) = -0.5d0*(qL(6)+qR(6))
    g(9) = 0d0
 
    ! Get the mean gradients

    do ipde = 1, npdes
       q(ipde) = 0.5d0 * (qL(ipde) + qR(ipde))
    end do
    do i = 1, 3
    do ipde = 1, npdes
       qp(ipde,i) = 0.5d0 * (qpL(ipde,i) + qpR(ipde,i))
    end do
    end do
    
    ! Get the mean modified rate of strain tensor

    do j = 1, 3
    do i = 1, 3
       smij(i,j) = 0.5d0 * (qp(i+1,j) + qp(j+1,i)) - (1d0/3d0)*(qp(1+1,1)+qp(2+1,2)+qp(3+1,3))
    end do
    end do

    ! Compute the viscosity

    mul = q(1)*nul

    ! Compute the viscous stresses

    do j = 1, 3
    do i = 1, 3
       tij(i,j) = 2d0 * mul * smij(i,j)
    end do
    end do

    ! Compute the thermal gradient

    ty = ggi * (qp(5,2)*q(1) - qp(1,2)*q(5)) / (q(1)*q(1))

    ! Find the diffusion of heat

    qy = -(mul / pr) * ty

    ! Turbulence modelling

    ! Set the length scale

    L = 0.25 * pi
   ! L = min(0.41 * dist(1), 0.25 *pi)

    ! Compute turbulent eddy viscosity

    sumsmij = (smij(1,1)**2 + smij(2,1)**2 + smij(1,2)**2 + smij(2,2)**2 *smij(3,1)**2 + smij(3,2)**2 + smij(3,3)**2)

    mut = q(1) * sqrt(L**2 * 2 * sumsmij)

    ! Calculate the Reynolds stress tensor
    
    do j = 1, 3
    do i = 1, 3
       rij(i,j) = 2 * mut * smij(i,j)
    end do
    end do
    
    ! Find the turbulent diffusion of heat
    
    qty = -(mut / pr) * ty

    ! Add on the viscous terms + turbulent stresses

    g(1) = g(1)
    g(2) = g(2) - tij(2,1) !- rij(2,1)
    g(3) = g(3) - tij(2,2) !- rij(2,2)
    g(4) = g(4) - tij(2,3) !- rij(2,3)
    g(5) = g(5) - q(2)*tij(2,1) - q(3)*tij(2,2) - q(4)*tij(2,3) + qy !+ qty - q(2)*rij(2,1) - q(3)*rij(2,2) - q(4)*rij(2,3)
    g(6) = g(6)
    g(7) = g(7)
    g(8) = g(8)
    g(9) = g(9)

    ! And dot with the normals

    do ipde = 1, npdes
       g(ipde) = g(ipde) * n(2)
    end do

  end function yflux

  !*******************************************************************
  !*******************************************************************

  function zflux(uL, qL, qpL, uR, qR, qpR, n, dist) result(h)

    implicit none

    ! Declare variables

    real(kind=RP) :: h(npdes), n(3)
    real(kind=RP) :: uL(npdes), qL(npdes), qpL(npdes,3)
    real(kind=RP) :: uR(npdes), qR(npdes), qpR(npdes,3)

    real(kind=RP) :: q(npdes), qp(npdes,3), smij(3,3), tij(3,3), qz, tz, mul

    real(kind=RP) :: HL, HR, RT, umix, vmix, wmix, Hmix, amix, smax, vall, smth

    integer(kind=IP) :: ipde, i, j

    ! Turbulence variables

    real(kind=RP) :: L, mut, rij(3,3), qtz, sumsmij, dist(1)

    ! Set the wave speeds in the relevant direction (Rusanov)

    HL = (uL(5) + qL(5)) / qL(1)
    HR = (uR(5) + qR(5)) / qR(1)

    RT = sqrt(qR(1) / qL(1))

    umix = (qL(2) + RT*qR(2)) / (1d0 + RT)
    vmix = (qL(3) + RT*qR(3)) / (1d0 + RT)
    wmix = (qL(4) + RT*qR(4)) / (1d0 + RT)

    Hmix = (HL    + RT*HR   ) / (1d0 + RT)
    
    vall = sqrt(umix*umix+vmix*vmix+wmix*wmix)
    amix = sqrt(gm1*(Hmix - 0.5d0*(vall*vall)))

    smax = abs(vall) + amix

    smth = 0.50d0

    ! Set the fluxes

    h(1) = 0.5d0*(qL(1)*qL(4)               + (qR(1)*qR(4)              )) - smth*smax*(uR(1)-uL(1))
    h(2) = 0.5d0*(qL(1)*qL(4)*qL(2)         + (qR(1)*qR(4)*qR(2)        )) - smth*smax*(uR(2)-uL(2))
    h(3) = 0.5d0*(qL(1)*qL(4)*qL(3)         + (qR(1)*qR(4)*qR(3)        )) - smth*smax*(uR(3)-uL(3))
    h(4) = 0.5d0*(qL(1)*qL(4)*qL(4) + qL(5) + (qR(1)*qR(4)*qR(4) + qR(5))) - smth*smax*(uR(4)-uL(4))
    h(5) = 0.5d0*(qL(1)*qL(4)*HL            + (qR(1)*qR(4)*HR           )) - smth*smax*(uR(5)-uL(5))
    
    !Wall distances
    
    h(6) = -(0.5d0*(qL(9)+qR(9)))
    h(7) = 0d0
    h(8) = 0d0
    h(9) = -0.5d0*(qL(6)+qR(6))
    
    ! Get the mean gradients

    do ipde = 1, npdes
       q(ipde) = 0.5d0 * (qL(ipde) + qR(ipde))
    end do
    do i = 1, 3
    do ipde = 1, npdes
       qp(ipde,i) = 0.5d0 * (qpL(ipde,i) + qpR(ipde,i))
    end do
    end do
    
    ! Get the mean modified rate of strain tensor

    do j = 1, 3
    do i = 1, 3
       smij(i,j) = 0.5d0 * (qp(i+1,j) + qp(j+1,i)) - (1d0/3d0)*(qp(1+1,1)+qp(2+1,2)+qp(3+1,3))
    end do
    end do

    ! Compute the viscosity

    mul = q(1)*nul

    ! Compute the viscous stresses

    do j = 1, 3
    do i = 1, 3
       tij(i,j) = 2d0 * mul * smij(i,j)
    end do
    end do

    ! Compute the thermal gradient

    tz = ggi * (qp(5,2)*q(1) - qp(1,2)*q(5)) / (q(1)*q(1))

    ! Find the diffusion of heat

    qz = -(mul / pr) * tz

    ! Turbulence modelling

    ! Set the length scale

    L = 0.25 * pi
    !L = min( 0.41 * dist(1), 0.25 * pi)

    ! Compute turbulent eddy viscosity

    sumsmij = (smij(1,1)**2 + smij(2,1)**2 + smij(1,2)**2 + smij(2,2)**2 *smij(3,1)**2 + smij(3,2)**2 + smij(3,3)**2)

    mut = q(1) * sqrt(L**2 * 2 * sumsmij)

    ! Calculate the Reynolds stress tensor
    
    do j = 1, 3
    do i = 1, 3
       rij(i,j) = 2 * mut * smij(i,j)
    end do
    end do
    
    ! Find the turbulent diffusion of heat
    
    qtz = -(mut / pr) * tz

    ! Add on the viscous terms

    h(1) = h(1)
    h(2) = h(2) - tij(3,1) !- rij(3,1)
    h(3) = h(3) - tij(3,2) !- rij(3,2)
    h(4) = h(4) - tij(3,3) !- rij(3,3)
    h(5) = h(5) - q(2)*tij(3,1) - q(3)*tij(3,2) - q(4)*tij(3,3) + qz !+ qtz -q(2)*rij(3,1) - q(3)*rij(3,2) - q(4)*rij(3,3)
    h(6) = h(6)
    h(7) = h(7)
    h(8) = h(8)
    h(9) = h(9)

    ! And dot with the normals

    do ipde = 1, npdes
       h(ipde) = h(ipde) * n(3)
    end do

  end function zflux

  !*******************************************************************
  !*******************************************************************

  function bc_type_a(uIn) result(uBn)

    implicit none

    ! Declare variables

    real   (kind=RP) :: uIn(npdes), uBn(npdes)

    real   (kind=RP) :: p0, t0, nx, ny, nz
    real   (kind=RP) :: ri, ui, vi, wi, pi

    real   (kind=RP) :: ci, h0, riemn, aqe, bqe, cqe, cb, mb, sb

    real   (kind=RP) :: rb, ub, vb, wb, pb, tb, eb

    ! Set the stagnation pressure and temperature

    p0 = 115000d0
    t0 = 300d0

    nx = 1d0
    ny = 0d0
    nz = 0d0

    ! Give them more convenient names

    ri = uIn(1)
    ui = uIn(2) / ri
    vi = uIn(3) / ri
    wi = uIn(4) / ri
    pi = gm1*(uIn(5) - 0.5d0*ri*(ui*ui + vi*vi + wi*wi))

    ! Compute the Riemann invariant according to the internal variables

    ci = sqrt(gam*pi/ri)
    h0 = (pi/ri) * ggi + 0.5d0*(ui*ui + vi*vi + wi*wi)
    riemn = sqrt(ui*ui + vi*vi + wi*wi) + 2d0*ci*gmi
    
    ! Solve the quadratic equation to find the new speed of sound

    aqe = 1d0 + 2d0*gmi
    bqe = -2d0*riemn
    cqe = 0.5d0*gm1*(riemn*riemn-2d0*h0)

    ! The new speed of sound is the largest of the two roots

    cb = -0.5d0*bqe/aqe + 0.5d0*sqrt(bqe*bqe-4d0*aqe*cqe)/aqe

    ! Update to get the flow speed and Mach number

    sb = riemn - 2d0*cb*gmi
    mb = sb / cb

    ! And recompute the primitive variables

    ub = sb * nx
    vb = sb * ny
    wb = sb * nz

    pb = p0 * (1d0 + 0.5d0*gm1*mb*mb)**(-ggi)
    tb = t0 * (pb/p0)**ggo

    rb = pb / (rgas * tb)

    eb = pb / (rb * gm1) + 0.5d0*(ub*ub + vb*vb + wb*wb)

    ! And finally construct the conserved variables

    uBn(1) = rb
    uBn(2) = rb*ub
    uBn(3) = rb*vb
    uBn(4) = rb*wb
    uBn(5) = rb*eb
    
    !Wall distances
 
    uBn(6) = -uIn(6)
    uBn(7) = uIn(7)
    uBn(8) = uIn(8)
    uBn(9) = uIn(9)

    ! Return cleanly

    return

  end function bc_type_a

  !*******************************************************************
  !*******************************************************************

  function bc_type_b(uIn) result(uBn)

    implicit none

    ! Declare variables

    real   (kind=RP) :: uIn(npdes), uBn(npdes)

    real   (kind=RP) :: ri, ui, vi, wi, pi

    real   (kind=RP) :: rb, ub, vb, wb, pb, eb
    

    ! Set the back pressure

    pb = 100000d0
    !pb = 108727d0
    
    ! Give them more convenient names

    ri = uIn(1)
    ui = uIn(2) / ri
    vi = uIn(3) / ri
    wi = uIn(4) / ri
    pi = gm1*(uIn(5) - 0.5d0*ri*(ui*ui + vi*vi + wi*wi))    

    ! Set the back variables

    rb = ri * pb / pi
    ub = ui
    vb = vi
    wb = wi
    eb = pb / (rb * gm1) + 0.5d0*(ub*ub + vb*vb + wb*wb)

    ! And finally construct the conserved variables

    uBn(1) = rb
    uBn(2) = rb*ub
    uBn(3) = rb*vb
    uBn(4) = rb*wb
    uBn(5) = rb*eb
    
    ! Wall distances

    uBn(6) = -uIn(6)
    uBn(7) = uIn(7)
    uBn(8) = uIn(8)
    uBn(9) = uIn(9)

    ! Return cleanly
   
    return

  end function bc_type_b

  !*******************************************************************
  !*******************************************************************

  function bc_type_c(uIn, nIn) result(uBn)

    implicit none

    ! Declare variables

    real   (kind=RP) :: uIn(npdes), nIn(3), uBn(npdes)

    real   (kind=RP) :: ri, ui, vi, wi, pi

    real   (kind=RP) :: rb, ub, vb, wb, pb, eb, ucdn

    ! Give them more convenient names

    ri = uIn(1)
    ui = uIn(2) / ri
    vi = uIn(3) / ri
    wi = uIn(4) / ri
    pi = gm1*(uIn(5) - 0.5d0*ri*(ui*ui + vi*vi + wi*wi))    

    ! Dot the velocity with the wall

    ucdn = ui*nIn(1) + vi*nIn(2) + wi*nIn(3)

    ! And eliminate the normal vectors

    ub = ui - 2d0*ucdn*nIn(1)
    vb = vi - 2d0*ucdn*nIn(2)
    wb = wi - 2d0*ucdn*nIn(3)

    ! Otherwise, set a zero gradient condition

    rb = ri
    pb = pi

    eb = pb / (rb * gm1) + 0.5d0*(ub*ub + vb*vb + wb*wb)

    ! And finally construct the conserved variables

    uBn(1) = rb
    uBn(2) = rb*ub
    uBn(3) = rb*vb
    uBn(4) = rb*wb
    uBn(5) = rb*eb

    ! Wall distances

    uBn(6) = -uIn(6)
    uBn(7) = uIn(7)
    uBn(8) = uIn(8)
    uBn(9) = uIn(9)

  end function bc_type_c

  !*******************************************************************
  !*******************************************************************

  function sources(q,f,vol) result(s)

    implicit none

    ! Declare variables

    real   (kind=RP) :: q(npdes), f(npdes), vol(1), s(npdes)

    ! Set the source terms (multiply by volume since this is an FV code)

    s(1) = 0d0  *vol(1)
    s(2) = f(2) *vol(1)
    s(3) = f(3) *vol(1)
    s(4) = f(4) *vol(1)
    s(5) = 0    *vol(1)
    
    !Wall distances
    s(6) = (1 + f(6))*vol(1)
    s(7) = -q(7)*vol(1)
    s(8) = -q(8)*vol(1)
    s(9) = -q(9)*vol(1)

    ! Return

    return

  end function sources

end module equations_module

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module equations_module                                   !
!                                                                    !
!********************************************************************!
!********************************************************************!

