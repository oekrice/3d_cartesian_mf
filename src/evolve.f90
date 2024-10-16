!*******************************************************************************

MODULE evolve
!*******************************************************************************
! Initialise the simulation. Read in/compute the initial condition (python code for that I think?), read in the parameters and all will be well
!*******************************************************************************
    USE shared_data
    USE mpi_tools
    !USE output
    !USE pressure
    IMPLICIT NONE

!*******************************************************************************
CONTAINS

SUBROUTINE timestep()


    !New timestep with a quasi-implicit bit. Using a predictor for the magnfield of the new bit but not the old
    CALL calculate_magnetic()

    CALL magnetic_boundary()

    !CALL check_solenoidal()
    !CALL calculate_jp()

    CALL calculate_current()

    CALL j_to_gridpts()
    CALL b_to_gridpts()

    CALL calculate_velocity()

    CALL add_boundary_flows()

    CALL calculate_electric()

    CALL MPI_Barrier(comm,ierr)  !Wait for t to be broadcast everywhere.

END SUBROUTINE timestep


SUBROUTINE calculate_magnetic()

    IMPLICIT NONE

    !Determine the interior points from the vector potential.
    !No boundary conditions here (do that in a bit)

    ! INTERIOR POINTS (DON'T USE INFORMATION FROM A THAT DOESN'T EXIST)

    bx(0:nx, 1:ny,1:nz) = (az(0:nx,1:ny,1:nz) - az(0:nx, 0:ny-1,1:nz))/dy - (ay(0:nx,1:ny,1:nz) - ay(0:nx,1:ny,0:nz-1))/dz

    by(1:nx, 0:ny,1:nz) = (ax(1:nx,0:ny,1:nz) - ax(1:nx, 0:ny,0:nz-1))/dz - (az(1:nx,0:ny,1:nz) - az(0:nx-1,0:ny,1:nz))/dx

    bz(1:nx, 1:ny,0:nz) = (ay(1:nx,1:ny,0:nz) - ay(0:nx-1, 1:ny,0:nz))/dx - (ax(1:nx,1:ny,0:nz) - ax(1:nx,0:ny-1,0:nz))/dy


END SUBROUTINE calculate_magnetic

SUBROUTINE magnetic_boundary()

    IMPLICIT NONE

    CALL bfield_mpi()
    CALL MPI_BARRIER(comm,ierr)

    !Boundary conditions on the magnetic field.
    !Uses the no-horizontal-current condition for the first ghost cells and then the solenoidal condition for the next lot.

    !LOWER BOUNDARY (Zero current)
    if (z_rank == 0) then
    by(0:nx+1,0:ny,0) = by(0:nx+1,0:ny,1) - dz*(bz(0:nx+1,1:ny+1,0) - bz(0:nx+1, 0:ny,0))/dy
    bx(0:nx, 0:ny+1,0) = bx(0:nx,0:ny+1,1) - dz*(bz(1:nx+1,0:ny+1,0) - bz(0:nx,0:ny+1,0))/dx
    end if

    !UPPER BOUNDARY (Zero Current)
    if (z_rank == z_procs-1) then
    by(0:nx+1,0:ny,nz+1) = by(0:nx+1,0:ny,nz) + dz*(bz(0:nx+1,1:ny+1,nz) - bz(0:nx+1, 0:ny,nz))/dy
    bx(0:nx, 0:ny+1,nz+1) = bx(0:nx,0:ny+1,nz) + dz*(bz(1:nx+1,0:ny+1,nz) - bz(0:nx,0:ny+1,nz))/dx
    end if

    !x boundaries (Zero current, and zero flux)
    if (x_rank == 0) then
    bz(0,0:ny+1,0:nz) = bz(1,0:ny+1,0:nz) - dx*(bx(0,0:ny+1,1:nz+1) - bx(0, 0:ny+1,0:nz))/dz
    by(0,0:ny,0:nz+1) = by(1, 0:ny,0:nz+1) - dx*(bx(0,1:ny+1,0:nz+1) - bx(0,0:ny,0:nz+1))/dy
    end if

    if (x_rank == x_procs-1) then
    bz(nx+1,0:ny+1,0:nz) = bz(nx,0:ny+1,0:nz) + dx*(bx(nx,0:ny+1,1:nz+1) - bx(nx, 0:ny+1,0:nz))/dz
    by(nx+1,0:ny,0:nz+1) = by(nx, 0:ny,0:nz+1) + dx*(bx(nx,1:ny+1,0:nz+1) - bx(nx,0:ny,0:nz+1))/dy
    end if

    !y boundaries (Zero current, and zero flux)
    if (y_rank == 0) then
    bz(0:nx+1,0,0:nz) = bz(0:nx+1, 1,0:nz) - dy*(by(0:nx+1,0,1:nz+1) - by(0:nx+1,0,0:nz))/dz
    bx(0:nx,0,0:nz+1) = bx(0:nx,1,0:nz+1) - dy*(by(1:nx+1,0,0:nz+1) - by(0:nx, 0,0:nz+1))/dx
    end if

    if (y_rank == y_procs-1) then
    bz(0:nx+1,ny+1,0:nz) = bz(0:nx+1, ny,0:nz) + dy*(by(0:nx+1,ny,1:nz+1) - by(0:nx+1,ny,0:nz))/dz
    bx(0:nx,ny+1,0:nz+1) = bx(0:nx,ny,0:nz+1) + dy*(by(1:nx+1,ny,0:nz+1) - by(0:nx, ny,0:nz+1))/dx
    end if

    CALL MPI_BARRIER(comm,ierr)

    !_________________________________________________
    !Solenoidal condition, using the above information
    !Isn't necessary to do more MPI here, I suppose. So I won't.

    !LOWER BOUNDARY (Zero current)
    bz(0:nx+1,0:ny+1,-1) = (1.0_num/(dx*dy))*(bz(0:nx+1,0:ny+1,0)*dx*dy - bx(-1:nx,0:ny+1,0)*dy*dz + bx(0:nx+1,0:ny+1,0)*dy*dz - by(0:nx+1,-1:ny,0)*dx*dz + by(0:nx+1,0:ny+1,0)*dx*dz)

    !UPPER BOUNDARY (Zero Current)
    bz(0:nx+1,0:ny+1,nz+1) = (1.0_num/(dx*dy))*(bz(0:nx+1,0:ny+1,nz)*dx*dy + bx(-1:nx,0:ny+1,nz+1)*dy*dz - bx(0:nx+1,0:ny+1,nz+1)*dy*dz + by(0:nx+1,-1:ny,nz+1)*dx*dz - by(0:nx+1,0:ny+1,nz+1)*dx*dz)

    !x boundaries (Zero current, and zero flux)

    bx(-1,0:ny+1,0:nz+1) = (1.0_num/(dy*dz))*(bx(0,0:ny+1,0:nz+1)*dy*dz - by(0,-1:ny,0:nz+1)*dx*dz + by(0,0:ny+1,0:nz+1)*dx*dz - bz(0,0:ny+1,-1:nz)*dx*dy + bz(0,0:ny+1,0:nz+1)*dx*dy)


    bx(nx+1,0:ny+1,0:nz+1) = (1.0_num/(dy*dz))*(bx(nx,0:ny+1,0:nz+1)*dy*dz + by(nx+1,-1:ny,0:nz+1)*dx*dz - by(nx+1,0:ny+1,0:nz+1)*dx*dz + bz(nx+1,0:ny+1,-1:nz)*dx*dy - bz(nx+1,0:ny+1,0:nz+1)*dx*dy)


    !y boundaries (Zero current, and zero flux)
    by(0:nx+1,-1,0:nz+1) = (1.0_num/(dx*dz))*(by(0:nx+1,0,0:nz+1)*dx*dz - bx(-1:nx,0,0:nz+1)*dy*dz + bx(0:nx+1,0,0:nz+1)*dy*dz - bz(0:nx+1,0,-1:nz)*dx*dy + bz(0:nx+1,0,0:nz+1)*dx*dy)


    by(0:nx+1,ny+1,0:nz+1) = (1.0_num/(dx*dz))*(by(0:nx+1,ny,0:nz+1)*dx*dz + bx(-1:nx,ny+1,0:nz+1)*dy*dz - bx(0:nx+1,ny+1,0:nz+1)*dy*dz + bz(0:nx+1,ny+1,-1:nz)*dx*dy - bz(0:nx+1,ny+1,0:nz+1)*dx*dy)

END SUBROUTINE magnetic_boundary

SUBROUTINE calculate_current()

    IMPLICIT NONE
    !Determine the current from the magnetic field (after boundary conditions etc.)

    jx(0:nx+1, 0:ny,0:nz) = (bz(0:nx+1,1:ny+1,0:nz) - bz(0:nx+1, 0:ny,0:nz))/dy - (by(0:nx+1,0:ny,1:nz+1) - by(0:nx+1,0:ny,0:nz))/dz

    jy(0:nx, 0:ny+1,0:nz) = (bx(0:nx,0:ny+1,1:nz+1) - bx(0:nx, 0:ny+1,0:nz))/dz - (bz(1:nx+1,0:ny+1,0:nz) - bz(0:nx,0:ny+1,0:nz))/dx

    jz(0:nx, 0:ny,0:nz+1) = (by(1:nx+1,0:ny,0:nz+1) - by(0:nx, 0:ny,0:nz+1))/dx - (bx(0:nx,1:ny+1,0:nz+1) - bx(0:nx,0:ny,0:nz+1))/dy

END SUBROUTINE calculate_current

SUBROUTINE j_to_gridpts
    !Averages the current field to raw grid points
    !Should only need to average in one direction
    IMPLICIT NONE
    jx1(0:nx,0:ny,0:nz) = 0.5_num*(jx(1:nx+1,0:ny,0:nz) + jx(0:nx,0:ny,0:nz))
    jy1(0:nx,0:ny,0:nz) = 0.5_num*(jy(0:nx,1:ny+1,0:nz) + jy(0:nx,0:ny,0:nz))
    jz1(0:nx,0:ny,0:nz) = 0.5_num*(jz(0:nx,0:ny,1:nz+1) + jz(0:nx,0:ny,0:nz))

END SUBROUTINE j_to_gridpts

SUBROUTINE b_to_gridpts
    !Averages the magnetic field to raw grid points
    !Need to average in two dimension
    IMPLICIT NONE
    bx1(0:nx,0:ny,0:nz) = 0.25_num*(bx(0:nx,0:ny,0:nz) + bx(0:nx,1:ny+1,0:nz) + bx(0:nx,0:ny,1:nz+1) + bx(0:nx,1:ny+1,1:nz+1))
    by1(0:nx,0:ny,0:nz) = 0.25_num*(by(0:nx,0:ny,0:nz) + by(1:nx+1,0:ny,0:nz) + by(0:nx,0:ny,1:nz+1) + by(1:nx+1,0:ny,1:nz+1))
    bz1(0:nx,0:ny,0:nz) = 0.25_num*(bz(0:nx,0:ny,0:nz) + bz(1:nx+1,0:ny,0:nz) + bz(0:nx,1:ny+1,0:nz) + bz(1:nx+1,1:ny+1,0:nz))

END SUBROUTINE b_to_gridpts

SUBROUTINE calculate_velocity
    !Calculates the magnetofrictional velocity
    IMPLICIT NONE
    b2 = bx1**2 + by1**2 + bz1**2 !B squared

    nu(:,:,:) = nu0

    if (abs(mf_delta) < 1e-10) then !No softening
        soft = b2
    else !Softening.
        soft = (b2 + mf_delta*exp(-b2/mf_delta))
    end if

    vx = nu*(jy1*bz1 - jz1*by1)/soft
    vy = nu*(jz1*bx1 - jx1*bz1)/soft
    vz = nu*(jx1*by1 - jy1*bx1)/soft

    if (z_down < 0) then
        vx(:,:,0) = 0.0_num; vy(:,:,0) = 0.0_num; vz(:,:,0) = 0.0_num
    end if

END SUBROUTINE calculate_velocity

SUBROUTINE add_boundary_flows()
    !Test script to add some velocity to the lower boundary
    IMPLICIT NONE
    real(num):: br, bl, kb    !Parameters for the boundary motions (taken from Pariat)
    real(num), dimension(:,:):: bzdy(0:nx+1,0:ny), bzdx(0:nx,0:ny+1)
    real(num), dimension(:,:):: bzdy0(0:nx,0:ny), bzdx0(0:nx,0:ny)
    real(num), dimension(:,:):: fact(0:nx+1,0:ny+1), fact0(0:nx,0:ny)

    if (z_down < 0) then
        br = 13.0_num; bl = 0.1_num; kb = 15.0_num

        bzdy = (bz(0:nx+1,1:ny+1,0) - bz(0:nx+1,0:ny,0)) / dy
        bzdx = (bz(1:nx+1,0:ny+1,0) - bz(0:nx,0:ny+1,0)) / dx

        bzdy0 = 0.5_num*(bzdy(1:nx+1,0:ny) + bzdy(0:nx,0:ny))
        bzdx0 = 0.5_num*(bzdx(0:nx,1:ny+1) + bzdx(0:nx,0:ny))

        fact = (kb*(br-bl))/(bz(0:nx+1,0:ny+1,0) + 1d-10)*tanh(kb*(bz(0:nx+1,0:ny+1,0)- bl)/(br-bl))
        fact0 = 0.25_num*(fact(0:nx,0:ny) + fact(1:nx+1,0:ny) + fact(0:nx,1:ny+1) + fact(1:nx+1, 1:ny+1))

        fact0(0,:) = 0.0_num; fact0(nx,:) = 0.0_num
        fact0(:,0) = 0.0_num; fact0(:,ny) = 0.0_num
        vx_surf(0:nx, 0:ny) = -fact0*bzdy0
        vy_surf(0:nx, 0:ny) = fact0*bzdx0

        vx(:,:,0) = vx(:,:,0) + shearfact*vx_surf
        vy(:,:,0) = vy(:,:,0) + shearfact*vy_surf
    end if

END SUBROUTINE add_boundary_flows


SUBROUTINE check_solenoidal()
    !Checks the solenoidal condition by calculating the divergence of all raw grid cells
    IMPLICIT NONE

    real(num), dimension(:,:,:):: div(0:nx+1,0:ny+1,0:nz+1)
    div = 0.0_num
    div(0:nx+1,0:ny+1,0:nz+1) = div(0:nx+1,0:ny+1,0:nz+1) + dx*dy*(bz(0:nx+1,0:ny+1,-1:nz) - bz(0:nx+1,0:ny+1,0:nz+1))
    div(0:nx+1,0:ny+1,0:nz+1) = div(0:nx+1,0:ny+1,0:nz+1) + dy*dz*(bx(-1:nx,0:ny+1,0:nz+1) - bx(0:nx+1,0:ny+1,0:nz+1))
    div(0:nx+1,0:ny+1,0:nz+1) = div(0:nx+1,0:ny+1,0:nz+1) + dx*dz*(by(0:nx+1,-1:ny,0:nz+1) - by(0:nx+1,0:ny+1,0:nz+1))

    print*, 'Max divergence', maxval(abs(div(0:nx+1,0:ny+1,0:nz+1)))

END SUBROUTINE check_solenoidal



SUBROUTINE calculate_electric()

    !Calculates the electric field - resistivity, magnetofriction and boundary effects
    IMPLICIT NONE

    ex = 0.0; ey = 0.0; ez = 0.0

    if (eta > 0) then
        !Determine the current from the magnetic field (after boundary conditions etc.)
        ex(1:nx, 0:ny,0:nz) = ex(1:nx, 0:ny,0:nz) + eta*jx(1:nx, 0:ny,0:nz)
        ey(0:nx, 1:ny,0:nz) = ey(0:nx, 1:ny,0:nz) + eta*jy(0:nx, 1:ny,0:nz)
        ez(0:nx, 0:ny,1:nz) = ez(0:nx, 0:ny,1:nz) + eta*jz(0:nx, 0:ny,1:nz)
    end if

    if (nu0 > 0) then  !There is some magnetofriction. Use existing velocity field to calculate electric field.
        ex1 = vz*by1 - vy*bz1
        ey1 = vx*bz1 - vz*bx1
        ez1 = vy*bx1 - vx*by1

        !Average to Ribs (interior only):
        ex(1:nx,0:ny,0:nz) = ex(1:nx,0:ny,0:nz)  + 0.5_num*(ex1(0:nx-1,0:ny,0:nz) + ex1(1:nx,0:ny,0:nz))
        ey(0:nx,1:ny,0:nz) = ey(0:nx,1:ny,0:nz)  + 0.5_num*(ey1(0:nx,0:ny-1,0:nz) + ey1(0:nx,1:ny,0:nz))
        ez(0:nx,0:ny,1:nz) = ez(0:nx,0:ny,1:nz)  + 0.5_num*(ez1(0:nx,0:ny,0:nz-1) + ez1(0:nx,0:ny,1:nz))

    end if




END SUBROUTINE calculate_electric























!*******************************************************************************
END MODULE evolve
!*******************************************************************************
