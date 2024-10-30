!*******************************************************************************

MODULE pressure
!*******************************************************************************
! Functions for the new pressure part
! Gives a function to add on to the velocity term, which is averaged to grid POINTS

!*******************************************************************************
    USE shared_data

    IMPLICIT NONE

!*******************************************************************************
CONTAINS

SUBROUTINE pressure_function()

    !Input the pressure function f(y) here - note y is the vertical coordinate like LARE
    !Should this be positive or negative?!?!
    implicit none
    integer:: k

    do k = 0, nz+1
        if (zstar < 0.01) then
            if (k == 0) print*, 'Base case, no pressure function'
            fz(:,:,k) = 0.0_num
        else
            if (decay_type == 0) then
                if (k == 0 .and. proc_num == 0) print*, 'Base case, no pressure function'
                    fz(:,:,k) = 0.0_num
            else if (decay_type == 1) then  !Exponential
                if (k == 0 .and. proc_num == 0) print*, 'Exponential Pressure Function'
                    fz(:,:,k) = a*exp(-zs(k)/b)
            else if (decay_type == 2) then !Smooth tanh
                if (k == 0 .and. proc_num == 0) print*, 'Smooth Tanh Pressure Function'
                    fz(:,:,k) = a*(1.0_num - b*tanh((zs(k)-zstar)/deltaz))
            else !Sharp tanh
                if (k == 0 .and. proc_num == 0) print*, 'Sharp Tanh Pressure Function'
                    fz(:,:,k) = a*(1.0_num - b*tanh((zs(k)-zstar)/deltaz))
            end if
        end if
    end do

END SUBROUTINE pressure_function

SUBROUTINE calculate_jp()
    !Calculates jp - the 'pressure current'
    implicit none

    jpx(0:nx+1, 0:ny,0:nz) = (fz(0:nx+1,1:ny+1,0:nz)*bz(0:nx+1,1:ny+1,0:nz) - fz(0:nx+1,0:ny,0:nz)*bz(0:nx+1, 0:ny,0:nz))/dy

    jpy(0:nx, 0:ny+1,0:nz) = -(fz(1:nx+1,0:ny+1,0:nz)*bz(1:nx+1,0:ny+1,0:nz) - fz(0:nx,0:ny+1,0:nz)*bz(0:nx,0:ny+1,0:nz))/dx

    !jpz1(0:nx, 0:ny) =  (fy(1:nx+1,0:ny)*by(1:nx+1,0:ny) - fy(0:nx,0:ny)*by(0:nx, 0:ny))/dx

END SUBROUTINE calculate_jp

SUBROUTINE calculate_pressure()
    !Does the cross product with b, averages to gridpoints and does the softening as for the velocity
    implicit none

    !REAL(num), DIMENSION(0:nx,0:ny):: gx2, gy2

    !Can get nu directly from the velocity calculation, so don't need to do this twice. But this function must come after calculate_velocity
    !gx2(0:nx,0:ny) = -nu0*(jpz1(0:nx,0:ny)*by1(0:nx,0:ny))/nu(0:nx,0:ny)
    !gy2(0:nx,0:ny) = nu0*(jpz1(0:nx,0:ny)*bx1(0:nx,0:ny))/nu(0:nx,0:ny)

    !vx(0:nx,0:ny) = vx(0:nx,0:ny) + gx2(0:nx, 0:ny)
    !vy(0:nx,0:ny) = vy(0:nx,0:ny) + gy2(0:nx, 0:ny)

    !vx(:,0) = 0.0; vy(:,0) = 0.0

END SUBROUTINE calculate_pressure





!*******************************************************************************
END MODULE pressure
!*******************************************************************************
