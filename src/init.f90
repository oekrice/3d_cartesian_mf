!*******************************************************************************
MODULE init
!*******************************************************************************
! Initialise the simulation. Read in/compute the initial condition (python code for that I think?), read in the parameters and all will be well
!*******************************************************************************
    USE shared_data
    USE mpi_tools
    USE netcdf
    USE pressure

    !USE output
    IMPLICIT NONE

!*******************************************************************************

CONTAINS

SUBROUTINE initialise()

    CALL read_parameters()

    CALL start_mpi()

    CALL MPI_BARRIER(comm, ierr)
    if (proc_num == -1) then
        print*, 'Run Number', int(run_number)
        print*, 'Outflow', voutfact
        print*, 'Shearing', shearfact
        print*, 'eta', eta
        print*, 'nu0', nu0
        print*, 'eta0', eta0
        print*, 'Data directory', data_directory
    end if

    CALL MPI_create_types

    CALL allocate_arrays()

    CALL establish_grid()

    CALL calculate_timestep()

    CALL set_outflow()

    CALL set_shearing()

    CALL pressure_function()

END SUBROUTINE initialise

SUBROUTINE read_parameters()
    !Reads in the parameters from the text file variables.txt
    REAL(num), DIMENSION(30):: variables
    CHARACTER(LEN=64):: input_value

    CHARACTER(LEN=64):: parameter_filename
    call get_command_argument(1, input_value)
    read(unit=input_value,fmt=*) run_number

    if (run_number < 10) then
      write (parameter_filename, "(A22, A2, I1, A4)") './parameters/variables', '00', int(run_number), '.txt'
    else if (run_number < 100) then
      write (parameter_filename, "(A22, A1, I2, A4)") './parameters/variables', '0', int(run_number), '.txt'
    else
      write (parameter_filename, "(A22, I3, A4)") './parameters/variables', int(run_number), '.txt'
    end if

    OPEN(1, FILE = parameter_filename)
    READ(1, *) variables
    CLOSE(1)

    run_number = variables(1)
    nx_global = int(variables(2))
    ny_global = int(variables(3))
    nz_global = int(variables(4))
    tmax = variables(5)

    nplots = int(variables(6))
    ndiags = int(variables(7))

    voutfact = variables(8)
    shearfact = variables(9)
    eta = variables(10)
    nu0 = variables(11)
    eta0 = variables(12)

    x0_global = variables(13)
    x1_global = variables(14)
    y0_global = variables(15)
    y1_global = variables(16)
    z0_global = variables(17)
    z1_global = variables(18)

    a = variables(19)
    b = variables(20)
    deltaz = variables(21)
    zstar = variables(22)

    hamilton_flag = int(variables(23))
    decay_type = int(variables(24))

END SUBROUTINE read_parameters

SUBROUTINE allocate_arrays()
    ! - Estalish the grid (not particularly difficult for a cartesian box). Try to use the same numbering system as LARE.
    ! - index 0 is always the cell on the boundary, or just OUTSIDE it. c cells are numbered from -1, s cells from -2
    ! - Arrays are local unless otherwise declared. Not sure if any global arrays are required, but perhaps

    !x0 = x0_global + x_rank*(x1_global-x0_global)/x_procs
    !x1 = x0_global + (x_rank + 1)*(x1_global-x0_global)/x_procs

    !y0 = y0_global + (y_rank)*(y1_global-y0_global)/y_procs
    !y1 = y0_global + (y_rank + 1)*(y1_global-y0_global)/y_procs

    !z0 = z0_global + z_rank*(z1_global-z0_global)/z_procs
    !z1 = z0_global + (z_rank + 1)*(z1_global-z0_global)/z_procs

    allocate(xs(-1:nx+1)); allocate(ys(-1:ny+1)); allocate(zs(-1:nz+1))
    allocate(xc(0:nx+1)); allocate(yc(0:ny+1)); allocate(zc(0:nz+1))

    !Allocate global arrays (things that could conceivably be plotted afterwards, I suppose)
    allocate(ax(0:nx+1,-1:ny+1,-1:nz+1)); allocate(ay(-1:nx+1,0:ny+1,-1:nz+1)); allocate(az(-1:nx+1,-1:ny+1,0:nz+1))
    allocate(bx(-1:nx+1,0:ny+1,0:nz+1)); allocate(by(0:nx+1,-1:ny+1,0:nz+1)); allocate(bz(0:nx+1,0:ny+1,-1:nz+1))
    allocate(jx(0:nx+1,-1:ny+1,-1:nz+1)); allocate(jy(-1:nx+1,0:ny+1,-1:nz+1)); allocate(jz(-1:nx+1,-1:ny+1,0:nz+1))
    allocate(ex(0:nx+1,-1:ny+1,-1:nz+1)); allocate(ey(-1:nx+1,0:ny+1,-1:nz+1)); allocate(ez(-1:nx+1,-1:ny+1,0:nz+1))

    !Fields averaged to gridpoints
    allocate(jx1(0:nx,0:ny,0:nz)); allocate(jy1(0:nx,0:ny,0:nz)); allocate(jz1(0:nx,0:ny,0:nz))
    allocate(bx1(0:nx,0:ny,0:nz)); allocate(by1(0:nx,0:ny,0:nz)); allocate(bz1(0:nx,0:ny,0:nz))
    allocate(ex1(0:nx,0:ny,0:nz)); allocate(ey1(0:nx,0:ny,0:nz)); allocate(ez1(0:nx,0:ny,0:nz))

    allocate(vx(0:nx,0:ny,0:nz)); allocate(vy(0:nx,0:ny,0:nz)); allocate(vz(0:nx,0:ny,0:nz))
    allocate(vx_surf(0:nx,0:ny)); allocate(vy_surf(0:nx,0:ny))

    allocate(lf(0:nx,0:ny,0:nz)); allocate(b2(0:nx,0:ny,0:nz)); allocate(nu(0:nx,0:ny,0:nz))
    allocate(soft(0:nx,0:ny,0:nz))

    allocate(bz_surf_reference(0:nx+1,0:ny+1))

    !Pressure arrays
    allocate(fz(0:nx+1,0:ny+1,-1:nz+1)) !Aligned with bz (I'm pretty sure)
    allocate(jpx(0:nx+1,-1:ny+1,-1:nz+1)); allocate(jpy(-1:nx+1,0:ny+1,-1:nz+1))
    allocate(jpx1(0:nx,0:ny,0:nz)); allocate(jpy1(0:nx,0:ny,0:nz))


    bx = 0.0_num; by = 0.0_num; bz = 0.0_num

END SUBROUTINE allocate_arrays

SUBROUTINE establish_grid()

    ! Imports the initial condition from the inits file.
    ! Will try to be smart with this, such that the vector potential is only read in in individual processes

    CHARACTER(LEN =64):: init_filename
    INTEGER:: ncid, vid

    if (run_number < 10) then
      write (init_filename, "(A12, A2, I1, A3)") './inits/init', '00', int(run_number), '.nc'
    else if (run_number < 100) then
      write (init_filename, "(A12, A1, I2, A3)") './inits/init', '0', int(run_number), '.nc'
    else
      write (init_filename, "(A12, I3, A3)") './inits/init', int(run_number), '.nc'
    end if

    call try(nf90_open(trim(init_filename), nf90_nowrite, ncid))

    call try(nf90_inq_varid(ncid, 'ax', vid))
    call try(nf90_get_var(ncid, vid, ax(1:nx,0:ny,0:nz), &
    start = (/x_rank*nx+1,y_rank*ny+1,z_rank*nz+1/),count = (/nx,ny+1,nz+1/)))

    call try(nf90_inq_varid(ncid, 'ay', vid))
    call try(nf90_get_var(ncid, vid, ay(0:nx,1:ny,0:nz), &
    start = (/x_rank*nx+1,y_rank*ny+1,z_rank*nz+1/),count = (/nx+1,ny,nz+1/)))

    call try(nf90_inq_varid(ncid, 'az', vid))
    call try(nf90_get_var(ncid, vid, az(0:nx,0:ny,1:nz), &
    start = (/x_rank*nx+1,y_rank*ny+1,z_rank*nz+1/),count = (/nx+1,ny+1,nz/)))

    call try(nf90_inq_varid(ncid, 'xs', vid))
    call try(nf90_get_var(ncid, vid, xs(0:nx), &
    start = (/x_rank*nx+1/),count = (/nx+1/)))

    call try(nf90_inq_varid(ncid, 'ys', vid))
    call try(nf90_get_var(ncid, vid, ys(0:ny), &
    start = (/y_rank*ny+1/),count = (/ny+1/)))

    call try(nf90_inq_varid(ncid, 'zs', vid))
    call try(nf90_get_var(ncid, vid, zs(0:nz), &
    start = (/z_rank*nz+1/),count = (/nz+1/)))

    call try(nf90_close(ncid))

    xs(-1) = 2*xs(0) - xs(1); xs(nx+1) = 2*xs(nx) - xs(nx-1)
    ys(-1) = 2*ys(0) - ys(1); ys(ny+1) = 2*ys(ny) - ys(ny-1)
    zs(-1) = 2*zs(0) - zs(1); zs(nz+1) = 2*zs(nz) - zs(nz-1)

    dx = sum(xs(0:nx+1) - xs(-1:nx))/(nx+2)
    dy = sum(ys(0:ny+1) - ys(-1:ny))/(ny+2)
    dz = sum(zs(0:nz+1) - zs(-1:nz))/(nz+2)

    xc = 0.5_num*(xs(-1:nx) + xs(0:nx+1))
    yc = 0.5_num*(ys(-1:ny) + ys(0:ny+1))
    zc = 0.5_num*(zs(-1:nz) + zs(0:nz+1))

    x0 = xs(0); x1 = xs(nx)
    y0 = ys(0); y1 = ys(ny)
    z0 = zs(0); z1 = zs(nz)

    volume = (x1-x0)*(y1-y0)*(z1-z0)
    volume_global = (x1_global-x0_global)*(y1_global-y0_global)*(z1_global-z0_global)

END SUBROUTINE establish_grid

subroutine try(status)
    ! Catch error in reading netcdf fild.
    INTEGER, INTENT(IN):: status

    if (status /= NF90_noerr) THEN
        PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
        call mpi_abort(comm, ierr)
    end if

end subroutine try

SUBROUTINE calculate_timestep()
    ! Calculates the timestep based on the parameters and cfl condition. Then adjusts so there are an integer number of steps in between snapshots, as that's neater
    REAL(num):: dt_ideal, nt_ft
    REAL(num):: plot_dt

    dt_ideal = 1e10 !Maximum time, will get smaller based on the conditions

    if (proc_num == 0) then
    if (nu0 > 0) dt_ideal = (min(dx,dy,dz))**2/(nu0*(1.0))
    if (nu0 > 0) print*, 'dt due to nu0', cfl*(min(dx,dy,dz))**2/(nu0*(1.0))
    if (voutfact > 0) dt_ideal = min(dt_ideal, dz/voutfact)
    if (voutfact > 0) print*, 'dt due to outflow',  cfl*dz/voutfact
    if (eta > 0 ) dt_ideal = min(dt_ideal, min(dx,dy,dz)**2/eta)
    if (eta > 0 ) print*, 'dt due to eta', cfl*min(dx,dy,dz)**2/eta
    if (shearfact > 0 ) dt_ideal = min(dt_ideal, min(dy,dx)/shearfact)
    if (shearfact > 0 ) print*, 'dt due to shearing', cfl*min(dy,dx)/(shearfact*15.0)
    if (eta0 > 0) dt_ideal = min(dt_ideal, min(dx,dy,dz)**2/eta0)
    if (eta0 > 0) print*, 'dt due to eta0',  cfl*min(dx,dy,dz)**2/eta0
    dt_ideal = cfl*dt_ideal
    print*, 'Ideal dt', dt_ideal

    !Adjust so there are an integer number of timesteps in between plots
    plot_dt = tmax/(ndiags-1)
    nt = int(ndiags-1)*(int((plot_dt-1d-6)/dt_ideal)+1)
    nt_ft = dble(nt)
    dt = tmax/float(nt)
    print*, 'Final dt', dt, ', total timesteps', nt, ', ', int(nt/(nplots-1)), 'per snapshot'
    end if

    call MPI_Bcast(dt, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(nt_ft, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    nt = int(nt_ft)

    call MPI_Barrier(comm,ierr)  !Wait for t to be broadcast everywhere.

END SUBROUTINE calculate_timestep

SUBROUTINE set_outflow()
    ! Set the outflow arrays vout and voutc, the sizes of which are declared in shared variables
    ! vouts is on gridpoints (like az), voutc is on x ribs (like ax). To allow for upwinding
    ! Need to take care with all this...
    INTEGER:: i, j, k
    REAL(num):: hfact

    allocate(voutx(0:nx+1,-1:ny+1,-1:nz+1))  !Outflow combining with ex
    allocate(vouty(-1:nx+1,0:ny+1,-1:nz+1))  !Outflow combining with ey

    do k = -1, nz+1
        hfact = (zs(k) - z0_global)/(z1_global - z0_global)   !Distance up the domain.
        do i = 0, nx+1
            do j = -1, ny+1
                voutx(i,j,k) = voutfact*hfact**2
            end do
        end do
        do i = -1, nx+1
            do j = 0, ny+1
                vouty(i,j,k) = voutfact*hfact**2
            end do
        end do
    end do

END SUBROUTINE set_outflow

SUBROUTINE set_shearing()
    ! Set the 1D array for the shearing velocity on the lower boundary. Is added to ex so is aligned with the x grid centres
    INTEGER:: i
    allocate(vz0(-1:nx+1))

    do i = -1, nx+1
        vz0(i) = shearfact*sin(xs(i)*PI/x1_global)
    end do

END SUBROUTINE set_shearing

!*******************************************************************************
END MODULE init
!*******************************************************************************
