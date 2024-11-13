!*******************************************************************************
PROGRAM main
!*******************************************************************************
! Main fortran file. Use this to initialise parameters, grid, produce/read in initial conditions and then run the code
!*******************************************************************************
    USE shared_data
    USE init
    USE mpi_tools
    USE evolve
    USE output
    !USE output, ONLY: save_snap, print_array, diagnostics
    IMPLICIT NONE

    ! Put some of the major variables in here - things that can be changed occasionally but not in a series of runs
    cfl  = 0.1
    mf_delta = 1e-4

    ! Import the parameters and set up the grid
    CALL initialise()

    if (.true.) then
    if (hamilton_flag < 0.5) then
        data_directory_root = '/extra/tmp/trcn27/mf3d/'
    else
        data_directory_root = '/nobackup/trcn27/mf3d0/'
    end if

    if (run_number < 10) then
        write (data_directory, "(A23,A2,I1,A1)") data_directory_root, '00',int(run_number), "/"
    else if (run_number < 100) then
        write (data_directory, "(A23,A1,I2,A1)") data_directory_root, '0', int(run_number), "/"
    else
        write (data_directory, "(A23,I3,A1)") data_directory_root, int(run_number), "/"
    end if

    if (proc_num == 0) print*, 'Initial condition set up in Fortran. Running...'
    do n = 0, nt-1  ! Actually run the code

        CALL timestep()  !Does everything except the actual timestep (for diagnostic reasons)

        if (MOD(n, (nt/int(nplots-1))) == 0) then   ! Save a snapshot (prints a message as well)
            CALL save_snap(int(n/(nt/int(nplots-1))))
            !if (proc_num == 0) print*, 'Snapshot ', int(n/(nt/int(nplots-1))), 'saved at time', t
        end if

        if (MOD(n, (nt/int(ndiags-1))) == 0) then   ! Save a snapshot (prints a message as well)
            CALL diagnostics(int(n/(nt/(ndiags-1))))
            !if (proc_num == 0) print*, 'Step', n, 'at time', t

            !print*, 'Max all currents', maxval(abs(jx(0:nx+1, 0:ny,0:nz))), maxval(abs(jy(0:nx, 0:ny+1,0:nz))), maxval(abs(jz(0:nx, 0:ny,0:nz+1)))

        end if

        ax = ax - dt*ex
        ay = ay - dt*ey
        az = az - dt*ez

        t = t + dt

    end do

    if (proc_num == 0) then
        print*, 'Open Flux', proc_num, z_rank, sum(abs(bz(1:nx,1:ny,nz)))
        print*, 'Max. currents', proc_num, sum(abs(jx(2:nx-2,2:ny-2,2:nz-1))), &
        sum(abs(jy(2:nx-2,2:ny-2,2:nz-1))), sum(abs(jz(2:nx-2,2:ny-2,2:nz-1)))
    end if
    !CALL diagnostics(int(n/(nt/(ndiags-1))))
    if (proc_num == 0) print*, 'Step', n, 'at time', t

    CALL save_snap(int(nplots-1.0))
    !CALL diagnostics(ndiags-1)
    end if
    if (proc_num == 0) print*, 'Fortran code completed sucessfully. Carry on.'
    CALL mpi_finalize(ierr)
    STOP

END PROGRAM main
