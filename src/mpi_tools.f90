!*******************************************************************************
MODULE mpi_tools
!*******************************************************************************
! Main fortran file. Use this to initialise parameters, grid, produce/read in initial conditions and then run the code. Theoretically. Want to run with intel debuggers on, otherwise definitely not going to pick up on all the mistakes.
!*******************************************************************************
    USE shared_data

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE start_mpi()

        IMPLICIT NONE
        integer:: proc_divide, dcount
        ! - Have already established the global rank.
        call mpi_init(ierr)  !Tells it to start using MPI

        call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr) !Number of processes globally.
        call mpi_comm_rank(MPI_COMM_WORLD, proc_num, ierr) !Returns the rank of current process

        ! - Establish grid decomposition into processes, which should be powers of 2
        proc_divide = nprocs
        x_procs = 1; y_procs = 1; z_procs = 1
        dcount = maxloc((/nx_global,ny_global,nz_global/),dim=1) - 1
        do while (.true.)
            if (mod(proc_divide, 2) .ne. 0) then
                exit
            else
                if (dcount == 2) x_procs = x_procs*2
                if (dcount == 1) y_procs = y_procs*2
                if (dcount == 0) z_procs = z_procs*2
                proc_divide = proc_divide/2
                dcount = mod(dcount + 1, 3)
            end if
        end do

        if (nprocs .ne. x_procs*y_procs*z_procs) then
            print*, 'MPI Decomposition failed'
            CALL MPI_abort(comm, ierr)
        end if

        z_rank = (z_procs*proc_num)/nprocs
        y_rank = mod(proc_num, x_procs*y_procs)/x_procs
        x_rank = mod(mod(proc_num, x_procs*y_procs), x_procs)

        nx = nx_global/x_procs; ny = ny_global/y_procs; nz = nz_global/z_procs
        t = 0.0

        x_up = -1; x_down = -1
        y_up = -1; y_down = -1
        z_up = -1; z_down = -1

        if (x_rank < x_procs-1) x_up = proc_num + 1
        if (x_rank > 0)      x_down = proc_num - 1


        if (y_rank < y_procs-1) y_up = proc_num + x_procs
        if (y_rank > 0) y_down = proc_num - x_procs

        if (z_rank < z_procs-1) z_up = proc_num + (x_procs*y_procs)
        if (z_rank > 0) z_down = proc_num - (x_procs*y_procs)

        !print*, proc_num, y_rank, y_procs, y_up, y_down

        return
    END SUBROUTINE start_mpi

    SUBROUTINE bfield_mpi

        IMPLICIT NONE

        integer:: j,k

        !MPI routines for passing the boundary data
        if (.true.) then
        !Send z data DOWN
        if (z_down >= 0) then
            do j = -1 ,ny+1
                call mpi_send(by(0:nx+1,j,1), nx+2, MPI_DOUBLE_PRECISION, z_down, j+1, comm, ierr)
            end do

            do j = 0 ,ny+1
                call mpi_send(bx(-1:nx+1,j,1), nx+3, MPI_DOUBLE_PRECISION, z_down, j + nx*ny + 1, comm, ierr)
            end do
        end if
        !Send z data UP
        if (z_up >= 0) then
            do j = -1 ,ny+1
                call mpi_send(by(0:nx+1,j,nz), nx+2, MPI_DOUBLE_PRECISION, z_up, j + 1, comm, ierr)
            end do
            do j = 0 ,ny+1
                call mpi_send(bx(-1:nx+1,j,nz), nx+3, MPI_DOUBLE_PRECISION, z_up, j + nx*ny + 1, comm, ierr)
            end do
        end if

        call MPI_BARRIER(comm, ierr)
        !Receive z data from ABOVE
        if (z_up >= 0) then
            do j = -1 ,ny+1
                call mpi_recv(by(0:nx+1,j,nz+1), nx+2, MPI_DOUBLE_PRECISION, z_up, j + 1, comm, MPI_STATUS_IGNORE, ierr)
            end do
            do j = 0 ,ny+1
                call mpi_recv(bx(-1:nx+1,j,nz+1), nx+3, MPI_DOUBLE_PRECISION, z_up, j + nx*ny + 1, comm, MPI_STATUS_IGNORE, ierr)
            end do
        end if
        !Receive z data from BELOW
        if (z_down >= 0) then
            do j = -1 ,ny+1
                call mpi_recv(by(0:nx+1,j,0), nx+2, MPI_DOUBLE_PRECISION, z_down, j + 1, comm, MPI_STATUS_IGNORE, ierr)
            end do
            do j = 0 ,ny+1
                call mpi_recv(bx(-1:nx+1,j,0), nx+3 , MPI_DOUBLE_PRECISION, z_down, j + nx*ny + 1, comm, MPI_STATUS_IGNORE, ierr)
            end do
        end if

        call MPI_BARRIER(comm, ierr)
        end if


        if (.true.) then
        !Send x data DOWN
        if (x_down >= 0) then
            do j = 0, ny+1
            do k = -1 ,nz+1
                call mpi_send(bz(1,j,k), 1, MPI_DOUBLE_PRECISION, x_down, k+2*ny*(j+1) + 1, comm, ierr)
            end do
            end do
            do j = -1, ny+1
            do k = 0 ,nz+1
                call mpi_send(by(1,j,k), 1, MPI_DOUBLE_PRECISION, x_down, k+2*ny*(j+1) + ny*nz + 1, comm, ierr)
            end do
            end do
        end if

        !Send x data UP

        if (x_up >= 0) then
            do j = 0, ny+1
            do k = -1 ,nz+1
                call mpi_send(bz(nx,j,k), 1, MPI_DOUBLE_PRECISION, x_up, k+2*ny*(j+1) + 1, comm, ierr)
            end do
            end do
            do j = -1, ny+1
            do k = 0 ,nz+1
                call mpi_send(by(nx,j,k), 1, MPI_DOUBLE_PRECISION, x_up,  k+2*ny*(j+1) + ny*nz + 1, comm, ierr)
            end do
            end do
        end if


        call MPI_BARRIER(comm, ierr)

        !Receive x from above
        if (x_up >= 0) then
            do j = 0, ny+1
            do k = -1 ,nz+1
                call mpi_recv(bz(nx+1,j,k), 1, MPI_DOUBLE_PRECISION, x_up, k+2*ny*(j+1) + 1, comm, MPI_STATUS_IGNORE,ierr)
            end do
            end do
            do j = -1, ny+1
            do k = 0 ,nz+1
                call mpi_recv(by(nx+1,j,k), 1, MPI_DOUBLE_PRECISION, x_up,  k+2*ny*(j+1) + ny*nz + 1, comm, MPI_STATUS_IGNORE, ierr)
            end do
            end do
        end if

        !Receive x from below

        if (x_down >= 0) then
            do j = 0, ny+1
            do k = -1 ,nz+1
                call mpi_recv(bz(0,j,k), 1, MPI_DOUBLE_PRECISION, x_down, k+2*ny*(j+1) + 1, comm, MPI_STATUS_IGNORE, ierr)
            end do
            end do
            do j = -1, ny+1
            do k = 0 ,nz+1
                call mpi_recv(by(0,j,k), 1, MPI_DOUBLE_PRECISION, x_down, k+2*ny*(j+1) + ny*nz + 1, comm, MPI_STATUS_IGNORE, ierr)
            end do
            end do
        end if

        call MPI_BARRIER(comm, ierr)

        end if

        if (.true.) then
        !Send y data DOWN
        if (y_down >= 0) then
            do k = -1 ,nz+1
                call mpi_send(bz(0:nx+1,1,k), nx+2, MPI_DOUBLE_PRECISION, y_down, k + 1, comm, ierr)
            end do
            do k = 0 ,nz+1
                call mpi_send(bx(-1:nx+1,1,k), nx+3, MPI_DOUBLE_PRECISION, y_down, k + nx*nz + 1, comm, ierr)
            end do
        end if

        !Send y data UP
        if (y_up >= 0) then
            do k = -1 ,nz+1
                call mpi_send(bz(0:nx+1,ny,k), nx+2, MPI_DOUBLE_PRECISION, y_up, k + 1, comm, ierr)
            end do
            do k = 0 ,nz+1
                call mpi_send(bx(-1:nx+1,ny,k), nx+3, MPI_DOUBLE_PRECISION, y_up, k + nx*nz + 1, comm, ierr)
            end do
        end if

        call MPI_BARRIER(comm, ierr)

        !Receive y from above
        if (y_up >= 0) then
            do k = -1 ,nz+1
                call mpi_recv(bz(0:nx+1,ny+1,k), nx+2, MPI_DOUBLE_PRECISION, y_up, k + 1, comm, MPI_STATUS_IGNORE,ierr)
            end do
            do k = 0 ,nz+1
                call mpi_recv(bx(-1:nx+1,ny+1,k), nx+3, MPI_DOUBLE_PRECISION, y_up, k + nx*nz + 1, comm,MPI_STATUS_IGNORE, ierr)
            end do
        end if

        !Receive y from below
        if (y_down >= 0) then
            do k = -1 ,nz+1
                call mpi_recv(bz(0:nx+1,0,k), nx+2, MPI_DOUBLE_PRECISION, y_down, k + 1, comm, MPI_STATUS_IGNORE,ierr)
            end do
            do k = 0 ,nz+1
                call mpi_recv(bx(-1:nx+1,0,k), nx+3, MPI_DOUBLE_PRECISION, y_down, k + nx*nz + 1, comm,MPI_STATUS_IGNORE, ierr)
            end do
        end if

        call MPI_BARRIER(comm, ierr)

        end if
        return
    END SUBROUTINE bfield_mpi

END MODULE mpi_tools





















