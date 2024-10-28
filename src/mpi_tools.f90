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
        integer:: i
        INTEGER:: mpi_dims(3), mpi_dims_best(3)
        LOGICAL:: mpi_periodic(3)
        REAL(num):: diff_best, diff, mean

        ! - Have already established the global rank.
        call mpi_init(ierr)  !Tells it to start using MPI

        call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr) !Number of processes globally.
        call mpi_comm_rank(MPI_COMM_WORLD, proc_num, ierr) !Returns the rank of current process

        ! Choose optimum division of procs that fits grid dimensions:
        diff_best = REAL(nprocs)
        mpi_dims_best = (/1,1,1/)
        DO i = 2, nprocs, 2
            IF (MOD(nprocs, i) == 0) THEN
                mpi_dims = (/0, 0, i/)
                ! Find optimum decomposition with i points in p:
                CALL MPI_DIMS_CREATE(nprocs, 3, mpi_dims, ierr)
                ! Check whether this is allowed:
                nx = nx_global / mpi_dims(1)
                ny = ny_global / mpi_dims(2)
                nz = nz_global / mpi_dims(3)
                IF (((nx * mpi_dims(1)) == nx_global) &
                    .AND. ((ny * mpi_dims(2) == ny_global)) &
                    .AND. ((nz * mpi_dims(3)) == nz_global)) THEN
                    mean = SUM(REAL(mpi_dims))/3.0_d
                    diff = MAXVAL(REAL(mpi_dims) - mean)
                    IF (diff < diff_best) THEN
                        diff_best = diff
                        mpi_dims_best = mpi_dims
                    END IF
                END IF
            END IF
        END DO
        IF (mpi_dims_best(1) * mpi_dims_best(2) * mpi_dims_best(3) == nprocs) THEN
            mpi_dims = mpi_dims_best
            nx = nx_global / mpi_dims(1)
            ny = ny_global / mpi_dims(2)
            nz = nz_global / mpi_dims(3)
        ELSE
            PRINT*,'ERROR: THIS NUMBER OF MPI PROCS DOES NOT FIT THE GRID'
            CALL MPI_abort(comm, ierr)
        END IF

        x_procs = mpi_dims(1); y_procs = mpi_dims(2); z_procs = mpi_dims(3);
        MPI_periodic = (/.false.,.false.,.false./)
        !Attempt to use the DUMFRIC way of establishing the communicator:
        CALL MPI_CART_CREATE(MPI_COMM_WORLD, 3, MPI_dims, MPI_periodic, .TRUE., &
        comm, ierr)

        !Redistribute ranks based on this new communicator
        CALL MPI_COMM_RANK(comm, proc_num, ierr)

        CALL MPI_CART_COORDS(comm, proc_num, 3, mpi_loc, ierr)
        CALL MPI_CART_SHIFT(comm, 0, 1, x_down, x_up, ierr)
        CALL MPI_CART_SHIFT(comm, 1, 1, y_down, y_up, ierr)
        CALL MPI_CART_SHIFT(comm, 2, 1, z_down, z_up, ierr)

        x_rank = mpi_loc(1); y_rank = mpi_loc(2); z_rank = mpi_loc(3)

        call MPI_BARRIER(comm, ierr)

        do i = 0, nprocs-1
            if (proc_num == i .and. .false.) then
                print*, proc_num, x_rank, y_rank, z_rank
                print*, 'x', x_down, x_up
                print*, 'y', y_down, y_up
                print*, 'z', z_down, z_up
                print*, '______________________________'
            end if
            call MPI_BARRIER(comm, ierr)
        end do

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
