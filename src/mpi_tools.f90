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

        CALL mpi_set_error_handler

        ! - Establish grid decomposition into processes, which should be powers of 2
        proc_divide = nprocs
        x_procs = 1; y_procs = 1; z_procs = 1
        dcount = maxloc((/nx_global,ny_global,nz_global/),dim=1) - 1
        do while (.true.)
            if (mod(proc_divide, 2) .ne. 0) then
                exit
            else
                if (dcount == 0) x_procs = x_procs*2
                if (dcount == 1) y_procs = y_procs*2
                if (dcount == 2) z_procs = z_procs*2
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

  SUBROUTINE mpi_set_error_handler

    INTEGER :: errhandler

    CALL MPI_COMM_CREATE_ERRHANDLER(mpi_error_handler, errhandler, ierr)
    CALL MPI_COMM_SET_ERRHANDLER(MPI_COMM_WORLD, errhandler, ierr)

  END SUBROUTINE mpi_set_error_handler



  SUBROUTINE mpi_error_handler(comm, error_code)

    INTEGER :: comm, error_code
    REAL :: tmp1, tmp2
    CHARACTER(LEN=29) :: errstring(0:MPI_ERR_LASTCODE)

    errstring(MPI_SUCCESS                  ) = 'MPI_SUCCESS                  '
    errstring(MPI_ERR_BUFFER               ) = 'MPI_ERR_BUFFER               '
    errstring(MPI_ERR_COUNT                ) = 'MPI_ERR_COUNT                '
    errstring(MPI_ERR_TYPE                 ) = 'MPI_ERR_TYPE                 '
    errstring(MPI_ERR_TAG                  ) = 'MPI_ERR_TAG                  '
    errstring(MPI_ERR_COMM                 ) = 'MPI_ERR_COMM                 '
    errstring(MPI_ERR_RANK                 ) = 'MPI_ERR_RANK                 '
    errstring(MPI_ERR_REQUEST              ) = 'MPI_ERR_REQUEST              '
    errstring(MPI_ERR_ROOT                 ) = 'MPI_ERR_ROOT                 '
    errstring(MPI_ERR_GROUP                ) = 'MPI_ERR_GROUP                '
    errstring(MPI_ERR_OP                   ) = 'MPI_ERR_OP                   '
    errstring(MPI_ERR_TOPOLOGY             ) = 'MPI_ERR_TOPOLOGY             '
    errstring(MPI_ERR_DIMS                 ) = 'MPI_ERR_DIMS                 '
    errstring(MPI_ERR_ARG                  ) = 'MPI_ERR_ARG                  '
    errstring(MPI_ERR_UNKNOWN              ) = 'MPI_ERR_UNKNOWN              '
    errstring(MPI_ERR_TRUNCATE             ) = 'MPI_ERR_TRUNCATE             '
    errstring(MPI_ERR_OTHER                ) = 'MPI_ERR_OTHER                '
    errstring(MPI_ERR_INTERN               ) = 'MPI_ERR_INTERN               '
    errstring(MPI_ERR_IN_STATUS            ) = 'MPI_ERR_IN_STATUS            '
    errstring(MPI_ERR_PENDING              ) = 'MPI_ERR_PENDING              '
    errstring(MPI_ERR_ACCESS               ) = 'MPI_ERR_ACCESS               '
    errstring(MPI_ERR_AMODE                ) = 'MPI_ERR_AMODE                '
    errstring(MPI_ERR_ASSERT               ) = 'MPI_ERR_ASSERT               '
    errstring(MPI_ERR_BAD_FILE             ) = 'MPI_ERR_BAD_FILE             '
    errstring(MPI_ERR_BASE                 ) = 'MPI_ERR_BASE                 '
    errstring(MPI_ERR_CONVERSION           ) = 'MPI_ERR_CONVERSION           '
    errstring(MPI_ERR_DISP                 ) = 'MPI_ERR_DISP                 '
    errstring(MPI_ERR_DUP_DATAREP          ) = 'MPI_ERR_DUP_DATAREP          '
    errstring(MPI_ERR_FILE_EXISTS          ) = 'MPI_ERR_FILE_EXISTS          '
    errstring(MPI_ERR_FILE_IN_USE          ) = 'MPI_ERR_FILE_IN_USE          '
    errstring(MPI_ERR_FILE                 ) = 'MPI_ERR_FILE                 '
    errstring(MPI_ERR_INFO_KEY             ) = 'MPI_ERR_INFO_KEY             '
    errstring(MPI_ERR_INFO_NOKEY           ) = 'MPI_ERR_INFO_NOKEY           '
    errstring(MPI_ERR_INFO_VALUE           ) = 'MPI_ERR_INFO_VALUE           '
    errstring(MPI_ERR_INFO                 ) = 'MPI_ERR_INFO                 '
    errstring(MPI_ERR_IO                   ) = 'MPI_ERR_IO                   '
    errstring(MPI_ERR_KEYVAL               ) = 'MPI_ERR_KEYVAL               '
    errstring(MPI_ERR_LOCKTYPE             ) = 'MPI_ERR_LOCKTYPE             '
    errstring(MPI_ERR_NAME                 ) = 'MPI_ERR_NAME                 '
    errstring(MPI_ERR_NO_MEM               ) = 'MPI_ERR_NO_MEM               '
    errstring(MPI_ERR_NOT_SAME             ) = 'MPI_ERR_NOT_SAME             '
    errstring(MPI_ERR_NO_SPACE             ) = 'MPI_ERR_NO_SPACE             '
    errstring(MPI_ERR_NO_SUCH_FILE         ) = 'MPI_ERR_NO_SUCH_FILE         '
    errstring(MPI_ERR_PORT                 ) = 'MPI_ERR_PORT                 '
    errstring(MPI_ERR_QUOTA                ) = 'MPI_ERR_QUOTA                '
    errstring(MPI_ERR_READ_ONLY            ) = 'MPI_ERR_READ_ONLY            '
    errstring(MPI_ERR_RMA_CONFLICT         ) = 'MPI_ERR_RMA_CONFLICT         '
    errstring(MPI_ERR_RMA_SYNC             ) = 'MPI_ERR_RMA_SYNC             '
    errstring(MPI_ERR_SERVICE              ) = 'MPI_ERR_SERVICE              '
    errstring(MPI_ERR_SIZE                 ) = 'MPI_ERR_SIZE                 '
    errstring(MPI_ERR_SPAWN                ) = 'MPI_ERR_SPAWN                '
    errstring(MPI_ERR_UNSUPPORTED_DATAREP  ) = 'MPI_ERR_UNSUPPORTED_DATAREP  '
    errstring(MPI_ERR_UNSUPPORTED_OPERATION) = 'MPI_ERR_UNSUPPORTED_OPERATION'
    errstring(MPI_ERR_WIN                  ) = 'MPI_ERR_WIN                  '
    errstring(MPI_ERR_LASTCODE             ) = 'MPI_ERR_LASTCODE             '

    PRINT*, 'Caught MPI error: ', TRIM(errstring(error_code))
    IF (comm == MPI_COMM_WORLD) THEN
      PRINT*, 'Communicator MPI_COMM_WORLD'
    ELSE
      PRINT*, 'Communicator ', comm, '(Not MPI_COMM_WORLD)'
    END IF

    ! Deliberately raise a divide-by-zero error
    tmp1 = 0.0
    tmp2 = 1.0 / tmp1

  END SUBROUTINE mpi_error_handler

  END MODULE mpi_tools
