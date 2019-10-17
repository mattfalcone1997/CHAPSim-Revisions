       SUBROUTINE start_mpi
       use mpi_info
       IMPLICIT NONE

       INTEGER(4)  :: NDIM11
       INTEGER(4)  :: DIMS(1)
       INTEGER(4)  :: ICOORDS(1)
       INTEGER(4)  :: NMINUS
       INTEGER(4)  :: NPLUS
       INTEGER(4)  :: IDIMS
       LOGICAL        PERIODS(1),REORDER

       CALL MPI_INIT(IERROR) !Interface with MPI module, starting it with the number
       if(IERROR==1) then
        call ERRHDL('mpi_init fails!',0) !If MPI fails write string to output file
       end if

       CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERROR)  !MYID gives the process its rank
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,SIZE,IERROR)  !Gives the nu,ber of processes

       NDIM11=1 !Mpich docs states this is the number of dimensions of cartesian grid - 1
       DIMS(1)=SIZE !DIMS(1) has the value of the number of processes

       PERIODS(1)=.TRUE.
       REORDER=.TRUE.

       !This section seems to indicate how arrays are split and ranks reordered
       !although this needs confirmation
       CALL MPI_CART_CREATE(MPI_COMM_WORLD,NDIM11,DIMS,PERIODS,  &
                             REORDER,ICOMM,IERROR)
      !Reordered Mpi ranks 'distributed'
       CALL MPI_COMM_RANK(ICOMM,MYID,IERROR)
       !Shifts source ranks to destination ranks - upwards by 1; NPLUS is the
       !rank of destination process
       CALL MPI_CART_SHIFT(ICOMM,0,1,NMINUS,NPLUS,IERROR)
       CALL MPI_CART_GET(ICOMM,1,IDIMS,PERIODS,ICOORDS,IERROR)

       NPTOT = SIZE
       NPSLV = SIZE-1 !Limit of the rank size - first process is rank 0 hence last is SIZE-1

       RETURN
       end subroutine start_mpi
