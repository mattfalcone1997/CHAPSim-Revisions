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

       CALL MPI_INIT(IERROR)
       if(IERROR==1) then
        call ERRHDL('mpi_init fails!',0)
       end if
       
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERROR)  
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,SIZE,IERROR)
       
       NDIM11=1
       DIMS(1)=SIZE

       PERIODS(1)=.TRUE.   
       REORDER=.TRUE.    
!    
       CALL MPI_CART_CREATE(MPI_COMM_WORLD,NDIM11,DIMS,PERIODS,  &
                             REORDER,ICOMM,IERROR)

       CALL MPI_COMM_RANK(ICOMM,MYID,IERROR)
       CALL MPI_CART_SHIFT(ICOMM,0,1,NMINUS,NPLUS,IERROR)
       CALL MPI_CART_GET(ICOMM,1,IDIMS,PERIODS,ICOORDS,IERROR)
       
       NPTOT = SIZE
       NPSLV = SIZE-1 
       
       RETURN
       end subroutine start_mpi
