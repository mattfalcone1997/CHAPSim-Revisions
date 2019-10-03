      SUBROUTINE BCAST_CPARAM
      use mpi_info
      use cparam 
      
      IMPLICIT NONE

      INTEGER(4)    :: POS0
      CHARACTER(1),ALLOCATABLE ::  buf(:)
      INTEGER(4)    :: NBUF
      
      NBUF = 10*4
      NBUF = 2*NBUF
      ALLOCATE (BUF(NBUF))
    
      POS0 = 0
      
      IF(MYID.EQ.0) THEN
         CALL mpi_pack(NCL1,    1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(NCL1_io, 1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(NCL2,    1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(NCL3,    1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         
         CALL mpi_pack(NND1,    1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(NND1_io, 1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(NND2,    1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(NND3,    1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)

         CALL mpi_bcast(buf, NBUF, MPI_PACKED, 0, ICOMM,IERROR)
         
      ELSE
         CALL mpi_bcast(buf, NBUF, MPI_PACKED, 0, ICOMM,IERROR)
         
         CALL mpi_unpack(buf, NBUF, pos0, NCL1, 1, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, NCL1_io, 1, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, NCL2, 1, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, NCL3, 1, MPI_INTEGER4, ICOMM,IERROR)
         
         CALL mpi_unpack(buf, NBUF, pos0, NND1, 1, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, NND1_io, 1, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, NND2, 1, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, NND3, 1, MPI_INTEGER4, ICOMM,IERROR)

      END IF
      
      DEALLOCATE(BUF)

      RETURN
      
      END SUBROUTINE BCAST_CPARAM

