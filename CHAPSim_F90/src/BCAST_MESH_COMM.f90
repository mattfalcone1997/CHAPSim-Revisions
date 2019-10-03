    SUBROUTINE BCAST_MESH_COMM
     
      use init_info
      use mesh_info
      
      IMPLICIT NONE

      INTEGER(4)    :: POS0
      CHARACTER(1),ALLOCATABLE ::  buf(:)
      INTEGER(4)    :: NBUF
      
      NBUF = 5*4 + 24*8
      NBUF = 2*NBUF
      ALLOCATE (BUF(NBUF))
      POS0 = 0
      
      IF(MYID.EQ.0) THEN
         CALL mpi_pack(NSST,   1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(IPWALL1,1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(IPWALL2,1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(JWLC1,  1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(JWGL1,  1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)         
         
         CALL mpi_pack(PI,     1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         
         CALL mpi_pack(HX,     1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(HZ,     1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(DX,     1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(DZ,     1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(ALX1,   1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(ALX2,   1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(ALX3,   1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         
         CALL mpi_pack(VL1313,    1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR) 
         CALL mpi_pack(VL1313_io, 1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(CVISC,     1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         
         
         CALL mpi_pack(TALP(0),4, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(TGAM(0),4, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(TROH(0),4, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(DXI,    1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(DZI,    1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(DXQI,   1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(DZQI,   1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(QDX1,   1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(QDX3,   1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         

         CALL mpi_bcast(buf, NBUF, MPI_PACKED, 0, ICOMM,IERROR)
         
      ELSE
         CALL mpi_bcast(buf, NBUF, MPI_PACKED, 0, ICOMM,IERROR)
         
         CALL mpi_unpack(buf, NBUF, pos0, NSST,   1, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, IPWALL1,1, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, IPWALL2,1, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, JWLC1,  1, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, JWGL1,  1, MPI_INTEGER4, ICOMM,IERROR)
         
         CALL mpi_unpack(buf, NBUF, pos0, PI,     1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         
         CALL mpi_unpack(buf, NBUF, pos0, HX,     1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, HZ,     1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, DX,     1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, DZ,     1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, ALX1,   1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, ALX2,   1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, ALX3,   1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         
         CALL mpi_unpack(buf, NBUF, pos0, VL1313, 1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, VL1313_io, 1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, CVISC,  1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         
         
         CALL mpi_unpack(buf, NBUF, pos0, TALP(0),4, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, TGAM(0),4, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, TROH(0),4, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, DXI,    1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, DZI,    1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, DXQI,   1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, DZQI,   1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, QDX1,   1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, QDX3,   1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
      
      END IF
      
      DEALLOCATE(BUF)
      
      RETURN
      
      END SUBROUTINE BCAST_MESH_COMM
