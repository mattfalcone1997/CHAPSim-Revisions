      SUBROUTINE BCAST_INI
     !==========================================================================
     !Broadcasts values from the ini file to all ranks
     !==========================================================================
      use init_info
      use mesh_info

      IMPLICIT NONE

      INTEGER(4)    :: POS0
      CHARACTER(1),ALLOCATABLE ::  buf(:)
      INTEGER(4)    :: NBUF

      NBUF = 14*4 + 18*8
      NBUF = 2*NBUF
      ALLOCATE (BUF(NBUF))
      POS0 = 0

      IF(MYID.EQ.0) THEN
         CALL mpi_pack(NREAD,   1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(NREAD_io,1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(NFLOW,   1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)

         CALL mpi_pack(N_START,1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(MULTIM, 1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(iswitch,1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(FLOWTP, 1, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)

         CALL mpi_pack(HX,     1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(HZ,     1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)

         CALL mpi_pack(REN,    1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(REINI,  1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(CFLGV,  1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)

         CALL mpi_pack(DT,     1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(VPER,   1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(SVPER,  1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)

         CALL mpi_pack(TSCN,   1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(TRST,   1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(TRST_io,1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(TSTOP,  1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(tbody,  1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)

         CALL mpi_pack(TTECCK, 1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(TSAVE1, 1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(TSTAV1, 1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)

         CALL mpi_pack(TLGRE,  1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(RATEM1, 1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(CFGV,   1, MPI_DOUBLE_PRECISION, buf, NBUF, pos0, ICOMM,IERROR)

         CALL mpi_pack(BCX,   2, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(BCZ,   2, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)
         CALL mpi_pack(BCX_io,2, MPI_INTEGER4, buf, NBUF, pos0, ICOMM,IERROR)

         CALL mpi_pack(IOFLOWflg,1, MPI_LOGICAL, buf, NBUF, pos0, ICOMM,IERROR)

         CALL mpi_bcast(buf, NBUF, MPI_PACKED, 0, ICOMM,IERROR)

      ELSE
         CALL mpi_bcast(buf, NBUF, MPI_PACKED, 0, ICOMM,IERROR)

         CALL mpi_unpack(buf, NBUF, pos0, NREAD,  1, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, NREAD_io,  1, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, NFLOW,  1, MPI_INTEGER4, ICOMM,IERROR)

         CALL mpi_unpack(buf, NBUF, pos0, N_START,1, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, MULTIM, 1, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, iswitch,1, MPI_INTEGER4, ICOMM,IERROR)

         CALL mpi_unpack(buf, NBUF, pos0, FLOWTP,1, MPI_INTEGER4, ICOMM,IERROR)

         CALL mpi_unpack(buf, NBUF, pos0, HX,     1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, HZ,     1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)

         CALL mpi_unpack(buf, NBUF, pos0, REN,    1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, REINI,  1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, CFLGV,   1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)

         CALL mpi_unpack(buf, NBUF, pos0, DT,     1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, VPER,   1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, SVPER,  1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)

         CALL mpi_unpack(buf, NBUF, pos0, TSCN,   1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, TRST,   1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, TRST_io,   1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, TSTOP,  1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, tbody,  1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)

         CALL mpi_unpack(buf, NBUF, pos0, TTECCK,  1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, TSAVE1,  1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, TSTAV1,  1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)

         CALL mpi_unpack(buf, NBUF, pos0, TLGRE,   1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, RATEM1,  1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, CFGV,  1, MPI_DOUBLE_PRECISION, ICOMM,IERROR)


         CALL mpi_unpack(buf, NBUF, pos0, BCX,  2, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, BCZ,  2, MPI_INTEGER4, ICOMM,IERROR)
         CALL mpi_unpack(buf, NBUF, pos0, BCX_io,  2, MPI_INTEGER4, ICOMM,IERROR)

         CALL mpi_unpack(buf, NBUF, pos0, IOFLOWflg,  1, MPI_LOGICAL, ICOMM,IERROR)

      END IF

      DEALLOCATE(BUF)

      RETURN

      END SUBROUTINE BCAST_INI

!***************************************************************************************************************
