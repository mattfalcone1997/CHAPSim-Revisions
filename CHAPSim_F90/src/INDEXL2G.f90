
      SUBROUTINE INDEXL2G
      use mesh_info
      IMPLICIT NONE

      INTEGER(4) :: J, JJ, K
      INTEGER(4),ALLOCATABLE :: CIP(:)

      ALLOCATE ( CIP(NCL2) )
      CIP    = 0

      DO J=1,N2DO(MYID)
         JCL2G(J)=JDSWT(MYID)-1+J
      END DO

      DO K=1,N3DO(MYID)
         KCL2G(K)=KDSWT(MYID)-1+K
      END DO

      DO J=1,N2DO(MYID)
         JJ = JCL2G(J)
         CIP(JJ)=MYID
      END DO
       !CIP is an array which indicates which rank each index belongs to
       !MPI_ALLREDUCE creates a single array from all the ranks then broadcasts
       !to all ranks hence creating a copy of the entire array on all ranks
      CALL MPI_BARRIER(ICOMM,IERROR) !All processes wait until all complete
      CALL MPI_ALLREDUCE(CIP(1),  CLGIP(1),  NCL2, MPI_INTEGER4,&
                MPI_SUM, ICOMM, IERROR)

      DEALLOCATE (CIP)

      RETURN

      END SUBROUTINE INDEXL2G
