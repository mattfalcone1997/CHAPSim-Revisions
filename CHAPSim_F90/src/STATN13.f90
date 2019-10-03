      SUBROUTINE STATN13(LO,LQ)
      use postprocess_info
      use mesh_info
      use flow_info
      IMPLICIT NONE

      INTEGER(4),INTENT(IN) :: LO
      INTEGER(4),INTENT(IN) :: LQ
      
      REAL(WP) :: UUJJ
      INTEGER(4) :: I, J, K, JJ
      REAL(WP),ALLOCATABLE :: UU(:)
      REAL(WP),ALLOCATABLE :: WORK_UU(:)
      
      ALLOCATE ( UU(NCL2)      )
      ALLOCATE ( WORK_UU(NCL2) )
      UU      = 0.0_WP
      WORK_UU = 0.0_WP
      
      DO J=1,N2DO(MYID)
         JJ=JCL2G(J)
         UUJJ=0.0_WP
         DO I=1,NCL1
            DO K=1,NCL3
               UUJJ=UUJJ+(Qtmp(I,J,K)-STA13(1,LQ,JJ))**LO
            ENDDO
         ENDDO
         UU(JJ)=UUJJ*VL1313
      ENDDO
      
      CALL MPI_BARRIER(ICOMM,IERROR)
      CALL MPI_ALLREDUCE(UU(1), WORK_UU(1), NCL2, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, ICOMM, IERROR)
     
      IF (MYID.EQ.0) THEN
         DO JJ=1,NCL2
            STA13(LO,LQ,JJ)=WORK_UU(JJ)
         ENDDO
      ENDIF
      
      DEALLOCATE ( UU      )
      DEALLOCATE ( WORK_UU )
      
      RETURN
      END
