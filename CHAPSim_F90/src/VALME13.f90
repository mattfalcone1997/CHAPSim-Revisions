      SUBROUTINE VALME13(LO,LQ)
      
      use mesh_info
      use init_info
      use flow_info
      use postprocess_info
      IMPLICIT NONE

      INTEGER(4),INTENT(IN) :: LO
      INTEGER(4),INTENT(IN) :: LQ
      
      REAL(WP)   :: UUJJ
      INTEGER(4) :: I, J, K, JJ   
      REAL(WP),ALLOCATABLE :: UU(:)
      REAL(WP),ALLOCATABLE :: WORK_UU(:)
      REAL(WP),ALLOCATABLE :: VMNO0(:)

      ALLOCATE ( UU(NCL2)      )
      ALLOCATE ( WORK_UU(NCL2) )
      ALLOCATE ( VMNO0(NCL2)   )
      UU      = 0.0_WP
      WORK_UU = 0.0_WP
      VMNO0   = 0.0_WP
      DO J=1,N2DO(MYID)
         JJ=JCL2G(J)
         UUJJ=0.0_WP
         DO I=1,NCL1
            DO K=1,NCL3
               UUJJ=Qtmp(I,J,K)+UUJJ
            ENDDO
        ENDDO
        UU(JJ)=UUJJ*VL1313
      ENDDO
   
!      if (iswitch.eq.2.and.myid.eq.0) then
!         DO J=1,1
!            JJ=JCL2G(J)
!            UUJJ=0.0_WP
!            DO I=1,NCL1
!               DO K=1,NCL3
!                  UUJJ=UUJJ+0.50_WP*(Qtmp(I,J,K)+Qtmp(I,J,isym(K)))
!               ENDDO
!            ENDDO
!            UU(NND2)=UUJJ*VL1313
!         ENDDO     
!      endif
       
      CALL MPI_BARRIER(ICOMM,IERROR)
      CALL MPI_ALLREDUCE(UU(1), WORK_UU(1), NCL2, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, ICOMM, IERROR)
     
      IF (MYID.EQ.0) THEN
          DO JJ=1,NCL2
             VMNO0(JJ)=WORK_UU(JJ)
          ENDDO
      ENDIF
      
!      CALL BRO_REAL( 0, VMNO0 ,M2,M2)
      CALL MPI_BCAST( VMNO0, NCL2, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR )
      
      DO JJ=1,NCL2
         STA13(LO,LQ,JJ)=VMNO0(JJ)
      ENDDO
      
      DEALLOCATE ( UU      )
      DEALLOCATE ( WORK_UU )
      DEALLOCATE ( VMNO0   )
      
      RETURN
      END
