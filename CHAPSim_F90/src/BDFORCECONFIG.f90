      SUBROUTINE BDFORCECONFIG
   
      use flow_info
      use mesh_info
      IMPLICIT NONE
      
      
      INTEGER(4)  :: J, JJ
      INTEGER(4)  :: I, IP
      INTEGER(4)  :: K
      
      REAL(WP)     :: UUJJ
      REAL(WP)     :: VELENTER
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
            IP=IPV(I)
            VELENTER = 0.0_WP
            DO K=1,NCL3
               VELENTER= ( Q(IP,J,K,1) + Q(I,J,K,1) ) * 0.50_WP
               UUJJ=VELENTER+UUJJ
            ENDDO
         ENDDO
         UU(JJ)=UUJJ*VL1313
      ENDDO
      
      CALL MPI_BARRIER(ICOMM,IERROR)
      CALL MPI_ALLREDUCE(UU(1), WORK_UU(1), NCL2, MPI_DOUBLE_PRECISION, &
            MPI_SUM, ICOMM, IERROR)
     
      IF (MYID.EQ.0) THEN
          DO JJ=1,NCL2
              VMNO0(JJ)=WORK_UU(JJ)
          ENDDO
      ENDIF
      
      CALL MPI_BCAST( VMNO0, NCL2, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR )
      
      DEALLOCATE ( UU      )
      DEALLOCATE ( WORK_UU )
      DEALLOCATE ( VMNO0   )
    
      RETURN
      END
