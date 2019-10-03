      SUBROUTINE CFL_tg
      use flow_info
      use init_info
      use mesh_info
      
      IMPLICIT NONE
      
      REAL(WP)      :: CFLMA
      REAL(WP)      :: CFLM_WORK
      REAL(WP)      :: QCF 
      INTEGER(4)  :: I, IP
      INTEGER(4)  :: J, JP, JJ
      INTEGER(4)  :: K, KP

      CFLMA=0.0_WP
      CFLM_WORK = 0.0_WP
      QCF = 0.0_WP
      
      DO K=1,NCL3
         KP=KPV(K)
         DO J=1,N2DO(MYID)
            JP=J+1
            JJ=JCL2G(J)
            DO I=1,NCL1
               IP=IPV(I)
               QCF=(DABS((Q(I,J,K,1)+Q(IP,J,K,1))*DXI)+        &
                    DABS((Q(I,J,K,2)+Q(I,JP,K,2))*DYCI(JJ))+        &
                    DABS((Q(I,J,K,3)+Q(I,J,KP,3))*DZI/rm(jj)))*0.50_WP  
               CFLMA=DMAX1(CFLMA,QCF)
            ENDDO
         ENDDO
      ENDDO
      
      CALL MPI_BARRIER(ICOMM,IERROR)
      CALL MPI_ALLREDUCE(CFLMA, CFLM_WORK, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, ICOMM, IERROR)
      CFLMM = CFLM_WORK
      
      RETURN
      END
      
!***********************************************************************

      SUBROUTINE CFL_io
      use flow_info
      use init_info
      use mesh_info
      
      IMPLICIT NONE
      
      REAL(WP)      :: CFLMA
      REAL(WP)      :: CFLM_WORK
      REAL(WP)      :: QCF 
      INTEGER(4)  :: I, IP
      INTEGER(4)  :: J, JP, JJ
      INTEGER(4)  :: K, KP

      CFLMA=0.0_WP
      CFLM_WORK = 0.0_WP
      QCF = 0.0_WP
      
      DO K=1,NCL3
         KP=KPV(K)
         DO J=1,N2DO(MYID)
            JP=J+1
            JJ=JCL2G(J)
            DO I=1,NCL1_io
               IP=IPV_io(I)
               QCF=(DABS((Q_io(I,J,K,1)+Q_io(IP,J,K,1))*DXI)+        &
                    DABS((Q_io(I,J,K,2)+Q_io(I,JP,K,2))*DYCI(JJ))+        &
                    DABS((Q_io(I,J,K,3)+Q_io(I,J,KP,3))*DZI/rm(jj)))*0.50_WP  
               CFLMA=DMAX1(CFLMA,QCF)
            ENDDO
         ENDDO
      ENDDO
      
      CALL MPI_BARRIER(ICOMM,IERROR)
      CALL MPI_ALLREDUCE(CFLMA, CFLM_WORK, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, ICOMM, IERROR)
      CFLMM_io = CFLM_WORK
      
      RETURN
      END

