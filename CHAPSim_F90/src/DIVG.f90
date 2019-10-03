      SUBROUTINE DIVG_tg(NS)
      use init_info
      use mesh_info
      use flow_info
      IMPLICIT NONE

      INTEGER(4),INTENT(IN) :: NS
      
      INTEGER(4)   :: K, KP
      INTEGER(4)   :: J, JP, JJ
      INTEGER(4)   :: I, IP
      REAL(WP)       :: DVIGVELO
      REAL(WP)     :: COE1
       
      RHSLLPHI = 0.0_WP   
      COE1 = DT*TALP(NS)
      DO K=1,NCL3
         KP=KPV(K)
         DO J=1,N2DO(MYID)
            JP=J+1
            JJ=JCL2G(J)
            DO I=1,NCL1
               IP=IPV(I)
               DVIGVELO=(Q(IP,J,K,1)-Q(I,J,K,1))*DXI*rm(jj)**2     &
                       +(Q(I,JP,K,2)-Q(I,J,K,2))*DYFI(JJ)*rm(jj)   & 
                       +(Q(I,J,KP,3)-Q(I,J,K,3))*DZI
               RHSLLPHI(I,J,K)=DVIGVELO/COE1
            ENDDO
         ENDDO
      ENDDO
      
      RETURN
      
      END

      SUBROUTINE DIVG_io(NS)
      use init_info
      use mesh_info
      use flow_info
      IMPLICIT NONE

      INTEGER(4),INTENT(IN) :: NS
      
      INTEGER(4)   :: K, KP
      INTEGER(4)   :: J, JP, JJ
      INTEGER(4)   :: I, IP
      REAL(WP)     :: DVIGVELO
      REAL(WP)     :: COE1
       
      RHSLLPHI_io = 0.0_WP   
      RHSLLPHI_io_tmp = 0.0_WP
      COE1 = DT*TALP(NS)
      
      DO J=1,N2DO(MYID)
         JP=J+1
         JJ=JCL2G(J)
         DO K=1,NCL3
            KP=KPV(K)
            DO I=1,NCL1_io
               IP=IPV_io(I)
               DVIGVELO=(Q_io(IP,J,K,1)-Q_io(I,J,K,1))*DXI*rm(jj)**2     &
                       +(Q_io(I,JP,K,2)-Q_io(I,J,K,2))*DYFI(JJ)*rm(jj)   & 
                       +(Q_io(I,J,KP,3)-Q_io(I,J,K,3))*DZI
               RHSLLPHI_io(I,J,K)=DVIGVELO/COE1
               RHSLLPHI_io_tmp(I,J,K) = RHSLLPHI_io(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      
      RETURN
      
      END
