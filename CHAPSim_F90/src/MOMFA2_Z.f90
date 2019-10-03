
      SUBROUTINE MOMFA2_Z_tg(N2I)
      use mesh_info
      use flow_info
      IMPLICIT NONE
                  
      INTEGER(4),INTENT(IN)    :: N2I
      
      REAL(WP)      :: AMKV(NCL3,N2DO(MYID))
      REAL(WP)      :: ACKV(NCL3,N2DO(MYID))
      REAL(WP)      :: APKV(NCL3,N2DO(MYID))
      REAL(WP)      :: FK  (NCL3,N2DO(MYID))
      
      REAL(WP)      :: RMC2(N2DO(MYID))
      
      INTEGER(4) :: I, J, K, JJ, JSZ
      
      AMKV = 0.0_WP
      ACKV = 0.0_WP
      APKV = 0.0_WP
      FK  = 0.0_WP
      
      RMC2= 0.0_WP
        
      DO J=N2I, N2DO(MYID)
          JJ = JCL2G(J)
          IF (n2i.eq.2) THEN
             RMC2(J)=1.0_WP/rc(jj)**2
          ELSE
             RMC2(J)=1.0_WP/rm(jj)**2
          END IF
      END DO

      DO I=1,NCL1
      
          DO K=1,NCL3
              DO J=N2I,N2DO(MYID)
                 ACKV(K,J)= 1.0_WP-COE*(-2.0_WP*DZQI)*RMC2(J)
                 APKV(K,J)=-COE*DZQI*RMC2(J)      
                 AMKV(K,J)=-COE*DZQI*RMC2(J)   
                 FK(K,J) = RHS(I,J,K)     
              ENDDO
         ENDDO

         JSZ = N2DO(MYID)-N2I+1
         CALL TDMAIJI_CYC(AMKV(1:NCL3,N2I:N2DO(MYID)),&
                          ACKV(1:NCL3,N2I:N2DO(MYID)),&
                          APKV(1:NCL3,N2I:N2DO(MYID)),&
                           FK(1:NCL3,N2I:N2DO(MYID)),&
                          1,NCL3,N2I,JSZ)
         DO K=1,NCL3
            DO J=N2I,N2DO(MYID)
               RHS(I,J,K)=FK(K,J)
            ENDDO
         ENDDO
         
      ENDDO
      
      RETURN
      END SUBROUTINE MOMFA2_Z_tg
      
      SUBROUTINE MOMFA2_Z_io(N2I)
      use mesh_info
      use flow_info
      IMPLICIT NONE
                  
      INTEGER(4),INTENT(IN)    :: N2I
      
      REAL(WP)      :: AMKV(NCL3,N2DO(MYID))
      REAL(WP)      :: ACKV(NCL3,N2DO(MYID))
      REAL(WP)      :: APKV(NCL3,N2DO(MYID))
      REAL(WP)      :: FK  (NCL3,N2DO(MYID))
      
      REAL(WP)      :: RMC2(N2DO(MYID))
      
      INTEGER(4) :: I, J, K, JJ, JSZ
      
      AMKV = 0.0_WP
      ACKV = 0.0_WP
      APKV = 0.0_WP
      FK  = 0.0_WP
      
      RMC2= 0.0_WP
         
      DO J=N2I, N2DO(MYID)
          JJ = JCL2G(J)
          IF (n2i.eq.2) THEN
             RMC2(J)=1.0_WP/rc(jj)**2
          ELSE
             RMC2(J)=1.0_WP/rm(jj)**2
          END IF
      END DO

      DO I=1,NCL1_io
      
          DO K=1,NCL3
              DO J=N2I,N2DO(MYID)
                 ACKV(K,J)= 1.0_WP-COE*(-2.0_WP*DZQI)*RMC2(J)
                 APKV(K,J)=-COE*DZQI*RMC2(J)       
                 AMKV(K,J)=-COE*DZQI*RMC2(J)    
                 FK(K,J) = RHS_io(I,J,K)     
              ENDDO
         ENDDO

         JSZ = N2DO(MYID)-N2I+1
         CALL TDMAIJI_CYC(AMKV(1:NCL3,N2I:N2DO(MYID)),&
                          ACKV(1:NCL3,N2I:N2DO(MYID)),&
                          APKV(1:NCL3,N2I:N2DO(MYID)),&
                           FK(1:NCL3,N2I:N2DO(MYID)),&
                          1,NCL3,N2I,JSZ)
         DO K=1,NCL3
            DO J=N2I,N2DO(MYID)
               RHS_io(I,J,K)=FK(K,J)
            ENDDO
         ENDDO
         
      ENDDO
      
      RETURN
      END SUBROUTINE MOMFA2_Z_io

