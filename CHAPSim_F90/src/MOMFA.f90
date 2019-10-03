 
      SUBROUTINE MOMFA(NS,IDR)
      use init_info
      use mesh_info
      use flow_info
      IMPLICIT NONE
 
      INTEGER(4),INTENT(IN)  :: NS
      INTEGER(4),INTENT(IN)  :: IDR
      
      INTEGER(4)  :: N2I
      INTEGER(4)  :: I, J, K
      
           
      COE=0.50_WP*TALP(NS)*DT*CVISC  
      N2I=1
      IF ((MYID.EQ.0).AND.(IDR.EQ.2)) THEN
          N2I=2
      ENDIF

      CALL MOMFA1_X_tg(N2I,IDR)
      CALL MOMFA2_Z_tg(N2I)
      CALL MOMFA3_Y_tg(IDR)
            
      DO K=1,NCL3
         DO J=N2I,N2DO(MYID)
            DO I=1,NCL1
               Q(I,J,K,IDR)=RHS(I,J,K)+Q(I,J,K,IDR)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
      
  
      SUBROUTINE MOMFA_io(NS,IDR)
      use init_info
      use mesh_info
      use flow_info
      IMPLICIT NONE
 
      INTEGER(4),INTENT(IN)  :: NS
      INTEGER(4),INTENT(IN)  :: IDR
      
      INTEGER(4)  :: N2I
      INTEGER(4)  :: I, J, K
      
          
      COE=0.50_WP*TALP(NS)*DT*CVISC
     
      N2I=1
      IF ((MYID.EQ.0).AND.(IDR.EQ.2)) THEN
          N2I=2
      ENDIF

      CALL MOMFA1_X_io(N2I,IDR)
      CALL MOMFA2_Z_io(N2I)
      CALL MOMFA3_Y_io(IDR)
         
      DO K=1,NCL3
         DO J=N2I,N2DO(MYID)
            DO I=1,NCL1_io
               Q_io(I,J,K,IDR)=RHS_io(I,J,K)+Q_io(I,J,K,IDR)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END

