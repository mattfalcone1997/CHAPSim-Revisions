
      SUBROUTINE MOMFA3_Y_tg(IDR)
      use mesh_info
      use flow_info
      IMPLICIT NONE      
      
      INTEGER(4),INTENT(IN) :: IDR
      
      INTEGER(4)  :: I, J, K 
      REAL(WP)     :: FJ  (NCL1,NND2) 
      REAL(WP)     :: AMJV(NCL1,NND2)
      REAL(WP)     :: ACJV(NCL1,NND2)
      REAL(WP)     :: APJV(NCL1,NND2)
      REAL(WP)     :: BCJ (NCL1,2)
      INTEGER(4)  :: NF 

      IF(IDR.EQ.2) THEN
         NF = NND2
      ELSE
         NF = NCL2 
      END IF 
      FJ  = 0.0_WP
      AMJV = 0.0_WP
      APJV = 0.0_WP
      ACJV = 0.0_WP
      BCJ = 0.0_WP
      F   = 0.0_WP
      CALL TRASP23L2G_RHS
       
      DO K=1,N3DO(MYID)  
         DO I=1,NCL1 
            DO J=1,NF
               FJ(I,J) = F(I,J,K)
               ACJV(I,J)= 1.0_WP-COE*ACVR(J,IDR)
               APJV(I,J)=-COE*APVR(J,IDR)
               AMJV(I,J)=-COE*AMVR(J,IDR)
           ENDDO
           BCJ(I,:) = 0.0_WP
         ENDDO
         
         CALL TDMAIJJ_nonCYC (AMJV(1:NCL1,1:NF),&
                              ACJV(1:NCL1,1:NF),&
                              APJV(1:NCL1,1:NF),&
                              FJ (1:NCL1,1:NF),&
                              BCJ(1:NCL1,1:2),&
                              1,NF,1,NCL1)
         
         DO I=1,NCL1
            DO J=1,NF
               F(I,J,K)=FJ(I,J)
            ENDDO
         ENDDO
      ENDDO
      
      CALL TRASP23G2L_RHS
      
      RETURN
      END SUBROUTINE MOMFA3_Y_tg

      SUBROUTINE MOMFA3_Y_io(IDR)
      use mesh_info
      use flow_info
      IMPLICIT NONE      
      
      INTEGER(4),INTENT(IN) :: IDR
      
      INTEGER(4)  :: I, J, K 
      REAL(WP)     :: FJ  (NCL1_io,NND2) 
      REAL(WP)     :: AMJV(NCL1_io,NND2)
      REAL(WP)     :: ACJV(NCL1_io,NND2)
      REAL(WP)     :: APJV(NCL1_io,NND2)
      REAL(WP)     :: BCJ (NCL1_io,2)
      INTEGER(4)  :: NF 

      IF(IDR.EQ.2) THEN
         NF = NND2
      ELSE
         NF = NCL2 
      END IF 
      FJ  = 0.0_WP
      AMJV = 0.0_WP
      APJV = 0.0_WP
      ACJV = 0.0_WP
      BCJ = 0.0_WP
      F   = 0.0_WP
 
      CALL TRASP23L2G_RHS_io
       
      DO K=1,N3DO(MYID)       
         DO I=1,NCL1_io 
            DO J=1,NF
               FJ(I,J) = F_io(I,J,K)
               ACJV(I,J)= 1.0_WP-COE*ACVR(J,IDR)
               APJV(I,J)=-COE*APVR(J,IDR)
               AMJV(I,J)=-COE*AMVR(J,IDR)
           ENDDO
           BCJ(I,:) = 0.0_WP
         ENDDO
         
         CALL TDMAIJJ_nonCYC (AMJV(1:NCL1_io,1:NF),&
                              ACJV(1:NCL1_io,1:NF),&
                              APJV(1:NCL1_io,1:NF),&
                              FJ (1:NCL1_io,1:NF),&
                              BCJ(1:NCL1_io,1:2),&
                              1,NF,1,NCL1_io)
         
         DO I=1,NCL1_io
            DO J=1,NF
               F_io(I,J,K)=FJ(I,J)
            ENDDO
         ENDDO
      ENDDO
      
      CALL TRASP23G2L_RHS_io
      
      RETURN
      END SUBROUTINE MOMFA3_Y_io

