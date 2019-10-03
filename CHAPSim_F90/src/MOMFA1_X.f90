      SUBROUTINE MOMFA1_X_tg(N2I,IDR)
      use mesh_info
      use flow_info
      IMPLICIT NONE
                  
      INTEGER(4),INTENT(IN)    :: N2I
      INTEGER(4),INTENT(IN)    :: IDR
      
      REAL(WP)  :: AMIV(NCL1,N2DO(MYID))
      REAL(WP)  :: ACIV(NCL1,N2DO(MYID))
      REAL(WP)  :: APIV(NCL1,N2DO(MYID))
      REAL(WP)  :: FI  (NCL1,N2DO(MYID))
      
      INTEGER(4) :: I, J, K, JSZ
      
      AMIV = 0.0_WP
      ACIV = 0.0_WP
      APIV = 0.0_WP
      FI  = 0.0_WP
             
      DO K=1,NCL3
         
         DO J=N2I,N2DO(MYID)
            DO I=1, NCL1
                ACIV(I,J)= 1.0_WP-COE*(-2.0_WP*DXQI)             
                APIV(I,J)=-COE*DXQI      
                AMIV(I,J)=-COE*DXQI      
                 FI(I,J)= RHS(I,J,K)     
            END DO
         END DO
         JSZ = N2DO(MYID)-N2I+1
         CALL TDMAIJI_CYC(AMIV(1:NCL1,N2I:N2DO(MYID)),&
                          ACIV(1:NCL1,N2I:N2DO(MYID)),&
                          APIV(1:NCL1,N2I:N2DO(MYID)),&
                          FI(1:NCL1,N2I:N2DO(MYID)),&
                          1,NCL1,N2I,JSZ)
         DO J=N2I,N2DO(MYID)
             DO I=1,NCL1
                RHS(I,J,K) = FI(I,J)
             END DO
         END DO
         
         IF(IOFLOWflg) THEN
            DO J=N2I,N2DO(MYID)
               BC_U_SSTAR(J,K,IDR) = RHS(NCL1,J,K)
            END DO
         END IF
         
      ENDDO
      
      RETURN
      END SUBROUTINE MOMFA1_X_tg
      
      SUBROUTINE MOMFA1_X_io(N2I,IDR)
      use mesh_info
      use flow_info
      IMPLICIT NONE
                  
      INTEGER(4),INTENT(IN)    :: N2I
      INTEGER(4),INTENT(IN)    :: IDR
      
      REAL(WP)  :: AMIV(NCL1_io,N2DO(MYID))
      REAL(WP)  :: ACIV(NCL1_io,N2DO(MYID))
      REAL(WP)  :: APIV(NCL1_io,N2DO(MYID))
      REAL(WP)  :: FI  (NCL1_io,N2DO(MYID))
      REAL(WP)  :: BCI(N2I:N2DO(MYID),2)
      
      INTEGER(4) :: I, J, K, JSZ
              
      DO K=1,NCL3
      
         AMIV = 0.0_WP
         ACIV = 0.0_WP
         APIV = 0.0_WP
         FI   = 0.0_WP
         BCI  = 0.0_WP
         
         DO J=N2I,N2DO(MYID)
            DO I=1, NCL1_io
                ACIV(I,J)= 1.0_WP-COE*(-2.0_WP*DXQI) 
                APIV(I,J)=     -COE*DXQI     
                AMIV(I,J)=     -COE*DXQI     
                 FI(I,J)=       RHS_io(I,J,K)
            END DO
            BCI(J,1) = BC_U_SSTAR(J,K,IDR)
            BCI(J,2) = BC_TDMA(1, J, K, IDR)
         END DO
         JSZ = N2DO(MYID)-N2I+1
         CALL TDMAIJI_nonCYC ( &
               AMIV(1:NCL1_io,N2I:N2DO(MYID)),  &
               ACIV(1:NCL1_io,N2I:N2DO(MYID)),  &
               APIV(1:NCL1_io,N2I:N2DO(MYID)),  &
               FI  (1:NCL1_io,N2I:N2DO(MYID)),  &
               BCI (          N2I:N2DO(MYID) ,1:2),1,NCL1_io,N2I,JSZ )
                              
                              
         DO J=N2I,N2DO(MYID)
             DO I=1,NCL1_io
                RHS_io(I,J,K) = FI(I,J)
             END DO
         END DO
         
      ENDDO
      
      RETURN
      END SUBROUTINE MOMFA1_X_io

