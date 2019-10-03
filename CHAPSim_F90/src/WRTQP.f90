
      SUBROUTINE WRTQP_tg
      
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      
      CHARACTER(7) :: PNTIM
      
      INTEGER(4) :: ITIM
      REAL(WP)    :: U_F0(NCL1,NCL2,  NCL3)
      REAL(WP)    :: V_F0(NCL1,NCL2+1,NCL3)
      REAL(WP)    :: W_F0(NCL1,NCL2,  NCL3)
      REAL(WP)    :: P_F0(NCL1,NCL2,  NCL3)
      
      INTEGER(4)  :: INN
      INTEGER(4)  :: INNN
      INTEGER(4)  :: I
      INTEGER(4)  :: J, JJ
      INTEGER(4)  :: K, KK
      INTEGER(4)  :: N2DOID
      REAL(WP)     :: QAUX (NCL1,0:N2DO(MYID)+1,NCL3,NDV,1:NPTOT)
      REAL(WP)     :: PRAUX(NCL1,0:N2DO(MYID)+1,NCL3,    1:NPTOT)
             
      INN=NCL1*(N2DO(MYID)+2)*NCL3*NDV
      CALL MPI_GATHER( Q, INN, MPI_DOUBLE_PRECISION, QAUX, INN,  &
            MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
      CALL MPI_BARRIER(ICOMM,IERROR)
     
      INNN=NCL1*(N2DO(MYID)+2)*NCL3
      CALL MPI_GATHER( PR, INNN, MPI_DOUBLE_PRECISION, PRAUX, INNN,  &
            MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
      CALL MPI_BARRIER(ICOMM,IERROR)
         
      IF (MYID.EQ.0) THEN
      
         DO KK=0,NPSLV
            N2DOID=JDEWT(KK)-JDSWT(KK)+1
            DO J=1,N2DOID
               JJ=JDSWT(KK)-1+J
               DO I=1,NCL1
                  DO K=1,NCL3
                     U_F0(I,JJ,K)=(QAUX(I,J,K,1,KK+1))
                     W_F0(I,JJ,K)=(QAUX(I,J,K,3,KK+1))
                     V_F0(I,JJ,K)=(QAUX(I,J,K,2,KK+1))
                     P_F0(I,JJ,K)=(PRAUX(I,J,K,KK+1))
                  ENDDO
               ENDDO
            ENDDO
         END DO
   
         DO I=1,NCL1
            DO K=1,NCL3
               V_F0(I,NND2,K)=(QAUX(I,N2DO(MYID)+1,K,2,NPTOT))  
            ENDDO
         ENDDO
         
         ITIM=INT(phyTIME_TG*MULTIM+0.3)
         WRITE(PNTIM,'(I7.7)') ITIM
!
    
         OPEN(10,FILE='DNS2_periodicxz_RU'//PNTIM//'U.D',FORM='UNFORMATTED')
         WRITE(10) NCL1,N2DO(MYID),NCL3,ITERG
         WRITE(10) phyTIME_TG,REN,DT
         WRITE(10) U_F0
         CLOSE(10)
!
         OPEN(10,FILE='DNS2_periodicxz_RU'//PNTIM//'W.D',FORM='UNFORMATTED')
         WRITE(10) NCL1,N2DO(MYID),NCL3,ITERG
         WRITE(10) phyTIME_TG,REN,DT
         WRITE(10) W_F0              
         CLOSE(10)
!
         OPEN(10,FILE='DNS2_periodicxz_RU'//PNTIM//'P.D',FORM='UNFORMATTED')
         WRITE(10) NCL1,N2DO(MYID),NCL3,ITERG
         WRITE(10) phyTIME_TG,REN,DT
         WRITE(10) P_F0
         CLOSE(10)
!
         OPEN(10,FILE='DNS2_periodicxz_RU'//PNTIM//'V.D',FORM='UNFORMATTED')
         WRITE(10) NCL1,N2DO(MYID),NCL3,ITERG
         WRITE(10) phyTIME_TG,REN,DT
         WRITE(10) V_F0
         CLOSE(10) 
         
      ENDIF
      
      RETURN
      END
      
      SUBROUTINE WRTQP_io
      
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      
      CHARACTER(7) :: PNTIM
      
      INTEGER(4) :: ITIM
      REAL(WP)    :: U_F0_io(0:NCL1_io+1,NCL2,  NCL3)
      REAL(WP)    :: V_F0_io(0:NCL1_io+1,NCL2+1,NCL3)
      REAL(WP)    :: W_F0_io(0:NCL1_io+1,NCL2,  NCL3)
      REAL(WP)    :: P_F0_io(0:NCL1_io+1,NCL2,  NCL3)
      
      INTEGER(4)  :: INN
      INTEGER(4)  :: INNN
      INTEGER(4)  :: I
      INTEGER(4)  :: J, JJ
      INTEGER(4)  :: K, KK
      INTEGER(4)  :: N2DOID
      REAL(WP)     :: QAUX_io (0:NCL1_io+1,0:N2DO(MYID)+1,NCL3,NDV,1:NPTOT)
      REAL(WP)     :: PRAUX_io(0:NCL1_io+1,0:N2DO(MYID)+1,NCL3,    1:NPTOT)
            
         INN=(NCL1_io+2)*(N2DO(MYID)+2)*NCL3*NDV
         CALL MPI_GATHER( Q_io, INN, MPI_DOUBLE_PRECISION, QAUX_io, INN,  &
               MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
         CALL MPI_BARRIER(ICOMM,IERROR)
    
         INNN=(NCL1_io+2)*(N2DO(MYID)+2)*NCL3
         CALL MPI_GATHER( PR_io, INNN, MPI_DOUBLE_PRECISION, PRAUX_io, INNN,  &
               MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
         CALL MPI_BARRIER(ICOMM,IERROR)
         
      IF (MYID.EQ.0) THEN
      
           DO KK=0,NPSLV
               N2DOID=JDEWT(KK)-JDSWT(KK)+1
               DO J=1,N2DOID
                  JJ=JDSWT(KK)-1+J
                  DO I=0,NCL1_io+1,1
                     DO K=1,NCL3
                        U_F0_io(I,JJ,K)= QAUX_io(I,J,K,1,KK+1)
                        W_F0_io(I,JJ,K)= QAUX_io(I,J,K,3,KK+1)
                        V_F0_io(I,JJ,K)= QAUX_io(I,J,K,2,KK+1)
                        P_F0_io(I,JJ,K)= PRAUX_io(I,J,K,KK+1)
                     ENDDO
                  ENDDO
               ENDDO
            END DO
   
            DO I=0,NCL1_io+1,1
               DO K=1,NCL3
                  V_F0_io(I,NND2,K)= QAUX_io(I,N2DO(MYID)+1,K,2,NPTOT) 
               ENDDO
            ENDDO

         ITIM=INT(phyTIME*MULTIM+0.3)
         WRITE(PNTIM,'(I7.7)') ITIM
!
    
         OPEN(10,FILE='DNS2_inoudomain_RU'//PNTIM//'U.D',FORM='UNFORMATTED')
         WRITE(10) NCL1_io,N2DO(MYID),NCL3,ITERG
         WRITE(10) phyTIME,REN,DT
         WRITE(10) U_F0_io
         CLOSE(10)
!
         OPEN(10,FILE='DNS2_inoudomain_RU'//PNTIM//'W.D',FORM='UNFORMATTED')
         WRITE(10) NCL1_io,N2DO(MYID),NCL3,ITERG
         WRITE(10) phyTIME,REN,DT
         WRITE(10) W_F0_io            
         CLOSE(10)
!
         OPEN(10,FILE='DNS2_inoudomain_RU'//PNTIM//'P.D',FORM='UNFORMATTED')
         WRITE(10) NCL1_io,N2DO(MYID),NCL3,ITERG
         WRITE(10) phyTIME,REN,DT
         WRITE(10) P_F0_io
         CLOSE(10)
!
         OPEN(10,FILE='DNS2_inoudomain_RU'//PNTIM//'V.D',FORM='UNFORMATTED')
         WRITE(10) NCL1_io,N2DO(MYID),NCL3,ITERG
         WRITE(10) phyTIME,REN,DT
         WRITE(10) V_F0_io
         CLOSE(10) 
         
      ENDIF
      
      RETURN
      END

