      MODULE TEC360_INFO
      USE WPRECISION
      USE cparam
      
      

      INTEGER(4),PARAMETER :: NI = 2
      INTEGER(4),PARAMETER :: NJ = 2
      INTEGER(4),PARAMETER :: NK = 1
      INTEGER(4) :: IID(NI)
      INTEGER(4) :: JID(NJ)
      INTEGER(4) :: KID(NK)
      INTEGER(4) :: NCOUNT(4) = 0

      REAL(WP),ALLOCATABLE     :: U_INTP(:, :, :)
      REAL(WP),ALLOCATABLE     :: V_INTP(:, :, :)
      REAL(WP),ALLOCATABLE     :: W_INTP(:, :, :)
      REAL(WP),ALLOCATABLE     :: P_INTP(:, :, :)
      
      REAL(WP),ALLOCATABLE     :: U_INTP_io(:, :, :)
      REAL(WP),ALLOCATABLE     :: V_INTP_io(:, :, :)
      REAL(WP),ALLOCATABLE     :: W_INTP_io(:, :, :)
      REAL(WP),ALLOCATABLE     :: P_INTP_io(:, :, :)
      
      
      END MODULE TEC360_INFO

      SUBROUTINE TEC360_NODE_EXTRAPOLATION
      use TEC360_INFO
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
      INTEGER(4)  :: I, IK,IP, IM
      INTEGER(4)  :: J, JJ,JP, JM
      INTEGER(4)  :: K, KK,KP, KM
      INTEGER(4)  :: N2DOID, NCL1TT
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
                     U_F0(I,JJ,K)= QAUX(I,J,K,1,KK+1)
                     W_F0(I,JJ,K)= QAUX(I,J,K,3,KK+1)
                     V_F0(I,JJ,K)= QAUX(I,J,K,2,KK+1)
                     P_F0(I,JJ,K)= PRAUX(I,J,K,KK+1)
                  ENDDO
               ENDDO
            ENDDO
         END DO
   
         DO I=1,NCL1
            DO K=1,NCL3
               V_F0(I,NND2,K)=(QAUX(I,N2DO(MYID)+1,K,2,NPTOT))  
            ENDDO
         ENDDO

         ! FOR U
         DO I=1,NCL1
            DO J=2,NCL2
               DO K=1,NCL3
                  JM = JMV(J)
                  KM = KMV(K)
                  U_INTP(I,J,K) = 0.250_WP * ( U_F0(I,JM,KM) + U_F0(I,JM,K) + U_F0(I,J,KM) + U_F0(I,J,K) )
               END DO
            END DO
         END DO
         ! I = NND1
         U_INTP(NND1,:,:) = U_INTP(1,:,:)
         ! K =1
         U_INTP(:,:,NND3) = U_INTP(:,:,1)
         ! J = 1 or NND2
         U_INTP(:,1,:) = 0.0_WP
         U_INTP(:,NND2,:) = 0.0_WP
         
         
         ! FOR V
         DO I=1,NCL1
            DO J=1,NND2
               DO K=1,NCL3
                  IM = IMV(I)
                  KM = KMV(K)
                  V_INTP(I,J,K) = 0.250_WP * ( V_F0(IM,J,KM) + V_F0(IM,J,K) + V_F0(I,J,KM) + V_F0(I,J,K) )
               END DO
            END DO
         END DO

         V_INTP(NND1,:,:) = V_INTP(1,:,:)
         V_INTP(:,:,NND3) = V_INTP(:,:,1)
         V_INTP(:,1,:) = 0.0_WP
         V_INTP(:,NND2,:) = 0.0_WP
         

         !FOR W
         DO I=1,NCL1
            DO J=2,NCL2
               DO K=1,NCL3
                  JM = JMV(J)
                  IM = IMV(I)
                  W_INTP(I,J,K) = 0.250_WP * ( W_F0(IM,JM,K) + W_F0(IM,J,K) + W_F0(I,JM,K) + W_F0(I,J,K) )
               END DO
            END DO
         END DO
         W_INTP(NND1,:,:) = W_INTP(1,:,:)
         W_INTP(:,:,NND3) = W_INTP(:,:,1)
         W_INTP(:,1,:) = 0.0_WP
         W_INTP(:,NND2,:) = 0.0_WP
         
         !FOR P
         DO I=1,NCL1
            DO J=2,NCL2
               DO K=1,NCL3
                  JM = JMV(J)
                  IM = IMV(I)
                  KM = KMV(K)
                  P_INTP(I,J,K) = 0.1250_WP * ( P_F0(IM,JM,KM) + P_F0(IM,JM,K) + &
                                              P_F0(IM,J ,KM) + P_F0(IM,J ,K) + &
                                              P_F0(I ,JM,KM) + P_F0(I ,JM,K) + &
                                              P_F0(I ,J ,KM) + P_F0(I ,J ,K) )
               END DO
            END DO
         END DO
         P_INTP(NND1,:,:) = P_INTP(1,:,:)
         P_INTP(:,:,NND3) = P_INTP(:,:,1)
         P_INTP(:,1,:) = P_INTP(:,2,:)
         P_INTP(:,NND2,:) = P_INTP(:,NCL2,:)
         
 
      ENDIF
      
      RETURN
      END SUBROUTINE TEC360_NODE_EXTRAPOLATION

!********************************************************************************************************      
      SUBROUTINE TEC360_NODE_EXTRAPOLATION_io
      use TEC360_INFO
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      
      CHARACTER(7) :: PNTIM
      
      INTEGER(4) :: ITIM
      INTEGER(4)  :: INN
      INTEGER(4)  :: INNN
      INTEGER(4)  :: I, IK,IP, IM
      INTEGER(4)  :: J, JJ,JP, JM
      INTEGER(4)  :: K, KK,KP, KM
      INTEGER(4)  :: N2DOID, NCL1TT
      
      REAL(WP)    :: U_F0_io(0:NCL1_io+1,NCL2,  NCL3)
      REAL(WP)    :: V_F0_io(0:NCL1_io+1,NCL2+1,NCL3)
      REAL(WP)    :: W_F0_io(0:NCL1_io+1,NCL2,  NCL3)
      REAL(WP)    :: P_F0_io(0:NCL1_io+1,NCL2,  NCL3)
      
      REAL(WP)     :: QAUX_io (0:NCL1_io+1,0:N2DO(MYID)+1,NCL3,NDV,1:NPTOT)
      REAL(WP)     :: PRAUX_io(0:NCL1_io+1,0:N2DO(MYID)+1,NCL3,    1:NPTOT)


     IF(.not.IOFLOWflg) RETURN

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

         ! FOR U
         DO I=1,NND1_io
            DO J=2,NCL2
               DO K=1,NCL3
                  JM = JMV(J)
                  KM = KMV(K)
                  U_INTP_io(I,J,K) = 0.250_WP * ( U_F0_io(I,JM,KM) + U_F0_io(I,JM,K) + U_F0_io(I,J,KM) + U_F0_io(I,J,K) )
               END DO
            END DO
         END DO
         ! I = NND1
         !U_INTP_io(NND1_io,1:NCL2,1:NCL3) = U_F0_io(NND1_io,1:NCL2,1:NCL3)
         ! K =1
         U_INTP_io(:,:,NND3) = U_INTP_io(:,:,1)
         ! J = 1 or NND2
         U_INTP_io(:,1,:) = 0.0_WP
         U_INTP_io(:,NND2,:) = 0.0_WP
         
         
         ! FOR V
         DO I=1,NND1_io
            DO J=1,NND2
               DO K=1,NCL3
                  IF(I==NND1_io) THEN 
                    IM = NCL1_io 
                  ELSE
                    IM = IMV_io(I)
                  END IF
                  KM = KMV(K)
                  V_INTP_io(I,J,K) = 0.250_WP * ( V_F0_io(IM,J,KM) + V_F0_io(IM,J,K) + V_F0_io(I,J,KM) + V_F0_io(I,J,K) )
               END DO
            END DO
         END DO
         ! I = NND1
         !V_INTP_io(NND1_io,1:NCL2,1:NCL3) = V_F0_io(NND1_io,1:NCL2,1:NCL3)
         ! K =1
         V_INTP_io(:,:,NND3) = V_INTP_io(:,:,1)
         ! J = 1 or NND2
         V_INTP_io(:,1,:) = 0.0_WP
         V_INTP_io(:,NND2,:) = 0.0_WP
         

         !FOR W
         DO I=1,NND1_io
            DO J=2,NCL2
               DO K=1,NCL3
                  JM = JMV(J)
                  IF(I==NND1_io) THEN 
                    IM = NCL1_io 
                  ELSE
                    IM = IMV_io(I)
                  END IF
                  W_INTP_io(I,J,K) = 0.250_WP * ( W_F0_io(IM,JM,K) + W_F0_io(IM,J,K) + W_F0_io(I,JM,K) + W_F0_io(I,J,K) )
               END DO
            END DO
         END DO
         ! I = NND1
         !W_INTP_io(NND1_io,1:NCL2,1:NCL3) = W_F0_io(NND1_io,1:NCL2,1:NCL3)
         ! K =1
         W_INTP_io(:,:,NND3) = W_INTP_io(:,:,1)
         ! J = 1 or NND2
         W_INTP_io(:,1,:) = 0.0_WP
         W_INTP_io(:,NND2,:) = 0.0_WP
         
         !FOR P
         DO I=1,NCL1_io
            DO J=2,NCL2
               DO K=1,NCL3
                  JM = JMV(J)
                  IF(I==NND1_io) THEN 
                    IM = NCL1_io 
                  ELSE
                    IM = IMV_io(I)
                  END IF
                  KM = KMV(K)
                  P_INTP_io(I,J,K) = 0.1250_WP * ( P_F0_io(IM,JM,KM) + P_F0_io(IM,JM,K) + &
                                                 P_F0_io(IM,J ,KM) + P_F0_io(IM,J ,K) + &
                                                 P_F0_io(I ,JM,KM) + P_F0_io(I ,JM,K) + &
                                                 P_F0_io(I ,J ,KM) + P_F0_io(I ,J ,K) )
               END DO
            END DO
         END DO
         ! I = NND1
         P_INTP_io(NND1_io,1:NCL2,1:NCL3) = P_INTP_io(NCL1_io,1:NCL2,1:NCL3)
         ! K =1
         P_INTP_io(:,:,NND3) = P_INTP_io(:,:,1)
         ! J = 1 or NND2
         P_INTP_io(:,1,:) = P_INTP_io(:,2,:)
         P_INTP_io(:,NND2,:) = P_INTP_io(:,NCL2,:)
         
      ENDIF
      
      RETURN
      END SUBROUTINE TEC360_NODE_EXTRAPOLATION_io

!**********************************************************************************************************************************
     SUBROUTINE  TEC360_WRITE
      use TEC360_INFO
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      INTEGER(4) :: N
     
      IF(MYID.EQ.0) THEN
          ALLOCATE ( U_INTP(NND1,NND2,NND3) )
          ALLOCATE ( V_INTP(NND1,NND2,NND3) )
          ALLOCATE ( W_INTP(NND1,NND2,NND3) )
          ALLOCATE ( P_INTP(NND1,NND2,NND3) )
      END IF
      CALL TEC360_NODE_EXTRAPOLATION
      
      
      IF(IOFLOWflg) THEN
          IF(MYID.EQ.0) THEN
              ALLOCATE ( U_INTP_io(NND1_io,NND2,NND3) )
              ALLOCATE ( V_INTP_io(NND1_io,NND2,NND3) )
              ALLOCATE ( W_INTP_io(NND1_io,NND2,NND3) )
              ALLOCATE ( P_INTP_io(NND1_io,NND2,NND3) )
          END IF
          CALL TEC360_NODE_EXTRAPOLATION_io
      END IF
      

      IF(MYID.EQ.0) THEN
               
!         !===================WRITE ALL DATA==================================         
!         IF((NCOUNT(4).EQ.0) .AND. (NREAD.NE.2).AND. (NREAD_io.NE.2)) THEN
!            CALL TEC360_ALL_NODES_FIRST
!         ELSE
!            CALL TEC360_ALL_NODES_OTHERS 
!         END IF
         
         !===================WRITE X-SLICE================================== 
         IF(IOFLOWflg) THEN 
             IF (NCOUNT(1).EQ.0)  THEN
                DO N=1,NI
                   IID(N) = N*(NND1_IO+1)/(2*NI) + 1
                END DO
             END IF 
              
             IF((NCOUNT(1).EQ.0) .AND. (NREAD.NE.2).AND. (NREAD_io.NE.2)) THEN
                DO N=1,NI
                   CALL TEC360_XSLICE_FIRST(N)
                END DO
             ELSE
                DO N=1,NI
                   CALL TEC360_XSLICE_OTHERS(N) 
                END DO
             END IF
         END IF
         
         !===================WRITE Y-SLICE==================================  
         IF (NCOUNT(2).EQ.0)  THEN
            DO N=1,NJ
               JID(N) = N*(NND2+1)/(2*NJ) + 1
            END DO
         END IF 
          
         IF((NCOUNT(2).EQ.0) .AND. (NREAD.NE.2).AND. (NREAD_io.NE.2)) THEN
            DO N=1,NJ
               CALL TEC360_YSLICE_FIRST(N)
            END DO
         ELSE
            DO N=1,NJ
               CALL TEC360_YSLICE_OTHERS(N) 
            END DO
         END IF
         
         
         !===================WRITE Z-SLICE==================================   
         IF (NCOUNT(3).EQ.0)  THEN
            DO N=1,NK
               KID(N) = N*(NND3+1)/(2*NK) + 1
            END DO
         END IF 
         
         IF((NCOUNT(3).EQ.0) .AND. (NREAD.NE.2).AND. (NREAD_io.NE.2)) THEN
            DO N=1,NK
               CALL TEC360_ZSLICE_FIRST(N)
            END DO
         ELSE
            DO N=1,NK
               CALL TEC360_ZSLICE_OTHERS(N) 
            END DO
         END IF
         
         !===================RELEASE MEMOERY==================================  
         DEALLOCATE ( U_INTP ) 
         DEALLOCATE ( V_INTP )
         DEALLOCATE ( W_INTP )
         DEALLOCATE ( P_INTP ) 
         IF(IOFLOWflg) THEN
             DEALLOCATE ( U_INTP_io ) 
             DEALLOCATE ( V_INTP_io )
             DEALLOCATE ( W_INTP_io )
             DEALLOCATE ( P_INTP_io )
         END IF            

     END IF
     
     END SUBROUTINE TEC360_WRITE

!**********************************************************************************************************************************
     SUBROUTINE TEC360_ALL_NODES_FIRST
      use TEC360_INFO
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      
      INTEGER(4)  :: I, J, K, N1
      INTEGER(4)  :: TECFLG = 202
           
        
        IF(IOFLOWflg) THEN 
           N1 = NND1+NND1_io
        ELSE
           N1 = NND1
        END IF
        
        OPEN(TECFLG,FILE='RESULT.UVWP.ALLPOINTS.plt')
        NCOUNT(4) = NCOUNT(4)+1
        WRITE(TECFLG,'(A)') 'TITLE = "DNS 3D FLOW"'
        WRITE(TECFLG,'(A)') 'VARIABLES = "X", "Y", "Z", "U", "V","W","P"'
        
        WRITE(TECFLG,'(A,1I6.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T=" ',ITERG,phyTIME, &
                      ' ", I=', N1, ', J=',NND2, ', K=',NND3, ', F=POINT'
                   
        DO K=1,NND3
           DO J=1,NND2
              DO I =1,NND1            
                 WRITE(TECFLG,'(7ES15.7)') XND(I),YND(J),ZND(K),U_INTP(I,J,K),V_INTP(I,J,K),W_INTP(I,J,K),P_INTP(I,J,K)
              END DO
              IF(IOFLOWflg) THEN
                  DO I =1,NND1_io            
                     WRITE(TECFLG,'(7ES15.7)') XND_io(I),YND(J),ZND(K), &
                      U_INTP_io(I,J,K),V_INTP_io(I,J,K),W_INTP_io(I,J,K),P_INTP_io(I,J,K)
                  END DO
              END IF
           END DO
        END DO
     
        IF(IOFLOWflg) THEN
            WRITE(TECFLG,'(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID,'// &
                                 'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK,'// &
                                 'F = POINT, S = GLOBAL'
            WRITE(TECFLG,'(I2.1)') 1     ! NUMBER OF LINE
            WRITE(TECFLG,'(I2.1)') 2     ! NUMBER OF POINTS IN LINE 1
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),YND(1)
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),YND(NND2)
        END IF
        
 
        CLOSE(TECFLG)
     
     
     END SUBROUTINE TEC360_ALL_NODES_FIRST
     
     
     
!**********************************************************************************************************************************
     SUBROUTINE TEC360_ALL_NODES_OTHERS
      use TEC360_INFO
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      
      INTEGER(4) :: I, J, K, N1
      INTEGER(4)  :: TECFLG = 202
      
        IF(IOFLOWflg) THEN 
           N1 = NND1+NND1_io
        ELSE
           N1 = NND1
        END IF
        
        OPEN(TECFLG,FILE='RESULT.UVWP.plt', position='append')
        NCOUNT(4) = NCOUNT(4)+1
        WRITE(TECFLG,'(A,1I6.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T=" ',ITERG,phyTIME, &
                      ' ", I=',N1, ', J=',NND2,', K=',NND3,', F=POINT'
        write(TECFLG,'(A)') 'VARSHARELIST=([1-3]=1)'
        
        DO K=1,NND3
           DO J=1,NND2
              DO I =1,NND1            
                 WRITE(TECFLG,'(7ES15.7)') U_INTP(I,J,K),V_INTP(I,J,K),W_INTP(I,J,K),P_INTP(I,J,K)
              END DO
              IF(IOFLOWflg) THEN
                  DO I =1,NND1_io            
                     WRITE(TECFLG,'(4ES15.7)') &
                      U_INTP_io(I,J,K),V_INTP_io(I,J,K),W_INTP_io(I,J,K),P_INTP_io(I,J,K)
                  END DO
              END IF
           END DO
        END DO
        
        IF(IOFLOWflg) THEN
            WRITE(TECFLG,'(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID,'// &
                                 'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK,'// &
                                 'F = POINT, S = GLOBAL'
            WRITE(TECFLG,'(I2.1)') 1
            WRITE(TECFLG,'(I2.1)') 2
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),YND(1)
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),YND(NND2)
            END IF
 
        CLOSE(TECFLG)
     
     END SUBROUTINE TEC360_ALL_NODES_OTHERS
     
 
!**********************************************************************************************************************************
     SUBROUTINE TEC360_XSLICE_FIRST(N)
      use TEC360_INFO
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      
      INTEGER(4) :: I, J, K, N1, N
      CHARACTER(4) :: PNTIM
      INTEGER(4)  :: TECFLG = 202
        
        IF(.NOT.IOFLOWflg) RETURN
        
        TECFLG = TECFLG + IID(N)
        
        WRITE(PNTIM,'(I4.4)') IID(N)
        
        OPEN(TECFLG,FILE='RESULT.UVWP.XSLICE.'//PNTIM//'.plt')
        NCOUNT(1) = NCOUNT(1)+1
        WRITE(TECFLG,'(A)') 'TITLE = "DNS FLOW X-SLICE"'
        WRITE(TECFLG,'(A)') 'VARIABLES = "X", "Y", "Z", "U", "V","W","P"'
        
        WRITE(TECFLG,'(A,1I6.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T=" ',ITERG, phyTIME, &
                      ' ", I=', 1, ', J=',NND2,', K=',NND3,', F=POINT'
                   
        DO K=1,NND3
           DO J=1,NND2
              DO I =IID(N),IID(N)            
                 WRITE(TECFLG,'(7ES15.7)') XND_io(I),YND(J),ZND(K), &
                  U_INTP_io(I,J,K),V_INTP_io(I,J,K),W_INTP_io(I,J,K),P_INTP_io(I,J,K)
              END DO
           END DO
        END DO
        
        IF(IOFLOWflg) THEN
            WRITE(TECFLG,'(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID,'// &
                                 'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK,'// &
                                 'F = POINT, S = GLOBAL'
            WRITE(TECFLG,'(I2.1)') 1
            WRITE(TECFLG,'(I2.1)') 2
            WRITE(TECFLG,'(2ES15.7)') ZND(1),YND(1)
            WRITE(TECFLG,'(2ES15.7)') ZND(1),YND(NND2)
        
        END IF

        CLOSE(TECFLG)
     
     
     END SUBROUTINE TEC360_XSLICE_FIRST
     
!**********************************************************************************************************************************
     SUBROUTINE TEC360_XSLICE_OTHERS(N)
      use TEC360_INFO
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      
      INTEGER(4) :: I, J, K, N1, N
      CHARACTER(4) :: PNTIM
      INTEGER(4)  :: TECFLG = 202
           
        IF(.NOT.IOFLOWflg) RETURN
        TECFLG = TECFLG + IID(N)
        
        WRITE(PNTIM,'(I4.4)') IID(N)
        
        OPEN(TECFLG,FILE='RESULT.UVWP.XSLICE.'//PNTIM//'.plt', position='append')
        NCOUNT(1) = NCOUNT(1)+1
        WRITE(TECFLG,'(A,1I6.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T=" ',ITERG,phyTIME, &
                      ' ", I=', 1, ', J=',NND2,', K=',NND3,', F=POINT'
        write(TECFLG,'(A)') 'VARSHARELIST=([1-3]=1)'

        DO K=1,NND3
           DO J=1,NND2
              IF(IOFLOWflg) THEN
                  DO I =IID(N),IID(N)            
                     WRITE(TECFLG,'(4ES15.7)') &
                      U_INTP_io(I,J,K),V_INTP_io(I,J,K),W_INTP_io(I,J,K),P_INTP_io(I,J,K)
                  END DO
              END IF
           END DO
        END DO
        
        IF(IOFLOWflg) THEN
            WRITE(TECFLG,'(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID,'// &
                                 'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK,'// &
                                 'F = POINT, S = GLOBAL'
            WRITE(TECFLG,'(I2.1)') 1
            WRITE(TECFLG,'(I2.1)') 2
            WRITE(TECFLG,'(2ES15.7)') ZND(1),YND(1)
            WRITE(TECFLG,'(2ES15.7)') ZND(1),YND(NND2)
        
        END IF

        CLOSE(TECFLG)
     
     
     END SUBROUTINE TEC360_XSLICE_OTHERS
     
     
     
!**********************************************************************************************************************************
     SUBROUTINE TEC360_YSLICE_FIRST(N)
      use TEC360_INFO
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      
      INTEGER(4) :: I, J, K, N1, N
      CHARACTER(4) :: PNTIM
      INTEGER(4)  :: TECFLG = 202
           
        TECFLG = TECFLG + JID(N)
        IF(IOFLOWflg) THEN 
           N1 = NND1+NND1_io
        ELSE
           N1 = NND1
        END IF
        WRITE(PNTIM,'(I4.4)') JID(N)
        
        OPEN(TECFLG,FILE='RESULT.UVWP.YSLICE.'//PNTIM//'.plt')
        NCOUNT(2) = NCOUNT(2)+1
        WRITE(TECFLG,'(A)') 'TITLE = "DNS FLOW Y-SLICE"'
        WRITE(TECFLG,'(A)') 'VARIABLES = "X", "Y", "Z", "U", "V","W","P"'
        
        WRITE(TECFLG,'(A,1I6.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T=" ',ITERG,phyTIME, &
                      ' ", I=', N1, ', J=',1,', K=',NND3,', F=POINT'
                   
        DO K=1,NND3
           DO J=JID(N),JID(N)
           
              DO I=1,NND1            
                 WRITE(TECFLG,'(7ES15.7)') XND(I),YND(J),ZND(K),U_INTP(I,J,K),V_INTP(I,J,K),W_INTP(I,J,K),P_INTP(I,J,K)
              END DO
              IF(IOFLOWflg) THEN
                  DO I =1,NND1_io            
                     WRITE(TECFLG,'(7ES15.7)') XND_io(I),YND(J),ZND(K), &
                      U_INTP_io(I,J,K),V_INTP_io(I,J,K),W_INTP_io(I,J,K),P_INTP_io(I,J,K)
                  END DO
              END IF
           END DO
        END DO
        
        IF(IOFLOWflg) THEN
            WRITE(TECFLG,'(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID,'// &
                                 'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK,'// &
                                 'F = POINT, S = GLOBAL'
            WRITE(TECFLG,'(I2.1)') 1
            WRITE(TECFLG,'(I2.1)') 2
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),ZND(1)
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),ZND(NND3)
        
        END IF

        CLOSE(TECFLG)
     
     
     END SUBROUTINE TEC360_YSLICE_FIRST
     
     !**********************************************************************************************************************************
     SUBROUTINE TEC360_YSLICE_OTHERS(N)
      use TEC360_INFO
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      
      INTEGER(4) :: I, J, K, N1, N
      CHARACTER(4) :: PNTIM
      INTEGER(4)  :: TECFLG = 202
           
        TECFLG = TECFLG + JID(N)
        IF(IOFLOWflg) THEN 
           N1 = NND1+NND1_io
        ELSE
           N1 = NND1
        END IF
        WRITE(PNTIM,'(I4.4)') JID(N)
        
        OPEN(TECFLG,FILE='RESULT.UVWP.YSLICE.'//PNTIM//'.plt', position='append')
        NCOUNT(2) = NCOUNT(2)+1
        WRITE(TECFLG,'(A,1I6.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T=" ',ITERG,phyTIME, &
                      ' ", I=', N1, ', J=',1,', K=',NND3,', F=POINT'
        write(TECFLG,'(A)') 'VARSHARELIST=([1-3]=1)'
        
        DO K=1,NND3
           DO J=JID(N),JID(N)
           
              DO I=1,NND1            
                 WRITE(TECFLG,'(4ES15.7)') U_INTP(I,J,K),V_INTP(I,J,K),W_INTP(I,J,K),P_INTP(I,J,K)
              END DO
              IF(IOFLOWflg) THEN
                  DO I =1,NND1_io            
                     WRITE(TECFLG,'(4ES15.7)') &
                      U_INTP_io(I,J,K),V_INTP_io(I,J,K),W_INTP_io(I,J,K),P_INTP_io(I,J,K)
                  END DO
              END IF
           END DO
        END DO
        
        IF(IOFLOWflg) THEN
            WRITE(TECFLG,'(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID,'// &
                                 'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK,'// &
                                 'F = POINT, S = GLOBAL'
            WRITE(TECFLG,'(I2.1)') 1
            WRITE(TECFLG,'(I2.1)') 2
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),ZND(1)
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),ZND(NND3)
        
        END IF

        CLOSE(TECFLG)
     
     
     END SUBROUTINE TEC360_YSLICE_OTHERS
         
!**********************************************************************************************************************************
     SUBROUTINE TEC360_ZSLICE_FIRST(N)
      use TEC360_INFO
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      
      INTEGER(4) :: I, J, K, N1,N
      CHARACTER(4) :: PNTIM
      INTEGER(4)  :: TECFLG = 202
           
        TECFLG = TECFLG + KID(N)
        IF(IOFLOWflg) THEN 
           N1 = NND1+NND1_io
        ELSE
           N1 = NND1
        END IF
        WRITE(PNTIM,'(I4.4)') KID(N)
        
        OPEN(TECFLG,FILE='RESULT.UVWP.ZSLICE.'//PNTIM//'.plt')
        NCOUNT(3) = NCOUNT(3)+1
        WRITE(TECFLG,'(A)') 'TITLE = "DNS FLOW Z-SLICE"'
        WRITE(TECFLG,'(A)') 'VARIABLES = "X", "Y", "Z", "U", "V","W","P"'
        
        WRITE(TECFLG,'(A,1I6.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T=" ',ITERG,phyTIME, &
                      ' ", I=', N1, ', J=',NND2,', K=',1,', F=POINT'
                   

        DO K=KID(N),KID(N)
           DO J=1,NND2
           
              DO I=1,NND1            
                 WRITE(TECFLG,'(7ES15.7)') XND(I),YND(J),ZND(K),U_INTP(I,J,K),V_INTP(I,J,K),W_INTP(I,J,K),P_INTP(I,J,K)
              END DO
              IF(IOFLOWflg) THEN
                      DO I =1,NND1_io            
                         WRITE(TECFLG,'(7ES15.7)') XND_io(I),YND(J),ZND(K), &
                          U_INTP_io(I,J,K),V_INTP_io(I,J,K),W_INTP_io(I,J,K),P_INTP_io(I,J,K)
                      END DO
              END IF
           END DO
        END DO
        
        IF(IOFLOWflg) THEN
            WRITE(TECFLG,'(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID,'// &
                                 'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK,'// &
                                 'F = POINT, S = GLOBAL'
            WRITE(TECFLG,'(I2.1)') 1
            WRITE(TECFLG,'(I2.1)') 2
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),YND(1)
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),YND(NND2)
        
        END IF

        CLOSE(TECFLG)
     
     
     END SUBROUTINE TEC360_ZSLICE_FIRST
     
     
     
!**********************************************************************************************************************************
     SUBROUTINE TEC360_ZSLICE_OTHERS(N)
      use TEC360_INFO
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      
      INTEGER(4) :: I, J, K,N1,N
      CHARACTER(4) :: PNTIM
      INTEGER(4)  :: TECFLG = 202
           
        TECFLG = TECFLG + KID(N)
        IF(IOFLOWflg) THEN 
           N1 = NND1+NND1_io
        ELSE
           N1 = NND1
        END IF
        WRITE(PNTIM,'(I4.4)') KID(N)
 
        OPEN(TECFLG,FILE='RESULT.UVWP.ZSLICE.'//PNTIM//'.plt', position='append')
        NCOUNT(3) = NCOUNT(3)+1
        WRITE(TECFLG,'(A,1I6.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T=" ',ITERG,phyTIME, &
                      ' ", I=', N1, ', J=',NND2,', K=',1,', F=POINT'
        write(TECFLG,'(A)') 'VARSHARELIST=([1-3]=1)'
        
        DO K=KID(N),KID(N)
           DO J=1,NND2
              DO I =1,NND1            
                 WRITE(TECFLG,'(4ES15.7)') U_INTP(I,J,K),V_INTP(I,J,K),W_INTP(I,J,K),P_INTP(I,J,K)
              END DO
              IF(IOFLOWflg) THEN
                      DO I =1,NND1_io            
                         WRITE(TECFLG,'(4ES15.7)') &
                          U_INTP_io(I,J,K),V_INTP_io(I,J,K),W_INTP_io(I,J,K),P_INTP_io(I,J,K)
                      END DO
              END IF
           END DO
        END DO
        IF(IOFLOWflg) THEN
            WRITE(TECFLG,'(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID,'// &
                                 'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK,'// &
                                 'F = POINT, S = GLOBAL'
            WRITE(TECFLG,'(I2.1)') 1
            WRITE(TECFLG,'(I2.1)') 2
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),YND(1)
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),YND(NND2)
        END IF
        CLOSE(TECFLG)
     
     END SUBROUTINE TEC360_ZSLICE_OTHERS




     
     
