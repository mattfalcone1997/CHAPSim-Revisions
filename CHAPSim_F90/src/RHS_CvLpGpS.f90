      SUBROUTINE RHS_CvLpGpS_tg(NS,IDR)   
      use init_info
      use mesh_info
      use flow_info
      use postprocess_info
      IMPLICIT NONE

      INTEGER(4),INTENT(IN) :: NS
      INTEGER(4),INTENT(IN) :: IDR

      REAL(WP) :: INTGRHSY
      REAL(WP) :: INTGRHSY_WORK
      REAL(WP) :: INTGU
      REAL(WP) :: INTGU_WORK
      REAL(WP) :: DPGRNS
      REAL(WP) :: INTGVOL
      REAL(WP) :: INTGVOL_WORK       
      REAL(WP) :: RHSCCC
      REAL(WP) :: SUCACJ
      REAL(WP) :: DERQ
      REAL(WP) :: RHSC, RHSL
      REAL(WP) :: VOLtmp
      REAL(WP) :: PGM      
      INTEGER(4) :: NII
      INTEGER(4) :: I, IC, IM, IP
      INTEGER(4) :: J, JC, JM, JP, JJ
      INTEGER(4) :: K, KC, KM, KP
      REAL(WP)    :: RMC2(N2DO(MYID))
      REAL(WP)    :: CONVH(NCL1,N2DO(0),NCL3) 
      REAL(WP)    :: COE1, COE2
       
      NII=1
      IF(IDR.EQ.2) THEN
      
         IF(MYID.EQ.0) THEN
            NII=2
         END IF
         
      END IF
 
      RMC2 = 0.0_WP
      IF(IDR.EQ.2)THEN
         DO JC=NII, N2DO(MYID)
             JJ=JCL2G(JC)
             RMC2(JC)=RC(JJ)**2
         END DO
      ELSE
         DO JC=NII, N2DO(MYID)
             JJ=JCL2G(JC)
             RMC2(JC)=RM(JJ)**2
         END DO
      END IF  
      
      CONVH = 0.0_WP
      IF(IDR.EQ.1) THEN
      
         DO I=1,NCL1
            DO J=1,N2DO(MYID)
               DO K=1,NCL3
                  CONVH(I,J,K) = Qtmp(I,J,K)
               END DO
            END DO
         END DO
      
      ELSE IF(IDR .EQ. 2) THEN
      
         DO I=1,NCL1
            DO J=1,N2DO(MYID)
               DO K=1,NCL3
                  CONVH(I,J,K) = DPH(I,J,K)
               END DO
            END DO
         END DO
      
      
      ELSE IF(IDR .EQ. 3) THEN
      
         DO I=1,NCL1
            DO J=1,N2DO(MYID)
               DO K=1,NCL3
                  CONVH(I,J,K) = RHSLLPHI(I,J,K)
               END DO
            END DO
         END DO
      
      ELSE
      END IF
             
      RHSC = 0.0_WP
      RHSL = 0.0_WP
      RHSCCC = 0.0_WP
      COE1 = TALP(NS)*CVISC
      DO KC=1,NCL3
         KM=KMV(KC)
         KP=KPV(KC)
         DO JC=NII,N2DO(MYID)
            JM=JC-1
            JP=JC+1
            JJ=JCL2G(JC)
            DO IC=1,NCL1
               IP=IPV(IC)
               IM=IMV(IC)
               RHSC  = TGAM(NS)*CONVH(IC,JC,KC)+TROH(NS)*CONVH0(IC,JC,KC,IDR)
               RHSL  = COE1*(                                           &
                       (       Q(IP,JC,KC,IDR)-                         &
                        2.0_WP*Q(IC,JC,KC,IDR)+                         &
                               Q(IM,JC,KC,IDR)   )*DXQI+                &
                       (       Q(IC,JC,KP,IDR)-                         &
                        2.0_WP*Q(IC,JC,KC,IDR)+                         &
                               Q(IC,JC,KM,IDR)   )*DZQI/RMC2(JC)+       &
                       (       Q(IC,JP,KC,IDR)*APVR(JJ,IDR)+             &
                               Q(IC,JC,KC,IDR)*ACVR(JJ,IDR)+             &
                               Q(IC,JM,KC,IDR)*AMVR(JJ,IDR) )            &
                       )

               RHSCCC=(RHSC+RHSL)*DT
               
               CONVH0(IC,JC,KC,IDR)=CONVH(IC,JC,KC)

               RHS(IC,JC,KC)=RHSCCC

            ENDDO
         ENDDO
      ENDDO
            
      PGM=0.0_WP
      COE2=TALP(NS)*DT
      IF (IDR.EQ.1) THEN 
          DO KC=1,NCL3
             DO JC=NII,N2DO(MYID)
                DO IC=1,NCL1
                   IM=IMV(IC)
                   PGM = (PR(IC,JC,KC)-PR(IM,JC,KC))*DXI*COE2 
                   RHS(IC,JC,KC)=RHS(IC,JC,KC)-PGM
                ENDDO
             ENDDO
          ENDDO
      ELSE IF (IDR.EQ.2) THEN 
          DO KC=1,NCL3
             DO JC=NII,N2DO(MYID)
                JM=JC-1
                JJ=JCL2G(JC)
                SUCACJ=DYCI(JJ)*rc(jj)
                DO IC=1,NCL1
                   PGM = (PR(IC,JC,KC)-PR(IC,JM,KC))*SUCACJ*COE2
                   RHS(IC,JC,KC)=RHS(IC,JC,KC)-PGM
                ENDDO
             ENDDO
         ENDDO      
      ELSE IF (IDR.EQ.3) THEN
          DO KC=1,NCL3
             KM=KMV(KC)
             DO JC=NII,N2DO(MYID)
                DO IC=1,NCL1
                   PGM = (PR(IC,JC,KC)-PR(IC,JC,KM))*DZI*COE2
                   RHS(IC,JC,KC)=RHS(IC,JC,KC)-PGM
                ENDDO
             ENDDO
          ENDDO  
      ELSE       
      ENDIF

      DPGRNS=0.0_WP
      DERQ=0.0_WP
      IF (IDR.EQ.NFLOW) THEN

          IF(FLOWTP==1) THEN         
             INTGRHSY=0.0_WP
             INTGU=0.0_WP
             INTGVOL=0.0_WP
          
             IF(iswitch.EQ.1)THEN
               VOLtmp = 1.0_WP
             ELSE IF(iswitch.EQ.2)THEN
               VOLtmp = DXI*DZI
             ELSE
          
             END IF 
                       
             DO IC=1,NCL1
                DO JC=1,N2DO(MYID)
                   DO KC=1,NCL3
                      JJ=JCL2G(JC)
                      INTGRHSY=INTGRHSY+RHS(IC,JC,KC)/DYFI(JJ)*rm(jj)/VOLtmp         
                      INTGVOL =INTGVOL +1.0_WP/DYFI(JJ)*rm(jj)/VOLtmp
                      INTGU   =INTGU   +Q(IC,JC,KC,NFLOW)/DYFI(JJ)
                   ENDDO
                ENDDO
             ENDDO
 
             CALL MPI_ALLREDUCE(INTGVOL,INTGVOL_WORK,1,MPI_DOUBLE_PRECISION, &
                             MPI_SUM, ICOMM, IERROR)   
             CALL MPI_ALLREDUCE(INTGRHSY,INTGRHSY_WORK,1,MPI_DOUBLE_PRECISION, &
                             MPI_SUM, ICOMM, IERROR)     
             CALL MPI_ALLREDUCE(INTGU,INTGU_WORK,1,MPI_DOUBLE_PRECISION,&
                             MPI_SUM, ICOMM, IERROR)     
          
          
             DPGRNS=INTGRHSY_WORK/INTGVOL_WORK
          ELSE IF(FLOWTP==2) THEN 
             DPGRNS = -0.50_WP*CFGV*COE2
          ELSE
          END IF   

          DERQ=DPCONS*COE2
          
          DO K=1,NCL3
             DO I=1,NCL1
                DO J=NII,N2DO(MYID)
                   RHS(I,J,K)=RHS(I,J,K) - ( DPGRNS - DERQ )
                ENDDO
             ENDDO
          ENDDO

      END IF 

      RETURN
      END SUBROUTINE RHS_CvLpGpS_tg
      
      SUBROUTINE RHS_CvLpGpS_io(NS,IDR)   !@
      use init_info
      use mesh_info
      use flow_info
      use postprocess_info
      IMPLICIT NONE

      INTEGER(4),INTENT(IN) :: NS
      INTEGER(4),INTENT(IN) :: IDR

      REAL(WP) :: INTGRHSY
      REAL(WP) :: INTGRHSY_WORK
      REAL(WP) :: INTGU
      REAL(WP) :: INTGU_WORK
      REAL(WP) :: DPGRNS
      REAL(WP) :: INTGVOL
      REAL(WP) :: INTGVOL_WORK       
      REAL(WP) :: RHSCCC
      REAL(WP) :: SUCACJ
      REAL(WP) :: DERQ
      REAL(WP) :: RHSC, RHSL
      REAL(WP) :: VOLtmp
      REAL(WP) :: PGM      
      INTEGER(4) :: NII
      INTEGER(4) :: I, IC, IM, IP
      INTEGER(4) :: J, JC, JM, JP, JJ
      INTEGER(4) :: K, KC, KM, KP
      REAL(WP)    :: RMC2(N2DO(MYID))
      REAL(WP)    :: CONVH_io(NCL1_io,N2DO(0),NCL3) 
      REAL(WP)    :: COE1, COE2 

      NII=1
      IF(IDR.EQ.2) THEN
      
         IF(MYID.EQ.0) THEN
            NII=2
         END IF
         
      END IF
 
      RMC2 = 0.0_WP
      IF(IDR.EQ.2)THEN
         DO JC=NII, N2DO(MYID)
             JJ=JCL2G(JC)
             RMC2(JC)=RC(JJ)**2
         END DO
      ELSE
         DO JC=NII, N2DO(MYID)
             JJ=JCL2G(JC)
             RMC2(JC)=RM(JJ)**2
         END DO
      END IF  
      
      
      CONVH_io = 0.0_WP
      IF(IDR.EQ.1) THEN
      
         DO I=1,NCL1_io
            DO J=1,N2DO(MYID)
               DO K=1,NCL3
                  CONVH_io(I,J,K) = Qtmp_io(I,J,K)
               END DO
            END DO
         END DO
      
      ELSE IF(IDR .EQ. 2) THEN
      
         DO I=1,NCL1_io
            DO J=1,N2DO(MYID)
               DO K=1,NCL3
                  CONVH_io(I,J,K) = DPH_io(I,J,K)
               END DO
            END DO
         END DO
      
      
      ELSE IF(IDR .EQ. 3) THEN
      
         DO I=1,NCL1_io
            DO J=1,N2DO(MYID)
               DO K=1,NCL3
                  CONVH_io(I,J,K) = RHSLLPHI_io(I,J,K)
               END DO
            END DO
         END DO
      
      ELSE
      END IF
             
      RHSC = 0.0_WP
      RHSL = 0.0_WP
      RHSCCC = 0.0_WP
      COE1 = TALP(NS)*CVISC
      DO KC=1,NCL3
         KM=KMV(KC)
         KP=KPV(KC)
         DO JC=NII,N2DO(MYID)
            JM=JC-1
            JP=JC+1
            JJ=JCL2G(JC)
            DO IC=1,NCL1_io
               IP=IPV_io(IC)
               IM=IMV_io(IC)
               RHSC  = TGAM(NS)*CONVH_io(IC,JC,KC)+TROH(NS)*CONVH0_io(IC,JC,KC,IDR)
               RHSL  = COE1*(                                    &
                       (     Q_io(IP,JC,KC,IDR)-                         &
                        2.0_WP*Q_io(IC,JC,KC,IDR)+                         &
                             Q_io(IM,JC,KC,IDR)   )*DXQI+              &
                       (     Q_io(IC,JC,KP,IDR)-                         &
                        2.0_WP*Q_io(IC,JC,KC,IDR)+                         &
                             Q_io(IC,JC,KM,IDR)   )*DZQI/RMC2(JC)+     &
                       (     Q_io(IC,JP,KC,IDR)*APVR(JJ,IDR)+             &
                             Q_io(IC,JC,KC,IDR)*ACVR(JJ,IDR)+             &
                             Q_io(IC,JM,KC,IDR)*AMVR(JJ,IDR) )            &
                       )

               RHSCCC=(RHSC+RHSL)*DT
               

               CONVH0_io(IC,JC,KC,IDR)=CONVH_io(IC,JC,KC)

               RHS_io(IC,JC,KC)=RHSCCC

            ENDDO
         ENDDO
      ENDDO

      PGM=0.0_WP
      COE2 = TALP(NS)*DT 
      IF (IDR.EQ.1) THEN 

          DO KC=1,NCL3
             DO JC=NII,N2DO(MYID)
                DO IC=1,NCL1_io
                   IM=IMV_io(IC)
                   PGM = (PR_io(IC,JC,KC)-PR_io(IM,JC,KC))*DXI*COE2 
                   RHS_io(IC,JC,KC)=RHS_io(IC,JC,KC)-PGM
                ENDDO
             ENDDO
          ENDDO
      ELSE IF (IDR.EQ.2) THEN 

          DO KC=1,NCL3
             DO JC=NII,N2DO(MYID)
                JM=JC-1
                JJ=JCL2G(JC)
                SUCACJ=DYCI(JJ)*rc(jj)
                DO IC=1,NCL1_io
                   PGM = (PR_io(IC,JC,KC)-PR_io(IC,JM,KC))*SUCACJ*COE2
                   RHS_io(IC,JC,KC)=RHS_io(IC,JC,KC)-PGM
                ENDDO
             ENDDO
         ENDDO      
      ELSE IF (IDR.EQ.3) THEN

          DO KC=1,NCL3
             KM=KMV(KC)
             DO JC=NII,N2DO(MYID)
                DO IC=1,NCL1_io
                   PGM = (PR_io(IC,JC,KC)-PR_io(IC,JC,KM))*DZI*COE2
                   RHS_io(IC,JC,KC)=RHS_io(IC,JC,KC)-PGM
                ENDDO
             ENDDO
          ENDDO  
      ELSE       
      ENDIF
    

      DPGRNS=0.0_WP
      DERQ=0.0_WP
      IF (IDR.EQ.NFLOW) THEN
 
          IF(FLOWTP ==1) THEN         
             INTGRHSY=0.0_WP
             INTGU=0.0_WP
             INTGVOL=0.0_WP
          
             IF(iswitch.EQ.1)THEN
               VOLtmp = 1.0_WP
             ELSE IF(iswitch.EQ.2)THEN
               VOLtmp = DXI*DZI
             ELSE
             END IF 
                       
             DO IC=1,NCL1_io
                DO JC=1,N2DO(MYID)
                   DO KC=1,NCL3
                      JJ=JCL2G(JC)
                      INTGRHSY=INTGRHSY+RHS_io(IC,JC,KC)/DYFI(JJ)*rm(jj)/VOLtmp         
                      INTGVOL =INTGVOL +1.0_WP/DYFI(JJ)*rm(jj)/VOLtmp
                      INTGU   =INTGU   +Q_io(IC,JC,KC,NFLOW)/DYFI(JJ)
                   ENDDO
                ENDDO
             ENDDO
 
             CALL MPI_ALLREDUCE(INTGVOL,INTGVOL_WORK,1,MPI_DOUBLE_PRECISION, &
                             MPI_SUM, ICOMM, IERROR)   
             CALL MPI_ALLREDUCE(INTGRHSY,INTGRHSY_WORK,1,MPI_DOUBLE_PRECISION, &
                             MPI_SUM, ICOMM, IERROR)     
             CALL MPI_ALLREDUCE(INTGU,INTGU_WORK,1,MPI_DOUBLE_PRECISION,&
                             MPI_SUM, ICOMM, IERROR)     

             DPGRNS=INTGRHSY_WORK/INTGVOL_WORK
          ELSE IF(FLOWTP==2) THEN
             DPGRNS = -0.50_WP*CFGV*COE2 
          ELSE
          END IF   

          DERQ=DPCONS*COE2
          
          DO K=1,NCL3
             DO I=1,NCL1_io
                DO J=NII,N2DO(MYID)
                   RHS_io(I,J,K)=RHS_io(I,J,K) - ( DPGRNS - DERQ )
                ENDDO
             ENDDO
          ENDDO

      END IF 

      RETURN
      END SUBROUTINE RHS_CvLpGpS_io

