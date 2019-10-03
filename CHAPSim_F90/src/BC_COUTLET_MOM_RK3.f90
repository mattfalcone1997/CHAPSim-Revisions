 !**********************************************************************
      SUBROUTINE BC_COUTLET_MOM_RK3(NS)
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      
      INTEGER(4),INTENT(IN) :: NS  
      
      INTEGER(4) :: J, K, IDR, N2I, I, JP, KP,JJ
      REAL(WP)    :: BC_RHS
      REAL(WP)    :: BC_CONV
      REAL(WP)    :: coeffou
      REAL(WP)    :: ACOEF
      
!>==========Convective outflow boundary condition====Method One==============================================         
      CALL CONVCTION_OUTLET_U
      
      N2I = 1
      DO IDR=1, NDV
          IF( (MYID.EQ.0) .AND. (IDR .EQ. 2)) N2I=2
          DO K=1, NCL3
                DO J=N2I, N2DO(myid)
                    BC_CONV=-1.0_WP*( Q_io(NCL1_io+1,J,K,IDR)- &
                               Q_io(NCL1_io,  J,K,IDR) )*DXI*U_OUTLET 
                    !BC_CONV= -1.0_WP*( Q_io(NCL1_io,J,K,IDR)- &
                     !                Q_io(NCL1_io-1,J,K,IDR) )*DXI*U_OUTLET 
                    BC_RHS = TGAM(NS)*BC_CONV + TROH(NS)*BC_CONV0(J,K,IDR)
                    BC_CONV0(J,K,IDR)=BC_CONV
                    
                    BC_TDMA(3,J,K,IDR) = DT * BC_RHS
                    Q_io(NCL1_io+1,J,K,IDR) = Q_io(NCL1_io+1,J,K,IDR) + &
                                              BC_TDMA(3,J,K,IDR)
                    DPH_io(NCL1_io+1, J, K)    = 0.0_WP
                    !DPH_io(NCL1_io+1, J, K)    = DPH_io(NCL1_io-1, J, K)
                END DO
          END DO
          
      END DO
      
      ! bottom wall, not necessary, given in intf
      IF(MYID==0) THEN
         DO K=1,NCL3
            Q_io(NCL1_io+1,1,K,2) = 0.0_WP
            Q_io(NCL1_io+1,0,K,2) = 0.0_WP
         END DO
      END IF
      ! top wall, not necessary, given in intf
      IF(MYID==NPSLV) THEN
         DO K=1,NCL3
            Q_io(NCL1_io+1,N2DO(MYID)+1,K,2) = 0.0_WP
         END DO
      END IF

!!>=========================CBC Method Two=========================================
!      CALL VELOUPDT_CBC_U
      
!      CALL CONVCTION_OUTLET_U
      
!     Calculate V, W for BC*********************
!      DO IDR=1,NDV
!          IF( (MYID.EQ.0) .AND. (IDR .EQ. 2)) N2I=2
!          DO K=1, NCL3
!                DO J=1, N2DO(MYID) !N2I, N2DO(myid)
!                    BC_CONV=-1.0_WP*( Q_io(NCL1_io,J,K,IDR)- &
!                                    Q_io(NCL1_io-1,J,K,IDR) )*DXI*U_OUTLET 
!                    BC_RHS = GAL*BC_CONV + ROL*BC_CONV0(J,K,IDR)
!                    BC_CONV0(J,K,IDR)=BC_CONV
                    
!                    BC_TDMA(3,J,K,IDR) = Q_io(NCL1_io+1,J,K,IDR) !last time step
!                    Q_io(NCL1_io+1,J,K,IDR) = Q_io(NCL1_io+1,J,K,IDR) + DT*BC_RHS
                                              
!                    DPH_io(NCL1_io+1, J, K)    = 0.0_WP
!                    PR_io (NCL1_io+1, J, K)    = PR_io(NCL1_io, J, K)  !not used!
!                END DO
!          END DO
          
!      END DO
      
!       bottom wall, not necessary, given in intf
!      IF(MYID==0) THEN
!         DO K=1,NCL3
!            Q_io(NCL1_io+1,1,K,2) = 0.0_WP
!            Q_io(NCL1_io+1,0,K,2) = 0.0_WP
!         END DO
!      END IF
!       top wall, not necessary, given in intf
!      IF(MYID==NPSLV) THEN
!         DO K=1,NCL3
!            Q_io(NCL1_io+1,N2DO(MYID)+1,K,2) = 0.0_WP
!         END DO
!      END IF
!      CALL INTFC_BC_QP_INOUT
      
!!     Calculate U***************     
!      CALL VELOUPDT_CBC_U
      
!      DO K=1, NCL3
!         KP=KPV(K)
!         DO J=1, N2DO(myid)
!            JP=J+1
!            JJ=JCL2G(J)
!            BC_CONV=   -1.0_WP* &
!                    ( (Q_io(NCL1_io,JP,K,2)-Q_io(NCL1_io,J,K,2))*DYFI(JJ) + &
!                      (Q_io(NCL1_io,J,KP,3)-Q_io(NCL1_io,J,K,3))*DZI ) !explicit ( U_{n+1}-U_n )/ dx
                      
!            BC_RHS = GAL*BC_CONV + ROL*BC_CONV0(J,K,1)
!            BC_CONV0(J,K,1)=BC_CONV
                    
!            BC_TDMA(3,J,K,1) = Q_io(NCL1_io+1,J,K,1) !last time step
!            Q_io(NCL1_io+1,J,K,1) = Q_io(NCL1_io+1,J,K,1) + DT*BC_RHS                        
!         END DO
!     END DO

!>========================================Check Mass Conservation===================================================      

!      CALL VELOUPDT_CBC_U
      
!      CALL MassConservationCheck(coeffou) 
      
!      write(*,*) 'coeffou',coeffou
      
!      Q_io(NCL1_io+1,:,:,:) = Q_io(NCL1_io+1,:,:,:)*coeffou
      
      
!      DO J = 1, N2DO(MYID)
!         DO K=1, NCL3
!             DO IDR=1,NDV
!                 BC_TDMA(3,J,K,IDR) = Q_io(NCL1_io+1, J, K , IDR)-BC_TDMA(3,J,K,IDR)
!             END DO
!         END DO
!      END DO 

      
      RETURN
      END SUBROUTINE BC_COUTLET_MOM_RK3

!**********************************************************************
 !**********************************************************************
      SUBROUTINE MassConservationCheck(coeffou)
!     Below is based on constant density in the domian.

      use mpi_info
      use init_info
      use mesh_info
      use flow_info
      IMPLICIT NONE
      
      REAL(WP) :: FluxINTG, FluxINTG0
      REAL(WP) :: FluxINTG_WORK, FluxINTG_WORK0
      REAL(WP) :: AreaINTG
      REAL(WP) :: AreaINTG_WORK
      REAL(WP) :: coeffou
      
      INTEGER(4) :: JC, KC, JJ, I, J, K
      
      REAL(WP) :: Uxmax, Uxmin, Uxmax1, Uxmin1
               
          CALL MPI_BARRIER(ICOMM,IERROR)
          
          FluxINTG =0.0_WP
          FluxINTG0=0.0_WP
          AreaINTG =0.0_WP            
          DO JC=1,N2DO(myid)
             DO KC=1,NCL3
                JJ=JCL2G(JC)
                AreaINTG    =AreaINTG   + 1.00_WP/DYFI(JJ)/DZI
                FluxINTG0   =FluxINTG0  + Q_io(0,        JC,KC,NFLOW)/DYFI(JJ)/DZI
                FluxINTG    =FluxINTG   + Q_io(NCL1_io+1,JC,KC,NFLOW)/DYFI(JJ)/DZI
             ENDDO
          ENDDO
          CALL MPI_ALLREDUCE(AreaINTG,AreaINTG_WORK,1,MPI_DOUBLE_PRECISION,&
                             MPI_SUM, ICOMM, IERROR)  
          CALL MPI_ALLREDUCE(FluxINTG,FluxINTG_WORK,1,MPI_DOUBLE_PRECISION,&
                             MPI_SUM, ICOMM, IERROR)
          CALL MPI_ALLREDUCE(FluxINTG0,FluxINTG_WORK0,1,MPI_DOUBLE_PRECISION,&
                             MPI_SUM, ICOMM, IERROR)   

          coeffou = FluxINTG_WORK0/FluxINTG_WORK
          !if(myid==0 .and. dabs(coeffou-1.0_WP)>1.0_WP ) write(*,*) 'coeffou2',ITERG, coeffou,FluxINTG_WORK0,FluxINTG_WORK
          if(myid==0 ) write(*,*) 'coeffou2',ITERG, coeffou,FluxINTG_WORK0,FluxINTG_WORK
          
      RETURN
      END SUBROUTINE MassConservationCheck
      
 !**********************************************************************
      SUBROUTINE CONVCTION_OUTLET_U
!     Below is based on constant density in the domian.

      use mpi_info
      use init_info
      use mesh_info
      use flow_info
      IMPLICIT NONE
      
      REAL(WP) :: FluxINTG
      REAL(WP) :: FluxINTG_WORK
      REAL(WP) :: AreaINTG
      REAL(WP) :: AreaINTG_WORK
      
      INTEGER(4) :: JC, KC, JJ, I, J, K
      
      REAL(WP) :: Uxmax, Uxmin, Uxmax1, Uxmin1
      
!>================Method One========================================          
!          FluxINTG=0.0_WP
!          AreaINTG=0.0_WP            
!          DO JC=1,N2DO(myid)
!             DO KC=1,NCL3
!                JJ=JCL2G(JC)
!                AreaINTG   =AreaINTG   + 1.00_WP/DYFI(JJ)/DZI
!                FluxINTG   =FluxINTG   + Q_io(0,JC,KC,NFLOW)/DYFI(JJ)/DZI
!             ENDDO
!          ENDDO
!          CALL MPI_ALLREDUCE(AreaINTG,AreaINTG_WORK,1,MPI_DOUBLE_PRECISION,&
!                             MPI_SUM, ICOMM, IERROR)  
!          CALL MPI_ALLREDUCE(FluxINTG,FluxINTG_WORK,1,MPI_DOUBLE_PRECISION,&
!                             MPI_SUM, ICOMM, IERROR)   
          
!          U_OUTLET = FluxINTG_WORK/AreaINTG_WORK
          
!          U_OUTLET = 0.10_WP * U_OUTLET
!          !write(*,*) 'U_OUTLET1',U_OUTLET
!>=============Method Two==============================================
          uxmax=-1.0E10
          uxmin= 1.0E10
          do k=1,NCL3
             do j=1,N2DO(MYID)
                if (Q_io(NCL1_io,J,K,NFLOW).gt.uxmax) uxmax=Q_io(NCL1_io,J,K,NFLOW)
                if (Q_io(NCL1_io,J,K,NFLOW).lt.uxmin) uxmin=Q_io(NCL1_io,J,K,NFLOW)
             enddo
          enddo
          call MPI_ALLREDUCE(uxmax,uxmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX,ICOMM, IERROR)
          call MPI_ALLREDUCE(uxmin,uxmin1,1,MPI_DOUBLE_PRECISION,MPI_MIN,ICOMM, IERROR)

          U_OUTLET = 0.5*(uxmax1+uxmin1)

      RETURN
      END SUBROUTINE CONVCTION_OUTLET_U
      
!**********************************************************************
      SUBROUTINE BC_CBC_TDMA(NS)
!     Check N2I = 1 or 2
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      
      INTEGER(4),INTENT(IN) :: NS  
      
      INTEGER(4) :: J, JJ, JM, JC, JP
      INTEGER(4) :: K,     KM, KC, KP
      INTEGER(4) :: IDR
      REAL(WP)    :: AMJV, ACJV, APJV
      REAL(WP)    :: AMKV, ACKV, APKV
      REAL(WP)    :: COE1
      
      
      COE1=0.50_WP*TALP(NS)*DT*CVISC
      
      CALL BC_INTFC_CBC(3)
!>    @note : With BC_INTFC_CBC(3) all known, and tri-matrix, 
!>            TRANSP23 is not necessary here.
      DO IDR = 1, NDV 
         DO K=1, NCL3
            DO JC=1,N2DO(myid)
               JM = JC-1
               JP = JC+1
               JJ = JCL2G(JC)
               ACJV = 1.0_WP-COE1*ACVR(JJ,IDR)
               APJV =     -COE1*APVR(JJ,IDR)
               AMJV =     -COE1*AMVR(JJ,IDR)
               BC_TDMA(2, JC, K, IDR) = AMJV*BC_TDMA(3,JM,K,IDR)+ &
                                        ACJV*BC_TDMA(3,JC,K,IDR)+ &
                                        APJV*BC_TDMA(3,JP,K,IDR)
            END DO
         END DO
      END DO
       
      !CALL BC_INTFC_CBC(2) not necessary

      DO IDR = 1, NDV 
         DO J=1,N2DO(myid)
            DO KC=1, NCL3
               KM = KMV(KC)
               KP = KPV(KC)
               ACKV = 1.0_WP-COE1*(-2.0_WP*DZQI)
               APKV =     -COE1*DZQI
               AMKV =     -COE1*DZQI
               BC_TDMA(1, J, KC, IDR) = AMKV*BC_TDMA(2,J,KM,IDR)+ &
                                        ACKV*BC_TDMA(2,J,KC,IDR)+ &
                                        APKV*BC_TDMA(2,J,KP,IDR)
            END DO
         END DO
      END DO
      
      !CALL BC_INTFC_CBC(1) not necessary
      

      RETURN
      END SUBROUTINE BC_CBC_TDMA

!**********************************************************************

!************************************************************************
      SUBROUTINE BC_INTFC_CBC(FASTEP)
      use mpi_info
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE

       INTEGER(4),INTENT(IN)  :: FASTEP
       
       INTEGER(4)  :: K
       INTEGER(4)  :: IDR
       INTEGER(4)  :: ITAG
       INTEGER(4)  :: IDESF
       INTEGER(4)  :: IDESB
       INTEGER(4)  :: TRC_STS(MPI_STATUS_SIZE)
       INTEGER(4)  :: NSZ

       BSEN_F_io = 0.0_WP
       BSEN_L_io = 0.0_WP
       BREC_F_io = 0.0_WP
       BREC_L_io = 0.0_WP
       NSZ      = 1*NCL3*NDV 
        
                       
       DO IDR=1,NDV
          DO K=1,NCL3
             BSEN_F_io(K,IDR)=BC_TDMA(FASTEP,1,         K,IDR) !y=local Floor
             BSEN_L_io(K,IDR)=BC_TDMA(FASTEP,N2DO(MYID),K,IDR) !y=local ROOF
          ENDDO
       ENDDO

       
       ITAG=0      
!>     @note rank=2,3,4...SIZE-1 (no b.c.)***************************************      
       IF ( (MYID.GT.0) .AND. (MYID.LT.NPSLV) ) THEN
          IDESF = MYID+1    !next myid, forward
          IDESB = MYID-1    !last myid, backward


          CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                MPI_DOUBLE_PRECISION,IDESB,ITAG,BREC_F_io(1,1),NSZ,&
                MPI_DOUBLE_PRECISION,IDESF,ITAG,ICOMM,TRC_STS,IERROR)
        
          CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                   &
                MPI_DOUBLE_PRECISION,IDESF,ITAG,BREC_L_io(1,1),NSZ, &
                MPI_DOUBLE_PRECISION,IDESB,ITAG,ICOMM,TRC_STS,IERROR)


!>        @note Constructing p,u,v,w of ghost cells y_id=0 and N2NO+1, received from adjacent partitions.      
         
          DO IDR=1,NDV
             DO K=1,NCL3
                BC_TDMA(FASTEP,0,           K,IDR)=BREC_L_io(K,IDR)
                BC_TDMA(FASTEP,N2DO(MYID)+1,K,IDR)=BREC_F_io(K,IDR)
             ENDDO
         ENDDO
 
       ENDIF
       
!>     @note rank=SIZE-1 (top wall)***************************************
       IF(MYID.EQ.NPSLV) THEN
          IDESF = 0  
          IDESB = MYID-1   
          
          CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                MPI_DOUBLE_PRECISION, IDESB, ITAG,BREC_F_io(1,1),NSZ,&
                MPI_DOUBLE_PRECISION, IDESF, ITAG,ICOMM,TRC_STS,IERROR)
                
          CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,               &
                MPI_DOUBLE_PRECISION, IDESF, ITAG,BREC_L_io(1,1),NSZ, &
                MPI_DOUBLE_PRECISION, IDESB, ITAG,ICOMM,TRC_STS,IERROR)
          
          DO IDR=1,3
              DO K=1,NCL3
                 BC_TDMA(FASTEP,0, K,IDR)=BREC_L_io(K,IDR)
              END DO
          ENDDO
          
          
          DO K=1,NCL3
                 BC_TDMA(FASTEP,N2DO(myid)+1, K,1) = -BC_TDMA(FASTEP,N2DO(myid), K,1)
                 BC_TDMA(FASTEP,N2DO(myid)+1, K,3) = -BC_TDMA(FASTEP,N2DO(myid), K,3)
                 BC_TDMA(FASTEP,N2DO(myid)+1, K,2) = 0.0_WP     
          ENDDO               
          
       ENDIF

!>     @note rank=0 (bottom wall)***************************************
       IF (MYID.EQ.0) THEN
          IDESF = MYID+1
          IDESB = NPSLV
         
          CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                 &
                MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F_io(1,1),NSZ,  &
                MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
          CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                 &
                MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L_io(1,1),NSZ,  &
                MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)
                
!>        @note bottom wall velocity B.C.
          DO IDR=1,3
               DO K=1,NCL3 
                  BC_TDMA(FASTEP,N2DO(myid)+1, K,IDR)=BREC_F_io(K,IDR)
               END DO
          ENDDO          
          
          DO K=1,NCL3                
              BC_TDMA(FASTEP,0,K,1)=-BC_TDMA(FASTEP,1,K,1)
              BC_TDMA(FASTEP,0,K,3)=-BC_TDMA(FASTEP,1,K,3)
              BC_TDMA(FASTEP,1,K,2)=0.0_WP
              BC_TDMA(FASTEP,0,K,2)=0.0_WP
         ENDDO
          
       ENDIF
       
      RETURN
      END SUBROUTINE BC_INTFC_CBC
      
!***********************************************************************  

!************************
      SUBROUTINE VELOUPDT_CBC_U
      use init_info
      use mesh_info
      use flow_info
      IMPLICIT NONE 
      
      INTEGER(4) :: I,J,K,JP,JJ,KP
      REAL(WP)    :: coeffou
      
      DO J = 1, N2DO(MYID)
         JP=J+1
         JJ=JCL2G(J)
         DO K=1, NCL3
             KP=KPV(K)
             BC_TDMA(3,J,K,1) = Q_io(NCL1_io+1, J, K, 1)
             
             Q_io(NCL1_io+1, J, K, 1)= -1.0_WP* &
             ( (Q_io(NCL1_io,JP,K,2)-Q_io(NCL1_io,J,K,2))*DYFI(JJ) + &
               (Q_io(NCL1_io,J,KP,3)-Q_io(NCL1_io,J,K,3))*DZI ) /DXI &
             + Q_io(NCL1_io,J,K,1)
         END DO
      END DO
      
      !CALL MassConservationCheck(coeffou)
      !write(*,*) 'coeffou',coeffou
      
      RETURN
      END SUBROUTINE VELOUPDT_CBC_U
