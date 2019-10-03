
       SUBROUTINE INIFIELD_FLOW_tg  
       use init_info
       use mesh_info
       use flow_info
       IMPLICIT NONE      
       
       
       INTEGER(4) :: I
       INTEGER(4) :: J, JJ
       INTEGER(4) :: K
       INTEGER(4) :: RDCOUNT, seed
!       INTEGER(4) :: IDUM
       REAL(WP)    :: VPERQ
       REAL(WP)    :: XRDM(3)
       REAL(WP)    :: XRD(3)
       REAL(WP), ALLOCATABLE    :: YCLtmp(:)
       
       INTEGER(4)  :: RANDOMTYPE=2  !=1 for real random, =2 for fixed random
       
       ALLOCATE ( YCLtmp(NCL2) )
       IF(ISWITCH.EQ.1) THEN
          YCLtmp(1:NCL2)=YCC(1:NCL2)
       ELSE IF(ISWITCH.EQ.2) THEN
          YCLtmp(1:NCL2)=RM(1:NCL2)
       ELSE
          CALL ERRHDL('NO SUCH GEO',myid)
       END IF    
       
       IF (RANDOMTYPE==1) THEN
       XRDM = 0.0_WP  
     
       DO J=1,N2DO(MYID)              
          JJ=JCL2G(J)              

          VPERQ=VPER                               
          IF((1.0_WP-DABS(YCLtmp(JJ))).LT.0.250_WP) VPERQ=VPER*SVPER
          
          DO I=1,NCL1
             DO K=1,NCL3             

                CALL RANDOM_NUMBER(XRDM)            
                XRD(1)=(-1.0_WP+2.0_WP*XRDM(1))
                XRD(2)=(-1.0_WP+2.0_WP*XRDM(2))
                XRD(3)=(-1.0_WP+2.0_WP*XRDM(3))           
                Q(I,J,K,1)=VPERQ*XRD(1)
                Q(I,J,K,2)=VPERQ*XRD(2)*RC(JJ)
                Q(I,J,K,3)=VPERQ*XRD(3)*RM(JJ)
                            
                PR(I,J,K)=1.0_WP                
             ENDDO
          ENDDO
       ENDDO
       END IF
       
       
       IF (RANDOMTYPE==2) THEN
       XRDM = 0.0_WP  
       RDCOUNT = 0
       DO J=1,N2DO(MYID)          
          JJ=JCL2G(J)              
          VPERQ=VPER                               
          IF((1.0_WP-DABS(YCLtmp(JJ))).LT.0.250_WP) VPERQ=VPER*SVPER
          
          DO I=1,NCL1
             DO K=1,NCL3
                RDCOUNT = RDCOUNT + 1
                seed = 1973+RDCOUNT          
                call random_initialize ( seed )
                call rvec_random ( -1.0_WP, 1.0_WP, 3, XRD )                
                Q(I,J,K,1)=VPERQ*XRD(1)
                Q(I,J,K,2)=VPERQ*XRD(2)*RC(JJ)
                Q(I,J,K,3)=VPERQ*XRD(3)*RM(JJ)
           
                PR(I,J,K)=1.0_WP                
             ENDDO
          ENDDO
       ENDDO
       END IF
       
       DEALLOCATE ( YCLtmp )
 
       RETURN
       END SUBROUTINE INIFIELD_FLOW_tg
       
     
       SUBROUTINE INIFIELD_FLOW_io  
       use init_info
       use mesh_info
       use flow_info
       IMPLICIT NONE      
       
       
       INTEGER(4) :: I
       INTEGER(4) :: J, JJ
       INTEGER(4) :: K
       INTEGER(4) :: RDCOUNT, seed
!       INTEGER(4) :: IDUM
       REAL(WP)    :: VPERQ
       REAL(WP)    :: XRDM(3)
       REAL(WP)    :: XRD(3)
       REAL(WP), ALLOCATABLE    :: YCLtmp(:)
       
       INTEGER(4)  :: RANDOMTYPE=2  !=1 for real random, =2 for fixed random
       
       ALLOCATE ( YCLtmp(NCL2) )
       IF(ISWITCH.EQ.1) THEN
          YCLtmp(1:NCL2)=YCC(1:NCL2)
       ELSE IF(ISWITCH.EQ.2) THEN
          YCLtmp(1:NCL2)=RM(1:NCL2)
       ELSE
          CALL ERRHDL('NO SUCH GEO',myid)
       END IF    
       
       IF (RANDOMTYPE==1) THEN
       XRDM = 0.0_WP  
   
       DO J=1,N2DO(MYID)               
          JJ=JCL2G(J)              

          VPERQ=VPER                               
          IF((1.0_WP-DABS(YCLtmp(JJ))).LT.0.250_WP) VPERQ=VPER*SVPER
          
          DO I=1,NCL1_io
             DO K=1,NCL3             

                CALL RANDOM_NUMBER(XRDM)
             
                XRD(1)=(-1.0_WP+2.0_WP*XRDM(1))
                XRD(2)=(-1.0_WP+2.0_WP*XRDM(2))
                XRD(3)=(-1.0_WP+2.0_WP*XRDM(3))
              
                Q_io(I,J,K,1)=VPERQ*XRD(1)
                Q_io(I,J,K,2)=VPERQ*XRD(2)*RC(JJ)
                Q_io(I,J,K,3)=VPERQ*XRD(3)*RM(JJ)
                            
                PR_io(I,J,K)=1.0_WP                
             ENDDO
          ENDDO
       ENDDO
       END IF
!******************************************************************
       
       IF (RANDOMTYPE==2) THEN
       
       XRD = 0.0_WP  
       RDCOUNT = 0
   
       DO J=1,N2DO(MYID)     
          JJ=JCL2G(J)              
          VPERQ=VPER                               
          IF((1.0_WP-DABS(YCLtmp(JJ))).LT.0.250_WP) VPERQ=VPER*SVPER
          
          DO I=1,NCL1_io
             DO K=1,NCL3
                RDCOUNT = RDCOUNT + 1
                seed = 1973+RDCOUNT+2048
                call random_initialize ( seed )
                call rvec_random ( -1.0_WP, 1.0_WP, 3, XRD )
               
                Q_io(I,J,K,1)=VPERQ*XRD(1)
                Q_io(I,J,K,2)=VPERQ*XRD(2)*RC(JJ)
                Q_io(I,J,K,3)=VPERQ*XRD(3)*RM(JJ)
            
                PR_io(I,J,K)=1.0_WP                
             ENDDO
          ENDDO
       ENDDO
       

       CALL BC_TINLET
       XRDM = 0.0_WP  
       RDCOUNT = 0
       DO J=1,N2DO(MYID)
           JJ=JCL2G(J) 
           VPERQ=VPER  
           IF((1.0_WP-DABS(YCLtmp(JJ))).LT.0.250_WP) VPERQ=VPER*SVPER
           DO K=1,NCL3
                RDCOUNT = RDCOUNT + 1
                seed = 1973+RDCOUNT+1024
                call random_initialize ( seed )
                call rvec_random ( -1.0_WP, 1.0_WP, 3, XRD )
              
                Q_io(NCL1_io+1,J,K,1)=VPERQ*XRD(1)
                Q_io(NCL1_io+1,J,K,2)=VPERQ*XRD(2)*RC(JJ)
                Q_io(NCL1_io+1,J,K,3)=VPERQ*XRD(3)*RM(JJ)
            
                PR_io(NCL1_io+1,J,K)=1.0_WP  
           END DO
       END DO
       
       
       END IF
       
       DEALLOCATE ( YCLtmp )
 
       RETURN
       END SUBROUTINE INIFIELD_FLOW_io

