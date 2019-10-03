      SUBROUTINE VALME13ORL_tg(LO,LQ)
      
      use mesh_info
      use flow_info
      use postprocess_info
      IMPLICIT NONE
      
      INTEGER(4),INTENT(IN) :: LO
      INTEGER(4),INTENT(IN) :: LQ
         
      REAL(WP)    :: UUJJ!, UUJJS
      REAL(WP)    :: VELENTER
      INTEGER(4) :: I, IP
      INTEGER(4) :: J, JP, JJ
      INTEGER(4) :: K, KP
      REAL(WP),ALLOCATABLE :: UU(:)
      REAL(WP),ALLOCATABLE :: WORK_UU(:)
      
      ALLOCATE ( UU(NCL2)      )
      ALLOCATE ( WORK_UU(NCL2) )
      UU      = 0.0_WP
      WORK_UU = 0.0_WP
       
             
       IF (LQ.EQ.1) THEN 
             
          DO J=1,N2DO(MYID)          
             JJ=JCL2G(J)            
             UUJJ    =0.0_WP
             VELENTER=0.0_WP            
             DO I=1,NCL1
                IP=IPV(I)        
                DO K=1,NCL3                                
                   VELENTER= ( Qtmp(IP,J,K) + Qtmp(I,J,K) ) * 0.50_WP     
                   UUJJ=VELENTER+UUJJ
                ENDDO
             ENDDO                 
             UU(JJ)  =UUJJ  *VL1313   
          ENDDO
          
       ELSE IF (LQ.EQ.3) THEN
       
          DO J=1,N2DO(MYID)         
             JJ=JCL2G(J)            
             UUJJ    =0.0_WP
             VELENTER=0.0_WP             
             DO I=1,NCL1       
                DO K=1,NCL3
                   KP=KPV(K)                              
                   VELENTER= ( Qtmp(I,J,KP) + Qtmp(I,J,K) ) * 0.50_WP     
                   UUJJ=VELENTER+UUJJ
                ENDDO
             ENDDO                  
             UU(JJ)  =UUJJ  *VL1313 
          ENDDO
       
       ELSE IF (LQ.EQ.2) THEN
      
          DO J=1,N2DO(MYID)          
             JJ=JCL2G(J)
             JP=J+1            
             UUJJ    =0.0_WP
             VELENTER=0.0_WP             
             DO I=1,NCL1      
                DO K=1,NCL3                               
                   VELENTER= ( Qtmp(I,JP,K) + Qtmp(I,J,K) ) * 0.50_WP     
                   UUJJ=VELENTER+UUJJ
                ENDDO
             ENDDO                  
             UU(JJ)  =UUJJ  *VL1313 
          ENDDO     
            
       ELSE
          CALL ERRHDL('DIRECTIONS SHOULD BE 1, 2 or 3, and no others.',myid)  
       END IF
       
 
       CALL MPI_BARRIER(ICOMM,IERROR)
       CALL MPI_ALLREDUCE(UU(1),  WORK_UU(1),  NCL2, MPI_DOUBLE_PRECISION,&
              MPI_SUM, ICOMM, IERROR)

       DO JJ=1,NCL2
          STA13(LO,LQ,JJ)= WORK_UU(JJ)  !VMNO0(JJ)
       ENDDO
       
      DEALLOCATE ( UU      )
      DEALLOCATE ( WORK_UU )
      
      RETURN
      
      END SUBROUTINE VALME13ORL_tg
      

      SUBROUTINE VALME13ORL_io(LO,LQ)
      
      use mesh_info
      use flow_info
      use postprocess_info
      IMPLICIT NONE
      
      INTEGER(4),INTENT(IN) :: LO
      INTEGER(4),INTENT(IN) :: LQ
         
      REAL(WP)    :: UUJJ!, UUJJS
      REAL(WP)    :: VELENTER
      INTEGER(4) :: I, IP
      INTEGER(4) :: J, JP, JJ
      INTEGER(4) :: K, KP
      REAL(WP),ALLOCATABLE :: UU(:)
      REAL(WP),ALLOCATABLE :: WORK_UU(:)
      
      ALLOCATE ( UU(NCL2)      )
      ALLOCATE ( WORK_UU(NCL2) )
      UU      = 0.0_WP
      WORK_UU = 0.0_WP
       
              
       IF (LQ.EQ.1) THEN 
             
          DO J=1,N2DO(MYID)          
             JJ=JCL2G(J)            
             UUJJ    =0.0_WP
             VELENTER=0.0_WP            
             DO I=1,NCL1_io
                IP=IPV_io(I)        
                DO K=1,NCL3                                
                   VELENTER= ( Qtmp_io(IP,J,K) + Qtmp_io(I,J,K) ) * 0.50_WP     
                   UUJJ=VELENTER+UUJJ
                ENDDO
             ENDDO                 
             UU(JJ)  =UUJJ  *VL1313_io 
          ENDDO
          
       ELSE IF (LQ.EQ.3) THEN
       
          DO J=1,N2DO(MYID)         
             JJ=JCL2G(J)            
             UUJJ    =0.0_WP
             VELENTER=0.0_WP             
             DO I=1,NCL1_io
                DO K=1,NCL3
                   KP=KPV(K)                                      
                   VELENTER= ( Qtmp_io(I,J,KP) + Qtmp_io(I,J,K) ) * 0.50_WP     
                   UUJJ=VELENTER+UUJJ
                ENDDO
             ENDDO                  
             UU(JJ)  =UUJJ  *VL1313_io 
          ENDDO
       
       ELSE IF (LQ.EQ.2) THEN
     
          DO J=1,N2DO(MYID)          
             JJ=JCL2G(J)
             JP=J+1            
             UUJJ    =0.0_WP
             VELENTER=0.0_WP             
             DO I=1,NCL1_io      
                DO K=1,NCL3                               
                   VELENTER= ( Qtmp_io(I,JP,K) + Qtmp_io(I,J,K) ) * 0.50_WP     
                   UUJJ=VELENTER+UUJJ
                ENDDO
             ENDDO                  
             UU(JJ)  =UUJJ  *VL1313_io 
          ENDDO     
            
       ELSE
          CALL ERRHDL('DIRECTIONS SHOULD BE 1, 2 or 3, and no others.',myid)  
       END IF
        
       CALL MPI_BARRIER(ICOMM,IERROR)
       CALL MPI_ALLREDUCE(UU(1),  WORK_UU(1),  NCL2, MPI_DOUBLE_PRECISION,&
              MPI_SUM, ICOMM, IERROR)

       DO JJ=1,NCL2
          STA13_io(LO,LQ,JJ)= WORK_UU(JJ) 
       ENDDO
       
      DEALLOCATE ( UU      )
      DEALLOCATE ( WORK_UU )
      
      RETURN
      
      END SUBROUTINE VALME13ORL_io

