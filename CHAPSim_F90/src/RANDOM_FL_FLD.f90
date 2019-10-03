       SUBROUTINE RANDOM_FL_FLD_TG
       use init_info
       use mesh_info
       use flow_info
       IMPLICIT NONE      
             
       INTEGER(4) :: I
       INTEGER(4) :: J, JJ
       INTEGER(4) :: K
       INTEGER(4) :: L
       
       IF(MYID.EQ.0) CALL CHKHDL('   (1-1) TG:Generating random velocity',myid)
       CALL INIFIELD_FLOW_tg        
       CALL VMAV_tg                  
       CALL INTFC_BC_QP_tg  
       IF(MYID.EQ.0) THEN
          CALL CHKHDL('         TG: Random velocity perturbation VMV(1:3)',myid)      
          WRITE(*,'(A,20X,3ES15.7)') '#',(VMV(L),L=1,3)
       END IF

       DO L=1,3
          CALL MEANVELOINI_tg(L)
          DO J=1,N2DO(MYID)
             JJ=JCL2G(J) 
             Q(:,J,:,L)= Qtmp(:,J,:) -STA13(1,L,JJ)
          ENDDO
       END DO   
       CALL VMAV_tg        
       CALL INTFC_BC_QP_tg
       IF(MYID.EQ.0) THEN
          CALL CHKHDL('   (2-1) TG:Random velocity perturbation after subtracting meanXZ VMV(1:3)',myid) 
          WRITE(*,'(A,20X,3ES15.7)') '#',(VMV(L),L=1,3)
       END IF

       DO J=1,N2DO(MYID)
          JJ=JCL2G(J)
          DO I=1,NCL1
             DO K=1,NCL3
                Q(I,J,K,NFLOW)= Q(I,J,K,NFLOW) + Vini(JJ)  
             ENDDO
          ENDDO
       ENDDO             
        
  
       DO J=1,N2DO(MYID)
          JJ=JCL2G(J)
          IF ( (JJ.EQ.1) .OR. (JJ.EQ.(NCL2+1)) ) THEN
              Q(:,J,:,2)= 0.0_WP
          ENDIF
       ENDDO        
       CALL WRTVPROF_tg
             
       CALL VMAV_tg 
       CALL INTFC_BC_QP_tg
       
       IF(MYID.EQ.0) THEN
          CALL CHKHDL('   (3-1) TG:Initial velocity with given profiles VMV(1:3)',myid) 
          WRITE(*,'(A,20X,3ES15.7)') '#',(VMV(L),L=1,3)
       END IF

       CALL DIVGCK    
       IF(MYID.EQ.0) THEN
          CALL CHKHDL('   (4-1) TG:Max divergence of initial flow field',myid) 
          WRITE(*,'(A,20X,1ES15.7)') '#',MAXDIVGV
       END IF
      
      RETURN
      
      END SUBROUTINE RANDOM_FL_FLD_TG
      
!*********************************************************************************************************************
      SUBROUTINE CALC_INITIALIZATION_tg
       use init_info
       use mesh_info
       use flow_info
       IMPLICIT NONE      
             
       INTEGER(4) :: I
       INTEGER(4) :: J, JJ
       INTEGER(4) :: K
       INTEGER(4) :: L
      

       IF(MYID.EQ.0) CALL CHKHDL('   (5-1) TG: Initial flow field, calc QP...',myid)       
       CALL DIVG_tg(0)      
       CALL FFT99_POIS3D_periodicxz
       CALL INTFC_BC_DPH_tg
       CALL VELOUPDT_tg(0)
       CALL PRCALC_tg(0)
       CALL INTFC_BC_QP_tg
       CALL VMAV_tg
       IF(MYID.EQ.0) THEN
          CALL CHKHDL('   (6-1) TG: Initial flow field with updated velocity, VMAX(1:3)',myid) 
          WRITE(*,'(A,20X,3ES15.7)') '#',(VMV(L),L=1,3)
       END IF
       CALL DIVGCK
       IF(MYID.EQ.0) THEN
          CALL CHKHDL('   (7-1) TG: MAX DIVERGENCE OF VELOCITY in final initial flow field',myid) 
          WRITE(*,'(A,20X,1ES15.7)') '#',MAXDIVGV
       END IF

       RETURN
       END SUBROUTINE CALC_INITIALIZATION_tg
       

       SUBROUTINE RANDOM_FL_FLD_io
       use init_info
       use mesh_info
       use flow_info
       IMPLICIT NONE      
             
       INTEGER(4) :: I
       INTEGER(4) :: J, JJ
       INTEGER(4) :: K
       INTEGER(4) :: L
 
       IF(.not.IOFLOWflg) RETURN
 

          IF(MYID.EQ.0) CALL CHKHDL('   (1-2) IO:Generating random velocity',myid)
          CALL INIFIELD_FLOW_io         
          CALL VMAV_io                
          CALL INTFC_BC_QP_io  
          IF(MYID.EQ.0) THEN
             CALL CHKHDL('         IO: Random velocity perturbation VMV(1:3)',myid)      
             WRITE(*,'(A,20X,3ES15.7)') '#',(VMV_io(L),L=1,3)
          END IF

          DO L=1,3
             CALL MEANVELOINI_io(L)
             DO J=1,N2DO(MYID)
                JJ=JCL2G(J)
                Q_io(:,J,:,L)= Qtmp_io(:,J,:) -STA13_io(1,L,JJ)
             ENDDO
          END DO      
          CALL VMAV_io          
          CALL INTFC_BC_QP_io
          IF(MYID.EQ.0) THEN
             CALL CHKHDL('   (2-2) IO:Random velocity perturbation after subtracting meanXZ VMV(1:3)',myid) 
             WRITE(*,'(A,20X,3ES15.7)') '#',(VMV_io(L),L=1,3)
          END IF
            
          DO J=1,N2DO(MYID)
             JJ=JCL2G(J)
             DO I=0,NCL1_io+1,1
                DO K=1,NCL3
                   Q_io(I,J,K,NFLOW)= Q_io(I,J,K,NFLOW) + Vini(JJ)  
                ENDDO
             ENDDO
          ENDDO          
          
    
          DO J=1,N2DO(MYID)
             JJ=JCL2G(J)
             IF ( (JJ.EQ.1) .OR. (JJ.EQ.(NCL2+1)) ) THEN
                 Q_io(:,J,:,2)= 0.0_WP
             ENDIF
          ENDDO        
          CALL WRTVPROF_io
              
          CALL VMAV_io
          CALL INTFC_BC_QP_io
          IF(MYID.EQ.0) THEN
             CALL CHKHDL('   (3-2) IO:INITIAL VELOCITY WITH PERB AND GIVEN PROFILE VMV(1:3)',myid) 
             WRITE(*,'(A,20X,3ES15.7)') '#',(VMV_io(L),L=1,3)
          END IF
          
          CALL DIVGCK_io    
          IF(MYID.EQ.0) THEN
             CALL CHKHDL('   (4-2) IO:MAX DIVERGENCE OF VELOCITY, Main field, INOUT B.C.',myid) 
             WRITE(*,'(A,20X,3ES15.7)') '#',MAXDIVGV_io(1),MAXDIVGV_io(2),MAXDIVGV_io(3)
          END IF

      
      RETURN
      
      END SUBROUTINE RANDOM_FL_FLD_io
      
!*********************************************************************************************************************
      SUBROUTINE CALC_INITIALIZATION_io
       use init_info
       use mesh_info
       use flow_info
       IMPLICIT NONE      
             
       INTEGER(4) :: I
       INTEGER(4) :: J, JJ
       INTEGER(4) :: K
       INTEGER(4) :: L
       
       IF(.not.IOFLOWflg) RETURN
      
          IF(MYID.EQ.0) CALL CHKHDL('   (5-2) IO:Initial flow field, calc QP...',myid) 
          
          CALL BC_TINLET
          CALL BC_COUTLET_MOM_RK3(0)
          CALL BC_CBC_TDMA(0)
          CALL INTFC_BC_QP_INOUT
           
          CALL DIVG_io(0)
          CALL FISHPACK_POIS3D_SIMPLE
          CALL INTFC_BC_DPH_io
          CALL VELOUPDT_io(0)
          CALL PRCALC_io(0)
          CALL INTFC_BC_QP_io
          CALL VELOUPDT_CBC_U
          CALL VMAV_io
          IF(MYID.EQ.0) THEN
             CALL CHKHDL('   (6-2) IO: Initial flow field with updated velocity, VMAX(1:3)',myid) 
             WRITE(*,'(A,20X,3ES15.7)') '#',(VMV_io(L),L=1,3)
          END IF
          CALL DIVGCK_io
          IF(MYID.EQ.0) THEN
             CALL CHKHDL('   (7-2) IO: Max divergence of the flow field. Main field, INOUT B.C.',myid) 
             WRITE(*,'(A,20X,3ES15.7)') '#',MAXDIVGV_io(1),MAXDIVGV_io(2),MAXDIVGV_io(3)
          END IF

       RETURN
       END SUBROUTINE CALC_INITIALIZATION_io

