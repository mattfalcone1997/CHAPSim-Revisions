
     SUBROUTINE FLOWSTART
       use init_info
       use flow_info 
       IMPLICIT NONE
       
       INTEGER(4) :: L
     
       IF (NREAD.EQ.0) THEN  ! 
          IF (MYID.EQ.0) CALL CHKHDL('12.TG:Flow initialization from random velocity field', myid)          
          CALL RANDOM_FL_FLD_TG
          CALL CALC_INITIALIZATION_tg
       ELSE IF (NREAD.EQ.1) THEN  !
          CALL ERRHDL('Removed this subroutine',myid)
          
       ELSE IF (NREAD.EQ.2) THEN     
          IF (MYID.EQ.0) CALL CHKHDL('12.TG:Flow initialization from last step using the same mesh.', myid) 
          CALL RESTART_tg 
                    
       ENDIF
       
       
       IF(IOFLOWflg) THEN
       
           IF (NREAD_io.EQ.0) THEN  ! 
              IF (MYID.EQ.0) CALL CHKHDL('13.IO:Flow initialization from random velocity field', myid)          
              CALL RANDOM_FL_FLD_io
              CALL CALC_INITIALIZATION_io
           ELSE IF (NREAD_io.EQ.1) THEN  !
              CALL ERRHDL('Removed this subroutine',myid)
              
           ELSE IF (NREAD_io.EQ.2) THEN     
              IF (MYID.EQ.0) CALL CHKHDL('13.IO:Flow initialization from last step using the same mesh.', myid) 
              CALL RESTART_io          
           ENDIF

       END IF

       CALL TEC360_WRITE
       

       CALL INTFC_BC_QP_tg
       IF(IOFLOWflg) CALL INTFC_BC_QP_io


       RETURN     
     END SUBROUTINE FLOWSTART
