      SUBROUTINE POSTPROCESS_tg
       use init_info
       use postprocess_info
       IMPLICIT NONE
       
       
          CALL DIVGCK
       
          IF (ITERG.LT.2.OR.DMOD(phyTIME,TSCN).LT.DT) CALL WRTSCNDAT
          
          IF ( (DMOD(phyTIME,TTECCK).LT.DT) .and. (.not.IOFLOWflg)) THEN
             CALL TEC360_WRITE
          END IF

          IF (DMOD(phyTIME,TSAVE).LT.DT.OR.DMOD(phyTIME,TSTAV1).LT.DT) THEN
             CALL WRTQP_tg
          ENDIF
          CALL MPI_BARRIER(ICOMM,IERROR)
       
       RETURN
      END SUBROUTINE POSTPROCESS_tg

      SUBROUTINE POSTPROCESS_io
       use init_info
       use postprocess_info
       IMPLICIT NONE
       
       
          CALL DIVGCK_io
       
          IF (ITERG.LT.2.OR.DMOD(phyTIME,TSCN).LT.DT) THEN
              CALL WRTSCNDAT_io
          END IF
          
          IF ( DMOD(phyTIME,TTECCK).LT.DT ) THEN
             CALL TEC360_WRITE
             CALL PP_WALLSTRESS_io
          END IF
          

          IF ( DMOD(phyTIME,TSAVE).LT.DT.OR.DMOD(phyTIME,TSTAV1).LT.DT) THEN
             CALL WRTQP_io
          ENDIF
          CALL MPI_BARRIER(ICOMM,IERROR)
       
       RETURN
      END SUBROUTINE POSTPROCESS_io
      
!*******************************************************************


