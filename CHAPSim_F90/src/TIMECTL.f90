
       SUBROUTINE TIMECTL
       use init_info
       IMPLICIT NONE
       
 
          IF (phyTIME.LT.TSTAV1) THEN
             DPCONS=RATEM1
             TSAVE=TSAVE1
          ELSE IF (phyTIME.GT.TSTAV1.AND.phyTIME.LE.tbody) THEN
             DPCONS=RATEM1
             TSAVE=TSAVE1*10
          ELSE IF (phyTIME.GT.tbody.AND.phyTIME.LE.tbody+1.0_WP) THEN
             DPCONS=RATEM1
             TSAVE=TSAVE1*2
          ELSE IF ((phyTIME.GT.(tbody+1.0_WP)) .AND.phyTIME.LE.TSTOP) THEN
             DPCONS=RATEM1
             TSAVE=TSAVE1
          ENDIF
          
          
     END SUBROUTINE TIMECTL
