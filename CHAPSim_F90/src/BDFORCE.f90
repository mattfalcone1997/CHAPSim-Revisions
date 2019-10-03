       SUBROUTINE BDFORCE
       use init_info
       IMPLICIT NONE
       
          if (dmod(phytime,(tbody-dt)).lt.dt) then
             call BDFORCECONFIG 
          endif
          if (dmod(phytime,(tbody+dt)).lt.dt) then 
             call BDFORCECONFIG 
          endif
          
       END SUBROUTINE BDFORCE
