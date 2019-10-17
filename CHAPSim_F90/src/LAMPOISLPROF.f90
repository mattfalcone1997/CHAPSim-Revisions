
      SUBROUTINE LAMPOISLPROF
      use init_info
      use mesh_info
      IMPLICIT NONE
      !=========================================================================
      !Not sure what this is doing
      !=========================================================================
      INTEGER(4) :: J



       if(iswitch.eq.1) then

          DO J=1,NCL2
             IF (YCC(J).LT.1.0_WP.AND.YCC(J).GT.-1.0_WP) THEN
                Vini(J)=(1.-YCC(J)**2) * (3.0_WP/2.0_WP)
             ELSE
                Vini(J)=0.0_WP
             ENDIF
          ENDDO

       else

          DO J=1,NCL2
             if (rm(j).lt.1.0_WP.and.rm(j).gt.0.0_WP) then
                 Vini(j)=(1.0_WP-rm(j)**2) * 2.0_WP
             else
                 Vini(j)=0.0_WP
             endif
          enddo

       endif


      RETURN

      END  SUBROUTINE LAMPOISLPROF
