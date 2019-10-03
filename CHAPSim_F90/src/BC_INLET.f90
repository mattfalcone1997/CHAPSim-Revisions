      SUBROUTINE BC_TINLET
      use mesh_info
      use flow_info
      IMPLICIT NONE
      
      INTEGER(4) :: J, K, IDR
      
      
      DO J = 0, N2DO(MYID)+1, 1
         DO K=1, NCL3
             DO IDR=1,NDV
                 Q_io  (0, J, K, IDR) = Q(NCL1, J, K , IDR) 
             END DO
             PR_io (0, J, K)    = PR(NCL1, J, K)
             DPH_io(0, J, K)    = 0.0_WP
             
             !Q_io  (1, J, K, 1) = Q(1, J, K , 1) 
             
         END DO
      END DO

      RETURN
      END SUBROUTINE BC_TINLET
