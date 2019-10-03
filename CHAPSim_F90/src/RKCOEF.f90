       SUBROUTINE RKCOEF
       use init_info
       IMPLICIT NONE

       NSST=3
       IF (NSST.EQ.3) THEN

           TGAM(0)= 1.0_WP
           TGAM(1)= 8.0_WP/15.0_WP
           TGAM(2)= 5.0_WP/12.0_WP
           TGAM(3)= 3.0_WP/4.0_WP
           
           TROH(0)= 0.0_WP
           TROH(1)= 0.0_WP
           TROH(2)=-17.0_WP/60.0_WP
           TROH(3)=-5.0_WP/12.0_WP
           
           TALP(0)= 1.0_WP
           TALP(1)= TGAM(1) + TROH(1)
           TALP(2)= TGAM(2) + TROH(2)
           TALP(3)= TGAM(3) + TROH(3)
           
       ELSE

           TGAM(1)= 1.50_WP
           TGAM(2)= 0.0_WP
           TGAM(3)= 0.0_WP
           
           TROH(1)= -0.50_WP
           TROH(2)=  0.0_WP
           TROH(3)=  0.0_WP
           
           TALP(0)= 1.0_WP
           TALP(1)= TGAM(1) + TROH(1)
           TALP(2)= TGAM(2) + TROH(2)
           TALP(3)= TGAM(3) + TROH(3)
       ENDIF


       END SUBROUTINE RKCOEF
