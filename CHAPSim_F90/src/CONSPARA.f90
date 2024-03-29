
       SUBROUTINE CONSPARA
       use init_info
       use mesh_info
       IMPLICIT NONE
       
       INTEGER(4) :: J
       
       
       PI=2.0_WP*(DASIN(1.0_WP))

       VL1313   =1.0_WP/DBLE(NCL1*NCL3)

       ALX1=HX
       ALX2=1.0_WP  
       if(iswitch.eq.1) then
          ALX3=HZ
       else 
          ALX3=2.0_WP*PI
       endif

       DX=  ALX1/DBLE(NCL1)
       DZ=  ALX3/DBLE(NCL3)
      
       DXI=1.0_WP/DX
       DZI=1.0_WP/DZ
      
       DXQI=DXI*DXI
       DZQI=DZI*DZI
       
       QDX1=DXI*0.250_WP
       QDX3=DZI*0.250_WP
                 
       CVISC=1.0_WP/REN        
       
       IF(IOFLOWflg) THEN
          VL1313_io=1.0_WP/DBLE(NCL1_io*NCL3)
       END IF          

       
       RETURN
       END SUBROUTINE CONSPARA
