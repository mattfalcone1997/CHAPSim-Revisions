
      SUBROUTINE LAPLACECOEF
      use mesh_info
      IMPLICIT NONE     

      INTEGER(4) :: JC, JM, JP
      INTEGER(4) :: IDR
      INTEGER(4) :: NFIL
      REAL(WP)    :: UCAJ 
     
      DO JC=1,NCL2      
         JP=JC+1
         AMPH(JC)= DYFI(JC)*DYCI(JC)
         APPH(JC)= DYFI(JC)*DYCI(JP)
      ENDDO
      

      ACPH(1:NCL2)= -( AMPH(1:NCL2)+APPH(1:NCL2) )
      
      ACPH(1)    = AMPH(1) + ACPH(1)
      AMPH(1)    = 0.0_WP
      
      ACPH(NCL2) = APPH(NCL2) + ACPH(NCL2)
      APPH(NCL2) = 0.0_WP
      

      DO IDR=1,3
         IF( (IDR.EQ.1) .OR. (IDR.EQ.3) ) THEN

            DO JC=2,NCL2-1
               JP=JC+1
               AMVR(JC,IDR)= DYFI(JC)*DYCI(JC)
               APVR(JC,IDR)= DYFI(JC)*DYCI(JP)
               ACVR(JC,IDR)=-( AMVR(JC,IDR) + APVR(JC,IDR) )
            ENDDO

            UCAJ=4.0_WP/( 2.0_WP/DYCI(2)+1.0_WP/DYFI(1) )
            AMVR(1,IDR)=   0.0_WP
            APVR(1,IDR)=   UCAJ*DYCI(2)
            ACVR(1,IDR)=  -UCAJ*( DYCI(2)+2.0_WP*DYCI(1) )

            UCAJ=4.0_WP/( 2.0_WP/DYCI(NCL2)+1.0_WP/DYFI(NCL2) )
            AMVR(NCL2,IDR)= UCAJ*DYCI(NCL2)
            APVR(NCL2,IDR)= 0.0_WP
            ACVR(NCL2,IDR)=-UCAJ*( 2.0_WP*DYCI(NND2)+DYCI(NCL2) )

         ELSE

            DO JC=2,NCL2
               JM=JC-1
               AMVR(JC,IDR)=  DYFI(JM)*DYCI(JC)
               APVR(JC,IDR)=  DYFI(JC)*DYCI(JC)
               ACVR(JC,IDR)= -( AMVR(JC,IDR)+APVR(JC,IDR) )
            ENDDO
            AMVR(1,IDR)=0.0_WP
            APVR(1,IDR)=0.0_WP
            ACVR(1,IDR)=1.0_WP
            AMVR(NND2,IDR)=0.0_WP
            APVR(NND2,IDR)=0.0_WP
            ACVR(NND2,IDR)=1.0_WP
            
         ENDIF
         
      ENDDO
      

      NFIL=16
      OPEN(NFIL,FILE='PHICOEF.dat')
      REWIND NFIL
      WRITE(NFIL,'(A)') &
               'AMPH(JC),  AMVR(JC,1), AMVR(JC,3), AMVR(JC,2),',  &
               'APPH(JC),  APVR(JC,1), APVR(JC,3), APVR(JC,2),',  &
               'ACPH(JC),  ACVR(JC,1), ACVR(JC,3), ACVR(JC,2)'
      DO JC=1,NCL2
          WRITE(NFIL,'(12ES13.5)') &
                AMPH(JC),  AMVR(JC,1), AMVR(JC,3), AMVR(JC,2),  &
                APPH(JC),  APVR(JC,1), APVR(JC,3), APVR(JC,2),  &
                ACPH(JC),  ACVR(JC,1), ACVR(JC,3), ACVR(JC,2)
      END DO
      CLOSE(NFIL)
  
      RETURN
      
      END SUBROUTINE LAPLACECOEF

!**********************************************************************
