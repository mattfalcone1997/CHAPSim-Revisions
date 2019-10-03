
      SUBROUTINE LAPLACECOEFpipe
      use mesh_info
      IMPLICIT NONE      

      INTEGER(4) :: JC, JM, JP
      INTEGER(4) :: NQ
      INTEGER(4) :: NFIL
      REAL(WP)    :: A22ICC
      REAL(WP)    :: A22ICP
      REAL(WP)    :: A22
      REAL(WP)    :: A22DP
      REAL(WP)    :: UCAJ11
      REAL(WP)    :: UGMMM
      REAL(WP)    :: AC2
!*********************************************************************        
       do jc=2,NCL2-1
          jp=jc+1
          A22ICC=rc(jc)*DYCI(jc)  
          A22ICP=rc(jp)*DYCI(jp)   
          AC2=-(A22ICC+A22ICP)
          UGMMM=rm(jc)*DYFI(JC)           
          AMPH(JC)=(A22ICC)*UGMMM 
          APPH(JC)=(A22ICP)*UGMMM 
          ACPH(JC)=AC2*UGMMM 
      enddo
      
      JC=1
      JP=JC+1
      A22ICP=rc(jp)*DBLE(jp-jc)*DYCI(jp)           
      UGMMM=rm(jc)*DYFI(JC)
      AMPH(JC)=0.0_WP
      APPH(JC)= (A22ICP)*UGMMM         
      ACPH(JC)=-(A22ICP)*UGMMM        
         
      JC=NCL2
      A22ICC=rc(jc)*DYCI(JC)
      UGMMM =rm(jc)*DYFI(JC)
      AMPH(JC)= (A22ICC)*UGMMM        
      APPH(JC)=0.0_WP
      ACPH(JC)=-(A22ICC)*UGMMM          
     
      do nq=1,3
         if(nq.eq.3) then                   
            do jc=2,NCL2-1
               jp=jc+1
               jm=jc-1     !@
               A22ICC=rc(jc)*DYCI(JC)        
               A22ICP=rc(jp)*DYCI(JP)         
               AMVR(JC,NQ)=  A22ICC*DYFI(JC)/rm(jm)            
               APVR(JC,NQ)=  A22ICP*DYFI(JC)/rm(jp)             
               ACVR(JC,NQ)=-(A22ICC+A22ICP)*DYFI(JC)/rm(jc)-1.0_WP/rm(jc)**2    
            enddo

            JC=1
            JP=JC+1                 
            UCAJ11=DYFI(JC)*rc(jp)*DYCI(JP)               
            AMVR(JC,NQ)=0.0_WP                   
            APVR(JC,NQ)= UCAJ11/rm(jp)               
            ACVR(JC,NQ)=-UCAJ11/rm(jc)-1.0_WP/rm(jc)**2       
!-----------
            JC=NCL2
            JP=NND2
            jm=jc-1 
            UGMMM=DYFI(JC)*rm(jc)       
                                       
            AMVR(JC,NQ)=UGMMM/rc(jc)*DYCI(jc)    
            APVR(JC,NQ)=0. !-UGMMM/rc(jp)/cac(jp)*2.0_WP                    
            ACVR(JC,NQ)=-UGMMM/rc(jc)*DYCI(jc)-UGMMM/rc(jp)*DYCI(jp)*2.0_WP    
         endif
         
         if(nq.eq.1) then
            do jc=2,NCL2-1
               jp=jc+1
               ugmmm=DYFI(JC)/rm(jc)        
               a22icc=+ugmmm*rc(jc)*DYCI(jc)        
               a22icp=+ugmmm*rc(jp)*DYCI(jp)         
               amvr(JC,NQ)=a22icc             
               apvr(JC,NQ)=a22icp             
               acvr(JC,NQ)=-(a22icc+a22icp)                    
           enddo
 
           jc=1
           jp=jc+1
           ugmmm=DYFI(jc)/rm(jc)            
           amvr(jc,nq)=0.0_WP                
           apvr(jc,nq)=ugmmm*(rc(jp)*DYCI(jp))        
           acvr(jc,nq)=-ugmmm*(rc(jp)*DYCI(jp))            
!-----------
           jc=NCL2
           jp=NND2
           ugmmm=DYFI(jc)/rm(jc)                                                                      
           amvr(jc,nq)=rc(jc)*DYCI(jc)*ugmmm      
           apvr(jc,nq)=0 !ugmmm*rc(jp)/cac(jp)*2.0_WP    
           acvr(jc,nq)=-ugmmm*rc(jc)*DYCI(jc)-ugmmm*rc(jp)*DYCI(jp)*2.0_WP            
         endif
         
         if(nq.eq.2) then
            DO jc=2,NCL2
               jm=jc-1
               a22  = DYCI(jc)
               a22dp= DYCI(jc)/(2.0_WP*rc(jc))
               apvr(jc,nq)=  a22*DYFI(jc)-a22dp
               acvr(jc,nq)=-(a22*DYFI(jc)+a22*DYFI(jm))
               amvr(jc,nq)=  a22*DYFI(jm)+a22dp
            ENDDO
                 
            do jc=1,NND2,NCL2
               apvr(jc,nq)=0.0_WP
               acvr(jc,nq)=0.0_WP
               amvr(jc,nq)=0.0_WP
            enddo
         endif
                                           
      enddo   ! end  nq   do

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
      
      END  SUBROUTINE LAPLACECOEFpipe
