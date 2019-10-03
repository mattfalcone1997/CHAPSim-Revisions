      SUBROUTINE PRCALC_tg(NS)
      use init_info
      use flow_info
      use mesh_info
      IMPLICIT NONE     

      INTEGER(4),INTENT(IN) :: NS
      REAL(WP) :: COE1
      REAL(WP) :: LLphi
      REAL(WP) :: ugmmm
      INTEGER(4) :: IC, IP, IM
      INTEGER(4) :: JC, JP, JM, JJ, jpp
      INTEGER(4) :: KC, KP, KM
     
      COE1=TALP(NS)*DT*0.50_WP*CVISC   
      
 
      IF(ISWITCH.EQ.1) THEN
      
      DO KC=1,NCL3
         KP=KPV(KC)
         KM=KMV(KC)
         DO JC=1,N2DO(MYID)
            JM=JC-1
            JP=JC+1
            JJ=JCL2G(JC)
            DO IC=1,NCL1                   
               LLphi = RHSLLPHI(IC,JC,KC)               
               PR(IC,JC,KC)=PR(IC,JC,KC)+DPH(IC,JC,KC)-COE1*LLphi
            ENDDO
         ENDDO
      ENDDO
     
      ELSE IF(ISWITCH.EQ.2) THEN
      
      DO KC=1,NCL3
         KP=KPV(KC)
         KM=KMV(KC)
         DO JC=1,N2DO(MYID)
            JM=JC-1
            JP=JC+1
            JJ=JCL2G(JC)
            jpp=jj+1
            ugmmm=DYFI(JJ)/RM(JJ)
            DO IC=1,NCL1
               IP=IPV(IC)
               IM=IMV(IC)               
               LLphi= DXQI*(DPH(IP,JC,KC)                         &
                        -2.0_WP*DPH(IC,JC,KC)                         &
                             +DPH(IM,JC,KC) )+                      &
                           ( (DPH(IC,JP,KC)                         &
                             -DPH(IC,JC,KC) )*rc(jpp)*DYCI(jpp)-     &
                           (  DPH(IC,JC,KC)                         &
                             -DPH(IC,JM,KC) )*rc(jj)*DYCI(jj))*ugmmm &              
                  +   DZQI*(DPH(IC,JC,KP)                         &
                        -2.0_WP*DPH(IC,JC,KC)                         &
                             +DPH(IC,JC,KM) )/(rm(jj)**2)   
               PR(IC,JC,KC)=PR(IC,JC,KC)+DPH(IC,JC,KC)-COE1*LLphi
            ENDDO
         ENDDO                                                                     
      ENDDO
      
      ELSE     
         CALL ERRHDL('NO SUCH GEO',myid)
      END IF
      
      RETURN
      END SUBROUTINE PRCALC_tg
      
      
      SUBROUTINE PRCALC_io(NS)
      use init_info
      use flow_info
      use mesh_info
      IMPLICIT NONE     

      INTEGER(4),INTENT(IN) :: NS  
      REAL(WP) :: COE1
      REAL(WP) :: LLphi
      REAL(WP) :: ugmmm
      INTEGER(4) :: IC, IP, IM
      INTEGER(4) :: JC, JP, JM, JJ, jpp
      INTEGER(4) :: KC, KP, KM   
      
      COE1=TALP(NS)*DT*0.50_WP*CVISC   
          
      IF(ISWITCH.EQ.1) THEN
      
      DO KC=1,NCL3
         KP=KPV(KC)
         KM=KMV(KC)
         DO JC=1,N2DO(MYID)
            JM=JC-1
            JP=JC+1
            JJ=JCL2G(JC)
            DO IC=1,NCL1_io            
               LLphi = RHSLLPHI_io(IC,JC,KC)               
               PR_io(IC,JC,KC)=PR_io(IC,JC,KC)+DPH_io(IC,JC,KC)-COE1*LLphi
            ENDDO
         ENDDO
      ENDDO
     
      ELSE IF(ISWITCH.EQ.2) THEN
      
      DO KC=1,NCL3
         KP=KPV(KC)
         KM=KMV(KC)
         DO JC=1,N2DO(MYID)
            JM=JC-1
            JP=JC+1
            JJ=JCL2G(JC)
            jpp=jj+1
            ugmmm=DYFI(JJ)/RM(JJ)
            DO IC=1,NCL1_io
               IP=IPV_io(IC)
               IM=IMV_io(IC)               
               LLphi= DXQI*(DPH_io(IP,JC,KC)                         &
                        -2.0_WP*DPH_io(IC,JC,KC)                         &
                             +DPH_io(IM,JC,KC) )+                      &
                           ( (DPH_io(IC,JP,KC)                         &
                             -DPH_io(IC,JC,KC) )*rc(jpp)*DYCI(jpp)-     &
                           (  DPH_io(IC,JC,KC)                         &
                             -DPH_io(IC,JM,KC) )*rc(jj)*DYCI(jj))*ugmmm &              
                  +   DZQI*(DPH_io(IC,JC,KP)                         &
                        -2.0_WP*DPH_io(IC,JC,KC)                         &
                             +DPH_io(IC,JC,KM) )/(rm(jj)**2)   
               PR_io(IC,JC,KC)=PR_io(IC,JC,KC)+DPH_io(IC,JC,KC)-COE1*LLphi
            ENDDO
         ENDDO                                                                     
      ENDDO
      
      ELSE     
         CALL ERRHDL('NO SUCH GEO',myid)
      END IF
      
      RETURN
      END SUBROUTINE PRCALC_io

