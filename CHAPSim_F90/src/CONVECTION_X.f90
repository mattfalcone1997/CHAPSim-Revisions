      SUBROUTINE CONVECTION_X_tg
      use mesh_info
      use flow_info
      IMPLICIT NONE    
      
      INTEGER(4)  :: IC, IM, IP
      INTEGER(4)  :: JC, JM, JP, JJ
      INTEGER(4)  :: KC, KM, KP
      REAL(WP)     :: QDX2
      REAL(WP)     :: H11, H12, H13
       
      Qtmp(:,:,:) = 0.0_WP

      DO KC=1,NCL3
         KM=KMV(KC)
         KP=KPV(KC)
         DO JC=1,N2DO(MYID)
            JM=JC-1
            JP=JC+1
            JJ=JCL2G(JC)
            QDX2 = DYFI(JJ) *0.250_WP 
            DO IC=1,NCL1
               IP=IPV(IC)
               IM=IMV(IC)
               H11=( (Q(IP,JC,KC,1)+Q(IC,JC,KC,1))*                  &
                     (Q(IP,JC,KC,1)+Q(IC,JC,KC,1))-                  &
                     (Q(IM,JC,KC,1)+Q(IC,JC,KC,1))*                  &
                     (Q(IM,JC,KC,1)+Q(IC,JC,KC,1)) )*QDX1
               H12=(((Q(IC,JP,KC,2)+Q(IM,JP,KC,2))*                  &
                     (Q(IC,JP,KC,1)+Q(IC,JC,KC,1))-                  &
                     (Q(IC,JC,KC,2)+Q(IM,JC,KC,2))*                  &
                     (Q(IC,JC,KC,1)+Q(IC,JM,KC,1)) )*QDX2 )/rm(jj)
                        
               H13=(((Q(IC,JC,KP,3)+Q(IM,JC,KP,3))*                   &
                     (Q(IC,JC,KP,1)+Q(IC,JC,KC,1))-                   &
                     (Q(IC,JC,KC,3)+Q(IM,JC,KC,3))*                   &
                     (Q(IC,JC,KC,1)+Q(IC,JC,KM,1)) )*QDX3/rm(jj) )/rm(jj)
                     
               Qtmp(IC,JC,KC)=-(H11+H12+H13)     
            ENDDO
         ENDDO
      ENDDO
        
      RETURN
      END SUBROUTINE CONVECTION_X_tg

      SUBROUTINE CONVECTION_X_io
      use mesh_info
      use flow_info
      IMPLICIT NONE    
      
      INTEGER(4)  :: IC, IM, IP
      INTEGER(4)  :: JC, JM, JP, JJ
      INTEGER(4)  :: KC, KM, KP
      REAL(WP)     :: QDX2 
      REAL(WP)     :: H11, H12, H13
       
      Qtmp_io(:,:,:)=0.0_WP

      DO KC=1,NCL3
         KM=KMV(KC)
         KP=KPV(KC)
         DO JC=1,N2DO(MYID)
            JM=JC-1
            JP=JC+1
            JJ=JCL2G(JC)
            QDX2 = DYFI(JJ) *0.250_WP 
            DO IC=1,NCL1_io
               IP=IPV_io(IC)
               IM=IMV_io(IC)
               H11=( (Q_io(IP,JC,KC,1)+Q_io(IC,JC,KC,1))*                  &
                     (Q_io(IP,JC,KC,1)+Q_io(IC,JC,KC,1))-                  &
                     (Q_io(IM,JC,KC,1)+Q_io(IC,JC,KC,1))*                  &
                     (Q_io(IM,JC,KC,1)+Q_io(IC,JC,KC,1)) )*QDX1
               H12=(((Q_io(IC,JP,KC,2)+Q_io(IM,JP,KC,2))*                  &
                     (Q_io(IC,JP,KC,1)+Q_io(IC,JC,KC,1))-                  &
                     (Q_io(IC,JC,KC,2)+Q_io(IM,JC,KC,2))*                  &
                     (Q_io(IC,JC,KC,1)+Q_io(IC,JM,KC,1)) )*QDX2 )/rm(jj)
                        
               H13=(((Q_io(IC,JC,KP,3)+Q_io(IM,JC,KP,3))*                   &
                     (Q_io(IC,JC,KP,1)+Q_io(IC,JC,KC,1))-                   &
                     (Q_io(IC,JC,KC,3)+Q_io(IM,JC,KC,3))*                   &
                     (Q_io(IC,JC,KC,1)+Q_io(IC,JC,KM,1)) )*QDX3/rm(jj) )/rm(jj)
                     
               Qtmp_io(IC,JC,KC)=-(H11+H12+H13) 
            ENDDO
         ENDDO
      ENDDO
        
      RETURN
      END SUBROUTINE CONVECTION_X_io
