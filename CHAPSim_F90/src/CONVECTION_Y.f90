
      SUBROUTINE CONVECTION_Y_tg
      use flow_info
      use mesh_info
      use init_info
      IMPLICIT NONE    
      
      INTEGER(4)  :: IC, IM, IP
      INTEGER(4)  :: JC, JM, JP, JJ, JJM, JJP
      INTEGER(4)  :: KC, KM, KP
      INTEGER(4)  :: imsy
      INTEGER(4)  :: NII
      REAL(WP)     :: QDX2
      REAL(WP)     :: H21, H22, H23    
      REAL(WP)     :: q2s1     
      REAL(WP)     :: h23n
      REAL(WP)     :: q1e, q1w     
      REAL(WP)     :: d11q2e
      
      DPH(:,:,:)=0.0_WP 
     
      NII=1
      IF (MYID.EQ.0) THEN
         NII=2
      ENDIF
     
      DO KC=1,NCL3
         KM=KMV(KC)
         KP=KPV(KC)
         DO JC=NII,N2DO(MYID)
            JM=JC-1
            JP=JC+1
            JJ=JCL2G(JC)
            QDX2=DYCI(JJ)*0.250_WP
            jjm=jmv(jj)
            jjp=jpv(jj)
            DO IC=1,NCL1
       
               IP=IPV(IC)
               IM=IMV(IC)
     
               H21=((Q(IP,JC,KC,1)+Q(IP,JM,KC,1))*                     &
                    (Q(IP,JC,KC,2)+Q(IC,JC,KC,2))-                     &
                    (Q(IC,JC,KC,1)+Q(IC,JM,KC,1))*                     &
                    (Q(IC,JC,KC,2)+Q(IM,JC,KC,2)) )*QDX1
               if (iswitch.eq.2.and.myid.eq.0.and.jc.eq.2) then
                  imsy = isym(kmv(kc))
                  q2s1 = (q(ic,2,kc,2) - q(ic,2,isym(kc),2))*0.50_WP/rc(2)
                  h22=( (q(ic,3,kc,2)/rc(3)+q(ic,2,kc,2)/rc(2))*       &
                        (q(ic,3,kc,2)+q(ic,2,kc,2)) -                  &
                        (q(ic,2,kc,2)/rc(2)+q2s1)*                     &
                        (q(ic,2,kc,2)+q(ic,1,kc,2))                    &
                       )*QDX2  
               else
                  H22=((Q(IC,JP,KC,2)/rc(jjp)+Q(IC,JC,KC,2)/rc(jj))*   &
                       (Q(IC,JP,KC,2)+Q(IC,JC,KC,2))-                  &
                       (Q(IC,JC,KC,2)/rc(jj)+Q(IC,JM,KC,2)/rc(jjm))*   &
                       (Q(IC,JC,KC,2)+Q(IC,JM,KC,2)) )*QDX2
              endif
              H23=( ( Q(IC,JC,KP,3)/rm(jj)+      &
                      Q(IC,JM,KP,3)/rm(jjm) )*   &
                    ( Q(IC,JC,KP,2)+Q(IC,JC,KC,2) )- &
                    ( Q(IC,JC,KC,3)/rm(jj)+      &
                      Q(IC,JM,KC,3)/rm(jjm) )*   &
                    ( Q(IC,JC,KC,2)+Q(IC,JC,KM,2) ) )*QDX3/rc(jj)        !@
              DPH(IC,JC,KC)=-(H21+H22+H23)   
            ENDDO
         ENDDO
      ENDDO


      if (iswitch.eq.2) then        
         do kc=1,NCL3
            km=kmv(kc)
            kp=kmv(kc)
            do jc=1,n2do(MYID)
               jm=jc-1
               jj=JCL2G(JC)
               jjm=jmv(jj)
               do ic=1,NCL1
                  q1e=q(ic,jc,kp,3)/rm(jj)+q(ic,jm,kp,3)/rm(jjm) 
                  q1w=q(ic,jc,kc,3)/rm(jj)+q(ic,jm,kc,3)/rm(jjm)
                  h23n=( (q1e+q1w)*.250_WP )**2
                  d11q2e= -(q1e-q1w)*DZI/rc(jj)
                  DPH(ic,jc,kc)=DPH(ic,jc,kc)+h23n+d11q2e/ren
               enddo
            enddo
         enddo
      endif
        
      RETURN
      END SUBROUTINE CONVECTION_Y_tg
      
      
      SUBROUTINE CONVECTION_Y_io
      use flow_info
      use mesh_info
      use init_info
      IMPLICIT NONE    
      
      INTEGER(4)  :: IC, IM, IP
      INTEGER(4)  :: JC, JM, JP, JJ, JJM, JJP
      INTEGER(4)  :: KC, KM, KP
      INTEGER(4)  :: imsy
      INTEGER(4)  :: NII
      REAL(WP)     :: QDX2
      REAL(WP)     :: H21, H22, H23    
      REAL(WP)     :: q2s1     
      REAL(WP)     :: h23n
      REAL(WP)     :: q1e, q1w     
      REAL(WP)     :: d11q2e
      
      DPH_io(:,:,:)=0.0_WP  

      NII=1
      IF (MYID.EQ.0) THEN
         NII=2
      ENDIF
     
      DO KC=1,NCL3
         KM=KMV(KC)
         KP=KPV(KC)
         DO JC=NII,N2DO(MYID)
            JM=JC-1
            JP=JC+1
            JJ=JCL2G(JC)
            QDX2=DYCI(JJ)*0.250_WP
            jjm=jmv(jj)
            jjp=jpv(jj)
            DO IC=1,NCL1_io
       
               IP=IPV_io(IC)
               IM=IMV_io(IC)
     
               H21=((Q_io(IP,JC,KC,1)+Q_io(IP,JM,KC,1))*                     &
                    (Q_io(IP,JC,KC,2)+Q_io(IC,JC,KC,2))-                     &
                    (Q_io(IC,JC,KC,1)+Q_io(IC,JM,KC,1))*                     &
                    (Q_io(IC,JC,KC,2)+Q_io(IM,JC,KC,2)) )*QDX1
               if (iswitch.eq.2.and.myid.eq.0.and.jc.eq.2) then
                  imsy = isym(kmv(kc))
                  q2s1 = (q_io(ic,2,kc,2) - q_io(ic,2,isym(kc),2))*0.50_WP/rc(2)
                  h22=( (q_io(ic,3,kc,2)/rc(3)+q_io(ic,2,kc,2)/rc(2))*       &
                        (q_io(ic,3,kc,2)+q_io(ic,2,kc,2)) -                  &
                        (q_io(ic,2,kc,2)/rc(2)+q2s1)*                     &
                        (q_io(ic,2,kc,2)+q_io(ic,1,kc,2))                    &
                       )*QDX2  
               else
                  H22=((Q_io(IC,JP,KC,2)/rc(jjp)+Q_io(IC,JC,KC,2)/rc(jj))*   &
                       (Q_io(IC,JP,KC,2)+Q_io(IC,JC,KC,2))-                  &
                       (Q_io(IC,JC,KC,2)/rc(jj)+Q_io(IC,JM,KC,2)/rc(jjm))*   &
                       (Q_io(IC,JC,KC,2)+Q_io(IC,JM,KC,2)) )*QDX2
              endif
              H23=( ( Q_io(IC,JC,KP,3)/rm(jj)+      &
                      Q_io(IC,JM,KP,3)/rm(jjm) )*   &
                    ( Q_io(IC,JC,KP,2)+Q_io(IC,JC,KC,2) )- &
                    ( Q_io(IC,JC,KC,3)/rm(jj)+      &
                      Q_io(IC,JM,KC,3)/rm(jjm) )*   &
                    ( Q_io(IC,JC,KC,2)+Q_io(IC,JC,KM,2) ) )*QDX3/rc(jj)        !@
              DPH_io(IC,JC,KC)=-(H21+H22+H23)   
            ENDDO
         ENDDO
      ENDDO


      if (iswitch.eq.2) then        
         do kc=1,NCL3
            km=kmv(kc)
            kp=kmv(kc)
            do jc=1,n2do(MYID)
               jm=jc-1
               jj=JCL2G(JC)
               jjm=jmv(jj)
               do ic=1,NCL1_io
                  q1e=q_io(ic,jc,kp,3)/rm(jj)+q_io(ic,jm,kp,3)/rm(jjm) 
                  q1w=q_io(ic,jc,kc,3)/rm(jj)+q_io(ic,jm,kc,3)/rm(jjm)
                  h23n=( (q1e+q1w)*.250_WP )**2
                  d11q2e= -(q1e-q1w)*DZI/rc(jj)
                  DPH_io(ic,jc,kc)=DPH_io(ic,jc,kc)+h23n+d11q2e/ren
               enddo
            enddo
         enddo
      endif
        
      RETURN
      END SUBROUTINE CONVECTION_Y_io

