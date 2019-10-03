      SUBROUTINE CONVECTION_Z_tg
      use init_info
      use mesh_info
      use flow_info
      IMPLICIT NONE    
      
      INTEGER(4)  :: IC, IM, IP
      INTEGER(4)  :: JC, JM, JP, JJ, JJM, JJP
      INTEGER(4)  :: KC, KM, KP
      INTEGER(4)  :: imsy
      REAL(WP)     :: QDX2 
      REAL(WP)     :: H31, H32, H33
      REAL(WP)     :: q2s1, q2s2      
      REAL(WP)     :: h32n
      REAL(WP)     :: q2e, q2w      
      REAL(WP)     :: d11q1e
      
      
      RHSLLPHI(:,:,:)=0.0_WP 
      
      DO KC=1,NCL3
         KM=KMV(KC)
         KP=KPV(KC)
         DO JC=1,N2DO(MYID)
            JM=JC-1
            JP=JC+1
            JJ=JCL2G(JC)
            QDX2=DYFI(JJ)*0.250_WP

            jjm=jmv(jj)           !@
            jjp=jpv(jj)           !@
            DO IC=1,NCL1
               IP=IPV(IC)
               IM=IMV(IC)                 
               H32=( ( Q(IC,JP,KC,2)+Q(IC,JP,KM,2) )*                &
                     ( Q(IC,JP,KC,3)/rm(jjp)+                        &
                       Q(IC,JC,KC,3)/rm(jj)   )-                     &       
                     ( Q(IC,JC,KC,2)+Q(IC,JC,KM,2) )*                &
                     ( Q(IC,JC,KC,3)/rm(jj)+                         &
                       Q(IC,JM,KC,3)/rm(jjm) ) )*QDX2   
                    
               H31=((Q(IP,JC,KC,1)+Q(IP,JC,KM,1))*                  &
                    (Q(IP,JC,KC,3)+Q(IC,JC,KC,3))-                  &
                    (Q(IC,JC,KC,1)+Q(IC,JC,KM,1))*                  &
                    (Q(IC,JC,KC,3)+Q(IM,JC,KC,3)) )*QDX1 
                    
               H33=((Q(IC,JC,KP,3)+Q(IC,JC,KC,3))*                   &
                    (Q(IC,JC,KP,3)+Q(IC,JC,KC,3))-                   &
                    (Q(IC,JC,KC,3)+Q(IC,JC,KM,3))*                   &
                    (Q(IC,JC,KC,3)+Q(IC,JC,KM,3)) )*QDX3/rm(jj)**2  

               RHSLLPHI(IC,JC,KC)=-(H31+H32+H33)
            ENDDO
         ENDDO
      ENDDO
            
      if (iswitch.eq.2) then 
          do kc=1,NCL3
             km=kmv(kc)
             kp=kpv(kc)
             DO jc=1,n2do(MYID)
                jm=jc-1
                jp=jc+1
                jj=JCL2G(JC)
                jjm=jmv(jj)           !@
                jjp=jpv(jj)           !@
                do ic=1,NCL1
                   if (jc.eq.1.and.myid.eq.0) then    
                      imsy = isym(kmv(kc))       
                      q2s1= (q(ic,jp,kc,2)-q(ic,jp,isym(kc),2))*0.50_WP/rc(jjp)
                      q2s2= (q(ic,jp,km,2)-q(ic,jp,imsy,    2))*0.50_WP/rc(jjp)
                      q2e =  q(ic,jp,kc,2)/rc(jjp) + q2s1
                      q2w =  q(ic,jp,km,2)/rc(jjp) + q2s2
                   else
                      q2e = q(ic,jp,kc,2)/rc(jjp)+q(ic,jc,kc,2)/rc(jj)
                      q2w = q(ic,jp,km,2)/rc(jjp)+q(ic,jc,km,2)/rc(jj)
                   end if
                   h32n= q(ic,jc,kc,3)*(q2e+q2w)*0.250_WP/rm(jj)
                   d11q1e=1.0_WP*(q2e-q2w)*DZI/rm(jj)
                   RHSLLPHI(ic,jc,kc)=-h32n+d11q1e/ren+RHSLLPHI(ic,jc,kc)              
                enddo
             enddo
          enddo 
     
      endif 
              
      RETURN
      END SUBROUTINE CONVECTION_Z_tg


      SUBROUTINE CONVECTION_Z_io
      use init_info
      use mesh_info
      use flow_info
      IMPLICIT NONE    
      
      INTEGER(4)  :: IC, IM, IP
      INTEGER(4)  :: JC, JM, JP, JJ, JJM, JJP
      INTEGER(4)  :: KC, KM, KP
      INTEGER(4)  :: imsy
      REAL(WP)     :: QDX2
      REAL(WP)     :: H31, H32, H33
      REAL(WP)     :: q2s1, q2s2      
      REAL(WP)     :: h32n
      REAL(WP)     :: q2e, q2w      
      REAL(WP)     :: d11q1e
      
      
      RHSLLPHI_io(:,:,:)=0.0_WP 
      
      DO KC=1,NCL3
         KM=KMV(KC)
         KP=KPV(KC)
         DO JC=1,N2DO(MYID)
            JM=JC-1
            JP=JC+1
            JJ=JCL2G(JC)
            QDX2=DYFI(JJ)*0.250_WP

            jjm=jmv(jj)           !@
            jjp=jpv(jj)           !@
            DO IC=1,NCL1_io
               IP=IPV_io(IC)
               IM=IMV_io(IC)                 
               H32=( ( Q_io(IC,JP,KC,2)+Q_io(IC,JP,KM,2) )*                &
                     ( Q_io(IC,JP,KC,3)/rm(jjp)+                        &
                       Q_io(IC,JC,KC,3)/rm(jj)   )-                     &       
                     ( Q_io(IC,JC,KC,2)+Q_io(IC,JC,KM,2) )*                &
                     ( Q_io(IC,JC,KC,3)/rm(jj)+                         &
                       Q_io(IC,JM,KC,3)/rm(jjm) ) )*QDX2   
                    
               H31=((Q_io(IP,JC,KC,1)+Q_io(IP,JC,KM,1))*                  &
                    (Q_io(IP,JC,KC,3)+Q_io(IC,JC,KC,3))-                  &
                    (Q_io(IC,JC,KC,1)+Q_io(IC,JC,KM,1))*                  &
                    (Q_io(IC,JC,KC,3)+Q_io(IM,JC,KC,3)) )*QDX1 
                    
               H33=((Q_io(IC,JC,KP,3)+Q_io(IC,JC,KC,3))*                   &
                    (Q_io(IC,JC,KP,3)+Q_io(IC,JC,KC,3))-                   &
                    (Q_io(IC,JC,KC,3)+Q_io(IC,JC,KM,3))*                   &
                    (Q_io(IC,JC,KC,3)+Q_io(IC,JC,KM,3)) )*QDX3/rm(jj)**2  

               RHSLLPHI_io(IC,JC,KC)=-(H31+H32+H33)
            ENDDO
         ENDDO
      ENDDO
            
      if (iswitch.eq.2) then 
          do kc=1,NCL3
             km=kmv(kc)
             kp=kpv(kc)
             DO jc=1,n2do(MYID)
                jm=jc-1
                jp=jc+1
                jj=JCL2G(JC)
                jjm=jmv(jj)           !@
                jjp=jpv(jj)           !@
                do ic=1,NCL1_io
                   if (jc.eq.1.and.myid.eq.0) then    
                      imsy = isym(kmv(kc))       
                      q2s1= (q_io(ic,jp,kc,2)-q_io(ic,jp,isym(kc),2))*0.50_WP/rc(jjp)
                      q2s2= (q_io(ic,jp,km,2)-q_io(ic,jp,imsy,    2))*0.50_WP/rc(jjp)
                      q2e =  q_io(ic,jp,kc,2)/rc(jjp) + q2s1
                      q2w =  q_io(ic,jp,km,2)/rc(jjp) + q2s2
                   else
                      q2e = q_io(ic,jp,kc,2)/rc(jjp)+q_io(ic,jc,kc,2)/rc(jj)
                      q2w = q_io(ic,jp,km,2)/rc(jjp)+q_io(ic,jc,km,2)/rc(jj)
                   end if
                   h32n= q_io(ic,jc,kc,3)*(q2e+q2w)*0.250_WP/rm(jj)
                   d11q1e=1.0_WP*(q2e-q2w)*DZI/rm(jj)
                   RHSLLPHI_io(ic,jc,kc)=-h32n+d11q1e/ren+RHSLLPHI_io(ic,jc,kc)              
                enddo
             enddo
          enddo 
     
      endif 
              
      RETURN
      END SUBROUTINE CONVECTION_Z_io
