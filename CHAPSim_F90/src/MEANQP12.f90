      SUBROUTINE MEANQP12
      use init_info
      use mesh_info
      use flow_info
      use postprocess_info
      IMPLICIT NONE     
      
      INTEGER(4)  :: I
      INTEGER(4)  :: J, JJ
      INTEGER(4)  :: K      
      INTEGER(4)  :: L
          
      if (iswitch.eq.2) then
         DO L=1,3    
            Qtmp = 0.0_WP
            DO I=1,NCL1
               DO J=1,N2DO(MYID)
                  JJ=JCL2G(J)
                  DO K=1,NCL3
                     if (l.eq.1) then
                        Qtmp(I,J,K) =Q(I,J,K,L)
                     elseif(l.eq.3) then
                        Qtmp(I,J,K) =Q(I,J,K,L)/rm(jj)
                     elseif(l.eq.2.and.myid.ne.0) then
                        Qtmp(I,J,K) =Q(I,J,K,L)/rc(jj)
                     elseif(l.eq.2.and.myid.eq.0.and.j.ne.1) then
                        Qtmp(I,J,K) =Q(I,J,K,L)/rc(jj)
                     elseif(l.eq.2.and.myid.eq.0.and.j.eq.1) then
                        Qtmp(i,1,k)=(q(i,2,k,2)-q(i,2,isym(k),2))*0.50_WP &
                                   /rc(2)
                     endif
                  ENDDO
               ENDDO
            ENDDO
            
            CALL VALME13(1,L)
            CALL STATN13(2,L)
            
         ENDDO
         
      else
      
         DO L=1,3
            Qtmp = 0.0_WP
            DO I=1,NCL1
               DO J=1,N2DO(MYID) 
                  JJ=JCL2G(J)
                  DO K=1,NCL3
                     Qtmp(I,J,K) =Q(I,J,K,L)
                  ENDDO
               ENDDO
            ENDDO

            CALL VALME13(1,L)
            CALL STATN13(2,L)
         ENDDO
         
      endif 
             
      RETURN
      END
