
       SUBROUTINE MEANVELOINI_tg(L)   
       use mesh_info
       use init_info
       use flow_info
       IMPLICIT NONE      
       
       INTEGER(4),INTENT(IN)  :: L
       
       INTEGER(4) :: I
       INTEGER(4) :: J, JJ
       INTEGER(4) :: K
       INTEGER(4) :: NII

       Qtmp = 0.0_WP
       IF(L.EQ.1) THEN
       
           Qtmp(:,:,:) = Q(:,:,:,L)
           
       ELSE IF(L.EQ.3) THEN
       
           DO J=1,N2DO(MYID)
              JJ=JCL2G(J)
              Qtmp(:,J,:) = Q(:,J,:,L)/RM(JJ)
           END DO  
       
       ELSE IF(L.EQ.2) THEN
           NII = 1
           IF(MYID.EQ.0 .and. ISWITCH.EQ.2 ) NII = 2
           
           DO J=NII,N2DO(MYID)
               JJ=JCL2G(J)
               Qtmp(:,J,:) = Q(:,J,:,L)/rc(jj)
           END DO
           
           IF(MYID.EQ.0 .and. ISWITCH.EQ.2 ) THEN
              J = 1
              do i=1,NCL1      
                 do k=1,NCL3
                    Qtmp(i,j,k)= 0.50_WP*( q(i,2,k,2)-q(i,2,isym(k),2) )/rc(2)
                 enddo
              end do
           END IF            
                 
       ELSE
         CALL ERRHDL('NO SUCH DIRCTION',myid) 
       END IF
       
       CALL VALME13ORL_tg(1,L)
       
       
       RETURN
       END SUBROUTINE MEANVELOINI_tg
       

       SUBROUTINE MEANVELOINI_io(L)   
       use mesh_info
       use init_info
       use flow_info
       IMPLICIT NONE      
       
       INTEGER(4),INTENT(IN)  :: L
       
       INTEGER(4) :: I
       INTEGER(4) :: J, JJ
       INTEGER(4) :: K
       INTEGER(4) :: NII

       Qtmp_io = 0.0_WP
       IF(L.EQ.1) THEN
       
           Qtmp_io(:,:,:) = Q_io(:,:,:,L)
           
       ELSE IF(L.EQ.3) THEN
       
           DO J=1,N2DO(MYID)
              JJ=JCL2G(J)
              Qtmp_io(:,J,:) = Q_io(:,J,:,L)/RM(JJ)
           END DO  
       
       ELSE IF(L.EQ.2) THEN
           NII = 1
           IF(MYID.EQ.0 .and. ISWITCH.EQ.2 ) NII = 2
           
           DO J=NII,N2DO(MYID)
               JJ=JCL2G(J)
               Qtmp_io(:,J,:) = Q_io(:,J,:,L)/rc(jj)
           END DO
           
           IF(MYID.EQ.0 .and. ISWITCH.EQ.2 ) THEN
              J = 1
              do i=1,NCL1      
                 do k=1,NCL3
                    Qtmp_io(i,j,k)= 0.50_WP*( q_io(i,2,k,2)-q_io(i,2,isym(k),2) )/rc(2)
                 enddo
              end do
           END IF            
                 
       ELSE
         CALL ERRHDL('NO SUCH DIRCTION',myid) 
       END IF
       
       CALL VALME13ORL_io(1,L)
       
       RETURN
       END SUBROUTINE MEANVELOINI_io

