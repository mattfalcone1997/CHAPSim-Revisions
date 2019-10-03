      SUBROUTINE VMAV_tg
      
      use mesh_info
      use flow_info
      IMPLICIT NONE
      
      INTEGER(4)  :: I, J, K, JJ
      REAL(WP)     :: VMA
      REAL(WP)     :: VMAX_WORK     
      
!>     @note Max. Q(i,j,k,1)       
       VMA = 0.0_WP
       DO K=1,NCL3
          DO J=1,N2DO(MYID)
             JJ=JCL2G(J) 
             DO I=1,NCL1              
                VMA=DMAX1(VMA,DABS(Q(I,J,K,1))) 
             ENDDO
          ENDDO
       ENDDO
     
       CALL MPI_BARRIER(ICOMM,IERROR)
       CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, ICOMM, IERROR)
       VMV(1)=VMAX_WORK
       
            
       VMA = 0.0_WP
       DO K=1,NCL3
          DO J=1,N2DO(MYID)
             JJ=JCL2G(J) 
             DO I=1,NCL1              
                VMA=DMAX1(VMA,DABS(Q(I,J,K,3)/rm(jj))) 
             ENDDO
          ENDDO
       ENDDO
     
       CALL MPI_BARRIER(ICOMM,IERROR)
       CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, ICOMM, IERROR)
       VMV(3)=VMAX_WORK
            
       VMA = 0.0_WP
       DO K=1,NCL3
          DO J=1,N2DO(MYID)
             JJ=JCL2G(J)
             IF(JJ==1) THEN
                DO I=1,NCL1              
                   VMA=DMAX1(VMA,DABS(Q(I,J,K,2)/rc(jj+1))) 
                ENDDO
             ELSE 
                DO I=1,NCL1              
                   VMA=DMAX1(VMA,DABS(Q(I,J,K,2)/rc(jj))) 
                ENDDO
             END IF
          ENDDO
       ENDDO
     
       CALL MPI_BARRIER(ICOMM,IERROR)
       CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, ICOMM, IERROR)
            
       VMV(2)=VMAX_WORK       
       
      RETURN
      
      END 
      

      SUBROUTINE VMAV_io
      use mesh_info
      use flow_info
      IMPLICIT NONE
      
      INTEGER(4)  :: I, J, K, JJ
      REAL(WP)     :: VMA
      REAL(WP)     :: VMAX_WORK     
      
    
       VMA = 0.0_WP
       DO K=1,NCL3
          DO J=1,N2DO(MYID)
             JJ=JCL2G(J) 
             DO I=0,NCL1_io+1,1             
                VMA=DMAX1(VMA,DABS(Q_io(I,J,K,1))) 
             ENDDO
          ENDDO
       ENDDO
     
       CALL MPI_BARRIER(ICOMM,IERROR)
       CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, ICOMM, IERROR)
       VMV_io(1)=VMAX_WORK
       
            
       VMA = 0.0_WP
       DO K=1,NCL3
          DO J=1,N2DO(MYID)
             JJ=JCL2G(J) 
             DO I=0,NCL1_io+1,1              
                VMA=DMAX1(VMA,DABS(Q_io(I,J,K,3)/rm(jj))) 
             ENDDO
          ENDDO
       ENDDO
     
       CALL MPI_BARRIER(ICOMM,IERROR)
       CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, ICOMM, IERROR)
       VMV_io(3)=VMAX_WORK
             
       VMA = 0.0_WP
       DO K=1,NCL3
          DO J=1,N2DO(MYID)
             JJ=JCL2G(J)
             IF(JJ==1) THEN
                DO I=0,NCL1_io+1,1              
                   VMA=DMAX1(VMA,DABS(Q_io(I,J,K,2)/rc(jj+1))) 
                ENDDO
             ELSE 
                DO I=0,NCL1_io+1,1              
                   VMA=DMAX1(VMA,DABS(Q_io(I,J,K,2)/rc(jj))) 
                ENDDO
             END IF
          ENDDO
       ENDDO
     
       CALL MPI_BARRIER(ICOMM,IERROR)
       CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, ICOMM, IERROR)
            
       VMV_io(2)=VMAX_WORK       
       
      RETURN
      
      END 

