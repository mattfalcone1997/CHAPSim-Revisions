      SUBROUTINE DIVGCK
      use mesh_info
      use flow_info
      IMPLICIT NONE
                     
      REAL(WP) :: DIVX
      REAL(WP) :: QMA
      REAL(WP) :: QMAX_WORK
      REAL(WP) :: DQCAP
      INTEGER(4) :: K, KP
      INTEGER(4) :: J, JP, JJ
      INTEGER(4) :: I, IP
          
       DIVX      =0.0_WP
       QMA       =0.0_WP
       QMAX_WORK = 0.0_WP
       DQCAP     = 0.0_WP
       DO K=1,NCL3
          KP=KPV(K)
          DO J=1,N2DO(MYID)
             JP=J+1
             JJ=JCL2G(J)
             DO I=1,NCL1
                IP=IPV(I)
                DQCAP=(Q(IP,J,K,1)-Q(I,J,K,1))*DXI*rm(jj)**2       &
                     +(Q(I,JP,K,2)-Q(I,J,K,2))*DYFI(JJ)*rm(jj)     &
                     +(Q(I,J,KP,3)-Q(I,J,K,3))*DZI
                DIVX=DMAX1(DABS(DQCAP),DIVX)
             ENDDO
          ENDDO
       ENDDO
      
       QMA=DIVX  

       CALL MPI_BARRIER(ICOMM,IERROR)
       CALL MPI_ALLREDUCE(QMA, QMAX_WORK, 1, MPI_DOUBLE_PRECISION,  &
             MPI_MAX, ICOMM, IERROR)
     
       MAXDIVGV=QMAX_WORK
      
      RETURN
      END

      SUBROUTINE DIVGCK_io
      use mesh_info
      use flow_info
      use init_info 
      IMPLICIT NONE
                     
      REAL(WP) :: DIVX1, DIVX2, DIVX3
      REAL(WP) :: QMA1, QMA2, QMA3
      REAL(WP) :: QMAX_WORK1, QMAX_WORK2, QMAX_WORK3
      REAL(WP) :: DQCAP
      INTEGER(4) :: K, KP
      INTEGER(4) :: J, JP, JJ
      INTEGER(4) :: I, IP
          
       DIVX1      = 0.0_WP
       QMA1       = 0.0_WP
       QMAX_WORK1 = 0.0_WP
       DIVX2      = 0.0_WP
       QMA2       = 0.0_WP
       QMAX_WORK2 = 0.0_WP
       DIVX3      = 0.0_WP
       QMA3       = 0.0_WP
       QMAX_WORK3 = 0.0_WP
       DQCAP      = 0.0_WP
       MAXDIVGV_io= 0.0_WP
       
       DO K=1,NCL3
          KP=KPV(K)
          DO J=1,N2DO(MYID)
             JP=J+1
             JJ=JCL2G(J)
             DO I=1,NCL1_io-1 !test
                IP=IPV_io(I)
                DQCAP=(Q_io(IP,J,K,1)-Q_io(I,J,K,1))*DXI*rm(jj)**2       &
                     +(Q_io(I,JP,K,2)-Q_io(I,J,K,2))*DYFI(JJ)*rm(jj)     &
                     +(Q_io(I,J,KP,3)-Q_io(I,J,K,3))*DZI
                DIVX1=DMAX1(DABS(DQCAP),DIVX1)

             ENDDO
             
                I = 0
                IP=IPV_io(I)
                DQCAP=(Q_io(IP,J,K,1)-Q_io(I,J,K,1))*DXI*rm(jj)**2       &
                     +(Q_io(I,JP,K,2)-Q_io(I,J,K,2))*DYFI(JJ)*rm(jj)     &
                     +(Q_io(I,J,KP,3)-Q_io(I,J,K,3))*DZI
                DIVX2=DMAX1(DABS(DQCAP),DIVX2)

                I = NCL1_io
                IP=IPV_io(I)
                DQCAP=(Q_io(IP,J,K,1)-Q_io(I,J,K,1))*DXI*rm(jj)**2       &
                     +(Q_io(I,JP,K,2)-Q_io(I,J,K,2))*DYFI(JJ)*rm(jj)     &
                     +(Q_io(I,J,KP,3)-Q_io(I,J,K,3))*DZI
                DIVX3=DMAX1(DABS(DQCAP),DIVX3)
             
          ENDDO
       ENDDO
      
       QMA1=DIVX1
       QMA2=DIVX2  
       QMA3=DIVX3

       CALL MPI_BARRIER(ICOMM,IERROR)
       CALL MPI_ALLREDUCE(QMA1, QMAX_WORK1, 1, MPI_DOUBLE_PRECISION,  &
             MPI_MAX, ICOMM, IERROR)
       CALL MPI_ALLREDUCE(QMA2, QMAX_WORK2, 1, MPI_DOUBLE_PRECISION,  &
             MPI_MAX, ICOMM, IERROR)
       CALL MPI_ALLREDUCE(QMA3, QMAX_WORK3, 1, MPI_DOUBLE_PRECISION,  &
             MPI_MAX, ICOMM, IERROR)
     
       MAXDIVGV_io(1)=QMAX_WORK1
       MAXDIVGV_io(2)=QMAX_WORK2
       MAXDIVGV_io(3)=QMAX_WORK3
      
      RETURN
      END
