      SUBROUTINE WRTVPROF_tg
      
      use mesh_info
      use init_info
      use flow_info
      IMPLICIT NONE
      
      INTEGER(4) :: FLLG=29
      INTEGER(4) :: J, JJ
      INTEGER(4) :: IP
      INTEGER(4) :: CID
      INTEGER(1) :: mm
      REAL(WP)     :: VELO(NCL2,3)
      REAL(WP)     :: Qttmp(3)

      IF(MYID.EQ.0) THEN
         VELO = 0.0_WP
         VELO(1:N2DO(MYID),1:3)=Q(1,1:N2DO(MYID),1,1:3)
         
         DO IP=1,NPSLV
            call mpi_ssend(mm, 1, MPI_INTEGER1, ip, ip, ICOMM,IERROR)
            DO JJ=JDSWT(1),NCL2
                IF(IP/=CLGIP(JJ)) CYCLE
                CALL MPI_RECV(CID,1,MPI_INTEGER4,IP,IP,ICOMM,STS,IERROR)
                CALL MPI_RECV(Qttmp(1),3,MPI_DOUBLE_PRECISION,IP,IP,ICOMM,STS,IERROR)
                VELO(CID,1:3)=Qttmp(1:3)
!                write(*,*) CID,VELO(CID,1:3)
            END DO 
         END DO
         
         OPEN(FLLG,FILE='Vini_POISEUILLE.dat')
          
         IF(iswitch.eq.1) THEN
             WRITE(FLLG,'(A)') '#  YCC, Uini, Uini+perb, Vini+perb, Wini+perb'
             DO JJ=1,NCL2
                WRITE(FLLG,'(I8,5ES15.7)') JJ,YCC(JJ),Vini(JJ), &
                          VELO(JJ,1:3)
             END DO
         ELSE if(iswitch.eq.2) THEN
             WRITE(FLLG,'(A)') '#  RM, Vini, Vini+perb, Uini, Wini'
             DO JJ=1,NCL2
                WRITE(FLLG,'(I8,5ES15.7)') JJ,RM(JJ),Vini(JJ), &
                          VELO(JJ,1:3)
             END DO
         ELSE
             CALL ERRHDL('NO SUCH GEO',myid)
         END IF
         
         CLOSE(FLLG)
      
      ELSE
         call mpi_recv(mm, 1, MPI_INTEGER1, 0, myid, ICOMM,STS,IERROR)
         DO J=1,N2DO(MYID)
            JJ=JCL2G(J)
            Qttmp(1:3)=Q(1,J,1,1:3)
            CALL MPI_SEND(JJ,1,MPI_INTEGER4,0,myid,ICOMM,IERROR)
            CALL MPI_SEND(Qttmp(1),3,MPI_DOUBLE_PRECISION,0,myid,ICOMM,IERROR)
         END DO
      
      END IF


      RETURN
    
      END  SUBROUTINE WRTVPROF_tg
      
      
      SUBROUTINE WRTVPROF_io
      
      use mesh_info
      use init_info
      use flow_info
      IMPLICIT NONE
      
      INTEGER(4) :: FLLG=29
      INTEGER(4) :: J, JJ
      INTEGER(4) :: IP
      INTEGER(4) :: CID
      INTEGER(1) :: mm
      REAL(WP)     :: VELO(NCL2,3)
      REAL(WP)     :: Qttmp(3)

!**********************************************************************
      IF(MYID.EQ.0) THEN
         VELO = 0.0_WP
         VELO(1:N2DO(MYID),1:3)=Q_io(1,1:N2DO(MYID),1,1:3)
         
         DO IP=1,NPSLV
            call mpi_ssend(mm, 1, MPI_INTEGER1, ip, ip, ICOMM,IERROR)
            DO JJ=JDSWT(1),NCL2
                IF(IP/=CLGIP(JJ)) CYCLE
                CALL MPI_RECV(CID,1,MPI_INTEGER4,IP,IP,ICOMM,STS,IERROR)
                CALL MPI_RECV(Qttmp(1),3,MPI_DOUBLE_PRECISION,IP,IP,ICOMM,STS,IERROR)
                VELO(CID,1:3)=Qttmp(1:3)
!                write(*,*) CID,VELO(CID,1:3)
            END DO 
         END DO
         
         OPEN(FLLG,FILE='Vini_POISEUILLE_io.dat')
          
         IF(iswitch.eq.1) THEN
             WRITE(FLLG,'(A)') '#  YCC, Uini, Uini+perb, Vini+perb, Wini+perb'
             DO JJ=1,NCL2
                WRITE(FLLG,'(I8,5ES15.7)') JJ,YCC(JJ),Vini(JJ), &
                          VELO(JJ,1:3)
             END DO
         ELSE if(iswitch.eq.2) THEN
             WRITE(FLLG,'(A)') '#  RM, Vini, Vini+perb, Uini, Wini'
             DO JJ=1,NCL2
                WRITE(FLLG,'(I8,5ES15.7)') JJ,RM(JJ),Vini(JJ), &
                          VELO(JJ,1:3)
             END DO
         ELSE
             CALL ERRHDL('NO SUCH GEO',myid)
         END IF
         
         CLOSE(FLLG)
      
      ELSE
         call mpi_recv(mm, 1, MPI_INTEGER1, 0, myid, ICOMM,STS,IERROR)
         DO J=1,N2DO(MYID)
            JJ=JCL2G(J)
            Qttmp(1:3)=Q_io(NCL1_io/2,J,NCL3/2,1:3)
            CALL MPI_SEND(JJ,1,MPI_INTEGER4,0,myid,ICOMM,IERROR)
            CALL MPI_SEND(Qttmp(1),3,MPI_DOUBLE_PRECISION,0,myid,ICOMM,IERROR)
         END DO
      
      END IF


      RETURN
    
      END  SUBROUTINE WRTVPROF_io

