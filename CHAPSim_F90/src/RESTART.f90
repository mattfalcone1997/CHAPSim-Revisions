      SUBROUTINE RESTART_tg
      use init_info
      use mesh_info
      use flow_info
      use postprocess_info
      IMPLICIT NONE
      
      CHARACTER(7) :: PNTIM
      CHARACTER(30) :: FIL1
      CHARACTER(30) :: FILE3
      
      INTEGER(4)   :: ITIM
      INTEGER(4)   :: I, J, K, JJ
      INTEGER(4)   :: N1ML
      INTEGER(4)   :: N2DOL
      INTEGER(4)   :: N3ML
      
      REAL(WP),ALLOCATABLE  :: U_F0(:,:,:) 
      REAL(WP),ALLOCATABLE  :: V_F0(:,:,:)
      REAL(WP),ALLOCATABLE  :: W_F0(:,:,:)
      REAL(WP),ALLOCATABLE  :: P_F0(:,:,:)
      REAL(WP)  :: TIMEL
      REAL(WP)  :: RENL
      REAL(WP)  :: DTL

      ITIM=INT(TRST*MULTIM+0.3)
      WRITE(PNTIM,'(I7.7)') ITIM
      
      ALLOCATE( U_F0(1:NCL1,1:NCL2,  1:NCL3) )
      ALLOCATE( W_F0(1:NCL1,1:NCL2,  1:NCL3) )
      ALLOCATE( V_F0(1:NCL1,1:NCL2+1,1:NCL3) )
      ALLOCATE( P_F0(1:NCL1,1:NCL2,  1:NCL3) )
      

      FILE3='DNS2_periodicxz_RU'//PNTIM//'U.D'
      FIL1=FILE3
      if(myid==0) WRITE(*,'(A,A)') '#  Start reading data from                    ', FIL1
      OPEN(10,FILE=FIL1, FORM='UNFORMATTED')
      READ(10) N1ML,N2DOL,N3ML,ITERG0_TG
      READ(10) TIMEL,RENL,DTL
      if(myid==0) WRITE(*,'(A,A)') '#    After reading controlling parameters from ',  FIL1
      READ(10) U_F0
      CLOSE(10)
      if(myid==0) WRITE(*,'(A,A)') '#    End reading parameters from               ',  FIL1
      IF(DABS(TIMEL-TRST)>TSAVE1) CALL ERRHDL('Warning: Restarting from an ealier time step!',myid)

      FILE3='DNS2_periodicxz_RU'//PNTIM//'W.D'
      FIL1=FILE3
      if(myid==0) WRITE(*,'(A,A)') '#  Start reading data from                     ', FIL1
      OPEN(10,FILE=FIL1,FORM='UNFORMATTED')
      READ(10) N1ML,N2DOL,N3ML,ITERG0_TG
      READ(10) TIMEL,RENL,DTL
      if(myid==0) WRITE(*,'(A,A)') '#    After reading controlling parameters from ',  FIL1
      READ(10) W_F0
      CLOSE(10)
      if(myid==0) WRITE(*,'(A,A)') '#    End reading parameters from               ',  FIL1
      IF(DABS(TIMEL-TRST)>TSAVE1) CALL ERRHDL('Warning: Restarting from an ealier time step!',myid)

      FILE3='DNS2_periodicxz_RU'//PNTIM//'P.D'
      FIL1=FILE3
      if(myid==0) WRITE(*,'(A,A)') '#  Start reading data from                     ', FIL1
      OPEN(10,FILE=FIL1,FORM='UNFORMATTED')
      READ(10) N1ML,N2DOL,N3ML,ITERG0_TG
      READ(10) TIMEL,RENL,DTL
      if(myid==0) WRITE(*,'(A,A)') '#    After reading controlling parameters from ',  FIL1
      READ(10) P_F0
      CLOSE(10)
      if(myid==0) WRITE(*,'(A,A)') '#    End reading parameters from               ',  FIL1
      IF(DABS(TIMEL-TRST)>TSAVE1) CALL ERRHDL('Warning: Restarting from an ealier time step!',myid)

      FILE3='DNS2_periodicxz_RU'//PNTIM//'V.D'
      FIL1=FILE3
      if(myid==0) WRITE(*,'(A,A)') '#  Start reading data from ', FIL1
      OPEN(10,FILE=FIL1,FORM='UNFORMATTED')
      READ(10) N1ML,N2DOL,N3ML,ITERG0_TG
      READ(10) TIMEL,RENL,DTL
      if(myid==0) WRITE(*,'(A,A)') '#    After reading controlling parameters from  ',  FIL1
      READ(10) V_F0
      CLOSE(10)
      if(myid==0) WRITE(*,'(A,A)') '#    End reading controlling parameters from    ',  FIL1
      IF(DABS(TIMEL-TRST)>TSAVE1) CALL ERRHDL('Warning: Restarting from an ealier time step!',myid)
 
      DO J=1,N2DO(MYID)
         JJ=JCL2G(J)
         DO I=1,NCL1
            DO K=1,NCL3
               Q(I,J,K,1)=(U_F0(I,JJ,K))
               Q(I,J,K,2)=(V_F0(I,JJ,K))
               Q(I,J,K,3)=(W_F0(I,JJ,K))
               PR(I,J,K)=(P_F0(I,JJ,K))
            ENDDO
         ENDDO
      ENDDO

      IF (MYID.EQ.NPSLV) THEN
         J = N2DO(MYID)+1
         DO I=1,NCL1
            DO K=1,NCL3
               Q(I,J,K,2)=(V_F0(I,NND2,K))
            ENDDO
         ENDDO
      ENDIF

      phyTIME_TG=TIMEL
      
      DEALLOCATE (U_F0,W_F0,P_F0,V_F0)
      CALL MPI_BARRIER(ICOMM,IERROR)
      
      RETURN
      END SUBROUTINE RESTART_tg
      

      SUBROUTINE RESTART_io
      use init_info
      use mesh_info
      use flow_info
      use postprocess_info
      IMPLICIT NONE
      
      CHARACTER(7) :: PNTIM
      CHARACTER(30) :: FIL1
      CHARACTER(30) :: FILE3
      
      INTEGER(4)   :: ITIM
      INTEGER(4)   :: I, J, K, JJ
      INTEGER(4)   :: N1ML
      INTEGER(4)   :: N2DOL
      INTEGER(4)   :: N3ML
      
      REAL(WP),ALLOCATABLE  :: U_F0_io(:,:,:) 
      REAL(WP),ALLOCATABLE  :: V_F0_io(:,:,:)
      REAL(WP),ALLOCATABLE  :: W_F0_io(:,:,:)
      REAL(WP),ALLOCATABLE  :: P_F0_io(:,:,:)
      REAL(WP)  :: TIMEL
      REAL(WP)  :: RENL
      REAL(WP)  :: DTL


      ITIM=INT(TRST_io*MULTIM+0.3)
      WRITE(PNTIM,'(I7.7)') ITIM
      
      ALLOCATE( U_F0_io(0:NCL1_io+1,1:NCL2,1:NCL3) )
      ALLOCATE( W_F0_io(0:NCL1_io+1,1:NCL2,1:NCL3) )
      ALLOCATE( V_F0_io(0:NCL1_io+1,1:NND2,1:NCL3) )
      ALLOCATE( P_F0_io(0:NCL1_io+1,1:NCL2,1:NCL3) )
      

      FILE3='DNS2_inoudomain_RU'//PNTIM//'U.D'
      FIL1=FILE3
      if(myid==0) WRITE(*,'(A,A)') '#  Start reading data from                    ', FIL1
      OPEN(10,FILE=FIL1, FORM='UNFORMATTED')
      READ(10) N1ML,N2DOL,N3ML,ITERG0  
      READ(10) TIMEL,RENL,DTL
      if(myid==0) WRITE(*,'(A,A)') '#    After reading controlling parameters from ',  FIL1
      READ(10) U_F0_io
      CLOSE(10)
      if(myid==0) WRITE(*,'(A,A)') '#    End reading parameters from               ',  FIL1
      IF(DABS(TIMEL-TRST_io)>TSAVE1) CALL ERRHDL('Warning: Restarting from an ealier time step!',myid)

      FILE3='DNS2_inoudomain_RU'//PNTIM//'W.D'
      FIL1=FILE3
      if(myid==0) WRITE(*,'(A,A)') '#  Start reading data from                     ', FIL1
      OPEN(10,FILE=FIL1,FORM='UNFORMATTED')
      READ(10) N1ML,N2DOL,N3ML,ITERG0
      READ(10) TIMEL,RENL,DTL
      if(myid==0) WRITE(*,'(A,A)') '#    After reading controlling parameters from ',  FIL1
      READ(10) W_F0_io
      CLOSE(10)
      if(myid==0) WRITE(*,'(A,A)') '#    End reading parameters from               ',  FIL1
      IF(DABS(TIMEL-TRST_io)>TSAVE1) CALL ERRHDL('Warning: Restarting from an ealier time step!',myid)

      FILE3='DNS2_inoudomain_RU'//PNTIM//'P.D'
      FIL1=FILE3
      if(myid==0) WRITE(*,'(A,A)') '#  Start reading data from                     ', FIL1
      OPEN(10,FILE=FIL1,FORM='UNFORMATTED')
      READ(10) N1ML,N2DOL,N3ML,ITERG0
      READ(10) TIMEL,RENL,DTL
      if(myid==0) WRITE(*,'(A,A)') '#    After reading controlling parameters from ',  FIL1
      READ(10) P_F0_io
      CLOSE(10)
      if(myid==0) WRITE(*,'(A,A)') '#    End reading parameters from               ',  FIL1
      IF(DABS(TIMEL-TRST_io)>TSAVE1) CALL ERRHDL('Warning: Restarting from an ealier time step!',myid)

      FILE3='DNS2_inoudomain_RU'//PNTIM//'V.D'
      FIL1=FILE3
      if(myid==0) WRITE(*,'(A,A)') '#  Start reading data from ', FIL1
      OPEN(10,FILE=FIL1,FORM='UNFORMATTED')
      READ(10) N1ML,N2DOL,N3ML,ITERG0
      READ(10) TIMEL,RENL,DTL
      if(myid==0) WRITE(*,'(A,A)') '#    After reading controlling parameters from  ',  FIL1
      READ(10) V_F0_io
      CLOSE(10)
      if(myid==0) WRITE(*,'(A,A)') '#    End reading controlling parameters from    ',  FIL1
      IF(DABS(TIMEL-TRST_io)>TSAVE1) CALL ERRHDL('Warning: Restarting from an ealier time step!',myid)
      
      DO J=1,N2DO(MYID)
         JJ=JCL2G(J)
         DO I=0,NCL1_io+1,1
            DO K=1,NCL3
               Q_io(I,J,K,1)=(U_F0_io(I,JJ,K))
               Q_io(I,J,K,2)=(V_F0_io(I,JJ,K))
               Q_io(I,J,K,3)=(W_F0_io(I,JJ,K))
               PR_io(I,J,K)=(P_F0_io(I,JJ,K))
            ENDDO
         ENDDO
      ENDDO

      IF (MYID.EQ.NPSLV) THEN
         J=N2DO(MYID)+1
         DO I=0,NCL1_io+1,1
            DO K=1,NCL3
               Q_io(I,J,K,2)=(V_F0_io(I,NND2,K))
            ENDDO
         ENDDO
      ENDIF

      phyTIME=TIMEL

      
      DEALLOCATE (U_F0_io,W_F0_io,P_F0_io,V_F0_io)
      CALL MPI_BARRIER(ICOMM,IERROR)
      
      RETURN
      END SUBROUTINE RESTART_io

