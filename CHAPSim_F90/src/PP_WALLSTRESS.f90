      SUBROUTINE PP_WALLSTRESS
      use init_info
      use mesh_info
      use flow_info
      use postprocess_info
      
      IMPLICIT NONE
      INTEGER(4)  :: I, K, J, JJ
      REAL(WP)     :: DQL1, DQU1
      REAL(WP)     :: DQLW, DQUW
      REAL(WP)     :: CFL11
      REAL(WP)     :: CFU11
  
      if (iswitch.eq.1) then
      
         CFLW1=0.0_WP
         CFUW1=0.0_WP
         CFAV =0.0_WP
         
         J=JWLC1
         JJ=JCL2G(J)
         IF (JJ.EQ.JWGL1) THEN
         
             DQL1=0.0_WP
             DO K=1,NCL3
                DO I=1,NCL1
                   DQLW=(Q(I,J,K,1)-0.0_WP)/(YCC(JJ)-YND(JJ)) 
                   DQL1= DQL1 + DQLW
                END DO
             END DO
             CFL11=DQL1*VL1313
             CFLW1=CFL11/REN
             
         ENDIF

         J=N2DO(MYID)
         JJ=JCL2G(J)
         IF (JJ.EQ.NCL2) THEN
             DQU1=0.0_WP
             DO K=1,NCL3
                DO I=1,NCL1
                   DQUW=(0.0_WP-Q(I,J,K,1))/(YCC(JJ)-YND(JJ+1))
                   DQU1= DQU1 + DQUW
                END DO
             END DO
             CFU11=DQU1*VL1313
             CFUW1=CFU11/REN
         ENDIF
         
         CALL MPI_BCAST(CFLW1,1, MPI_DOUBLE_PRECISION,IPWALL1,ICOMM, &
              IERROR)  
         CALL MPI_BCAST(CFUW1,1, MPI_DOUBLE_PRECISION,IPWALL2,ICOMM, &
              IERROR)
         CFAV=( DABS(CFUW1) + DABS(CFLW1) ) / 2.0_WP
         
      elseif(iswitch.eq.2) then 

         cfav=0.0_WP

         j=n2do(MYID)
         jj=JCL2G(J)
         if (jj.eq.NCL2) then
            DQU1=0.0_WP
            DO K=1,NCL3
               DO I=1,NCL1
                  DQUW=-(Q(I,J,K,1))/(rm(jj)-rc(jj+1))
                  DQU1= DQU1 + DQUW
               END DO
            END DO
            CFU11=DQU1*VL1313
            CFUW1=CFU11/REN
         ENDIF
         CALL MPI_BCAST(CFUW1,1, MPI_DOUBLE_PRECISION,IPWALL2,ICOMM, &
              IERROR)      
         CFAV=DABS(CFUW1)
      endif
    
      RETURN
      END
      

      MODULE WALLSTRESS_info
      USE WPRECISION
      
      REAL(WP),ALLOCATABLE     :: Cf_LW_io(:)
      REAL(WP),ALLOCATABLE     :: Cf_UW_io(:)
      REAL(WP),ALLOCATABLE     :: Cf_LW_tg(:)
      REAL(WP),ALLOCATABLE     :: Cf_UW_tg(:)
      INTEGER(4)               :: NCOUNT1 = 0
      INTEGER(4)               :: TECFLG 
      END MODULE WALLSTRESS_info
      SUBROUTINE WALLSTRESS_WRT
      USE WALLSTRESS_info
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      
      INTEGER(4) :: I
      
      
            TECFLG = 122
            
            IF( (NCOUNT1.EQ.1) .AND. (NREAD.NE.2) .AND. (NREAD_io.NE.2) ) THEN 
                OPEN(TECFLG,FILE='RESULT.CF.plt')
                WRITE(TECFLG,'(A)') 'TITLE = "Cf on the wall"'
                WRITE(TECFLG,'(A)') 'VARIABLES = "X", "Cf_LW", "Cf_UW" '
                WRITE(TECFLG,'(A,1I6.1,1ES13.5,A)') 'ZONE T=" ',ITERG,phyTIME, '  "'
            ELSE
                OPEN(TECFLG,FILE='RESULT.CF.plt', position='append')
                WRITE(TECFLG,'(A,1I6.1,1ES13.5,A)') 'ZONE T=" ',ITERG,phyTIME, '  "'
            END IF

            DO I=1,NCL1
               WRITE(TECFLG,'(3ES15.7)') XND(I),Cf_LW_tg(I),Cf_UW_tg(I)
            END DO
            DO I=1,NCL1_io
               WRITE(TECFLG,'(3ES15.7)') XND_io(I),Cf_LW_io(I),Cf_UW_io(I)
            END DO

      RETURN
      END SUBROUTINE WALLSTRESS_WRT
!*******************************************************************
      SUBROUTINE PP_WALLSTRESS_io
      USE WALLSTRESS_info
      use init_info
      use mesh_info
      use flow_info
      use postprocess_info
      
      IMPLICIT NONE
      INTEGER(4)  :: I, K, J, JJ
      REAL(WP)     :: DQL1, DQU1
      REAL(WP)     :: DQLW, DQUW
      REAL(WP)     :: CFL11
      REAL(WP)     :: CFU11
!***********************************************************

      NCOUNT1 = NCOUNT1 + 1
      
      IF(NCOUNT1 == 1) THEN
        ALLOCATE ( Cf_LW_io(NCL1_io) )
        ALLOCATE ( Cf_UW_io(NCL1_io) )
        ALLOCATE ( Cf_LW_tg(NCL1) )
        ALLOCATE ( Cf_UW_tg(NCL1) )
      END IF
      
      if (iswitch.eq.1) then
!-------------------------------------------------------
!        LOWER WALLS SHEAR
         J=JWLC1
         JJ=JCL2G(J)
         IF (JJ.EQ.JWGL1) THEN
             Cf_LW_tg(:)=0.0_WP
             DO I=1,NCL1
                DO K=1,NCL3
                   DQLW=(Q(I,J,K,1)-0.0_WP)/(YCC(JJ)-YND(JJ)) 
                   Cf_LW_tg(I)= Cf_LW_tg(I) + DQLW
                END DO
             END DO
             Cf_LW_tg(:) =2.0_WP * DABS( Cf_LW_tg(:) )/NCL3/REN
         
             Cf_LW_io(:)=0.0_WP
             DO I=1,NCL1_io
                DO K=1,NCL3
                   DQLW=(Q_io(I,J,K,1)-0.0_WP)/(YCC(JJ)-YND(JJ)) 
                   Cf_LW_io(I)= Cf_LW_io(I) + DQLW
                END DO
             END DO
             Cf_LW_io(:) =2.0_WP * DABS( Cf_LW_io(:) )/NCL3/REN
         ENDIF

         J=N2DO(MYID)
         JJ=JCL2G(J)
         IF (JJ.EQ.NCL2) THEN
             Cf_UW_tg = 0.0_WP
             DO I=1,NCL1
                DO K=1,NCL3
                   DQUW=(0.0_WP-Q(I,J,K,1))/(YCC(JJ)-YND(JJ+1))
                   Cf_UW_tg(I)= Cf_UW_tg(I) + DQUW
                END DO
             END DO
             Cf_UW_tg(:) =2.0_WP * DABS( Cf_UW_tg(:) )/NCL3/REN
             
             Cf_UW_io = 0.0_WP
             DO I=1,NCL1_io
                DO K=1,NCL3
                   DQUW=(0.0_WP-Q_io(I,J,K,1))/(YCC(JJ)-YND(JJ+1))
                   Cf_UW_io(I)= Cf_UW_io(I) + DQUW
                END DO
             END DO
             Cf_UW_io(:) =2.0_WP * DABS( Cf_UW_io(:) )/NCL3/REN
             
         ENDIF
         
         CALL MPI_BCAST(Cf_LW_tg,NCL1, MPI_DOUBLE_PRECISION,IPWALL1,ICOMM, &
              IERROR)  
         CALL MPI_BCAST(Cf_UW_tg,NCL1, MPI_DOUBLE_PRECISION,IPWALL2,ICOMM, &
              IERROR)
              
         CALL MPI_BCAST(Cf_LW_io,NCL1_io, MPI_DOUBLE_PRECISION,IPWALL1,ICOMM, &
              IERROR)  
         CALL MPI_BCAST(Cf_UW_io,NCL1_io, MPI_DOUBLE_PRECISION,IPWALL2,ICOMM, &
              IERROR)
         
      elseif(iswitch.eq.2) then 

         cfav=0.0_WP

         j=n2do(MYID)
         jj=JCL2G(J)
         if (jj.eq.NCL2) then
            DQU1=0.0_WP
            DO K=1,NCL3
               DO I=1,NCL1
                  DQUW=-(Q(I,J,K,1))/(rm(jj)-rc(jj+1))
                  DQU1= DQU1 + DQUW
               END DO
            END DO
            CFU11=DQU1*VL1313
            CFUW1=CFU11/REN
         ENDIF
         CALL MPI_BCAST(CFUW1,1, MPI_DOUBLE_PRECISION,IPWALL2,ICOMM, &
              IERROR)      
         CFAV=DABS(CFUW1)
      endif
      
      IF(MYID == 0) CALL WALLSTRESS_WRT
      
      RETURN
      END

