       SUBROUTINE SOLVE
       use init_info
       use flow_info
       IMPLICIT NONE

       REAL(WP)    :: ENDTIME(2)
       REAL(WP)    :: STARTTIME
       REAL(WP)    :: REN_TEMP

       INTEGER(4) :: NTSTF
       INTEGER(4) :: NS

       phyTIME    = 0.0_WP
       phyTIME_TG = 0.0_WP
       ITERG0 = 0
       ITERL  = 0

       CALL FLOWSTART

       REN_TEMP=REN
       CONVH0 = 0.0_WP
       IF(IOFLOWflg) CONVH0_io = 0.0_WP
       IF(MYID.EQ.0) THEN
          CALL WRTSCNINI_tg
          IF(IOFLOWflg) CALL WRTSCNINI_io
       END IF

       CALL MPI_BARRIER(ICOMM,IERROR)

       NTSTF = IDINT(TSTOP/DT)*1000
       DO ITERG=ITERG0+1,NTSTF
          STARTTIME=MPI_WTIME()

          ITERL = ITERL + 1
          phyTIME=phyTIME+DT
          phyTIME_TG = phyTIME_TG + DT
          IF(phyTIME .GT. TSTOP) CALL MPI_FINALIZE(IERROR)

          IF ( phyTIME_TG.LT.TLGRE ) THEN
              REN=REINI
              CVISC=1.0_WP/REN
          ELSE
              REN=REN_TEMP
              CVISC=1.0_WP/REN
          ENDIF

          IF (ITERG.GT.1) THEN
             CALL CFL_tg
             IF(IOFLOWflg) THEN
                CALL CFL_io
                CFLMM = DMAX1(CFLMM,CFLMM_io)
             END IF
             DT = DMIN1(DT,CFLGV/CFLMM)
          ENDIF

          !CALL BDFORCE

          CALL TIMECTL

          DO NS=1,NSST
             CALL SOLVERRK3_MOM_tg(NS)
          END DO
          ENDTIME(1)=MPI_WTIME()
          CPUTIME_tg=ENDTIME(1)-STARTTIME
          
          !Solved the inflow first then the spatially developing
          !one if it is enabled
          IF(IOFLOWflg) THEN
             CALL BC_TINLET
             DO NS=1,NSST
                CALL SOLVERRK3_MOM_io(NS)
             END DO
          ENDIF

          ENDTIME(2)=MPI_WTIME()
          CPUTIME_io=ENDTIME(2)-ENDTIME(1)
          CPUTIME   =ENDTIME(2)-STARTTIME

          CALL POSTPROCESS_tg
          IF(IOFLOWflg) CALL POSTPROCESS_io
!
       END DO


       RETURN

       END SUBROUTINE SOLVE
