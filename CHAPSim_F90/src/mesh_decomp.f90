       SUBROUTINE mesh_Ydecomp_master
       use mesh_info
       IMPLICIT NONE

       INTEGER(4)  :: NTOT1
       INTEGER(4)  :: NLOCAL
       INTEGER(4)  :: DEFICIT
       INTEGER(4)  :: IP

          ALLOCATE ( JDSWT(0:NPSLV) )
          ALLOCATE ( JDEWT(0:NPSLV) )
          ALLOCATE ( N2DO (0:NPSLV) )
          !initialising arrays
          JDSWT=0
          JDEWT=0

          NTOT1=NCL2  !Nodes in the y direction

          DO IP=0,NPSLV !Look over all ranks
             !Number of nodes in the y direction / number of ranks
             NLOCAL  = NTOT1 / NPTOT
             !^integer division therefore remainder is calculated
             DEFICIT = MOD(NTOT1,NPTOT)
             !Calculates array start point
             JDSWT(IP) = IP * NLOCAL + 1
             !Adds the deficit onto ranks below the value of the deficit
             JDSWT(IP) = JDSWT(IP) + MIN(IP,DEFICIT)
             !Recalculates NLOCAl
             IF (IP .LT. DEFICIT) THEN
                NLOCAL = NLOCAL + 1
             ENDIF
            !Calculation of end point for each process
             JDEWT(IP) = JDSWT(IP) + NLOCAL - 1
!
             IF ( (JDEWT(IP) .GT. NTOT1) .OR. ( IP .EQ. NPSLV ) ) THEN
                JDEWT(IP) = NTOT1
             END IF
          END DO

          call CHKHDL('   Ydecomp, rank, Start, End, Interval',MYID)
          DO IP=0,NPSLV
              N2DO(IP) = JDEWT(IP)-JDSWT(IP)+1 !Calculation of array size
               WRITE(*,'(A, 4I8.1)') '# MYID   0        Ydecomp,',ip, JDSWT(IP), JDEWT(IP),N2DO(IP)
          END DO

       end subroutine mesh_Ydecomp_master


       SUBROUTINE mesh_Zdecomp_master
       use mesh_info
       IMPLICIT NONE
       !=======================================================================
       !Similar reasoning to above but in the spanwise rather than wall-normal
       !direction
       !=======================================================================
       INTEGER(4)  :: NTOT1
       INTEGER(4)  :: NLOCAL
       INTEGER(4)  :: DEFICIT
       INTEGER(4)  :: IP

          ALLOCATE ( KDSWT(0:NPSLV) )
          ALLOCATE ( KDEWT(0:NPSLV) )
          ALLOCATE ( N3DO (0:NPSLV) )
          KDSWT=0
          KDEWT=0

          NTOT1=NCL3

          DO IP=0,NPSLV

             NLOCAL  = NTOT1 / NPTOT
             DEFICIT = MOD(NTOT1,NPTOT)

             KDSWT(IP) = IP * NLOCAL + 1
             KDSWT(IP) = KDSWT(IP) + MIN(IP,DEFICIT)
!
             IF (IP .LT. DEFICIT) THEN
                NLOCAL = NLOCAL + 1
             ENDIF
!
             KDEWT(IP) = KDSWT(IP) + NLOCAL - 1
!
             IF ( (KDEWT(IP) .GT. NTOT1) .OR. ( IP .EQ. NPSLV ) ) THEN
                KDEWT(IP) = NTOT1
             END IF
          END DO

          call CHKHDL('   Zdecomp, rank, Start, End, Interval',MYID)
          DO IP=0,NPSLV
              N3DO(IP) = KDEWT(IP)-KDSWT(IP)+1
              WRITE(*,'(A, 4I8.1)') '# MYID   0        Zdecomp,',ip, KDSWT(IP), KDEWT(IP),N3DO(IP)
          END DO


       end subroutine mesh_Zdecomp_master

!***********************************************************************
       SUBROUTINE mesh_decomp_bcast
       use mesh_info
       IMPLICIT NONE
       !========================================================================
       !Broadcasts mesh values and rank info calculated in the previous step to
       !all ranks
       !========================================================================
       IF(MYID.NE.0) THEN
           ALLOCATE ( JDSWT(0:NPSLV) )
           ALLOCATE ( JDEWT(0:NPSLV) )
           ALLOCATE ( KDSWT(0:NPSLV) )
           ALLOCATE ( KDEWT(0:NPSLV) )
           ALLOCATE ( N2DO (0:NPSLV) )
           ALLOCATE ( N3DO (0:NPSLV) )
           KDSWT=0
           KDEWT=0
           JDSWT=0
           JDEWT=0
           N2DO =0
           N3DO =0
       END IF

           CALL MPI_BCAST(JDSWT(0),NPTOT,MPI_INTEGER4,0,ICOMM,IERROR)
           CALL MPI_BCAST(JDEWT(0),NPTOT,MPI_INTEGER4,0,ICOMM,IERROR)

           CALL MPI_BCAST(KDSWT(0),NPTOT,MPI_INTEGER4,0,ICOMM,IERROR)
           CALL MPI_BCAST(KDEWT(0),NPTOT,MPI_INTEGER4,0,ICOMM,IERROR)

           CALL MPI_BCAST(N2DO(0),NPTOT,MPI_INTEGER4,0,ICOMM,IERROR)
           CALL MPI_BCAST(N3DO(0),NPTOT,MPI_INTEGER4,0,ICOMM,IERROR)


       END SUBROUTINE mesh_decomp_bcast
