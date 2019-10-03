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
          JDSWT=0
          JDEWT=0

          NTOT1=NCL2  

          DO IP=0,NPSLV
          
             NLOCAL  = NTOT1 / NPTOT  
             DEFICIT = MOD(NTOT1,NPTOT)
          
             JDSWT(IP) = IP * NLOCAL + 1
             JDSWT(IP) = JDSWT(IP) + MIN(IP,DEFICIT)
!
             IF (IP .LT. DEFICIT) THEN
                NLOCAL = NLOCAL + 1
             ENDIF
!
             JDEWT(IP) = JDSWT(IP) + NLOCAL - 1
!
             IF ( (JDEWT(IP) .GT. NTOT1) .OR. ( IP .EQ. NPSLV ) ) THEN
                JDEWT(IP) = NTOT1 
             END IF
          END DO

          call CHKHDL('   Ydecomp, rank, Start, End, Interval',MYID)
          DO IP=0,NPSLV
              N2DO(IP) = JDEWT(IP)-JDSWT(IP)+1
               WRITE(*,'(A, 4I8.1)') '# MYID   0        Ydecomp,',ip, JDSWT(IP), JDEWT(IP),N2DO(IP)
          END DO
       
       end subroutine mesh_Ydecomp_master
       
    
       SUBROUTINE mesh_Zdecomp_master
       use mesh_info
       IMPLICIT NONE      

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
       
