      SUBROUTINE TRASP23L2G_RHS_io   !I23=1
      use flow_info
      use mesh_info
      IMPLICIT NONE
      
      REAL(WP)  :: SENBLK( (NCL1_io), N2DO(0), N3DO(0) )
      REAL(WP)  :: RECBLK( (NCL1_io), N2DO(0), N3DO(0) )
      INTEGER(4)  :: TRAS_STS(MPI_STATUS_SIZE)
      INTEGER(4)  :: NTOL 
      INTEGER(4)  :: JL, JG
      INTEGER(4)  :: IG
      INTEGER(4)  :: KL, KG
      INTEGER(4)  :: IP      
      INTEGER(4)  :: ISEND
      INTEGER(4)  :: IRECV
      INTEGER(4)  :: ITAG
      INTEGER(4)  :: NNTR
      
      NNTR=NCL1_io
      
      DO IP=0,NPSLV
          IF(MYID.EQ.IP) THEN
                   
            DO JL=1,N2DO(MYID)
                JG=JCL2G(JL)   
                DO KL=1,N3DO(MYID)
                   KG=KDSWT(MYID)-1+KL
                   DO IG=1,NNTR
                      F_io(IG,JG,KL) = RHS_io(IG,JL,KG)
                   END DO
                END DO 
            END DO 
               
          ELSE
          
            NTOL = NNTR * N2DO(0) * N3DO(0)
            DO KL=1,N3DO(IP)
               KG=KDSWT(IP)-1+KL
               DO IG=1,NNTR
                  DO JL=1,N2DO(MYID)
                     SENBLK(IG,JL,KL) = RHS_io(IG,JL,KG)
                  END DO
               END DO
            END DO
            
            ISEND = IP
            IRECV = IP
            ITAG=0
            CALL MPI_SENDRECV(SENBLK(1,1,1),NTOL,MPI_DOUBLE_PRECISION, &
                 IRECV,ITAG,RECBLK(1,1,1),NTOL,MPI_DOUBLE_PRECISION,   &
                 ISEND,ITAG,ICOMM,TRAS_STS,IERROR)
                 
            DO JL=1,N2DO(IP)
               JG=JDSWT(IP)+JL-1
               DO IG=1,NNTR
                  DO KL=1,N3DO(MYID)
                     F_io(IG,JG,KL) = RECBLK(IG,JL,KL)
                  END DO
               END DO 
            END DO 
                 
          END IF
       END DO
       
       RETURN
       END SUBROUTINE TRASP23L2G_RHS_io
       
!>*********************************************************************       
      SUBROUTINE TRASP23G2L_RHS_io  !I23=-1
     
      use mesh_info
      use flow_info
      IMPLICIT NONE
      
      REAL(WP)  :: SENBLK((NCL1_io),N2DO(0),N3DO(0))
      REAL(WP)  :: RECBLK((NCL1_io),N2DO(0),N3DO(0))
      INTEGER(4)  :: TRAS_STS(MPI_STATUS_SIZE)
      INTEGER(4)  :: NTOL 
      INTEGER(4)  :: JL, JG
      INTEGER(4)  :: IG
      INTEGER(4)  :: KL, KG
      INTEGER(4)  :: IP      
      INTEGER(4)  :: ISEND
      INTEGER(4)  :: IRECV
      INTEGER(4)  :: ITAG
      INTEGER(4)  :: NNTR
      

          
      NNTR=NCL1_io
      DO IP=0,NPSLV
          IF(MYID.EQ.IP) THEN
                   
            DO KL=1,N3DO(MYID)
               KG=KDSWT(MYID)-1+KL 
                DO JL=1,N2DO(MYID)
                   JG=JDSWT(MYID)-1+JL
                   DO IG=1,NNTR
                      RHS_io(IG,JL,KG) = F_io(IG,JG,KL)
                   END DO
                END DO 
            END DO 
               
          ELSE
            NTOL = NNTR * N3DO(0) * N2DO(0)
            DO JL=1,N2DO(IP)
               JG=JDSWT(IP)-1+JL
               DO IG=1,NNTR
                  DO KL=1,N3DO(MYID)
                     SENBLK(IG,JL,KL) = F_io(IG,JG,KL)
                  END DO
               END DO
            END DO
            ISEND = IP
            IRECV = IP
            ITAG=0
            CALL MPI_SENDRECV(SENBLK(1,1,1),NTOL,MPI_DOUBLE_PRECISION, &
                 IRECV,ITAG,RECBLK(1,1,1),NTOL,MPI_DOUBLE_PRECISION,   &
                 ISEND,ITAG,ICOMM,TRAS_STS,IERROR)
                 
            
            DO KL=1,N3DO(IP)
               KG=KDSWT(IP)+KL-1
               DO IG=1,NNTR
                  DO JL=1,N2DO(MYID)
                     RHS_io(IG,JL,KG) = RECBLK(IG,JL,KL)
                  END DO
               END DO 
            END DO 
                 
          END IF
       END DO
       
       RETURN
       END SUBROUTINE TRASP23G2L_RHS_io
