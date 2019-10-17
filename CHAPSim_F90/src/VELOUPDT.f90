      SUBROUTINE VELOUPDT_tg(NS)

      use init_info
      use mesh_info
      use flow_info
      IMPLICIT NONE

      INTEGER(4),INTENT(IN) :: NS

      INTEGER(4) :: IC, IM
      INTEGER(4) :: JC, JM, JJ
      INTEGER(4) :: KC, KM
      INTEGER(4) :: NII
      REAL(WP)    :: DFX11, DFX22,DFX33
      REAL(WP)    :: COE1

      COE1 = DT*TALP(NS)

      DO KC=1,NCL3
         DO JC=1,N2DO(MYID)
            DO IC=1,NCL1
               IM=IMV(IC)
               DFX11=(DPH(IC,JC,KC)-DPH(IM,JC,KC))*DXI
               !Possibly where x and z velocities are modified
               Q(IC,JC,KC,1)=Q(IC,JC,KC,1)-DFX11*COE1
            ENDDO
         ENDDO
      ENDDO


      DO KC=1,NCL3
         KM=KMV(KC)
         DO JC=1,N2DO(MYID)
             DO IC=1,NCL1
                DFX33=(DPH(IC,JC,KC)-DPH(IC,JC,KM))*DZI
                Q(IC,JC,KC,3)=Q(IC,JC,KC,3)-DFX33*COE1
             ENDDO
         ENDDO
      ENDDO


      NII=1
      IF (MYID.EQ.0) THEN

         NII=2
         DO KC=1,NCL3
            DO IC=1,NCL1
               Q(IC,1,KC,2)=0.0_WP
            ENDDO
         ENDDO
      ENDIF

      IF (MYID.EQ.NPSLV) THEN
         DO KC=1,NCL3
            DO IC=1,NCL1
               Q(IC,N2DO(MYID)+1,KC,2)=0.0_WP
            ENDDO
         ENDDO
      ENDIF

      DO KC=1,NCL3
         DO JC=NII,N2DO(MYID)
            JJ=JCL2G(JC)
            DO IC=1,NCL1
               JM=JC-1
               DFX22=(DPH(IC,JC,KC)-DPH(IC,JM,KC))*DYCI(JJ)*rc(jj)  !
               Q(IC,JC,KC,2)=Q(IC,JC,KC,2)-DFX22*COE1
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE VELOUPDT_tg

      SUBROUTINE VELOUPDT_io(NS)

      use init_info
      use mesh_info
      use flow_info
      IMPLICIT NONE

      INTEGER(4),INTENT(IN) :: NS

      INTEGER(4) :: IC, IM
      INTEGER(4) :: JC, JM, JJ
      INTEGER(4) :: KC, KM
      INTEGER(4) :: NII
      REAL(WP)    :: DFX11, DFX22,DFX33
      REAL(WP)    :: COE1
      !A modification of the x velocity near the wall
      COE1 = DT*TALP(NS)
      DO KC=1,NCL3
         DO JC=1,N2DO(MYID)
            DO IC=1,NCL1_io
               IM=IMV_io(IC)
               DFX11=(DPH_io(IC,JC,KC)-DPH_io(IM,JC,KC))*DXI
               Q_io(IC,JC,KC,1)=Q_io(IC,JC,KC,1)-DFX11*COE1
            ENDDO
         ENDDO
      ENDDO

      DO KC=1,NCL3
         KM=KMV(KC)
         DO JC=1,N2DO(MYID)
             DO IC=1,NCL1_io
                DFX33=(DPH_io(IC,JC,KC)-DPH_io(IC,JC,KM))*DZI
                Q_io(IC,JC,KC,3)=Q_io(IC,JC,KC,3)-DFX33*COE1
             ENDDO
         ENDDO
      ENDDO

      NII=1
      IF (MYID.EQ.0) THEN
         NII=2
         DO KC=1,NCL3
            DO IC=1,NCL1_io
               Q_io(IC,1,KC,2)=0.0_WP
            ENDDO
         ENDDO
      ENDIF

      IF (MYID.EQ.NPSLV) THEN
         DO KC=1,NCL3
            DO IC=1,NCL1_io
               Q_io(IC,N2DO(MYID)+1,KC,2)=0.0_WP
            ENDDO
         ENDDO
      ENDIF

      DO KC=1,NCL3
         DO JC=NII,N2DO(MYID)
            JJ=JCL2G(JC)
            DO IC=1,NCL1_io
               JM=JC-1
               DFX22=(DPH_io(IC,JC,KC)-DPH_io(IC,JM,KC))*DYCI(JJ)*rc(jj)  !
               Q_io(IC,JC,KC,2)=Q_io(IC,JC,KC,2)-DFX22*COE1
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE VELOUPDT_io
