
      SUBROUTINE INTFC_BC_QP_tg
        !====================================================
        !Updating information from MPI halo
        !====================================================
      use mesh_info
      use flow_info
      IMPLICIT NONE


       INTEGER(4)  :: I, IK
       INTEGER(4)  :: K
       INTEGER(4)  :: N
       INTEGER(4)  :: ITAG
       INTEGER(4)  :: IDESF
       INTEGER(4)  :: IDESB
       INTEGER(4)  :: TRC_STS(MPI_STATUS_SIZE)
       INTEGER(4)  :: NSZ

       BSEN_F = 0.0_WP
       BSEN_L = 0.0_WP
       BREC_F = 0.0_WP
       BREC_L = 0.0_WP
       NSZ =NCL1*NCL3*(NDV+1)

       DO I=1,NCL1
          DO K=1,NCL3
             IK=(I-1)*NCL3+K
             BSEN_F(IK,1)=PR(I,1,K)     !y=local Floor
             BSEN_L(IK,1)=PR(I,N2DO(MYID),K)  !y=local ROOF
          ENDDO
       ENDDO

       DO N=1,NDV
          DO I=1,NCL1
              DO K=1,NCL3
                 IK=(I-1)*NCL3+K
                 BSEN_F(IK,N+1)=Q(I,1,K,N)        !y=local Floor
                 BSEN_L(IK,N+1)=Q(I,N2DO(MYID),K,N)     !y=local ROOF
              ENDDO
          ENDDO
       ENDDO


       ITAG=0

       IF ( (MYID.GT.0) .AND. (MYID.LT.NPSLV) ) THEN
          IDESF = MYID+1    !next myid
          IDESB = MYID-1    !last myid


          CALL  MPI_SENDRECV( BSEN_F(1,1),NSZ,                &
                MPI_DOUBLE_PRECISION,IDESB,ITAG,BREC_F(1,1),NSZ,&
                MPI_DOUBLE_PRECISION,IDESF,ITAG,ICOMM,TRC_STS,IERROR)

          CALL  MPI_SENDRECV( BSEN_L(1,1),NSZ,                   &
                MPI_DOUBLE_PRECISION,IDESF,ITAG,BREC_L(1,1),NSZ, &
                MPI_DOUBLE_PRECISION,IDESB,ITAG,ICOMM,TRC_STS,IERROR)

          DO I=1,NCL1
             DO K=1,NCL3
                IK=(I-1)*NCL3+K
                PR(I,0,K)=BREC_L(IK,1)
                PR(I,N2DO(MYID)+1,K)=BREC_F(IK,1)
             ENDDO
          ENDDO


          DO N=1,NDV
             DO I=1,NCL1
                DO K=1,NCL3
                   IK=(I-1)*NCL3+K
                   Q(I,0,K,N)=BREC_L(IK,N+1)
                   Q(I,N2DO(MYID)+1,K,N)=BREC_F(IK,N+1)
                ENDDO
             ENDDO
         ENDDO

       ENDIF

       IF(MYID.EQ.NPSLV) THEN
          IDESF = 0
          IDESB = MYID-1

          CALL  MPI_SENDRECV( BSEN_F(1,1),NSZ,                &
                MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F(1,1),NSZ,&
                MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
          CALL  MPI_SENDRECV( BSEN_L(1,1),NSZ,               &
                MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L(1,1),NSZ, &
                MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)

          DO I=1,NCL1
             DO K=1,NCL3
                IK=(I-1)*NCL3+K
                PR(I,0,K)=BREC_L(IK,1)
                PR(I,N2DO(MYID)+1,K)=PR(I,N2DO(MYID),K)
             ENDDO
          ENDDO

          DO I=1,NCL1
             DO K=1,NCL3
                 IK=(I-1)*NCL3+K
                 DO N=1,3
                    Q(I,0,K,N)=BREC_L(IK,N+1)
                 ENDDO

                 DO N=1,3,2
                    Q(I,N2DO(MYID)+1,K,N)=-Q(I,N2DO(MYID),K,N)
                 ENDDO
                 Q(I,N2DO(MYID)+1,K,2)=0.0_WP
                ENDDO
             ENDDO

       ENDIF


       IF (MYID.EQ.0) THEN
          IDESF = MYID+1
          IDESB = NPSLV

          CALL  MPI_SENDRECV( BSEN_F(1,1),NSZ,                 &
                MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F(1,1),NSZ,  &
                MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
          CALL  MPI_SENDRECV( BSEN_L(1,1),NSZ,                 &
                MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L(1,1),NSZ,  &
                MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)

          DO I=1,NCL1
             DO K=1,NCL3
                IK=(I-1)*NCL3+K
                PR(I,0,K)=PR(I,1,K)
                PR(I,N2DO(MYID)+1,K)=BREC_F(IK,1)
             ENDDO
          ENDDO

          DO I=1,NCL1
             DO K=1,NCL3
                 IK=(I-1)*NCL3+K
                 DO N=1,3
                    Q(I,N2DO(MYID)+1,K,N)=BREC_F(IK,N+1)
                 ENDDO
                 DO N=1,3,2
                    Q(I,0,K,N)=-Q(I,1,K,N) !Where it enforces the lower wall BC
                 ENDDO
                 Q(I,1,K,2)=0.0_WP !Sets y velocity near the wall to zero
                 Q(I,0,K,2)=0.0_WP
               ENDDO
          ENDDO

       ENDIF


      RETURN
      END SUBROUTINE INTFC_BC_QP_tg

      SUBROUTINE INTFC_BC_QP_io
      use mesh_info
      use flow_info
      IMPLICIT NONE


       INTEGER(4)  :: I, IK
       INTEGER(4)  :: K
       INTEGER(4)  :: N
       INTEGER(4)  :: ITAG
       INTEGER(4)  :: IDESF
       INTEGER(4)  :: IDESB
       INTEGER(4)  :: TRC_STS(MPI_STATUS_SIZE)
       INTEGER(4)  :: NSZ

!********************************************************************
       BSEN_F_io = 0.0_WP
       BSEN_L_io = 0.0_WP
       BREC_F_io = 0.0_WP
       BREC_L_io = 0.0_WP
       NSZ =(NCL1_io+2)*NCL3*(NDV+1)

       DO I=0,NCL1_io+1,1
          DO K=1,NCL3
             IK=I*NCL3+K
             BSEN_F_io(IK,1)=PR_io(I,1,K)     !y=local Floor
             BSEN_L_io(IK,1)=PR_io(I,N2DO(MYID),K)  !y=local ROOF
          ENDDO
       ENDDO

       DO N=1,NDV
          DO I=0,NCL1_io+1,1
              DO K=1,NCL3
                 IK=I*NCL3+K
                 BSEN_F_io(IK,N+1)=Q_io(I,1,K,N)        !y=local Floor
                 BSEN_L_io(IK,N+1)=Q_io(I,N2DO(MYID),K,N)     !y=local ROOF
              ENDDO
          ENDDO
       ENDDO


       ITAG=0

       IF ( (MYID.GT.0) .AND. (MYID.LT.NPSLV) ) THEN
          IDESF = MYID+1    !next myid
          IDESB = MYID-1    !last myid


          CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                MPI_DOUBLE_PRECISION,IDESB,ITAG,BREC_F_io(1,1),NSZ,&
                MPI_DOUBLE_PRECISION,IDESF,ITAG,ICOMM,TRC_STS,IERROR)


          CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                   &
                MPI_DOUBLE_PRECISION,IDESF,ITAG,BREC_L_io(1,1),NSZ, &
                MPI_DOUBLE_PRECISION,IDESB,ITAG,ICOMM,TRC_STS,IERROR)


          DO I=0,NCL1_io+1,1
             DO K=1,NCL3
                IK=I*NCL3+K
                PR_io(I,0,K)=BREC_L_io(IK,1)
                PR_io(I,N2DO(MYID)+1,K)=BREC_F_io(IK,1)
             ENDDO
          ENDDO


          DO N=1,NDV
             DO I=0,NCL1_io+1,1
                DO K=1,NCL3
                   IK=I*NCL3+K
                   Q_io(I,0,K,N)=BREC_L_io(IK,N+1)
                   Q_io(I,N2DO(MYID)+1,K,N)=BREC_F_io(IK,N+1)
                ENDDO
             ENDDO
         ENDDO

       ENDIF

       IF(MYID.EQ.NPSLV) THEN
          IDESF = 0
          IDESB = MYID-1

          CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F_io(1,1),NSZ,&
                MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
          CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,               &
                MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L_io(1,1),NSZ, &
                MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)


          DO I=0,NCL1_io+1,1
             DO K=1,NCL3
                IK=I*NCL3+K
                PR_io(I,0,K)=BREC_L_io(IK,1)
                PR_io(I,N2DO(MYID)+1,K)=PR_io(I,N2DO(MYID),K)
             ENDDO
          ENDDO

          DO I=0,NCL1_io+1,1
             DO K=1,NCL3
                 IK=I*NCL3+K
                 DO N=1,3
                    Q_io(I,0,K,N)=BREC_L_io(IK,N+1)
                 ENDDO

                 DO N=1,3,2
                    Q_io(I,N2DO(MYID)+1,K,N)=-Q_io(I,N2DO(MYID),K,N)
                 ENDDO
                 Q_io(I,N2DO(MYID)+1,K,2)=0.0_WP
                
                ENDDO
             ENDDO

       ENDIF


       IF (MYID.EQ.0) THEN
          IDESF = MYID+1
          IDESB = NPSLV

          CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                 &
                MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F_io(1,1),NSZ,  &
                MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
          CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                 &
                MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L_io(1,1),NSZ,  &
                MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)

          DO I=0,NCL1_io+1,1
             DO K=1,NCL3
                IK=I*NCL3+K
                PR_io(I,0,K)=PR_io(I,1,K)
                PR_io(I,N2DO(MYID)+1,K)=BREC_F_io(IK,1)
             ENDDO
          ENDDO

          DO I=0,NCL1_io+1,1
             DO K=1,NCL3
                 IK=I*NCL3+K
                 DO N=1,3
                    Q_io(I,N2DO(MYID)+1,K,N)=BREC_F_io(IK,N+1)
                 ENDDO
                 DO N=1,3,2
                    Q_io(I,0,K,N)=-Q_io(I,1,K,N) !Enforcing the no slip through gradient through 0
                 ENDDO
                 Q_io(I,1,K,2)=0.0_WP
                 Q_io(I,0,K,2)=0.0_WP
               ENDDO
          ENDDO

       ENDIF


      RETURN
      END SUBROUTINE INTFC_BC_QP_io

      SUBROUTINE INTFC_BC_QP_INOUT
      use mesh_info
      use flow_info
      IMPLICIT NONE


       INTEGER(4)  :: I, IK
       INTEGER(4)  :: K
       INTEGER(4)  :: N
       INTEGER(4)  :: ITAG
       INTEGER(4)  :: IDESF
       INTEGER(4)  :: IDESB
       INTEGER(4)  :: TRC_STS(MPI_STATUS_SIZE)
       INTEGER(4)  :: NSZ

       BSEN_F_io = 0.0_WP
       BSEN_L_io = 0.0_WP
       BREC_F_io = 0.0_WP
       BREC_L_io = 0.0_WP

       NSZ =2*NCL3*(NDV+1)


       DO K=1,NCL3
          BSEN_F_io(K,1)=PR_io(0,1,K)           !y=local Floor
          BSEN_L_io(K,1)=PR_io(0,N2DO(MYID),K)  !y=local ROOF
       ENDDO

       DO K=1,NCL3
          IK = NCL3+K
          BSEN_F_io(IK,1)=PR_io(NCL1_io+1,1,K)           !y=local Floor
          BSEN_L_io(IK,1)=PR_io(NCL1_io+1,N2DO(MYID),K)  !y=local ROOF
       END DO

       DO N=1,NDV
          DO K=1,NCL3
              BSEN_F_io(K,N+1)=Q_io(0,1,K,N)              !y=local Floor
              BSEN_L_io(K,N+1)=Q_io(0,N2DO(MYID),K,N)     !y=local ROOF
          ENDDO
       ENDDO

       DO N=1,NDV
          DO K=1,NCL3
              IK = NCL3+K
              BSEN_F_io(IK,N+1)=Q_io(NCL1_io+1,1,K,N)           !y=local Floor
              BSEN_L_io(IK,N+1)=Q_io(NCL1_io+1,N2DO(MYID),K,N)  !y=local ROOF
          ENDDO
       ENDDO


       ITAG=0

       IF ( (MYID.GT.0) .AND. (MYID.LT.NPSLV) ) THEN
          IDESF = MYID+1    !next myid
          IDESB = MYID-1    !last myid


          CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                MPI_DOUBLE_PRECISION,IDESB,ITAG,BREC_F_io(1,1),NSZ,&
                MPI_DOUBLE_PRECISION,IDESF,ITAG,ICOMM,TRC_STS,IERROR)


          CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                   &
                MPI_DOUBLE_PRECISION,IDESF,ITAG,BREC_L_io(1,1),NSZ, &
                MPI_DOUBLE_PRECISION,IDESB,ITAG,ICOMM,TRC_STS,IERROR)


          DO K=1,NCL3
             PR_io(0,0,K)           =BREC_L_io(K,1)
             PR_io(0,N2DO(MYID)+1,K)=BREC_F_io(K,1)
          ENDDO
          DO K=1,NCL3
             IK = NCL3+K
             PR_io(NCL1_io+1,0,K)           =BREC_L_io(IK,1)
             PR_io(NCL1_io+1,N2DO(MYID)+1,K)=BREC_F_io(IK,1)
          ENDDO

          DO N=1,NDV
             DO K=1,NCL3
                Q_io(0,0,K,N)            = BREC_L_io(K,N+1)
                Q_io(0,N2DO(MYID)+1,K,N) = BREC_F_io(K,N+1)
             ENDDO
          ENDDO

          DO N=1,NDV
             DO K=1,NCL3
                IK = NCL3+K
                Q_io(NCL1_io+1,0,K,N)            = BREC_L_io(IK,N+1)
                Q_io(NCL1_io+1,N2DO(MYID)+1,K,N) = BREC_F_io(IK,N+1)
             ENDDO
          ENDDO


       ENDIF

       IF(MYID.EQ.NPSLV) THEN
          IDESF = 0
          IDESB = MYID-1

          CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F_io(1,1),NSZ,&
                MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
          CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,               &
                MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L_io(1,1),NSZ, &
                MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)


          DO K=1,NCL3
             PR_io(0,0,K)           =BREC_L_io(K,1)
             PR_io(0,N2DO(MYID)+1,K)=PR_io(0,N2DO(MYID),K)
          ENDDO
          DO K=1,NCL3
             IK = NCL3+K
             PR_io(NCL1_io+1,0,K)           =BREC_L_io(IK,1)
             PR_io(NCL1_io+1,N2DO(MYID)+1,K)=PR_io(NCL1_io+1,N2DO(MYID),K)
          ENDDO

          DO N=1,NDV
             DO K=1,NCL3
                Q_io(0,0,K,N)  = BREC_L_io(K,N+1)
             ENDDO
          ENDDO
          DO K=1,NCL3
             Q_io(0,N2DO(MYID)+1,K,1) = -1.0_WP*Q_io(0,N2DO(MYID),K,1)
             Q_io(0,N2DO(MYID)+1,K,3) = -1.0_WP*Q_io(0,N2DO(MYID),K,3)
             Q_io(0,N2DO(MYID)+1,K,2) =  0.0_WP
          END DO

          DO N=1,NDV
             DO K=1,NCL3
                IK = NCL3+K
                Q_io(NCL1_io+1,0,K,N) = BREC_L_io(IK,N+1)
             ENDDO
          ENDDO

          DO K=1,NCL3
             Q_io(NCL1_io+1,N2DO(MYID)+1,K,1) = -1.0_WP*Q_io(NCL1_io+1,N2DO(MYID),K,1)
             Q_io(NCL1_io+1,N2DO(MYID)+1,K,3) = -1.0_WP*Q_io(NCL1_io+1,N2DO(MYID),K,3)
             Q_io(NCL1_io+1,N2DO(MYID)+1,K,2) =  0.0_WP
          END DO


       ENDIF

       IF (MYID.EQ.0) THEN
          IDESF = MYID+1
          IDESB = NPSLV

          CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                 &
                MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F_io(1,1),NSZ,  &
                MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
          CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                 &
                MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L_io(1,1),NSZ,  &
                MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)

          DO K=1,NCL3
             PR_io(0,0,K)           =PR_io(0,1,K)
             PR_io(0,N2DO(MYID)+1,K)=BREC_F_io(K,1)
          ENDDO
          DO K=1,NCL3
             IK = NCL3+K
             PR_io(NCL1_io+1,0,K)           =PR_io(NCL1_io+1,1,K)
             PR_io(NCL1_io+1,N2DO(MYID)+1,K)=BREC_F_io(IK,1)
          ENDDO

          DO N=1,NDV
             DO K=1,NCL3
                Q_io(0,N2DO(MYID)+1,K,N)  = BREC_F_io(K,N+1)
             ENDDO
          ENDDO
          DO K=1,NCL3
             Q_io(0,0,K,1) = -1.0_WP*Q_io(0,1,K,1)
             Q_io(0,0,K,3) = -1.0_WP*Q_io(0,1,K,3)
             Q_io(0,0,K,2) =  0.0_WP
             Q_io(0,1,K,2) =  0.0_WP
          END DO

          DO N=1,NDV
             DO K=1,NCL3
                IK = NCL3+K
                Q_io(NCL1_io+1,N2DO(MYID)+1,K,N) = BREC_F_io(IK,N+1)
             ENDDO
          ENDDO

          DO K=1,NCL3
             Q_io(NCL1_io+1,0,K,1) = -1.0_WP*Q_io(NCL1_io+1,1,K,1)
             Q_io(NCL1_io+1,0,K,3) = -1.0_WP*Q_io(NCL1_io+1,1,K,3)
             Q_io(NCL1_io+1,0,K,2) =  0.0_WP
             Q_io(NCL1_io+1,1,K,2) =  0.0_WP
          END DO


       ENDIF

      RETURN
      END SUBROUTINE INTFC_BC_QP_INOUT
