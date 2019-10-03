      SUBROUTINE INTFC_BC_DPH_tg
      use mesh_info
      use flow_info
      IMPLICIT NONE
       
       INTEGER(4)  :: I, IK
       INTEGER(4)  :: K
       INTEGER(4)  :: ITAG
       INTEGER(4)  :: IDESF
       INTEGER(4)  :: IDESB
       INTEGER(4)  :: TRC_STS(MPI_STATUS_SIZE)
       INTEGER(4)  :: NSZ

       BSEN_F = 0.0_WP
       BSEN_L = 0.0_WP
       BREC_F = 0.0_WP
       BREC_L = 0.0_WP
       
       NSZ=NCL1*NCL3    
       DO I=1,NCL1
          DO K=1,NCL3
             IK=(I-1)*NCL3+K
             BSEN_F(IK,1)=DPH(I,1,K)           !y=local Floor     
             BSEN_L(IK,1)=DPH(I,N2DO(MYID),K)  !y=local ROOF                                           
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
                DPH(I,0,K)=BREC_L(IK,1)
                DPH(I,N2DO(MYID)+1,K)=BREC_F(IK,1)
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
                DPH(I,0,K)=BREC_L(IK,1)
                DPH(I,N2DO(MYID)+1,K)=DPH(I,N2DO(MYID),K)      
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
                DPH(I,0,K)=DPH(I,1,K)
                DPH(I,N2DO(MYID)+1,K)=BREC_F(IK,1)
             ENDDO
          ENDDO
          
       ENDIF
      
      
      RETURN
      END SUBROUTINE INTFC_BC_DPH_tg
      

      SUBROUTINE INTFC_BC_DPH_io
      use mesh_info
      use flow_info
      IMPLICIT NONE
       
       INTEGER(4)  :: I, IK
       INTEGER(4)  :: K, J
       INTEGER(4)  :: ITAG
       INTEGER(4)  :: IDESF
       INTEGER(4)  :: IDESB
       INTEGER(4)  :: TRC_STS(MPI_STATUS_SIZE)
       INTEGER(4)  :: NSZ

       BSEN_F_io = 0.0_WP
       BSEN_L_io = 0.0_WP
       BREC_F_io = 0.0_WP
       BREC_L_io = 0.0_WP
       
       NSZ=(NCL1_io+2)*NCL3    
       DO I=0,NCL1_io+1,1
          DO K=1,NCL3
             IK=(I)*NCL3+K
             BSEN_F_io(IK,1)=DPH_io(I,1,K)   
             BSEN_L_io(IK,1)=DPH_io(I,N2DO(MYID),K)                                      
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
                IK=(I)*NCL3+K
                DPH_io(I,0,K)=BREC_L_io(IK,1)
                DPH_io(I,N2DO(MYID)+1,K)=BREC_F_io(IK,1)
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
                IK=(I)*NCL3+K
                DPH_io(I,0,K)=BREC_L_io(IK,1)
                DPH_io(I,N2DO(MYID)+1,K)=DPH_io(I,N2DO(MYID),K)      
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
                IK=(I)*NCL3+K
                DPH_io(I,0,K)=DPH_io(I,1,K)
                DPH_io(I,N2DO(MYID)+1,K)=BREC_F_io(IK,1)
             ENDDO
          ENDDO
          
       ENDIF
      
      
      RETURN
      END SUBROUTINE INTFC_BC_DPH_io

