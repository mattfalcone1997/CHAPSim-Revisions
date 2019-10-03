      SUBROUTINE COORDJR
      use mesh_info
      use init_info
      IMPLICIT NONE     

      INTEGER(4)  :: I, IMY
      INTEGER(4)  :: J, JJ, JC      
      INTEGER(4)  :: K
      INTEGER(4)  :: NFIL
      REAL(WP)      :: X2 
      REAL(WP)      :: STR2tmp
      REAL(WP),ALLOCATABLE      :: YNDtmp(:) 


      do jc=1,NND2
         rc(jc)=1.0_WP
         rm(jc)=1.0_WP
      enddo
      
      IF(ISTR2.EQ.1) THEN       

         IF(iswitch.EQ.1) THEN 
            YND(1)=-1.0_WP 
            STR2tmp = STR2*0.50_WP 
            DO J=1,NND2            
               X2=(J-1)/DBLE(NND2-1)  
               YND(J)=DTANH(STR2*(X2-0.50_WP))/DTANH(STR2tmp)
            ENDDO         
         ELSE  
            rc(1)=0.0_WP
            STR2tmp = STR2
            DO J=1,NND2
               X2=(J-1)/DBLE(NND2-1)  
               rc(j)=DTANH(STR2*(X2))/DTANH(STR2tmp)   
            ENDDO 
         END IF       

      ENDIF
    
      IF(iswitch.NE.1) THEN
         DO J=1,NCL2
            rm(J)=(rc(j)+rc(j+1))*0.50_WP      !@
         END DO
         rm(NCL2+1)=1.0_WP
      END IF

      DO J=1,NCL2
         YCC(J)=(YND(J)+YND(J+1))*0.50_WP
      ENDDO

      DO  I=1,NND1
          XND(I) = DBLE(I-1)*DX
      ENDDO
      
      IF(IOFLOWflg) THEN
          DO  I=1,NND1_io
              XND_io(I) = DBLE(I-1)*DX+DX*DBLE(NCL1)
          ENDDO
      END IF

      DO K=1,NND3
         ZND(K)=DBLE(K-1)/DBLE(NCL3)*ALX3
      ENDDO

      ALLOCATE ( YNDtmp(NND2) )
      YNDtmp = 0.0_WP 
                  
      IF(iswitch.EQ.1) THEN
         YNDtmp(1:NND2)=YND(1:NND2)
      ELSE
         YNDtmp(1:NND2)=RC(1:NND2)
      END IF 
           
      DO J=1,NCL2
         DYFI(J)=1.0_WP/ ( YNDtmp(J+1)-YNDtmp(J) ) 
      ENDDO   
      
      DYCI(1) = 1.0_WP/( YNDtmp(2) -YNDtmp(1)    )
      DYCI(NND2)= 1.0_WP/( YNDtmp(NND2)-YNDtmp(NND2-1) )       
      DO J=2,NND2-1
         DYCI(J)= 1.0_WP/ ( ( YNDtmp(J+1)-YNDtmp(J-1) )*0.50_WP )  
      END DO
      
      DEALLOCATE ( YNDtmp )
       
      JWGL1=1   
      DO IMY=0,NPSLV   
         DO JC=1,N2DO(MYID)
            JJ=JDSWT(IMY)-1+JC     
            IF(JJ.EQ.JWGL1) THEN  
               JWLC1=JC           
               IPWALL1=IMY       
            ENDIF
         ENDDO
      ENDDO
      IPWALL2=NPSLV   
      
      NFIL=15  
      
      OPEN(NFIL,FILE='COORDX1.dat')
      REWIND NFIL
      IF(IOFLOWflg) THEN
         WRITE(NFIL,*) NND1_io + NND1
      ELSE
         WRITE(NFIL,*) NND1
      END IF
      WRITE(NFIL,*) (XND(I),I=1,NND1)   
      IF(IOFLOWflg) WRITE(NFIL,*) (XND_io(I),I=1,NND1_io)    
      CLOSE(NFIL)

      NFIL=NFIL + 1
      OPEN(NFIL,FILE='COORDY2.dat')
      REWIND NFIL
      WRITE(NFIL,*) NND2  
      IF(iswitch .EQ. 1) THEN    
         WRITE(NFIL,*) (YND(J),J=1,NND2)   
      ELSE 
         WRITE(NFIL,*) (RC(J),J=1,NND2)
      END IF
      CLOSE(NFIL)
      
      NFIL=NFIL + 2
      OPEN(NFIL,FILE='COORDZ3.dat')
      REWIND NFIL
      WRITE(NFIL,*) NND3      
      WRITE(NFIL,*) (ZND(K),K=1,NND3)     
      CLOSE(NFIL)
      
      NFIL=NFIL + 3      
      OPEN (NFIL,FILE='METRIC_FIN.dat')
      IF(iswitch .EQ. 1) THEN
         WRITE(NFIL,'(A)') '# J, YND, YCC, DYC, DYF'     
         DO J=1,NCL2
            WRITE(NFIL,'(I10, 4ES15.7)') J,YND(J),YCC(J),1.0_WP/DYCI(J),1.0_WP/DYFI(J)  
         ENDDO
      ELSE
         WRITE(NFIL,'(A)') '# J, RC, RM, DYC, DYF'
         DO J=1,NCL2
            WRITE(NFIL,'(I10, 4ES15.7)') J,RC(J),RM(J),1.0_WP/DYCI(J),1.0_WP/DYFI(J)  
         ENDDO
      END IF    
      CLOSE(NFIL)
      
     
      RETURN
      END SUBROUTINE COORDJR
