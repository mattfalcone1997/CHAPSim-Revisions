!********************************MODULE**********************************
      MODULE FFT99_info
      USE WPRECISION
      
         REAL(WP),ALLOCATABLE  :: AK1(:)
         REAL(WP),ALLOCATABLE  :: AK3(:)
         
         REAL(WP),ALLOCATABLE  :: TRIGXX1(:)
         REAL(WP),ALLOCATABLE  :: TRIGXX2(:)
         REAL(WP),ALLOCATABLE  :: TRIGXX3(:)
        
         REAL(WP),ALLOCATABLE  :: TRIGXC1(:)
         REAL(WP),ALLOCATABLE  :: TRIGXC2(:)
         REAL(WP),ALLOCATABLE  :: TRIGXC3(:)
      
         INTEGER(4)   :: IFXX1(13)
         INTEGER(4)   :: IFXX2(13)
         INTEGER(4)   :: IFXX3(13)
         INTEGER(4)   :: IFXC3(13)
         INTEGER(4)   :: IFXC1(13)
         INTEGER(4)   :: IFXC2(13)
         
         REAL(WP),ALLOCATABLE     :: XR  (:,:)
         REAL(WP),ALLOCATABLE     :: WORK(:,:)
         COMPLEX(8),ALLOCATABLE  :: XA  (:,:)
         COMPLEX(8),ALLOCATABLE  :: WOR (:,:)
       
         
         INTEGER(4)   :: LHFP
         
         INTEGER(4)   :: L,M,N, ML, NL, MLmax, NLmax
         
         REAL(WP)       :: DXQI0, DZQI0
         
         REAL(WP),ALLOCATABLE  :: AMJP(:,:,:)
         REAL(WP),ALLOCATABLE  :: ACJP(:,:,:)
         REAL(WP),ALLOCATABLE  :: APJP(:,:,:)
         
         REAL(WP),ALLOCATABLE :: RHSLLPHIRe(:,:,:)
         REAL(WP),ALLOCATABLE :: RHSLLPHIIm(:,:,:)
         
         REAL(WP),ALLOCATABLE :: FJ(:,:)
         REAL(WP),ALLOCATABLE :: BCJ(:,:)
      
      END MODULE FFT99_info
!***********************************************************************
!***********************************************************************
      SUBROUTINE FFT99_POIS3D_INIT
      USE FFT99_info
      USE mesh_info
      use mpi_info
      IMPLICIT NONE
      
      L = NCL1
      M = NCL2
      N = NCL3
      
      ML = N2DO(myid)
      NL = N3DO(myid)
      
      MLmax = N2DO(0)
      NLmax = N3DO(0)
      
      LHFP=L/2+1
      
      DXQI0 = DXQI
      DZQI0 = DZQI
      
      ALLOCATE ( AK1(L) )
      ALLOCATE ( AK3(N) )
         
      ALLOCATE ( TRIGXX1(3*L/2+1) )
      ALLOCATE ( TRIGXX2(3*M/2+1) )
      ALLOCATE ( TRIGXX3(3*N/2+1) )
        
      ALLOCATE ( TRIGXC1(2*L) )
      ALLOCATE ( TRIGXC2(2*M) )
      ALLOCATE ( TRIGXC3(2*N) )
      
      AK1 = 0.0_WP
      AK3 = 0.0_WP
      TRIGXX1 = 0.0_WP
      TRIGXX2 = 0.0_WP
      TRIGXX3 = 0.0_WP
      TRIGXC1 = 0.0_WP
      TRIGXC2 = 0.0_WP
      TRIGXC3 = 0.0_WP
      
      ALLOCATE (  XR  (L+2,N)  )
      ALLOCATE (  WORK(L+1,N)  )
      ALLOCATE (  XA  (N,L)  )
      ALLOCATE (  WOR (N,L)  )
      
      ALLOCATE ( AMJP(LHFP,M,NLmax ) )
      ALLOCATE ( ACJP(LHFP,M,NLmax ) )
      ALLOCATE ( APJP(LHFP,M,NLmax ) )
      AMJP = 0.0_WP
      ACJP = 0.0_WP
      APJP = 0.0_WP
      
      ALLOCATE ( RHSLLPHIRe(LHFP,MLmax,N) )
      ALLOCATE ( RHSLLPHIIm(LHFP,MLmax,N) )
      RHSLLPHIRe = 0.0_WP
      RHSLLPHIIm = 0.0_WP
      
      
      ALLOCATE ( FJ (LHFP,M) )
      ALLOCATE ( BCJ(LHFP,2)  )

      FJ  = 0.0_WP
      BCJ = 0.0_WP

      CALL FFT99_ROOT
      CALL TDMA_PHI_COEF
      
      RETURN      
      END SUBROUTINE FFT99_POIS3D_INIT
      
!***********************************************************************
!***********************************************************************
      SUBROUTINE FFT99_ROOT
      USE FFT99_info
      IMPLICIT NONE
      
      INTEGER(4)  :: I, K
      INTEGER(4)  :: N1MH, N3MH 
      INTEGER(4)  :: N1MP, N3MP
      
      REAL(WP),ALLOCATABLE :: AP(:)
      REAL(WP),ALLOCATABLE :: AN(:)
      REAL(WP)  :: PI
      
!>    @note: set up coefficients IFXX1, TRIGXX1 in x and z directions
      CALL FFTFAX(L,IFXX1,TRIGXX1)  
      CALL FFTFAX(N,IFXX3,TRIGXX3)
!>    @note: set up coefficients IFXC1, TRIGXC1 in x and z directions
      CALL CFTFAX(N,IFXC3,TRIGXC3)
      CALL CFTFAX(L,IFXC1,TRIGXC1)
      
 
      N1MH=L/2
      N3MH=N/2
      N1MP=N1MH+1
      N3MP=N3MH+1
      
      PI=2.0_WP*(DASIN(1.0_WP))

!>    @note calculate 2-2cos(2pi/nk) for z direction.
      ALLOCATE( AN(N) )
      AN = 0.0_WP     
      DO K=1,N3MH
         AN(K)=(K-1)*2.0_WP*PI
      ENDDO    
      DO K=N3MP,N
         AN(K)=-(N-K+1)*2.0_WP*PI
      ENDDO
      
      DO K=1,N
         AK3(K)=2.0_WP*(1.0_WP-DCOS(AN(K)/N))*DZQI0
      END DO 
      DEALLOCATE(AN)
           
!>    @note calculate 2-2cos(2pi/nk) for x direction.
      ALLOCATE( AP(L) ) 
      AP = 0.0_WP      
      DO I=1,N1MH
         AP(I)=(I-1)*2.0_WP*PI
      ENDDO      
      DO I=N1MP,L
         AP(I)=-(L-I+1)*2.0_WP*PI
      ENDDO    
  
      DO I=1,L
         AK1(I)=2.0_WP*(1.0_WP-DCOS(AP(I)/L))*DXQI0
      END DO   
      DEALLOCATE(AP)
     
      RETURN
      END SUBROUTINE FFT99_ROOT
!***********************************************************************
      SUBROUTINE TDMA_PHI_COEF
      use mesh_info
      use FFT99_info
      IMPLICIT NONE
      
      INTEGER(4)  :: I, J, K, KK
      
      ! For TDMA of Phi
      DO K=1,NL
         KK = KCL2G(K)
!>       @note calcuate FJ(I,J), ACJP(I,J), APJP(I,J), AMJP(I,J).         
         DO I=1,LHFP 
            DO J=1,M
               !ACC = 1.0_WP/(ACPH(J)+(-AK1(I)*(rm(j)**2)-AK3(KK))) 
               APJP(I,J,K)=APPH(J)!*ACC
               AMJP(I,J,K)=AMPH(J)!*ACC
               ACJP(I,J,K)= ACPH(J)+(-AK1(I)*(rm(j)**2)-AK3(KK))               
            ENDDO
         ENDDO
         
         IF (KK.EQ.1) THEN
             !write(*,*) '111,',AMJP(1,1,K),ACJP(1,1,K), APJP(1,1,K)
             J=1
             I=1
             ACJP(1,1,K)=ACPH(J)+(-AK1(I)*(rm(j)**2)-AK3(KK))
             APJP(1,1,K)=0.0_WP
             AMJP(1,1,K)=0.0_WP
             !write(*,*) '111,',AMJP(1,1,K),ACJP(1,1,K), APJP(1,1,K)
         ENDIF
      END DO
      
      RETURN
      
      END SUBROUTINE TDMA_PHI_COEF
      
!***********************************************************************
!***********************************************************************
!*********************************************************************** 
      SUBROUTINE FFT99_POIS3D_periodicxz
      USE mpi_info
      USE FFT99_info
      use mesh_info
      use flow_info
      implicit none
      
      integer(4) :: I, J, K, KK
      REAL(WP)    :: UN3M
      
      RHSLLPHIRe=0.0_WP
      RHSLLPHIIm=0.0_WP

!     FFT for x direction, and then CFFT for z direction.      
      UN3M=1.0_WP/DBLE(N)
      DO J=1,ML   
                       
          DO K=1,N
             XR(1,  K)  =RHSLLPHI(L,J,K)    !1
             XR(L+2,K)  =RHSLLPHI(1,J,K)   !NCL1+2
          ENDDO
          DO I=1,L
             DO K=1,N
                XR (I+1,K)=RHSLLPHI(I,J,K) !2 to NCL1+1
                WORK(I,K )=0.0_WP 
             ENDDO
          ENDDO
                          
          CALL FFT99(XR,WORK,TRIGXX1,IFXX1,1,L+2,L,N,-1)
             
          DO K=1,N
             DO I=1,LHFP
                XA(K,I)    =DCMPLX(XR(2*I-1,K),XR(2*I,K))
                XR(2*I-1,K)=0.0_WP
                XR(2*I,  K)=0.0_WP
             ENDDO
          ENDDO

          CALL CFFT99(XA,WOR,TRIGXC3,IFXC3,1,N,N,LHFP,-1)
            
          DO I=1,LHFP
             DO K=1,N
                RHSLLPHIRe(I,J,K)=DREAL(XA(K,I)*UN3M)
                RHSLLPHIIm(I,J,K)=DIMAG(XA(K,I)*UN3M)
             END DO
          END DO
            
      END DO
!     Correction for b.c.
      IF (MYID.EQ.0) THEN
          RHSLLPHIRe(1,1,1)=0.0_WP
          RHSLLPHIIm(1,1,1)=0.0_WP
      ENDIF 

!     TDMA for the real part of FFTxz
      FJ  = 0.0_WP
      BCJ = 0.0_WP
      F   = 0.0_WP     
      CALL TRASP23L2G_PHIRe
      DO K=1,NL
         KK=KCL2G(K)
         DO I=1,LHFP 
            DO J=1,M
               FJ(I,J)=F(I,J,K)         
           ENDDO
           BCJ(I,:) = 0.0_WP
         ENDDO
         
         IF (KK.EQ.1) THEN
             FJ(1,1)=0.0_WP
         ENDIF
         
         CALL TDMAIJJ_nonCYC(AMJP(:,:,K),ACJP(:,:,K),APJP(:,:,K), &
                             FJ,BCJ,1,M,1,LHFP)
         DO I=1,LHFP
            DO J=1,M
               F(I,J,K)=FJ(I,J)
            ENDDO
         ENDDO
         
      ENDDO
      CALL TRASP23G2L_PHIRe
      
!     TDMA for the imaginary part of FFTxz
      FJ  = 0.0_WP
      BCJ = 0.0_WP
      F   = 0.0_WP      
      CALL TRASP23L2G_PHIIm
      DO K=1,NL
         KK=KCL2G(K)        
         DO I=1,LHFP 
            DO J=1,M
               FJ(I,J)=F(I,J,K)           
            ENDDO
            BCJ(I,:) = 0.0_WP
         ENDDO
         
         IF (KK.EQ.1) THEN
             FJ(1,1)=0.0_WP
         ENDIF

         CALL TDMAIJJ_nonCYC(AMJP(:,:,K),ACJP(:,:,K),APJP(:,:,K), &
               FJ,BCJ,1,M,1,LHFP)

         DO I=1,LHFP
            DO J=1,M
               F(I,J,K)=FJ(I,J)
            ENDDO
         ENDDO
         
      ENDDO
      CALL TRASP23G2L_PHIIm
      
!     backward CFFT for z direction, and then backward FFT for x direction.  
      DPH = 0.0_WP 
      DO J=1,ML
           
          DO K=1,N
             DO I=1,LHFP
                XA(K,I)=DCMPLX(RHSLLPHIRe(I,J,K),RHSLLPHIIm(I,J,K))
             END DO
          END DO  

          CALL CFFT99(XA,WOR,TRIGXC3,IFXC3,1,N,N,LHFP,+1)
            
          DO  I=1,LHFP
              DO  K=1,N
                  XR(2*I-1,K)=DREAL(XA(K,I))
                  XR(2*I,  K)=DIMAG(XA(K,I))
              END DO
          END DO 
  
          CALL FFT99(XR,WORK,TRIGXX1,IFXX1,1,L+2,L,N,+1)
          
          DO I=1,L
             DO K=1,N
                DPH (I,J,K)=XR(I+1,K)
                WORK(I,K  )=0.0_WP 
             END DO
          END DO
      
      END DO
      
      RETURN
      END SUBROUTINE FFT99_POIS3D_periodicxz
!*********************************************************************** 

