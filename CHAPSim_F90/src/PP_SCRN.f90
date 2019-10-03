      SUBROUTINE PP_SCRN
      use init_info
      use mesh_info
      use flow_info
      use postprocess_info
      IMPLICIT NONE     
      
      INTEGER(4)  :: J, JJ
      INTEGER(4)  :: L
      REAL(WP)     :: UMX1,UMX1_WORK
      REAL(WP)     :: UME1,UME1_WORK
      REAL(WP)     :: UR1(NCL2)
      REAL(WP)     :: VR1(0:NCL2)
      REAL(WP)     :: WR1(NCL2)
     
      L=1
      ENE11=0.0_WP
      UME1=0.0_WP
      UMX1=0.0_WP
      URMSX1=0.0_WP
      DO J=1,N2DO(MYID)  !@
         JJ=JCL2G(J)
         UME1=UME1+rm(jj)/DYFI(JJ)*STA13(1,L,JJ)/(rc(NND2)**2)
         UMX1=DMAX1(UMX1,DABS(STA13(1,L,JJ))) 
      ENDDO

      CALL MPI_BARRIER(ICOMM,IERROR)
      CALL MPI_ALLREDUCE(UMX1, UMX1_WORK, 1, MPI_DOUBLE_PRECISION, &
                   MPI_MAX, ICOMM, IERROR)
      CALL MPI_ALLREDUCE(UME1, UME1_WORK, 1, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, ICOMM, IERROR)
      UMX1U(L)=UMX1_WORK
      if (iswitch.eq.1) then
          UME1U(L)=UME1_WORK/2.0_WP
      else
          UME1U(L)=UME1_WORK*2.0_WP
      endif      
         
      DO J=1,NCL2
         UR1(J) = DSQRT(STA13(2,1,J))
         WR1(J) = DSQRT(STA13(2,3,J))
         ENE11=ENE11+1.0_WP/DYFI(J)*(UR1(J)**2+WR1(J)**2)
         URMSX1=DMAX1(URMSX1,UR1(J))
      END DO
   
      DO J=1,NCL2
         VR1(J) = DSQRT(STA13(2,2,J))
         ENE11=ENE11+rm(j)/DYCI(J)*VR1(J)**2
      END DO
   
      if (iswitch.eq.1) then
         URMSC=(UR1(NCL2/2)+UR1(NCL2/2+1))/2.0_WP
         VRMSC=VR1(NCL2/2)
         WRMSC=(WR1(NCL2/2)+WR1(NCL2/2+1))/2.0_WP
         ENE11=ENE11/2.0_WP  
         ENE11=ENE11/2.0_WP
      else
         URMSC=UR1(1)
         VRMSC=VR1(1)
         WRMSC=WR1(1)
         ENE11=ENE11/2.0_WP  
         ENE11=ENE11/2.0_WP
      endif
      
      RETURN
      END
