      SUBROUTINE  WRTSCNINI_tg
      use init_info
      use mesh_info
      use wrt_info
      IMPLICIT NONE
      
      WRITE(*,'(A)') '#**********Mesh and flow details based on initial flow field******'
      WRITE(*,'(A,I6,2(3X,A,F9.5))') '#   NCL1_TG=',NCL1,    'HX=', ALX1,                     'DX=   ', DX
      IF(IOFLOWflg) &
      WRITE(*,'(A,I6,2(3X,A,F9.5))') '#   NCL1_IO=',NCL1_io, 'HX=', ALX1*REAL(NCL1_io/NCL1), 'DX=   ', DX
      WRITE(*,'(A,I6,2(3X,A,F9.5))') '#   NCL3=   ',NCL3,    'HZ=', ALX3,                     'DZ=   ', DZ
      WRITE(*,'(A,I6,3(3X,A,F9.5))') '#   NCL2=   ',NCL2,    'HY=', ALX2,                     'DYMIN=', &
                      1.0_WP/DYCI(1), 'DYMAX=', 1.0_WP/DYCI(NCL2/2)
      WRITE(*,'(A,I12)'            ) '#   MESH_TG=',NCL1*NCL2*NCL3
      IF(IOFLOWflg) &
      WRITE(*,'(A,I12)'            ) '#   MESH_IO=',NCL1_io*NCL2*NCL3
      WRITE(*,'(A)') '#******************************************************************'
      
      IF(IOFLOWflg)  THEN
         CALL date_and_time(DATE=date,TIME=time)
         fllog=date(1:4)//'.'//date(5:8)//'.'//time(1:4)//'.log'
         OPEN(logflg,FILE='history.periodicxz.'//fllog)
      ELSE
         logflg = 6
      END IF
      
      
      IF (iswitch.eq.1) then
          WRITE(logflg,1000) 
 1000     FORMAT('# ',2X,4HSTEP,5X,4HPhyT,5X,5HCpuTT,2X,6HDivMax,&
           5X,3HCFL,3X,5HUmean,4X,4HUmax,5X,3HReM,4X,3HReX,&
           4X,6HEnergy,5X,5HU'max,6X,3HU'c,8X,3HV'c,8X,3HW'c,&
           8X,6HRE_tau,7X,2HCf,3X)
      ELSE IF(iswitch.eq.2) then
          WRITE(logflg,1004) 
 1004    FORMAT(1X,'# ',1X,4HSTEP,7X,1HT,8X,6HC_TIME,3X,6HDIVMAX,    &
         5X,3HCFL,4X,5HUMEAN,3X,4HUMAX,4X,3HREM,4X,3HREX,            &
         5X,6HENERGY,2X,6HUz'MAX,3X,4HUz'C,3X,4HUr'C,4X,4HUa'C,      &
         4X,6HRE_TAU,7X,2HCf,3X)
      ELSE
!         !do nothing!
      END IF
      
      RETURN
      END  SUBROUTINE WRTSCNINI_tg 


!**********************************************************************
      SUBROUTINE WRTSCNDAT
      
      use flow_info
      use init_info
      use mesh_info
      use postprocess_info
      use wrt_info
      IMPLICIT NONE 

       INTEGER(4) :: FLFLG1
       INTEGER(4) :: J
       
       REAL(WP)    :: Ret  
       REAL(WP)    :: Cf
       
       
       FLFLG1 = 16
       IF(ITERL.EQ.1) THEN
          if (iswitch.eq.1) then
             OPEN(FLFLG1,FILE='SECONDSTAT_Channel.dat')
          else
             OPEN(FLFLG1,FILE='SECONDSTAT_Pipe.dat')
          endif
       ELSE
          if (iswitch.eq.1) then
             OPEN(FLFLG1,FILE='SECONDSTAT_Channel.dat', position='append')
          else
             OPEN(FLFLG1,FILE='SECONDSTAT_Pipe.dat', position='append')
          endif
       END IF
          
          
       CALL PP_WALLSTRESS
       CALL MEANQP12
       CALL VMAV_tg           
       CALL PP_SCRN
          
       IF (MYID.EQ.0) THEN
           Ret    =REN*DSQRT(DABS(CFUW1))
           Cf     =2.0_WP*(DABS(CFUW1))   
           
           WRITE(logflg,1100) ITERG,phyTIME_TG,CPUTIME_tg,MAXDIVGV,                 &
                       CFLMM*DT,UME1U(1),UMX1U(1),                     &
                       IDINT(REN*UME1U(1)),                           &
                       IDINT(REN*UMX1U(1)),                           &
                       ENE11,URMSX1,URMSC,VRMSC,                       &
                       WRMSC,Ret,Cf
           
 1100     FORMAT(1I8,1F12.5,1F7.3,1ES11.3,1F6.2,2F9.5,2I7,5ES11.3,&
                  2ES13.5)
             
           IF ( (ITERG.EQ.2) .OR.                                &
                (DMOD(phyTIME,TSCN).LT.DT) .AND.                     &
                (DMOD(phyTIME,TSCN*100.0_WP).LT.DT) ) THEN
             
                if(iswitch.eq.1) then
                    WRITE(FLFLG1,*) 'ypj,ysj,uprim,vprim,wprim,Prim,Um,Vm,Wm'
                    DO J=1,NCL2
                        WRITE(FLFLG1,'(8F12.5)') YCC(J),YND(J),      &
                        DSQRT(STA13(2,1,J)), DSQRT(STA13(2,2,J)),   &
                        DSQRT(STA13(2,3,J)),                        &
                        STA13(1,1,J),STA13(1,2,J),STA13(1,3,J)
                    ENDDO
                elseif(iswitch.eq.2) then
                    WRITE(FLFLG1,*) 'rmj,rcj,uzprim,urprim,uaprim,Prim,',  &
                     'Uzm,Urm,Uam'
                    DO J=1,NCL2
                        WRITE(FLFLG1,'(8F12.5)') rm(J),rc(J),        &
                        DSQRT(STA13(2,1,J)), DSQRT(STA13(2,2,J)),    &
                        DSQRT(STA13(2,3,J)),                         &
                        STA13(1,1,J),STA13(1,2,J),STA13(1,3,J)
                    ENDDO
                endif
                
           ENDIF
       ENDIF
       
       CLOSE(FLFLG1)
       
       RETURN
      END SUBROUTINE WRTSCNDAT 
      
      !***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani, Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> Print screen data
!
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
! 06/02/2014- Subroutine structured is optimized by Wei Wang

!**********************************************************************
      SUBROUTINE  WRTSCNINI_io
      use init_info
      use mesh_info
      use wrt_info
      IMPLICIT NONE
      
      !CALL date_and_time(DATE=date,TIME=time)
      !fllog=date(1:4)//'.'//date(5:8)//'.'//time(1:4)//'.log'
      !OPEN(logflg_io,FILE='history.inloutflow.'//fllog)
      logflg_io = 6 
      IF (iswitch.eq.1) then
          WRITE(logflg_io,1000) 
 1000     FORMAT('# ',2X,4HStep,5X,4HPhyT,5X,5HCpuTL,2X,5HCpuTT,     &
           2X,3HCFL,3X,8HDivMaxMF,3X,8HDivMaxIn,3X,8HDivMaxOu,3X,1HU,10X,1HV,10X,1HW,10X,1HP,10X,2HCf)
      ELSE IF(iswitch.eq.2) then
          WRITE(logflg_io,1004) 
 1004    FORMAT(1X,'# ',1X,4HSTEP,7X,1HT,8X,6HC_TIME,3X,6HDIVMAX,    &
         5X,3HCFL,4X,5HUMEAN,3X,4HUMAX,4X,3HREM,4X,3HREX,            &
         5X,6HENERGY,2X,6HUz'MAX,3X,4HUz'C,3X,4HUr'C,4X,4HUa'C,      &
         4X,6HRE_TAU,3X,2HCf,3X)
      ELSE
!         !do nothing!1I8,1F12.5,2F7.3,1F6.2, 6ES11.3
      END IF
      
      RETURN
      END  SUBROUTINE WRTSCNINI_io


!**********************************************************************
      SUBROUTINE WRTSCNDAT_io
      
      use flow_info
      use init_info
      use mesh_info
      use postprocess_info
      use wrt_info
      IMPLICIT NONE 

       INTEGER(4) :: FLFLG1
       INTEGER(4) :: J
       
       REAL(WP)    :: U, V, W, P  
       REAL(WP)    :: Cf
       
       

       !CALL PP_WALLSTRESS
       !CALL MEANQP12
       !CALL VMAV           
       !CALL PP_SCRN
       logflg_io = 6    
       IF (MYID.EQ.0) THEN
           !Ret    =REN*DSQRT(DABS(CFUW1))
           Cf     = 2.0_WP*(Q_io(NCL1_io/2,1,NCL3/3,1)-0.0_WP)/(YCC(1)-YND(1))/REN
           U      = Q_io(NCL1_io/2,N2DO(myid),NCL3/2,1)
           V      = Q_io(NCL1_io/2,N2DO(myid),NCL3/2,2)
           W      = Q_io(NCL1_io/2,N2DO(myid),NCL3/2,3)
           P      = PR_io(NCL1_io/2,N2DO(myid),NCL3/2)
           WRITE(logflg_io,'(1I8,1F12.5,2F7.3,1F6.2, 8ES11.3)') &
                       ITERG,phyTIME,CPUTIME_io,CPUTIME,CFLMM*DT, MAXDIVGV_io(1),&
                       MAXDIVGV_io(2),MAXDIVGV_io(3),U, V, W, P, Cf


       ENDIF
       
       
       RETURN
      END SUBROUTINE WRTSCNDAT_io



