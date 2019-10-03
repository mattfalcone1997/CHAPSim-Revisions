

      FUNCTION DGAL(AL,ZITAF)
      USE WPRECISION
      IMPLICIT NONE
      REAL(WP),INTENT(IN) :: AL
      REAL(WP),INTENT(IN) :: ZITAF
      REAL(WP) :: DGAL
      
      REAL(WP) :: DG1
      REAL(WP) :: DG2
      REAL(WP) :: DG3
      REAL(WP) :: ALZITF

      
      ALZITF=AL*ZITAF
      DG1=1.0_WP/(DSINH(ALZITF)*DCOSH(ALZITF))
      DG2=-ALZITF/DSINH(ALZITF)**2
      DG3=-ALZITF/DCOSH(ALZITF)**2
      DGAL=DG1+DG2+DG3
      
      RETURN
      END
      
      
      
!************************************************
      FUNCTION GAL(AL,ZITAF)
      USE WPRECISION
      IMPLICIT NONE     
      REAL(WP),INTENT(IN) :: AL
      REAL(WP),INTENT(IN) :: ZITAF
      REAL(WP) :: GAL      
      REAL(WP) :: ALZITF
      
      ALZITF=AL*ZITAF
      GAL=AL/(DCOSH(ALZITF)*DSINH(ALZITF))
      
      RETURN
      END
      
!!>>>>>>>>>>>>>>>>>>    Begin of Function   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!!
    logical function ISNAN1(a)
    USE WPRECISION
    IMPLICIT NONE
       REAL(WP) :: a
      if( ((a + 1.0).eq.a) .or. (a/=a) ) then
        isnan1 = .true.
      else
        isnan1 = .false.
      end if 

       
       return
    end function ISNAN1
!!<<<<<<<<<<<<<<<<<<    End of Function     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!



!!>>>>>>>>>>>>>>>>>>    Begin of Function   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!!
    logical function ISINF1(a)
    USE WPRECISION
    IMPLICIT NONE
    
       REAL(WP) :: a
       REAL(WP) :: b
       b = 0.0_WP
       b = HUGE(b)

      if( dabs(a) .GT. dabs(b) ) then
        ISINF1 = .true.
      else
        ISINF1 = .false.
      end if 

       
       return
    end function ISINF1
!!<<<<<<<<<<<<<<<<<<    End of Function     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!

    REAL(WP) function logbaseb(a,b)
    USE WPRECISION
    IMPLICIT NONE
    
    REAL(WP) :: a
    REAL(WP) :: b

    logbaseb = log10(a)/log10(b)
     
    return
    end function logbaseb
