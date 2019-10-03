      SUBROUTINE TDMAIJI_CYC(A,B,C,R,IS,ISZ,JS,JSZ)  
      USE WPRECISION
      IMPLICIT NONE                                  

      INTEGER(4),INTENT(IN)   :: IS
      INTEGER(4),INTENT(IN)   :: ISZ
      INTEGER(4),INTENT(IN)   :: JS
      INTEGER(4),INTENT(IN)   :: JSZ
      REAL(WP),INTENT(IN)      :: A(IS:ISZ+IS-1,JS:JSZ+JS-1)
      REAL(WP),INTENT(IN)      :: B(IS:ISZ+IS-1,JS:JSZ+JS-1)
      REAL(WP),INTENT(IN)      :: C(IS:ISZ+IS-1,JS:JSZ+JS-1)
      REAL(WP),INTENT(INOUT)   :: R(IS:ISZ+IS-1,JS:JSZ+JS-1)
      
      INTEGER(4)  :: I, IEND
      INTEGER(4)  :: J, JEND
      REAL(WP)     :: PP
      REAL(WP)     :: QQ1, QQ2
      REAL(WP)     :: H0(IS:ISZ+IS-1,JS:JSZ+JS-1)
      REAL(WP)     :: G1(IS:ISZ+IS-1,JS:JSZ+JS-1)
      REAL(WP)     :: G2(IS:ISZ+IS-1,JS:JSZ+JS-1)
      
      H0 = 0.0_WP
      G1 = 0.0_WP
      G2 = 0.0_WP
      
      IEND = IS+ISZ-1
      JEND = JS+JSZ-1
      

       I = IS
       DO J=JS,JEND
          H0(I,J) =  C(I,J)/B(I,J)
          G1(I,J) =  R(I,J)/B(I,J)
          G2(I,J) = -A(I,J)/B(I,J)
       END DO
      
       DO I=IS+1,IEND-2
          DO J=JS,JEND
             PP = 1.0_WP / ( B(I,J) - H0(I-1,j)*A(I,J) )
             QQ1 = G1(I-1,J) * A(I,J)
             QQ2 = G2(I-1,J) * A(I,J)
             H0(I,J) = C(I,J) * PP
             G1(I,J) = ( R(I,J) - QQ1 ) * PP
             G2(I,J) = -QQ2 * PP
          END DO
       END DO
       
      I = IEND - 1
      DO J=JS, JEND
         PP = 1.0_WP / ( B(I,J) - H0(I-1,j)*A(I,J) )
         QQ1 = G1(I-1,J) * A(I,J)
         QQ2 = G2(I-1,J) * A(I,J)
         G1(I,J) = (  R(I,J) - QQ1 ) * PP
         G2(I,J) = ( -C(I,J) - QQ2 ) * PP
      END DO 
!            
      DO I = IEND-2, IS, -1
         DO J = JS, JEND
            G1(I,J) = -H0(I,J) * G1(I+1,J) + G1(I,J)
            G2(I,J) = -H0(I,J) * G2(I+1,J) + G2(I,J)
         END DO
      END DO

      I = IEND
      DO J=JS, JEND
         PP = R(I,J) - C(I,J) * G1(IS,J) - A(I,J) * G1(I-1,J)
         QQ1= B(I,J) + C(I,J) * G2(IS,J) + A(I,J) * G2(I-1,J)
         R(I,J) = PP / QQ1
      END DO      

      DO I = IEND-1, IS, -1
         DO J = JS, JEND
            R(I,J) = G1(I,J) + G2(I,J) * R(IEND,J)
         END DO
      END DO      
                                                           
      RETURN                                                            
      END SUBROUTINE TDMAIJI_CYC
