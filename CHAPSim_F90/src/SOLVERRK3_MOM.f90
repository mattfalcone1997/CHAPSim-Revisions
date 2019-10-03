      SUBROUTINE SOLVERRK3_MOM_tg(NS)
      use cparam
      USE WPRECISION
      IMPLICIT NONE

      INTEGER(4),INTENT(IN) :: NS
      INTEGER(4) :: IDR
      INTEGER(4) :: NII


      CALL CONVECTION_X_tg
      CALL CONVECTION_Y_tg
      CALL CONVECTION_Z_tg
      DO IDR=1,3
         CALL RHS_CvLpGpS_tg(NS,IDR)
         CALL MOMFA(NS,IDR)
      END DO

      CALL INTFC_BC_QP_tg
      CALL DIVG_tg(NS)

      CALL FFT99_POIS3D_periodicxz

      CALL INTFC_BC_DPH_tg
      CALL VELOUPDT_tg(NS)

      CALL PRCALC_tg(NS)

      CALL INTFC_BC_QP_tg

      RETURN
      END SUBROUTINE SOLVERRK3_MOM_tg


      SUBROUTINE SOLVERRK3_MOM_io(NS)
      use cparam
      USE WPRECISION
      IMPLICIT NONE

      INTEGER(4),INTENT(IN) :: NS
      INTEGER(4) :: IDR
      INTEGER(4) :: NII


      CALL BC_COUTLET_MOM_RK3(NS)
      CALL INTFC_BC_QP_INOUT
      CALL BC_CBC_TDMA(NS)

      CALL CONVECTION_X_io
      CALL CONVECTION_Y_io
      CALL CONVECTION_Z_io
      DO IDR=1,3
         CALL RHS_CvLpGpS_io(NS,IDR)
         CALL MOMFA_io(NS,IDR)
      END DO

      CALL INTFC_BC_QP_io
      CALL DIVG_io(NS) !poisson equation source term


      CALL FISHPACK_POIS3D_SIMPLE !solving poisson equation
      CALL INTFC_BC_DPH_io
      CALL VELOUPDT_io(NS)

      CALL PRCALC_io(NS)

      CALL INTFC_BC_QP_io 

      RETURN
      END SUBROUTINE SOLVERRK3_MOM_io
