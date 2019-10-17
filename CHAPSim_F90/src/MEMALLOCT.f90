      SUBROUTINE MEM_ALLOCAT
      use mesh_info
      use init_info
      use flow_info
      use postprocess_info
      IMPLICIT NONE

      ALLOCATE ( JCL2G( N2DO(0) ) )
      ALLOCATE ( KCL2G( N3DO(0)) )
      ALLOCATE ( CLGIP( NCL2) )
      JCL2G = 0
      KCL2G = 0
      CLGIP = 0

      ALLOCATE ( IPV( NCL1 ) )
      ALLOCATE ( IMV( NCL1 ) )
      ALLOCATE ( KPV(NCL3) )
      ALLOCATE ( KMV(NCL3) )
      ALLOCATE ( jmv(NCL2) )
      ALLOCATE ( jpv(NCL2) )
      ALLOCATE ( isym(NCL3) )
      IPV = 0
      IMV = 0
      KPV = 0
      KMV = 0
      JMV = 0
      JPV = 0
      ISYM= 0

      ALLOCATE ( AMPH(NCL2) )
      ALLOCATE ( ACPH(NCL2) )
      ALLOCATE ( APPH(NCL2) )

      ALLOCATE ( AMVR(NND2,3) )
      ALLOCATE ( ACVR(NND2,3) )
      ALLOCATE ( APVR(NND2,3) )

      ALLOCATE ( XND(NND1) )
      ALLOCATE ( YND(NND2) )
      ALLOCATE ( ZND(NND3) )

      ALLOCATE ( rm(NND2) )
      ALLOCATE ( rc(NND2) )

      ALLOCATE ( YCC(NCL2) )

      ALLOCATE ( DYFI(NCL2) )
      ALLOCATE ( DYCI(NND2) )

      ALLOCATE ( Vini(NCL2) )
      Vini = 0.0_WP

      AMPH  = 0.0_WP
      ACPH  = 0.0_WP
      APPH  = 0.0_WP
      AMVR  = 0.0_WP
      ACVR  = 0.0_WP
      APVR  = 0.0_WP
      XND   = 0.0_WP

      YND   = 0.0_WP
      ZND   = 0.0_WP
      rm    = 0.0_WP
      rc    = 0.0_WP
      YCC   = 0.0_WP
      DYFI  = 0.0_WP
      DYCI  = 0.0_WP
      !This implies that the MPI is parallelised in the wall normal direction
      !for turbulence generator
      !Choice of ND2O(0) potentially as that is the largest array hence ensures
      !there is sufficient mememory allocation
      ALLOCATE ( Q   (NCL1,0:N2DO(0)+1,NCL3,NDV) )
      ALLOCATE ( PR  (NCL1,0:N2DO(0)+1,NCL3    ) )
      ALLOCATE ( Qtmp(NCL1,0:N2DO(0)+1,NCL3    ) )
      ALLOCATE ( DPH (NCL1,0:N2DO(0)+1,NCL3    ) )

      ALLOCATE ( F   ( NCL1,NND2,N3DO(0) )     )

      ALLOCATE ( CONVH0   (NCL1,N2DO(0),NCL3,NDV)  )
      ALLOCATE ( RHS      (NCL1,N2DO(0),NCL3)      )
      ALLOCATE ( RHSLLPHI (NCL1,N2DO(0),NCL3)      )

      Q       = 0.0_WP
      PR      = 0.0_WP
      Qtmp    = 0.0_WP
      DPH     = 0.0_WP
      F       = 0.0_WP

      CONVH0  = 0.0_WP
      RHSLLPHI= 0.0_WP
      RHS     = 0.0_WP

      ALLOCATE ( STA13(6,4,NCL2) )
      STA13   = 0.0_WP

      ALLOCATE ( BSEN_F(NCL1*NCL3,NDV+1) )
      ALLOCATE ( BSEN_L(NCL1*NCL3,NDV+1) )
      ALLOCATE ( BREC_F(NCL1*NCL3,NDV+1) )
      ALLOCATE ( BREC_L(NCL1*NCL3,NDV+1) )
      BSEN_F = 0.0_WP
      BSEN_L = 0.0_WP
      BREC_F = 0.0_WP
      BREC_L = 0.0_WP

      UMX1U = 0.0_WP
      UME1U = 0.0_WP


      IF(IOFLOWflg) THEN !If there is a statially developing region
         ALLOCATE ( IPV_io(0:NCL1_io   ) )
         ALLOCATE ( IMV_io(1:NCL1_io+1 ) )
         IPV_io = 0
         IMV_io = 0
         ALLOCATE ( XND_io( NND1_io ) )
         XND_io= 0.0_WP

         ALLOCATE ( Q_io   (0:NCL1_io+1,0:N2DO(0)+1,NCL3,NDV) )
         ALLOCATE ( PR_io  (0:NCL1_io+1,0:N2DO(0)+1,NCL3    ) )
         ALLOCATE ( Qtmp_io(0:NCL1_io+1,0:N2DO(0)+1,NCL3    ) )
         ALLOCATE ( DPH_io (0:NCL1_io+1,0:N2DO(0)+1,NCL3    ) )

         ALLOCATE ( F_io   ( NCL1_io,NND2,N3DO(0) )     )

         ALLOCATE ( CONVH0_io   (NCL1_io,N2DO(0),NCL3,NDV)  )
         ALLOCATE ( RHS_io      (NCL1_io,N2DO(0),NCL3)      )
         ALLOCATE ( RHSLLPHI_io (NCL1_io,N2DO(0),NCL3)      )
         ALLOCATE ( RHSLLPHI_io_tmp (NCL1_io,N2DO(0),NCL3)      )

         Q_io       = 0.0_WP
         PR_io      = 0.0_WP
         Qtmp_io    = 0.0_WP
         DPH_io     = 0.0_WP
         F_io       = 0.0_WP

         CONVH0_io  = 0.0_WP
         RHSLLPHI_io= 0.0_WP
         RHS_io     = 0.0_WP

         RHSLLPHI_io_tmp = 0.0_WP

         ALLOCATE ( BSEN_F_io((NCL1_io+2)*NCL3,NDV+1) )
         ALLOCATE ( BSEN_L_io((NCL1_io+2)*NCL3,NDV+1) )
         ALLOCATE ( BREC_F_io((NCL1_io+2)*NCL3,NDV+1) )
         ALLOCATE ( BREC_L_io((NCL1_io+2)*NCL3,NDV+1) )
         BSEN_F_io = 0.0_WP
         BSEN_L_io = 0.0_WP
         BREC_F_io = 0.0_WP
         BREC_L_io = 0.0_WP

         ALLOCATE ( STA13_io(6,4,NCL2) )
         STA13_io   = 0.0_WP

          ALLOCATE ( BC_CONV0(      N2DO(0),NCL3,NDV  ) )
          ALLOCATE ( BC_TDMA (3,0:N2DO(0)+1,NCL3,NDV  ) )
          BC_CONV0 = 0.0_WP
          BC_TDMA  = 0.0_WP
          ALLOCATE ( BC_U_SSTAR(      N2DO(0),NCL3,NDV  ) )
          BC_U_SSTAR = 0.0_WP

      END IF


      END SUBROUTINE MEM_ALLOCAT


!**********************************************************************
      SUBROUTINE MEM_DEALLOCAT
      use mesh_info
      use init_info
      use flow_info
      use postprocess_info
      IMPLICIT NONE


      DEALLOCATE ( JDSWT )
      DEALLOCATE ( JDEWT )
      DEALLOCATE ( KDSWT )
      DEALLOCATE ( KDEWT )

      DEALLOCATE ( JCL2G )
      DEALLOCATE ( KCL2G )

      DEALLOCATE ( IPV )
      DEALLOCATE ( IMV )

      DEALLOCATE ( KPV )
      DEALLOCATE ( KMV )
      DEALLOCATE ( jmv )
      DEALLOCATE ( jpv )
      DEALLOCATE ( isym )

      DEALLOCATE ( AMPH )
      DEALLOCATE ( ACPH )
      DEALLOCATE ( APPH )

      DEALLOCATE ( AMVR )
      DEALLOCATE ( ACVR )
      DEALLOCATE ( APVR )

      DEALLOCATE ( YND )

      DEALLOCATE ( rm )
      DEALLOCATE ( rc )

      DEALLOCATE ( YCC )

      DEALLOCATE ( DYFI )
      DEALLOCATE ( DYCI )

      DEALLOCATE ( RHS  )
      DEALLOCATE ( Q    )
      DEALLOCATE ( PR   )
      DEALLOCATE ( F    )
      DEALLOCATE ( Qtmp )

      DEALLOCATE ( CONVH0 )
      DEALLOCATE ( DPH    )

      DEALLOCATE ( RHSLLPHI   )



      DEALLOCATE ( BSEN_F )
      DEALLOCATE ( BSEN_L )
      DEALLOCATE ( BREC_F )
      DEALLOCATE ( BREC_L )

      DEALLOCATE ( Vini )

      DEALLOCATE ( STA13 )

      IF(IOFLOWflg) THEN
          DEALLOCATE ( IPV_io )
          DEALLOCATE ( IMV_io )
          DEALLOCATE ( RHS_io  )
          DEALLOCATE ( Q_io    )
          DEALLOCATE ( PR_io   )
          DEALLOCATE ( F_io    )
          DEALLOCATE ( Qtmp_io )

          DEALLOCATE ( CONVH0_io )
          DEALLOCATE ( DPH_io    )

          DEALLOCATE ( RHSLLPHI_io   )
          DEALLOCATE ( BSEN_F_io )
          DEALLOCATE ( BSEN_L_io )
          DEALLOCATE ( BREC_F_io )
          DEALLOCATE ( BREC_L_io )

          DEALLOCATE ( STA13_io )

      END IF

      RETURN
      END SUBROUTINE MEM_DEALLOCAT
