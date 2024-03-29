      SUBROUTINE BCAST_MESH_COEF
     
      use mesh_info
      use init_info
      IMPLICIT NONE

      CALL MPI_BCAST( CLGIP, NCL2, MPI_INTEGER4,0,ICOMM,IERROR )
      CALL MPI_BCAST( IMV, NCL1, MPI_INTEGER4,0,ICOMM,IERROR )
      CALL MPI_BCAST( IPV, NCL1, MPI_INTEGER4,0,ICOMM,IERROR )
      CALL MPI_BCAST( KMV, NCL3, MPI_INTEGER4,0,ICOMM,IERROR )
      CALL MPI_BCAST( KPV, NCL3, MPI_INTEGER4,0,ICOMM,IERROR )      
      
      CALL MPI_BCAST( isym, NCL3, MPI_INTEGER4,0,ICOMM,IERROR )
      CALL MPI_BCAST( jmv,  NCL2, MPI_INTEGER4,0,ICOMM,IERROR )
      CALL MPI_BCAST( jpv,  NCL2, MPI_INTEGER4,0,ICOMM,IERROR )

      CALL MPI_BCAST( YND, NND2, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR )
      
      CALL MPI_BCAST( YCC, NCL2, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR )
      
      CALL MPI_BCAST( RC,  NND2, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR )
      CALL MPI_BCAST( RM,  NND2, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR )
      
      CALL MPI_BCAST( ACPH, NCL2, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR )
      CALL MPI_BCAST( AMPH, NCL2, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR )
      CALL MPI_BCAST( APPH, NCL2, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR )
      
      CALL MPI_BCAST( DYCI, NND2, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR )
      CALL MPI_BCAST( DYFI, NCL2, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR )
            
      CALL MPI_BCAST( Vini, NCL2, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR )

      CALL MPI_BCAST( ACVR, NND2*3, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR )
      CALL MPI_BCAST( AMVR, NND2*3, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR )
      CALL MPI_BCAST( APVR, NND2*3, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR ) 
      
      IF(IOFLOWflg) THEN
         CALL MPI_BCAST( IMV_io(1), NCL1_io+1, MPI_INTEGER4,0,ICOMM,IERROR )
         CALL MPI_BCAST( IPV_io(0), NCL1_io+1, MPI_INTEGER4,0,ICOMM,IERROR )
      END IF
                             
      RETURN
      END SUBROUTINE BCAST_MESH_COEF
