
      SUBROUTINE INDEXSET
      use mesh_info
      IMPLICIT NONE
      !=========================================================================
      !Linking nodes
      !=========================================================================
      INTEGER(4) :: IC
      INTEGER(4) :: JC
      INTEGER(4) :: K, KC
      !Find neighbouring nodes in the x direction (streamwise)
      DO IC=1,NCL1
         IPV(IC)=IC+1
         IMV(IC)=IC-1
         IF(IC.EQ.NCL1) IPV(IC)=1 !If the last streamwise cell
         IF(IC.EQ.1)    IMV(IC)=NCL1 !If the first - could relate to periodicity
      END DO

      !Same for developing flow
      IF(IOFLOWflg) THEN
         DO IC=0,NCL1_io
            IPV_io(IC)=IC+1
         END DO
         DO IC=1,NCL1_io+1
            IMV_io(IC)=IC-1
         END DO
      END IF
      !Find neighbouring nodes in the z direction (Spanwise)
      DO KC=1,NCL3
         KPV(KC)=KC+1
         KMV(KC)=KC-1
         IF(KC.EQ.NCL3) KPV(KC)=1
         IF(KC.EQ.1)    KMV(KC)=NCL3
      END DO

      do k=1,NCL3
         isym(k) = k + NCL3/2
         if(isym(k).gt.NCL3) isym(k) = isym(k) - NCL3
      enddo

      do jc=1,NCL2
         jmv(jc)=jc-1
         jpv(jc)=jc+1
         if(jc.eq.1)    jmv(jc)=jc !First and last seem to link to itself?
         if(jc.eq.NCL2) jpv(jc)=jc
      END DO

      RETURN

      END SUBROUTINE INDEXSET
