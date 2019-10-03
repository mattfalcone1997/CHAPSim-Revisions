!************************************************************************
    MODULE WPRECISION
    
      INTEGER,PARAMETER  :: WP=KIND(0.0D0) !WORKING PRECESION
      
    END MODULE WPRECISION
!************************************************************************

    MODULE cparam
           
      INTEGER(4),PARAMETER ::  NDV=3 
!>    @param NDV is three directions
  
      INTEGER(4) ::  NCL1
      INTEGER(4) ::  NCL2
      INTEGER(4) ::  NCL3
      INTEGER(4) ::  NCL1_io

      INTEGER(4) ::  NND1
      INTEGER(4) ::  NND2
      INTEGER(4) ::  NND3 
      INTEGER(4) ::  NND1_io                  
      
      LOGICAL    ::  IOFLOWflg
      
      END MODULE cparam
!********************************MODULE**********************************
      MODULE wrt_info
      character(8)  :: date
      character(10) :: time
      
      character(18)   :: fllog
      integer(4)      :: logflg = 1023
      integer(4)      :: logflg_io = 1022

      
      END MODULE wrt_info
!********************************MODULE**********************************

!********************************MODULE**********************************
!********************************MODULE**********************************
      MODULE mpi_info
        include 'mpif.h'
     
        INTEGER(4) ::  MYID
        INTEGER(4) ::  IERROR
        INTEGER(4) ::  SIZE
        INTEGER(4) ::  ICOMM
        integer(4) ::  STS(MPI_STATUS_SIZE)
     
        INTEGER(4) ::  NPSLV
        INTEGER(4) ::  NPTOT
       
      END MODULE mpi_info
     
!********************************MODULE**********************************
!********************************MODULE**********************************
      MODULE mesh_info
        use cparam 
        use mpi_info
        use WPRECISION
        
        INTEGER(4),ALLOCATABLE  :: N2DO(:)
        INTEGER(4),ALLOCATABLE  :: N3DO(:)
        
        INTEGER(4),ALLOCATABLE  :: JDSWT(:)
        INTEGER(4),ALLOCATABLE  :: JDEWT(:)
        INTEGER(4),ALLOCATABLE  :: KDSWT(:)
        INTEGER(4),ALLOCATABLE  :: KDEWT(:)
        INTEGER(4),ALLOCATABLE  :: JCL2G(:)
        INTEGER(4),ALLOCATABLE  :: KCL2G(:)
        INTEGER(4),ALLOCATABLE  :: CLGIP(:)
        
        INTEGER(4),ALLOCATABLE  :: IPV(:)
        INTEGER(4),ALLOCATABLE  :: IMV(:)
        INTEGER(4),ALLOCATABLE  :: IPV_io(:)
        INTEGER(4),ALLOCATABLE  :: IMV_io(:)
        INTEGER(4),ALLOCATABLE  :: KPV(:)
        INTEGER(4),ALLOCATABLE  :: KMV(:)
        INTEGER(4),ALLOCATABLE  :: jmv(:)
        INTEGER(4),ALLOCATABLE  :: jpv(:)
        INTEGER(4),ALLOCATABLE  :: isym(:)
        
        
        REAL(WP)      :: HX,HZ
        REAL(WP)      :: ALX1,ALX2,ALX3
        REAL(WP)      :: DX, DXI, DXQI, QDX1
        REAL(WP)      :: DZ, DZI, DZQI, QDX3 
        REAL(WP)      :: VL1313, VL1313_io
        
        REAL(WP),ALLOCATABLE  :: AMPH(:)
        REAL(WP),ALLOCATABLE  :: ACPH(:)
        REAL(WP),ALLOCATABLE  :: APPH(:)
        
        REAL(WP),ALLOCATABLE  :: AMVR(:,:)
        REAL(WP),ALLOCATABLE  :: ACVR(:,:)
        REAL(WP),ALLOCATABLE  :: APVR(:,:)
        
        REAL(WP),ALLOCATABLE  :: XND(:)
        REAL(WP),ALLOCATABLE  :: XND_io(:)
        REAL(WP),ALLOCATABLE  :: YND(:)
        REAL(WP),ALLOCATABLE  :: ZND(:)
            
        REAL(WP),ALLOCATABLE  :: rm(:)
        REAL(WP),ALLOCATABLE  :: rc(:) 
      
        REAL(WP),ALLOCATABLE  :: YCC(:)
        
        REAL(WP),ALLOCATABLE  :: DYFI(:)
        REAL(WP),ALLOCATABLE  :: DYCI(:)
        
        
      END MODULE mesh_info


!********************************MODULE**********************************
      MODULE BC_info
      use cparam
      use WPRECISION
      REAL(WP)  :: U_OUTLET
      REAL(WP), ALLOCATABLE :: BC_CONV0(:,:,:)
      REAL(WP), ALLOCATABLE :: BC_TDMA (:,:,:, :)
      REAL(WP), ALLOCATABLE :: BC_U_SSTAR(:,:,:)
      
      END MODULE BC_info


!********************************MODULE**********************************
      MODULE init_info
        use cparam
        use mpi_info
        use WPRECISION
        INTEGER(4)  :: ITERG0
        INTEGER(4)  :: ITERG0_TG
        INTEGER(4)  :: ITERL
        INTEGER(4)  :: ITERG
        INTEGER(4)  :: NSST
        INTEGER(4)  :: NFLOW
        INTEGER(4)  :: N_START
        INTEGER(4)  :: MULTIM
        INTEGER(4)  :: NREAD, NREAD_io
        INTEGER(4)  :: IPWALL1,IPWALL2
        INTEGER(4)  :: iswitch
        INTEGER(4)  :: JWLC1
        INTEGER(4)  :: JWGL1
        INTEGER(4)   :: ISTR2
        
        INTEGER(4)   :: FLOWTP
        
        REAL(WP)      :: PI
        REAL(WP)      :: DT
        
        REAL(WP)      :: CPUTIME
        REAL(WP)      :: CPUTIME_tg
        REAL(WP)      :: CPUTIME_io
        REAL(WP)      :: CFLGV
        REAL(WP)      :: phyTIME
        REAL(WP)      :: phyTIME_TG
        REAL(WP)      :: TSCN
        REAL(WP)      :: TRST
        REAL(WP)      :: TRST_io
        REAL(WP)      :: TSTOP
        REAL(WP)      :: tbody
        REAL(WP)      :: TLGRE

        REAL(WP)      :: TSAVE      
        REAL(WP)      :: TSAVE1  
        REAL(WP)      :: TSTAV1

        REAL(WP)      :: TTECCK

        REAL(WP)      :: RATEM1
        REAL(WP)      :: CFGV
        REAL(WP)      :: DPCONS
   
        REAL(WP)      :: REINI
        REAL(WP)      :: REN
        REAL(WP)      :: CVISC
 
        REAL(WP)      :: VPER
        REAL(WP)      :: SVPER
        
        REAL(WP)      :: STR2
        REAL(WP)      :: TGAM(0:3), TROH(0:3), TALP(0:3)
        
        REAL(WP),ALLOCATABLE :: Vini(:)
        
        
        INTEGER(4)    :: BCX(2)
        INTEGER(4)    :: BCZ(2)
        
        INTEGER(4)    :: BCX_io(2)
        
      
      END MODULE init_info
      
!********************************MODULE**********************************
      MODULE flow_info
        use cparam
        use mpi_info
        use BC_info
        use WPRECISION
        REAL(WP)      :: MAXDIVGV
        REAL(WP)      :: CFLMM
        REAL(WP)      :: VMV(3)
        REAL(WP)      :: COE
        
        REAL(WP),ALLOCATABLE :: RHS (:,:,:)
        REAL(WP),ALLOCATABLE :: Q   (:,:,:,:)
        REAL(WP),ALLOCATABLE :: PR  (:,:,:)
        REAL(WP),ALLOCATABLE :: F   (:,:,:)
        REAL(WP),ALLOCATABLE :: Qtmp(:,:,:)
        
        REAL(WP),ALLOCATABLE :: CONVH0(:,:,:,:)
        REAL(WP),ALLOCATABLE :: DPH   (:,:,:)
        
        REAL(WP),ALLOCATABLE :: RHSLLPHI  (:,:,:)    
        
        REAL(WP),ALLOCATABLE :: BSEN_F(:,:)
        REAL(WP),ALLOCATABLE :: BSEN_L(:,:)
        REAL(WP),ALLOCATABLE :: BREC_F(:,:)
        REAL(WP),ALLOCATABLE :: BREC_L(:,:)
        
        REAL(WP),ALLOCATABLE :: STA13(:,:,:)
        
        
        REAL(WP)      :: MAXDIVGV_io(3)
        REAL(WP)      :: VMV_io(3)
        REAL(WP)      :: CFLMM_io
        
        REAL(WP),ALLOCATABLE :: RHS_io (:,:,:)
        REAL(WP),ALLOCATABLE :: Q_io   (:,:,:,:)
        REAL(WP),ALLOCATABLE :: PR_io  (:,:,:)
        REAL(WP),ALLOCATABLE :: F_io   (:,:,:)
        REAL(WP),ALLOCATABLE :: Qtmp_io(:,:,:)
        
        REAL(WP),ALLOCATABLE :: CONVH0_io(:,:,:,:)
        REAL(WP),ALLOCATABLE :: DPH_io   (:,:,:)
        
        REAL(WP),ALLOCATABLE :: RHSLLPHI_io  (:,:,:)
        REAL(WP),ALLOCATABLE :: RHSLLPHI_io_tmp  (:,:,:)    
        
        REAL(WP),ALLOCATABLE :: BSEN_F_io(:,:)
        REAL(WP),ALLOCATABLE :: BSEN_L_io(:,:)
        REAL(WP),ALLOCATABLE :: BREC_F_io(:,:)
        REAL(WP),ALLOCATABLE :: BREC_L_io(:,:)

        REAL(WP),ALLOCATABLE :: STA13_io(:,:,:)
        
      END MODULE flow_info
           
!********************************MODULE**********************************      
      MODULE postprocess_info
        use cparam
        use mpi_info
        use WPRECISION
        INTEGER(4)  :: NAV11
        
        
        REAL(WP)      :: CFUW1
        REAL(WP)      :: CFLW1
        REAL(WP)      :: CFAV
        REAL(WP)      :: URMSX1
        REAL(WP)      :: ENE11
        REAL(WP)      :: URMSC,VRMSC,WRMSC
        
        REAL(WP)      :: UMX1U(9)
        REAL(WP)      :: UME1U(9)

      END MODULE postprocess_info 
      
!********************************MODULE**********************************      
