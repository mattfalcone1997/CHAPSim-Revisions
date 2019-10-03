#################################################################################################
# Makefile for DNSpipe, by Wei Wang, Dec 2013.                                                                                              #
# Usage:                                                                                        #
#       make all        to make all files with -O2                                              #
# For debugging run:                                                                            #
# mpirun -np 4 valgrind --leak-check=full --track-origins=yes \                                 #
#                       --log-file=valgrind_output.txt ./dgdes* < solver_input > solver_output  #
#                                                                                               #
#################################################################################################
FMAKE=make
SHELL  = /bin/sh

VERSION = F90
OUTPUT = CHAPSim-$(VERSION)

SRCDIR= ../src
FFT99DIR= ../libs/fft99
FISHPACK= ../libs/fishpack4.1

FC = mpif90 #-vt       

FOPT = -fdefault-real-8 -fdefault-double-8 # this is for gfortran compiler
#FOPT = -r8 # this is for intel compiler

FOPT1 =  -O2
FOPT2 =  -O2 
FOPT3 =  -O2 

#-------------------------------------------------------------------------------------#
# CAUTION: After this line, NOTHING should be really changed!!!                       #
#-------------------------------------------------------------------------------------#

OBJ1= \
        modules.o\
        MEMALLOCT.o\
        IncompSim_F90.o\
        BCAST_CPARAM.o\
		BCAST_INI.o\
		BCAST_MESH_COEF.o\
		BCAST_MESH_COMM.o\
		BC_COUTLET_MOM_RK3.o\
		BC_INLET.o\
		BDFORCECONFIG.o\
		BDFORCE.o\
		CFL.o\
		CONSPARA.o\
		CONVECTION_X.o\
		CONVECTION_Y.o\
		CONVECTION_Z.o\
		COORDJR.o\
		DIVGCK.o\
		DIVG.o\
		FLOWSTART.o\
		functions.o\
		INDEXL2G.o\
		INDEXSET.o\
		INIFIELD_FLOW.o\
		INTFC_BC_DPH.o\
		INTFC_BC_QP.o\
		LAMPOISLPROF.o\
		LAPLACECOEF.o\
		LAPLACECOEFpipe.o\
		MEANQP12.o\
		MEANVELOINI.o\
		mesh_decomp.o\
		MOMFA1_X.o\
		MOMFA2_Z.o\
		MOMFA3_Y.o\
		MOMFA.o\
		possion3d_xz_io.o\
		possion3d_xz_periodic.o\
		POSTPROESS.o\
		PP_SCRN.o\
		PP_WALLSTRESS.o\
		PRCALC.o\
		RANDOM_FL_FLD.o\
		randomgen.o\
		READINI.o\
		RESTART.o\
		RHS_CvLpGpS.o\
		RKCOEF.o\
		SOLVE.o\
		SOLVERRK3_MOM.o\
		start_mpi.o\
		STATN13.o\
		TDMAIJI_CYC.o\
		TDMAIJI_nonCYC.o\
		TDMAIJJ_nonCYC.o\
		TEC360_DATA_CHECK.o\
		TIMECTL.o\
		TRASP23_PHIIm.o\
		TRASP23_PHIRe.o\
		TRASP23_RHS.o\
		TRASP23_RHS_io.o\
		TRASP23_RHSLLPHI_io.o\
		VALME13.o\
		VALME13ORL.o\
		VELOUPDT.o\
		VMAV.o\
		WRTHDL.o\
		WRTQP.o\
		WRTSCN.o\
		WRTVPROF.o\
             
OBJ2= \
        CFFT99.o\
        CFTFAX.o\
        CFTRIG.o\
        FACT.o\
        FAX.o\
        FFT99A.o\
        FFT99B.o\
        FFT99.o\
        FFTFAX.o\
        FFTRIG.o\
        VPASSM.o\
        
OBJ3= \
        fftpack.o\
   
default:
	@cd bin; make $(OUTPUT) -f ../Makefile

$(OUTPUT): $(OBJ1) $(OBJ2) $(OBJ3)
	@echo -n "Linking... "
	$(FC) -o $(OUTPUT) *.o $(FOPT1) $(FOPT2) $(FOPT3) $(FOPT)
	@echo -e "Done.\a"
                	
$(OBJ1): %.o: $(SRCDIR)/%.f90
	$(FC) -c   $(FOPT1) $(FOPT)$<
	
$(OBJ2): %.o: $(FFT99DIR)/%.f90
	$(FC) -c  $(FOPT2) $(FOPT1) $(FOPT)$<

$(OBJ3): %.o: $(FISHPACK)/%.f
	$(FC) -c  $(FOPT3) $(FOPT)$<

clean:
	rm -f bin/*

#clear
#@echo "All binary files and executeable program deleted."

all:
	@make clean -f Makefile

	@cd bin; make $(OUTPUT) -f ../Makefile
	
	@echo "Code Compiling Is Successful."
