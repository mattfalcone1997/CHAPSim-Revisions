!/***************************************************************************************************
!>Document of the Code "IncompSim_F90"
!>
!>INTRODUCTION
!>  "IncompSim_F90" is a branch of the code "IncompSim", which simulates steady/transient fully
!>   developed/spacial developping channel/pipe flows.
!>  "IncompSim" is a finite-difference based Direct Numerical Simulation code. A brief
!>   flowchart of the code development and branches is given below:
!>                        IncompSim_V1
!>       ______________________|____________________________
!>       |                     |                            |
!>   IncompSim_LES       IncompSim_Pipe                IncompSim_IBM
!>                             |
!>                       IncompSim_F90
!>              _______________|__________________
!               |                                |
!>        IncompSim_TH                    IncompSim_BL
!>
!>CODE DEVELOPERS
!>  -Main Developers
!>    * Dr Mehdi Seddighi
!>
!>  -Co-Developers
!>    * Dr Wei Wang
!>
!>  -Other Developers
!>    * Kui He
!>
!>RESEARCH GROUP
!>  - Heat, Flow and Turbulence Research Group, (http://www.sheffield.ac.uk/heft)
!>    Department of Mechanical Engineering,
!>    The University of Sheffield
!>
!>PI
!>  - Professor Shuisheng He
!>
!>CONTACT
!>  - Dr Mehdi Seddighi, seddighi@sheffield.ac.uk
!>  - Dr Wei Wang, wei.wang@sheffield.ac.uk
!>  - Kui He, mep11kh@sheffield.ac.uk
!>  - Professor Shuisheng He, s.he@sheffield.ac.uk

!>CODE DEVELOPMENT HISTORY
!>  -IncompSim_V1
!>     * A DNS code, discretized on a rectangular staggered mesh with spatial derivatives
!>       formulated using a second-order central finite-difference method, to solve N-S
!>       equations for a channel flow. For the temporal discretization, an explicit low
!>       storage, third-order Runge–Kutta scheme and a second order implicit Crank–Nicholson
!>       scheme are used for the nonlinear and the viscous terms, respectively. These are
!>       combined with the fractional-step method to enforce continuity. The Poisson equation
!>       is solved using FFT. The code is parallelised using the Message-Passing Interface (MPI).
!>       (Mehdi Seddighi, 2011)
!>
!>  -IncompSim_Pipe
!>     * Extended to include cylindrical Coordinates for pipe flow. (Kui He, 2013)
!>
!>  -IncompSim_F90
!>     * Transferred the code from F77 to F90 and optimized code structure/memory (Wei Wang, 2013)
!>     * Developed the developing flow solver; included a separate turbulence generator.
!>      (Wei Wang, 2014)
!>
!>REFERENCES
!>  - Seddighi, Mehdi. (2011) Study of turbulence and wall shear stress in unsteady flow over
!>    smooth and rough wall surfaces. PhD thesis, University of Aberdeen.
!>  - He, Shuisheng. & Seddighi, Mehdi. (2013) Turbulence in transient channel flow. J. Fluid Mech.
!>    715, 60–102.
!>
!>RESEARCH GROUP
!>  - Heat, Flow and Turbulence Research Group, (http://www.sheffield.ac.uk/heft)
!>    Department of Mechanical Engineering,
!>    The University of Sheffield
!>
!>PI
!>  - Professor Shuisheng He
!>
!>CONTACT
!>  - Dr Mehdi Seddighi, seddighi@sheffield.ac.uk
!>  - Dr Wei Wang, wei.wang@sheffield.ac.uk
!>  - Kui He, mep11kh@sheffield.ac.uk
!>  - Professor Shuisheng He, s.he@sheffield.ac.uk
!>
!>TERMS AND CONDITIONS
!>
!>  Copyright (C) 2015  The University Of Sheffield
!>
!>  This program is free software; you can redistribute it and/or modify it under the terms of
!>  the GNU General Public License as published by the Free Software Foundation; either version
!>  2 of the License, or (at your option) any later version. This program is distributed in the
!>  hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!>  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
!>  more details. You should have received a copy of the GNU General Public License along with
!>  this program (http://www.gnu.org/licenses/gpl-2.0.html); if not, write to the Free Software
!>  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. Also add
!>  information on how to contact you by electronic and paper mail.
!>
!>  Head of Heat, Flow and Turbulence Research Group in the department of Mechanical Engineering,
!>  The University of Sheffield, Shuisheng He, hereby disclaims all copyright interest in the CFD
!>  code 'IncompSim'.
!>
!>  Shuisheng He,  1 May, 2015
!>  Head of Heat, Flow and Turbulence Research Group in Department of Mechanical Engineering,
!>  The University of Sheffield
!>
!>*/*************************************************************************************************

       PROGRAM IncomSim_F90
       use init_info
       IMPLICIT NONE

       call START_MPI !initailises MPI
       IF (MYID.EQ.0) CALL CHKHDL('****DNS SOVLER STARTING******',myid)

       IF (MYID.eq.0) CALL CHKHDL('1. READ INI FILE...',myid)
       IF (MYID.EQ.0) THEN
           CALL READINI !Reads input file
           CALL CHKHDL('2. MESH DECOMPOSITION...',myid)
           CALL mesh_Ydecomp_master
           CALL mesh_Zdecomp_master
       END IF
       CALL mesh_decomp_bcast !Broadcasts mesh info
       CALL BCAST_INI !Broadcasts ini file info
       CALL BCAST_CPARAM  !Broadcasts mesh parameter e.g. nodes in each direction
       CALL MEM_ALLOCAT !Allocates memory

       CALL INDEXL2G !Creates an array on ranks indiciting which index in on each rank

       IF (MYID.EQ.0) THEN
           CALL CHKHDL('3. SET UP RK COEFFICIENTS...',myid)
           CALL RKCOEF !Presumably setup Runge Kutta coefficients

           CALL CHKHDL('4. SET UP +1/-1 INDEX...',myid)
           CALL INDEXSET  !Links nodes

           CALL CHKHDL('5. SET UP CONST MESH INFO IN X, Z DIRECTIONS...',myid)
           CALL CONSPARA

           CALL CHKHDL('6. SET UP MESH INFO IN Y...',myid)
           CALL COORDJR

           CALL CHKHDL('7. SET UP INITIAL VELOCITY PROFILE...',myid)
           CALL LAMPOISLPROF

           CALL CHKHDL('8. SET UP COEFFICIENTS FOR LAPLACE EQ...',myid)
           if(iswitch.eq.1) then
             CALL LAPLACECOEF
           else
             call LAPLACECOEFpipe
           endif

       END IF

       IF (MYID.EQ.0) CALL CHKHDL('9. BCAST PARAMETERS AND MESHES...',myid)
       CALL BCAST_MESH_COMM !broadcast constants such as mesh length and spacing
       CALL BCAST_MESH_COEF !broadcast coefficients related to the mesh such as BCFI

       IF (MYID.EQ.0) CALL CHKHDL('10.INITILIZATING FFT',myid)
       IF(IOFLOWflg) THEN
          CALL FFT99_POIS3D_INIT
          CALL FISHPACK_POIS3D_INIT
       ELSE
          CALL FFT99_POIS3D_INIT
       END IF

       IF (MYID.EQ.0) CALL CHKHDL('11.THE MAIN SOLVER STARTS...',myid)
       CALL SOLVE

       CALL MEM_DEALLOCAT
       CALL MPI_FINALIZE(IERROR)
       STOP

       END
