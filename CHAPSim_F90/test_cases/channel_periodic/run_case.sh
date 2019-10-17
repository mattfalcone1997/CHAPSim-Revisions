#!/bin/bash
nohup mpirun -np 2 ../../bin/CHAPSim-F90 <readdata.ini >OUTPUT_$(date +%Y-%m-%d_%H.%M).log &
