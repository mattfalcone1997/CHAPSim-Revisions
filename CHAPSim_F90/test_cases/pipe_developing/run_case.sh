#!/bin/bash
nohup mpirun -np 3 ../../bin/Incomp3D-* <readdata.ini >OUTPUT_$(date +%Y-%m-%d_%H.%M).log &
