#! /bin/bash
export OMP_NUM_THREADS=8
export STACKSIZE=1000000
/Users/lacquema/Oracle.env/Oracle/Algorithm/bin/astrom_mcmco <<!
1 # Simulation continuation
dump.dat
10000000 # Dump frequency
exit
!