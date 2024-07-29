#! /bin/bash
export OMP_NUM_THREADS=8
export STACKSIZE=1000000
./astrom_mcmcop <<!
2
0 0 0 0
1
1d-7
0
0.00 pc
0.00 Msun


10000000
1
0 0 
Distribution