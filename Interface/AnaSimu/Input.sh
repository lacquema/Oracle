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
0. mj


10000000
0
0. mj 0. 0. 0. 0. 0. 0. 
30.00 300000.00
0.20 2000.00
0.00 0.10
exit!