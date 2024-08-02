#! /bin/bash
export OMP_NUM_THREADS=8
export STACKSIZE=1000000
/home/lacquema/OrbitsFits/Oracle/Interface/NewSimu/../../Algorithm/bin/astrom_mcmcop <<!
2 ! New simulation LM and MCMC
1 0 1 0 ! Type of data (RelAstro? AbsRV? RelRV? AbsAstro?)
1 ! Number of orbits
data.txt ! Data file
1d-7 ! Precision
1 1 1 ! Format data (1:DDMMYYYY/2:JD 1:(DEC,RA)/2:(SEP,PA) CorrCoeff?
1 ! Jitter?
0.00 pc ! Distance
0.0 ! First guess of center mass (ms)
adjustment.dat ! Result file
dump.dat ! Dump file
10000000 ! Dump frequency
0 ! Number of masses prior
0 ! Reference of time
0.00 0.00 ! Initial VO and Jitter
0. mj 0. 0. 0. 0. 0. 0.  ! First guess of orbit parameters (m[mj] a[AU] e i[deg] Om[deg] om[deg] tp[MJD])
30.00 300000.00 ! Range of permited period
0.20 2000.00 ! Range of permited half major axis
0.00 0.10 ! Range of permited eccentricity
exit
!