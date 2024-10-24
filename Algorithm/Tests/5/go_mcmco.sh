#! /bin/bash
export OMP_NUM_THREADS=8
export STACKSIZE=1000000
/Users/lacquema/Oracle.env/Oracle/Algorithm/bin/astrom_mcmco <<!
2 # New simulation LM and MCMC
1 0 0 0 # Type of data (RelAstro? AbsRV? RelRV? AbsAstro?)
1 # Number of orbits
data.mod
1d-7 # Precision
1 1 0 # Format data (1=DDMMYYYY/2=JD 1=(DEC,RA)/2=(SEP,PA) CorrCoeff?
0.00 pc # Distance
0.0 # First guess of center mass [ms]
adjustment.dat
dump.dat
10000000 # Dump frequency
0 # Number of masses prior
0 # Reference of time
0. mj 0. 0. 0. 0. 0. 0.  # First guess of orbit 1 parameters (m[mj] a[AU] e i[deg] Om[deg] om[deg] tp[MJD])
0.00 0.00 # Range of permited period
0.00 0.00 # Range of permited semi-major axis
0.00 0.00 # Range of permited eccentricity
!