#! /bin/bash
export OMP_NUM_THREADS=8
export STACKSIZE=1000000
/Users/lacquema/Oracle.env/Oracle/Interface/NewSimu/../../Algorithm/bin/astrom_mcmco <<!
2 # New simulation LM and MCMC
1 0 0 0 # Type of data (RelAstro? AbsRV? RelRV? AbsAstro?)
1 # Number of orbits
ggtau_Ab12_hci.txt
1d-7 # Precision
1 1 1 # Format data (1=DDMMYYYY/2=JD 1=(DEC,RA)/2=(SEP,PA) CorrCoeff?
145.00 pc # Distance
1.431882 # First guess of center mass [ms]
adjustment.dat
dump.dat
10000000 # Dump frequency
0 # Number of masses prior
56583 # Reference of time
400 mj 4 0.5 40 110 290 57754  # First guess of orbit 1 parameters (m[mj] a[AU] e i[deg] Om[deg] om[deg] tp[MJD])
0.00 10000.00 # Range of permited period
0.00 10.00 # Range of permited semi-major axis
0.00 0.90 # Range of permited eccentricity
!