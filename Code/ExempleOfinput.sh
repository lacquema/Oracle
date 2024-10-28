#! /bin/bash
export OMP_NUM_THREADS=8
export STACKSIZE=1000000
./astrom_mcmcop <<!
2        ! Action
1 1 0 0 ! Relative astrometry + radial velocities
2        ! Number of planets
hd196885_hci.dat
hd196885_rv.dat
1d-7     !    Precision
3        ! Data format hci; 3 = ipl JD Dec RA dDec dRA Corr
2        ! Data format rv; JD Rv dRV
1        ! Jitter
34.00475870 pc     ! distance (pc)
1.245     ! Stellar mass (Msun, first guess)
solhd196885
dumphd196885.dat
10000000
1   ! 1 mass priors
1. 0. 0.  ! Prior on central mass
1         ! 1 = Gaussian prior
1.245 0.045 ms ! Mean & standard deviation
50000.
-32.072180 1d-4 ! Initial V0 and jitter
2.227986 mj 2.472655 0.677069 130.661408 65.163924 -168.653531 49060.437997 ! m a e i Om om tp (Aa)
517.168534 mj 20.074605 0.432331 119.239667 79.513589 -125.118806 45776.684398 ! m a e i Om om tp (AB)
30. 300000.     ! Gamme de périodes
0.2 2000. ! Gamme de demi-grands axes
0.  0.99  ! Gamme d''excentricites
exit
!

