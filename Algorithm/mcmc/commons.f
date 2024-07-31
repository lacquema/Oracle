c
C...  Definition of constants and structures common to all versions of the code.
C...  To be included into the module 'DATA' in every version 

        INTEGER*4, PARAMETER :: NDATATYP = 4 ! # of data types
        INTEGER*4, PARAMETER :: NDATMAX = 100
        INTEGER*4, PARAMETER :: NCH = 8 ! Nombre de chaines
        INTEGER*4, PARAMETER :: NINFER = 50000
        INTEGER*4, PARAMETER :: NMAX = NINFER*NCH
        INTEGER*8, PARAMETER :: NDIAG = 500000_8
        INTEGER*8, PARAMETER :: SD0 = 6
        INTEGER*8, PARAMETER :: SDL = 35
        INTEGER*4, PARAMETER :: SD = SDL
        INTEGER*8, PARAMETER :: VERBFREQ = 100000_8
        INTEGER*8, PARAMETER :: WRITFREQ = 10000000_8
        REAL*8, PARAMETER :: BETAMAX = 1.d10
        REAL*8, PARAMETER :: PI = 3.14159265358979323846d0     
        REAL*8, parameter :: AU = 1.495978707d8  ! (km)      
        REAL*8, PARAMETER :: DAY = 86400.d0
        REAL*8, PARAMETER :: YJD = 365.25636042d0 ! Sideral year
        REAL*8, PARAMETER :: YEAR = YJD*DAY       
        REAL*8, PARAMETER :: DR = 180.d0/PI
        REAL*8, PARAMETER :: DPI = 2.d0*PI
        REAL*8, PARAMETER :: LDPI = LOG(DPI)
        REAL*8, PARAMETER :: RCONV = 1.01d0
        REAL*8, PARAMETER :: KG = 0.01720209895d0 ! Gauss' constant
        REAL*8, PARAMETER :: KG2 = KG*KG   ! KG ~ 2*Pi/YJD
        REAL*8, PARAMETER :: SMAS = KG2 ! (mass unit)^(-1) = Solar mass
c...      Time Unit = day, length unit = ua, mass unit = 1/KG^2 Msun => G=1
        REAL*8, PARAMETER :: MPS = DAY/AU
c...  (Velocity unit)^(-1) = 1 km/s expressed in AU per DAY
        REAL*8, PARAMETER :: MJUP = SMAS/1047.355d0
        REAL*8, PARAMETER :: MTINY = 1d-6*MJUP
        REAL*8, PARAMETER :: MHUGE = 100.d0*SMAS
        REAL*8, PARAMETER :: TSCONV = 1d3
        REAL*8, PARAMETER :: RINF = 1.1d0
        REAL*8, PARAMETER :: TSINF = 1d2
        REAL*8, PARAMETER :: MINPROB = 1.d-1
        REAL*8, PARAMETER :: SIGV0 = 0.1d0*MPS ! Jitter norm = 0.1 km/s
        CHARACTER*80, PARAMETER :: LOGFILE = 'logfile.dat'
        INTEGER*4 :: FREQDUMP     ! Saving frequency
        REAL*8 :: EPS
        CHARACTER*1 :: CRR
        INTEGER*4, DIMENSION(NDATATYP) :: DATATYP ! Type data line
        LOGICAL, DIMENSION(NDATATYP) :: ISDATA ! types of data present
        LOGICAL :: RADVEL
        INTEGER*4 :: JITNUM     ! 1 if Jitter to be considered
        INTEGER*4, PARAMETER :: NFIL = 3 ! # of files
        INTEGER*4, PARAMETER :: NLIM = 3
        INTEGER*4, PARAMETER :: NEL = 7 ! # of elements to fit per planet
        
        TYPE CHAIN
          REAL*8, DIMENSION(:), ALLOCATABLE ::
     &       MOY,                ! Means of chain
     &       BETA,               ! Scaling vector
     &       P                   ! Parameters of last model            
          REAL*8, DIMENSION(:,:), ALLOCATABLE ::
     &       COV                 ! Covariance matrix
          REAL*8 :: CHI2         ! Chi2 of last model
        END TYPE CHAIN               

        TYPE PLANET
          INTEGER*4 :: NDATAS,NDATVR ! # of data of each kind 
          REAL*8, DIMENSION(:), ALLOCATABLE :: X,Y,TAS ! Astrometric data
          REAL*8, DIMENSION(:), ALLOCATABLE :: SIGX,SIGY,SIGXM2,SIGYM2,
     &                     SIGXM1,SIGYM1,RHOXY,UFAC ! Errors in ast. data
          REAL*8, DIMENSION(:), ALLOCATABLE ::
     &                      TVR,V,SIGV,SIGV2,SIGVM2 ! Pl. RV data
          REAL*8 ::  MU,DMU,        ! mass = m_i
     &               MDYN,DMDYN,    ! Dynamical mass = sum(m_j,j=0..i)
     &               A,DA,            ! Semi-major axis                    
     &               PER,DPER,  ! Period
     &               EXC,EXQ,DEXC, ! Eccentricity
     &               INC,DINC,TI,   ! Inclination & tan(i/2)
     &               O,DOO,CO,SO,   ! Long. of ascending node
     &               OM,DOM,COM,SOM,! Arg. of periastron
     &               W,DW,CW,SW,    ! W = OM+O
     &               PHI,CP,SP,DPHI,! PHI = OM-O
     &               TP,DTP         ! Time of periastron      
        END TYPE PLANET

        TYPE STARDATA
          INTEGER*4 :: NDATAS,NDATVR ! # of data of each kind 
          REAL*8, DIMENSION(:), ALLOCATABLE :: X,Y,TAS ! Abs. Astrometric data
          REAL*8, DIMENSION(:), ALLOCATABLE :: SIGX,SIGY,SIGXM2,SIGYM2,
     &                     SIGXM1,SIGYM1,RHOXY,UFAC ! Errors in ast. data
          REAL*8, DIMENSION(:), ALLOCATABLE ::
     &                      TVR,V,SIGV,SIGV2,SIGVM2   ! Stel. RV data
          REAL*8 ::  MASS,DMASS,    ! Mass & error
     &               V0,DV0,        ! Velocity offset & error
     &               SIGJV,DSIGJV,  ! Vec. Jitter & error
     &               JDOFFSET       ! Julian day offset
          REAL*8 :: DIST,PARX       ! Distance, parallax
          REAL*8 :: T0              ! Refrence time for mean logitude
        END TYPE STARDATA
        
        TYPE PRIOR             ! Priors for masses 
          INTEGER*4 :: TYP     ! Type of prior : 0 = logarithmic (default)
c    1 = Gaussian, 2 = Gauss-Normal, 3 = linear, 4 = fixed, 5 = sin        
c    6 = haeviside
          REAL*8, DIMENSION(2) :: BOUND ! Bounds
          REAL*8 :: MEAN,SDEV           ! Mean & Standard deviation
          REAL*8, DIMENSION(:), ALLOCATABLE :: ACOF,BCOF
c           Prior on sum(acof(i)*m(i))=sum(bcof(i)*mdyn(i))  
        END TYPE PRIOR
                
        TYPE(PLANET), DIMENSION(:), ALLOCATABLE :: PLA
        TYPE(STARDATA) :: STAR 
        REAL*8, DIMENSION(:), ALLOCATABLE :: PSTART
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: PSAV
        TYPE(PRIOR), DIMENSION(:), ALLOCATABLE :: MPRIOR
        TYPE(PRIOR) :: VPRIOR(NLIM)
        INTEGER*4 :: NPLA       ! # of planets        
        INTEGER*4 :: NPAR      ! # of parameters to fit
        INTEGER*4 :: NFREE     ! # of degrees of freedom
        INTEGER*4 :: NPRIOR    ! # of Prior prescriptions for masses
        LOGICAL :: MULTIPLA    ! True if more than 1 planet 
        INTEGER*4 :: DATEINPUT ! 1 = JJMMYYYY, 2 = JD-offset
        INTEGER*4 :: XYFORMAT  ! 1 = (DEC,RA), 2 = (SEP,PA)
        INTEGER*4 :: CORRNUM   ! 1 = Correlation XY, 0 = no
        LOGICAL :: CORR        ! True if correlation
        
