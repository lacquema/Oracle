c-------------------------------------------------------------------------
c     Programme de test des priors sur les masses : evalue directement
c     USER_PRIOR_MAP (MAP), USER_PRIOR_RATIO (MCMC) et FMAP_PRIOR (L-M)
c     pour un prior et une liste de masses lus sur l'entree standard :
c        type, moyenne, ecart-type, borne 1, borne 2
c        nombre de masses, puis les masses (unites internes)
c     Sorties (lignes "TEST") :
c        MAP  i  terme ajoute au MAP pour la masse i
c        EX   i  masse apres l'appel de USER_PRIOR_MAP (doit etre inchangee)
c        LM   i  F de FMAP_PRIOR (-2 ln prior)
c        RAT  i j  rapport du MCMC pour passer de la masse i a la masse j
c-------------------------------------------------------------------------
        MODULE DATA

        PUBLIC

        INCLUDE 'commons.f'

        END MODULE DATA

c-------------------------------------------------------------------------

        PROGRAM PRIORS

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    I,J,NM          ! Indexes, number of masses
        REAL*8, DIMENSION(:), ALLOCATABLE :: M ! Masses
        REAL*8 ::       FMAP,           ! MAP term
     &                  EX,EXT,         ! Masses given to the routines
     &                  RAPQ,           ! Prior ratio
     &                  F,DF,D2F        ! Output of FMAP_PRIOR

        NPLA = 1
        NPRIOR = 0
        ALLOCATE(MPRIOR(0:0))
        ALLOCATE(MPRIOR(0)%ACOF(0:NPLA))
        ALLOCATE(MPRIOR(0)%BCOF(0:NPLA))
        MPRIOR(0)%ACOF = 0.d0
        MPRIOR(0)%BCOF = 0.d0
        MPRIOR(0)%MUNIT = 1.d0
        READ(5,*)MPRIOR(0)%TYP,MPRIOR(0)%MEAN,MPRIOR(0)%SDEV,
     &           MPRIOR(0)%BOUND(1:2)
        READ(5,*)NM
        ALLOCATE(M(NM))
        READ(5,*)M(1:NM)

 1      FORMAT('TEST ',A,1X,I4,1X,ES24.16)
 2      FORMAT('TEST ',A,2(1X,I4),1X,ES24.16)
        DO I = 1,NM
           FMAP = 0.d0
           EX = M(I)
           CALL USER_PRIOR_MAP(0,FMAP,EX)
           WRITE(SD,1)'MAP',I,FMAP
           WRITE(SD,1)'EX',I,EX
           EX = M(I)
           CALL FMAP_PRIOR(0,EX,F,DF,D2F)
           WRITE(SD,1)'LM',I,F
        END DO
        DO I = 1,NM
           DO J = 1,NM
              EX = M(I)
              EXT = M(J)
              CALL USER_PRIOR_RATIO(0,RAPQ,EX,EXT)
              WRITE(SD,2)'RAT',I,J,RAPQ
           END DO
        END DO

        END
