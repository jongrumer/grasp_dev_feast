************************************************************************
*                                                                      *
      SUBROUTINE LODRWFI(NAME,NTESTG)
*                                                                      *
*   This subroutine loads  radial wavefunctions from the  .rwf  file   *
*   and performs some related setup.                                   *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, INTRPQ, ORTHSC.                *
*                                                                      *
*   Written by Per Jonsson                              June 1996      *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      LOGICAL LDBPR

      PARAMETER (NLMAX = 40)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590) 
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)

      CHARACTER*2 NHII, NH(NNNW)
      CHARACTER*6 G92RWF
      CHARACTER*128 NAME

      DIMENSION NAK(NNNW),NP(NNNW),E(NNNW),GAMA(NNNW)
*
      POINTER (PNTRPA,PA(*)),(PNTRQA,QA(*)),(PNTRRA,RA(*))
*
      POINTER (PNTRPFII,PFII(NNNP,*)),(PNTRQFII,QFII(NNNP,*))
*
      COMMON/DEBUGR/LDBPR(30)
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
*
      COMMON/SBDAT/NAKINVII(NNNW),NSHLII(NLMAX),NSHLPII(NLMAX,NLMAX),   
     :             NAKINVFF(NNNW),NSHLFF(NLMAX),NSHLPFF(NLMAX,NLMAX),   
     :             NSHLPPII(NLMAX,NNNW),NSHLPPFF(NLMAX,NNNW),
     :             NINII(NLMAX),NINFF(NLMAX),IKAPPA(NLMAX),KAMAX
*
*   Common relevant for the final state
*
      COMMON/ORB1II/EII(NNNW),GAMAII(NNNW)
     :      /ORB2II/NCFII,NWII
     :      /ORB4II/NPII(NNNW),NAKII(NNNW)
     :      /ORB10II/NHII(NNNW)
     :      /WAVEII/PZII(NNNW),PNTRPFII,PNTRQFII,MFII(NNNW)
*
      NTESTL = 00
      NTEST = MAX0(NTESTL,NTESTG)
      NTEST = 0

*
*   Write entry message
*
      PRINT *, 'Loading Radial WaveFunction File for initial state...'
*
*   Open the radial wave function file
*
      J = INDEX(NAME,' ')
      OPEN (UNIT = 21,FILE=NAME(1:J-1)//'.w',FORM='UNFORMATTED',
     :     STATUS='OLD')
*
*   Save NAK, NP and NH
*
       
      DO K = 1,NWII
        NAK(K) = NAKII(K)
        NP(K)  = NPII(K)
        NH(K)  = NHII(K)
      ENDDO
*
*   Allocate storage to orbital arrays
*
      CALL ALLOC (PNTRPFII,NNNP*NWII,8)
      CALL ALLOC (PNTRQFII,NNNP*NWII,8)
*
      CON = Z/C
      CON = CON*CON
*
      DO J = 1,NWII
        DO I = 1,NNNP
          PFII(I,J) = 0.0D 00
          QFII(I,J) = 0.0D 00
        ENDDO 
*
        K = ABS (NAK(J))
        IF (NPARM .EQ. 0) THEN
          FKK = DBLE (K*K)
          IF (FKK .GE. CON) THEN
            GAMA(J) = SQRT (FKK-CON)
          ELSE
            PRINT *, 'LODRWF: Imaginary gamma parameter'
            PRINT *, ' for ',NP(J),NH(J),' orbital; the'
            PRINT *, ' point model for the nucleus'
            PRINT *, ' is inappropriate for Z > ',C,'.'
            STOP
          ENDIF
        ENDIF
*
      ENDDO
*
*   Read orbital information from Read Orbitals File;
*
      NWIN = 0
      READ (21,IOSTAT = IOS) G92RWF
      IF ((IOS .NE. 0) .OR. (G92RWF .NE. 'G92RWF')) THEN
        PRINT *, 'This is not a Radial WaveFunction File;'
        CLOSE (21)
      ENDIF

      IF (NTEST.GE.100) THEN
        WRITE(*,*) '******************'
        WRITE(*,*) ' Entering lodrwfi'
        WRITE(*,*) '******************'
      ENDIF

    3 READ (21,IOSTAT = IOS) NPY,NAKY,EY,MY
      IF (IOS .EQ. 0) THEN
        CALL ALLOC (PNTRPA,MY,8)
        CALL ALLOC (PNTRQA,MY,8)
        CALL ALLOC (PNTRRA,MY,8)
        READ (21) PZY,(PA(I),I = 1,MY),(QA(I),I = 1,MY)
        READ (21) (RA(I),I = 1,MY)
*
*    Orbital order as defined in kapdata
*
        JJ = 0
        DO K = 1,KAMAX
          IF (K.GT.1) JJ = NSHLII(K-1) + JJ
          DO J = 1,NSHLII(K)
            KK = NSHLPII(K,J)
            IF ((NPY .EQ. NP(KK)) .AND.
     :        (NAKY .EQ. NAK(KK)) ) THEN
              JJJ = JJ + J
              PZII(JJJ) = PZY
              EII(JJJ) = EY
              NAKII(JJJ) = NAK(KK)
              NPII(JJJ) = NP(KK)
              NHII(JJJ) = NH(KK)
              GAMAII(JJJ) = GAMA(KK)
              CALL INTRPQI (PA,QA,MY,RA,JJJ,DNORM)
              IF (NTEST.GE.100) THEN
                WRITE (*,301) NPII(JJJ),NHII(JJJ),EII(JJJ),DNORM
              ENDIF
              IF (NTEST.GT.1000) THEN
                WRITE(*,*) 'PF              QF             RA'
                DO KKK = 1,MFII(JJJ)
                  WRITE(*,*) PFII(KKK,JJJ),QFII(KKK,JJJ),RA(KKK)
                ENDDO
              ENDIF 
              NWIN = NWIN+1
            ENDIF
          ENDDO
        ENDDO
        CALL DALLOC (PNTRPA)
        CALL DALLOC (PNTRQA)
        CALL DALLOC (PNTRRA)
        GOTO 3
      ENDIF
      IF (LDBPR(3)) WRITE (99,*) ' orbitals renormalised;'
*
*   Stop with an error message if all orbitals are not known
*
      IF (NWIN .LT. NWII) THEN
         PRINT *, 'LODRWF: All required orbitals not'
         PRINT *, ' found.'
         IERR = 1
         GOTO 5
      ENDIF
*
      PRINT *, ' ... load complete;'
*
    5 CLOSE(21)
      IF (NTEST.GE.100) THEN
        WRITE(*,*) 'Sorted order should be the same as from kapdat'
        DO J = 1,NWII 
          WRITE (*,301) NPII(J),NHII(J),EII(J),DNORM
        ENDDO
        WRITE(*,*)
        WRITE(*,*) '*****************'
        WRITE(*,*) ' Leaving lodrwfi'
        WRITE(*,*) '*****************'
      ENDIF

      RETURN
*
  301 FORMAT (2X,I2,A2,4X,1P,1D22.15,4X,1D22.15)
*
      END
