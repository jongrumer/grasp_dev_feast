************************************************************************
*                                                                      *
      SUBROUTINE LODRWFF(NAME,NTESTG)
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

      PARAMETER (NLMAX = 20)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)

      CHARACTER*2 NHFF, NH(NNNW)
      CHARACTER*6 G92RWF
      CHARACTER*24 NAME

      DIMENSION NAK(NNNW),NP(NNNW),E(NNNW),GAMA(NNNW)
*
      POINTER (PNTRPA,PA(*)),(PNTRQA,QA(*)),(PNTRRA,RA(*))
*
      POINTER (PNTRPFFF,PFFF(NNNP,*)),(PNTRQFFF,QFFF(NNNP,*))
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
      COMMON/ORB1FF/EFF(NNNW),GAMAFF(NNNW)
     :      /ORB2FF/NCFFF,NWFF
     :      /ORB4FF/NPFF(NNNW),NAKFF(NNNW)
     :      /ORB10FF/NHFF(NNNW)
     :      /WAVEFF/PZFF(NNNW),PNTRPFFF,PNTRQFFF,MFFF(NNNW)
*
      NTESTL = 00
      NTEST = MAX0(NTESTL,NTESTG)
      NTEST = 0

*
*   Write entry message
*
      PRINT *, 'Loading Radial WaveFunction File for final state...'
*
*   Open the radial wave function file
*
      J = INDEX(NAME,' ')
      OPEN (UNIT = 21,FILE=NAME(1:J-1)//'.w',FORM='UNFORMATTED',
     :     STATUS='OLD')
*
*   Save NAK, NP and NH
*
      DO K = 1,NWFF
        NAK(K) = NAKFF(K)
        NP(K)  = NPFF(K)
        NH(K)  = NHFF(K)
      ENDDO
*
*   Allocate storage to orbital arrays
*
      CALL ALLOC (PNTRPFFF,NNNP*NWFF,8)
      CALL ALLOC (PNTRQFFF,NNNP*NWFF,8)
*
      CON = Z/C
      CON = CON*CON
*
      DO J = 1,NWFF
        DO I = 1,NNNP
          PFFF(I,J) = 0.0D 00
          QFFF(I,J) = 0.0D 00
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
        WRITE(*,*) ' Entering lodrwff'
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
          IF (K.GT.1) JJ = NSHLFF(K-1) + JJ
          DO J = 1,NSHLFF(K)
            KK = NSHLPFF(K,J)
            IF ((NPY .EQ. NP(KK)) .AND.
     :        (NAKY .EQ. NAK(KK)) ) THEN
              JJJ = JJ + J
              PZFF(JJJ) = PZY
              EFF(JJJ) = EY
              NAKFF(JJJ) = NAK(KK)
              NPFF(JJJ) = NP(KK)
              NHFF(JJJ) = NH(KK)
              GAMAFF(JJJ) = GAMA(KK)
              CALL INTRPQF (PA,QA,MY,RA,JJJ,DNORM)
              IF (NTEST.GE.100) THEN
                WRITE (*,301) NPFF(JJJ),NHFF(JJJ),EFF(JJJ),DNORM
              ENDIF
              IF (NTEST.GT.1000) THEN
                WRITE(*,*) 'PF              QF             RA'
                DO KKK = 1,MFFF(JJJ)
                  WRITE(*,*) PFFF(KKK,JJJ),QFFF(KKK,JJJ),RA(KKK)
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
      IF (NWIN .LT. NWFF) THEN
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
        DO J = 1,NWFF 
          WRITE (*,301) NPFF(J),NHFF(J),EFF(J),DNORM
        ENDDO
        WRITE(*,*)
        WRITE(*,*) '*****************'
        WRITE(*,*) ' Leaving lodrwff'
        WRITE(*,*) '*****************'
      ENDIF

      RETURN
*
  301 FORMAT (2X,I2,A2,4X,1P,1D22.15,4X,1D22.15)
*
      END
