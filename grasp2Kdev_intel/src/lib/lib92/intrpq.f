************************************************************************
*                                                                      *
      SUBROUTINE INTRPQ (PA,QA,MA,RA,J,DNORM)
*                                                                      *
*   This  subprogram  interpolates  the  arrays  PA(1:MA), QA(1:MA),   *
*   tabulated on grid RA(1:MA) into the COMMON arrays PF(1:MF(J),J),   *
*   QF(1:MF(J),J). (Aitken's  algorithm is used. See F B Hildebrand,   *
*   Introduction  to  Numerical  Analysis, 2nd ed., McGraw-Hill, New   *
*   York, NY, 1974.) The orbital is renormalized.                      *
*                                                                      *
*   SUBROUTINEs called: RINT.                                          *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 14 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      LOGICAL LDBPR,SET,USED
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      DIMENSION PA(*),QA(*),RA(*)
*     DIMENSION PA(1),QA(1),RA(1)
      DIMENSION USED(NNNP)
*
      COMMON/DEBUGR/LDBPR(30)
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
*   MXORD is the maximum order of the interpolation
*
      PARAMETER (MXORD = 13)
      DIMENSION X(MXORD),DX(MXORD),POLYP((MXORD*(MXORD+1))/2)
     :                            ,POLYQ((MXORD*(MXORD+1))/2)
*
*   This function dispenses with the need for a two-dimensional
*   array
*
      ILOC (IND1,IND2) = (IND1*(IND1-1))/2+IND2
*
*   Initialization
*
      RAMA = RA(MA)
      RN = R(N)
*
*   This is always true in GRASP
*
      PF(1,J) = 0.0D 00
      QF(1,J) = 0.0D 00
*
*   Checks
*
      IF (RAMA .GT. RN) THEN
         WRITE (*,300) RN,RAMA
         PRINT*, "N =",N,"MA =",MA
         STOP
      ENDIF
*
*   Determine end of grid
*
      I = N
    1 I = I-1
      IF (R(I) .LE. RAMA) THEN
         MFJ = I
      ELSE
         GOTO 1
      ENDIF
      MF(J) = MFJ
*
*   Overall initialization for interpolation
*
      NRSTLO = 0
      KOUNT = 0
*
*   Perform interpolation
*
      DO 7 I = 2,MFJ
*
*   Initialization for interpolation
*
         XBAR = R(I)
         IROW = 0
         PESTL = 0.0D 00
         QESTL = 0.0D 00
*
*   Determine the nearest two grid points bounding the present
*   grid point
*
    2    K = NRSTLO+1
         IF (RA(K) .LT. XBAR) THEN
            NRSTLO = K
            GOTO 2
         ELSE
            NRSTHI = K
         ENDIF
*
*   Clear relevant piece of use-indicator array
*
         LLO = MAX (NRSTLO-MXORD, 1)
         LHI = MIN (NRSTHI+MXORD,MA)
         DO 3 K = LLO,LHI
            USED(K) = .FALSE.
    3    CONTINUE
*
*   Determine next-nearest grid point
*
    4    IROW = IROW+1
         LLO = MAX (NRSTLO-IROW+1, 1)
         LHI = MIN (NRSTHI+IROW-1,MA)
         SET = .FALSE.
         DO 5 K = LLO,LHI
            IF (.NOT. USED(K)) THEN
               IF (.NOT. SET) THEN
                  DIFF = RA(K)-XBAR
                  LOCNXT = K
                  SET = .TRUE.
               ELSE
                  DIFFT = RA(K)-XBAR
                  IF (ABS (DIFFT) .LT. ABS (DIFF)) THEN
                     DIFF = DIFFT
                     LOCNXT = K
                  ENDIF
               ENDIF
            ENDIF
    5    CONTINUE
         USED(LOCNXT) = .TRUE.
         X(IROW) = RA(LOCNXT)
         DX(IROW) = DIFF
*
*   Fill table for this row
*
         DO 6 K = 1,IROW
            ILIROK = ILOC (IROW,K)
            IF (K .EQ. 1) THEN
               POLYP(ILIROK) = PA(LOCNXT)
               POLYQ(ILIROK) = QA(LOCNXT)
            ELSE
               ILDIAG = ILOC (K-1,K-1)
               ILOTHR = ILOC (IROW,K-1)
               DXKMN1 = DX(K-1)
               DXIROW = DX(IROW)
               FACTOR = 1.0D 00 / ( X(IROW)-X(K-1) )
               POLYP(ILIROK) =  ( POLYP(ILDIAG)*DXIROW
     :                           -POLYP(ILOTHR)*DXKMN1 )
     :                         * FACTOR
               POLYQ(ILIROK) =  ( POLYQ(ILDIAG)*DXIROW
     :                           -POLYQ(ILOTHR)*DXKMN1 )
     :                         * FACTOR
            ENDIF
    6    CONTINUE
*
*   Check for convergence
*
         ILDIAG = ILOC (IROW,IROW)
         PESTT = POLYP(ILDIAG)
         QESTT = POLYQ(ILDIAG)
         IF ((PESTT .EQ. 0.0D 00) .OR.
     :       (QESTT .EQ. 0.0D 00)) THEN
            IF (IROW .LT. MXORD) THEN
               GOTO 4
            ELSE
               PF(I,J) = PESTT
               QF(I,J) = QESTT
            ENDIF
         ELSE
            DPBP = ABS ((PESTT-PESTL)/PESTT)
            DQBQ = ABS ((QESTT-QESTL)/QESTT)
            IF ((DQBQ .LT. ACCY) .AND. (DPBP .LT. ACCY)) THEN
               PF(I,J) = PESTT
               QF(I,J) = QESTT
            ELSE
               PESTL = PESTT
               QESTL = QESTT
               IF (IROW .LT. MXORD) THEN
                  GOTO 4
               ELSE
                  PF(I,J) = PESTT
                  QF(I,J) = QESTT
                  KOUNT = KOUNT+1
               ENDIF
            ENDIF
         ENDIF
*
    7 CONTINUE
*
*   Ensure that all points of the array are defined by setting the
*   tail to zero
*
      MFJP1 = MFJ+1
      DO 8 I = MFJP1,N
         PF(I,J) = 0.0D 00
         QF(I,J) = 0.0D 00
    8 CONTINUE
*
      IF (LDBPR(3) .AND. (KOUNT .GT. 0)) WRITE (99,301) ACCY,KOUNT,MFJ
*
*   Normalization
*
      DNORM = RINT (J,J,0)
      DNFAC = 1.0D 00/SQRT (DNORM)
      DO 9 I = 1,MFJ
         PF(I,J) = PF(I,J)*DNFAC
         QF(I,J) = QF(I,J)*DNFAC
    9 CONTINUE
*
      RETURN
*
  300 FORMAT (/'INTRPQ: Grid of insufficient extent:'
     :        /' Present grid has R(N) = ',1P,1D19.12,' Bohr radii'
     :        /'          Require R(N) = ',   1D19.12,' Bohr radii')
  301 FORMAT (/'INTRPQ: Interpolation procedure not converged to',
     :         1P,1D19.12,' for ',1I3,' of ',1I3,' tabulation points')
*
      END
