************************************************************************
*                                                                      *
      SUBROUTINE INTERP (XARR,YARR,NARR,XVAL,YVAL,ACCY)
*                                                                      *
*   This routine returns  YVAL  given a value  XVAL by interpolating   *
*   using a pair of arrays XARR(1:NARR), YARR(1:NARR), that tabulate   *
*   a  function.  ACCY  is the  desired  accuracy of the estimate: a   *
*   warning message is issued if this is not achieved.  A warning is   *
*   also issued when the routine is extrapolating.  Aitken's algori-   *
*   thm is used. See, for  instance, F B Hildebrand, Introduction to   *
*   Numerical Analysis,  2nd ed., McGraw-Hill, New York, NY, 1974.     *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
      LOGICAL SET,USED
*
      DIMENSION XARR(NARR),YARR(NARR)
*
*   MXORD is the maximum order of the interpolation
*
      PARAMETER (MXORD = 11)
      DIMENSION DX(MXORD),X(MXORD),
     :          EST(MXORD),
     :          USED(2*MXORD+2),
     :          POLY((MXORD*(MXORD+1))/2)
*
*   This function dispenses with the need for a two-dimensional
*   array for the Aitken lower triangle matrix
*
      ILOC (IND1,IND2) = (IND1*(IND1-1))/2+IND2
*
*   Determine the nearest two XARR entries bounding XVAL
*
      IF     (XVAL .LT. XARR(   1)) THEN
         NRSTLO =    1
         NRSTHI =    1
         WRITE (*,300)
      ELSEIF (XVAL .GT. XARR(NARR)) THEN
         NRSTLO = NARR
         NRSTHI = NARR
         WRITE (*,300)
      ELSE
         K = 0
    1    K = K+1
         IF (XARR(K) .LT. XVAL) THEN
            NRSTLO = K
            GOTO 1
         ELSE
            NRSTHI = K
         ENDIF
      ENDIF
*
*   Clear relevant piece of use-indicator array
*
      LLO = MAX (NRSTLO-MXORD,   1)
      LHI = MIN (NRSTHI+MXORD,NARR)
      LLR = LLO-1
      DO 2 K = LLO,LHI
         USED(K-LLR) = .FALSE.
    2 CONTINUE
*
*   Determine next-nearest XARR entry
*
      DO 5 IROW = 1,MXORD
         LLO = MAX (NRSTLO-IROW+1,   1)
         LHI = MIN (NRSTHI+IROW-1,NARR)
         SET = .FALSE.
         DO 3 K = LLO,LHI
            IF (.NOT. USED(K-LLR)) THEN
               IF (.NOT. SET) THEN
                  DIFF = XARR(K)-XVAL
                  LOCNXT = K
                  SET = .TRUE.
               ELSE
                  DIFFT = XARR(K)-XVAL
                  IF (ABS (DIFFT) .LT. ABS (DIFF)) THEN
                     DIFF = DIFFT
                     LOCNXT = K
                  ENDIF
               ENDIF
            ENDIF
    3    CONTINUE
         USED(LOCNXT-LLR) = .TRUE.
         X(IROW) = XARR(LOCNXT)
         DX(IROW) = DIFF
*
*   Fill table for this row
*
         DO 4 K = 1,IROW
            ILIROK = ILOC (IROW,K)
            IF (K .EQ. 1) THEN
               POLY(ILIROK) = YARR(LOCNXT)
            ELSE
               ILDIAG = ILOC (K-1,K-1)
               ILOTHR = ILOC (IROW,K-1)
               POLY(ILIROK) =  ( POLY(ILDIAG)*DX(IROW)
     :                          -POLY(ILOTHR)*DX(K -1) )
     :                       / ( X(IROW)-X(K-1) )
            ENDIF
    4    CONTINUE
*
*   Pick off the diagonal element
*
         ILDIAG = ILOC (IROW,IROW)
         EST(IROW) = POLY(ILDIAG)
*
    5 CONTINUE
*
*   Now the estimate vector is filled in, so obtain the
*   best estimate
*
      DEBEB = ABS ((EST(2)-EST(1))/EST(2))
      IBEST = 2
      DO 6 IROW = 3,MXORD
         DEBE = ABS ((EST(IROW)-EST(IROW-1))/EST(IROW))
         IF (DEBE .LT. DEBEB) THEN
            DEBEB = DEBE
            IBEST = IROW
         ENDIF
    6 CONTINUE
      YVAL = EST(IBEST)
*
      IF (DEBEB .GT. ACCY) THEN
         WRITE (*,301) DEBEB,ACCY
      ENDIF
*
      RETURN
*
  300 FORMAT ('INTERP: Extrapolating, not interpolating.')
  301 FORMAT ('INTERP: Accuracy of interpolation (',1P,1D10.3,') is',
     :        ' below input criterion (',1D10.3,').')
*
      END
