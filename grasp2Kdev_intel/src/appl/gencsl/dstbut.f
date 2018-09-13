************************************************************************
*                                                                      *
      SUBROUTINE DSTBUT (IOCCS,NCORE,NORB,NLAST)
*                                                                      *
*   This  subroutine distributes the electrons in the peel subshells   *
*   without regard to occupation number.                               *
*                                                                      *
*   Written by Farid A Parpia and Wasantha P Wijesundera, at Oxford    *
*                                                                      *
*                                            Last update 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
      DIMENSION IOCCS(NNNW)
*
      IF (NLAST .NE. NORB) THEN
         IOCCS(NLAST) = IOCCS(NLAST)-1
         NLAST = NLAST+1
         IOCCS(NLAST) = IOCCS(NLAST)+1
      ELSE
         NOCC = IOCCS(NLAST)
         IOCCS(NLAST) = 0
         DO 1 ILOC = NORB-1,NCORE+1,-1
            IF (IOCCS(ILOC) .NE. 0) THEN
               NLAST = ILOC
               GOTO 2
            ENDIF
    1    CONTINUE
    2    IOCCS(NLAST) = IOCCS(NLAST)-1
         NLAST = NLAST+1
         IOCCS(NLAST) = IOCCS(NLAST)+1+NOCC
      ENDIF
*
      RETURN
      END
