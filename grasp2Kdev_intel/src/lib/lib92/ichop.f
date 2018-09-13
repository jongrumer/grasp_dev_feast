************************************************************************
*                                                                      *
      FUNCTION ICHOP (ISUBSH,ICSF)
*                                                                      *
*   ICHOP is -1 if subshell ISUBSH is empty in CSF  ICSF,  +1 if the   *
*   subshell is full, and 0 if it is open.                             *
*                                                                      *
*   Call(s) to: [LIB92]: IQ.                                           *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 30 Oct 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
*
* cff  Since ICHOP is always called from within a do-loop over the
*      appropriate range, testing seems redundant
*     IF ((ISUBSH .GE. 1) .AND. (ISUBSH .LE. NW)) THEN
*        IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCF)) THEN
            IOCC = IQ (ISUBSH,ICSF)
            IFULL = NKJ(ISUBSH)+1
            IF (IOCC .EQ. 0) THEN
               ICHOP = -1
            ELSEIF (IOCC .EQ. IFULL) THEN
               ICHOP = 1
            ELSE
               ICHOP = 0
            ENDIF
*        ELSE
*           PRINT *, 'ICHOP: Argument ICSF is out of range.'
*           STOP
*        ENDIF
*     ELSE
*        PRINT *, 'ICHOP: Argument ISUBSH is out of range.'
*        STOP
*     ENDIF
*
      RETURN
      END
