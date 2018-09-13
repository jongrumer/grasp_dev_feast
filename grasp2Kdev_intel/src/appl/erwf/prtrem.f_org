************************************************************************
*                                                                      *
      SUBROUTINE PRTREM (ALL)
*                                                                      *
*   Prints a list of subshells that remain to be estimated.  ALL  is   *
*   .TRUE.  if all subshells have been estimated and  .FALSE. other-   *
*   wise.                                                              *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL ALL,SET
      CHARACTER*80 RECORD
      CHARACTER*2 CNUM,NH
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
      COMMON/LEFT/SET(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
      COMMON/iounit/istdi,istdo,istde
*
*   Determine if there are any subshell radial wavefunctions that
*   remain to be estimated
*
      DO 1 I = 1,NW
         IF (.NOT. SET(I)) THEN
            ALL = .FALSE.
            GOTO 2
         ENDIF
    1 CONTINUE
      ALL = .TRUE.
      GOTO 4
*
*   Print a list of subshell radial wavefunctions that remain to
*   be estimated; this list is no more than 80 characters wide
*
    2 WRITE(istde,*) 'The following subshell radial wavefunctions '
     &             , 'remain to be estimated:'
      IEND = 0
      DO 3 I = 1,NW
         IF (.NOT. SET(I)) THEN
            IF (IEND .GT. 75) THEN
               WRITE(istde,*) RECORD(1:IEND)
               IEND = 0
            ENDIF
            IF (IEND .GT. 0) THEN
               IBEG = IEND+1
               IEND = IBEG
               RECORD(IBEG:IEND) = ' '
            ENDIF
            IBEG = IEND+1
            CALL CONVRT (NP(I),CNUM,LENTH)
            IEND = IBEG+LENTH-1
            RECORD(IBEG:IEND) = CNUM(1:LENTH)
            IF (NAK(I) .LT. 0) THEN
               LENTH = 1
            ELSE
               LENTH = 2
            ENDIF
            IBEG = IEND+1
            IEND = IBEG+LENTH-1
            RECORD(IBEG:IEND) = NH(I)(1:LENTH)
         ENDIF
    3 CONTINUE
      IF (IEND .GT. 1) WRITE(istde,*) RECORD(1:IEND)
*
    4 RETURN
      END
