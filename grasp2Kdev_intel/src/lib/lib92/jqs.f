************************************************************************
*                                                                      *
      FUNCTION JQS (IWHICH,ISUBSH,ICSF)
*                                                                      *
*   JQS is a subshell quantum number for subshell ISUBSH in configu-   *
*   ration state function  ICSF:  the seniority if IWHICH is 1;  the   *
*   quantum number w if IWHICH is 2, and 2J+1 if IWHICH is 3.          *
*                                                                      *
*   Call(s) to: [LIB92]: IUNPCK.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
*                                                                      *
************************************************************************
*
! Packed array declared I*4 and a common added (iounit)
! XHH 1997.02.12
Cww      INTEGER PNJCUP,PNTRIQ
      POINTER (PNJCUP,JCUPDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      include 'parameters.def'
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
      INTEGER*4 JQSA
      POINTER (PNTJQS,JQSA(NNNWP,3,*))
*
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /STAT/PNTJQS,PNJCUP
      COMMON/iounit/istdi,istdo,istde
*
c     IF ((IWHICH .GE. 1) .AND. (IWHICH .LE. 3)) THEN
c        IF ((ISUBSH .GE. 1) .AND. (ISUBSH .LE. NW)) THEN
c           IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCF)) THEN
*cff           JQS = IUNPCK (JQSA(1,IWHICH,ICSF),ISUBSH)
               JQS = IBITS( JQSA((ISUBSH-1)/4+1,IWHICH,ICSF),
     :                      8*MOD(ISUBSH-1,4), 8)
c           ELSE
c              WRITE(istde,*) 'JQS: Argument ICSF is out of range.'
c              STOP
c           ENDIF
c        ELSE
c           WRITE(istde,*) 'JQS: Argument ISUBSH is out of range.'
c           STOP
c        ENDIF
c     ELSE
c        WRITE(istde,*) 'JQS: Argument IWHICH is out of range.'
c        STOP
c     ENDIF
*
      RETURN
      END
