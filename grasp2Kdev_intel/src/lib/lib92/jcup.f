************************************************************************
*                                                                      *
      FUNCTION JCUP (LOC,ICSF)
*                                                                      *
*   JCUP is the 2J+1 value of the LOCth nontrivial intermediate ang-   *
*   ular momentum in CSF  ICSF.                                        *
*                                                                      *
*   Call(s) to: [LIB92]: IUNPCK.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
*                                                                      *
************************************************************************
*
! Packed array declared I*4 and a common (iounit) added
! XHH 1997.02.12
Cww      INTEGER PNTJQS,PNTRIQ
      POINTER (PNTJQS,JQSDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      include 'parameters.def'
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
      INTEGER*4 JCUPA
      POINTER (PNJCUP,JCUPA(NNNWP,*))
*
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /STAT/PNTJQS,PNJCUP
      COMMON/iounit/istdi,istdo,istde
*
      IF ((LOC .GE. 1) .AND. (LOC .LE. NW-1)) THEN
         IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCF)) THEN
*cff        JCUP = IUNPCK (JCUPA(1,ICSF),LOC)
            JCUP = IBITS (JCUPA((LOC-1)/4+1,ICSF),8*MOD(LOC-1,4), 8)
         ELSE
            WRITE(istde,*) 'JCUP: Argument ICSF is out of range.'
            STOP
         ENDIF
      ELSE
         WRITE(istde,*) 'JCUP: Argument LOC is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
