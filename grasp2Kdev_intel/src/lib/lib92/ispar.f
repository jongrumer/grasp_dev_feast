************************************************************************
*                                                                      *
      FUNCTION ISPAR (ICSF)
*                                                                      *
*   ISPAR is the value of P for CSF number ICSF.                       *
*                                                                      *
*   Call(s) to: [LIB92]: IUNPCK.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
*                                                                      *
************************************************************************
*
! The packed integer array jcupa declared as I*4 
! COMMON/iounit/ added
! XHH 1997.02.12
Cww      INTEGER PNTJQS,PNTRIQ
      POINTER (PNTJQS,JQSDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
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
      IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCF)) THEN
*cff     ISPAR = IUNPCK (JCUPA(1,ICSF),NNNW)
*        .. note this bit extraction does not preserve sign
         ISPAR = IBITS (JCUPA((NNNW-1)/4+1,ICSF), 8*MOD(NNNW-1,4), 8)
         IF (ISPAR .gt. 127) ISPAR = ISPAR -256
         ISPAR = SIGN (1,ISPAR)
      ELSE
         WRITE(istde,*) 'ISPAR: Argument ICSF is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
