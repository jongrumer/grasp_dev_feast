************************************************************************
*                                                                      *
      FUNCTION IQ (ISUBSH,ICSF)
*                                                                      *
*   IQ is the occupation of subshell ISUBSH in CSF  ICSF.              *
*                                                                      *
*   Call(s) to: [LIB92]: IUNPCK.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 30 Oct 1992   *
*                                                                      *
************************************************************************
*
! IQA declared I*4 and COMMON/iounit/ added
! XHH 1997.02.12
      include 'parameters.def'
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
! XHH added
      COMMON/iounit/istdi,istdo,istde
      INTEGER*4 IQA

      POINTER (PNTRIQ,IQA(NNNWP,*))
*
      COMMON/ORB2/NCF,NW,PNTRIQ
*
c     IF ((ISUBSH .GE. 1) .AND. (ISUBSH .LE. NW)) THEN
c        IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCF)) THEN
*cff        IQ = IUNPCK (IQA(1,ICSF),ISUBSH)
            IQ = IBITS (IQA((ISUBSH-1)/4+1,ICSF),8*MOD(ISUBSH-1,4),8)
c        ELSE
c           WRITE(istde,*) 'IQ: Argument ICSF is out of range.'
c           STOP
c        ENDIF
c     ELSE
c        WRITE(istde,*) 'IQ: Argument ISUBSH is out of range.'
c        WRITE(istde,*) 'ISUBSH=',ISUBSH, '  NW=',NW
c        STOP
c     ENDIF
*
      RETURN
      END
