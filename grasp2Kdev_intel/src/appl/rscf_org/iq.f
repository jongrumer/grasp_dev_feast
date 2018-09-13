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
*   The only (real) difference from iq.f in lib92/ is the removal
*   of the check for ICSF.
*
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
      IF ((ISUBSH .GE. 1) .AND. (ISUBSH .LE. NW)) THEN
         IQ = IBITS (IQA((ISUBSH-1)/4+1,ICSF),8*MOD(ISUBSH-1,4),8)
      ELSE
         WRITE(istde,*) 'IQ: Argument ISUBSH is out of range.'
         WRITE(istde,*) 'ISUBSH=',ISUBSH, '  NW=',NW
         STOP
      ENDIF
*
      RETURN
      END
