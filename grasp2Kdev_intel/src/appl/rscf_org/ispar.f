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
*   The only (real) difference from ispar.f in lib92/ is the removal
*   of the check for ICSF.
*
      POINTER (PNTJQS,JQSDUMMY)
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
      INTEGER*4 JCUPA
      POINTER (PNJCUP,JCUPA(NNNWP,1))
*
      COMMON/STAT/PNTJQS,PNJCUP
      COMMON/iounit/istdi,istdo,istde
*        .. note this bit extraction does not preserve sign
      ISPAR = IBITS (JCUPA((NNNW-1)/4+1,ICSF), 8*MOD(NNNW-1,4), 8)
      IF (ISPAR .gt. 127) ISPAR = ISPAR -256
      ISPAR = SIGN (1,ISPAR)
*
      RETURN
      END
