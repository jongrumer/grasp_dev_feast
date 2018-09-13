************************************************************************
*                                                                      *
      FUNCTION ITJPO (ICSF)
*                                                                      *
*   ITJPO is the value of 2J+1 for CSF number ICSF.                    *
*                                                                      *
*   Call(s) to: [LIB92]: IUNPCK.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
*                                                                      *
************************************************************************
*   The only (real) difference from itjpo.f in lib92/ is the removal
*   of the check for ICSF.
*
      POINTER (PNTJQS,JQSDUMMY)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)

      INTEGER*4 JCUPA
      POINTER (PNJCUP,JCUPA(NNNWP,1))

      COMMON/STAT/PNTJQS,PNJCUP
      COMMON/iounit/istdi,istdo,istde
*        bit extraction will not give the proper sign.
      ITJPO = IBITS (JCUPA((NNNW-1)/4+1,ICSF),8*MOD(NNNW-1,4),8)
      IF (ITJPO .gt. 127) ITJPO = 256-ITJPO
*
      RETURN
      END
