************************************************************************
*                                                                      *
      FUNCTION ITJPOR (ICSF)
*                                                                      *
*   ITJPOR is the value of 2J+1 for CSF number ICSF.                   *
*                                                                      *
*   Call(s) to: [LIB92]: IUNPCK.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PNJQSR,PNTIQR
      POINTER (PNJQSR,JQSRDUMMY)
      POINTER (PNTIQR,IQRDUMMY)
*
      include 'parameters.def'
Cff      INTEGER NNNW
Cff      PARAMETER (NNNW = 120)
Cff      INTEGER NNNWP
Cff      PARAMETER (NNNWP = 30)
*
      INTEGER*4 JCUPAR
      POINTER (PJCUPR,JCUPAR(NNNWP,1))
*
      COMMON/ORB2R/NCFR,NWR,PNTIQR
     :      /STATR/PNJQSR,PJCUPR
*
      IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCFR)) THEN
*cff        ITJPOR = IUNPCK (JCUPAR(1,ICSF),NNNW)
*        bit extraction will not give the proper sign.
         ITJPOR = IBITS (JCUPAR((NNNW-1)/4+1,ICSF),8*MOD(NNNW-1,4),8)
         IF (ITJPOR .gt. 127) ITJPOR = 256-ITJPOR
      ELSE
         PRINT *, 'ITJPOR: Argument ICSF is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
