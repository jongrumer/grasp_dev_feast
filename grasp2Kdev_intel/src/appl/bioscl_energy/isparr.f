************************************************************************
*                                                                      *
      FUNCTION ISPARR (ICSF)
*                                                                      *
*   ISPARR is the value of P for CSF number ICSF.                      *
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
*        ISPARR = IUNPCK (JCUPAR(1,ICSF),NNNW)
*        .. note this bit extraction does not preserve sign
         ISPARR = IBITS (JCUPAR((NNNW-1)/4+1,ICSF), 8*MOD(NNNW-1,4), 8)
         IF (ISPARR .gt. 127) ISPARR = ISPARR -256
         ISPARR = SIGN (1,ISPARR)
      ELSE
         PRINT *, 'ISPARR: Argument ICSF is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
