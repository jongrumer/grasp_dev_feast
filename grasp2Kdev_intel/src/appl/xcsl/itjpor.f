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
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
      POINTER (PJCUPR,JCUPAR(NNNWP,*))
*
      COMMON/ORB2R/NCFR,NWR,PNTIQR
     :      /STATR/PNJQSR,PJCUPR
*
      IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCFR)) THEN
         ITJPOR = IUNPCK (JCUPAR(1,ICSF),NNNW)
         ITJPOR = ABS (ITJPOR)
      ELSE
         PRINT *, 'ITJPOR: Argument ICSF is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
