************************************************************************
*                                                                      *
      FUNCTION JCUPR (LOC,ICSF)
*                                                                      *
*   JCUPR is the 2J+1 value of the LOCth nontrivial intermediate ang-  *
*   ular momentum in CSF  ICSF.                                        *
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
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
      POINTER (PJCUPR,JCUPAR(NNNWP,*))
*
      COMMON/ORB2R/NCFR,NWR,PNTIQR
     :      /STATR/PNJQSR,PJCUPR
*
      IF ((LOC .GE. 1) .AND. (LOC .LE. NWR-1)) THEN
         IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCFR)) THEN
            JCUPR = IUNPCK (JCUPAR(1,ICSF),LOC)
         ELSE
            PRINT *, 'JCUPR: Argument ICSF is out of range.'
            STOP
         ENDIF
      ELSE
         PRINT *, 'JCUPR: Argument LOC is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
