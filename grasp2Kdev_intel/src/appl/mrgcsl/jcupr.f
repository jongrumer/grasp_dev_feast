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
! Packed array declared I*4, common/iounit/ added
! XHH 1997.02.12
Cww      INTEGER PNJQSR,PNTIQR
      POINTER (PNJQSR,JQSRDUMMY)
      POINTER (PNTIQR,IQRDUMMY)
*
      include 'parameters.def'
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
      INTEGER*4 JCUPAR
      POINTER (PJCUPR,JCUPAR(NNNWP,*))
*
      COMMON/ORB2R/NCFR,NWR,PNTIQR
     :      /STATR/PNJQSR,PJCUPR
      COMMON/iounit/istdi,istdo,istde
*
      IF ((LOC .GE. 1) .AND. (LOC .LE. NWR-1)) THEN
         IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCFR)) THEN
            JCUPR = IUNPCK (JCUPAR(1,ICSF),LOC)
         ELSE
            WRITE(istde,*) 'JCUPR: Argument ICSF is out of range.'
            STOP
         ENDIF
      ELSE
         WRITE(istde,*) 'JCUPR: Argument LOC is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
