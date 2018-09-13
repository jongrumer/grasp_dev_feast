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
! Packed array declared I*4, common/iounit/ added
! XHH 1997.02.12
Cww      INTEGER PNJQSR,PNTIQR
      POINTER (PNJQSR,JQSRDUMMY)
      POINTER (PNTIQR,IQRDUMMY)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
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
      IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCFR)) THEN
         ISPARR = IUNPCK (JCUPAR(1,ICSF),NNNW)
         ISPARR = SIGN (1,ISPARR)
      ELSE
         WRITE(istde,*) 'ISPARR: Argument ICSF is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
