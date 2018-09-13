************************************************************************
*                                                                      *
      FUNCTION IQR (ISUBSH,ICSF)
*                                                                      *
*   IQR is the occupation of subshell ISUBSH in CSF  ICSF.             *
*                                                                      *
*   Call(s) to: [LIB92]: IUNPCK.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
      include 'parameters.def'
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
      POINTER (PNTIQR,IQAR(NNNWP,*))
*
      COMMON/ORB2R/NCFR,NWR,PNTIQR
*
      IF ((ISUBSH .GE. 1) .AND. (ISUBSH .LE. NWR)) THEN
         IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCFR)) THEN
            IQR = IUNPCK (IQAR(1,ICSF),ISUBSH)
         ELSE
            PRINT *, 'IQR: Argument ICSF is out of range.'
            STOP
         ENDIF
      ELSE
         PRINT *, 'IQR: Argument ISUBSH is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
