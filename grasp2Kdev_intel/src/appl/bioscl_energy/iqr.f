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
Cff      INTEGER NNNWP

      include 'parameters.def'
Cff      PARAMETER (NNNWP = 30)
*
      INTEGER*4 IQAR

      POINTER (PNTIQR,IQAR(NNNWP,1))
*
      COMMON/ORB2R/NCFR,NWR,PNTIQR
*
      IF ((ISUBSH .GE. 1) .AND. (ISUBSH .LE. NWR)) THEN
         IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCFR)) THEN
*cff        IQR = IUNPCK (IQAR(1,ICSF),ISUBSH)
            IQR = IBITS (IQAR((ISUBSH-1)/4+1,ICSF),8*MOD(ISUBSH-1,4),8)
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
