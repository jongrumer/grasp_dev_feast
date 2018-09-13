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
! Packed array declared I*4, common/iounit/ added
! XHH 1997.02.12
      include 'parameters.def'
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
      INTEGER*4 IQAR
      POINTER (PNTIQR,IQAR(NNNWP,*))
*
      COMMON/ORB2R/NCFR,NWR,PNTIQR
      COMMON/iounit/istdi,istdo,istde
*
      IF ((ISUBSH .GE. 1) .AND. (ISUBSH .LE. NWR)) THEN
         IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCFR)) THEN
            IQR = IUNPCK (IQAR(1,ICSF),ISUBSH)
         ELSE
            WRITE(istde,*) 'IQR: Argument ICSF is out of range.'
            STOP
         ENDIF
      ELSE
         WRITE(istde,*) 'IQR: Argument ISUBSH is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
