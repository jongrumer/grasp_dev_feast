************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***        *****   *******  **   **   *****    *****   **            ***
***       **   **  **   **  ***  **  **   **  **   **  **            ***
***       **       **       ** * **  **       **       **            ***
***       **       ****     ** * **  **        *****   **            ***
***       **  ***  **       **  ***  **            **  **            ***
***       **   **  **   **  **  ***  **   **  **   **  **   **       ***
***        *****   *******  **   **   *****    *****   *******       ***
***                                                                  ***
***    Package to generate a list of configuration symmetries for    ***
***    use by the remaining GRASP92 packages.                        ***
***                                                                  ***
***    This program is a derivative of JJCAS (F. A. Parpia, W. P.    ***
***    Wijesundera, and  I. P. Grant, Computer Physics Communica-    ***
***    tions, in preparation).                                       ***
***                                                                  ***
***                            GRASP92                               ***
***          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
***                                                                  ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      PROGRAM GENCSL
*                                                                      *
*   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
*               [GENCSL]: GENNRL, GENRL.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL GETYN,YES
      CHARACTER*256 FILNAM
      CHARACTER*11 DEFNAM
      CHARACTER*9 FORM
      CHARACTER*3 STATUS
      include 'parameters.def'
CGG      INTEGER NNNW, NNNWM1, NNNWM2
CGG      PARAMETER (NNNW = 120)
CGG      PARAMETER (NNNWM1 = 119)
CGG      PARAMETER (NNNWM2 = 118)

      EXTERNAL CONSTS
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
      COMMON/iounit/istdi,istdo,istde
*
*   Check critical plant substitutions
*
      IF ((NNNW-1) .NE. NNNWM1) THEN
         WRITE(istde,*) ' Values of plants NW and NWM1 are not '
     &,                 'consistent.'
         WRITE(istde,*) ' Reprocess and recompile package.'
         STOP
      ENDIF
*
      IF ((NNNW-2) .NE. NNNWM2) THEN
         WRITE(istde,*) ' Values of plants NW and NWM2 are not '
     &,                 'consistent.'
         WRITE(istde,*) ' Reprocess and recompile package.'
         STOP
      ENDIF
*
*   Set up the  .csl file; it is formatted; open it on unit 21
*
      DEFNAM = 'rcsl.out'
      FORM = 'FORMATTED'
      STATUS = 'NEW'
*
      FILNAM = DEFNAM
*
    1 CALL OPENFL (21,FILNAM,FORM,STATUS,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(istde,*) 'Error when opening rcsl.out'
         STOP
      ENDIF
*
*   Which convention?
*
      WRITE(istde,*) '''Nonrelativistic'' scheme?'
      YES = GETYN ()
      IF (YES) THEN
         CALL GENNRL
      ELSE
         CALL GENRL
      ENDIF
*
*   CLOSE  .csl
*
      CLOSE (21)
*
      STOP
      END
