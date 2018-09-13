************************************************************************
*                                                                      *
      SUBROUTINE CHKPLT
*                                                                      *
*   This code checks for  consistent substitutions of plants between   *
*   RIS3 and the LIB92 subprograms.                                     *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT, LODPLT.                               *
*                                                                      *
*   Written by Farid A Parpia               Last update: 24 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL LPLANT
      CHARACTER*256 RECORD
      include 'parameters.def'
CGG      INTEGER KEYORB
CGG      PARAMETER (KEYORB = 121)
*
CGG      INTEGER NNNP
CGG      PARAMETER (NNNP = 590)
CGG      INTEGER NNN1
CGG      PARAMETER (NNN1 = 600)
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
      COMMON/LIB92P/LPLANT,NPLANT(4)
*
*   Load COMMON/LIB92P/
*
      CALL LODPLT
*
*   Consistently DOUBLEPRECISION or REAL?
*
      IF (LPLANT .NEQV. .TRUE.) THEN
         IF (LPLANT) THEN
            PRINT *, 'Plant DB was set to .TRUE. in LIB92,'
            PRINT *, ' but to .FALSE. in RIS3.'
         ELSE
            PRINT *, 'Plant DB was set to .FALSE. in LIB92,'
            PRINT *, ' but to .TRUE. in RIS3.'
         ENDIF
         STOP
      ENDIF
*
*   Consistent numerical plants?
*
      IF (NPLANT(1) .NE. KEYORB) THEN
         CALL CONVRT (NPLANT(1),RECORD,LENTH)
         PRINT *, 'CHKPLT: Plant KEYORB has been set to '
         PRINT *, ' '//RECORD(1:LENTH)//' in LIB92, but to KEYORB'
         PRINT *, ' in RIS3.'
         STOP
      ENDIF
*
      IF (NPLANT(2) .NE. NNNP) THEN
         CALL CONVRT (NPLANT(2),RECORD,LENTH)
         PRINT *, 'CHKPLT: Plant NP has been set to '
         PRINT *, ' '//RECORD(1:LENTH)//' in LIB92, but to NNNP'
         PRINT *, ' in RIS3.'
         STOP
      ENDIF
*
      IF (NNN1 .NE. NNNP+10) THEN
         CALL CONVRT (NNNP+10,RECORD,LENTH)
         PRINT *, 'CHKPLT: Plant N1 should be set'
         PRINT *, ' to '//RECORD(1:LENTH)//'.'
         STOP
      ENDIF
*
      IF (NPLANT(3) .NE. NNNW) THEN
         CALL CONVRT (NPLANT(3),RECORD,LENTH)
         PRINT *, 'CHKPLT: Plant NW has been set to '
         PRINT *, ' '//RECORD(1:LENTH)//' in LIB92, but to NNNW'
         PRINT *, ' in RIS3.'
         STOP
      ENDIF
*
      IF (NPLANT(4) .NE. NNNWP) THEN
         CALL CONVRT (NPLANT(4),RECORD,LENTH)
         PRINT *, 'CHKPLT: Plant NWP has been set to '
         PRINT *, ' '//RECORD(1:LENTH)//' in LIB92, but to NNNWP'
         PRINT *, ' in RIS3.'
         STOP
      ENDIF
*
      IF (MOD (NNNW,4) .EQ. 0) THEN
         NWP = NNNW/4
      ELSE
         NWP = NNNW/4+1
      ENDIF
*
      IF (NNNWP .NE. NWP) THEN
         CALL CONVRT (NWP,RECORD,LENTH)
         PRINT *, 'CHKPLT: Plant NWP should be set'
         PRINT *, ' to '//RECORD(1:LENTH)//'.'
         STOP
      ENDIF
*
      RETURN
      END
