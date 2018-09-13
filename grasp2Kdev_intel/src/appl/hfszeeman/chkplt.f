************************************************************************
*                                                                      *
      SUBROUTINE CHKPLT
*                                                                      *
*   This code checks for  consistent substitutions of plants between   *
*   HFSZEEMAN05 and the LIB92 subprograms.                             *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT, LODPLT.                               *
*                                                                      *
*   Written by Farid A Parpia               Last update: 24 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL LPLANT
      CHARACTER*256 RECORD
CPJ      INTEGER KEYORB
      include 'parameters.def'
CPJ      PARAMETER (KEYORB =121)
*
CPJ      INTEGER NNNP
CPJ      PARAMETER (NNNP = 590)
CPJ      INTEGER NNN1
CPJ      PARAMETER (NNN1 = 600)
CPJ      INTEGER NNNW
CPJ      PARAMETER (NNNW = 120)
CPJ      INTEGER NNNWP
CPJ      PARAMETER (NNNWP = 30)
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
            PRINT *, ' but to .FALSE. in HFSZEEMAN.'
         ELSE
            PRINT *, 'Plant DB was set to .FALSE. in LIB92,'
            PRINT *, ' but to .TRUE. in HFSZEEMAN.'
         ENDIF
         STOP
      ENDIF
*
*   Consistent numerical plants?
*
      PRINT *, 'NPLANTS: ', NPLANT(1), NPLANT(2), NPLANT(3), NPLANT(4)
      IF (NPLANT(1) .NE. KEYORB) THEN
         CALL CONVRT (NPLANT(1),RECORD,LENTH)
         PRINT *, 'CHKPLT: Plant KEYORB has been set to '
         PRINT *, ' '//RECORD(1:LENTH)//' in LIB92, but to KEYORB'
         PRINT *, ' in HFSZEEMAN.'
         STOP
      ENDIF
*
      IF (NPLANT(2) .NE. NNNP) THEN
         CALL CONVRT (NPLANT(2),RECORD,LENTH)
         PRINT *, 'CHKPLT: Plant NP has been set to '
         PRINT *, ' '//RECORD(1:LENTH)//' in LIB92, but to NNNP'
         PRINT *, ' in HFSZEEMAN.'
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
         PRINT *, ' in HFSZEEMAN.'
         STOP
      ENDIF
*
      IF (NPLANT(4) .NE. NNNWP) THEN
         CALL CONVRT (NPLANT(4),RECORD,LENTH)
         PRINT *, 'CHKPLT: Plant NWP has been set to '
         PRINT *, ' '//RECORD(1:LENTH)//' in LIB92, but to NNNWP'
         PRINT *, ' in HFSZEEMAN.'
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
