************************************************************************
*                                                                      *
      SUBROUTINE CHKPLT
*                                                                      *
*   This code checks for  consistent substitutions of plants between   *
*   ERWF and the LIB92 subprograms.                                    *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT, LODPLT.                               *
*                                                                      *
*   Written by Farid A Parpia               Last update: 15 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL LPLANT
      CHARACTER*256 RECORD
*
      include 'parameters.def'
CGG      INTEGER NNNP
CGG      PARAMETER (NNNP = 590)
CGG      INTEGER NNN1
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
      COMMON/LIB92P/LPLANT,NPLANT(4)

      COMMON/iounit/istdi,istdo,istde
*
*   Load COMMON/LIB92P/
*
      CALL LODPLT
*
*   Consistently REAL*8          or REAL?
*
      IF (LPLANT .NEQV. .TRUE.) THEN
         IF (LPLANT) THEN
            WRITE(istde,*) 'Plant DB was set to .TRUE. in LIB92,'
     &,                    ' but to .FALSE. in ERWF.'
         ELSE
            WRITE(istde,*) 'Plant DB was set to .FALSE. in LIB92,'
     &,                    ' but to .TRUE. in ERWF'
         ENDIF
         STOP
      ENDIF
*
*   Consistent numerical plants?
*
      IF (NPLANT(2) .NE. NNNP) THEN
         CALL CONVRT (NPLANT(2),RECORD,LENTH)
         WRITE(istde,*) 'CHKPLT: Plant NP has been set to '
     &                  //RECORD(1:LENTH)//' in LIB92, but to ', NNNP
     &,                 ' in ERWF.'
         STOP
      ENDIF
*
      IF (NNN1 .NE. NNNP+10) THEN
         CALL CONVRT (NNNP+10,RECORD,LENTH)
         WRITE(istde,*) 'CHKPLT: Plant N1 should be set'
     &,                 ' to '//RECORD(1:LENTH)//'.'
         STOP
      ENDIF
*
      IF (NPLANT(3) .NE. NNNW) THEN
         CALL CONVRT (NPLANT(3),RECORD,LENTH)
         WRITE(istde,*) 'CHKPLT: Plant NW has been set to '
     &                  //RECORD(1:LENTH)//' in LIB92, but to ',NNNW
     &,                 ' in ERWF.'
         STOP
      ENDIF
*
      IF (NPLANT(4) .NE. NNNWP) THEN
         CALL CONVRT (NPLANT(4),RECORD,LENTH)
         WRITE(istde,*) 'CHKPLT: Plant NWP has been set to '
     &                  //RECORD(1:LENTH)//' in LIB92, but to ',NNNWP
     &,                 ' in ERWF.'
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
         WRITE(istde,*) 'CHKPLT: Plant NWP should be set'
     &,                 ' to '//RECORD(1:LENTH)//'.'
         STOP
      ENDIF
*
      RETURN
      END
