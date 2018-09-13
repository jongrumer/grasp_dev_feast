************************************************************************
*                                                                      *
      SUBROUTINE CHKPLT (from)
*                                                                      *
*   This code checks for  consistent substitutions of plants between   *
*   the application (from) and the LIB92 subprograms.                  *
*                                                                      *
*   Call(s) to: [LIB92]: LODPLT.                                       *
*                                                                      *
*   Written by Farid A Parpia               Last update: 09 Dec 1992   *
*   Modified by Xinghong He                 Last update: 15 Jul 1998   *
*                                                                      *
************************************************************************
      CHARACTER*(*) from
      CHARACTER*(*), PARAMETER:: myname = 'CHKPLT'
      include 'parameters.def'
CGG      INTEGER KEYORB, NNNP, NNN1, NNNW, NNNWP
CGG      PARAMETER (KEYORB = 121)
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
CGG      PARAMETER (NNNWP = 30)
CGG      PARAMETER (NNNWP = 54)

      LOGICAL LPLANT
      COMMON/LIB92P/LPLANT,NPLANT(4)

      COMMON/iounit/istdi,istdo,istde
!-----------------------------------------------------------------------
      lenfrom = LEN_TRIM (from)
*
*   Load COMMON/LIB92P/
*
      CALL LODPLT
*
*   Consistent numerical plants?
*
      IF (NPLANT(1) .NE. KEYORB) THEN
         WRITE (istde,*) myname,': Plant KEYORB has been set to ',
     &                   NPLANT(1), ' in LIB92, but to ',
     &                   KEYORB, ' in ', from(1:lenfrom)
         STOP
      ENDIF

      IF (NPLANT(2) .NE. NNNP) THEN
         WRITE (istde,*) myname,': Plant NP has been set to ',
     &                   NPLANT(2), ' in LIB92, but to ',
     &                   NNNP, ' in ', from(1:lenfrom)
         STOP
      ENDIF

      IF (NNN1 .NE. NNNP+10) THEN
         WRITE (istde,*) myname,': Plant N1 should be set to ', NNNP+10
         STOP
      ENDIF

      IF (NPLANT(3) .NE. NNNW) THEN
         WRITE (istde,*) myname,': Plant NW has been set to ',
     &                   NPLANT(3),' in LIB92, but to ',
     &                   NNNW, ' in ', from(1:lenfrom)
         STOP
      ENDIF

      IF (NPLANT(4) .NE. NNNWP) THEN
         WRITE (istde,*) myname,': Plant NWP has been set to ',
     &                   NPLANT(4),' in LIB92, but to ',
     &                   NNNWP, ' in ', from(1:lenfrom)
         STOP
      ENDIF

      IF (MOD (NNNW,4) .EQ. 0) THEN
         NWP = NNNW/4
      ELSE
         NWP = NNNW/4+1
      ENDIF

      IF (NNNWP .NE. NWP) THEN
         WRITE (istde,*) myname,': Plant NWP should be set to', NWP
         STOP
      ENDIF

      RETURN
      END
