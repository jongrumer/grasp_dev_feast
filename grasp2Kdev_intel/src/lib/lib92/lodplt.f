************************************************************************
*                                                                      *
      SUBROUTINE LODPLT
*                                                                      *
      include 'parameters.def'
CGG      INTEGER NNNP
CGG      PARAMETER (NNNP = 590)
CGG      INTEGER NNN1
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
CGG      INTEGER KEYORB
CGG      PARAMETER (KEYORB = 121)
*   This subroutine loads   COMMON/LIB92P/    with the values of the   *
*   plants substituted in the LIB92 subprograms.                       *
*                                                                      *
*   Written by Farid A Parpia               Last update: 09 Dec 1992   *
*                                                                      *
************************************************************************
      LOGICAL LPLANT
*
      COMMON/LIB92P/LPLANT,NPLANT(4)
*
*   Load COMMON/LIB92P/
*
      LPLANT = .TRUE.
*
      NPLANT(1) = KEYORB
      NPLANT(2) = NNNP
*
*   Consistency check
*
      IF (NNN1 .NE. NNNP+10) THEN
CGG         PRINT *, 'LODPLT: Plant N1 should be equal to ',NP+10
         PRINT *, 'LODPLT: Plant N1 should be equal to ',NNNP+10
     &, ' in the LIB92 subprograms.'
         STOP
      ENDIF
*
      NPLANT(3) = NNNW
      NPLANT(4) = NNNWP
*
      RETURN
      END
