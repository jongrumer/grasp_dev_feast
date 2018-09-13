************************************************************************
*                                                                      *
      SUBROUTINE COUNT (FR,MTPFR,NNCFF,SGN)
*                                                                      *
*   This subroutine counts the nodes in the radial function FR using   *
*   the criteria  given by C Froese Fischer, Comp Phys Rep, 3 (1986)   *
*   314-315 . The  function FR is assumed defined on the first MTPFR   *
*   points of the radial grid. The sign of the function at the first   *
*   oscillation is also determined.                                    *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
*
      DIMENSION FR(NNNP)
      DIMENSION LCEXT(NNNP)
*
      COMMON/COUN/THRESH
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
*
*   (1) Find all extrema in FR
*   (2) Find the maximum amplitudes of FR
*
      NEXT = 1
      EXT = 0.0D 00
      LCEXT(1) = 1
      EMX = 0.0D 00
      DO 1 I = 2,MTPFR
         ABFRI = ABS (FR(I))
         TEST = ABS ( SIGN (1.0D 00,FR(I))
     :               +SIGN (1.0D 00,FR(I-1)))
         IF (TEST .LE. ACCY) THEN
            NEXT = NEXT+1
            LCEXT(NEXT) = 0
            EXT = 0.0D 00
         ENDIF
         IF (ABFRI .GT. EXT) THEN
            EXT = ABFRI
            LCEXT(NEXT) = I
         ENDIF
         IF (ABFRI .GT. EMX) THEN
            EMX = ABFRI
         ENDIF
    1 CONTINUE
*
*   Eliminate oscillations with amplitude less than THRESH times
*   the maximum
*
      LOC = 0
      THRESE = THRESH*EMX
    4 CONTINUE
      LOC = LOC+1
      IF (LOC .LE. NEXT) THEN
         IF (LCEXT(LOC) .EQ. 0) THEN
            ABLCL = 0.0D 00
         ELSE
            ABLCL = ABS (FR(LCEXT(LOC)))
         ENDIF
         IF (ABLCL .LT. THRESE) THEN
            NEXT = NEXT-1
            NSTPS = NEXT-LOC
            DO 5 I = 0,NSTPS
               LCEXT(LOC+I) = LCEXT(LOC+I+1)
    5       CONTINUE
            LOC = LOC-1
         ENDIF
         GOTO 4
      ENDIF
*
*   Count changes of sign using the remaining oscillations
*
      NNCFF = 0
      DO 6 I = 2,NEXT
         TEST = ABS (SIGN (1.0D 00,FR(LCEXT(I  ))) +
     :               SIGN (1.0D 00,FR(LCEXT(I-1))) )
         IF (TEST .LE. ACCY) THEN
            NNCFF = NNCFF+1
         ENDIF
    6 CONTINUE
*
*   Determine the position of the first oscillation, and the
*   sign of the function at this location
*
      SGN = SIGN (1.0D 00,FR(LCEXT(1)))
*
      RETURN
      END
