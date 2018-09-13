************************************************************************
*                                                                      *
      SUBROUTINE SETPOT (J,JP)
*                                                                      *
*   This  subroutine  sets  up the  arrays TF and TG for use by  the   *
*   subprograms IN, OUT, and SBSTEP.                                   *
*                                                                      *
*   Arguments:                                                         *
*                                                                      *
*      J:  (Input) Index of orbital                                    *
*      JP: (Output) Join point: point where TG changes sign            *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 16 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      COMMON/iounit/istdi,istdo,istde

      LOGICAL NOTSET
*
      COMMON/DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INT3/TF(NNNP),TG(NNNP),XU(NNNP),XV(NNNP)
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /POTE/YP(NNNP),XP(NNNP),XQ(NNNP)
*
      NOTSET = .TRUE.
*
*   Define constants
*
      DMHB2C = -H/(2.0D 00*C)
      ENERGY = E(J)
      ENEFAC = 2.0D 00*C*C-ENERGY
*
*   Set up arrays TF and TG
*
*   Since TF(1) and TG(1) are never used, set them
*   to some arbitrary value
*
      JP = 0
      TF(1) = 0.0D 00
      TG(1) = 0.0D 00
      DO 1 I = 2,N
         RPI = RP(I)
         YPRPOR = YP(I)*RPOR(I)
         TF(I) = DMHB2C*(ENEFAC*RPI+YPRPOR)
         TG(I) = DMHB2C*(ENERGY*RPI-YPRPOR)
         IF (NOTSET) THEN
            IF (ABS ( SIGN (1.0D 00,TG(I))
     :               +SIGN (1.0D 00,TG(1))) .LT. ACCY) THEN
               JP = I
               NOTSET = .FALSE.
            ENDIF
         ENDIF
    1 CONTINUE
*
*   Trap for inappropriate grid
*
      IF (JP .EQ. 0) THEN
         WRITE(istde,*) 'SETPOT: Grid of insufficient extent.'
         STOP
      ENDIF
*
      RETURN
      END
